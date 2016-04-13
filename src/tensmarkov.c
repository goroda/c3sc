#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "tensmarkov.h"
#include "util.h"
#include "simulate.h"
#include "dynamics.h"
#include "cost.h"


/** \struct MCNode
 *  \brief Markov Chain Node (tree like structure)
 *  \var MCNode::d
 *  dimension of node value
 *  \var MCNode::x
 *  node value
 *  \var MCNode::N
 *  number of neighbors
 *  \var MCNode::neighbors
 *  neighboring nodes
 *  \var MCNode::pself
 *  transition probability to oneself
 *  \var MCNode::p
 *  transition probabilities to neighbors
 *  \var MCNode::gpself
 *  Gradient associated with transitionining to oneself
 *  \var MCNode::gp
 *  gradients associated with transitioning to neighbors
 */
struct MCNode
{
    size_t d;
    size_t du; // dimension of gradients
    int xalloc;
    double * x;

    
    //size_t N;
    //struct MCNode ** neighbors;
    struct MCNList * neigh;
    double pself;
    int galloc;
    double * gpself;
//    double * p; // transition probabilities
//    double ** gp;

};

/**********************************************************//**
    Allocate node

    \param[in] d - dimension of node

    \return allocate node
**************************************************************/
struct MCNode * mcnode_alloc(size_t d)
{
    struct MCNode * mcn;
    mcn = malloc(sizeof(struct MCNode));
    assert (mcn != NULL);
    
    mcn->d = d;
    mcn->du = 0;
    mcn->xalloc = 0;
    mcn->x = NULL;
    mcn->pself = 1.0;
    mcn->galloc = 0;
    mcn->gpself = NULL;
    mcn->neigh = NULL;
    
    return mcn;
}

/**********************************************************//**
    Free a node
**************************************************************/
void mcnode_free(struct MCNode * mcn)
{
    if (mcn != NULL){
        if (mcn->xalloc == 1){
            free(mcn->x); mcn->x = NULL;
        }
        if (mcn->galloc == 1){
            free(mcn->gpself);
        }
        mcnlist_free(mcn->neigh); mcn->neigh = NULL;
        free(mcn); mcn = NULL;
        
        /* //        free_dd(mcn->N,mcn->gp); mcn->gp = NULL; */
        /* if (mcn->N > 0){ */
        /*     mcnode_free_array(mcn->neighbors,mcn->N); */
        /*     free(mcn->neighbors); mcn->neighbors = NULL; */
        /* } */
    }
}

/**********************************************************//**
    Initialize a node by referencing an array x of size d, 
    node has no neighbors and only transitions to itself
**************************************************************/
struct MCNode * mcnode_init(size_t d, double * x)
{
    struct MCNode * mm = mcnode_alloc(d);
    assert (mm != NULL);
    mm->x = x;
    mm->pself = 1.0;
    return mm;
}

/**********************************************************//**
    Set pself;
**************************************************************/
void mcnode_set_pself(struct MCNode * mcn,double pself)
{
    assert (mcn != NULL);
    mcn->pself = pself;
}

/**********************************************************//**
    Add gradient for mcnode
**************************************************************/
void mcnode_set_gradient(struct MCNode * mcn,size_t du,double * grad)
{
    assert (mcn != NULL);
    mcn->du = du;
    mcn->gpself = grad;
}

/**********************************************************//**
    Initialize gradient of prob of transitioning to oneself
**************************************************************/
void mcnode_init_self_grad(struct MCNode * mcn, size_t du)
{
    mcn->du = du;
    mcn->galloc = 1;
    mcn->gpself = calloc_double(du);
}

/**********************************************************//**
    Add Neighbor
**************************************************************/
void mcnode_prepend_neigh(struct MCNode * mcn, size_t dir,
                          double val, double p, double * grad)
{
    mcnlist_prepend(&(mcn->neigh),dir,val,p,grad);
}

size_t mcnode_get_d(const struct MCNode * mcn){ return mcn->d; }
double * mcnode_get_xref(const struct MCNode * mcn){ return mcn->x;}
double mcnode_get_pself(const struct MCNode * mcn){ return mcn->pself; }
size_t mcnode_get_du(const struct MCNode * mcn){ return mcn->du;}
double * mcnode_get_gpself(const struct MCNode * mcn){ return mcn->gpself;}
struct MCNList * mcnode_get_neigh(const struct MCNode * mcn){ return mcn->neigh;}
size_t mcnode_get_n(const struct MCNode * mcn){ return mcnlist_length(mcn->neigh);}

//////////////////////////////////

struct MCNList
{
    double p; // probability to transition to this node
    int galloc;
    double * gradp; // gradient of transition probability
    size_t dir; // direction of neighbor
    double val;
    struct MCNList * next;
};

double mcnlist_get_p(const struct MCNList * mcn){return mcn->p;}
double * mcnlist_get_gradp(const struct MCNList * mcn){return mcn->gradp;}
size_t mcnlist_get_dir(const struct MCNList * mcn){return mcn->dir;}
double mcnlist_get_val(const struct MCNList * mcn){return mcn->val;}
struct MCNList * mcnlist_get_next(const struct MCNList * mcn){return mcn->next;}
size_t mcnlist_length(const struct MCNList * mcn)
{
    size_t N = 0;
    const struct MCNList * temp = mcn;
    while (temp != NULL){
        N++;
        temp = temp->next;
    }
    return N;
}


/**********************************************************//**
    Allocate memory for a list

    \return allocated list
**************************************************************/
struct MCNList * mcnlist_alloc()
{
    struct MCNList * mcnl = malloc(sizeof(struct MCNList));
    if (mcnl == NULL){
        fprintf(stderr,"Error allocating memory for MCNList\n");
        exit(1);
    }
    mcnl->p = 0;
    mcnl->galloc = 0;
    mcnl->gradp = NULL;
    mcnl->next = NULL;
    return mcnl;
}

/**********************************************************//**
    Free a list
**************************************************************/
void mcnlist_free(struct MCNList * mcnl)
{
    if (mcnl != NULL){
        if (mcnl->galloc == 1){
            free(mcnl->gradp); mcnl->gradp = NULL;
        }
        //mcnode_free(mcnl->node); mcnl->node = NULL;
        mcnlist_free(mcnl->next); mcnl->next = NULL;
        free(mcnl); mcnl = NULL;
    }
}

/**********************************************************//**
    Add a neighbor to the front of the list
**************************************************************/
struct MCNList * mcnlist_prepend(struct MCNList ** mcn, size_t dir,
                                 double val, double p, double * grad)
{
    struct MCNList * node = mcnlist_alloc();
    node->p = p;
    node->gradp = grad;
    node->dir = dir;
    node->val = val;
    node->next = *mcn;
    *mcn = node;
    return node;
}

/**********************************************************//**
    Prepend an empty element
**************************************************************/
struct MCNList * mcnlist_prepend_empty(struct MCNList ** mcn, size_t dir,
                                       double val, size_t dgrad)
{
    struct MCNList * node = mcnlist_alloc();
    node->dir = dir;
    node->val = val;
//    node->p = p;
    if (dgrad > 0){
        node->galloc = 1;
        node->gradp = calloc_double(dgrad);
    }
    node->next = *mcn;
    *mcn = node;
    return node;
}

/**********************************************************//**
    Compute the expectation of a function around a node

    \param[in]     mc   - node
    \param[in]     f    - function(nvals,times,states,evas,arg)
    \param[in]     arg  - additional function arguments
    \param[in,out] grad - If (!NULL) return gradient

    \return average
**************************************************************/
double mcnode_expectation(
    const struct MCNode * mc,
    //void (*f)(size_t,double*,double**x,double*,void*),
    double (*f)(double,double*,void*),
    void * arg,double * grad)
{

    if (grad!= NULL){
        for (size_t ii = 0; ii < mc->du; ii++){
            grad[ii] = 0.0;
        }
    }

    double * point2 = calloc_double(mc->d);
    memmove(point2,mc->x,mc->d * sizeof(double));

    double avg = 0.0;
    double time = 0.0;
    double eval;
    // take care of self transition
    if (mc->pself > 0.0){
        eval = f(time,point2,arg);
        avg += mc->pself * eval;
        if (grad != NULL){
            cblas_daxpy(mc->du,eval,mc->gpself,1,grad,1);
        }
    }
    
    // take care of transitions to neighbors
    if (mc->neigh != NULL){
        struct MCNList * temp = mc->neigh;
        while (temp != NULL){
            point2[temp->dir] = temp->val;
            //dprint(mc->d,point2);
            eval = f(time,point2,arg);
            avg += temp->p * eval;
            if (grad != NULL){
                cblas_daxpy(mc->du,eval,temp->gradp,1,grad,1);  
            } 
            point2[temp->dir] = mc->x[temp->dir];
            temp = temp->next;
        }
    }
    free(point2); point2 = NULL;
    return avg;
}

/**********************************************************//**
    Print a node

    \param[in] n          - state
    \param[in] fp         - stream to print to
    \param[in] prec       - precision to print
**************************************************************/
void mcnode_print(struct MCNode * n, FILE * fp, int prec)
{
    if (n == NULL){
        fprintf(fp,"NULL");
        return;
    }
    else if (n->x == NULL){
        fprintf(fp,"NULL");
        return;
    }

    fprintf(fp,"Node is ");
    for (size_t ii = 0; ii < n->d; ii++){
        fprintf(fp,"%3.*G ",prec,n->x[ii]);
    }
    fprintf(fp,"\n");

    size_t N = 0;
    if (n->neigh != NULL){
        struct MCNList * temp = n->neigh;
        while (temp != NULL){
            N++;
            temp = temp->next;
        }
    }
    fprintf(fp,"Number of neighbors are %zu\n",N);
    /* for (size_t ii = 0; ii < n->N; ii++){ */
    /*     fprintf(fp,"\t"); */
    /*     mcnode_print(n->neighbors[ii],fp,prec); */
    /* } */
}

/** \struct MCA
 *  \brief Markov Chain Approximation struct
 *  \var MCA::d
 *  dimension of state space
 *  \var MCA::dw
 *  dimension of noise
 *  \var MCA::h
 *  distance between neighbors in each direction
 *  \var MCA::minh
 *  smallest distance in dimensions
 *  \var MCA::dyn
 *  dynamics
 *  \var MCA::bound
 *  boundary
 */
struct MCA
{
    size_t d;
    size_t du;
    size_t dw; 
    double * h;
    double minh;

    struct Dyn * dyn;
    struct Boundary * bound;

    double * drift;
    double * gdrift;
    double * diff;
    double * gdiff;
    
};

/**********************************************************//**
    Allocate MCA

    \param[in] d  - dimension of state space
    \param[in] du - dimension of control
    \param[in] dw - dimension of noise
    \param[in] h  - spacing of nodes in each direction

    \return mca
**************************************************************/
struct MCA * mca_alloc(size_t d, size_t du, size_t dw, double * h)
{

    struct MCA * mca = malloc(sizeof(struct MCA));
    assert (mca != NULL);
    
    mca->d = d;
    mca->du = du;
    mca->dw =dw;
    mca->h = calloc_double(d);
    mca->minh = h[0];
    for (size_t ii = 0; ii < d; ii++){
        mca->h[ii] = h[ii];
        if (h[ii] < mca->minh){
            mca->minh = h[ii];
        }
    }
    mca->dyn = NULL;
    mca->bound = NULL;

    mca->drift = calloc_double(d);
    mca->gdrift = calloc_double(d * du);
    mca->diff = calloc_double(d*dw);
    mca->gdiff = calloc_double(d * dw * du);
    
    return mca;
}

/**********************************************************//**
    Free memory allocated to an MCA
**************************************************************/
void mca_free(struct MCA * mca)
{
    if (mca != NULL){
        free(mca->h); mca->h = NULL;
        free(mca->drift); mca->drift = NULL;
        free(mca->diff); mca->diff = NULL;

        free(mca->gdrift); mca->gdrift = NULL;
        free(mca->gdiff); mca->gdiff = NULL;

        free(mca); mca = NULL;
    }
}

/**********************************************************//**
    Free memory allocated to an MCA along with boundary 
    and dynamics
**************************************************************/
void mca_free_deep(struct MCA * mca)
{
    if (mca != NULL){
        boundary_free(mca->bound); mca->bound = NULL;
        dyn_free_deep(mca->dyn); mca->dyn = NULL;
        mca_free(mca); mca = NULL;
        free(mca); mca = NULL;
    }
}

/**********************************************************//**
    Add a reference to dynamics to the MCA 
**************************************************************/
void mca_attach_dyn(struct MCA * mca, struct Dyn * dyn)
{
    assert (mca != NULL);
    mca->dyn = dyn;
}

/**********************************************************//**
    Get a reference to dynamics to the MCA 
**************************************************************/
struct Dyn * mca_get_dyn(struct MCA * mca)
{
    assert (mca != NULL);
    return mca->dyn;
}

/**********************************************************//**
    Add a reference to Boundary to the MCA 
**************************************************************/
void mca_attach_bound(struct MCA * mca,struct Boundary * bound)
{
    assert (mca != NULL);
    mca->bound = bound;
}

/**********************************************************//**
    Get size of state
**************************************************************/
size_t mca_get_dx(struct MCA * mca)
{
    return mca->d;
}

/**********************************************************//**
    Get size of control
**************************************************************/
size_t mca_get_du(struct MCA * mca)
{
    return mca->du;
}

/**********************************************************//**
    Generate a node that is outside the boundays
    
    \param[in] mca  - markov chain approximation 
    \param[in] time - time (not used)
    \param[in] x    - current state
    
    \return node
**************************************************************/
struct MCNode *
mca_outbound_node(struct MCA * mca, double time, double * x)
{
    (void)(time);
    struct MCNode * mcn = mcnode_init(mca->d,x);
    size_t du = mca_get_du(mca);
    mcnode_init_self_grad(mcn,du);
    return mcn;
}

/**********************************************************//**
    Generate a node that is inside the boundarys and its 
    transitions
    
    \param[in]     mca  - markov chain approximation 
    \param[in]     time - time
    \param[in]     x    - current state
    \param[in]     u    - current control
    \param[in,out] dt   - output transition time
    \param[in,out] gdt  - allocated space for gradient of dt
    \param[in]     grad - flag for specifying whether or not to
                          compute gradients of transitions with
                          respect to control (NULL for no)

    \return node

    \note 
    No funny boundary conditions
**************************************************************/
struct MCNode *
mca_inbound_node(struct MCA * mca, double time, double * x,
                 double * u, double *dt,double * gdt, double *grad)
{

    assert(mca->dyn != NULL);
    assert (mca->dw == mca->d);

    struct MCNode * mcn = mcnode_init(mca->d,x);
    
    //double * probs = calloc_double(2*mca->d);
    double * dQ = NULL;

    int res;
    size_t du = dyn_get_du(mca->dyn);
    size_t ngrad = 0;
    if (grad != NULL){
        ngrad = du;
        mcnode_init_self_grad(mcn,du);
            
        //gprob = malloc_dd(2*mca->d);
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       mca->gdrift,mca->diff,mca->gdiff);
        dQ = calloc_double(du);

    }
    else{
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       NULL,mca->diff,NULL);
    }

    assert(res == 0);
    
//    printf("got through first part\n");
    double diff,t;
    double Qh = 0.0;
    double minh2 = pow(mca->minh,2);
    double hrel;
    for (size_t ii = 0; ii < mca->d; ii++){
        // prepend the - and + nodes.
        // mcn->neigh points to plus node
        // mcn->neigh->next points to minus node
        mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]-mca->h[ii],ngrad);
        mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]+mca->h[ii],ngrad);
        struct MCNList * node_plus = mcn->neigh;
        struct MCNList * node_minus = mcn->neigh->next;
            
        hrel = (minh2/mca->h[ii]);

        diff = pow(mca->diff[ii*mca->d+ii],2);
        // general case (also needs further modifications)
        //diff = cblas_ddot(mca->dw,
        //                  mca->diff+ii,mca->d,
        //                  mca->diff+ii,mca->d);
        t = pow(mca->minh/mca->h[ii],2);
        
        node_minus->p = t * diff/2.0;
        node_plus->p  = t * diff/2.0;

        if (grad != NULL){
            cblas_daxpy(du,t*2.0,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,dQ,1);
            cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                        node_plus->gradp,1);
            cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                        node_minus->gradp,1);
        }
        
        Qh += t*diff;
        
        if (mca->drift[ii] > 0){
            node_plus->p += hrel * mca->drift[ii];
            Qh += hrel*mca->drift[ii];
            if (grad != NULL){
                cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,node_plus->gradp,1);
                cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,dQ,1);
            }

        }
        else{
            node_minus->p += hrel * (-mca->drift[ii]);
            Qh += hrel * (-mca->drift[ii]);
            if (grad != NULL){
                cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,node_minus->gradp,1);
                cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,dQ,1);
            }
        }
        // printf("%G\n",mcn->neigh->p);
//        dprint(du,node_minus)
    }

//    printf("got through second part \n");
//    printf("Qh = %G\n",Qh);
    *dt = minh2 / Qh;

    double pself = 1.0;
    if (grad != NULL){
        double Qh2 = pow(Qh,2);
        for (size_t jj = 0; jj < du; jj++){
//            printf("dQ[%zu]=%3.15G\n",jj,dQ[jj]);
            gdt[jj] = -minh2 * dQ[jj] / Qh2;
        }
        struct MCNList * temp = mcn->neigh;
        while (temp != NULL){
            for (size_t jj = 0; jj < du; jj++){
                temp->gradp[jj] = (Qh * temp->gradp[jj] - dQ[jj]*temp->p)/Qh2;
            }
            temp->p /= Qh;
            pself -= temp->p;
            
            cblas_daxpy(du,-1.0,temp->gradp,1,mcn->gpself,1);
            
            temp = temp->next; 
        }
        free(dQ); dQ = NULL;
    }
    else{
        struct MCNList * temp = mcn->neigh;
        while (temp != NULL){
            temp->p /= Qh;
            pself -= temp->p;
            temp = temp->next; 
        }
    }
//    printf("got through last part\n");
    if (pself < 0){
        pself = 0.0;
    }
    mcnode_set_pself(mcn,pself);
//    printf("return\n");
    return mcn;
}

void mca_upwind_dim(struct MCA * mca, struct MCNode * mcn,
                    double * x, size_t ngrad,
                    double minh2, double * Qh, double * dQ, size_t ii,
                    double * grad)
{
    size_t du = mca_get_du(mca);
    // prepend the - and + nodes.
    // mcn->neigh points to plus node
    // mcn->neigh->next points to minus node
    mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]-mca->h[ii],ngrad);
    mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]+mca->h[ii],ngrad);
    struct MCNList * node_plus = mcn->neigh;
    struct MCNList * node_minus = mcn->neigh->next;
            
    double hrel = (minh2/mca->h[ii]);

    double diff = pow(mca->diff[ii*mca->d+ii],2);
    // general case (also needs further modifications)
    //diff = cblas_ddot(mca->dw,
    //                  mca->diff+ii,mca->d,
    //                  mca->diff+ii,mca->d);
    double t = pow(mca->minh/mca->h[ii],2);
        
    node_minus->p = t * diff/2.0;
    node_plus->p  = t * diff/2.0;

    if (grad != NULL){
        cblas_daxpy(du,t*2.0,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,dQ,1);
        cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                    node_plus->gradp,1);
        cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                    node_minus->gradp,1);
    }
        
    *Qh += t*diff;
        
    if (mca->drift[ii] > 0){
        node_plus->p += hrel * mca->drift[ii];
        *Qh += hrel*mca->drift[ii];
        if (grad != NULL){
            cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,node_plus->gradp,1);
            cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,dQ,1);
        }

    }
    else{
        node_minus->p += hrel * (-mca->drift[ii]);
        *Qh += hrel * (-mca->drift[ii]);
        if (grad != NULL){
            cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,node_minus->gradp,1);
            cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,dQ,1);
        }
    }    

}

/**********************************************************//**
    Generate a node that is part of a reflective boundary and its 
    transitions
    
    \param[in]     mca  - markov chain approximation 
    \param[in]     time - time
    \param[in]     x    - current state
    \param[in]     u    - current control
    \param[in,out] dt   - output transition time
    \param[in,out] gdt  - allocated space for gradient of dt
    \param[in]     grad - flag for specifying whether or not to
                          compute gradients of transitions with
                          respect to control (NULL for no)
    \param[in]     bi   - boundary info

    \return node

    \note 
    No funny boundary conditions
**************************************************************/
struct MCNode *
mca_reflect_node(struct MCA * mca, double time, double * x,
                 double * u, double *dt,double * gdt, double *grad,
                 struct BoundInfo * bi)
{


    assert(mca->dyn != NULL);
    assert (mca->dw == mca->d);

    struct MCNode * mcn = mcnode_init(mca->d,x);
    
    //double * probs = calloc_double(2*mca->d);
    double * dQ = NULL;

    int res;
    size_t du = dyn_get_du(mca->dyn);
    size_t ngrad = 0;
    if (grad != NULL){
        ngrad = du;
        mcnode_init_self_grad(mcn,du);
            
        //gprob = malloc_dd(2*mca->d);
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       mca->gdrift,mca->diff,mca->gdiff);
        dQ = calloc_double(du);

    }
    else{
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       NULL,mca->diff,NULL);
    }

    assert(res == 0);
    
//    printf("got through first part\n");
    /* printf("this node is reflecting "); */
    /* dprint(mca->d,x); */
    double Qh = 0.0;
    double minh2 = pow(mca->minh,2);
    double hrel;
    for (size_t ii = 0; ii < mca->d; ii++){
        //printf("ii=%zu\n",ii);
        int on_bound = bound_info_onbound_dim(bi,ii);

        if (on_bound == 0){
            mca_upwind_dim(mca,mcn,x,ngrad,minh2,&Qh,dQ,ii,grad);
        }
        else{
            // general case needs modification
            double diff = pow(mca->diff[ii*mca->d+ii],2);
            double t = pow(mca->minh/mca->h[ii],2);
            Qh += t*diff;
            hrel = (minh2/mca->h[ii]);
            int reflect = bound_info_reflect(bi);
            assert (reflect == 1);
            int reflect_dir = bound_info_reflect_dim_dir(bi,ii);
            //printf("before\n");
            if (reflect_dir == -1){ // just do right side poition
                mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]+mca->h[ii],ngrad);
                struct MCNList * node_plus = mcn->neigh;
                node_plus->p  = t * diff/2.0;

                if (grad != NULL){
                    cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                                dQ,1);
                    cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                                node_plus->gradp,1);
                }
                if (mca->drift[ii] > 0){
                    node_plus->p += hrel * mca->drift[ii];
                    Qh += hrel*mca->drift[ii];
                    if (grad != NULL){
                        cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,
                                    node_plus->gradp,1);
                        cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,dQ,1);
                    }
                }
                else{
                    Qh += hrel * (-mca->drift[ii]);
                    if (grad != NULL){
                        cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,dQ,1);
                    }
                }
            }
            else{ // just do left side portion
                mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]-mca->h[ii],ngrad);
                struct MCNList * node_minus = mcn->neigh;
                node_minus->p  = t * diff/2.0;
                if (grad != NULL){
                    cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                                dQ,1);
                    cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                                node_minus->gradp,1);
                }
                if (mca->drift[ii] < 0){
                    node_minus->p += hrel * (-mca->drift[ii]);
                    Qh += hrel * (-mca->drift[ii]);
                    if (grad != NULL){
                        cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,
                                    node_minus->gradp,1);
                        cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,dQ,1);
                    }
                }
                else{
                    Qh += hrel*mca->drift[ii];
                    if (grad != NULL){
                        cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,dQ,1);
                    }
                }
            }
        }
    }

    *dt = minh2 / Qh;
    double pself = 1.0;
    if (grad != NULL){
        double Qh2 = pow(Qh,2);
        for (size_t jj = 0; jj < du; jj++){
            gdt[jj] = -minh2 * dQ[jj] / Qh2;
        }
        struct MCNList * temp = mcn->neigh;
        while (temp != NULL){
            for (size_t jj = 0; jj < du; jj++){
                temp->gradp[jj] = (Qh * temp->gradp[jj] - dQ[jj]*temp->p)/Qh2;
            }
            temp->p /= Qh;
            pself -= temp->p;
            
            cblas_daxpy(du,-1.0,temp->gradp,1,mcn->gpself,1);
            
            temp = temp->next; 
        }
        free(dQ); dQ = NULL;
    }
    else{
        struct MCNList * temp = mcn->neigh;
        while (temp != NULL){
            temp->p /= Qh;
            pself -= temp->p;
            temp = temp->next; 
        }
    }
    if (pself < 0){
        pself = 0.0;
    }
    mcnode_set_pself(mcn,pself);
    return mcn;
}

/**********************************************************//**
    Generate a node that is part of a periodic boundary and its 
    transitions
    
    \param[in]     mca  - markov chain approximation 
    \param[in]     time - time
    \param[in]     x    - current state
    \param[in]     u    - current control
    \param[in,out] dt   - output transition time
    \param[in,out] gdt  - allocated space for gradient of dt
    \param[in]     grad - flag for specifying whether or not to
                          compute gradients of transitions with
                          respect to control (NULL for no)
    \param[in]     bi   - boundary info

    \return node

    \note 
    No funny boundary conditions
**************************************************************/
struct MCNode *
mca_period_node(struct MCA * mca, double time, double * x,
                 double * u, double *dt,double * gdt, double *grad,
                 struct BoundInfo * bi)
{


    assert(mca->dyn != NULL);
    assert (mca->dw == mca->d);

    struct MCNode * mcn = mcnode_init(mca->d,x);
    
    //double * probs = calloc_double(2*mca->d);
    double * dQ = NULL;

    int res;
    size_t du = dyn_get_du(mca->dyn);
    size_t ngrad = 0;
    if (grad != NULL){
        ngrad = du;
        mcnode_init_self_grad(mcn,du);
            
        //gprob = malloc_dd(2*mca->d);
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       mca->gdrift,mca->diff,mca->gdiff);
        dQ = calloc_double(du);
    }
    else{
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       NULL,mca->diff,NULL);
    }

    assert(res == 0);
    
//    printf("got through first part\n");
    /* printf("this node is reflecting "); */
    /* dprint(mca->d,x); */
    double Qh = 0.0;
    double minh2 = pow(mca->minh,2);
    double hrel;
    for (size_t ii = 0; ii < mca->d; ii++){
        //printf("ii=%zu\n",ii);
        int on_bound = bound_info_onbound_dim(bi,ii);

        if (on_bound == 0){
            mca_upwind_dim(mca,mcn,x,ngrad,minh2,&Qh,dQ,ii,grad);
        }
        else{
            // general case needs modification
            double diff = pow(mca->diff[ii*mca->d+ii],2);
            double t = pow(mca->minh/mca->h[ii],2);
            Qh += t*diff;
            hrel = (minh2/mca->h[ii]);
            int period = bound_info_period(bi);
            assert (period == 1);
            int period_dir = bound_info_period_dim_dir(bi,ii);
            double xmap = bound_info_period_xmap(bi,ii);
            //printf("before\n");
            if (period_dir == -1){ // left maps to something
                mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]+mca->h[ii],ngrad);
                mcnlist_prepend_empty(&(mcn->neigh),ii,xmap-mca->h[ii],ngrad);
            }
            else{ // just do left side portion
                mcnlist_prepend_empty(&(mcn->neigh),ii,xmap+mca->h[ii],ngrad);
                mcnlist_prepend_empty(&(mcn->neigh),ii,x[ii]-mca->h[ii],ngrad);
            }

            struct MCNList * node_plus = mcn->neigh;
            struct MCNList * node_minus = mcn->neigh->next;
            node_minus->p = t * diff/2.0;
            node_plus->p  = t * diff/2.0;

            if (grad != NULL){
                cblas_daxpy(du,t*2.0,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,dQ,1);
                cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                            node_plus->gradp,1);
                cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                            node_minus->gradp,1);
            }
        
            if (mca->drift[ii] > 0){
                node_plus->p += hrel * mca->drift[ii];
                Qh += hrel*mca->drift[ii];
                if (grad != NULL){
                    cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,node_plus->gradp,1);
                    cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,dQ,1);
                }

            }
            else{
                node_minus->p += hrel * (-mca->drift[ii]);
                Qh += hrel * (-mca->drift[ii]);
                if (grad != NULL){
                    cblas_daxpy(du,-hrel,mca->gdrift+ii,
                                mca->d,node_minus->gradp,1);
                    cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,dQ,1);
                }
            }
        }
    }
    *dt = minh2 / Qh;
    double pself = 1.0;
    if (grad != NULL){
        double Qh2 = pow(Qh,2);
        for (size_t jj = 0; jj < du; jj++){
            gdt[jj] = -minh2 * dQ[jj] / Qh2;
        }
        struct MCNList * temp = mcn->neigh;
        while (temp != NULL){
            for (size_t jj = 0; jj < du; jj++){
                temp->gradp[jj] = (Qh * temp->gradp[jj] - dQ[jj]*temp->p)/Qh2;
            }
            temp->p /= Qh;
            pself -= temp->p;
            
            cblas_daxpy(du,-1.0,temp->gradp,1,mcn->gpself,1);
            
            temp = temp->next; 
        }
        free(dQ); dQ = NULL;
    }
    else{
        struct MCNList * temp = mcn->neigh;
        while (temp != NULL){
            temp->p /= Qh;
            pself -= temp->p;
            temp = temp->next; 
        }
    }
    if (pself < 0){
        pself = 0.0;
    }
    mcnode_set_pself(mcn,pself);
    return mcn;
}

/**********************************************************//**
    Generate a node of the Markov chain
    
    \param[in]     mca    - markov chain approximation 
    \param[in]     time   - time
    \param[in]     x      - current state
    \param[in]     u      - current control
    \param[in,out] dt     - time step
    \param[in,out] gdt    - gradient of time step
    \param[in,out] bi     - boundary info
    \param[in]     grad   - flag for specifying whether or not to
                            compute gradients of transitions with 
                            respect to control (NULL for no)
**************************************************************/
struct MCNode *
mca_get_node(struct MCA * mca, double time, double * x,
             double * u, double * dt, double * gdt,
             struct BoundInfo ** bi,
             double *grad)
{

    *bi = boundary_type(mca->bound,time,x);
    struct MCNode * mcn = NULL;
    

//    printf("get node\n"); dprint(mca->d, x);
    int onbound = bound_info_onbound(*bi);
//    printf("onbound = %d\n",onbound);
    /* if (x[1] < -2.0){ */
    /*     printf("what!\n"); */
    /*     if (onbound == 0){ */
    /*         printf("get node\n"); dprint(mca->d, x); */
    /*         printf("onbound = %d\n",onbound); */
    /*         exit(1); */
    /*     } */
    /* } */
    if (onbound == 0){
        mcn = mca_inbound_node(mca,time,x,u,dt,gdt,grad);
    }
    else{
        int absorb = bound_info_absorb(*bi);
        //      printf("absorb = %d\n",absorb);
        if (absorb == 1){
            *dt = 0;
            if (gdt != NULL){
                for (size_t ii = 0; ii < mca->du; ii++){
                    gdt[ii] = 0.0;
                }
            }
            mcn = mca_outbound_node(mca,time,x);            
        }
        else{
            int reflect = bound_info_reflect(*bi);
            if (reflect != 0){
//                assert (1 == 0);
                mcn = mca_reflect_node(mca,time,x,u,dt,gdt,grad,*bi);
            }
            else{
                int period = bound_info_period(*bi);
                if (period != 0){
                    /* assert (1 == 0); */
                    mcn = mca_period_node(mca,time,x,u,dt,gdt,grad,*bi);
                }
                else{
                    fprintf(stderr,"Not sure what to do with this node\n");
                    printf("onbound = %d\n",onbound);
                    dprint(mca->d,x);
                    exit(1);                
                }
            }
        }
    }
    // printf("got node\n");
    //bound_info_free(bi); bi = NULL;
    return mcn;
}

/**********************************************************//**
    Compute expectation of a function
    
    \param[in]     mca  - markov chain approximation 
    \param[in]     time - time
    \param[in]     x    - current state
    \param[in]     u    - current control
    \param[in,out] dt   - time step
    \param[in]     f    - function(npts,times,states,out,arg)
    \param[in,out] arg  - function arguments
    \param[in,out] bi   - boundary information
    \param[in,out] grad - gradient
**************************************************************/
double
mca_expectation(struct MCA * mca, double time,
                double * x, double * u, double * dt, double * gdt,
//                void (*f)(size_t,double *,double **,double*,void*),
                double (*f)(double,double *,void*),
                void * arg, struct BoundInfo ** bi, double * grad)
{

//    printf("compute expectation\n");
    struct MCNode * mcn = mca_get_node(mca,time,x,u,dt,gdt,bi,grad);

//    printf("before onbound=%zu ",bound_info_onbound(*bi)); dprint(mca->d, x);
    double val = mcnode_expectation(mcn,f,arg,grad);
//    printf("after\n");
    mcnode_free(mcn); mcn = NULL;
//    printf("computed\n");
    return val;
}

/* /\**********************************************************\//\** */
/*     Sample a transition in the markov chain */
    
/*     \param[in]     mca   - markov chain approximation  */
/*     \param[in]     time  - time */
/*     \param[in]     x     - current state */
/*     \param[in]     u     - current control */
/*     \param[in]     noise - [0,1] uniform sample */
/*     \param[in,out] dt    - output transition time */
/*     \param[in,out] newx  - space for new state */
/* **************************************************************\/ */
/* void mca_step(struct MCA * mca, double time, double * x, double * u, */
/*          double noise, double * dt, double * newx) */
/* { */

/*     struct MCNode * node = mca_inbound_node(mca,time,x,u,dt,NULL,NULL);//,NULL); */
    
/*     mcnode_sample_neighbor(node,noise,newx); */
/*     mcnode_free(node); node = NULL; */
/* } */

///////////////////////////////////////////////////////////////


/* void mcnode_init_neighbors_interior_dim(struct MCNode * mcn, */
/*                                         size_t dim, double h) */
/* { */
/*     /\* mcn->neighbors[2*dim] = mcnode_init(mcn->d,mcn->x); *\/ */
/*     /\* mcn->neighbors[2*dim+1] = mcnode_init(mcn->d,mcn->x); *\/ */
/*     /\* mcn->neighbors[2*dim]->x[dim] = mcn->x[dim]-h; *\/ */
/*     /\* mcn->neighbors[2*dim+1]->x[dim] = mcn->x[dim]+h; *\/ */
/* } */

/* void mcnode_init_neighbors_periodr_dim(struct MCNode * mcn, */
/*                                        size_t dim, double h, */
/*                                        double xmap) */
/* { */
/*     mcn->neighbors[2*dim] = mcnode_init(mcn->d,mcn->x); */
/*     mcn->neighbors[2*dim+1] = mcnode_init(mcn->d,mcn->x); */
/*     mcn->neighbors[2*dim]->x[dim] = mcn->x[dim]-h; */
/*     mcn->neighbors[2*dim+1]->x[dim] = xmap+h; */
/* } */

/* void mcnode_init_neighbors_periodl_dim(struct MCNode * mcn, */
/*                                        size_t dim, double h, */
/*                                        double xmap) */
/* { */
/*     mcn->neighbors[2*dim] = mcnode_init(mcn->d,mcn->x); */
/*     mcn->neighbors[2*dim+1] = mcnode_init(mcn->d,mcn->x); */
/*     mcn->neighbors[2*dim]->x[dim] = xmap-h; */
/*     mcn->neighbors[2*dim+1]->x[dim] = mcn->x[dim]+h;     */
/* } */

/* /\**********************************************************\//\** */
/*     Add neighbors to a node that are a distance of  */
/*     *h* away in each coordinate direction */

/*     \param[in,out] mcn   - node */
/*     \param[in]     h     - spacing of nodes in each direction */
/*     \param[in]     pself - new probability of self transition */
/*     \param[in]     probs - transition probs to neighbors */
/*     \param[in]     bi    - boundary info */
/* **************************************************************\/ */
/* void mcnode_add_neighbors_hspace(struct MCNode * mcn,  */
/*                                  const double * h, */
/*                                  double pself, */
/*                                  const double * probs, */
/*                                  const struct BoundInfo * bi) */
/* { */
/*     assert (mcn != NULL); */
/*     assert (mcn->N == 0); */
/*     mcn->N = 2*mcn->d; */
/*     mcn->neighbors = mcnode_alloc_array(2*mcn->d); */
/*     for (size_t ii = 0; ii < mcn->d; ii++ ){ */
/*         int period_dir = 0; */
/*         if (bi == NULL){ */
/*             mcnode_init_neighbors_interior_dim(mcn,ii,h[ii]); */
/*         } */
/*         else{ */
/*             period_dir = bound_info_period_dim_dir(bi,ii); */
/*             if (period_dir == 1){ // to the right is periodic */
/*                 double xmap = bound_info_period_xmap(bi,ii); */
/*                 mcnode_init_neighbors_periodr_dim(mcn,ii, */
/*                                                   h[ii],xmap); */
/*             } */
/*             else if (period_dir == -1){ */
/*                 double xmap = bound_info_period_xmap(bi,ii); */
/*                 mcnode_init_neighbors_periodl_dim(mcn,ii, */
/*                                                   h[ii],xmap); */
/*             } */
/*             else{ // in interior */
/*                 mcnode_init_neighbors_interior_dim(mcn,ii,h[ii]); */
/*             } */
/*         } */
/*     } */
/*     mcn->pself = pself; */
/*     mcn->p = calloc_double(2*mcn->d); */
/*     memmove(mcn->p,probs,2*mcn->d * sizeof(double)); */
/* } */



/* /\**********************************************************\//\** */
/*     Add gradients of probability transitions to neighbors with  */
/*     respect to control   */

/*     \param[in,out] mcn    - node */
/*     \param[in]     N      - number of gradients to neighbors */
/*     \param[in]     du     - size of gradient for each transition */
/*     \param[in]     gpself - gradient of transition to oneself */
/*     \param[in]     gprob  - gradients of transition to neighbors */
/* **************************************************************\/ */
/* void mcnode_add_gradients(struct MCNode * mcn, size_t N, */
/*                           size_t du, */
/*                           const double *gpself, */
/*                           double ** gprob) */
/* { */
/*     assert (mcn != NULL); */
/*     assert (mcn->N == N); */
/*     mcn->du = du; */
    
/*     free_dd(mcn->N,mcn->gp); mcn->gp = NULL; */
/*     mcn->gp = malloc_dd(mcn->N); */
/*     for (size_t ii = 0; ii < N; ii++){ */
/*         mcn->gp[ii] = calloc_double(du); */
/*         memmove(mcn->gp[ii],gprob[ii],du*sizeof(double)); */
/*     } */

/*     free(mcn->gpself); mcn->gpself = NULL; */
/*     mcn->gpself = calloc_double(du); */
/*     memmove(mcn->gpself,gpself, du*sizeof(double)); */
/* } */



/**********************************************************//**
    Sample a transition to neighbors or oneself

    \param[in]     node     - node
    \param[in]     noise    - sample from [0,1]
    \param[in,out] neighbor - (node->d) sized allocated array
**************************************************************/
/* void mcnode_sample_neighbor(struct MCNode * node, double noise, */
/*                             double * neighbor) */
/* { */
/*     assert (node != NULL); */
/*     assert (noise >= 0.0); */
/*     assert (noise <= 1.0); */
    
/*     double * tprobs = calloc_double(node->N + 1); */
/*     tprobs[0] = node->pself; */
/*     if (node->N > 0){ */
/*         memmove(tprobs+1,node->p, node->N*sizeof(double)); */
/*     } */

/*     size_t ind = c3sc_sample_discrete_rv(node->N+1,tprobs,noise); */
/*     if (ind == 0){ */
/*         memmove(neighbor,node->x,node->d*sizeof(double)); */
/*     } */
/*     else{ */
/*         //printf("whats up\n"); */
/*         double * nx = mcnode_get_xref(node->neighbors[ind-1]); */
/*         //dprint(node->d,nx); */
/*         memmove(neighbor,nx,node->d*sizeof(double)); */
/*     } */

/* } */

/**********************************************************//**
    Allocate node array

    \param[in] N - Number of ndoes

    \return Array of N nodes set to NULL;
**************************************************************/
struct MCNode ** mcnode_alloc_array(size_t N)
{
    struct MCNode ** mcn;
    mcn = malloc(N * sizeof(struct MCNode *));
    assert (mcn != NULL);
    for (size_t ii = 0; ii < N; ii++){
        mcn[ii] = NULL;
    }
    
    return mcn;
}

/**********************************************************//**
    Free an array of N nodes
**************************************************************/
void mcnode_free_array(struct MCNode ** mcn, size_t N)
{
    if (mcn != NULL){
        for (size_t ii = 0; ii < N; ii++){
            mcnode_free(mcn[ii]); mcn[ii] = NULL;
        }
        //free(mcn); mcn = NULL;
    }
}

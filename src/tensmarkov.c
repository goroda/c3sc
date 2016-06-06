#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "tensmarkov.h"
#include "util.h"
#include "dynamics.h"
#include "cost.h"

/* struct Node */
/* { */
/*     size_t d; */
/*     double * x; */
/*     double c; // cost */
/*     double * pert; // (- +) in each direction //2d */
/*     double * cost; */
/* }; */

struct Node * node_alloc(size_t d, const double * x)
{
    struct Node * n = malloc(sizeof(struct Node));
    assert (n != NULL);
    n->d = d;
    n->x = calloc_double(d);
    memmove(n->x,x,d*sizeof(double));
    n->pert = NULL;
    n->use = calloc_int(2*d);
    n->cost = NULL;
    return n;
}

void node_free(struct Node * node)
{
    if (node != NULL){
        free(node->x); node->x = NULL;
        free(node->pert); node->pert = NULL;
        free(node->use); node->use = NULL;
        free(node->cost); node->cost = NULL;
        free(node); node = NULL;
    }
}

struct Node * node_init(size_t d, const double * x, const double * h,
                        struct Cost * cost, struct BoundInfo * bi)
{
    int onbound = bound_info_onbound(bi);
    struct Node * n = node_alloc(d,x);
    /* printf("allocated onbound=%d\n",onbound); */
    if (onbound == 0){
        n->pert = calloc_double(2*d);
        n->cost = calloc_double(2*d);
        for (size_t ii = 0; ii < d; ii++){
            n->pert[ii*2]   = x[ii] - h[ii];
            n->pert[ii*2+1] = x[ii] + h[ii];
            n->use[ii*2]   = 1;
            n->use[ii*2+1] = 1;
        }
        int res = cost_eval_neigh(cost,0.0,n->x,&(n->c),n->pert,n->cost);
        assert (res == 0);
    }
    else{
        int absorb = bound_info_absorb(bi);
        if (absorb == 1){
            n->pert = NULL;
            n->cost = NULL;
            int res = cost_eval(cost,0.0,n->x,&(n->c));
            assert (res == 0);
        }
        else{
            n->pert = calloc_double(2*d);
            n->cost = calloc_double(2*d);
            int reflect = bound_info_reflect(bi);
            if (reflect != 0){
                for (size_t ii = 0; ii < d; ii++){
                    int on_bound = bound_info_onbound_dim(bi,ii);
                    if (on_bound == 0){
                        n->pert[ii*2]   = x[ii] - h[ii];
                        n->pert[ii*2+1] = x[ii] + h[ii];
                        n->use[2*ii]   = 1;
                        n->use[2*ii+1] = 1;
                    }
                    else{
                        int reflect_dir = bound_info_reflect_dim_dir(bi,ii);
                        //printf("before\n");
                        if (reflect_dir == -1){ // just do right side poition
                            n->pert[ii*2] = x[ii];
                            n->pert[ii*2+1] = x[ii]+h[ii];
                            n->use[2*ii+1] = 1;
                        }
                        else{
                            n->pert[ii*2] = x[ii] - h[ii];
                            n->pert[ii*2+1] = x[ii];
                            n->use[2*ii] = 1;
                        }
                    }
                }
                int res = cost_eval_neigh(cost,0.0,n->x,&(n->c),n->pert,n->cost);
                assert (res == 0);
            }
            else{
                int period = bound_info_period(bi);
                if (period != 0){
                    for (size_t ii = 0; ii < d; ii++){
                        // use both no matter what
                        n->use[2*ii]   = 1;
                        n->use[2*ii+1] = 1;
                        
                        int on_bound = bound_info_onbound_dim(bi,ii);
                        if (on_bound == 0){
                            n->pert[ii*2]   = x[ii] - h[ii];
                            n->pert[ii*2+1] = x[ii] + h[ii];
                        }
                        else{
                            int period_dir = bound_info_period_dim_dir(bi,ii);
                            double xmap = bound_info_period_xmap(bi,ii);
                            if (period_dir == -1){ // left maps to something
                                n->pert[ii*2]   = xmap  - h[ii];
                                n->pert[ii*2+1] = x[ii] + h[ii];
                            }
                            else{ // just do left side portion
                                n->pert[ii*2]   = x[ii] - h[ii];
                                n->pert[ii*2+1] = xmap + h[ii];
                            }
                        }
                    }
                    int res = cost_eval_neigh(cost,0.0,n->x,&(n->c),n->pert,n->cost);
                    assert (res == 0);
                }
                else{
                    fprintf(stderr,"Not sure what to do with this node\n");
                    printf("onbound = %d\n",onbound);
                    dprint(d,x);
                    exit(1);                
                }
            }
        }
    }
    return n;
}

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
    struct Node * node;

    size_t du;
    double p;
    double * gradp;
    double * pn; // probability transitions to neighbors;
    double * gradpn; // gradient of probability transitions to neighbors
    
};

/**********************************************************//**
    Allocate node

    \param[in] node - node
    \param[in] du   - number of controls

    \return allocate node
**************************************************************/
struct MCNode * mcnode_alloc(struct Node * node, size_t du)
{
    assert (node != NULL);
    struct MCNode * mcn;
    mcn = malloc(sizeof(struct MCNode));
    assert (mcn != NULL);
    
    mcn->node = node;

    /* printf("du = %zu\n",du); */
    /* printf("d = %zu\n",node->d); */
    mcn->du = du;
    mcn->gradp = calloc_double(du);
    mcn->pn = calloc_double(2*node->d);
    mcn->gradpn = calloc_double(2* du * node->d);
    
    return mcn;
}

/**********************************************************//**
    Free a node
**************************************************************/
void mcnode_free(struct MCNode * mcn)
{
    if (mcn != NULL){
        free(mcn->gradp); mcn->gradp = NULL;
        free(mcn->pn); mcn->pn = NULL;
        free(mcn->gradpn); mcn->gradpn = NULL;
        free(mcn); mcn = NULL;
    }
}

size_t mcnode_get_d(const struct MCNode * mcn){ return mcn->node->d; }
double * mcnode_get_xref(const struct MCNode * mcn){ return mcn->node->x;}
double mcnode_get_pself(const struct MCNode * mcn){ return mcn->p; }
size_t mcnode_get_du(const struct MCNode * mcn){ return mcn->du;}
double * mcnode_get_gpself(const struct MCNode * mcn){ return mcn->gradp;}

//////////////////////////////////

/**********************************************************//**
    Compute the expectation of a function around a node

    \param[in]     mc   - node
    \param[in,out] grad - If (!NULL) return gradient

    \return average
**************************************************************/
double mcnode_expectation2(
    const struct MCNode * mc,
    double * grad)
{

    assert (mc != NULL);
    assert (mc->node != NULL);
    
    if (grad!= NULL){
        for (size_t ii = 0; ii < mc->du; ii++){
            grad[ii] = 0.0;
        }
    }


    struct Node * node = mc->node;
    /* printf("lets go in exp !\n"); */
    /* iprint(2*node->d,node->use); */
    double avg = 0.0;
    for (size_t ii = 0; ii < 2*node->d; ii++){
        /* printf("ii = %zu\n",ii); */
        if (node->use[ii] == 1){
            /* printf("prob = %G, cost = %G\n",mc->pn[ii],node->cost[ii]); */
            /* assert(isinf(mc->pn[ii]) == 0); */
            avg += mc->pn[ii] * node->cost[ii];
            if (grad != NULL){
                cblas_daxpy(mc->du,node->cost[ii],mc->gradpn + ii*mc->du,1,grad,1);
            }
        }
    }
    /* printf("got it in exp \n"); */
    if (mc->p > 0.0){
        avg += mc->p * node->c; // self eval;
        if (grad != NULL){
            cblas_daxpy(mc->du,node->c,mc->gradp,1,grad,1);
        }
    }
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
    else if (n->node == NULL){
        fprintf(fp,"NULL");
        return;
    }

    fprintf(fp,"Node is ");
    for (size_t ii = 0; ii < n->node->d; ii++){
        fprintf(fp,"%3.*G ",prec,n->node->x[ii]);
    }
    fprintf(fp,"\n");

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
    
    // workspace
    mca->drift = calloc_double(d);
    mca->gdrift = calloc_double(d * du);
    mca->diff = calloc_double(d*dw);
    mca->gdiff = calloc_double(d * dw * du);
    
    mca->dyn = NULL;
    mca->bound = NULL;

    return mca;
}

/**********************************************************//**
    Copy a markov model
**************************************************************/
struct MCA * mca_copy_deep(struct MCA * mca)
{
    if (mca == NULL){
        return NULL;
    }

    struct MCA * newm = mca_alloc(mca->d,mca->du,mca->dw,mca->h);
    
    newm->dyn = dyn_copy_deep(mca->dyn);
    newm->bound = boundary_copy_deep(mca->bound);

    return newm;
    
}

/**********************************************************//**
    Set a new discretization level for MCA
**************************************************************/
void mca_set_newh(struct MCA * mca, const double * h)
{
    free(mca->h); mca->h = NULL;
    mca->h = calloc_double(mca->d);
    mca->minh = h[0];
    for (size_t ii = 0; ii < mca->d; ii++){
        mca->h[ii] = h[ii];
        if (h[ii] < mca->minh){
            mca->minh = h[ii];
        }
    }
}

/**********************************************************//**
    Generate a new markov model by either cutting h in half (*inhalf=1*)
    or multiplying it by 2 (*inhalf=0*)
**************************************************************/
struct MCA * mca_interp2(struct MCA * old, int inhalf)
{

    double * h = calloc_double(old->d);
    if (inhalf == 1){
        for (size_t ii = 0;ii < old->d; ii++){
            h[ii] = old->h[ii] / 2.0;
        }
    }
    else{
        for (size_t ii = 0; ii < old->d; ii++){
            h[ii] = h[ii] * 2.0;
        }
    }

    struct MCA * mnew = mca_copy_deep(old);
    mca_set_newh(mnew,h);
    free(h); h = NULL;
    return NULL;
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
    Get a reference to the boundary
**************************************************************/
struct Boundary * mca_get_boundary(struct MCA * mca)
{
    assert (mca != NULL);
    return mca->bound;
}

/**********************************************************//**
    Get distances between nodes
**************************************************************/
double * mca_get_h(struct MCA * mca)
{
    assert (mca != NULL);
    return mca->h;
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
    Generate a node of the Markov chain
    
    \param[in]     mca    - markov chain approximation 
    \param[in]     time   - time
    \param[in]     node   - graph node information
    \param[in]     u      - current control
    \param[in,out] dt     - time step
    \param[in,out] gdt    - gradient of time step
    \param[in]     grad   - flag for specifying whether or not to
                            compute gradients of transitions with 
                            respect to control (NULL for no)
**************************************************************/
struct MCNode *
mca_get_node2(struct MCA * mca, double time, struct Node * node,
              const double * u, double * dt, double * gdt,
              double *grad)
{

    assert(mca->dyn != NULL);
    assert (mca->dw == mca->d);

    size_t du = dyn_get_du(mca->dyn);

    /* printf("alloc node\n"); */
    struct MCNode * mcn = mcnode_alloc(node,du);
    /* printf("allocated\n"); */

    int res;
    double * dQ = NULL;
    size_t ngrad = 0;
    if (grad != NULL){
        ngrad = du;
        //gprob = malloc_dd(2*mca->d);
        res = dyn_eval(mca->dyn,time,node->x,u,mca->drift,
                       mca->gdrift,mca->diff,mca->gdiff);
        dQ = calloc_double(du);
    }
    else{
        res = dyn_eval(mca->dyn,time,node->x,u,mca->drift,
                       NULL,mca->diff,NULL);
    }


    /* printf("drift = "); dprint(mca->d,mca->drift); */
    
    double Qh = 0.0;
    double minh2 = pow(mca->minh,2);
    double hrel;
    for (size_t ii = 0; ii < mca->d; ii++){
        hrel = (minh2/mca->h[ii]);
        double diff = pow(mca->diff[ii*mca->d+ii],2);
        double t = pow(mca->minh/mca->h[ii],2);
        /* printf("ii=%zu, hrel = %G\n",ii,hrel); */
        
        // minus node
        if (node->use[2*ii] == 1){
            mcn->pn[2*ii] = t * diff/2.0;
            Qh += t/2.0*diff;
            if (grad != NULL){
                cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,dQ,1);
                cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                            mcn->gradpn + (2*ii)*du,1);
            }
            if (mca->drift[ii] <= 0){
                /* printf("drift less than zero\n"); */
                mcn->pn[2*ii] += hrel * (-mca->drift[ii]);
                Qh += hrel * (-mca->drift[ii]);
                if (grad != NULL){

                    cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,mcn->gradpn + (2*ii)*du,1);
                    cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,dQ,1);
                    /* printf("computed gradpn = "); */
                    /* dprint(du,mcn->gradpn + 2*ii*du); */
                }
            }
        }

        // plus node
        if (node->use[2*ii+1] == 1){
            mcn->pn[2*ii+1] = t * diff/2.0;
            Qh += t/2.0*diff;
            if (grad != NULL){
                cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,dQ,1);
                cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                            mcn->gradpn + (2*ii+1)*du,1);
            }
            if (mca->drift[ii] > 0){
                /* printf("drift greather than zero\n"); */
                mcn->pn[2*ii+1] += hrel * mca->drift[ii];
                Qh += hrel*mca->drift[ii];
                if (grad != NULL){
                    cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,mcn->gradpn + (2*ii+1)*du,1);
                    cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,dQ,1);
                    /* printf("computed gradpn = "); */
                    /* dprint(du,mcn->gradpn + (2*ii+1)*du); */
                }
            }
        }
    }

    if (fabs(Qh) < 1e-14){
        *dt = 0;
        mcn->p = 1.0;
        free(dQ); dQ = NULL;
        return mcn;
    }
    // dela with normalization
    *dt = minh2 / Qh;
    mcn->p = 1.0;
    if (grad != NULL){
        double Qh2 = pow(Qh,2);
        for (size_t jj = 0; jj < du; jj++){
            gdt[jj] = -minh2 * dQ[jj] / Qh2;
        }
        for (size_t ii = 0; ii < 2*mca->d; ii++){
            if (node->use[ii] == 1){
                for (size_t jj = 0; jj < du; jj++){
                    mcn->gradpn[ii*du+jj] = (Qh * mcn->gradpn[ii*du+jj] - dQ[jj]*mcn->pn[ii])/Qh2;
                }
                mcn->pn[ii] /= Qh;
                mcn->p -= mcn->pn[ii];
                cblas_daxpy(du,-1.0,mcn->gradpn+ii*du,1,mcn->gradp,1);
            }
        }
        free(dQ); dQ = NULL;
    }
    else{
        for (size_t ii = 0; ii < 2*mca->d; ii++){
            if (node->use[ii] == 1){
                mcn->pn[ii] /= Qh;
                mcn->p -= mcn->pn[ii];
            }
        }
    }

    /* printf("u = "); dprint(du,u); */
    /* printf("transition probabilities are "); */
    /* dprint(2*mca->d,mcn->pn); */
    /* printf("grad == NULL = %d\n", grad == NULL); */
    /* printf("gradients are\n"); */
    /* for (size_t ii = 0; ii < 2*mca->d; ii++){ */
    /*     dprint(du,mcn->gradpn+ii*du); */
    /* } */

    return mcn;
}

/**********************************************************//**
    Compute expectation of a cost function
    
    \param[in]     mca  - markov chain approximation 
    \param[in]     time - time
    \param[in]     node - current state with cost info for neighbors
    \param[in]     u    - current control
    \param[in,out] dt   - time step
    \param[in,out] grad - gradient
    \param[in,out] info - 0 if successfull 1 if error

    \return expectation of f
**************************************************************/
double
mca_expectation_cost2(struct MCA * mca, double time,
                      struct Node * node,
                      const double * u,
                      double * dt, double * gdt,
                      double * grad,
                      int * info)
{

//    printf("compute expectation\n");
    double * x = node->x;
    if (grad != NULL){
        int res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                           mca->gdrift,mca->diff,mca->gdiff);
        assert (res == 0);
    }
    else{
        int res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                           NULL,mca->diff,NULL);
        assert (res == 0);
    }

    /* printf("get mcn\n"); */
    struct MCNode * mcn = mca_get_node2(mca,time,node,u,dt,gdt,grad);
    /* printf("got it\n"); */

    double val = 0.0;
    if (mcn == NULL){
        fprintf(stderr, "Warning: MCNode obtained in mca_expectation is NULL\n");
        fprintf(stderr, "\tx = ");
        for (size_t ii = 0; ii < mca->d; ii++ ){
            fprintf(stderr,"%G ",x[ii]);
        }
        fprintf(stderr,"\n\tControl = ");
        for (size_t ii = 0; ii < mca->du; ii++){
            fprintf(stderr,"%G ",u[ii]);
        }
        fprintf(stderr,"\n");
        fprintf(stderr,"\tBoundary info is ?????\n");
        *info = 1;
    }
    else{
        /* printf("compute mcnode expectation\n"); */
        val = mcnode_expectation2(mcn,grad);
        /* printf("computed\n"); */
        mcnode_free(mcn); mcn = NULL;
        *info = 0;
    }
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
/* struct MCNode ** mcnode_alloc_array(size_t N) */
/* { */
/*     struct MCNode ** mcn; */
/*     mcn = malloc(N * sizeof(struct MCNode *)); */
/*     assert (mcn != NULL); */
/*     for (size_t ii = 0; ii < N; ii++){ */
/*         mcn[ii] = NULL; */
/*     } */
    
/*     return mcn; */
/* } */

/**********************************************************//**
    Free an array of N nodes
**************************************************************/
/* void mcnode_free_array(struct MCNode ** mcn, size_t N) */
/* { */
/*     if (mcn != NULL){ */
/*         for (size_t ii = 0; ii < N; ii++){ */
/*             mcnode_free(mcn[ii]); mcn[ii] = NULL; */
/*         } */
/*         //free(mcn); mcn = NULL; */
/*     } */
/* } */

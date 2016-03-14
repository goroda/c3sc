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
    double * x;
    
    size_t N;
    struct MCNode ** neighbors;
    double pself;
    double * p; // transition probabilities

    size_t du; // dimension of gradients
    double ** gp;
    double * gpself;
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
    mcn->x = NULL;
    
    mcn->N = 0;
    mcn->neighbors = NULL;
    mcn->pself = 1.0;
    mcn->p = NULL;

    mcn->du = 0;
    mcn->gp = NULL;
    mcn->gpself = NULL;
    return mcn;
}

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

/**********************************************************//**
    Free a node
**************************************************************/
void mcnode_free(struct MCNode * mcn)
{
    if (mcn != NULL){
        free(mcn->x); mcn->x = NULL;
        free(mcn->p); mcn->p = NULL;
        free_dd(mcn->N,mcn->gp);
        free(mcn->gpself);
        if (mcn->N > 0){
            mcnode_free_array(mcn->neighbors,mcn->N);
            free(mcn->neighbors); mcn->neighbors = NULL;
        }
        free(mcn); mcn = NULL;
    }
}

/**********************************************************//**
    Initialize a node by copying an an array x of size d, 
    node has no neighbors and only transitions to itself
**************************************************************/
struct MCNode * mcnode_init(size_t d, const double * x)
{
    struct MCNode * mm = mcnode_alloc(d);
    assert (mm != NULL);
    mm->x = calloc_double(d);
    memmove(mm->x,x,d*sizeof(double));
    mm->N = 0;
    mm->neighbors = NULL;
    mm->pself = 1.0;
    mm->p = NULL;

    return mm;
}

/**********************************************************//**
    Get the dimension
**************************************************************/
size_t mcnode_get_d(const struct MCNode * mcn)
{
    return mcn->d;
}


/**********************************************************//**
    Get a reference to the node state
**************************************************************/
double * mcnode_get_xref(const struct MCNode * mcn)
{
    return mcn->x;
}

/**********************************************************//**
    Get the probability of self transition
**************************************************************/
double mcnode_get_pself(const struct MCNode * mcn)
{
    return mcn->pself;
}

/**********************************************************//**
    Get the number of neighbors
**************************************************************/
size_t mcnode_get_N(const struct MCNode * mcn)
{
    return mcn->N;
}

/**********************************************************//**
    Get a reference to the neighbors
**************************************************************/
struct MCNode **  mcnode_get_neighbors(const struct MCNode * mcn)
{
    return mcn->neighbors;
}

/**********************************************************//**
    Get a reference to the transition probabilities to neighbors
**************************************************************/
double *  mcnode_get_pref(const struct MCNode * mcn)
{
    return mcn->p;
}

/**********************************************************//**
    Get the dimension of the gradients
**************************************************************/
size_t mcnode_get_du(const struct MCNode * mcn)
{
    return mcn->du;
}

/**********************************************************//**
    Get the derivative of self transition probability
**************************************************************/
double * mcnode_get_gpself(const struct MCNode * mcn)
{
    return mcn->gpself;
}

/**********************************************************//**
    Get the derivative of transition probabilities to neighbors
**************************************************************/
double ** mcnode_get_gp(const struct MCNode * mcn)
{
    return mcn->gp;
}



/**********************************************************//**
    Add neighbors to a node that are a distance of 
    *h* away in each coordinate direction

    \param[in,out] mcn   - node
    \param[in]     h     - spacing of nodes in each direction
    \param[in]     pself - new probability of self transition
    \param[in]     probs - transition probs to neighbors
    \param[in]     bi    - boundary info
**************************************************************/
void mcnode_add_neighbors_hspace(struct MCNode * mcn, 
                                 const double * h, double pself,
                                 const double * probs,
                                 const struct BoundInfo * bi)
{
    assert (mcn != NULL);
    assert (mcn->N == 0);
    mcn->N = 2*mcn->d;
    mcn->neighbors = mcnode_alloc_array(2*mcn->d);
    for (size_t ii = 0; ii < mcn->d; ii++ ){
        mcn->neighbors[2*ii] = mcnode_init(mcn->d,mcn->x);
        mcn->neighbors[2*ii+1] = mcnode_init(mcn->d,mcn->x);

        int period_dir = 0;
        //int bigood = 0;
        if (bi == NULL){
            mcn->neighbors[2*ii]->x[ii] = mcn->x[ii]-h[ii];
            mcn->neighbors[2*ii+1]->x[ii] = mcn->x[ii]+h[ii];
        }
        else{
            //bigood = 1;
            period_dir = bound_info_period_dim_dir(bi,ii);
            if (period_dir == 1){ // to the right is periodic

                //printf("periodic right"); 
                //dprint(mcn->d,mcn->x);
                double xmap = bound_info_period_xmap(bi,ii);

                //printf("xmap = %G\n",xmap);
                mcn->neighbors[2*ii]->x[ii] = mcn->x[ii]-h[ii];
                mcn->neighbors[2*ii+1]->x[ii] = xmap+h[ii];
            }
            else if (period_dir == -1){
                //printf("periodic left "); 
                //dprint(mcn->d,mcn->x);
                double xmap = bound_info_period_xmap(bi,ii);
                //printf("xmap = %G\n",xmap);
                mcn->neighbors[2*ii]->x[ii] = xmap-h[ii];
                mcn->neighbors[2*ii+1]->x[ii] = mcn->x[ii]+h[ii];
            }
            else{
                mcn->neighbors[2*ii]->x[ii] = mcn->x[ii]-h[ii];
                mcn->neighbors[2*ii+1]->x[ii] = mcn->x[ii]+h[ii];
            }
        }
        /* double bound = 2.0; */
        /* if (fabs(mcn->neighbors[2*ii]->x[2]) > 2.0){ */
        /*     printf("going left ii=%zu\n",ii); */
        /*     printf("period dir = %d\n",period_dir); */
        /*     printf("bigood =%d\n",bigood); */
        /*     printf("absorb=%d\n",bound_info_absorb(bi)); */
        /*     printf("period=%d\n",bound_info_period(bi)); */
        /*     printf("onbound=%d\n",bound_info_onbound(bi)); */
        /*     mcnode_print(mcn,stdout,5); */
        /* } */
        /* if (fabs(mcn->neighbors[2*ii+1]->x[2]) > 2.0){ */
        /*     printf("going right=%zu\n",ii); */
        /*     mcnode_print(mcn,stdout,5); */
        /* } */
    }

    mcn->pself = pself;
    mcn->p = calloc_double(2*mcn->d);
    memmove(mcn->p,probs,2*mcn->d * sizeof(double));
}


void mcnode_init_self_grad(struct MCNode * mcn, size_t du)
{
    mcn->gpself = calloc_double(du);
}

/**********************************************************//**
    Add gradients of probability transitions to neighbors with 
    respect to control  

    \param[in,out] mcn    - node
    \param[in]     N      - number of gradients to neighbors
    \param[in]     du     - size of gradient for each transition
    \param[in]     gpself - gradient of transition to oneself
    \param[in]     gprob  - gradients of transition to neighbors
**************************************************************/
void mcnode_add_gradients(struct MCNode * mcn, size_t N,
                          size_t du,
                          const double *gpself,
                          double ** gprob)
{
    assert (mcn != NULL);
    assert (mcn->N == N);
    mcn->du = du;
    
    free_dd(mcn->N,mcn->gp);
    mcn->gp = malloc_dd(mcn->N);
    for (size_t ii = 0; ii < N; ii++){
        mcn->gp[ii] = calloc_double(du);
        memmove(mcn->gp[ii],gprob[ii],du*sizeof(double));
    }

    free(mcn->gpself);
    mcn->gpself = calloc_double(du);
    memmove(mcn->gpself,gpself, du*sizeof(double));
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
    void (*f)(size_t,double*,double**x,double*,void*),
    void * arg,double * grad)
{

    if (grad!= NULL){
        //printf("num neighbors = %zu\n",mc->N);
        if (mc->N > 0){
            assert(mc->gp != NULL);
        }
        assert(mc->gpself != NULL);
        for (size_t ii = 0; ii < mc->du; ii++){
            grad[ii] = 0.0;
        }
    }

    double avg = 0.0;
    double * evals = NULL;
    if (mc->pself > 0.0)
    {
        evals = calloc_double(mc->N+1);
        double * t = calloc_double(mc->N+1);
        double ** x = malloc_dd(mc->N+1);
        x[0] = mcnode_get_xref(mc);
        for (size_t ii = 1; ii < mc->N+1; ii++)
        {
            x[ii] = mcnode_get_xref(mc->neighbors[ii-1]);
            /* if (fabs(x[ii][2]) > M_PI){ */
            /*     mcnode_print(mc,stdout,3); */
            /* } */
        }

        f(mc->N+1,t,x,evals,arg);
        avg = cblas_ddot(mc->N,mc->p,1,evals+1,1);
        avg += mc->pself*evals[0];
        if (grad != NULL){
            for (size_t ii = 1; ii < mc->N+1; ii++ ){
                cblas_daxpy(mc->du,evals[ii],mc->gp[ii-1],1,grad,1);
            }
            cblas_daxpy(mc->du,evals[0],mc->gpself,1,grad,1);
        }

        free(t); t = NULL;
        free(x); x = NULL;
        free(evals); evals = NULL;
    }
    else{
        evals = calloc_double(mc->N);
        double * t = calloc_double(mc->N);
        double ** x = malloc_dd(mc->N);
        for (size_t ii = 0; ii < mc->N; ii++)
        {
            x[ii] = mcnode_get_xref(mc->neighbors[ii]);
        }
        f(mc->N,t,x,evals,arg);
        avg = cblas_ddot(mc->N,mc->p,1,evals,1);

        if (grad != NULL){
            for (size_t ii = 0; ii < mc->N; ii++ ){
                cblas_daxpy(mc->du,evals[ii],mc->gp[ii],1,grad,1);
            }
        }

        free(t); t = NULL;
        free(x); x = NULL;
        free(evals); evals = NULL;
    }
    return avg;
}


/**********************************************************//**
    Sample a transition to neighbors or oneself

    \param[in]     node     - node
    \param[in]     noise    - sample from [0,1]
    \param[in,out] neighbor - (node->d) sized allocated array
**************************************************************/
void mcnode_sample_neighbor(struct MCNode * node, double noise,
                            double * neighbor)
{

    assert (node != NULL);
    assert (noise >= 0.0);
    assert (noise <= 1.0);
    
    double * tprobs = calloc_double(node->N + 1);
    tprobs[0] = node->pself;
    if (node->N > 0){
        memmove(tprobs+1,node->p, node->N*sizeof(double));
    }

    size_t ind = c3sc_sample_discrete_rv(node->N+1,tprobs,noise);
    if (ind == 0){
        memmove(neighbor,node->x,node->d*sizeof(double));
    }
    else{
        //printf("whats up\n");
        double * nx = mcnode_get_xref(node->neighbors[ind-1]);
        //dprint(node->d,nx);
        memmove(neighbor,nx,node->d*sizeof(double));
    }

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
    
    fprintf(fp,"Number of neighbors are %zu\n",n->N);
    for (size_t ii = 0; ii < n->N; ii++){
        fprintf(fp,"\t");
        mcnode_print(n->neighbors[ii],fp,prec);
    }
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
**************************************************************/
struct MCNode *
mca_inbound_node(struct MCA * mca, double time, double * x,
                 double * u, double *dt,double * gdt, double *grad,
                 struct BoundInfo * bi)
{

    assert(mca->dyn != NULL);
    double * probs = calloc_double(2*mca->d);
    double ** gprob = NULL;
    double * gpself = NULL;
    double * dQ = NULL;

    int res;
    size_t du = dyn_get_du(mca->dyn);;
    if (grad != NULL){
        gprob = malloc_dd(2*mca->d);
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       mca->gdrift,mca->diff,mca->gdiff);
        dQ = calloc_double(du);

    }
    else{
        res = dyn_eval(mca->dyn,time,x,u,mca->drift,
                       NULL,mca->diff,NULL);
    }

    assert(res == 0);

    double diff,t;
    double Qh = 0.0;
    double minh2 = pow(mca->minh,2);
    double hrel;
    for (size_t ii = 0; ii < mca->d; ii++){

        hrel = (minh2/mca->h[ii]);
        //printf("hrel = %G\n",hrel);
        assert (mca->dw == mca->d);
        diff = pow(mca->diff[ii*mca->d+ii],2);

        // general case (also needs further modifications)
        //diff = cblas_ddot(mca->dw,
        //                  mca->diff+ii,mca->d,
        //                  mca->diff+ii,mca->d);

        t = pow(mca->minh/mca->h[ii],2);
        probs[ii*2] = t * diff/2.0;
        probs[ii*2+1] = t * diff/2.0;

        if (gprob != NULL){
            gprob[2*ii] = calloc_double(du);
            gprob[2*ii+1] = calloc_double(du);
            //printf("blah t = %G ii =%zu\n",t,ii);
            //dprint2d_col(2,2,mca->gdiff);
            cblas_daxpy(du,t*2.0,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,dQ,1);

            cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                        gprob[2*ii],1);
            cblas_daxpy(du,t,mca->gdiff+ii*mca->d+ii,mca->d*mca->dw,
                        gprob[2*ii+1],1);
            //dprint(mca->du,gprob[2*ii]);
//            printf("grad after diff = "); dprint(mca->du,gprob[2*ii]);
        }
        
        Qh += t*diff;
        
        if (mca->drift[ii] > 0){
            probs[ii*2+1] += hrel * mca->drift[ii];
            Qh += probs[ii*2+1];
            if (gprob != NULL){
                cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,gprob[2*ii+1],1);
                cblas_daxpy(du,hrel,mca->gdrift+ii,mca->d,dQ,1);
            }

        }
        else{
            probs[2*ii] += hrel * (-mca->drift[ii]);
            if (gprob != NULL){
                cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,gprob[2*ii],1);
                cblas_daxpy(du,-hrel,mca->gdrift+ii,mca->d,dQ,1);

            }
            Qh += probs[ii*2];
        }

    }

    *dt = minh2 / Qh;

    if (gprob != NULL){
        double Qh2 = pow(Qh,2);
        for (size_t jj = 0; jj < du; jj++){
            gdt[jj] = -minh2 * dQ[jj] / Qh2;
        }
        for (size_t ii = 0;ii < mca->d; ii++){
            for (size_t jj = 0; jj < du; jj++){
                gprob[2*ii][jj] = (Qh * gprob[2*ii][jj] - dQ[jj]*probs[ii*2])/Qh2;
                gprob[2*ii+1][jj]=(Qh * gprob[2*ii+1][jj] - dQ[jj]*probs[ii*2+1])/Qh2;
            }
        }

    }

    double pself = 1.0;
    for (size_t ii = 0; ii < mca->d ; ii++){
        probs[ii*2] /= Qh;
        probs[ii*2+1] /= Qh;
        pself -= probs[2*ii];
        pself -= probs[2*ii+1];

    }

    
    struct MCNode * mcn = mcnode_init(mca->d,x);
    mcnode_add_neighbors_hspace(mcn,mca->h,pself,probs,bi);

    if (gprob!= NULL){
        gpself = calloc_double(du);
        for (size_t ii = 0; ii < mca->d; ii++){
            cblas_daxpy(du,-1.0,gprob[2*ii],1,gpself,1);
            cblas_daxpy(du,-1.0,gprob[2*ii+1],1,gpself,1);
        }
        mcnode_add_gradients(mcn,mca->d*2,mca->du,gpself,gprob);
        free_dd(2*mca->d,gprob); gprob = NULL;
        free(gpself); gpself = NULL;
        free(dQ); dQ = NULL;
    }


    free(probs); probs = NULL;
    
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

    //printf("get node\n");
    int onbound = bound_info_onbound(*bi);
    if (onbound == 0){
        mcn = mca_inbound_node(mca,time,x,u,dt,gdt,grad,*bi);
    }
    else{
        int absorb = bound_info_absorb(*bi);
        if (absorb == 1){
            *dt = 0;
            mcn = mca_outbound_node(mca,time,x);            
        }
        else{
            int period = bound_info_period(*bi);
            if (period != 0){
                mcn = mca_inbound_node(mca,time,x,u,dt,gdt,grad,*bi);
            }
            else{
                fprintf(stderr,"Not sure what to do with this node\n");
                printf("onbound = %d\n",onbound);
                dprint(mca->d,x);
                exit(1);                
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
                void (*f)(size_t,double *,double **,double*,void*),
                void * arg, struct BoundInfo ** bi, double * grad)
{

    struct MCNode * mcn = mca_get_node(mca,time,x,u,dt,gdt,bi,grad);
    double val = mcnode_expectation(mcn,f,arg,grad);
    mcnode_free(mcn); mcn = NULL;
    return val;
}

/**********************************************************//**
    Sample a transition in the markov chain
    
    \param[in]     mca   - markov chain approximation 
    \param[in]     time  - time
    \param[in]     x     - current state
    \param[in]     u     - current control
    \param[in]     noise - [0,1] uniform sample
    \param[in,out] dt    - output transition time
    \param[in,out] newx  - space for new state
**************************************************************/
void mca_step(struct MCA * mca, double time, double * x, double * u,
         double noise, double * dt, double * newx)
{

    struct MCNode * node = mca_inbound_node(mca,time,x,u,dt,NULL,NULL,NULL);
    
    mcnode_sample_neighbor(node,noise,newx);
    mcnode_free(node); node = NULL;
}

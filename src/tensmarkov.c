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
struct MCNode * mcnode_init(size_t d, double * x)
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
    Get a reference to the node state
**************************************************************/
double * mcnode_getx_ref(struct MCNode * mcn)
{
    return mcn->x;
}

/**********************************************************//**
    Add neighbors to a node that are a distance of 
    *h* away in each coordinate direction

    \param[in,out] mcn   - node
    \param[in]     h     - spacing of nodes in each direction
    \param[in]     pself - new probability of self transition
    \param[in]     probs - transition probs to neighbors
**************************************************************/
void mcnode_add_neighbors_hspace(struct MCNode * mcn, 
                                 double * h, double pself,
                                 double * probs)
{
    assert (mcn != NULL);
    assert (mcn->N == 0);
    mcn->N = 2*mcn->d;
    mcn->neighbors = mcnode_alloc_array(2*mcn->d);
    for (size_t ii = 0; ii < mcn->d; ii++ ){
        mcn->neighbors[2*ii] = mcnode_init(mcn->d,mcn->x);
        mcn->neighbors[2*ii+1] = mcnode_init(mcn->d,mcn->x);
        mcn->neighbors[2*ii]->x[ii] = mcn->x[ii]-h[ii];
        mcn->neighbors[2*ii+1]->x[ii] = mcn->x[ii]+h[ii];
    }

    mcn->pself = pself;
    mcn->p = calloc_double(2*mcn->d);
    memmove(mcn->p,probs,2*mcn->d * sizeof(double));
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
                          size_t du, double *gpself,
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
    struct MCNode * mc,
    void (*f)(size_t,double*,double**x,double*,void*),
    void * arg,double * grad)
{

    if (grad!= NULL){
        assert (mc->gp != NULL);
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
        x[0] = mcnode_getx_ref(mc);
        for (size_t ii = 1; ii < mc->N+1; ii++)
        {
            x[ii] = mcnode_getx_ref(mc->neighbors[ii-1]);
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
            x[ii] = mcnode_getx_ref(mc->neighbors[ii]);
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
        double * nx = mcnode_getx_ref(node->neighbors[ind-1]);
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
    Determine type of node for a given state and tie
**************************************************************/
enum NodeType
mca_node_type(struct MCA * mca, double time, double * x)
{
    assert (mca->bound != NULL);
    assert (mca->dyn   != NULL);
    
    int location = boundary_type(mca->bound,time,x);
    enum NodeType ntype;
    if (location == 0){
        ntype = INBOUNDS;
    }
    else if (location == 1){
        ntype = OUTBOUNDS;
    }
    else if (location == -1){
        ntype = ONBOUNDS;
    }
    else{
        fprintf(stderr,"Cannot determine node type\n");
        fprintf(stderr,"Check to make sure boundary checking function\n");
        fprintf(stderr,"Checks for all possibilities.\nNode is \n");
        for (size_t ii = 0; ii < mca->d; ii++)
        {
            fprintf(stderr,"%G ",x[ii]);
        }
        fprintf(stderr,"\ntime = %G\n",time);
        exit(1);
    }
    return ntype;
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
                 double * u, double *dt,double * gdt, double *grad)
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
    mcnode_add_neighbors_hspace(mcn,mca->h,pself,probs);
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
    \param[in,out] ntype  - node type
    \param[in]     grad   - flag for specifying whether or not to
                            compute gradients of transitions with 
                            respect to control (NULL for no)
**************************************************************/
struct MCNode *
mca_get_node(struct MCA * mca, double time, double * x,
             double * u, double * dt, double * gdt,
             enum NodeType * ntype,
             double *grad)
{
    *ntype = mca_node_type(mca,time,x);
    struct MCNode * mcn = NULL;
    if (*ntype == INBOUNDS){
        mcn = mca_inbound_node(mca,time,x,u,dt,gdt,grad);
    }
    else if (*ntype == OUTBOUNDS){
        *dt = 0;
        mcn = mca_outbound_node(mca,time,x);
    }
    else if (*ntype == ONBOUNDS){
        fprintf(stderr,"Not sure what to do with ONBOUND nodes yet\n");
        exit(1);
    }

    return mcn;
}

/**********************************************************//**
    Compute expectation of a function
    
    \param[in]     mca    - markov chain approximation 
    \param[in]     time   - time
    \param[in]     x      - current state
    \param[in]     u      - current control
    \param[in,out] dt     - time step
    \param[in]     f      - function(npts,times,states,out,arg)
    \param[in,out] arg    - function arguments
    \param[in,out] absorb - absorption
    \param[in,out] grad   - gradient
**************************************************************/
double
mca_expectation(struct MCA * mca, double time,
                double * x, double * u, double * dt, double * gdt,
                void (*f)(size_t,double *,double **,double*,void*),
                void * arg,int * absorb, double * grad)
{
    enum NodeType ntype;

    struct MCNode * mcn = mca_get_node(mca,time,x,u,dt,gdt,&ntype,grad);

    /* if (grad){ */
    /*     for (size_t ii = 0; ii < mca->d; ii++){ */
    /*         printf("minus grad "); */
    /*         dprint(1,mcn->gp[2*ii]); */
    /*         printf("plus grad "); */
    /*         dprint(1,mcn->gp[2*ii+1]); */
    /*     } */
    /*     printf("grad self = "); */
    /*     dprint(1,mcn->gpself); */
    /*     printf("probs = "); */
    /*     dprint(2*mcn->d,mcn->p); */
    /*     printf("self prob = %G\n",mcn->pself); */
    /* } */
    double val = 0.0;
    if (ntype == INBOUNDS){
        val = mcnode_expectation(mcn,f,arg,grad);
        *absorb = 0;
    }
    else if (ntype == OUTBOUNDS){
        *absorb = 1;
    }
    else if (ntype == ONBOUNDS){
        *absorb = 1;
        fprintf(stderr,"Not sure what to do with ONBOUND nodes yet\n");
        exit(1);
    }
    else{
        fprintf(stderr,"Node type %d is unknown\n",ntype);
        exit(1);
    }
    
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
void
mca_step(struct MCA * mca, double time, double * x, double * u,
         double noise, double * dt, double * newx)
{

    struct MCNode * node = mca_inbound_node(mca,time,x,u,dt,NULL,NULL);
    
    //mcnode_print(node,stdout,5);
    //fprintf(stdout,"\n");
    
    //printf("begin sample\n");
    mcnode_sample_neighbor(node,noise,newx);
    //printf("ended sample\n");
    //dprint(node->d,newx);
    mcnode_free(node); node = NULL;
    //printf("freed node out out \n");
}


// assumes only transitions to neighbors
/* double tensor_mm_cost(struct TensorMM * mm, */
/*                       double t, double * x, */
/*                       double * uu, */
/*                       struct Cost * cost, */
/*                       double * dt, */
/*                       double * probs) */
/* { */
    
/*     int res = tensor_mm_tprob(mm,t,x,uu,probs,dt); */
/*     assert (res == 0); */
/*     double evals[2]; */
/*     double pt[2]; */
/*     double out = 0.0; */

/*     //evaluate transition to oneself */
/*     res = cost_eval(cost,t,x,&out); */
/*     if (res != 0){ */
/*         fprintf(stderr, "point x is out of bounds\n \t x = "); */
/*         for (size_t jj = 0; jj < mm->d; jj++){ */
/*             fprintf(stderr,"%G ", x[jj]); */
/*         } */
/*         fprintf(stderr,"\n"); */
/*         exit(1); */
/*     } */
/* //    assert (res == 0); */

/*     //printf("probs are "); */
/*     //dprint(2*mm->d+1,probs); */
/*     out *= probs[2*mm->d]; */

/*     for (size_t ii = 0; ii < mm->d; ii++){ */
/*         pt[0] = x[ii]-mm->h; */
/*         pt[1] = x[ii]+mm->h; */
/*         res = cost_eval_neigh(cost,t,x,ii,pt,evals); */
/*         if (res != 0){ */
/*             fprintf(stderr, "point x is out of bounds\n\t x="); */
/*             for (size_t jj = 0; jj < mm->d; jj++){ */
/*                 fprintf(stderr,"%G ", x[ii]); */
/*             } */
/*             fprintf(stderr,"\n"); */
/*             exit(1); */
/*         } */
/*         //      assert (res == 0); */
/*         out += probs[ii*mm->d]*evals[0]; */
/*         out += probs[ii*mm->d+1]*evals[1]; */
/*     } */

/*     return out; */
/* } */

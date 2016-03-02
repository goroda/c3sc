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
 */
struct MCNode
{
    size_t d;
    double * x;
    
    size_t N;
    struct MCNode ** neighbors;
    double pself;
    double * p; // transition probabilities
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
        //printf("free node\n");
        free(mcn->x); mcn->x = NULL;
        free(mcn->p); mcn->p = NULL;
        //printf("free array\n");
        if (mcn->N > 0){
            //printf("free nodes neighbors N = %zu\n",mcn->N);
            mcnode_free_array(mcn->neighbors,mcn->N);
            free(mcn->neighbors); mcn->neighbors = NULL;
            //printf("freed array\n");
        }
        //printf("freed neighbors\n");
        free(mcn); mcn = NULL;
        //printf("done freeing node\n");
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
    Sample a transition to neighbors or oneself

    \param[in]     node     - node
    \param[in]     noise    - sample from [0,1]
    \param[in,out] neighbor - (node->d) sized allocated array
**************************************************************/
void mcnode_sample_neighbor(struct MCNode * node, double noise, double * neighbor)
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
    size_t dw; 
    double * h;
    double minh;

    struct Dyn * dyn;
    struct Boundary * bound;

    double * drift;
    double * diff;
    
};

/**********************************************************//**
    Allocate MCA

    \param[in] d  - dimension of state space
    \param[in] dw - dimension of noise
    \param[in] h  - spacing of nodes in each direction

    \return mca
**************************************************************/
struct MCA * mca_alloc(size_t d, size_t dw, double * h)
{

    struct MCA * mca = malloc(sizeof(struct MCA));
    assert (mca != NULL);
    
    mca->d = d;
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
    mca->diff = calloc_double(d*dw);
    
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
    Add a reference to Boundary to the MCA 
**************************************************************/
void mca_attach_bound(struct MCA * mca, struct Boundary * bound)
{
    assert (mca != NULL);
    mca->bound = bound;
}

/**********************************************************//**
    Determine type of node for a given state and tie
**************************************************************/
int mca_node_type(struct MCA * mca, double time, double * x)
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
    Generate a node and its transitions
    
    \param[in]     mca  - markov chain approximation 
    \param[in]     time - time
    \param[in]     x    - current state
    \param[in]     u    - current control
    \param[in,out] dt   - output transition time
    
    \return node
**************************************************************/
struct MCNode *
mca_inbound_node(struct MCA * mca, double time, double * x, double * u, double *dt)
{

    assert(mca->dyn != NULL);
    int res = dyn_eval(mca->dyn,time,x,u,
                       mca->drift, mca->diff);
    assert(res == 0);
    
    double * probs = calloc_double(2*mca->d);
    double diff,t;
    double Qh = 0.0;
    for (size_t ii = 0; ii < mca->d; ii++){
        diff = cblas_ddot(mca->dw,
                          mca->diff+ii,mca->d,
                          mca->diff+ii,mca->d);

        t = pow(mca->minh/mca->h[ii],2);
        probs[ii*2] = t * diff/2.0;
        probs[ii*2+1] = t * diff/2.0;
        Qh += t*diff;
        if (mca->drift[ii] > 0){
            probs[ii*2+1] += ((pow(mca->minh,2)/mca->h[ii]) * mca->drift[ii]);
            Qh += probs[ii*2+1];
        }
        else{
            probs[ii*2] += ((pow(mca->minh,2)/mca->h[ii]) *(-mca->drift[ii]));
            Qh += probs[ii*2];
        }
    }

    *dt = pow(mca->minh,2) / Qh;
    
    double pself = 1.0;
    for (size_t ii = 0; ii < mca->d ; ii++){
        probs[ii*2] /= Qh;
        probs[ii*2+1] /= Qh;
        pself -= probs[ii*2];
        pself -= probs[ii*2+1];
    }

    struct MCNode * mcn = mcnode_init(mca->d,x);
    mcnode_add_neighbors_hspace(mcn,mca->h,pself,probs);

    free(probs); probs = NULL;
    
    return mcn;
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
mca_step(struct MCA * mca, double time, double * x, double * u,double noise,
         double * dt, double * newx)
{

    struct MCNode * node = mca_inbound_node(mca,time,x,u,dt);
    
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

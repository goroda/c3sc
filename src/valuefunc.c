#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "c3.h"
#include "util.h"

/** \struct ValueF
\brief Value Function
\var ValueF::d        dimension of state space
\var ValueF::cost     function train cost
\var ValueF::N        size of grid
\var ValueF::cores    core evaluations at grid
\var ValueF::grid     grid of values
*/
struct ValueF
{
    size_t d;
    struct FunctionTrain * cost;
    size_t * N;


    // private memory allocations for efficient doing of things
    double ** cores;
    double ** center;
    double ** fprod;
    double ** bprod;
    double * workspace;

};

struct ValueF * valuef_create(size_t d)
{
    struct ValueF * vf = malloc(sizeof(struct ValueF));
    assert (vf != NULL);

    vf->d = d;
    vf->N = calloc_size_t(d);

    vf->cores     = NULL;
    vf->center    = NULL;
    vf->fprod     = NULL;
    vf->bprod     = NULL;
    vf->workspace = NULL;

    return vf;
}

/**********************************************************//**
    Free memory allocated to the value function

    \param[in,out] vf - structure whose memory to free 
**************************************************************/
void valuef_destroy(struct ValueF * vf)
{
    if (vf != NULL){
        function_train_free(vf->cost); vf->cost = NULL;
        free(vf->N); vf->N = NULL;

        free_dd(vf->d,vf->cores); vf->cores = NULL;

        /* freeing private stuff */
        /* free_dd(vf->d,vf->center); */
        /* free(vf->center); vf->center = NULL; */

        free_dd(vf->d,vf->fprod); vf->fprod = NULL;
        /* free(vf->fprod); vf->fprod = NULL; */

        free_dd(vf->d,vf->bprod); vf->bprod = NULL;
        /* free(vf->bprod); vf->bprod = NULL; */

        free(vf->workspace); vf->workspace = NULL;

        free(vf); vf = NULL;
    }
}

void valuef_precompute_cores(struct ValueF * cost)
{
    assert (cost != NULL);
    assert (cost->cost != NULL);
    size_t nrows, ncols, nvals;
    size_t * ranks = function_train_get_ranks(cost->cost);
    if (cost->cores == NULL){
        cost->cores = malloc_dd(cost->d);
        for (size_t jj = 0; jj < cost->d; jj++){
            cost->cores[jj] = NULL;
        }
    }
    for (size_t ii = 0; ii < cost->d; ii++){
        free(cost->cores[ii]); cost->cores[ii] = NULL;
        nrows = ranks[ii];
        ncols = ranks[ii+1];
        nvals = cost->N[ii];
        struct GenericFunction ** arr = cost->cost->cores[ii]->funcs;
        cost->cores[ii] = calloc_double(nrows * ncols * nvals);
        for (size_t jj = 0; jj < nvals; jj++){
            generic_function_1darray_eval2_ind(nrows * ncols,arr,jj,
                                               cost->cores[ii] + jj * nrows * ncols);
        }
    }
}

/**********************************************************//**
    Copy a value function                                                           
**************************************************************/
struct ValueF * valuef_copy(struct ValueF * vf)
{
    assert (vf != NULL);

    struct ValueF * vf2 = valuef_create(vf->d);
    vf2->cost = function_train_copy(vf->cost);
    memmove(vf2->N,vf->N,vf->d * sizeof(size_t));
    valuef_precompute_cores(vf2);
    return vf2;
}

/**********************************************************//**
    Get the value function ranks
**************************************************************/
size_t * valuef_get_ranks(struct ValueF * vf)
{
    assert (vf != NULL);
    assert (vf->cost != NULL);
    return function_train_get_ranks(vf->cost);
}

/**********************************************************//**
    Get the norm of the value function

    \param[in] vf - value function
    
    \return L2 Norm
**************************************************************/
double valuef_norm(struct ValueF * vf)
{
    assert (vf != NULL);
    assert (vf->cost != NULL);
    return function_train_norm2(vf->cost);
}

/**********************************************************//**
    Get the norm between value functions

    \param[in] vf  - value function
    \param[in] vf2 - second value function
    
    \return L2 difference 
**************************************************************/
double valuef_norm2diff(struct ValueF * vf, struct ValueF * vf2)
{
    assert (vf != NULL);
    assert (vf2 != NULL);
    return function_train_norm2diff(vf->cost,vf2->cost);
}

/**********************************************************//**
    Evaluate the value function

    \param[in] vf - value function to evaluate 
    \param[in] x  - location at which to Evaluate
    
    \return evaluation
**************************************************************/
double valuef_eval(struct ValueF * vf, const double * x)
{
    assert (vf != NULL);
    assert (vf->cost != NULL);
    return function_train_eval(vf->cost,x);
}



/**********************************************************//**
   Evaluate the cost of a fiber and associated neighbors

   \param[in]     vf             - value function
   \param[in]     fixed_ind      - (d,) indices of fixed points
   \param[in]     dim_vary       - dimension along which fiber exists
   \param[in]     neighbors      - (2*(d-1),) array of neighbors (2 in each dim that is not varying)
   \param[in]     neighbors_vary - (2*N,) array of neighbors in dim_vary direction
   \param[in,out] out            - (N*(2d+1),) array of values at neighbors of each element including fibers

   \returns  0 if successful
   
   \note
   out ordered as an array of N arrays
   each arary is (2d+1) ordered by cost in (- +) neighbor in each direction
   fiber followed by the value at the fiber
**************************************************************/
int valuef_eval_fiber_ind_nn(struct ValueF * vf, const size_t * fixed_ind,
                              size_t dim_vary, const size_t * neighbors,
                              const size_t * neighbors_vary, double * out)
{
    assert (vf != NULL);
    assert (vf->cost != NULL);
    assert (vf->cores != NULL);
    size_t maxrank = function_train_get_maxrank(vf->cost);
    size_t dim = vf->d;

    /* if (vf->center == NULL){ */
    /*     vf->center = malloc_dd(dim); */
    /*     for (size_t ii = 0; ii < dim; ii++){ */
    /*         vf->center[ii] = calloc_double(maxrank); */
    /*     } */
    /* } */
    if (vf->fprod == NULL){
        assert (vf->bprod == NULL); // both should be null

        size_t nmax = vf->N[0];
        for (size_t ii = 0; ii < dim; ii++){
            if (vf->N[ii] > nmax){
                nmax = vf->N[ii];
            }
        }

        vf->fprod = malloc_dd(dim);
        vf->bprod = malloc_dd(dim);
        for (size_t ii = 0; ii < dim; ii++){
            vf->fprod[ii] = calloc_double(maxrank*nmax);
            vf->bprod[ii] = calloc_double(maxrank*nmax);
        }
    }
    if (vf->workspace == NULL){
        vf->workspace = calloc_double(maxrank);
    }

    double * space1 = vf->workspace;
    /* double ** center = vf->center; */
    double ** bprod = vf->bprod;
    double ** fprod = vf->fprod;

    // get the cores for front to back until dim_vary
    size_t nrows, ncols, nvals;
    size_t * ranks = function_train_get_ranks(vf->cost);
    for (size_t ii = 0; ii < dim_vary; ii++){
        /* printf("not here\n"); */
        nrows = ranks[ii];
        ncols = ranks[ii+1];

        double * center = vf->cores[ii] + fixed_ind[ii]*nrows*ncols;
        if (ii == 0){
            memmove(fprod[ii],center,nrows*ncols*sizeof(double));
        }
        else{
            cblas_dgemv(CblasColMajor,CblasTrans, nrows, ncols,
                        1.0, center, nrows,
                        fprod[ii-1], 1, 0.0, fprod[ii], 1);
        }
    }
    /* printf("got front to back\n"); */

    // for back to front until dim_vary
    for (size_t ii = vf->d-1; ii > dim_vary; ii--){
        nrows = ranks[ii];
        ncols = ranks[ii+1];
        nvals = vf->N[ii];
        double * center = vf->cores[ii] + fixed_ind[ii]*nrows*ncols;
        if (ii == (vf->d-1)){
            memmove(bprod[ii],center,nrows*ncols*sizeof(double));
        }
        else{
            cblas_dgemv(CblasColMajor,CblasNoTrans,
                        nrows, ncols, 1.0,
                        center, nrows,
                        bprod[ii+1], 1, 0.0, bprod[ii], 1);
        }
    }
    /* printf("got back to front\n"); */

    // dim_vary for forward 
    nrows = ranks[dim_vary];
    ncols = ranks[dim_vary+1];
    nvals = vf->N[dim_vary];
    if (dim_vary == 0){
        for (size_t jj = 0; jj < nvals; jj++){
            memmove(fprod[dim_vary]+jj*ncols,vf->cores[dim_vary] + jj * nrows * ncols,
                    nrows * ncols * sizeof(double));
        }   
    }
    else{
        for (size_t jj = 0; jj < nvals; jj++){
            cblas_dgemv(CblasColMajor,CblasTrans, nrows, ncols,
                        1.0, vf->cores[dim_vary] + jj * nrows * ncols, 
                        nrows, fprod[dim_vary-1], 1, 0.0, fprod[dim_vary]+jj*ncols, 1);
        }
    }

    // for backward
    if (dim_vary == (vf->d-1)){
        for (size_t jj = 0; jj < nvals; jj++){
            memmove(bprod[dim_vary]+jj*nrows,vf->cores[dim_vary] + jj * nrows * ncols,
                    nrows * ncols * sizeof(double));
        }
    }
    else{
        for (size_t jj = 0; jj < nvals; jj++){
            cblas_dgemv(CblasColMajor,CblasNoTrans, nrows, ncols,
                        1.0, vf->cores[dim_vary] + jj * nrows * ncols, 
                        nrows, bprod[dim_vary+1], 1, 0.0, bprod[dim_vary]+jj*nrows, 1);
        }
    }
    
    /* printf("got dim_vary\n"); */

    // after dim_vary for front-to back
    for (size_t ii = dim_vary+1; ii < vf->d; ii++){
        nrows = ranks[ii];
        ncols = ranks[ii+1];

        double * center = vf->cores[ii] + fixed_ind[ii]*nrows*ncols;
        for (size_t jj = 0; jj < nvals; jj++){
            cblas_dgemv(CblasColMajor,CblasTrans, nrows, ncols,
                        1.0, center, nrows,
                        fprod[ii-1]+jj*nrows, 1, 0.0, fprod[ii]+jj*ncols, 1);            
        }
    }

    // after dim_vary for back to front 
    for (int ii = dim_vary-1; ii >= 0;  ii--){
        nrows = ranks[ii];
        ncols = ranks[ii+1];

        double * center = vf->cores[ii] + fixed_ind[ii]*nrows*ncols;
        for (size_t jj = 0; jj < nvals; jj++){
            cblas_dgemv(CblasColMajor,CblasNoTrans, nrows, ncols,
                        1.0, center, nrows,
                        bprod[ii+1]+jj*ncols, 1, 0.0, bprod[ii]+jj*nrows, 1);            
        }

    }

    /* /\* printf("precomputed everything\n"); *\/ */
    /* memmove(out_fiber,bprod[0],nvals*sizeof(double)); */

    size_t stride = 2 * dim + 1;
    for (size_t ii = 0; ii < nvals; ii++){
        out[ii*stride + 2*dim_vary] = bprod[0][neighbors_vary[2*ii]];
        out[ii*stride + 2*dim_vary+1] = bprod[0][neighbors_vary[2*ii+1]];
        out[ii*stride + 2*dim] = bprod[0][ii];
    }

    // now have everything and only need to assemble
    for (size_t ii = 0; ii < dim_vary; ii++){
        nrows = ranks[ii];
        ncols = ranks[ii+1];
        
        if (ii == 0){
            for (size_t jj = 0; jj < nvals; jj++){
                out[jj*stride + 2*ii] = cblas_ddot(ncols,vf->cores[ii]+neighbors[2*ii]*ncols,1,
                                                   bprod[ii+1]+jj*ncols,1);
                out[jj*stride + 2*ii+1] = cblas_ddot(ncols,vf->cores[ii]+neighbors[2*ii+1]*ncols,1,
                                                     bprod[ii+1]+jj*ncols,1);
            }
        }
        else{
            for (size_t jj = 0; jj < nvals; jj++){
                cblas_dgemv(CblasColMajor,CblasNoTrans, nrows, ncols,
                            1.0, vf->cores[ii] + neighbors[2*ii]*nrows*ncols, nrows,
                            bprod[ii+1]+jj*ncols, 1, 0.0, space1, 1);
                out[jj*stride + 2*ii] = cblas_ddot(nrows,space1,1,fprod[ii-1],1);

                cblas_dgemv(CblasColMajor,CblasNoTrans, nrows, ncols,
                            1.0, vf->cores[ii] + neighbors[2*ii+1]*nrows*ncols, nrows,
                            bprod[ii+1]+jj*ncols, 1, 0.0, space1, 1);
                out[jj*stride + 2*ii+1] = cblas_ddot(nrows,space1,1,fprod[ii-1],1);
            }
        }
    }
    
    for (size_t ii = dim_vary+1; ii < dim; ii++){
        /* printf("ii = %zu\n",ii); */
        nrows = ranks[ii];
        ncols = ranks[ii+1];

        if (ii == dim-1){
            /* printf("here!\n"); */
            for (size_t jj = 0; jj < nvals; jj++){
                /* printf("jj=%zu, neighbor=%zu\n",jj,neighbors[2*(ii-1)]); */
                out[jj*stride + 2*ii] = cblas_ddot(nrows,vf->cores[ii]+neighbors[2*(ii-1)]*nrows,1,
                                                   fprod[ii-1]+jj*nrows,1);
                /* printf("\tgot first\n"); */
                out[jj*stride + 2*ii+1] = cblas_ddot(nrows,vf->cores[ii]+neighbors[2*(ii-1)+1]*nrows,1,
                                                     fprod[ii-1]+jj*nrows,1);
            }
        }
        else{
            for (size_t jj = 0; jj < nvals; jj++){
                cblas_dgemv(CblasColMajor,CblasTrans, nrows, ncols,
                            1.0, vf->cores[ii] + neighbors[2*(ii-1)] * nrows * ncols, 
                            nrows, fprod[ii-1] + jj * nrows, 1, 0.0, space1, 1);
                out[jj*stride + 2*ii] = cblas_ddot(ncols,space1,1,bprod[ii+1],1);

                cblas_dgemv(CblasColMajor,CblasTrans, nrows, ncols,
                            1.0, vf->cores[ii] + neighbors[2*(ii-1)+1] * nrows * ncols, 
                            nrows, fprod[ii-1] + jj * nrows, 1, 0.0, space1, 1);
                out[jj*stride + 2*ii+1] = cblas_ddot(ncols,space1,1,bprod[ii+1],1);
            }   
        }
        /* printf("last ind = %zu\n",(nvals-1)*stride+2*ii+1); */
    }

    return 0;
}

/**********************************************************//**
   Create a value function through interpolation

   \param[in] d       - state space dimension
   \param[in] f       - function
   \param[in] args    - function arguments
   \param[in] N       - number of nodes in each dimension of discretization
   \param[in] grid    - discretization nodes
   \param[in] start   - starting nodes for cross 
   \param[in] aargs   - interpolation arguments
   \param[in] verbose - verbosity level

   \return Value function
**************************************************************/
struct ValueF * valuef_interp(size_t d, int (*f)(size_t,const double *,double*,void*),void * args,
                              const size_t * N, double ** grid, double ** start,
                              struct ApproxArgs * aargs, int verbose)
{
    struct ValueF * vf = valuef_create(d);
    memmove(vf->N,N,d * sizeof(size_t));

    struct C3Approx * c3a = c3approx_create(CROSS,d);
    struct LinElemExpAopts ** aopts = malloc(d * sizeof(struct LinElemExpAopts *));
    struct OneApproxOpts ** qmopts = malloc(d * sizeof(struct OneApproxOpts *));

    double cross_tol = approx_args_get_cross_tol(aargs);
    double round_tol = approx_args_get_round_tol(aargs);
    size_t kickrank  = approx_args_get_kickrank(aargs);
    size_t maxrank   = approx_args_get_maxrank(aargs);
    size_t startrank = approx_args_get_startrank(aargs);
    for (size_t ii = 0; ii < d; ii++){
        aopts[ii] = lin_elem_exp_aopts_alloc(N[ii],grid[ii]);
        qmopts[ii] = one_approx_opts_alloc(LINELM,aopts[ii]);
        c3approx_set_approx_opts_dim(c3a,ii,qmopts[ii]);
    }
    c3approx_init_cross(c3a,startrank,verbose,start);
    c3approx_set_cross_tol(c3a,cross_tol);
    c3approx_set_cross_maxiter(c3a,10);
    c3approx_set_round_tol(c3a,round_tol);
    c3approx_set_adapt_kickrank(c3a,kickrank);
    size_t minN = N[0];
    for (size_t ii = 0; ii < d; ii++){
        if (N[ii] < minN){ minN = N[ii];}
    }
    if (maxrank < minN){
        c3approx_set_adapt_maxrank_all(c3a,maxrank);
    }
    else{
        c3approx_set_adapt_maxrank_all(c3a,minN);
    }

    struct Fwrap * fw = fwrap_create(d,"general-vec");
    fwrap_set_fvec(fw,f,args);

    int adapt = approx_args_get_adapt(aargs);;
    vf->cost = c3approx_do_cross(c3a,fw,adapt);
    valuef_precompute_cores(vf);
    /* iprint_sz(4,vf->cost->ranks); */
    

    for (size_t ii = 0; ii < d; ii++){
        one_approx_opts_free_deep(&(qmopts[ii]));
    }
    free(qmopts); qmopts = NULL;
    free(aopts); aopts = NULL;
    c3approx_destroy(c3a); c3a = NULL;
    fwrap_destroy(fw);

    return vf;
}

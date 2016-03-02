#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "c3.h"

#include "util.h"
#include "cost.h"


/**********************************************************//**
    Allocate the cost
**************************************************************/
struct Cost * cost_alloc()
{
    struct Cost * cost = NULL;
    cost = malloc(sizeof(struct Cost));
    if (cost == NULL){
        fprintf(stderr,"Could not allocate cost\n");
        exit(1);
    }
    cost->discrete = 0;
    cost->ndisc = 0;
    cost->bds = NULL;
    cost->cost = NULL;
    return cost;
}


/**********************************************************//**
    Initialize the cost of a discrete cost function

    \param[in,out] cost - allocated cost
    \param[in]     dim  - dimension of cost function
    \param[in]     N    - discretization level (same each dim)
    \param[in]     ft   - discrete cost function
**************************************************************/
void cost_init_discrete(struct Cost * cost,size_t dim,
                        size_t N, struct FunctionTrain * ft)
{
    assert (cost != NULL);
    assert (cost->bds == NULL);
    assert (cost->cost == NULL);
    cost->discrete = 1;
    cost->ndisc = N;
    cost->bds = bounding_box_init_std(dim); // [-1,1] bounds
    cost->cost = function_train_copy(ft);
}


/**********************************************************//**
    Free a cost function
**************************************************************/
void cost_free(struct Cost * c)
{
    if (c != NULL){
        bounding_box_free(c->bds); c->bds = NULL;
        function_train_free(c->cost); c->cost = NULL;
        free(c); c= NULL;
    }
}

/**********************************************************//**
    Set cost function to an approximation of some input
    function

    \param[in,out] c       - cost structure
    \param[in]     f       - function
    \param[in]     args    - function arguments
    \param[in]     verbose - verbosity level for approximation

    \note
    c->bds and c->cost should be NULL
**************************************************************/
void cost_approx(struct Cost * c,
                 double (*f)(double *, void *),
                 void * args, int verbose)
{
    assert (c != NULL);
    assert (c->bds != NULL);
    assert (c->cost == NULL);
    
    size_t n = c->ndisc;
    size_t d = c->bds->dim;
    double lb = -1.0;
    double ub = 1.0;
    
    struct LinElemExpAopts aopts = {n,0};
    struct FtApproxArgs * fapp = NULL;
    struct FiberOptArgs * fopt = NULL;
    fapp = ft_approx_args_create_le(d,&aopts);
    
    double * xnodes = linspace(lb,ub,n);
    struct c3Vector c3v = {n,xnodes};
    fopt = fiber_opt_args_bf_same(d,&c3v);

    struct FtCrossArgs fca;
    ft_cross_args_init(&fca);
    fca.dim = d;
    fca.ranks = calloc_size_t(d+1);
    for (size_t ii = 0; ii < d+1; ii++){
        fca.ranks[ii] = 3;
    }
    fca.ranks[0] = 1;
    fca.ranks[d] = 1;
    fca.epsilon = 1e-7;
    fca.maxiter = 10;
    fca.epsround = 1e-10;
    fca.kickrank = 2;
    fca.maxiteradapt = 5;
    fca.verbose = verbose;
    fca.optargs = fopt;

    double ** start = malloc_dd(d);

    for (size_t ii = 0; ii < d; ii++){
        start[ii] = calloc_double(3);
        start[ii][0] = xnodes[0];
        start[ii][1] = xnodes[4];
        start[ii][2] = xnodes[n-1];
    }

    c->cost = function_train_cross(f,args,c->bds,
                                   start,&fca,fapp);

    ft_approx_args_free(fapp); fapp = NULL;
    free(fca.ranks); fca.ranks = NULL;
    free(xnodes); xnodes = NULL;
    free_dd(d,start);
    fiber_opt_args_free(fopt); fopt = NULL;
}


/**********************************************************//**
    Approximate a new function that takes

    \note 
    This function doesn't make sense here 
**************************************************************/
struct FunctionTrain *
cost_approx_new(struct Cost * c,
                double (*f)(double *, void *),
                void * args, int verbose)
{
    assert (c != NULL);
    assert (c->bds != NULL);
    
    size_t n = c->ndisc;
    size_t d = c->bds->dim;
    double lb = -1.0;
    double ub = 1.0;
    
    struct LinElemExpAopts aopts = {n,0};
    struct FtApproxArgs * fapp =
        ft_approx_args_create_le(d,&aopts);
    
    double * xnodes = linspace(lb,ub,n);
    
    struct c3Vector c3v = {n,xnodes};
    struct FiberOptArgs * fopt =
        fiber_opt_args_bf_same(d,&c3v);

    struct FtCrossArgs fca;
    ft_cross_args_init(&fca);   
    fca.dim = d;
    fca.ranks = calloc_size_t(d+1);
    for (size_t ii = 0; ii < d+1; ii++){
        fca.ranks[ii] = 3;
    }
    fca.ranks[0] = 1;
    fca.ranks[d] = 1;
    fca.epsilon = 1e-10;
    fca.maxiter = 10;
    fca.epsround = 1e-10;
    fca.kickrank = 5;
    fca.maxiteradapt = 5;
    fca.verbose = verbose;
    fca.optargs = fopt;

    double ** start = malloc_dd(d);

    size_t mid = n/2;
    //printf("mid = %zu\n",mid);
    for (size_t ii = 0; ii < d; ii++){
        start[ii] = calloc_double(3);
        start[ii][0] = xnodes[0];
        start[ii][1] = xnodes[mid];
        start[ii][2] = xnodes[n-1];
    }

    struct FunctionTrain * cost = 
        function_train_cross(f,args,c->bds,start,&fca,fapp);

    ft_approx_args_free(fapp); fapp = NULL;
    free(fca.ranks); fca.ranks = NULL;
    free(xnodes); xnodes = NULL;
    free_dd(d,start);
    fiber_opt_args_free(fopt); fopt = NULL;

    return cost;
}


/**********************************************************//**
    Update the Cost function with cost (function train)

    \param[in,out] c    - cost funciton
    \param[in]     cost - new compressed cost

    \note
    This function free's the memory of the old cost function
    and refers to the new one by reference
**************************************************************/
void cost_update_ref(struct Cost * c,
                     struct FunctionTrain * cost)
{
    
    function_train_free(c->cost); c->cost = NULL;
    c->cost = cost;
}

/**********************************************************//**
    Evaluate a cost function

    \param[in] cost - cost funciton
    \param[in] time - time at which to evaluate
    \param[in] x    - location in space at which to evaluate
    \param[in] eval - pointer to evaluation location

    \return res - 0 if everything is ok 
                 !0 if x is out of expected bounds
**************************************************************/
int cost_eval(struct Cost * cost,
              double time,
              double * x,
              double * eval)
{

    (void)(time);
//    printf("cost bounds are\n");
//    dprint(cost->bds->dim, cost->bds->lb);
//    dprint(cost->bds->dim, cost->bds->ub);
//    printf("x is \n");
//    dprint(cost->bds->dim,x);
    int res = c3sc_check_bounds(cost->bds->dim,
                                cost->bds->lb,
                                cost->bds->ub,
                                x);
//    printf("res = %d\n",res);

    if (res != 0){
        return res;
    }

    *eval = function_train_eval(cost->cost,x);

    return 0;
    
}


/**********************************************************//**
    Evaluate a cost function at two neighboring points
    around x

    \param[in]     cost  - cost funciton
    \param[in]     time  - time at which to evaluate
    \param[in]     x     - location at which to evaluate
    \param[in]     ii    - dimension at which to perturb x
    \param[in]     pt    - perturbed values of x[ii] 
    \param[in,out] evals - space allocated for evaluation

    \return res - 0 if everything is ok 
                 !0 if x is out of expected bounds
**************************************************************/
int cost_eval_neigh(struct Cost * cost,
                    double time,
                    double * x,
                    size_t ii,
                    double pt[2],
                    double evals[2])
{
    // need special function train eval;

    // slow for now
    size_t d = cost->bds->dim;
    double * x2 = calloc_double(d);
    memmove(x2,x,d*sizeof(double));
    x2[ii] = pt[0];
    int res = cost_eval(cost,time,x2,evals);
    if (res != 0){
        free(x2); x2 = NULL;
        return res;
    }
    x2[ii] = pt[1];
    res = cost_eval(cost,time,x2,evals+1);
    free(x2); x2 = NULL;
    return res;
    
}

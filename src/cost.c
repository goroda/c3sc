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
struct Cost * cost_alloc(size_t d,double * lb, double * ub)
{
    struct Cost * cost = NULL;
    cost = malloc(sizeof(struct Cost));
    if (cost == NULL){
        fprintf(stderr,"Could not allocate cost\n");
        exit(1);
    }
    cost->d = d;
    cost->bds = bounding_box_vec(d,lb,ub);
    cost->cost = NULL;
    
    cost->N = NULL;
    cost->x = NULL;

    return cost;
}

/**********************************************************//**
    Free a cost function
**************************************************************/
void cost_free(struct Cost * c)
{
    if (c != NULL){
        bounding_box_free(c->bds); c->bds = NULL;
        function_train_free(c->cost); c->cost = NULL;
        free_dd(c->d,c->x); c->x = NULL;
        free(c->N); c->N = NULL;
        free(c); c= NULL;
    }
}

/**********************************************************//**
    Return a reference to cost funciton lower bounds
**************************************************************/
double * cost_get_lb(struct Cost * c)
{
    assert (c != NULL);
    if (c->bds == NULL){
        return NULL;
    }
    return bounding_box_get_lb(c->bds);
}

/**********************************************************//**
    Return a reference to cost funciton upper
**************************************************************/
double * cost_get_ub(struct Cost * c)
{
    assert (c != NULL);
    if (c->bds == NULL){
        return NULL;
    }
    return bounding_box_get_ub(c->bds);
}

/**********************************************************//**
    Get cost function ranks
**************************************************************/
size_t * cost_get_ranks(struct Cost * c)
{
    assert (c != NULL);
    return c->cost->ranks;
}

/**********************************************************//**
    Initialize the cost of a discrete cost function

    \param[in,out] cost - allocated cost
    \param[in]     N    - number of nodes in each dimension
    \param[in]     x    - nodes in each dimension
**************************************************************/
void cost_init_discrete(struct Cost * cost,size_t * N,double ** x)
{
    assert (cost != NULL);
    
    cost->x = malloc_dd(cost->d);
    cost->N = calloc_size_t(cost->d);
    for (size_t ii = 0; ii < cost->d; ii++){
        cost->N[ii] = N[ii];
        cost->x[ii] = calloc_double(N[ii]);
        memmove(cost->x[ii],x[ii],N[ii]*sizeof(double));
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
    assert (c->N != NULL);
    assert (c->x != NULL);

    
    struct LinElemExpAopts ** aopts = NULL;
    aopts = malloc(c->d*sizeof(struct LinElemExpAopts *));
    assert (aopts != NULL);

    for (size_t ii = 0; ii < c->d;ii++){
        aopts[ii] = lin_elem_exp_aopts_alloc(c->N[ii], c->x[ii]);
    }

    struct FtApproxArgs * fapp = NULL;
    struct FiberOptArgs * fopt = NULL;
    fapp = ft_approx_args_create_le2(c->d,aopts); 
    
    struct c3Vector ** c3v = c3vector_alloc_array(c->d);
    for (size_t ii = 0; ii < c->d; ii++){
        c3v[ii] = c3vector_alloc(c->N[ii],c->x[ii]);
    }
    fopt = fiber_opt_args_bf(c->d,c3v);

    struct FtCrossArgs fca;
    ft_cross_args_init(&fca);
    fca.dim = c->d;
    fca.ranks = calloc_size_t(c->d+1);
    for (size_t ii = 0; ii < c->d+1; ii++){
        fca.ranks[ii] = 3;
    }
    fca.ranks[0] = 1;
    fca.ranks[c->d] = 1;
    fca.epsilon = 1e-7;
    fca.maxiter = 10;
    fca.epsround = 1e-10;
    fca.kickrank = 2;
    fca.maxiteradapt = 5;
    fca.verbose = verbose;
    fca.optargs = fopt;

    double ** start = malloc_dd(c->d);

    for (size_t ii = 0; ii < c->d; ii++){
        assert (c->N[ii] > 5);
        size_t mid = c->N[ii]/2;
        //printf("mid = %zu\n",mid);
        start[ii] = calloc_double(3);
        start[ii][0] = c->x[ii][1];
        start[ii][1] = c->x[ii][mid];
        start[ii][2] = c->x[ii][c->N[ii]-2];
        //printf("start is "); dprint(3,start[ii]);
    }

    c->cost = function_train_cross(f,args,c->bds,
                                   start,&fca,fapp);
    
    for (size_t ii = 0; ii < c->d; ii++){
        lin_elem_exp_aopts_free(aopts[ii]); aopts[ii] = NULL;
    }
    free(aopts); aopts = NULL;
    c3vector_free_array(c3v,c->d); c3v = NULL;
    ft_approx_args_free(fapp); fapp = NULL;
    free(fca.ranks); fca.ranks = NULL;
    free_dd(c->d,start);
    fiber_opt_args_free(fopt); fopt = NULL;
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
        printf("point is not in bounds for cost evaluation \n");
        dprint(cost->bds->dim,x);
        return res;
    }

    *eval = function_train_eval(cost->cost,x);

    return 0;
    
}

/**********************************************************//**
    Black box (bb) interface to cost function evaluation

    \param[in]     N    - number of points at which to evaluate
    \param[in]     t    - times of points
    \param[in]     x    - locations of points
    \param[in,out] out  - allocated evaluation space
    \param[in]     args - pointer to cost function structure
**************************************************************/
void cost_eval_bb(size_t N, double * t, double ** x,
                  double * out, void * args)
{
    struct Cost * c = args;
    int res;
    for (size_t ii = 0; ii < N; ii++)
    {
        res = cost_eval(c,t[ii],x[ii],out+ii);
        assert (res == 0);
    }
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

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
    Copy a cost function
**************************************************************/
struct Cost * cost_copy_deep(struct Cost * old)
{
    if (old == NULL){
        return NULL;
    }

    struct Cost * newc = cost_alloc(old->d,old->bds->lb,old->bds->ub);
    if (old->N != NULL){
        newc->N = calloc_size_t(newc->d);
        memmove(newc->N,old->N,newc->d*sizeof(size_t));
    }
    if (old->x != NULL){
        newc->x = malloc_dd(newc->d);
        for (size_t ii = 0; ii < newc->d; ii++){
            newc->x[ii] = calloc_double(newc->N[ii]);
            memmove(newc->x[ii],old->x[ii],newc->N[ii]*sizeof(double));
        }
    }
    
    return newc;
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

    \param[in,out] c         - cost structure
    \param[in]     f         - function
    \param[in]     args      - function arguments
    \param[in]     verbose   - verbosity level for approximation
    \param[in]     cross_tol - cross approximation tolerance
    \param[in]     round_tol - rounding tolerance
    \param[in]     kickrank  - rank increase level
    \note
    c->cost should be NULL
**************************************************************/
void cost_approx(struct Cost * c,
                 double (*f)(double *, void *),
                 void * args, int verbose,double cross_tol,
                 double round_tol, size_t kickrank)
{
    assert (c != NULL);
    assert (c->bds != NULL);
    assert (c->cost == NULL);
    assert (c->N != NULL);
    assert (c->x != NULL);

    struct C3Approx * c3a = c3approx_create(CROSS,c->d,c->bds->lb,c->bds->ub);
    c3approx_init_lin_elem(c3a);
    c3approx_set_lin_elem_fixed(c3a,c->N,c->x);

    size_t init_rank = 5;
    c3approx_init_cross(c3a,init_rank,verbose);
    c3approx_set_fiber_opt_brute_force(c3a,c->N,c->x);
    
    c3approx_set_cross_tol(c3a,cross_tol);
    c3approx_set_cross_maxiter(c3a,10);
    c3approx_set_round_tol(c3a,round_tol);
    c3approx_set_adapt_kickrank(c3a,kickrank);
    size_t minN = c->N[0];
    for (size_t ii = 0; ii < c->d; ii++){
        if (c->N[ii] < minN){ minN = c->N[ii];}
    }
    c3approx_set_adapt_maxrank_all(c3a,minN);

    // determine starting nodes
//    printf("x[0] = "); dprint(c->N[0],c->x[0]);
    double ** start = malloc_dd(c->d);
    size_t * nstart = calloc_size_t(c->d);
    for (size_t ii = 0; ii < c->d; ii++){
        nstart[ii] = 5;
        assert (c->N[ii] > 5);
        size_t mid = c->N[ii]/2;
        start[ii] = calloc_double(nstart[ii]);
        start[ii][3] = c->x[ii][0];
        start[ii][1] = c->x[ii][1];
        start[ii][0] = c->x[ii][mid];
        start[ii][2] = c->x[ii][c->N[ii]-2];
        start[ii][4] = c->x[ii][c->N[ii]-1];
        //printf("start mid = %G\n",start[ii][0]); 
    }
    c3approx_set_start(c3a,nstart,start);

    c->cost = c3approx_do_cross(c3a,f,args);
    
    free(nstart); nstart = NULL;
    free_dd(c->d,start); start = NULL;
    c3approx_destroy(c3a); c3a = NULL;
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

    double * xuse = NULL;
    if (res != 0){
        xuse = calloc_double(cost->bds->dim);
        for (size_t ii = 0; ii < cost->bds->dim; ii++){
            if (x[ii] < cost->bds->lb[ii]){
                xuse[ii] = cost->bds->lb[ii];
            }
            else if (x[ii] > cost->bds->ub[ii]){
                xuse[ii] = cost->bds->ub[ii];
            }
            else{
                xuse[ii] = x[ii];
            }
        }
    }
    else{
        xuse = x;
    }
//    printf("res = %d\n",res);

    /* if (res != 0){ */
    /*     printf("point is not in bounds for cost evaluation \n"); */
    /*     dprint(cost->bds->dim,x); */
    /*     return res; */
    /* } */
    
    *eval = function_train_eval(cost->cost,xuse);
    if (res != 0){
        free(xuse); xuse = NULL;
    }

    return 0;
    
}

/**********************************************************//**
    Black box (bb) interface to cost function evaluation

    \param[in]     t    - times of points
    \param[in]     x    - locations of points
    \param[in]     args - pointer to cost function structure

    \return value
**************************************************************/
double cost_eval_bb(double t,double * x,void * args)
{
    struct Cost * c = args;
    double out;
    int res = cost_eval(c,t,x,&out);
    assert (res == 0);
    return out;
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

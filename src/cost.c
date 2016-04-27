#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "c3.h"
#include "util.h"
#include "cost.h"

/** \struct Cost
 *  \brief Cost function
 *  \var Cost::d
 *  dimension of state space
 *  \var Cost::bds
 *  bounding box
 *  \var Cost::cost
 *  function train cost
 *  \var Cost::N
 *  Number of nodes of discretization in each dimension
 *  \var Cost::x
 *  discretization in each dimension
 *  \var Cost::Nobs
 *  Number of nodes that describes any obstacles
 *  \var Cost::xobs
 *  discretization of obstacles
 *  \var Cost::fm
 *  function monitor
 */
struct Cost 
{
    size_t d;
    struct BoundingBox * bds;
    struct FunctionTrain * cost;

    size_t * N;
    double ** x;

    size_t * Nobs;
    double ** xobs;

    struct FunctionMonitor * fm; // for storing evaluations
};
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
    cost->Nobs = NULL;
    cost->xobs = NULL;
    cost->fm = NULL;

    return cost;
}

/**********************************************************//**
    Save a cost function
**************************************************************/
int cost_save(struct Cost * cost, char *filename)
{
    assert (cost->cost != NULL);
    int res = function_train_save(cost->cost,filename);
    if (res == 1){
        // 0 should be good according to unix but haven't updated
        // function train code!!
        return 0; 
    }
    return 1;
}

/**********************************************************//**
    Load a cost function
**************************************************************/
int cost_load(struct Cost * cost, char * filename)
{
    assert (cost != NULL);
    function_train_free(cost->cost); cost->cost = NULL;
    cost->cost = function_train_load(filename);
    if (cost->cost != NULL){
        return 0;
    }
    return 1;
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
    newc->cost = function_train_copy(old->cost);
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
    if (old->Nobs != NULL){
        newc->Nobs = calloc_size_t(newc->d);
        memmove(newc->Nobs,old->Nobs,newc->d*sizeof(size_t));
    }
    if (old->xobs != NULL){
        newc->xobs = malloc_dd(newc->d);
        for (size_t ii = 0; ii < newc->d; ii++){
            newc->xobs[ii] = calloc_double(newc->Nobs[ii]);
            memmove(newc->xobs[ii],old->xobs[ii],newc->Nobs[ii]*sizeof(double));
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
        free_dd(c->d,c->xobs); c->xobs = NULL;
        free(c->Nobs); c->Nobs = NULL;
        function_monitor_free(c->fm); c->fm = NULL;
        free(c); c = NULL;
    }
}

/**********************************************************//**
    Return the dimension of a cost function
**************************************************************/
size_t cost_get_d(struct Cost * c)
{
    assert (c != NULL);
    return c->d;
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
   Get spacing between nodes
**************************************************************/
void cost_get_h(const struct Cost * c, double * h)
{
    assert (c != NULL);
    for (size_t ii = 0; ii < c->d; ii++){
        h[ii] = c->x[ii][1] - c->x[ii][0];
    }
}

/**********************************************************//**
    Get cost function ranks
**************************************************************/
size_t * cost_get_ranks(struct Cost * c)
{
    assert (c != NULL);
    assert (c->cost != NULL);
    return c->cost->ranks;
}

/**********************************************************//**
    Get cost norm
**************************************************************/
double cost_norm2(struct Cost * c)
{
    assert (c != NULL);
    return function_train_norm2(c->cost);
}

/**********************************************************//**
    Get norm2 difference between cost functions
**************************************************************/
double cost_norm2_diff(struct Cost * c, struct Cost * nc)
{
    assert (c != NULL);
    assert (nc != NULL);
    return function_train_norm2diff(c->cost,nc->cost);
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
    Add obstacle nodes                                                        

    \param[in,out] cost - allocated cost
    \param[in]     lb   - lower bounds of obstacles
    \param[in]     ub   - upper bounds of obstacles
    \param[in]     N    - number of nodes with which to represent obstacle
**************************************************************/
void cost_add_nodes(struct Cost * cost, double *lb, double *ub, size_t N)
{
    assert (cost != NULL);
    assert (cost->xobs == NULL);
    assert (cost->Nobs == NULL);
    cost->xobs = malloc_dd(cost->d);
    cost->Nobs = calloc_size_t(cost->d);
    for (size_t ii = 0; ii < cost->d; ii++){
        cost->Nobs[ii] = N;
        cost->xobs[ii] = linspace(lb[ii],ub[ii],N);
    }
}

/**********************************************************//**
    Get the total discretization could be quite large and result
    in overflow!!
**************************************************************/
size_t cost_get_size(const struct Cost * cost)
{
    assert (cost->N != NULL);
    assert (cost->x != NULL);
    double ** xuse = malloc_dd(cost->d);
    size_t * Nuse = calloc_size_t(cost->d);
    for (size_t ii = 0; ii < cost->d; ii++){
        if (cost->Nobs != NULL){
            xuse[ii] = c3sc_combine_and_sort(cost->N[ii],cost->x[ii],
                                             cost->Nobs[ii],cost->xobs[ii],
                                             Nuse+ii);
        }
        else{
            xuse[ii] = calloc_double(cost->N[ii]);
            memmove(xuse[ii],cost->x[ii],cost->N[ii]*sizeof(double));
            Nuse[ii] = cost->N[ii];
        }
    }
    size_t ntot = 1;
    for (size_t ii = 0; ii < cost->d; ii++){
        ntot *= Nuse[ii];
    }
    free_dd(cost->d,xuse); xuse = NULL;
    free(Nuse); Nuse = NULL;
    return ntot;
}

/**********************************************************//**
    Interpolate a new cost function from an old one

    \param[in,out] cnew - new cost function with allocated nodes and such
    \param[in]     cold - old cost function

    \note
    c->cost should be NULL
**************************************************************/
void cost_interpolate_new(struct Cost * cnew, struct Cost * cold)
{
    assert (cold->cost != NULL);
    assert (cnew->N != NULL);
    assert (cnew->x != NULL);
    function_train_free(cnew->cost); cnew->cost = NULL;
    double ** xuse = malloc_dd(cnew->d);
    size_t * Nuse = calloc_size_t(cnew->d);
    for (size_t ii = 0; ii < cnew->d; ii++){
        if (cnew->Nobs != NULL){
            xuse[ii] = c3sc_combine_and_sort(cnew->N[ii],cnew->x[ii],
                                             cnew->Nobs[ii],cnew->xobs[ii],
                                             Nuse+ii);
            
        }
        else{
            xuse[ii] = calloc_double(cnew->N[ii]);
            memmove(xuse[ii],cnew->x[ii],cnew->N[ii]*sizeof(double));
            Nuse[ii] = cnew->N[ii];
        }
    }

    cnew->cost = function_train_create_nodal(cold->cost,Nuse,xuse);
    free_dd(cnew->d,xuse); xuse = NULL;
    free(Nuse); Nuse = NULL;
}

/**********************************************************//**
    Divide the number of nodes in half
**************************************************************/
void cost_interp_inhalf(struct Cost * cnew, int inhalf)
{
 
    assert (cnew->N != NULL);
    assert (cnew->x != NULL);

    double * lb = cost_get_lb(cnew);
    double * ub = cost_get_ub(cnew);
    double ** xuse = malloc_dd(cnew->d);
    size_t * Nuse = calloc_size_t(cnew->d);
    for (size_t ii = 0; ii < cnew->d; ii++){
        if (inhalf == 1){
            cnew->N[ii] = cnew->N[ii]/2;
        }
        else{
            cnew->N[ii] = cnew->N[ii]*2;
        }
        free(cnew->x[ii]);
        cnew->x[ii] = linspace(lb[ii],ub[ii],cnew->N[ii]);

        if (cnew->Nobs != NULL){
            xuse[ii] = c3sc_combine_and_sort(cnew->N[ii],cnew->x[ii],
                                             cnew->Nobs[ii],cnew->xobs[ii],
                                             Nuse+ii);
        }
        else{
            xuse[ii] = calloc_double(cnew->N[ii]);
            memmove(xuse[ii],cnew->x[ii],cnew->N[ii]*sizeof(double));
            Nuse[ii] = cnew->N[ii];
        }

    }

    if (cnew->cost != NULL){
        struct FunctionTrain * newcost = 
            function_train_create_nodal(cnew->cost,Nuse,xuse);
        function_train_free(cnew->cost);
        cnew->cost = function_train_copy(newcost);
        function_train_free(newcost); newcost = NULL;
    }

    free_dd(cnew->d,xuse); xuse = NULL;
    free(Nuse); Nuse = NULL;
}


/**********************************************************//**
    Set cost function to an approximation of some input
    function

    \param[in,out] c       - cost structure
    \param[in]     f       - function
    \param[in]     args    - function arguments
    \param[in]     verbose - verbosity level for approximation
    \param[in]     aargs   - approximation arguments
    \note
    c->cost should be NULL
**************************************************************/
void cost_approx(struct Cost * c,
                 double (*f)(double *, void *),
                 void * args, int verbose,
                 const struct ApproxArgs * aargs)
{
    assert (c != NULL);
    assert (c->bds != NULL);
    assert (c->N != NULL);
    assert (c->x != NULL);
    if (c->cost != NULL){
        function_train_free(c->cost); c->cost = NULL;
    }

    double cross_tol = approx_args_get_cross_tol(aargs);
    double round_tol = approx_args_get_round_tol(aargs);
    size_t kickrank  = approx_args_get_kickrank(aargs);

    struct C3Approx * c3a = c3approx_create(CROSS,c->d,c->bds->lb,c->bds->ub);
    double ** xuse = malloc_dd(c->d);
    size_t * Nuse = calloc_size_t(c->d);
    for (size_t ii = 0; ii < c->d; ii++){
        if (c->Nobs != NULL){
            xuse[ii] = c3sc_combine_and_sort(c->N[ii],c->x[ii],
                                             c->Nobs[ii],c->xobs[ii],Nuse+ii);
            
        }
        else{
            xuse[ii] = calloc_double(c->N[ii]);
            memmove(xuse[ii],c->x[ii],c->N[ii]*sizeof(double));
            Nuse[ii] = c->N[ii];
        }
    }
    
    c3approx_init_lin_elem(c3a);
    c3approx_set_lin_elem_fixed(c3a,Nuse,xuse);

    size_t init_rank = 5;
    c3approx_init_cross(c3a,init_rank,verbose);
    c3approx_set_fiber_opt_brute_force(c3a,Nuse,xuse);
    
    c3approx_set_cross_tol(c3a,cross_tol);
    c3approx_set_cross_maxiter(c3a,10);
    c3approx_set_round_tol(c3a,round_tol);
    c3approx_set_adapt_kickrank(c3a,kickrank);
    size_t minN = Nuse[0];
    for (size_t ii = 0; ii < c->d; ii++){
        if (c->N[ii] < minN){ minN = Nuse[ii];}
    }
    c3approx_set_adapt_maxrank_all(c3a,minN);

    // determine starting nodes
    double ** start = malloc_dd(c->d);
    size_t * nstart = calloc_size_t(c->d);
    for (size_t ii = 0; ii < c->d; ii++){
        /* printf("x[%zu] = ",ii); dprint(Nuse[ii],xuse[ii]); */
        nstart[ii] = 5;
        assert (Nuse[ii] > 5);
        size_t mid = Nuse[ii]/2;
        if (c->Nobs != NULL){
            mid = c->Nobs[ii]/2;
        }
        start[ii] = calloc_double(nstart[ii]);
        if (c->Nobs != NULL){
            start[ii][0] = c->xobs[ii][mid];
        }
        else{
            start[ii][0] = xuse[ii][mid];
        }
        start[ii][1] = xuse[ii][1];
        start[ii][2] = xuse[ii][Nuse[ii]-2];
        start[ii][3] = xuse[ii][0];
        start[ii][4] = xuse[ii][Nuse[ii]-1];
        /* printf("start[%zu] = ",ii); dprint(nstart[ii],start[ii]); */
//        printf("start mid = %G\n",start[ii][0]); 
    }
    c3approx_set_start(c3a,nstart,start);

    c->cost = c3approx_do_cross(c3a,f,args);
    
    free(nstart); nstart = NULL;
    free_dd(c->d,start); start = NULL;
    free_dd(c->d,xuse); xuse = NULL;
    free(Nuse); Nuse = NULL;
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
    
    *eval = function_train_eval(cost->cost,xuse);
    if (res != 0){
        free(xuse); xuse = NULL;
    }

    return 0;
}

/**********************************************************//**
    Evaluate a cost function

    \param[in] x    - location in space at which to evaluate
    \param[in] eval - pointer to evaluation location

**************************************************************/
double cost_eval_to_wrap(double * x,void * arg)
{
    struct Cost * cost = arg;
    double out;
    double t = 0.0;
    int res = cost_eval(cost,t,x,&out);
    assert (res == 0);
    return out;
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
    (void)(t);
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
                    double * eval,
                    double * points,
                    double * evals)
{
    *eval = function_train_eval_co_perturb(cost->cost,x,points,evals);
    return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "c3.h"
#include "util.h"
#include "cost.h"

/** \struct Cost
\brief Cost function
\var Cost::d        dimension of state space
\var Cost::bds      bounding box
\var Cost::cost     function train cost
\var Cost::grid     grid of values
\var Cost::hashgrid map from grid values to indices
\var Cost::fm       function monitor
\var Cost::Nobs     number of nodes in each grid that are obstacles
\var Cost::obs      array of indices of grid which are obstacles
*/
struct Cost 
{
    size_t d;
    struct BoundingBox * bds;
    struct FunctionTrain * cost;
    double ** cores;

    struct c3Vector ** grid;
    struct HashGrid ** hashgrid;
    int cost_eval_grid;
    
    size_t * Nobs;
    size_t ** obs;

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
    cost->cores = malloc_dd(d);
    cost->cost_eval_grid = 1;
    cost->grid = NULL;
    cost->hashgrid = NULL;
    cost->Nobs = NULL;
    cost->obs  = NULL;


    cost->fm = NULL;

    return cost;
}

/**********************************************************//**
    Save a cost function
**************************************************************/
int cost_save(const struct Cost * cost, char *filename)
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
struct Cost * cost_copy_deep(const struct Cost * old)
{
    if (old == NULL){
        return NULL;
    }

    double * lb = bounding_box_get_lb(old->bds);
    double * ub = bounding_box_get_ub(old->bds);
    struct Cost * newc = cost_alloc(old->d,lb,ub);
    newc->cost = function_train_copy(old->cost);
    if (old->grid != NULL){
        newc->grid = c3vector_array_copy(newc->d,old->grid);
        newc->hashgrid = hash_grid_create_ndgrid(10000,newc->d,newc->grid);
    }
    if (old->Nobs != NULL){
        newc->Nobs = calloc_size_t(old->d);
        memmove(newc->Nobs,old->Nobs, old->d*sizeof(size_t));
        newc->obs = malloc(old->d * sizeof(size_t *));
        for (size_t ii = 0; ii < newc->d; ii++){
            if (newc->Nobs[ii] > 0){
                newc->obs[ii]  = calloc_size_t(newc->Nobs[ii]);
                memmove(newc->obs[ii],old->obs[ii],newc->Nobs[ii]*sizeof(size_t));
            }
            else{
                newc->obs[ii] = NULL;
            }
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
        c3vector_array_free(c->d,c->grid);
        hash_grid_free_ndgrid(c->d,c->hashgrid);
        function_monitor_free(c->fm); c->fm = NULL;
        free_dd(c->d,c->cores);
        if (c->Nobs != NULL){
            free(c->Nobs); c->Nobs = NULL;
        }
        if (c->obs != NULL){
            for (size_t ii = 0; ii < c->d; ii++){
                free(c->obs[ii]); c->obs[ii] = NULL;
            }
            free(c->obs);
        }
        free(c); c = NULL;
    }
}

/**********************************************************//**
    Return the dimension of a cost function
**************************************************************/
size_t cost_get_d(const struct Cost * c)
{
    assert (c != NULL);
    return c->d;
}

/**********************************************************//**
    Get the total discretization could be quite large and result
    in overflow!!
**************************************************************/
size_t cost_get_size(const struct Cost * cost)
{
    /* assert (cost->N != NULL); */
    assert (cost->grid != NULL);
    size_t ntot = 1;
    for (size_t ii = 0; ii < cost->d; ii++){
        ntot *= cost->grid[ii]->size;
    }
//    printf("ntot = %zu\n",ntot);
    
    return ntot;
}

/**********************************************************//**
    Return a reference to cost funciton lower bounds
**************************************************************/
double * cost_get_lb(const struct Cost * c)
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
double * cost_get_ub(const struct Cost * c)
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
        h[ii] = c->grid[ii]->elem[1] - c->grid[ii]->elem[0];
    }
}

/**********************************************************//**
    Get cost function ranks
**************************************************************/
size_t * cost_get_ranks(const struct Cost * c)
{
    assert (c != NULL);
    assert (c->cost != NULL);
    return c->cost->ranks;
}

/**********************************************************//**
    Get cost norm
**************************************************************/
double cost_norm2(const struct Cost * c)
{
    assert (c != NULL);
    return function_train_norm2(c->cost);
}

/**********************************************************//**
    Get norm2 difference between cost functions
**************************************************************/
double cost_norm2_diff(const struct Cost * c, const struct Cost * nc)
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
void cost_init_grid(struct Cost * cost,size_t * N,double ** x)
{
    assert (cost != NULL);
    assert (cost->grid == NULL);
    assert (cost->hashgrid == NULL);
    cost->grid = c3vector_array_alloc(cost->d);
    for (size_t ii = 0; ii < cost->d; ii++){
        /* dprint(N[ii],x[ii]); */
        cost->grid[ii] = c3vector_alloc(N[ii],x[ii]);
    }
    cost->hashgrid = hash_grid_create_ndgrid(10000,cost->d,cost->grid);
}

/**********************************************************//**
    Add obstacle nodes                                                        

    \param[in,out] cost - allocated cost
    \param[in]     lb   - lower bounds of obstacles
    \param[in]     ub   - upper bounds of obstacles
**************************************************************/
void cost_add_nodes(struct Cost * cost, double *lb, double *ub)
{
    assert (cost != NULL);
    assert (cost->obs == NULL);
    assert (cost->Nobs == NULL);
    cost->obs = malloc(cost->d * sizeof(size_t *));
    assert (cost->obs != NULL);
    cost->Nobs = calloc_size_t(cost->d);
    size_t Nobs = 4;
    for (size_t ii = 0; ii < cost->d; ii++){
        /* size_t newsize = cost->grid[ii]->size + 4; */
        size_t newsize;
        double * obs = linspace(lb[ii],ub[ii],Nobs);
        double * newg = c3sc_combine_and_sort(cost->grid[ii]->size,
                                              cost->grid[ii]->elem,
                                              Nobs,
                                              obs,
                                              &newsize);
        c3vector_free(cost->grid[ii]);
        cost->grid[ii] = c3vector_alloc(newsize,newg);
        printf("grid is\n \t ");
        dprint(cost->grid[ii]->size,cost->grid[ii]->elem);
        free(newg); newg = NULL;
        free(obs); obs = NULL;
    }

    for (size_t ii = 0; ii < cost->d; ii++){

        for (size_t jj = 0; jj < cost->grid[ii]->size; jj++){
            if ((cost->grid[ii]->elem[jj] > lb[ii]-1e-13) && (cost->grid[ii]->elem[jj] < ub[ii]+1e-13)){
                cost->Nobs[ii]++;
            }
        }
        cost->obs[ii] = calloc_size_t(cost->Nobs[ii]);
        size_t ind = 0;
        for (size_t jj = 0; jj < cost->grid[ii]->size; jj++){
            if ((cost->grid[ii]->elem[jj] > lb[ii]-1e-13) && (cost->grid[ii]->elem[jj] < ub[ii]+1e-13)){
                cost->obs[ii][ind] = jj;
                ind++;
            }
        }
    }
    
    size_t allnonzero = 1;
    for (size_t ii = 0; ii < cost->d; ii++){
        if (cost->Nobs[ii] == 0){
            allnonzero = 0;
            break;
        }
    }
    if (allnonzero == 0){
        fprintf(stderr,"There are no valid nodes in the grid that represent this obstacle!\n");
        /* printf("Boundary[%zu] in (%G,%G)\n",ii,lb[ii],ub[ii]); */
        /* printf("\t grid is "); dprint(cost->grid[ii]->size,cost->grid[ii]->elem); */
        exit(1);
    }
    

}

/**********************************************************//**
    Interpolate a new cost function from an old one

    \param[in,out] cnew - new cost function with allocated nodes and such
    \param[in]     cold - old cost function

    \note
    c->cost should be NULL
**************************************************************/
void cost_interpolate_new(struct Cost * cnew, const struct Cost * cold)
{
    assert (cold->cost != NULL);
    assert (cnew->grid != NULL);

    function_train_free(cnew->cost); cnew->cost = NULL;
    double ** xuse = malloc_dd(cnew->d);
    size_t * Nuse = calloc_size_t(cnew->d);
    for (size_t ii = 0; ii < cnew->d; ii++){
        xuse[ii] = cnew->grid[ii]->elem;
        Nuse[ii] = cnew->grid[ii]->size;
    }

    cnew->cost = function_train_create_nodal(cold->cost,Nuse,xuse);
    free_dd(cnew->d,xuse); xuse = NULL;
    free(Nuse); Nuse = NULL;
}

/**********************************************************//**
    Divide the number of nodes in half
**************************************************************/
/* void cost_interp_inhalf(struct Cost * cnew, int inhalf) */
/* { */
 
/*     assert (cnew->N != NULL); */
/*     assert (cnew->x != NULL); */

/*     double * lb = cost_get_lb(cnew); */
/*     double * ub = cost_get_ub(cnew); */
/*     double ** xuse = malloc_dd(cnew->d); */
/*     size_t * Nuse = calloc_size_t(cnew->d); */
/*     for (size_t ii = 0; ii < cnew->d; ii++){ */
/*         if (inhalf == 1){ */
/*             cnew->N[ii] = cnew->N[ii]/2; */
/*         } */
/*         else{ */
/*             cnew->N[ii] = cnew->N[ii]*2; */
/*         } */
/*         free(cnew->x[ii]); */
/*         cnew->x[ii] = linspace(lb[ii],ub[ii],cnew->N[ii]); */

/*         if (cnew->Nobs != NULL){ */
/*             xuse[ii] = c3sc_combine_and_sort(cnew->N[ii],cnew->x[ii], */
/*                                              cnew->Nobs[ii],cnew->xobs[ii], */
/*                                              Nuse+ii); */
/*         } */
/*         else{ */
/*             xuse[ii] = calloc_double(cnew->N[ii]); */
/*             memmove(xuse[ii],cnew->x[ii],cnew->N[ii]*sizeof(double)); */
/*             Nuse[ii] = cnew->N[ii]; */
/*         } */

/*     } */

/*     if (cnew->cost != NULL){ */
/*         struct FunctionTrain * newcost =  */
/*             function_train_create_nodal(cnew->cost,Nuse,xuse); */
/*         function_train_free(cnew->cost); */
/*         cnew->cost = function_train_copy(newcost); */
/*         function_train_free(newcost); newcost = NULL; */
/*     } */

/*     free_dd(cnew->d,xuse); xuse = NULL; */
/*     free(Nuse); Nuse = NULL; */
/* } */

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
                 double (*f)(const double *, void *),
                 void * args, int verbose,
                 const struct ApproxArgs * aargs)
{
    assert (c != NULL);
    assert (c->bds != NULL);
    assert (c->grid != NULL);
    if (c->cost != NULL){
        function_train_free(c->cost); c->cost = NULL;
    }

    double cross_tol = approx_args_get_cross_tol(aargs);
    double round_tol = approx_args_get_round_tol(aargs);
    size_t kickrank  = approx_args_get_kickrank(aargs);
    size_t maxrank   = approx_args_get_maxrank(aargs);
    size_t startrank = approx_args_get_startrank(aargs);

    struct C3Approx * c3a = c3approx_create(CROSS,c->d);
    struct LinElemExpAopts ** aopts =
        malloc(c->d * sizeof(struct LinElemExpAopts *));
    struct OneApproxOpts ** qmopts =
        malloc(c->d * sizeof(struct OneApproxOpts *));

    size_t init_rank = startrank;
    // determine starting nodes
    double ** start = malloc_dd(c->d);
    for (size_t ii = 0; ii < c->d; ii++){
        start[ii] = calloc_double(init_rank);

        for (size_t jj = 0; jj < init_rank; jj++){
            start[ii][jj] = c->grid[ii]->elem[jj];
        }
        
        if (c->Nobs != NULL){ // add the first element
            start[ii][0] = c->grid[ii]->elem[c->obs[ii][0]];
        }
        /* printf("ii=%zu ",ii); */
        /* dprint(init_rank,start[ii]); */
        /* printf("set it\n"); */
        /* dprint(c->grid[ii]->size,c->grid[ii]->elem); */

        aopts[ii] = lin_elem_exp_aopts_alloc(c->grid[ii]->size,c->grid[ii]->elem);

        qmopts[ii] = one_approx_opts_alloc(LINELM,aopts[ii]);
        c3approx_set_approx_opts_dim(c3a,ii,qmopts[ii]);
        /* printf("did it\n"); */
    }
    /* printf("really?!\n"); */
    c3approx_init_cross(c3a,init_rank,verbose,start);
    /* printf("why\n"); */
    c3approx_set_cross_tol(c3a,cross_tol);
    c3approx_set_cross_maxiter(c3a,30);
    /* c3approx_set_cross_maxiter(c3a,10); */
    c3approx_set_round_tol(c3a,round_tol);
    c3approx_set_adapt_kickrank(c3a,kickrank);
    /* printf("set all things\n"); */
    size_t minN = c->grid[0]->size;
    for (size_t ii = 0; ii < c->d; ii++){
        if (c->grid[ii]->size < minN){ minN = c->grid[ii]->size;}
    }
    if (maxrank < minN){
        c3approx_set_adapt_maxrank_all(c3a,maxrank);
    }
    else{
        c3approx_set_adapt_maxrank_all(c3a,minN);
    }
    /* printf("wrap function\n"); */
    struct Fwrap * fw = fwrap_create(c->d,"general");
    fwrap_set_f(fw,f,args);

    /* printf("start approximation\n"); */
    int adapt = 0;
    c->cost = c3approx_do_cross(c3a,fw,adapt);
    /* printf("approx end\n"); */
    

    free_dd(c->d,start); start = NULL;
    for (size_t ii = 0; ii < c->d; ii++){
        one_approx_opts_free_deep(&(qmopts[ii]));
    }
    free(qmopts);
    free(aopts);;
    c3approx_destroy(c3a); c3a = NULL;
    fwrap_destroy(fw);
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
              const double * x,
              double * eval)
{


    (void)(time);
    assert (cost->cost != NULL);
    int dothis = 0;
    if (dothis){
        assert (cost->hashgrid != NULL);
        size_t * ind = calloc_size_t(cost->d);
        /* printf("get ind\n"); */
        int success = hash_grid_ndgrid_get_ind(cost->hashgrid,cost->d,x,ind);
        if (success == 0){
            /* printf("obtain evaluation\n"); */
            *eval = function_train_eval_ind(cost->cost,ind);
            /* printf("evaluation is %G\n",*eval); */
            /* iprint_sz(cost->d,ind); */
            free(ind); ind = NULL;
            return 0;
        }
        free(ind); ind = NULL;
        printf("x that is not indexed is\n");
        dprint(cost->d,x);
        assert(1 == 0);
    }
    
    double * lb = bounding_box_get_lb(cost->bds);
    double * ub = bounding_box_get_ub(cost->bds);
    int res = c3sc_check_bounds(cost->d,lb,ub,x);
    
    /* printf("x = "); dprint(3,x); */
    double * xuse = NULL;
    if (res != 0){
        xuse = calloc_double(cost->d);
        for (size_t ii = 0; ii < cost->d; ii++){
            if (x[ii] < lb[ii]){
                xuse[ii] = lb[ii];
            }
            else if (x[ii] > ub[ii]){
                xuse[ii] = ub[ii];
            }
            else{
                xuse[ii] = x[ii];
            }
        }
    }

    /* for (size_t ii = 0; ii < cost->d; ii++){ */
    /*     int okd = 0; */
    /*     for (size_t jj = 0; jj < cost->N[ii]; jj++){ */
    /*         if (fabs(xuse[ii]-cost->x[ii][jj]) < 1e-15){ */
    /*             okd = 1; */
    /*             break; */
    /*         } */
    /*     } */
    /*     if (okd == 0){ */
    /*         fprintf(stderr,"Evaluation point is not a proper candidate\n"); */
    /*         printf("pt = "); dprint(cost->d,xuse); */
    /*         dprint(cost->N[ii], cost->x[ii]); */
    /*         exit(1); */
    /*     } */
    /* } */

    if (xuse != NULL){
        *eval = function_train_eval(cost->cost,xuse);
        free(xuse); xuse = NULL;
    }
    else{
        *eval = function_train_eval(cost->cost,x);
    }

    return 0;
}

/**********************************************************//**
    Evaluate a cost function

    \param[in] x   - location in space at which to evaluate
    \param[in] arg - pointer to cost function

**************************************************************/
double cost_eval_to_wrap(const double * x,void * arg)
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
double cost_eval_bb(double t,const double * x,void * args)
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

    \param[in]     cost   - cost funciton
    \param[in]     time   - time at which to evaluate
    \param[in]     x      - location at which to evaluate
    \param[in]     eval   - evaluation at x
    \param[in]     points - perturbed values of x[ii] - + in each dimension
    \param[in,out] evals  - space allocated to evaluation of each point

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
    (void)(time);
    assert (cost->cost != NULL);
    assert (cost->hashgrid != NULL);
    int dothis = 0;
    if (dothis){
        size_t * ind = calloc_size_t(cost->d);
        size_t * pert = calloc_size_t(cost->d*2);
        /* printf("get ind\n"); */
        int success = hash_grid_ndgrid_get_ind(cost->hashgrid,cost->d,x,ind);
        /* printf("got ind\n"); */
        /* dprint(cost->d*2, points); */

        if (success == 0){
            /* int exist; */
            for (size_t jj = 0; jj < cost->d; jj++){
                
                /* printf("jj = %zu\n",jj); */
                /* printf("grid = "); dprint(cost->grid[jj]->size,cost->grid[jj]->elem); */
                /* printf("hashgrid is\n"); */
                /* hash_grid_print(cost->hashgrid[jj],stdout); */
                /* printf("point is %3.15G\n",points[jj*cost->d]); */
                /* pert[jj*2] = hash_grid_get_ind(cost->hashgrid[jj],points[jj*2],&exist); */
                /* if (exist == 0){ */
                /*     success = 1; */
                /*     break; */
                /* } */

            
                /* pert[jj*2+1] = hash_grid_get_ind(cost->hashgrid[jj],points[jj*2+1],&exist); */
                /* if (exist == 0){ */
                /*     success = 1; */
                /*     break; */
                /* } */
                if (points[jj*2] < x[jj]){
                    pert[jj*2] = ind[jj]-1;
                }
                else{
                    pert[jj*2] = ind[jj];
                }
                if (points[jj*2+1] > x[jj]){
                    pert[jj*2+1] = ind[jj]+1;
                }
                else{
                    pert[jj*2+1] = ind[jj];
                }

            }
            if (success == 0){
                /* printf("sometimes here\n"); */
                *eval = function_train_eval_co_perturb_ind(cost->cost,ind,pert,evals);
            
                /* printf("x should be "); */
                /* dprint(cost->d, x); */
                /* printf("x is\n"); */
                /* for (size_t ii = 0; ii < cost->d; ii++){ */
                /*     printf("%G ",cost->grid[ii]->elem[ind[ii]]); */
                /* } */
                /* printf("\n"); */

                /* printf("pert should be "); */
                /* dprint(2*cost->d, points); */
                /* printf("pert is\n "); */
                /* for (size_t ii = 0; ii < cost->d; ii++){ */
                /*     printf("(%G,%G)  ",cost->grid[ii]->elem[pert[2*ii]], cost->grid[ii]->elem[pert[2*ii+1]]); */
                /* } */
                /* printf("\n"); */
                /* exit(1); */
            


                /* double eval2; */
                /* double * evals = calloc_double(cost->d*2); */
                /* *eval = function_train_eval_co_perturb(cost->cost,x,points,evals); */
                /* free(ind); ind = NULL; */
                /* free(pert); pert = NULL; */
                return 0;
            }

            /* return 0; */
        }
        free(ind); ind = NULL;
        free(pert); pert = NULL;
    }
    /* assert(1 == 0); */
    /* exit(1); */

    /* printf("in cost eval\n"); */
    /* printf("\t x = "); dprint(3,x); */
    /* printf("\t neigh = "); dprint(6,points); */
    *eval = function_train_eval_co_perturb(cost->cost,x,points,evals);
    
    /* printf("\n\n ------------------\n evaluate neighbor!\n"); */
    /* printf("x = "); dprint(cost->d, x); */
    /* for (size_t ii = 0; ii < cost->d; ii++){ */
    /*     dprint(cost->grid[ii]->size,cost->grid[ii]->elem); */
    /*     printf("perturb (-,+)=(%G,%G)\n",points[2*ii],points[2*ii+1]); */
    /* } */
    /*     int okd = 0; */
    /*     for (size_t jj = 0; jj < cost->N[ii]; jj++){ */
    /*         if (fabs(x[ii]-cost->x[ii][jj]) < 1e-15){ */
    /*             okd = 1; */
    /*             break; */
    /*         } */
    /*     } */
    /*     if (okd == 0){ */
    /*         fprintf(stderr,"Evaluation point is not a proper candidate\n"); */
    /*         printf("pt = "); dprint(cost->d,x); */
    /*         dprint(cost->N[ii], cost->x[ii]); */
    /*         exit(1); */
    /*     } */
    /*     okd = 0; */
    /*     for (size_t jj = 0; jj < cost->N[ii]; jj++){ */
    /*         if (fabs(points[2*ii]-cost->x[ii][jj]) < 1e-15){ */
    /*             okd = 1; */
    /*             break; */
    /*         } */
    /*     } */
    /*     if (okd == 0){ */
    /*         fprintf(stderr,"Evaluation point is not a proper candidate\n"); */
    /*         printf("pt = "); dprint(cost->d,x); */
    /*         dprint(cost->N[ii], cost->x[ii]); */
    /*         exit(1); */
    /*     } */
    /*     okd = 0; */
    /*     for (size_t jj = 0; jj < cost->N[ii]; jj++){ */
    /*         if (fabs(points[2*ii+1]-cost->x[ii][jj]) < 1e-15){ */
    /*             okd = 1; */
    /*             break; */
    /*         } */
    /*     } */
    /*     if (okd == 0){ */
    /*         fprintf(stderr,"Evaluation point is not a proper candidate\n"); */
    /*         printf("pt = "); dprint(cost->d,x); */
    /*         dprint(cost->N[ii], cost->x[ii]); */
    /*         exit(1); */
    /*     } */
    /* } */
    
    /* *eval = function_train_eval_co_perturb(cost->cost,x,points,evals); */
    /* printf("done evaluating neighbor\n ----------------------\n\n\n"); */
    return 0;
}

void cost_eval_fiber_ind(struct Cost * cost, const size_t * fixed_ind, 
                         size_t N, const size_t * ind, size_t fiber_dim,
                         double * out)
{
    assert (cost != NULL);
    assert (cost->cost != NULL);
    function_train_eval_fiber_ind(cost->cost, fixed_ind, N, ind, fiber_dim, out);
}


void cost_precompute_cores(struct Cost * cost)
{
    assert (cost != NULL);
    assert (cost->cost != NULL);
    size_t nrows, ncols, nvals;
    size_t * ranks = function_train_get_ranks(cost->cost);
    for (size_t ii = 0; ii < cost->d; ii++){
        free(cost->cores[ii]);
        nrows = ranks[ii];
        ncols = ranks[ii+1];
        nvals = cost->grid[ii]->size;
        struct GenericFunction ** arr = cost->cost->cores[ii]->funcs;
        cost->cores[ii] = calloc_double(nrows * ncols * nvals);
        for (size_t jj = 0; jj < nvals; jj++){
            generic_function_1darray_eval2_ind(nrows * ncols,arr,jj,
                                               cost->cores[ii] + jj * nrows * ncols);
        }
    }
}


void cost_eval_fiber_ind_nn(struct Cost * cost, const size_t * fixed_ind,
                            size_t N, const size_t * ind, size_t fiber_dim,
                            size_t * neighbors, double ** out)
{
    
}

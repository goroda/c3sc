#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "c3.h"
#include "control.h"
#include "simulate.h"
#include "util.h"


/** \struct Policy
 *  \brief Policy
 *  \var Policy::dx
 *  dimension of state space
 *  \var Policy::du
 *  dimension of control
 *  \var Policy::lbx
 *  Lower bounds for state space
 *  \var Policy::ubx
 *  Upper bounds for state space
 *  \var Policy::feedback
 *  Function mapping time,state to control
 *  \var Policy::farg
 *   arguments to feedback function
 */
struct Policy
{
    size_t dx;
    size_t du;
    
    struct BoundingBox * bds;
    int (*feedback)(double, double *, double *,void*);
    void * farg;

    size_t * N;
    double ** x;
    struct FT1DArray * ftc;

};


/**********************************************************//**
    Allocate memory for a policy

    \param[in]     dx  - dimension of state space
    \param[in]     du  - dimension of control space

    \return allocated policy
**************************************************************/
struct Policy * policy_alloc(size_t dx, size_t du)
{
    struct Policy * pol = NULL;
    pol = malloc(sizeof(struct Policy));
    if (pol == NULL){
        fprintf(stderr, "Cannot allocate Policy \n");
        exit(1);
    }
    pol->dx = dx;
    pol->du = du;
    pol->bds = NULL;
    pol->feedback = NULL;
    pol->farg = NULL;

    pol->N = NULL;
    pol->x = NULL;
    pol->ftc = NULL;
    return pol;
}

/**********************************************************//**
    Free a policy
**************************************************************/
void policy_free(struct Policy * pol)
{
    if (pol != NULL){
        bounding_box_free(pol->bds); pol->bds = NULL;
        free(pol->N); pol->N = NULL;
        free_dd(pol->du,pol->x); pol->x = NULL;
        ft1d_array_free(pol->ftc); pol->ftc = NULL;
        free(pol); pol = NULL;
    }
}

/**********************************************************//**
    Set bounds of the policy

    \param[in,out] pol - allocated policy
    \param[in]     lbx - lower bounds of state space
    \param[in]     ubx - upper bounds of state space
**************************************************************/
void policy_set_bounds(struct Policy * pol,double * lbx, double * ubx)
{
    assert (pol != NULL);
    assert (lbx != NULL);
    assert (ubx != NULL);
    pol->bds = bounding_box_vec(pol->dx,lbx,ubx);
}

/**********************************************************//**
    Get policy ranks
**************************************************************/
size_t * policy_get_ranks(struct Policy * pol,size_t ind)
{
    assert (pol != NULL);
    assert (pol->ftc != NULL);
    return pol->ftc->ft[ind]->ranks;
}

/**********************************************************//**
    Initialize the cost of a discrete cost function

    \param[in,out] pol - allocated policy
    \param[in]     N   - number of nodes in each dimension
    \param[in]     x   - nodes in each dimension
**************************************************************/
void policy_init_discrete(struct Policy * pol,size_t * N,double ** x)
{
    assert (pol != NULL);
    
    pol->x = malloc_dd(pol->dx);
    pol->N = calloc_size_t(pol->dx);
    for (size_t ii = 0; ii < pol->dx; ii++){
        pol->N[ii] = N[ii];
        pol->x[ii] = calloc_double(N[ii]);
        memmove(pol->x[ii],x[ii],N[ii]*sizeof(double));
    }
}

/**********************************************************//**
    Set policy to an approximation of some input
    function

    \param[in,out] c       - policy structure
    \param[in]     f       - function(x,ind,args)
    \param[in]     args    - function arguments
    \param[in]     verbose - verbosity level for approximation

    \note
    c->bds and c->ftc should be NULL
**************************************************************/
void policy_approx(struct Policy * c,
                   double (*f)(double *,size_t, void *),
                   void * args, int verbose)
{
    assert (c != NULL);
    assert (c->bds != NULL);
    assert (c->ftc == NULL);
    assert (c->N != NULL);
    assert (c->x != NULL);

    
    struct LinElemExpAopts ** aopts = NULL;
    aopts = malloc(c->dx*sizeof(struct LinElemExpAopts *));
    assert (aopts != NULL);

    for (size_t ii = 0; ii < c->dx;ii++){
        aopts[ii] = lin_elem_exp_aopts_alloc(c->N[ii], c->x[ii]);
    }

    struct FtApproxArgs * fapp = NULL;
    struct FiberOptArgs * fopt = NULL;
    fapp = ft_approx_args_create_le2(c->dx,aopts); 
    
    struct c3Vector ** c3v = c3vector_alloc_array(c->dx);
    size_t minN = c->N[0];
    for (size_t ii = 0; ii < c->dx; ii++){
        c3v[ii] = c3vector_alloc(c->N[ii],c->x[ii]);
        if (c->N[ii] < minN){
            minN = c->N[ii];
        }
    }
    fopt = fiber_opt_args_bf(c->dx,c3v);

    struct FtCrossArgs fca;
    ft_cross_args_init(&fca);
    fca.dim = c->dx;
    fca.ranks = calloc_size_t(c->dx+1);
    for (size_t ii = 0; ii < c->dx+1; ii++){
        fca.ranks[ii] = 3;
    }
    fca.ranks[0] = 1;
    fca.ranks[c->dx] = 1;
    fca.epsilon = 1e-7;
    fca.maxiter = 10;
    fca.epsround = 1e-10;
    fca.kickrank = 5;
    fca.maxiteradapt = (minN-3)/(fca.kickrank);
    assert( (3 + fca.kickrank*fca.maxiteradapt) <= minN);
    fca.verbose = verbose;
    fca.optargs = fopt;

    double ** start = malloc_dd(c->dx);
    
    //printf("c->dx = %zu, c->du = \n",c->dx,c->du);
    for (size_t ii = 0; ii < c->dx; ii++){
        assert (c->N[ii] > 5);
        size_t mid = c->N[ii]/2;
        //printf("mid = %zu\n",mid);
        start[ii] = calloc_double(3);
        start[ii][0] = c->x[ii][1];
        start[ii][1] = c->x[ii][mid];
        start[ii][2] = c->x[ii][c->N[ii]-2];
        //printf("start is "); dprint(3,start[ii]);
    }
    
    c->ftc = ft1d_array_cross(f,args,c->du,c->bds,
                              start,&fca,fapp);

    //exit(1);    
    for (size_t ii = 0; ii < c->dx; ii++){
        lin_elem_exp_aopts_free(aopts[ii]); aopts[ii] = NULL;
    }
    free(aopts); aopts = NULL;
    c3vector_free_array(c3v,c->dx); c3v = NULL;
    ft_approx_args_free(fapp); fapp = NULL;
    free(fca.ranks); fca.ranks = NULL;
    free_dd(c->dx,start);
    fiber_opt_args_free(fopt); fopt = NULL;
}

void policy_add_feedback(struct Policy * pol,
                         int (*f)(double,double*,double*,void*),
                         void * farg)
{
    assert (pol != NULL);
    pol->feedback = f;
    pol->farg = farg;
}

int policy_eval(struct Policy * pol, double time,
                double * x, 
                struct Control ** u)
{
    assert (pol != NULL);

    double * lbx = bounding_box_get_lb(pol->bds);
    double * ubx = bounding_box_get_ub(pol->bds);
    int res = c3sc_check_bounds(pol->dx,lbx,ubx,x);
    if (res != 0){
        return res;
    }
    
    assert (*u == NULL);
    (*u) = control_alloc();
    (*u)->d = pol->du;
    (*u)->u = calloc_double(pol->du);
    
    if (pol->ftc != NULL){
        ft1d_array_eval2(pol->ftc,x,(*u)->u);
    }
    else{
        assert (pol->feedback != NULL);        
        res = pol->feedback(time,x,(*u)->u,pol->farg);
    }



    return res;
}

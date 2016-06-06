#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "cdyn/src/simulate.h"

#include "util.h"
#include "cost.h"
#include "dp.h"
#include "tensmarkov.h"


/** \struct C3SC
 *  \brief Stochastic control problem
 *  \var C3SC::type
 *  problem type
 *  \var C3SC::dx
 *  dimension of state
 *  \var C3SC::du
 *  dimension of control
 *  \var C3SC::dw
 *  dimension of noise
 *  \var C3SC::dpih
 *  infinite horizon problem
 *  \var C3SC::dpfh
 *  finite horizon problem
 */
struct C3SC
{
    enum SCTYPE type;
    size_t dx;
    size_t du;
    size_t dw;
    double * lbx;
    double * ubx;

    struct DPih * dpih;
    struct DPfh * dpfh;

    // auxilary stuff
    struct Boundary * bound;
    struct Dyn * dyn;
    struct MCA * mca;
    struct Cost * cost;
//    struct ImplictPolicy * pol;
    struct c3Opt * opt;
};

/**********************************************************//**
    Create a stochastic control problem

    \param[in] type  - IH (infinite horizon) or FH (finite horizon)
    \param[in] dx    - dimension of state space
    \param[in] du    - dimension of control space
    \param[in] dw    - dimension of noise

    \return stochastic control problem
**************************************************************/
struct C3SC * c3sc_create(enum SCTYPE type, size_t dx, size_t du, size_t dw)
{
    struct C3SC * sc = malloc(sizeof(struct C3SC));
    if (sc == NULL){
        fprintf(stderr,"Cannot create c3sc\n");
        exit(1);
    }
    
    sc->type = type;
    sc->dx = dx;
    sc->du = du;
    sc->dw = dw;
    sc->lbx = NULL;
    sc->ubx = NULL;
    
    sc->bound = NULL;
    sc->dyn = NULL;
    sc->mca = NULL;
    sc->cost = NULL;
//    sc->pol = NULL;
    sc->opt = NULL;
    return sc;
}

/**********************************************************//**
    Free memory allocated to stochastic control problem
**************************************************************/
void c3sc_destroy(struct C3SC * sc)
{
    if (sc != NULL){
        assert (sc->type != FH);
        free(sc->lbx); sc->lbx = NULL;
        free(sc->ubx); sc->ubx = NULL;
        dpih_free_deep(sc->dpih); sc->dpih = NULL;
        free(sc); sc = NULL;
    }
}

/**********************************************************//**
    Add bounds on state and initialize boundary conditions
**************************************************************/
void
c3sc_set_state_bounds(struct C3SC * sc, double * lb, double * ub)
{
    assert (sc != NULL);
    if (lb != NULL){
        sc->lbx = calloc_double(sc->dx);
        memmove(sc->lbx,lb,sc->dx * sizeof(double));
    }
    if (ub != NULL){
        sc->ubx = calloc_double(sc->dx);
        memmove(sc->ubx,ub,sc->dx * sizeof(double));
    }

    sc->bound = boundary_alloc(sc->dx,sc->lbx,sc->ubx);
}

/**********************************************************//**
    Set the external boundary type of a particular dimension
**************************************************************/
void c3sc_set_external_boundary(struct C3SC * sc, size_t dim,
                                char * type)
{
    boundary_external_set_type(sc->bound,dim,type);
}

/**********************************************************//**
    Add a rectangular obstacle
**************************************************************/
void c3sc_add_obstacle(struct C3SC * sc, double * center, double * lengths)
{
    boundary_add_obstacle(sc->bound,center,lengths);
    
}

/**********************************************************//**
    Add dynamics of stochastic cotnrol

    \param[in,out] sc    - stochastic control problem
    \param[in]     b     - drift dynamics
    \param[in]     bargs - drift dynamics arguments
    \param[in]     s     - diffusion dynamics
    \param[in]     sargs - diffusion dynamics arguments
**************************************************************/
void c3sc_add_dynamics(struct C3SC * sc,  
                       int (*b)(double,const double*,const double*,
                                double*,double*,void*),
                       void * bargs,
                       int (*s)(double,const double*,const double*,
                                double*,double*,void*),
                       void * sargs)
{
    assert (sc != NULL);
    struct Drift * drift = drift_alloc(sc->dx,sc->du);
    drift_add_func(drift,b,bargs);
    struct Diff * diff = diff_alloc(sc->dx,sc->du,sc->dw);
    diff_add_func(diff,s,sargs);
    sc->dyn = dyn_alloc(drift,diff);
}

/**********************************************************//**
    Initialize the markov chain approximation method 
    (cost and mca structure)

    \param[in,out] sc - stochastic control problem
    \param[in]     N  - number of nodes in each dimension
**************************************************************/
void c3sc_init_mca(struct C3SC * sc, size_t * N)
{
    assert (sc != NULL);
    if (sc->dyn == NULL){
        fprintf(stderr,"Must add dynamics to c3sc before initializing\n");
        fprintf(stderr,"the Markov Chain approximation algorithm\n");
        exit(1);
    }
    if (sc->bound == NULL){
        fprintf(stderr,"Must add boundary cond. to c3sc before initializing\n");
        fprintf(stderr,"the Markov Chain approximation algorithm\n");
        exit(1);
    }

    sc->cost = cost_alloc(sc->dx,sc->lbx,sc->ubx);
    double ** x = malloc_dd(sc->dx);
    double * h = calloc_double(sc->dx);
    for (size_t ii = 0; ii < sc->dx; ii++){
        x[ii] = linspace(sc->lbx[ii],sc->ubx[ii],N[ii]);
        h[ii] = x[ii][1] - x[ii][0];
        /* printf("x[%zu] = ",ii);dprint(N[ii],x[ii]); */
        /* printf("(%G,%G)\n",sc->lbx[ii],sc->ubx[ii]); */
    }
    cost_init_grid(sc->cost,N,x);

    size_t nobs = boundary_get_nobs(sc->bound);
    for (size_t jj = 0; jj < nobs; jj++){
        double * lb = boundary_obstacle_get_lb(sc->bound,jj);
        double * ub = boundary_obstacle_get_ub(sc->bound,jj);
        cost_add_nodes(sc->cost,lb,ub);
    }

    sc->mca = mca_alloc(sc->dx,sc->du,sc->dw,h);
    mca_attach_dyn(sc->mca,sc->dyn);
    mca_attach_bound(sc->mca,sc->bound);
    
    free(h); h = NULL; 
    free_dd(sc->dx,x); x = NULL;
}

/**********************************************************//**
    Attach (by reference )the optimization options and 
    algorithms to the problem

    \param[in,out] sc  - stochastic control problem
    \param[in]     opt - optimization structure
**************************************************************/
void c3sc_attach_opt(struct C3SC * sc, struct c3Opt * opt)
{
    assert(sc != NULL);
    assert(opt != NULL);
    sc->opt = opt;
}

/**********************************************************//**
    Attach (by reference )the optimization options and 
    algorithms to the problem

    \param[in,out] sc   - stochastic control problem
    \param[in]     beta - discount factor
    \param[in]     s    - stage cost
    \param[in]     b    - boundary cost
    \param[in]     o    - obstacle cost (NULL if none)
**************************************************************/
void c3sc_init_dp(struct C3SC * sc, double beta,
                  int (*s)(double,const double*,const double*,double*,double*),
                  int (*b)(double,const double*,double*),
                  int (*o)(const double*,double*))
{
    assert(sc != NULL);

    if (sc->dyn == NULL){
        fprintf(stderr,"Must add dynamics to c3sc before initializing\n");
        fprintf(stderr,"the dynamic program\n");
        exit(1);
    }
    if (sc->bound == NULL){
        fprintf(stderr,"Must add boundary cond. to c3sc before initializing\n");
        fprintf(stderr,"the dynamic program\n");
        exit(1);
    }
    if (sc->mca == NULL){
        fprintf(stderr,"Must initialize markov chain approximation method\n");
        fprintf(stderr,"before initializing the dynamic program\n");
        exit(1);
    }
    if (sc->opt == NULL){
        fprintf(stderr,"Must add optimizer to c3sc before initializing\n");
        fprintf(stderr,"the dynamic program\n");
        exit(1);
    }
    
    assert (sc->type == IH);
    sc->dpih = dpih_alloc(beta,s,b,o);
    dpih_attach_mca(sc->dpih, &(sc->mca));
    dpih_attach_cost(sc->dpih, &(sc->cost));
    dpih_attach_opt(sc->dpih, &(sc->opt));
}


/**********************************************************//**
    Load a cost function
**************************************************************/
int c3sc_cost_load(struct C3SC * sc, char * filename )
{
    assert (sc != NULL);
    assert (sc->cost != NULL);
    int load_success = cost_load(sc->cost,filename);
    return load_success;
}

/**********************************************************//**
    Approximate a cost function
**************************************************************/
int c3sc_cost_approx(struct C3SC * sc, double (*f)(const double *, void *),
                     void * args, int verbose,
                     const struct ApproxArgs * aargs)
{
    assert (sc != NULL);
    assert (sc->cost != NULL);
    cost_approx(sc->cost,f,args,verbose,aargs);
    return 0;
}


/**********************************************************//**
    Update the cost function to be that corresponding 
    to the currently stored implicit policy
**************************************************************/
void c3sc_pol_solve(struct C3SC * sc, size_t max_solve_iter,
                    double solve_tol,
                    int verbose,
                    const struct ApproxArgs * aargs)
{
    struct ImplicitPolicy * pol = c3sc_create_implicit_policy(sc);
    dpih_iter_pol_solve(sc->dpih,pol,max_solve_iter,solve_tol,verbose,aargs);
    sc->cost = dpih_get_cost(sc->dpih);
    implicit_policy_free(pol);
}

/**********************************************************//**
   Take a step of value iteration        
**************************************************************/
double c3sc_iter_vi(struct C3SC * sc,int verbose,const struct ApproxArgs * aargs,
                    struct C3SCDiagnostic * diag)
{
    struct Cost * newcost = dpih_iter_vi(sc->dpih,verbose-1,aargs,diag);
    struct Cost * oldcost = dpih_get_cost(sc->dpih);
    double diff = cost_norm2_diff(newcost,oldcost);
    dpih_attach_cost_ow(sc->dpih,newcost);
    sc->cost = dpih_get_cost(sc->dpih);
    cost_free(newcost); newcost = NULL;
    return diff;
}

/**********************************************************//**
    Get a reference to the dynamic programming problem
**************************************************************/
void * c3sc_get_dp(const struct C3SC * sc)
{
    assert (sc != NULL);
    assert (sc->type == IH);
    return sc->dpih;
}

/**********************************************************//**
    Get a reference to the cost
**************************************************************/
struct Cost * c3sc_get_cost(const struct C3SC * sc)
{
    assert (sc != NULL);
    return sc->cost;
}

/**********************************************************//**
    Get the control size
**************************************************************/
size_t c3sc_get_du(const struct C3SC * sc)
{
    assert (sc != NULL);
    return sc->du;
}

/** \struct ImplicitPolicy
 *  \brief Stores information needed to evaluate implicit control
 *  \var ImplicitPolicy::dpalloc
 *  flag to indicate whether need to free dp on cleanup
 *  \var ImplicitPolicy::dp
 *  infinite horizon dp problem
 *  \var ImplicitPolicy::du
 *  size of control
 *  \var ImplicitPolicy::fm
 *  store previously calculated controls
 *  \var ImplicitPolicy::transform
 *  transformation function in case need to perform state transformation before evaluating policy
 */
struct ImplicitPolicy
{
    int dpalloc;
    struct DPih * dp;
    size_t du;
    struct HashtableCpair ** fm;

    size_t dx;
    void (*transform)(size_t,const double*,double*);
};

/**********************************************************//**
   Allocate memory for policy
**************************************************************/
struct ImplicitPolicy * implicit_policy_alloc()
{
    struct ImplicitPolicy * ip = malloc(sizeof(struct ImplicitPolicy));
    if (ip == NULL){
        fprintf(stderr, "Failure allocating memory for implicit policy\n");
        exit(1);
    }
    ip->dpalloc = 0;
    ip->dp = NULL;
    ip->fm = NULL;

    ip->dx = 0;
    ip->transform = NULL;
    return ip;
}

/**********************************************************//**
    Create an implicit policy
**************************************************************/
struct ImplicitPolicy * c3sc_create_implicit_policy(struct C3SC * sc)
{
    assert (sc != NULL);
    assert (sc->type == IH);
    struct ImplicitPolicy * ip = implicit_policy_alloc();
    struct DPih * dp = dpih_copy_deep(sc->dpih);
    implicit_policy_set_dp(ip,dp);


    ip->du = c3sc_get_du(sc);
    ip->fm = malloc(ip->du * sizeof(struct HashtableCpair *));
    if (ip->fm == NULL){
        fprintf(stderr,"Cannot allocate memory for implicit policy\n");
        exit(1);
    }
    size_t d = dpih_get_d(dp);
    for (size_t ii = 0; ii < ip->du; ii++)
    {
        ip->fm[ii] = create_hashtable_cp(d*1000);
    }
    /* printf("created implicit policy!\n"); */

    return ip;
}

/**********************************************************//**
    Add transform to the policy
**************************************************************/
void implicit_policy_add_transform(struct ImplicitPolicy * pol,size_t dx, void (*f)(size_t, const double *, double *))
{
    pol->dx = dx;
    pol->transform = f;
}


/**********************************************************//**
   Free policy
**************************************************************/
void implicit_policy_free(struct ImplicitPolicy * ip)
{
    if (ip != NULL){
        if (ip->dpalloc == 1){
            dpih_free_deep(ip->dp); ip->dp = NULL;
        }
        for (size_t ii = 0; ii < ip->du; ii++){
            free_hashtable_cp(ip->fm[ii]); ip->fm[ii] = 0;
        }
        free(ip->fm); ip->fm = NULL;
        free(ip); ip = NULL;
    }
}

/**********************************************************//**
   SetDP
**************************************************************/
void implicit_policy_set_dp(struct ImplicitPolicy * ip, struct DPih * dp)
{
    if (ip->dpalloc == 1){
        if (ip->dp != NULL){
            dpih_free_deep(ip->dp); ip->dp = NULL;
        }
    }
    ip->dpalloc = 1;
    ip->dp = dp;
}

/**********************************************************//**
   Evaluate a policy
   
   \param[in]     ip  - policy to evaluate
   \param[in]     t   - time at which to evaluate
   \param[in]     xin - location at which to evaluate
   \param[in,out] u   - resulting control

   \return 0 - success, else failure
**************************************************************/
int implicit_policy_eval(struct ImplicitPolicy * ip,double t,
                         const double * xin, double * u)
{
    (void)(t);

    /* printf("here\n"); */
    struct DPX2 dpx;
    dpx.dp = ip->dp;

    struct MCA * mca = dpih_get_mca(dpx.dp);
    struct Cost * cost = dpih_get_cost(dpx.dp);
    size_t dx = mca_get_dx(mca);
    size_t du = mca_get_du(mca);
    struct Boundary * bound = mca_get_boundary(mca);
    double * h = mca_get_h(mca);
    
    struct BoundInfo * bi = NULL; 
    struct Node * node = NULL; 

    double * x = (double *)xin;
    if (ip->transform == NULL){
        bi = boundary_type(bound,t,xin);
        node = node_init(dx,xin,h,cost,bi);
    }
    else{
        x = calloc_double(ip->dx);
        ip->transform(ip->dx,xin,x);
        bi = boundary_type(bound,t,x);
        node = node_init(dx,x,h,cost,bi);
    }
    dpx.node = node;
    dpx.bi = bi;

    int serialize = 1;
    
    if (serialize == 1){
        /* if (dx > 100) assert(1 == 0); */
        /* size_t ind[100]; */
        /* for (size_t jj = 0; jj < dx; jj++){ */
        /*     ind =  */
        /* } */
        
        char * ser = serialize_darray_to_text(dx,(double *)dpx.node->x);
        char * sval = lookup_key(ip->fm[0],ser);
        /* if (1 == 0){ */
        if (sval != NULL){
            /* printf("sval exists\n"); */
            u[0] = deserialize_double_from_text(sval);
            for (size_t ii = 1; ii < du; ii++){
                char * sval2 = lookup_key(ip->fm[ii],ser);
                u[ii] = deserialize_double_from_text(sval2);
                free(sval2); sval2 = NULL;
            }
        }
        else{
            /* printf("sval doesnt exists\n"); */
            struct c3Opt * opt = dpih_get_opt(dpx.dp);
            assert (opt != NULL);
            assert (cost != NULL);
            assert (mca != NULL);
            double val = 0.0;
            /* c3opt_add_objective(opt,dpih_rhs_opt_bb,&dpx); */
            c3opt_add_objective(opt,dpih_bellman_evalu,&dpx);
            /* printf("minimize \n"); */
            int res = c3opt_minimize(opt,u,&val);
            /* printf("minimized!\n"); */
            if (res < -1){
                printf("max iter reached in optimization res=%d\n",res);

                printf("x = ");
                dprint(dx,dpx.node->x);
                dprint(du,u);
                for (size_t ii = 0; ii < du; ii++){
                    u[ii] = 0.0;
                }
                val = 0.0;
                printf("restart with verbose\n");
                c3opt_set_verbose(opt,2);
                int res2 = c3opt_minimize(opt,u,&val);
                printf("res2 = %d\n",res2);
            }

            for (size_t ii = 0; ii < du; ii++){
                char * sval2 = serialize_double_to_text(u[ii]);
                struct Cpair * cp = cpair_create(ser,sval2);
                add_cpair(ip->fm[ii],cp);
                cpair_free(cp); cp = NULL;
                free(sval2);
            }
        }
        free(ser); ser = NULL;
        free( sval); sval = NULL;
    }
    else{
        struct c3Opt * opt = dpih_get_opt(dpx.dp);
        assert (opt != NULL);
        assert (cost != NULL);
        assert (mca != NULL);
        double val = 0.0;
        /* c3opt_add_objective(opt,dpih_rhs_opt_bb,&dpx); */
        c3opt_add_objective(opt,dpih_bellman_evalu,&dpx);
        /* printf("minimize \n"); */
        int res = c3opt_minimize(opt,u,&val);
        /* printf("minimized!\n"); */
        if (res < -1){
            printf("max iter reached in optimization res=%d\n",res);

            printf("x = ");
            dprint(dx,dpx.node->x);
            dprint(du,u);
            for (size_t ii = 0; ii < du; ii++){
                u[ii] = 0.0;
            }
            val = 0.0;
            printf("restart with verbose\n");
            c3opt_set_verbose(opt,2);
            int res2 = c3opt_minimize(opt,u,&val);
            printf("res2 = %d\n",res2);
        }
    }
//    printf("we are done!\n");
    if (ip->transform != NULL){
        free(x); x = NULL;
//        free(dpx.x); dpx.x = NULL;
    }
    bound_info_free(bi); bi = NULL;
    node_free(node); node = NULL;
    
    int res_pol = 0;
    return res_pol;
}


/**********************************************************//**
   Interface for integrator
**************************************************************/
int implicit_policy_controller(double t,const double * x, double * u, void * args)
{

    struct ImplicitPolicy * ip = args;
    return implicit_policy_eval(ip,t,x,u);
}



/////////////////////////////////////////////////////////////////////////
// diagnostic stuff
////////////////////////////////////////////////////////////////////////
struct C3SCDiagnostic
{
    size_t itervi;
    size_t iterpi;
    struct Trajectory * traj_vi;
};

struct C3SCDiagnostic * c3sc_diagnostic_init()
{
    struct C3SCDiagnostic * diag = malloc(sizeof(struct C3SCDiagnostic));
    if (diag == NULL){
        fprintf(stderr,"Failure allocating diagnostic arguments\n");
        exit(1);
    }
    diag->itervi = 0;
    diag->iterpi = 0;
    diag->traj_vi = NULL;

    return diag;
}

void c3sc_diagnostic_free(struct C3SCDiagnostic * diag)
{
    if (diag != NULL){
        trajectory_free(diag->traj_vi); diag->traj_vi = NULL;
        free(diag); diag = NULL;
    }
}

int c3sc_diagnostic_save(struct C3SCDiagnostic * diag, 
                         char * filename, size_t nprec)
{
    FILE * fp = fopen(filename,"w");
    if (fp == NULL){
        return 1;
    }
    trajectory_print(diag->traj_vi,fp,nprec);
    fclose(fp);
    return 0;
}

void c3sc_diagnostic_vi_update(struct C3SCDiagnostic * diag,
                               struct Cost * newcost,
                               struct DPih * dp,
                               struct FunctionMonitor * fm)
{
    struct Cost * old = dpih_get_cost(dp);
    size_t d = cost_get_d(old);

    double norm = cost_norm2(newcost);
    double diff = cost_norm2_diff(old,newcost);
    size_t * ranks = cost_get_ranks(newcost);
    size_t ntot = cost_get_size(newcost);
    size_t nevals = nstored_hashtable_cp(fm->evals);
    double frac = (double) nevals / ntot;

    size_t ndim = d+1+3;
    double * pt = calloc_double(ndim);
    pt[0] = norm;
    pt[1] = diff;
    pt[2] = frac;
    for (size_t ii = 0; ii < d+1; ii++){
        pt[ii+3] = (double) ranks[ii];
    }
    
    diag->itervi += 1;
    trajectory_add(&(diag->traj_vi),ndim,0,(double)(diag->itervi),pt,NULL);

    free(pt); pt = NULL;

}


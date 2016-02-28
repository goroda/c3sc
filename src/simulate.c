#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "simulate.h"
#include "control.h"
#include "integrate.h"
#include "tensmarkov.h"

/**********************************************************//**
    Evaluate a linear transformation

    \param[in]  lt  - linear transformation
    \param[in]  x   - (d,)
    \param[out] out - must be preallocated to (d,)
**************************************************************/
void lin_transform_eval(struct LinTransform * lt, double * x, 
                        double * out)
{
    for (size_t ii = 0; ii < lt->d; ii++){
        out[ii] = lt->slope[ii]*x[ii] + lt->offset[ii];
    }
}

/**********************************************************//**
    Evaluate the slope of the transformation in dimension ii

    \param[in]  lt - linear transformation
    \param[in]  ii - dimension 

    \return Slope lt->slope[ii]
**************************************************************/
double lin_transform_get_slopei(struct LinTransform * lt,
                                size_t ii)
{
    assert (ii < lt->d);
    return lt->slope[ii];
}

/**********************************************************//**
    Allocate State structure and set everything to 0/null

    \return allocated state
**************************************************************/
struct State * state_alloc()
{
    struct State * ss = NULL;
    ss = malloc(sizeof(struct State));
    if (ss == NULL){
        fprintf(stderr, "Cannot Allocate System State\n");
        exit(1);
    }
    ss->d = 0;
    ss->t = 0.0;
    ss->x = NULL;
    ss->spec = NULL;
    return ss;
}

/**********************************************************//**
    Free memory allocated for state
**************************************************************/
void state_free(struct State *s)
{
    if (s != NULL){
        free(s->x); s->x = NULL;
        free(s); s = NULL;
    }
}

/**********************************************************//**
    Initialize a state by copying \f$ x \f$ 
**************************************************************/
void state_init(struct State * s,size_t d,double time,double * x)
{
    assert (s != NULL);
    assert (x != NULL);
    s->d = d;
    s->t = time;
    s->x = calloc_double(d);
    memmove(s->x,x,d*sizeof(double));
    s->spec = NULL;
}

/**********************************************************//**
    Initialize a state to zero
**************************************************************/
void state_init_zero(struct State * s,size_t d, double time)
{
    s->d = d;
    s->t = time;
    s->x = calloc_double(d);
    s->spec = NULL;
}

/**********************************************************//**
    Copy a state
**************************************************************/
struct State * state_copy(struct State * ss)
{
    if (ss == NULL){
        return NULL;
    }
    struct State * s = state_alloc();
    state_init(s,ss->d,ss->t,ss->x);
    return s;
}

/**********************************************************//**
    Update a state by copying
**************************************************************/
void state_up(struct State * s,size_t d,double time,double * x)
{
    s->d = d;
    s->t = time;
    assert (s->x != NULL);
    memmove(s->x,x,d*sizeof(double));
    s->spec = NULL;
}

/***********************************************************//**
    Return dimension of state
**************************************************************/
size_t state_getd(struct State * s)
{
    return s->d;
}

/**********************************************************//**
    Return time of state
**************************************************************/
double state_gett(struct State * s)
{
    return s->t;
}

/**********************************************************//**
    Return state by reference
**************************************************************/
double * state_getx_ref(struct State * s)
{
    return s->x;
}

/**********************************************************//**
    Print a state

    \param[in] s          - state
    \param[in] fp         - stream to print to
    \param[in] print_time - 1 print time, 0 dont
    \param[in] prec       - precision to print
**************************************************************/
void state_print(struct State * s, FILE * fp, int print_time, 
                 int prec)
{
    if (s == NULL){
        fprintf(fp,"NULL");
        return;
    }
    else if (s->x == NULL){
        fprintf(fp,"NULL");
        return;
    }

    if (print_time == 0){
        for (size_t ii = 0; ii < s->d; ii++){
            fprintf(fp,"%3.*G ",prec,s->x[ii]);
        }
    }
    else{
        fprintf(fp,"%3.*G ",prec,s->t);
        for (size_t ii = 0; ii < s->d; ii++){
            fprintf(fp,"%3.*G ",prec,s->x[ii]);
        }
    }
}

//////////////////////////////////////////////////////////////

/**********************************************************//**
    Allocate control
**************************************************************/
struct Control * control_alloc()
{
    struct Control * ss = NULL;
    ss = malloc(sizeof(struct Control));
    if (ss == NULL){
        fprintf(stderr, "Cannot Allocate Control\n");
        exit(1);
    }
    ss->d = 0;
    ss->u = NULL;

    return ss;
}

/**********************************************************//**
    Free control
**************************************************************/
void control_free(struct Control *s)
{
    if (s != NULL){
        free(s->u); s->u = NULL;
        free(s); s = NULL;
    }
}

/**********************************************************//**
    Initialize control

    \param[in] s - control
    \param[in] d - dimension of control
    \param[in] u - control value to copy
**************************************************************/
void control_init(struct Control * s, size_t d, double * u)
{
    assert (u != NULL);
    s->d = d;
    s->u = calloc_double(d);
    memmove(s->u,u,d*sizeof(double));
}

/**********************************************************//**
    Copy a control
**************************************************************/
struct Control * control_copy(struct Control * ss)
{
    struct Control * s = control_alloc();
    control_init(s,ss->d,ss->u);
    return s;
}

/**********************************************************//**
    Update a control 

    \note
    Overwrites previous structure elements
**************************************************************/
void control_up(struct Control * s,size_t d, double * u)
{
    assert (s->u != NULL);
    s->d = d;
    memmove(s->u,u,d*sizeof(double));
}

/**********************************************************//**
    Get dimension of control
**************************************************************/
size_t control_getd(struct Control * s)
{
    return s->d;
}

/**********************************************************//**
    Get a reference to control value
**************************************************************/
double * control_getu_ref(struct Control * s)
{
    return s->u;
}

/**********************************************************//**
    Print control

    \param[in] c    - control
    \param[in] fp   - stream to print to
    \param[in] prec - precision to print
**************************************************************/
void control_print(struct Control * s, FILE * fp, int prec)
{
    if (s == NULL){
        fprintf(fp,"NULL");
        return;
    }
    else if (s->u == NULL){
        fprintf(fp,"NULL");
        return;
    }
    for (size_t ii = 0; ii < s->d; ii++){
        fprintf(fp,"%3.*G ",prec,s->u[ii]);
    }
}

////////////////////////////////////////////////////////

/**********************************************************//**
    Allocate Trajectory
**************************************************************/
struct Trajectory * trajectory_alloc()
{
    struct Trajectory * t = malloc(sizeof(struct Trajectory ));
    if (t == NULL){
        fprintf(stderr, "Cannot allocate trajectory\n");
        exit(1);
    }
    
    t->s = NULL;
    t->u = NULL;
    t->next = NULL;
    return t;
}

/**********************************************************//**
    Free Trajectory
**************************************************************/
void trajectory_free(struct Trajectory * traj)
{
    if (traj != NULL){
        state_free(traj->s); traj->s = NULL;
        control_free(traj->u); traj->u = NULL;
        trajectory_free(traj->next); traj->next = NULL;
        free(traj); traj = NULL;
    }
}

/**********************************************************//**
    Add a state and control to the trajectory

    \param[in,out] traj - trajectory
    \param[in]     s    - state
    \param[in]     u    - control
**************************************************************/
int trajectory_add(struct Trajectory ** traj,
                   struct State * s,
                   struct Control * u)
{
    if (*traj == NULL){
        (*traj) = trajectory_alloc();
        (*traj)->s = state_copy(s);
        (*traj)->u = control_copy(u);
        (*traj)->next = NULL;
    }
    else{
        trajectory_add(&((*traj)->next),s,u);
    }
    return 0;
}

/**********************************************************//**
    Add a state and control to the trajectory by reference

    \param[in,out] traj - trajectory
    \param[in]     s    - state
    \param[in]     u    - control
**************************************************************/
int trajectory_add_ref(struct Trajectory ** traj,
                       struct State * s,
                       struct Control * u)
{
    if (*traj == NULL){
        (*traj) = trajectory_alloc();
        (*traj)->s = s;
        (*traj)->u = u;
        (*traj)->next = NULL;
    }
    else{
        trajectory_add_ref(&((*traj)->next),s,u);
    }
    return 0;
}

/**********************************************************//**
    Print the trajectory

    \param[in,out] traj - trajectory
    \param[in]     fp   - stream to print to
    \param[in]     u    - precision with which to print
**************************************************************/
void trajectory_print(struct Trajectory * traj,
                      FILE *fp,
                      int prec)
{
    
    if (traj == NULL){
//        fprintf(fp,"NULL\n");
        fprintf(fp,"\n");
    }
    else{
        //fprintf(fp, "State: ");
        state_print(traj->s,fp,1,prec);
        //fprintf(fp, "Control: ");
        control_print(traj->u,fp,prec);
        fprintf(fp,"\n");
        trajectory_print(traj->next,fp,prec);
    }
}

/**********************************************************//**
    Get a reference to the last state
**************************************************************/
struct State * trajectory_last_state(struct Trajectory * traj)
{
    if (traj == NULL){
        return NULL;
    }
    else{
        while (traj->next != NULL){
            traj = traj->next;
        }
        return traj->s;
    }
}

/**********************************************************//**
    Take a step of a trajectory

    \param[in,out] traj   - trajectory
    \param[in]     pol    - policy
    \param[in]     dyn    - dynamics
    \param[in]     dt     - time step
                            Used when method=
                            "euler" or 
                            "euler-maruyama"
    \param[in]     method - integration algorithm
    \param[in]     space  - free space to use for computation
    \param[in]     noise  - noise arguments
    \param[in]     args   - additional arguments

    \note
    Let d denote dimension of state \n
    Let dw denote dimension of noise \n
    Method "euler" - space (d), \n
                     noise NULL, \n
                     args NULL  \n
    Method "euler-maruyama" - space(d + d*dw), \n
                              noise (dw) gaussian samples, \n 
                              args NULL \n
    Method "markov-chain" (kushner 2001) \n
                          - space (2d+1) \n
                          - noise (1) uniform (0,1) sample \n
                          - args TensorMM
**************************************************************/
int trajectory_step(struct Trajectory * traj,
                    struct Policy * pol, 
                    struct Dyn * dyn,
                    double dt, char * method, 
                    double * space,
                    void * noise, void * args)
{
    
    if (traj == NULL){
        fprintf(stderr,"Warning: cannot advance trajectory starting\n");
        fprintf(stderr,"         from NULL Trajectory\n");
        return 1;
    }
    else if (pol == NULL){
        fprintf(stderr,"Warning: cannot advance trajectory starting\n");
        fprintf(stderr,"         from NULL Policy\n");
        return 1;
    }
    else if (dyn == NULL){
        fprintf(stderr,"Warning: cannot advance trajectory starting\n");
        fprintf(stderr,"         from NULL Dyn\n");
        return 1;
    }
    assert (dt > 0);

    struct State * current_state = trajectory_last_state(traj);
    if (current_state == NULL){
        return 1;
    }

    size_t d = state_getd(current_state);
    double t = state_gett(current_state);
    double * x = state_getx_ref(current_state);
    if (x == NULL){
        return 1;
    }

    struct Control * u = NULL;
    int res = policy_eval(pol,t,x,&u);
    if (res != 0){
        control_free(u);
        return res;
    }
    struct State * s = NULL;
    if (strcmp(method,"euler") == 0){
        s = euler_step(current_state,u,dt,
                       dyn,space);
    }
    else if (strcmp(method,"euler-maruyama") == 0){
        // args is a realization of the noise
        // should be generated with variance dt
        s = euler_maruyama_step(current_state,noise,u,dt,
                                dyn,space,space+d);
    }
    else if (strcmp(method,"markov-chain") == 0){
        struct TensorMM * mm = args;
        double no = *(double *)noise;
        s = tensor_mm_step(mm,current_state,no,u,
                           &dt,space);
    }
    else{
        return 1;
    }

    res = trajectory_add_ref(&traj,s,u);
    return res;

}

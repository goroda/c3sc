#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c3.h"
#include "simulate.h"
#include "control.h"
#include "integrate.h"

void lin_transform_eval(struct LinTransform * lt, double * x, 
                        double * out)
{
    for (size_t ii = 0; ii < lt->d; ii++){
        out[ii] = lt->slope[ii]*x[ii] + lt->offset[ii];
    }
}

double lin_transform_get_slopei(struct LinTransform * lt, size_t ii)
{
    return lt->slope[ii];
}

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

void state_free(struct State *s)
{
    if (s != NULL){
        free(s->x); s->x = NULL;
        free(s); s = NULL;
    }
}

void state_init(struct State * s,size_t d, double time, double * x)
{
    assert (x != NULL);
    s->d = d;
    s->t = time;
    s->x = calloc_double(d);
    memmove(s->x,x,d*sizeof(double));
    s->spec = NULL;
}

void state_init_zero(struct State * s,size_t d, double time)
{
    s->d = d;
    s->t = time;
    s->x = calloc_double(d);
    s->spec = NULL;
}

struct State * state_copy(struct State * ss)
{
    struct State * s = state_alloc();
    state_init(s,ss->d,ss->t,ss->x);
    return s;
}


void state_up(struct State * s,size_t d, double time, double * x)
{
    s->d = d;
    s->t = time;
    assert (s->x != NULL);
    memmove(s->x,x,d*sizeof(double));
    s->spec = NULL;
}

size_t state_getd(struct State * s)
{
    return s->d;
}

double state_gett(struct State * s)
{
    return s->t;
}

double * state_getx_ref(struct State * s)
{
    return s->x;
}

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

void control_free(struct Control *s)
{
    if (s != NULL){
        free(s->u); s->u = NULL;
        free(s); s = NULL;
    }
}

void control_init(struct Control * s, size_t d, double * u)
{
    assert (u != NULL);
    s->d = d;
    s->u = calloc_double(d);
    memmove(s->u,u,d*sizeof(double));
}

struct Control * control_copy(struct Control * ss)
{
    struct Control * s = control_alloc();
    control_init(s,ss->d,ss->u);
    return s;
}

void control_up(struct Control * s,size_t d, double * u)
{
    assert (s->u != NULL);
    s->d = d;
    memmove(s->u,u,d*sizeof(double));
}

size_t control_getd(struct Control * s)
{
    return s->d;
}

double * control_getu_ref(struct Control * s)
{
    return s->u;
}

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

void trajectory_free(struct Trajectory * traj)
{
    if (traj != NULL){
        state_free(traj->s); traj->s = NULL;
        control_free(traj->u); traj->u = NULL;
        trajectory_free(traj->next); traj->next = NULL;
        free(traj); traj = NULL;
    }
}

int trajectory_add(struct Trajectory ** traj, struct State * s,
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

int trajectory_add_ref(struct Trajectory ** traj, struct State * s,
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

void trajectory_print(struct Trajectory * traj, FILE *fp, int prec)
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


int trajectory_step(struct Trajectory * traj, struct Policy * pol, 
                    struct Dyn * dyn, double dt, char * method, 
                    double * space, void * args)
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
        s = euler_maruyama_step(current_state,args,u,dt,
                                dyn,space,space+d);
    }
    else{
        return 1;
    }

    res = trajectory_add_ref(&traj,s,u);
    return res;

}

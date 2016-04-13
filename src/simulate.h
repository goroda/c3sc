#ifndef SIMULATE_H
#define SIMULATE_H

#include <stdlib.h>
#include <stdio.h>

#include "c3sc_elements.h"
#include "control.h"

struct State * state_alloc();
void state_free(struct State *);

void state_init(struct State *,size_t,double,double *);
void state_init_zero(struct State *,size_t, double);
struct State * state_copy(struct State *);
void state_up(struct State *,size_t,double,double *);

size_t state_getd(struct State *);
double state_gett(struct State *);
double * state_getx_ref(struct State *);
void state_print(struct State *, FILE *, int, int);

//////////////////////////////////////////////


struct Control * control_alloc();
void control_free(struct Control *);
void control_init(struct Control *,size_t, double *);
struct Control * control_copy(struct Control *);
void control_up(struct Control *,size_t, double *);
size_t control_getd(struct Control *);
double * control_getu_ref(struct Control *);

//////////////////////////////////////////

/* struct Trajectory * trajectory_alloc(); */
/* void trajectory_free(struct Trajectory *); */
/* int trajectory_add(struct Trajectory **, struct State *,  */
/*                    struct Control *); */
/* int trajectory_add_ref(struct Trajectory **, struct State *, */
/*                        struct Control *); */
/* struct State * trajectory_last_state(struct Trajectory *); */
/* int trajectory_step(struct Trajectory *, struct Policy *,  */
/*                     struct Dyn *, double, char *,  */
/*                     double *, void *, void *); */
/* void trajectory_print(struct Trajectory *, FILE *, int); */

#endif  

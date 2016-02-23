#ifndef CONTROL_H
#define CONTROL_H


#include "c3sc_elements.h"

struct Policy * policy_alloc();
void policy_init(struct Policy *, size_t, size_t, 
                 double *, double *);
void policy_free(struct Policy *);

void policy_add_feedback(struct Policy *, 
                         int (*f)(double,double*,double*));
int policy_eval(struct Policy *, double, double *, struct Control **);

#endif









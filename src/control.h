#ifndef CONTROL_H
#define CONTROL_H


#include "c3sc_elements.h"

struct Policy;
struct Policy * policy_alloc(size_t,size_t);
void policy_free(struct Policy *);

void policy_set_bounds(struct Policy *,double *, double *);
size_t * policy_get_ranks(struct Policy *,size_t);
void policy_init_discrete(struct Policy *,size_t *,double **);
void policy_approx(struct Policy *,
                   double (*)(double *,size_t, void *),
                   void *, int);

void 
policy_add_feedback(struct Policy *, 
                    int (*f)(double,double*,double*,void*),
                    void*);
int policy_eval(struct Policy *, double, double *, struct Control **);

#endif









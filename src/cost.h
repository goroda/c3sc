#ifndef C3SC_COST_H
#define C3SC_COST_H

#include "c3.h"

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
 */
struct Cost 
{
    size_t d;
    struct BoundingBox * bds;
    struct FunctionTrain * cost;

    size_t * N;
    double ** x;
};

struct Cost * cost_alloc(size_t, double *, double *);
struct Cost * cost_copy_deep(struct Cost *);
void cost_free(struct Cost *);
double * cost_get_lb(struct Cost *);
double * cost_get_ub(struct Cost *);
size_t * cost_get_ranks(struct Cost *);

void cost_init_discrete(struct Cost *,size_t *,double **);
void cost_approx(struct Cost *,
                 double (*)(double *, void *),
                 void *, int,double,double,size_t);

int cost_eval(struct Cost *, double, double *, double *);
//void cost_eval_bb(size_t,double *,double **,double *,void *);
double cost_eval_bb(double,double *,void *);
int cost_eval_neigh(struct Cost *, double, double *,
                    size_t, double[],double []);


#endif

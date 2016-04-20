#ifndef C3SC_COST_H
#define C3SC_COST_H

#include "c3.h"
#include "util.h"

struct Cost;
struct Cost * cost_alloc(size_t, double *, double *);
struct Cost * cost_copy_deep(struct Cost *);
int cost_save(struct Cost *, char *);
int cost_load(struct Cost *, char *);
void cost_free(struct Cost *);
size_t cost_get_d(struct Cost *);
double * cost_get_lb(struct Cost *);
double * cost_get_ub(struct Cost *);
size_t * cost_get_ranks(struct Cost *);
double cost_norm2(struct Cost *);
double cost_norm2_diff(struct Cost *,struct Cost *);

void cost_init_discrete(struct Cost *,size_t *,double **);
void cost_add_nodes(struct Cost *, double *, double *, size_t);
size_t cost_get_size(const struct Cost *);
void cost_interpolate_new(struct Cost *, struct Cost *);
void cost_approx(struct Cost *,
                 double (*)(double *, void *),
                 void *,int,const struct ApproxArgs *);

int cost_eval(struct Cost *, double, double *, double *);
//void cost_eval_bb(size_t,double *,double **,double *,void *);
double cost_eval_bb(double,double *,void *);
int cost_eval_neigh(struct Cost *, double, double *,
                    size_t, double[],double []);


#endif

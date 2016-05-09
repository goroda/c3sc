#ifndef C3SC_COST_H
#define C3SC_COST_H

#include "c3.h"
#include "util.h"

struct Cost;
struct Cost * cost_alloc(size_t, double *, double *);
struct Cost * cost_copy_deep(const struct Cost *);
int cost_save(const struct Cost *, char *);
int cost_load(struct Cost *, char *);
void cost_free(struct Cost *);
size_t cost_get_d(const struct Cost *);
size_t cost_get_size(const struct Cost *);
double * cost_get_lb(const struct Cost *);
double * cost_get_ub(const struct Cost *);
void cost_get_h(const struct Cost *, double *);
size_t * cost_get_ranks(const struct Cost *);
double cost_norm2(const struct Cost *);
double cost_norm2_diff(const struct Cost *,const struct Cost *);

void cost_init_grid(struct Cost *,size_t *,double **);
/* void cost_init_discrete(struct Cost *,size_t *,double **); */
void cost_add_nodes(struct Cost *, double *, double *);

void cost_interpolate_new(struct Cost *, const struct Cost *);
void cost_interp_inhalf(struct Cost *, int);
void cost_approx(struct Cost *,
                 double (*)(const double *, void *),
                 void *,int,const struct ApproxArgs *);

int cost_eval(struct Cost *, double, const double *, double *);
//void cost_eval_bb(size_t,double *,double **,double *,void *);
double cost_eval_bb(double,const double *,void *);
/* int cost_eval_neigh(struct Cost *, double, double *, */
/*                     size_t, double[],double []); */
int cost_eval_neigh(struct Cost *,
                    double,
                    double *,
                    double *,
                    double *,
                    double *);

#endif

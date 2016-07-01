#ifndef BELLMAN_H
#define BELLMAN_H

#include <stdlib.h>

double bellmanrhs(size_t dx, size_t du, double stage_cost, const double * stage_grad, 
                      double discount, const double * prob, const double * prob_grad, 
                      double dt, const double * dtgrad, const double * cost, double * grad);

struct MCAparam;
struct MCAparam * mca_param_create(size_t, size_t);
void mca_add_grid_refs(struct MCAparam *, size_t *, double **,
                       double, double *);
void mca_param_destroy(struct MCAparam *);

struct DPparam;
struct DPparam * dp_param_create(size_t, size_t, size_t, double );
void dp_param_destroy(struct DPparam *);
void dp_param_add_drift(struct DPparam *,
                        int (*)(double,const double*,const double*,
                                double*,double*,void*),
                        void *);
void dp_param_add_diff(struct DPparam *,
                       int (*)(double,const double*,const double*,
                                double*,double*,void*),
                       void *);
void dp_param_add_stagecost(struct DPparam *,
                            int (*)(double,const double*,
                                    const double*,double*,double*));
void dp_param_add_boundcost(struct DPparam *, int (*)(double,const double*,double*));
void dp_param_add_obscost(struct DPparam *, int (*)(const double*,double*));


struct ControlParams;
struct ControlParams *
control_params_create(size_t, size_t, struct DPparam *,
                      struct MCAparam *, struct c3Opt *);
void control_params_add_state_info(struct ControlParams *,
                                   double, const double *, int,
                                   const double *);
int control_params_get_last_res(const struct ControlParams *);
void control_params_destroy(struct ControlParams *);

///////////////////////////////////////////////////////////////
double bellman_control(size_t, double *, double *, void *);
int bellman_optimal(size_t, double *, double *, void *);
#endif

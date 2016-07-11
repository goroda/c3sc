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
void dp_param_add_boundary(struct DPparam *, struct Boundary *);
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


///////////////////////////////////////////////////////////////
struct VIparam;
struct VIparam * vi_param_create(double);
void vi_param_destroy(struct VIparam *);
void vi_param_add_cp(struct VIparam *, struct ControlParams *);
void vi_param_add_value(struct VIparam * vi, struct ValueF *);
int bellman_vi(size_t, const double *, double *, void *);


///////////////////////////////////////////////////////////////
struct PIparam;
struct PIparam * pi_param_create(double, struct ValueF *);
void pi_param_destroy(struct PIparam *);
void pi_param_add_cp(struct PIparam * pi, struct ControlParams *);
void pi_param_add_value(struct PIparam *, struct ValueF *);
int bellman_pi(size_t, const double *, double *, void *);


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

struct Diag;
void diag_destroy(struct Diag **);
struct Diag * diag_create(size_t iter, int type, double norm,
                          double abs_diff, size_t dim, size_t * ranks,
                          double frac);
void diag_append(struct Diag ** diag,
                 size_t iter, int type, double norm,
                 double abs_diff, size_t dim, size_t * ranks,
                 double frac);
void diag_print(struct Diag * head, FILE * fp);
int diag_save(struct Diag * head, char * filename);

struct C3Control *
c3control_create(size_t, size_t, size_t,
                     double *, double *,
                     size_t *, double);

void c3control_destroy(struct C3Control *);
size_t * c3control_get_ngrid(struct C3Control *);
double ** c3control_get_xgrid(struct C3Control *);
void c3control_add_policy_sim(struct C3Control *, struct ValueF *, 
                              struct c3Opt * opt_sim,
                              void (*transform)(size_t, const double *, double*));

void c3control_set_external_boundary(struct C3Control *, size_t,
                                     char *);
void c3control_add_obstacle(struct C3Control *, double *, double *);

void c3control_add_drift(struct C3Control * c3c, 
                             int (*b)(double,const double*,const double*,
                                          double*,double*,void*),
                             void * args);
void c3control_add_diff(struct C3Control * c3c, 
                            int (*s)(double,const double*,const double*,
                                         double*,double*,void*),
                            void * sargs);
void c3control_add_stagecost(struct C3Control * c3c,
                             int (*stagecost)(double,const double*,
                                                  const double*,double*,double*));
void c3control_add_boundcost(struct C3Control * c3c,
                                 int (*boundcost)(double,const double*,double*));

void c3control_add_obscost(struct C3Control * c3c, int (*obscost)(const double*,double*));

int c3control_controller(double,const double *, double *, void *);

struct ValueF * c3control_step_vi(struct C3Control * c3c, struct ValueF * vf,
                                  struct ApproxArgs * apargs,
                                  struct c3Opt * opt,
                                  int verbose,
                                  size_t * nevals);

struct ValueF * c3control_step_pi(struct C3Control * c3c, struct ValueF * vf,
                                  struct PIparam * poli,
                                  struct ApproxArgs * apargs,
                                  struct c3Opt * opt,
                                  int verbose,
                                  size_t * nevals_iter);
struct ValueF *
c3control_init_value(struct C3Control * c3c,
                     int (*f)(size_t,const double *,double*,void*),void * args,
                     struct ApproxArgs * aargs, int verbose);

struct ValueF * c3control_vi_solve(struct C3Control * c3c,
                                   size_t maxiter, double abs_conv_tol,
                                   struct ValueF * vo,
                                   struct ApproxArgs * apargs,
                                   struct c3Opt * opt,
                                   int verbose, 
                                   struct Diag ** diag);

struct ValueF * c3control_pi_solve(struct C3Control * c3c,
                                   size_t maxiter, double abs_conv_tol,
                                   struct ValueF * policy,
                                   struct ApproxArgs * apargs,
                                   struct c3Opt * opt,
                                   int verbose, 
                                   struct Diag **);
#endif

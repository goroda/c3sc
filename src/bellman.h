// Copyright (c) 2015-2016, Massachusetts Institute of Technology
// Copyright (c) 2018, University of Michigan
//
// This file is part of the C3 for Stochastic Optimal Control (C3SC) toolbox
// Author: Alex A. Gorodetsky 
// Contact: goroda@umich.edu
// Website: https://www.alexgorodetsky.com

// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice, 
//    this list of conditions and the following disclaimer in the documentation 
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its contributors 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//Code


#ifndef BELLMAN_H
#define BELLMAN_H

#include <stdlib.h>
#include "c3/lib_linalg.h"
#include "c3/lib_optimization.h"

#include "util.h"

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
                      struct MCAparam *, struct Workspace *, struct c3Opt *);
void control_params_add_time_and_states(struct ControlParams *, double, size_t, const double *);
int control_params_get_last_res(const struct ControlParams *);
void control_params_destroy(struct ControlParams *);

///////////////////////////////////////////////////////////////
double bellman_control(size_t, const double *, double *, void *);
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

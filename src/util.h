// Copyright (c) 2015-2016, Massachusetts Institute of Technology
//
// This file is part of the C3 for Stochastic Optimal Control (C3SC) toolbox
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

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


#ifndef C3SC_UTIL_H
#define C3SC_UTIL_H

/* #include "c3.h" */

#include "c3/array.h"
#include "c3/lib_funcs.h"
#include "hashgrid.h"

struct ApproxArgs;
struct ApproxArgs * approx_args_init(void);
void approx_args_free(struct ApproxArgs *);
void approx_args_set_function_class(struct ApproxArgs *, enum function_class);
enum function_class approx_args_get_function_class(const struct ApproxArgs *);

void approx_args_set_cross_tol(struct ApproxArgs *, double);
double approx_args_get_cross_tol(const struct ApproxArgs *);
void approx_args_set_round_tol(struct ApproxArgs *, double);
double approx_args_get_round_tol(const struct ApproxArgs *);
void approx_args_set_kickrank(struct ApproxArgs *, size_t);
size_t approx_args_get_kickrank(const struct ApproxArgs *);
void approx_args_set_maxrank(struct ApproxArgs *, size_t);
size_t approx_args_get_maxrank(const struct ApproxArgs *);
void approx_args_set_startrank(struct ApproxArgs *, size_t);
size_t approx_args_get_startrank(const struct ApproxArgs *);
void approx_args_set_adapt(struct ApproxArgs *, int);
int approx_args_get_adapt(const struct ApproxArgs *);

int c3sc_check_bounds(size_t,double*,double*,const double*);
size_t c3sc_sample_discrete_rv(size_t,double *,double);
double * c3sc_combine_and_sort(size_t,double *,size_t,double *,size_t *);

struct HashGrid;
struct HashGrid * hash_grid_create(size_t);
struct HashGrid * hash_grid_create_grid(size_t,const struct c3Vector *);
struct HashGrid ** hash_grid_create_ndgrid(size_t, size_t, struct c3Vector **);
void hash_grid_free_ndgrid(size_t d, struct HashGrid **);
int hash_grid_add_element(struct HashGrid *,size_t,double);
void hash_grid_print(struct HashGrid *, FILE *);
size_t hash_grid_get_ind(struct HashGrid *,double,int *);
int hash_grid_ndgrid_get_ind(struct HashGrid **, size_t, const double *, size_t *);
void hash_grid_free(struct HashGrid *);

size_t uniform_stride(size_t, size_t);


struct Workspace;
struct Workspace * workspace_alloc(size_t,size_t,size_t,size_t);
void workspace_free(struct Workspace *);

void workspace_increment_vi_iter(struct Workspace * w);
size_t workspace_get_vi_iter(const struct Workspace * w);
struct HTable * workspace_get_vi_htable(const struct Workspace * w);

void workspace_increment_pi_iter(struct Workspace * w);
size_t workspace_get_pi_iter(const struct Workspace * w);
void workspace_increment_pi_subiter(struct Workspace * w);
size_t workspace_get_pi_subiter(const struct Workspace * w);

struct HTable * workspace_get_pi_htable(const struct Workspace *);
struct HTable * workspace_get_pi_prob_htable(const struct Workspace * w);


void workspace_set_active(struct Workspace * w, size_t active);
size_t workspace_get_active(struct Workspace * w);
double * workspace_get_drift(struct Workspace * w, size_t node);
double * workspace_get_grad_drift(struct Workspace * w, size_t node);
double * workspace_get_diff(struct Workspace * w, size_t node);
double * workspace_get_grad_diff(struct Workspace * w, size_t node);
double * workspace_get_dt(struct Workspace * w, size_t node);
double * workspace_get_grad_dt(struct Workspace * w, size_t node);
double * workspace_get_prob(struct Workspace * w, size_t node);
double * workspace_get_grad_prob(struct Workspace * w, size_t node);
double * workspace_get_grad_stage(struct Workspace * w, size_t node);
double * workspace_get_control_size_extra(struct Workspace * w, size_t node);
double * workspace_get_u(struct Workspace * w, size_t node);
double * workspace_get_costs(struct Workspace * w, size_t node);
int *    workspace_get_absorbed(struct Workspace * w, size_t node);
size_t * workspace_get_ind_to_serialize(struct Workspace * w);
size_t * workspace_get_absorbed_no(struct Workspace *);
size_t * workspace_get_absorbed_yes(struct Workspace *);
#endif

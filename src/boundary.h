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


#ifndef C3SC_BOUNDARY_H
#define C3SC_BOUNDARY_H

enum EBTYPE {
    EB_NONE=0,
    ABSORB=1,
    PERIODIC=2,
    REFLECT=3
};

struct Boundary;
struct Boundary * boundary_alloc(size_t,double *,double*);
double * boundary_obstacle_get_lb(struct Boundary *, size_t);
double * boundary_obstacle_get_ub(struct Boundary *, size_t);

size_t boundary_get_nobs(struct Boundary *);
struct Boundary * boundary_copy_deep(struct Boundary *);
void boundary_free(struct Boundary *);
void boundary_external_set_type(struct Boundary *,size_t,char *);
double outer_bound_dim(const struct Boundary *, size_t, 
                       double, int *);
enum EBTYPE boundary_type_dim(const struct Boundary *,size_t, int);
int boundary_in_obstacle(const struct Boundary *, const double *);

struct BoundInfo * boundary_type(const struct Boundary *,double,const double *);

enum BOUNDRESULT {
    IN,
    LEFT,
    RIGHT,
};

struct BoundInfo;
struct BoundInfo * bound_info_alloc(size_t);
void bound_info_free(struct BoundInfo *);
void boundary_add_obstacle(struct Boundary *, double *, double *);
int bound_info_set_dim(struct BoundInfo *, enum BOUNDRESULT,
                       enum EBTYPE, size_t);
int bound_info_set_xmap_dim(struct BoundInfo *, double, size_t);
int bound_info_onbound(const struct BoundInfo *);
int bound_info_onbound_dim(const struct BoundInfo *,size_t);
int bound_info_absorb(const struct BoundInfo *);
int bound_info_period(const struct BoundInfo *);
int bound_info_period_dim_dir(const struct BoundInfo *,size_t);
double bound_info_period_xmap(const struct BoundInfo *,size_t);
int bound_info_reflect(const struct BoundInfo *);
int bound_info_reflect_dim_dir(const struct BoundInfo *,size_t);

int bound_info_get_in_obstacle(const struct BoundInfo *);
#endif

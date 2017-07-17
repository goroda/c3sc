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


#ifndef NODEUTIL_H
#define NODEUTIL_H

#include <stdlib.h>
#include "boundary.h"
#include "valuefunc.h"
int transition_assemble_old(size_t dx, size_t du, size_t dw, double h, const double * hvec,
                        const double * drift, const double * grad_drift,
                        const double * ddiff, const double * grad_ddiff,
                        double * prob, double * grad_prob,
                        double * dt,   double * grad_dt, double * space);

int transition_assemble(size_t dx, size_t du, size_t dw, double h, const double * hvec,
                        const double * drift, const double * grad_drift,
                        const double * ddiff, const double * grad_ddiff,
                        double * prob, double * grad_prob,
                        double * dt,   double * grad_dt, double * space);

int process_fibers(size_t d, size_t nvals, const double * x,
                   int * boundv, int * neighbors, 
                   struct Boundary * bound);

int process_fibers_neighbor(size_t d, const size_t * fixed_ind, size_t dim_vary, 
                            const double * x, int * absorbed, size_t * neighbors_vary,
                            size_t * neighbors_fixed, const size_t * ngrid, 
                            const struct Boundary * bound); 

int convert_fiber_to_ind(size_t d, size_t N, const double * x, 
                         const size_t * Ngrid, double ** xgrid,
                         size_t * fixed_ind, size_t * dim_vary);

int mca_get_neighbor_costs(size_t d,size_t N,const double * x,struct Boundary * bound,
                           struct ValueF * vf, const size_t * ngrid, double ** xgrid,
                           size_t *, size_t *,
                           int * absorbed, double * out);

int mca_get_neighbor_node_costs(size_t d, const double * x,
                                struct Boundary * bound,
                                struct ValueF * vf, 
                                const size_t * ngrid, 
                                double ** xgrid,
                                int * absorbed, double * out);

#endif

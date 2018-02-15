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


#ifndef DYNAMICS_H
#define DYNAMICS_H

struct Drift;
struct Diff;
struct Dyn;

struct Drift * drift_alloc(size_t, size_t);
struct Drift * drift_copy(struct Drift *);

void drift_free(struct Drift *);
void drift_add_func(struct Drift *,
                    int (*)(double,const double*,const double*,
                            double*,double*,void*),
                    void *);
size_t drift_get_dx(struct Drift *);
int drift_eval(struct Drift *,double,const double*,
               const double*,double*, double*);

////////////////////////////////////
struct Diff * diff_alloc(size_t, size_t, size_t);
struct Diff * diff_copy(struct Diff *);
void diff_free(struct Diff *);
void diff_add_func(struct Diff *,
                   int (*)(double,const double*,const double*,
                           double*,double*,void*),
                   void *);
int diff_eval(struct Diff *,double,const double*,const double*,
              double*,double*);
size_t diff_get_dw(struct Diff *);


///////////////////////////////////////////////////////////

struct Dyn * dyn_alloc(struct Drift *, struct Diff *);
struct Dyn * dyn_copy_deep(struct Dyn *);
void dyn_free(struct Dyn *);
void dyn_free_deep(struct Dyn *);

void dyn_init_ref(struct Dyn*,struct Drift *, struct Diff *);
size_t dyn_get_dx(struct Dyn *);
size_t dyn_get_dw(struct Dyn *);
size_t dyn_get_du(struct Dyn *);
int dyn_eval(struct Dyn *,double,const double*,const double*,
             double*,double*,double*,double*);

////////////////////////////////////////////////////////////

#endif

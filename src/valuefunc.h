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


#ifndef VALUE_H
#define VALUE_H

#include <stdlib.h>

#include "util.h"

struct ValueF;
void valuef_destroy(struct ValueF *);
struct ValueF * valuef_copy(struct ValueF *);
struct CrossIndex ** valuef_get_isl(const struct ValueF *);
int valuef_save(struct ValueF *, char *);
struct ValueF * valuef_load(char *, size_t *, double **);
int valuef_savetxt(struct ValueF *, char *);
struct ValueF * valuef_loadtxt(char *,size_t *,double **);

size_t * valuef_get_ranks(struct ValueF *);
double valuef_norm(struct ValueF *);
double valuef_norm2diff(struct ValueF *, struct ValueF *);
double valuef_eval(struct ValueF *, const double *);


int
valuef_eval_fiber_ind_nn(struct ValueF *, const size_t *,
                         size_t, const size_t *, const size_t *,
                         double *);


/* struct ValueF *  */
/* valuef_interp(size_t,int (*)(size_t,const double *,double*,void*),void *, */
/*               const size_t *, double **, double **, */
/*               struct ApproxArgs *, int); */

struct ValueF * 
valuef_interp(size_t,int (*)(size_t,const double *,double*,void*),void *,
              const size_t *, double **, struct ValueF*,
              struct ApproxArgs *, int);
#endif

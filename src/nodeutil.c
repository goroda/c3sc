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


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "c3.h"

#include "boundary.h"
#include "valuefunc.h"

/**********************************************************//**
    Assemble transition probabilities

    \param[in]     dx         - dimension of state space
    \param[in]     du         - dimension of control space
    \param[in]     dw         - dimension of stochastic space
    \param[in]     h          - minimum step-size in each direction
    \param[in]     hvec       - step-size in each direction
    \param[in]     drift      - drift evaluation
    \param[in]     grad_drift - gradient of drift
    \param[in]     ddiff      - diagonal diffusion
    \param[in]     grad_ddiff - gradient of diffusion
    \param[in,out] prob       - vector of probabililities (2d+1)
    \param[in,out] grad_prob  - gradient of probabilities (du*(2d+1))
    \param[in,out] dt         - interpolation interval
    \param[in,out] grad_dt    - gradient of interpolation interval
    \param[in,out] space      - workspace needed if gradients requested (du,) array

    \return 0 if everything is ok, 1 if computed dt is infinity, 2 if derivative ambiguous

    \note
    The probabilities are ordered as ("left", "right") for each dimension
    followed by the transition probability to one-self

    gradient of drift is dimension (dx,du) and stored as x varying first
    (df1/du1 ..... df1/dun
     .                .
     dfm/du1  ..... dfm/dun)
    
**************************************************************/
int transition_assemble(size_t dx, size_t du, size_t dw, double h, double * hvec,
                        const double * drift, const double * grad_drift,
                        const double * ddiff, const double * grad_ddiff,
                        double * prob, double * grad_prob,
                        double * dt, double * grad_dt,
                        double * space)
{

    if (space != NULL){
        assert (grad_drift != NULL);
        assert (grad_ddiff != NULL);
        for (size_t ii = 0; ii < du; ii++){
            space[ii] = 0.0;
        }
    }

    double normalization = 0.0;
    double h2 = h * h;
    double t2, diff, t;
    int res = 0;
    for (size_t ii = 0; ii < dx; ii++){
        /* printf("ii = %zu\n",ii); */
        t = h2/hvec[ii];
        diff = ddiff[ii*dx+ii] * ddiff[ii*dx+ii];
        t2 = h2 / hvec[ii] / hvec[ii];

        prob[2*ii] = t2 * diff / 2.0;
        prob[2*ii+1] = t2 * diff / 2.0;
        if (drift[ii] < -1e-14){ // add transition to the left
            prob[2*ii] -= t * drift[ii];
        }
        else if (drift[ii] > 1e-14){ // add transition to the right
            prob[2*ii+1] += t * drift[ii];
        }

        /* printf("ii=%zu (-,+) = (%G,%G)\n",ii,prob[2*ii],prob[2*ii+1]); */
        normalization += prob[2*ii];
        normalization += prob[2*ii+1];

        if (grad_prob != NULL){
            for (size_t jj = 0; jj < du; jj++){
                grad_prob[2*ii*du+jj] = 0.0;
                grad_prob[(2*ii+1)*du+jj] = 0.0;
            }
            // add diffusions
            cblas_daxpy(du,t2,grad_ddiff+ii*dx+ii,dx*dw,grad_prob + 2*ii*du,1);
            cblas_daxpy(du,t2,grad_ddiff+ii*dx+ii,dx*dw,grad_prob + (2*ii+1)*du,1);
            
            // now drifts
            if (drift[ii] < -1e-14){
                cblas_daxpy(du,-t,grad_drift+ii,dx,grad_prob + (2*ii)*du,1);
            }
            else if (drift[ii] > 1e-14){
                cblas_daxpy(du,t,grad_drift+ii,dx,grad_prob + (2*ii+1)*du,1);
            }
            else{
                /* printf("warning! not sure what to do here!\n"); */
                for (size_t jj = 0; jj < du; jj++){
                    /* printf("jj = %zu\n",jj); */
                    if (grad_drift[jj*dx+ii] < 0){ // wants to go left with changing u
                        cblas_daxpy(du,-t,grad_drift+ii,dx,grad_prob + (2*ii)*du,1);
                    }
                    else if (grad_drift[jj*dx+ii] > 0){ // wants to go right with changing u
                        cblas_daxpy(du,t,grad_drift+ii,dx,grad_prob + (2*ii+1)*du,1);
                    }
                    else{ // do nothing because doesn't want to change!
                        res = 2;
                        /* printf("warning! not sure what to do here transition assemble!\n"); */
                    }

                }
                // don't know what to do about gradients if it is zero
            }
            /* printf("\t (g-,g+) = (%G,%G)\n",grad_prob[2*ii],grad_prob[2*ii+1]); */

            // now the normalizations
            cblas_daxpy(du,1,grad_prob+(2*ii  )*du,1,space,1);
            cblas_daxpy(du,1,grad_prob+(2*ii+1)*du,1,space,1);
        }
    }

    if (normalization < 1e-14){ // basically stationarry
        return 1;
    }
    
    *dt = h2 / normalization;

    prob[2*dx] = 1.0;    
    if (grad_prob != NULL){
        double norm2 = normalization * normalization;
        for (size_t jj = 0; jj < du; jj++){
            grad_prob[2*dx*du+jj] = 0.0;            
            grad_dt[jj] = -space[jj]*h2 / norm2;
        }

        for (size_t ii = 0; ii < 2*dx; ii++){

            for (size_t jj = 0; jj < du; jj++){
                grad_prob[ii*du+jj] = (normalization * grad_prob[ii*du+jj] - 
                                       space[jj]*prob[ii]) / norm2;
            }
            /* printf("\t (ii,g) = (%zu,%G)\n",ii,grad_prob[ii]); */
            prob[ii] /= normalization;
            prob[2*dx] -= prob[ii];
            cblas_daxpy(du,-1.0,grad_prob+ii*du,1,grad_prob + 2*dx*du,1);
        }
    }
    else{
        for (size_t ii = 0; ii < dx; ii++){
            prob[2*ii]   /= normalization;
            prob[2*ii+1] /= normalization;
            prob[2*dx]   -= prob[2*ii];
            prob[2*dx]   -= prob[2*ii+1];
        }
    }

    return res;
}

size_t convert_x_to_ind(double x, size_t N, double * grid)
{
    double tol = 1e-14;
    // slow O(N)
    for (size_t ii = 0; ii < N; ii++){
        double diff = fabs(x- grid[ii]);
        if (diff < tol){
            return ii;
        }
    }
    return N;
}

/**********************************************************//**
    Convert double-based fibers obtained during cross
    to index based

    \param[in]     d         - dimension of state space
    \param[in]     N         - number of evaluations
    \param[in]     x         - evaluation locations
    \param[in]     Ngrid     - size of each dimensions discretization
    \param[in]     xgrid     - grid-values
    \param[in,out] fixed_ind - (d,) vector of fixed indices
    \param[in,out] dim_vary  - dimension that is varying

    \return 0 if everything is ok, 
            1 if evaluation location is not on the grid
            2 if number of evaluations does not match up with discretization size
**************************************************************/
int convert_fiber_to_ind(size_t d, size_t N, const double * x,
                         const size_t * Ngrid, double ** xgrid,
                         size_t * fixed_ind, size_t * dim_vary)
{
    
    // convert first one
    for (size_t ii = 0; ii < d; ii++){
        fixed_ind[ii] = convert_x_to_ind(x[ii],Ngrid[ii],xgrid[ii]);
        if (fixed_ind[ii] == Ngrid[ii]){
            return 1;
        }
    }

    // now convert the second one to see which index is different
    size_t temp;
    *dim_vary = d;
    for (size_t ii = 0; ii < d; ii++){
        temp = convert_x_to_ind(x[ii+d],Ngrid[ii],xgrid[ii]);
        if (temp != fixed_ind[ii]){
            *dim_vary = ii;
            break;
        }
    }
    if (*dim_vary == d){
        return 1;
    }
    else if (N != Ngrid[*dim_vary]){
        return 2;
    }

    return 0;
}


/**********************************************************//**
    Convert double-based fibers obtained during cross
    to index based

    \param[in]     d               - dimension of state space
    \param[in]     fixed_ind       - (d,) vector of fixed indices
    \param[in]     dim_vary        - dimension that varies
    \param[in]     x               - values of the fiber
    \param[in,out] absorbed        - (0,no), (1,yes), (-1,obstacle)
    \param[in,out] neighbors_vary  - (2*ngrid[dim_vary],)  neighbor indices in varying dimension
    \param[in,out] neighbors_fixed - (2*(d-1),) neighbor indices for fixed dimension
    \param[in]     ngrid           - size of grid in each dimension
    \param[in]     bound           - boundary information

    \return 0 if everything is ok, 
**************************************************************/
int process_fibers_neighbor(size_t d, const size_t * fixed_ind, size_t dim_vary, 
                            const double * x, int * absorbed, size_t * neighbors_vary,
                            size_t * neighbors_fixed, const size_t * ngrid, 
                            const struct Boundary * bound)
{
    assert (bound != NULL);
    /* printf("process fibers N=%zu\n",ngrid[dim_vary]); */
    for (size_t ii = 0; ii < ngrid[dim_vary]; ii++){
        /* dprint(d,x+ii*d); */
        /* printf("ii = %zu\n",ii); */
        int inobs = boundary_in_obstacle(bound,x+ii*d);
        /* printf("inobs = %d\n",inobs); */
        if (inobs == 0){
            absorbed[ii] = 0;
        }
        else{
            /* printf("x is in an obstacle!! = "); dprint(d,x+ii*d); */
            /* assert(1 == 0); */
            absorbed[ii] = -1;
        }
    }

    /* printf("kk\n"); */
    size_t onind = 0;
    for (size_t ii = 0; ii < d; ii++){
        if (ii != dim_vary){
            if ( fixed_ind[ii] == 0){ // left boundary
                enum EBTYPE b = boundary_type_dim(bound,ii,0);
                if (b == ABSORB){
                    neighbors_fixed[onind  ] = fixed_ind[ii];
                    neighbors_fixed[onind+1] = fixed_ind[ii];
                    for (size_t jj = 0; jj < ngrid[dim_vary]; jj++){
                        absorbed[jj] = 1;
                    }
                }
                else if (b == REFLECT){
                    neighbors_fixed[onind  ] = fixed_ind[ii];
                    neighbors_fixed[onind+1] = fixed_ind[ii]+1;
                }
                else if (b == PERIODIC){
                    neighbors_fixed[onind  ] = ngrid[ii]-2;
                    neighbors_fixed[onind+1] = fixed_ind[ii]+1;
                }
                else{
                    fprintf(stderr,"No boundary specified!\n");
                    assert(1 == 0);
                }
            } 
            else if (fixed_ind[ii] == (ngrid[ii]-1)){ // right boundary
                enum EBTYPE b = boundary_type_dim(bound,ii,1);
                if (b == ABSORB){
                    neighbors_fixed[onind  ] = fixed_ind[ii];
                    neighbors_fixed[onind+1] = fixed_ind[ii];
                    for (size_t jj = 0; jj < ngrid[dim_vary]; jj++){
                        absorbed[jj] = 1;
                    }
                }
                else if (b == REFLECT){
                    neighbors_fixed[onind  ] = fixed_ind[ii]-1;
                    neighbors_fixed[onind+1] = fixed_ind[ii];
                }
                else if (b == PERIODIC){
                    neighbors_fixed[onind  ] = fixed_ind[ii]-1;
                    neighbors_fixed[onind+1] = 1;
                }
                else{
                    fprintf(stderr,"No boundary specified!\n");
                    assert(1 == 0);
                }
                // handle boundary
            }
            else{
                neighbors_fixed[onind  ] = fixed_ind[ii]-1;
                neighbors_fixed[onind+1] = fixed_ind[ii]+1;
            }
            onind += 2;
        }
    }
    
    /* printf("now\n"); */
    // first point
    enum EBTYPE b = boundary_type_dim(bound,dim_vary,0);
    if (b == ABSORB){
        neighbors_vary[0] = 0;
        neighbors_vary[1] = 0;
        absorbed[0] = 1;
    }
    else if (b == REFLECT){
        neighbors_vary[0] = 0;
        neighbors_vary[1] = 1;
        absorbed[0] = 0;
    }
    else if (b == PERIODIC){
        neighbors_vary[0] = ngrid[dim_vary]-2;
        neighbors_vary[1] = 1;
        absorbed[0] = 0;
    }
    else{
        fprintf(stderr,"Should not be here!\n");
        assert(1 == 0);
    }
    
    /* printf("on last point\n"); */
    // last point
    b = boundary_type_dim(bound,dim_vary,1);
    if (b == ABSORB){
        neighbors_vary[2*(ngrid[dim_vary]-1)] = ngrid[dim_vary]-1;
        neighbors_vary[2*(ngrid[dim_vary]-1)+1] = ngrid[dim_vary]-1;
        absorbed[ngrid[dim_vary]-1] = 1;
    }
    else if (b == REFLECT){
        neighbors_vary[2*(ngrid[dim_vary]-1)] =  ngrid[dim_vary]-2;
        neighbors_vary[2*(ngrid[dim_vary]-1)+1] =  ngrid[dim_vary]-1;
        absorbed[ngrid[dim_vary]-1] = 0;
    }
    else if (b == PERIODIC){
        neighbors_vary[2*(ngrid[dim_vary]-1)] = ngrid[dim_vary]-2;
        neighbors_vary[2*(ngrid[dim_vary]-1)+1] = 1;
        absorbed[ngrid[dim_vary]-1] = 0;
    }
    else{
        fprintf(stderr,"Should not be here!\n");
        assert(1 == 0);
    }
    
    /* printf("lets go\n"); */
    for (size_t ii = 1; ii < ngrid[dim_vary]-1; ii++){
        if (absorbed[ii] == 0){
            neighbors_vary[2*ii] = ii-1;
            neighbors_vary[2*ii+1] = ii+1;
        }
        else{
            neighbors_vary[2*ii] = ii;
            neighbors_vary[2*ii+1] = ii;
        }
    }
    /* printf("finished\n"); */
    return 0;
}

/**********************************************************//**
    Process a fiber of nodes obtaining the boundary information
    for each node and the costs of all nodes neded for the MCA method
    
    \param[in]     d               - dimension of state space
    \param[in]     N               - number of elements in the fiber
    \param[in]     x               - evaluation locations
    \param[in]     bound           - boundary information
    \param[in]     vf              - value function
    \param[in]     ngrid           - size of grid in each dimension
    \param[in]     xgrid           - size of grid in each dimension
    \param[in,out] fixed_ind       - (d,) fixed indices with some fill for varying one
    \param[in,out] dim_vary        - dimension that is varying
    \param[in,out] absorbed        - (0,no), (1,yes), (-1,obstacle)
    \param[in,out] out             - (N*(2d+1),) array of values at neighbors of each element including fibers

    \return 0 if everything is ok, 
**************************************************************/
int mca_get_neighbor_costs(size_t d,size_t N,const double * x,
                           struct Boundary * bound,
                           struct ValueF * vf, 
                           const size_t * ngrid, double ** xgrid,
                           size_t * fixed_ind , size_t * dim_vary,
                           int * absorbed, double * out)
{
    // BellmanParam must have

    // convert fibers to indices
    /* size_t * fixed_ind = calloc_size_t(d); */
    /* size_t dim_vary; */
    /* printf("convert fiber\n"); */
    int res = convert_fiber_to_ind(d,N,x,ngrid,xgrid,fixed_ind,dim_vary);
    if (res != 0){
        printf("\n======================================\n");
        printf("Error calling convert fiber to _ind!!\n\n");
        printf("grid is \n");
        for (size_t ii = 0; ii < d; ii++){
            dprint(ngrid[ii],xgrid[ii]);
        }
        printf("\n evaluate at %zu points\n",N);
        for (size_t ii = 0; ii < N; ii++){
            printf("first point = "); dprint(d,x+ii*d);
        }
        printf("\n======================================\n");
    }
    assert (res == 0);
    assert (*dim_vary < d);
    assert (N == ngrid[*dim_vary]);

    /* printf("converted fiber, dim_vary = %zu \n",*dim_vary); */
    /* printf("fixed_ind = "); iprint_sz(d,fixed_ind); */
    
    
    /* int * absorbed = calloc_int(N); */
    size_t * neighbors_vary = calloc_size_t(2*N);
    size_t * neighbors_fixed = calloc_size_t(2*(d-1));
    /* printf("process neighbors\n"); */
    res = process_fibers_neighbor(d,fixed_ind,*dim_vary,x,absorbed,
                                  neighbors_vary,neighbors_fixed,
                                  ngrid,bound);
    assert (res == 0);
    /* printf("done processing \n"); iprint(N,absorbed); */
    /* printf("processed neighbors\n"); */
    // evaluate cost associated with all neighbors
    /* printf("evaluated neighbors\n"); */
    res = valuef_eval_fiber_ind_nn(vf, fixed_ind, *dim_vary,
                                   neighbors_fixed, neighbors_vary,
                                   out);

    assert (res == 0);
    /* free(fixed_ind); fixed_ind = NULL; */
    /* free(absorbed); absorbed = NULL; */
    free(neighbors_vary); neighbors_vary = NULL;
    free(neighbors_fixed); neighbors_fixed = NULL;
    /* printf("done with calls\n"); */
    
    return 0;
}

/**********************************************************//**
    Process node for use within an implicit policy
**************************************************************/
int mca_get_neighbor_node_costs(size_t d, const double * x,
                                struct Boundary * bound,
                                struct ValueF * vf, 
                                const size_t * ngrid, double ** xgrid,
                                int * absorbed, double * out)
{

    int inobs = boundary_in_obstacle(bound,x);
    if (inobs == 1){
        *absorbed = -1;
        double val = valuef_eval(vf,x);
        for (size_t ii = 0; ii < 2*d+1; ii++){
            out[ii] = val;
        }
        return 0;
    }

    *absorbed = 0;
    double * xtemp = calloc_double(d);
    for (size_t ii = 0; ii < d; ii++){
        xtemp[ii] = x[ii];
    }

    for (size_t ii = 0; ii < d; ii++){
        double lb = xgrid[ii][0];
        double ub = xgrid[ii][ngrid[ii]-1];
        double h = xgrid[ii][1] - xgrid[ii][0];
        if ( ((x[ii] + h) < ub) && (x[ii]-h > lb)){ // standard case
            xtemp[ii] = x[ii]-h;
            out[2*ii] = valuef_eval(vf,xtemp);
            xtemp[ii] = x[ii]+h;
            out[2*ii+1] = valuef_eval(vf,xtemp);
            xtemp[ii] = x[ii];
        }
        else if ((x[ii]-h) <= lb){ // left boundary is hit
            xtemp[ii] = x[ii]+h;
            out[2*ii+1] = valuef_eval(vf,xtemp);

            enum EBTYPE b = boundary_type_dim(bound,ii,0);
            if ( (b == ABSORB) || (b == REFLECT) ){
                xtemp[ii] = lb;
                out[2*ii] = valuef_eval(vf,xtemp);
            }
            else if (b == PERIODIC){
                if (x[ii] > lb){
                    double diff = x[ii] - lb;
                    double left = h - diff; // assumes x[ii] > lb
                    xtemp[ii] = ub-left;
                    out[2*ii] = valuef_eval(vf,xtemp);
                }
                else{ // lb > x[ii]
                    double diff = lb-x[ii];
                    double te = ub - diff;
                    xtemp[ii] = te-h;
                    out[2*ii] = valuef_eval(vf,xtemp);
                }
            }
            else{
                fprintf(stderr,"No boundary specified!\n");
                assert(1 == 0);
            }            
            xtemp[ii] = x[ii];
        }
        else{ // right boundary is hit
            xtemp[ii] = x[ii]-h;
            out[2*ii] = valuef_eval(vf,xtemp);

            enum EBTYPE b = boundary_type_dim(bound,ii,1);
            if ( (b == ABSORB) || (b == REFLECT) ){
                xtemp[ii] = ub;
                out[2*ii+1] = valuef_eval(vf,xtemp);
            }
            else if (b == PERIODIC){
                /* printf ("hit it %G\n",x[ii]); */
                if (x[ii] < ub){
                    double diff = ub - x[ii];
                    double left = h - diff; // assumes x[ii] > lb
                    xtemp[ii] = lb+left;
                    /* printf("xtemp is %G \n",xtemp[ii]); */
                    out[2*ii+1] = valuef_eval(vf,xtemp);
                }
                else{ // lb > x[ii]
                    double diff = x[ii] - ub;
                    double te = lb + diff;
                    xtemp[ii] = te+h;
                    out[2*ii+1] = valuef_eval(vf,xtemp);
                }
            }
            else{
                fprintf(stderr,"No boundary specified!\n");
                assert(1 == 0);
            }            
            xtemp[ii] = x[ii];
        }
    }

    free(xtemp); xtemp = NULL;
    return 0;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// I don't think below is used
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/**********************************************************//**
    Process fibers to figure out boundary information

    \param[in]     d         - dimension of state space
    \param[in]     nvals     - number of evaluations
    \param[in]     x         - (Nd,) evaluation locations
    \param[in,out] boundv    - (N,) boundary condition specifies 
    \param[in,out] neighbors - (2Nd,) neighbor specifiers
    \param[in]     bound     - structure holding bc information

    \return 0 if everything is ok, 
            
    \note
    boundv[ii] = -1 if absorbed in an obstacle
    boundv[ii] = 0  if not absorbed
    boundv[ii] = 1  if absored in a boundary
**************************************************************/
int process_fibers(size_t d, 
                   size_t nvals,
                   const double * x, 
                   int * boundv, 
                   int * neighbors, 
                   struct Boundary * bound)
{
    double time = 0.0;
    for (size_t ii = 0; ii < nvals; ii++){
        struct BoundInfo * bi = boundary_type(bound,time,x+ii*d);
        boundv[ii] = bound_info_onbound(bi);
        if (boundv[ii] == 0){ // use all neighbors
            for (size_t jj = 0; jj < d; jj++){
                neighbors[ii*(2*d) + 2*jj] = 1;
                neighbors[ii*(2*d) + 2*jj+1] = 1;
            }
        }
        else{
            if (bound_info_absorb(bi) == 1){ // absorbed
                /* printf("x = "); dprint(d,x+ii*d); */
                int in_obstacle = bound_info_get_in_obstacle(bi);
                if (in_obstacle >= 0){
                    /* printf("in obstacle\n"); */
                    boundv[ii] = -1;
                }
                for (size_t jj = 0; jj < d; jj++){ // use none of the neighbors
                    neighbors[ii*(2*d) + 2*jj] = 0;
                    neighbors[ii*(2*d) + 2*jj+1] = 0;
                }   
            }
            else{
                if (bound_info_reflect(bi) == 0){ // reflected back
                    for (size_t jj = 0; jj < d; jj++){
                        int on_bound = bound_info_onbound_dim(bi,ii);
                        if (on_bound == 0){
                            neighbors[ii*(2*d) + 2*jj] = 1;
                            neighbors[ii*(2*d) + 2*jj+1] = 1;
                        }
                        else{
                            int reflect_dir = bound_info_reflect_dim_dir(bi,ii);
                            //printf("before\n");
                            if (reflect_dir == -1){ // just do right side poition
                                neighbors[ii*(2*d) + 2*jj] = 0;
                                neighbors[ii*(2*d) + 2*jj+1] = 1;
                            }
                            else{
                                neighbors[ii*(2*d) + 2*jj] = 1;
                                neighbors[ii*(2*d) + 2*jj+1] = 0;
                            }
                        }
                    }
                }
                else{
                    fprintf(stderr, "Have not dealt with periodic boundaries yet\n");
                    exit(1);
                }
            }
        }
        bound_info_free(bi); bi = NULL;
    }
    return 0;
}


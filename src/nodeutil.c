#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "c3.h"

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
    \param[in,out] grad_prob  - gradient of probabilities
    \param[in,out] dt         - interpolation interval
    \param[in,out] grad_dt    - gradient of interpolation interval
    \param[in,out] space      - workspace needed if gradients requested (du,) array

    \return 0 if everything is ok, 1 if computed dt is infinity

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
    for (size_t ii = 0; ii < dx; ii++){
        t = h2/hvec[ii];
        diff = ddiff[ii*dx+ii] * ddiff[ii*dx+ii];
        t2 = h2 / hvec[ii] / hvec[ii];

        prob[2*ii] = t2 * diff / 2.0;
        prob[2*ii+1] = t2 * diff / 2.0;
        if (drift[ii] < 0){ // add transition to the left
            prob[2*ii] -= t * drift[ii];
        }
        else if (drift[ii] > 0){ // add transition to the right
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
            if (drift[ii] < 0){
                cblas_daxpy(du,-t,grad_drift+ii,dx,grad_prob + (2*ii)*du,1);
            }
            else if (drift[ii] > 0){
                cblas_daxpy(du,t,grad_drift+ii,dx,grad_prob + (2*ii+1)*du,1);
            }
            else{
                /* printf("warning! not sure what to do here!\n"); */
                for (size_t jj = 0; jj < du; jj++){
                    if (grad_drift[jj*dx+ii] < 0){ // wants to go left with changing u
                        cblas_daxpy(du,-t,grad_drift+ii,dx,grad_prob + (2*ii)*du,1);
                    }
                    else if (grad_drift[jj*dx+ii] > 0){ // wants to go right with changing u
                        cblas_daxpy(du,t,grad_drift+ii,dx,grad_prob + (2*ii+1)*du,1);
                    }
                    else{ // do nothing because doesn't want to change!
                        /* printf("warning! not sure what to do here!\n"); */
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

    return 0;
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
int convert_fiber_to_ind(size_t d, size_t N, const double ** x, const size_t * Ngrid, double ** xgrid,
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

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
    \param[in,out] absorbed        - (0,no), (1,yes), (2,obstacle)
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

    for (size_t ii = 0; ii < ngrid[dim_vary]; ii++){
        int inobs = boundary_in_obstacle(bound,x+ii*d);
        if (inobs == 0){
            absorbed[ii] = 0;
        }
        else{
            absorbed[ii] = -1;
        }
    }

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
    }
    else if (b == PERIODIC){
        neighbors_vary[0] = ngrid[dim_vary]-2;
        neighbors_vary[1] = 1;
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
    }
    else if (b == PERIODIC){
        neighbors_vary[2*(ngrid[dim_vary]-1)] = ngrid[dim_vary]-2;
        neighbors_vary[2*(ngrid[dim_vary]-1)+1] = 1;
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

int mca_get_neighbor_costs(size_t d,size_t N,const double * x,struct Boundary * bound,
                           struct ValueF * vf, const size_t * ngrid, double ** xgrid,
                           int * absorbed,
                           double * out)
{
    // BellmanParam must have

    // convert fibers to indices
    size_t * fixed_ind = calloc_size_t(d);
    size_t dim_vary;
    int res = convert_fiber_to_ind(d,N,x,ngrid,xgrid,fixed_ind,&dim_vary);
    assert (res == 0);

    /* int * absorbed = calloc_int(N); */
    size_t * neighbors_vary = calloc_size_t(2*N);
    size_t * neighbors_fixed = calloc_size_t(2*(d-1));
    res = process_fibers_neighbor(d,fixed_ind,dim_vary,x,absorbed,
                                  neighbors_vary,neighbors_fixed,
                                  ngrid,bound);
    assert (res == 0);

    // evaluate cost associated with all neighbors
    res = valuef_eval_fiber_ind_nn(vf, fixed_ind, dim_vary,
                                   neighbors_fixed, neighbors_vary,
                                   out);

    free(fixed_ind); fixed_ind = NULL;
    /* free(absorbed); absorbed = NULL; */
    free(neighbors_vary); neighbors_vary = NULL;
    free(neighbors_fixed); neighbors_fixed = NULL;
    assert (res == 0);
    
    return 0;
}


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
int process_fibers(size_t d, size_t nvals, const double * x, int * boundv, int * neighbors, 
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

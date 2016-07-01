// Copyright (c) 2016, Massachusetts Institute of Technology
// Authors: Alex Gorodetsky
// Email  : goroda@mit.edu

//Code

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CuTest.h"	

#include "c3.h"
#include "nodeutil.h"
#include "valuefunc.h"
#include "bellman.h"
#include "c3sc.h"

//////////////////////////////////////////////////////
// Dynamics 
//////////////////////////////////////////////////////

int f1(double t, const double * x, const double * u, double * out,
       double * jac, void * args)
{
    (void)(t);
    (void)(args);
    out[0] = sin(x[1]+ x[0]);
    out[1] = pow(x[0],2) * u[0];
    if (jac != NULL){
        //df1/du
        jac[0] = 0.0;
        jac[1] = pow(x[0],2);
    }
    return 0;
}

int s1(double t,const double * x,const double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(args);
    /* printf("x in diffusion = "); dprint(2,x); */
    out[0] = 1e-2;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 1e-2;
    if (grad != NULL){
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
    }
    return 0;
}

// 3 u
// 3 x
int f2(double t, const double * x, const double * u, double * out,
       double * jac, void * args)
{
    (void)(t);
    (void)(args);
    out[0] = x[0]*pow(x[2],2)*u[0] * cos(u[1]);
    out[1] = -x[1]  * u[2];
    out[2] = x[0]*x[1]*u[0] + 2 * u[1];
    if (jac != NULL){
        
        jac[0] = x[0]*pow(x[2],2) * cos(u[1]); 
        jac[1] = 0.0;                           
        jac[2] =  x[0]*x[1];                    
        
        jac[3] = x[0]*pow(x[2],2)* u[0] * (- sin(u[1]));
        jac[4] = 0.0;
        jac[5] = 2.0;

        jac[6] = 0.0;
        jac[7] = -x[1];
        jac[8] = 0.0;
    }
    return 0;
}

int s2(double t,const double * x,const double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(args);
    for (size_t ii = 0; ii < 9; ii++){
        out[ii] = 0.0;
    }
    for (size_t jj = 0; jj < 3; jj++){
        out[jj*3+jj] = 1e-2;
    }

    if (grad != NULL){
        for (size_t ii = 0; ii < 27; ii++){
            grad[ii] = 0.0;
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////
// Cost functions
////////////////////////////////////////////////////////////////
int stagecost2d(double t,const double * x,const double * u, double * out, 
              double * grad)
{
    (void)(t);

    *out = 0.0;

    // states
    *out += 20.0 * x[0]*x[0];
    *out += 50.0 * x[1]*x[1];

    // control
    *out += 0.1 * u[0] * u[0];
    if (grad!= NULL){
        grad[0] = 0.1 * 2.0 * u[0];
    }

    return 0;
}

int stagecost3d(double t,const double * x,const double * u, double * out, 
              double * grad)
{
    (void)(t);

    *out = 0.0;

    // states
    *out += 20.0 * x[0]*x[0];
    *out += 50.0 * x[1]*x[1];
    *out += 70.0 * x[2]*x[2];

    // control
    *out += 0.1 * u[0] * u[0];
    *out += 0.5 * u[1] * u[1];
    *out += 3.0 * u[2] * u[2];
    if (grad!= NULL){
        grad[0] = 0.1 * 2.0 * u[0];
        grad[1] = 0.5 * 2.0 * u[1];
        grad[2] = 3.0 * 2.0 * u[2];
    }

    return 0;
}

int boundcost(double t,const  double * x, double * out)
{

    (void)(t);
    (void)(x);
    *out = 0.0;
    *out = 100.0;
    return 0;
}

int ocost(const double * x,double * out)
{
    (void)(x);
    /* printf("got ocost!!\n"); */
    /* dprint(2,x); */
    *out = 0.0;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

void Test_tprob_probsum(CuTest * tc)
{
    printf("Testing Function: transition_assemble (1)\n");

    size_t dx = 2;
    size_t du = 1;
    size_t dw = 2;
    double h[2] = {1e-1, 1e-2};
    double hmin = h[1];

    double pt[2] = {-2.0,-0.3};
    double u[1] = {0.0};
    double t = 0.2;
    
    double drift[2];
    double grad_drift[2];
    double diff[4];
    double grad_diff[4];

    int res = f1(t,pt,u,drift,grad_drift,NULL);
    CuAssertIntEquals(tc,0,res);

    res = s1(t,pt,u,diff,grad_diff,NULL);
    CuAssertIntEquals(tc,0,res);

    double prob[5];
    double dt;
    res = transition_assemble(dx,du,dw,hmin,h,
                              drift,NULL,diff,NULL,prob,NULL,&dt,NULL,NULL);
    CuAssertIntEquals(tc,0,res);

    double psum = 0.0;
    for (size_t ii = 0; ii < 2*dx+1; ii++){
        psum += prob[ii];
        CuAssertIntEquals(tc,1,prob[ii]>-1e-15);
    }
    CuAssertDblEquals(tc,1,psum,1e-15);
}

void Test_tprob_grad(CuTest * tc)
{
    printf("Testing Function: transition_assemble (2)\n");

    int res;
    size_t dx = 2;
    size_t du = 1;
    size_t dw = 2;
    double h[2] = {1e-1, 1e-2};
    double hmin = h[1];

    double drift[2];
    double grad_drift[2];
    double diff[4];
    double grad_diff[4];

    double space[1];
    double prob[5];
    double prob2[5];
    double grad_prob[5*1];
    double dt,dt2;
    double dt_grad[1];

    size_t nrand = 1000;
    for (size_t kk = 0; kk < nrand; kk++){
        /* printf("kk = %zu\n",kk); */
        /* double pt[2] = {1.0,2.0}; */
        double pt[2] = {randu()*3.0-1.5,randu()*3.0-1.5};
        
        double t = 0.2;

        double u[1] = {randu()*3.0-1.5};
        /* double u[1] = {0.0}; */
    
        res =f1(t,pt,u,drift,grad_drift,NULL);
        CuAssertIntEquals(tc,0,res);
        res = s1(t,pt,u,diff,grad_diff,NULL);
        CuAssertIntEquals(tc,0,res);
        res = transition_assemble(dx,du,dw,hmin,h,drift,grad_drift,
                                  diff,grad_diff,prob,grad_prob,&dt,dt_grad,space);
        CuAssertIntEquals(tc,0,res);

        /* printf("drift = (%G,%G)\n",drift[0],drift[1]); */
        /* printf("gdrift = (%G,%G)\n",grad_drift[0],grad_drift[1]); */
        double delta = 1e-10;
        u[0] += delta;
        res =f1(t,pt,u,drift,grad_drift,NULL);
        CuAssertIntEquals(tc,0,res);
        res = s1(t,pt,u,diff,grad_diff,NULL);
        CuAssertIntEquals(tc,0,res);
        u[0] -= delta;
        res = transition_assemble(dx,du,dw,hmin,h,drift,NULL,diff,NULL,prob2,NULL,
                                  &dt2,NULL,NULL);
        CuAssertIntEquals(tc,0,res);

        /* printf("prob transitions = "); dprint(5,prob); */
        for (size_t ii = 0; ii < 2*dx+1; ii++){
            /* printf("ii=(%zu/%zu)\n",ii,2*dx); */
            double fdgrad = (prob2[ii] - prob[ii])/delta;
            CuAssertDblEquals(tc,fdgrad,grad_prob[ii],1e-5);
            /* printf("(%G,%G) \n",grad_prob[ii],fdgrad); */
        }
        double fddtgrad = (dt2-dt)/delta;
        CuAssertDblEquals(tc,fddtgrad,dt_grad[0],1e-5);
            
    }
}

void Test_tprob_grad2(CuTest * tc)
{
    printf("Testing Function: transition_assemble (3)\n");

    int res;
    size_t dx = 3;
    size_t du = 3;
    size_t dw = 3;
    double h[3] = {1e-1, 1e-2, 1e0};
    double hmin = h[1];

    double drift[3];
    double grad_drift[9];
    double diff[9];
    double grad_diff[27];

    double space[3];
    double prob[7];
    double prob2[7];
    double grad_prob[7*3];
    double dt,dt2;
    double dt_grad[3];

    size_t nrand = 1000;
    for (size_t kk = 0; kk < nrand; kk++){
        /* printf("kk = %zu\n",kk); */
        /* double pt[2] = {1.0,2.0}; */
        double pt[3] = {randu()*3.0-1.5,
                        randu()*3.0-1.5,
                        randu()*3.0-1.5};
        
        double t = 0.2;

        double u[3] = {randu()*3.0-1.5, 
                       randu()*3.0-1.5,
                       randu()*3.0-1.5
                      };

        /* double u[1] = {0.0}; */
    
        res = f2(t,pt,u,drift,grad_drift,NULL);
        CuAssertIntEquals(tc,0,res);
        res = s2(t,pt,u,diff,grad_diff,NULL);
        CuAssertIntEquals(tc,0,res);
        res = transition_assemble(dx,du,dw,hmin,h,drift,grad_drift,
                                  diff,grad_diff,prob,grad_prob,&dt,dt_grad,
                                  space);
        CuAssertIntEquals(tc,0,res);

        for (size_t jj = 0; jj < 2*dx+1; jj++){
            /* printf("jj = (%zu/%zu)\n",jj,2*dx); */
            for (size_t ii = 0; ii < du; ii++){
                /* printf("\t ii = (%zu/%zu)\n",ii,du-1); */
                double delta = 1e-8;
                u[ii] += delta;                
                res = f2(t,pt,u,drift,grad_drift,NULL);
                CuAssertIntEquals(tc,0,res);
                res = s2(t,pt,u,diff,grad_diff,NULL);
                CuAssertIntEquals(tc,0,res);
                u[ii] -= delta;

                res = transition_assemble(dx,du,dw,hmin,h,
                                          drift,NULL,diff,NULL,prob2,NULL,
                                          &dt2,NULL,NULL);
                CuAssertIntEquals(tc,0,res);
                /* printf("ii=(%zu/%zu)\n",ii,2*dx); */
                double fdgrad = (prob2[jj] - prob[jj])/delta;
                CuAssertDblEquals(tc,fdgrad,grad_prob[jj*du+ii],1e-4);
                /* printf("\t\t (%G,%G) \n",grad_prob[jj*du+ii],fdgrad); */
                double fddtgrad = (dt2-dt)/delta;
                CuAssertDblEquals(tc,fddtgrad,dt_grad[ii],1e-5);

            }
        }
    }
}


CuSuite * TProbGetSuite()
{
    // tests assembly of transition probabilities
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_tprob_probsum);
    SUITE_ADD_TEST(suite, Test_tprob_grad);
    SUITE_ADD_TEST(suite, Test_tprob_grad2);

    return suite;
}


double quad(const double * x, void * arg)
{
    (void) (arg);
    double out = x[0]*x[0] + x[0]*x[1] + x[2]*x[2];
    return out;
}

void Test_valuef_neighbor_eval(CuTest * tc)
{
    printf("Testing Functions associated with evaluating the cost of a fibers neighbors \n");
    printf("                  valuef_eval_fiber_ind_nn\n");

    size_t dim = 3;

    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-8);
    approx_args_set_round_tol(aargs,1e-7);
    approx_args_set_kickrank(aargs,10);
    approx_args_set_adapt(aargs,0);
    size_t start_rank = 20;
    approx_args_set_startrank(aargs,start_rank);
    approx_args_set_maxrank(aargs,start_rank);

    //initialize
    const size_t N[3] = {30, 43, 24};
    double * xgrid[3];
    double * start[3];
    for (size_t ii = 0; ii < 3; ii++){
        xgrid[ii] = linspace(-1.0,2.0,N[ii]);
        start[ii] = calloc_double(start_rank);
        size_t stride = uniform_stride(N[ii],start_rank);
        for (size_t jj = 0; jj < start_rank; jj++){
            start[ii][jj] = xgrid[ii][stride*jj];
        }
        /* printf("start[%zu] = ",ii); dprint(start_rank,start[ii]); */
    }

    //interpolate
    struct ValueF * vf = valuef_interp(3,quad,NULL,N,xgrid,start,aargs,0);

    // check neighbor cost computation function
    size_t fixed_ind[3] = { 3, 5, 9};
    size_t neighbors[4] = { 1, 4, 4, 6};
    double pt[3];

    size_t dim_vary = 1;
    double * fiber_vals = calloc_double(N[dim_vary] * (2*dim+1));
    size_t * neighbor_vary = calloc_size_t(N[dim_vary]*2);
    neighbor_vary[0] = 0;
    neighbor_vary[1] = 1;
    for (size_t ii = 1; ii < N[dim_vary]; ii++){
        neighbor_vary[2*ii] = ii-1;
        neighbor_vary[2*ii+1] = ii+1;
    }
    neighbor_vary[2*N[dim_vary]-1] = N[dim_vary]-1;

    valuef_eval_fiber_ind_nn(vf, fixed_ind, dim_vary,
                             neighbors, neighbor_vary,
                             fiber_vals);

    /* printf("size think = %zu\n",nevals); */
    /* dprint(nevals,vals); */
    for (size_t jj = 0; jj < N[dim_vary]; jj++){

        for (size_t kk = 0; kk < dim; kk++){
            pt[kk] = xgrid[kk][fixed_ind[kk]];
        }
        pt[dim_vary] = xgrid[dim_vary][jj];

        double val_should = valuef_eval(vf,pt);
        double val_is = fiber_vals[jj*(2*dim+1) + 2*dim];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
        if (jj == 0){
            pt[1] = xgrid[1][jj];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[1] = xgrid[1][jj+1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }
        else if (jj < (N[dim_vary]-1)){
            pt[1] = xgrid[1][jj-1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[1] = xgrid[1][jj+1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }
        else{
            pt[1] = xgrid[1][jj-1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[1] = xgrid[1][jj];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }

        size_t ii = 0;
        // first point
        pt[0] = xgrid[0][neighbors[0]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][neighbors[1]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        ii = 2;
        // second point
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][neighbors[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][neighbors[3]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
    }
    free(neighbor_vary); neighbor_vary = NULL;
    free(fiber_vals); fiber_vals = NULL;

    /* printf("second\n"); */
    dim_vary = 0;
    fiber_vals = calloc_double((2*dim+1)*N[dim_vary]);
    neighbor_vary = calloc_size_t(N[dim_vary]*2);
    neighbor_vary[0] = 0;
    neighbor_vary[1] = 1;
    for (size_t ii = 1; ii < N[dim_vary]; ii++){
        neighbor_vary[2*ii] = ii-1;
        neighbor_vary[2*ii+1] = ii+1;
    }
    neighbor_vary[2*N[dim_vary]-1] = N[dim_vary]-1;

    valuef_eval_fiber_ind_nn(vf, fixed_ind, dim_vary,
                             neighbors, neighbor_vary,
                             fiber_vals);
    
    /* printf("got it\n"); */
    
    /* printf("size think = %zu\n",nevals); */
    /* dprint(nevals,vals); */
    for (size_t jj = 0; jj < N[dim_vary]; jj++){
        // no neighbors
        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][fixed_ind[2]];
        
        double val_should = valuef_eval(vf,pt);
        double val_is = fiber_vals[jj*(2*dim+1) + 2*dim];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        if (jj == 0){
            pt[0] = xgrid[0][jj];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[0] = xgrid[0][jj+1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }
        else if (jj < (N[dim_vary]-1)){
            pt[0] = xgrid[0][jj-1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[0] = xgrid[0][jj+1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }
        else{
            pt[0] = xgrid[0][jj-1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[0] = xgrid[0][jj];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }


        size_t ii = 1;
        // first point
        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][neighbors[0]];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][neighbors[1]];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii+1];

        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        ii = 2;
        // second point
        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][neighbors[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][neighbors[3]];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
    }
    free(neighbor_vary); neighbor_vary = NULL;
    free(fiber_vals); fiber_vals = NULL;

    dim_vary = 2;
    fiber_vals = calloc_double((2*dim+1)*N[dim_vary]);
    neighbor_vary = calloc_size_t(N[dim_vary]*2);
    neighbor_vary[0] = 0;
    neighbor_vary[1] = 1;
    for (size_t ii = 1; ii < N[dim_vary]; ii++){
        neighbor_vary[2*ii] = ii-1;
        neighbor_vary[2*ii+1] = ii+1;
    }
    neighbor_vary[2*N[dim_vary]-1] = N[dim_vary]-1;

    valuef_eval_fiber_ind_nn(vf, fixed_ind, dim_vary,
                             neighbors, neighbor_vary,
                             fiber_vals);

    for (size_t jj = 0; jj < N[dim_vary]; jj++){

        // no neighbors
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][jj];
        
        double val_should = valuef_eval(vf,pt);
        double val_is = fiber_vals[jj*(2*dim+1) + 2*dim];
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
        
        if (jj == 0){
            pt[2] = xgrid[2][jj];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[2] = xgrid[2][jj+1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }
        else if (jj < (N[dim_vary]-1)){
            pt[2] = xgrid[2][jj-1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[2] = xgrid[2][jj+1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }
        else{
            pt[2] = xgrid[2][jj-1];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
            
            pt[2] = xgrid[2][jj];
            val_should = valuef_eval(vf,pt);
            val_is = fiber_vals[jj*(2*dim+1) + 2*(dim_vary)+1];        
            CuAssertDblEquals(tc,val_should,val_is,1e-14);
        }

        // first point
        size_t ii = 0;

        pt[0] = xgrid[0][neighbors[0]];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][jj];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][neighbors[1]];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][jj];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        ii = 1;
        // second point
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][neighbors[2]];
        pt[2] = xgrid[2][jj];
        
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        // second point
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][neighbors[3]];
        pt[2] = xgrid[2][jj];
            
        val_should = valuef_eval(vf,pt);
        val_is = fiber_vals[jj*(2*dim+1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
    }
    free(fiber_vals); fiber_vals = NULL;
    free(neighbor_vary); neighbor_vary = NULL;


    /* printf("done sdaskd\n"); */

    valuef_destroy(vf); vf = NULL;
    for (size_t ii = 0; ii < dim; ii++){
        free(xgrid[ii]); xgrid[ii] = NULL;
        free(start[ii]); start[ii] = NULL;
    }
    approx_args_free(aargs); aargs = NULL;

}

void Test_valuef_fiber_to_ind(CuTest * tc)
{
    printf("Testing Function: fiber_to_ind\n");

    size_t dim = 3;

    //initialize
    const size_t N[3] = {30, 43, 24};
    double * xgrid[3];
    for (size_t ii = 0; ii < 3; ii++){
        xgrid[ii] = linspace(-1.0,2.0,N[ii]);
    }

    size_t fixed_ind_true[3] = { 10, 12, 13};
    for (size_t ll = 0; ll < dim; ll++){
        /* printf("ll = %zu\n",ll); */
        size_t dim_vary_true = ll;
        double * x = calloc_double(N[dim_vary_true] * dim);
        for (size_t ii = 0; ii < N[dim_vary_true]; ii++){
            x[ii*dim] = xgrid[0][fixed_ind_true[0]];
            x[ii*dim+1] = xgrid[1][fixed_ind_true[1]];
            x[ii*dim+2] = xgrid[2][fixed_ind_true[2]];
            x[ii*dim+dim_vary_true] = xgrid[dim_vary_true][ii];
        }
    
        size_t fixed_ind[3];
        size_t dim_vary;
        int res = convert_fiber_to_ind(dim,N[dim_vary_true],x,N,xgrid,fixed_ind,&dim_vary);
        CuAssertIntEquals(tc,0,res);
        CuAssertIntEquals(tc,dim_vary_true,dim_vary);
        for (size_t ii = 0; ii < dim; ii++){
            if (ii != dim_vary_true){
                CuAssertIntEquals(tc,fixed_ind_true[ii],fixed_ind[ii]);
            }
        }
        free(x); x = NULL;
    }
    for (size_t ii = 0; ii < dim; ii++){
        free(xgrid[ii]); xgrid[ii] = 0;
    }
}

CuSuite * ValueFGetSuite()
{
    //printf("----------------------------\n");
    // tests value function

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_valuef_neighbor_eval);
    SUITE_ADD_TEST(suite, Test_valuef_fiber_to_ind);

    return suite;
}

void Test_bound_nodes(CuTest * tc)
{
    printf("Testing function: process_fibers\n");
    size_t d = 3;
    double lb[3] = {-1.0, -2.0, -3.0};
    double ub[3] = { 1.0, 2.0, 3.0 };

    double center[3] = {0.0,0.0,0.0};
    double lengths[3] = {0.2,0.2,0.8};
    
    struct Boundary * bound = boundary_alloc(d,lb,ub);
    boundary_add_obstacle(bound,center,lengths);

    size_t N[3] = {20, 40, 60};
    double * grid[3];
    size_t nvals = 1;
    for (size_t ii = 0; ii < d; ii++){
        grid[ii] = linspace(lb[ii],ub[ii],N[ii]);
        nvals *= N[ii];
    }

    double * x = calloc_double(d * nvals);
    
    size_t onval = 0;
    int * bound_true = calloc_int(nvals);
    for (size_t ii = 0; ii < N[0]; ii++){
        for (size_t jj = 0; jj < N[1]; jj++){
            for (size_t kk = 0; kk < N[2]; kk++){
                x[onval*d+0] = grid[0][ii];
                x[onval*d+1] = grid[1][jj];
                x[onval*d+2] = grid[2][kk];
                if ( (fabs(grid[0][ii]) < (lengths[0]/2.0)) &&
                     (fabs(grid[1][jj]) < (lengths[1]/2.0)) &&
                     (fabs(grid[2][kk]) < (lengths[2]/2.0)))
                {
                    /* printf("here\n"); */
                    bound_true[onval] = -1;
                }
                else if ( (ii == 0) || (ii == (N[0]-1)) ||
                          (jj == 0) || (jj == (N[1]-1)) ||
                          (kk == 0) || (kk == (N[2]-1)))
                {
                    /* printf("there\n"); */
                    bound_true[onval] = 1;
                }
                onval++;
            }
        }
    }
    

    int * neighbors = calloc_int(2*d*nvals);
    int * boundaries = calloc_int(nvals);
    int res = process_fibers(d,nvals,x,boundaries,neighbors,bound);
    CuAssertIntEquals(tc,0,res);
    for (size_t ii = 0; ii < nvals; ii++){
        /* printf("ii = %zu ",ii); dprint(3,x + ii*d); */
        CuAssertIntEquals(tc,bound_true[ii],boundaries[ii]);
        if (bound_true[ii] == 1){
            CuAssertIntEquals(tc,0,neighbors[ii*2*d]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+1]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+2]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+3]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+4]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+5]);
        }
        else if (bound_true[ii] == -1){
            CuAssertIntEquals(tc,0,neighbors[ii*2*d]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+1]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+2]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+3]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+4]);
            CuAssertIntEquals(tc,0,neighbors[ii*2*d+5]);
        }
        else{
            CuAssertIntEquals(tc,1,neighbors[ii*2*d]);
            CuAssertIntEquals(tc,1,neighbors[ii*2*d+1]);
            CuAssertIntEquals(tc,1,neighbors[ii*2*d+2]);
            CuAssertIntEquals(tc,1,neighbors[ii*2*d+3]);
            CuAssertIntEquals(tc,1,neighbors[ii*2*d+4]);
            CuAssertIntEquals(tc,1,neighbors[ii*2*d+5]);
        }
    }

    free(x); x = NULL;
    free(neighbors); neighbors = NULL;
    free(boundaries); boundaries = NULL;
    free(bound_true); bound_true = NULL;
    for (size_t ii = 0; ii < d; ii++){
        free(grid[ii]); grid[ii] = NULL;
    }
    boundary_free(bound); bound = NULL;
}

void Test_process_fibers_neighbor(CuTest * tc)
{
    printf("Testing Function: process_fibers_neighbor\n");

    size_t dim = 3;
    double lb[3] = {-1.0, -2.0, -3.0};
    double ub[3] = { 1.0, 2.0, 3.0 };

    double center[3] = {0.0,0.0,0.0};
    double lengths[3] = {0.8,0.8,0.8};
    
    struct Boundary * bound = boundary_alloc(dim,lb,ub);
    boundary_add_obstacle(bound,center,lengths);

    //initialize
    const size_t N[3] = {30, 43, 24};
    double * xgrid[3];
    for (size_t ii = 0; ii < 3; ii++){
        xgrid[ii] = linspace(lb[ii],ub[ii],N[ii]);
    }

    size_t fixed_ind_true[3];// = { 15, 21, 12};
    for (size_t aa = 0; aa < N[0]; aa++){
        fixed_ind_true[0] = aa;
        for (size_t bb = 0; bb < N[1]; bb++){
            fixed_ind_true[1] = bb;
            for (size_t cc = 0; cc < N[2]; cc++){
                fixed_ind_true[2] = cc;
                for (size_t ll = 1; ll < dim; ll++){
                    /* size_t ll = 0; */
                    /* printf("ll = %zu\n",ll); */
                    size_t dim_vary_true = ll;


                    int * absorbed = calloc_int(N[dim_vary_true]);
                    size_t * neighbors_fixed = calloc_size_t(2*(dim-1));
                    size_t * neighbors_vary = calloc_size_t(2*N[dim_vary_true]);
                    int * absorbed_true = calloc_int(N[dim_vary_true]);

                    double * x = calloc_double(N[dim_vary_true] * dim);
                    for (size_t ii = 0; ii < N[dim_vary_true]; ii++){
                        x[ii*dim] = xgrid[0][fixed_ind_true[0]];
                        x[ii*dim+1] = xgrid[1][fixed_ind_true[1]];
                        x[ii*dim+2] = xgrid[2][fixed_ind_true[2]];
                        x[ii*dim+dim_vary_true] = xgrid[dim_vary_true][ii];

                        if ( (fabs(x[ii*dim]) < (lengths[0]/2.0)) &&
                             (fabs(x[ii*dim+1]) < (lengths[1]/2.0)) &&
                             (fabs(x[ii*dim+2]) < (lengths[2]/2.0)))
                        {
                            /* printf("here\n"); */
                            absorbed_true[ii] = -1;
                        }
                        else if ( (x[ii*dim] <= lb[0]+1e-12) || 
                                  (x[ii*dim] >= ub[0]-1e-12) ||
                                  (x[ii*dim+1] <= lb[1]+1e-12) || 
                                  (x[ii*dim+1] >= ub[1]-1e-12) ||
                                  (x[ii*dim+2] <= lb[2]+1e-12) || 
                                  (x[ii*dim+2] >= ub[2]-1e-12))
                        {
                            /* printf("there\n"); */
                            /* dprint(3,x+ii*dim); */
                            absorbed_true[ii] = 1;
                        }
                    }
    
                    size_t fixed_ind[3];
                    size_t dim_vary;
                    /* printf("there!\n"); */
                    int res = convert_fiber_to_ind(dim,N[dim_vary_true],x,N,xgrid,
                                                   fixed_ind,&dim_vary);
                    CuAssertIntEquals(tc,0,res);
                    for (size_t ii = 0; ii < dim; ii++){
                        if (ii != dim_vary_true){
                            CuAssertIntEquals(tc,fixed_ind_true[ii],fixed_ind[ii]);
                        }
                    }        

                    res = process_fibers_neighbor(dim,fixed_ind,dim_vary,x,absorbed,
                                                  neighbors_vary,neighbors_fixed,N,bound);
        
                    for (size_t ii = 0; ii < N[dim_vary_true]; ii++){
                        /* printf("x = "); dprint(3,x+ii*dim); */
                        CuAssertIntEquals(tc,absorbed_true[ii],absorbed[ii]);
                    }

                    free(absorbed_true); absorbed_true = NULL;
                    free(absorbed); absorbed = NULL;
                    free(neighbors_fixed); neighbors_fixed = NULL;
                    free(neighbors_vary); neighbors_vary = NULL;
                    free(x); x = NULL;
                }
            }
        }
    }


    boundary_free(bound); bound = NULL;
    for (size_t ii = 0; ii < dim; ii++){
        free(xgrid[ii]); xgrid[ii] = 0;
    }
}

void Test_bellman_grad1(CuTest * tc)
{
    printf("Testing Function: bellmanrhs gradient\n");

    size_t dx = 2;
    size_t du = 1;
    size_t dw = 2;
    double h[2] = {1e-1, 1e-2};
    double hmin = h[1];

    double drift[2];
    double grad_drift[2];
    double diff[4];
    double grad_diff[4];

    double space[1];
    double prob[5];
    double prob2[5];
    double grad_prob[5*1];
    double dt,dt2;
    double dt_grad[1];

    double grad_stage[1];
    double grad_bellman[1];

    double discount = 0.1;

    size_t nrand = 100;
    for (size_t kk = 0; kk < nrand; kk++){
        /* printf("kk = %zu\n",kk); */
        /* double pt[2] = {1.0,2.0}; */
        double t = 0.2;
        double pt[2] = {randu()*3.0-1.5,randu()*3.0-1.5};
        double u[1] = {randu()*3.0-1.5};
        double cost[5] = {randu(),randu(),randu(),randu(),randu()};

        double stage_cost;
        int res = stagecost2d(t,pt,u,&stage_cost,grad_stage);
        CuAssertIntEquals(tc,0,res);
    
        res =f1(t,pt,u,drift,grad_drift,NULL);
        CuAssertIntEquals(tc,0,res);
        res = s1(t,pt,u,diff,grad_diff,NULL);
        CuAssertIntEquals(tc,0,res);
        res = transition_assemble(dx,du,dw,hmin,h,drift,grad_drift,
                                  diff,grad_diff,prob,grad_prob,&dt,dt_grad,space);
        CuAssertIntEquals(tc,0,res);

        double val = bellmanrhs(dx,du,stage_cost,grad_stage,
                                discount, prob, grad_prob, dt, dt_grad,
                                cost, grad_bellman);


        /* printf("drift = (%G,%G)\n",drift[0],drift[1]); */
        /* printf("gdrift = (%G,%G)\n",grad_drift[0],grad_drift[1]); */
        double delta = 1e-9;
        u[0] += delta;
        res =f1(t,pt,u,drift,grad_drift,NULL);
        CuAssertIntEquals(tc,0,res);
        res = s1(t,pt,u,diff,grad_diff,NULL);
        CuAssertIntEquals(tc,0,res);
        res = transition_assemble(dx,du,dw,hmin,h,drift,NULL,diff,NULL,prob2,NULL,
                                  &dt2,NULL,NULL);

        double stage_cost2;
        res = stagecost2d(t,pt,u,&stage_cost2,NULL);
        CuAssertIntEquals(tc,0,res);
        double newval = bellmanrhs(dx,du,stage_cost2,NULL,
                                   discount, prob2,NULL, dt2, NULL,
                                   cost, NULL);

        /* printf("val = %G, newval=%G\n",val,newval); */
        double fdgrad = (newval - val)/delta;
        CuAssertDblEquals(tc,fdgrad,grad_bellman[0],1e-5);

        u[0] -= delta;

    }
}

void Test_bellman_grad3d(CuTest * tc)
{
    printf("Testing Function: bellmanrhs gradient with 3d control\n");

    size_t dx = 3;
    size_t du = 3;
    size_t dw = 3;
    double h[3] = {1e-1,1e-2,2e-1};
    double hmin = h[1];

    double drift[3];
    double grad_drift[9];
    double diff[9];
    double grad_diff[27];

    double space[3];
    double prob[7];
    double prob2[7];
    double prob3[7];
    double grad_prob[7*3];
    double dt,dt2;
    double dt_grad[3];

    double grad_stage[3];
    double grad_bellman[3] = {0.0,0.0,0.0};

    double discount = 0.1;

    size_t nrand = 1000;
    for (size_t kk = 0; kk < nrand; kk++){
        /* printf("kk = %zu\n",kk); */
        /* double pt[2] = {1.0,2.0}; */
        double t = 0.2;
        double pt[3] = {randu()*3.0-1.5,
                        randu()*3.0-1.5,
                        randu()*3.0-1.5,
                       };
        double u[3] = {randu()*3.0-1.5,
                       randu()*3.0-1.5,
                       randu()*3.0-1.5};

        if (kk == 0){
            pt[0] = -0.951613;
            pt[1] = -1.305556;
            pt[2] = -2.615385;
            u[0] = -5.0;
            u[1] = -5.0;
            u[2] = 0.0;
        }

        double v[3] = {u[0],u[1],u[2]};
        
        double cost[7] = {randu(),
                          randu(),
                          randu(),
                          randu(),
                          randu(),
                          randu(),
                          randu()};

        double stage_cost;
        int res = stagecost3d(t,pt,u,&stage_cost,grad_stage);
        CuAssertIntEquals(tc,0,res);
    
        res = f2(t,pt,u,drift,grad_drift,NULL);
        CuAssertIntEquals(tc,0,res);
        res = s2(t,pt,u,diff,grad_diff,NULL);
        CuAssertIntEquals(tc,0,res);
        res = transition_assemble(dx,du,dw,hmin,h,drift,grad_drift,
                                  diff,grad_diff,prob,grad_prob,&dt,dt_grad,space);
        CuAssertIntEquals(tc,0,res==1); // res should not be 1

        
        bellmanrhs(dx,du,stage_cost,grad_stage,
                   discount, prob, grad_prob, dt, dt_grad,
                   cost, grad_bellman);
        /* printf("grad_bellman = "); dprint(3,grad_bellman); */

        /* printf("drift = (%G,%G,%G)\n",drift[0],drift[1],drift[2]); */
        /* printf("gdrift = (%G,%G,%G)\n",grad_drift[0],grad_drift[1],grad_drift[2]); */
        double delta = 1e-7;
        for (size_t zz = 0; zz < du; zz++){
            v[zz] = u[zz]+delta;
            res =f2(t,pt,v,drift,grad_drift,NULL);
            CuAssertIntEquals(tc,0,res);
            res = s2(t,pt,v,diff,grad_diff,NULL);
            CuAssertIntEquals(tc,0,res);
            res = transition_assemble(dx,du,dw,hmin,h,drift,NULL,diff,NULL,prob2,NULL,
                                      &dt2,NULL,NULL);

            double stage_cost2;
            res = stagecost3d(t,pt,v,&stage_cost2,NULL);
            CuAssertIntEquals(tc,0,res);
            double newval = bellmanrhs(dx,du,stage_cost2,NULL,
                                       discount, prob2,NULL, dt2, NULL,
                                       cost, NULL);

            v[zz] = u[zz];//-delta;
            res =f2(t,pt,v,drift,grad_drift,NULL);
            CuAssertIntEquals(tc,0,res);
            res = s2(t,pt,v,diff,grad_diff,NULL);
            CuAssertIntEquals(tc,0,res);
            res = transition_assemble(dx,du,dw,hmin,h,drift,NULL,diff,NULL,prob3,NULL,
                                      &dt2,NULL,NULL);

            double stage_cost3;
            res = stagecost3d(t,pt,v,&stage_cost3,NULL);
            CuAssertIntEquals(tc,0,res);
            double newval3 = bellmanrhs(dx,du,stage_cost3,NULL,
                                        discount, prob3,NULL, dt2, NULL,
                                        cost, NULL);

            /* double gstage = (stage_cost2 - stage_cost3)/2.0/delta; */
            double gstage = (stage_cost2 - stage_cost3)/delta;
            /* printf("\t grad stage should =  %G, grad is %G\n",gstage,grad_stage[zz]); */
            CuAssertDblEquals(tc,gstage,grad_stage[zz],1e-3);
            for (size_t ll = 0; ll < 2*dx+1; ll++){
                /* double gprob = (prob2[ll]-prob3[ll])/2.0/delta; */
                double gprob = (prob2[ll]-prob3[ll])/delta;
                double gprob_is = grad_prob[ll*du + zz];
                /* printf("\t\t grad prob should =  %G, grad is %G\n",gprob,gprob_is); */
                CuAssertDblEquals(tc,gprob,gprob_is,1e-3);
            }
            /* printf("val = %G, newval=%G\n",val,newval); */
            /* double fdgrad = (newval - newval3)/2.0/delta; */
            double fdgrad = (newval - newval3)/delta;
            /* printf("\t should = %G, is = %G\n",fdgrad,grad_bellman[zz]); */
            double diffe = fabs(fdgrad - grad_bellman[zz]);
            if (fabs(fdgrad) > 1){
                diffe /= fabs(fdgrad);
            }
            CuAssertDblEquals(tc,0,diffe,1e-3);

            v[zz] = u[zz];
        }

    }
}

double quad2d(const double * x, void * arg)
{
    (void)(arg);
    return x[0]*x[0] + x[1]*x[1];
}

void Test_bellman_control1d(CuTest * tc)
{
    printf("Testing Function: bellman_control (1d) \n");

    size_t dx = 2;
    size_t du = 1;
    size_t dw = 2;

    double lb[2] = {-1.0, -2.0};
    double ub[2] = {2.0, 3.0 };
    double center[2] = {0.0,0.0};
    double lengths[2] = {0.8,0.8};
    struct Boundary * bound = boundary_alloc(dx,lb,ub);
    boundary_add_obstacle(bound,center,lengths);
    

    size_t ngrid[2] = {63, 37};
    double * xgrid[2];
    xgrid[0] = linspace(lb[0],ub[0],ngrid[0]);
    xgrid[1] = linspace(lb[1],ub[1],ngrid[1]);
    double h[2] = {xgrid[0][1]-xgrid[0][0], xgrid[1][1]-xgrid[1][0]};
    double hmin = fmin(h[0],h[1]);
    struct MCAparam * mca = mca_param_create(dx,du);
    mca_add_grid_refs(mca,ngrid,xgrid,hmin,h);
    

    double discount = 0.1;
    struct DPparam * dp = dp_param_create(dx,du,dw,discount);
    dp_param_add_drift(dp,f1,NULL);
    dp_param_add_diff(dp,s1,NULL);
    dp_param_add_stagecost(dp,stagecost2d);
    dp_param_add_boundcost(dp,boundcost);
    dp_param_add_obscost(dp,ocost);


    // cross approximation tolerances
    size_t start_rank = 10;
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-8);
    approx_args_set_round_tol(aargs,1e-7);
    approx_args_set_kickrank(aargs,10);
    approx_args_set_adapt(aargs,0);
    approx_args_set_startrank(aargs,start_rank);
    approx_args_set_maxrank(aargs,start_rank);

    double * start[2];
    for (size_t ii = 0; ii < dx; ii++){
        start[ii] = calloc_double(start_rank);
        size_t stride = uniform_stride(ngrid[ii],start_rank);
        for (size_t jj = 0; jj < start_rank; jj++){
            start[ii][jj] = xgrid[ii][stride*jj];
        }
    }
    struct ValueF * vf = valuef_interp(dx,quad2d,NULL,ngrid,xgrid,start,aargs,0);


    double lbu = -5.0;
    double ubu = 5.0;
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,&lbu);
    c3opt_add_ub(opt,&ubu);
    c3opt_set_absxtol(opt,1e-15);
    c3opt_set_relftol(opt,1e-15);
    c3opt_set_gtol(opt,1e-30);
    c3opt_set_verbose(opt,0);

    struct ControlParams * cp = control_params_create(dx,dw,dp,mca,opt);
    for (size_t kk = 0; kk < ngrid[0]; kk++){
        for (size_t ll = 0; ll < ngrid[1]; ll++){
            size_t ind[2] = {kk,ll};
            for (size_t vary_ind = 0; vary_ind < dx; vary_ind++){

                /* size_t vary_ind = 0;; */
                double * x = calloc_double(dx*ngrid[vary_ind]);
                int * absorbed = calloc_int(ngrid[vary_ind]);
                double * costs = calloc_double(ngrid[vary_ind]*(2*dx+1));
                for (size_t ii = 0; ii < ngrid[vary_ind]; ii++){
                    for (size_t jj = 0; jj < dx; jj++){
                        x[ii*dx+jj] = xgrid[jj][ind[jj]];
                    }
                    x[ii*dx+vary_ind] = xgrid[vary_ind][ii];
                }

                /* printf("Hi\n"); */
                int res = mca_get_neighbor_costs(dx,ngrid[vary_ind],x,
                                                 bound,vf,ngrid,xgrid,
                                                 absorbed,costs);
                CuAssertIntEquals(tc,0,res);


                double time = 0.0;
                double val = 0.0;
                double val2 = 0.0;
                double val3 = 0.0;
                double delta = 1e-5;
                double grad[1];
                size_t nopts = 100;
                double * uopts = linspace(lbu,ubu,nopts);
                for (size_t ii = 0; ii < ngrid[vary_ind]; ii++){
                    /* printf("ii=%zu\n",ii); */
                    control_params_add_state_info(cp,time,x+ii*dx,absorbed[ii],costs+ii*(2*dx+1));

                    double umin = 1000;
                    double minim = 0.0;
                    for (size_t jj = 0; jj < nopts; jj++){
                        double u[1] = {uopts[jj]};
                        double v[1] = {u[0]+delta};
        
                        /* printf("control\n"); */
                        val = bellman_control(du,u,grad,cp);
                        if (val < umin){
                            umin = val;
                            minim = u[0];
                        }
                        
                        val2 = bellman_control(du,v,NULL,cp);
                        v[0] = u[0]-delta;
                        val3 = bellman_control(du,v,NULL,cp);

                        double diff_should = (val2-val3)/2.0/delta;
                        double diff_is = grad[0];
                
                        double err = (diff_is-diff_should);
                        if (fabs(diff_should) > 1e-13){
                            err /= diff_should;
                        }
                        CuAssertDblEquals(tc,0,err,1e-2);
                        /* printf("val=%G\n",val); */
                    }
                    if (absorbed[ii] == 0){ // this breaks!!! ofcourse...
                        double uopt[1];
                        double valopt;
                        int res2 = bellman_optimal(du,uopt,&valopt,cp);
                        CuAssertIntEquals(tc,0,res2);
                        /* printf("valopt = %3.10G, optbf=%3.10G\n",valopt,umin); */
                        /* printf("locminopt = %3.10G\n",uopt[0]); */
                        /* printf("locmin bf = %3.10G\n",minim); */
                        CuAssertIntEquals(tc,1,valopt < umin+1e-10);
                    }
                    /*     struct c3Opt * opt = c3opt_alloc(BFGS,du); */
                    /*     c3opt_add_lb(opt,&lbu); */
                    /*     c3opt_add_ub(opt,&ubu); */
                    /*     c3opt_set_absxtol(opt,1e-14); */
                    /*     c3opt_set_relftol(opt,1e-14); */
                    /*     c3opt_set_gtol(opt,1e-14); */
                    /*     c3opt_set_verbose(opt,0); */
                    /*     c3opt_add_objective(opt,&bellman_control,cp); */
                    /*     double ustart[1] = {-1.0}; */
                    /*     double valopt; */
                    /*     c3opt_minimize(opt,ustart,&valopt); */
                    /*     if (valopt > umin+1e-13){ */
                    /*         printf("x = "); dprint(dx,x+ii*dx); */
                    /*         printf("umin=%3.10G, valopt=%3.10G\n",umin,valopt); */
                    /*         printf("minimizer=%3.10G, optmin=%3.10G\n",minim,ustart[0]); */
                    /*         printf("not minimuzed\n"); */
                    /*         exit(1); */
                    /*     } */
                    /*     c3opt_free(opt); opt = NULL; */
                    /* } */
                }

                free(uopts); uopts = NULL;
                free(x); x = NULL;
                free(absorbed); absorbed = NULL;
                free(costs); costs = NULL;
            }
        }
    }
    /* printf("weird\n"); */
    /* double drift[2]; */
    /* double grad_drift[2]; */
    /* double diff[4]; */
    /* double grad_diff[4]; */

    control_params_destroy(cp); cp = NULL;
    boundary_free(bound); bound = NULL;
    mca_param_destroy(mca); mca = NULL;
    free(xgrid[0]); xgrid[0] = NULL;
    free(xgrid[1]); xgrid[1] = NULL;
    free(start[0]); start[0] = NULL;
    free(start[1]); start[1] = NULL;
    dp_param_destroy(dp); dp = NULL;
    valuef_destroy(vf); vf = NULL;
    approx_args_free(aargs); aargs = NULL;
    c3opt_free(opt); opt = NULL;
}

double quad3d(const double * x, void * arg)
{
    (void)(arg);
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

void Test_bellman_control3d(CuTest * tc)
{
    printf("Testing Function: bellman_control (3d) \n");

    size_t dx = 3;
    size_t du = 3;
    size_t dw = 3;

    double lb[3] = {-1.0, -2.0,-3.0};
    double ub[3] = {2.0, 3.0,1.0 };
    double center[3] = {0.0,0.0,0.0};
    double lengths[3] = {0.8,0.8,0.8};
    struct Boundary * bound = boundary_alloc(dx,lb,ub);
    boundary_add_obstacle(bound,center,lengths);
    

    size_t ngrid[3] = {13, 24, 41};
    double * xgrid[3];
    xgrid[0] = linspace(lb[0],ub[0],ngrid[0]);
    xgrid[1] = linspace(lb[1],ub[1],ngrid[1]);
    xgrid[2] = linspace(lb[2],ub[2],ngrid[2]);
    double h[3] = {xgrid[0][1]-xgrid[0][0],
                   xgrid[1][1]-xgrid[1][0],
                   xgrid[2][1]-xgrid[2][0]};
    double hmin = fmin(h[0],h[1]);
    hmin = fmin(hmin,h[2]);
    
    struct MCAparam * mca = mca_param_create(dx,du);
    mca_add_grid_refs(mca,ngrid,xgrid,hmin,h);
    
    double discount = 0.1;
    struct DPparam * dp = dp_param_create(dx,du,dw,discount);
    dp_param_add_drift(dp,f2,NULL);
    dp_param_add_diff(dp,s2,NULL);
    dp_param_add_stagecost(dp,stagecost3d);
    dp_param_add_boundcost(dp,boundcost);
    dp_param_add_obscost(dp,ocost);


    // cross approximation tolerances
    size_t start_rank = 10;
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-8);
    approx_args_set_round_tol(aargs,1e-7);
    approx_args_set_kickrank(aargs,10);
    approx_args_set_adapt(aargs,0);
    approx_args_set_startrank(aargs,start_rank);
    approx_args_set_maxrank(aargs,start_rank);

    double * start[3];
    for (size_t ii = 0; ii < dx; ii++){
        start[ii] = calloc_double(start_rank);
        size_t stride = uniform_stride(ngrid[ii],start_rank);
        for (size_t jj = 0; jj < start_rank; jj++){
            start[ii][jj] = xgrid[ii][stride*jj];
        }
    }
    struct ValueF * vf = valuef_interp(dx,quad3d,NULL,ngrid,xgrid,start,aargs,0);

    double lbu = -5.0;
    double ubu = 5.0;
    double lbarr[3] = {lbu,lbu,lbu};
    double ubarr[3] = {ubu,ubu,ubu};
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbarr);
    c3opt_add_ub(opt,ubarr);
    c3opt_set_absxtol(opt,1e-15);
    c3opt_set_relftol(opt,1e-15);
    c3opt_set_gtol(opt,0.0);
    c3opt_set_verbose(opt,0);
    c3opt_ls_set_beta(opt,0.8);

    struct ControlParams * cp = control_params_create(dx,dw,dp,mca,opt);
    for (size_t kk = 0; kk < ngrid[0]; kk = kk + 5){
        for (size_t ll = 0; ll < ngrid[1]; ll = ll + 5){
            for (size_t zz = 0; zz < ngrid[2]; zz = zz + 5){
                printf("(%zu/%zu), (%zu/%zu), (%zu/%zu)\n",
                       kk+1,ngrid[0],ll+1,ngrid[1],zz+1,ngrid[2]);
                size_t ind[3] = {kk,ll,zz};
                for (size_t vary_ind = 0; vary_ind < dx; vary_ind++){

                    /* size_t vary_ind = 0;; */
                    double * x = calloc_double(dx*ngrid[vary_ind]);
                    int * absorbed = calloc_int(ngrid[vary_ind]);
                    double * costs = calloc_double(ngrid[vary_ind]*(2*dx+1));
                    for (size_t ii = 0; ii < ngrid[vary_ind]; ii++){
                        for (size_t jj = 0; jj < dx; jj++){
                            x[ii*dx+jj] = xgrid[jj][ind[jj]];
                        }
                        x[ii*dx+vary_ind] = xgrid[vary_ind][ii];
                    }

                    /* printf("Hi\n"); */
                    int res = mca_get_neighbor_costs(dx,ngrid[vary_ind],x,
                                                     bound,vf,ngrid,xgrid,
                                                     absorbed,costs);
                    CuAssertIntEquals(tc,0,res);


                    double time = 0.0;
                    double val = 0.0;
                    double val2 = 0.0;
                    double val3 = 0.0;
                    double delta = 1e-7;
                    double grad[3];
                    size_t nopts = 10;
                    /* double lbu = -5.0; */
                    /* double ubu = 5.0; */
                    double * uopts = linspace(lbu+5e-1,ubu-5e-1,nopts);
                    for (size_t ii = 0; ii < ngrid[vary_ind]; ii++){
                        /* printf("ii=%zu\n",ii); */
                        control_params_add_state_info(cp,time,x+ii*dx,absorbed[ii],
                                                      costs+ii*(2*dx+1));

                        double umin = 1000;
                        double minim[3] = {0.0, 0.0, 0.0};
                        for (size_t jj = 0; jj < nopts; jj++){
                            for (size_t ff = 0; ff < nopts; ff++ ){
                                for (size_t qq = 0; qq < nopts; qq++){
                                    
                                    double u[3] = {uopts[jj],uopts[ff],uopts[qq]};
                                                                /* printf("control\n"); */
                                    val = bellman_control(du,u,grad,cp);
                                    if (val < umin){
                                        umin = val;
                                        minim[0] = u[0];
                                        minim[1] = u[1];
                                        minim[2] = u[2];
                                    }


                                    double v[3] = {u[0],u[1],u[2]};
                                    for (size_t rr = 0; rr < du; rr++){
                                        v[rr] = u[rr]+delta;
                                        val2 = bellman_control(du,v,NULL,cp);
                                        v[rr] = u[rr]-delta;
                                        val3 = bellman_control(du,v,NULL,cp);
                                        v[rr] = u[rr];
                                        double diff_should = (val2-val3)/2.0/delta;
                                        /* double diff_should = (val2-val3)/delta; */
                                        double diff_is = grad[rr];
                                        double err = (diff_is-diff_should);


                                        /* if (fabs(diff_should) > 1e-13){ */
                                        /*     err /= diff_should; */
                                        /* } */
                                        if (fabs(err) > 1e-2){ // make sure its due to weird deriv
                                                               // due to 0 drift

                                            int last_res = control_params_get_last_res(cp);
                                            /* printf("last_res = %d\n",last_res); */
                                            /* printf("u = "); dprint(3,u); */
                                            /* printf("err = %G\n",err); */
                                            /* printf("diff_should = %G\n",diff_should); */
                                            /* printf("diff_is = %G\n",diff_is); */
                                            /* printf("rr = %zu\n", rr); */
                                            /* printf("x = "); dprint(dx,x+ii*dx); */
                                            /* printf("absorbed = %d\n",absorbed[ii]); */

                                            /* double drift[3]; */
                                            /* f2(time,x+ii*dx,u,drift,NULL,NULL); */
                                            /* printf("drift = "); */
                                            /* for (size_t llf = 0; llf < dx; llf++){ */
                                            /*     printf("%3.10G ",drift[llf]); */
                                            /* } */
                                            /* printf("\n"); */
                                            CuAssertIntEquals(tc,2,last_res);
                                        }
                                        /* CuAssertDblEquals(tc,0,err,1e-2); */
                                    }
                                }
                            }
                        }
                        /* if (absorbed[ii] == 0){ // this breaks!!! ofcourse... */
                        /*     double uopt[3]; */
                        /*     double valopt; */
                        /*     /\* printf("\n\n\n optimize\n"); *\/ */
                        /*     /\* printf("\toptimize\n"); *\/ */
                        /*     /\* c3opt_set_verbose(opt,2); *\/ */
                        /*     int res2 = bellman_optimal(du,uopt,&valopt,cp); */
                        /*     CuAssertIntEquals(tc,0,res2); */
                        /*     if (valopt > umin){ */
                        /*         int last_res = control_params_get_last_res(cp); */
                        /*         printf("last_res = %d\n",last_res); */
                        /*         printf("valopt = %3.10G, optbf=%3.10G\n",valopt,umin); */
                        /*         printf("\t locminopt = %3.10G,%3.10G,%3.10G\n", */
                        /*                uopt[0],uopt[1],uopt[2]); */
                        /*         printf("\t locbf     = %3.10G,%3.10G,%3.10G\n", */
                        /*                minim[0],minim[1],minim[2]); */
                        /*     } */
                        /*     CuAssertIntEquals(tc,1,valopt < umin+1e-10); */
                        /* } */
                    }

                    free(uopts); uopts = NULL;
                    free(x); x = NULL;
                    free(absorbed); absorbed = NULL;
                    free(costs); costs = NULL;
                }
            }
        }
    }
    /* printf("weird\n"); */
    /* double drift[2]; */
    /* double grad_drift[2]; */
    /* double diff[4]; */
    /* double grad_diff[4]; */

    control_params_destroy(cp); cp = NULL;
    boundary_free(bound); bound = NULL;
    mca_param_destroy(mca); mca = NULL;
    free(xgrid[0]); xgrid[0] = NULL;
    free(xgrid[1]); xgrid[1] = NULL;
    free(xgrid[2]); xgrid[2] = NULL;
    free(start[0]); start[0] = NULL;
    free(start[1]); start[1] = NULL;
    free(start[2]); start[2] = NULL;
    dp_param_destroy(dp); dp = NULL;
    valuef_destroy(vf); vf = NULL;
    approx_args_free(aargs); aargs = NULL;
    c3opt_free(opt); opt = NULL;
}


CuSuite * BellmanGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    /* SUITE_ADD_TEST(suite, Test_process_fibers_neighbor); */
    /* SUITE_ADD_TEST(suite, Test_bound_nodes); */
    /* SUITE_ADD_TEST(suite, Test_bellman_grad1); */
    /* SUITE_ADD_TEST(suite, Test_bellman_grad3d); */
    /* SUITE_ADD_TEST(suite, Test_bellman_control1d); */
    SUITE_ADD_TEST(suite, Test_bellman_control3d);

    return suite;
}

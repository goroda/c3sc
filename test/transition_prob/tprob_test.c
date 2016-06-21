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
#include "c3sc.h"

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

/* void fd_grad(size_t d, double * u, double h, double * out, */
/*              double (*f)(double * u, void *), void * arg ) */
/* { */
/*     double * utemp = calloc_double(d); */
/*     memmove(utemp,u,d*sizeof(double)); */
/*     double val1 = f(u,arg); */
/*     for (size_t ii = 0; ii < d; ii++){ */
/*         utemp[ii] += h */
/*         double val2 = f(utemp,arg); */
/*         utemp[ii] -= h; */
/*         out[ii] = (val2 - val1)/h; */
/*     } */
/*     free (utemp); */
/* } */

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
        //df1/du
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

void Test_cost_neighbor_eval(CuTest * tc)
{
    printf("Testing Functions associated with evaluating neighbors \n");

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
    size_t nevals = N[dim_vary]*2*(dim-1);
    double * vals = calloc_double(nevals);
    double * fiber_vals = calloc_double(N[dim_vary]);
    cost_eval_fiber_ind_nn(vf,fixed_ind,dim_vary,neighbors,vals,fiber_vals);

    /* printf("size think = %zu\n",nevals); */
    /* dprint(nevals,vals); */
    for (size_t jj = 0; jj < N[dim_vary]; jj++){
       
        // no neighbors
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][fixed_ind[2]];
        
        double val_should = valuef_eval(vf,pt);
        double val_is = fiber_vals[jj];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        size_t ii = 0;
        // first point
        pt[0] = xgrid[0][neighbors[0]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][neighbors[1]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        ii = 1;
        // second point
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][neighbors[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][jj];
        pt[2] = xgrid[2][neighbors[3]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
    }
    free(vals); vals = NULL;
    free(fiber_vals); fiber_vals = NULL;


    dim_vary = 0;
    nevals = N[dim_vary]*2*(dim-1);
    vals = calloc_double(nevals);
    fiber_vals = calloc_double(N[dim_vary]);
    cost_eval_fiber_ind_nn(vf,fixed_ind,dim_vary,neighbors,vals,fiber_vals);
    /* printf("size think = %zu\n",nevals); */
    /* dprint(nevals,vals); */
    for (size_t jj = 0; jj < N[dim_vary]; jj++){

        // no neighbors
        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][fixed_ind[2]];
        
        double val_should = valuef_eval(vf,pt);
        double val_is = fiber_vals[jj];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        size_t ii = 0;
        // first point
        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][neighbors[0]];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][neighbors[1]];
        pt[2] = xgrid[2][fixed_ind[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        ii = 1;
        // second point
        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][neighbors[2]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][jj];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][neighbors[3]];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
    }
    free(vals); vals = NULL;
    free(fiber_vals); fiber_vals = NULL;

    dim_vary = 2;
    nevals = N[dim_vary]*2*(dim-1);
    vals = calloc_double(nevals);
    fiber_vals = calloc_double(N[dim_vary]);
    cost_eval_fiber_ind_nn(vf,fixed_ind,dim_vary,neighbors,vals,fiber_vals);
    /* printf("size think = %zu\n",nevals); */
    /* dprint(nevals,vals); */
    for (size_t jj = 0; jj < N[dim_vary]; jj++){

        // no neighbors
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][jj];
        
        double val_should = valuef_eval(vf,pt);
        double val_is = fiber_vals[jj];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        // first point
        size_t ii = 0;

        pt[0] = xgrid[0][neighbors[0]];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][jj];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        pt[0] = xgrid[0][neighbors[1]];
        pt[1] = xgrid[1][fixed_ind[1]];
        pt[2] = xgrid[2][jj];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        ii = 1;
        // second point
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][neighbors[2]];
        pt[2] = xgrid[2][jj];
        
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);

        // second point
        pt[0] = xgrid[0][fixed_ind[0]];
        pt[1] = xgrid[1][neighbors[3]];
        pt[2] = xgrid[2][jj];
            
        val_should = valuef_eval(vf,pt);
        val_is = vals[jj*2*(dim-1) + 2*ii+1];
        /* printf("f(%G,%G,%G) = %G\n",pt[0],pt[1],pt[2],val_should); */
        /* printf("\t val=%G\n",val_is); */
        CuAssertDblEquals(tc,val_should,val_is,1e-14);
    }
    free(vals); vals = NULL;
    free(fiber_vals); fiber_vals = NULL;


    /* printf("done sdaskd\n"); */

    valuef_destroy(vf); vf = NULL;
    for (size_t ii = 0; ii < dim; ii++){
        free(xgrid[ii]); xgrid[ii] = NULL;
        free(start[ii]); start[ii] = NULL;
    }
    approx_args_free(aargs); aargs = NULL;

}

CuSuite * ValueFGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_cost_neighbor_eval);

    return suite;
}

void Test_bellman_grad1(CuTest * tc)
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

void Test_bound_nodes(CuTest * tc)
{
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

CuSuite * BellmanGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    /* SUITE_ADD_TEST(suite, Test_bellman_grad1); */
    SUITE_ADD_TEST(suite, Test_bound_nodes);

    return suite;
}

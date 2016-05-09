// Copyright (c) 2016, Massachusetts Institute of Technology
// Authors: Alex Gorodetsky
// Email  : goroda@mit.edu

//Code

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CuTest.h"	
#include "tensmarkov.h"

#define dglob 5

void Test_mcnode_alloc(CuTest * tc)
{
    printf("Testing Function: mcnode_alloc,mcnode_free \n");

    struct MCNode * mcn = NULL;
    mcn = mcnode_alloc(dglob);
    CuAssertIntEquals(tc,1,mcn != NULL);
    mcnode_free(mcn); mcn = NULL;
    CuAssertIntEquals(tc,1,mcn == NULL);
}

void Test_mcnode_init(CuTest * tc)
{
    printf("Testing Function: mcnode_init \n");

    struct MCNode * mcn = NULL;
    double x[dglob] = {1.0,2.0,3.0,4.0,5.0};
    mcn = mcnode_init(dglob,x);

    size_t dg = mcnode_get_d(mcn);
    CuAssertIntEquals(tc,dglob,dg);

    double * xt = mcnode_get_xref(mcn);
    double diff = norm2diff(x,xt,5);
    CuAssertDblEquals(tc,0.0,diff,1e-15);

    double pself = mcnode_get_pself(mcn);
    CuAssertDblEquals(tc,1.0,pself,1e-15);

    mcnode_free(mcn); mcn = NULL;
}

void Test_mcnode_prepend(CuTest * tc)
{
    printf("Testing Function: mcnlist_prepend (with no BondaryInfo)\n");


    struct MCNode * mcn = NULL;
    double x[dglob] = {1.0,2.0,3.0,4.0,5.0};
    double h[dglob] = {1.0,0.5,0.25,0.25,0.2};
    double pselfin = 0.2;
    double probs[dglob] = {0.2,0.2,0.2,0.1,0.1};

    mcn = mcnode_init(dglob,x);
    mcnode_set_pself(mcn,pselfin);

    for (size_t ii = 0; ii < dglob;ii++){
        mcnode_prepend_neigh(mcn,ii,x[ii]-h[ii],probs[ii],NULL);
    }

    double pself = mcnode_get_pself(mcn);
    CuAssertDblEquals(tc,pselfin,pself,1e-15);

    struct MCNList * temp = mcnode_get_neigh(mcn);
    int on = 4;
    while (temp != NULL){
        size_t dir = mcnlist_get_dir(temp);
        double xval = mcnlist_get_val(temp);
        CuAssertIntEquals(tc,on,dir);
        CuAssertDblEquals(tc,x[on]-h[on],xval,1e-15);
        temp = mcnlist_get_next(temp);
        on--;
    }
    
    mcnode_free(mcn); mcn = NULL;
}

double expfunc(double t, const double * x, void * args)
{
    (void)(t);
    (void)(args);
    double out = x[0] * x[1] + x[3] * sin(x[1]);
    return out;
}

void Test_mcnode_expectation(CuTest * tc)
{
    printf("Testing Function: mcnode_expectation\n");

    struct MCNode * mcn = NULL;
    double x[5] = {1.0,2.0,3.0,4.0,5.0};
    double h[5] = {1.0,0.5,0.25,0.25,0.2};
    double pselfin = 0.2;
    double probs[5] = {0.2,0.2,0.2,0.1,0.1};
    size_t du = 4;
    double ** gprob = malloc_dd(dglob);
    for (size_t ii = 0; ii < dglob; ii++){
        gprob[ii] = calloc_double(du);
        for (size_t jj = 0; jj < du; jj++){
            gprob[ii][jj] = randu();
        }
    }
    double * gprobself = calloc_double(du);
    for (size_t jj = 0; jj < du; jj++){
        gprobself[jj] = randu();
    }

    mcn = mcnode_init(dglob,x);
    mcnode_set_pself(mcn,pselfin);
    mcnode_set_gradient(mcn,du,gprobself);

    for (size_t ii = 0; ii < dglob; ii++){
        mcnode_prepend_neigh(mcn,ii,x[ii]-h[ii],probs[ii],gprob[ii]);
    }
    
    double avg = mcnode_expectation(mcn,expfunc,NULL,NULL);
    
    double xtemp[5];
    memmove(xtemp,x,5*sizeof(double));
    double as = 0.0;
    for (size_t ii = 0; ii < 5; ii++){
        xtemp[ii] = x[ii] - h[ii];
        as += probs[ii] * expfunc(0,xtemp,NULL);
        xtemp[ii] = x[ii];
    }
    as += pselfin * expfunc(0,x,NULL);
            
    CuAssertDblEquals(tc,as,avg,1e-14);
    
    free_dd(dglob,gprob);
    free(gprobself);
    mcnode_free(mcn); mcn = NULL;
}

int f1(double t, const double * x, const double * u, double * out,
       double * jac, void * args)
{
    (void)(t);
    (void)(args);
    out[0] = x[1];
    out[1] = u[0];
    if (jac != NULL){
        //df1/du
        jac[0] = 0.0;
        jac[1] = 1.0;
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
    out[0] = 1.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 1.0;
    if (grad != NULL){
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
    }
    return 0;
}

void Test_mca_get_node_absorb(CuTest * tc)
{
    printf("Testing Function: mca_get_node w/ absorbing conditions(1)\n");

    size_t d = 2;
    size_t du = 1;
    size_t dw = 2;
    double h[2] = {1e-1, 1e-2};
    double lb[2] = {-2.0, -2.0};
    double ub[2] = {2.0, 2.0};

    struct Boundary * bound = boundary_alloc(d,lb,ub);
    struct Drift * drift = drift_alloc(d,du);
    drift_add_func(drift,f1,NULL);
    struct Diff * diff = diff_alloc(d,du,dw);
    diff_add_func(diff,s1,NULL);
    struct Dyn  * dyn = dyn_alloc(drift,diff);

    struct MCA * mca = mca_alloc(d,du,dw,h);
    mca_attach_dyn(mca,dyn);
    mca_attach_bound(mca,bound);

    struct BoundInfo * bi = NULL;
    
    // first check corner point that is absorbing
    double pt[2] = {-2.0,-2.0};
    double u[1] = {0.0};
    double time = 0.0;

    double dt; 
    struct MCNode * mcn = mca_get_node(mca,time,pt,u,&dt,NULL,&bi,NULL);
    int onbound = bound_info_onbound(bi);
    CuAssertIntEquals(tc,1,onbound);
    int absorb = bound_info_absorb(bi);
    CuAssertIntEquals(tc,1,absorb);
    int period = bound_info_period(bi);
    CuAssertIntEquals(tc,0,period);
    int reflect = bound_info_reflect(bi);
    CuAssertIntEquals(tc,0,reflect);
    int N = mcnode_get_n(mcn);
    CuAssertIntEquals(tc,0,N);
    double p = mcnode_get_pself(mcn);
    CuAssertDblEquals(tc,1.0,p,1e-15);
    mcnode_free(mcn);
    bound_info_free(bi);

    // check non-corner but absorbing point
    pt[0] = 0.4; pt[1] = 2.0;
    mcn = mca_get_node(mca,time,pt,u,&dt,NULL,&bi,NULL);
    onbound = bound_info_onbound(bi);
    CuAssertIntEquals(tc,1,onbound);
    absorb = bound_info_absorb(bi);
    CuAssertIntEquals(tc,1,absorb);
    period = bound_info_period(bi);
    CuAssertIntEquals(tc,0,period);
    reflect = bound_info_reflect(bi);
    CuAssertIntEquals(tc,0,reflect);
    N = mcnode_get_n(mcn);
    CuAssertIntEquals(tc,0,N);
    p = mcnode_get_pself(mcn);
    CuAssertDblEquals(tc,1.0,p,1e-15);
    mcnode_free(mcn);
    bound_info_free(bi);

    // check non-absorbing point
    pt[0] = 0.4; pt[1] = 0.8;
    double evald[2]; 
    drift_eval(drift,time,pt,u,evald,NULL);
    double evaldiff[4];
    diff_eval(diff,time,pt,u,evaldiff,NULL);

    
    double Qh = h[1]*(h[1]/h[0]*fabs(evald[0]) + fabs(evald[1]));
    Qh +=  evaldiff[0]*evaldiff[0]*h[1]*h[1]/h[0]/h[0] + 
           evaldiff[3]*evaldiff[3];
//    printf("Qh should = %G\n",Qh);

    double dtshould = h[1]*h[1] / Qh;
    mcn = mca_get_node(mca,time,pt,u,&dt,NULL,&bi,NULL);
    CuAssertDblEquals(tc,dtshould,dt,1e-14);
    onbound = bound_info_onbound(bi);
    CuAssertIntEquals(tc,0,onbound);
    absorb = bound_info_absorb(bi);
    CuAssertIntEquals(tc,0,absorb);
    period = bound_info_period(bi);
    CuAssertIntEquals(tc,0,period);
    reflect = bound_info_reflect(bi);
    CuAssertIntEquals(tc,0,reflect);
    N = mcnode_get_n(mcn);
    CuAssertIntEquals(tc,4,N);
   
    // left,right, up, down
    double p1=0.0,p2=0.0,p3=0.0,p4=0.0; 
    if (evald[0] > 0){
        p2 = h[1]*h[1]/h[0]*evald[0] + pow(evaldiff[0]*h[1]/h[0],2)/2.0;
        p1 = pow(evaldiff[0]*h[1]/h[0],2)/2.0;
    }
    else{
        p1 = -h[1]*h[1]/h[0]*evald[0] + pow(evaldiff[0]*h[1]/h[0],2)/2.0;
        p2 = pow(evaldiff[0]*h[1]/h[0],2)/2.0;
    }
    
    if (evald[1] > 0){
        p3 = h[1]*evald[1] + pow(evaldiff[3],2)/2.0;
        p4 = pow(evaldiff[3],2)/2.0;
    }
    else{
        p4 = -h[1]*evald[1] + pow(evaldiff[3],2)/2.0;
        p3 = pow(evaldiff[3],2)/2.0;
    }

    p1 /= Qh;
    p2 /= Qh;
    p3 /= Qh;
    p4 /= Qh;

    struct MCNList * temp = mcnode_get_neigh(mcn);
    while (temp != NULL){
        size_t dir = mcnlist_get_dir(temp);
        double val = mcnlist_get_val(temp);
        double pp = mcnlist_get_p(temp);
        if (dir == 0){
            if (val > pt[dir]){
                CuAssertDblEquals(tc,val,pt[0]+h[0],1e-15);
                CuAssertDblEquals(tc,p2,pp,1e-15);
            }
            else{
                CuAssertDblEquals(tc,val,pt[0]-h[0],1e-15);
                CuAssertDblEquals(tc,p1,pp,1e-15);
            }
        }
        if (dir == 1){
            if (val > pt[dir]){
                CuAssertDblEquals(tc,val,pt[1]+h[1],1e-15);
                CuAssertDblEquals(tc,p3,pp,1e-15);
            }
            else{
                CuAssertDblEquals(tc,val,pt[1]-h[1],1e-15);
                CuAssertDblEquals(tc,p4,pp,1e-15);
            }
        }
        temp = mcnlist_get_next(temp);
    }

    /* p = mcnode_get_pself(mcn); */
    /* CuAssertDblEquals(tc,1.0,p,1e-15); */
    mcnode_free(mcn);
    bound_info_free(bi);

    mca_free_deep(mca);
}

void Test_mca_get_node_absorb_grad(CuTest * tc)
{
    printf("Testing Function: mca_get_node w/ absorbing conditions (2)\n");

    size_t d = 2;
    size_t du = 1;
    size_t dw = 2;
    double h[2] = {1e-1, 1e-2};
    double lb[2] = {-2.0, -2.0};
    double ub[2] = {2.0, 2.0};

    struct Boundary * bound = boundary_alloc(d,lb,ub);
    struct Drift * drift = drift_alloc(d,du);
    drift_add_func(drift,f1,NULL);
    struct Diff * diff = diff_alloc(d,du,dw);
    diff_add_func(diff,s1,NULL);
    struct Dyn  * dyn = dyn_alloc(drift,diff);

    struct MCA * mca = mca_alloc(d,du,dw,h);
    mca_attach_dyn(mca,dyn);
    mca_attach_bound(mca,bound);

    struct BoundInfo * bi = NULL;
    
    // first check corner point that is absorbing
    double pt[2] = {-2.0,-2.0};
    double u[1] = {0.0};
    double time = 0.0;

    double dt; 
    double gdt[1];
    double grad[1];
    struct MCNode * mcn = mca_get_node(mca,time,pt,u,&dt,gdt,&bi,grad);
    double * gself = mcnode_get_gpself(mcn);
    CuAssertDblEquals(tc,0.0,gself[0],1e-15);
    mcnode_free(mcn);
    bound_info_free(bi);
    
//    printf("doesnt get here!\n");
    // check non-corner but absorbing point
    pt[0] = 0.4; pt[1] = 2.0;
    mcn = mca_get_node(mca,time,pt,u,&dt,gdt,&bi,grad);
    gself = mcnode_get_gpself(mcn);
    CuAssertDblEquals(tc,0.0,gself[0],1e-15);
    mcnode_free(mcn);
    bound_info_free(bi);

    // check non-absorbing point
    pt[0] = 0.4; pt[1] = 0.8;
    double evald[2];
    double evaldg[2];
    drift_eval(drift,time,pt,u,evald,evaldg);
    double evaldiff[4];
    diff_eval(diff,time,pt,u,evaldiff,NULL);

    double Qh = h[1]*(h[1]/h[0]*fabs(evald[0]) + fabs(evald[1]));
    Qh +=  evaldiff[0]*evaldiff[0]*h[1]*h[1]/h[0]/h[0] +
           evaldiff[3]*evaldiff[3];

    double dQ;
    if (evald[0] > 0){
        dQ = h[1]*h[1]/h[0] * evaldg[0];
    }
    else{
        dQ = -h[1]*h[1]/h[0] * evaldg[0];
    }
    if (evald[1] > 0){
        dQ = h[1] * evaldg[1];
    }
    else{
        dQ = -h[1] * evaldg[1];
    }

    //printf("dQshould = %G\n",dQ);
    double dtshould = h[1]*h[1] / Qh;
    //printf("dtshould = %G\n",dtshould);
    double gdtshould = - h[1]*h[1]/ Qh/Qh * dQ;
    mcn = mca_get_node(mca,time,pt,u,&dt,gdt,&bi,grad);
    gself = mcnode_get_gpself(mcn);
    CuAssertDblEquals(tc,dtshould,dt,1e-14);
    CuAssertDblEquals(tc,gdtshould,gdt[0],1e-14);

    // left,right, up, down
    double p1=0.0,p2=0.0,p3=0.0,p4=0.0; 
    if (evald[0] > 0){
        p2 = h[1]*h[1]/h[0]*evald[0] + pow(evaldiff[0]*h[1]/h[0],2)/2.0;
        p1 = pow(evaldiff[0]*h[1]/h[0],2)/2.0;
    }
    else{
        p1 = -h[1]*h[1]/h[0]*evald[0] + pow(evaldiff[0]*h[1]/h[0],2)/2.0;
        p2 = pow(evaldiff[0]*h[1]/h[0],2)/2.0;
    }
    
    if (evald[1] > 0){
        p3 = h[1]*evald[1] + pow(evaldiff[3],2)/2.0;
        p4 = pow(evaldiff[3],2)/2.0;
    }
    else{
        p4 = -h[1]*evald[1] + pow(evaldiff[3],2)/2.0;
        p3 = pow(evaldiff[3],2)/2.0;
    }
   
    //gradients of left,right,up,down
    double dp1=0.0,dp2=0.0,dp3=0.0,dp4=0.0;
    if (evald[0] > 0){
        dp2 = h[1]*h[1]/h[0]*evaldg[0];
        dp1 = 0.0; 
    }
    else{
        dp1 = -h[1]*h[1]/h[0]*evaldg[0];
        dp2 = 0.0;
    }
    if (evald[1] > 0){
        dp3 = h[1]*evaldg[1];
        dp4 = 0.0;
    }
    else{
        dp4 = -h[1]*evaldg[1];
        dp3 = 0.0;
    }
    
    dp1 = (Qh*dp1 - p1*dQ)/Qh/Qh;
    dp2 = (Qh*dp2 - p2*dQ)/Qh/Qh;
    dp3 = (Qh*dp3 - p3*dQ)/Qh/Qh;
    dp4 = (Qh*dp4 - p4*dQ)/Qh/Qh;

    struct MCNList * temp = mcnode_get_neigh(mcn);
    while (temp != NULL){
        size_t dir = mcnlist_get_dir(temp);
        double val = mcnlist_get_val(temp);
        double * gradp = mcnlist_get_gradp(temp);
        if (dir == 0){
            if (val > pt[dir]){
                CuAssertDblEquals(tc,val,pt[0]+h[0],1e-15);
                CuAssertDblEquals(tc,dp2,gradp[0],1e-15);
            }
            else{
                CuAssertDblEquals(tc,val,pt[0]-h[0],1e-15);
                CuAssertDblEquals(tc,dp1,gradp[0],1e-15);
            }
        }
        if (dir == 1){
            if (val > pt[dir]){
                CuAssertDblEquals(tc,val,pt[1]+h[1],1e-15);
                CuAssertDblEquals(tc,dp3,gradp[0],1e-15);
            }
            else{
                CuAssertDblEquals(tc,val,pt[1]-h[1],1e-15);
                CuAssertDblEquals(tc,dp4,gradp[0],1e-15);
            }
        }
        temp = mcnlist_get_next(temp);
    }

    mcnode_free(mcn);
    bound_info_free(bi);

    mca_free_deep(mca);
}

void Test_mca_get_node_reflect(CuTest * tc)
{
    printf("Testing Function: mca_get_node w/ reflecting conditions(1)\n");

    size_t d = 2;
    size_t du = 1;
    size_t dw = 2;
    double h[2] = {1e-1, 1e-2};
    double lb[2] = {-2.0, -2.0};
    double ub[2] = {2.0, 2.0};

    struct Boundary * bound = boundary_alloc(d,lb,ub);
    boundary_external_set_type(bound,0,"reflect");
    boundary_external_set_type(bound,1,"reflect");

    struct Drift * drift = drift_alloc(d,du);
    drift_add_func(drift,f1,NULL);
    struct Diff * diff = diff_alloc(d,du,dw);
    diff_add_func(diff,s1,NULL);
    struct Dyn  * dyn = dyn_alloc(drift,diff);

    struct MCA * mca = mca_alloc(d,du,dw,h);
    mca_attach_dyn(mca,dyn);
    mca_attach_bound(mca,bound);

    struct BoundInfo * bi = NULL;
    
    // first check that corner point that is reflecting
    double pt[2] = {-2.0,-2.0};
    double u[1] = {0.0};
    double time = 0.0;
    double evald[2]; 
    drift_eval(drift,time,pt,u,evald,NULL);
    double evaldiff[4];
    diff_eval(diff,time,pt,u,evaldiff,NULL);

    double dt;
    struct MCNode * mcn = mca_get_node(mca,time,pt,u,&dt,NULL,&bi,NULL);
    int onbound = bound_info_onbound(bi);
    CuAssertIntEquals(tc,1,onbound);
    int absorb = bound_info_absorb(bi);
    CuAssertIntEquals(tc,0,absorb);
    int period = bound_info_period(bi);
    CuAssertIntEquals(tc,0,period);
    int reflect = bound_info_reflect(bi);
    CuAssertIntEquals(tc,1,reflect);
    int N = mcnode_get_n(mcn);
    CuAssertIntEquals(tc,2,N);
    struct MCNList * temp = mcnode_get_neigh(mcn);
    double sum_prob = 0.0;
    while (temp != NULL){
        /* size_t dir = mcnlist_get_dir(temp); */
        /* double val = mcnlist_get_val(temp); */
        double p = mcnlist_get_p(temp);
//        printf("new pt[%zu] = %G p_trans=%G\n",dir,val,p);
        sum_prob += p;
        temp = mcnlist_get_next(temp);
    }
    /* double p1; // right; */
    /* double p2; // down; */

    /* double Qh = h[1]*(h[1]/h[0]*fabs(evald[0]) + fabs(evald[1])); */
    /* Qh +=  evaldiff[0]*evaldiff[0]*h[1]*h[1]/h[0]/h[0] +  */
    /*        evaldiff[3]*evaldiff[3]; */

    double p = mcnode_get_pself(mcn);
    CuAssertIntEquals(tc,1,p>=0);
    sum_prob += p;
    CuAssertDblEquals(tc,1.0,sum_prob,1e-14);
//    CuAssertDblEquals(tc,1.0,p,1e-15);
    mcnode_free(mcn);
    bound_info_free(bi);

    // check derivative
    double gradc[1];
    double gdt[1];
    double gradc2[1];
    double gdt2[1];
    u[0] = 2.0;
    struct MCNode * mcn1 = mca_get_node(mca,time,pt,u,&dt,gdt,&bi,gradc);
    double v1 = mcnode_expectation(mcn1,expfunc,NULL,gradc);
    bound_info_free(bi);
    u[0] = u[0]+1e-4;
    double dt2;
    struct MCNode * mcn2 = mca_get_node(mca,time,pt,u,&dt2,gdt2,&bi,gradc2);
    double v2 = mcnode_expectation(mcn2,expfunc,NULL,gradc2);
    double grad_diff = (v2-v1)/1e-4;
    double gdt_diff = (dt2-dt)/1e-4;
    /* printf("%G \n",grad_diff-gradc[0]); */
    /* printf("%G \n",gdt_diff-gdt[0]); */
    CuAssertDblEquals(tc,grad_diff,gradc[0],1e-3);
    CuAssertDblEquals(tc,gdt_diff,gdt[0],1e-3);
    bound_info_free(bi);
    mcnode_free(mcn1);
    mcnode_free(mcn2);
    

    // check non-corner but reflecting point
    pt[0] = 0.4; pt[1] = 2.0;
    mcn = mca_get_node(mca,time,pt,u,&dt,NULL,&bi,NULL);
    onbound = bound_info_onbound(bi);
    CuAssertIntEquals(tc,1,onbound);
    absorb = bound_info_absorb(bi);
    CuAssertIntEquals(tc,0,absorb);
    period = bound_info_period(bi);
    CuAssertIntEquals(tc,0,period);
    reflect = bound_info_reflect(bi);
    CuAssertIntEquals(tc,1,reflect);
    N = mcnode_get_n(mcn);
    CuAssertIntEquals(tc,3,N);
    p = mcnode_get_pself(mcn);
    CuAssertIntEquals(tc,1,p>=0);
    /* CuAssertDblEquals(tc,1.0,p,1e-15); */
    mcnode_free(mcn);
    bound_info_free(bi);

    // check non-areflecting point
    pt[0] = 0.4; pt[1] = 0.8;
    mcn = mca_get_node(mca,time,pt,u,&dt,NULL,&bi,NULL);
    onbound = bound_info_onbound(bi);
    CuAssertIntEquals(tc,0,onbound);
    absorb = bound_info_absorb(bi);
    CuAssertIntEquals(tc,0,absorb);
    period = bound_info_period(bi);
    CuAssertIntEquals(tc,0,period);
    reflect = bound_info_reflect(bi);
    CuAssertIntEquals(tc,0,reflect);
    N = mcnode_get_n(mcn);
    CuAssertIntEquals(tc,4,N);
    p = mcnode_get_pself(mcn);
//    printf("p = %G\n",p);
    CuAssertIntEquals(tc,1,p>=0);
    /* CuAssertDblEquals(tc,1.0,p,1e-15); */
    mcnode_free(mcn);
    bound_info_free(bi);

    mca_free_deep(mca);
}


CuSuite * MCAGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_mcnode_alloc);
    SUITE_ADD_TEST(suite, Test_mcnode_init);
    SUITE_ADD_TEST(suite, Test_mcnode_prepend);
    SUITE_ADD_TEST(suite, Test_mcnode_expectation);
    SUITE_ADD_TEST(suite, Test_mca_get_node_absorb);
    SUITE_ADD_TEST(suite, Test_mca_get_node_absorb_grad);
    SUITE_ADD_TEST(suite, Test_mca_get_node_reflect);
    return suite;
}

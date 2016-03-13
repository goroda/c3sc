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


static const size_t dglob = 5;

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
    double x[5] = {1.0,2.0,3.0,4.0,5.0};
    mcn = mcnode_init(dglob,x);

    size_t dg = mcnode_get_d(mcn);
    CuAssertIntEquals(tc,dglob,dg);

    double * xt = mcnode_get_xref(mcn);
    double diff = norm2diff(x,xt,5);
    CuAssertDblEquals(tc,0.0,diff,1e-15);

    double pself = mcnode_get_pself(mcn);
    CuAssertDblEquals(tc,1.0,pself,1e-15);

    size_t N = mcnode_get_N(mcn);
    CuAssertIntEquals(tc,0,N);

    struct MCNode ** neigh = mcnode_get_neighbors(mcn);
    CuAssertIntEquals(tc,1,neigh==NULL);

    double * p = mcnode_get_pref(mcn);
    CuAssertIntEquals(tc,1,p==NULL);

    mcnode_free(mcn); mcn = NULL;
}

void Test_mcnode_hspace(CuTest * tc)
{
    printf("Testing Function: mcnode_add_neighbors_hspace (with no BondaryInfo)\n");


    struct MCNode * mcn = NULL;
    double x[5] = {1.0,2.0,3.0,4.0,5.0};
    mcn = mcnode_init(dglob,x);

    double h[5] = {1.0,0.5,0.25,0.25,0.2};
    double pselfin = 0.2;
    double probs[10] = {0.2,0.2,0.2,0.1,0.1,
                        0.0,0.0,0.0,0.0,0.0};

    mcnode_add_neighbors_hspace(mcn,h,pselfin,probs,NULL);

    double pself = mcnode_get_pself(mcn);
    CuAssertDblEquals(tc,pselfin,pself,1e-15);

    size_t N = mcnode_get_N(mcn);
    CuAssertIntEquals(tc,2*dglob,N);

    struct MCNode ** neigh = mcnode_get_neighbors(mcn);
    for (size_t ii = 0; ii < dglob; ii++){
        double * y = mcnode_get_xref(neigh[2*ii]);
        double * z = mcnode_get_xref(neigh[2*ii+1]);
        for (size_t jj = 0; jj < dglob;jj++){
            if (jj == ii){
                CuAssertDblEquals(tc,x[ii]-h[ii],y[ii],1e-15);
                CuAssertDblEquals(tc,x[ii]+h[ii],z[ii],1e-15);
            }
            else{
                CuAssertDblEquals(tc,x[jj],y[jj],1e-15);
                CuAssertDblEquals(tc,x[jj],z[jj],1e-15);
            }
        }
    }
    double * p = mcnode_get_pref(mcn);
    double diff = norm2diff(p,probs,2*dglob);
    CuAssertDblEquals(tc,0.0,diff,1e-15);

    mcnode_free(mcn); mcn = NULL;
}

void Test_mcnode_hspace2(CuTest * tc)
{
    printf("Testing Function: mcnode_add_neighbors_hspace (with absorb BondaryInfo)\n");


    struct MCNode * mcn = NULL;
    double x[5] = {1.0,2.0,3.0,4.0,5.0};
    mcn = mcnode_init(dglob,x);

    double h[5] = {1.0,0.5,0.25,0.25,0.2};
    double pselfin = 0.2;
    double probs[10] = {0.2,0.2,0.2,0.1,0.1,
                        0.0,0.0,0.0,0.0,0.0};

    struct BoundInfo * bi = bound_info_alloc(dglob);
    int res;
    for (size_t ii = 0; ii < dglob;ii++){
        res = bound_info_set_dim(bi,IN,NONE,ii);
        CuAssertIntEquals(tc,0,res);
    }

    double mapval = 0.12345;
    size_t mapind = 2;
    res = bound_info_set_dim(bi,LEFT,PERIODIC,mapind);
    CuAssertIntEquals(tc,1,res);

    bound_info_set_xmap_dim(bi,mapval,mapind);
    
    mcnode_add_neighbors_hspace(mcn,h,pselfin,probs,bi);

    double pself = mcnode_get_pself(mcn);
    CuAssertDblEquals(tc,pselfin,pself,1e-15);

    size_t N = mcnode_get_N(mcn);
    CuAssertIntEquals(tc,2*dglob,N);

    //mcnode_print(mcn,stdout,3);
    struct MCNode ** neigh = mcnode_get_neighbors(mcn);
    for (size_t ii = 0; ii < dglob; ii++){
        double * y = mcnode_get_xref(neigh[2*ii]);
        double * z = mcnode_get_xref(neigh[2*ii+1]);
        for (size_t jj = 0; jj < dglob;jj++){
            if (jj == ii){
                if (ii != mapind){
                    CuAssertDblEquals(tc,x[ii]-h[ii],y[ii],1e-15);
                    CuAssertDblEquals(tc,x[ii]+h[ii],z[ii],1e-15);
                }
                else{
                    CuAssertDblEquals(tc,mapval-h[ii],y[ii],1e-15);
                    CuAssertDblEquals(tc,x[ii]+h[ii],z[ii],1e-15);
                }
            }
            else{
                CuAssertDblEquals(tc,x[jj],y[jj],1e-15);
                CuAssertDblEquals(tc,x[jj],z[jj],1e-15);
            }
        }
    }
    double * p = mcnode_get_pref(mcn);
    double diff = norm2diff(p,probs,2*dglob);
    CuAssertDblEquals(tc,0.0,diff,1e-15);

    bound_info_free(bi); bi = NULL;
    mcnode_free(mcn); mcn = NULL;
}

void Test_mcnode_gradients(CuTest * tc)
{
    printf("Testing Function: mcnode_add_gradients\n");

    struct MCNode * mcn = NULL;
    double x[5] = {1.0,2.0,3.0,4.0,5.0};
    mcn = mcnode_init(dglob,x);

    double h[5] = {1.0,0.5,0.25,0.25,0.2};
    double pselfin = 0.2;
    double probs[10] = {0.2,0.2,0.2,0.1,0.1,
                        0.0,0.0,0.0,0.0,0.0};

    size_t du = 4;
    double ** gprob = malloc_dd(2*dglob);
    for (size_t ii = 0; ii < 2*dglob; ii++){
        gprob[ii] = calloc_double(du);
        for (size_t jj = 0; jj < du; jj++){
            gprob[ii][jj] = randu();
        }
    }
    double * gprobself = calloc_double(du);
    for (size_t jj = 0; jj < du; jj++){
        gprobself[jj] = randu();
    }
    

    mcnode_add_neighbors_hspace(mcn,h,pselfin,probs,NULL);
    mcnode_add_gradients(mcn,2*dglob,du,gprobself,gprob);

    size_t ddu = mcnode_get_du(mcn);
    CuAssertIntEquals(tc,du,ddu);
    double * gpself = mcnode_get_gpself(mcn);
    CuAssertIntEquals(tc,1,gpself!=NULL);

    double ** gpn = mcnode_get_gp(mcn);
    CuAssertIntEquals(tc,1,gpn!=NULL);

    CuAssertIntEquals(tc,1,0);
    mcnode_free(mcn); mcn = NULL;
}

void expfunc(size_t Nvals, double * times, double ** states, double * out, void * arg)
{
    (void)(times);
    (void)(arg);

    for (size_t ii = 0; ii < Nvals; ii++){
        double * x = states[ii];
        out[ii] = x[1]*x[2] + sin(x[3]);
    }
}

void Test_mcnode_expectation(CuTest * tc)
{
    printf("Testing Function: mcnode_expectation\n");

    struct MCNode * mcn = NULL;
    double x[5] = {1.0,2.0,3.0,4.0,5.0};
    mcn = mcnode_init(dglob,x);

    double h[5] = {1.0,0.5,0.25,0.25,0.2};
    double pselfin = 0.2;
    double probs[10] = {0.2,0.2,0.2,0.1,0.1,
                        0.0,0.0,0.0,0.0,0.0};

    size_t du = 4;
    double ** gprob = malloc_dd(2*dglob);
    for (size_t ii = 0; ii < 2*dglob; ii++){
        gprob[ii] = calloc_double(du);
        for (size_t jj = 0; jj < du; jj++){
            gprob[ii][jj] = randu();
        }
    }
    double * gprobself = calloc_double(du);
    for (size_t jj = 0; jj < du; jj++){
        gprobself[jj] = randu();
    }
    
    mcnode_add_neighbors_hspace(mcn,h,pselfin,probs,NULL);
    mcnode_add_gradients(mcn,2*dglob,du,gprobself,gprob);
    
    double avg = mcnode_expectation(mcn,expfunc,NULL,NULL);
    double as = 0.0;
//    as += 
    
    CuAssertIntEquals(tc,1,0);
    mcnode_free(mcn); mcn = NULL;
}


CuSuite * MCAGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_mcnode_alloc);
    SUITE_ADD_TEST(suite, Test_mcnode_init);
    SUITE_ADD_TEST(suite, Test_mcnode_hspace);
    SUITE_ADD_TEST(suite, Test_mcnode_hspace2);
    SUITE_ADD_TEST(suite, Test_mcnode_gradients);
    return suite;
}

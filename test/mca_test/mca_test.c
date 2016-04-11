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

    struct MCNList * neigh = mcnode_get_neigh(mcn);
    for (size_t ii = 0; ii < dglob;ii++){
        mcnlist_prepend(&neigh,ii,x[ii]-h[ii],probs[ii],NULL);
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
    }
    
    mcnode_free(mcn); mcn = NULL;
}

double expfunc(double t, double * x, void * args)
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

CuSuite * MCAGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_mcnode_alloc);
    SUITE_ADD_TEST(suite, Test_mcnode_init);
    SUITE_ADD_TEST(suite, Test_mcnode_prepend);
    SUITE_ADD_TEST(suite, Test_mcnode_expectation);
    return suite;
}

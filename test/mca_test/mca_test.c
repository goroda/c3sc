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


size_t dglob = 5;

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

    // struct MCNode ** neigh = mcnode_get_neighbors(mcn);
    //CuAssertIntEquals(tc,1,neigh==NULL);

    double * p = mcnode_get_pref(mcn);
    double diff = norm2diff(p,probs,2*dglob);
    CuAssertDblEquals(tc,0.0,diff,1e-15);

    mcnode_free(mcn); mcn = NULL;
}


CuSuite * MCAGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_mcnode_alloc);
    SUITE_ADD_TEST(suite, Test_mcnode_init);
    SUITE_ADD_TEST(suite, Test_mcnode_hspace);
    return suite;
}

// Copyright (c) 2016, Massachusetts Institute of Technology
// Authors: Alex Gorodetsky
// Email  : goroda@mit.edu

//Code

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CuTest.h"	
#include "util.h"


void Test_hashgrid(CuTest * tc)
{
    printf("Testing Function: hashgrid \n");

    size_t size = 10000;
    
    size_t N = 20;
    double * x = linspace(-10.0,10.0,N);

    struct HashGrid * hg = hash_grid_create(size);
    
    for (size_t ii = 0; ii < N; ii++){
        printf("hashing ii=%zu, x=%G\n",ii,x[ii]);
        int exists = hash_grid_add_element(hg,ii,x[ii]);

        /* int exp; */
        /* frexp(x[ii],&exp); */
        /* printf("exp should be %d\n",exp); */
        CuAssertIntEquals(tc,0,exists);
    }

    for (size_t ii = 0; ii < N; ii++){
        size_t ind = hash_grid_get_ind(hg,x[ii]);
        printf("ii=%zu, ind=%zu\n",ii,ind);
//        CuAssertIntEquals(tc,ii,ind);
    }

    hash_grid_free(hg);
    free(x);
    CuAssertIntEquals(tc,0,0);
}

CuSuite * HashGridGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_hashgrid);
    return suite;
}

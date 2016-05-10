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
    
    size_t N = 100000;
    double * x = linspace(-10.0,10.0,N);

    struct HashGrid * hg = hash_grid_create(size);
    
    for (size_t ii = 0; ii < N; ii++){
        /* printf("hashing ii=%zu, x=%G\n",ii,x[ii]); */
        int exists = hash_grid_add_element(hg,ii,x[ii]);

        /* int exp; */
        /* frexp(x[ii],&exp); */
        /* printf("exp should be %d\n",exp); */
        CuAssertIntEquals(tc,0,exists);
    }

    for (size_t ii = 0; ii < N; ii++){
        int exists;
        size_t ind = hash_grid_get_ind(hg,x[ii],&exists);
        CuAssertIntEquals(tc,1,exists);
        CuAssertIntEquals(tc,ii,ind);
    }

    hash_grid_free(hg);
    free(x);
}

void Test_hashgrid_grid(CuTest * tc)
{
    printf("Testing Function: hashgrid_grid \n");

    size_t size = 10000;
    
    size_t N = 100000;
    double * x = linspace(-10.0,10.0,N);
    struct c3Vector * c3v = c3vector_alloc(N,x);
    struct HashGrid * hg = hash_grid_create_grid(size,c3v);

    for (size_t ii = 0; ii < N; ii++){
        int exists;
        size_t ind = hash_grid_get_ind(hg,x[ii],&exists);
        CuAssertIntEquals(tc,1,exists);
        CuAssertIntEquals(tc,ii,ind);
    }

    hash_grid_free(hg);
    free(x);
    c3vector_free(c3v);
}

void Test_hashgrid_ndgrid(CuTest * tc)
{
    printf("Testing Function: hashgrid_ndgrid \n");

    size_t size = 10000;
    size_t dim = 3;
    size_t N[3] = {1000,2000,3000};
    

    struct c3Vector ** grid = c3vector_array_alloc(dim);
    for (size_t ii = 0; ii < dim; ii++){
        double * x = linspace(-10.0,10.0,N[ii]);
        grid[ii] = c3vector_alloc(N[ii],x);
        free(x); x = NULL;
    }
    struct HashGrid ** hg = hash_grid_create_ndgrid(size,3,grid);

    for (size_t jj = 0; jj < dim; jj++){
        double * x = linspace(-10.0,10.0,N[jj]);
        for (size_t ii = 0; ii < N[jj]; ii++){
            int exists;
            size_t ind = hash_grid_get_ind(hg[jj],x[ii],&exists);
            CuAssertIntEquals(tc,1,exists);
            CuAssertIntEquals(tc,ii,ind);
        }
        free(x); x = NULL;
    }

    hash_grid_free_ndgrid(dim,hg);
    c3vector_array_free(dim,grid);
}

CuSuite * HashGridGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_hashgrid);
    SUITE_ADD_TEST(suite, Test_hashgrid_grid);
    SUITE_ADD_TEST(suite, Test_hashgrid_ndgrid);
    return suite;
}

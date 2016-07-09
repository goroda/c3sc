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

#include "hashgrid.h"

/* void Test_hashgrid(CuTest * tc) */
/* { */
/*     printf("Testing Function: hashgrid \n"); */

/*     size_t size = 10000; */
    
/*     size_t N = 100000; */
/*     double * x = linspace(-10.0,10.0,N); */

/*     struct HashGrid * hg = hash_grid_create(size); */
    
/*     for (size_t ii = 0; ii < N; ii++){ */
/*         /\* printf("hashing ii=%zu, x=%G\n",ii,x[ii]); *\/ */
/*         int exists = hash_grid_add_element(hg,ii,x[ii]); */

/*         /\* int exp; *\/ */
/*         /\* frexp(x[ii],&exp); *\/ */
/*         /\* printf("exp should be %d\n",exp); *\/ */
/*         CuAssertIntEquals(tc,0,exists); */
/*     } */

/*     for (size_t ii = 0; ii < N; ii++){ */
/*         int exists; */
/*         size_t ind = hash_grid_get_ind(hg,x[ii],&exists); */
/*         CuAssertIntEquals(tc,1,exists); */
/*         CuAssertIntEquals(tc,ii,ind); */
/*     } */

/*     hash_grid_free(hg); */
/*     free(x); */
/* } */

/* void Test_hashgrid_grid(CuTest * tc) */
/* { */
/*     printf("Testing Function: hashgrid_grid \n"); */

/*     size_t size = 10000; */
    
/*     size_t N = 100000; */
/*     double * x = linspace(-10.0,10.0,N); */
/*     struct c3Vector * c3v = c3vector_alloc(N,x); */
/*     struct HashGrid * hg = hash_grid_create_grid(size,c3v); */

/*     for (size_t ii = 0; ii < N; ii++){ */
/*         int exists; */
/*         size_t ind = hash_grid_get_ind(hg,x[ii],&exists); */
/*         CuAssertIntEquals(tc,1,exists); */
/*         CuAssertIntEquals(tc,ii,ind); */
/*     } */

/*     hash_grid_free(hg); */
/*     free(x); */
/*     c3vector_free(c3v); */
/* } */

/* void Test_hashgrid_ndgrid(CuTest * tc) */
/* { */
/*     printf("Testing Function: hashgrid_ndgrid \n"); */

/*     size_t size = 10000; */
/*     size_t dim = 3; */
/*     size_t N[3] = {1000,2000,3000}; */
    

/*     struct c3Vector ** grid = c3vector_array_alloc(dim); */
/*     for (size_t ii = 0; ii < dim; ii++){ */
/*         double * x = linspace(-10.0,10.0,N[ii]); */
/*         grid[ii] = c3vector_alloc(N[ii],x); */
/*         free(x); x = NULL; */
/*     } */
/*     struct HashGrid ** hg = hash_grid_create_ndgrid(size,3,grid); */

/*     for (size_t jj = 0; jj < dim; jj++){ */
/*         double * x = linspace(-10.0,10.0,N[jj]); */
/*         for (size_t ii = 0; ii < N[jj]; ii++){ */
/*             int exists; */
/*             size_t ind = hash_grid_get_ind(hg[jj],x[ii],&exists); */
/*             CuAssertIntEquals(tc,1,exists); */
/*             CuAssertIntEquals(tc,ii,ind); */
/*         } */
/*         free(x); x = NULL; */
/*     } */

/*     hash_grid_free_ndgrid(dim,hg); */
/*     c3vector_array_free(dim,grid); */
/* } */

CuSuite * HashGridGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    /* SUITE_ADD_TEST(suite, Test_hashgrid); */
    /* SUITE_ADD_TEST(suite, Test_hashgrid_grid); */
    /* SUITE_ADD_TEST(suite, Test_hashgrid_ndgrid); */
    return suite;
}


void Test_htable(CuTest * tc)
{
    printf("Testing htable\n");

    size_t size = 100000;
    struct HTable * hg = htable_create(size);
    CuAssertIntEquals(tc,0,hg==NULL);
    
    size_t nvals = 600;
    size_t sval = 10;
    double ** x = malloc(nvals * sizeof(double *));
    for (size_t ii = 0; ii < nvals; ii++){
        /* printf("ii = %zu\n",ii); */
        size_t val = ii+1;
        char * blah = size_t_to_char(val);
        x[ii] = calloc(sval,sizeof(double));
        for (size_t jj = 0; jj < sval; jj++){
            x[ii][jj] = randu()*10.0 - 5.0;
        }
        
        /* printf("key=%s\n",blah); */
        int res = htable_add_element(hg,blah,x[ii],sval*sizeof(double));
        CuAssertIntEquals(tc,0,res);
        /* free(blah); blah = NULL; */
    }

    /* printf("\n\n\n\n"); */
    for (size_t ii = 0; ii < nvals; ii++){
        size_t val = ii+1;
        char * blah = size_t_to_char(val);
        double * data_get = NULL;
        size_t nbytes;

        data_get = htable_get_element(hg,blah,&nbytes);
        CuAssertIntEquals(tc,sval*sizeof(double),nbytes);
        CuAssertIntEquals(tc,0,data_get==NULL);

        for (size_t jj = 0; jj < sval; jj++){
            /* printf("data[%zu]=(%G,%G)\n ",jj,data_get[jj],x[ii][jj]); */
            CuAssertDblEquals(tc,x[ii][jj],data_get[jj],1e-15);
        }
        free(blah);
    }

    for (size_t ii = 0; ii < nvals; ii++){
        free(x[ii]); x[ii] = NULL;
    }
    free(x); x = NULL;
    
    htable_destroy(hg); hg = NULL;

}

void Test_htable2(CuTest * tc)
{
    printf("Testing htable2 - with size_t array keys \n");

    size_t size = 10000;
    struct HTable * hg = htable_create(size);
    CuAssertIntEquals(tc,0,hg==NULL);
    
    size_t nvals = 8000;
    size_t sval = 10;
    double ** x = malloc(nvals * sizeof(double *));
    for (size_t ii = 0; ii < nvals; ii++){
        /* printf("ii = %zu\n",ii); */
        size_t val[5] = {ii+1,ii+2,ii+3,ii+4,ii+5};
        char * blah = size_t_a_to_char(val,5);

        
        x[ii] = calloc(sval,sizeof(double));
        for (size_t jj = 0; jj < sval; jj++){
            x[ii][jj] = randu()*10.0 - 5.0;
        }
        
        /* printf("key=%s\n",blah); */
        int res = htable_add_element(hg,blah,x[ii],sval*sizeof(double));
        CuAssertIntEquals(tc,0,res);
        /* free(blah); blah = NULL; */
    }

    /* printf("\n\n\n\n"); */
    for (size_t ii = 0; ii < nvals; ii++){

        
        size_t val[5] = {ii+1,ii+2,ii+3,ii+4,ii+5};
        char * blah = size_t_a_to_char(val,5);
        double * data_get = NULL;
        size_t nbytes;

        data_get = htable_get_element(hg,blah,&nbytes);
        CuAssertIntEquals(tc,sval*sizeof(double),nbytes);
        CuAssertIntEquals(tc,0,data_get==NULL);

        for (size_t jj = 0; jj < sval; jj++){
            /* printf("data[%zu]=(%G,%G)\n ",jj,data_get[jj],x[ii][jj]); */
            CuAssertDblEquals(tc,x[ii][jj],data_get[jj],1e-15);
        }
        free(blah);
    }

    for (size_t ii = 0; ii < nvals; ii++){
        free(x[ii]); x[ii] = NULL;
    }
    free(x); x = NULL;
    
    htable_destroy(hg); hg = NULL;

}

CuSuite * HTableGetSuite()
{
    //printf("----------------------------\n");

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_htable);
    SUITE_ADD_TEST(suite, Test_htable2);

    return suite;
}

// Copyright (C) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

//Code

#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"

CuSuite * HashGridGetSuite();


void RunAllTests(void) {
    
    printf("Running test suite for: lib_linalg\n");
    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * mca = HashGridGetSuite();

    CuSuiteAddSuite(suite, mca);

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);
    
    CuSuiteDelete(mca);
    
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}

// Copyright (C) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

//Code

#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"

CuSuite * TProbGetSuite();
CuSuite * ValueFGetSuite();

void RunAllTests(void) {
    
    printf("Running test suite for: lib_linalg\n");
    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * mca = TProbGetSuite();
    CuSuite * val = ValueFGetSuite();

    /* CuSuiteAddSuite(suite, mca); */
    CuSuiteAddSuite(suite, val);

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);
    
    CuSuiteDelete(mca);
    CuSuiteDelete(val);
    
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}

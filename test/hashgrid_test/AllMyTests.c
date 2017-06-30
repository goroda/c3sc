// Copyright (C) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

//Code

#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"

CuSuite * HashGridGetSuite(void);
CuSuite * HTableGetSuite(void);


void RunAllTests(void) {
    
    printf("Running test suite for: hashing\n");
    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * mca = HashGridGetSuite();
    CuSuite * tab = HTableGetSuite();

    CuSuiteAddSuite(suite, mca);
    CuSuiteAddSuite(suite, tab);

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);
    
    CuSuiteDelete(mca);
    CuSuiteDelete(tab);
    
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}

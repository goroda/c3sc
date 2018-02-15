// Copyright (C) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

//Code

#include <stdio.h>
#include <stdlib.h>
#include "CuTest.h"

CuSuite * TProbGetSuite(void);
CuSuite * ValueFGetSuite(void);
CuSuite * BellmanGetSuite(void);
CuSuite * DPAlgsGetSuite(void);

void RunAllTests(void) {
    
    printf("Running test suite\n");
    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * mca = TProbGetSuite();
    CuSuite * val = ValueFGetSuite();
    CuSuite * bel = BellmanGetSuite();
    CuSuite * dp = DPAlgsGetSuite();

    /* CuSuiteAddSuite(suite, mca); */
    /* CuSuiteAddSuite(suite, val); */
    /* CuSuiteAddSuite(suite, bel); */
    CuSuiteAddSuite(suite, dp);

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);
    
    CuSuiteDelete(mca);
    CuSuiteDelete(val);
    CuSuiteDelete(bel);
    CuSuiteDelete(dp);
    
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}

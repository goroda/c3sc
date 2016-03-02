#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "c3.h"

#include "boundary.h"



/** \struct Boundary
 *  \brief Structure to handle any boundary conditions
 *  \var Boundary::d
 *  dimension of state space
 *  \var Boundary::dirs
 *  (d) sized array with locations on boundaries (-1->left,1->right)
 *  \var Boundary::bcheck
 *  user supplied function for checking
 *  \var Boundary::args
 *  arguments to bcheck
 */
struct Boundary
{
    size_t d;
    int * dirs;
    int clean;
    
    int (*bcheck)(double,double *, void *, int *);
    void * args;

};

/**********************************************************//**
    Allocate Boundary

    \param[in] d    - dimension of state space
    \param[in] b    - function(time,x,args,dirs)
    \param[in] arg  - function arguments

    \return boundary
**************************************************************/
struct Boundary *
boundary_alloc(size_t d, int (*b)(double, double *,void *,int*), void *arg)
{
    struct Boundary * bound;
    bound = malloc(sizeof(struct Boundary));
    assert (bound != NULL);

    bound->d = d;
    bound->dirs = calloc_int(d);
    bound->clean = 1;
    
    bound->bcheck = b;
    bound->args = arg;
    return bound;
}

void boundary_free(struct Boundary * bound)
{
    if (bound != NULL){
        free(bound->dirs); bound->dirs = NULL;
        free(bound); bound = NULL;
    }
}

int boundary_type(struct Boundary * bound, double time, double * x)
{
    bound->clean = 1;
    int outbounds = bound->bcheck(time,x,bound->args,bound->dirs);
    bound->clean = 0;
    return outbounds;
}


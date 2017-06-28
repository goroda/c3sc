// Copyright (c) 2015-2016, Massachusetts Institute of Technology
//
// This file is part of the C3 for Stochastic Optimal Control (C3SC) toolbox
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice, 
//    this list of conditions and the following disclaimer in the documentation 
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its contributors 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//Code

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* #include "c3.h" */
#include "c3/array.h"
#include "boundary.h"

/** \struct ExternalBoundary
 *  \brief Structure to handle external boundary conditions
  *  \var ExternalBoundary::left
 *  left boundary
 *  \var ExternalBoundary::right
 *  right boundary
 *  \var ExternalBoundary::type
 *  type of boundary
 */
struct ExternalBoundary
{
    double left;
    double right;
    enum EBTYPE type;
};

/**********************************************************//**
    Allocate External boundary

    \param[in] left  - left edge
    \param[in] right - right edge
    \param[in] type  - "absorb","periodic","reflect"

    \return boundary
**************************************************************/
struct ExternalBoundary * 
external_boundary_alloc(double left, double right, char * type)
{
    struct ExternalBoundary * db = malloc(sizeof(struct ExternalBoundary));
    if (db == NULL){
        fprintf(stderr,"Mem error allocating ExternalBoundary\n");
        exit(1);
    }
    db->left = left;
    db->right = right;
    if (strcmp(type,"absorb") == 0){
        db->type = ABSORB;
    }
    else if (strcmp(type,"periodic") == 0){
        db->type = PERIODIC;
    }
    else if (strcmp(type,"reflect") == 0){
        db->type = REFLECT;
    }
    else{
        fprintf(stderr, "External boundary of type %s is unknown\n",type);
        exit(1);
    }
    return db;
}

/**********************************************************//**
    Copy External boundary
**************************************************************/
struct ExternalBoundary *
external_boundary_copy(struct ExternalBoundary * old)
{
    if (old == NULL){
        return NULL;
    }
    struct ExternalBoundary * db = malloc(sizeof(struct ExternalBoundary));
    if (db == NULL){
        fprintf(stderr,"Mem error allocating ExternalBoundary\n");
        exit(1);
    }
    db->left = old->left;
    db->right = old->right;
    db->type = old->type;

    return db;
}

/**********************************************************//**
    Free External boundary
**************************************************************/
void external_boundary_free(struct ExternalBoundary * db)
{
    if (db != NULL){
        free(db); db = NULL;
    }
}

/**********************************************************//**
    Allocate an array of external boundaries

    \param[in] N - number of arrays

    \return array with each element set to NULL
**************************************************************/
struct ExternalBoundary ** 
external_boundary_alloc_array(size_t N)
{
    struct ExternalBoundary ** db = malloc(N*sizeof(struct ExternalBoundary *));
    if (db == NULL){
        fprintf(stderr,"Mem error allocating ExternalBoundary array\n");
        exit(1);
    }
    for (size_t ii = 0; ii < N; ii++){
        db[ii] = NULL;
    }
    return db;
}


/**********************************************************//**
    Free external boundary array
**************************************************************/
void external_boundary_free_array(struct ExternalBoundary ** db, size_t N)
{
    if (db != NULL){
        for (size_t ii = 0; ii < N; ii++){
            external_boundary_free(db[ii]); db[ii] = NULL;
        }
        free(db); db = NULL;
    }
}

/**********************************************************//**
    Set type of external boundary
**************************************************************/
void external_boundary_set_type(struct ExternalBoundary * db, char * type)
{
    if (strcmp(type,"absorb") == 0){
        db->type = ABSORB;
    }
    else if (strcmp(type,"periodic") == 0){
        db->type = PERIODIC;
    }
    else if (strcmp(type,"reflect") == 0){
        db->type = REFLECT;
    }
    else{
        fprintf(stderr, "External boundary of type %s is unknown\n",type);
        exit(1);
    }
}

/**********************************************************//**
    Get type of external boundary
**************************************************************/
enum EBTYPE external_boundary_get_type(struct ExternalBoundary *bd)
{
    assert (bd != NULL);
    return bd->type;
}

/**********************************************************//**
    Get left boundary
**************************************************************/
double external_boundary_get_left(struct ExternalBoundary *bd)
{
    assert (bd != NULL);
    return bd->left;
}

/**********************************************************//**
    Get right boundary
**************************************************************/
double external_boundary_get_right(struct ExternalBoundary *bd)
{
    assert (bd != NULL);
    return bd->right;
}

/**********************************************************//**
    Check if point violates left boundary 

    \return 1 - violates
            0 - doesnt
**************************************************************/
int external_boundary_check_left(struct ExternalBoundary * bd, double x)
{
    assert (bd != NULL);
    if (x <= bd->left){
        return 1;
    }
    return 0;
}

/**********************************************************//**
    Check if point violates right

    \return 1 - violates
            0 - doesnt
**************************************************************/
int external_boundary_check_right(struct ExternalBoundary * bd, double x)
{
    assert (bd != NULL);
    if (x >= bd->right){
        return 1;
    }
    return 0;
}

struct BoundRect
{
    size_t dim;
    double * center;
    double * lb;
    double * ub;
};

struct BoundRect * bound_rect_init(size_t dim, double * center, double * lengths)
{
    struct BoundRect * br;
    br = malloc(sizeof(struct BoundRect));
    if (br == NULL){
        fprintf(stderr, "Error allocating BoundRect\n");
        exit(1);
    }
    br->dim = dim;
    br->center = calloc_double(dim);
    memmove(br->center,center,dim*sizeof(double));
    br->lb = calloc_double(dim);
    br->ub = calloc_double(dim);
    for (size_t ii = 0; ii < dim; ii++){
        br->lb[ii] = center[ii] - lengths[ii]/2.0;
        br->ub[ii] = center[ii] + lengths[ii]/2.0;
    }
    return br;
}

struct BoundRect * bound_rect_copy(struct BoundRect * old)
{
    if (old == NULL){
        return NULL;
    }

    struct BoundRect * br;
    br = malloc(sizeof(struct BoundRect));
    if (br == NULL){
        fprintf(stderr, "Error allocating BoundRect\n");
        exit(1);
    }
    br->dim = old->dim;
    br->center = calloc_double(old->dim);
    br->lb = calloc_double(old->dim);
    br->ub = calloc_double(old->dim);
    memmove(br->center,old->center,old->dim*sizeof(double));
    memmove(br->lb,old->lb,old->dim*sizeof(double));
    memmove(br->ub,old->ub,old->dim*sizeof(double));

    return br;
}

struct BoundRect ** bound_rect_alloc_array(size_t N)
{
    struct BoundRect ** br;
    br = malloc(N * sizeof(struct BoundRect *));
    if (br == NULL){
        fprintf(stderr,"Error allocating an array of BoundRect\n");
        exit(1);
    }
    for (size_t ii = 0; ii < N; ii++)
    {
        br[ii] = NULL;
    }
    return br;
}

void bound_rect_free(struct BoundRect * br)
{
    if (br != NULL){
        free(br->center); br->center = NULL;
        free(br->lb); br->lb = NULL;
        free(br->ub); br->ub = NULL;
        free(br); br = NULL;
    }
}

void bound_rect_free_array(struct BoundRect ** br, size_t N)
{
    if (br != NULL){
        for (size_t ii = 0; ii < N; ii++){
            bound_rect_free(br[ii]); br[ii] = NULL;
        }
        free(br); br = NULL;
    }
}

int bound_rect_inside(struct BoundRect * br,const double * x)
{
    int inside_obstacle = 1;
    /* printf("br dim = %zu\n",br->dim); */
    /* printf("x in bound = "); dprint(br->dim,x); */
    assert (br->dim > 0);
    for (size_t ii = 0; ii < br->dim; ii++){
        /* printf("ii=%zu, lb=%G,ub=%G\n",ii,br->lb[ii],br->ub[ii]); */

        if ((x[ii] < br->lb[ii]) || (x[ii] > br->ub[ii])){
            inside_obstacle = 0;
            break;
        }
    }
    return inside_obstacle;
}

/** \struct Boundary
 *  \brief Structure to handle any boundary conditions
 *  \var Boundary::d
 *  dimension of state space
 *  \var Boundary::eb
 *  external boundaries
 */
struct Boundary
{
    size_t d;

    struct ExternalBoundary ** eb;

    size_t n;
    size_t nalloc;
    struct BoundRect ** br;
    double * costs;
};

/**********************************************************//**
    Allocate Boundary and initialize each dimension to absorbing

    \param[in] d  - dimension of state space
    \param[in] lb - lower bounds of state space
    \param[in] ub - upper bounds of state space

    \return boundary
**************************************************************/
struct Boundary *
boundary_alloc(size_t d,double * lb, double * ub)
{
    struct Boundary * bound;
    bound = malloc(sizeof(struct Boundary));
    assert (bound != NULL);

    bound->d = d;    
    bound->eb = external_boundary_alloc_array(d);
    /* printf("bounds = \n");  */
    /* dprint(d,lb); */
    /* dprint(d,ub); */
    for (size_t ii = 0; ii < d; ii++){
        bound->eb[ii] = external_boundary_alloc(lb[ii],ub[ii],"absorb");
        /* printf("%G,%G\n",external_boundary_get_left(bound->eb[ii]), */
        /*        external_boundary_get_right(bound->eb[ii])); */
    }

    bound->n = 0;
    bound->nalloc = 10;
    bound->br = bound_rect_alloc_array(bound->nalloc);
    
    return bound;
}

/**********************************************************//**
    Copy a boundary
**************************************************************/
struct Boundary * boundary_copy_deep(struct Boundary * old)
{
    if (old == NULL){
        return old;
    }

    struct Boundary * bound;
    bound = malloc(sizeof(struct Boundary));
    assert (bound != NULL);

    bound->d = old->d;    
    bound->eb = external_boundary_alloc_array(old->d);
    for (size_t ii = 0; ii < old->d; ii++){
        bound->eb[ii] = external_boundary_copy(old->eb[ii]);
    }
    bound->n = old->n;
    bound->nalloc = old->nalloc;
    bound->br = bound_rect_alloc_array(bound->nalloc);
    for (size_t ii = 0; ii < old->n; ii++){
        bound->br[ii] = bound_rect_copy(old->br[ii]);
    }
    return bound;
}

/**********************************************************//**
    Free memory allocated to boundary
**************************************************************/
void boundary_free(struct Boundary * bound)
{
    if (bound != NULL){
        external_boundary_free_array(bound->eb, bound->d);
        bound->eb = NULL;
        bound_rect_free_array(bound->br,bound->nalloc);
        bound->br = NULL;
        free(bound); bound = NULL;
    }
}

/**********************************************************//**
    Get number of obstacles
**************************************************************/
size_t boundary_get_nobs(struct Boundary * bound)
{
    assert (bound != NULL);
    return bound->n;
}

/**********************************************************//**
    Get lower bound of obstacle
**************************************************************/
double * boundary_obstacle_get_lb(struct Boundary * bound, size_t ind)
{
    assert (bound != NULL);
    return bound->br[ind]->lb;
}

/**********************************************************//**
    Get upper bound of obstacle
**************************************************************/
double * boundary_obstacle_get_ub(struct Boundary * bound, size_t ind)
{
    assert (bound != NULL);
    return bound->br[ind]->ub;
}

/**********************************************************//**
    Add a rectangular obstacle
**************************************************************/
void boundary_add_obstacle(struct Boundary * bound, double * center, double * lengths)
{
    if (bound->n == bound->nalloc){
        fprintf(stderr,"Not enough space allocated for obstacles in boundary\n");
        exit(1);
    }
    /* printf("Adding boundary obstacle n = %zu!!\n",bound->n); */
    assert (bound->d > 0);
    bound->br[bound->n] = bound_rect_init(bound->d,center,lengths);
    /* printf("br dim = %zu\n",bound->br[bound->n]->dim); */
    bound->n = bound->n+1;
}

/**********************************************************//**
    Set external boundary type for the specified dimension
**************************************************************/
void boundary_external_set_type(struct Boundary * b,size_t dim,char * type)
{
    external_boundary_set_type(b->eb[dim],type);
}

struct BoundInfo
{
    size_t d;
    enum BOUNDRESULT * br;
    enum EBTYPE * type;
    double * xmap; // equivalence mappigns for periodic boundaries
    int absorb_overall; // if shared edge
    int in_obstacle;
};

/**********************************************************//**
    Allocate boundary info of a certain dimension
**************************************************************/
struct BoundInfo * bound_info_alloc(size_t d)
{
    struct BoundInfo * bi = malloc(sizeof(struct BoundInfo));
    if (bi == NULL){
        fprintf(stderr,"Mem Error allocating boundary info\n");
        exit(1);
    }
    bi->d = d;
    bi->br = malloc(d * sizeof(enum BOUNDRESULT));
    if (bi->br == NULL){
        fprintf(stderr,"Mem Error allocating boundary info\n");
        exit(1);
    }

    bi->xmap = calloc_double(d);
    bi->type = malloc(d * sizeof(enum EBTYPE));
    if (bi->type == NULL){
        fprintf(stderr,"Mem Error allocating boundary info\n");
        exit(1);
    }
    bi->absorb_overall = 0;
    bi->in_obstacle = -1;
    return bi;
}

/**********************************************************//**
    Free boundary info 
**************************************************************/
void bound_info_free(struct BoundInfo * bi)
{
    if (bi != NULL){
        free(bi->br); bi->br = NULL;
        free(bi->type); bi->type = NULL;
        free(bi->xmap); bi->xmap = NULL;
        free(bi); bi = NULL;
    }
}

/**********************************************************//**
    Set a boundary info dimension to a given value

    \return 
            -1 if don't know what to do with type
            0 if successful
            1 if periodic and need further information
**************************************************************/
int bound_info_set_dim(struct BoundInfo * bi, enum BOUNDRESULT br,
                       enum EBTYPE type, size_t dim)
{
    bi->br[dim] = br;
    bi->type[dim] = type;
    if (type == ABSORB){
        bi->absorb_overall = 1;
    }
    else if (type == PERIODIC){
        return 1;
    }
    else if ((type != EB_NONE) & (type != REFLECT)){
        return -1;
    }

    return 0;
}
/**********************************************************//**
    Set a mapping for x
**************************************************************/
int bound_info_set_xmap_dim(struct BoundInfo * bi, double x, size_t dim)

{
    bi->xmap[dim] = x;
    return 0;
}

double outer_bound_dim(const struct Boundary * bound, size_t dim, 
                       double x, int *map)
{
    
    *map = 0;
    if (external_boundary_check_left(bound->eb[dim],x)){
        enum EBTYPE type = external_boundary_get_type(bound->eb[dim]);
        if (type == PERIODIC){
            *map = 1;
            return external_boundary_get_right(bound->eb[dim]);
        }
    }
    else if (external_boundary_check_right(bound->eb[dim],x)){
        enum EBTYPE type = external_boundary_get_type(bound->eb[dim]);
        if (type == PERIODIC){
            *map = 2;
            return external_boundary_get_left(bound->eb[dim]);
        }
    }
    return x;
}

/**********************************************************//**
    Get external boundary type for a particular dimension.
    So far can't have different types on left and right but 
    in the future one might want to
**************************************************************/
enum EBTYPE boundary_type_dim(const struct Boundary * bound, size_t dim, int right)
{
    enum EBTYPE type;
    if (right == 0){ // left boundary
        type = external_boundary_get_type(bound->eb[dim]);
    }
    else{
        type = external_boundary_get_type(bound->eb[dim]);
    }
    return type;
}

/**********************************************************//**
    Return boundary info
**************************************************************/
struct BoundInfo * boundary_type(const struct Boundary * bound,double time,const double * x)
{
    (void)(time);

    struct BoundInfo * bi = bound_info_alloc(bound->d);
    /* if (fabs(x[2])>=3.14159){ */
    /*     printf("bound check above, x= "); */
    /*     dprint(3,x); */
    /*}*/
    enum EBTYPE type;
    int res;
    for (size_t ii = 0; ii < bound->d; ii++){
        bi->br[ii] = IN;
        if (external_boundary_check_left(bound->eb[ii],x[ii])){
            type = external_boundary_get_type(bound->eb[ii]);
            res = bound_info_set_dim(bi,LEFT,type,ii);
            assert (res > -1);
            if (res == 1){
                bound_info_set_xmap_dim(bi,external_boundary_get_right(bound->eb[ii]),ii);
            }
        }
        else if (external_boundary_check_right(bound->eb[ii],x[ii])){
            type = external_boundary_get_type(bound->eb[ii]);
            res = bound_info_set_dim(bi,RIGHT,type,ii);
            assert (res > -1);
            if (res == 1){
                bound_info_set_xmap_dim(bi,external_boundary_get_left(bound->eb[ii]),ii);
            }
        }
    }

    for (size_t ii = 0; ii < bound->n; ii++){

        int inobs = bound_rect_inside(bound->br[ii],x);
//        dprint(6,x);
        if (inobs == 1){
            //          printf("in obstacle! \n");
            bi->in_obstacle = (int) ii;
            bi->absorb_overall = 1;
            break;
        }
    }
    return bi;
}

/**********************************************************//**
    Return 1 if in an obstacle
    Return 0 if not
**************************************************************/
int boundary_in_obstacle(const struct Boundary * bound, const double * x)
{
    /* printf("n obs = %zu\n",bound->n); */
    /* printf("x = "); dprint(bound->d,x); */
    for (size_t ii = 0; ii < bound->n; ii++){
        int inobs = bound_rect_inside(bound->br[ii],x);
        if (inobs == 1){
            /* printf("IN OBSTACLE\n"); */
            return 1;
        }
    }
    return 0;
}

/**********************************************************//**
    Return 0 if not on boundary
**************************************************************/
int bound_info_onbound(const struct BoundInfo * bi)
{
    if (bi->in_obstacle > -1){
        return 1;
    }
    for (size_t ii = 0; ii < bi->d; ii++){
        if (bi->br[ii] != IN){
            return 1;
        }
    }
    return 0;
}

/**********************************************************//**
    Return 0 if dimension dim is not on boundary or 
    1 if in an obstacle or this dimension is on a boundary
**************************************************************/
int bound_info_onbound_dim(const struct BoundInfo * bi,size_t dim)
{
    if (bi->in_obstacle > -1){
        return 1;
    }
    if (bi->br[dim] != IN){
        return 1;
    }
    
    return 0;
}

/**********************************************************//**
    Return 0 if not on boundary
**************************************************************/
int bound_info_absorb(const struct BoundInfo * bi)
{
    if (bi->absorb_overall == 1){
        return 1;
    }
    return 0;
}

/**********************************************************//**
    Return 0 if not a periodic boundary
**************************************************************/
int bound_info_period(const struct BoundInfo * bi)
{
    for (size_t ii = 0; ii < bi->d; ii++){
        if (bi->br[ii] != IN){
            if (bi->type[ii] == PERIODIC){
                return 1;
            }
        }
    }
    return 0;
}

/**********************************************************//**
    Return 0 if periodic boundary, -1 if on left boundary, 1 if on right
**************************************************************/
int bound_info_period_dim_dir(const struct BoundInfo * bi,size_t dim)
{
    if (bi->type[dim] == PERIODIC){
        if (bi->br[dim] == LEFT){
            return -1;
        }
        else{
            return 1;
        }
    }
    return 0;
}

/**********************************************************//**
    Return mapping for periodic boundary conditions
**************************************************************/
double bound_info_period_xmap(const struct BoundInfo * bi, size_t dim)
{
    return bi->xmap[dim];
}

/**********************************************************//**
    Return 0 if not a reflective boundary
**************************************************************/
int bound_info_reflect(const struct BoundInfo * bi)
{
    for (size_t ii = 0; ii < bi->d; ii++){
        if (bi->br[ii] != IN){
            if (bi->type[ii] == REFLECT){
                return 1;
            }
        }
    }
    return 0;
}

/**********************************************************//**
    Return 0 if not reflect boundary, -1 if on left boundary, 1 if on right
**************************************************************/
int bound_info_reflect_dim_dir(const struct BoundInfo * bi,size_t dim)
{
    if (bi->type[dim] == REFLECT){
        if (bi->br[dim] == LEFT){
            return -1;
        }
        else{
            return 1;
        }
    }
    return 0;
}

/**********************************************************//**
    Return the in_bound parameter
**************************************************************/
int bound_info_get_in_obstacle(const struct BoundInfo * bi)
{
    return bi->in_obstacle;
}

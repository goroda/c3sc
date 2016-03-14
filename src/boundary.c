#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "c3.h"

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
    \param[in] type  - "absorb" or "periodic"

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
    else{
        fprintf(stderr, "External boundary of type %s is unknown\n",type);
        exit(1);
    }
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
    for (size_t ii = 0; ii < br->dim; ii++){
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
    Add a rectangular obstacle
**************************************************************/
void boundary_add_obstacle(struct Boundary * bound, double * center, double * lengths)
{
    if (bound->n == bound->nalloc){
        fprintf(stderr,"Not enough space allocated for obstacles in boundary\n");
        exit(1);
    }
    bound->br[bound->n] = bound_rect_init(bound->d,center,lengths);
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
    else if (type != NONE){
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
        if (inobs == 1){
            bi->in_obstacle = (int) ii;
            bi->absorb_overall = 1;
            break;
        }
    }
    return bi;
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
    Return 0 if not on boundary
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
    Return 0 if not on boundary, -1 if on left boundary, 1 if on right
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
    Return the in_bound parameter
**************************************************************/
int bound_info_get_in_obstacle(const struct BoundInfo * bi)
{
    return bi->in_obstacle;
}

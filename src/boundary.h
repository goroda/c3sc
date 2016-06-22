#ifndef C3SC_BOUNDARY_H
#define C3SC_BOUNDARY_H

enum EBTYPE {
    NONE=0,
    ABSORB=1,
    PERIODIC=2,
    REFLECT=3
};

struct Boundary;
struct Boundary * boundary_alloc(size_t,double *,double*);
double * boundary_obstacle_get_lb(struct Boundary *, size_t);
double * boundary_obstacle_get_ub(struct Boundary *, size_t);

size_t boundary_get_nobs(struct Boundary *);
struct Boundary * boundary_copy_deep(struct Boundary *);
void boundary_free(struct Boundary *);
void boundary_external_set_type(struct Boundary *,size_t,char *);
double outer_bound_dim(const struct Boundary *, size_t, 
                       double, int *);
enum EBTYPE boundary_type_dim(const struct Boundary *,size_t, int);
int boundary_in_obstacle(const struct Boundary *, const double *);

struct BoundInfo * boundary_type(const struct Boundary *,double,const double *);

enum BOUNDRESULT {
    IN,
    LEFT,
    RIGHT,
};

struct BoundInfo;
struct BoundInfo * bound_info_alloc(size_t);
void bound_info_free(struct BoundInfo *);
void boundary_add_obstacle(struct Boundary *, double *, double *);
int bound_info_set_dim(struct BoundInfo *, enum BOUNDRESULT,
                       enum EBTYPE, size_t);
int bound_info_set_xmap_dim(struct BoundInfo *, double, size_t);
int bound_info_onbound(const struct BoundInfo *);
int bound_info_onbound_dim(const struct BoundInfo *,size_t);
int bound_info_absorb(const struct BoundInfo *);
int bound_info_period(const struct BoundInfo *);
int bound_info_period_dim_dir(const struct BoundInfo *,size_t);
double bound_info_period_xmap(const struct BoundInfo *,size_t);
int bound_info_reflect(const struct BoundInfo *);
int bound_info_reflect_dim_dir(const struct BoundInfo *,size_t);

int bound_info_get_in_obstacle(const struct BoundInfo *);
#endif

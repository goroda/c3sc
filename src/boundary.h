#ifndef C3SC_BOUNDARY_H
#define C3SC_BOUNDARY_H

enum EBTYPE {
    NONE=0,
    ABSORB=1, 
    PERIODIC=2
};

struct Boundary;
struct Boundary *
boundary_alloc(size_t,double *,double*);
void boundary_free(struct Boundary *);
void boundary_external_set_type(struct Boundary *,size_t,char *);
struct BoundInfo * boundary_type(const struct Boundary *,double,const double *);




enum BOUNDRESULT {
    IN,
    LEFT,
    RIGHT
};
struct BoundInfo;
struct BoundInfo * bound_info_alloc(size_t);
void bound_info_free(struct BoundInfo *);

int bound_info_set_dim(struct BoundInfo *, enum BOUNDRESULT,
                       enum EBTYPE, size_t);
int bound_info_set_xmap_dim(struct BoundInfo *, double, size_t);
int bound_info_onbound(const struct BoundInfo *);
int bound_info_absorb(const struct BoundInfo *);
int bound_info_period(const struct BoundInfo *);
int bound_info_period_dim_dir(const struct BoundInfo *,size_t);
double bound_info_period_xmap(const struct BoundInfo *,size_t);
#endif

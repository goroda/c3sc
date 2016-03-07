#ifndef C3SC_BOUNDARY_H
#define C3SC_BOUNDARY_H

struct Boundary;
struct Boundary *
boundary_alloc(size_t,double *,double*);
void boundary_free(struct Boundary *);
void boundary_external_set_type(struct Boundary *,size_t,char *);
struct BoundInfo * boundary_type(struct Boundary *,double,double *);


struct BoundInfo;
enum BOUNDRESULT {
    IN,
    LEFT,
    RIGHT
};
void bound_info_free(struct BoundInfo *);
int bound_info_onbound(struct BoundInfo *);
int bound_info_absorb(struct BoundInfo *);
int bound_info_period(struct BoundInfo *);
int bound_info_period_dim_dir(struct BoundInfo *,size_t);
double bound_info_period_xmap(struct BoundInfo *,size_t);
#endif

#ifndef DP_H
#define DP_H

#include "c3sc_elements.h"

#include "boundary.h"
#include "cost.h"
#include "tensmarkov.h"
#include "control.h"

enum SCTYPE {IH,FH};
struct C3SC;
typedef struct C3SC* c3sc;

struct C3SC * c3sc_create(enum SCTYPE, size_t, size_t, size_t);
void c3sc_destroy(struct C3SC *);
void c3sc_set_state_bounds(struct C3SC *,double *,double *);
void c3sc_add_dynamics(struct C3SC *,  
                       int (*)(double,double*,double*,
                                double*,double*,void*),
                       void *,
                       int (*)(double,double*,double*,
                                double*,double*,void*),
                       void *);
void c3sc_add_boundary(struct C3SC *,
                       int (*)(double,double *,void *,int*), 
                       void *);
void c3sc_init_mca(struct C3SC *, size_t *);
void c3sc_attach_opt(struct C3SC *, struct c3Opt *);
void c3sc_init_dp(struct C3SC *, double,
                  int (*)(double,double*,double*,double*,double*),
                  int (*)(double,double*,double*));

void * c3sc_get_dp(struct C3SC *);
////////////////////////////////////////////////////
struct DPih;
struct DPih * 
dpih_alloc(double,
           int (*)(double,double*,double*,double*,double*),
           int (*)(double,double*,double*));

void dpih_free(struct DPih *);
void dpih_free_deep(struct DPih *);
void dpih_attach_mca(struct DPih *, struct MCA *);
void dpih_attach_cost(struct DPih *, struct Cost *);
void dpih_attach_policy(struct DPih *, struct Policy *);
void dpih_attach_opt(struct DPih *, struct c3Opt *);

struct Cost * dpih_get_cost(struct DPih *);
struct Dyn * dpih_get_dyn(struct DPih *);

double dpih_rhs(struct DPih *,double *,double *,double *);
struct Cost * dpih_iter_pol(struct DPih *,int);
struct Cost * dpih_iter_vi(struct DPih *,int);
struct Policy * dpih_iter_vi_pol(struct DPih *,int);

struct DPfh;
#endif

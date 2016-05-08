#ifndef DP_H
#define DP_H

#include "util.h"
#include "cost.h"
#include "dynamics.h"
#include "tensmarkov.h"


// prototypes of structures
enum SCTYPE {IH,FH};
struct C3SC;
struct ImplicitPolicy;
struct DPih;
struct C3SCDiagnostic;

typedef struct C3SC* c3sc;

struct C3SC * c3sc_create(enum SCTYPE, size_t, size_t, size_t);
void c3sc_destroy(struct C3SC *);
void c3sc_set_state_bounds(struct C3SC *,double *,double *);
void c3sc_set_external_boundary(struct C3SC *,size_t,char *);
void c3sc_add_obstacle(struct C3SC *, double *, double *);
void c3sc_add_dynamics(struct C3SC *,  
                       int (*)(double,const double*,const double*,
                               double*,double*,void*),
                       void *,
                       int (*)(double,const double*,const double*,
                                double*,double*,void*),
                       void *);
void c3sc_init_mca(struct C3SC *, size_t *);
void c3sc_attach_opt(struct C3SC *, struct c3Opt *);
void c3sc_init_dp(struct C3SC *, double,
                  int (*)(double,const double*,const double*,double*,double*),
                  int (*)(double,const double*,double*),
                  int (*)(const double*,double*));
int c3sc_cost_load(struct C3SC * sc, char *);
int c3sc_cost_approx(struct C3SC *, double (*)(const double *, void *),
                     void *, int, const struct ApproxArgs *);
void c3sc_pol_solve(struct C3SC *,size_t, double,
                    int, const struct ApproxArgs *);
double c3sc_iter_vi(struct C3SC *,int,
                        const struct ApproxArgs *,struct C3SCDiagnostic *);
void * c3sc_get_dp(const struct C3SC *);
struct Cost * c3sc_get_cost(const struct C3SC *);
size_t c3sc_get_du(const struct C3SC *);

////////////////////////////////////////////////////


struct DPih * 
dpih_alloc(double,
           int (*)(double,const double*,const double*,double*,double*),
           int (*)(double,const double*,double*),
           int (*)(const double*,double*));
struct DPih * dpih_copy_deep(struct DPih *);
void dpih_free(struct DPih *);
void dpih_free_deep(struct DPih *);
void dpih_attach_mca(struct DPih *, struct MCA **);
void dpih_attach_cost(struct DPih *, struct Cost **);
void dpih_attach_cost_ow(struct DPih *, struct Cost *);
void dpih_attach_opt(struct DPih *, struct c3Opt **);

struct Cost * dpih_get_cost(struct DPih *);
struct Dyn * dpih_get_dyn(struct DPih *);
size_t dpih_get_d(const struct DPih *);

double dpih_rhs(struct DPih *,const double *,const double *,double *);
struct Cost * dpih_iter_vi(struct DPih *,int,const struct ApproxArgs *,
                           struct C3SCDiagnostic *);
struct Cost * dpih_iter_pol(struct DPih *, struct ImplicitPolicy *,
                            int, const struct ApproxArgs *);
void
dpih_iter_pol_solve(struct DPih *, struct ImplicitPolicy *,
                    size_t, double, int,const struct ApproxArgs *);

struct DPfh;

struct ImplicitPolicy * implicit_policy_alloc();
struct ImplicitPolicy * c3sc_create_implicit_policy(struct C3SC *);
void implicit_policy_add_transform(struct ImplicitPolicy *,size_t, void (*)(size_t, const double *, double *));
void implicit_policy_free(struct ImplicitPolicy *);
void implicit_policy_set_dp(struct ImplicitPolicy *, struct DPih *);
int implicit_policy_eval(struct ImplicitPolicy *,double,const double *, double *);
int implicit_policy_controller(double, const double *, double *, void *);


///////////////////////////////////////////////////////////////////////

struct C3SCDiagnostic * c3sc_diagnostic_init();
void c3sc_diagnostic_free(struct C3SCDiagnostic *);
int c3sc_diagnostic_save(struct C3SCDiagnostic *, 
                             char *, size_t);
void c3sc_diagnostic_vi_update(struct C3SCDiagnostic *,
                                   struct Cost *,
                                   struct DPih *,
                                   struct FunctionMonitor *);

#endif

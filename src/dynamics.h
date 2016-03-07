#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "c3sc_elements.h"
struct Drift * drift_alloc(size_t, size_t);
void drift_free(struct Drift *);
void drift_add_func(struct Drift *,
                    int (*)(double,double*,double*,
                             double*,double*,void*),
                    void *);
size_t drift_get_dx(struct Drift *);
int drift_eval(struct Drift *,double,double*,double*,
               double*,double*);

////////////////////////////////////
struct Diff * diff_alloc(size_t, size_t, size_t);
void diff_free(struct Diff *);
void diff_add_func(struct Diff *,
                   int (*)(double,double*,double*,
                           double*,double*,void*),
                   void *);
int diff_eval(struct Diff *,double,double*,double*,
              double*,double*);
size_t diff_get_dw(struct Diff *);


///////////////////////////////////////////////////////////

struct Dyn * dyn_alloc(struct Drift *, struct Diff *);
void dyn_free(struct Dyn *);
void dyn_free_deep(struct Dyn *);

void dyn_init_ref(struct Dyn*,struct Drift *, struct Diff *);
size_t dyn_get_dx(struct Dyn *);
size_t dyn_get_dw(struct Dyn *);
size_t dyn_get_du(struct Dyn *);
int dyn_eval(struct Dyn *,double,double*,double*,double*,
             double*,double*,double*);



////////////////////////////////////////////////////////////

#endif

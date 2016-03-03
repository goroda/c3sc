#ifndef DP_H
#define DP_H

#include "c3sc_elements.h"

#include "boundary.h"
#include "cost.h"
#include "tensmarkov.h"
#include "control.h"


struct DPih;
struct DPih * 
dpih_alloc(double,
           int (*)(double,double*,double*,double*,double*),
           int (*)(double,double*,double*));

void dpih_free(struct DPih *);
void dpih_attach_mca(struct DPih *, struct MCA *);
void dpih_attach_cost(struct DPih *, struct Cost *);
void dpih_attach_policy(struct DPih *, struct Policy *);

double dpih_rhs(struct DPih *,double *,double *);
struct Cost * dpih_iter_pol(struct DPih *,int);
#endif

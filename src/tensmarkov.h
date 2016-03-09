#ifndef TENS_MARKOV_H
#define TENS_MARKOV_H

#include "c3sc_elements.h"

#include "dynamics.h"
#include "boundary.h"
#include "cost.h"


enum NodeType {INBOUNDS,OUTBOUNDS,ONBOUNDS};

struct MCNode;
struct MCNode * mcnode_alloc(size_t);
struct MCNode ** mcnode_alloc_array(size_t);
void mcnode_free_array(struct MCNode **, size_t);
void mcnode_free(struct MCNode *);
struct MCNode * mcnode_init(size_t, double *);
size_t mcnode_get_d(struct MCNode *);
double * mcnode_get_xref(struct MCNode *);
double mcnode_get_pself(struct MCNode *);
size_t mcnode_get_N(struct MCNode *);
struct MCNode **  mcnode_get_neighbors(struct MCNode *);
double * mcnode_get_pref(struct MCNode *);

void mcnode_add_neighbors_hspace(struct MCNode *, 
                                 double *, double,
                                 double *, struct BoundInfo *);
void mcnode_add_gradients(struct MCNode *, size_t,
                          size_t, double *,
                          double **);
double mcnode_expectation(
    struct MCNode *,
    void (*)(size_t,double*,double**x,double*,void*),
    void *,double*);
void mcnode_sample_neighbor(struct MCNode *, double, double *);
void mcnode_print(struct MCNode *, FILE *, int);

/////////////////////////////////////////////////
struct MCA;
struct MCA * mca_alloc(size_t, size_t, size_t,double *);
void mca_free(struct MCA *);
void mca_free_deep(struct MCA *);
void mca_attach_dyn(struct MCA *, struct Dyn *);
struct Dyn * mca_get_dyn(struct MCA *);
void mca_attach_bound(struct MCA *, struct Boundary *);
size_t mca_get_dx(struct MCA *);
size_t mca_get_du(struct MCA *);
double
mca_expectation(struct MCA *,double,double *,double *,double *,
                double *,
                void(*)(size_t,double*,double**,double*,void*),
                void *,int*,double*);

struct MCNode *
mca_inbound_node(struct MCA*,double,double*,
                 double*,double*,double*,double*,
                 struct BoundInfo *);
struct MCNode *
mca_outbound_node(struct MCA *, double, double *);
struct MCNode *
mca_get_node(struct MCA *,double,double *,
             double *,double *,double*,
             enum NodeType *,double*);

void mca_step(struct MCA *, double, double *, double *,double,
              double *, double *);

/* double tensor_mm_cost(struct TensorMM *, */
/*                       double, double *, */
/*                       double *, */
/*                       struct Cost *, */
/*                       double *, */
/*                       double *); */
#endif

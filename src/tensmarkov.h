#ifndef TENS_MARKOV_H
#define TENS_MARKOV_H

#include "dynamics.h"
#include "boundary.h"
#include "cost.h"

struct Node
{
    size_t d;
    double * x;
    double c; // cost
    double * pert; // (- +) in each direction //2d
    int * use;
    double * cost;
};
struct Node * node_init(size_t, const double *, const double *,
                        struct Cost *, struct BoundInfo *);
void node_free(struct Node *);

enum NodeType {INBOUNDS,OUTBOUNDS,ONBOUNDS};

struct MCNode;
struct MCNode * mcnode_alloc(struct Node *, size_t);
void mcnode_free(struct MCNode *);
size_t mcnode_get_d(const struct MCNode *);
double * mcnode_get_xref(const struct MCNode *);
double mcnode_get_pself(const struct MCNode *);
size_t mcnode_get_du(const struct MCNode *);
double * mcnode_get_gpself(const struct MCNode *);

void mcnode_print(struct MCNode *, FILE *, int);

/////////////////////////////////////////////////
struct MCA;
struct MCA * mca_alloc(size_t, size_t, size_t,double *);
struct MCA * mca_copy_deep(struct MCA *);
void mca_set_newh(struct MCA *, const double *);
void mca_free(struct MCA *);
void mca_free_deep(struct MCA *);
void mca_attach_dyn(struct MCA *, struct Dyn *);
struct Dyn * mca_get_dyn(struct MCA *);
struct Boundary * mca_get_boundary(struct MCA *);
double * mca_get_h(struct MCA *);
void mca_attach_bound(struct MCA *, struct Boundary *);
size_t mca_get_dx(struct MCA *);
size_t mca_get_du(struct MCA *);

double
mca_expectation_cost2(struct MCA *, double,
                      struct Node *,
                      const double *,
                      double *, double *,
                      double *,
                      int *);
/* void mca_step(struct MCA *, double, double *, double *,double, */
/*               double *, double *); */

/* double tensor_mm_cost(struct TensorMM *, */
/*                       double, double *, */
/*                       double *, */
/*                       struct Cost *, */
/*                       double *, */
/*                       double *); */
#endif

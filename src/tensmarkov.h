#ifndef TENS_MARKOV_H
#define TENS_MARKOV_H

#include "c3sc_elements.h"

#include "dynamics.h"
#include "boundary.h"
#include "cost.h"


enum NodeType {INBOUNDS,OUTBOUNDS,ONBOUNDS};

struct MCNode;
struct MCNode * mcnode_alloc(size_t);
/* struct MCNode ** mcnode_alloc_array(size_t); */
/* void mcnode_free_array(struct MCNode **, size_t); */
void mcnode_free(struct MCNode *);
struct MCNode * mcnode_init(size_t,double *);
void mcnode_set_pself(struct MCNode *,double);
void mcnode_set_gradient(struct MCNode *,size_t,double*);
void mcnode_init_self_grad(struct MCNode *, size_t);
size_t mcnode_get_d(const struct MCNode *);
double * mcnode_get_xref(const struct MCNode *);
double mcnode_get_pself(const struct MCNode *);
size_t mcnode_get_du(const struct MCNode *);
double * mcnode_get_gpself(const struct MCNode *);
struct MCNList * mcnode_get_neigh(const struct MCNode *);
void mcnode_prepend_neigh(struct MCNode *, size_t,
                          double, double, double *);

struct MCNList;
struct MCNList * mcnlist_alloc();
double mcnlist_get_p(const struct MCNList *);
double * mcnlist_get_gradp(const struct MCNList *);
size_t mcnlist_get_dir(const struct MCNList *);
double mcnlist_get_val(const struct MCNList *);
struct MCNList * mcnlist_get_next(const struct MCNList *);

void mcnlist_free(struct MCNList *);
struct MCNList * mcnlist_prepend(struct MCNList **, size_t,
                                 double, double, double *);

/* struct MCNode **  mcnode_get_neighbors(const struct MCNode *); */
/* double * mcnode_get_pref(const struct MCNode *); */

/* void mcnode_add_neighbors_hspace(struct MCNode *,  */
/*                                  const double *, */
/*                                  const double, */
/*                                  const double *, */
/*                                  const struct BoundInfo *); */
/* void mcnode_add_gradients(struct MCNode *, size_t, */
/*                           size_t, const double *, */
/*                           double **); */
double mcnode_expectation(
    const struct MCNode *,
    //void (*)(size_t,double*,double**x,double*,void*),
    double (*)(double,double*,void*),
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
//                void(*)(size_t,double*,double**,double*,void*),
                double(*)(double,double*,void*),
                void *,struct BoundInfo **,double*);

struct MCNode *
mca_inbound_node(struct MCA*,double,double*,
                 double*,double*,double*,double*);

struct MCNode *
mca_outbound_node(struct MCA *, double, double *);
struct MCNode *
mca_get_node(struct MCA *,double,double *,
             double *,double *,double*,
             struct BoundInfo **,double*);

void mca_step(struct MCA *, double, double *, double *,double,
              double *, double *);

/* double tensor_mm_cost(struct TensorMM *, */
/*                       double, double *, */
/*                       double *, */
/*                       struct Cost *, */
/*                       double *, */
/*                       double *); */
#endif

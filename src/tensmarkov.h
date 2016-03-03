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
double * mcnode_getx_ref(struct MCNode *);
void mcnode_add_neighbors_hspace(struct MCNode *, 
                                 double *, double,
                                 double *);
double mcnode_expectation(
    struct MCNode *,
    void (*)(size_t,double*,double**x,double*,void*),
    void *);
void mcnode_sample_neighbor(struct MCNode *, double, double *);
void mcnode_print(struct MCNode *, FILE *, int);


/////////////////////////////////////////////////
struct MCA;
struct MCA * mca_alloc(size_t, size_t, double *);
void mca_free(struct MCA *);
void mca_attach_dyn(struct MCA *, struct Dyn *);
void mca_attach_bound(struct MCA *, struct Boundary *);
enum NodeType mca_node_type(struct MCA *, double, double *);
double
mca_expectation(struct MCA *,double,double *, double *, double *,
                void (*)(size_t, double *, double **, double *, void *),
                void *, int*);

struct MCNode *
mca_inbound_node(struct MCA *,double,double *,double *,double *);
struct MCNode *
mca_outbound_node(struct MCA *, double, double *);
struct MCNode *
mca_get_node(struct MCA *,double,double *,
             double *,double *, enum NodeType *);

void mca_step(struct MCA *, double, double *, double *,double,
              double *, double *);

/* double tensor_mm_cost(struct TensorMM *, */
/*                       double, double *, */
/*                       double *, */
/*                       struct Cost *, */
/*                       double *, */
/*                       double *); */
#endif

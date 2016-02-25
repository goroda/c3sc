#ifndef c3sc_ELEMENTS_H
#define c3sc_ELEMENTS_H


struct LinTransform
{
    size_t d;
    double * slope;
    double * offset;
};
void lin_transform_eval(struct LinTransform *, double *, double *);
double lin_transform_get_slopei(struct LinTransform *, size_t);

struct Control
{
    size_t d;
    double * u;
};

struct State
{
    size_t d; // dimension
    double t; // time
    double * x; // state
    void * spec; // special arguments
};

struct Trajectory
{
    struct State * s;
    struct Control * u;
    struct Trajectory * next;
};

struct Policy
{
    size_t dx;
    size_t du;
    double * lbx;
    double * ubx;
    int (*feedback)(double, double *, double *);

    int trans;
    struct LinTransform * lt;
    double * space;
};


struct Drift
{
    size_t dx;
    size_t du;
    double * lbx;
    double * ubx;
    double * lbu;
    double * ubu;

    int (*b)(double,double *, double *, double *, void *);
    void * bargs;
};

struct Diff
{
    size_t dw;
    size_t dx;
    size_t du;
    double * lbx;
    double * ubx;
    double * lbu;
    double * ubu;

    int (*s)(double,double *, double *, double *, void *);
    void * sargs;
};

struct Dyn
{
    int trans;
    struct Drift * drift;
    struct Diff  * diff;

// slope and offset for going from
// [-1,1] -> [a,b]
// the drift and dynamics above are all defined on
// [a,b]
    struct LinTransform * lt;
    double * space;
};




#endif

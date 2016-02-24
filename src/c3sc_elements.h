#ifndef c3sc_ELEMENTS_H
#define c3sc_ELEMENTS_H

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
    struct Drift * drift;
    struct Diff  * diff;
};


// for normalizing states
struct StdDyn
{
    size_t d;
    double * slope; //slope
    double * off; //offset
    double * space;
    struct Dyn * dyn;
};

#endif

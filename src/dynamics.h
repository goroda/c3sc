#ifndef DYNAMICS_H
#define DYNAMICS_H


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

void drift_init(struct Drift *,size_t,size_t,
                double*,double*,double*,double*);
int drift_eval(struct Drift *,double, double *, double *,double*);

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
void diff_init(struct Diff *,size_t,size_t,size_t,
                double*,double*,double*,double*);
int diff_eval(struct Diff *, double, double *, double *,double*);

struct Dyn
{
    struct Drift * b;
    struct Diff  * s;
};

int dyn_eval(struct Dyn *, double, double *, double *,double *, double *);
#endif

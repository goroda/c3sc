#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include "c3.h"
#include "c3sc.h"

static char * program_name;

void print_code_usage (FILE *, int) __attribute__ ((noreturn));
void print_code_usage (FILE * stream, int exit_code)
{

    fprintf(stream, "Usage: %s options \n", program_name);
    fprintf(stream,
            " -h --help      Display this usage information.\n"
            " -d --directory Output directory (defaults to .)\n"
            " -v --verbose   Output words (default 0)\n"
        );
    exit (exit_code);
}

int f1(double t,double *x,double *u,double *out, double *jac, 
       void * args)
{
    (void)(t);
    (void)(args);
    (void)(jac);
    out[0] = x[1] + sin(4.0*x[0]);
    out[1] = u[0];
    return 0;
}

int s1(double t,double *x,double *u,double *out,double*jac,void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(jac);
    (void)(args);
    
    out[0] = sin(3.0*x[0]);
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = cos(8.0*x[1]);
    
    return 0;
}

int polfunc(double t, double * x, double * u)
{
    (void)(t);
    if (x[1] > 0){
        u[0] = -1.0;
    }
    else if (x[1] < 0){
        u[0] = 1.0;
    }
    else{
        u[0] = 0.0;
    }
    //u[0] = 0.0;
    return 0;
}

int main(int argc, char * argv[])
{
    int next_option;
    const char * const short_options = "hd:v:";
    const struct option long_options[] = {
        { "help"     , 0, NULL, 'h' },
        { "directory", 1, NULL, 'd' },
        { "verbose"  , 1, NULL, 'v' },
        { NULL       , 0, NULL, 0   }
    };
    program_name = argv[0];

    char * dirout = ".";
    int verbose = 0;

    do {
        next_option = getopt_long (argc, argv, short_options, long_options, NULL);
        switch (next_option)
        {
            case 'h': 
                print_code_usage(stdout, 0);
            case 'd':
                dirout = optarg;
                break;
            case 'v':
                verbose = strtol(optarg,NULL,10);
                break;
            case '?': // The user specified an invalid option 
                print_code_usage (stderr, 1);
            case -1: // Done with options. 
                break;
            default: // Something unexpected
                abort();
        }
    } while (next_option != -1);
    
    size_t dx = 2;
    size_t dw = 2;
    size_t du = 1;
    
    struct Drift drift;
    drift_init(&drift,dx,du,NULL,NULL,NULL,NULL);
    drift.b = f1;

    struct Diff diff;
    diff_init(&diff,dw,dx,du,NULL,NULL,NULL,NULL);
    diff.s = s1;

    struct Dyn dyn;
    dyn_init_ref(&dyn,&drift,&diff);

    double t0 = 0.0;
    double xs[2] = {0.2,0.5};
    double us[1] = {0.0};

    struct State * state = state_alloc();
    state_init(state,dx,t0,xs);

    struct Control * control = control_alloc();    
    control_init(control,du,us);

    struct Trajectory * traj = NULL;
    trajectory_add(&traj,state,control);

    struct Policy * pol = policy_alloc();
    policy_init(pol,dx,du,NULL,NULL);

    size_t nsteps = 100;
    double space[2 + 4];
    double dt = 1e-2;
    double noise[2];
    int res;
    for (size_t ii = 0; ii < nsteps; ii++){
        noise[0] = randn()*sqrt(dt);
        noise[1] = randn()*sqrt(dt);
        /* res = trajectory_step(traj,pol,&dyn,dt,"euler", */
        /*                       space,NULL,NULL); */
        
        res = trajectory_step(traj,pol,&dyn,dt,
                              "euler-maruyama",
                              space,noise,NULL);
        if (res != 0){
            break;
        }
    }

    if (verbose == 1){
        trajectory_print(traj,stdout,4);
    }

    char filename[256];
    sprintf(filename,"%s/%s.dat",dirout,"traj");
    FILE * fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);


    state_free(state);
    control_free(control);
    trajectory_free(traj);
    policy_free(pol);
    
    return 0;
}

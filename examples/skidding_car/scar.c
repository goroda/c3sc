#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include "c3.h"
#include "cdyn/src/simulate.h"
#include "cdyn/src/integrate.h"

#include "nodeutil.h"
#include "valuefunc.h"
#include "bellman.h"
#include "c3sc.h"

static char * program_name;

void print_code_usage (FILE *, int) __attribute__ ((noreturn));
void print_code_usage (FILE * stream, int exit_code)
{

    fprintf(stream, "Usage: %s options \n", program_name);
    fprintf(stream,
            " -h --help      Display this usage information.\n"
            " -d --directory Output directory (defaults to .)\n"
            " -n --nodes     Number of nodes (defaults to 20)\n"
            " -s --steps     Number of iterations (default 100)\n"
            " -v --verbose   Output words (default 0)\n"
            "                1 - output main file stuff\n"
            "                >1 - also output approximation info\n"
        );
    exit (exit_code);
}

int f1(double t, const double * x, const double * u, double * out,
       double * jac, void * args)
{
    (void)(args);
    (void)(t);
    
    // states 
    // x[0] : position in x coordinate
    // x[1] : position in y coordinate
    // x[2] : orientation (0, 2pi)
    // x[3] : speed (2, 5)

    // inputs
    // u[0] : steering angle (-15 pi / 180,  15 pi /180)
    // u[1] : acceleration (-1, 1)

    double L = 0.2;
    double vc = 8.0;
    double alpha = 2.0;


    double pre = (1.0 / (1.0 + (x[3]/vc))) * (x[3] / L);
    out[0] = x[3] * cos(x[2]);
    out[1] = x[3] * sin(x[2]);
    out[2] = pre * tan (u[0]);
    out[3] = alpha * u[1];

    if (jac != NULL){
        //df1/du
        jac[0] = 0.0;    
        jac[1] = 0.0;    
        jac[2] = pre * pow(cos(u[0]),-2);
        jac[3] = 0;

        jac[4] = 0.0;
        jac[5] = 0.0;
        jac[6] = 0.0;   
        jac[7] = alpha;
    }

    /* printf("x = "); dprint(4,x); */
    /* printf("u = "); dprint(2,u); */
    /* printf("out = "); dprint(4,out); */
    /* exit(1); */
    return 0;
}

int s1(double t,const double * x,const double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(args);

    for (size_t ii = 0; ii < 16; ii++){
        out[ii] = 0.0;
    }
    double vpos = 1e-2;
    double vorient = 1e-2;
    double vspeed = 1e-2;

    out[0] = vpos;
    out[5] = vpos;
    out[10] = vorient;
    out[15] = vspeed;

    if (grad != NULL){
        for (size_t ii = 0; ii < 16*2; ii++){
            grad[ii] = 0.0;
        }
    }
    return 0;
}

int stagecost(double t,const double * x,const double * u, double * out, 
              double * grad)
{
    (void)(t);

    *out = 0.0;
    *out = 1.0;
    (void)(x);
    (void)(u);
    
    if (grad!= NULL){
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
    }
    return 0;
}

int boundcost(double t, const double * x, double * out)
{
    /* printf("not here!!\n"); */
    (void)(t);
    (void)(x);
    *out = 0.0;
    *out = 10.0;
    return 0;
}

int ocost(const double * x,double * out)
{
    /* printf("absorbed \n"); */
    /* dprint(6,x); */
    (void)(x);
    *out = 0.0;
    return 0;
}

int startcost(size_t N, const double * x, double * out, void * args)
{
    (void)(args);
    (void)(x);
    for (size_t ii = 0; ii < N; ii++){
        out[ii] = 20.0;
    }
    return 0;
}

void state_transform(size_t ndim, const double * x, double * y)
{
    (void)(ndim);

    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    y[3] = x[3];
    /* int trans = 0; */
    if ((x[2] < M_PI ) && (x[2] > -M_PI)){
        y[2] = x[2];
    }
    else if (x[2] > M_PI){
        /* trans = 1; */
        while (y[2] > 2 * M_PI){
            y[2] -= 2.0*M_PI;
        }
        if (y[2] > M_PI){
            y[2] = y[2] - 2.0*M_PI;
        }
    }
    else if (x[2] < -M_PI){
        /* trans = 1; */
        while (y[2] < - 2.0* M_PI){
            y[2] += 2.0 * M_PI;
        }
        if (y[2] < -M_PI){
            y[2] += 2.0 * M_PI;
        }
    }
    /* if (trans == 1){ */
    /*     printf("out here x = (%G,%G,%G)\n",x[0],x[1],x[2]); */
    /*     printf("out here y = (%G,%G,%G)\n",y[0],y[1],y[2]); */
        
    /* } */
}

int main(int argc, char * argv[])
{
    int next_option;
    const char * const short_options = "hd:n:s:v:";
    const struct option long_options[] = {
        { "help"     , 0, NULL, 'h' },
        { "directory", 1, NULL, 'd' },
        { "nodes"    , 1, NULL, 'n' },
        { "steps"    , 1, NULL, 's' },
        { "verbose"  , 1, NULL, 'v' },
        { NULL       , 0, NULL, 0   }
    };
    program_name = argv[0];

    char * dirout = ".";
    int verbose = 0;
    size_t N = 40;
    size_t niter = 100;
    do {
        next_option = getopt_long (argc, argv, short_options, long_options, NULL);
        switch (next_option)
        {
            case 'h': 
                print_code_usage(stdout, 0);
            case 'd':
                dirout = optarg;
                break;
            case 'n':
                N = strtoul(optarg,NULL,10);
                break;
            case 's':
                niter = strtoul(optarg,NULL,10);
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

    size_t dx = 4;
    size_t dw = 4;
    size_t du = 2;
    double lb[4] = {-4.0, -4.0, -M_PI, 2.0};
    double ub[4] = {4.0, 4.0, M_PI,  5.0};
    size_t Narr[4] = {N,N,N,N};
    double beta = 0.0;

    double lbu[2] = {-15.0*M_PI/180.0, -1.0};
    double ubu[2] = {15.0*M_PI/180.0, 1.0};
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    c3opt_set_absxtol(opt,1e-10);
    c3opt_set_relftol(opt,1e-10);
    c3opt_set_gtol(opt,0.0);
    c3opt_ls_set_maxiter(opt,10);
    c3opt_ls_set_alpha(opt,0.1);
    c3opt_ls_set_beta(opt,0.2);
    c3opt_set_verbose(opt,0);

    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-5);
    approx_args_set_round_tol(aargs,1e-5);
    approx_args_set_startrank(aargs,5);
    approx_args_set_adapt(aargs,1);
    approx_args_set_kickrank(aargs,5);
    approx_args_set_maxrank(aargs,5);


    // setup problem
    struct C3Control * c3c = c3control_create(dx,du,dw,lb,ub,Narr,beta);
    c3control_add_drift(c3c,f1,NULL);
    c3control_add_diff(c3c,s1,NULL);
    c3control_add_stagecost(c3c,stagecost);
    c3control_add_boundcost(c3c,boundcost);
    c3control_add_obscost(c3c,ocost);

    c3control_set_external_boundary(c3c,2,"periodic");
    c3control_set_external_boundary(c3c,3,"reflect");

    // possible obstacle
    double w = 0.5;
    double center[4] = {0.0,0.0,0.0,0.0};
    double width[4] = {w,w,ub[2]-lb[2],ub[3]-lb[3]};
    c3control_add_obstacle(c3c,center,width);

    char filename[256];
    sprintf(filename,"%s/%s.c3",dirout,"cost");
    double ** xgrid = c3control_get_xgrid(c3c);
    struct ValueF * cost = valuef_load(filename,Narr,xgrid);
    if (cost == NULL){
        cost = c3control_init_value(c3c,startcost,NULL,aargs,0);
    }

    printf("\n\n\n\n\n\n\n");
    printf("Start Solver Iterations\n");
    printf("\n\n\n\n\n\n\n\n");
    size_t maxiter_vi = niter;
    double abs_conv_vi = 1e-3;
    size_t maxiter_pi = 10;
    double abs_conv_pi = 1e-2;
    struct Diag * diag = NULL;
    char filename_diag[256];
    sprintf(filename_diag,"%s/%s.dat",dirout,"diagnostic");
    printf("filename = %s\n",filename_diag);
    for (size_t ii = 0; ii < maxiter_vi; ii++){

        struct ValueF * next = c3control_pi_solve(c3c,maxiter_pi,abs_conv_pi,
                                                  cost,aargs,opt,verbose,&diag);

        valuef_destroy(cost); cost = NULL;
        cost = c3control_vi_solve(c3c,1,abs_conv_vi,next,aargs,opt,verbose,&diag);
        valuef_destroy(next); next = NULL;

        sprintf(filename,"%s/%s.c3",dirout,"cost");
        int saved = valuef_save(cost,filename);
        assert (saved == 0);

        saved = diag_save(diag,filename_diag);
        assert (saved == 0);

        size_t * ranks = valuef_get_ranks(cost);
        if (verbose != 0){
            printf("ii=%zu ranks=",ii);
            iprint_sz(5,ranks);
        }

    }

    c3control_add_policy_sim(c3c,cost,opt,state_transform);

/* //    printf("created policy\n"); */
    /* char odename[256] = "rkf45"; */
    char odename[256] = "forward-euler";
    struct Integrator * ode_sys =
        integrator_create_controlled(4,2,f1,NULL,c3control_controller,c3c);
    integrator_set_type(ode_sys,odename);
    /* integrator_set_adaptive_opts(ode_sys,1e-5,1e-2,1e-7); */
    integrator_set_dt(ode_sys,1e-3);
    integrator_set_verbose(ode_sys,0);
//    printf("initialized integrator\n");

    double time = 0.0;
    double state[4] = {1.0, 1.0, -M_PI/2.0, 3.0};
    double con[2] = {0.0, 0.0};

    struct Trajectory * traj = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj,4,2,time,state,con);
    printf("initialized trajectory\n");


    double dt = 1e-2;
    double final_time = 1e0;
    int res;
    while (time < final_time){
        /* printf("time = %G\n",time); */
        res = trajectory_step(traj,ode_sys,dt);
        if (res != 0){
            break;
        }
//        assert(res == 0);
        time = time + dt;
    }

    if (verbose == 1){
        trajectory_print(traj,stdout,4);
    }
    
    sprintf(filename,"%s/%s.dat",dirout,"traj");
    FILE * fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);
    integrator_destroy(ode_sys); ode_sys = NULL;
    trajectory_free(traj); traj = NULL;

    
    valuef_destroy(cost); cost = NULL;
    c3control_destroy(c3c); c3c = NULL;
    diag_destroy(&diag); diag = NULL;
    c3opt_free(opt); opt = NULL;
    approx_args_free(aargs); aargs = NULL;

    return 0;
}

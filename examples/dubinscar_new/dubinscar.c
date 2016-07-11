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
            " -n --nodes     Number of nodes (defaults to 10)\n"
            " -s --steps     Number of iterations (default 100)\n"
            " -v --verbose   Output words (default 0)\n"
            "                1 - output main file stuff\n"
            "                >1 - also output approximation info\n"
        );
    exit (exit_code);
}

int f1(double t, const double * x, const double * u,
       double * out,
       double * jac, 
       void * args)
{
    (void)(t);
    (void)(args);
        
    out[0] = cos(x[2]);
    out[1] = sin(x[2]);
    out[2] = u[0];
    
    if (jac != NULL){
        jac[0] = 0.0;
        jac[1] = 0.0;
        jac[2] = 1.0;
    }
    
    return 0;
}

int s1(double t,const double * x,const double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(args);

    double val1 = 1e0;
    double val2 = 1e-2;
    out[0] = val1; out[3] = 0e0; out[6] = 0.0;
    out[1] = 0.0; out[4] = val1; out[7] = 0.0;
    out[2] = 0.0; out[5] = 0e0; out[8] = val2;

    if (grad != NULL){
        for (size_t ii = 0; ii < 3*3; ii++){
            grad[ii] = 0.0;
        }
    }
    return 0;
}

int stagecost(double t, const double * x, const double * u, double * out, 
              double * grad)
{
    (void)(t);
    (void)(u);
    (void)(x);
    
    *out = 0.0;
    //*out += pow(x[0],2) + pow(x[1],2);
    *out = 1.0;

    if (grad!= NULL){
        grad[0] = 0.0; 
    }

    return 0;
}

// cost of going into the boundaries
int boundcost(double t, const double * x, double * out)
{

    (void)(t);
    (void)(x);
    *out = 0.0;
    *out = 10.0;
    return 0;
}

// cost of hitting obstacle
int ocost(const double * x,double * out)
{
    /* printf("absorbed\n"); */
    /* dprint(3,x); */
    (void)(x);
    
    *out = 0.0;
    return 0;
}

// starting cost
int startcost(size_t N, const double * x, double * out, void * args)
{
    (void)(args);
    (void)(x);
    //double out = 10.0*x[0]*x[0] + 10.0*x[1]*x[1];
    for (size_t ii = 0; ii < N; ii++){
        out[ii] = 20.0;
    }
    return 0;
}

// needed for simulation
int f1sym(double t, const double * x, const double * u,
          double * out,
          double * jac, void * args)
{

    (void)(args);

    if (fabs(x[0]) < 0.1){
        if (fabs(x[1]) < 0.1){
            return 0;
        }
    }
    out[0] = cos(x[2]);
    out[1] = sin(x[2]);
    out[2] = u[0];

    (void)(t);

    /* printf("t=%G,u=%G\n",t,u[0]); */
    /* printf("x = "); dprint(3,x); */
    /* printf("out = "); dprint(3,out); */
    /* printf("----\n"); */
    
    if (jac != NULL){
        jac[0] = 0.0;
        jac[1] = 0.0;
        jac[2] = 1.0;
    }
    
    return 0;
}

void state_transform(size_t ndim, const double * x, double * y)
{
    (void)(ndim);

    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
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

void print_cost(char * filename, struct ValueF * cost, size_t N1, size_t N2, 
                double * lb, double * ub, double angle)
{

    FILE * fp2 =  fopen(filename, "w");
    if (fp2 == NULL){
        fprintf(stderr, "cat: can't open %s\n", filename);
        exit(0);
    }

    fprintf(fp2,"x y f\n");
    double * xtest = linspace(lb[0],ub[0],N1);
    double * ytest = linspace(lb[1],ub[1],N2);

    double pt3[3];
    double v2;
    for (size_t zz = 0; zz < N1; zz++){
        for (size_t jj = 0; jj < N2; jj++){
            pt3[0] = xtest[zz]; pt3[1] = ytest[jj];
            pt3[2] = angle;
            v2 = valuef_eval(cost,pt3);
            fprintf(fp2, "%3.5f %3.5f %3.5f \n",
                    xtest[zz],ytest[jj],v2);
        }
        fprintf(fp2,"\n");
    }
    free(xtest); xtest = NULL;
    free(ytest); ytest = NULL;

    fclose(fp2);
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
    size_t N = 30;
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

    size_t dx = 3;
    size_t dw = 3;
    size_t du = 1;
    double lb[3] = {-4.0, -4.0,-M_PI};
    double ub[3] = {4.0, 4.0, M_PI};
    size_t Narr[3] = {N, N, N};
    double beta = 0.0;
    
    struct c3Opt * opt = c3opt_alloc(BRUTEFORCE,du);
    size_t dopts = 3;
    double uopts[3] = {-1.0,0.0,1.0};
    c3opt_set_brute_force_vals(opt,dopts,uopts);
        
    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-5);
    approx_args_set_round_tol(aargs,1e-5);
    approx_args_set_kickrank(aargs,5);
    approx_args_set_maxrank(aargs,15);
    approx_args_set_startrank(aargs,15);


    // setup problem
    struct C3Control * c3c = c3control_create(dx,du,dw,lb,ub,Narr,beta);
    c3control_add_drift(c3c,f1,NULL);
    c3control_add_diff(c3c,s1,NULL);
    c3control_add_stagecost(c3c,stagecost);
    c3control_add_boundcost(c3c,boundcost);
    c3control_add_obscost(c3c,ocost);

    c3control_set_external_boundary(c3c,2,"periodic");

    // possible obstacle
    double w = 0.5;
    double center[3] = {0.0,0.0,0.0};
    double width[3] = {w,w,2.0*M_PI};
    c3control_add_obstacle(c3c,center,width);

    char filename[256];
    sprintf(filename,"%s/%s.c3",dirout,"cost");
    double ** xgrid = c3control_get_xgrid(c3c);
    struct ValueF * cost = valuef_load(filename,Narr,xgrid);
    if (cost == NULL){
        cost = c3control_init_value(c3c,startcost,NULL,aargs,0);
    }

    size_t maxiter_vi = niter+1;
    double abs_conv_vi = 1e-3;
    size_t maxiter_pi = 10;
    double abs_conv_pi = 1e-2;
    struct Diag * diag = NULL;
    char filename_diag[256] = "diagnostic.dat";

    printf("\n\n\n\n\n\n\n\n\n\n");
    printf("Start Solver Iterations\n");
    printf("\n\n\n\n\n\n\n\n\n\n");
    for (size_t ii = 0; ii < maxiter_vi; ii++){
    /* for (size_t ii = 0; ii < 0; ii++){ */
        struct ValueF * next = c3control_pi_solve(c3c,maxiter_pi,abs_conv_pi,
                                                  cost,aargs,opt,verbose,&diag);

        valuef_destroy(cost); cost = NULL;
        cost = c3control_vi_solve(c3c,1,abs_conv_vi,next,aargs,opt,verbose,&diag);
        valuef_destroy(next); next = NULL;

        sprintf(filename,"%s/%s.c3",dirout,"cost");
        int saved = valuef_save(cost,filename);
        assert (saved == 0);

        if (verbose != 0){
            printf("ii=%zu ranks=",ii);
            size_t * ranks = valuef_get_ranks(cost);
            iprint_sz(4,ranks);
        }
        
        saved = diag_save(diag,filename_diag);
        assert (saved == 0);
    }

    sprintf(filename,"%s/%s_angle=1.6.dat",dirout,"cost");
    print_cost(filename,cost,30,30,lb,ub,1.6);

    sprintf(filename,"%s/%s_angle=-1.6.dat",dirout,"cost");
    print_cost(filename,cost,30,30,lb,ub,-1.6);

    sprintf(filename,"%s/%s_angle=pi_3.dat",dirout,"cost");
    print_cost(filename,cost,30,30,lb,ub,M_PI/3.0);


    c3control_add_policy_sim(c3c,cost,opt,state_transform);

    printf("created policy\n");
    char odename[256] = "rk4";
    struct Integrator * ode_sys = NULL;
    ode_sys = integrator_create_controlled(
        3,1,f1sym, NULL,c3control_controller,c3c);
    integrator_set_type(ode_sys,odename);
    integrator_set_dt(ode_sys,1e-2);
    integrator_set_verbose(ode_sys,0);
    printf("initialized integrator\n");
    

    // Initialize trajectories for filter and for observations
    double time = 0.0;
    double state[3] = {-1.0, 0.0, 3*M_PI/4.0};
    double con[1] = {0.0};
    
    struct Trajectory * traj = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj,3,1,time,state,con);
    printf("initialized trajectory\n");

    double final_time = 1e1;
    double dt = 1e-2;
    int res;
    while (time < final_time){
        /* printf("time = %G\n",time); */
        res = trajectory_step(traj,ode_sys,dt);
        double * ls = trajectory_get_last_state(traj);
        if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){
            break;
        }
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

    /* double final_time = 1e1; */
    /* double dt = 1e-2; */
    /* int res; */
    double state2[3] = {3.0, 2.0, -M_PI/2};
    /* double state2[3] = {-1.0, -1.0, M_PI/3.0}; */
    double con2[1] = {0.0};
    time = 0.0;
    struct Trajectory * traj2 = NULL;
    trajectory_add(&traj2,3,1,time,state2,con2);
//    final_time = 1e1;
    time = 0.0;
    while (time < final_time){
        res = trajectory_step(traj2,ode_sys,dt);
        double * ls = trajectory_get_last_state(traj2);
        if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){
            break;
        }
        if (res != 0){
            break;
        }
        time = time + dt;
    }
    if (verbose == 1){
        trajectory_print(traj2,stdout,4);
    }
    sprintf(filename,"%s/%s.dat",dirout,"traj2");
    fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj2,fp,4);
    fclose(fp);

    double state3[3] = {-3.0, -1.0, 0.3};
    double con3[1] = {0.0};
    time = 0.0;
    struct Trajectory * traj3 = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj3,3,1,time,state3,con3);
    printf("initialized trajectory\n");
    final_time = 1e1;
    time = 0.0;
    while (time < final_time){
        res = trajectory_step(traj3,ode_sys,dt);
        double * ls = trajectory_get_last_state(traj3);
        if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){
            break;
        }
        if (res != 0){
            break;
        }
        time = time + dt;
    }
    sprintf(filename,"%s/%s.dat",dirout,"traj3");
    fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj3,fp,4);
    fclose(fp);

    integrator_destroy(ode_sys);
    trajectory_free(traj);
    trajectory_free(traj2);

    valuef_destroy(cost); cost = NULL;
    c3control_destroy(c3c); c3c = NULL;
    diag_destroy(&diag); diag = NULL;
    c3opt_free(opt); opt = NULL;
    approx_args_free(aargs); aargs = NULL;

    return 0;
}

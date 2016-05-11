#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include "c3.h"
#include "cdyn/src/simulate.h"
#include "cdyn/src/integrate.h"

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
       double * jac, void * args)
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
    *out = 100.0;
    return 0;
}

// cost of hitting obstacle
int ocost(const double * x,double * out)
{
    //dprint(3,x);
    (void)(x);
    
    *out = 0.0;
    return 0;
}

// starting cost
double startcost(const double * x, void * args)
{
    (void)(args);
    (void)(x);
    //double out = 10.0*x[0]*x[0] + 10.0*x[1]*x[1];
    double out = 20.0;
    return out;
}

// needed for simulation
int f1sym(double t, const double * x, const double * u,
          double * out,
          double * jac, void * args)
{
    (void)(t);
    (void)(args);

    if (fabs(x[0]) < 0.1){
        if (fabs(x[1]) < 0.1){
            return 0;
        }
    }
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

void state_transform(size_t ndim, const double * x, double * y)
{
    (void)(ndim);
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    if ((x[2] < M_PI ) && (x[2] > -M_PI)){
        y[2] = x[2];
    }
    else if (x[2] > M_PI){
        while (y[2] > 2 * M_PI){
            y[2] -= 2.0*M_PI;
        }
        if (y[2] > M_PI){
            y[2] = y[2] - 2.0*M_PI;
        }
    }
    else if (x[2] < -M_PI){
        while (y[2] < - 2.0* M_PI){
            y[2] += 2.0 * M_PI;
        }
        if (y[2] < -M_PI){
            y[2] += 2.0 * M_PI;
        }
    }
}

void print_cost(FILE * fp2, struct Cost * cost, size_t N1, size_t N2, double * lb, double * ub, double angle)
{

    fprintf(fp2,"x y f\n");
    double * xtest = linspace(lb[0],ub[0],N1);
    double * ytest = linspace(lb[1],ub[1],N2);

    double pt3[3];
    double v2;
    for (size_t zz = 0; zz < N1; zz++){
        for (size_t jj = 0; jj < N2; jj++){
            pt3[0] = xtest[zz]; pt3[1] = ytest[jj];
            pt3[2] = angle;
            cost_eval(cost,0.0,pt3,&v2);
            fprintf(fp2, "%3.5f %3.5f %3.5f \n",
                    xtest[zz],ytest[jj],v2);
        }
        fprintf(fp2,"\n");
    }
    free(xtest); xtest = NULL;
    free(ytest); ytest = NULL;
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

    char filename[256];
    
    size_t dx = 3;
    size_t dw = 3;
    size_t du = 1;
    double lb[3] = {-4.0, -4.0,-M_PI};
    double ub[3] = {4.0, 4.0, M_PI};
    size_t Narr[3] = {N, N, N};
    
    struct c3Opt * opt = c3opt_alloc(BRUTEFORCE,du);
    size_t dopts = 3;
    double uopts[3] = {-1.0,0.0,1.0};
    c3opt_set_brute_force_vals(opt,dopts,uopts);
        
    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-5);
    approx_args_set_round_tol(aargs,1e-5);
    approx_args_set_kickrank(aargs,5);

    double beta = 0.0;
    // setup problem
    c3sc sc = c3sc_create(IH,dx,du,dw);
    c3sc_set_state_bounds(sc,lb,ub);
    c3sc_set_external_boundary(sc,0,"reflect");
    c3sc_set_external_boundary(sc,1,"reflect");
    c3sc_set_external_boundary(sc,2,"periodic");
    // possible obstacle
    double w = 0.4;
    double center[3] = {0.0,0.0,0.0};
    double width[3] = {w,w,2.0*M_PI};
    c3sc_add_obstacle(sc,center,width);
    c3sc_add_dynamics(sc,f1,NULL,s1,NULL);
    c3sc_init_mca(sc,Narr);
    c3sc_attach_opt(sc,opt);
    c3sc_init_dp(sc,beta,stagecost,boundcost,ocost);

    sprintf(filename,"%s/%s.dat",dirout,"saved_cost.dat");
    int load_success = c3sc_cost_load(sc,filename);
    if (load_success != 0){
        c3sc_cost_approx(sc,startcost,NULL,0,aargs);
    }
    struct Cost * cost_init = c3sc_get_cost(sc);
    int saved = cost_save(cost_init,"saved_cost.dat");
    assert (saved == 0);
    printf("save cost\n");

    double solve_tol = 1e-4;
    size_t npol = 10;

    struct C3SCDiagnostic * diag = c3sc_diagnostic_init();
    FILE *fp;

    printf("\n\n\n\n\n\n\n\n\n\n");
    printf("Start Solver Iterations\n");
    printf("\n\n\n\n\n\n\n\n\n\n");
    for (size_t ii = 0; ii < niter+1; ii++){
    /* for (size_t ii = 0; ii < 0; ii++){ */

        if (ii > 0){
            c3sc_pol_solve(sc,npol,solve_tol,verbose,aargs);
        }

        double diff = c3sc_iter_vi(sc,verbose,aargs,diag);

        struct Cost * cost = c3sc_get_cost(sc);
        size_t * ranks = cost_get_ranks(cost);
        double normval = cost_norm2(cost);
        if (verbose != 0){
            printf("ii=%zu diff=%G,norm=%G diff/norm=%G, ranks=",ii,diff,normval,diff/normval);
            iprint_sz(4,ranks);
        }
        
//        printf("get cost\n");
        sprintf(filename,"%s/%s.dat",dirout,"saved_cost.dat");
        saved = cost_save(cost,filename);
        assert (saved == 0);

        sprintf(filename,"%s/%s.dat",dirout,"diagnostic");
        int dres = c3sc_diagnostic_save(diag,filename,4);
        assert (dres == 0);
        //      printf("save cost\n");
        /* if (verbose != 0){ */
        /*     printf("ii=%zu diff = %G\n",ii,diff); */
        /* } */
        if (diff < 1e-2){
            break;
        }
    }

    sprintf(filename,"%s/%s_angle=1.6.dat",dirout,"cost");
    fp =  fopen(filename, "w");
    if (fp == NULL){
        fprintf(stderr, "cat: can't open %s\n", filename);
        return 0;
    }
    size_t N1 = N, N2 = N;
    struct Cost * cost = c3sc_get_cost(sc);
    print_cost(fp,cost,N1,N2,lb,ub,1.6);
    fclose(fp);


    sprintf(filename,"%s/%s_angle=-1.6.dat",dirout,"cost");
    fp =  fopen(filename, "w");
    if (fp == NULL){
        fprintf(stderr, "cat: can't open %s\n", filename);
        return 0;
    }
    cost = c3sc_get_cost(sc);
    print_cost(fp,cost,N1,N2,lb,ub,-1.6);
    fclose(fp);

    sprintf(filename,"%s/%s_angle=pi_3.dat",dirout,"cost");
    fp =  fopen(filename, "w");
    if (fp == NULL){
        fprintf(stderr, "cat: can't open %s\n", filename);
        return 0;
    }
    cost = c3sc_get_cost(sc);
    print_cost(fp,cost,N1,N2,lb,ub,M_PI/3.0);
    fclose(fp);

    
    sprintf(filename,"%s/%s.dat",dirout,"diagnostic");
    int dres = c3sc_diagnostic_save(diag,filename,4);
    assert (dres == 0);

    struct ImplicitPolicy * pol = c3sc_create_implicit_policy(sc);
    implicit_policy_add_transform(pol,dx,state_transform);
    printf("created policy\n");
    char odename[256] = "rk4";
    struct Integrator * ode_sys = NULL;
    ode_sys = integrator_create_controlled(
        3,1,f1sym, NULL,implicit_policy_controller,pol);
    integrator_set_type(ode_sys,odename);
    integrator_set_dt(ode_sys,1e-2);
    integrator_set_verbose(ode_sys,0);
    printf("initialized integrator\n");
    
    // Initialize trajectories for filter and for observations
    double time = 0.0;
    double state[3] = {3.0, 0.0, -M_PI/2.0};
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
    fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);

    double state2[3] = {-1.0, 0.0, 3.0*M_PI/4.0};
    double con2[1] = {0.0};
    time = 0.0;
    struct Trajectory * traj2 = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj2,3,1,time,state2,con2);
    printf("initialized trajectory\n");
    final_time = 1e1;
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
        trajectory_print(traj,stdout,4);
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

    printf("cost ranks are ");
    cost = c3sc_get_cost(sc);
    size_t * ranks = cost_get_ranks(cost);
    iprint_sz(dx+1,ranks);
    
    integrator_destroy(ode_sys);
    trajectory_free(traj);
    trajectory_free(traj2);
    c3sc_destroy(sc);
    return 0;
}

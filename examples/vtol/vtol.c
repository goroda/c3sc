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
    double eta = 0.01;
    double tau = u[1];
    //0 1  2 3  4     5
    //x dx y dy theta dtheta    
    out[0] = x[1];
    out[1] = -u[0]*sin(x[4]) + eta * tau * cos (x[4]);
    out[2] = x[3];
    out[3] = u[0]*cos(x[4]) + eta * tau * sin (x[4]) - 9.81;
    out[4] = x[5];
    out[5] = u[1];

    if (jac != NULL){
        jac[0] = 0.0;
        jac[1] = -sin(x[4]);
        jac[2] = 0.0;
        jac[3] = cos(x[4]);
        jac[4] = 0.0;
        jac[5] = 0.0;

        jac[6] = 0.0;
        jac[7] = eta * cos(x[4]);
        jac[8] = 0.0;
        jac[9] = eta * sin(x[4]);
        jac[10] = 0.0;
        jac[11] = 1.0;
    }
    
    return 0;
}

int s1(double t,const double * x,
       const double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(args);


    for (size_t ii = 0; ii < 36; ii++){
        out[ii] = 0.0;
    }
    
    double val1 = 1e0;
    double val2 = 1e0;
    out[0] = val1;
    out[7] = val2;
    out[14] = val1;
    out[21] = val2;
    out[28] = val1;
    out[35] = val2;
   
    if (grad != NULL){
        for (size_t ii = 0; ii < 6*6; ii++){
            grad[ii] = 0.0;
        }
    }
    return 0;
}

int stagecost(double t,const double * x,
              const double * u,
              double * out, 
              double * grad)
{
    (void)(t);
    (void)(u);
    (void)(x);
    
    *out = 0.0;
    //*out += pow(x[0],2) + pow(x[1],2);
    *out = 1.0 + u[0]*u[0] + u[1]*u[1];

    if (grad!= NULL){
        grad[0] = 2.0 * u[0];
        grad[1] = 2.0 * u[1];
    }

    return 0;
}

// cost of going into the boundaries
int boundcost(double t,const  double * x, double * out)
{

    (void)(t);
    (void)(x);
    *out = 100.0;
    if (fabs(x[2] - 0.0) < 1e-13){
        *out = 0.0;
        for (size_t ii = 0; ii < 6; ii++){
            *out += x[ii]*x[ii];
        }
    }
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
    double out = 0.0;
    for (size_t ii = 0; ii < 6; ii++){
        out += x[ii]*x[ii];
    }

    return out;
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

/* void print_cost(FILE * fp2, struct Cost * cost, size_t N1, size_t N2, double * lb, double * ub, double angle) */
/* { */

/*     fprintf(fp2,"x y f\n"); */
/*     double * xtest = linspace(lb[0],ub[0],N1); */
/*     double * ytest = linspace(lb[1],ub[1],N2); */

/*     double pt3[3]; */
/*     double v2; */
/*     for (size_t zz = 0; zz < N1; zz++){ */
/*         for (size_t jj = 0; jj < N2; jj++){ */
/*             pt3[0] = xtest[zz]; pt3[1] = ytest[jj]; */
/*             pt3[2] = angle; */
/*             cost_eval(cost,0.0,pt3,&v2); */
/*             fprintf(fp2, "%3.5f %3.5f %3.5f \n", */
/*                     xtest[zz],ytest[jj],v2); */
/*         } */
/*         fprintf(fp2,"\n"); */
/*     } */
/*     free(xtest); xtest = NULL; */
/*     free(ytest); ytest = NULL; */
/* } */

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

    size_t dx = 6;
    size_t dw = 6;
    size_t du = 2;
    double lb[6] = {-4.0,-8.0,0.0,-1.0,-M_PI,-5.0};
    double ub[6] = {4.0,8.0,2.0,1.0,M_PI,5.0};
    size_t Narr[6] = {N, N, N, N, N, N};
    
    /* double lbu[3] = {-0.1, -0.1, -0.01}; */
    /* double ubu[3] = {0.1, 0.1, 0.01}; */
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    /* c3opt_add_lb(opt,lbu); */
    /* c3opt_add_ub(opt,ubu); */
    c3opt_set_absxtol(opt,1e-10);
    c3opt_set_relftol(opt,1e-10);
    c3opt_set_gtol(opt,1e-10);
    c3opt_set_verbose(opt,0);
        
    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-5);
    approx_args_set_round_tol(aargs,1e-4);
    approx_args_set_kickrank(aargs,15);

    double beta = 0.0;
    // setup problem
    c3sc sc = c3sc_create(IH,dx,du,dw);
    c3sc_set_state_bounds(sc,lb,ub);
    c3sc_set_external_boundary(sc,4,"periodic");
    
    // possible obstacle
    /* double w = 0.4; */
    /* double center[3] = {0.0,0.0,0.0}; */
    /* double width[3] = {w,w,2.0*M_PI}; */
    /* c3sc_add_obstacle(sc,center,width); */
    c3sc_add_dynamics(sc,f1,NULL,s1,NULL);
    c3sc_init_mca(sc,Narr);
    c3sc_attach_opt(sc,opt);
    c3sc_init_dp(sc,beta,stagecost,boundcost,ocost);
    int load_success = c3sc_cost_load(sc,"saved_cost.dat");
    printf("approximate initial costs\n");
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
    char filename[256];
    FILE *fp;

    printf("\n\n\n\n\n\n\n\n\n\n");
    printf("Start Solver Iterations\n");
    printf("\n\n\n\n\n\n\n\n\n\n");
    for (size_t ii = 0; ii < niter; ii++){
    /* for (size_t ii = 0; ii < 0; ii++){ */

        if (ii > 5){
            c3sc_pol_solve(sc,npol,solve_tol,verbose-1,aargs);
        }

        double diff = c3sc_iter_vi(sc,verbose-1,aargs,diag);

        struct Cost * cost = c3sc_get_cost(sc);
        size_t * ranks = cost_get_ranks(cost);
        iprint_sz(7,ranks);
//        printf("get cost\n");
        saved = cost_save(cost,"saved_cost.dat");
        assert (saved == 0);
        //      printf("save cost\n");

        sprintf(filename,"%s/%s.dat",dirout,"diagnostic");
        int dres = c3sc_diagnostic_save(diag,filename,4);
        assert (dres == 0);
        if (verbose != 0){
            printf("ii=%zu diff = %G\n",ii,diff);
        }
        if (diff < 1e-2){
            break;
        }
    }


    struct ImplicitPolicy * pol = c3sc_create_implicit_policy(sc);
    implicit_policy_add_transform(pol,dx,state_transform);
    printf("created policy\n");
    char odename[256] = "forward-euler";
    struct Integrator * ode_sys = NULL;
    ode_sys = integrator_create_controlled(
        6,2,f1, NULL,implicit_policy_controller,pol);
    integrator_set_type(ode_sys,odename);
    integrator_set_dt(ode_sys,1e-3);
    integrator_set_verbose(ode_sys,0);
    printf("initialized integrator\n");
    
    // Initialize trajectories for filter and for observations
    double time = 0.0;
    double state[6] = {2.0,0.0,1.0,0.0,0.1,0.0};
    double con[2] = {0.0};
    
    struct Trajectory * traj = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj,6,2,time,state,con);
    printf("initialized trajectory\n");

    double final_time = 1e0;
    double dt = 1e-2;
    int res;
    while (time < final_time){
        printf("time = %G\n",time);
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
    fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);

    /* double state2[3] = {-1.0, 0.0, 3.0*M_PI/4.0}; */
    /* double con2[1] = {0.0}; */
    /* time = 0.0; */
    /* struct Trajectory * traj2 = NULL; */
    /* printf("add trajectory\n"); */
    /* trajectory_add(&traj2,3,1,time,state2,con2); */
    /* printf("initialized trajectory\n"); */
    /* final_time = 1e1; */
    /* time = 0.0; */
    /* while (time < final_time){ */
    /*     res = trajectory_step(traj2,ode_sys,dt); */
    /*     double * ls = trajectory_get_last_state(traj2); */
    /*     if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){ */
    /*         break; */
    /*     } */
    /*     if (res != 0){ */
    /*         break; */
    /*     } */
    /*     time = time + dt; */
    /* } */
    /* if (verbose == 1){ */
    /*     trajectory_print(traj,stdout,4); */
    /* } */
    /* sprintf(filename,"%s/%s.dat",dirout,"traj2"); */
    /* fp = fopen(filename,"w"); */
    /* assert (fp != NULL); */
    /* trajectory_print(traj2,fp,4); */
    /* fclose(fp); */
    /* trajectory_free(traj2); */

    /* double state3[3] = {-3.0, -1.0, 0.3}; */
    /* double con3[1] = {0.0}; */
    /* time = 0.0; */
    /* struct Trajectory * traj3 = NULL; */
    /* printf("add trajectory\n"); */
    /* trajectory_add(&traj3,3,1,time,state3,con3); */
    /* printf("initialized trajectory\n"); */
    /* final_time = 1e1; */
    /* time = 0.0; */
    /* while (time < final_time){ */
    /*     res = trajectory_step(traj3,ode_sys,dt); */
    /*     double * ls = trajectory_get_last_state(traj3); */
    /*     if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){ */
    /*         break; */
    /*     } */
    /*     if (res != 0){ */
    /*         break; */
    /*     } */
    /*     time = time + dt; */
    /* } */
    /* sprintf(filename,"%s/%s.dat",dirout,"traj3"); */
    /* fp = fopen(filename,"w"); */
    /* assert (fp != NULL); */
    /* trajectory_print(traj3,fp,4); */
    /* fclose(fp); */

    printf("cost ranks are ");
    struct Cost * lcost = c3sc_get_cost(sc);
    size_t * ranks = cost_get_ranks(lcost);
    iprint_sz(dx+1,ranks);
    
    integrator_destroy(ode_sys);
    trajectory_free(traj);

    c3sc_destroy(sc);
    return 0;
}

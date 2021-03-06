#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include <time.h>

#include "c3/c3.h"
#include "cdyn/simulate.h"
#include "cdyn/integrate.h"

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
            " -a --lbu       Lower bound of control (default -1.0)\n"
            " -b --ubu       Upper bound of control (default 1.0)\n"
            " -r --diff1     Size of diffusion for x (default 0.01)\n"
            " -f --diff2     Size of diffusion for y (default 0.01)\n"
            " -s --steps     Number of iterations (default 100)\n"
            " -t --type      0 for absorb 1 for reflect (default 0)\n"
            " -e --epsilon   rounding tol (default is 1e-7)\n"
            " -v --verbose   Output words (default 0)\n"
            "                1 - output main file stuff\n"
            "                >1 - also output approximation info\n"
        );
    exit (exit_code);
}

static size_t dim;
void print_cost(char * filename, struct ValueF * cost,
                size_t N1, size_t N2, double * lb, double * ub)
{
    /* if (dim == 2){ */
        FILE * fp2 =  fopen(filename, "w");
        if (fp2 == NULL){
            fprintf(stderr, "cat: can't open %s\n", filename);
            exit(0);
        }

        fprintf(fp2,"x y f\n");
        double * xtest = linspace(lb[0],ub[0],N1);
        double * ytest = linspace(lb[1],ub[1],N2);

        double * pt = calloc_double(dim);
        double v2;
        for (size_t zz = 0; zz < N1; zz++){
            for (size_t jj = 0; jj < N2; jj++){
                pt[0] = xtest[zz]; pt[1] = ytest[jj];
                v2 = valuef_eval(cost,pt);
                fprintf(fp2, "%3.5f %3.5f %3.5f \n",
                        xtest[zz],ytest[jj],v2);
            }
            fprintf(fp2,"\n");
        }
        free(xtest); xtest = NULL;
        free(ytest); ytest = NULL;
        free(pt);
        fclose(fp2);
    /* } */
    /* else{ */
    /*     fprintf(stdout,"Not printing cost function, it is not 2D dynamics\n"); */
    /* } */
}


int f1(double t, const double * x, const double * u, double * out,
       double * jac, void * args)
{
    (void)(t);
    (void)(args);

    /* double a = 0.423; */
    /* double b = 2.0; */
    /* double c = 4.0; */
    double a = 0.1;
    double b = 0.1;
    double c = 14.0;
    out[0] = -x[1] - x[2];
    out[1] = x[0] + a*x[1] + u[0];
    out[2] = b + x[2]*(x[0] - c);

    if (jac != NULL){
        jac[0] = 0;
        jac[1] = 1;
        jac[2] = 0;
    }

    return 0;
}

int s1(double t,const double * x,const double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);

    double * val = args;

    for (size_t ii = 0; ii < dim*dim; ii++){
        out[ii] = 0.0;
    }
    for (size_t ii = 0; ii < dim-1; ii++){
        out[ii*dim + ii] = val[0];
    }
    size_t ii = dim-1;
    out[ii*dim + ii] = val[1];

    if (grad != NULL){
        for (ii = 0; ii < dim*dim; ii++){
            grad[ii] = 0.0;
        }
    }
    return 0;
}

int stagecost(double t,const double * x,const double * u, double * out, 
              double * grad)
{
    (void)(t);
    (void)(x);
    (void)(u);

    double scale = 1e2;
    *out = 0.0;
    for (size_t ii = 0; ii < dim; ii++){
        *out = *out + scale*x[ii]*x[ii];
    }

    double rho = 1e0;
    *out += rho * u[0] * u[0];
    if (grad!= NULL){
        grad[0] = 2.0 * rho * u[0];
    }
    return 0;
}

int boundcost(double t,const  double * x, double * out)
{

    (void)(t);
    (void)(x);
    *out = 0.0;
    *out = 1000.0;
    return 0;
}

static int obs = 0;
int ocost(const double * x,double * out)
{
    (void)(x);
    *out = 0.0;
    if (obs == 0){
        printf("Hit the obstacle!\n");
        obs=1;
    }
    return 0;
}

int startcost(size_t N, const double * x, double * out, void * args)
{
    (void)(args);
    double scale = 10.0;
    for (size_t ii = 0; ii < N; ii++){
        out[ii] = 0.0;
        for (size_t jj = 0; jj < dim; jj++){
            out[ii] += scale*(x[dim*ii+jj]*x[dim*ii+jj]);
        }
    }
    /* printf("done\n"); */
    return 0;
}

int main(int argc, char * argv[])
{
    srand(time(NULL));
    int next_option;
    const char * const short_options = "hd:n:a:b:r:f:s:e:v:";
    const struct option long_options[] = {
        { "help"     , 0, NULL, 'h' },
        { "directory", 1, NULL, 'd' },
        { "nodes"    , 1, NULL, 'n' },
        { "lbu"      , 1, NULL, 'a' },
        { "ubu"      , 1, NULL, 'b' },
        { "diff1"    , 1, NULL, 'r' },
        { "diff2"    , 1, NULL, 'f' },
        { "steps"    , 1, NULL, 's' },
        { "epsilon"  , 1, NULL, 'e' },
        { "verbose"  , 1, NULL, 'v' },
        { NULL       , 0, NULL, 0   }
    };
    program_name = argv[0];

    char * dirout = ".";
    int verbose = 0;
    size_t N = 20;
    size_t niter = 100;
    double lbu[1] = {-4.0};
    double ubu[1] = {4.0};
    double ss[3] = {1e0,1e0,1e0};
    double roundtol = 1e-8;
    dim = 3;
    do {
        next_option = getopt_long (argc, argv, short_options, long_options,
                                   NULL);
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
            case 'a':
                lbu[0] = strtod(optarg,NULL);
                break;
            case 'b':
                ubu[0] = strtod(optarg,NULL);
                break;
            case 'r':
                ss[0] = strtod(optarg,NULL);
                break;
            case 'f':
                ss[1] = strtod(optarg,NULL);
                break;
            case 's':
                niter = strtoul(optarg,NULL,10);
                break;
            case 'e':
                roundtol = strtod(optarg,NULL);
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

    size_t dx = dim;
    size_t dw = dim;
    size_t du = 1;

    double * lb = calloc_double(dx);
    double * ub = calloc_double(dx);
    size_t * Narr = calloc_size_t(dx);
    for (size_t ii = 0; ii < dx; ii++){
        lb[ii] = -1.0;
        ub[ii] = 1.0;
        Narr[ii] = N;
    }
    double beta = 0.1;

    // optimization arguments
    /* struct c3Opt * opt = c3opt_alloc(BRUTEFORCE,du); */
    /* size_t dopts = 3; */
    /* double uopts[3] = {lbu[0],0.0,ubu[0]}; */
    /* c3opt_set_brute_force_vals(opt,dopts,uopts); */
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    c3opt_set_absxtol(opt,1e-8);
    c3opt_set_relftol(opt,1e-8);
    c3opt_set_gtol(opt,1e-25);
    c3opt_set_verbose(opt,0);
    
    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-10);
    approx_args_set_round_tol(aargs,roundtol);
    approx_args_set_kickrank(aargs,2);
    approx_args_set_adapt(aargs,1);
    approx_args_set_startrank(aargs,2);
    approx_args_set_maxrank(aargs,N);

    // setup problem
    struct C3Control * c3c = c3control_create(dx,du,dw,lb,ub,Narr,beta);
    c3control_add_drift(c3c,f1,NULL);
    c3control_add_diff(c3c,s1,ss);
    c3control_add_stagecost(c3c,stagecost);
    c3control_add_boundcost(c3c,boundcost);
    c3control_add_obscost(c3c,ocost);
    

    for (size_t ii = 0; ii < dx; ii++){
        c3control_set_external_boundary(c3c,ii,"reflect");
    }


    double * center = calloc_double(dx);
    double * width = calloc_double(dx);
    for (size_t ii = 0; ii < dx; ii++){
        width[ii] = 0.4;
    }
    /* c3control_add_obstacle(c3c,center,width); */
    
    char filename[256];
    sprintf(filename,"%s/%s_%d.c3",dirout,"costfunc",0);
    double ** xgrid = c3control_get_xgrid(c3c);
    /* struct ValueF * cost = NULL; */
    struct ValueF * cost = valuef_load(filename,Narr,xgrid);
    if (cost == NULL){
        printf("NOT HERE\n");
        cost = c3control_init_value(c3c,startcost,NULL,aargs,0);
    }

    int saved = 0;
    /* sprintf(filename,"%s/%s_%d.c3",dirout,"costfunc",0); */
    /* int saved = valuef_save(cost,filename); */
    /* assert (saved == 0); */
    
    size_t maxiter_vi = niter;
    double abs_conv_vi = 1e-6;
    size_t maxiter_pi = 10;
    double abs_conv_pi = 1e-8;
    struct Diag * diag = NULL;
    char filename_diag[256];
    sprintf(filename_diag,"%s/dim%zu_%s",dirout,dx,"diagnostic.dat");
    printf("\n\n\nmaxiter_vi = %zu\n", maxiter_vi);
    for (size_t ii = 0; ii < maxiter_vi; ii++){
        if (verbose != 0){
            printf("\n\n\n\n\n\n\n");
        }

        // comment the next call and uncomment the copy for pure VI
        struct ValueF * next = c3control_pi_solve(c3c,maxiter_pi,
                                                  abs_conv_pi,
                                                  cost,aargs,opt,
                                                  verbose,&diag);

        /* struct ValueF * next = valuef_copy(cost); */

        valuef_destroy(cost); cost = NULL;
        cost = c3control_vi_solve(c3c,1,abs_conv_vi,next,aargs,opt,verbose,&diag);
        double err = valuef_norm2diff(cost,next);
        valuef_destroy(next); next = NULL;

        saved = diag_save(diag,filename_diag);
        assert (saved == 0);

        /* if (verbose != 0){ */
        printf("ii=%zu norm=%3.5E, err=%3.5E ranks=",ii,valuef_norm(cost),err);
        size_t * ranks = valuef_get_ranks(cost);
        iprint_sz(dx+1,ranks);

        sprintf(filename,"%s/%s_%zu.dat",dirout,"costfunc",ii+1);
        print_cost(filename,cost,Narr[0],Narr[1],lb,ub);


        sprintf(filename,"%s/%s_%d.c3",dirout,"costfunc",0);
        saved = valuef_save(cost,filename);
        assert (saved == 0);

        /* } */
        if (err < abs_conv_vi){
            if (ii > 10){
                printf("Converged!");
                /* exit(1); */
                break;
            }
        }
    }

    if (verbose != 0){
        size_t * ranks = valuef_get_ranks(cost);
        iprint_sz(dx+1,ranks);
    }

    c3control_add_policy_sim(c3c,cost,opt,NULL);
    
    printf("created policy\n");
    char odename[256] = "rk4";
    /* char odename[256] = "forward-euler"; */
    struct Integrator * ode_sys = NULL;
    ode_sys = integrator_create_controlled(dx,1,f1,NULL,
                                           c3control_controller,c3c);
    integrator_set_type(ode_sys,odename);
    integrator_set_dt(ode_sys,1e-3);
    integrator_set_verbose(ode_sys,0);

    // Initialize trajectories for state
    
    double time = 0.0;
    double * state = calloc_double(dx);
    for (size_t ii = 0; ii < dx; ii++){
        /* state[ii] = 0.5; */
        /* state[ii] = 0.8; */
        state[ii] = randu()*2.0-1.0;
    }

    /* double state[2] = {-0.5, -0.5}; */
    double con[1] = {0.0};
    // Integrate
    struct Trajectory * traj = NULL;
    trajectory_add(&traj,dx,1,time,state,con);
    double final_time = 1e1;
    double dt = 1e-2;
    int res;
    printf("Integrating Trajectory\n");
    while (time < final_time){
        res = trajectory_step(traj,ode_sys,dt);
        assert(res == 0);
        time = time + dt;
    }
    printf("Saving Trajectory\n");
    sprintf(filename,"%s/%s.dat",dirout,"traj1");
    FILE * fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);

    // cleanup integrator stuff
    integrator_destroy(ode_sys); ode_sys = NULL;
    trajectory_free(traj); traj = NULL;
    free(state); state = NULL;
    

    free(center); center = NULL;
    free(width); width = NULL;
    free(lb); lb = NULL;
    free(ub); ub = NULL;
    free(Narr); Narr = NULL;
    valuef_destroy(cost); cost = NULL;
    c3control_destroy(c3c); c3c = NULL;
    diag_destroy(&diag); diag = NULL;
    c3opt_free(opt); opt = NULL;
    approx_args_free(aargs); aargs = NULL;
    
    return 0;
}

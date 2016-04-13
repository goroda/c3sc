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
            " -n --nodes     Number of nodes (defaults to 10)\n"
            " -s --steps     Number of iterations (default 100)\n"
            " -v --verbose   Output words (default 0)\n"
            "                1 - output main file stuff\n"
            "                >1 - also output approximation info\n"
        );
    exit (exit_code);
}

int f1(double t, const double * x, const double * u, double * out, double * jac,
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

int s1(double t,double * x,double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(args);

    double val = 1e0;
    out[0] = val; out[3] = 0e0; out[6] = 0.0;
    out[1] = 0.0; out[4] = val; out[7] = 0.0;
    out[2] = 0.0; out[5] = 0e0; out[8] = val;

    if (grad != NULL){
        for (size_t ii = 0; ii < 3*3; ii++){
            grad[ii] = 0.0;
        }
    }
    return 0;
}

int trajboundcheck(double time, double * x, void * args, int * dirs)
{
    
    // two boundaries
    // outside of the box([-2,2]^2)
    // and the inside box [-0.1,0.1]^2
    (void)(time);
    (void)(args);
    double bound1 = 2.0;
    int out = 0;
    for (size_t ii = 0; ii < 2; ii++){
        dirs[ii] = 0;
        if (x[ii] > bound1){
            return 1;
        }
        else if (x[ii] < -bound1){
            return 1;
        }
        else if (fabs(x[ii] - bound1) < 1e-15){
            return 1;
        }
        else if (fabs(x[ii] + bound1) < 1e-15){
            return 1;
        }
    }
    
    
    return out;
}

int stagecost(double t, double * x, double * u, double * out, 
              double * grad)
{
    (void)(t);
    (void)(u);
    
    *out = 0.0;
    *out += pow(x[0],2) + pow(x[1],2);

    if (grad!= NULL){
        grad[0] = 0.0; 
    }

    return 0;
}

int boundcost(double t, double * x, double * out)
{

    (void)(t);
    (void)(x);
    *out = 0.0;
    *out = 3.0;
    return 0;
}

double startcost(double * x, void * args)
{
    (void)(args);
    (void)(x);
    return 0.2;

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
    size_t N = 10;
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
    double lb[3] = {-2.0, -2.0,-M_PI};
    double ub[3] = {2.0, 2.0, M_PI};
    size_t Narr[3] = {N, N, N};

    //double * xt = linspace(-M_PI,M_PI,N);
   
    double lbu[1] = {-1.0};
    double ubu[1] = {1.0};
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    c3opt_set_relftol(opt,1e-4);
    c3opt_set_verbose(opt,0);
    
    double beta = 1.0;

    // setup problem
    c3sc sc = c3sc_create(IH,dx,du,dw);
    c3sc_set_state_bounds(sc,lb,ub);
    c3sc_set_external_boundary(sc,2,"periodic");
    c3sc_add_dynamics(sc,f1,NULL,s1,NULL);
    c3sc_init_mca(sc,Narr);
    c3sc_attach_opt(sc,opt);
    c3sc_init_dp(sc,beta,stagecost,boundcost,NULL);

    struct DPih * dp = c3sc_get_dp(sc);
    struct Cost * cost = dpih_get_cost(dp);
    cost_approx(cost,startcost,NULL,verbose-1);
    for (size_t ii = 0; ii < niter+1; ii++){
        struct Cost * newcost = dpih_iter_vi(dp,verbose-1);
        cost_free(cost);
        cost = newcost;
        dpih_attach_cost(dp,cost);
        if (verbose != 0){
            printf("ii=%zu\n",ii);
        }
    }

    struct Policy * pol = dpih_iter_vi_pol(dp,verbose-1);

    double t0 = 0.0;
    double xs[3] = {0.5,0.0,-M_PI/2.0};
    double us[1] = {0.0};

    struct State * state = state_alloc();
    state_init(state,dx,t0,xs);

    struct Control * control = control_alloc();
    control_init(control,du,us);

    struct Trajectory * traj = NULL;
    trajectory_add(&traj,state,control);

    struct Dyn * dyn = dpih_get_dyn(dp);

    size_t nsteps = 1000;
    double space[3 + 9];
    double dt = 1e-2;
    int dirs[3];
    int res;
    for (size_t ii = 0; ii < nsteps; ii++){
        res = trajectory_step(traj,pol,dyn,dt,"euler",
                               space,NULL,NULL);
        if (res != 0){
            break;
        }
        struct State * scheck = trajectory_last_state(traj);
        double * xcheck = state_getx_ref(scheck);
        double tcheck = state_gett(scheck);

        res = trajboundcheck(tcheck,xcheck,NULL,dirs);
        /* res = trajectory_step(traj,pol,dyn,dt, */
        /*                       "euler-maruyama", */
        /*                       space,noise,NULL); */
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

    
    printf("cost ranks are ");
    size_t * ranks = cost_get_ranks(cost);
    iprint_sz(dx+1,ranks);

    printf("policy ranks are \n");
    for (size_t ii = 0; ii < du; ii++){
        ranks = policy_get_ranks(pol,ii);
        iprint_sz(dx+1,ranks);
    }

/*     state_free(state); state = NULL; */
/*     control_free(control); control = NULL; */
/*     trajectory_free(traj); traj = NULL; */

    policy_free(pol); 

    
    c3sc_destroy(sc);
    return 0;
}

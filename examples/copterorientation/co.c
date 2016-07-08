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

    double lb[6] = {-M_PI/2.0, -M_PI/3.0, -M_PI/3.0, -2.0, -2.0, -0.5};
    double ub[6] = {M_PI/2.0, M_PI/3.0, M_PI/3.0, 2.0, 2.0, 0.5};
    if ((x[0] < lb[0]) || (x[0] > ub[0]) ){
        printf(" state[0] out of bounds, %G ,t = %G \n",x[0],t);
        return 1;
    }
    if ((x[1] < lb[1]) || (x[1] > ub[1]) ){
        printf(" state[1] out of bounds, %G, t = %G\n",x[1],t);
        return 1;
    }
    if ((x[2] < lb[2]) || (x[2] > ub[2]) ){
        printf(" state[2] out of bounds, %G, t = %G\n",x[2],t);
        dprint(6,(double *)x);
        return 1;
    }
    if ((x[3] < lb[3]) || (x[3] > ub[3]) ){
        printf(" state[3] out of bounds, %G, t = %G\n",x[3],t);
        return 1;
    }
    if ((x[4] < lb[4]) || (x[4] > ub[4]) ){
        printf(" state[4] out of bounds, %G, t = %G\n",x[4],t);
        return 1;
    }
    if ((x[5] < lb[5]-1e-13) || (x[5] > ub[5]+1e-13) ){
        printf(" state[5] out of bounds, %3.15G < %3.15G\n",x[5],lb[5]);
        return 1;
    }

    double lbu[3] = {-0.1, -0.1, -0.01};
    double ubu[3] = {0.1, 0.1, 0.01};
    if ((u[0] < lbu[0]) || (u[0] > ubu[0]) ){
        return 1;
    }
    if ((u[1] < lbu[1]) || (u[1] > ubu[1]) ){
        return 1;
    }
    if ((u[2] < lbu[2]) || (u[2] > ubu[2]) ){
        return 1;
    }
    
    double cphi = cos(x[2]);
    double sphi = sin(x[2]);
    
    double cth = cos(x[1]);
    double sth = sin(x[1]);

    double Jx = 5.78e-3;
    double Jy = 5.84e-3;
    double Jz = 10.35e-3;
       
    double qr = x[4]*x[5];
    double pr = x[3]*x[5];
    double pq = x[3]*x[4];
    
    out[0] = x[5]*cphi/cth + x[4] * sphi/cth;
    out[1] = x[4]*cphi - x[5]*sphi;
    out[2] = x[3] + x[5]*cphi*sth/cth + x[4]*sth*sphi/cth;
    out[3] = (u[0] + Jy * qr - Jz * qr)/Jx;
    out[4] = (u[1] - Jx * pr + Jz * pr)/Jy;
    out[5] = (u[2] + Jx * pq - Jy * pq)/Jz;

    if (jac != NULL){
        //df1/du
        jac[0] = 0.0;    jac[6]  = 0.0;    jac[12] = 0.0;
        jac[1] = 0.0;    jac[7]  = 0.0;    jac[13] = 0.0;
        jac[2] = 0.0;    jac[8]  = 0.0;    jac[14] = 0.0;
        jac[3] = 1.0/Jx; jac[9]  = 0.0;    jac[15] = 0.0;
        jac[4] = 0.0;    jac[10] = 1.0/Jy; jac[16] = 0.0;
        jac[5] = 0.0;    jac[11] = 0.0;    jac[17] = 1.0/Jz;
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

    for (size_t ii = 0; ii < 36; ii++){
        out[ii] = 0.0;
    }
    double vpos = 1e-1;
    double vspeed = 1e-1;
    double vspeed_last = 1e-1;

    out[0] = vpos;
    out[7] = vpos;
    out[14] = vpos;
    out[21] = vspeed;
    out[28] = vspeed;
    out[35] = vspeed_last;

    if (grad != NULL){
        for (size_t ii = 0; ii < 36*3; ii++){
            grad[ii] = 0.0;
        }
    }
    return 0;
}

int stagecost(double t,const double * x,const double * u, double * out, 
              double * grad)
{
    (void)(t);
    /* (void)(u); */
    /* (void)(x); */
    *out = 0.0;

    // states
    /* *out = 5.0*0.41*pow(x[0],2) + 0.91*pow(x[1],2) + 0.91 * pow(x[2],2); */
    /* *out = *out + 0.04 * pow(x[3],2) + 0.04 * pow(x[4],2) + 0.04 * pow(x[5],2); */

    // controls
//    double c = 37.18;
    double wu1 = 1e-1;
    double wu2 = 1e-1;
    double wu3 = 1e-1;
    *out = *out + wu1 * pow(u[0],2) + wu2 * pow(u[1],2) +
                                        wu3 * pow(u[2],2);
    /* *out += 1e0; */
    *out += 100.0*pow(x[0],2);
    *out += 20.0*pow(x[1],2);
    *out += 10.0*pow(x[2],2);
    *out += pow(x[3],2);
    *out += pow(x[4],2);
    *out += pow(x[5],2);

    //*out = 5.0;
    
    if (grad!= NULL){
        grad[0] = 2.0 * wu1 * u[0];
        grad[1] = 2.0 * wu2 * u[1];
        grad[2] = 2.0 * wu3 * u[2];
        /* grad[0] = 0.0; */
        /* grad[1] = 0.0; */
        /* grad[2] = 0.0; */
    }
    return 0;
}

int boundcost(double t, const double * x, double * out)
{
    /* printf("not here!!\n"); */
    (void)(t);
    (void)(x);
    *out = 0.0;
    *out = 1e2;
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

double startcost(const double * x, void * args)
{
    (void)(args);
    (void)(x);
    return 1e3;
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

    size_t dx = 6;
    size_t dw = 6;
    size_t du = 3;
    double lb[6] = {-M_PI/2.0, -M_PI/3.0, -M_PI/3.0, -2.0, -2.0, -0.5};
    double ub[6] = {M_PI/2.0, M_PI/3.0, M_PI/3.0, 2.0, 2.0, 0.5};
    size_t Narr[6] = {N, N, N, N, N, N};

    /* double uopts[8*3] =  */
    /*     { -0.1, -0.1, -0.01,  */
    /*       -0.1, -0.1,  0.01, */
    /*       -0.1,  0.1, -0.01, */
    /*       -0.1,  0.1,  0.01, */
    /*       0.1, -0.1, -0.01, */
    /*        0.1, -0.1,  0.01, */
    /*        0.1,  0.1, -0.01, */
    /*        0.1,  0.1,  0.01 */
    /*     }; */
    /* double uopts[27*3] =  */
    /*     { -0.1, -0.1, -0.01,  */
    /*       -0.1, -0.1,  0.01, */
    /*       -0.1, -0.1,  0.00, */
    /*       -0.1,  0.1, -0.01, */
    /*       -0.1,  0.1,  0.01, */
    /*       -0.1,  0.1,  0.00, */
    /*       -0.1,  0.0, -0.01, */
    /*       -0.1,  0.0,  0.01, */
    /*       -0.1,  0.0,  0.00, */
    /*        0.1, -0.1, -0.01,  */
    /*        0.1, -0.1,  0.01, */
    /*        0.1, -0.1,  0.00, */
    /*        0.1,  0.1, -0.01, */
    /*        0.1,  0.1,  0.01, */
    /*        0.1,  0.1,  0.00, */
    /*        0.1,  0.0, -0.01, */
    /*        0.1,  0.0,  0.01, */
    /*        0.1,  0.0,  0.00, */
    /*        0.0, -0.1, -0.01,  */
    /*        0.0, -0.1,  0.01, */
    /*        0.0, -0.1,  0.00, */
    /*        0.0,  0.1, -0.01, */
    /*        0.0,  0.1,  0.01, */
    /*        0.0,  0.1,  0.00, */
    /*        0.0,  0.0, -0.01, */
    /*        0.0,  0.0,  0.01, */
    /*        0.0,  0.0,  0.00, */
    /*     }; */

    /* struct c3Opt * opt = c3opt_alloc(BRUTEFORCE,du); */
    /* c3opt_set_brute_force_vals(opt,27,uopts); */

    double lbu[3] = {-0.1, -0.1, -0.01};
    double ubu[3] = {0.1, 0.1, 0.01};
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    c3opt_set_absxtol(opt,1e-15);
    c3opt_set_relftol(opt,1e-15);
    c3opt_set_gtol(opt,1e-15);
    
    c3opt_ls_set_maxiter(opt,200);
    c3opt_ls_set_alpha(opt,0.1);
    c3opt_ls_set_beta(opt,0.2);

    c3opt_set_verbose(opt,0);

    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-5);
    approx_args_set_round_tol(aargs,1e-10);
    approx_args_set_kickrank(aargs,15);
    approx_args_set_maxrank(aargs,10);
    approx_args_set_startrank(aargs,10);

    double beta = 1.0;

    /* double * xt = linspace(lb[4],ub[4],N); */
    /* dprint(N,xt); */
    /* free(xt); */
    //exit(1);
    // setup problem
    c3sc sc = c3sc_create(IH,dx,du,dw);
    c3sc_set_state_bounds(sc,lb,ub);
    printf("bounds = \n");
    for (size_t ii = 0; ii < 6; ii++){
        c3sc_set_external_boundary(sc,ii,"reflect");
        printf("(%3.14G,%3.14G)\n",lb[ii],ub[ii]);
        /* double * xt = linspace(lb[4],ub[4],N); */
        /* dprint(N,xt); */
    }
    /* double center[6] = {0.0,0.0,0.0,0.0,0.0}; */
    /* double width[6] = {ub[0]-lb[0],ub[1]-lb[1],ub[2]-lb[2],0.4,0.4,0.2}; */
    /* double width[6] = {0.2, 0.2, 0.2, 0.25, 0.25, 0.2}; */
    /* dprint(6,width); */
    /* c3sc_add_obstacle(sc,center,width); */
    c3sc_add_dynamics(sc,f1,NULL,s1,NULL);
    c3sc_init_mca(sc,Narr);
    c3sc_attach_opt(sc,opt);
    c3sc_init_dp(sc,beta,stagecost,boundcost,ocost);
    int load_success = c3sc_cost_load(sc,"saved_cost.dat");
    if (load_success != 0){
        c3sc_cost_approx(sc,startcost,NULL,0,aargs);
    }

    double solve_tol = 1e-5;
    size_t npol = 5;

    struct C3SCDiagnostic * diag = c3sc_diagnostic_init();
    char filename[256];
    FILE *fp;//, *fp2;

    printf("\n\n\n\n\n\n\n");
    printf("Start Solver Iterations\n");
    printf("\n\n\n\n\n\n\n\n");
    for (size_t ii = 0; ii < niter; ii++){

        if (ii > 2){
            c3sc_pol_solve(sc,npol,solve_tol,verbose,aargs);
        }
        double diff = c3sc_iter_vi(sc,verbose,aargs,diag);

        struct Cost * cost = c3sc_get_cost(sc);
        int saved = cost_save(cost,"saved_cost.dat");
        assert (saved == 0);

        size_t * ranks = cost_get_ranks(cost);
        double normval = cost_norm2(cost);
        if (verbose != 0){
            printf("ii=%zu diff=%G,norm=%G ranks=",ii,diff,normval);
            iprint_sz(7,ranks);
        }
        sprintf(filename,"%s/%s.dat",dirout,"diagnostic");
        int dres = c3sc_diagnostic_save(diag,filename,4);
        assert (dres == 0);
        if (diff < 1e-2){
            break;
        }
    }
    /* exit(1); */
    struct Cost * cost = c3sc_get_cost(sc);
    int saved = cost_save(cost,"saved_cost.dat");
    assert (saved == 0);
    /* exit(1); */

    sprintf(filename,"%s/%s.dat",dirout,"diagnostic");
    int dres = c3sc_diagnostic_save(diag,filename,4);
    assert (dres == 0);

    struct ImplicitPolicy * pol = c3sc_create_implicit_policy(sc);
//    printf("created policy\n");
    /* char odename[256] = "rkf45"; */
    char odename[256] = "forward-euler";
    struct Integrator * ode_sys =
        integrator_create_controlled(6,3,f1,NULL,implicit_policy_controller,pol);
    integrator_set_type(ode_sys,odename);
    /* integrator_set_adaptive_opts(ode_sys,1e-5,1e-2,1e-7); */
    integrator_set_dt(ode_sys,1e-3);
    integrator_set_verbose(ode_sys,0);
//    printf("initialized integrator\n");

    double time = 0.0;
    /* double state[6] = {0.3, 0.2, 0.1, 0.1, 0.1, 0.02}; */
    double state[6] = {0.3, 0.2, 0.1, 0.0, 0.0, 0.0};
//    double state[6] = {0.1, 0.6, -0.2, 0.4, 0.2, 0.04};
    double con[3] = {0.0, 0.0, 0.0};

    struct Trajectory * traj = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj,6,3,time,state,con);
    printf("initialized trajectory\n");

    double final_time = 1.0;
    double dt = 1e-2;
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
    fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);
    
    printf("cost ranks are ");
    size_t * ranks = cost_get_ranks(cost);
    iprint_sz(dx+1,ranks);

    c3sc_destroy(sc);
    return 0;
}

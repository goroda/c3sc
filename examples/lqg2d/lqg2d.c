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
            " -a --lbu       Lower bound of control (default -1.0)\n"
            " -b --ubu       Upper bound of control (default 1.0)\n"
            " -r --diff1     Size of diffusion for x (default 1.0)\n"
            " -f --diff2     Size of diffusion for y (default 1.0)\n"
            " -s --steps     Number of iterations (default 100)\n"
            " -v --verbose   Output words (default 0)\n"
            "                1 - output main file stuff\n"
            "                >1 - also output approximation info\n"
        );
    exit (exit_code);
}

void print_cost(FILE * fp2, struct Cost * cost, size_t N1, size_t N2, double * lb, double * ub)
{

    fprintf(fp2,"x y f\n");
    double * xtest = linspace(lb[0],ub[0],N1);
    double * ytest = linspace(lb[1],ub[1],N2);

    double pt3[2];
    double v2;
    for (size_t zz = 0; zz < N1; zz++){
        for (size_t jj = 0; jj < N2; jj++){
            pt3[0] = xtest[zz]; pt3[1] = ytest[jj];
            cost_eval(cost,0.0,pt3,&v2);
            fprintf(fp2, "%3.5f %3.5f %3.5f \n",
                    xtest[zz],ytest[jj],v2);
        }
        fprintf(fp2,"\n");
    }
    free(xtest); xtest = NULL;
    free(ytest); ytest = NULL;
}

/* void print_policy(FILE * fp2, struct Policy *pol, size_t N1, size_t N2, */
/*                   double * lb, double * ub) */
/* { */

/*     fprintf(fp2,"x y u1 u2\n"); */
/*     double * xtest = linspace(lb[0],ub[0],N1); */
/*     double * ytest = linspace(lb[1],ub[1],N2); */

/*     double pt3[2]; */
/*     for (size_t zz = 0; zz < N1; zz++){ */
/*         for (size_t jj = 0; jj < N2; jj++){ */

/*             pt3[0] = xtest[zz]; pt3[1] = ytest[jj]; */
/*             /\* printf("\n\n\n\n\n\n"); *\/ */
/*             /\* printf("pt = "); dprint(2,pt3); *\/ */
/*             struct Control * u = NULL; */
/*             policy_eval(pol,0.0,pt3,&u); */
/* //            assert (res == 0); */
/*             /\* printf("done computing policy\n"); *\/ */
/*             fprintf(fp2, "%3.5f %3.5f %3.5f %3.5f\n", */
/*                     xtest[zz],ytest[jj],u->u[0],0.0);//u->u[1]); */
/*             control_free(u); */
/*         } */
/*         fprintf(fp2,"\n"); */
/*     } */
/*     free(xtest); xtest = NULL; */
/*     free(ytest); ytest = NULL; */
/* } */

void print_policy_implict(FILE * fp2, struct ImplicitPolicy * pol,
                          size_t N1, size_t N2,
                          double * lb, double * ub)
{

    fprintf(fp2,"x y u1 u2\n");
    double * xtest = linspace(lb[0],ub[0],N1);
    double * ytest = linspace(lb[1],ub[1],N2);

    double pt3[2];
//    double u[2];
    double u[1];
    for (size_t zz = 0; zz < N1; zz++){
        for (size_t jj = 0; jj < N2; jj++){

            pt3[0] = xtest[zz]; pt3[1] = ytest[jj];
            /* printf("\n\n\n\n\n\n"); */
            /* printf("pt = "); dprint(2,pt3); */

            int res = implicit_policy_eval(pol,0.0,pt3,u);
            assert (res >-1);
            /* printf("done computing policy\n"); */
            /* fprintf(fp2, "%3.5f %3.5f %3.5f %3.5f\n", */
            /*         xtest[zz],ytest[jj],u[0],u[1]); */
            fprintf(fp2, "%3.5f %3.5f %3.5f\n",
                    xtest[zz],ytest[jj],u[0]);

        }
        fprintf(fp2,"\n");
    }
    free(xtest); xtest = NULL;
    free(ytest); ytest = NULL;
}

int f1(double t, const double * x, const double * u, double * out,
       double * jac, void * args)
{
    (void)(t);
    (void)(args);
    
    out[0] = x[1];
    out[1] = u[0];

    if (jac != NULL){
        //df1/du
        jac[0] = 0.0;
        jac[1] = 1.0;
    }

    /* out[0] = 3*x[1] + u[1]; */
    /* out[1] = u[0]; */

    /* if (jac != NULL){ */
    /*     jac[0] = 0.0; jac[2] = 1.0; */
    /*     jac[1] = 1.0; jac[3] = 0.0; */
    /* } */
    
    return 0;
}

int s1(double t,double * x,double * u,double * out, double * grad,
       void * args)
{
    (void)(t);
    (void)(x);
    (void)(u);

    double * val = args;
    
    out[0] = val[0];
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = val[1];

    if (grad != NULL){
        /* for (size_t ii = 0; ii < 2*2*2;ii++){ */
        /*     grad[ii] = 0.0; */
        /* } */
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
    }
    return 0;
}

int stagecost(double t, double * x, double * u, double * out, 
              double * grad)
{
    (void)(t);
    *out = 0.0;

    /* *out += pow(x[0],2) + 4.0*pow(x[1],2) + pow(u[0],2) + 2.0*pow(u[1],2); */

    /* if (grad!= NULL){ */
    /*     grad[0] = 2 * u[0]; */
    /*     grad[1] = 4 * u[1]; */
    /* } */

    *out += pow(x[0],2) + pow(x[1],2) + pow(u[0],2);
    if (grad!= NULL){
        grad[0] = 2 * u[0];
    }
    return 0;
}

int boundcost(double t, double * x, double * out)
{

    (void)(t);
    (void)(x);
    *out = 0.0;
    *out = 100.0;
    return 0;
}

int ocost(double * x,double * out)
{
    (void)(x);
    //printf("got ocost!!\n");
    *out = 0.0;
    return 0;
}

double startcost(double * x, void * args)
{
    (void)(args);
    (void)(x);
    if ((fabs(x[0]) <= 2e-1) && (fabs(x[1]) < 2e-1)){
        return 0.2;
    }
    else{
        return 0.2;
    }
}

int main(int argc, char * argv[])
{
    int next_option;
    const char * const short_options = "hd:n:a:b:r:f:s:v:";
    const struct option long_options[] = {
        { "help"     , 0, NULL, 'h' },
        { "directory", 1, NULL, 'd' },
        { "nodes"    , 1, NULL, 'n' },
        { "lbu"      , 1, NULL, 'a' },
        { "ubu"      , 1, NULL, 'b' },
        { "diff1"    , 1, NULL, 'r' },
        { "diff2"    , 1, NULL, 'f' },
        { "steps"    , 1, NULL, 's' },
        { "verbose"  , 1, NULL, 'v' },
        { NULL       , 0, NULL, 0   }
    };
    program_name = argv[0];

    char * dirout = ".";
    int verbose = 0;
    size_t N = 20;
    size_t niter = 100;
    double lbu[1] = {-1.0};
    double ubu[1] = {1.0};
    double ss[2] = {1.0,1.0};
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
    double lb[2] = {-2.0, -2.0};
    double ub[2] = {2.0, 2.0};
    size_t Narr[2] = {N, N};

    // optimization arguments
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    c3opt_set_absxtol(opt,1e-8);
    c3opt_set_relftol(opt,1e-8);
    c3opt_set_gtol(opt,1e-10);
    c3opt_set_verbose(opt,0);
    
    // cross approximation tolerances
    double cross_tol = 1e-10;
    double round_tol = 1e-10;
    size_t kickrank = 5;

    // setup problem
    double beta = 0.1; 
    c3sc sc = c3sc_create(IH,dx,du,dw);
    c3sc_set_state_bounds(sc,lb,ub);
    /* c3sc_set_external_boundary(sc,0,"reflect"); */
    /* c3sc_set_external_boundary(sc,1,"reflect"); */
    double center[2] = {0.0,0.0};
    double width[2] = {0.1,0.1};
    c3sc_add_obstacle(sc,center,width);
    c3sc_add_dynamics(sc,f1,NULL,s1,ss);
    c3sc_init_mca(sc,Narr);
    c3sc_attach_opt(sc,opt);
    c3sc_init_dp(sc,beta,stagecost,boundcost,ocost);

    size_t N1 = N, N2 = N;
    struct DPih * dp = c3sc_get_dp(sc);
    struct Cost * cost = dpih_get_cost(dp);
    cost_approx(cost,startcost,NULL,verbose-1,cross_tol,round_tol,kickrank);


    
    // keep track of 1) Norm 2) Diff 3) rank
    struct Trajectory * diagnostics = NULL;
    double track[3];
    size_t npol = 1000;
//    struct Cost * tcost = NULL;
    for (size_t ii = 0; ii < niter+1; ii++){

        FILE *fp2;
        char filename[256];
        sprintf(filename,"%s/%s_%zu.dat",dirout,"costfunc",ii);
        fp2 =  fopen(filename, "w");
        if (fp2 == NULL){
            fprintf(stderr, "cat: can't open %s\n", filename);
            return 0;
        }

        print_cost(fp2,cost,N1,N2,lb,ub);
        fclose(fp2);

        if (ii > 14){
            struct ImplicitPolicy * pol = c3sc_create_implicit_policy(sc);
            dpih_iter_pol_solve(dp,pol,npol,cross_tol,
                                verbose,cross_tol,
                                round_tol,kickrank);
            printf("done here!\n");
            implicit_policy_free(pol);
        }
        cost = dpih_get_cost(dp);
        printf("lets continue\n");
        
        //struct Cost * newwhcost = dpih_iter_pol(dp,verbose-1);
        struct Cost * newcost = dpih_iter_vi(dp,verbose-1,
                                             cross_tol,round_tol,kickrank);
        printf("got new cost\n");
        track[0] = cost_norm2(newcost);
        track[1] = cost_norm2_diff(newcost,cost);
        size_t * ranks = cost_get_ranks(newcost);
        track[2] = ranks[1];
        trajectory_add(&diagnostics,3,0,(double)(ii+1),track,NULL);

        cost_free(cost);
        cost = newcost;
        dpih_attach_cost(dp,cost);
        printf("attached new cost\n");

    /*     delta = dpih_pi_iter_approx(&prob,verbose); */
        if (verbose != 0){
            printf("ii=%zu diff = %G, norm=%G\n",ii,track[1],track[0]);
        }
        if (track[1] < 1e-5){
            break;
        }
    }
//    printf("here!?\n");
    FILE *fp2;
    char filename[256];
    sprintf(filename,"%s/absorb_ulb%3.2f_uub%3.2f_s1%3.2f_s2%3.2f_%s_final.dat",
            dirout,lbu[0],ubu[0],ss[0],ss[1],"cost");
    fp2 =  fopen(filename, "w");
    if (fp2 == NULL){
        fprintf(stderr, "cat: can't open %s\n", filename);
        return 0;
    }
    print_cost(fp2,cost,N1,N2,lb,ub);
    fclose(fp2);

    sprintf(filename,"%s/%s.dat",dirout,"diagnostic");
    FILE * fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(diagnostics,fp,4);
    fclose(fp);

    struct ImplicitPolicy * pol = c3sc_create_implicit_policy(sc);
    printf("created policy\n");
    char odename[256] = "rk4";
    struct Integrator * ode_sys =
        integrator_create_controlled(2,1,f1,NULL,implicit_policy_controller,pol);
    integrator_set_type(ode_sys,odename);
    integrator_set_dt(ode_sys,1e-2);
//    integrator_set_adaptive_opts(ode_sys,dtmin,dtmax,tol);
    integrator_set_verbose(ode_sys,0);
    printf("initialized integrator\n");
    // Initialize trajectories for filter and for observations
    double time = 0.0;
    double state[2] = {0.5, 0.5};
    double con[1] = {0.0};

    struct Trajectory * traj = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj,2,1,time,state,con);
    printf("initialized trajectory\n");

    double final_time = 5e0;
    double dt = 1e-2;
    int res;
    while (time < final_time){
        // printf("time = %G\n",time);
        res = trajectory_step(traj,ode_sys,dt);
        assert(res == 0);
        time = time + dt;
    }

    /* if (verbose == 1){ */
    /*     trajectory_print(traj,stdout,4); */
    /* } */

    sprintf(filename,"%s/%s.dat",dirout,"traj");
    fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);

    printf("cost ranks are ");
    size_t * ranks = cost_get_ranks(cost);
    iprint_sz(dx+1,ranks);

    integrator_destroy(ode_sys); ode_sys = NULL;
    trajectory_free(traj); traj = NULL;
    trajectory_free(diagnostics); diagnostics = NULL;
    implicit_policy_free(pol); pol = NULL;
    c3sc_destroy(sc); sc= NULL;
    return 0;
}

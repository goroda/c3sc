#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include "c3/c3.h"
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
       double * grad, void * args)
{
    (void)(t);
    (void)(args);

    double m = 0.05;
    double g = 9.81;
    double rho = 1.292;
    double S_w = 0.1;
    double S_e = 0.025;
    double I = 6e-3;
    double l = 0.35;
    double l_w = -0.03;
    double l_e = 0.04;
    /* double l_t = 0.05; */

    /* double aspectRatio = 8.0/3.0;// #Cory 08 */
    /* double chordLength = 98.0/1000.0;// #Cory 08 */
    /* double surfaceArea = .1022;// #k */

    //for ease of notation
    //Chord and distance lengths
    /* double aspectRatioElevator = 2.0; */
        
//        #Spans
    /* double b = l_w * aspectRatio; */
    /* double b_e = l_e * aspectRatioElevator; */


    double lb[7] = {-4.0, -1.0, -M_PI/2.0, -2.0*M_PI/9.0, 0.0, -5.0, -10.0};
    double ub[7] = {0.0, 1.0, M_PI/2.0, 2.0*M_PI/9.0, 7.0, 5.0, 10.0};

    if ((x[0] < lb[0]) || (x[0] > ub[0]) ){
        /* printf(" state[0] out of bounds, %G\n",x[0]); */
        /* return 1; */
    }
    if ((x[1] < lb[1]) || (x[1] > ub[1]) ){
        printf(" state[1] out of bounds, %G\n",x[1]);
        /* return 1; */
    }
    if ((x[2] < lb[2]) || (x[2] > ub[2]) ){
        /* printf(" state[2] out of bounds, %G\n",x[2]); */
        /* dprint(7,(double *)x); */
        /* return 1; */
    }
    if ((x[3] < lb[3]) || (x[3] > ub[3]) ){
        /* printf(" state[3] out of bounds, %G\n",x[3]); */
        /* return 1; */
    }
    if ((x[4] < lb[4]) || (x[4] > ub[4]) ){
        /* printf(" state[4] out of bounds, %G\n",x[4]); */
        /* return 1; */
    }
    if ((x[5] < lb[5]) || (x[5] > ub[5]) ){
        /* printf(" state[5] out of bounds, %G\n",x[5]); */
        /* return 1; */
    }
    if ((x[6] < lb[6]) || (x[6] > ub[6]) ){
        printf(" state[5] out of bounds, %G\n",x[5]);
        /* return 1; */
    }

    // \vec{x} = [x, z, \theta, \phi, \dot{x}, \dot{z}, \dot{\theta}]
    //         = [0, 1, 2,      3,    4,       5,       6]

    double c_t = cos(x[2]);
    double c_tp = cos(x[2]+x[3]);
    double c_p = cos(x[3]);
    double s_t = sin(x[2]);
    double s_tp =sin(x[2]+x[3]);

    /* double x_w[2] ={x[0] - l_w * c_t,  */
    /*                 x[1] - l_w * s_t}; */
    double dx_w[2] = {x[4] +  l_w*x[6]*s_t,
                      x[5] -  l_w*x[6]*c_t};
    double dx_w_norm_sq = dx_w[0]*dx_w[0] + dx_w[1]*dx_w[1];

    /* double x_e[2] = {x[0] -  l*c_t -  l_e*c_tp,  */
    /*                  x[1] -  l*s_t -  l_e*s_tp}; */
   
    double dx_e[2] = {x[4] +  l*x[6]*s_t + l_e*(x[6]+u[0])*s_tp,
                      x[5] -  l*x[6]*c_t -l_e*(x[6]+u[0])*c_tp};
    double dx_e_norm_sq = dx_e[0]*dx_e[0] + dx_e[1]*dx_e[1];
        
    double alpha_w = x[2]-atan2(dx_w[1],dx_w[0]);
    double alpha_e = x[2]+x[3]-atan2(dx_e[1],dx_e[0]);

    double f_w = rho * S_w * dx_w_norm_sq * sin(alpha_w);
    double f_e = rho * S_e * dx_e_norm_sq * sin(alpha_e);

    out[0] = x[4];
    out[1] = x[5];
    out[2] = x[6];
    out[3] = u[0];
    out[4] = (-f_w * s_t - f_e * s_tp) / m;
    out[5] = (f_w*c_t + f_e*c_tp - m * g)/ m;
    out[6] = (-f_w * l_w - f_e*( l*c_p+ l_e))/ I;

    if (grad != NULL){

        double ddxe_du =  l_e * s_tp;
        double ddze_du = - l_e * c_tp;
        double ddze_ddx_e_du = (dx_e[0]*ddze_du - dx_e[1]*ddxe_du) / (dx_e[0] * dx_e[0]);
        
        double dtaninv = 1.0/(1.0 + pow(dx_e[1]/dx_e[0],2.0)) * ddze_ddx_e_du;
        double        dalphae_du = -dtaninv;

        double dsae_du = cos(alpha_e)*dalphae_du;
        double dabsxe_du = 2.0*dx_e[0]*ddxe_du + 2.0*dx_e[1]*ddze_du;
        double sumdx_e2 = dx_e[0]*dx_e[0] + dx_e[1] * dx_e[1];
        double dfe_du =  rho * S_e *(dabsxe_du*sin(alpha_e) + sumdx_e2*dsae_du);
        
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = -s_tp / m * dfe_du;
        grad[5] = c_tp / m * dfe_du;
        grad[6] = -(l*c_p + l_e) / I * dfe_du;
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

    for (size_t ii = 0; ii < 7*7; ii++){
        out[ii] = 0.0;
    }
    for (size_t ii = 0; ii < 7; ii++){
        out[ii*7+ii] = 1e-9;
    }
    
    if (grad != NULL){
        for (size_t ii = 0; ii < 7*7; ii++){
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

    //*out += 10.0;

    // states
    *out += 20.0 * x[0]*x[0];
    /* *out += 100 * x[0]*x[0]; */
    *out += 50.0 * x[1]*x[1]; // original

    /* *out += 50.0 * x[1]*x[1]; */


    *out += 10.0 * x[2]*x[2]; // original
    /* *out +=  1.0 * x[2]*x[2];  */
    *out +=  1.0 * x[3]*x[3];
    *out +=  1.0 * x[4]*x[4]; // original

    /* *out += 10.0 * x[4]*x[4]; */
    //*out +=  4.0 * x[4]*x[4]; // addition

    *out +=  1.0 * x[5]*x[5];

//    *out += 5.0 * x[5]*x[5]; // addition

    *out +=  1.0 * x[6]*x[6];

    // control
    *out += 0.1 * u[0] * u[0];
    if (grad!= NULL){
        grad[0] = 0.1 * 2.0 * u[0];
    }

    return 0;
}

int boundcost(double t,const double * x, double * out)
{
    /* printf("not here!!\n"); */
    (void)(t);
    *out = 0.0;
    *out += 600.0    * x[0]*x[0]; // original
    /* *out += 400 * x[0]*x[0]; */
    /* *out += 400.0    * x[0]*x[0]; */
    *out += 400.0 * x[1]*x[1];// original
    /* *out += 600.0 * x[1] * x[1]; //addition */
    /* *out += 800.0    * x[1]*x[1]; */
    *out += 1.0/9.0  * x[2]*x[2];
    *out += 5.0 * (x[2]-M_PI/2.0) * (x[2]-M_PI/2.0);
    *out += 1.0/9.0  * x[3]*x[3];
    *out +=  1.0     * x[4]*x[4]; // original
    /* *out +=  100.0     * x[4]*x[4]; */
    *out +=  1.0     * (x[5]+1.5)*(x[5]+1.5); // original 
    /* *out +=  100.0     * (x[5]+1.5)*(x[5]+1.5); // original  */
    *out += 1.0/9.0  * (x[6]+0.5)*(x[6]+0.5);
    
    return 0;
}

int startcost(size_t N, const double * xin, double * out, void * args)
{
    (void)(args);
    for (size_t ii = 0; ii < N; ii++){
        double x[7] = {xin[ii*7+0],xin[ii*7+1],xin[ii*7+2],xin[ii*6+3],
                       xin[ii*7+4],xin[ii*7+5],xin[ii*7+6]};
        out[ii] = 0.0;
        out[ii] += 600.0    * x[0]*x[0];
        out[ii] += 400.0    * x[1]*x[1];
        out[ii] += 1.0/9.0  * x[2]*x[2];
        out[ii] += 1.0/9.0  * x[3]*x[3];
        out[ii] +=  1.0     * x[4]*x[4];
        out[ii] +=  1.0     * x[5]*x[5];
        out[ii] += 1.0/9.0  * x[6]*x[6];
    }
    return 0.0;
}

static int inobs;
int ocost(const double * x,double * out)
{
    /* dprint(6,x); */
    /* printf("here\n"); */
    (void)(x);
    inobs = 1;
    *out = 0.0;
    return 0;
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
    size_t N = 20;
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

    size_t dx = 7;
    size_t dw = 7;
    size_t du = 1;
    double lb[7] = {-4.0, -1.0, -M_PI/2.0, -2.0*M_PI/9.0, 0.0, -5.0, -10.0};
    double ub[7] = {0.0, 1.0, M_PI/2.0, 2.0*M_PI/9.0, 7.0, 5.0, 10.0};
    size_t Narr[7] = {N, N, N, N, N, N, N};
    double beta = 1.0;

    double lbu[1] = {-2.0*M_PI};
    double ubu[1] = {2.0*M_PI};
    struct c3Opt * opt = c3opt_alloc(BRUTEFORCE,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);

    size_t nopts = 40;
    double * uopts = linspace(lbu[0],ubu[0],nopts);
    c3opt_set_brute_force_vals(opt,nopts,uopts);

    /* struct c3Opt * opt = c3opt_alloc(BFGS,du); */

    /* c3opt_set_absxtol(opt,1e-10); */
    /* c3opt_set_relftol(opt,1e-10); */
    /* c3opt_set_gtol(opt,1e-30); */
    
    /* c3opt_ls_set_maxiter(opt,10); */
    /* c3opt_ls_set_alpha(opt,0.3); */
    /* c3opt_ls_set_beta(opt,0.2); */
    c3opt_set_verbose(opt,0);

    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-10);
    approx_args_set_round_tol(aargs,1e-10);
    approx_args_set_kickrank(aargs,10);
    approx_args_set_adapt(aargs,1);

    approx_args_set_startrank(aargs,5);
    approx_args_set_maxrank(aargs,10);

    // setup problem
    struct C3Control * c3c = c3control_create(dx,du,dw,lb,ub,Narr,beta);
    c3control_add_drift(c3c,f1,NULL);
    c3control_add_diff(c3c,s1,NULL);
    c3control_add_stagecost(c3c,stagecost);
    c3control_add_boundcost(c3c,boundcost);
    c3control_add_obscost(c3c,ocost);
    
    // x y aoa phi dx dy daoa
    double center[7];
    double width[7];
    center[0] = 0.0; width[0] = 0.1;
    center[1] = 0.0, width[1] = 0.1;
    center[2] = 0.0; width[2] = ub[2]-lb[2];
    center[3] = 0.0; width[3] = ub[3]-lb[3];
    center[4] = 0.0; width[4] = 0.5;
    center[5] = -2.0; width[5] = 0.5;
    center[6] = 0.0; width[6] = ub[6] - lb[6];

    c3control_add_obstacle(c3c,center,width);
    
    /* c3control_set_external_boundary(c3c,0,"reflect"); */
    /* c3control_set_external_boundary(c3c,1,"reflect"); */
    /* c3control_set_external_boundary(c3c,2,"reflect"); */
    /* c3control_set_external_boundary(c3c,3,"reflect"); */
    /* c3control_set_external_boundary(c3c,4,"reflect"); */
    /* c3control_set_external_boundary(c3c,5,"reflect"); */
    /* c3control_set_external_boundary(c3c,6,"reflect"); */

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

        inobs = 0;
        struct ValueF * next = c3control_pi_solve(c3c,maxiter_pi,abs_conv_pi,
                                                  cost,aargs,opt,verbose,&diag);
        //assert (inobs == 1);
        inobs = 0;
        valuef_destroy(cost); cost = NULL;
        cost = c3control_vi_solve(c3c,1,abs_conv_vi,next,aargs,opt,verbose,&diag);
        valuef_destroy(next); next = NULL;
        //assert (inobs == 1);

        sprintf(filename,"%s/%s.c3",dirout,"cost");
        int saved = valuef_save(cost,filename);
        assert (saved == 0);

        saved = diag_save(diag,filename_diag);
        assert (saved == 0);

        size_t * ranks = valuef_get_ranks(cost);
        if (verbose != 0){
            printf("ii=%zu ranks=",ii);
            iprint_sz(dx+1,ranks);
        }

    }

    printf("beginning simulations\n");
    c3control_add_policy_sim(c3c,cost,opt,NULL);

    char odename[256] = "forward-euler";
    struct Integrator * ode_sys =
        integrator_create_controlled(7,1,f1,NULL,c3control_controller,c3c);
    integrator_set_type(ode_sys,odename);
    integrator_set_dt(ode_sys,1e-4);
    integrator_set_verbose(ode_sys,0);

    double time = 0.0;
    /* double state[7] = {-3.5, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0}; */
    /* double state[7] = {-3.0, 0.2, 0.0, 0.0, 4.5, 0.0, 0.0}; */
    double state[7] = {-3.5, 0.0, 0.0, 0.0, 6.5, 0.0, 0.0};
    double con[1] = {0.0};
    struct Trajectory * traj = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj,7,1,time,state,con);
    printf("initialized trajectory\n");
    /* double final_time = 2.5; */
    double final_time = 1.2;
    double dt = 1e-2;
    int res;
    while (time < final_time){
        /* printf("time = %G\n",time); */
        res = trajectory_step(traj,ode_sys,dt);
        double * ls = trajectory_get_last_state(traj);
        if ( (ls[0] > 0.1) || (fabs(ls[1]) > 1.0)){
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
    
    integrator_destroy(ode_sys); ode_sys = NULL;    
    valuef_destroy(cost); cost = NULL;
    c3control_destroy(c3c); c3c = NULL;
    diag_destroy(&diag); diag = NULL;
    c3opt_free(opt); opt = NULL;
    approx_args_free(aargs); aargs = NULL;
    return 0;
}

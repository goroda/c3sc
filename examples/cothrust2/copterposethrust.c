#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

#include "c3.h"
#include "cdyn/src/simulate.h"
#include "cdyn/src/integrate.h"

#include "c3sc/src/nodeutil.h"
#include "c3sc/src/valuefunc.h"
#include "c3sc/src/bellman.h"
#include "c3sc/src/c3sc.h"

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


static int order[6] = {0,1,2,3,4,5};
/* static int order[6] = {0,4,5,1,2,3};  */

int f1(double t, const double * xin, const double * u, double * out,
       double * jac, void * args)
{
    (void)(t);
    (void)(args);

    double lbs[6] = {-3.5, -3.5, -2, -5, -5, -5};
    double ubs[6] = { 3.5,  3.5,  2,  5,  5,  5};

    /* printf("xin = "); dprint(6,xin); */

    double lb[6];
    double ub[6];
    double x[6];
    for (size_t ii = 0; ii < 6; ii++){
        lb[ii] = lbs[ii];
        ub[ii] = ubs[ii];
        x[ii] = xin[order[ii]];
    }
    /* printf("lb = "); dprint(6,lb); */
    /* printf("ub = "); dprint(6,ub); */
    /* printf("x = "); dprint(6,x); */

    if ((x[0] < lb[0]) || (x[0] > ub[0]) ){
        printf(" state[0] out of bounds, %G\n",x[0]);
        return 1;
    }
    if ((x[1] < lb[1]) || (x[1] > ub[1]) ){
        printf(" state[1] out of bounds, %G\n",x[1]);
        return 1;
    }
    if ((x[2] < lb[2]) || (x[2] > ub[2]) ){
        printf(" state[2] out of bounds, %G\n",x[2]);
        dprint(6,(double *)x);
        return 1;
    }
    if ((x[3] < lb[3]) || (x[3] > ub[3]) ){
        printf(" state[3] out of bounds, %G\n",x[3]);
        return 1;
    }
    if ((x[4] < lb[4]) || (x[4] > ub[4]) ){
        printf(" state[4] out of bounds, %G\n",x[4]);
        return 1;
    }
    if ((x[5] < lb[5]) || (x[5] > ub[5]) ){
        printf(" state[5] out of bounds, %G\n",x[5]);
        return 1;
    }

    double lbu[3] = {-1.5, -0.4, -0.4};
    double ubu[3] = {1.5, 0.4, 0.4};

    if ((u[0] < lbu[0]) || (u[0] > ubu[0]) ){
        return 1;
    }
    if ((u[1] < lbu[1]) || (u[1] > ubu[1]) ){
        return 1;
    }
    if ((u[2] < lbu[2]) || (u[2] > ubu[2]) ){
        return 1;
    }

    double m = 1.227;
    double g = 9.81;
    double mg=m*g;
    
    double cphi = cos(u[1]);
    double sphi = sin(u[1]);
    
    double cth = cos(u[2]);
    double sth = sin(u[2]);
       
    out[0] = x[3];
    out[1] = x[4];
    out[2] = x[5];
    out[3] = cphi*sth*(u[0]-mg)/m;
    out[4] = -sphi*(u[0]-mg)/m;
    out[5] = g + cth*cphi*(u[0]-mg)/m;

    if (jac != NULL){
        //df1/du
        jac[0] = 0.0;      jac[6]  = 0.0;                     jac[12] = 0.0;
        jac[1] = 0.0;      jac[7]  = 0.0;                     jac[13] = 0.0;
        jac[2] = 0.0;      jac[8]  = 0.0;                     jac[14] = 0.0;
        jac[3] = cphi*sth/m; jac[9]  = -sphi*sth*(u[0]-mg)/m; jac[15] = cphi*cth*(u[0]-mg)/m;
        jac[4] = -sphi/m;    jac[10] = -cphi*(u[0]-mg)/m;     jac[16] = 0.0;
        jac[5] = cth*cphi/m; jac[11] = -cth*sphi*(u[0]-mg)/m; jac[17] = -sth*cphi*(u[0]-mg)/m;
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

    /*double vpos = 1e-2; //orig
    double vspeed = 1e-1; //orig
    double vspeed_last = 1e-1; //orig
    */
    double vposX = 1e-1;
    double vposZ = 2e-1; //2e-2 (a)
    double vspeed = 12e-1; //15e-1 (a)

    out[0] = vposX;
    out[7] = vposX;
    out[14] = vposZ;
    out[21] = vspeed;
    out[28] = vspeed;
    out[35] = vspeed;

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
    (void)(u);
    (void)(x);
    *out = 0.0;
 
    // controls
    double wu1 = 2;//0.0;//1/400 (a);//0.5;//1/144;//1/170;0.5; //0.5. 10@n8, 10@n16,10@n32
    double wu2 = 1;//0.03; //0.5 // 0.2 (a);   //2
    double wu3 = 6;//0.2; //0.5 // 0.2 (a);	//2
    //*out = *out + 10.0 + wu1 * pow(u[0],2) + wu2 * pow(u[1],2) + wu3 * pow(u[2],2); (a)
    *out = *out + 60.0 + wu1 * pow(u[0],2) + wu2 * pow(u[1],2) + wu3 * pow(u[2],2); //30 //40

    /* *out = *out + 8.0*pow(x[3]-1.0,2.0); */
    *out = *out + 8.0*pow(x[2],2.0);
    *out = *out + 6.0*pow(x[1],2.0);
    *out = *out + 8.0*pow(x[0],2.0);
    
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

int boundcost(double t,const double * x, double * out)
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
    /* dprint(6,x); */
    /* printf("here\n"); */
    (void)(x);
    *out = 0.0;
    return 0;
}

int startcost(size_t N, const double * x, double * out,void * args)
{
    (void)(args);
    (void)(x);
    for (size_t ii = 0; ii < N; ii++){
      out[ii] = 100.0;
    }
    
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

    size_t dx = 6;
    size_t dw = 6;
    size_t du = 3;
    double lbs[6] = {-3.5, -3.5, -2, -5, -5, -5};
    double ubs[6] = {3.5, 3.5, 2, 5, 5, 5};
    size_t Narr[6] = {N, N, N, N, N, N};

    double lb[6];
    double ub[6];
    for (size_t ii = 0; ii < 6; ii++){
        lb[ii] = lbs[order[ii]];
        ub[ii] = ubs[order[ii]];
    }

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

    double lbu[3] = {-1.5, -0.4, -0.4};
    double ubu[3] = {1.5, 0.4, 0.4};

    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    c3opt_set_absxtol(opt,1e-10); //-14
    c3opt_set_relftol(opt,1e-10); //-14
    c3opt_set_gtol(opt,1e-30); //-30
    
    c3opt_ls_set_maxiter(opt,200);
    c3opt_ls_set_alpha(opt,0.1);
    c3opt_ls_set_beta(opt,0.2);

    c3opt_set_verbose(opt,0);


    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-5); //decreased to e-7 at some point
    approx_args_set_round_tol(aargs,1e-5); //decrtease to e-6 at some point
    approx_args_set_kickrank(aargs,5);

    approx_args_set_startrank(aargs,10);
    approx_args_set_maxrank(aargs,10); //5 better since matches guess (?), 10 gives nice results


    double beta = 1.0;

    /* double * xt = linspace(lb[4],ub[4],N); */
    /* dprint(N,xt); */
    /* free(xt); */
    //exit(1);

    // setup problem
    struct C3Control * c3c = c3control_create(dx,du,dw,lb,ub,Narr,beta);
    c3control_add_drift(c3c,f1,NULL);
    c3control_add_diff(c3c,s1,NULL);
    c3control_add_stagecost(c3c,stagecost);
    c3control_add_boundcost(c3c,boundcost);
    c3control_add_obscost(c3c,ocost);

    printf("bounds = \n");
    for (size_t ii = 0; ii < 6; ii++){
        c3control_set_external_boundary(c3c,ii,"reflect");
        printf("(%3.14G,%3.14G)\n",lb[ii],ub[ii]);
    }
    double center[6] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    //double width[6] = {2,2,1,ub[3]-lb[3],ub[4]-lb[4],1};
    double width[6] =  {0.4, 0.4, 0.4, 0.5 , 0.3, 0.3};

    double c[6];
    double w[6];
    for (size_t ii = 0; ii < 6; ii++){
        c[order[ii]] = center[ii];
        w[order[ii]] = width[ii];
    }

    dprint(6,width);
    c3control_add_obstacle(c3c,c,w);
    

    char filename[256];
    sprintf(filename,"%s/%s.c3",dirout,"saved_cost");
    double ** xgrid = c3control_get_xgrid(c3c);
    struct ValueF * cost = valuef_load(filename,Narr,xgrid);
    if (cost == NULL){
        cost = c3control_init_value(c3c,startcost,NULL,aargs,0);
    }


    size_t maxiter_vi = niter+1;
    double abs_conv_vi = 1e-3;
    size_t npol = 10;
    size_t maxiter_pi = npol;
    double abs_conv_pi = 1e-2;
    struct Diag * diag = NULL;
    char filename_diag[256];
    sprintf(filename_diag,"%s/%s.dat",dirout,"diagnostic");
    if (niter!=0)
    {
    for (size_t ii = 0; ii < maxiter_vi; ii++){
        struct ValueF * next = c3control_pi_solve(c3c,maxiter_pi,abs_conv_pi,
                                                  cost,aargs,opt,verbose,&diag);

        valuef_destroy(cost); cost = NULL;
        cost = c3control_vi_solve(c3c,1,abs_conv_vi,next,aargs,opt,verbose,&diag);
        valuef_destroy(next); next = NULL;

        sprintf(filename,"%s/%s.c3",dirout,"saved_cost");
        int saved = valuef_save(cost,filename);
        assert (saved == 0);

        sprintf(filename,"%s/%s_%zu.c3",dirout,"saved_cost",ii);
        saved = valuef_save(cost,filename);
        assert (saved == 0);

        saved = diag_save(diag,filename_diag);
        assert (saved == 0);

        size_t * ranks = valuef_get_ranks(cost);
        double normval = valuef_norm(cost);
        if (verbose != 0){
            printf("ii=%zu norm=%G ranks=",ii,normval);
            iprint_sz(7,ranks);
        }

    }
    }
    c3control_add_policy_sim(c3c,cost,opt,NULL);

//    printf("created policy\n");
    /* char odename[256] = "rkf45"; */
    /* char odename[256] = "rk4"; */
    char odename[256] = "forward-euler";
    struct Integrator * ode_sys =
        integrator_create_controlled(6,3,f1,NULL,c3control_controller,c3c);
    integrator_set_type(ode_sys,odename);
    /* integrator_set_adaptive_opts(ode_sys,1e-5,1e-2,1e-7); */
    integrator_set_dt(ode_sys,1e-4);
    integrator_set_verbose(ode_sys,0);
//    printf("initialized integrator\n");

    double time = 0.0;
	//double state[6] = {-2.5, -0.5, 0.0, 0.0, 0.0, 0.0};
	//double state[6] = {-3.0, -0.5, 0.0, 0.0, 0.0, 0.0};
	//double state[6] = {-2.5, -0.5, -0.5, 0.0, 0.0, 0.0};

	/* double state[6] = {-2.0, -0.02, 0.02, 0.02, -0.02, 0.01}; */
        //double state[6] = {0.3, 0.3, 0.3, 0.0, 0.0, 0.0};
        //double state[6] = {0.1, 0.1, 0.1, 0.0, 0.0, 0.0};
	//double state[6] = {2.5, 0.5, 0.5, 0.0, 0.0, 0.0};
        //double state[6] = {-1.0, -0.02, 0.02, 0.02, -0.02, 0.01}; 
        //double state[6] = {-1.0, -0.5, 1.5, 0.02, -0.02, 0.01};
        //double state[6] = {-1.0, -0.5, -1, 0.02, -0.02, 0.01};
        //double state[6] = {-2.0, 0.005, 0.005, -0.005, -0.02, -0.0005};
	//double state[6] = {-1.0, 0.1, 1.0, 0.0, 0.0, 0.0};

	//double state[6] = {-1.0,2.0,1.0,-0.1,0.5,1.2};
	double state[6] = {-3.0, -1.0, -0.5, 0.0, 0.0, 0.0}; 
	// double state[6] = {-3.0, -0.5, 0.5, 0.0, 0.0, 0.0}; 
        
    // these three work
    
    double con[3] = {0.0, 0.0, 0.0};

    double ss[6];
    for (size_t ii = 0; ii < 6; ii++){
        ss[order[ii]] = state[ii];
    }

    struct Trajectory * traj = NULL;
    printf("add trajectory\n");
    trajectory_add(&traj,6,3,time,ss,con);
    printf("initialized trajectory\n");

    double final_time = 1.9e0;
    double dt = 1e-3;
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
    
    printf("cost ranks are ");
    size_t * ranks = valuef_get_ranks(cost);
    iprint_sz(dx+1,ranks);

    integrator_destroy(ode_sys); ode_sys = NULL;
    trajectory_free(traj); traj = NULL;


    valuef_destroy(cost); cost = NULL;
    c3control_destroy(c3c); c3c = NULL;
    diag_destroy(&diag); diag = NULL;
    c3opt_free(opt); opt = NULL;
    approx_args_free(aargs); aargs = NULL;

    return 0;
}

//optim log
/*
-y cheaper
-ychaper
*/

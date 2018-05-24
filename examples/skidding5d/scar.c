#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <getopt.h>

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
            " -s --steps     Number of iterations (default 100)\n"
            " -v --verbose   Output words (default 0)\n"
            "                1 - output main file stuff\n"
            "                >1 - also output approximation info\n"
        );
    exit (exit_code);
}

static int order[5] = { 0, 1, 2, 3, 4}; // x y orientation omega v
/* static int order[5] = { 0, 4, 3, 2, 1}; // x y orientation omega v */

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


    /* int order[4] = { 0, 1, 2, 3}; // x y o v */
    /* double posx = x[order[0]]; */
    /* double posy = x[order[1]]; */

    double orient = x[order[2]];
    double angvel = x[order[3]];
    double speed = x[order[4]];

    double steering = u[0];
    
    double m = 1460.0;
    double cf = 17000.0;
    double ct = 20000.0;
    double a = 1.2;
    double b = 1.5;
    double In = 2170.0;
    double s = 27.0;

    double co = cos(orient);
    double so = sin(orient);

    // forces
    double ff = cf * ( (speed + a * angvel)/s + steering);
    double ft = ct * (speed - b * angvel)/s;

    // derivative of ff
    double dff = cf;

    out[order[0]] = s * co - speed * so;
    out[order[1]] = s * so + speed * co;
    out[order[2]] = angvel;
    out[order[3]] = (a * ff - b * ft)/ In;
    out[order[4]] = -s * angvel + (ff + ft)/m;

    if (jac != NULL){
        //df1/du
        jac[order[0]] = 0.0;    
        jac[order[1]] = 0.0;    
        jac[order[2]] = 0.0;
        jac[order[3]] = a * dff / In;
        jac[order[4]] = dff / m;

        /* jac[order[0]] = 0.0; */
        /* jac[order[1]] = 0.0; */
        /* jac[order[2]] = 0.0;    */
        /* jac[order[3]] = alpha; */
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

    for (size_t ii = 0; ii < 5*5; ii++){
        out[ii] = 0.0;
    }
    double vpos = 1e-5;
    double vorient = 1e-5;
    double vspeed = 1e-5;

    out[0] = vpos;
    out[6] = vpos;
    out[12] = vorient;
    out[28] = vspeed;
    out[24] = vspeed;

    // note 1 control!!
    if (grad != NULL){
        for (size_t ii = 0; ii < 25*1; ii++){
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
    /* *out = 10.0 + 10.0*pow(x[order[0]],2) + 10.0*pow(x[order[1]],2);// + 5.0*pow(x[order[4]],2); */
    *out = 1.0 + 0.02 * pow(x[order[0]],2) + 0.02 * pow(x[order[1]],2);
    *out = *out + pow(x[order[3]],2) + pow(x[order[4]],2);
    (void)(u);
    
    if (grad!= NULL){
        grad[0] = 0.0;
        /* grad[1] = 0.0; */
    }
    return 0;
}

int boundcost(double t, const double * x, double * out)
{
    /* printf("not here!!\n"); */
    (void)(t);
    *out = 0.0;
    *out = 0.1 * pow(x[order[0]],2) + 0.1 * pow(x[order[1]],2);
    *out = *out + 0.1 * pow(x[order[3]],2) + 0.1*pow(x[order[4]],2);;
    return 0;
}

int ocost(const double * x,double * out)
{
    /* printf("absorbed "); */
    /* dprint(5,x); */
    (void)(x);
    *out = 0.0;
    return 0;
}

int startcost(size_t N, const double * x, double * out, void * args)
{
    (void)(args);
    (void)(x);
    for (size_t ii = 0; ii < N; ii++){
        /* out[ii] = pow(x[ii*4+(size_t)order[0]],2) + pow(x[ii*4+(size_t)order[1]],2); */
                 /* + pow(x[ii*4+2],2) + pow(x[ii*4+3],2); */
        out[ii] = 0.1;
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
    y[4] = x[4];
    /* int trans = 0; */
    if ((x[order[2]] < M_PI ) && (x[order[2]] > -M_PI)){
        y[order[2]] = x[order[2]];
    }
    else if (x[order[2]] > M_PI){
        /* trans = 1; */
        while (y[order[2]] > 2 * M_PI){
            y[order[2]] -= 2.0*M_PI;
        }
        if (y[order[2]] > M_PI){
            y[order[2]] = y[order[2]] - 2.0*M_PI;
        }
    }
    else if (x[order[2]] < -M_PI){
        /* trans = 1; */
        while (y[order[2]] < - 2.0* M_PI){
            y[order[2]] += 2.0 * M_PI;
        }
        if (y[order[2]] < -M_PI){
            y[order[2]] += 2.0 * M_PI;
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

    size_t dx = 5;
    size_t dw = 5;
    size_t du = 1;

    // lower and upper bounds of standard ordering
    double lbs[5] = {-500.0, -500.0, -M_PI, -0.5,  -10.0};
    double ubs[5] = {500.0, 500.0, M_PI, 0.5,  10.0};
    
    double lb[5] = {lbs[order[0]],lbs[order[1]],lbs[order[2]],lbs[order[3]],lbs[order[4]]};
    double ub[5] = {ubs[order[0]],ubs[order[1]],ubs[order[2]],ubs[order[3]],ubs[order[4]]};
    size_t Narr[5] = {N,N,N,N,N};
    double beta = 1.0;

    double lbu[1] = {-5.0*M_PI/180.0};//, -1.0};
    double ubu[1] = {5.0*M_PI/180.0};//, 1.0};

    struct c3Opt * opt = c3opt_alloc(BRUTEFORCE,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    size_t dopts = 20;
    double * uopts = linspace(lbu[0],ubu[0],dopts);
    c3opt_set_brute_force_vals(opt,dopts,uopts);

    /* c3opt_set_absxtol(opt,1e-8); */
    /* c3opt_set_relftol(opt,1e-8); */
    /* c3opt_set_gtol(opt,1e-13); */
    /* c3opt_ls_set_maxiter(opt,10); */
    /* c3opt_ls_set_alpha(opt,0.1); */
    /* c3opt_ls_set_beta(opt,0.2); */
    /* c3opt_set_verbose(opt,0); */

    // cross approximation tolerances
    struct ApproxArgs * aargs = approx_args_init();
    approx_args_set_cross_tol(aargs,1e-8);
    approx_args_set_round_tol(aargs,1e-5);
    approx_args_set_adapt(aargs,1);
    approx_args_set_kickrank(aargs,10);
    approx_args_set_startrank(aargs,5);
    approx_args_set_maxrank(aargs,15);


    // setup problem
    struct C3Control * c3c = c3control_create(dx,du,dw,lb,ub,Narr,beta);
    c3control_add_drift(c3c,f1,NULL);
    c3control_add_diff(c3c,s1,NULL);
    c3control_add_stagecost(c3c,stagecost);
    c3control_add_boundcost(c3c,boundcost);
    c3control_add_obscost(c3c,ocost);

    c3control_set_external_boundary(c3c,(size_t)order[0],"reflect");
    c3control_set_external_boundary(c3c,(size_t)order[1],"reflect");
    c3control_set_external_boundary(c3c,(size_t)order[2],"periodic");
    /* c3control_set_external_boundary(c3c,(size_t)order[3],"reflect"); */
    /* c3control_set_external_boundary(c3c,(size_t)order[4],"reflect"); */

    // possible obstacle
    double w = 40.0;
    double center[5] = {0.0, 0.0, 0.0, 0.0,0.0};
    center[order[3]] = (ubs[order[3]]+lbs[order[3]])/2.0;
    center[order[4]] = (ubs[order[4]]+lbs[order[4]])/2.0;
    
    double width[5];
    width[order[0]] = w;
    width[order[1]] = w;
    width[order[2]] = (ubs[order[2]] - lbs[order[2]]);
    width[order[3]] = (ubs[order[3]] - lbs[order[3]]);
    width[order[4]] = (ubs[order[4]] - lbs[order[4]]);
    
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
    size_t maxiter_pi = 40;
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
            iprint_sz(dx+1,ranks);
        }

    }

    printf("beginning simulations\n");
    c3control_add_policy_sim(c3c,cost,opt,state_transform);

/* //    printf("created policy\n"); */
    /* char odename[256] = "rk4"; */
    char odename[256] = "forward-euler";
    struct Integrator * ode_sys =
        integrator_create_controlled(dx,du,f1,NULL,c3control_controller,c3c);
    integrator_set_type(ode_sys,odename);
    integrator_set_dt(ode_sys,1e-2);
    integrator_set_verbose(ode_sys,0);
//    printf("initialized integrator\n");



    /* double state_st[5] = {30.0, 3.0, -M_PI / 2.0, 0.0, 3.0}; */
    /* double state_st[5] = {30.0, 30.0, - 3*M_PI / 4.0, 0.0, 0.0}; */

/////////////////////////////////////////////////////////////////////////////
    double state_st[5] = {30.0, 40.0, 3.0*M_PI/4.0, 0.0, 0.0}; // doesnt work
    double con[1] = {0.0};
    double state[5] = {state_st[order[0]], state_st[order[1]], state_st[order[2]],
                       state_st[order[3]], state_st[order[4]]};

    double time = 0.0;    
    struct Trajectory * traj = NULL;
    trajectory_add(&traj,dx,du,time,state,con);
    double dt = 1e-2;
    double final_time = 4e0;
    int res;
    while (time < final_time){
        res = trajectory_step(traj,ode_sys,dt);
        double * ls = trajectory_get_last_state(traj);
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
    sprintf(filename,"%s/%s.dat",dirout,"traj");
    FILE * fp = fopen(filename,"w");
    assert (fp != NULL);
    trajectory_print(traj,fp,4);
    fclose(fp);
    trajectory_free(traj); traj = NULL;

/////////////////////////////////////////////////////////////////////////////
    /* state_st[0] = 20.0; */
    /* state_st[1] = 0.0; */
    /* state_st[2] = -3.0*M_PI/4.0; */
    /* state_st[3] = 0.0; */
    /* state_st[4] = 0.0; */
    /* for (size_t ii = 0; ii < dx; ii++){ */
    /*     state[ii] = state_st[order[ii]]; */
    /* } */
    /* con[0] = 0; */

    /* traj = NULL; */
    /* time = 0.0; */
    /* trajectory_add(&traj,5,1,time,state,con); */
    /* while (time < final_time){ */
    /*     res = trajectory_step(traj,ode_sys,dt); */
    /*     double * ls = trajectory_get_last_state(traj); */
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
    /* trajectory_print(traj,fp,4); */
    /* fclose(fp); */
    /* trajectory_free(traj); traj = NULL; */

/////////////////////////////////////////////////////////////////////////////

/*     state_st[0] = -30.0; */
/*     state_st[1] = 20.0; */
/*     state_st[2] = -M_PI/3.0; */
/*     state_st[3] = 0.0; */
/*     state_st[4] = 0.0; */
/*     for (size_t ii = 0; ii < dx; ii++){ */
/*         state[ii] = state_st[order[ii]]; */
/*     } */
/*     con[0] = 0; */

/*     traj = NULL; */
/*     time = 0.0; */
/*     trajectory_add(&traj,5,1,time,state,con); */
/*     while (time < final_time){ */
/*         res = trajectory_step(traj,ode_sys,dt); */
/*         double * ls = trajectory_get_last_state(traj); */
/*         if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){ */
/*             break; */
/*         } */
/*         if (res != 0){ */
/*             break; */
/*         } */
/*         time = time + dt; */
/*     } */
/*     if (verbose == 1){ */
/*         trajectory_print(traj,stdout,4); */
/*     } */
/*     sprintf(filename,"%s/%s.dat",dirout,"traj3"); */
/*     fp = fopen(filename,"w"); */
/*     assert (fp != NULL); */
/*     trajectory_print(traj,fp,4); */
/*     fclose(fp); */
/*     trajectory_free(traj); traj = NULL; */


/* ///////////////////////////////////////////////////////////////////////////// */

/*     state_st[0] = 10.0; */
/*     state_st[1] = 20.0; */
/*     state_st[2] = M_PI; */
/*     state_st[3] = 0.0; */
/*     state_st[4] = 0.0; */
/*     for (size_t ii = 0; ii < dx; ii++){ */
/*         state[ii] = state_st[order[ii]]; */
/*     } */
/*     con[0] = 0; */

/*     traj = NULL; */
/*     time = 0.0; */
/*     trajectory_add(&traj,5,1,time,state,con); */
/*     while (time < final_time){ */
/*         res = trajectory_step(traj,ode_sys,dt); */
/*         double * ls = trajectory_get_last_state(traj); */
/*         if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){ */
/*             break; */
/*         } */
/*         if (res != 0){ */
/*             break; */
/*         } */
/*         time = time + dt; */
/*     } */
/*     if (verbose == 1){ */
/*         trajectory_print(traj,stdout,4); */
/*     } */
/*     sprintf(filename,"%s/%s.dat",dirout,"traj4"); */
/*     fp = fopen(filename,"w"); */
/*     assert (fp != NULL); */
/*     trajectory_print(traj,fp,4); */
/*     fclose(fp); */
/*     trajectory_free(traj); traj = NULL; */


/* ///////////////////////////////////////////////////////////////////////////// */
/*     state_st[0] = 30.0; */
/*     state_st[1] = 30.0; */
/*     state_st[2] = -M_PI/1.5; */
/*     state_st[3] = 0.0; */
/*     state_st[4] = 0.0; */
/*     for (size_t ii = 0; ii < dx; ii++){ */
/*         state[ii] = state_st[order[ii]]; */
/*     } */
/*     con[0] = 0; */

/*     traj = NULL; */
/*     time = 0.0; */
/*     trajectory_add(&traj,5,1,time,state,con); */
/*     while (time < final_time){ */
/*         res = trajectory_step(traj,ode_sys,dt); */
/*         double * ls = trajectory_get_last_state(traj); */
/*         if ((fabs(ls[0]) < w/2.0) && (fabs(ls[1]) < w/2.0) ){ */
/*             break; */
/*         } */
/*         if (res != 0){ */
/*             break; */
/*         } */
/*         time = time + dt; */
/*     } */
/*     if (verbose == 1){ */
/*         trajectory_print(traj,stdout,4); */
/*     } */
/*     sprintf(filename,"%s/%s.dat",dirout,"traj5"); */
/*     fp = fopen(filename,"w"); */
/*     assert (fp != NULL); */
/*     trajectory_print(traj,fp,4); */
/*     fclose(fp); */
/*     trajectory_free(traj); traj = NULL; */


    integrator_destroy(ode_sys); ode_sys = NULL;    
    valuef_destroy(cost); cost = NULL;
    c3control_destroy(c3c); c3c = NULL;
    diag_destroy(&diag); diag = NULL;
    c3opt_free(opt); opt = NULL;
    approx_args_free(aargs); aargs = NULL;

    return 0;
}

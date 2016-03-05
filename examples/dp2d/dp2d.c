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

int f1(double t, double * x, double * u, double * out, double * jac,
       void * args)
{
    (void)(t);
    (void)(args);
    
    out[0] = x[1];// + sin(4.0*x[0]);
    out[1] = u[0];

    if (jac != NULL){
        //df1/du
        jac[0] = 0.0;
        jac[1] = 1.0;
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
    
    out[0] = 1.0;//sin(3.0*x[0]);
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 1.0;//cos(8.0*x[1]);


    if (grad != NULL){
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
    }
    return 0;
}

int polfunc(double t, double * x, double * u)
{
    (void)(t);
    if (x[1] > 0){
        u[0] = -1.0;
    }
    else if (x[1] < 0){
        u[0] = 1.0;
    }
    else{
        u[0] = 0.0;
    }
    // u[0] = 0.0;
    return 0;
}

int outbounds(double time, double * x, void * args, int * dirs)
{
    
    // two boundaries
    // outside of the box([-2,2]^2)
    // and the inside box [-0.1,0.1]^2
    (void)(time);
    (void)(args);
    
    double bound = 2.0;
    int out = 0;
    for (size_t ii = 0; ii < 2; ii++){
        dirs[ii] = 0;
        if (x[ii] > bound){
            return 1;
        }
        else if (x[ii] < -bound){
            return 1;
        }
        else if (fabs(x[ii] - bound) < 1e-15){
            // for onbound behaviour other than absorbing
            // uncomment the next two lines
            // dirs[ii] = 1;
            //out = -1;
            
            return 1;
        }
        else if (fabs(x[ii] + bound) < 1e-15){
            // for onbound behaviour other than absorbing
            // uncomment the next two lines
            //dirs[ii] = -1;
            //out = -1;
            
            return 1;
        }
    }

    if ( (fabs(x[0]) <= 2e-1) && (fabs(x[1]) <= 2e-1)){
//        printf("here!!!!! (%G,%G)\n",x[0],x[1]);
        return 1;
    }
    
    return out;
}

int stagecost(double t, double * x, double * u, double * out, 
              double * grad)
{
    (void)(t);
    *out = 0.0;
    *out += pow(x[0],2) + pow(x[1],2) + pow(u[0],2);

    if (grad!= NULL){
        grad[0] = 2 * u[0];
    }

    return 0;
}

int boundcost(double t, double * x, double * out)
{

    (void)(t);
    *out = 0.0;

    if ( (fabs(x[0]) <= 2e-1) && (fabs(x[1]) <= 2e-1)){
        //printf("here !\n");
        *out = 0.0;
    }
    else{
        *out = 5.0;
    }

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

    size_t dx = 2;
    size_t dw = 2;
    size_t du = 1;
    
    struct Drift drift;
    drift_init(&drift,dx,du,NULL,NULL,NULL,NULL);
    drift.b = f1;

    struct Diff diff;
    diff_init(&diff,dw,dx,du,NULL,NULL,NULL,NULL);
    diff.s = s1;

    struct Dyn dyn;
    dyn_init_ref(&dyn,&drift,&diff);

    struct Boundary * bound = boundary_alloc(dx,outbounds,NULL);
    
    //double t3[2];
    struct Policy * pol = policy_alloc();
    policy_init(pol,dx,du,NULL,NULL);
    policy_add_feedback(pol,polfunc);
    //policy_add_transform_ref(pol,&lt,t3);

    double lb[2] = {-2.0, -2.0};
    double ub[2] = {2.0, 2.0};

    size_t Narr[2] = {N, N};
    double * x = linspace(lb[0],ub[0],Narr[0]);
    double * y = linspace(lb[1],ub[1],Narr[1]);
    if (verbose > 0){
        printf("x = "); dprint(Narr[0],x);
        printf("y = "); dprint(Narr[1],y);
    }
    double * z[2]; 
    z[0] = x;
    z[1] = y;
    struct Cost * cost = cost_alloc(dx,lb,ub);
    cost_init_discrete(cost,Narr,z);
    cost_approx(cost,startcost,NULL,verbose-1);

    double h[2] = {x[1]-x[0],y[1]-y[0]};
    free(x); x = NULL;
    free(y); y = NULL;
    struct MCA * mm = mca_alloc(dx,du,dw,h);
    mca_attach_dyn(mm,&dyn);
    mca_attach_bound(mm,bound);

    double pt[2] = {1.0,1.0};
    double u[1] = {-1.0};
    double dt;
    int absorb;
    double val = mca_expectation(mm,0.0,pt,u,&dt,NULL,cost_eval_bb,cost,&absorb,NULL);
//    assert (absorb == 0);
    printf("value is %G\n", val);
    double grad[1];
    double gdt[1];
    val = mca_expectation(mm,0.0,pt,u,&dt,gdt,cost_eval_bb,cost,&absorb,grad);
    printf("value is %G, grad=%G, dtG=%G\n", val,grad[0],gdt[0]);
    // below should be true if startcost is constant
    //assert(fabs(val - startcost(pt,NULL)) < 1e-15);
    
    double beta = 1.0;
    struct DPih * dp = dpih_alloc(beta,stagecost,boundcost);
    dpih_attach_mca(dp, mm);
    dpih_attach_cost(dp, cost);
    dpih_attach_policy(dp, pol);


    //double delta;
    for (size_t ii = 0; ii < niter+1; ii++){

        FILE *fp2;
        char filename[256];
        sprintf(filename,"%s/%s_%zu.dat",dirout,"costfunc",ii);
        fp2 =  fopen(filename, "w");
        if (fp2 == NULL){
            fprintf(stderr, "cat: can't open %s\n", filename);
            return 0;
        }

        fprintf(fp2,"x y f\n");
        size_t N1 = 50;
        size_t N2 = 50;
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
        fclose(fp2);
        
        struct Cost * newcost = dpih_iter_pol(dp,verbose-1);
        cost_free(cost);
        cost = newcost;
        dpih_attach_cost(dp,cost);


    /*     delta = dpih_pi_iter_approx(&prob,verbose); */
        if (verbose != 0){
            printf("ii=%zu\n",ii);
        }
    }

    double pt2[2] = {0.00,0.00};
    double value;
    cost_eval(cost,0.0,pt2,&value);
    printf("value = %G -- should be %G\n",value,0.0);

    pt2[0] = 0.6;
    pt2[1] = -0.4;
    cost_eval(cost,0.0,pt2,&value);
    printf("value = %G \n",value);

    u[0] = 0.5;
    val = mca_expectation(mm,0.0,pt,u,&dt,gdt,cost_eval_bb,cost,&absorb,grad);
    printf("value is %3.15G, grad=%3.15G,gdt=%G\n", val,grad[0],gdt[0]);
    double dhu = 1e-4;
    u[0] = 0.5+dhu;
    double val2;
    double dt2;
    val2 = mca_expectation(mm,0.0,pt,u,&dt2,gdt,cost_eval_bb,cost,&absorb,grad);
    printf("value is %3.15G, grad=%3.15G,gdt=%G\n", val,grad[0],gdt[0]);
    printf("grad should be about %3.15G\n",(val2-val)/dhu);
    printf("gdt should be about %3.15G\n",(dt2-dt)/dhu);

    boundary_free(bound);
    policy_free(pol);
    //printf("free mca\n");
    mca_free(mm);
    //printf("free cost\n");
    cost_free(cost);
    //printf("free dpih\n");
    dpih_free(dp);
    
    return 0;
}

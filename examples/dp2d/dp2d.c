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

int f1(double t, double * x, double * u, double * out, double * jac,
       void * args)
{
    (void)(t);
    (void)(args);
    
    out[0] = x[1] + sin(4.0*x[0]);
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
    
    out[0] = sin(3.0*x[0]);
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = cos(8.0*x[1]);

    if (grad != NULL){
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;//cos(8.0*x[1]);
    }
    return 0;
}

int polfunc(double t, double * x, double * u, void * arg)
{
    (void)(t);
    (void)(arg);
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
    double lb[2] = {-2.0, -2.0};
    double ub[2] = {2.0, 2.0};
    size_t Narr[2] = {N, N};

    double lbu[1] = {-1.0};
    double ubu[1] = {1.0};
    struct c3Opt * opt = c3opt_alloc(BFGS,du);
    c3opt_add_lb(opt,lbu);
    c3opt_add_ub(opt,ubu);
    
    double beta = 1.0;

    // setup problem
    c3sc sc = c3sc_create(IH,dx,du,dw);
    c3sc_set_state_bounds(sc,lb,ub);
    c3sc_add_dynamics(sc,f1,NULL,s1,NULL);
    c3sc_init_mca(sc,Narr);
    c3sc_attach_opt(sc,opt);
    c3sc_init_dp(sc,beta,stagecost,boundcost,NULL);

    size_t N1 = 50, N2 = 50;
    struct DPih * dp = c3sc_get_dp(sc);
    struct Cost * cost = dpih_get_cost(dp);
    cost_approx(cost,startcost,NULL,verbose-1);
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
        
        //struct Cost * newwhcost = dpih_iter_pol(dp,verbose-1);
        struct Cost * newcost = dpih_iter_vi(dp,verbose-1);
        cost_free(cost);
        cost = newcost;
        dpih_attach_cost(dp,cost);


    /*     delta = dpih_pi_iter_approx(&prob,verbose); */
        if (verbose != 0){
            printf("ii=%zu\n",ii);
        }
    }
    
    c3sc_destroy(sc);
    return 0;
}

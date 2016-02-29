#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "c3.h"

#include "util.h"
#include "cost.h"

int boundary_ob(struct Boundary * bound, double * x)
{

    int outbounds = bound->bcheck(x,bound->args);
    // 0 inbounds
    return outbounds;
    
}

double dpih_rhs(struct DPih * dp, double t, double * x,
                double * u)
{

    // first check if x is in bounds
    int bound = boundary_ob(dp->bound,x);

    double out; 
    if (bound == 0){
        // inbounds, can do standard kushner update

        NEED TO ADD SCALING
        // stagecost
        int res = dp->stagecost(t,x,u,&out);
        assert (res == 0);

        // additional cost
        double dt;
        double probs[1000];
        out += tensor_mm_cost(dp->mm,t,x,u,&dt,probs);
    }
    else{

        // check direction of movement
        int dirin = 0;
        if (dirin)
        { //direction of movement is into the domain so ok
            
        }
        else{
            // direction of movement is out of the domain so
            // need to use boundary
            int res = dp->boundcost(t,x,&out);
            assert (res == 0);
        }
    }

    return out;
}


double dpih_pi_iter(struct DPih * dp)
{

    return 1.0;
    
}

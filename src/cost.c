#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "util.h"
#include "cost.h"

struct Cost * cost_alloc()
{
    struct Cost * cost = NULL;
    cost = malloc(sizeof(struct Cost));
    if (cost == NULL){
        fprintf(stderr,"Could not allocate cost\n");
        exit(1);
    }
    return cost;
}

void cost_free(struct Cost * c)
{
    if (cost != NULL){
        bounding_box_free(c->bds); c->bds = NULL;
        function_train_free(c->cost); c->cost = NULL;
    }
}

int cost_eval(struct Cost * cost, double time, double * x,
              double * eval)
{

    int res = c3sc_check_bounds(cost->bds->dim,
                                cost->bds->lb,
                                cost->bds->ub,
                                x);

    if (res != 0){
        return res;
    }

    *eval = function_train_eval(cost->cost,x);

    return 0;
    
}


int cost_eval_neigh(struct Cost * cost, double time,
                    double * x, size_t ii, double pt[2],
                    double evals[2]);
{
    // special function train eval;

    return 0;
    
}

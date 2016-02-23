#include <stdlib.h>

#include "util.h"

int c3sc_check_bounds(size_t dx, double * lbx, double * ubx, double * x)
{
    if ((lbx == NULL) || (ubx == NULL)){
        return 0;
    }
    for (size_t ii = 0; ii < dx; ii++){
        if (x[ii] < lbx[ii]){
            return - (ii+1);
        }
        else if (x[ii] > ubx[ii]){
            return ii+1;
        }
    }
    return 0;
}

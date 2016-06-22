#ifndef BELLMAN_H
#define BELLMAN_H

#include <stdlib.h>

double bellmanrhs(size_t dx, size_t du, double stage_cost, const double * stage_grad, 
                      double discount, const double * prob, const double * prob_grad, 
                      double dt, const double * dtgrad, const double * cost, double * grad);
#endif

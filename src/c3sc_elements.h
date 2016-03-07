#ifndef c3sc_ELEMENTS_H
#define c3sc_ELEMENTS_H

/** \struct State
 *  \brief Defines a state
 *  \var State::d
 *  dimension
 *  \var State::t
 *  time
 *  \var State::x
 *  state
 *  \var State::spec
 *  some additional arguments
 */
struct State
{
    size_t d; 
    double t; 
    double * x; 
    void * spec; 
};

/** \struct Control
 *  \brief Defines a control
 *  \var Control::d
 *  dimension
 *  \var Control::u
 *  control value
 */
struct Control
{
    size_t d;
    double * u;
};

/** \struct Trajectory
 *  \brief Linked List defining a trajectory
 *  \var Trajectory::s
 *  State
 *  \var Trajectory::u
 *  Control
 *  \var Trajectory::next
 *  Pointer to next element
 */
struct Trajectory
{
    struct State * s;
    struct Control * u;
    struct Trajectory * next;
};


/** \struct Policy
 *  \brief Policy
 *  \var Policy::dx
 *  dimension of state space
 *  \var Policy::du
 *  dimension of control
 *  \var Policy::lbx
 *  Lower bounds for state space
 *  \var Policy::ubx
 *  Upper bounds for state space
 *  \var Policy::feedback
 *  Function mapping time,state to control
 */
struct Policy
{
    size_t dx;
    size_t du;
    double * lbx;
    double * ubx;
    int (*feedback)(double, double *, double *);
};


struct Drift;
struct Diff;

/** \struct Dyn
 *  \brief Stochastic Differential Equation Dynamics
 *  \var Dyn::drift
 *  drift dynamics
 *  \var Dyn::diff
 *  diffusion dynamics
 *  \var Dyn::s
 *  RHS of diffusion term of stochastic differential equation
 *  \var Dyn::sargs
 *  Additional arguments to dynamics
 */
struct Dyn
{
    struct Drift * drift;
    struct Diff  * diff;
};




#endif

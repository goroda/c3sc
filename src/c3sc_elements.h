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


struct Drift;
struct Diff;
struct Dyn;


#endif

#pragma once
#include "../riemann/state.h"

struct Cons{
    double mass;
    double momentum;
    double energy;
};

//im converting primitive variable to conserved variables to use the conservative laws in calculations
Cons prim_to_cons(const State& W);
State cons_to_prim(const Cons& U);
Cons euler_flux(const State& W);    //returns eulerh's flux for state W
double sound_speed(const State& W); //computing speed of sound
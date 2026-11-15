#include "../header/euler.h"
#include <cmath>

//converting primitive variables (u, rho, p) to conserved varibles for conservative laws
Cons prim_to_cons(const State& W){
    Cons U;
    U.mass = W.rho;
    U.momentum = W.rho * W.u;
    U.energy   = W.p / (GAMMA - 1.0) + 0.5 * W.rho * W.u * W.u;     
    return U;
}

//conservative -> primitive
State cons_to_prim(const Cons& U){
    State W;
    W.rho = U.mass;
    W.u   = U.momentum / U.mass;
    W.p   = (GAMMA - 1.0) * (U.energy - 0.5 * U.momentum * U.momentum / U.mass);   
    return W;
}

Cons euler_flux(const State& W){
    const double E = W.p / (GAMMA - 1.0) + 0.5 * W.rho * W.u * W.u;
    Cons F;
    F.mass = W.rho * W.u;                    // mass carried by the flow
    F.momentum = W.rho * W.u * W.u + W.p;        // momentum carried by the flow + the push of pressure
    F.energy   = W.u * (E + W.p);                // energy carried by the flow + work done by the pressure
    return F;
}

double sound_speed(const State& W){
    return std::sqrt(GAMMA * W.p / W.rho);
}

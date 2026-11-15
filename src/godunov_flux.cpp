#include "../header/godunov_flux.h"
#include "../header/riemann_at.h"

Cons godunov_flux(const State& WL, const State& WR){
    State W0 = riemann_at(WL, WR, 0.0);      // for x/t = 0
    return euler_flux(W0);                   // physical flux of that state
}

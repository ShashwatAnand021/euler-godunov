#include "../header/godunov.h"
#include "../header/godunov_flux.h"
#include <cmath>
#include <algorithm>

void apply_transmissive(std::vector<Cons>& U){
    int last = (int)(U.size()) - 1;      // index of the right ghost cell
    U[0] = U[1];                //a copy of the first cell itself
    U[last] = U[last - 1];      //copy of last cell
}

double compute_dt(const std::vector<Cons>& U, double dx, double cfl, double t, double t_end){
    double smax = 0.0;
    for(const Cons& c : U){                         // "for each cell c in U" (ghosts included)
        State W = cons_to_prim(c);                              // c as primitives
        double speed = sound_speed(W) + std::fabs(W.u);                         // |u| + a
        smax = std::max(smax, speed);
    }
    double dt = (cfl*dx) / smax;
    if(t + dt > t_end) dt = t_end - t;                               // don't step past t_end
    return dt;
}

void compute_fluxes(const std::vector<Cons>& U, std::vector<Cons>& F){
    int M = static_cast<int>(U.size()) - 2;         // number of real cells
    for(int i = 0; i <= M; i++){                    // M+1 walls
        State WL = cons_to_prim(U[i]);                             // the cell on the left of wall i
        State WR = cons_to_prim(U[i + 1]);                          // the cell on the right of wall i
        F[i] = godunov_flux(WL, WR);
    }
}

void update(std::vector<Cons>& U, const std::vector<Cons>& F, double dt, double dx){
    int M = (int)(U.size()) - 2;
    double dt_over_dx = dt/dx;
    for(int i = 1; i <= M; i++){                    // avoiding ghost cells
        // cell i sits between wall i-1 (left) and wall i (right)
        U[i].mass     += dt_over_dx * (F[i - 1].mass     - F[i].mass);
        U[i].momentum += dt_over_dx * (F[i - 1].momentum - F[i].momentum);
        U[i].energy   += dt_over_dx * (F[i - 1].energy   - F[i].energy);
    }
}
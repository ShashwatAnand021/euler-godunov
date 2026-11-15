#include <iostream>
#include <fstream>
#include <vector>
#include "../header/euler.h"
#include "../header/godunov.h"
#include "../header/riemann_at.h"

int main(){

    // ---- problem: Toro Test 1 ----
    State WL{1.0, 1.0, 0.75};              // { p, rho, u }
    State WR{0.1, 0.125, 0.0};
    double x0 = 0.3, t_end = 0.2;
    int M = 100;
    double x_min = 0.0, x_max = 1.0;
    double dx  = (x_max - x_min) / M; 
    double cfl = 0.9;

    std::vector<Cons> U(M+2), F(M+1);

    // ---- initial condition ----
    Cons UL = prim_to_cons(WL), UR = prim_to_cons(WR);     // the two states as conserved variables
    for(int i = 1; i <= M; i++){
        double x = x_min + (i - 0.5) * dx;                 // centre of cell i
        U[i] = (x < x0) ? UL : UR;
    }

    // ---- time loop ----
    double t = 0.0;
    int step = 0;
    while(t < t_end - 1e-12){
        apply_transmissive(U);                                // 1. ghosts
        double c  = (step < 5) ? 0.2 * cfl : cfl;             // 2. smaller CFL for the first 5 steps
        double dt = compute_dt(U, dx, c, t, t_end);
        compute_fluxes(U, F);                                 // 3. fluxes
        update(U, F, dt, dx);                                 // 4. updating
        t += dt;
        step++;
    }

    // ---- output ----
    std::ofstream out("godunov.csv");
    out << "x,rho,u,p,rho_exact,u_exact,p_exact\n";
    for(int i = 1; i <= M; i++){
        double x = x_min + (i - 0.5) * dx;
        State W = cons_to_prim(U[i]);                         // numerical solution
        State E = riemann_at(WL, WR, (x - x0) / t_end);       // exact solution at this x, from the solver
        out << x << "," << W.rho << "," << W.u << "," << W.p << ","
            << E.rho << "," << E.u << "," << E.p << "\n";
    }
    std::cout << "steps = " << step << "\n";
}
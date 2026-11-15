#include "../header/riemann_at.h"
#include "../riemann/solver.h"
#include "../riemann/newtonRaphson.h"
#include "../riemann/rho_star.h"
#include "../riemann/sample_solution.h"

State riemann_at(const State& WL, const State& WR, double S){

    double p_star = calc_p_star(WL, WR);                         // Newton-Raphson for p*

    auto [f_l, fd_l] = calc_f(p_star, WL);
    auto [f_r, fd_r] = calc_f(p_star, WR);
    double u_star = 0.5 * (WL.u + WR.u) + 0.5 * (f_r - f_l);     // u*

    double rho_star_l = calc_rho_star(WL, p_star);               // rho_l* and rho_r*
    double rho_star_r = calc_rho_star(WR, p_star);

    return sample_solution(S, WL, WR, p_star, u_star, rho_star_l, rho_star_r);
}

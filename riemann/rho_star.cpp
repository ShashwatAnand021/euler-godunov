#include"state.h"
#include"rho_star.h"
#include<cmath>

double calc_rho_star(const State& W, const double p_star){
    if (p_star > W.p) {
        //shock
        double ratio = p_star / W.p;
        double gamma_term = (GAMMA - 1.0) / (GAMMA + 1.0);
        return W.rho * (ratio + gamma_term) / (gamma_term * ratio + 1.0);
    } else {
        //rarefaction
        return W.rho * std::pow(p_star / W.p, 1.0 / GAMMA);
    }
}

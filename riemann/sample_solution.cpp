#include "sample_solution.h"
#include <cmath>

double get_sound_speed(double p, double rho) {
    return std::sqrt(GAMMA * p / rho);
}

State sample_solution(double S, const State& WL, const State& WR, 
                      double p_star, double u_star, double rho_star_l, double rho_star_r) {
    
    State out;
    double a_L = get_sound_speed(WL.p, WL.rho);
    double a_R = get_sound_speed(WR.p, WR.rho);
    double a_star_L = get_sound_speed(p_star, rho_star_l);
    double a_star_R = get_sound_speed(p_star, rho_star_r);

    if (S <= u_star) {
        if (p_star > WL.p) { 
            double S_L = WL.u - a_L * std::sqrt(((GAMMA + 1.0) / (2.0 * GAMMA)) * (p_star / WL.p) + ((GAMMA - 1.0) / (2.0 * GAMMA)));
            if (S < S_L) {
                out = WL; 
            } else {
                out = {p_star, rho_star_l, u_star}; 
            }
        } else { 
            double S_HL = WL.u - a_L;
            double S_TL = u_star - a_star_L;
            
            if (S < S_HL) {
                out = WL; 
            } else if (S > S_TL) {
                out = {p_star, rho_star_l, u_star}; 
            } else {
                out.u = (2.0 / (GAMMA + 1.0)) * (a_L + ((GAMMA - 1.0) / 2.0) * WL.u + S);
                double a_fan = (2.0 / (GAMMA + 1.0)) * (a_L + ((GAMMA - 1.0) / 2.0) * (WL.u - S));
                out.p = WL.p * std::pow(a_fan / a_L, 2.0 * GAMMA / (GAMMA - 1.0));
                out.rho = WL.rho * std::pow(a_fan / a_L, 2.0 / (GAMMA - 1.0));
            }
        }
    } 
    else {
        if (p_star > WR.p) { 
            double S_R = WR.u + a_R * std::sqrt(((GAMMA + 1.0) / (2.0 * GAMMA)) * (p_star / WR.p) + ((GAMMA - 1.0) / (2.0 * GAMMA)));
            if (S > S_R) {
                out = WR; 
            } else {
                out = {p_star, rho_star_r, u_star}; 
            }
        } else { 
            double S_HR = WR.u + a_R;
            double S_TR = u_star + a_star_R;
            
            if (S > S_HR) {
                out = WR; 
            } else if (S < S_TR) {
                out = {p_star, rho_star_r, u_star}; 
            } else {
                out.u = (2.0 / (GAMMA + 1.0)) * (-a_R + ((GAMMA - 1.0) / 2.0) * WR.u + S);
                double a_fan = (2.0 / (GAMMA + 1.0)) * (a_R - ((GAMMA - 1.0) / 2.0) * (WR.u - S));
                out.p = WR.p * std::pow(a_fan / a_R, 2.0 * GAMMA / (GAMMA - 1.0));
                out.rho = WR.rho * std::pow(a_fan / a_R, 2.0 / (GAMMA - 1.0));
            }
        }
    }
    return out;
}
#include"state.h"
#include"newtonRaphson.h"
#include"solver.h"
#include<iostream>
#include<cmath>

const int MAX_ITER = 100;
const double TOL = 1e-6;

double calc_p_star(const State& WL, const State& WR){

    double p_old = 0.5 * (WL.p + WR.p);     //initial guess - AM of the initial pressures
    double p_new = p_old;

    for(int i=0;i<MAX_ITER;i++){

        auto [f_l, fd_l] = calc_f(p_old, WL);
        auto [f_r, fd_r] = calc_f(p_old, WR);

        double f = f_l + f_r + (WR.u - WL.u);
        double fd = fd_l + fd_r;

        p_new = p_old - f/fd;       //the recurring relation

        if(p_new < 0.0) p_new = TOL;    //since pressure can't be -ve

        double change = std::abs(p_old - p_new)/p_old;
        if(change < TOL) {
            return p_new;
        }
        p_old = p_new;
    }

    std::cout << "NR didn't converge!\n";
    return p_new;

}
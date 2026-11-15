#include"state.h"
#include"solver.h"
#include<cmath>

std::pair<double,double> calc_f (double p, const State& Wk){
    double f, fd;

    double A = 2/((GAMMA + 1)*Wk.rho);
    double B = ((GAMMA - 1)*Wk.p)/(GAMMA + 1);
    double a = std::sqrt((GAMMA*Wk.p)/Wk.rho);

    //calclating f and fd

    if(p > Wk.p){
        //shock branch
        f = (p - Wk.p) * std::sqrt((A / (p + B)));
        fd = std::sqrt(A/(B + p)) * (1 - (p - Wk.p)/(2*(B + p)));
    }else{
        //rarefaction branch
        f = ((2.0*a) / (GAMMA - 1)) * ((std::pow(p/Wk.p, (GAMMA-1)/(2*GAMMA))) - 1);
        fd = (1.0/(Wk.rho * a)) * std::pow((p/Wk.p),(-(GAMMA + 1)/(2*GAMMA)));
    }

    return {f, fd};
}

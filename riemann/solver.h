#pragma once
#include"state.h"
#include<utility>

//f(p) = fl(p) + fr(p) + ul- ur = 0 
// calc. p* = p using the Newton Raphson method
// p_n+1 = p_n - f_p_n/f'_p_n

std::pair<double,double> calc_f (double p, const State& Wk);

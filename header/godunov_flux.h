#pragma once
#include"../riemann/state.h"
#include"euler.h"

Cons godunov_flux(const State& WL, const State& WR);
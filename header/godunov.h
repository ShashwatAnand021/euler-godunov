#pragma once
#include"euler.h"
#include<vector>

void apply_transmissive(std::vector<Cons>& U);
double compute_dt(const std::vector<Cons>& U, double dx, double cfl, double t, double t_end);
void compute_fluxes(const std::vector<Cons>& U, std::vector<Cons>& F);
void update(std::vector<Cons>& U, const std::vector<Cons>& F, double dt, double dx);
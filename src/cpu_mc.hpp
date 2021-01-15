#ifndef CPU_MC_H
#define CPU_MC_H
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

double monte_carlo_call_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T) {
  double S_adjust = S * exp(T*(r-0.5*v*v));
  double S_cur = 0.0;
  double payoff_sum = 0.0;
  std::mt19937 gen(42);
  std::normal_distribution<double> d(0, 1);

  for (int i=0; i<num_sims; i++) {
    double gauss_bm = d(gen);
    S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
    payoff_sum += std::max(S_cur - K, 0.0);

  }

  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}

double monte_carlo_put_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T) {
  double S_adjust = S * exp(T*(r-0.5*v*v));
  double S_cur = 0.0;
  double payoff_sum = 0.0;
  std::mt19937 gen(42);
  std::normal_distribution<double> d(0, 1);

  for (int i=0; i<num_sims; i++) {
    double gauss_bm = d(gen);
    S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
    payoff_sum += std::max(K - S_cur, 0.0);
  }

  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}



#endif

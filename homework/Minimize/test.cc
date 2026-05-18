#include <ios>
#include <iostream>
#include <iomanip>
#include "core/vector.h"
#include "minimize.h"

using namespace pp;
int main(){
  size_t evals = 0;
  auto Rosenbrock = [&](const vector& x){
    double val1 = 1.0 - x[0];
    double val2 = x[1] - x[0]*x[0];
    evals += 1;
    return val1*val1 + 100*val2*val2;
  };
  auto Himmelblau = [&](const vector& x){
    double val1 = (x[0]*x[0]+x[1]-11.0);
    double val2 = (x[0]+x[1]*x[1]-7.0);
    evals += 1;
    return val1*val1 + val2*val2;
  };
  NewtonForward nf; BacktrackingSearch s; FindMinimum mini_f(nf, s);
  NewtonCentral nc; FindMinimum mini_c(nc, s);
  vector x0 = {5.0, -1.0};

  auto res_f = mini_f.minimize(Rosenbrock, x0, 1e-8); int evals_f = evals; evals = 0;
  auto res_c = mini_c.minimize(Rosenbrock, x0, 1e-8); int evals_c = evals; evals = 0;
  std::cout << std::scientific << std::setprecision(10) << "---Rosenbrock---\n"
    << "Minima: (1, 1)\n"
    << "---Forward difference---\n"
    << "Result: " << res_f.first
    <<"\nFunction evals: " << evals_f
    << "\n---Central difference---\n"
    << "Result: " << res_c.first
    << "\nFunction evals: " << evals_c << "\n";
  x0 = {3, 2};
  res_f = mini_f.minimize(Himmelblau, x0, 1e-8); evals_f = evals; evals = 0;
  res_c = mini_c.minimize(Himmelblau, x0, 1e-8); evals_c = evals;
  std::cout << "\n---Himmelblau---\n"
    << "Minima: (3, 2)\n"
    << "---Forward difference---\n"
    << "Result: " << res_f.first
    << "\nFunction evals: " << evals_f
    << "\n---Central difference---\n"
    << "Result: " << res_c.first
    << "\nFunction evals: " << evals_c << "\n";
  return 0;
}

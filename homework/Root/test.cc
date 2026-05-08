#include <iostream>
#include "core/vector.h"
#include "root.h"

using namespace pp;
int main(){
  NewtonUpdate n(2);
  BroydenUpdate b(2);
  BacktrackingSearch s1;
  QuadraticInterpSearch s2(1,1);
  FindRoot solver1(n, s1);
  FindRoot solver2(n, s2);
  FindRoot solver3(b, s1);
  FindRoot solver4(b, s2);
  int evals = 0;
  auto Rosenbrock = [&](const vector& x){
    vector res(x);
    res[0] = -2.0*(1.0 - x[0]) - 400*x[0]*(x[1]-x[0]*x[0]);
    res[1] = 200*(x[1] - x[0]*x[0]);
    evals += 1;
    return res;
  };
  auto Himmelblau = [&](const vector& x){
    vector res(x);
    double val1 = 2.0*(x[0]*x[0]+x[1]-11.0);
    double val2 = 2.0*(x[0]+x[1]*x[1]-7.0);
    res[0] = val1*2.0*x[0] + val2;
    res[1] = val1 + val2*2.0*x[1];
    evals += 1;
    return res;
  };
  vector x0{5, -1};
  size_t max_steps = 4000;
  double acc = 1e-10, alpha_min=1e-6;
  auto res1 = solver1.solve(Rosenbrock, x0, acc, alpha_min, max_steps); int evals1 = evals;

  evals = 0;
  auto res2 = solver2.solve(Rosenbrock, x0, acc, alpha_min, max_steps); int evals2 = evals;

  evals = 0;
  auto res3 = solver3.solve(Rosenbrock, x0, acc, alpha_min, max_steps); int evals3 = evals;

  evals = 0; b.reset();
  auto res4 = solver4.solve(Rosenbrock, x0, acc, alpha_min, max_steps); int evals4 = evals;
  
  std::cout << "---Convergense test for Rosenbrock with start guess (5, -1)---\n"
    << "Newton's method with backtracking linesearch\n" << "Result: " << res1
    << "\nFunction evals: " << evals1 << "\n" 
    << "Newton's method with quadratic interpolation\n" << "Result: " << res2
    << "\nFunction evals: " << evals2 << "\n"
    << "Broyden's method with backtracking linesearch\n" << "Result: " << res3
    << "\nFunction evals: " << evals3 << "\n"
    << "Broyden's method with quadratic interpolation\n" << "Result: " << res4
    << "\nFunction evals: " << evals4 << "\n\n";

  evals = 0;
  res1 = solver1.solve(Himmelblau, x0, acc, alpha_min, max_steps); evals1 = evals;
  evals = 0;
  res2 = solver2.solve(Himmelblau, x0, acc, alpha_min, max_steps); evals2 = evals;
  evals = 0;
  b.reset();
  res3 = solver3.solve(Himmelblau, x0, acc, alpha_min, max_steps); evals3 = evals;
  evals = 0;
  b.reset();
  res4 = solver4.solve(Himmelblau, x0, acc, alpha_min, max_steps); evals4 = evals;
  std::cout << "---Convergense test for Himmelblau with start guess (5, -1)---\n"
    << "Newton's method with backtracking linesearch\n" << "Result: " << res1
    << "\nFunction evals: " << evals1 << "\n" 
    << "Newton's method with quadratic interpolation\n" << "Result: " << res2
    << "\nFunction evals: " << evals2 << "\n"
    << "Broyden's method with backtracking linesearch\n" << "Result: " << res3
    << "\nFunction evals: " << evals3 << "\n"
    << "Broyden's method with quadratic interpolation\n" << "Result: " << res4
    << "\nFunction evals: " << evals4 << "\n";
  return 0;
}

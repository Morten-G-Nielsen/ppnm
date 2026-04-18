#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "quad.h"

int main(){
  double acc = 1e-6, eps=1e-6;
  pp::NewtonCotes nc;
  auto f1 = [](double x){return std::sqrt(x);};
  auto f2 = [](double x){return 1/std::sqrt(x);};
  auto f3 = [](double x){return std::sqrt(1 - (x*x));};
  auto f4 = [](double x){return std::log(x)/std::sqrt(x);};
  double sol1 = 2.0/3.0, sol2 = 2.0, sol3 = M_PI/4.0, sol4 = -4.0;
  double a = 0, b = 1;
  
  pp::Integrator quad(nc, acc, eps);
  pp::QuadResult res = quad.integrate(f1, a, b);
  std::cout << std::fixed << std::setprecision(10);
  std::cout << "---Newton-Cotes---\n" << "Integral of sqrt(x) from 0 to 1\n";
  std::cout << "Calculated integral:  " << res.integral << "\nTrue solution:        "<< sol1;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol1) << "\nFunction evals:       " << res.evals << "\n\n";

  res = quad.integrate(f2, a, b);
  std::cout <<"Integral of 1/sqrt(x) from 0 to 1\n";
  std::cout << "Calculated integral:  " << res.integral << "\nTrue solution:        "<< sol2;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol2) << "\nFunction evals:       " << res.evals << "\n\n";

  res = quad.integrate(f3, a, b);
  std::cout << "Integral of sqrt(1 - x^2) from 0 to 1\n";
  std::cout << "Calculated integral:  " << res.integral << "\nTrue solution:        "<< sol3;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol3) << "\nFunction evals:       " << res.evals << "\n\n";

  res = quad.integrate(f4, a, b);
  std::cout << "Integral of ln(x)/sqrt(x) from 0 to 1\n";
  std::cout << "Calculated integral:  " << res.integral << "\nTrue solution:        "<< sol4;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol4) << "\nFunction evals:       " << res.evals << "\n\n";

  pp::ClenshawCurtis cc;
  pp::Integrator quad2(cc, acc, eps);

  res = quad2.integrate(f2, a, b);
  std::cout << "---Clenshaw-Curtis---\n" << "Integral of 1/sqrt(x) from 0 to 1\n";
  std::cout << "Calculated integral:  " << res.integral << "\nTrue solution:        "<< sol2;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol2) << "\nFunction evals:       " << res.evals << "\n\n";

  res = quad2.integrate(f4, a, b);
  std::cout << "Integral of ln(x)/sqrt(x) from 0 to 1\n";
  std::cout <<"Calculated integral: "<< res.integral << "\nTrue solution:       "<< sol4;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol4) << "\nFunction evals:       " << res.evals << "\n\n";

  pp::InfiniteRule inf(cc);
  pp::Integrator quad3(inf, acc, eps);
  auto f5 = [](double x){return std::exp(-(x*x));};
  auto f6 = [](double x){return 1/(1+x*x);};
  double sol5 = std::sqrt(M_PI), sol6 = M_PI/2.0;
  a = -INFINITY; b = INFINITY;

  res = quad3.integrate(f5, a, b);
  std::cout << "Integral of exp(-x^2) from -inf to inf\n";
  std::cout << "Calculated integral:  " << res.integral << "\nTrue solution:        "<< sol5;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol5) << "\nFunction evals:       " << res.evals << "\n\n";

  res = quad3.integrate(f6, 0, b);
  std::cout << "Integral of 1/(1 + x^2) from 0 to inf\n";
  std::cout << "Calculated integral:  " << res.integral << "\nTrue solution:        "<< sol6;
  std::cout << "\nEstimated error:      " << res.error << "\nTrue error:           " << std::abs(res.integral-sol6) << "\nFunction evals:       " << res.evals << "\n\n";
  return 0;
}

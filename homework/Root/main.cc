#include <ios>
#include <iostream>
#include <iomanip>
#include <functional>
#include <string>
#include "core/vector.h"
#include "ode.h"
#include "root.h"

using namespace pp;
int main(int argc, char* argv[]){
  double rmin = 1e-7;
  double rmax = 20.0;
  double acc = 1e-12;
  double eps = 1e-12;
  double h = 1e-3;
  int id = 1;
  for(size_t i = 0; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "-rmin" && i+1<argc) rmin = std::stod(argv[++i]);
    if(arg == "-rmax" && i+1<argc) rmax = std::stod(argv[++i]);
    if(arg == "-acc" && i+1<argc) acc = std::stod(argv[++i]);
    if(arg == "-eps" && i+1<argc) eps = std::stod(argv[++i]);
    if(arg == "-h" && i+1<argc) h = std::stod(argv[++i]);
    if(arg == "-id" && i+1<argc) id = std::stoi(argv[++i]);
  }
  vector yinit{rmin-rmin*rmin, 1.0-2.0*rmin};
  double E_exact = -0.5;

  auto func = [=](const vector& E_vec){
    double E = E_vec[0];
    auto f = [&](double x, const vector& y){
      vector dy(2);
      dy[0] = y[1];
      dy[1] = -2.0*(1.0/x + E)*y[0];
      return dy;
    };
    rk_45 rule;
    auto res = driver(rule, f, rmin, rmax, yinit, h, acc, eps);
    auto xvals = res.x_values;
    auto yvals = res.y_values;
    vector result = {yvals[xvals.size()-1][0]};
    return result;
  };
  NewtonUpdate n(1);
  BacktrackingSearch l;
  FindRoot root(n, l);

  vector x0 = {-1.0};
  auto res = root.solve(func, x0, 1e-12, 1e-3);
  double E0 = res[0];
  if(id == 1){
    std::cout << std::scientific << std::setprecision(15) << std::abs(E0 - E_exact) << "\n";
  } else if(id == 2){
    auto f = [&](double x, vector y){
      vector dy(2);
      dy[0] = y[1];
      dy[1] = -2.0*(1.0/x + E0)*y[0];
      return dy;
    };
    rk_45 stepper;
    auto ress = driver(stepper, f, rmin, rmax, yinit, h, acc, eps);
    auto xvals = ress.x_values;
    auto yvals = ress.y_values;
    for(size_t i = 0; i<xvals.size(); i++){
      std::cout << xvals[i];
      for(size_t j = 0; j<2; j++){
        std::cout << " " << yvals[i][j];
      }
      std::cout << "\n";
    }
  }
  return 0;
}

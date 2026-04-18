#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "quad.h"

double erf(pp::Integrator& quad, double z){
  if(z<0){
    return -erf(-z);
  } else if(z>=0 && z<=1){
    auto f = [](double x){return std::exp(-(x*x));};
    pp::QuadResult res = quad.integrate(f, 0, z);
    return 2/std::sqrt(M_PI)*res.integral;
  } else {
    auto f = [=](double x){return std::exp(-std::pow(z+(1-x)/x, 2))/x/x;};
    pp::QuadResult res = quad.integrate(f, 0, 1);
    return 1 - 2/std::sqrt(M_PI)*res.integral;
  }
}

int main(int argc, char* argv[]){
  double z = 0;
  double acc = 1e-6, eps = 1e-6;

  for(int i = 0; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "-z" && i+1<argc) z = std::stod(argv[++i]);
    if(arg == "-acc" && i+1<argc) acc = std::stod(argv[++i]);
    if(arg == "-eps" && i+1<argc) eps = std::stod(argv[++i]);
  }
  pp::ClenshawCurtis cc; pp::InfiniteRule inf(cc);
  pp::Integrator quad(inf, acc, eps);
  std::cout << std::setprecision(20) << z << " " << erf(quad, z) << "\n";
  return 0;
}

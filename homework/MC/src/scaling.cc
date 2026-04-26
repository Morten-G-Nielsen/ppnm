#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include "mc.h"

using namespace pp;
int main (int argc, char* argv[]) {
  size_t N = 100;
  for(size_t i = 0; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "-N" && i+1<argc) N = static_cast<size_t>(std::stod(argv[++i]));
  }

  std::vector<double> a{-1.0, -1.0};
  std::vector<double> b{1.0, 1.0};
  double exact = M_PI;

  auto f = [](const std::vector<double>& x){
    double val = x[0]*x[0] + x[1]*x[1];
    return (val <= 1) ? 1.0:0.0;
  }; 

  LCG gen(42);
  MCPlain rule(gen);
  MCIntegrator I(rule);
  auto res = I.integrate(f, a, b, N);
  std::cout << N << " " << std::abs(res.integral - exact) << " " << res.error << "\n";
  return 0;
}

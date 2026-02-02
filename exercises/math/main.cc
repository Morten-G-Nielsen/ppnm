#include <iostream>
#include <cmath>
#include <numbers>
#include <cstdio>
#include "sfuns.h"
int main(){
  std::cout << "sqrt(2)=" << std::sqrt(2.0) << "\n";
  std::cout << "2^(1/5)=" << std::exp2(1.0/5) << "\n";
  std::cout << "e^(pi)=" << std::exp(std::numbers::pi) << "\n";
  std::cout << "pi^(e)=" << std::pow(std::numbers::pi, std::exp(1)) << "\n\n";

  for(double x=1;x<=10;x+=1)
    std::cout << "fgamma("<< x <<")=" << sfuns::fgamma(x)
      << " tgamma("<< x <<")=" << std::tgamma(x) << "\n";
  double x = 100;
  std::cout << "\nlngamma(100)=" << sfuns::lngamma(x) << " lgamma(100)=" << std::lgamma(x) <<"\n";  
  return 0;
}

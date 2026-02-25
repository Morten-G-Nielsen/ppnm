#include <iostream>
#include <cmath>
#include <complex>
#include "sfuns.h"

const double PI = 3.1415926535897932; 
const std::complex<double> I(0,1);
int main(){
  std::cout << "sqrt(2) = " << std::sqrt(2.0) << "\n";
  std::cout << "2^(1/5) = " << std::exp2(1.0/5) << "\n";
  std::cout << "e^(pi) = " << std::exp(PI) << "\n";
  std::cout << "e^i = " <<  std::exp(I) << "\n";
  std::cout << "pi^(e) = " << std::pow(PI, std::exp(1)) << "\n";
  std::cout << "pi^i = " << std::pow(PI, I) << "\n";
  std::cout << "i^i = " << std::pow(I, I) << "\n";
  std::cout << "log(i) = " << std::log(I) << "\n\n";

  for(double x=1;x<=10;x+=1)
    std::cout << "fgamma("<< x <<")=" << sfuns::fgamma(x)
      << " tgamma("<< x <<")=" << std::tgamma(x) << "\n";
  double x = 100;
  std::cout << "\nlngamma(100)=" << sfuns::lngamma(x) << " lgamma(100)=" << std::lgamma(x) <<"\n";  
  return 0;
}

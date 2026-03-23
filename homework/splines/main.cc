#include <iostream>
#include <cmath>
#include "core/vector.h"
#include "interp.h"

int main(int argc, char* argv[]){
  bool l=false, q=false, c=false;
  int a = 11;
  pp::vector x(a), y(a);
  for(int i=0; i<a; i++){
    x[i] = 1*i;
    y[i] = std::cos(x[i]) + std::cos(2*x[i] + 0.5);
  }

  for(int i = 0; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "-l" && i<argc) l = true;
    if(arg == "-q" && i<argc) q = true;
    if(arg == "-c" && i<argc) c = true;
  }
  int n = 1000;
  double dx = double(a-1)/(n-1);
  if(l == true){
    pp::lspline spline(x,y);
    for(int i = 0; i<n; i++){
      double x = dx*i;
      double y = spline.eval(x);
      double z = spline.deriv(x);
      double w = spline.integ(x);
      std::cout << x << " " << y << " " << z << " " << w << "\n";
    }
  }
  else if(q == true){
    pp::qspline spline(x,y);
    for(int i = 0; i<n; i++){
      double x = dx*i;
      double y = spline.eval(x);
      double z = spline.deriv(x);
      double w = spline.integ(x);
      std::cout << x << " " << y << " " << z << " " << w << "\n";
    }
  }
  else if(c == true){
    pp::cspline spline(x,y);
    for(int i = 0; i<n; i++){
      double x = dx*i;
      double y = spline.eval(x);
      double z = spline.deriv(x);
      double w = spline.integ(x);
      std::cout << x << " " << y << " " << z << " " << w << "\n";
    }
  }
  else {
    for(int i = 0; i<n; i++){
      double x = dx*i;
      double y = std::cos(x) + std::cos(2*x+0.5);
      double deriv = -std::sin(x) -2*std::sin(2*x+0.5);
      double integ = std::sin(x) + std::sin(2*x+0.5)/2;
      std::cout << x << " " << y << " " << deriv << " " << integ << "\n";
    }
  }
  return 0;
}

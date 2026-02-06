#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

bool approx(double a, double b, double acc=1e-9, double eps=1e-9){
  double diff = std::abs(a-b);
  if(diff <= acc) return true;
  double max_ab = std::max(std::abs(a), std::abs(b));
  return diff <= eps*max_ab;
}

int main() {
  float f = 1.0f;
  while((float) (1.0f+f) != 1.0f){f/=2.0f;} f*=2.0f;
  double d = 1.0d;
  while((1.0d+d) != 1.0d){d/=2.0d;} d*=2.0d;
  long double l = 1.0l;
  while((1.0l+l) != 1.0l){l/=2.0l;} l*=2.0l;

  std::cout << "float eps = "<< f << "\n";
  std::cout << "double eps = "<< d << "\n";
  std::cout << "long double eps = "<< l << "\n\n";

  std::cout << "true float eps = "<< std::numeric_limits<float>::epsilon() << "\n";
  std::cout << "true double eps = "<< std::numeric_limits<double>::epsilon() << "\n";
  std::cout << "true long double eps = "<< std::numeric_limits<long double>::epsilon() << "\n\n";
 
  double epsilon = std::pow(2,-52);
  double tiny = epsilon/2;
  double a=1+tiny+tiny;
  double b=tiny+tiny+1;

  std::cout << "a==b: "<< (a==b ? "true":"false") << "\n";
  std::cout << "a>1: "<< (a>1 ? "true":"false") << "\n";
  std::cout << "b>1: "<< (b>1 ? "true":"false") << "\n";

  std::cout << std::fixed << std::setprecision(17);
  std::cout << "a = " << a << "\n";
  std::cout << "b = "<< b << "\n";
  std::cout << "tiny = "<< tiny << "\n\n";

  double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
  double d2 = 8*0.1;
  std::cout << "d1 == d2: " << (d1 == d2 ? "true":"false") << "\n";
  std::cout << "d1 = " << d1 << "\n";
  std::cout << "d2 = " << d2 << "\n";

  std::cout << "approx(d1, d2, acc=1e-9, eps=1e-9): "<< (approx(d1, d2) ? "true":"false") << "\n";
  return 0;
}

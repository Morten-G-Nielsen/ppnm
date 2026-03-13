#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include "core/vector.h"
#include "ols.h"

int main(){
  std::vector<std::function<double(double)>> fs{[](double z){return 1;},[](double z){return z;}};
  pp::vector x{1,  2,  3, 4, 6, 9,   10,  13,  15};
  pp::vector y{117,100,88,72,53,29.5,25.2,15.2,11.1};
  pp::vector dy{6,5,4,4,4,3,3,2,2};
  for(int i = 0; i<y.size(); i++){
    dy[i] = dy[i]/y[i];
    y[i] = std::log(y[i]);
  }
  auto [coeff, sigma] = pp::lsfit(fs, x, y, dy);
  pp::vector c(fs.size());
  for(int i = 0; i<sigma.row_count(); i++) c[i] = std::sqrt(sigma(i,i));
  std::cout << "a="<<coeff[0]<<"\n"<<"b="<<coeff[1]<<"\n";
  std::cout << "da="<<c[0]<<"\n"<<"db="<<c[1]<<"\n";
  return 0;
}

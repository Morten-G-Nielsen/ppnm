#include <iostream>
#include <vector>
#include <cmath>
#include "core/vector.h"
#include "minimize.h"

using namespace pp;
int main(){
  std::vector<double> energy, signal, error;
  double x, y, z;
  while(std::cin >> x >> y >> z){
    energy.push_back(x); signal.push_back(y); error.push_back(z);
  }
  auto F = [](double E, double m, double gamma, double A){
    return A / ((E-m)*(E-m) + gamma*gamma/4.0);
  };
  auto D = [&](const vector& x){
    double sum_sq = 0;
    if(x[1] < 0.5) sum_sq += 100*std::pow(0.5 - x[1], 2);
    for(size_t i = 0; i<energy.size(); i++){
      sum_sq += std::pow((F(energy[i], x[0], x[1], x[2]) - signal[i])/error[i], 2);
    }
    return sum_sq;
  };

  NewtonCentral n; BacktrackingSearch s; FindMinimum min(n, s);
  vector x0{126, 2, 10};
  auto result = min.minimize(D, x0, 1e-6);
  auto res = result.first;

  double E_min = energy.front(), E_max = energy.back();
  double step = (E_max - E_min)/999.0;
  for(int i = 0; i<1000; i++){
    double E = E_min + step*i;
    std::cout << E << " " << F(E, res[0], res[1], res[2]) << "\n";
  }
  std::cerr << "---Fitting parametres---"
    << "\nMass m: " << res[0] << " +/- "<< result.second[0] 
    << "\nResonance width Gamma: " << res[1] << " +/- " << result.second[1]
    << "\nScale factor A: " << res[2] << " +/- " << result.second[2] << "\n\n";
  return 0;
}

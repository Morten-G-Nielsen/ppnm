#include <iostream>
#include <string>
#include <cmath>
#include "core/vector.h"
#include "ode.h"
#include "ann.h"

using namespace pp;
int main(int argc, char* argv[]){
  int id = 1;
  for(int i = 0; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "-id") id = std::stoi(argv[++i]);
  }

  if(id == 1){
    NeuralNetwork nn_interp(30);
    Interpolator interpolator(nn_interp);

    int n = 50;
    vector x_data(n), y_data(n);
    for(size_t i = 0; i<n; i++){
      x_data[i] = -3.0 + i*6.0/(n-1);
      y_data[i] = std::cos(5.0*x_data[i]-1.0)*std::exp(-x_data[i]*x_data[i]);
    }
    interpolator.train(x_data, y_data);

    int N = 200;
    for(size_t i = 0; i<N; i++){
      double x = -3.0 + i*6.0/(N-1);
      double y = std::cos(5.0*x-1.0)*std::exp(-x*x);
      std::cout << x << " " << y 
        << " " << interpolator.get_value(x)
        << " " << interpolator.get_derivative(x)
        << " " << interpolator.get_second_derivative(x)
        << " " << interpolator.get_definite_integral(-3.0, x) << "\n";
    }
  } else if(id == 2){
    NeuralNetwork nn_ode(20);
    double a=0.0, b=4.0;
    double y = 1.0, y_prime = -0.1;
    ODE_Solver solver(nn_ode, a, b, a, y, y_prime, 100, 100);

    auto pendulum_ann = [](double d2y, double dy, double y, double x){
      return d2y + std::sin(y);
    };
    auto pendulum_rk23 = [](double x, vector y){
      vector dydx(2);
      dydx[0] = y[1];
      dydx[1] = -std::sin(y[0]);
      return dydx;
    };
    solver.train(pendulum_ann, 1e-6);
    rk_45 stepper;
    vector yinit = {y, y_prime};
    auto res_rk23 = driver(stepper, pendulum_rk23, a, b, yinit, 1e-3, 1e-6, 1e-6); 
    auto x_data = res_rk23.x_values;
    auto y_data = res_rk23.y_values;

    for(size_t i = 0; i<100; i++){
      double x = i*4.0/99.0;
      if(i<x_data.size()){
        std::cout << x << " "
          << solver.get_solution(x) << " "
          << x_data[i] << " "
          << y_data[i][0] << " " 
          << y_data[i][1] << "\n";
      } else{
        std::cout << x << " " << solver.get_solution(x) << "\n";
      }
    }
  }
  return 0;
}

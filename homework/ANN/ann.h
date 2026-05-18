#pragma once
#include <functional>
#include "core/vector.h"

namespace pp{

struct NeuralNetwork{
  size_t n;
  std::vector<std::vector<double>> params;

  NeuralNetwork(size_t n);

  vector get_params()const; 
  void set_params(const vector& p);
  double response(double x)const;
  double derivative(double x)const;
  double second_derivative(double x)const;
  double integrate(double x)const;
};

struct Interpolator{
  NeuralNetwork& nn;
  double x_min, x_max;
  double y_min, y_max;
  double dx_domain, dy_domain;

  Interpolator(NeuralNetwork& n);

  double scale_x(double x)const;
  double scale_y(double y)const;
  double unscale_y(double y_scaled)const;

  void train(const vector& x_data, const vector& y_data, double acc=1e-6);
  double get_value(double x)const;
  double get_derivative(double x)const;
  double get_second_derivative(double x)const;
  double get_definite_integral(double a, double b)const;
};

struct ODE_Solver{
  NeuralNetwork& nn;
  double a, b, c;
  double yc, yc_prime;
  double alpha, beta;

  double dx_domain;
  double scale_factor_d1;
  double scale_factor_d2;

  ODE_Solver(NeuralNetwork& n, double start, double end, double anchor, double target_y,
      double target_dy, double alpha=1000.0, double beta=1000.0);

  double scale_x(double x)const;

  void train(const std::function<double(double y_d2, double y_d1, double y, double x)>& phi,
      double acc=1e-3);
  double get_solution(double x)const;
};
}

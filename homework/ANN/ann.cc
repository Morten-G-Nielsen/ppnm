#include <cmath>
#include <vector>
#include "core/vector.h"
#include "quad.h"
#include "minimize.h"
#include "ann.h"

namespace pp{

NeuralNetwork::NeuralNetwork(size_t n) : n(n) {
  params.assign(n, std::vector<double>(3));
  for(size_t i = 0; i<n; i++){
    params[i][0] = 0.1;
    params[i][1] = -0.5 + 1.0*double(i)/(n-1);
    params[i][2] = 0.3;
  }
}

vector NeuralNetwork::get_params()const{
  vector p(3*n);
  size_t idx = 0;
  for(size_t i = 0; i<n; i++){
    p[idx++] = params[i][0];
    p[idx++] = params[i][1];
    p[idx++] = params[i][2];
  }
  return p;
}

void NeuralNetwork::set_params(const vector& p){
  size_t idx = 0;
  for(size_t i = 0; i<n; i++){
    params[i][0] = p[idx++];
    params[i][1] = p[idx++];
    params[i][2] = p[idx++];
  }
}

double NeuralNetwork::response(double x)const{
  double sum = 0.0;
  for(const auto& p : params){
    double w = p[0], a = p[1], b = p[2];
    double z = (x - a)/b;

    sum += w*z*std::exp(-z*z);
  }
  return sum;
}

double NeuralNetwork::derivative(double x)const{
  double sum = 0.0;
  for(const auto& p : params){
    double w = p[0], a = p[1], b = p[2];
    double z = (x - a)/b;

    double df_dz = (1.0 - 2.0*z*z)*std::exp(-z*z);
    sum += (w/b) * df_dz;
  }
  return sum;
}

double NeuralNetwork::second_derivative(double x)const{
  double sum = 0.0;
  for(const auto& p : params){
    double w = p[0], a = p[1], b = p[2];
    double z = (x - a)/b;

    double d2f_dz = (4*z*z - 6)*z*std::exp(-z*z);
    sum += w*d2f_dz/b/b;
  }
  return sum;
}

double NeuralNetwork::integrate(double x)const{
  double sum = 0.0;
  for(const auto& p : params){
    double w = p[0], a = p[1], b = p[2];
    double z = (x - a)/b;

    double F_integral = -0.5*std::exp(-z*z);
    sum += w*b*F_integral;
  }
  return sum;
}

Interpolator::Interpolator(NeuralNetwork& network) : nn(network),
  x_min(INFINITY), x_max(-INFINITY), y_min(INFINITY), y_max(-INFINITY) {}

double Interpolator::scale_x(double x)const{return 2.0*(x-x_min)/dx_domain - 1.0;}
double Interpolator::scale_y(double y)const{return 2.0*(y-y_min)/dy_domain - 1.0;}
double Interpolator::unscale_y(double y_scaled)const{return (y_scaled+1.0)/2.0*dy_domain+y_min;}

void Interpolator::train(const vector& x_data, const vector& y_data, double acc){
  x_min = INFINITY; x_max = -INFINITY;
  y_min = INFINITY; y_max = -INFINITY;

  for(size_t i = 0; i<x_data.size(); i++){
    double x_val = x_data[i], y_val = y_data[i];
    if(x_min>x_val) x_min = x_val;
    if(x_max<x_val) x_max = x_val;
    if(y_min>y_val) y_min = y_val;
    if(y_max<y_val) y_max = y_val;
  }
  dx_domain = x_max - x_min;
  dy_domain = y_max - y_min;

  auto cost = [&](const vector& p){
    nn.set_params(p);
    double sum_sq = 0.0;

    for(size_t i = 0; i<x_data.size(); i++){
      double xs = scale_x(x_data[i]);
      double ys = scale_y(y_data[i]);
      sum_sq += std::pow(nn.response(xs) - ys, 2);
    }
    return sum_sq;
  };
  vector p0 = nn.get_params();
  BFGS optimizer(p0.size());
  BacktrackingSearch line_search;
  FindMinimum solver(optimizer, line_search);
  auto res = solver.minimize(cost, p0, acc);
  nn.set_params(res.first);
}
double Interpolator::get_value(double x)const{
  return unscale_y(nn.response(scale_x(x)));
}
double Interpolator::get_derivative(double x)const{
  return nn.derivative(scale_x(x))*(dy_domain/dx_domain);
}
double Interpolator::get_second_derivative(double x)const{
  return nn.second_derivative(scale_x(x))*(dy_domain/(dx_domain*dx_domain));
}
double Interpolator::get_definite_integral(double a, double b)const{
  double start = nn.integrate(scale_x(a));
  double end = nn.integrate(scale_x(b));
  return (end-start)*(dx_domain*dy_domain/4.0);
}

ODE_Solver::ODE_Solver(NeuralNetwork& n, double start, double end, double anchor,
    double target_y, double target_dy, double alpha, double beta) : nn(n),
  a(start), b(end), c(anchor), yc(target_y), yc_prime(target_dy), alpha(alpha), beta(beta) {
    dx_domain = b - a;
    scale_factor_d1 = 2.0/dx_domain;
    scale_factor_d2 = 4.0/(dx_domain*dx_domain);
  }

double ODE_Solver::scale_x(double x)const{return 2.0*(x-a)/dx_domain - 1.0;}

void ODE_Solver::train(
    const std::function<double(double y_d2, double y_d1, double y, double x)>& phi,
    double acc){
  NewtonCotes nc;
  Integrator quad(nc);
  auto f = [&](double x){
    double y = nn.response(scale_x(x));
    double dy_dx = nn.derivative(scale_x(x))*scale_factor_d1;
    double d2y_dx2 = nn.second_derivative(scale_x(x))*scale_factor_d2;
    return std::pow(phi(d2y_dx2, dy_dx, y, x), 2);
  };
  auto cost = [&](const vector& p){
    nn.set_params(p);
    auto quad_res = quad.integrate(f, a, b);

    double xs_c = scale_x(c);
    double boundary__pos_error = std::pow(nn.response(xs_c) - yc, 2);
    double boundary_deriv_error = std::pow(nn.derivative(xs_c)*scale_factor_d1 - yc_prime, 2);
    
    return quad_res.integral + alpha*boundary__pos_error + beta*boundary_deriv_error;
  };
  vector p0 = nn.get_params();
  BFGS optimizer(p0.size());
  BacktrackingSearch line_search;
  FindMinimum solver(optimizer, line_search);

  auto result = solver.minimize(cost, p0, acc);
  nn.set_params(result.first);
}
double ODE_Solver::get_solution(double x)const{
  return nn.response(scale_x(x));
}
}

#include <cmath>
#include <functional>
#include "quad.h"

namespace pp{
QuadResult Integrator::integrate(const std::function<double(double)>& f, double a, double b){
  auto trans_f = rule.get_function_transform(f, a, b);
  auto [start, end] = rule.get_limits(a, b);

  std::vector<double> fs;
  return recursive_step(trans_f, start, end, acc, 0, fs);
}
QuadResult Integrator::recursive_step(const std::function<double(double)>& f,
    double a, double b, double local_acc, int depth, std::span<const double> view){

  std::vector<double> fs;

  auto [Q, q, evals] = rule.estimate(f, a, b, view, fs);
  double err = std::abs(Q - q);

  if(err < local_acc + eps*std::abs(Q) || depth >= max_depth){
    return {Q, err, evals};
  }

  double mid = a + (b-a) / 2.0;
  local_acc /= std::sqrt(2);

  auto left_view = rule.get_left_points(fs);
  auto right_view = rule.get_right_points(fs);

  auto left = recursive_step(f, a, mid, local_acc, depth+1, left_view);
  auto right = recursive_step(f, mid, b, local_acc, depth+1, right_view);

  return{
    left.integral+right.integral,
    std::sqrt(left.error*left.error+right.error*right.error),
    evals + left.evals + right.evals
  };
}

EstimateResult NewtonCotes::estimate(const std::function<double(double)>& f,
    double a, double b, std::span<const double> view, std::vector<double>& fs){
  double h = b - a;
  double f1, f2, f3, f4;
  int nfev = 0;

  if(view.empty()){
    f2 = f(a + 2*h/6);
    f3 = f(a + 4*h/6);
    nfev += 2;
  } else {
    f2 = view[0];
    f3 = view[1];
  }
  f1 = f(a + h/6);
  f4 = f(a + 5*h/6);
  nfev += 2;

  fs = {f1, f2, f3, f4};
  double Q = (2*f1+f2+f3+2*f4)/6*h;
  double q = (f1+f2+f3+f4)/4*h;
  return {Q, q, nfev};
}
std::span<const double> NewtonCotes::get_left_points(std::span<const double> fs) const {
  return std::span(fs).subspan(0, 2);
}
std::span<const double> NewtonCotes::get_right_points(std::span<const double> fs) const {
  return std::span(fs).subspan(2,2);
}

EstimateResult ClenshawCurtis::estimate(const std::function<double(double)>& g,
    double a, double b, std::span<const double> view, std::vector<double>& gs){
  double h = b - a;
  double g1, g2, g3, g4;
  int nfev = 0;

  if(view.empty()){
    g2 = g(a + 2*h/6);
    g3 = g(a + 4*h/6);
    nfev += 2;
  } else {
    g2 = view[0];
    g3 = view[1];
  }
  g1 = g(a + h/6);
  g4 = g(a + 5*h/6);
  nfev += 2;

  gs = {g1, g2, g3, g4};
  double Q = (2*g1+g2+g3+2*g4)/6*h;
  double q = (g1+g2+g3+g4)/4*h;
  return {Q, q, nfev};
}
std::pair<double, double> ClenshawCurtis::get_limits(double a, double b) const {
  return {0, M_PI};
}
std::function<double(double)> ClenshawCurtis::get_function_transform(const std::function<double(double)>& f, double a, double b) const {
  auto g = [=](double theta){
    double x = (a + b)/2.0 + (b - a)/2.0*std::cos(theta);
    double weight = (b - a)/2.0*std::sin(theta);
    return f(x)*weight;
  };
  return g;
}
std::span<const double> ClenshawCurtis::get_left_points(std::span<const double> gs) const {
  return std::span(gs).subspan(0,2);
}
std::span<const double> ClenshawCurtis::get_right_points(std::span<const double> gs) const {
  return std::span(gs).subspan(2,2);
}

EstimateResult InfiniteRule::estimate(const std::function<double(double)>& f,
    double a, double b, std::span<const double> view, std::vector<double>& fs){
  return base_rule.estimate(f, a, b, view, fs);
}
std::span<const double> InfiniteRule::get_left_points(std::span<const double> fs) const {
  return base_rule.get_left_points(fs);
}
std::span<const double> InfiniteRule::get_right_points(std::span<const double> fs) const {
  return base_rule.get_right_points(fs);
}
std::pair<double, double> InfiniteRule::get_limits(double a, double b) const {
  if(std::isinf(b) && std::isinf(a)){
    return base_rule.get_limits(-1.0, 1.0);
  } else if(std::isinf(b)) {
    return base_rule.get_limits(0.0, 1.0);
  } else if(std::isinf(a)) {
    return base_rule.get_limits(-1.0, 0.0);
  } else return base_rule.get_limits(a, b);
}
std::function<double(double)> InfiniteRule::get_function_transform(const std::function<double(double)>& f, double a, double b) const {
  if(std::isinf(b) && std::isinf(a)){
    auto g = [=](double t){
     double x = t/(1-t*t); 
     double jacobian = (1+t*t)/((1-t*t)*(1-t*t));
     return f(x)*jacobian;
    };
    return base_rule.get_function_transform(g, -1.0, 1.0);
  } else if(std::isinf(b)){
    auto g = [=](double t){
      double x = a + (1-t)/t;
      double jacobian = 1/t/t;
      return f(x)*jacobian;
    };
    return base_rule.get_function_transform(g, 0.0, 1.0);
  } else if(std::isinf(a)){
    auto g = [=](double t){
      double x = b - (1-t)/t;
      double jacobian = 1/t/t;
      return f(x)*jacobian;
    };
    return base_rule.get_function_transform(g, 0.0, 1.0);
  } else return base_rule.get_function_transform(f, a, b);
}
}

#pragma once
#include <span>
#include <vector>
#include <functional>

namespace pp{
struct QuadResult{
  double integral;
  double error;
  int evals;
};
struct EstimateResult{
  double Q;
  double q;
  int nfev;
};
struct QuadRule{
  virtual EstimateResult estimate(const std::function<double(double)>& f,
      double a, double b, std::span<const double> view, std::vector<double>& fs) = 0;

  virtual std::span<const double> get_left_points(std::span<const double> fs) const = 0;
  virtual std::span<const double> get_right_points(std::span<const double> fs) const = 0;
  virtual std::pair<double, double> get_limits(double a, double b) const{
    return {a, b};
  }
  virtual std::function<double(double)> get_function_transform(const std::function<double(double)>& f,
      double a, double b) const {
    return f;
  }

  virtual ~QuadRule()=default;
};

struct NewtonCotes : public QuadRule{
  EstimateResult estimate(const std::function<double(double)>& f,
      double a, double b, std::span<const double> view, std::vector<double>& fs) override;

  std::span<const double> get_left_points(std::span<const double> fs) const override;
  std::span<const double> get_right_points(std::span<const double> fs) const override;
};
struct ClenshawCurtis : public QuadRule{
  EstimateResult estimate(const std::function<double(double)>& f,
      double a, double b, std::span<const double> view, std::vector<double>& fs) override;

  std::span<const double> get_left_points(std::span<const double> fs) const override;
  std::span<const double> get_right_points(std::span<const double> fs) const override;
  std::pair<double, double> get_limits(double a, double b) const override;
  std::function<double(double)> get_function_transform(const std::function<double(double)>& f, double a, double b) const override;
};
struct InfiniteRule : public QuadRule{
  QuadRule& base_rule;

  InfiniteRule(QuadRule& rule) : base_rule(rule){};

  EstimateResult estimate(const std::function<double(double)>& f,
      double a, double b, std::span<const double> view, std::vector<double>& fs) override;
  std::span<const double> get_left_points(std::span<const double> fs) const override;
  std::span<const double> get_right_points(std::span<const double> fs) const override;
  std::pair<double, double> get_limits(double a, double b) const override;
  std::function<double(double)> get_function_transform(const std::function<double(double)>& f, double a, double b) const override;
};

struct Integrator{
  QuadRule& rule;
  double acc = 1e-3;
  double eps = 1e-3;
  int max_depth = 25;

  QuadResult integrate(const std::function<double(double)>& f, double a, double b);
  QuadResult recursive_step(const std::function<double(double)>& f, double a, double b,
      double local_acc, int depth, std::span<const double> fs);
};
}

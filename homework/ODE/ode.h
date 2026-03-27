#pragma once
#include <vector>
#include <functional>
#include "core/vector.h"

namespace pp{
struct rk_step_res{
  vector y_new;
  vector error;
};

struct ode_res{
  std::vector<double> x_values;
  std::vector<vector> y_values;
};

struct rk_stepper{
  virtual rk_step_res step(
      const std::function<vector(double, vector)>& f,
      double x,
      const vector& y,
      double h) = 0;
  virtual ~rk_stepper()=default;
};

struct rk_12 : public rk_stepper{
  rk_step_res step(
      const std::function<vector(double, vector)>& f,
      double x, const vector& y, double h) override;
};

struct rk_45 : public rk_stepper{
  rk_step_res step(
      const std::function<vector(double, vector)>& f,
      double x, const vector& y, double h) override;
};

ode_res driver(
    rk_stepper& stepper,
    std::function<vector(double, vector)> f,
    double a, double b,
    vector y_start,
    double h,
    double acc,
    double eps
    );
}

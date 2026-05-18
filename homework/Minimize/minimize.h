#pragma once
#include <functional>
#include "core/vector.h"
#include "core/matrix.h"


namespace pp{
struct UpdateRule{
  virtual ~UpdateRule() = default;
  virtual vector get_update(const std::function<double(const vector&)>& f,
      vector& x, double fx, vector& grad, matrix& Hessian) = 0;
};
struct SearchRule{
  virtual ~SearchRule() = default;
  virtual double get_step(const std::function<double(const vector&)>& f,
      const vector& x, double fx, const vector& dx, double alpha_min=1e-3) = 0;
};

struct BacktrackingSearch : SearchRule{
  double get_step(const std::function<double(const vector&)>& f,
      const vector& x, double fx, const vector& dx, double alpha_min) override;
};

struct NewtonForward : UpdateRule{
  vector get_update(const std::function<double(const vector&)>& f,
      vector& x, double fx, vector& grad, matrix& Hessian) override;
};

struct NewtonCentral : UpdateRule{
  vector get_update(const std::function<double(const vector&)>& f,
      vector& x, double fx, vector& grad, matrix& Hessian) override;
};

struct BFGS : UpdateRule{
  vector last_x;
  vector last_grad;
  bool is_first;
  BFGS(size_t dim) : last_x(dim), last_grad(dim), is_first(true) {}

  vector get_update(const std::function<double(const vector&)>& f,
      vector& x, double fx, vector& grad, matrix& Hessian) override;
};

struct FindMinimum{
  UpdateRule& update_rule;
  SearchRule& search_rule;

  FindMinimum(UpdateRule& u, SearchRule& s) : update_rule(u), search_rule(s) {}

  std::pair<vector, vector> minimize(const std::function<double(const vector&)>& f,
      const vector& x0, double acc, size_t max_iter = 1000);
};
}

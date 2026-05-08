#pragma once
#include <functional>
#include "core/vector.h"
#include "core/matrix.h"

namespace pp{
struct UpdateRule{
  virtual ~UpdateRule() = default;
  virtual vector get_update(const std::function<vector(const vector&)>& f,
      const vector& x, const vector& fx) = 0;

  void fill_jacobian(const std::function<vector(const vector&)>& f,
      const vector& x, const vector& fx, matrix& J);
};

struct SearchRule{
  virtual ~SearchRule() = default;
  virtual double get_step(const std::function<vector(const vector&)>& f,
      const vector& x, const vector& fx, const vector& dx, double alpha_min) = 0;
};

struct NewtonUpdate : UpdateRule{
  matrix J;
  NewtonUpdate(size_t dim) : J(dim, dim) {}

  vector get_update(const std::function<vector(const vector&)>& f,
      const vector& x, const vector& fx) override;
};

struct BroydenUpdate : UpdateRule{
  matrix B;
  vector last_x;
  vector last_fx;
  bool is_first;
  BroydenUpdate(size_t dim) : B(dim, dim), last_x(dim), last_fx(dim), is_first(true) {}

  vector get_update(const std::function<vector(const vector&)>& f,
      const vector& x, const vector& fx) override;

  void reset(){is_first = true;}
};

struct BacktrackingSearch : SearchRule{
  double get_step(const std::function<vector(const vector&)>& f,
      const vector& x, const vector& fx, const vector& dx, double alpha_min) override;
};

struct QuadraticInterpSearch : SearchRule{
  double lambda_trial;
  double phi_trial;

  QuadraticInterpSearch(double l, double p) : lambda_trial(l), phi_trial(p) {}

  double get_step(const std::function<vector(const vector&)>& f,
      const vector& x, const vector& fx, const vector& dx, double alpha_min) override;
};

struct FindRoot{
  UpdateRule& update_rule;
  SearchRule& search_rule;

  FindRoot(UpdateRule& u, SearchRule& s) : update_rule(u), search_rule(s) {}

  vector solve(const std::function<vector(const vector&)>& f, vector x0, double acc,
      double alpha_min, size_t max_iter=100);
};
}

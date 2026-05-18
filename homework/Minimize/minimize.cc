#include <cmath>
#include <functional>
#include "core/vector.h"
#include "core/matrix.h"
#include "minimize.h"
#include "qr.h"

namespace pp{

double BacktrackingSearch::get_step(const std::function<double(const vector&)>& f,
    const vector& x, double fx, const vector& dx, double alpha_min){
  double alpha = 1.0;
  vector z = x;
  while(alpha >= alpha_min){
    z = x + dx*alpha;
    double fz = f(z);

    if(fz < fx){
      return alpha;
    }
    alpha/=2;
  }
  return alpha;
}

vector NewtonForward::get_update(const std::function<double(const vector&)>& f,
    vector& x, double fx, vector& grad, matrix& Hessian){
  size_t n = x.size();
  for(size_t i = 0; i<n; i++){
    double old_x = x[i];
    double dx = (1.0 + std::abs(x[i]))*std::pow(2, -26);
    x[i] += dx;
    grad[i] = (f(x) - fx)/dx;
    x[i] = old_x;
  }

  double eps13 = std::pow(2, -13);
  vector fx_steps(n);
  for(size_t k = 0; k<n; k++){
    double old_xk = x[k];
    x[k] += (1.0 + std::abs(old_xk))*eps13;
    fx_steps[k] = f(x);
    x[k] = old_xk;
  }

  for(size_t j = 0; j<n; j++){
    double old_xj = x[j];
    double dxj = (1.0 + std::abs(x[j]))*eps13;
    x[j] += 2.0*dxj;
    double fx_plus_plus_j = f(x);
    Hessian(j, j) = (fx_plus_plus_j - 2.0*fx_steps[j] + fx)/(dxj*dxj) + 1e-6;
    x[j] = old_xj;

    for(size_t i = j+1; i<n; i++){
      double old_xi = x[i];
      double dxi = (1.0 + std::abs(old_xi))*eps13; 

      x[i] += dxi;
      x[j] += dxj;
      double fx_plus_ij = f(x);
      Hessian(i, j) = (fx_plus_ij - fx_steps[i] - fx_steps[j] + fx)/(dxi*dxj);
      Hessian(j, i) = Hessian(i, j);
      x[j] = old_xj;
      x[i] = old_xi;
    }
  }
  QR qr(Hessian);
  return qr.solve(-1.0*grad);
}

vector NewtonCentral::get_update(const std::function<double(const vector&)>& f,
    vector& x, double fx, vector& grad, matrix& Hessian){
  size_t n = x.size();
  double eps13 = std::pow(2, -17);

  vector f_plus(n), f_minus(n);
  for(size_t i = 0; i<n; i++){
    double old_x = x[i];
    double dx = (1.0 + std::abs(old_x))*eps13;

    x[i] += dx; f_plus[i] = f(x);
    x[i] = old_x - dx; f_minus[i] = f(x);
    x[i] = old_x;

    grad[i] = (f_plus[i] - f_minus[i]) / (2.0*dx);
    Hessian(i, i) = (f_plus[i] -2.0*fx + f_minus[i])/(dx*dx) + 1e-6;
  }

  for(size_t j = 0; j<n; j++){
    double old_xj = x[j];
    double dxj = (1.0 + std::abs(old_xj))*eps13;
    for(size_t i = j+1; i<n; i++){
      double old_xi = x[i];
      double dxi = (1.0 + std::abs(x[i]))*eps13;

      x[i] += dxi;
      x[j] += dxj;
      double f_pp = f(x);

      x[j] = old_xj - dxj;
      double f_pm = f(x);

      x[i] = old_xi - dxi;
      double f_mm = f(x);

      x[j] = old_xj + dxj;
      double f_mp = f(x);

      Hessian(i, j) = (f_pp - f_pm - f_mp + f_mm)/(4.0*dxi*dxj);
      Hessian(j, i) = Hessian(i, j);

      x[i] = old_xi;
      x[j] = old_xj;
    }
  }
  QR qr(Hessian);
  return qr.solve(-1.0*grad);
}

vector BFGS::get_update(const std::function<double(const vector&)>& f,
    vector& x, double fx, vector& grad, matrix& H){
  size_t n = x.size();
  double eps13 = std::pow(2, -17);
  double f_p = 0.0;
  double f_m = 0.0;
  for(size_t i = 0; i<n; i++){
    double old_x = x[i];
    double dx = (1 + std::abs(old_x))*eps13;

    x[i] += dx; f_p = f(x);
    x[i] = old_x - dx; f_m = f(x);
    x[i] = old_x;

    grad[i] = (f_p - f_m)/(2.0*dx);
  }

  if(is_first){
    is_first = false;
    last_x = x;
    last_grad = grad;
    return -(H*grad);
  }

  vector s = x - last_x;
  vector y = grad - last_grad;
  double rho = 1/dot(y, s);

  vector Hy = H*y;
  double yHy = dot(y, Hy);

  for(size_t i = 0; i<n; i++){
    for(size_t j = 0; j<n; j++){
      H(i,j) = H(i,j) - rho*(s[i]*Hy[j] + Hy[i]*s[j]) + (rho*rho*yHy+rho)*(s[i]*s[j]);
    }
  }
  last_x = x;
  last_grad = grad;
  return -(H*grad);
}

std::pair<vector, vector> FindMinimum::minimize(
    const std::function<double(const vector&)>& f,
    const vector& x0, double acc, size_t max_iter){
  size_t n = x0.size();
  vector x = x0;
  vector err(x);
  vector grad(n);
  matrix H = matrix::identity(n);
  for(size_t i = 0; i<max_iter; i++){
    double fx = f(x);
    vector dx = update_rule.get_update(f, x, fx, grad, H);

    if(grad.norm() < acc){
      QR qr(H);
      matrix Hinv = qr.inverse();
      for(size_t j = 0; j<n; j++) err[j] = std::sqrt(2.0*Hinv(j,j));
      return {x, err};
    }

    double alpha = search_rule.get_step(f, x, fx, dx);
    x += dx*alpha;
  }
  QR qr(H);
  matrix Hinv = qr.inverse();
  for(size_t i = 0; i<n; i++){
    err[i] = std::sqrt(2.0*Hinv(i,i));
  }
  return {x, err};
}
}

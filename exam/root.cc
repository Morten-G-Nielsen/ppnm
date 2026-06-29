#include <cmath>
#include <functional>
#include <stdexcept>
#include "core/vector.h"
#include "core/matrix.h"
#include "qr.h"
#include "root.h"

namespace pp{

void UpdateRule::fill_jacobian(const std::function<vector(const vector&)>& f,
    const vector& x, const vector& fx, matrix& J){
  size_t dim = x.size();
  vector x_plus = x;
  vector x_minus = x;

  for(size_t j = 0; j<dim; j++){
    double h = std::max(x[j], 1.0)*std::pow(2.0, -26);
    x_plus = x; x_plus[j] += h;
    x_minus = x; x_minus[j] -= h;

    vector f_plus = f(x_plus);
    vector f_minus = f(x_minus);

    for(size_t i = 0; i<dim; i++){
      J(i, j) = (f_plus[i] - f_minus[i])/(2.0*h);
    }
  }
}

vector NewtonUpdate::get_update(const std::function<vector(const vector&)>& f,
       const vector& x, const vector& fx){
  fill_jacobian(f, x, fx, J); 
  QR qr(J);
  return qr.solve(-fx);
}

vector GoodBroydenUpdate::get_update(const std::function<vector(const vector&)>& f,
    const vector&x, const vector& fx){
  size_t dim = x.size();
  if(is_first){
    fill_jacobian(f, x, fx, B);
    is_first = false;
    last_x = x;
    last_fx = fx;

    QR qr(B);
    B = qr.inverse();
    return B*(-1.0*fx);
  }
  vector dx = x - last_x;
  vector df = fx - last_fx;

  if(force_recalculate){
    fill_jacobian(f, x, fx, B);
    try{
      QR qr(B);
      B = qr.inverse();
    }
    catch(std::runtime_error){
      B = B.identity(B.col_count());
    }
    force_recalculate = false;
    return B*(-1.0*fx);
  }

  vector u = dx - B*df;
  vector w(dim);
  for(size_t i = 0; i<dim; i++){
    double sum = 0.0;
    for(size_t j = 0; j<dim; j++){
      sum += dx[j]*B(j,i);
    }
    w[i] = sum;
  }

  double denom = dot(w, df);

  for(size_t i = 0; i<dim; i++){
    for(size_t j = 0; j<dim; j++){
      B(i,j) += (u[i] * w[j])/denom;
    }
  }
  last_x = x;
  last_fx = fx;
  return B*(-1.0*fx);
}

vector BadBroydenUpdate::get_update(const std::function<vector(const vector&)>& f,
    const vector& x, const vector& fx){
  size_t dim = x.size();
  if(is_first){
    fill_jacobian(f, x, fx, B);
    is_first = false;
    last_x = x;
    last_fx = fx;

    QR qr(B);
    B = qr.inverse();
    return B*(-1.0*fx);
  }
  vector dx = x - last_x;
  vector df = fx - last_fx;

  if(force_recalculate){
    fill_jacobian(f, x, fx, B);
    try {
      QR qr(B);
      B = qr.inverse();
    }
    catch(std::runtime_error){
      B.identity(B.col_count());
    }

    force_recalculate = false;
    return B*(-1.0*fx);
  }

  vector u = dx - B*df;
  double denom = df.norm_sq();
  for(size_t i = 0; i<dim; i++){
    for(size_t j = 0; j<dim; j++){
      B(i,j) += (u[i]*df[j])/denom;
    }
  }
  last_x = x;
  last_fx = fx;
  return B*(-1.0*fx);
}

double BacktrackingSearch::get_step(const std::function<vector(const vector&)>& f,
    const vector& x, const vector& fx, const vector& dx, double alpha_min){
  double alpha = 1.0;
  while(true){
    vector z = x + dx*alpha;
    vector fz = f(z);
    if(fz.norm() < fx.norm()) break;
    if(alpha < alpha_min) break;
    alpha/=2;
  }
  return alpha;
};

double QuadraticInterpSearch::get_step(const std::function<vector(const vector&)>& f,
    const vector& x, const vector& fx, const vector& dx, double alpha_min){

  double alpha = 1.0;
  double phi0 = 0.5*fx.norm_sq();
  double phi_prime0 = -fx.norm_sq();
  while(true){
    vector z = x + dx*alpha;
    vector fz = f(z);
    double phi_alpha = 0.5*fz.norm_sq();
    if(phi_alpha < phi0){
      return alpha;
    }
    if (alpha < alpha_min){
      return alpha;
    }
    double c = (phi_alpha - phi0 - phi_prime0*alpha)/(alpha*alpha);
    double alpha_next;
    if(c <= 0){
      alpha_next = alpha*0.5;
    } else{
      alpha_next = -phi_prime0/(2.0*c);
    }
    if(alpha_next < 0.1*alpha) alpha_next = 0.1*alpha;
    if(alpha_next > 0.5*alpha) alpha_next = 0.5*alpha;
    alpha = alpha_next;
  }
}

vector FindRoot::solve(const std::function<vector(const vector&)>& f, vector x0, double acc,
    double alpha_min, size_t max_iter, bool recalc){
  vector x = x0;
  for(size_t i = 0; i<max_iter; i++){
    vector fx = f(x);
    if(fx.norm() < acc) return x;

    vector dx = update_rule.get_update(f, x, fx);
    double alpha = search_rule.get_step(f, x, fx, dx, alpha_min);

    if(alpha < 1.0/128.0 && recalc) update_rule.set_force_recalculate(true);

    x += (dx*alpha);
  }
  return x;
}

size_t FindRoot::solve_count_iters(const std::function<vector(const vector&)>& f, vector x0,
    double acc, double alpha_min, size_t max_iter, bool recalc){
  vector x = x0;
  for(size_t i = 0; i<max_iter; i++){
    vector fx = f(x);
    if(fx.norm() < acc) return i;

    vector dx = update_rule.get_update(f, x, fx);
    double alpha = search_rule.get_step(f, x, fx, dx, alpha_min);

    if(alpha < 1.0/128.0 && recalc) update_rule.set_force_recalculate(true);

    x += (dx*alpha);
  }
  return max_iter+20;
}

}

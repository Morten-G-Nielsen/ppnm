#include <iostream>
#include <cmath>
#include <functional>
#include "core/vector.h"
#include "ode.h"

namespace pp{
rk_step_res rk_12::step(
    const std::function<vector(double, vector)>& f,
    double x, const vector& y, double h){
  vector k0 = f(x, y);
  vector k1 = f(x+h/2, y+k0*(h/2));
  vector yh = y + k1*h;
  vector dy = (k1-k0)*h;
  return {yh, dy};
}
rk_step_res rk_45::step(
    const std::function<vector(double, vector)>& f,
    double x, const vector& y, double h){
  vector k0 = f(x,y);
  vector k1 = f(x + h/4.0, y + k0/4.0*h);
  vector k2 = f(x + 3.0/8.0*h, y + (3.0/32.0*k0 + 9.0/32.0*k1)*h);
  vector k3 = f(x + 12.0/13.0*h, y + (1932.0/2197.0*k0 - 7200.0/2197.0*k1 + 7296.0/2197.0*k2)*h);
  vector k4 = f(x + h, y + (439.0/216.0*k0 - 8.0*k1 + 3680.0/513.0*k2 - 845.0/4104.0*k3)*h);
  vector k5 = f(x + h/2.0, y + (-8.0/27.0*k0 + 2*k1 - 3544.0/2565.0*k2 + 1859.0/4104.0*k3 - 11.0/40.0*k4)*h);
  vector y_next = y + (16.0/135.0*k0 + 6656.0/12825.0*k2 + 28561.0/56430.0*k3 - 9.0/50.0*k4 + 2.0/55.0*k5)*h;
  vector y_lower = y + (25.0/216.0*k0 + 1408.0/2565.0*k2 + 2197.0/4104.0*k3 - k4/5.0)*h;
  vector dy = y_next - y_lower;
  return {y_next, dy};
}

ode_res driver(
    rk_stepper& stepper,
    std::function<vector(double, vector)> f,
    double a, double b,
    vector y_start,
    double h,
    double acc,
    double eps){
    double x = a; vector y(y_start);
    std::vector<double> xlist; xlist.push_back(x);
    std::vector<vector> ylist; ylist.push_back(y);
    while(true){
      if(x>=b) return {xlist, ylist};
      if(x+h>b) h = b-x;
      auto [yh, dy] = stepper.step(f, x, y, h);
      double tol = (acc + eps*yh.norm())*std::sqrt(h/(b-a));
      double err = dy.norm();
      if(err<=tol){
        x+=h; y=yh;
        xlist.push_back(x);
        ylist.push_back(y);
      }
      double factor = std::pow(tol/err, 0.25)*0.95;
      if(err>0 && factor < 2) h*= factor;
      else h*=2;
    };
}
}

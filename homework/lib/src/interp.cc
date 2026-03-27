#include <cassert>
#include "core/vector.h"
#include "interp.h"

namespace pp{
int binsearch(const vector& x, double z){
  assert(z>=x[0] && z<=x[x.size()-1]);
  int i = 0;
  int j = x.size() - 1;
  while(j-i>1){
    int mid = (i+j)/2;
    if(z>x[mid]) i = mid; else j = mid;
  }
  return i;
}
lspline::lspline(const vector& x, const vector& y) : n(x.size()), x(x), y(y), cumint(n) {
  // Calculates a cumalitive integral
  for(int j = 0; j<n-1; j++){
    double dx = x[j+1]-x[j];
    double y_mean = (y[j+1]+y[j])/2;
    cumint[j+1] = cumint[j] + dx * y_mean; 
  };
};
int lspline::find_bin(double z)const{
  // Remembers the last index to skip a binarysearch if z is in the reagion
  if(z >= x[last_index] && z <= x[last_index+1]) return last_index;
  last_index = binsearch(x, z);
  return last_index;
};
double lspline::eval(double z)const{
  int i = find_bin(z);
  double dx = x[i+1]-x[i];
  double dy = y[i+1]-y[i];
  return y[i] + dy/dx * (z-x[i]);
};
double lspline::deriv(double z)const{
  int i = find_bin(z);
  double dx = x[i+1]-x[i];
  double dy = y[i+1]-y[i];
  return dy/dx;
};
double lspline::integ(double z)const{
  int i = find_bin(z);

  double total = cumint[i];
  double dx = x[i+1]-x[i];
  double dy = y[i+1] - y[i];
  double h = z - x[i];

  total += y[i]*h + 0.5*dy/dx * h*h;
  return total;
};

qspline::qspline(const vector& x, const vector& y)
  :n(x.size()),x(x),y(y),b(n-1),c(n-1),cumint(n) {
  vector dx(n-1), p(n-1);
  for(int i = 0; i<n-1; i++){dx[i] = x[i+1]-x[i]; p[i] = (y[i+1]-y[i])/dx[i];};
  for(int i = 0; i<n-2; i++) c[i+1] = 1/dx[i+1]*(p[i+1]-p[i] - c[i]*dx[i]);
  c[n-2]/=2;
  for(int i = n-3; i>=0; i--) c[i] = 1/dx[i] * (p[i+1] - p[i] - c[i+1]*dx[i+1]);
  for(int i = 0; i<n-1; i++) b[i] = p[i] - c[i]*dx[i];
  
  for(int i = 0; i<n-1; i++){
    cumint[i+1] = cumint[i] + y[i]*dx[i]+b[i]*dx[i]*dx[i]/2 + c[i]*dx[i]*dx[i]*dx[i]/3;
  };
};
int qspline::find_bin(double z)const{
  if(z>=x[last_index] && z<=x[last_index+1]) return last_index;
  last_index = binsearch(x, z);
  return last_index;
};
double qspline::eval(double z)const{
  int i = find_bin(z);
  double h = z - x[i];
  return y[i] + b[i]*h + c[i]*h*h;
};
double qspline::deriv(double z)const{
  int i = find_bin(z);
  return b[i] + 2*c[i]*(z-x[i]);
};
double qspline::integ(double z)const{
  int i = find_bin(z);

  double total = cumint[i];
  double h = z - x[i];
  total += y[i]*h + b[i]*h*h/2 + c[i]*h*h*h/3; 
  return total;
};
cspline::cspline(const vector& x, const vector& y)
  : n(x.size()), x(x), y(y), b(n), c(n-1), d(n-1), cumint(n) {
  vector h(n-1), p(n-1);
  for(int i = 0; i<n-1; i++){h[i] = x[i+1]-x[i]; p[i] = (y[i+1]-y[i])/h[i];}

  vector D(n), Q(n-1), B(n);
  D[0] = 2;
  Q[0] = 1;
  B[0] = 3*p[0];
  for(int i = 0; i<n-2; i++){
    D[i+1] = 2*h[i]/h[i+1]+2;
    Q[i+1] = h[i]/h[i+1];
    B[i+1] = 3*(p[i]+p[i+1]*h[i]/h[i+1]);
  }
  D[n-1] = 2;
  B[n-1] = 3*p[n-2];

  for(int i = 1; i<n; i++){
    D[i] -= Q[i-1]/D[i-1];
    B[i] -= B[i-1]/D[i-1];
  }

  b[n-1] = B[n-1]/D[n-1];
  for(int i = n-2; i>=0; i--){
    b[i] = (B[i]-Q[i]*b[i+1])/D[i];
  }
  for(int i = 0; i<n-1; i++){
    double dx = h[i];
    c[i] = (-2*b[i]-b[i+1]+3*p[i])/dx;
    d[i] = (b[i]+b[i+1]-2*p[i])/(dx*dx);
    cumint[i+1] = cumint[i] + dx*(y[i]+dx*(b[i]/2+dx*(c[i]/3 + dx*d[i]/4)));
  }
};
int cspline::find_bin(double z)const{
  if(z>=x[last_index] && z<=x[last_index+1]) return last_index;
  last_index = binsearch(x, z);
  return last_index;
};
double cspline::eval(double z)const{
  int i = find_bin(z);
  double h = z - x[i];
  return y[i] + b[i]*h + c[i]*h*h + d[i]*h*h*h;
};
double cspline::deriv(double z)const{
  int i = find_bin(z);
  double h = z - x[i];
  return b[i] + 2*c[i]*h + 3*d[i]*h*h;
};
double cspline::integ(double z)const{
  int i = find_bin(z);
  double h = z - x[i];
  double total = cumint[i];
  total += h*(y[i]+h*(b[i]/2+h*(c[i]/3 + h*d[i]/4)));
  return total;
};
}

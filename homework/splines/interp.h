#pragma once
#include <cassert>
#include "core/vector.h"

namespace pp{

struct lspline{
  const int n;
  vector x,y,cumint;
  mutable int last_index = 0;

  lspline(const vector& x, const vector& y);
  int find_bin(double z)const;
  double eval(double z)const;
  double deriv(double z)const;
  double integ(double z)const;
};
struct qspline{
  const int n;
  vector x,y,b,c,cumint;
  mutable int last_index = 0;

  qspline(const vector& x, const vector& y);
  int find_bin(double z)const;
  double eval(double z)const;
  double deriv(double z)const;
  double integ(double z)const;
};

struct cspline{
  const int n;
  vector x,y,b,c,d,cumint;
  mutable int last_index = 0;

  cspline(const vector& x, const vector& y);
  int find_bin(double z)const;
  double eval(double z)const;
  double deriv(double z)const;
  double integ(double z)const;
};

}

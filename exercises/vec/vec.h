#include <iostream>
#include <string>
#include <vector>
#pragma once
struct vec{
  double x,y,z;

  //constructors
  vec(double x, double y, double z) : x(x), y(y), z(z) {}
  vec() : vec(0, 0, 0) {}
  vec(const vec&) = default;
  vec(vec&&) = default;
  ~vec() = default;

  //assigment
  vec& operator=(const vec&) = default;
  vec& operator=(vec&&) = default;
  
  // arithmetic
  vec& operator+=(const vec&);
  vec& operator-=(const vec&);
  vec& operator*=(double);
  vec& operator/=(double);

  void print(const std::string& s="") const;
  friend std::ostream& operator<<(std::ostream&, const vec&);
};
// non-member
vec operator-(const vec&);
vec operator-(const vec&, const vec&);
vec operator+(const vec&, const vec&);
vec operator*(const vec&, double);
vec operator*(double, const vec&);
vec operator/(const vec&, double);

double dot(const vec&, const vec&);
double norm(const vec&);
vec cross(const vec&, const vec&);

std::vector<vec> GramSchmidt(const std::vector<vec>& input);

bool approx(const vec&, const vec&, double acc=1e-6, double eps=1e-6);

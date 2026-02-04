#include <iostream>
#include <cmath>
#include <string>
#include "vec.h"

// member operations
vec& vec::operator+=(const vec& other) {
  x += other.x;
  y += other.y;
  z += other.z;
  return *this;
}

vec& vec::operator-=(const vec& other) {
  x -=other.x;
  y -=other.y;
  z -=other.z;
  return *this;
}

vec& vec::operator*=(double s) {
  x *= s;
  y *= s;
  z *= s;
  return *this;
}

vec& vec::operator/=(double s) {
  x /= s;
  y /= s;
  z /= s;
  return *this;
}

void vec::print(const std::string& s) const {
  std::cout << s << x << " " << y << " " << z << std::endl;
}

std::ostream& operator<<(std::ostream& os, const vec& a) {
  os << "{ " << a.x << ", " << a.y << ", " << a.z << " } ";
  return os;
}

// non-member operations
vec operator-(const vec& v) {
  return vec(-v.x, -v.y, -v.z);
}

vec operator-(const vec& a, const vec& b) {
  vec res = a;
  res -= b;
  return res;
}

vec operator+(const vec& a, const vec& b) {
  vec res = a;
  res += b;
  return res;
}

vec operator*(const vec& a, double s) {
  vec res = a;
  res *= s;
  return res;
}

vec operator*(double s, const vec& a) {
  return a * s;
}

vec operator/(const vec& a, double s) {
  vec res = a;
  res /= s;
  return res;
}

double dot(const vec& a, const vec& b) {
  return a.x*b.x+a.y*b.y+a.z*b.z;
}

double norm(const vec& v) {
  return std::sqrt(dot(v, v));
}

vec cross(const vec& a, const vec& b) {
  return vec(a.y*b.z - a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

bool approx(double a, double b, double acc=1e-6, double eps=1e-6){
  if(std::fabs(a-b)<acc) return true;
  if(std::fabs(a-b)<eps*(std::fabs(a)+std::fabs(b))) return true;
  return false;
}

bool approx(const vec& a, const vec& b, double acc, double eps){
  if(!approx(a.x, b.x, acc, eps)) return false;
  if(!approx(a.y, b.y, acc, eps)) return false;
  if(!approx(a.z, b.z, acc, eps)) return false;
  return true;
}

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>

#pragma once
template<typename T>
struct vec{
  T x,y,z;

  //constructors
  vec(T x, T y, T z) : x(x), y(y), z(z) {};
  vec() : vec(0, 0, 0) {};
  vec(const vec&) = default;
  vec(vec&&) = default;
  ~vec() = default;

  //assigments 
  vec& operator=(const vec&) = default;
  vec& operator=(vec&&) = default;

  template <typename U>
    vec(const vec<U>& other) 
        : x(static_cast<T>(other.x)), 
          y(static_cast<T>(other.y)), 
          z(static_cast<T>(other.z)) {}

  template <typename U>
    vec& operator=(const vec<U>& other) {
        x = static_cast<T>(other.x);
        y = static_cast<T>(other.y);
        z = static_cast<T>(other.z);
        return *this;
    }

  // arithmetic
  template <typename U>
  vec& operator+=(const vec<U>& v){
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }
  template <typename U>
  vec& operator-=(const vec<U>& v){
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  template <typename S>
  vec& operator*=(S s){
    x*=s; y*=s; z*=s;
    return *this;
  }
  template <typename S>
  vec& operator/=(S s){
    x/=s; y/=s; z/=s;
    return *this;
  }

  void print(const std::string& s="") const {
    std::cout << s << x << " " << y << " " << z << std::endl;
  }
  friend std::ostream& operator<<(std::ostream& os, const vec& v){
    return os << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
  }
};
// non-member
template <typename T>
vec<T> operator-(const vec<T>& v){
  return vec(-v.x, -v.y, -v.z);
}

template <typename T, typename U>
auto operator-(vec<T> v, const vec<U>& u) -> vec<decltype(v.x-u.x)>{
  v -= u;
  return v;
}

template <typename T, typename U>
auto operator+(vec<T> v, const vec<U>& u) -> vec<decltype(v.x+u.x)>{
  v += u;
  return v;
}

template <typename T, typename S>
auto operator*(vec<T> v, S s){
  v *= s;
  return v;
}

template <typename T, typename S>
auto operator*(S s, const vec<T>& v){
  return v*s;
}

template <typename T, typename S>
auto operator/(vec<T> v, S s){
  v/=s;
  return v;
}

template <typename T>
auto smart_conj(const T& val){
  if constexpr (std::is_same_v<T, std::complex<float>> ||
                std::is_same_v<T, std::complex<double>> ||
                std::is_same_v<T, std::complex<long double>>) {
    return std::conj(val);
  }
  else return val;
}

template <typename T, typename U>
auto dot(const vec<T>& v, const vec<U>& u) -> decltype(v.x*u.x){
    return smart_conj(v.x)*u.x + smart_conj(v.y)*u.y + smart_conj(v.z)*u.z;
}

template <typename T>
auto norm(const vec<T>& v){
  return std::sqrt(std::abs(dot(v, v)));
}

template <typename T, typename U>
auto cross(const vec<T>& v, const vec<U>& u){
  return vec(v.y*u.z-v.z*u.y, v.z*u.x-v.x*u.z, v.x*u.y-v.y*u.x);
}

template <typename T>
std::vector<vec<T>> GramSchmidt(const std::vector<vec<T>>& input){
  std::vector<vec<T>> basis;
  for(const auto& v : input){
    vec<T> w = v;
    for(const auto u: basis){
      auto proj_scale = dot(u, w)/dot(u,u);
      w -= u*proj_scale;
    }
    auto w_norm = norm(w);
    if(w_norm > 1e-9) basis.push_back(w/w_norm);
  }
  return basis;
}

template <typename T, typename U>
bool approx(const vec<T>& v, const vec<U>& u, double acc=1e-6, double eps=1e-6){
  if(std::abs(v.x - u.x) < acc 
    && std::abs(v.y - u.y) < acc
    && std::abs(v.z - u.z) < acc){
    return true;
  }
  else if(std::abs(v.x-u.x)/std::max(std::abs(v.x), std::abs(u.x)) <= eps 
      && std::abs(v.y-u.y)/std::max(std::abs(v.y), std::abs(u.y)) <= eps 
      && std::abs(v.z-u.z)/std::max(std::abs(v.z), std::abs(u.z)) <= eps){
    return true;
  }
  else return false;
}

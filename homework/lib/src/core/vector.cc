#include <iostream>
#include <iomanip>
#include <random>
#include "core/vector.h"

namespace pp{
vector::vector(int n) : data(n, 0.0), size_(n), stride(1){
  start = data.data();
}
vector::vector(std::initializer_list<double> list) : data(list), size_(list.size()), stride(1){
  start = data.data();
}
vector::vector(double* ptr, int n, int s) : start(ptr), size_(n), stride(s) {};

vector::vector(const pp::vector& other)
  : data(other.size_), size_(other.size_), stride(1) 
  {
    double* dst = data.data();
    for(int i = 0; i<size_; i++){
      dst[i] = other[i];
    }
    this->start = dst;
  };

vector& vector::operator+=(const vector& other){
  if(size_ != other.size_){
    throw std::runtime_error("Vector size mismatch in +=");
  }
  for(int i = 0; i<size(); i++){
    (*this)[i] += other[i];
  }
  return *this;
}
vector& vector::operator-=(const vector& other){
  if(size_ != other.size_){
    throw std::runtime_error("Vector size mismatch in -=");
  }
  for(int i = 0; i<size_; i++){
    (*this)[i] -= other[i];
  }
  return *this;
}
vector& vector::operator*=(const double s){
  for(int i = 0; i<size_; i++){
    (*this)[i] *= s;
  }
  return *this;
}
vector& vector::operator/=(const double s){
  for(int i = 0; i<size_; i++){
    (*this)[i] /= s;
  }
  return *this;
}
double vector::norm()const{
  double sum = 0;
  double* ptr = start;
  for(int i = 0; i<size_; ++i){
    double val = *ptr;
    sum += val*val;
    ptr += stride;
  }
  return std::sqrt(sum);
}
void vector::copy_from(const vector& other){
  if(this->size_ != other.size_){
    throw std::runtime_error("Vector size mismatch in copy_from");
  }
  double* dst = this->start;
  double* src = other.start;
  int s_dst = this->stride;
  int s_src = other.stride;
  for(int i = 0; i<size_; ++i){
    *dst = *src;
    dst += s_dst;
    src += s_src;
  }
}
void vector::point_to(const vector& other){
  this->data.clear();
  this->start = other.start;
  this->size_ = other.size_;
  this->stride = other.stride;
}
vector& vector::operator=(const vector& other){
  if(this == &other) return *this;

  if(this->size_ == other.size_ && this->stride == 1){
    this->copy_from(other);
  }
  else{
    this->size_ = other.size_;
    this->stride = 1;
    this->data = other.data;

    if(other.data.empty() && other.size_ > 0){
      this->data.assign(other.size_, 0.0);
      this->start = this->data.data();
      this->copy_from(other);
    }
    else this->start = this->data.data();
  }
  return *this;
}
void vector::fill_random(double min, double max){
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);
  for(int i = 0; i<size(); i++){
    (*this)[i] = dist(gen);
  }
}
std::ostream& operator<<(std::ostream& os, const vector& v){
  os << "[";
  for(int i = 0; i<v.size(); ++i){
    os << std::setw(10) << v[i];
    if(i<v.size()-1) os << ", ";
  }
  os << "]";
  return os;
}
vector operator+(vector v, const vector& u){
  v+=u;
  return v;
}
vector operator-(vector v, const vector& u){
  v-=u;
  return v;
}
vector operator-(vector v){
  v*=(-1);
  return v;
}
vector operator*(vector v, const double s){
  v*=s;
  return v;
}
vector operator*(const double s, vector v){
  v*=s;
  return v;
}
vector operator/(vector v, const double s){
  v/=s;
  return v;
}
double dot(const vector& v, const vector& u){
  int n = v.size();
  if(n != u.size()){
    throw std::runtime_error("Vector size mismatch in dot");
  }
  double sum = 0;
  double* ptr_v = v.start;
  double* ptr_u = u.start;
  int stride_v = v.stride;
  int stride_u = u.stride;
  for(int i = 0; i<v.size(); ++i){
    sum += (*ptr_v)*(*ptr_u);
    ptr_v += stride_v;
    ptr_u += stride_u;
  }
  return sum;
}
}

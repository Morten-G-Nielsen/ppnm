#pragma once
#include <iostream>
#include <vector>
#include <initializer_list>

namespace pp{
  struct vector {
    std::vector<double> data;

    double* start = nullptr;
    int size_ = 0;
    int stride = 1;

    vector(int n);
    vector(std::initializer_list<double> list);
    vector(double* ptr, int n, int s=1);
    vector(const vector& other);

    inline double& operator[](int i){
      if(stride == 1){
        return start[i];
      }
      return start[i*stride];
    }
    inline double operator[](int i) const {
      if(stride == 1){
        return start[i];
      }
      return start[i*stride];
    }
    inline void sub_scaled(const vector& other, double s){
      double* dst = this->start;
      double* src = other.start;
      int dst_stride = this->stride;
      int src_stride = other.stride;
      for(int k=0;k<size_; ++k){
        *dst -= (*src) * s;
        dst += dst_stride;
        src += src_stride;
      }
    }
    int size() const {return size_;} 

    vector& operator+=(const vector& other);
    vector& operator-=(const vector& other);
    vector& operator*=(const double s);
    vector& operator/=(const double s);

    vector& operator=(const vector& other); // Deep copy assigement
    void copy_from(const vector& other); // Data move
    void point_to(const vector& other); // Pointer rerouting used for the "View"

    double norm() const;
    void fill_random(double min=-1, double max=1);
    friend std::ostream& operator<<(std::ostream& os, const vector& v);
  };

  vector operator+(const vector& v, const vector& u);
  vector operator-(vector v, const vector& u);
  vector operator-(vector v);
  vector operator*(vector v, const double s);
  vector operator*(const double s, vector v);
  vector operator/(vector v, const double s);
  double dot(const vector& v, const vector& u);
  bool approx(double x, double y, double acc=1e-6, double eps=1e-6);
  bool approx(const vector& v, const vector& u, double acc=1e-6, double eps=1e-6);
}

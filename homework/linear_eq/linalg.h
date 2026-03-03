#pragma once
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <initializer_list>
#include <random>

namespace pp {
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

    void copy_from(const vector& other); // Deep copy of values
    vector& operator=(const vector& other); // Shallow copy of pointers/view

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

  struct matrix {
    std::vector<double> data;
    int rows_;
    int cols_;
    bool transposed = false;
    matrix(int r, int c);
    matrix(const matrix& other);

    int row_count() const {return rows_;}
    int col_count() const {return cols_;}

    double& operator()(int r, int c){
      if(transposed) return data[c*rows_+r];
      return data[r*cols_ + c];
    }
    const double& operator()(int r, int c)const{
      if(transposed) return data[c*cols_+r];
      return data[r*cols_ + c];
    }
    inline void lazy_transpose(){
      std::swap(rows_, cols_);
      transposed = !transposed;
    }

    vector row(int r);
    vector row(int r)const;
    vector column(int c);
    vector column(int c)const;

    matrix transpose() const;
    static matrix identity(int n);
    void fill_random(double min=-1, double max = 1);

    matrix& operator+=(const matrix& other);
    matrix& operator-=(const matrix& other);
    matrix& operator*=(const double s);
    matrix& operator/=(const double s);
    matrix& operator*=(const matrix& other);

    friend std::ostream& operator<<(std::ostream& os, const matrix& A);
  };

  matrix operator+(matrix A, const matrix& B);
  matrix operator-(matrix A, const matrix& B);
  matrix operator*(matrix A, const double s);
  matrix operator*(const double s, matrix A);
  matrix operator*(matrix A, const matrix& B);
  vector operator*(const matrix& A, const vector& v);
  matrix operator/(matrix A, const double s);

struct QR{
  matrix Q;
  matrix R;
  bool is_valid = false;
  QR(const matrix& other);

  vector solve(const vector& b)const;
  double det()const;
  matrix inverse()const;
};
}

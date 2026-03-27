#pragma once
#include <vector>
#include "core/vector.h"

namespace pp{
  struct matrix {
    std::vector<double> data;
    int rows_;
    int cols_;
    bool transposed = false;
    matrix(int r, int c);
    matrix(const matrix& other);

    int row_count() const {return transposed ? cols_ : rows_;}
    int col_count() const {return transposed ? rows_ : cols_;}

    double& operator()(int r, int c){
      if(transposed) return data[c*cols_+r];
      return data[r*cols_ + c];
    }
    const double& operator()(int r, int c)const{
      if(transposed) return data[c*cols_+r];
      return data[r*cols_ + c];
    }
    inline void lazy_transpose(){
      this->transposed = !this->transposed;
    }

    vector row(int r);
    vector row(int r)const;
    vector column(int c);
    vector column(int c)const;

    matrix transpose() const;
    static matrix identity(int n);
    void fill_random(double min=-1, double max = 1);
    bool is_symetric(double tol = 1e-12)const;

    matrix& operator=(const matrix& other);
    matrix& operator=(matrix&& other)noexcept;
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
}

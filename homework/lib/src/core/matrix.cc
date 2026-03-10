#include <iostream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include "../../include/core/vector.h"
#include "../../include/core/matrix.h"

namespace pp{ 
matrix::matrix(int r, int c): data(r*c), rows_(r), cols_(c) {};
matrix::matrix(const matrix& other)
  : data(other.data), rows_(other.rows_), cols_(other.cols_), transposed(other.transposed)
  {};
vector matrix::row(int r){
  if(transposed){
    return vector(&data[r], rows_, cols_);
  }
  return vector(&data[r * cols_], cols_, 1);
}
vector matrix::row(int r)const{
  if(transposed){
    return vector(const_cast<double*>(&data[r]), rows_, cols_);
  }
  return vector(const_cast<double*>(&data[r*cols_]),cols_,1);
}
vector matrix::column(int c){
  if(transposed){
    return vector(&data[c*cols_], cols_, 1);
  }
  return vector(&data[c], rows_, cols_);
}
vector matrix::column(int c)const{
  if(transposed){
    return vector(const_cast<double*>(&data[c*cols_]), cols_, 1);
  }
  return vector(const_cast<double*>(&data[c]), rows_, cols_);
}
matrix matrix::transpose()const{
  matrix R(this->col_count(), this->row_count());
  for(int i = 0; i<this->row_count(); i++){
    for(int j = 0; j<this->col_count(); j++){
      R(j,i) = (*this)(i, j);
    }
  }
  return R;
}
matrix matrix::identity(int n){
 matrix I(n,n);
 for(int i = 0; i<n; i++){
   I(i,i) = 1;
 }
 return I;
}
void matrix::fill_random(double min, double max){
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<double> dist(min, max);
  for(int i = 0; i< row_count(); i++){
    for(int j = 0; j<col_count(); j++){
      (*this)(i,j) = dist(gen);
    }
  }
}
bool matrix::is_symetric(double tol)const{
  if(rows_ != cols_) return false;
  for(int i = 0; i<rows_; ++i){
    for(int j = i+1; j<cols_; ++j){
      if(std::abs((*this)(i,j)-(*this)(j,i)) > tol){
        return false;
      }
    }
  }
  return true;
}
matrix& matrix::operator=(const matrix& other){
  if(this == &other) return *this;
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->transposed = other.transposed;
  this->data = other.data;
  return *this;
}
matrix& matrix::operator=(matrix&& other)noexcept{
  if(this == &other) return *this;
  this->data = std::move(other.data);
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->transposed = other.transposed;
  return *this;
}
matrix& matrix::operator+=(const matrix& other){
  for(int i = 0; i<rows_; i++) (*this).row(i) += other.row(i);
  return *this;
}
matrix& matrix::operator-=(const matrix& other){
  for(int i = 0; i<rows_; i++) (*this).row(i) -= other.row(i);
  return *this;
}
matrix& matrix::operator*=(const double s){
  for(int i = 0; i<rows_; i++) (*this).row(i) *= s;
  return *this;
}
matrix& matrix::operator/=(const double s){
  for(int i = 0; i<rows_; i++) (*this).row(i) /= s;
  return *this;
}
matrix& matrix::operator*=(const matrix& other){
  if(this->col_count() != other.row_count()) throw std::runtime_error("Dimension mismatch");

  matrix res(this->row_count(), other.col_count());
  matrix B_T = other.transpose();

  vector row_A(0);

  for(int i = 0; i<rows_; ++i){
      row_A.point_to(this->row(i));
      for(int j = 0; j<other.col_count(); ++j){
        res(i,j) += dot(row_A, B_T.row(j));
      }
    }
  *this = std::move(res);
  return *this;
}
std::ostream& operator<<(std::ostream& os, const matrix& A){
  for(int i = 0; i<A.row_count(); i++){
    for(int j = 0; j<A.col_count(); j++){
      os << std::setw(12) << std::fixed << std::setprecision(4) << A(i,j);
    }
    os << "\n";
  }
  return os;
}

matrix operator+(matrix A, const matrix& B){
  if(A.row_count() != B.row_count()){
    throw std::runtime_error("Matrix row mismatch in +");
  }
  if(A.col_count() != B.col_count()){
    throw std::runtime_error("Matrix column mismatch in +");
  }
  for(int i = 0; i<A.row_count(); i++) A.row(i) += B.row(i);
  return A;
}
matrix operator-(matrix A, const matrix& B){
  if(A.row_count() != B.row_count()){
    throw std::runtime_error("Matrix row mismatch in -");
  }
  if(A.col_count() != B.col_count()){
    throw std::runtime_error("Matrix column mismatch in -");
  }
  for(int i = 0; i<A.row_count(); i++) A.row(i) -= B.row(i);
  return A;
}
matrix operator/(matrix A, const double s){
  for(int i = 0; i<A.row_count(); i++) A.row(i)/=s;
  return A;
}
matrix operator*(matrix A, const double s){
  for(int i = 0; i<A.row_count(); i++) A.row(i)*=s;
  return A;
}
matrix operator*(const double s, matrix A){
  for(int i = 0; i<A.row_count(); i++) A.row(i)*=s;
  return A;
}
vector operator*(const matrix& A, const vector& v){
  vector res(A.row_count());
  for(int i = 0; i<A.row_count(); i++){
    res[i] = dot(A.row(i), v);
  }
  return res;
}
matrix operator*(matrix A, const matrix& B){
  A*=B;
  return A;
}
}

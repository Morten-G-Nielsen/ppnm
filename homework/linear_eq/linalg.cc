#include <iostream>
#include <iomanip>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <random>
#include "linalg.h"

namespace pp{
// Vector implementation
vector::vector(int n) : data(n, 0.0), size_(n), stride(1){
  start = data.data();
}
vector::vector(std::initializer_list<double> list) : data(list), size_(list.size()), stride(1){
  start = data.data();
}
vector::vector(double* ptr, int n, int s) : start(ptr), size_(n), stride(s) {};

vector::vector(const pp::vector& other)
  : size_(other.size_), stride(1)
  {
    this->data.reserve(other.size_);
    for(int i = 0; i<other.size_; i++){
      this->data.push_back(other[i]);
    }
    this->start = this->data.data();
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
  if(size_ != other.size_){
    throw std::runtime_error("Vector size mismatch in copy_from");
  }
  for(int i = 0; i<size_; i++){
    (*this)[i] = other[i];
  }
}
vector& vector::operator=(const vector& other){
  this->start = other.start;
  this->size_ = other.size_;
  this->stride = other.stride;
  this->data = other.data;
  if(!this->data.empty()){
    this->start = this->data.data();
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

// Matrix implementation
matrix::matrix(int r, int c): data(r*c), rows_(r), cols_(c) {};
matrix::matrix(const matrix& other)
  : data(other.data), rows_(other.rows_), cols_(other.cols_) 
{
  this->data = other.data;
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->transposed = other.transposed;
};
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
    return vector(&data[c*rows_], rows_, 1);
  }
  return vector(&data[c], rows_, cols_);
}
vector matrix::column(int c)const{
  if(transposed){
    return vector(const_cast<double*>(&data[c*rows_]), rows_, 1);
  }
  return vector(const_cast<double*>(&data[c]), rows_, cols_);
}
matrix matrix::transpose()const{
  matrix R(cols_, rows_);
  for(int i = 0; i<rows_; i++){
    for(int j = 0; j<cols_; j++){
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
  if(this->cols_ != other.rows_) throw std::runtime_error("Dimension mismatch");

  matrix res(this->rows_, other.cols_);
  matrix B_T = other.transpose();

  for(int i = 0; i<rows_; ++i){
      vector row_A = this->row(i);
      for(int j = 0; j<other.cols_; ++j){
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

// QR implementation
QR::QR(const matrix& other)
  : Q(other.transpose()), R(other.col_count(), other.col_count()){
    int n = other.col_count();
    for(int i = 0; i<n; i++){
      vector qi = Q.row(i);
      double norm = qi.norm();
      this->R(i,i) = norm;
      if(this->R(i,i) < 1e-12) throw std::runtime_error("Linearly dependent columns");
      qi/=norm;
      for(int j = i+1; j<n; j++){
        vector qj = Q.row(j);
        double dot_val = dot(qi, qj);
        this->R(i,j) = dot_val;
        qj.sub_scaled(qi, dot_val);
      }
    }
    is_valid = true;
    this->Q.lazy_transpose();
  }
vector QR::solve(const vector& b)const{
  if(!is_valid) throw std::runtime_error("QR not initilized");
  int n = R.col_count();

  vector x(n);
  for(int i = 0; i<n; i++){
    double sum = dot(Q.column(i), b);
    x[i] = sum;
  }

  for(int i = n-1; i>=0; i--){
    double sum = 0;
    const double* r_ptr = &R(i,0);
    double* x_ptr = x.start;

    for(int j = i+1; j<n; j++){
      sum += r_ptr[j]*x_ptr[j];
    }
    x_ptr[i] = (x_ptr[i] - sum)/r_ptr[i];
  }
  return x;
}
double QR::det()const{
  double prod = 1;
  for(int i = 0; i<this->R.row_count(); i++){
    prod *= this->R(i,i);
  }
  return prod;
}
matrix QR::inverse()const{
  if(this->Q.row_count() != this->Q.col_count())
    throw std::runtime_error("Inverse needs square matrix");
  int n = this->R.col_count();
  matrix inv(n,n);

  vector x_buffer(n);

  for(int k = 0; k<n; ++k){
    x_buffer.copy_from(Q.row(k));

    double* x_ptr = x_buffer.start;
    for(int i = n-1; i>=0; --i){
      double sum = 0;
      const double* r_ptr = &R(i, 0);

      for(int j = i+1; j<n; ++j){
        sum += r_ptr[j] * x_ptr[j];
      }
      x_ptr[i] = (x_ptr[i]-sum)/r_ptr[i];
    }
    inv.column(k).copy_from(x_buffer);
  }
  return inv;
}
} 

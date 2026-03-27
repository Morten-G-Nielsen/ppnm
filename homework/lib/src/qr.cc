#include "core/vector.h"
#include "core/matrix.h"
#include "qr.h"

namespace pp{
QR::QR(const matrix& other)
  : Q(other.transpose()), R(other.col_count(), other.col_count()){
    int n = other.col_count();
    vector qi(0);
    vector qj(0);
    for(int i = 0; i<n; i++){
      qi.point_to(Q.row(i));
      double norm = qi.norm();
      this->R(i,i) = norm;
      if(this->R(i,i) < 1e-12) throw std::runtime_error("Linearly dependent columns");
      qi/=norm;
      for(int j = i+1; j<n; j++){
        qj.point_to(Q.row(j));
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

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

    vector(int n) : data(n), start(data.data()), size_(n), stride(1) {};
    vector(std::initializer_list<double> list)
      : data(list), start(data.data()), size_(list.size()), stride(1) {};

    vector(double* ptr, int n, int s=1) : start(ptr), size_(n), stride(s) {}
    vector(const vector& other)
      : data(other.data),
        size_(other.size_),
        stride(other.stride)
      {
        if(other.data.empty()) start = other.start;
        else start = data.data();
      }
    vector& operator=(const vector& other){
      if(size_ != other.size_) throw std::invalid_argument("size mismatch");
      for(int i = 0; i < size_; i++){
        (*this)[i] = other[i];
      }
      return *this;
    }
    
    int size() const {return size_;} 
    double& operator[](int i) {return start[i*stride];}
    const double& operator[](int i) const {return start[i*stride];}

    vector& operator+=(const vector& other){
      for(int i = 0; i<size(); i++)(*this)[i] += other[i];
      return *this;
    }
    vector& operator-=(const vector& other){
      for(int i = 0; i<size(); i++)(*this)[i] -= other[i];
      return *this;
    }
    vector& operator*=(const double& s){
      for(int i = 0; i<size(); i++)(*this)[i] *= s;
      return *this;
    }
    vector& operator/=(double s){
      for(int i = 0; i<size(); i++)(*this)[i] /= s;
      return *this;
    }
    double norm() const {
      double res = 0;
      for(int i = 0; i<size(); i++) res += (*this)[i] * (*this)[i];
      return std::sqrt(res);
    }
    void copy_from(const vector& other){
      for(int i = 0; i<size_; i++)(*this)[i] = other[i];
    }
    void fill_random(double min=-1, double max = 1){
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dist(min, max);
      for(int i = 0; i<size(); i++){
        (*this)[i] = dist(gen);
      }
    }
    friend std::ostream& operator<<(std::ostream& os, const vector& v){
      os << "{ ";
      for(int n=0; n<v.size(); n++){
        os << v[n] << " ";
      }
      return os << "}\n";
    }
  };

  vector operator+(const vector& v, const vector& u){
    vector res(v.size());
    for(int i=0;i<v.size(); i++) res[i] = v[i]+u[i];
    return res;
  }
  vector operator-(const vector& v, vector& u){
    vector res(v.size());
    for(int i = 0; i<v.size(); i++) res[i] = v[i]-u[i];
    return res;
  }
  vector operator-(const vector& v){
    vector res(v.size());
    for(int i=0;i<v.size(); i++) res[i] = -v[i];
    return res;
  }
  vector operator*(const vector& v, const double s){
    vector res(v.size());
    for(int i=0; i<v.size(); i++) res[i] = v[i]*s;
    return res;
  }
  vector operator*(const double s, const vector& v){
    vector res(v.size());
    res = v*s;
    return res;
  }
  vector operator/(const vector& v, const double s){
    vector res(v.size());
    res = v * (1/s);
    return res;
  }
  double dot(const vector& v, const vector& u){
    double res = 0;
    for(int i = 0; i<v.size(); i++) res+=v[i]*u[i];
    return res;
  }
  bool approx(double x, double y, double acc=1e-6, double eps=1e-6){
    if(std::abs(x-y) < acc) return true;
    if(std::abs(x-y) < eps*std::max(std::abs(x), std::abs(y))) return true;
    return false;
  }
  bool approx(const vector& v, const vector& u, double acc=1e-6, double eps=1e-6){
    if(v.size() != u.size()) return false;
    for(int i = 0; i<v.size(); i++){
      if(!approx(v[i], u[i])) return false;
    }
    return true;
  }

  struct matrix {
    std::vector<double> data;
    int rows_;
    int cols_;
    matrix(int r, int c) : data(r*c), rows_(r), cols_(c) {}

    int row_count() const {return rows_;}
    int col_count() const {return cols_;}

    double& operator()(int r, int c){return data[r*cols_ + c];}
    const double& operator()(int r, int c)const{return data[r*cols_ + c];}

    vector row(int r){return vector(&data[r*cols_], cols_, 1);}
    vector row(int r)const{
      return vector(const_cast<double*>(&data[r*cols_]), cols_, 1);
    }
    vector column(int c){return vector(&data[c], rows_, cols_);}
    vector column(int c)const{
      return vector(const_cast<double*>(&data[c]), rows_, cols_);
    }

    matrix transpose() const{
      matrix R(col_count(), row_count());
      for(int row = 0; row<row_count(); row++){
        for(int col = 0; col<col_count(); col++){
          R(col,row)=(*this)(row,col);
        }
      }
    return R;
    }
    static matrix identity(int n){
      matrix I(n,n);
      for(int i = 0; i<n; i++){
        I(i,i) = 1;
      }
      return I;
    }
    void fill_random(double min=-1, double max = 1){
      std::random_device rd;
      std::mt19937 gen(rd());

      std::uniform_real_distribution<double> dist(min, max);
      for(int i = 0; i< row_count(); i++){
        for(int j = 0; j<col_count(); j++){
          (*this)(i,j) = dist(gen);
        }
      }
    }
    matrix& operator+=(const matrix& other){
      for(int i = 0; i<row_count(); i++)(*this).row(i)+=other.row(i);
      return *this;
    }
    matrix& operator-=(const matrix& other){
      for(int i = 0; i<row_count(); i++)(*this).row(i)-=other.row(i);
      return *this;
    }
    matrix& operator*=(const double s){
      for(int i = 0; i<row_count(); i++)(*this).row(i)*=s;
      return *this;
    }
    matrix& operator/=(const double s){
      for(int i = 0; i<row_count(); i++)(*this).row(i)/=s;
      return *this;
    }
    matrix& operator*=(const matrix& other){
      if(cols_ != other.rows_) throw std::invalid_argument("size mismatch");
      matrix res(rows_, other.cols_);
      
      for(int i = 0; i<rows_; i++){
        for(int j = 0; j<other.cols_; j++){
          double sum = 0.0;
          for(int k = 0; k<cols_; k++){
            sum += (*this)(i,k) * other(k,j);
          }
          res(i,j) = sum;
        }
      }
    *this = res;
    return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const matrix& A){
      for(int i=0; i<A.row_count(); i++){
        os << A.row(i);
      }
      return os;
    }
  };

  matrix operator+(matrix A, const matrix& B){
    for(int i = 0; i<A.row_count(); i++) A.row(i)+=B.row(i);
    return A;
  }
  matrix operator-(matrix A, const matrix& B){
    for(int i = 0; i<A.row_count(); i++) A.row(i)-=B.row(i);
    return A;
  }
  matrix operator*(matrix A, const double s){
    A*=s;
    return A;
  }
  matrix operator*(const double s, matrix A){
    A*=s;
    return A;
  }
  matrix operator*(matrix A, const matrix& B){
    A*=B;
    return A;
  }
  vector operator*(const matrix& A, const vector& v){
    vector u(v.size());
    for(int i = 0; i<A.row_count(); i++){
      u[i] = dot(A.row(i),v);
    }
    return u;
  }
  matrix operator/(matrix A, const double s){
    A/=s;
    return A;
  }
}

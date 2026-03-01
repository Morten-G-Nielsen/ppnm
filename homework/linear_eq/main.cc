#include <iostream>
#include <string>
#include "linalg.h"

struct QR{
  pp::matrix Q;
  pp::matrix R;
  bool is_valid = false;
  QR(const pp::matrix A) : Q(A.row_count(), A.col_count()), R(A.col_count(), A.col_count()){
    int n = A.col_count();
    for(int i = 0; i<n; i++){
      auto ai = A.column(i);
      this->R(i,i) = ai.norm();
      if(this->R(i,i) < 1e-12) throw std::runtime_error("Linearly dependent columns");

      this->Q.column(i).copy_from(ai);
      this->Q.column(i)/=this->R(i,i);
      for(int j = i+1; j<n; j++){
        auto aj = A.column(j);
        this->R(i,j) = pp::dot(this->Q.column(i), aj);
        aj -= this->Q.column(i)*this->R(i,j);
      }
    }
    is_valid = true;
  }
   
  pp::vector solve(const pp::vector& b) const{
    if(!is_valid) throw std::runtime_error("QR not initilized");
    int n = R.col_count();
    int m = Q.row_count();

    pp::vector d(n);
    for(int i = 0; i < n; i++){
      double sum = 0;
      for(int j = 0; j < m; j++){
        sum += Q(j,i)*b[j];
      }
      d[i] = sum;
    }
    pp::vector x(n);
    for(int i = n-1; i>=0; i--){
      double sum = 0;
      for(int j = i+1; j<n; j++){
        sum += R(i,j)*x[j];
      }
      x[i] = (d[i] - sum)/R(i,i);
    }
    return x;
  }
  double det()const{
      double prod = 1;
      for(int i = 0; i < R.row_count(); i++){
        prod *= R(i,i);
      }
      return prod;
  }
  pp::matrix inverse()const{
    if(Q.row_count() != Q.col_count())
      throw std::runtime_error("Inverse needs square matrix");

    int n = Q.row_count();
    pp::matrix I = pp::matrix::identity(n);
    pp::matrix A_inv(n, n);

    for(int j = 0; j<n; j++){
      pp::vector col_ = this->solve(I.column(j));
      A_inv.column(j).copy_from(col_);
    }
    return A_inv;
  }
};

int main(int argc, char** argv){
  if(argc<2){
    int n1 = 5;
    int m = 3;
    pp::matrix A1(n1,m); A1.fill_random();
    QR sys1(A1);
    auto Q = sys1.Q;
    auto R = sys1.R;
    std::cout << "QR-decomp test\n"<< "A1\n"<< A1;
    std::cout << "R\n"<< R;
    std::cout << "Q^T*Q\n"<< Q.transpose()*Q<<"\n";
    std::cout << "Q*R-A\n"<< Q*R - A1 << "\n";

    int n2 = 4;
    pp::matrix A2(n2, n2); A2.fill_random();
    pp::vector b(n2); b.fill_random();
    QR sys2(A2);
    auto x = sys2.solve(b);
    std::cout << "Solve test\n"<< "A2\n"<< A2 << "b\n" << b;
    std::cout << "A2*x - b\n" << A2*x-b << "\n";
    auto B = sys2.inverse();
    std::cout << "Inverse test\n"<<"A2*B\n" << A2*B << "\n";
  }
  else {
    int N = 1;
    for(int i = 0; i<argc; i++){
      std::string arg = argv[i];
      if(arg == "-size" && i+1<argc) N = std::stoi(argv[++i]);
    }
    pp::matrix A(N, N); A.fill_random();
    QR sys(A);
  }
  return 0;
}

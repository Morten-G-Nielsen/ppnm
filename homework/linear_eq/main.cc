#include <iostream>
#include <string>
#include "linalg.h"

int main(int argc, char** argv){
  if(argc<2){
    int n1 = 2;
    int m = 2;
    pp::matrix A1(n1,m);
    A1.row(0).copy_from(pp::vector{3,1});
    A1.row(1).copy_from(pp::vector{4,2});
    std::cout << A1 << "\n";
    std::cout << A1.column(0) <<"\n";
    pp::QR sys1(A1);
    auto Q = sys1.Q;
    auto R = sys1.R;
    std::cout << "QR-decomp test\n"<< "A1\n"<< A1;
    std::cout << "R\n"<< R;
    std::cout << "Q^T*Q\n"<< Q.transpose()*Q<<"\n";
    std::cout << "Q*R-A\n"<< Q*R - A1 << "\n";

    int n2 = 4;
    pp::matrix A2(n2, n2); A2.fill_random();
    pp::vector b(n2); b.fill_random();
    pp::QR sys2(A2);
    auto x = sys2.solve(b);
    std::cout << "Solve test\n"<< "A2\n"<< A2 << "b\n" << b;
    std::cout << "\nA2*x - b\n" << A2*x-b << "\n";
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
  }
  return 0;
}

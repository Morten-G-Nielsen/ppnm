#include <iostream>
#include "core/matrix.h"
#include "evd.h"

int main(){
  int n = 3;
  pp::matrix A(n,n);
  A.fill_random();
  for(int i = 0; i<n; i++){
    for(int j = i+1; j<n; j++){
      A(i,j)=A(j,i);
    }
  }

  pp::EVD sys(A, 1e-12);
  auto V = sys.eigenvectores;
  auto D = V.transpose()*A*V;
  std::cout << "A\n"<<A<<"\n";
  std::cout << "V^TAV = D\n"<<D<<"\n";
  std::cout << "V\n"<<V<<"\n";
  std::cout << "VDV^T\n"<<V*D*V.transpose()<<"\n";
  std::cout << "V^TV\n"<<V.transpose()*V <<"\n";
  std::cout << "VV^T\n"<<V*V.transpose() << "\n";
  double sum1 = 0;
  double sum2 = 0;
  for(int i = 0; i<A.row_count(); i++){
    sum1 += A(i,i);
    sum2 += D(i,i);
  }
  std::cout << "Tr A = " << sum1<<"\n";
  std::cout << "Tr D = " << sum2<<"\n";
  return 0;
}

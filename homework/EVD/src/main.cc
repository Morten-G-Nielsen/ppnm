#include <iostream>
#include <cmath>
#include <string>
#include "vector.h"
#include "matrix.h"
#include "evd.h"

int main(int argc, char* argv[]){
  int N = 1;
  double rmax = 10;
  double dr = 0.1;
  double eps = 1e-12;
  bool wf = false;
  for(int i = 0; i<argc; ++i){
    std::string arg = argv[i];
    if(arg == "-rmax" && i+1<argc) rmax = std::stod(argv[++i]);
    if(arg == "-dr" && i+1<argc) dr = std::stod(argv[++i]);
    if(arg == "-wf" && i+1<argc) wf = true;
    if(arg == "-eps" && i+1<argc) eps = std::stod(argv[++i]);
    if(arg == "-n" && i+1<argc) N = std::stoi(argv[++i]);
  }
  // For timing
  if(N!=1){
    pp::matrix A(N,N);
    A.fill_random();
    for(int i = 0; i<N; i++){
      for(int j = i+1; j<N; j++){
        A(i,j)=A(j,i);
      }
    }
    pp::EVD sys(A, eps);
  }
  else {
    int npoints = int(rmax/dr)-1;
    pp::vector r(npoints);
    for(int i = 0; i<npoints; ++i)r[i] = dr*(i+1);
    pp::matrix H(npoints, npoints);
    for(int i = 0; i<npoints-1; ++i){
      H(i,i) = -2*(-0.5/dr/dr);
      H(i,i+1) = H(i+1,i)= 1*(-0.5/dr/dr);
    }
    H(npoints-1, npoints-1)=-2*(-0.5/dr/dr);
    for(int i=0; i<npoints;++i)H(i,i) += -1/r[i];

    pp::EVD sys(H, eps);
    auto d = sys.eigenvalues;
    auto V = sys.eigenvectores;
    // For radial wavefunction
    if(wf){
      auto V = sys.eigenvectores;
      for(int i = 0; i<npoints; i++){
        double r = (i+1)*dr;
        std::cout << r << " " << V(i,0)/std::sqrt(dr) << " " << V(i,1)/std::sqrt(dr) << " " << V(i,2)/std::sqrt(dr) << "\n";
      }
    }
    // For ground state energy
    else {
      std::cout << d[0] << "\n";
    }
  }
  return 0;
}

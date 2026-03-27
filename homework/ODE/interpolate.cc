#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "core/vector.h"
#include "interp.h"

int main(int argc, char* argv[]){
  std::string infile, outfile;
  for(int i = 1; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "--input" && i+1<argc) infile = argv[++i];
    else if (arg == "--output" && i+1<argc) outfile = argv[++i];
  }
  std::ifstream input(infile);
  std::ofstream output(outfile);
  double x;
  std::vector<double> phi, u;
  if(input.is_open() && output.is_open()){
    int t = 0;
     while(input >> x){
       if(t%3 == 0) phi.push_back(x); 
       if((t+2)%3 == 0) u.push_back(x);
       t += 1;
     } 
  }
  int n = u.size();
  pp::vector xvec(n), yvec(n);
  for(int i = 0; i<xvec.size(); i++){
    xvec[i] = phi.data()[i];
    yvec[i] = u.data()[i];
  }
  int steps = 1000;
  double xmax = xvec[n-1];
  double dx = (xmax-1)/(steps-1);
  pp::cspline spline(xvec, yvec);
  for(int i = 0; i<steps; i++){
    double x = dx*i;
    double y = spline.eval(x);
    output << x << " " << y << "\n";
  }
  return 0;
}

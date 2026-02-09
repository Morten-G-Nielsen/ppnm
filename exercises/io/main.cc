#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

int main (int argc, char *argv[]) {
  std::vector<double> numbers;
  std::string infile, outfile;

  for(int i=1; i<argc; i++){
    std::string arg = argv[i];
    if(arg=="-n" && i+1<argc){
      numbers.push_back(std::stod(argv[++i]));
    }
    else if(arg=="--input" && i+1<argc){
      infile=argv[++i];
    }
    else if(arg=="--output" && i+1<argc){
      outfile=argv[++i];
    } 
    else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return 1;
    }
  }

  std::ifstream myinput(infile);
  std::ofstream myoutput(outfile);

  double x;
  if(myinput.is_open() && myoutput.is_open()) {
    while(myinput >> x){
      myoutput << x << " " << std::sin(x) << " " << std::cos(x) << "\n";
    }
  }
  else if (!infile.empty() && !myinput){
    std::cerr << "Error opening input file: " << infile << "\n";
    return 1;
  }
  if(!numbers.empty()){
    for(auto n: numbers) std::cout << n << " " << std::sin(n) << " " << std::cos(n) << "\n";
  }
  else while(std::cin >> x) std::cout << x << " " << std::sin(x) << " " << std::cos(x) << "\n";

  return 0;
}

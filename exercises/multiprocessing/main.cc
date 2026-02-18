#include <iostream>
#include <vector>
#include <thread>
#include <string>

struct datum {long int start,end; double sum;};

void harmonic(datum& p){
  long int start=p.start, end = p.end;
  double sum = 0;
  for(long int i=start; i<end; i++) sum+= 1.0/i;
  p.sum = sum;
}

int main(int argc, char** argv){
  long int nterms=(long int)1e9, nthreads=1;
  for(int i=0; i<argc; i++){
    std::string arg = argv[i];
    if(arg=="-terms" && i+1<argc)nterms = (long int)std::stod(argv[++i]);
    if(arg=="-threads" && i+1<argc)nthreads = std::stoi(argv[++i]);
  }
  std::cout << "threads: " << nthreads << "\n";
  std::cout << "terms: " << (float)nterms << "\n";
  std::vector< std::thread > threads;
  threads.reserve(nthreads);
  std::vector<datum> data(nthreads);
  for(int i=0; i<nthreads; ++i){
    data[i].start = 1 + (nterms/nthreads)*i;
    data[i].end = 1 + (nterms/nthreads)*(i+1);
    threads.emplace_back(harmonic, std::ref(data[i]));
  }
  for(std::thread &thread : threads)thread.join();
  double total = 0;
  for(datum &d : data)total+=d.sum;
  std::cout << "total sum= " << total << "\n";
}

#include <cctype>
#include <vector>
#include <functional>
#include <iostream>
#include "core/vector.h"
#include "ode.h"

std::function<pp::vector(double, pp::vector)> get_system(int id, double epsilon){
  if(id == 1){
    return [](double x, pp::vector y){
      pp::vector dydx(2);
      dydx[0] = y[1];
      dydx[1] = -0.2*y[1] -y[0];
      return dydx;
    };
  } 
  if(id == 2){
    return [](double x, pp::vector y){
      pp::vector dydx(2);
      dydx[0] = 1.5*y[0] - y[0]*y[1];
      dydx[1] = -3*y[1] + y[0]*y[1];
      return dydx;
    };
  }
  if(id == 3){
    return [=](double x, pp::vector y){
      pp::vector dydx(2);
      dydx[0] = y[1];
      dydx[1] = 1 - y[0] + epsilon*y[0]*y[0];
      return dydx;
    };
  }
  if(id == 4){
    return [](double x, pp::vector y){
      pp::vector dydx(12);
      double G = 1;
      auto pos = [&](int i){return pp::vector{y[i*2], y[i*2+1]};};
      auto vel = [&](int i){return pp::vector{y[i*2+6], y[i*2+7]};};

      for(int i = 0; i<3; i++){
        pp::vector acc(2);
        for(int j = 0; j<3; j++){
          if(i == j) continue;
          pp::vector r_ij = pos(j) - pos(i);
          double dist = r_ij.norm();
          double force = G / (dist*dist*dist+1e-6);
          acc[0] += (force*r_ij[0]);
          acc[1] += (force*r_ij[1]);
        }
        dydx[i*2] = y[i*2+6];
        dydx[i*2+1] = y[i*2+7];
        dydx[i*2+6] = acc[0];
        dydx[i*2+7] = acc[1];
      }
      return dydx;
    };
  } else return nullptr;
}

int main(int argc, char* argv[]){
  int id = 0;
  double acc = 1e-3, eps = 1e-3, h = 0.1;
  double epsilon = 0, a = 0, b = 10;
  std::vector<double> yinit;
  for(int i = 1; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "-id" && i+1<argc) id = std::stoi(argv[++i]);
    if(arg == "-acc" && i+1<argc) acc = std::stod(argv[++i]);
    if(arg == "-eps" && i+1<argc) eps = std::stod(argv[++i]);
    if(arg == "-h" && i+1<argc) h = std::stod(argv[++i]);
    if(arg == "-a" && i+1<argc) a = std::stod(argv[++i]);
    if(arg == "-b" && i+1<argc) b = std::stod(argv[++i]);
    if(arg == "-epsilon" && i+1<argc) epsilon = std::stod(argv[++i]);
    if(arg == "-y" && i+1<argc){
      while(i+1<argc){
        std::string next_val = argv[i+1];
        if(next_val[0] == '-' && next_val.size() > 1 && !std::isdigit(next_val[1])){
          break;
        }
        yinit.push_back(std::stod(argv[++i]));
      }
    };
  }
  pp::vector y_start(yinit.size());
  for(int i = 0; i<yinit.size(); i++) y_start[i] = yinit.data()[i];
  pp::rk_45 stepper;
  auto res = pp::driver(stepper, get_system(id, epsilon), a, b, y_start, h, acc, eps);
  auto xdata = res.x_values;
  auto ydata = res.y_values;
  for(int i = 0; i<xdata.size(); i++){
    std::cout << xdata[i];
    for(int j = 0; j<ydata[i].size(); j++){
      std::cout << " " << ydata[i][j];
    }
    std::cout << "\n";
  }
  return 0;
}

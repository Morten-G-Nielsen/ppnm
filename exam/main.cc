#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>
#include "core/vector.h"
#include "root.h"


using namespace pp;
int main (int argc, char *argv[]) {
  int id = 0;
  for(size_t i = 0; i<argc; i++){
    std::string arg = argv[i];
    if(arg == "-id" && i+1<argc) id = std::stoi(argv[++i]);
  } 
  if(id == 0){
    auto chained_system = [](const vector& x){
      size_t n = x.size();
      vector fx(n);
      for(size_t i = 0; i<n-1; i++){
        fx[i] = x[i]*x[i] - x[i+1];
      }
      fx[n-1] = x[n-1]*x[n-1] - 1.0;
      return fx;
    };
    for(size_t N = 2; N<=202; N += 5){
      vector x0(N); for(size_t i = 0; i<N; i++) x0[i] = 2.0;

      BacktrackingSearch s;
      GoodBroydenUpdate u1(N);
      BadBroydenUpdate u2(N);
      NewtonUpdate u3(N);

      FindRoot solver1(u1, s);
      FindRoot solver2(u2, s);
      FindRoot solver3(u3, s);

      auto start = std::chrono::high_resolution_clock::now();
      vector root1 = solver1.solve(chained_system, x0, 1e-7, 1e-3);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> duration1 = end - start;

      start = std::chrono::high_resolution_clock::now();
      vector root2 = solver2.solve(chained_system, x0, 1e-7, 1e-3);
      end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> duration2 = end - start;

      start = std::chrono::high_resolution_clock::now();
      vector root3 = solver3.solve(chained_system, x0, 1e-7, 1e-3);
      end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> duration3 = end - start;

      bool found1 = (chained_system(root1).norm() < 1e-7);
      bool found2 = (chained_system(root2).norm() < 1e-7);
      bool found3 = (chained_system(root3).norm() < 1e-7);

      std::cout << N << " ";
      if(found1){
        std::cout << duration1.count() << " ";
      } else {
        std::cout << "NaN ";
      }
      if(found2){
        std::cout << duration2.count() << " ";
      } else {
        std::cout << "NaN ";
      }
      if(found3){
        std::cout << duration3.count() << "\n";
      } else {
        std::cout << "NaN\n";
      }
    }
  } else if(id == 1){
    size_t evals = 0;
    auto chained_system = [&](const vector& x){
      size_t n = x.size();
      vector fx(n);
      for(size_t i = 0; i<n-1; i++){
        fx[i] = x[i]*x[i] - x[i+1];
      }
      fx[n-1] = x[n-1]*x[n-1] - 1.0;
      evals += 1;
      return fx;
    };
    for(size_t N = 2; N<=202; N += 5){
      vector x0(N); for(size_t i = 0; i<N; i++) x0[i] = 2.0;

      BacktrackingSearch s;
      GoodBroydenUpdate u1(N);
      BadBroydenUpdate u2(N);
      NewtonUpdate u3(N);

      FindRoot solver1(u1, s);
      FindRoot solver2(u2, s);
      FindRoot solver3(u3, s);

      vector root1 = solver1.solve(chained_system, x0, 1e-7, 1e-3);
      size_t evals1 = evals; evals = 0;
      vector root2 = solver2.solve(chained_system, x0, 1e-7, 1e-3);
      size_t evals2 = evals; evals = 0;
      vector root3 = solver3.solve(chained_system, x0, 1e-7, 1e-3);
      size_t evals3 = evals; evals = 0;

      bool found1 = (chained_system(root1).norm() < 1e-7);
      bool found2 = (chained_system(root2).norm() < 1e-7);
      bool found3 = (chained_system(root3).norm() < 1e-7);

      std::cout << N << " ";
      if(found1){
        std::cout << evals1 << " ";
      } else {
        std::cout << "NaN ";
      }
      if(found2){
        std::cout << evals2 << " ";
      } else {
        std::cout << "NaN ";
      }
      if(found3){
        std::cout << evals3 << "\n";
      } else {
        std::cout << "NaN\n";
      }
    }
  } else if(id == 2){
    auto stress_system = [](const vector& x){
      return vector{x[0]*x[0]+x[1]*x[1]-4, std::sin(x[0]+x[1]*x[1])-0.5};
    };

    BacktrackingSearch s;
    BadBroydenUpdate u(2);
    FindRoot solver(u, s);

    for(double x_start = -3.0; x_start<=3.05; x_start += 0.12){
      for(double y_start = -3.0; y_start<=3.05; y_start += 0.12){
        vector x0 = {x_start, y_start};

        size_t iterations = solver.solve_count_iters(stress_system, x0, 1e-7, 1e-3,100,false);

        std::cout << x_start << " " << y_start << " " << iterations << "\n";
      }
      std::cout << "\n";
    }
  } else if(id == 3){
    auto stress_system = [](const vector& x){
      return vector{x[0]*x[0]+x[1]*x[1]-4, std::sin(x[0]+x[1]*x[1])-0.5};
    };

    BacktrackingSearch s;
    GoodBroydenUpdate u(2);
    FindRoot solver(u, s);

    for(double x_start = -3.0; x_start<=3.05; x_start += 0.12){
      for(double y_start = -3.0; y_start<=3.05; y_start += 0.12){
        vector x0 = {x_start, y_start};

        size_t iterations = solver.solve_count_iters(stress_system, x0, 1e-7, 1e-3, 100, false);

        std::cout << x_start << " " << y_start << " " << iterations << "\n";
      }
      std::cout << "\n";
    }
  } else if(id == 4){
    auto rosenbrock_root = [](const vector& x){
      return vector{10.0*(x[1]-x[0]*x[0]), 1.0-x[0]};
    };
    std::vector<vector> backtrack_mem;
    std::vector<vector> quadratic_mem;
    
    BacktrackingSearch s1;
    QuadraticInterpSearch s2;
    GoodBroydenUpdate u(2);

    FindRoot solver1(u, s1);
    FindRoot solver2(u, s2);

    vector x = {-2.0, 2.0};
    backtrack_mem.push_back(x);
    for(size_t i = 0; i<100; i++){
      x = solver1.solve(rosenbrock_root, x, 1e-7, 1e-3, 1);
      backtrack_mem.push_back(x);
      if(rosenbrock_root(x).norm() < 1e-7) break;
    }
    u.reset();
    x = {-2.0, 2.0};
    quadratic_mem.push_back(x);
    for(size_t i = 0; i<100; i++){
      x = solver2.solve(rosenbrock_root, x, 1e-7, 1e-3, 1);
      quadratic_mem.push_back(x);
      if(rosenbrock_root(x).norm() < 1e-7) break;
    }

    size_t max_points = std::max(backtrack_mem.size(), quadratic_mem.size());
    for(size_t i = 0; i<max_points; i++){

      if(i < backtrack_mem.size()){
        std::cout << backtrack_mem[i][0] << " " << backtrack_mem[i][1] << " ";
      } else {
        std::cout << "NaN NaN ";
      }

      if(i < quadratic_mem.size()){
        std::cout << quadratic_mem[i][0] << " " << quadratic_mem[i][1] << "\n";
      } else {
        std::cout << "NaN NaN\n";
      }
    }
  } else if(id == 5){
    size_t evals = 0;
    auto chain = [&](const vector& x){
      size_t n = x.size();
      size_t N = n/2;
      vector f(n);

      double k = 100.0, g=9.81, m=1.0;
      double L0 = 11/(N+1.0);

      double x_left = -5.0, y_left = 0.0;
      double x_right = 5.0, y_right = 0.0;

      for(size_t i = 0; i<N; i++){
        double xi = x[2*i];
        double yi = x[2*i+1];

        double x_prev = (i == 0) ? x_left : x[2*i-2];
        double y_prev = (i == 0) ? y_left : x[2*i-1];

        double x_next = (i == N-1) ? x_right : x[2*i+2];
        double y_next = (i == N-1) ? y_right : x[2*i+3];

        double dx_left = xi - x_prev;
        double dy_left = yi - y_prev;
        double dx_right = x_next - xi;
        double dy_right = y_next - yi;

        double r_left = std::sqrt(dx_left*dx_left + dy_left*dy_left);
        double r_right = std::sqrt(dx_right*dx_right + dy_right*dy_right);

        f[2*i] = -k*(1-L0/r_left)*dx_left + k*(1-L0/r_right)*dx_right;
        f[2*i+1] = -k*(1-L0/r_left)*dy_left + k*(1-L0/r_right)*dy_right - m*g;
      }
      evals += 1;
      return f;
    };
    size_t N = 100;
    vector x0(2*N);
    for(size_t i = 0; i<N; i++){
        x0[2*i] = -5.0 + 10.0/(N + 1.0)*(i+1.0);
        x0[2*i+1] = -1.0;
    }
    QuadraticInterpSearch s;
    GoodBroydenUpdate u(2*N);
    FindRoot solver(u, s);
    vector root = solver.solve(chain, x0, 1e-7, 1e-3);

    for(size_t i = 0; i<N; i++){
      std::cout << root[2*i] << " " << root[2*i+1] << "\n";
    }
    std::cout << "\nNorm of force: " << chain(root).norm() << "\n";
  }
  return 0;
}

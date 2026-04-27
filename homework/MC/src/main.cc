#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>
#include "mc.h"

using namespace pp;
struct TestValues{
  double mean_integral = 0;
  double mean_error = 0;
  double mean_true_error = 0;
  double std_dev_integral = 0;
  double std_dev_error = 0;
  double exact;
  size_t trials;
};
TestValues run_test(MCRule& rule, RNG& gen, int seed,
    const std::function<double(const std::vector<double>&)>& f,
    const std::vector<double>& a, const std::vector<double>& b,
    size_t N, size_t trials, double exact){

    std::vector<double> integrals;
    std::vector<double> errors;
    double sum_integrals = 0;
    double sum_errors = 0;
    double sum_true_error = 0;
    for(size_t i = 0; i<trials; i++){
      gen.set_seed(seed + i);
      MCIntegrator integrator(rule);
      auto res = integrator.integrate(f, a, b, N);
      integrals.push_back(res.integral);
      errors.push_back(res.error);
      sum_integrals += res.integral;
      sum_errors += res.error;
      sum_true_error += std::abs(res.integral - exact);
    }
    double mean_integral = sum_integrals/trials;
    double mean_error = sum_errors/trials;
    double mean_true_error = sum_true_error/trials;
    double sum_sq_diff = 0;
    double sum_sq_diff_error = 0;
    for(const double r : integrals){
      sum_sq_diff += (r - mean_integral)*(r - mean_integral);
    }
    for(const double e : errors){
      sum_sq_diff_error += (e - mean_error)*(e - mean_error);
    }
    double std_dev_integral = std::sqrt(sum_sq_diff/(trials - 1));
    double std_dev_error = std::sqrt(sum_sq_diff_error/(trials - 1));
    return {mean_integral, mean_error, mean_true_error,
        std_dev_integral, std_dev_error, exact, trials};
}

void print_res(TestValues& res){
    std::cout << std::setprecision(10) << "Trials:                   " << res.trials
      << "\nExact result:             " << res.exact
      << "\nMean integral:            " << res.mean_integral << " +/- " << res.std_dev_integral
      << "\nMean estimated error:     " << res.mean_error << " +/- " << res.std_dev_error
      << "\nMean true error:          " << res.mean_true_error
      << "\nTrueError/EstimatedError: " << res.mean_true_error/res.mean_error << "\n";
}

int main(int argc, char* argv[]) {
    size_t N = 100;
    int seed = 12345;
    int id = 1;
    for(int i = 0; i<argc; i++){
      std::string arg = argv[i];
      if(arg == "-N" && i+1<argc) N = static_cast<size_t>(std::stod(argv[++i]));
      if(arg == "-seed" && i+1<argc) seed = std::stoi(argv[++i]);
      if(arg == "-id" && i+1<argc) id = std::stoi(argv[++i]);
    }

    if(id == 1){
      std::vector<double> a = {0.0, 0.0, 0.0};
      std::vector<double> b = {M_PI, M_PI, M_PI};

      double PI3 = M_PI*M_PI*M_PI;
      double exact = std::pow(std::tgamma(1.0/4), 4)/(4*PI3);
      auto f = [=](const std::vector<double>& x) {
        return 1/(1-std::cos(x[0])*std::cos(x[1])*std::cos(x[2]))/PI3;
      };
      auto g = [](const std::vector<double>& x){
        double val = x[0]*x[0] + x[1]*x[1]/4 + x[2]*x[2]/9;
        return (val <= 1) ? 1.0 : 0.0;
      };
      std::vector<double> a_eli = {-1.0, -2.0, -3.0}, b_eli = {1.0, 2.0, 3.0};
      double volume = 8*M_PI;

      size_t trials = 20;

      LCG gen1(seed);
      MCPlain rule1(gen1);
      auto res = run_test(rule1, gen1, seed, g, a_eli, b_eli, N, trials, volume);
      auto res1 = run_test(rule1, gen1, seed, f, a, b, N, trials, exact);

      StdRNG gen2(seed);
      MCPlain rule2(gen2);
      auto res2 = run_test(rule2, gen2, seed, f, a, b, N, trials, exact);

      Halton hal;
      Halton hal2(10000);
      MCQuasi rule3(hal, hal2);
      MCIntegrator integrator(rule3);
      auto res3 = integrator.integrate(f, a, b, N);

      Lattice lat;
      MCQuasi rule4(hal, lat);
      integrator.set_rule(rule4);
      auto res4 = integrator.integrate(f, a, b, N);

      gen1.set_seed(seed);
      MCStrat rule5(gen1);
      auto res5 = run_test(rule5, gen1, seed, f, a, b, N, trials, exact);
    
      std::cout << "---Plain Monte Carlo on elipsoid---\n";
      print_res(res);
      std::cout << "\n\n---Comparison between methods on singular integral---\n";
      std::cout << "---Plain Monte Carlo with LCG---\n";
      print_res(res1);
      std::cout << "---Plain Monte Carlo with StdRNG---\n";
      print_res(res2);
      std::cout << "---Quasi Monte Carlo with Halton---\n" << std::setprecision(10)
        << "Exact result:             " << exact
        << "\nIntegral:                 " << res3.integral
        << "\nEstimated error:          " << res3.error 
        << "\nTrue error:               " << std::abs(res3.integral - exact)
        << "\nTrueError/EstimatedError: " << std::abs(res3.integral - exact)/res3.error 
        << "\n";
      std::cout << "---Quasi Monte Carlo with Halton and Lattice---\n"
        << std::setprecision(10)
        << "Exact result:             " << exact
        << "\nIntegral:                 " << res4.integral
        << "\nEstimated error:          " << res4.error 
        << "\nTrue error:               " << std::abs(res4.integral - exact)
        << "\nTrueError/EstimatedError: " << std::abs(res4.integral - exact)/res4.error 
        << "\n";
      std::cout << "---Stratified Monte Carlo with LCG---\n";
      print_res(res5);
    } else if(id == 2){
      std::vector<double> a{0.0, 0.0}, b{1.0, 1.0};
      double exact = 2.0/3.0;
      auto f = [](const std::vector<double>& x){
        return x[0]*x[0] + x[1]*x[1];
      };
      size_t trials = 20;

      LCG gen1(seed);
      MCPlain rule1(gen1);
      auto res1 = run_test(rule1, gen1, seed, f, a, b, N, trials, exact);

      Halton hal;
      Lattice lat;
      MCQuasi rule2(hal, lat);
      MCIntegrator integrator(rule2);
      auto res2 = integrator.integrate(f, a, b, N);

      gen1.set_seed(seed);
      MCStrat rule3(gen1);
      auto res3 = run_test(rule3, gen1, seed, f, a, b, N, trials, exact);
      std::cout << N << " " << res1.mean_true_error << " " 
        << std::abs(res2.integral - exact)
        << " " << res3.mean_true_error << "\n";
    }
    return 0;
}

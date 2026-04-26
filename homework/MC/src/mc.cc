#include <cmath>
#include <vector>
#include "mc.h"

namespace pp{
double LCG::next(){
  // Learned this from Gemeni where a overflow is used instead of modolus
  state = a*state + c;
  return static_cast<double>(state) / static_cast<double>(~0ULL);
}
void LCG::set_seed(unsigned int seed){
  state = seed;
}

double StdRNG::next(){
  return dist(gen);
}
void StdRNG::set_seed(unsigned int seed){
  gen.seed(seed);
}

const std::vector<int>& Halton::get_prime(size_t max_dim){
  static std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19};
  int candidate = primes.back() + 2;
  while(primes.size() < max_dim){
    bool is_prime = true;
    for(int p : primes){
      if(p*p > candidate) break;
      if(candidate % p == 0){is_prime=false; break;}
    }
    if(is_prime) primes.push_back(candidate);
    candidate += 2;
  }
  return primes;
}

double Halton::get(int index, size_t dim){
  const auto& primes = get_prime(dim);
  int base = primes[dim];
  double q = 0.0;
  double bk = 1.0/base;
  while(index > 0){
    q += bk*(index % base);
    index /= base;
    bk /= base;
  }
  return q;
}

MCResult MCPlain::integrate(const std::function<double(const std::vector<double>&)>& f,
    const std::vector<double>& a,
    const std::vector<double>& b,
    size_t N) const {
    if(N <= 0) return {0.0, 0.0}; 
    double sum = 0;
    double sum_sq = 0;
    double V = 1;
    size_t dim = a.size();
    for(size_t i = 0; i<dim; i++) V*=std::abs(b[i] - a[i]);

    std::vector<double> x(dim);

    for(size_t i = 0; i<N; i++){
      for(size_t j = 0; j<dim; j++){ 
        x[j] = (a[j]+rng.next()*(b[j]-a[j]));
      }
      double fx = f(x);
      sum += fx;
      sum_sq += fx*fx;
    }
    double mean = sum/N;
    double var = std::max(0.0, (sum_sq/N)-mean*mean);
    return {mean*V, std::sqrt(var/N)*V};
    }


MCResult MCHalton::integrate(const std::function<double(const std::vector<double>&)>& f,
    const std::vector<double>& a,
    const std::vector<double>& b,
    size_t N) const {
  std::pair<double, double> sum = {0.0, 0.0};
  double V = 1;
  size_t dim = a.size();
  for(size_t i = 0; i<dim; i++) V*=std::abs(b[i]-a[i]);
  std::vector<double> x(dim);
  std::vector<double> y(dim);
  Halton gen;

  for(size_t i = 0; i<N; i++){
    for(size_t j = 0; j<dim; j++){
      x[j] = a[j]+gen.get(i+1, j)*(b[j]-a[j]);
      y[j] = a[j]+gen.get(i+100000, j)*(b[j]-a[j]);
    }
    sum.first += f(x);
    sum.second += f(y);
  }
  double integral1 = sum.first/N*V;
  double integral2 = sum.second/N*V;
  double err = std::abs(integral1-integral2);
  return {(integral1+integral2)/2.0, err};
}

MCResult MCStrat::integrate(const std::function<double(const std::vector<double>&)>& f,
    const std::vector<double>& a,
    const std::vector<double>& b,
    size_t N) const {
    std::vector<double> mut_a = a;
    std::vector<double> mut_b = b;
    std::vector<double> points(a.size());
    std::vector<std::pair<std::vector<double>, double>> samples;
    MCPlain plain(rng);
    return recursive(f, mut_a, mut_b, N, plain, 0, points, samples);
}

MCResult MCStrat::recursive(const std::function<double(const std::vector<double>&)>& f,
    std::vector<double>& a,
    std::vector<double>& b,
    size_t N, MCRule& plain, size_t depth, std::vector<double>& points,
    std::vector<std::pair<std::vector<double>, double>>& samples) const {

  size_t dim = a.size();
  size_t nmin = 16*dim;
  if(N < 3*nmin || depth>=15){
    // Plain MC with N
    return plain.integrate(f, a, b, N);
  }


  // Finds the sample points and stores them
  if(samples.size() < nmin) samples.resize(nmin);
  for(size_t i = 0; i<nmin; i++){
    for(size_t j = 0; j<dim; j++){
      points[j] = (a[j]+rng.next()*(b[j]-a[j]));
    }
    double val = f(points);
    samples[i] = {points, val};
  }

  size_t best_dim = 0;
  double max_diff = -1.0;
  for(size_t d = 0; d<dim; d++){
    double mid = (a[d] + b[d])/2;
    double sum_left = 0.0, sum_right = 0.0;
    size_t count_left = 0, count_right = 0;

    for(size_t i = 0; i<nmin; i++){
      if(samples[i].first[d] < mid){
        sum_left += samples[i].second;
        count_left++;
      } else {
        sum_right += samples[i].second;
        count_right++;
      }
    }

    double mean_left = (count_left > 0) ? (sum_left/count_left) : 0.0;
    double mean_right = (count_right > 0) ? (sum_right/count_right) : 0.0;

    double diff = std::abs(mean_left - mean_right);
    if(diff > max_diff){
      max_diff = diff;
      best_dim = d;
    }
  }

  // I have tried for multilple days to get the variance spliting of points
  // to work, but it always lead to the method being extremely unstable
  // where it performed worse then the plain MC and where you could not trust 
  // that increasing the point by a factor of 1000 would be better the by 100.
  // So i have given up and just split them evenly.
  size_t NRemaning = N - nmin;
  size_t NLeft = NRemaning / 2;
  size_t NRight = NRemaning - NLeft;

  double a_original = a[best_dim];
  double b_original = b[best_dim];
  double mid = (a[best_dim] + b[best_dim])/2.0;

  b[best_dim] = mid;
  MCResult res_left = recursive(f, a, b, NLeft, plain, depth+1, points, samples);
  b[best_dim] = b_original;
  a[best_dim] = mid;
  MCResult res_right = recursive(f, a, b, NRight, plain, depth+1, points, samples);
  a[best_dim] = a_original;

  return {res_left.integral + res_right.integral,
    std::sqrt(res_left.error*res_left.error + res_right.error*res_right.error)};
}
}

#pragma once
#include <random>
#include <vector>
#include <functional>

namespace pp{
struct RNG{
  virtual ~RNG() = default;
  virtual double next() = 0;
  virtual void set_seed(unsigned int) = 0;
};

struct LCG : public RNG{
  unsigned long long state;
  static constexpr unsigned long long a = 6364136223846793005ULL;
  static constexpr unsigned long long c = 1442695040888963407ULL;

  LCG(unsigned long long seed) : state(seed) {}
  double next() override;
  void set_seed(unsigned int) override;
};

struct StdRNG : public RNG{
  std::mt19937 gen;
  std::uniform_real_distribution<double> dist{0.0, 1.0};
  
  StdRNG(unsigned int seed) : gen(seed) {}
  double next() override;
  void set_seed(unsigned int) override;
};

static const std::vector<int>& get_prime(size_t max_dim);

struct QRNG{
  int offset;
  QRNG(int offset) : offset(offset) {}
  QRNG() : offset(0) {}

  virtual ~QRNG() = default;
  virtual double get(int index, size_t dim) = 0;
};

struct Halton : QRNG{
  Halton(int offset) : QRNG(offset) {}
  Halton() : QRNG() {}

  double get(int index, size_t dim) override;
};

struct Lattice : QRNG{
  Lattice(int offset) : QRNG(offset) {}
  Lattice() : QRNG() {}

  double get(int index, size_t dim) override;
};

struct MCResult{
  double integral;
  double error;
}; 

struct MCRule{
  virtual ~MCRule() = default;

  virtual MCResult integrate(const std::function<double(const std::vector<double>&)>& f,
      const std::vector<double>& a,
      const std::vector<double>& b,
      size_t N) const = 0;
};

struct MCPlain : public MCRule{
  RNG& rng;
  MCPlain(RNG& gen) : rng(gen) {}

  MCResult integrate(const std::function<double(const std::vector<double>&)>& f,
      const std::vector<double>& a,
      const std::vector<double>& b,
      size_t N) const override;
};

struct MCQuasi : public MCRule{
  QRNG& gen1;
  QRNG& gen2;
  MCQuasi(QRNG& gen1, QRNG& gen2) : gen1(gen1), gen2(gen2) {}

  MCResult integrate(const std::function<double(const std::vector<double>&)>& f,
      const std::vector<double>& a, 
      const std::vector<double>& b,
      size_t N) const override;
};

struct MCStrat : public MCRule{
  RNG& rng;
  MCStrat(RNG& gen) : rng(gen) {}

  MCResult recursive(const std::function<double(const std::vector<double>&)>& f,
      std::vector<double>& a,
      std::vector<double>& b,
      size_t N, MCRule& plain,
      size_t depth, std::vector<double>& points,
      std::vector<std::pair<std::vector<double>, double>>& samples) const;

  MCResult integrate(const std::function<double(const std::vector<double>&)>& f,
      const std::vector<double>& a,
      const std::vector<double>& b,
      size_t N) const override;
};

struct MCIntegrator{
  MCRule* rule;
  MCIntegrator(MCRule& r) : rule(&r) {}

  void set_rule(MCRule& r){
    rule = &r;
  }

  MCResult integrate(const std::function<double(const std::vector<double>&)>& f,
      const std::vector<double>& a,
      const std::vector<double>& b,
      size_t N){
    return rule->integrate(f, a, b, N);
  }
};
}

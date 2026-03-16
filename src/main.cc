#include "Vector.h"
#include "LennardJones.h"
#include "Optimizers.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <random>
#include <string>

struct Result {
  int n, m;
  double u_tot, u_at;
  bool operator<(const Result& other) const {
    if (std::abs(u_at - other.u_at) > 1e-8) return u_at < other.u_at;
    return u_tot < other.u_tot;
  }
};

void print_help(const char* prog) {
  std::cout << "Usage: " << prog << " [n_start [n_end [m_start [m_end]]]]\n"
            << "Example: " << prog << " 5 5 3 3 (runs only n=5, m=3)\n"
            << "Default: n=3..20, m=1..6\n";
}

int main(int argc, char* argv[]) {
  int n_start = 3, n_end = 20;
  int m_start = 1, m_end = 6;

  if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
    print_help(argv[0]);
    return 0;
  }

  if (argc >= 3) {
    n_start = std::stoi(argv[1]);
    n_end = std::stoi(argv[2]);
  }
  if (argc >= 5) {
    m_start = std::stoi(argv[3]);
    m_end = std::stoi(argv[4]);
  }

  std::vector<Result> results;
  std::mt19937 rng(42); // Fixed seed for reproducibility
  std::uniform_real_distribution<double> dist(-3.0, 3.0);

  for (int n = n_start; n <= n_end; ++n) {
    for (int m = m_start; m <= m_end; ++m) {
      LennardJonesSystem system(n, m);
      Vector x(n * m);
      for (int i = 0; i < n * m; ++i) {
        x[i] = dist(rng);
      }

      std::cout << "--- Solving for n=" << n << ", m=" << m << " ---" << std::endl;
      std::cout << "U initial: " << system.compute_total_energy(x) << std::endl;

      // Iterative refinement as in original code
      x = DynamicOptimizer::optimize(system, x);
      x = PowellOptimizer::optimize(system, x);
      x = DynamicOptimizer::optimize(system, x);
      x = PowellOptimizer::optimize(system, x);

      double u_tot = system.compute_total_energy(x);
      double u_at = system.compute_average_energy(x);
      std::cout << "Final Energy: U_tot=" << u_tot << ", U_at=" << u_at << std::endl;

      results.push_back({n, m, u_tot, u_at});
    }
  }

  FILE* out = fopen("output.txt", "w");
  std::sort(results.begin(), results.end());
  for (const auto& res : results) {
    fprintf(out, "%d %d %.4lf %.4lf\n", res.n, res.m, res.u_tot, res.u_at);
  }
  fclose(out);

  return 0;
}

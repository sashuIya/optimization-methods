#ifndef LENNARD_JONES_HPP
#define LENNARD_JONES_HPP

#include "Vector.h"
#include <vector>
#include <cmath>

class LennardJonesSystem {
public:
  LennardJonesSystem(int n, int m) : n_(n), m_(m) {}

  double compute_total_energy(const Vector& x) const {
    double total = 0.0;
    for (int i = 0; i < n_; ++i) {
      for (int j = i + 1; j < n_; ++j) {
        total += compute_pair_energy(x, i, j);
      }
    }
    return total;
  }

  double compute_average_energy(const Vector& x) const {
    return compute_total_energy(x) / n_;
  }

  Vector compute_gradient(const Vector& x) const {
    Vector grad(n_ * m_, 0.0);
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        if (i == j) continue;
        
        double d2 = squared_dist(x, i, j);
        double d7 = std::pow(d2, 3.5);
        double d3 = std::pow(d2, 1.5);

        // Gradient of U_ij = (1/r_ij)^12 - (1/r_ij)^6
        // U_ij = (d2)^-6 - (d2)^-3
        // dU/dx_ik = -6(d2)^-7 * (d(d2)/dx_ik) - (-3)(d2)^-4 * (d(d2)/dx_ik)
        // d(d2)/dx_ik = 2(x_ik - x_jk)
        
        // Original code: result[i][j] += -2.0 * (points[i][j] - points[index][j]) * (6.0 - 3.0 * pow(d, 3.0)) * pow(d, -7.0) / n;
        // Here d = sqrt(d2).
        // d^-7 = (d2)^-3.5
        // d^3 = (d2)^1.5
        
        double d = std::sqrt(d2);
        double factor = -2.0 * (6.0 - 3.0 * std::pow(d, 3.0)) * std::pow(d, -7.0) / n_;

        for (int k = 0; k < m_; ++k) {
          grad[i * m_ + k] += (x[i * m_ + k] - x[j * m_ + k]) * factor;
        }
      }
    }
    return grad;
  }

  int getN() const { return n_; }
  int getM() const { return m_; }

private:
  int n_, m_;

  double squared_dist(const Vector& x, int i, int j) const {
    double res = 0.0;
    for (int k = 0; k < m_; ++k) {
      double diff = x[i * m_ + k] - x[j * m_ + k];
      res += diff * diff;
    }
    return res;
  }

  double compute_pair_energy(const Vector& x, int i, int j) const {
    double d2 = squared_dist(x, i, j);
    return std::pow(d2, -6.0) - std::pow(d2, -3.0);
  }
};

#endif // LENNARD_JONES_HPP

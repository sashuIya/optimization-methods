#ifndef OPTIMIZATION_METHODS_SRC_LENNARD_JONES_H_
#define OPTIMIZATION_METHODS_SRC_LENNARD_JONES_H_

#include <cmath>
#include <vector>

#include "Vector.h"

class LennardJonesSystem {
 public:
  LennardJonesSystem(int n, int m) : n_(n), m_(m) {}

  double ComputeTotalEnergy(const Vector& x) const {
    double total = 0.0;
    for (int i = 0; i < n_; ++i) {
      for (int j = i + 1; j < n_; ++j) {
        total += ComputePairEnergy(x, i, j);
      }
    }
    return total;
  }

  double ComputeAverageEnergy(const Vector& x) const {
    return ComputeTotalEnergy(x) / n_;
  }

  Vector ComputeGradient(const Vector& x) const {
    Vector grad(n_ * m_, 0.0);
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        if (i == j) continue;

        double d2 = SquaredDist(x, i, j);
        double d = std::sqrt(d2);
        double factor =
            -2.0 * (6.0 - 3.0 * std::pow(d, 3.0)) * std::pow(d, -7.0) / n_;

        for (int k = 0; k < m_; ++k) {
          grad[i * m_ + k] += (x[i * m_ + k] - x[j * m_ + k]) * factor;
        }
      }
    }
    return grad;
  }

  int GetN() const { return n_; }
  int GetM() const { return m_; }

 private:
  int n_, m_;

  double SquaredDist(const Vector& x, int i, int j) const {
    double res = 0.0;
    for (int k = 0; k < m_; ++k) {
      double diff = x[i * m_ + k] - x[j * m_ + k];
      res += diff * diff;
    }
    return res;
  }

  double ComputePairEnergy(const Vector& x, int i, int j) const {
    double d2 = SquaredDist(x, i, j);
    return std::pow(d2, -6.0) - std::pow(d2, -3.0);
  }
};

#endif  // OPTIMIZATION_METHODS_SRC_LENNARD_JONES_H_

#ifndef OPTIMIZERS_HPP
#define OPTIMIZERS_HPP

#include "Vector.hpp"
#include "LennardJones.hpp"
#include <cmath>
#include <algorithm>
#include <functional>

class LineSearch {
public:
  static double gold_method(std::function<double(double)> f, double a, double c, double eps = 1e-8) {
    double r = (3.0 - std::sqrt(5.0)) / 2.0;
    double b = a + r * (c - a);
    double y_b = f(b);

    while (std::abs(c - a) > eps) {
      double d = c + r * (a - c);
      double y_d = f(d);

      if (y_d < y_b) {
        a = b;
        b = d;
        y_b = y_d;
      } else {
        c = a;
        a = d;
      }
    }
    return c;
  }

  static void get_lambda_segment(std::function<double(double)> f, double d, double& a, double& c) {
    double left = 0.0, middle = 0.0, right = 0.0;
    double y_curr = f(0.0);
    double y_next = y_curr;
    double r = (3. - std::sqrt(5.0)) / 2.;
    r = (1. - r) / r;

    do {
      y_curr = y_next;
      left = middle;
      middle = right;
      right += d;
      y_next = f(right);
      d *= r;
    } while (y_next < y_curr);

    a = left;
    c = right;
  }

  static double get_best_lambda(std::function<double(double)> f, double delta = 1e-6) {
    double a, c;
    get_lambda_segment(f, delta, a, c);
    double lambda1 = gold_method(f, a, c);

    get_lambda_segment(f, -delta, a, c);
    double lambda2 = gold_method(f, a, c);

    double best_lambda = 0.0;
    double best_U = f(0.0);
    
    double u1 = f(lambda1);
    if (u1 < best_U) {
      best_U = u1;
      best_lambda = lambda1;
    }

    double u2 = f(lambda2);
    if (u2 < best_U) {
      best_lambda = lambda2;
    }

    return best_lambda;
  }
};

class PowellOptimizer {
public:
  static Vector optimize(const LennardJonesSystem& system, Vector x, double eps = 1e-8) {
    int dim = x.size();
    std::vector<Vector> s(dim);
    for (int i = 0; i < dim; ++i) {
      s[i] = Vector(dim, 0.0);
      s[i][i] = 1.0;
    }

    double u_prev;
    int reset_count = 0;
    do {
      u_prev = system.compute_average_energy(x);

      for (int k = 0; k < 3 * dim; ++k) {
        Vector prev_x = x;
        for (int i = 0; i < dim; ++i) {
          auto f = [&](double lambda) { return system.compute_total_energy(x + lambda * s[i]); };
          double lambda = LineSearch::get_best_lambda(f);
          x += lambda * s[i];
        }

        for (int i = 0; i < dim - 1; ++i) s[i] = s[i + 1];
        s[dim - 1] = Vector::normed_residual(x, prev_x, eps);

        auto f = [&](double lambda) { return system.compute_total_energy(x + lambda * s[dim - 1]); };
        double lambda = LineSearch::get_best_lambda(f);
        x += lambda * s[dim - 1];

        if (++reset_count > dim) {
          for (int i = 0; i < dim; ++i) {
            s[i] = Vector(dim, 0.0);
            s[i][i] = 1.0;
          }
          reset_count = 0;
        }
      }
    } while (std::abs((system.compute_average_energy(x) - u_prev) / system.compute_average_energy(x)) > eps);

    return x;
  }
};

class DynamicOptimizer {
public:
  static Vector optimize(const LennardJonesSystem& system, Vector x, double eps = 1e-8) {
    int dim = x.size();
    double dt = 0.085;
    double beta = 0.998;
    Vector v(dim, 0.0);

    double u_prev;
    do {
      u_prev = system.compute_average_energy(x);
      Vector grad = system.compute_gradient(x);

      double g_norm = grad.norm();
      if (g_norm > 4.0) {
        grad /= (100.0 * g_norm);
      }

      v = v * beta - dt * grad;
      x += dt * v;

    } while (std::abs((system.compute_average_energy(x) - u_prev) / system.compute_average_energy(x)) > eps);

    return x;
  }
};

#endif // OPTIMIZERS_HPP

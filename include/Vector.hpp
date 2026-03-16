#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <cmath>

class Vector {
public:
  Vector() = default;
  explicit Vector(int size, double value = 0.0) : data_(size, value) {}
  Vector(const std::vector<double>& data) : data_(data) {}

  int size() const { return data_.size(); }
  double& operator[](int index) { return data_[index]; }
  const double& operator[](int index) const { return data_[index]; }

  Vector& operator+=(const Vector& other) {
    for (int i = 0; i < size(); ++i) data_[i] += other.data_[i];
    return *this;
  }

  Vector& operator-=(const Vector& other) {
    for (int i = 0; i < size(); ++i) data_[i] -= other.data_[i];
    return *this;
  }

  Vector& operator*=(double scalar) {
    for (int i = 0; i < size(); ++i) data_[i] *= scalar;
    return *this;
  }

  Vector& operator/=(double scalar) {
    for (int i = 0; i < size(); ++i) data_[i] /= scalar;
    return *this;
  }

  friend Vector operator+(Vector lhs, const Vector& rhs) {
    lhs += rhs;
    return lhs;
  }

  friend Vector operator-(Vector lhs, const Vector& rhs) {
    lhs -= rhs;
    return lhs;
  }

  friend Vector operator*(Vector lhs, double scalar) {
    lhs *= scalar;
    return lhs;
  }

  friend Vector operator*(double scalar, Vector rhs) {
    rhs *= scalar;
    return rhs;
  }

  friend Vector operator/(Vector lhs, double scalar) {
    lhs /= scalar;
    return lhs;
  }

  double norm_sq() const {
    double res = 0.0;
    for (double x : data_) res += x * x;
    return res;
  }

  double norm() const { return std::sqrt(norm_sq()); }

  static Vector normed_residual(const Vector& a, const Vector& b, double eps = 1e-8) {
    Vector res = a - b;
    double n = res.norm();
    if (n > eps) res /= n;
    return res;
  }

private:
  std::vector<double> data_;
};

#endif // VECTOR_HPP

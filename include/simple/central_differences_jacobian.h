#ifndef SIMPLE_CENTRAL_DIFFERENCES_JACOBIAN_H_
#define SIMPLE_CENTRAL_DIFFERENCES_JACOBIAN_H_

#include <Eigen/Dense>

namespace simple {

template <typename Data>
struct CentralDifferencesJacobian {
  using argument_type = Eigen::Vector3i;
  using value_type = Eigen::Matrix2d;

  Eigen::Vector2d spacing;
  Data data;

  CentralDifferencesJacobian(const Eigen::Vector2d& spacing_, Data data_)
  : spacing(spacing_), data(data_) {}

  int dimension() const { return 3; }
  value_type zero() const { return Eigen::Matrix2d::Zero(); }
  bool in_domain(const argument_type &x) const { return data.in_domain(x); }

  bool operator()(const argument_type &x, Eigen::Ref<value_type> value) const {
    value = zero();
    if (!in_domain(x))
      return false;

    auto value_left = data.zero();
    auto value_right = data.zero();
    Eigen::Vector3i y = x;
    for (int j = 0; j < 2; ++j) {
      int distance = 2;
      y[j] = x[j] - 1;
      if (!data.in_domain(y)) {
        y[j] = x[j];
        distance = 1;
      }
      if (!data(y, value_left)) {
        return false;
      }

      y[j] = x[j] + 1;
      if (!data.in_domain(y)) {
        y[j] = x[j];
        distance = 1;
      }
      if (!data(y, value_right)) {
        return false;
      }

      value.col(j) = (value_right - value_left) / (distance * spacing[j]);
      y[j] = x[j];
    }

    return true;
  }
};

template <typename Data>
CentralDifferencesJacobian<Data> central_differences_jacobian(const Eigen::Vector2d& spacing, Data data) {
  return CentralDifferencesJacobian<Data>(spacing, data);
}

}

#endif //SIMPLE_CENTRAL_DIFFERENCES_JACOBIAN_H_

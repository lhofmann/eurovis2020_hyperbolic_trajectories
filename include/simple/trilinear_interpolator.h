#ifndef SIMPLE_TRILINEAR_INTERPOLATOR_H_
#define SIMPLE_TRILINEAR_INTERPOLATOR_H_

#include <Eigen/Dense>

namespace simple {

template <typename Data>
struct TrilinearInterpolator {
  typedef Eigen::Vector3d argument_type;
  typedef typename Data::value_type value_type;

  Eigen::Vector3i dimensions;
  Eigen::Vector3d grid_min, grid_max, spacing;
  Data data;

  TrilinearInterpolator(
      const Eigen::Vector3i& dimensions_,
      const Eigen::Vector3d& grid_min_,
      const Eigen::Vector3d& grid_max_,
      Data data_
  )
  : dimensions(dimensions_),
    grid_min(grid_min_),
    grid_max(grid_max_),
    data(data_)
  {
    spacing = (grid_max - grid_min).cwiseQuotient((dimensions.array() - 1).template cast<double>().matrix());
  }

  int dimension() const { return 3; }
  value_type zero() const { return data.zero(); }

  bool in_domain(const argument_type &x) const {
    return (x.array() >= grid_min.array()).all() && (x.array() <= grid_max.array()).all();
  }

  bool operator()(const argument_type &x, Eigen::Ref<value_type> value) const {
    value = zero();
    if (!in_domain(x))
      return false;

    Eigen::Vector3d x_local;
    Eigen::Vector3i cell_index;

    x_local = (x - grid_min).cwiseQuotient(spacing);
    cell_index = x_local.array().floor().template cast<int>();
    x_local -= cell_index.template cast<double>();

    auto node_value = data.zero();
    for (unsigned int i = 0; i < 8; i++) {
      double coeff = 1.;
      Eigen::Vector3i node_index = cell_index;

      for (unsigned int j = 0; j < 3; j++) {
        if ((i & (1U << j)) != 0) {
          coeff *= x_local[j];
          node_index[j]++;
        } else {
          coeff *= (1. - x_local[j]);
        }
      }

      if (data(node_index, node_value)) {
        value += coeff * node_value;
      }
    }

    return true;
  }
};

template <typename Data>
TrilinearInterpolator<Data> trilinear_interpolator(
    const Eigen::Vector3i& dimensions,
    const Eigen::Vector3d& grid_min,
    const Eigen::Vector3d& grid_max,
    Data data
) {
  return TrilinearInterpolator<Data>(dimensions, grid_min, grid_max, data);
}

}


#endif //SIMPLE_TRILINEAR_INTERPOLATOR_H_

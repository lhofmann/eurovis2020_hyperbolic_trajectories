#ifndef SIMPLE_VTK_GRID_3D_H_
#define SIMPLE_VTK_GRID_3D_H_

#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vector>
#include <cassert>

namespace simple {

struct vtk_grid_3d {
  using argument_type = Eigen::Vector3i;
  using value_type = Eigen::Vector2d;

  Eigen::Vector3i dimensions;
  vtkSmartPointer<vtkDataArray> data;

  vtk_grid_3d(const Eigen::Vector3i& dimensions_, vtkDataArray* data_)
  : dimensions(dimensions_), data(data_) {}

  int dimension() const { return 3; }
  value_type zero() const { return Eigen::Vector2d::Zero(); }
  bool in_domain(const argument_type &x) const {
    return (x.array() >= 0).all() && (x.array() < dimensions.array()).all();
  }

  bool operator()(const argument_type &x, Eigen::Ref<value_type> value) const {
    value = zero();
    if (!in_domain(x))
      return false;

    long index = x[0] + dimensions[0] * (x[1] + dimensions[1] * x[2]);
    std::vector<double> tuple(data->GetNumberOfComponents());
    data->GetTuple(index, &tuple[0]);

    assert(tuple.size() >= 2);
    value[0] = tuple[0];
    value[1] = tuple[1];

    return true;
  }
};

}

#endif //SIMPLE_VTK_GRID_3D_H_

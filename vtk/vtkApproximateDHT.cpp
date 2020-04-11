#include "vtkApproximateDHT.h"

#include "hyperbolic_trajectories/approximate_dht.h"
#include "hyperbolic_trajectories/continuous_svd.h"
#include "simple/central_differences_jacobian.h"
#include "simple/system_adapter.h"
#include "simple/trilinear_interpolator.h"
#include "simple/vtk_grid_3d.h"
#include "util/default_request_data_object.h"

#include <vtkCommand.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkSmartPointer.h>
#include <vtkDataObject.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>

#include <vtkImageData.h>
#include <vtkPolyData.h>

vtkStandardNewMacro(vtkApproximateDHT);

vtkApproximateDHT::vtkApproximateDHT() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

vtkApproximateDHT::~vtkApproximateDHT() = default;

int vtkApproximateDHT::FillInputPortInformation(int port, vtkInformation *info) {
  if (port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  } else if (port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  }
  return 1;
}

int vtkApproximateDHT::FillOutputPortInformation(int port, vtkInformation *info) {
  if (port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  }
  return 1;
}

int vtkApproximateDHT::RequestDataObject(
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector ** vtkNotUsed(inputVector),
    vtkInformationVector *outputVector) {
  util::default_request_data_object<vtkPolyData>(this, outputVector, 0);
  return 1;
}

int vtkApproximateDHT::RequestInformation(
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {
  return 1;
}

int vtkApproximateDHT::RequestUpdateExtent(
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {
  return 1;
}

int vtkApproximateDHT::RequestData(
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {
  vtkImageData *input = vtkImageData::GetData(inputVector[0], 0);
  vtkPolyData *input_seeds = vtkPolyData::GetData(inputVector[1], 0);
  vtkPolyData *output = vtkPolyData::GetData(outputVector, 0);

  if (!input || !input_seeds || !output)
    return 0;

  vtkDataArray* array = this->GetInputArrayToProcess(0, input);
  if (!array || array->GetNumberOfComponents() < 2)
    return 0;

  Eigen::Vector3i dimensions;
  Eigen::Vector3d grid_min, grid_max, spacing;

  input->GetDimensions(dimensions.data());
  input->GetSpacing(spacing.data());
  double bounds[6];
  input->GetBounds(bounds);
  grid_min = Eigen::Vector3d(bounds[0], bounds[2], bounds[4]);
  grid_max = Eigen::Vector3d(bounds[1], bounds[3], bounds[5]);

  auto grid = simple::vtk_grid_3d(dimensions, array);
  auto vector_field = simple::trilinear_interpolator(dimensions, grid_min, grid_max, grid);
  auto system = simple::system_adapter(vector_field);

  auto jacobian_grid = simple::central_differences_jacobian(spacing.template head<2>(), grid);
  auto jacobian_field = simple::trilinear_interpolator(dimensions, grid_min, grid_max, jacobian_grid);
  auto jacobian = simple::system_adapter(jacobian_field);

  Eigen::MatrixXd seeds(3, input_seeds->GetNumberOfPoints());
  {
    double point[3];
    for (vtkIdType i = 0; i < input_seeds->GetNumberOfPoints(); ++i) {    
      input_seeds->GetPoint(i, point);
      seeds(0, i) = point[0];
      seeds(1, i) = point[1];
      seeds(2, i) = point[2];
    }
  }

  output->ShallowCopy(input_seeds);

  if (input_seeds->GetNumberOfLines() <= 0 || !input_seeds->GetLines()) {
    return 0;
  }

  // get input lines
  std::vector<Eigen::MatrixXd> lines_coordinates;
  std::vector<std::vector<vtkIdType>> lines_ids;
  std::vector<vtkIdType> lines_cell_ids;
  {
    auto ids = vtkSmartPointer<vtkIdList>::New();
    for (vtkIdType i = 0; i < input_seeds->GetNumberOfCells(); ++i) {
      if (input_seeds->GetCellType(i) != VTK_POLY_LINE)
        continue;

      input_seeds->GetCellPoints(i, ids);
      std::vector<vtkIdType> point_ids;
      vtkIdType *p = ids->GetPointer(0);
      std::copy(p, p + ids->GetNumberOfIds(), std::back_inserter(point_ids));

      if (point_ids.size() <= 3) {
        // skip too short lines
        continue;
      }

      Eigen::MatrixXd coordinates(3, ids->GetNumberOfIds());
      for (size_t i = 0; i < point_ids.size(); ++i) {
        coordinates.col(i) = seeds.col(point_ids[i]);
      }

      lines_ids.emplace_back(point_ids);
      lines_coordinates.emplace_back(coordinates);
      lines_cell_ids.push_back(i);
    }
  }

  Eigen::VectorXd dht_confidence = Eigen::VectorXd::Zero(lines_coordinates.size());
  std::vector<Eigen::MatrixXd> dht_positions(lines_coordinates.size());
  std::vector<std::vector<Eigen::MatrixXd>> dht_frames(lines_coordinates.size());
  {
    for (size_t i = 0; i < lines_coordinates.size(); ++i) {
      Eigen::MatrixXd position = lines_coordinates[i].topRows(2);
      Eigen::VectorXd time = lines_coordinates[i].row(2);
      const long n = position.cols();

      bool valid = true;
      for (long j = 0; j < n - 1; ++j) {
        if (time[j + 1] <= time[j]) {
          valid = false;
          break;
        }
      }
      if (!valid) {
        continue;
      }

      bool diverged = false;
      int iteration = 0;
      for (; iteration < NumberOfIterations; ++iteration) {
        std::vector<Eigen::MatrixXd> localized(n);

        bool valid = true;
        for (long j = 0; j < n; ++j) {
          localized[j] = jacobian.zero();
          valid = valid && jacobian(position.col(j), time[j], localized[j]);
          if (!valid) {
            break;
          }
        }
        if (!valid) {
          break;
        }

        std::vector<Eigen::MatrixXd> B, RT;
        std::vector<Eigen::VectorXd> sigma;
        hyperbolic_trajectories::continuous_svd(
            localized, time, AbsTol, RelTol, MaxDtFactor, InitialDtFactor, B, sigma, RT);

        Eigen::VectorXd log_multipliers = sigma.back();
        double decimals = (log_multipliers).array().abs().maxCoeff();

        dht_confidence[i] = decimals;

        Eigen::MatrixXd dht_position(2, n);
        std::vector<Eigen::MatrixXd> dht_frame(n);
        diverged = false;

        // largest time interval such that precision suffices
        const double time_block = (time[n - 1] - time[0]) / decimals * SplitPrecisionThreshold;
        const int blocks = std::ceil((time[n - 1] - time[0]) / time_block);

        // first pass: non-overlapping blocks
        for (long start = 0; start < n; ++start) {
          long end = start;
          while ((end < n - 1) && ((time[end] - time[start]) < time_block)) {
            ++end;
          }
          const long block_n = end - start + 1;

          Eigen::MatrixXd block_position(2, block_n);
          Eigen::VectorXd block_time(block_n);
          std::vector<Eigen::MatrixXd> block_localized(block_n);
          std::vector<Eigen::MatrixXd> block_B(block_n), block_RT(block_n);
          std::vector<Eigen::VectorXd> block_sigma(block_n);
          for (long j = 0; j < block_n; ++j) {
            block_position.col(j) = position.col(start + j);
            block_time[j] = time[start + j];
            block_localized[j] = localized[start + j];
          }
          hyperbolic_trajectories::continuous_svd(
              block_localized,
              block_time,
              AbsTol,
              RelTol,
              MaxDtFactor,
              InitialDtFactor,
              block_B,
              block_sigma,
              block_RT);

          Eigen::MatrixXd block_dht_position;
          std::vector<Eigen::MatrixXd> block_dht_frame;
          int iterations = 0;
          double residual = 0.;
          hyperbolic_trajectories::approximate_dht(
              system,
              block_position, block_time, block_localized, block_B, block_sigma, block_RT,
              FixedPointMaxIterations, FixedPointTolerance,
              iterations, residual, block_dht_position, block_dht_frame);

          for (long j = 0; j < block_n; ++j) {
            dht_position.col(start + j) = block_dht_position.col(j);
            dht_frame[start + j] = block_dht_frame[j];
          }
          if (residual > FixedPointConvergenceTolerance) {
            diverged = true;
          }
          start = end;
        }

        if (blocks > 1) {
          // second pass: blocks with half-overlap fix the "gaps"
          long start = 0;
          while ((start < n - 1) && ((time[start] - time[0]) < time_block / 2.)) {
            ++start;
          }
          double time_center = time[0];
          for (; start < n; ++start) {
            long end = start;
            while ((end < n - 1) && ((time[end] - time[start]) < time_block)) {
              ++end;
            }
            // grow final block to the left
            if (end == n - 1) {
              while ((start > 0) && ((time[end] - time[start]) < time_block)) {
                --start;
              }
            }
            time_center += time_block;
            const long block_n = end - start + 1;

            Eigen::MatrixXd block_position(2, block_n);
            Eigen::VectorXd block_time(block_n);
            std::vector<Eigen::MatrixXd> block_localized(block_n);
            std::vector<Eigen::MatrixXd> block_B(block_n), block_RT(block_n);
            std::vector<Eigen::VectorXd> block_sigma(block_n);
            for (long j = 0; j < block_n; ++j) {
              block_position.col(j) = position.col(start + j);
              block_time[j] = time[start + j];
              block_localized[j] = localized[start + j];
            }
            hyperbolic_trajectories::continuous_svd(
                block_localized,
                block_time,
                AbsTol,
                RelTol,
                MaxDtFactor,
                InitialDtFactor,
                block_B,
                block_sigma,
                block_RT);

            Eigen::MatrixXd block_dht_position;
            std::vector<Eigen::MatrixXd> block_dht_frame;
            int iterations = 0;
            double residual = 0.;
            hyperbolic_trajectories::approximate_dht(
                system,
                block_position, block_time, block_localized, block_B, block_sigma, block_RT,
                FixedPointMaxIterations, FixedPointTolerance,
                iterations, residual, block_dht_position, block_dht_frame);
            if (residual > FixedPointConvergenceTolerance) {
              diverged = true;
            }

            for (long j = 0; j < block_n; ++j) {
              double alpha = std::abs(time_center - block_time[j]) / (time_block / 3.);
              double weight = (alpha < M_PI_2) ? std::cos(alpha * M_PI_2) : 0.;
              weight = (weight < 0) ? 0. : weight;

              dht_position.col(start + j) =
                  weight * block_dht_position.col(j) + (1. - weight) * dht_position.col(start + j);
              dht_frame[start + j] = weight * block_dht_frame[j] + (1. - weight) * dht_frame[start + j];
            }
            start = end;
          }
        }

        dht_positions[i] = dht_position;
        dht_frames[i] = dht_frame;

        double residual = (position - dht_position).norm();
        if (residual <= IterationTolerance)
          break;
        std::swap(position, dht_position);
      }
      if (diverged) {
        dht_confidence[i] = -1;
      }
    }
  }

  // output
  auto output_coordinates = vtkSmartPointer<vtkDoubleArray>::New();
  output_coordinates->SetName("Coordinates");
  output_coordinates->SetNumberOfComponents(3);
  output_coordinates->SetNumberOfTuples(input_seeds->GetNumberOfPoints());
  std::fill(
    output_coordinates->GetPointer(0), 
    output_coordinates->GetPointer(0) + 3 * input_seeds->GetNumberOfPoints(), 
    std::numeric_limits<double>::quiet_NaN());

  auto output_frame = vtkSmartPointer<vtkDoubleArray>::New();
  output_frame->SetName("Frame");
  output_frame->SetNumberOfComponents(4);
  output_frame->SetNumberOfTuples(input_seeds->GetNumberOfPoints());
  std::fill(
    output_frame->GetPointer(0), 
    output_frame->GetPointer(0) + 4 * input_seeds->GetNumberOfPoints(), 
    std::numeric_limits<double>::quiet_NaN());

  auto output_confidence = vtkSmartPointer<vtkDoubleArray>::New();
  output_confidence->SetName("HyperbolicityStrength");
  output_confidence->SetNumberOfComponents(1);
  output_confidence->SetNumberOfValues(input_seeds->GetNumberOfCells());
  std::fill(
    output_confidence->GetPointer(0), 
    output_confidence->GetPointer(0) + input_seeds->GetNumberOfCells(), 
    std::numeric_limits<double>::quiet_NaN());

  auto output_points = vtkSmartPointer<vtkPoints>::New();
  output_points->SetDataTypeToDouble();
  output_points->SetNumberOfPoints(input_seeds->GetNumberOfPoints());
  for (vtkIdType i = 0; i < input_seeds->GetNumberOfPoints(); ++i) {
    output_points->SetPoint(i, input_seeds->GetPoint(i));
  }

  for (size_t i = 0; i < lines_ids.size(); ++i) {
    Eigen::VectorXd time = lines_coordinates[i].row(2);
    if ((dht_positions[i].cols() != lines_ids[i].size()) || (dht_frames[i].size() != lines_ids[i].size())) {
      continue;
    }
    output_confidence->SetValue(lines_cell_ids[i], dht_confidence[i]);

    if (dht_confidence[i] < 0 || !dht_positions[i].array().isFinite().all()) {
      continue;
    }

    for (size_t j = 0; j < lines_ids[i].size(); ++j) {
      Eigen::VectorXd position = Eigen::VectorXd::Zero(3);
      position.head(2) = dht_positions[i].col(j);
      position[2] = time[j];
      output_points->SetPoint(lines_ids[i][j],
                              position[0], position[1], position[2]);
      output_coordinates->SetTuple(lines_ids[i][j], position.data());
      output_frame->SetTuple(lines_ids[i][j], dht_frames[i][j].data());
    }
  }

  output->SetPoints(output_points);
  output->GetPointData()->AddArray(output_coordinates);
  output->GetPointData()->AddArray(output_frame);
  output->GetCellData()->AddArray(output_confidence);

  return 1;
}

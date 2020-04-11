#ifndef HYPERBOLIC_TRAJECTORIES_CONTINUOUS_SVD_H
#define HYPERBOLIC_TRAJECTORIES_CONTINUOUS_SVD_H

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/StdVector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/iterator/times_iterator.hpp>
#include <algorithm>
#include <cassert>
#include <vector>

namespace hyperbolic_trajectories {

/**
 * @brief Continuous SVD method
 *
 * Algorithm 3.2 of
 * Dieci, L. and Elia, C., 2008.
 * SVD algorithms to approximate spectra of dynamical systems.
 * Mathematics and Computers in Simulation, 79(4), pp.1235-1254.
 */
void continuous_svd(
    const std::vector<Eigen::MatrixXd>& A,
    const Eigen::VectorXd& time,
    double abs_tol, double rel_tol, double max_dt_factor, double initial_dt_factor,
    std::vector<Eigen::MatrixXd>& B,
    std::vector<Eigen::VectorXd>& sigma,
    std::vector<Eigen::MatrixXd>& RT
) {
  assert(A.size() == time.size());
  const int D = A.front().rows();
  const size_t n = time.size();
  B.resize(n);
  sigma.resize(n);
  RT.resize(n);

  if (A.size() < 3) {
    for (size_t i = 0; i < n; ++i) {
      B[i] = RT[i] = Eigen::MatrixXd::Identity(D, D);
      sigma[i] = Eigen::VectorXd::Zero(D);
    }
    return;
  }

  B[0] = RT[0] = Eigen::MatrixXd::Identity(D, D);
  sigma[0] = Eigen::VectorXd::Zero(D);

  using namespace boost::numeric::odeint;

  const double max_dt = (time.tail(time.size() - 1) - time.head(time.size() - 1)).mean() * max_dt_factor;
  const double initial_dt = max_dt * initial_dt_factor;
  std::vector<double> X0(D * D);
  Eigen::Map<Eigen::MatrixXd>(X0.data(), D, D) = Eigen::MatrixXd::Identity(D, D);

  std::vector<double> times(time.data(), time.data() + time.size());

  auto A_interpolated = [&](double t) -> Eigen::MatrixXd {
    // 1D interpolation of A along times
    auto it = std::lower_bound(times.begin(), times.end(), t);
    size_t i = std::distance(times.begin(), it);
    if (i >= times.size() - 1)
      i = times.size() - 1;
    if (i < 1)
      i = 1;
    double alpha = (t - times[i - 1]) / (times[i] - times[i - 1]);
    return (alpha * A[i] + (1. - alpha) * A[i - 1]);
  };

  auto system = [&](const std::vector<double>& x, std::vector<double>& dxdt, double t) {
    auto X = Eigen::Map<const Eigen::MatrixXd>(x.data(), D, D);
    auto Y = Eigen::Map<Eigen::MatrixXd>(dxdt.data(), D, D);
    Y = A_interpolated(t) * X;
  };

  auto stepper = runge_kutta_dopri5<std::vector<double>>();
  auto integrator = make_dense_output(abs_tol, rel_tol, max_dt, stepper);

  size_t i = 0;
  Eigen::VectorXd singular = Eigen::VectorXd::Constant(D, 1.0);
  auto it_range = make_times_range(integrator, system, X0, times.begin(), times.end(), initial_dt);
  for (auto it = it_range.first; it != it_range.second; ++it)
  {
    if (i > 0) {
      auto X = Eigen::Map<const Eigen::MatrixXd>(it->data(), D, D);
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
      B[i] = svd.matrixU();
      sigma[i] = svd.singularValues().array().log().matrix();
      RT[i] = svd.matrixV().transpose();

      singular = svd.singularValues();
      if (std::abs(singular[0] - singular[D - 1]) > 0.1) {
        break;
      }
    }
    ++i;
  }
  if (i >= n - 1) {
    return;
  }

  size_t initial_svd_index = i;
  times = std::vector<double>(time.data() + i, time.data() + time.size());

  // log(nv_1), ..., log(nv_D-1), nv_D, U, V
  std::vector<double> svd_X0(D + 2 * D * D);
  {
    auto log_nv = Eigen::Map<Eigen::VectorXd>(svd_X0.data(), D);
    auto U = Eigen::Map<Eigen::MatrixXd>(svd_X0.data() + D, D, D);
    auto V = Eigen::Map<Eigen::MatrixXd>(svd_X0.data() + D + D * D, D, D);
    U = B[i];
    V = RT[i].transpose();
    for (int j = 0; j < D - 1; ++j) {
      log_nv[j] = std::log(singular[j + 1] / singular[j]);
    }
    log_nv[D - 1] = std::log(singular[D - 1]);
  }

  auto svd_system = [D, &A_interpolated](const std::vector<double>& x, std::vector<double>& dxdt, double t) {
    auto d_nv = Eigen::Map<Eigen::VectorXd>(dxdt.data(), D);
    auto d_U = Eigen::Map<Eigen::MatrixXd>(dxdt.data() + D, D, D);
    auto d_V = Eigen::Map<Eigen::MatrixXd>(dxdt.data() + D + D * D, D, D);

    auto log_nv = Eigen::Map<const Eigen::VectorXd>(x.data(), D);
    auto U = Eigen::Map<const Eigen::MatrixXd>(x.data() + D, D, D);
    auto V = Eigen::Map<const Eigen::MatrixXd>(x.data() + D + D * D, D, D);

    Eigen::VectorXd nv = log_nv.head(D - 1).array().exp();

    Eigen::MatrixXd C = U.transpose() * A_interpolated(t) * U;

    for (int j = 0; j < D - 1; ++j) {
      d_nv[j] = C(j + 1, j + 1) - C(j, j);
    }
    d_nv[D - 1] = C(D - 1, D - 1);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(D, D);
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(D, D);
    for (int i = 0; i < D; ++i) {
      for (int j = i + 1; j < D; ++j) {
        double nv_ij = 1.;
        for (int k = i; k <= j - 1; ++k) {
          nv_ij *= (k < D - 1) ? nv[k] : log_nv[D - 1];
        }
        H(i, j) = (C(i, j) * nv_ij + C(j, i)) / (nv_ij * nv_ij - 1.);
        K(i, j) = (C(i, j) + C(j, i)) / (nv_ij * nv_ij - 1.) * nv_ij;
        H(j, i) = -H(i, j);
        K(j, i) = -K(i, j);
      }
    }

    d_U = U * H;
    d_V = V * K;
  };


  auto svd_stepper = runge_kutta_dopri5<std::vector<double>>();
  auto svd_integrator = make_dense_output(abs_tol, rel_tol, max_dt, svd_stepper);

  integrate_times(svd_integrator, svd_system, svd_X0, times, initial_dt, [&](const std::vector<double>& state, double t) {
    // modify state of integrator stored in svd_X0
    auto U = Eigen::Map<Eigen::MatrixXd>(svd_X0.data() + D, D, D);
    auto V = Eigen::Map<Eigen::MatrixXd>(svd_X0.data() + D + D * D, D, D);
    auto U_qr = U.householderQr();
    auto V_qr = V.householderQr();
    U = U_qr.householderQ();
    V = V_qr.householderQ();
    // multiply with signs of R's diagonal for consistent orientation
    U *= U_qr.matrixQR().diagonal().array().sign().matrix().asDiagonal();
    V *= V_qr.matrixQR().diagonal().array().sign().matrix().asDiagonal();

    if (i > initial_svd_index) {
      auto log_nv = Eigen::Map<const Eigen::VectorXd>(state.data(), D);
      Eigen::MatrixXd U = Eigen::Map<const Eigen::MatrixXd>(state.data() + D, D, D);
      Eigen::MatrixXd V = Eigen::Map<const Eigen::MatrixXd>(state.data() + D + D * D, D, D);
      auto U_qr = U.householderQr();
      auto V_qr = V.householderQr();
      U = U_qr.householderQ();
      V = V_qr.householderQ();
      // multiply with signs of R's diagonal for consistent orientation
      U *= U_qr.matrixQR().diagonal().array().sign().matrix().asDiagonal();
      V *= V_qr.matrixQR().diagonal().array().sign().matrix().asDiagonal();

      B[i] = U;
      sigma[i] = Eigen::VectorXd(D);
      sigma[i][D - 1] = log_nv[D - 1];
      for (int j = D - 2; j >= 0; --j) {
        sigma[i][j] = sigma[i][j + 1] - log_nv[j];
      }
      RT[i] = V.transpose();
    }
    ++i;
  });

  assert(i == n);
}

}


#endif //HYPERBOLIC_TRAJECTORIES_CONTINUOUS_SVD_H

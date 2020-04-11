#ifndef HYPERBOLIC_TRAJECTORIES_APPROXIMATE_DHT_H
#define HYPERBOLIC_TRAJECTORIES_APPROXIMATE_DHT_H

#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace hyperbolic_trajectories {

template <typename SystemV>
void approximate_dht(
    SystemV velocity,
    const Eigen::MatrixXd& position,
    const Eigen::VectorXd& time,
    const std::vector<Eigen::MatrixXd>& localized,
    const std::vector<Eigen::MatrixXd>& B,
    const std::vector<Eigen::VectorXd>& sigma,
    const std::vector<Eigen::MatrixXd>& RT,
    int fixed_point_max_iterations,
    double fixed_point_tolerance,
    int& iterations,
    double& residual,
    Eigen::MatrixXd& dht_X,
    std::vector<Eigen::MatrixXd>& dht_frame
) {
  assert(time.size() == position.cols());
  const long n = position.cols();
  assert(localized.size() == n);
  assert(B.size() == n);
  assert(sigma.size() == n);
  assert(RT.size() == n);

  dht_X.resizeLike(position);
  dht_frame.resize(n);

  Eigen::MatrixXd position_dot;
  position_dot.resizeLike(position);

  if (n > 1) {
    position_dot.col(0) = (position.col(1) - position.col(0)) / (time[1] - time[0]);
    position_dot.col(n - 1) = (position.col(n - 1) - position.col(n - 2)) / (time[n - 1] - time[n - 2]);
  }
  for (long i = 1; i < n - 1; ++i) {
    position_dot.col(i) = (position.col(i + 1) - position.col(i - 1)) / (time[i + 1] - time[i - 1]);
  }

  std::vector<Eigen::MatrixXd> T(n), T_inv(n);

  T.front() = RT.back();
  T.back() = B.back().transpose();
  T_inv.front() = T.front().transpose();
  T_inv.back() = T.back().transpose();

  Eigen::VectorXd D = sigma.back() / (time[n - 1] - time[0]);

  for (size_t i = 1; i < n - 1; ++i) {
    Eigen::MatrixXd Y_inv = RT[i].transpose() * (-sigma[i]).array().exp().matrix().asDiagonal() * B[i].transpose();
    Eigen::MatrixXd Y = B[i] * sigma[i].array().exp().matrix().asDiagonal() * RT[i];

    T[i] = ((time[i] - time[0]) * D).array().exp().matrix().asDiagonal() * RT.back() * Y_inv;
    T_inv[i] = Y * RT.back().transpose() * (-(time[i] - time[0]) * D).array().exp().matrix().asDiagonal();
  }
  dht_frame = T_inv;

  auto integrate = [&D, &time](const Eigen::MatrixXd& h) -> Eigen::MatrixXd {
    Eigen::MatrixXd w = Eigen::MatrixXd::Zero(h.rows(), h.cols());

    for (long i = 0; i < h.cols(); ++i) {
      for (long c = 0; c < h.rows(); ++c) {
        double d = D[c];
        long from = (d < 0) ? 0 : i;
        long to = (d < 0) ? i : (h.cols() - 1);
        for (long j = from; j < to; ++j) {
          w(c, i) += (time[j + 1] - time[j]) / 2. *
              (std::exp(d * (time[i] - time[j + 1])) * h(c, j + 1)
                  + std::exp(d * (time[i] - time[j])) * h(c, j));
        }
        if (d >= 0) {
          w(c, i) *= -1.;
        }
      }
    }

    return w;
  };

  Eigen::MatrixXd w, h;
  w.resizeLike(position);
  h.resizeLike(position);
  w.setZero();
  int iter = 0;
  for (; iter < fixed_point_max_iterations; ++iter) {
    for (long i = 0; i < n; ++i) {
      Eigen::VectorXd T_inv_w = T_inv[i] * w.col(i);
      Eigen::VectorXd position_transformed = position.col(i) + T_inv_w;
      Eigen::VectorXd v = velocity.zero();
      velocity(position_transformed, time[i], v);
      h.col(i) = T[i] * (v - localized[i] * T_inv_w - position_dot.col(i));
    }

    Eigen::MatrixXd w_next = integrate(h);
    double res = (w - w_next).norm();
    residual = res;
    // Krasnoselskii-Mann iteration with alpha=.5
    w = .5 * w + .5 * w_next;
    if (res <= fixed_point_tolerance) {
      break;
    }
  }
  iterations = iter + 1;

  for (long i = 0; i < n; ++i) {
    dht_X.col(i) = position.col(i) + T_inv[i] * w.col(i);
  }
}

}


#endif //HYPERBOLIC_TRAJECTORIES_APPROXIMATE_DHT_H

#include "kalman_filter.h"
#include "tools.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

//void
//KalmanFilter::init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
//                   MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
//{
//  x_ = x_in;
//  P_ = P_in;
//  F_ = F_in;
//  H_ = H_in;
//  R_ = R_in;
//  Q_ = Q_in;
//}

void
KalmanFilter::predict()
{
  // Predicting the state
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void
KalmanFilter::update(const VectorXd& z)
{
  // Updating the state by using Kalman Filter equations
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  // New estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void
KalmanFilter::updateEKF(const VectorXd& z)
{
  // Update the state by using Extended Kalman Filter equations
  VectorXd y = z - h();
  // theta must be normalized into the (-PI,+PI) range
  if (y[1] > M_PI) y[1] = y[1] - floor(y[1] / M_PI) * M_PI;
  else if (y[1] < -M_PI) y[1] = y[1] - ceil(y[1] / M_PI) * M_PI;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  // New estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

Eigen::VectorXd
KalmanFilter::h() const
{
  Eigen::VectorXd Hxprime(3);
  Hxprime[0] = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
  Hxprime[1] = std::atan2(x_[1], x_[0]); // atan2 handles zero x_[0]
  if (Hxprime[0] > Tools::EPS)
    Hxprime[2] = (x_[0] * x_[2] + x_[1] * x_[3]) / Hxprime[0];
  else
    Hxprime[2] = 0;

  return Hxprime;
}

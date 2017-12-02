#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"

class FusionEKF
{
public:
  FusionEKF();
  virtual ~FusionEKF() = default;

  // Running the whole flow of the Kalman Filter from here
  void processMeasurement(const MeasurementPackage &measurement_pack);

  // Kalman Filter update and prediction math lives in here
  KalmanFilter ekf_;

private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;

  const float Noise_ax;
  const float Noise_ay;
};

#endif /* FusionEKF_H_ */

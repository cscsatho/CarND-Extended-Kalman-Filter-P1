#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF(bool verbose, bool nolidar, bool noradar)
 : is_initialized_(false), verbose_(verbose), nolidar_(nolidar), noradar_(noradar), previous_timestamp_(0), 
   Noise_ax(9.), Noise_ay(9.)
{
  // Measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0     ,
              0     , 0.0225;

  // Measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0     , 0   ,
              0   , 0.0009, 0   ,
              0   , 0     , 0.09;

  // Initializing the rest of FusionEKF
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_      = MatrixXd(3, 4);
//  Hj_ << 0, 0, 0, 0,
//         0, 0, 0, 0,
//         0, 0, 0, 0;

  // The initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // State covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0   , 0,
             0, 1, 0   , 0,
             0, 0, 1000, 0,
             0, 0, 0   , 1000;

  ekf_.x_ = VectorXd(4);
//  ekf_.x_ << 1, 1, 1, 1;

  // The initial transition matrix Q_
  ekf_.Q_ = MatrixXd(4, 4);
//  ekf_.Q_ << 0, 0, 0, 0,
//             0, 0, 0, 0,
//             0, 0, 0, 0,
//             0, 0, 0, 0;

  // No need to initialize ekf_.H_ and ekf_.R_
}

void
FusionEKF::processMeasurement(const MeasurementPackage& measurement_pack)
{
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  const float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1.0E6; // dt - expressed in usec

  if (nolidar_ && measurement_pack.sensor_type_ == MeasurementPackage::LASER) return;
  if (noradar_ && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) return;

  if (!is_initialized_ || fabs(dt) > 60)
  {
    // Initializing the state ekf_.x_ with the first measurement
    // Creating the covariance matrix

    // first measurement
    if (verbose_) cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Convert radar from polar to cartesian coordinates and initialize state
      // [0]=ro, [1]=theta, [2]=ro_dot
      ekf_.x_ << cos(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0],
                 sin(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0],
                 cos(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[2],
                 sin(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[2];
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      // Initialize state
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // Done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  previous_timestamp_ = measurement_pack.timestamp_;

  // Modifying the F matrix so that the time is integrated
  ekf_.F_(0, 2) = ekf_.F_(1, 3) = dt;

  // Setting the process covariance matrix Q
  const float dt2 = dt * dt;
  const float dt3 = dt2 * dt / 2.0;
  const float dt4 = dt3 * dt / 2.0;
  ekf_.Q_ << dt4 * Noise_ax, 0., dt3 * Noise_ax, 0.,
             0., dt4 * Noise_ay, 0., dt3 * Noise_ay,
             0., 0., dt2 * Noise_ax, 0.,
             0., 0., 0., dt2 * Noise_ay;
  ekf_.Q_(2, 0) = ekf_.Q_(0, 2);
  ekf_.Q_(3, 1) = ekf_.Q_(1, 3);

  // Calling the Kalman Filter predict() function
  ekf_.predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Useing the sensor type to perform the update step
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    Tools::CalculateJacobian(ekf_.x_, Hj_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.updateEKF(measurement_pack.raw_measurements_);
  }
  else // if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
  {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.update(measurement_pack.raw_measurements_);
  }

  // print the output
  if (verbose_)
  {
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
  }
}

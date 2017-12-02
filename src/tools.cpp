#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

const float Tools::EPS = 1E-5;
//const int Tools::MaxWnd = 300;

VectorXd
Tools::CalculateRMSE(const vector<VectorXd>& estimations,
                     const vector<VectorXd>& ground_truth)
{
  // calculating the RMSE here - with sliding window of maximal Tools::MaxWnd entries

   VectorXd rmse(4);
   rmse << 0, 0, 0, 0;

   // check the validity of the following inputs:
   // the estimation vector size should not be zero
   if (estimations.empty())
   {
       std::cout << "estimations vector is empty" << std::endl;
       return rmse;
   }

   // the estimation vector size should equal ground truth vector size
   if (estimations.size() != ground_truth.size())
   {
       std::cout << "estimations vector and ground_truth vector sizes do not match" << std::endl;
       return rmse;
   }

   // accumulate squared residuals
   //for (int i = estimations.size() <= MaxWnd ? 0 : estimations.size() - MaxWnd; i < estimations.size(); ++i)
   for (int i = 0; i < estimations.size(); ++i)
   {
//       if (i < 5) std::cout << "("<<i<<")est=" << estimations[i] << std::endl; // FIXME
//       if (i < 5) std::cout << "("<<i<<")gnd=" << ground_truth[i] << std::endl; // FIXME
       VectorXd diff = estimations[i] - ground_truth[i];
//       if (i < 5) std::cout << "("<<i<<")diff=" << diff << std::endl; // FIXME
       diff = diff.array() * diff.array();
       rmse += diff; //diff.array() * diff.array();
   }

   // calculate the mean
   rmse /= estimations.size();

   // calculate the squared root
   rmse = rmse.array().sqrt();

   // return the result
   return rmse;
}

void
Tools::CalculateJacobian(const VectorXd& x_state, MatrixXd& Hj)
{
  // calculating a Jacobian here
  // Hj must be of type MatrixXd(3, 4)
  assert (Hj.rows() == 3 && Hj.cols() == 4);

  // px = x_state[0], py = x_state[1]
  // vx = x_state[2], vy = x_state[3]

  // Pre-compute a set of terms to avoid repeated calculation
  float c1 = x_state[0] * x_state[0] + x_state[1] * x_state[1];
  float c2 = sqrt(c1);
  float c3 = (c1 * c2);

  // Check division by zero
  if (fabs(c1) < EPS)
  {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return;
  }

  // Computing the Jacobian matrix
  Hj << (x_state[0] / c2), (x_state[1] / c2), 0, 0,
       -(x_state[1] / c1), (x_state[0] / c1), 0, 0,
        x_state[1] * (x_state[1] * x_state[2] - x_state[0] * x_state[3]) / c3,
        x_state[0] * (x_state[0] * x_state[3] - x_state[1] * x_state[2]) / c3,
        x_state[0] / c2, x_state[1] / c2;
}

#include <math.h>       /* atan2 */
#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  float phi = atan2(x_(1), x_(0)); // atan2() returns values between -pi and pi.
  float rho_dot;
  // avoid devide by zero
  if (fabs(rho) < 0.0001) {
    rho_dot = 5; // initial speed
  } else {
    rho_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
  }

  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  
  VectorXd y = z - z_pred;
  
  cout << "Using updateEKF funciton" << endl;
  cout << "+++++++++++++++++++++++++" << endl;
  cout << "Angle before regulization is : " << y(1) << endl;
  cout << "+++++++++++++++++++++++++" << endl;
  
  // adjust y data so that they are between -pi and pi
  double pi = 3.14;
  double ang = y(1);
  // data type need to match?
  if (ang >= pi)
  {
	  ang = ang - (2*pi);
  }
  else if (ang <= -3.14)
  {
	  ang = ang + (2*pi);
  }
  y(1) = ang;
  cout << "------------------------" << endl;
  cout << "Angle for the update is : " << y(1) << endl;
  cout << "-------------------------" << endl;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

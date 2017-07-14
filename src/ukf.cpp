#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.281;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.265;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //set initially
  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //set spreading parameter
  lambda_ = 3 - n_aug_;

  //set NIS laser
  NIS_laser_ = 0;

  //set NIS radar
  NIS_radar_ = 0;

  // initial matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // initial vector for weights_
  weights_ = VectorXd::Zero(2 * n_aug_ + 1);

  // Update covariance matrix
  P_ << 0.504, 0,    0,     0,     0,
        0,     0.45, 0,     0,     0,
        0,     0,    3.145, 0,     0,
        0,     0,    0,     0.034, 0,
        0,     0,    0,     0,     4.208;
}

UKF::~UKF()
{
  // Do nothing
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    // first measurement
    cout << "UKF Initialization: " << endl;

    //set vector for weights_
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++) //2n+1 weights_
    {
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
    }

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float px = ro * cos(theta);
      float py = ro * sin(theta);
      x_(0) = px;
      x_(1) = py;
      cout << "RADAR x: " << x_[0] << " " << x_[1] << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      */
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      x_(0) = px;
      x_(1) = py;
      cout << "LASER x: " << x_[0] << " " << x_[1] << endl;
    }

    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  // This fix does not work with dataset obj_pose-laser-radar-synthetic-input.txt
  // because in dataset have delta_t is always 0.05
  while (delta_t > 0.1)
  {
    const double dt = 0.05;
    Prediction(dt);
    delta_t -= dt;
  }
  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    if (use_radar_== true)
    {
      // Radar updates
      UpdateRadar(meas_package);
    }
  }
  else
  {
    if (use_laser_== true)
    {
      // Laser updates
      UpdateLidar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /*****************************************************************************
   * Generating Sigma Points
   ****************************************************************************/

  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  /*****************************************************************************
   * Sigma Point Prediction
   ****************************************************************************/

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v/yawd * (sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw+yawd*delta_t));
    }
    else
    {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  /*****************************************************************************
   * Predicted Mean and Covariance
   ****************************************************************************/

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //iterate over sigma points
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //iterate over sigma points
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = angle_normalize(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //set measurement dimension, lidar can measure px, py
  int n_z = 2;

  /*****************************************************************************
   * Predict Lidar Measurements
   ****************************************************************************/

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //2n+1 simga points
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  //create matrix for cross correlation Tc
  //use for "Update State" step
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //2n+1 simga points
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = angle_normalize(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd::Zero(n_z,n_z);
  R << std_a_*std_a_, 0,
       0, std_a_*std_a_;
  S = S + R;

  /*****************************************************************************
   * Update State
   ****************************************************************************/

  //create vector for incoming lidar measurement
  VectorXd z = VectorXd::Zero(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //calculate the lidar NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  /*****************************************************************************
   * Predict Radar Measurements
   ****************************************************************************/

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //2n+1 simga points
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    if (p_x == 0 && p_y == 0)
    {
      Zsig(0,i) = 0;
      Zsig(1,i) = 0;
      Zsig(2,i) = 0;
    }
    else
    {
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
      Zsig(1,i) = atan2(p_y,p_x);                                 //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  //create matrix for cross correlation Tc
  //use for "Update State" step
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) //2n+1 simga points
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = angle_normalize(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = angle_normalize(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd::Zero(n_z,n_z);
  R << std_radr_*std_radr_,     0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S = S + R;

  /*****************************************************************************
   * Update State
   ****************************************************************************/

  //create vector for incoming radar measurement
  VectorXd z = VectorXd::Zero(n_z);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1],
       meas_package.raw_measurements_[2];

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = angle_normalize(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //calculate the radar NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

double UKF::angle_normalize(double angle)
{
  angle = fmod(angle, (2.* M_PI));
  if (angle > M_PI)
  {
    angle -= 2.* M_PI;
  }
  else if (angle < -M_PI)
  {
    angle += 2.* M_PI;
  }

  return angle;
}

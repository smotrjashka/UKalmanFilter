#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */


    is_initialized_ = false;

    Xsig_pred_ = MatrixXd(5,15);

    n_x_ = 5;

    n_aug_ = 7;

    lambda_ = 3 - n_aug_;

    weights_ = VectorXd(15);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */


    if (!is_initialized_)
    {
        if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            use_radar_ = false;
            x_[0] = meas_package.raw_measurements_[0];
            x_[1] = meas_package.raw_measurements_[1];
            x_[2] = 0.2;
            x_[3] = 0;
            x_[4] = 0;
            P_<< std_laspx_*std_laspx_, 0, 0, 0, 0,
                    0, std_laspy_*std_laspy_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            use_laser_ = false;
            double rho = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            double rhodot = meas_package.raw_measurements_[2];
            x_[0] = rho*cos(phi);
            x_[1] = rho*sin(phi);
            x_[2] = rhodot;
            x_[3] = phi;
            x_[4] = 0;
            P_<< std::pow(std_radr_, 2), 0, 0, 0, 0,
                    0, std::pow(std_radr_, 2), 0, 0, 0,
                    0, 0, std::pow(std_radrd_, 2), 0, 0,
                    0, 0, 0, 0.09, 0,
                    0, 0, 0, 0, 0.09;
        }
        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }
    double delta_t = (meas_package.timestamp_-time_us_)/1e6;
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        UpdateRadar(meas_package);
    }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

    VectorXd x_aug = VectorXd(7);
    MatrixXd Xsig_aug = MatrixXd(n_aug_,2*n_aug_+1);
    MatrixXd P_aug = MatrixXd(7,7);

    x_aug.head(5) = x_;
    x_aug[5] = 0;
    x_aug[6] = 0;

    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std::pow(std_a_, 2);
    P_aug(6,6) = std::pow(std_yawdd_, 2);

    MatrixXd L = P_aug.llt().matrixL();
    Xsig_aug.fill(0.0);
    int last = 2*n_aug_+1;
    for (int i=0; i<last; ++i)
    {
        Xsig_aug.col(i) = x_aug;
    }
    Xsig_aug.block<7,7>(0,1) += sqrt(lambda_+n_aug_)*L;
    Xsig_aug.block<7,7>(0,n_aug_+1) -= sqrt(lambda_+n_aug_)*L;

    for (int i=0; i<2*n_aug_+1; ++i)
    {
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double  v  = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yaw_dot = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yaw_dotdot = Xsig_aug(6,i);

        double px_pred, py_pred;

        if (fabs(yaw_dot) > 0.001)
        {
            px_pred = p_x + v/yaw_dot * (sin(yaw+yaw_dot*delta_t) - sin(yaw));
            py_pred = p_y - v/yaw_dot * (cos(yaw+yaw_dot*delta_t) - cos(yaw));
        }
        else
        {
            px_pred = p_x + v*cos(yaw)*delta_t;
            py_pred = p_y + v*sin(yaw)*delta_t;
        }
        double v_pred = v;
        double yaw_pred = yaw + yaw_dot*delta_t;
        double yaw_dot_pred = yaw_dot;

        px_pred += 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_pred += 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_pred += nu_a*delta_t;
        yaw_pred += 0.5*nu_yaw_dotdot*delta_t*delta_t;
        yaw_dot_pred += nu_yaw_dotdot*delta_t;

        Xsig_pred_(0,i) = px_pred;
        Xsig_pred_(1,i) = py_pred;
        Xsig_pred_(2,i) = v_pred;
        Xsig_pred_(3,i) = yaw_pred;
        Xsig_pred_(4,i) = yaw_dot_pred;
    }

    weights_(0) = lambda_/(lambda_+n_aug_);
    for (int i=1; i<2*n_aug_+1; ++i)
    {
        weights_(i) = 1/(2*lambda_+2*n_aug_);
    }
    x_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; ++i)
    {
        x_ += weights_(i)*Xsig_pred_.col(i);
    }
    P_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; ++i)
    {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        while (x_diff(3)>M_PI) {
            x_diff(3)-=2.*M_PI;
        }
        while (x_diff(3)<-M_PI) {
            x_diff(3)+=2.*M_PI;
        }

        P_ = P_ + weights_(i)*(std::pow(x_diff, 2)).transpose();
    }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

    VectorXd z = meas_package.raw_measurements_;
    int n_z = 2;
    int last = 2*n_aug_+1;
    MatrixXd Zsig = MatrixXd(n_z, last);
    MatrixXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);
    MatrixXd Tc = MatrixXd(n_x_,n_z);

    Zsig.fill(0.0);
    for (int i=0; i<last; i++)
    {
        Zsig(0,i) = Xsig_pred_(0,i);
        Zsig(1,i) = Xsig_pred_(1,i);
    }
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; ++i)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    S.fill(0.0);
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {

        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

    }

    MatrixXd R = MatrixXd(n_z,n_z);
    R <<  std::pow(std_laspx_, 2), 0,
            0, std::pow(std_laspy_, 2);
    S = S + R;

    MatrixXd K = Tc * S.inverse();
    VectorXd z_diff = z - z_pred;

    x_ = x_ + K*z_diff;
    P_ = P_ - K*S*K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

    VectorXd z = meas_package.raw_measurements_;
    int n_z = 3;
    int last = 2 * n_aug_ + 1;
    MatrixXd Zsig = MatrixXd(n_z, last);
    MatrixXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);
    MatrixXd Tc = MatrixXd(n_x_,n_z);
    Zsig.fill(0.0);
    for (int i=0; i<last; i++)
    {
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
        Zsig(1,i) = atan2(p_y,p_x);
        Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
    }
    z_pred.fill(0.0);
    for (int i=0; i < last; ++i){
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    S.fill(0.0);
    Tc.fill(0.0);
    for (int i = 0; i < last; ++i) {

        VectorXd z_diff = Zsig.col(i) - z_pred;

        while (z_diff(1)> M_PI) {
            z_diff(1) -= 2. * M_PI;
        }
        while (z_diff(1)<-M_PI) {
            z_diff(1) += 2. * M_PI;
        }

        S = S + weights_(i) * z_diff * z_diff.transpose();
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        while (x_diff(3)> M_PI){
            x_diff(3)-=2. * M_PI;
        }
        while (x_diff(3)<-M_PI) {
            x_diff(3) += 2. * M_PI;
        }

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd R = MatrixXd(n_z,n_z);
    R <<  std::pow(std_radr_, 2), 0, 0,
            0, std::pow(std_radphi_, 2), 0,
            0, 0, std::pow(std_radrd_, 2);
    S = S + R;

    MatrixXd K = Tc * S.inverse();
    VectorXd z_diff = z - z_pred;

    while(z_diff(1)>M_PI) {
        z_diff(1) -= 2. * M_PI;
    }
    while(z_diff(1)<-M_PI) {
        z_diff(1) += 2. * M_PI;
    }

    x_ = x_ + K*z_diff;
    P_ = P_ - K*S*K.transpose();

}
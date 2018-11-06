#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double EPS = 0.0001;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
	std_a_ = 1;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.8;
	
	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
	
	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
 // initially set to false, set to true in first call of ProcessMeasurement
	is_initialized_ = false;

	// time when the state is true, in us
	time_us_ = 0.0;

	// state dimension
	n_x_ = 5;

	// Augmented state dimension
	n_aug_ = 7;

	//create augmented mean vector
	x_aug_ = VectorXd(7);

	//create augmented state covariance
	P_aug_ = MatrixXd(7, 7);

	//create sigma point matrix
	Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	// Sigma point spreading parameter
	lambda_ = 3 - n_x_;

	// predicted sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	//create vector for weights
	weights_ = VectorXd(2 * n_aug_ + 1);
	// double weight_0 = lambda_/(lambda_+n_aug_);
	weights_(0) = lambda_/(lambda_+n_aug_);
	for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
		// double weight = 0.5/(n_aug_+lambda_);
		weights_(i) = 0.5/(n_aug_+lambda_);
	}

	H_laser_ = MatrixXd(2, 5);
	H_laser_ << 1, 0, 0, 0, 0,
			  	0, 1, 0, 0, 0;
	


	// the current NIS for radar
	NIS_radar_ = 0.0;

	// the current NIS for laser
	NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/
	if (!is_initialized_) {
		
		P_ << 	1,    0,    0,    0,   0,
				0,    1,    0,    0,    0,
				0,    0,    1,    0,    0,
				0,    0,    0,    1,    0,
				0,    0,    0,    0,    1;
		time_us_ = meas_package.timestamp_;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			float ro = meas_package.raw_measurements_(0);
			float phi = meas_package.raw_measurements_(1);
			float ro_dot = meas_package.raw_measurements_(2);

			x_ << ro*cos(phi), ro*sin(phi), 1, 1, 0.1;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 1, 1, 0.1;
		}

		is_initialized_ = true;
	}

	float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;

    Prediction(dt);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar(meas_package);
	}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/

	//create sigma point matrix
	MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

	// calculate square root of P
	// P.llt().matrixL() produces the lower triangular matrix L of the matrix P such that P = L*L^.
	MatrixXd A = P_.llt().matrixL();

	//calculate sigma points ...
	//set sigma points as columns of matrix Xsig
	Xsig.col(0)  = x_;
	for (int i = 0; i < n_x_; i++) {
		Xsig.col(i+1)      = x_ + sqrt(lambda_+n_x_) * A.col(i);
		Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
	}

	/****************************
	 Create augmented variable
	 ***************************/
	//create augmented mean state
	// x.head(n) = y, where n is the number of elements from first element, and y is an input vector of that size.
	x_aug_.head(5) = x_;
	x_aug_(5) = 0;
	x_aug_(6) = 0;

	//create augmented covariance matrix
	// P_aug = [P, 0,
	//          0, Q]
	// Q = [std_a^2,           0,
	//            0, std_yawdd^2]

	P_aug_.topLeftCorner(5,5) = P_;
	P_aug_(5,5) = std_a_*std_a_;
	P_aug_(6,6) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug_.llt().matrixL();
	
	//create augmented sigma points
	Xsig_aug_.col(0) = x_aug_;
	for(int i =0; i<n_aug_; ++i) {
		Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
		Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
	}

	
	//predict sigma points
	for(size_t i = 0; i < 2*n_aug_+1; i++) {
		const double p_x      = Xsig_aug_(0, i);
		const double p_y      = Xsig_aug_(1, i);
		const double v        = Xsig_aug_(2, i);
		const double yaw      = Xsig_aug_(3, i);
		const double yawd     = Xsig_aug_(4, i);
		const double nu_a     = Xsig_aug_(5, i);
		const double nu_yawdd = Xsig_aug_(6, i);

		double px_p, py_p;

		// avoid division by zero
		if (fabs(yawd) > EPS) {
			px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw)- cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}
		
		double v_p = v;
		double yaw_p = yaw+ yawd*delta_t;
		double yawd_p = yawd;

		// add noise
		px_p = px_p + 0.5 * nu_a * delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5 * nu_a * delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a *delta_t;

		yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd * delta_t;

		// write predicted sigma points into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;

	}
	// set weights
	

	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	
	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) {
			x_diff(3) -= 2.*M_PI;
			// std::cout << "- x_diff(3): " << x_diff(3) << "\n";
		} 
		while (x_diff(3)< -M_PI) {
			x_diff(3) += 2.*M_PI;
			// std::cout << "+ x_diff(3): " << x_diff(3) << "\n";
		}

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = meas_package.raw_measurements_ - z_pred;
	MatrixXd Ht = H_laser_.transpose();

	MatrixXd R = MatrixXd(2, 2);
	R << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;
	MatrixXd S = H_laser_ * P_ * Ht + R;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;


/*
	//extract measurement as VectorXd
	VectorXd z = meas_package.raw_measurements_;

	//set measurement dimension, lidar can only measure p_x and p_y
	int n_z = 2;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);

		// measurement model
		Zsig(0, i) = p_x;
		Zsig(1, i) = p_y;
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	//cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

	
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * z_diff * z_diff.transpose();
		
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//calculate cross correlation matrix Tc
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S = S + R;
*/
	

	/*****************************************************************************
	 *  UKF Update for Lidar
	 ****************************************************************************/
/*	
	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
*/

	//calculate NIS
	NIS_laser_ = y.transpose() * Si * y;
	cout << " Lidar measurement : " << NIS_laser_ << "\n";
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/

	//extract measurement as VectorXd
	VectorXd z = meas_package.raw_measurements_;

	//set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		const double p_x = Xsig_pred_(0, i);
		const double p_y = Xsig_pred_(1, i);
		const double v   = Xsig_pred_(2, i);
		const double yaw = Xsig_pred_(3, i);

		const double v1 = cos(yaw)*v;
		const double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                       //r
		if(p_y==0 && p_x==0) {
			Zsig(1, i) = 0;
		} else {
			Zsig(1, i) = atan2(p_y, p_x);                               //phi
		}
		
		Zsig(2, i) = (p_x*v1 + p_y*v2) / std::max(EPS, sqrt(p_x*p_x + p_y*p_y));   //r_dot
  	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R <<  std_radr_*std_radr_,                       0,                     0,
							0, std_radphi_*std_radphi_,                     0,
							0,                       0, std_radrd_*std_radrd_;
	S = S + R;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	/*****************************************************************************
	 *  UKF Update for Radar
	 ****************************************************************************/
	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	//calculate NIS
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
	cout << " Radar measurement : " << NIS_laser_ << "\n";
}

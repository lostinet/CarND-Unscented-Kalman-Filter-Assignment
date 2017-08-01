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
  // this need to be adjusted in order to get Kalman filter working.
  // 30 means around -60m/s^2 to 60m/s^2 or -60rad/s^2 to 60m/s^2 in 95% of the time.
  // an perspective upper limit for sport car is estimated at 12m/s^2, which equals to 6?
  // this shall be tuned if the UKF works.
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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
    // initially set to false, make it "true" after calling ProcessMeasurement
    is_initialized_ = false;
    
    // initialize dimension of nx for UKF
    int n_x_ = 5;
    
    // initialize augmented dimension n_aug_ for UKF
    int n_aug_ = 7;
    
    //initialize lambda
    lambda = 3 - n_x_;
    
    //initial weights vectors
    weights = VectorXd(2 * n_aug_ + 1);
    
    // initialize Xsig_pred
    Xsig_pred = MatrixXd(n_x_,2 * n_aug_ + 1);
    
    // time stamp for initialization.
    init_stamp_= 0.0;
    
    //initialize radar NIS term
    NIS_radar_ = 0.0;
    
    // initialize laser NIS term
    NIS_laser_ = 0.0;
    
  }

UKF::~UKF() {}

/**************************************************************************************/

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

/*initialization*/
    

        if (!is_initialized_){
        
        //first data
        x_ << 1,1,1,1,1;
        
        P_ << 1,0,0,0,0,
              0,1,0,0,0,
              0,0,1,0,0,
              0,0,0,1,0,
              0,0,0,0,1;
            
        //init timestamp
        init_stamp_ = meas_package.timestamp_;
        
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
            // convert radar from polar to cartesian coordinates and initialize state
            float rho = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
          float rho_dot = meas_package.raw_measurements_(2); // no use
            
            // get measurement from radar
            x_(0) = rho * cos(phi);
            x_(1) = rho * sin(phi);
            //converthe state to ukf state model
        }
        if(meas_package.sensor_type_ == MeasurementPackage::LASER){
            //get measurement from laser
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }
            //capture the timestap for next iteration
        is_initialized_ = true;
        return;
            
        }
    
    
 /* Prediction */
    if(meas_package.timestamp_ != init_stamp_){
    // compute the time elapsed  between current and previous measurements
    float dt = (meas_package.timestamp_ - init_stamp_)/1000000.0;
    init_stamp_ = meas_package.timestamp_;
    
    // call Predict function
    Prediction(dt);
    }
    
/* Call for Update */
        
     //call update function for laser data
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
        UpdateLidar(meas_package);
    }
    
    // call update function for radar data
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
        UpdateRadar(meas_package);
    }
    

    
}




void UKF::Prediction(double delta_t) {


        
/* generate sigma points */
    
    // create sigma point matrix
    MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
    
    // calculate square root P
    MatrixXd A = P_.llt().matrixL();
    
    // Lambda setup for non-augmented sigma points
    lambda = 3 - n_x_;
    
    // set first column of sigma pint matrix
    Xsig.col(0) = x_;
    
    // set remaining sigma points
    for(int i=0;i<n_x_;i++){
        
        Xsig.col(i+1)      =  x_ + sqrt(lambda + n_x_) * A.col(i);
        Xsig.col(i+1+n_x_) =  x_ - sqrt(lambda + n_x_) * A.col(i);
        
    }

    
/*Augmented sigma points*/
        
    // create augmented state vector
    VectorXd x_aug_ = VectorXd(n_aug_);
    
    // create augmented state covariance
    MatrixXd P_aug_ = MatrixXd(n_aug_,n_aug_);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_,2 * n_aug_ + 1);
    
    // set lambda for augmented sigma points
    lambda = 3 - n_aug_;
    
    // create augmented mean state vector
    x_aug_.head(5) = x_;
    x_aug_(5) = 0;
    x_aug_(6) = 0;
    
    //create augmented coariance matrix
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(5,5) = P_;
    P_aug_(5,5) = std_a_ * std_a_;
    P_aug_(6,6) = std_yawdd_ * std_yawdd_;
    
    //create square root matrix
    MatrixXd L = P_aug_.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug.col(0) = x_aug_;
    
    for (int i=0;i<n_aug_;i++)
    {
        Xsig_aug.col(i+1)        = x_aug_ + sqrt(lambda + n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug_ - sqrt(lambda + n_aug_) * L.col(i);
        
    }
    

    
/*predict sigma points*/
    
    for (int i=0;i<2*n_aug_+1;i++){
        
        double p_x =      Xsig_aug(0,i);
        double p_y =      Xsig_aug(1,i);
        double v   =      Xsig_aug(2,i);
        double yaw =      Xsig_aug(3,i);
        double yawd =     Xsig_aug(4,i);
        double nu_a =     Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p,py_p;
        
        //avoid division by zero. calculate the predicted value of px and py
        if(fabs(yawd)>0.001)
        {
            px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        }
        else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }
        // predicted velocity, yaw and yaw rate
        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;
        
        //add noise section
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a * delta_t;
        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;
        

        // assign sigma point prediction
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
        
    }
    
    
/*Predict Mean and Covariance Matrix*/
        
        //set weights
    
        double weight_0 = lambda / (lambda + n_aug_);
        weights(0) = weight_0;
    
        for (int i=1;i<2 * n_aug_ + 1;i++){
            double weight = 0.5 / ( n_aug_ + lambda);
            weights(i) = weight;
        }
        
        //predict state mean
        x_.fill(0.0);
        for(int i=0;i<2*n_aug_+1;i++){
            //iterate over sigma points to calculate mean
            x_ = x_ + weights(i) * Xsig_pred.col(i);
        }
        
        //Predicted state covariance matrix
        P_.fill(0.0);
        for(int i=0;i<2*n_aug_+1;i++){
            //state diff
            VectorXd x_diff = Xsig_pred.col(i) - x_;
            
            //angle normalization
            while(x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
            while(x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
            
            // calculate covariance
            P_ = P_ + weights(i) * x_diff * x_diff.transpose();
            
        }

    
}


    
    /**************************************************************************************/
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
    // get measurement as z
    VectorXd z_ = meas_package.raw_measurements_;
    
    // set z_ dimension
    int n_z_ = 2;
    
    // create matrix for sigma point in measurement space.
    MatrixXd Zsig = MatrixXd(n_z_,2 * n_aug_ + 1);
    
    // transform sigma points into measurement space
    for(int i=0;i<2 * n_aug_ +1;i++){
        
        
        Zsig(0,i) = Xsig_pred(0,i);
        Zsig(1,i) = Xsig_pred(1,i);
        
    }
    
    // mean predicted measurement, no R considered
    VectorXd z_pred = VectorXd(n_z_);
    z_pred.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_,n_z_);
    S.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
        
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //calculate covariance S
        S = S + weights(i) * z_diff * z_diff.transpose();
    
    }
    
    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z_,n_z_);
    
    // define R
    R << std_laspx_*std_laspx_,0,
        0,std_laspy_*std_laspy_;
    
    // update S with R
    S = S + R;
    
/* Update State */
    
    //create matrix for cress correlation Tc
    MatrixXd Tc = MatrixXd(n_x_,n_z_);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        
        //state diff
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        
        // Calculate Tc
        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
        
    }
    
    //Kalman gain K
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = z_ - z_pred;
    
    //Calculate NIS
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
    
    // update state mean and covariance matrix
    
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
    

}

/**************************************************************************************/

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    VectorXd z_ = meas_package.raw_measurements_;
    
    // set z_ dimension
    int n_z_ = 3;
    
    // create matrix for sigma point in measurement space.
    MatrixXd Zsig = MatrixXd(n_z_,2 * n_aug_ + 1);
    
    // transform sigma points into measurement space
    for(int i=0;i<2 * n_aug_ +1;i++){
        
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v   = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);
        
        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;
        
        //measurement model
        Zsig(0,i) = sqrt(p_x * p_x + p_y * p_y);//r
        Zsig(1,i) = atan2(p_y,p_x);//phi
        Zsig(2,i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);//r_dot
        
    }
    
    // mean predicted measurement, no R considered
    VectorXd z_pred = VectorXd(n_z_);
    z_pred.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_,n_z_);
    S.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
    
    
    //angle normalization
    while(z_diff(1) >  M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;
    
    //calculate covariance S
    S = S + weights(i) * z_diff * z_diff.transpose();
        
    }
    
    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z_,n_z_);
    
    // define R
    R << std_radr_*std_radr_,0,0,
    0,std_radphi_* std_radphi_,0,
    0,0,std_radrd_*std_radrd_;
    
    //update S with R
    S = S + R;
        
    
    
    /* Update State */
    
    //create matrix for cress correlation Tc
    MatrixXd Tc = MatrixXd(n_x_,n_z_);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        while(z_diff(1)>  M_PI) z_diff(1) -= 2.*M_PI;
        while(z_diff(1)< -M_PI) z_diff(1) += 2.*M_PI;
        
        
        //state diff
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        
        //angle normalization
        while(z_diff(1)>  M_PI) z_diff(1) -= 2.*M_PI;
        while(z_diff(1)< -M_PI) z_diff(1) += 2.*M_PI;
        
        //calculate Tc
        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
        
    }
    
    //Kalman gain K
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = z_ - z_pred;
    
    //angle normalization
    while(z_diff(1)>  M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1)< -M_PI) z_diff(1) += 2.*M_PI;
    
    //Calculate NIS
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
    
    // update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    
    
}

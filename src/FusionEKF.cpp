#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  
  //measurement covariance matrix - radar
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // set the state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * We'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      cout << "Kalman Filter Initialization: Radar " << endl;
      
    // Set the state with the initial location and velocity
    double rho = measurement_pack.raw_measurements_[0]; // range
    double phi = measurement_pack.raw_measurements_[1]; // bearing
    
    // Coordinates convertion from polar to cartesian
    double px = rho * cos(phi);
    if ( px < 0.0001 ) {
        px = 0.0001;
    }
    double py = rho * sin(phi);
    if ( py < 0.0001 ) {
       py = 0.0001;
    }

    double vx = 0;
    double vy = 0;
     
    // Initialise the variables
    ekf_.x_ << px, py, vx, vy;
   
    previous_timestamp_ = measurement_pack.timestamp_;      

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    
    // Initialize state.
    cout << "Kalman Filter Initialization: Laser " << endl;

    // set the state with the initial location and zero velocity
    ekf_.x_ << measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1], 
              0, 
              0;
   
    previous_timestamp_ = measurement_pack.timestamp_;
      
    }
    
    
    // Done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  // Compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
    
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  // Create the state transition matrix F_
  ekf_.F_ = MatrixXd(4,4);
  // Modify the F matrix so that the time is integrated
  ekf_.F_ <<  1,0,dt,0,
              0,1,0,dt,
              0,0,1,0,
              0,0,0,1;
  
 
  //Set the process covariance matrix Q
  float noise_ax = 9;
  float noise_ay = 9;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  //Call the Kalman Filter predict() function
  ekf_.Predict();
    

  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar update equations

    // Update/Initialise required H_in (Measurement matrix) and R_in (Measurement covariance matrix)
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    
    // Call the Kalman Filter update() function for radar
    //      with the most recent raw measurements_
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser update equations
    
    // Update/Initialise required H_in (Measurement matrix) and R_in (Measurement covariance matrix)
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    
    // Call the Kalman Filter update() function for laser
    //      with the most recent raw measurements_
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}

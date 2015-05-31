#include "unscented_KF.h"

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file unscented_KF.cpp
 *  @brief Implements a UKF filter with quaternions to filter angular positions
 *  @version 1.0
 *  @date 12-11-2013
 */

UKF::UKF(const MatrixXf &initial_state,
         const MatrixXf &initial_cov,
         const MatrixXf &process_noise,
         const MatrixXf &measurement_noise,
         float kappa)  {

  m_time_step = 1.f;

  b_first_time = true;
  m_particles.reset();
  m_q_particles.reset();
  _propagateFunc    = NULL;
  _measurementFunc  = NULL;

  if ((initial_state.cols () != 1)
      || (initial_cov.rows () != initial_state.rows ())
      || (initial_cov.rows () != initial_cov.cols ())
      || (process_noise.rows () != process_noise.cols ())
      || (measurement_noise.rows () != measurement_noise.cols ())){

    THROW_ERR("UKF : Matrices don't have compatible shapes");
  }

  // Take extended state vector into account
  int dim_ext = initial_state.rows ()
      + process_noise.rows ()
      + measurement_noise.rows ();

  m_dim = initial_state.rows ();

  // Initialize matrices --> extended state taking noises into account
  k_state_pre.setZero (dim_ext, 1);
  k_state_pre.block (0,0, initial_state.rows (), 1) = initial_state;

  k_cov_pre.setZero (dim_ext, dim_ext);
  k_cov_pre.block (0,0, m_dim, m_dim) = initial_cov;
  k_cov_pre.block (m_dim,m_dim, process_noise.cols (), process_noise.cols ()) = process_noise;
  k_cov_pre.block (m_dim + process_noise.cols (),
                   m_dim + process_noise.cols (),
                   measurement_noise.cols (),
                   measurement_noise.cols ()) = measurement_noise;

  // Set _post values
  k_state_post  = k_state_pre;
  k_cov_post    = k_cov_pre;

  k_expected_cov.setIdentity (dim_ext, dim_ext);

  k_process_noise     = process_noise;
  k_measurement_noise = measurement_noise;

  k_cov_cross_pred_meas.setZero (dim_ext, dim_ext);
  k_cov_cross_pred_ref.setZero (dim_ext, dim_ext);

  k_gain.setZero (dim_ext, dim_ext);
  k_innovation.setZero (dim_ext, 1);

  m_kappa = kappa;

  // Standard "vector space" UKF, no quaternions
  b_use_quaternions = false;
}


UKF::UKF(const MatrixXf &initial_state,
         const MatrixXf &initial_q_state,
         const MatrixXf &initial_cov,
         const MatrixXf &initial_q_cov,
         const MatrixXf &process_noise,
         const MatrixXf &process_q_noise,
         const MatrixXf &measurement_noise,
         const MatrixXf &measurement_q_noise,
         float kappa,
         float kappa_q)  {

  m_time_step = 1.f;

  b_first_time = true;
  m_particles.reset();
  m_q_particles.reset();

  _propagateFunc = NULL;
  _measurementFunc = NULL;
  _propagateQFunc = NULL;
  _measurementQFunc = NULL;

  if ((initial_state.cols () != 1)
      || (initial_cov.rows () != initial_state.rows ())
      || (initial_cov.rows () != initial_cov.cols ())
      || (process_noise.rows () != process_noise.cols ())
      || (measurement_noise.rows () != measurement_noise.cols ())
      || (initial_q_cov.rows () != initial_q_state.rows ())
      || (initial_q_cov.rows () != initial_q_cov.cols ())
      || (process_q_noise.rows () != process_q_noise.cols ())
      || (measurement_q_noise.rows () != measurement_q_noise.cols ())){
    THROW_ERR("UKF : Matrices don't have compatible shapes");
  }

  // Take extended state vector into account
  int dim_ext = initial_state.rows ()
      + process_noise.rows ()
      + measurement_noise.rows ();

  m_dim = initial_state.rows ();

  //Initialize matrices
  //  --> extended state taking noises into account

  k_state_pre.setZero (dim_ext, 1);
  k_state_pre.block (0,0, initial_state.rows (), 1) = initial_state;

  k_cov_pre.setZero (dim_ext, dim_ext);
  k_cov_pre.block (0,0, m_dim, m_dim) = initial_cov;
  k_cov_pre.block (m_dim,m_dim, process_noise.cols (), process_noise.cols ()) = process_noise;
  k_cov_pre.block (m_dim + process_noise.cols (),
                   m_dim + process_noise.cols (),
                   measurement_noise.cols (),
                   measurement_noise.cols ()) = measurement_noise;

  // Initialize quaternions matrices
  // TODO : check initial dimensions
  m_dim_q = initial_q_state.rows ();

  k_state_q_pre = initial_q_state;

  k_cov_q_pre = initial_q_cov;

  // Set _post values
  k_state_post = k_state_pre;
  k_cov_post = k_cov_pre;

  k_expected_cov.setIdentity (dim_ext, dim_ext);

  k_process_noise = process_noise;
  k_measurement_noise = measurement_noise;

  k_cov_cross_pred_meas.setZero (dim_ext, dim_ext);
  k_cov_cross_pred_ref.setZero (dim_ext, dim_ext);

  k_gain.setZero (dim_ext, dim_ext);
  k_innovation.setZero(dim_ext, 1);

  // Set _post quaternion values
  k_state_q_post = k_state_q_pre;
  k_cov_q_post = k_cov_q_pre;

  k_process_q_noise = process_q_noise;
  k_measurement_q_noise = measurement_q_noise;

  m_kappa    = kappa;
  m_kappa_q  = kappa_q;

  b_use_quaternions = true;
}

UKF::~UKF() {
    // Nothing to do, smart pointers should go out of scope
}


/*!
 * \brief UKF::computeKalmanGain
 * Computes the Kalman gain needed in the update step
 */
void UKF::computeKalmanGain(const MatrixXf &cross_correlation,
                            const MatrixXf &covariance_predicted,
                            MatrixXf &kalman_gain) {
  /*!
   * Kalman gain is P_cross_correlation * (P_predicted)^-1
   */

  /* Compute gain */
  MatrixXf innov_cov_inverse;

  innov_cov_inverse = covariance_predicted.inverse();

#ifdef DEBUG_LINUX
  cout << "UKF : covariance_predicted : \n" << covariance_predicted << endl;
  cout << "UKF : pseudo-inverse : \n" << innov_cov_inverse << endl;
#endif

  // Compute Kalman gain
  kalman_gain = cross_correlation * innov_cov_inverse;
}

void UKF::getStatePre(MatrixXf &state_pre) const {

  if (!b_use_quaternions) {
    state_pre = k_state_pre.block(0,0, m_dim, 1); // Only get the estimated variables, no noise elements
  }
  else {
    state_pre.resize (m_dim + k_state_q_pre.rows (), 1);
    state_pre.block(0,0,m_dim, 1) = k_state_pre.block(0,0, m_dim, 1);
    state_pre.block(m_dim,0,k_state_q_pre.rows (), 1) = k_state_q_pre;
  }
}

void UKF::getStatePost(MatrixXf &state_post) const {
  if (!b_use_quaternions) {
    state_post = k_state_post.block(0,0, m_dim, 1);
  } else {
    state_post.resize (m_dim + k_state_q_post.rows (), 1);
    state_post.block(0,0,m_dim, 1) = k_state_post.block(0,0, m_dim, 1);
    state_post.block(m_dim,0,k_state_q_post.rows (), 1) = k_state_q_post;
  }

#ifdef DEBUG_LINUX
  printf("UKF : Post-update state : \n");
  printEigenMatrix (state_post);
#endif
}

void UKF::predict() {
  /*
   * Compute sigma points from mean and cov
   * (Unscented Transform http://en.wikipedia.org/wiki/Unscented_transform)
   */
  if (!m_particles) {
    m_particles.reset( new SigmaPoints(k_state_post, k_cov_post, m_kappa));

    // If a propagation function has not been defined, declare it
    if (_propagateFunc != NULL) {
      m_particles->setPropagationFunction (_propagateFunc);
    }

    // .. Same for our measurement function
    if (_measurementFunc != NULL) {
      m_particles->setMeasurementFunction (_measurementFunc);
    }
  } else {
    /*
     *  Update the sigma set if already present
     */

    // Add the measurement noise and the process noise to the cov matrix :
    k_cov_post.block (m_dim,m_dim,
                      k_process_noise.cols (),
                      k_process_noise.cols ()) = k_process_noise;

    k_cov_post.block (m_dim + k_process_noise.cols (),
                      m_dim + k_process_noise.cols (),
                      k_measurement_noise.cols (),
                      k_measurement_noise.cols ()) = k_measurement_noise;

    // Add the correlations between state errors and process noise :
    // TODO
    m_particles->setState (k_state_post,
                          k_cov_post);
  }

  // Propagate sigma points according to motion model
  propagateSigmaPoints();


  // Fill in predicted state from propagated sigma-points
  m_particles->getPredictedState(k_state_pre,
                                k_cov_pre,
                                k_cov_cross_pred_ref);

  // Add process noise to the predicted covariance (independent additive noise ?)
  k_cov_pre.block(0,0, k_process_noise.rows (), k_process_noise.rows ()) += k_process_noise;


  /*
   * Compute sigma quaternions from mean and cov
   * (Unscented Transform http://en.wikipedia.org/wiki/Unscented_transform)
   */
  if (b_use_quaternions) {
    if (!m_q_particles) {

#ifdef DEBUG_LINUX
      printf("UKF : Allocate q_sigma points\n");
#endif

      cout << k_state_q_post << endl << endl << k_cov_q_post << endl << endl;

      m_q_particles.reset( new SigmaQPoints(k_state_q_post.block(0,0, 3,1),
                                            k_cov_q_post.block(0,0, 3,3),
                                            m_kappa_q) );

      m_q_particles->setProcessNoise (k_process_q_noise);

      // If a propagation function has not been defined, declare it
      if (_propagateQFunc != NULL) {
        printf("UKF : Set propagation Q function\n");
        m_q_particles->setPropagationFunction (_propagateQFunc);
      }

      // .. Same for our measurement function
      if (_measurementQFunc != NULL) {
        printf("UKF : Set measurement Q function\n");
        m_q_particles->setMeasurementFunction (_measurementQFunc);
      }

    }
    else {
      // Rmk : process noise is already added to the covariance estimation inside the
      // sigma quaternions
      m_q_particles->setState (k_state_q_post,
                              k_cov_q_post);
    }

    // Propagate quaternions
    propagateSigmaQPoints ();

    // Get predicted state
    m_q_particles->getStatePre(k_state_q_pre,
                              k_cov_q_pre);
  }

  if (b_first_time) {
    b_first_time = false;
  }

}

/*!
 * \brief propagateSigmaPoints according to the recorded motion model
 */
void UKF::propagateSigmaPoints() {
  /*!
   * Steps :
   * - new angular position is incremented knowing angular speed (supposed constant)
   * - new position is computed knowing axial and angular speed
   */

  if (b_first_time) {
    printf("UKF : set propagation function and propagate sigma points\n");
    if (_propagateFunc == NULL) {
      THROW_ERR("UKF : Propagation function must be defined before propagating sigma points");
    }
    m_particles->setPropagationFunction (_propagateFunc);
    m_particles->propagateSigmaPoints ();
  } else {
    m_particles->propagateSigmaPoints ();
  }
}

/*!
 * \brief propagateSigmaQPoints according to the recorded motion model
 */
void UKF::propagateSigmaQPoints() {
  if (b_first_time) {
    printf("UKF : set propagation function and propagate sigma Q points\n");
    if (_propagateQFunc == NULL) {
      THROW_ERR("UKF : Quaternions propagation function must be defined before propagating sigma points");
    }
    m_q_particles->setPropagationFunction (_propagateQFunc);
    m_q_particles->propagateSigmaQPoints ();
  } else {
    m_q_particles->propagateSigmaQPoints ();
  }

  m_q_particles->computeQMeanAndCovariance(20, 0.001f);
  // TODO: set the number of iterations as a parameter
}


void UKF::setState (const MatrixXf &new_state, const MatrixXf &new_cov) {
  k_state_post.setZero ();
  k_state_post.block(0,0, new_state.cols (), 1)  = new_state;
  k_state_pre   = k_state_post;

  k_cov_post  = new_cov;
  k_cov_pre   = k_cov_post;
}

/*!
 * \brief UKF::setMeasurementNoise
 * \param measurement_noise
 * \warning model_noise must have been set beforehand !
 */
void UKF::setMeasurementNoise(const MatrixXf &measurement_noise) {
  k_measurement_noise = measurement_noise;

  // FIXME : dimension problem if k_cov is too small !
  k_cov_pre.block (m_dim + k_process_noise.cols (),
                   m_dim + k_process_noise.cols (),
                   measurement_noise.cols (),
                   measurement_noise.cols ()) = measurement_noise;

  k_cov_post = k_cov_pre;
}


void UKF::setMeasurementFunction (void (*_meas_function)(const MatrixXf &, MatrixXf &)) {
  // Store the function pointer if no particles are yet allocated
  // else define particles propgation function
  if(m_particles == NULL) {
    _measurementFunc = _meas_function;
  } else {
    m_particles->setMeasurementFunction (_meas_function);
  }
}

void UKF::setMeasurementQFunction (void (*_meas_function)(const Quaternionf &, Quaternionf &)) {
  // Store the function pointer if no particles are yet allocated
  // else define particles propgation function
  if(m_q_particles == NULL) {
    _measurementQFunc = _meas_function;
  } else {
    m_q_particles->setMeasurementFunction (_meas_function);
  }
}

void UKF::setProcessNoise(const MatrixXf &process_noise) {
  k_process_noise = process_noise;
  // FIXME : dimension problem if k_cov is too small !
  k_cov_pre.block (m_dim,
                   m_dim,
                   process_noise.cols (),
                   process_noise.cols ()) = process_noise;

  k_cov_post = k_cov_pre;
}

void UKF::setProcessQNoise(const MatrixXf &process_q_noise) {
  // TODO: foolproof check in dimensions ?
  k_process_q_noise = process_q_noise;
}

void UKF::setMeasurementQNoise (const MatrixXf &measurement_q_noise) {
  // TODO: foolproof check in dimensions ?
  k_measurement_noise = measurement_q_noise;
}

void UKF::setPropagationFunction (void (*_prop_function)(const MatrixXf &, MatrixXf &)) {
  // Store the function pointer if no particles are yet allocated
  // else define particles propgation function
  if(m_particles == NULL) {
    _propagateFunc = _prop_function;
  } else {
    m_particles->setPropagationFunction (_prop_function);
  }
}

void UKF::setPropagationQFunction (void (*_prop_function) (const Quaternionf &, Quaternionf &)) {

  if (m_q_particles == NULL) {
    // Sigma Q points not allocated yet --> store the function pointer
    _propagateQFunc = _prop_function;

  } else {
    // Sigma Q points allocated, define their propagation function
    m_q_particles->setPropagationFunction(_prop_function);
  }
}


/*!
 * \brief UKF::update
 * \param new_measure : a 6x1 float Eigen matrix !
 */
void UKF::update (const MatrixXf &new_measure) {
  // Failsafe stupid tests
  if (new_measure.rows () != m_dim) {
    THROW_ERR("UKF : Wrong measure vector dimension\n");
  } else if (_measurementFunc == NULL) {
    THROW_ERR("UKF : Measurement function is not defined");
  }

  if (b_use_quaternions) {
    THROW_ERR("UKF : Filter including quaternions, you cannot just update vector values\n");
  }

  updateParticles (new_measure);
}


/*!
 * \brief UKF::update
 * \param new_measure : a 6x1 float Eigen matrix !
 */
void UKF::update (const MatrixXf &vec_measure, const MatrixXf &angle_measure) {

  // Failsafe stupid tests
  if (vec_measure.rows () != m_dim) {
    THROW_ERR("UKF : Wrong measure vector dimension\n");
  } else if (_measurementFunc == NULL) {
    THROW_ERR("UKF : Measurement function is not defined");
  }

  if (!b_use_quaternions) {
    THROW_ERR("UKF : Vector-only filter, you cannot update quaternions\n");
  }

  // Update Sigma points
  updateParticles (vec_measure);

  // Update Sigma quaternions
  updateQParticles (angle_measure);
}


void UKF::updateParticles (const MatrixXf &new_measure) {

  // Project on measurement space & Get expected measure
  m_particles->measureSigmaPoints ();

  m_particles->getMeasuredState (k_state_meas,
                                k_cov_meas,
                                k_cov_cross_pred_meas);

  // Compute innovation & Kalman gain:
  k_measure_extended.setZero (new_measure.rows () +
                              k_measurement_noise.cols () +
                              k_process_noise.cols (),
                              1);

  k_measure_extended.block(0,0, new_measure.rows (), 1) = new_measure;

  k_innovation = k_measure_extended - k_state_meas;

  k_cov_extended = k_cov_meas;

  k_cov_extended.block (0,
                        0,
                        k_measurement_noise.cols (),
                        k_measurement_noise.cols ()) += k_measurement_noise;

  computeKalmanGain (k_cov_cross_pred_meas,
                     k_cov_extended,
                     k_gain);

  // Compute a posteriori estimate
  k_state_post  = k_state_pre +
      k_gain * k_innovation;

  k_cov_post    = k_cov_pre -
      k_gain * k_cov_extended * k_gain.transpose(); // State and measurement space are the same, for now..

#ifdef DEBUG
  cout << "\nCorrected state : \n" << k_state_post.block(0,0,6,1) << endl << "***********************" << endl;
#endif

#ifdef DEBUG_LINUX
  printf("k_innovation : \n");
  printEigenMatrix(k_innovation);

  printf("k_gain : \n");
  printEigenMatrix(k_gain);

  printf("k_cov_pre : \n");
  printEigenMatrix(k_cov_pre);

  printf("k_cov_extended : \n");
  printEigenMatrix(k_cov_extended);

  printf("k_cov_cross_pred_meas : \n");
  printEigenMatrix(k_cov_cross_pred_meas);

  printf("k_cov_post : \n");
  printEigenMatrix(k_cov_post);

  printf("k_cov_meas : \n");
  printEigenMatrix(k_cov_meas);
#endif
}


void UKF::updateQParticles(const MatrixXf &new_angular_measure) {

  // Compute innovation & Kalman gain
  k_innovation = new_angular_measure - k_state_q_pre;


  computeKalmanGain(k_cov_q_pre,  // TODO: implement a proper cross-correlation measurement ?
                    k_cov_q_pre + k_measurement_q_noise,
                    k_gain);


  // Compute a posteriori estimate :
  k_state_q_post = k_state_q_pre +
      k_gain * k_innovation;

  k_cov_q_post    = k_cov_q_pre -
      k_gain * (k_cov_q_pre + k_measurement_q_noise) * k_gain.transpose();

#ifdef DEBUG_LINUX
  cout << "UKF : Gain :" << endl << k_gain << endl << endl;

  cout << "UKF : Measure :" << endl << new_angular_measure << endl << endl;

  cout << "UKF : Post-state : " << endl << k_state_q_post << endl << endl<<  k_cov_q_post << endl << endl;
#endif


}

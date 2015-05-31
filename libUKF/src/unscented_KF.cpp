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
         void (*_meas_function)(const MatrixXf &, MatrixXf &),
         void (*_prop_function)(const MatrixXf &, MatrixXf &),
         float kappa)  {

    m_time_step = 1.f;

    m_q_particles.reset();

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

    // Allocate particles (we already know the dimension)
    m_particles.reset( new SigmaPoints(k_state_post, k_cov_post, m_kappa));
    m_particles->setMeasurementFunction (_meas_function);
    m_particles->setPropagationFunction (_prop_function);

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
         void (*_meas_function)(const MatrixXf &, MatrixXf &),
         void (*_prop_function)(const MatrixXf &, MatrixXf &),
         void (*_meas_Qfunction)(const Quaternionf &, Quaternionf &),
         void (*_prop_Qfunction)(const Quaternionf &, Quaternionf &),
         float kappa,
         float kappa_q)  {

    m_time_step = 1.f;

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

    // Allocate particles and set measurement and propagation functions
    m_particles.reset( new SigmaPoints(k_state_post, k_cov_post, m_kappa));
    m_particles->setMeasurementFunction (_meas_function);
    m_particles->setPropagationFunction (_prop_function);

    // Allocate quaternion particles and set measurement and propagation functions
    m_q_particles.reset( new SigmaQPoints(k_state_q_post.block(0,0, 3,1),
                                          k_cov_q_post.block(0,0, 3,3),
                                          m_kappa_q) );

    m_q_particles->setProcessNoise (k_process_q_noise);
    m_q_particles->setMeasurementFunction(_meas_Qfunction);
    m_q_particles->setPropagationFunction(_prop_Qfunction);


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
}

void UKF::predict() {
    /*
   * Compute sigma points from mean and cov
   * (Unscented Transform http://en.wikipedia.org/wiki/Unscented_transform)
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
        // Rmk : process noise is already added to the covariance estimation inside the
        // sigma quaternions
        m_q_particles->setState (k_state_q_post,
                                 k_cov_q_post);

        // Propagate quaternions
        propagateSigmaQPoints ();

        // Get predicted state
        m_q_particles->getStatePre(k_state_q_pre,
                                   k_cov_q_pre);
    }
}

/*!
 * \brief propagateSigmaPoints according to the recorded motion model
 */
void UKF::propagateSigmaPoints() {
    m_particles->propagateSigmaPoints ();
}

/*!
 * \brief propagateSigmaQPoints according to the recorded motion model
 */
void UKF::propagateSigmaQPoints() {
    m_q_particles->propagateSigmaQPoints ();

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
}

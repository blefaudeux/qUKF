#include "motion_filtering.h"

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file motion_filtering.cpp
 *  @brief Defines an instance of UKF filtering
 *  @version 1.0
 *  @date 12-11-2013
 */

using namespace qukf;

Vec3<float> pointPropagation_speed(Vec3<float> const &points_in) {
    // Constant acceleration model... :
    MatrixXf prop_matrix;
    int dim = points_in.rows ();
    prop_matrix.setIdentity(dim,dim);
    prop_matrix(0,3) = 1.f;
    prop_matrix(1,4) = 1.f;
    prop_matrix(2,5) = 1.f;

    return prop_matrix * points_in;
}

void pointPropagation_angularSpeed(const Quaternionf &q_in, Quaternionf &q_out) {
    // Constant velocity model up to now...
    // a bit lame..
    q_out = q_in;

    /*
  // Rotate the positions and speeds from the known iterative rotation :
  Matrix3f rot_mat;
  Vector3f rot_vec;
  Vector3f temp_vec;
  i = 0;

  while (i < sigma_points.size()) {
    rot_vec = sigma_points[i].block<3,1>(6,0) * this->time_step; // Angular offset

    rot_mat = AngleAxisf(rot_vec(0), Vector3f::UnitZ())
        * AngleAxisf(rot_vec(1), Vector3f::UnitY())
        * AngleAxisf(rot_vec(2), Vector3f::UnitZ());

    temp_vec = sigma_points_propagated[i].segment(0,3);
    sigma_points_propagated[i].segment(0,3) = rot_mat * temp_vec;

    temp_vec = sigma_points_propagated[i].segment(3,3);
    sigma_points_propagated[i].segment(3,3) = rot_mat * temp_vec;

    ++i;
  }
  */
}

Vec<float,3> meas_function( Vec<float,3> const &vec_in) {
    return vec_in;
}

void meas_q_function(const Quaternionf &q_in, Quaternionf &q_measured) {
    q_measured = q_in;
}


/*!
 * \brief MotionEstimation::MotionEstimation
 * \param pos
 * \param speed
 */
MotionEstimation::MotionEstimation(const float *variable,
                                   const float ukf_measure_noise = 1.f,
                                   const float ukf_process_noise = 1.f,
                                   const float ukf_kappa = 0.5f) {

    // Constructor for a new pose estimator
    // Contains 2 types of variables :
    // - initial position
    // - initial speed

    _measure.setZero();
    _measure(0) = variable[0];
    _measure(1) = variable[1];
    _measure(2) = variable[2];
    _measure(3) = variable[0];
    _measure(4) = variable[1];
    _measure(5) = variable[2];

    _initial_cov.setIdentity(3, 3);
    _model_noise.setIdentity(3, 3);
    _measurement_noise.setIdentity(3, 3);

    _initial_cov *= 1.f;
    _model_noise *= ukf_process_noise;
    _measurement_noise *= ukf_measure_noise;

    // Allocate UKF and set propagation function
    UKF<float, 3,3>::MeasurementFunc meas = meas_function;
    UKF<float, 3,3>::PropagationFunc prop = pointPropagation_speed;

    m_filter.reset( new UKF<float, 3, 3>(_measure,
                                         _initial_cov, _model_noise,
                                         _measurement_noise,
                                         meas, prop, ukf_kappa));

    m_lastMeasure = _measure;
    _filter_angular_speed = false;
}

MotionEstimation::~MotionEstimation () {
}

void MotionEstimation::setMeasurementSettings(const float ukf_measure_noise,
                                              const float ukf_measure_q_noise,
                                              const float ukf_process_noise,
                                              const float ukf_process_q_noise) {

    // Vector space noise matrices
    _model_noise.setIdentity(3,3);
    _measurement_noise.setIdentity(3,3);

    _model_noise.block(0,0,3,3)       *= ukf_process_noise;
    _measurement_noise.block(0,0,3,3) *= ukf_measure_noise;

    m_filter->setProcessNoise (_model_noise);
    m_filter->setMeasurementNoise (_measurement_noise);

    // Quaternion space noise matrices
    _model_q_noise.setIdentity(3,3);
    _measurement_q_noise.setIdentity(3,3);

    _model_q_noise *= ukf_process_q_noise;
    _measurement_q_noise *= ukf_measure_q_noise;

    m_filter->setProcessQNoise (_model_q_noise);
    m_filter->setMeasurementQNoise (_measurement_q_noise);
}

void MotionEstimation::update(VectorXf const & measure)
{
    m_filter->update(measure);
}

void MotionEstimation::update(VectorXf const & speed, VectorXf const & angular_speed)
{

    // Get latest updated state :
    VecX<float> lastest_state;
    m_filter->getStatePost (lastest_state);

    m_lastMeasure(0) = speed[0];
    m_lastMeasure(1) = speed[1];
    m_lastMeasure(2) = speed[2];

    m_lastMeasure(3) = angular_speed[0];
    m_lastMeasure(4) = angular_speed[1];
    m_lastMeasure(5) = angular_speed[2];

    m_filter->update(m_lastMeasure.head(3), m_lastMeasure.tail(3));
}

void MotionEstimation::predict()
{
    m_filter->predict();
}

void MotionEstimation::getLatestState(VectorXf & state) const
{
    m_filter->getStatePost (state);
}

void MotionEstimation::getPropagatedState(VectorXf & state) const
{
    m_filter->getStatePost( state );
}

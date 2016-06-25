#include "motion_filtering.h"

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file motion_filtering.cpp
 *  @brief Defines an instance of UKF filtering
 *  @version 1.0
 *  @date 12-11-2012
 */

Vec2f pointPropagation_speed(Vec2f const & points_in) {
    // TODO: Set a constant speed model
    return points_in;
}

Vec2f meas_function( Vec2f const &vec_in) {
    return vec_in;
}

MotionEstimation::MotionEstimation(const float *variable,
                                   const float measure_noise = 1.f,
                                   const float process_noise = 1.f,
                                   const float ukf_kappa = 0.5f) {

    m_measure = { variable[0], variable[1]};

    m_initial_cov.setIdentity(2, 2);
    m_model_noise.setIdentity(2, 2);
    m_measurement_noise.setIdentity(2, 2);

    m_model_noise *= process_noise;
    m_measurement_noise *= measure_noise;

    // Allocate UKF and set propagation function
    UKF<float, 2,2>::MeasurementFunc meas = meas_function;
    UKF<float, 2,2>::PropagationFunc prop = pointPropagation_speed;

    m_filter.reset( new UKF<float, 2, 2>(m_measure,
                                         m_initial_cov, m_model_noise,
                                         m_measurement_noise,
                                         meas, prop, ukf_kappa));

    m_lastMeasure = m_measure;
    _filter_angular_speed = false;
}


void MotionEstimation::setMeasurementSettings(const float ukf_measure_noise,
                                              const float ukf_measure_q_noise,
                                              const float ukf_process_noise,
                                              const float ukf_process_q_noise) {

    // Vector space noise matrices
    m_model_noise.setIdentity(2,2);
    m_measurement_noise.setIdentity(2,2);

    m_model_noise.block(0,0,2,2)       *= ukf_process_noise;
    m_measurement_noise.block(0,0,2,2) *= ukf_measure_noise;

    m_filter->setProcessNoise (m_model_noise);
    m_filter->setMeasurementNoise (m_measurement_noise);
}

void MotionEstimation::update(const Vec2f &measure)
{
    m_filter->update(measure);
}

void MotionEstimation::predict()
{
    m_filter->predict();
}

void MotionEstimation::getStatePost(Vec2f & state) const
{
    m_filter->getStatePost( state );
}

#ifndef QKALMANFILTER_H
#define QKALMANFILTER_H

#include "sigma_point.h"
#include "sigma_q_point.h"
#include <memory>
#include <functional>

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file unscented_KF.
 *  @brief Implements a UKF filter with quaternions to filter angular positions
 *  @version 1.1
 */

namespace qukf {

    using namespace Eigen;
    using namespace std;

    template <typename T, size_t DimState, size_t DimMeas>
    class UKF
    {
        public:
            typedef Vec<T, DimState> VecState;
            typedef Vec<T, DimState*2> VecStateExt;
            typedef Vec<T, DimMeas> VecMeas;

            typedef typename SigmaPoints<T,DimState, DimMeas>::MeasurementFunc MeasurementFunc;
            typedef typename SigmaPoints<T,DimState, DimMeas>::PropagationFunc PropagationFunc;

        public:

            UKF( Vec<T, DimState> const & initial_state,
                 MatSquare<T, DimState> const & initial_cov,
                 MatSquare<T, DimState> const & model_noise,
                 MatSquare<T, DimMeas> const & measurement_noise,
                 MeasurementFunc & meas_function,
                 PropagationFunc & prop_function,
                 T alpha)
            {

                // Take extended state vector into account / See Julier and Uhlmann original paper
                int const dim_ext = 2*DimState + DimMeas;

                // Initialize matrices
                k_state_pre = initial_state;
                k_cov_pre = initial_cov;

                k_state_post  = k_state_pre;
                k_cov_post    = k_cov_pre;

                k_process_noise     = model_noise;
                k_measurement_noise = measurement_noise;

                k_cov_cross_pred_meas.fill(T(0.));
                k_cov_cross_pred_ref.fill(T(0.));
                k_gain.fill(T(0.));
                k_innovation.fill(T(0.));

                // Allocate particles
                m_kappa = alpha;
                m_particles.reset( new SigmaPoints<float, DimState*2, DimMeas>(k_state_post, k_cov_post, m_kappa));
                m_particles->setMeasurementFunction( meas_function );
                m_particles->setPropagationFunction( prop_function );

                // Standard "vector space" UKF, no quaternions
                b_use_quaternions = false;
                m_q_particles.reset();
            }

            void predict() {
                // Add the measurement noise and the process noise to the cov matrix :
                k_cov_post.block (DimState,DimState,
                                  k_process_noise.cols (),
                                  k_process_noise.cols ()) = k_process_noise;

                k_cov_post.block (DimState + k_process_noise.cols (),
                                  DimState + k_process_noise.cols (),
                                  k_measurement_noise.cols (),
                                  k_measurement_noise.cols ()) = k_measurement_noise;

                // Add the correlations between state errors and process noise :
                // TODO
                m_particles->setState (k_state_post, k_cov_post);

                // Propagate sigma points according to motion model
                propagateSigmaPoints();

                // Fill in predicted state from propagated sigma-points
                m_particles->getPredictedState(k_state_pre, k_cov_pre, k_cov_cross_pred_ref);

                // Add process noise to the predicted covariance (independent additive noise ?)
                k_cov_pre.block(0,0, k_process_noise.rows (), k_process_noise.rows ()) += k_process_noise;

                if (b_use_quaternions) {
                    m_q_particles->setState (k_state_q_post, k_cov_q_post);

                    // Propagate quaternions
                    propagateSigmaQPoints ();

                    // Get predicted state
                    m_q_particles->getStatePre(k_state_q_pre, k_cov_q_pre);
                }
            }

            void update(const MatX<T> &new_measure) {
                // Failsafe stupid tests
                if (new_measure.rows () != DimState) {
                    THROW_ERR("UKF : Wrong measure vector dimension\n");
                }

                if (b_use_quaternions) {
                    THROW_ERR("UKF : Filter including quaternions, you cannot just update vector values\n");
                }

                updateParticles (new_measure);
            }

            void update (const VecMeas &vec_measure, const MatX<T> &angle_measure) {

                if (!b_use_quaternions) {
                    THROW_ERR("UKF : Vector-only filter, you cannot update quaternions\n");
                }

                // Update Sigma points
                updateParticles (vec_measure);

                // Update Sigma quaternions
                updateQParticles (angle_measure);
            }

            // Gets -----
            void getStatePre(MatX<T> &state_pre) const {
                if (!b_use_quaternions) {
                    state_pre = k_state_pre.block(0,0, DimState, 1); // Only get the estimated variables, no noise elements
                }
                else {
                    state_pre.resize(DimState + k_state_q_pre.rows (), 1);
                    state_pre.block(0,0,DimState, 1) = k_state_pre.block(0,0, DimState, 1);
                    state_pre.block(DimState,0,k_state_q_pre.rows (), 1) = k_state_q_pre;
                }
            }

            void getStatePost(MatX<T> &state_post) const {
                if (!b_use_quaternions) {
                    state_post = k_state_post.block(0,0, DimState, 1);
                } else {
                    state_post.resize (DimState + k_state_q_post.rows (), 1);
                    state_post.block(0,0,DimState, 1) = k_state_post.block(0,0, DimState, 1);
                    state_post.block(DimState,0,k_state_q_post.rows (), 1) = k_state_q_post;
                }
            }

            // Sets -----
            void setState(const MatX<T> &new_state, const MatX<T> &new_cov) {
                k_state_post.setZero ();
                k_state_post.block(0,0, new_state.cols (), 1)  = new_state;
                k_state_pre   = k_state_post;

                k_cov_post  = new_cov;
                k_cov_pre   = k_cov_post;
            }

            void setMeasurementNoise(const MatX<T> &measurement_noise) {
                k_measurement_noise = measurement_noise;

                // FIXME : dimension problem if k_cov is too small !
                k_cov_pre.block (DimState + k_process_noise.cols (),
                                 DimState + k_process_noise.cols (),
                                 measurement_noise.cols (),
                                 measurement_noise.cols ()) = measurement_noise;

                k_cov_post = k_cov_pre;
            }

            void setMeasurementQNoise(const MatX<T> &measurement_q_noise) {
                k_measurement_noise = measurement_q_noise;
            }

            void setProcessNoise(const MatX<T> &process_noise) {
                k_process_noise = process_noise;
                // FIXME : dimension problem if k_cov is too small !
                k_cov_pre.block (DimState,
                                 DimState,
                                 process_noise.cols (),
                                 process_noise.cols ()) = process_noise;

                k_cov_post = k_cov_pre;
            }

            void setProcessQNoise(const MatX<T> &process_q_noise) {
                k_process_q_noise = process_q_noise;
            }

            void setPropagationQFunction (void (*_prop_function)
                                          (const Quaternionf &, Quaternionf &)) {
                if(m_q_particles == NULL) {
                    _propagateQFunc = _prop_function;
                } else {
                    m_q_particles->setPropagationFunction (_prop_function);
                }
            }

            void setMeasurementQFunction(void (*_meas_function)
                                         (const Quaternionf &, Quaternionf &)) {
                if(m_q_particles == NULL) {
                    _measurementQFunc = _meas_function;
                } else {
                    m_q_particles->setMeasurementFunction (_meas_function);
                }
            }

        private :

            void (*_propagateQFunc)(const Quaternionf &, Quaternionf &);

            void (*_measurementQFunc)(const Quaternionf &, Quaternionf &);


            void computeKalmanGain(MatX<T> const & cross_correlation,
                                   MatX<T> const & covariance_predicted,
                                   MatX<T> & kalman_gain)
            {
                kalman_gain = cross_correlation * covariance_predicted.inverse();
            }

            void computeSigmaQPoints() { m_q_particles->computeSigmaQPoints(); }

            void propagateSigmaPoints() { m_particles->propagateSigmaPoints (); }

            void propagateSigmaQPoints() {
                m_q_particles->propagateSigmaQPoints ();
                m_q_particles->computeQMeanAndCovariance(20, 0.001f);
                // TODO: set the number of iterations as a parameter
            }

            void updateParticles(const MatX<T> &new_measure) {
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

                k_cov_extended.block (0, 0, k_measurement_noise.cols (),
                                      k_measurement_noise.cols ()) += k_measurement_noise;

                computeKalmanGain (k_cov_cross_pred_meas,
                                   k_cov_extended,
                                   k_gain);

                // Compute a posteriori estimate
                k_state_post  = k_state_pre + k_gain * k_innovation;

                // State and measurement space are the same, for now..
                k_cov_post    = k_cov_pre -
                        k_gain * k_cov_extended * k_gain.transpose();
            }

            void updateQParticles(const MatX<T> &new_angular_measure) {
                // Compute innovation & Kalman gain
                k_innovation = new_angular_measure - k_state_q_pre;

                // TODO: implement a proper cross-correlation measurement ?
                computeKalmanGain(k_cov_q_pre, k_cov_q_pre + k_measurement_q_noise, k_gain);

                // Compute a posteriori estimate :
                k_state_q_post = k_state_q_pre + k_gain * k_innovation;
                k_cov_q_post    = k_cov_q_pre - k_gain * (k_cov_q_pre + k_measurement_q_noise) * k_gain.transpose();
            }

        private :
            bool  b_first_time;
            bool  b_use_quaternions;

            int   DimState_q;

            T m_kappa, m_kappa_q;

            unique_ptr<SigmaPoints<float, DimState, DimMeas>> m_particles;
            unique_ptr<SigmaQPoints> m_q_particles;

            MatSquare<T, DimState>  k_process_noise;
            MatSquare<T, DimMeas>  k_measurement_noise;

            MatX<T>  k_gain;
            MatX<T>  k_innovation;

            VecState  k_state_pre;
            VecState  k_state_post;
            VecState  k_state_meas;

            MatSquare<T, DimState>  k_cov_pre;
            MatSquare<T, DimState>  k_cov_post;
            MatSquare<T, DimState>  k_cov_meas;
            MatX<T>  k_cov_extended;

            Mat<T, DimMeas, DimState>  k_cov_cross_pred_meas;
            Mat<T, DimMeas, DimState>  k_cov_cross_pred_ref;

            MatX<T>  k_cross_corr;
            MatX<T>  k_measure_extended;

            /* Same matrices needed for quaternion stuff */
            Vec3<T>  k_state_q_pre;
            Vec3<T>  k_state_q_post;

            Mat3<T>  k_cov_q_pre;
            Mat3<T>  k_cov_q_post;
            Mat3<T>  k_process_q_noise;
            Mat3<T>  k_measurement_q_noise;
    };
}

#endif // QKALMANFILTER_H

#ifndef QKALMANFILTER_H
#define QKALMANFILTER_H

#include "sigma_point.h"
#include "sigma_q_point.h"
#include <memory>

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file unscented_KF.
 *  @brief Implements a UKF filter with quaternions to filter angular positions
 *  @version 1.0
 *  @date 12-11-2013
 */

using namespace Eigen;

/*
 * Quaternion tools ...
 */
void QuaternionToEuler(const Quaternion<float> &q_input,
                       Matrix <float, 3,1> &output);

void EulerToQuaternion(const Matrix <float, 3,1> &input,
                       Quaternion <float> &q_out);


/*
 * The UKF filter class...
 */
class UKF
{

    public:

        UKF(const MatrixXf &initial_state,
            const MatrixXf &initial_cov,
            const MatrixXf &model_noise,
            const MatrixXf &measurement_noise,
            void (*_meas_function)(const MatrixXf &, MatrixXf &),
            void (*_prop_function)(const MatrixXf &, MatrixXf &),
            float alpha);


        UKF(const MatrixXf &initial_state,
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
            float kappa_q);

        ~UKF();

        // --------------------------------//
        // -- Intrinsic filtering logic -- //
        void predict();

        /*!
     * \brief update (in case of a vector-only filter)
     * \param new_measure
     */
        void update(const MatrixXf &new_measure);

        /*!
     * \brief update (in case of a vector + quaternions filter)
     * \param angle_measure
     * \param vec_measure
     */
        void update (const MatrixXf &vec_measure, const MatrixXf &angle_measure);

        // -----------------//
        // -- Get values -- //

        /*!
     * \brief getStatePre : get propagated state prior to measurement correction
     * \param state_pre
     */
        void getStatePre(MatrixXf &state_pre) const;

        /*!
     * \brief getStatePost : get state posterior to measurement correction
     * \param state_post
     */
        void getStatePost(MatrixXf &state_post) const;


        // ---------------------//
        // -- Set parameters -- //

        /*!
     * \brief setState : set initial filter position ( (0) by default)
     * \param new_state
     * \param new_cov
     */
        void setState(const MatrixXf &new_state, const MatrixXf &new_cov);

        /*!
     * \brief setMeasurementNoise : set the constant (additive) measure noise
     * \param measurement_noise
     */
        void setMeasurementNoise(const MatrixXf &measurement_noise);

        /*!
     * \brief setMeasurementQNoise
     * \param measurement_q_noise
     */
        void setMeasurementQNoise(const MatrixXf &measurement_q_noise);

        /*!
     * \brief setProcessNoise : set the constant (additive) process noise
     * \param process_noise
     */
        void setProcessNoise(const MatrixXf &process_noise);

        /*!
     * \brief setProcessNoise : set the constant (additive) process noise in quaternion space
     * \param process_noise
     */
        void setProcessQNoise(const MatrixXf &process_q_noise);

        /*!
     * \brief setPropagationQFunction
     *  Define the function used to propagate the quaternions
     */
        void setPropagationQFunction (void (*_prop_function) (const Quaternionf &, Quaternionf &));

        /*!
     * \brief setMeasurementQFunction :
     *  Define  the function used to propagate the quaternion sigma points
     */
        void setMeasurementQFunction(void (*_meas_function)(const Quaternionf &, Quaternionf &));

    private :
        /*
     * Sigma_points : in vector, measure and quaternions space
     */
        void (*_propagateFunc)(const MatrixXf &, MatrixXf &);

        void (*_measurementFunc)(const MatrixXf &, MatrixXf &);

        void (*_propagateQFunc)(const Quaternionf &, Quaternionf &);

        void (*_measurementQFunc)(const Quaternionf &, Quaternionf &);

        /*!
     *  Methods (private) :
     */

        /*!
     * \brief Computes the Kalman gain needed in the update step
     */
        void computeKalmanGain(const MatrixXf &cross_correlation,
                               const MatrixXf &covariance_predicted,
                               MatrixXf &kalman_gain);

        /*!
     * \brief Compute the new sigma points in quaternions space
     */
        void computeSigmaQPoints();

        /*!
     * \brief Propagate sigma points from the vector space
     */
        void propagateSigmaPoints();

        /*!
     * \brief Propagate Quaternions sigma points
     */
        void propagateSigmaQPoints();

        /*!
     * \brief updateParticles
     * \param new_measure
     */
        void updateParticles(const MatrixXf &new_measure);

        /*!
     * \brief updateParticles
     * \param new_measure
     */
        void updateQParticles(const MatrixXf &new_angular_measure);

    private :
        bool  b_first_time;
        bool  b_use_quaternions;

        int   m_dim;
        int   m_dim_q;

        float m_time_step;
        float m_kappa, m_kappa_q;

        std::unique_ptr<SigmaPoints> m_particles;
        std::unique_ptr<SigmaQPoints> m_q_particles;

        /*
     *  State vector includes (excepting quaternions)
     *  x/y/z           - position
     *  x°/y°/z°        - linear speed
     *  w_x / w_y / w_z - rotation rate
     *  theta/phi/psi   - angular position (possibly absent, depending on quaternions)
     */
        MatrixXf  k_process_noise;
        MatrixXf  k_measurement_noise;

        MatrixXf  k_gain;
        MatrixXf  k_innovation;
        MatrixXf  k_expected_cov;

        MatrixXf  k_state_pre;
        MatrixXf  k_state_post;
        MatrixXf  k_state_meas;

        MatrixXf  k_cov_pre;
        MatrixXf  k_cov_post;
        MatrixXf  k_cov_meas;
        MatrixXf  k_cov_extended;
        MatrixXf  k_cov_cross_pred_meas;
        MatrixXf  k_cov_cross_pred_ref;


        MatrixXf  k_cross_corr;
        MatrixXf  k_measure_extended;

        /* Same matrices needed for quaternion stuff */
        Vector3f  k_state_q_pre;
        Vector3f  k_state_q_post;
        Matrix3f  k_cov_q_pre;
        Matrix3f  k_cov_q_post;
        Matrix3f  k_process_q_noise;
        Matrix3f  k_measurement_q_noise;
};

#endif // QKALMANFILTER_H

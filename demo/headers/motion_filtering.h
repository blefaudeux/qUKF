#ifndef MOTION_FILTERING_H
#define MOTION_FILTERING_H

#include "unscented_KF.h"

using namespace qukf;

class MotionEstimation
{
    public:
        MotionEstimation(const float *variable,
                         const float ukf_measure_noise,
                         const float ukf_process_noise,
                         const float ukf_kappa);

        MotionEstimation(const float *speed,
                         const float *angular_speed,
                         const float ukf_measure_noise,
                         const float ukf_measure_q_noise,
                         const float ukf_process_noise,
                         const float ukf_process_q_noise,
                         const float ukf_kappa,
                         const float ukf_kappa_q);

        ~MotionEstimation();

        void getLatestState(VectorXf & state) const;

        void getPropagatedState(VectorXf & state_out) const;

        void predict();

        void setMeasurementSettings(const float ukf_measure_noise,
                                    const float ukf_measure_q_noise,
                                    const float ukf_process_noise,
                                    const float ukf_process_q_noise);

        void update(VectorXf const & measure);

        void update(VectorXf const & speed, VectorXf const & angular_speed);

    public:
        float dt;


    private :
        bool _filter_angular_speed;

        std::vector<float> past_positions;
        std::vector<float[3]> past_speed;

        std::unique_ptr<UKF<float,3,3>> m_filter;

        UKF<float,3,3>::VecMeas _measure;
        VectorXf m_lastMeasure;
        MatrixXf _initial_cov;
        MatrixXf _model_noise;
        MatrixXf _measurement_noise;

        MatrixXf _initial_q_cov;
        MatrixXf _model_q_noise;
        MatrixXf _measurement_q_noise;


};

#endif // POSE_ESTIMATION_H

#ifndef MOTION_FILTERING_H
#define MOTION_FILTERING_H

#include "unscented_KF.h"

using namespace std;
using namespace qukf;


class MotionEstimation
{
    public:

        MotionEstimation(const float *variable,
                         const float measure_noise,
                         const float process_noise,
                         const float ukf_kappa);

        void getStatePost(Vec2f & state) const;

        void predict();

        void setMeasurementSettings(const float ukf_measure_noise,
                                    const float ukf_measure_q_noise,
                                    const float ukf_process_noise,
                                    const float ukf_process_q_noise);

        void update(Vec2f const & measure);


    public:
        float dt;


    private :
        bool _filter_angular_speed;

        vector<float> past_positions;
        vector<float[3]> past_speed;

        unique_ptr<UKF<float,2,2>> m_filter;

        UKF<float,2,2>::VecMeas m_measure;
        Vec2f m_lastMeasure;
        Mat2f m_initial_cov;
        Mat2f m_model_noise;
        Mat2f m_measurement_noise;

};

#endif // POSE_ESTIMATION_H

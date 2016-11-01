#ifndef SIGMAQPOINT_H
#define SIGMAQPOINT_H


/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file sigma_q_point.h
 *  @version 1.0
 *  @date 12-11-2013
 */

#include "eigen_tools.h"
#include <vector>


/*!
 * \brief The SigmaQPoints class
 *
 * Represents sigma points from a 3D angular set,
 * in the quaternion space
 * Dimensions are fixed to 3 right now, do not use elsewere..
 *
 */

namespace qukf
{
    using namespace Eigen;
    using namespace std;

    class SigmaQPoints
    {

        public:
            SigmaQPoints(const Matrix<float, 3, 1> &angles_mean,
                         const Matrix<float, 3, 3> &angles_cov,
                         float kappa);


            void computeQMeanAndCovariance(int max_iterations,
                                           float max_err_norm );


            void projectSigmaQPoints();


            void  propagateSigmaQPoints ();

            void getStatePre(Matrix<float, 3, 1> &mean,
                             Matrix<float, 3, 3> &cov) const;

            void setPropagationFunction(void (*_prop_function)(const Quaternionf &, Quaternionf &));

            void setMeasurementFunction(void (*meas_function)(const Quaternionf &, Quaternionf &));

            void setProcessNoise(const MatrixXf &process_noise);

            void setState(const MatrixXf &mean, const MatrixXf &cov) ;

        private:

            void computeSigmaQPoints();

            void quaternionToEuler(const Quaternion <float> &q_input,
                                   Vector3f &output) const;

            void eulerToQuaternion(const Vector3f &input,
                                   Quaternion <float> &q_out) const;

            void angleAxisNormalize(const AngleAxis<float> &input, AngleAxis <float> &output);

            void averageQuaternionsAngleAxis(const vector <Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                             const vector <float>       &q_weight,
                                             Quaternionf  &q_avg);

            void averageQuaternionsSlerp(const vector <Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                         const vector <float> &q_weight,
                                         Quaternionf &q_avg);

            void averageQuaternionsIterative(const vector <Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                             const vector <float> &q_weight,
                                             int    n_iterations,
                                             float  max_err,
                                             Quaternionf &q_avg);

            void pickMeanQuaternion(const vector<Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                    const vector<float> &q_weight,
                                    Quaternionf *mean_quaternion);

            float distQuaternions(const Quaternionf &q_1,
                                  const Quaternionf &q_2);

        private :
            int _dim;
            float _kappa;

            Quaternionf _q_mean_reference;
            Quaternionf _q_mean_predicted;
            Quaternionf _q_mean_measure;
            Quaternionf _q_cov;

            Vec3<float>    _mean_reference;
            Vector3f    _mean_predicted;
            Vector3f    _mean_measure;

            Matrix3f    _cov_reference;
            Matrix3f    _cov_predicted;
            Matrix3f    _cov_measure;

            vector <float> _weight_mean;
            vector <float> _weight_cov;

//            Implementation is different from sigma_point.cpp,
//            we don't introduce noisy propagation with and extended state vector,
//            we consider a static process noise which we add to covariance
//            prior to generating sigma points
            Matrix3f  _process_noise;

            vector <Quaternionf, aligned_allocator<Quaternionf> >  _q_sigma_pt;
            vector <Vector3f, aligned_allocator<Vector3f> >        _sigma_pt;

            void (*_propagateFunc)(const Quaternionf &, Quaternionf &);
            void (*_measurementFunc)(const Quaternionf &, Quaternionf &);

    };
}



#endif // SIGMAQPOINT_H

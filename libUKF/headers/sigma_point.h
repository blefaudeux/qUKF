#ifndef SIGMAPOINT_H
#define SIGMAPOINT_H

#include "statistic_tools.h"
#include "eigen_tools.h"
#include <cassert>
#include <functional>

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file sigma_point.h
 *  @brief Implements a class of sigma points, to be used with UKF
 *  @version 1.0
 *  @date 12-11-2013
 */


using namespace Eigen;
using namespace std;
using namespace qukf;

/*
 * Vector space sigma point
 */

template <typename T>
class SigmaPoints {
    public:
        SigmaPoints( MatX<T> const &mean,
                     MatX<T> const  &cov,
                     float kappa) {

            _dim   = mean.rows ();

            _mean_reference   = mean;

            _cov_reference        = cov;
            _cov_measure          = cov; // Useless ?
            _cov_cross_pred_meas  = cov; // Useless ?
            _cov_cross_pred_state = cov; // Useless ?

            _kappa  = kappa;

            // Create sigma points distribution (Unscented Transform)
            computeSigmaPoints ();
        }


        /*!
     * \brief computeSigmaPoints
     * Update the set of sigma points, knowing mean, cov, and dispersion parameters
     * // Unscented Transform //
     */
        void computeSigmaPoints() {
            int n_sigma_points = 2 * _dim + 1;

            _point_reference.resize(n_sigma_points);
            _point_predicted.resize(n_sigma_points);
            _point_measure.resize(n_sigma_points);
            _weight_mean.resize (n_sigma_points);
            _weight_cov.resize (n_sigma_points);

            // Mean
            _point_reference[0]  = _mean_reference;
            _weight_mean[0] = _kappa / (_dim + _kappa);
            _weight_cov[0]  = _kappa / (_dim + _kappa);;

            // Compute "square root matrix" for covariance to get sigma points
            LLT <MatX<T>> lltOfCov(_cov_reference);
            MatX<T> L = lltOfCov.matrixL ();

            // Distributed points..
            for (int i=1; i<=_dim; ++i) {
                // Sigma point position
                _point_reference[i  ]       = _mean_reference +  L.col (i-1) * sqrt(_dim + _kappa);
                _point_reference[_dim + i]  = _mean_reference -  L.col (i-1) * sqrt(_dim + _kappa);

                // Weights : same weights for everyone right now..
                _weight_mean[i]       = 1.f / (2.f * (_kappa + _dim));
                _weight_mean[_dim + i] = _weight_mean[i];

                _weight_cov[i]        = 1.f / (2.f * (_kappa + _dim));
                _weight_cov[_dim + i]  = _weight_cov[i];
            }
        }

        /*!
     * \brief getState
     * \param _mean
     * \param _cov
     */
        void getState(MatX<T> &mean,
                      MatX<T> &cov) const {
            mean  = _mean_reference;
            cov   = _cov_reference;
        }

        /*!
     * \brief getMeasuredState
     * \param mean
     * \param cov
     * \param cov_cross
     */
        void getMeasuredState(MatX<T> &mean,
                              MatX<T> &cov_measure,
                              MatX<T> &cov_cross) const {
            mean        = _mean_measure;
            cov_measure = _cov_measure;
            cov_cross   = _cov_cross_pred_meas;
        }

        /*!
     * \brief getPredictedState
     * \param mean
     * \param cov
     * \param cov_cross
     */
        void getPredictedState(MatX<T> &mean,
                               MatX<T> &cov,
                               MatX<T> &cov_cross) const {
            mean      = _mean_predicted;
            cov       = _cov_predicted;
            cov_cross = _cov_cross_pred_state;
        }

        /*!
     * \brief measureSigmaPoints : project the sigma points on the measurement space
     * populates the _point_measurement_space vector
     */
        void measureSigmaPoints() {
            if (!_measurementFunc) {
                THROW_ERR("Sigma points : measurement function undefined");
            }

            _weight = 0.f;

            // Compute the projection of the sigma points onto the measurement space
            _point_measure.resize (_point_predicted.size ());

            for (unsigned int i=0; i<_point_predicted.size (); ++i) {
                _measurementFunc(_point_predicted[i], _point_measure[i]);
            }

            // Compute the mean of the measured points :
            _mean_measure.setZero (_dim, 1);
            for (unsigned int i=0; i<_point_measure.size (); ++i) {
                _mean_measure += _point_measure[i] * _weight_mean[i];
                _weight += _weight_mean[i];
            }

            if (_weight != 0.f) {
                _mean_measure /= _weight;
            }

            // Compute the intrinsic covariance of the measured points :
            _cov_measure.setZero (_dim, _dim);
            _weight = 0.f;

            for (unsigned int i=0; i<_point_measure.size (); ++i) {
                _cov_measure += _weight_cov[i] * ((_point_measure[i] - _mean_measure)
                                                  * (_point_measure[i].transpose() - _mean_measure.transpose()));

                _weight += _weight_cov[i];
            }

            if (_weight != 0.f) {
                _cov_measure /= _weight;
            }

            // Compute the crossed covariance between measurement space and intrisic space
            _cov_cross_pred_meas.setZero (_dim, _dim);

            for (unsigned int i=0; i<_point_predicted.size (); ++i) {
                _cov_cross_pred_meas += _weight_cov[i] * ((_point_measure[i] - _mean_measure)
                                                          * (_point_predicted[i].transpose() - _mean_predicted.transpose()));
            }

            if (_weight != 0.f) {
                _cov_cross_pred_meas /= _weight;
            }
        }


        /*!
     * \brief propagateSigmaPoints
     */
        void propagateSigmaPoints() {
            if (!_propagateFunc) {
                THROW_ERR("SPoint : Invalid propagation function");
            }

            // Propagate existing set of points (supposed to be representative)
            _point_predicted.resize(_point_reference.size());

            unsigned  int i=0;
            while (i < this->_point_reference.size()) {
                _propagateFunc(_point_reference[i], _point_predicted[i]);
                ++i;
            }

            // Update statistics
            updateMean(_point_predicted,
                       _mean_predicted);

            updateCov (_point_predicted,
                       _mean_predicted,
                       _point_predicted,
                       _mean_predicted,
                       _cov_predicted);

            updateCov (_point_predicted,
                       _mean_predicted,
                       _point_reference,
                       _mean_reference,
                       _cov_cross_pred_state);
        }

        void setState(const MatX<T> &mean, const MatX<T> &cov) {
            _mean_reference = mean;
            _cov_reference  = cov;

            computeSigmaPoints ();
        }

        void setPropagationFunction(std::function<void(const MatX<T> &, MatX<T> &)> & _prop_function) {
            _propagateFunc = _prop_function;
        }

        void setMeasurementFunction(std::function<void(const MatX<T> &, MatX<T> &)> &_meas_function) {
            _measurementFunc = _meas_function;
        }

    private:

        std::function<void(const MatX<T> &, MatX<T> &)> _propagateFunc;
        std::function<void(const MatX<T> &, MatX<T> &)> _measurementFunc;


        void updateMean(const std::vector < MatX<T> > &sigma_points,
                        MatX<T> &mean) {
            mean.setZero (_dim, 1);

            _weight = 0.f;

            unsigned int i = 0;

            while (i < sigma_points.size ()) {
                mean += sigma_points[i] * _weight_mean[i];

                _weight += _weight_mean[i];
                ++i;
            }

            // Normalize weights
            if (_weight != 0.f) {
                mean /= _weight;
            }
        }

        void updateCov(const vector < MatX<T> > &sigma_points_1, const MatX<T> &mean_1,
                       const vector < MatX<T> > &sigma_points_2, const MatX<T> &mean_2,
                       MatX<T> &cov) {
            // Deal with mistakes...
            if (mean_1.cols() != mean_2.cols()) {
                THROW_ERR("SPoints : vectors have different sizes");
            }

            // Compute inter-state covariance
            cov.setZero (_dim, _dim);

            _weight = 0.f;

            for (unsigned int i=0; i<_point_predicted.size (); ++i) {
                cov += _weight_cov[i] * ((sigma_points_1[i] - mean_1)
                                         * (sigma_points_2[i].transpose() - mean_2.transpose()));

                _weight += _weight_cov[i];
            }

            // Normalize weights
            if (_weight != 0.f) {
                cov /= _weight;
            }
        }

    private :
        int     _dim; // Move to template
        float   _kappa; // Sigma distribution parameters
        float   _weight;


        MatX<T> _mean_reference;
        MatX<T> _mean_predicted;
        MatX<T> _mean_measure;

        MatX<T> _cov_reference;          // Covariance of the initial state
        MatX<T> _cov_predicted;          // Covariance of the predicted state
        MatX<T> _cov_measure;            // Covariance within measured points
        MatX<T> _cov_cross_pred_meas;    // Covariance between propagated and measured points  (~Measurement noise)
        MatX<T> _cov_cross_pred_state;   // Covariance between propagated and initial points   (~Process noise)

        /*
     *  Sigma points
     */
        vector < MatX<T> > _point_predicted;
        vector < MatX<T> > _point_reference;
        vector < MatX<T> > _point_measure;

        vector <float> _weight_mean;
        vector <float> _weight_cov;
};

#endif // SIGMAPOINT_H

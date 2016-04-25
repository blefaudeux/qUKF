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

template <typename T, size_t DimState, size_t DimMeas>
class SigmaPoints {
        typedef Vec<T, DimState> VecState;
        typedef Vec<T, DimMeas> VecMeas;
        typedef MatSquare<T, DimState> MatCovState;
        typedef MatSquare<T, DimMeas>  MatCovMeas;

    public:
        typedef function<VecMeas(VecMeas const &)> MeasurementFunc;
        typedef function<VecState(VecState const &)> PropagationFunc;

    public:
        SigmaPoints( VecState const &mean, Mat<T, DimState, DimMeas> const  &cov, float kappa)
        {
            m_mean_ref   = mean;

            m_cov_ref        = cov;
            m_cov_meas          = cov; // Useless ?
            m_cov_x_pred_meas  = cov; // Useless ?
            m_cov_x_pred_state = cov; // Useless ?

            m_kappa  = kappa;

            // Create sigma points distribution (Unscented Transform)
            computeSigmaPoints ();
        }

        void computeSigmaPoints()
        {
            int const n_sigma_points = 2 * DimState + 1;

            m_point_ref.setZero(DimState, n_sigma_points);
            m_point_pred.resize(n_sigma_points);
            m_point_meas.resize(n_sigma_points);
            m_weight_mean.resize (n_sigma_points);
            m_weight_cov.resize (n_sigma_points);

            // Mean
            m_point_ref.col(0)  = m_mean_ref;
            m_weight_mean[0] = m_kappa / (DimState + m_kappa);
            m_weight_cov[0]  = m_kappa / (DimState + m_kappa);;

            // Compute "square root matrix" for covariance to get sigma points
            LLT <MatX<T>> lltOfCov(m_cov_ref);
            MatX<T> L = lltOfCov.matrixL ();

            // Distributed points..
            // TODO : Rewrite this matricially !
            for (int i=1; i<=DimState; ++i) {
                // Sigma point position
                m_point_ref.col(i)              = m_mean_ref +  L.col(i-1) * sqrt(DimState + m_kappa);
                m_point_ref.col(DimState + i)   = m_mean_ref -  L.col(i-1) * sqrt(DimState + m_kappa);

                // Weights : same weights for everyone right now..
                m_weight_mean[i]       = 1.f / (2.f * (m_kappa + DimState));
                m_weight_mean[DimState + i] = m_weight_mean[i];

                m_weight_cov[i]        = 1.f / (2.f * (m_kappa + DimState));
                m_weight_cov[DimState + i]  = m_weight_cov[i];
            }
        }

        void getState(VecState &mean, MatCovState &cov) const
        {
            mean  = m_mean_ref;
            cov   = m_cov_ref;
        }

        void getMeasuredState(VecMeas &mean, MatCovMeas &cov_measure, MatCovState &cov_cross) const
        {
            mean        = m_mean_meas;
            cov_measure = m_cov_meas;
            cov_cross   = m_cov_x_pred_meas;
        }

        void getPredictedState(VecState &mean, MatCovState &cov, MatCovState &cov_cross) const
        {
            mean      = m_mean_pred;
            cov       = m_cov_predicted;
            cov_cross = m_cov_x_pred_state;
        }

        void measureSigmaPoints()
        {
            // Compute the projection of the sigma points onto the measurement space
            m_point_meas.resize (m_point_pred.size ());

            for (unsigned int i=0; i<m_point_pred.cols(); ++i) {
                m_point_meas.col(i) = _measurementFunc(m_point_pred.col(i));
            }

            // Compute the mean of the measured points :
            updateMean(m_point_meas, m_mean_meas);

            // Compute the intrinsic covariance of the measured points :
            auto const zm_meas = m_point_meas.colwise() - m_mean_meas;
            m_cov_meas = m_weight_cov.normalized().asDiagonal() * zm_meas * zm_meas.transpose();

            // Compute the crossed covariance between measurement space and intrisic space
            auto const zm_pred = m_point_pred.colwise() - m_mean_pred;
            m_cov_x_pred_meas = m_weight_cov.normalized().asDiagonal() * zm_meas * zm_pred.transpose();
        }

        void propagateSigmaPoints() {
            // Propagate existing set of points (supposed to be representative)
            m_point_pred.resize(m_point_ref.size());

            unsigned  int i=0;
            while (i < m_point_ref.cols()) {
                m_point_pred.col(i) = _propagateFunc(m_point_ref.col(i++));
            }

            // Update statistics
            updateMean(m_point_pred, m_mean_pred);

            updateCov (m_point_pred, m_mean_pred,
                       m_point_pred, m_mean_pred,
                       m_cov_predicted);

            updateCov (m_point_pred, m_mean_pred,
                       m_point_ref, m_mean_ref,
                       m_cov_x_pred_state);
        }

        void setState(VecState const & mean, MatCovState const & cov)
        {
            m_mean_ref = mean;
            m_cov_ref  = cov;

            computeSigmaPoints ();
        }

        void setPropagationFunction(PropagationFunc & _prop_function)
        {
            _propagateFunc = _prop_function;
        }

        void setMeasurementFunction(MeasurementFunc & _meas_function)
        {
            _measurementFunc = _meas_function;
        }

    private:

        void updateMean( MatR<T, DimState> const & sp, VecState &mean)
        {
            mean = (sp * m_weight_mean.asDiagonal()).rowwise().mean();
        }

        void updateCov(MatR<T, DimState> const & sigma_points_1, VecState const & mean_1,
                       MatR<T, DimState> const & sigma_points_2, VecState const & mean_2,
                       MatSquare<T, DimState> &cov)
        {
            auto sgp1_zm = sigma_points_1.colwise() - mean_1;
            auto sgp2_zm = sigma_points_2.colwise() - mean_2;

            cov = m_weight_cov.normalized().asDiagonal() * sgp1_zm * sgp2_zm.transpose();
        }

    private :
        PropagationFunc _propagateFunc;
        MeasurementFunc _measurementFunc;

        float   m_kappa; // Sigma distribution parameters

        VecState    m_mean_ref;
        VecState    m_mean_pred;
        VecMeas     m_mean_meas;

        MatCovState                 m_cov_ref;         // Covariance of the initial state
        MatCovState                 m_cov_predicted;   // Covariance of the predicted state
        MatCovMeas                  m_cov_meas;        // Covariance within measured points
        MatCovState                 m_cov_x_pred_state; // Covariance between propagated and initial points   (~Process noise)
        Mat<T, DimState, DimMeas>      m_cov_x_pred_meas;  // Covariance between propagated and measured points  (~Measurement noise)

        MatR<T, DimState>    m_point_pred;
        MatR<T, DimState>    m_point_ref;
        MatR<T, DimMeas>     m_point_meas;

        VecState m_weight_mean;
        VecState m_weight_cov;
};

#endif // SIGMAPOINT_H

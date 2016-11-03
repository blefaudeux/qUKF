#ifndef SIGMAPOINT_H
#define SIGMAPOINT_H

#include "eigen_tools.h"
#include <cassert>
#include <functional>

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file sigma_point.h
 *  @brief Implements a class of sigma points, to be used with UKF
 *  The structure is not generic, the noise vector being of the state dimension
 *  @version 2.0
 *  @date 12-11-2013
 */


namespace qukf {

    using namespace Eigen;
    using namespace std;

    /*
     * Vector space sigma point
     */

    template <typename T, size_t DimState, size_t DimMeas>
    class SigmaPoints
    {
            typedef Vec<T, DimState> VecState;
            typedef Vec<T, DimMeas> VecMeas;

            typedef MatSquare<T, DimState> MatState;
            typedef MatSquare<T, DimMeas> MatMeas;

            // Extended state definitions
            // See the seminal paper from Julier & Uhlmann
            static size_t const DimExtMeas = DimMeas*2;
            static size_t const DimExtState = DimState * 2;

            typedef Vec<T, DimExtMeas> VecExtMeas;
            typedef Vec<T, DimExtState> VecExtState;

            typedef MatSquare<T, DimState> MatExtState;
            typedef MatSquare<T, DimMeas>  MatExtMeas;

        public:
            typedef function<VecExtMeas(VecExtMeas const &)> MeasurementFunc;
            typedef function<VecExtState(VecExtState const &)> PropagationFunc;

        public:
            SigmaPoints( VecState const &mean,
                         Mat<T, DimState, DimMeas> const  &cov,
                         Mat<T, DimState, DimState> const  &process_noise,
                         Mat<T, DimState, DimMeas> const  &process_cross_noise,
                         float kappa)
            {
                m_mean_ref   = mean;

                m_cov_ref = m_cov_meas = m_cov_pred = cov;

                m_process_noise = process_noise;
                m_process_cross_noise = process_cross_noise;

                m_kappa  = kappa;

                // Create sigma points distribution (Unscented Transform)
                computeSigmaPoints ();
            }

            void computeSigmaPoints()
            {
                // The sigma points are generated over an extended state,
                // as proposed by Julier and Uhlmann

                int const n_sigma_points = 4 * DimState + 1;

                m_point_ref.setZero(DimState * 2, n_sigma_points);
                m_point_pred.setZero(DimState * 2, n_sigma_points);

                // Mean
                m_point_ref.col(0).head(DimState)  = m_mean_ref;
                m_weight_mean(0) = m_kappa / (DimState + m_kappa);
                m_weight_cov(0)  = m_kappa / (DimState + m_kappa);

                // Compute "square root matrix" for covariance to get sigma points
                // TODO: Keep the main cov matrix constant, just update the "real" state part
                int const dimExt = 2 * DimState;
                Mat<T, dimExt, dimExt> covExt;
                covExt.topLeftCorner(DimState,DimState) = m_cov_ref;
                covExt.block(DimState,DimState, DimState, DimState) = m_process_noise;

                covExt.block(0,DimState, DimState, DimState) = m_process_cross_noise;
                covExt.block(DimState,0, DimState, DimState) = m_process_cross_noise;

                LLT <MatX<T>> lltOfCov(covExt);
                MatX<T> L = lltOfCov.matrixL ();

                // Distributed points.. | This can be tuned if needed
                Vec<T, DimState * 2> meanExt;

                meanExt.head(DimState) = m_mean_ref;
                meanExt.fill(T(0.)); // TODO: introduce the possibility to have noise bias

                for (unsigned i=1; i <= 2*DimState; ++i)
                {
                    m_point_ref.col(i)              = meanExt +  L.col(i-1) * sqrt(DimState + m_kappa);
                    m_point_ref.col(DimState + i)   = meanExt -  L.col(i-1) * sqrt(DimState + m_kappa);

                    m_weight_mean[i]       = 1.f / (2.f * (m_kappa + DimState));
                    m_weight_mean[DimState + i] = m_weight_mean[i];

                    m_weight_cov[i]        = 1.f / (2.f * (m_kappa + DimState));
                    m_weight_cov[DimState + i]  = m_weight_cov[i];
                }
            }

            void getState(VecState &mean, MatExtState &cov) const
            {
                mean  = m_mean_ref;
                cov   = m_cov_ref;
            }

            void getMeasuredState(VecState &mean, MatMeas &cov_measure, MatMeas &cov_cross) const
            {
                mean        = m_mean_meas.top(DimState);
                cov_measure = m_cov_meas.topLeftCorner(DimMeas, DimMeas);
                cov_cross   = m_cov_meas.topRightCorner(DimMeas, DimMeas);
            }

            void getPredictedState(VecMeas &mean, MatExtState &cov, MatExtState &cov_cross) const
            {
                mean      = m_mean_pred.top(DimMeas);
                cov       = m_cov_pred;
                cov_cross = m_cov_meas;
            }

            void measureSigmaPoints()
            {
                // Compute the projection of the sigma points onto the measurement space
                for (unsigned i=0; i<m_point_pred.cols(); ++i) {
                    m_point_meas.col(i).head(DimState) = _measurementFunc(m_point_pred.col(i).head(DimState));
                    m_point_meas.col(i).tail(DimState) = _measurementFunc(m_point_pred.col(i).tail(DimState));
                }

                // Compute the mean of the measured points :
                updateMean(m_point_meas, m_mean_meas);

                // Compute the intrinsic covariance of the measured points :
                Vec<T, DimState * 2> measExt;
                measExt.setZero();
                measExt.head(DimState) = m_mean_meas;

                auto const zm_meas = (m_point_meas.colwise() - measExt).eval();
                m_cov_meas = m_weight_cov.normalized().asDiagonal() * zm_meas * zm_meas.transpose();

                // Compute the crossed covariance between measurement space and intrisic space
                Vec<T, DimState * 2> predExt;
                predExt.setZero();
                predExt.head(DimState) = m_mean_pred;

                auto const zm_pred = (m_point_pred.colwise() - predExt).eval();
                m_cov_meas = ((m_weight_cov.normalized().asDiagonal() * zm_meas * zm_pred.transpose()).topLeftCorner<DimState>()).eval();
            }

            void propagateSigmaPoints()
            {
                // Propagate existing set of points (supposed to be representative)
                for(unsigned i=0; i < m_point_ref.cols(); ++i )
                {
                    m_point_pred.col(i) = _propagateFunc(m_point_ref.col(i));
                }

                // Update statistics
                updateMean(m_point_pred, m_mean_pred);

                updateCov (m_point_pred, m_mean_pred,
                           m_point_pred, m_mean_pred,
                           m_cov_pred);

                updateCov (m_point_pred, m_mean_pred,
                           m_point_ref, m_mean_ref,
                           m_cov_x_pred_state);
            }

            void setState(VecState const & mean, MatExtState const & cov)
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
                mean = (sp * m_weight_mean.asDiagonal()).rowwise().mean().head(mean.rows());
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

            // State
            VecExtState m_mean_ref;
            VecExtState m_mean_pred;
            VecExtMeas  m_mean_meas;

            MatExtState m_cov_ref;
            MatExtState m_cov_pred;
            MatExtMeas  m_cov_meas;

            Mat<T, DimState, DimState>  m_process_noise;
            Mat<T, DimState, DimMeas>   m_process_cross_noise;

            // Sigma points machinery
            MatR<T, DimExtState>    m_point_pred;
            MatR<T, DimExtState>    m_point_ref;
            MatR<T, DimMeas>        m_point_meas;

            Vec<T, DimState * 2> m_weight_mean;
            Vec<T, DimState * 2> m_weight_cov;
    };
}



#endif // SIGMAPOINT_H

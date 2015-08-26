#ifndef SIGMAPOINT_H
#define SIGMAPOINT_H

#include "statistic_tools.h"
#include "eigen_tools.h"

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file sigma_point.h
 *  @version 1.0
 *  @date 12-11-2013
 */


using namespace Eigen;
using namespace std;

/*
 * Vector space sigma point
 */
class SigmaPoints {
  private :
    int _dim;
    float _kappa; // Sigma distribution parameters
    float _weight;
    /*
     * State estimate
     */
    MatrixXf _mean_reference;
    MatrixXf _mean_predicted;
    MatrixXf _mean_measure;

    MatrixXf _cov_reference;          // Covariance of the initial state
    MatrixXf _cov_predicted;          // Covariance of the predicted state
    MatrixXf _cov_measure;            // Covariance within measured points
    MatrixXf _cov_cross_pred_meas;    // Covariance between propagated and measured points  (~Measurement noise)
    MatrixXf _cov_cross_pred_state;   // Covariance between propagated and initial points   (~Process noise)

    /*
     *  Sigma points
     */
    vector < MatrixXf > _point_predicted;
    vector < MatrixXf > _point_reference;
    vector < MatrixXf > _point_measure;

    vector <float> _weight_mean;
    vector <float> _weight_cov;

    /*
     *  Functions (pointers)
     */
    void (*_propagateFunc)(const MatrixXf &, MatrixXf &);
    void (*_measurementFunc)(const MatrixXf &, MatrixXf &);

    /*!
     * \brief updateMean :
     * Update the -mean field, considering the current state of sigma points
     */
    void updateMean(const std::vector < MatrixXf > &sigma_points,
                    MatrixXf &mean);

    /*!
     * \brief updateCov :
     * Update the -cov matrix, considering the current state of sigma points
     */
    void updateCov(const vector < MatrixXf > &sigma_points_1, const MatrixXf &mean_1,
                   const vector < MatrixXf > &sigma_points_2, const MatrixXf &mean_2,
                   MatrixXf &cov);


  public :
    /*
     *  Constructors
     */
    SigmaPoints(const MatrixXf &_mean,
                const MatrixXf &_cov,
                float kappa);



    /*
     *  Methods
     */

    /*!
     * \brief computeSigmaPoints
     * Update the set of sigma points, knowing mean, cov, and dispersion parameters
     * // Unscented Transform //
     */
    void computeSigmaPoints();

    /*!
     * \brief getState
     * \param _mean
     * \param _cov
     */
    void getState(MatrixXf &mean,
                  MatrixXf &cov) const;

    /*!
     * \brief getMeasuredState
     * \param mean
     * \param cov
     * \param cov_cross
     */
    void getMeasuredState(MatrixXf &mean,
                          MatrixXf &cov_measure,
                          MatrixXf &cov_cross) const;

    /*!
     * \brief getPredictedState
     * \param mean
     * \param cov
     * \param cov_cross
     */
    void getPredictedState(MatrixXf &mean,
                           MatrixXf &cov,
                           MatrixXf &cov_cross) const;

    /*!
     * \brief measureSigmaPoints : project the sigma points on the measurement space
     * populates the _point_measurement_space vector
     */
    void measureSigmaPoints();

    /*!
     * \brief propagateSigmaPoints
     */
    void propagateSigmaPoints();

    /*!
     * \brief setState
     * \param _mean
     * \param _cov
     */
    void setState(const MatrixXf &mean,
                  const MatrixXf &cov);

    /*!
     * \brief setPropagationFunction :
     *  Define  the function used to propagate the sigma points
     *  The propagation function parameters apply on the state vector
     */
    void setPropagationFunction(void (*_prop_function)(const MatrixXf &, MatrixXf &));

    /*!
     * \brief setMeasurementFunction :
     *  Define  the function used to project the sigma points on measurement space
     *  The measurement function parameters apply on the state vector
     */
    void setMeasurementFunction(void (*meas_function)(const MatrixXf &, MatrixXf &));
};

#endif // SIGMAPOINT_H

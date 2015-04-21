#include "sigma_point.h"

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file sigma_point.cpp
 *  @brief Implements a class of sigma points, to be used with UKF
 *  @version 1.0
 *  @date 12-11-2013
 */

SigmaPoints::SigmaPoints(const MatrixXf &mean,
                         const MatrixXf &cov,
                         float kappa) {

  /*
   * Copy distribution parameters
   */

  if ( (mean.rows () != cov.cols ()) ||
       cov.cols () != cov.rows ()) {
    THROW_ERR("SPoint : Could not create sigma points, matrix sizes do not match");
  }

  _dim   = mean.rows ();

#ifdef DEBUG
  printf("Sigma points : dimensions %li x %li\n", mean.rows (), mean.cols ());
#endif

  _mean_reference   = mean;

  _cov_reference        = cov;
  _cov_measure          = cov; // Useless ?
  _cov_cross_pred_meas  = cov; // Useless ?
  _cov_cross_pred_state = cov; // Useless ?

  _kappa  = kappa;

  /*
   * Create sigma points distribution (Unscented Transform)
   */
  computeSigmaPoints ();

  _propagateFunc = NULL;
  _measurementFunc = NULL;
}


/*!
 * \brief SigmaPoint::computeSigmaPoints
 * Compute the distribution of sigma points from mean, covariance and
 *  distribution parameters (alpha, beta)
 */
void SigmaPoints::computeSigmaPoints () {

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
  LLT <MatrixXf> lltOfCov(_cov_reference);
  MatrixXf L = lltOfCov.matrixL ();

#ifdef DEBUG
  cout << "Square root matrix : \n" << L << endl;
#endif

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
 * \brief SigmaPoints::getState
 * \param _mean
 * \param _cov
 */
void SigmaPoints::getState(MatrixXf &mean,
                           MatrixXf &cov) const {
  mean  = _mean_reference;
  cov   = _cov_reference;
}


/*!
 * \brief SigmaPoints::getPredictedState
 * \param _mean
 * \param _cov
 */
void SigmaPoints::getPredictedState(MatrixXf &mean,
                                    MatrixXf &cov,
                                    MatrixXf &cov_cross) const {
  mean      = _mean_predicted;
  cov       = _cov_predicted;
  cov_cross = _cov_cross_pred_state;
}


/*!
 * \brief SigmaPoints::getMeasuredState
 * \param _mean
 * \param _cov
 */
void SigmaPoints::getMeasuredState(MatrixXf &mean,
                                   MatrixXf &cov_measure,
                                   MatrixXf &cov_cross) const {
  mean        = _mean_measure;
  cov_measure = _cov_measure;
  cov_cross   = _cov_cross_pred_meas;
}


/*!
 * \brief SigmaPoints::measureSigmaPoints
 *  Project the propagated sigma points onto the measurement space
 */
void SigmaPoints::measureSigmaPoints () {
  if (_measurementFunc == NULL) {
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

//#ifdef DEBUG
//  printf("Measure sigma_points : \n");
//  printEigenVector(_point_measure);

//  printf("Mean measure : \n");
//  printEigenMatrix(_mean_measure);
//#endif

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
 * \brief SigmaPoints::propagateSigmaPoints
 */
void SigmaPoints::propagateSigmaPoints() {
  if (_propagateFunc == NULL) {
    THROW_ERR("SPoint : Invalid propagation function");
  }

  // Propagate existing set of points (supposed to be representative)
  _point_predicted.resize(_point_reference.size());

  unsigned  int i=0;
  while (i < this->_point_reference.size()) {
    _propagateFunc(_point_reference[i], _point_predicted[i]);
    ++i;
  }

#ifdef DEBUG
  cout << "\nInitial mean :\n" << _mean_reference.block(0,0,6,1) << endl;
#endif

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

#ifdef DEBUG
  cout << "\nPredicted mean :\n" << _mean_predicted.block(0,0,6,1) << endl;
#endif
}

/*!
 * \brief SigmaPoints::setState
 * \param _mean
 * \param _cov
 */
void SigmaPoints::setState(const MatrixXf &mean,
                           const MatrixXf &cov) {
  _mean_reference = mean;
  _cov_reference  = cov;

  computeSigmaPoints ();
}


/*!
 * \brief SigmaPoints::setMeasurementFunction
 */
void SigmaPoints::setMeasurementFunction(void (*meas_function)(const MatrixXf &, MatrixXf &)) {
  _measurementFunc = meas_function;
}


/*!
 * \brief SigmaPoints::setPropagationFunction
 */
void SigmaPoints::setPropagationFunction(void (*_prop_function)(const MatrixXf &, MatrixXf &)) {
  _propagateFunc = _prop_function;
}


/*!
 * \brief SigmaPoint::updateMean
 *  Computes mean from a given sigma points set
 */
void SigmaPoints::updateMean(const vector < MatrixXf > &sigma_points,
                             MatrixXf &mean) {
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


/*!
 * \brief SigmaPoint::updateCov
 * Computes covariance between two given sigma points sets
 */
void SigmaPoints::updateCov(const vector < MatrixXf > &sigma_points_1,
                            const MatrixXf &mean_1,
                            const vector < MatrixXf > &sigma_points_2,
                            const MatrixXf &mean_2,
                            MatrixXf &cov) {

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

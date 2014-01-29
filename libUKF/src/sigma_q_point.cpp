#include "sigma_q_point.h"

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file sigma_q_point.cpp
 *  @brief Implements a class of sigma points, to be used with UKF. Points are quaternions in this case !
 *  @version 1.0
 *  @date 12-11-2013
 */

SigmaQPoints::SigmaQPoints(const Matrix<float, 3, 1> &angles_mean,
                           const Matrix<float, 3, 3> &angles_cov,
                           float kappa) {

  _dim = 3;
  _kappa  = kappa;

  // Convert angles to quaternion space
  eulerToQuaternion(angles_mean, _q_mean_reference);
  _mean_reference = angles_mean;
  _cov_reference = angles_cov;

  // Deal with fucking rdm number generator
  _st_tool = new statisticTools();

  // Reset a few matrices..
  _process_noise.setZero (3, 3);
}


SigmaQPoints::~SigmaQPoints () {
  delete _st_tool;

  // No need to handle Eigen matrices..
}

void SigmaQPoints::angleAxisNormalize(const AngleAxis <float> &input, AngleAxis <float> &output) {
  float norm = sqrt(pow(input.axis ()(0), 2.f) +
                    pow(input.axis ()(1), 2.f) +
                    pow(input.axis ()(2), 2.f));

  norm *= input.angle ();

  if (norm != 0.f) {
    output.axis () = input.axis () / norm;
    output.angle () = input.angle ();
  } else  {
    output.axis () = input.axis ();
    output.angle () = 0.f;
  }
}


void SigmaQPoints::averageQuaternionsSlerp(const vector <Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                           const vector <float> &q_weight,
                                           Quaternionf &q_avg) {

  /*
   * !!! TODO !!!
   *  Use the interpolation method on the error vector (see paper) to iteratively converge towards the
   *  Real quaternions mean !!
   *
   */

  // We use the slerp interpolation to get the average quaternion best describing
  // a set of quaternions
  if (q_list.size () != q_weight.size ()) {
    THROW_ERR("Error averaging quaternions, \n weights and quaternion lists have different size");
  }

  if (q_list.size () == 0) {
    printf("Empty list of quaternions\n");
    return;
  }

  float         current_weight = 0.f;
  unsigned int  n_avg = 0, iq;
  Quaternionf   q_avg_temp;

  q_avg = q_list[0];
  current_weight = q_weight[0];

  for (iq = 1; iq < q_list.size (); ++iq) {
    if (q_weight[iq] != 0.f) {
      q_avg_temp = q_avg.slerp( q_weight[iq] /
                                (q_weight[iq] + current_weight),
                                q_list[iq]);
      q_avg = q_avg_temp;

      current_weight += q_weight[iq];
      n_avg ++;
    }
  }
}

void SigmaQPoints::averageQuaternionsAngleAxis(const vector <Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                               const vector <float>       &q_weight,
                                               Quaternionf  &q_avg) {

  if (q_list.size () != q_weight.size ()) {
    THROW_ERR("Error averagin quaternions, weights and quaternions lists have different sizes !");
  }

  AngleAxis <float> vec_mean;

  vec_mean = AngleAxisf(q_list[0]);

  for (unsigned int i = 0; i< q_list.size (); ++i) {
    // TODO: Finish angle/axis implementation of quaternion pool averaging
  }

  //  // in
  //  quat_with_weight inputs[N];

  //  // result net torque as a rotation vector
  //  vec3 accum(0,0,0);

  //  foreach(q in inputs)
  //  {
  //      // add the torque contributation of this rotation
  //      accum += q.getAngle() * q.weight * q.getUnitAxis();
  //  }

  //  // convert angle/axis back to quaternion
  //  quat result(accum.magnitude(), accum.asUnitVector());

}

void SigmaQPoints::averageQuaternionsIterative(const vector <Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                               const vector <float>       &q_weight,
                                               int          n_iterations,
                                               float        max_err,
                                               Quaternionf  &q_avg) {

  /*!
   * Steps :
   * - start with an initial estimate (first quaternion of the list for example)
   *
   * - compute the mean error vector between all quaternions of the list,
   *  and the current mean
   *
   * - update the current mean quaternion with the mean error vector
   *
   * - stop when error vector is below a threshold, or iteration number is above max
   *
   * \remarks : Use the Angle Axis representation from Eigen
   **/

  // Method is from Xavier Pennec (INRIA)

  // See http://en.wikipedia.org/wiki/Generalized_quaternion_interpolation

  if (q_list.size () != q_weight.size ()) {
    THROW_ERR("Error averagin quaternions, weights and quaternions lists have different sizes !");
  }

  vector < AngleAxis<float> > err_aa_vector;
  AngleAxis<float>            err_aa_mean;
  Quaternionf                 err_q_mean;
  unsigned int i = 0, iq = 0;

  err_aa_vector.resize(q_list.size());

  err_aa_mean.angle() = M_PI; // Unused initial value, must just be above "max_err_norm"

  q_avg = q_list[0];

  while ( (i < n_iterations) &&
          (err_aa_mean.angle() >= max_err)) {

    /* Update error quaternions (in AngleAxis form) :
      * they are defined as the quaternions to go from
      * sigma points to current mean estimation
      */
    for (iq = 0; iq < q_list.size(); ++iq) {
      err_aa_vector[iq] = AngleAxisf(q_list[iq] * q_avg.conjugate ());
    }

    /*! Compute mean error quaternion (mean of the axis)
       * - Computation of the mean is initially in the AngleAxis space
       * - We finally update the quaternion mean once barycentric mean is obtained
       * from the angle axis
       */

    err_aa_mean.axis ().setZero (3,1);
    err_aa_mean.angle () = 0.f;

    for (iq= 0; iq < err_aa_vector.size(); ++iq) {
      err_aa_mean.axis() += err_aa_vector[iq].axis() *
          err_aa_vector[iq].angle() *
          _weight_mean[iq];
    }

    // Refactor an angle-axis representation
    float norm = sqrt(pow(err_aa_mean.axis ()(0), 2.f) +
                      pow(err_aa_mean.axis ()(1), 2.f) +
                      pow(err_aa_mean.axis ()(2), 2.f));

    if (norm != 0.f) {
      err_aa_mean.axis() /= norm;
      err_aa_mean.angle() = norm;
    } else {
      err_aa_mean.Identity ();
    }

    /* Switch back to the quaternion space */
    err_q_mean = Quaternionf(err_aa_mean);

    /* Update mean quaternion using mean error vector */
    q_avg = err_q_mean * q_avg;
    ++i;
  }
}



/*!
 * \brief Compute mean values from all the quaternion sigma points...
 * Not that simple knowing quaternions algebra, we'll use gradient descent
 * Covariance in quaternion space is another output of the calculus
 */

void SigmaQPoints::computeQMeanAndCovariance(int max_iterations,
                                             float max_err_norm) {

  unsigned int iq = 0;


  averageQuaternionsIterative(_q_sigma_pt,
                              _weight_mean,
                              max_iterations,
                              max_err_norm,
                              _q_mean_predicted);

  quaternionToEuler (_q_mean_predicted, _mean_predicted);

#ifdef DEBUG
  cout << "SigmaQPoints : predicted mean\n" << _mean_predicted << endl;
#endif

  /* --- Update covariance of the sigma points : ---
  * We use the last error vector to compute the covariance of the set.
  * classically, covariance is defined as :
  *
  * sig(x,y) = E( (x - E(x))(y - E(y))
  *
  * That's about the same for us as regards the sigma points set, we
  * use the last error vector as the best "esperance" estimation
  *
  */

  // Compute the covariance matrix in vector space
  Quaternionf q_mismatch;
  Vector3f    euler_mismatch;
  float       weight_overall = 0.f;
  _cov_predicted.setZero(3,3);

  for (iq = 0; iq < _q_sigma_pt.size(); ++iq) {
    // Computes the covariance using the diff between mean estimate and sigma points

//    // -------------------------------------------------
//    // Quaternion representing the mismatch between "mean" and this sigma point
//    q_mismatch =  _q_mean_predicted.inverse() * _q_sigma_pt[iq]; // FIXME !

//    // Switch to Euler (3D) representation
//    quaternionToEuler(q_mismatch, euler_mismatch);
//    // -------------------------------------------------

    // -------------------------------------------------
    // Quaternion representing the mismatch between "mean" and this sigma point
    q_mismatch = _q_sigma_pt[iq]; // FIXME !

    // Switch to Euler (3D) representation
    quaternionToEuler(q_mismatch, euler_mismatch);
    euler_mismatch -= _mean_predicted;
    // -------------------------------------------------

    weight_overall += _weight_cov[iq];
    _cov_predicted += _weight_cov[iq] * (euler_mismatch * euler_mismatch.transpose());

  }

  if (weight_overall != 0.f) {
    _cov_predicted /= weight_overall;
  }


#ifdef DEBUG_LINUX
  cout << "SigmaQPoints : Mismatch vector : \n"  << _mean_predicted - _mean_reference << endl;
  cout << "SigmaQPoints : Covariance : \n"  << _cov_predicted << endl;
#endif


}

void SigmaQPoints::computeSigmaQPoints () {
  _sigma_pt.resize (2*_dim+1);
  _q_sigma_pt.resize (2*_dim+1);
  _weight_mean.resize (2*_dim + 1);
  _weight_cov.resize (2*_dim + 1);

  // Mean
  _sigma_pt[0]  = _mean_reference;
  eulerToQuaternion (_mean_reference, _q_sigma_pt[0]);
  _weight_mean[0] = _kappa / (_dim + _kappa);
  _weight_cov[0]  = _weight_mean[0];

  /*
   * Compute "square root matrix" from covariance to get sigma points
   *
   * Rmq : we add the process noise BEFORE the computation of the sigma points,
   * this ensures that it is automatically taken into account
   * _process_noise should not be used elsewhere !
   */

  LLT <Matrix <float, 3,3> > lltOfCov(_cov_reference + _process_noise);
  Matrix <float, 3,3> L = lltOfCov.matrixL ();

  // Distributed points..
  for (int i=1; i<=_dim; ++i) {
    // Sigma point position using covariance square root
    _sigma_pt[i  ]       = _mean_reference +  L.col (i-1) * sqrt(_dim + _kappa);
    _sigma_pt[_dim + i]  = _mean_reference -  L.col (i-1) * sqrt(_dim + _kappa);

    // Equivalent in Quaternion space :
    eulerToQuaternion (_sigma_pt[i  ] , _q_sigma_pt[ i]);
    eulerToQuaternion (_sigma_pt[_dim + i] , _q_sigma_pt[_dim + i]);

    // Weights : same weights for everyone right now..
    _weight_mean[i]       = 1.f / (2.f * (_kappa + _dim));
    _weight_mean[_dim + i] = _weight_mean[i];

    _weight_cov[i]        = 1.f / (2.f * (_kappa + _dim));
    _weight_cov[_dim + i]  = _weight_cov[i];
  }
}

float SigmaQPoints::distQuaternions(const Quaternionf &q_1,
                                    const Quaternionf &q_2) {

  return AngleAxis<float> (q_1 * q_2.conjugate ()).angle();
}


void SigmaQPoints::eulerToQuaternion(const Vector3f &input,
                                     Quaternion <float> &q_out) const {

  /*
   * Conversion from Euler angles to Quaternion
   *
   * see : http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
   */

  // Create rotation matrix, then quaternion from rot matrix ?
  //  Matrix <float, 3,3> rot_mat;

  /*
   \cos\theta \cos\psi  & -\cos\phi \sin\psi + \sin\phi \sin\theta \cos\psi   & \sin\phi \sin\psi + \cos\phi \sin\theta \cos\psi
   \cos\theta \sin\psi  & \cos\phi \cos\psi + \sin\phi \sin\theta \sin\psi    & -\sin\phi \cos\psi + \cos\phi \sin\theta \sin\psi
   -\sin\theta          & \sin\phi \cos\theta                                 & \cos\phi \cos\theta
  */

  // order {phi, theta, psi}

  float phi = input(0), theta = input(1), psi = input(2);

  /*
  rot_mat(0,0) = cos(input(1)) * cos(input(2));
  rot_mat(0,1) = - cos(input(0)) * sin(input(2)) + sin(phi) * sin(theta) * cos(psi);
  rot_mat(0,2) = sin(phi) * sin(psi) + cos(phi) * sin(theta) * cos(psi);

  rot_mat(1,0) = cos(theta) * sin (psi);
  rot_mat(1,1) = cos(phi) * cos(psi) + sin(phi) * sin(theta) * sin(psi);
  rot_mat(1,2) = -sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi);

  rot_mat(2,0) = -sin(theta);
  rot_mat(2,1) = sin(phi) * cos(theta);
  rot_mat(2,2) = cos(phi) * cos(theta);
  */

  /*
  rot_mat = AngleAxisf(input(0,0), Vector3f::UnitZ())
      * AngleAxisf(input(1,0), Vector3f::UnitY())
      * AngleAxisf(input(2,0), Vector3f::UnitZ());
  */

  // Quaternion from rot matrix :
  //  q_out = Quaternion <float> (rot_mat.transpose());


  q_out.w () = cos(phi/2.f) * cos(theta/2.f) * cos(psi/2.f)
      + sin(phi/2.f) * sin(theta/2.f) * sin(psi/2.f);

  q_out.x () = sin(phi/2.f) * cos(theta/2.f) * cos(psi/2.f)
      - cos(phi/2.f) * sin(theta/2.f) * sin(psi/2.f);

  q_out.y () = cos(phi/2.f) * sin(theta/2.f) * cos(psi/2.f)
      + sin(phi/2.f) * cos(theta/2.f) * sin(psi/2.f);

  q_out.z () = cos(phi/2.f) * cos(theta/2.f) * sin(psi/2.f)
      - sin(phi/2.f) * sin(theta/2.f) * cos(psi/2.f);

}
void SigmaQPoints::getStatePre(Matrix<float, 3,1> &mean,
                               Matrix<float, 3,3> &cov) const {

  mean = _mean_predicted;
  cov = _cov_predicted;
}


void SigmaQPoints::pickMeanQuaternion(const vector<Quaternionf, aligned_allocator<Quaternionf> > &q_list,
                                      const vector<float> &q_weight,
                                      Quaternionf *mean_quaternion) {

  if (q_list.size () != q_weight.size ()) {
    THROW_ERR("Error averaging quaternions, \n weights and quaternion lists have different size");
  }

  float dist_min = 10e100, dist = 0;

  unsigned int i,j, bi = 0;

  for (i = 0; i<q_list.size (); ++i) {
    dist = 0.f;
    for (j = 0; j<q_list.size (); ++j) {
      dist += distQuaternions (q_list[i], q_list[j]) * q_weight[j];
    }

    if (dist < dist_min) {
      dist_min = dist;
      bi = i;
    }
  }

  *mean_quaternion = q_list[bi];
}

/*!
 * \brief propagate quaternions
 */
void SigmaQPoints::propagateSigmaQPoints () {
  computeSigmaQPoints();

  vector<Quaternionf, aligned_allocator<Quaternionf> > q_points = _q_sigma_pt;

  for (unsigned int i=0; i<_q_sigma_pt.size(); ++i) {
    _propagateFunc(q_points[i], _q_sigma_pt[i]);
  }
}

void SigmaQPoints::quaternionToEuler(const Quaternion <float> &q_input,
                                     Vector3f &output) const {

  /*
   * Conversion from Quaternions to Euler angles
   *
   * see : http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
   *
   * and
   *  " Euler Angles, Quaternions, and Transformation Matrices" - D. Henderson - NASA ;-)
   */
  Vector4f q;
  q = q_input.coeffs (); // /!\ Coeffs are stored in the order [x, y, z, w] and not [w, x, y,z] !!

  output(0) = atan2(2.f * ( q(3) * q(0) + q(1) * q(2)),
                    (1.f - 2.f * (q(0) * q(0)
                                  + q(1) * q(1))));

  output(1) = asin( 2.f * (q(3) * q(1) - q(2) * q(0)));

  output(2) = atan2( 2.f * (q(3) * q(2) + q(0) * q(1)),
                     1.f - 2.f * (q(1) * q(1) + q(2) * q(2)));
}

/*!
 * \brief SigmaPoints::setMeasurementFunction
 */
void SigmaQPoints::setMeasurementFunction(void (*meas_function)(const Quaternionf &, Quaternionf &)) {
  _measurementFunc = meas_function;
}

void SigmaQPoints::setProcessNoise(const MatrixXf &process_noise) {
  printf("Quaternions : set process noise %d x %d\n",
         process_noise.rows (),
         process_noise.cols ());

  _process_noise = process_noise;
}


/*!
 * \brief SigmaPoints::setPropagationFunction
 */
void SigmaQPoints::setPropagationFunction(void (*_prop_function)(const Quaternionf &, Quaternionf &)) {
  _propagateFunc = _prop_function;
}


/*!
 * \brief SigmaQPoints::setState
 * \param _mean
 * \param _cov
 */
void SigmaQPoints::setState(const MatrixXf &mean,
                            const MatrixXf &cov) {

  _mean_reference = mean;
  _cov_reference  = cov;
}

/*!
 * \brief Project sigma points onto the measurement space
 */
void SigmaQPoints::projectSigmaQPoints() {
  _sigma_pt.clear();
  unsigned int i = 0;

  // Get angles back from the quaternions and get the tail values for sigma_meas_points
  i = 0;

  while (i < this->_q_sigma_pt.size()) {
    quaternionToEuler(_q_sigma_pt[i], _sigma_pt[i]);
    ++i;
  }

}

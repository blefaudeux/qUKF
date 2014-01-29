#include "statistic_tools.h"

/*
 *  @license GPL
 *  @author Benjamin Lefaudeux (blefaudeux at github)
 *  @file statistic_tools.cpp
 *  @brief Some tools to be used here and there..
 *  @version 1.0
 *  @date 12-11-2013
 */

statisticTools::statisticTools() {
  // Initialize random number generator
  rng = new boost::mt19937;
  rand_gen = new boost::random::normal_distribution<> (0.0, 1.0);
}

/*!
 * \brief randN : generate random number from a normal distribution
 * std is 1.f
 */
float statisticTools::normal() {
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > rand_var(*rng, *rand_gen);

  return rand_var();
}


void statisticTools::multivariate_normal( MatrixXf &covar_mat,
                                          VectorXf &mean_vec,
                                          VectorXf &samples) {

  // Initialize random number generator
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > rand_var(*rng, *rand_gen);

  // TODO : check that random number are indeed random from time to time...


  // Allocate matrices and compute weighted samples
  int size = samples.size();

  SelfAdjointEigenSolver<Eigen::MatrixXf> eigenSolver(covar_mat);
  MatrixXf rot = eigenSolver.eigenvectors();
  VectorXf scl = eigenSolver.eigenvalues();
  VectorXf samples_temp(size);

  for (int ii=0; ii<size; ++ii) {
    scl(ii,0) = sqrt(scl(ii,0));
  }

  for (int ii=0;ii<size;++ii) {
    samples_temp(ii,0) = rand_var()*scl(ii,0);
  }

  samples = mean_vec + rot * samples_temp;
}



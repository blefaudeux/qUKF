#ifndef STATISTICTOOLS_H
#define STATISTICTOOLS_H

/*!
 *  A few tools to compute or generate statistically-relevant data..
 * \author : Benjamin Lefaudeux
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ostream>
#include <vector>

#include <Eigen/Eigen>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace Eigen;

class statisticTools {
  public :
    statisticTools();

    /*!
     * \brief normal : normal distribution, zero mean, unit variance
     * \return random variable
     */
    float normal();


    /*!
     * \brief multivariate_normal
     * \param covar_mat : covariance for the distribution
     * \param mean_mat  : mean values for the distribution
     * \param samples   : filled with random variables from multivariate normal distribution
     *
     * \warning : covar_mat and mean_mat must have appropriate sizes !
     */
    void multivariate_normal(MatrixXf &covar_mat,
                             VectorXf &mean_vec,
                             VectorXf &samples);

  private :
    boost::mt19937 *rng;
    boost::random::normal_distribution<> *rand_gen;
};


#endif // STATISTICTOOLS_H

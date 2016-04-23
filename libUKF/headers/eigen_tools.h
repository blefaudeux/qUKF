#ifndef EIGEN_TOOLS_H
#define EIGEN_TOOLS_H

#include <Eigen/Eigen>
#include <vector>
#include <stdio.h>
#include "def.h"
#include <iostream>
#include <Eigen/StdVector>
#include <math.h>

using namespace std;

namespace qukf {
    template <typename T, size_t Rows>
    using Vec = Eigen::Matrix<T, Rows, 1>;

    template <typename T>
    using Vec2 = Eigen::Matrix<T, 2, 1>;

    template <typename T>
    using Vec3 = Eigen::Matrix<T, 3, 1>;

    template <typename T, size_t Rows, size_t Cols>
    using Mat = Eigen::Matrix<T, Rows, Cols>;

    template <typename T>
    using Mat3 = Eigen::Matrix<T, 3, 3>;

    template <typename T, size_t Rows_Cols>
    using MatSquare = Eigen::Matrix<T, Rows_Cols, Rows_Cols>;

    template <typename T>
    using MatX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    template <typename T, size_t Rows, size_t Cols>
    using Mat = Eigen::Matrix<T, Rows, Cols>;

    using MatXf = Eigen::MatrixXf;

    float pseudo_inv(const Eigen::MatrixXf *mat_in,
                     Eigen::MatrixXf *mat_out);

    float pseudo_inv_6(const Eigen::Matrix <float, 6,6> *mat_in,
                    Eigen::Matrix <float, 6,6> *mat_out);


    float pseudo_inv_12(const Eigen::Matrix <float, 12, 12> *mat_in,
                    Eigen::Matrix <float, 12, 12> *mat_out);

}

#endif // EIGEN_TOOLS_H

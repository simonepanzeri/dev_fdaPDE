//
// Created by simonepanzeri on 24/11/2021.
//

#ifndef DEV_FDAPDE_FDAPDE_H
#define DEV_FDAPDE_FDAPDE_H

#include <iostream>

#include <cstdlib>
#include <cmath>
#include <numeric>
#include <limits>
#include <string>
#include <type_traits>
#include <iomanip>
#include <utility>
#include <memory>
#include <stack>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <unordered_set>
#include <tuple>


#include "eigen_library/Eigen/Sparse"
//#include "eigen_library/Eigen/SparseLU"
#include "eigen_library/Eigen/Dense"
#include "eigen_library/Eigen/StdVector"
#include "eigen_library/Eigen/IterativeLinearSolvers"

typedef double Real;
typedef int UInt;

typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<UInt,Eigen::Dynamic,Eigen::Dynamic> MatrixXi;
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> VectorXr;
typedef Eigen::Matrix<UInt,Eigen::Dynamic,1> VectorXi;
typedef Eigen::Matrix<VectorXr,Eigen::Dynamic,Eigen::Dynamic> MatrixXv;
typedef Eigen::SparseMatrix<Real> SpMat;
typedef Eigen::SparseVector<Real> SpVec;
typedef Eigen::Triplet<Real> coeff;

template <bool ... b>
struct multi_bool_type
{};

typedef multi_bool_type<true> t_type;
typedef multi_bool_type<false> f_type;
typedef multi_bool_type<true, true> tt_type;
typedef multi_bool_type<false, true> ft_type;
typedef multi_bool_type<true, false> tf_type;
typedef multi_bool_type<false, false> ff_type;

#endif //DEV_FDAPDE_FDAPDE_H

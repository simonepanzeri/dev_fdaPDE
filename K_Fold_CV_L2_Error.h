//
// Created by simon on 01/12/2021.
//

#ifndef DEV_FDAPDE_K_FOLD_CV_L2_ERROR_H
#define DEV_FDAPDE_K_FOLD_CV_L2_ERROR_H

#include "FdaPDE.h"

// This file implements the cross validation error based on the L2 norm useful for the Density Estimation problem

//! @brief A class to compute the L2 error during cross-validation.
template<UInt ORDER, UInt mydim, UInt ndim>
class KfoldCV_L2_error{
private:
    // A member to acess  data problem methods
    const DataProblem<ORDER, mydim, ndim>& dataProblem_;

public:
    //! A constructor.
    KfoldCV_L2_error(const DataProblem<ORDER, mydim, ndim>& dp): dataProblem_(dp) {};
    //! A call operator to compute the L2 error.
    Real operator()(const SpMat& Psi, const VectorXr& g) {
        Real integral = dataProblem_.FEintegrate_exponential(2.*g);
        Real test = (Psi*g).array().exp().sum();

        return (integral - 2./Psi.rows() *test);
    }

};

#endif //DEV_FDAPDE_K_FOLD_CV_L2_ERROR_H

//
// Created by simonepanzeri on 30/11/2021.
//

#ifndef DEV_FDAPDE_FUNCTIONAL_PROBLEM_IMP_H
#define DEV_FDAPDE_FUNCTIONAL_PROBLEM_IMP_H

#include "Kronecker_Product.h"

template<UInt ORDER, UInt mydim, UInt ndim>
std::pair<Real, VectorXr>
FunctionalProblem<ORDER, mydim, ndim>::computeIntegrals(const VectorXr& g) const{

    using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1>>;

    // Initialization
    Real int1 = 0.;
    VectorXr int2 = VectorXr::Zero(dataProblem_.getNumNodes());
    for(UInt triangle = 0; triangle < dataProblem_.getNumElements(); triangle++){

        Element<EL_NNODES, mydim, ndim> tri_activated = dataProblem_.getElement(triangle);
// (1) -------------------------------------------------
        Eigen::Matrix<Real,EL_NNODES,1> sub_g;
        for (UInt i = 0; i < EL_NNODES; i++){
            sub_g[i] = g[tri_activated[i].getId()];
        }
// (2) -------------------------------------------------
        Eigen::Matrix<Real, Integrator::NNODES, 1> expg = (dataProblem_.getPsiQuad()*sub_g).array().exp();

        Eigen::Matrix<Real, EL_NNODES, 1> sub_int2;

        int1+=expg.dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0])) * tri_activated.getMeasure();
        sub_int2 = dataProblem_.getPsiQuad().transpose() * expg.cwiseProduct(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0])) * tri_activated.getMeasure();

        for (UInt i = 0; i < EL_NNODES; i++){
            int2[tri_activated[i].getId()] += sub_int2[i];
        }
    }

    return std::pair<Real, VectorXr> (int1, int2);
}

template<UInt ORDER, UInt mydim, UInt ndim>
std::tuple<Real, VectorXr, Real, Real>
FunctionalProblem<ORDER, mydim, ndim>::computeFunctional_g(const VectorXr& g, Real lambda, const SpMat& Psi) const{

    Real int1;
    VectorXr int2;
    std::tie(int1,int2) = computeIntegrals(g);

    const UInt n = Psi.rows();
    const Real llik = -(Psi*g).sum() + n*int1;
    const Real pen = g.dot(dataProblem_.getP()*g);

    VectorXr grad1 = - VectorXr::Constant(n,1).transpose()*Psi;
    VectorXr grad2 =  n*int2;
    VectorXr grad3 = 2*g.transpose()*dataProblem_.getP();

    VectorXr grad = grad1 + grad2 + lambda*grad3;

    return std::make_tuple(llik+lambda*pen, grad, llik, pen);

}


template<UInt ORDER, UInt mydim, UInt ndim>
std::pair<Real, Real>
FunctionalProblem<ORDER, mydim, ndim>::computeLlikPen_f(const VectorXr& f) const{

    Real llik = - (dataProblem_.getGlobalPsi()*f).array().log().sum() + dataProblem_.dataSize()*dataProblem_.FEintegrate(f);
    VectorXr tmp = f.array().log();
    Real pen = tmp.dot(dataProblem_.getP()*tmp);

    return std::pair<Real, Real>(llik,pen);
}

template<UInt ORDER, UInt mydim, UInt ndim>
std::pair<Real, VectorXr>
FunctionalProblem_time<ORDER, mydim, ndim>::computeIntegrals(const VectorXr& g) const{
    // Kronecker product of the Gauss quadrature rules weights
    VectorXr weights_kronecker;
    weights_kronecker.resize(Integrator::NNODES*IntegratorP5::NNODES);
    UInt k=0;
    for (UInt i = 0;  i < Integrator::NNODES; i++) {
        for (UInt j = 0;  j < IntegratorP5::NNODES; j++){
            weights_kronecker(k) = Integrator::WEIGHTS[i]*IntegratorP5::WEIGHTS[j];
            ++k;
        }
    }
    weights_kronecker.transpose();
    // Initialization
    Real int1 = 0.;
    VectorXr int2 = VectorXr::Zero(dataProblem_time_.getNumNodes()*dataProblem_time_.getSplineNumber());
    const MatrixXr& PsiQuad = dataProblem_time_.getPsiQuad();
    UInt global_idx = 0; //index that keeps track of the first B-spline basis function active in the current time-interval
    for (int time_step = 0; time_step < dataProblem_time_.getNumNodes_time()-1; ++time_step) {
        MatrixXr PhiQuad = dataProblem_time_.fillPhiQuad(time_step);
        MatrixXr Psi_kronecker_Phi = kroneckerProduct_Matrix(PsiQuad,PhiQuad);
        for(UInt triangle = 0; triangle < dataProblem_time_.getNumElements(); triangle++) {
            Element<EL_NNODES, mydim, ndim> tri_activated = dataProblem_time_.getElement(triangle);
// (1) -------------------------------------------------
            VectorXr sub_g;
            UInt k=0; //index for sub_g
            for (UInt i = 0; i < PsiQuad.cols(); ++i) {
                UInt global_location = tri_activated[i].getId()*dataProblem_time_.getSplineNumber();
                for (UInt j = global_idx; j < global_idx + PhiQuad.cols(); ++j) {
                    sub_g.resize(k+1);
                    sub_g(k)=g[global_location + j];
                    ++k;
                }
            }
            VectorXr expg = (Psi_kronecker_Phi*sub_g).array().exp();
            int1 += expg.dot(weights_kronecker) * tri_activated.getMeasure() * (dataProblem_time_.getMesh_time()[time_step+1]-dataProblem_time_.getMesh_time()[time_step])/2;
// (2) -------------------------------------------------
            VectorXr sub_int2;
            sub_int2 = Psi_kronecker_Phi.transpose() *
                       expg.cwiseProduct(weights_kronecker) * tri_activated.getMeasure() * (dataProblem_time_.getMesh_time()[time_step+1]-dataProblem_time_.getMesh_time()[time_step])/2;
            k=0;
            for (UInt i = 0; i < PsiQuad.cols(); ++i) {
                UInt global_location = tri_activated[i].getId()*dataProblem_time_.getSplineNumber();
                for (UInt j = global_idx; j < global_idx + PhiQuad.cols(); ++j) {
                    int2[global_location + j] += sub_int2[k];
                    ++k;
                }
            }
        }
        if(time_step>=dataProblem_time_.getSplineDegree()) ++global_idx;
    }
    return std::pair<Real, VectorXr> (int1, int2);
}

template<UInt ORDER, UInt mydim, UInt ndim>
std::tuple<Real, VectorXr, Real, Real, Real>
FunctionalProblem_time<ORDER, mydim, ndim>::computeFunctional_g(const VectorXr& g, Real lambda_S, Real lambda_T,
                                                                const SpMat& Upsilon) const {
    Real int1 = 0;
    VectorXr int2;
    std::tie(int1,int2) = computeIntegrals(g);

    const UInt n = Upsilon.rows(); // dataProblem_time_.dataSize()
    const Real llik = -(Upsilon*g).sum() + n * int1;

    const SpMat K1 = dataProblem_time_.getPen_s();
    //const SpMat K1 = kroneckerProduct(dataProblem_time_.getP().sparseView(), dataProblem_time_.getTimeMass());
    const SpMat K2 = dataProblem_time_.getPen_t();
    //const SpMat K2 = kroneckerProduct(dataProblem_time_.getMass(), dataProblem_time_.getPt());

    const Real pen_S = g.dot(K1 * g);
    const Real pen_T = g.dot(K2 * g);

    VectorXr grad1 = - VectorXr::Constant(n,1).transpose()*Upsilon;
    VectorXr grad2 = n * int2;
    VectorXr grad3_S = 2*g.transpose() * K1;
    VectorXr grad3_T = 2*g.transpose() * K2;

    VectorXr grad = grad1 + grad2 + lambda_S * grad3_S + lambda_T * grad3_T;

    return std::make_tuple(llik + lambda_S * pen_S + lambda_T * pen_T, grad, llik, pen_S, pen_T);
}

template<UInt ORDER, UInt mydim, UInt ndim>
std::tuple<Real, Real, Real>
FunctionalProblem_time<ORDER, mydim, ndim>::computeLlikPen_f(const VectorXr& f) const {
    Real llik = (dataProblem_time_.getUpsilon()*f).array().log().sum() + dataProblem_time_.dataSize() * dataProblem_time_.FEintegrate(f);
    VectorXr tmp = f.array().log();
    const SpMat K1 = dataProblem_time_.getPen_s();
    //const SpMat K1 = kroneckerProduct(dataProblem_time_.getP().sparseView(), dataProblem_time_.getTimeMass());
    const SpMat K2 = dataProblem_time_.getPen_t();
    //const SpMat K2 = kroneckerProduct(dataProblem_time_.getMass(), dataProblem_time_.getPt());
    Real pen_S = tmp.dot(K1 * tmp);
    Real pen_T = tmp.dot(K2 * tmp);

    return std::make_tuple(llik, pen_S, pen_T);
}

#endif //DEV_FDAPDE_FUNCTIONAL_PROBLEM_IMP_H

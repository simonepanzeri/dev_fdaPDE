//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_DATA_PROBLEM_H
#define DEV_FDAPDE_DATA_PROBLEM_H

#include "FdaPDE.h"
#include "DE_Data.h"
#include "Kronecker_Product.h"
#include "Spline.h"
#include "Projection.h"
#include "Finite_Element.h"
#include "Matrix_Assembler.h"
#include "Integration.h"

// This file contains data information for the Density Estimation problem

//! @brief A class to store common data for the problem.
template<UInt ORDER, UInt mydim, UInt ndim>
class DataProblem{
protected:
    using Integrator = typename DensityIntegratorHelper::Integrator<mydim>;
    static constexpr UInt EL_NNODES = how_many_nodes(ORDER,mydim);
    DEData<ndim> deData_;
    MeshHandler<ORDER, mydim, ndim> mesh_;
    SpMat R0_, R1_, GlobalPsi_;
    MatrixXr P_;
    Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES> PsiQuad_;

    //! A method to compute the finite element matrices.
    void fillFEMatrices();
    //! A method to compute the matrix which evaluates the basis function at the quadrature EL_NNODES.
    void fillPsiQuad();
    //! A method to compute the GlobalPsi (evaluation on the data nodes)
    virtual void fillGlobalPsi(); //!***

public:
    //! A constructor: it delegates DEData and MeshHandler constructors.
    /* THE FOLLOWING INPUT PARAMETERS
     * const RNumericMatrix& points, const RIntegerMatrix& sides,
     * const RIntegerMatrix& elements, const RIntegerMatrix& neighbors
     * ARE USED, INSTEAD OF SEXP Rmesh, SO THAT THEY CAN BE USED IN THE CONSTRUCTOR OF MeshHandler in Mesh.h
     */
    DataProblem(const std::vector<Point<ndim>>& data, const UInt& order,
                const VectorXr& fvec, Real heatStep, UInt heatIter,
                const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print,
                UInt search, const RNumericMatrix& points, const RIntegerMatrix& sides,
                const RIntegerMatrix& elements, const RIntegerMatrix& neighbors, bool isTime = 0);

    //! A method to compute the integral of a function.
    Real FEintegrate(const VectorXr& f) const {return (R0_*f).sum();}
    //! A method to compute the integral of the square of a function.
    Real FEintegrate_square(const VectorXr& f) const {return f.dot(R0_*f);}
    //! A method to compute the integral of the exponential of a function.
    Real FEintegrate_exponential(const VectorXr& g) const;
    //! A method to compute the matrix which evaluates the basis function at the data points.
    SpMat computePsi(const std::vector<UInt>& indices) const;

    // Getters
    //! A method to access the data. It calls the same method of DEData class.
    const std::vector<Point<ndim> >& data() const {return deData_.data();}
    //! A method returning a datum. It calls the same method of DEData class.
    const Point<ndim>& data(UInt i) const {return deData_.data(i);}

    //! A method returning the number of observations. It calls the same method of DEData class.
    UInt dataSize() const {return deData_.dataSize();}
    //! A method returning the the input order. It calls the same method of DEData class.
    UInt getOrder() const {return deData_.getOrder();}
    //! A method returning the initial coefficients for the density. It calls the same method of DEData class.
    VectorXr getFvec() const {return deData_.getFvec();}
    //! A method returning a bool which says if there is a user's initial density. It calls the same method of DEData class.
    bool isFvecEmpty() const {return deData_.isFvecEmpty();}
    //! A method returning the heat diffusion process alpha parameter. It calls the same method of DEData class.
    Real getHeatStep() const {return deData_.getHeatStep();}
    //! A method returning the number of iterations for the heat diffusion process. It calls the same method of DEData class.
    UInt getHeatIter() const {return deData_.getHeatIter();}
    //! A method returning the penalization parameters. It calls the same method of DEData class.
    Real getLambda(UInt i) const {return deData_.getLambda(i);}
    //! A method returning the number of lambdas. It calls the same method of DEData class.
    UInt getNlambda()  const {return deData_.getNlambda();}
    //! A method returning the number of folds for CV. It calls the same method of DEData class.
    UInt getNfolds()  const {return deData_.getNfolds();}
    //! A method returning the number of iterations to use in the optimization algorithm. It calls the same method of DEData class.
    UInt getNsimulations() const {return deData_.getNsimulations();}
    //! A method returning the number of parameters for fixed step methods. It calls the same method of DEData class.
    UInt getNstepProposals() const {return deData_.getNstepProposals();}
    //! A method returning a parameter for fixed step methods. It calls the same method of DEData class.
    Real getStepProposals(UInt i) const {return deData_.getStepProposals(i);}
    //! A method returning the tolerance for optimization algorithm first termination criteria. It calls the same method of DEData class.
    Real getTol1() const {return deData_.getTol1();}
    //! A method returning the tolerance for optimization algorithm second termination criteria. It calls the same method of DEData class.
    Real getTol2() const {return deData_.getTol2();}
    //! A method returning the boolean print member. It calls the same method of DEData class.
    bool Print() const {return deData_.Print();}
    //! A method returning the integer that specifies the search algorithm type.
    UInt getSearch() const {return deData_.getSearch();}

    //getter for mesh
    //! A method returning the mesh.
    const MeshHandler<ORDER, mydim, ndim>& getMesh() const {return mesh_;}
    //getter for specific mesh features
    //! A method returning the number of mesh EL_NNODES. It calls the same method of MeshHandler class.
    UInt getNumNodes() const {return mesh_.num_nodes();}
    //! A method returning the number of mesh elements. It calls the same method of MeshHandler class.
    UInt getNumElements() const {return mesh_.num_elements();}
    //! A method returning a node. It calls the same method of MeshHandler class.
    Point<ndim> getPoint(Id id) const {return mesh_.getPoint(id);}
    //! A method returning an element. It calls the same method of MeshHandler class.
    Element<EL_NNODES,mydim,ndim> getElement(Id id) const {return mesh_.getElement(id);}
    //! A method returning the element in which the point in input is located. It calls the same method of MeshHandler class.
    Element<EL_NNODES,mydim,ndim> findLocation(const Point<ndim>& point) const {return mesh_.findLocation(point);}

    //getter for matrices
    //! A method returning the P matrix.
    MatrixXr getP() const {return P_;}
    //! A method returning the PsiQuad_ matrix.
    const Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES>& getPsiQuad() const {return PsiQuad_;}
    //! A method returning the GlobalPsi_ matrix.
    SpMat getGlobalPsi() const {return GlobalPsi_;}
};

//######################################################################################################################
//############################### CHILD CLASS CREATED TO DEAL WITH DE_SPACE_TIME PROBLEM ###############################
//############## contains time mesh, data, and all the numerical matrices needed for the implementation ################
//######################################################################################################################

template<UInt ORDER, UInt mydim, UInt ndim>
class DataProblem_time : public DataProblem<ORDER,mydim,ndim>{
private:
    static const UInt SPLINE_DEGREE=3;
    static const UInt ORDER_DERIVATIVE=2;
    using IntegratorP5 = IntegratorGaussP5;
    Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline;
    //! A DEData_time object to store time data.
    DEData_time<ndim> deData_time_;
    //! A vector containing the time mesh.
    std::vector<Real> mesh_time_;
    //! Matrix of the evaluations of the spline basis functions in the time locations and mass matrix.
    SpMat GlobalPhi_, K0_;
    //! Time penalty matrix.
    SpMat Pt_;
    //! Kronecker product between GlobalPsi_ and GlobalPhi_.
    SpMat Upsilon_;

    void fillGlobalPsi(void) override;
    void fillGlobalPhi(void);
    void fillTimeMass(void);
    void fillTimeSecondDerivative(void);

public:
    //! A constructor:
    DataProblem_time(const std::vector<Point<ndim>>& data, const std::vector<Real>& data_time,
                     const UInt& order, const VectorXr& fvec, Real heatStep, UInt heatIter,
                     const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                     const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print,
                     UInt search, const RNumericMatrix& points, const RIntegerMatrix& sides,
                     const RIntegerMatrix& elements, const RIntegerMatrix& neighbors,
                     const std::vector<Real>& mesh_time, bool isTime = 1);

    //! A method filling the current PhiQuad matrix needed for the discretization of the exponential integral.
    Eigen::Matrix<Real, IntegratorP5::NNODES, SPLINE_DEGREE+1> fillPhiQuad(UInt time_node) const {
        Eigen::Matrix<Real, IntegratorP5::NNODES, SPLINE_DEGREE+1> phi;
        Real t_a = mesh_time_[time_node], t_b = mesh_time_[time_node+1];
        std::array<Real,IntegratorP5::NNODES> ref_nodes;
        for(UInt k = 0; k < IntegratorP5::NNODES; ++k)
            ref_nodes[k]=((t_b-t_a)*IntegratorP5::NODES[k]+t_a+t_b)/2;
        unsigned n_basis = time_node;
        for(UInt j = 0; j < SPLINE_DEGREE+1; j++){
            for(UInt i = 0; i < IntegratorP5::NNODES; i++){
                if(n_basis >= 0){
                    phi(i,j) = spline.BasisFunction(time_node+j, ref_nodes[i]);
                }
                else
                    phi(i,j) = 0.;
            }
            --n_basis;
        }
        return phi;
    }

    //! A method computing the Upsilon matrix (kronecker product between GlobalPsi_ and GlobalPhi_, calculated
    //! considering for each spatial location only the proper rows of GlobalPhi_ corresponding to the time instants
    //! when that location is observed).
    SpMat computeUpsilon(const SpMat& psi, const SpMat& phi, const std::map<UInt, std::set<UInt>>& data_noD) const;

    // Getters
    //! A method to access the data. It calls the same method of DEData class.
    const std::vector<Real>& data_time() const {return deData_time_.data();}
    //! A method returning a datum.
    Real data_time(UInt i) const {return data_time()[i];}
    //! A method returning the time penalty matrix.
    const SpMat & getPt(void) const {return Pt_;}
    //! A method returning the matrix of the evaluations of the spline basis functions.
    const SpMat & getGlobalPhi(void) const {return GlobalPhi_;}
    //! A method returning the time mass matrix.
    const SpMat & getTimeMass(void) const {return K0_;}
    //! A method returning the spline degree.
    UInt getSplineDegree(void) const {return SPLINE_DEGREE; }

    //getter for time mesh
    //! A method returning the time mesh.
    const std::vector<Real>& getMesh_time() const {return mesh_time_;}
    //! A method returning the number of time mesh nodes.
    UInt getNumNodes_time() const {return mesh_time_.size();}
    //! *******
};

#include "Data_Problem_imp.h"

#endif //DEV_FDAPDE_DATA_PROBLEM_H

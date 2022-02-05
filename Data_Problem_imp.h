//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_DATA_PROBLEM_IMP_H
#define DEV_FDAPDE_DATA_PROBLEM_IMP_H

#include "Integration.h"
#include "Matrix_Assembler.h"

#include <fstream>

template<UInt ORDER, UInt mydim, UInt ndim>
DataProblem<ORDER, mydim, ndim>::DataProblem(const std::vector<Point<ndim>>& data, const UInt& order,
                                             const VectorXr& fvec, Real heatStep, UInt heatIter,
                                             const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                                             const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print,
                                             UInt search, const RNumericMatrix& points, const RIntegerMatrix& sides,
                                             const RIntegerMatrix& elements, const RIntegerMatrix& neighbors, bool isTime) :
    deData_(data, order, fvec, heatStep, heatIter, lambda, nfolds, nsim, stepProposals, tol1, tol2, print, search),
    mesh_(points, sides, elements, neighbors, search) {

    std::vector<Point<ndim>>& data_ = deData_.data();

    // PROJECTION
    if(mydim == 2 && ndim == 3){
        std::cout << "##### DATA PROJECTION #####\n";

        projection<ORDER, mydim, ndim> projection(mesh_, data_);
        data_ = projection.computeProjection();
    }

    // REMOVE POINTS NOT IN THE DOMAIN
    if(!isTime) {
        for (auto it = data_.begin(); it != data_.end();) {
            Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(data_[it - data_.begin()]);
            if (tri_activated.getId() == Identifier::NVAL) {
                it = data_.erase(it);
                std::cout << "WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n";
            } else {
                it++;
            }
            std::cout << it - data_.begin() << std::endl;
        }
    }

    // FILL MATRICES
    fillFEMatrices();
    fillPsiQuad();

    if(!isTime)
        fillGlobalPsi();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillGlobalPsi() {
    std::vector <UInt> v(deData_.dataSize());
    std::iota(v.begin(), v.end(), 0);
    GlobalPsi_ = computePsi(v);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillFEMatrices(){
    //fill R0 and R1
    FiniteElement<ORDER, mydim, ndim> fe;
    typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
    typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
    Assembler::operKernel(mass, mesh_, fe, R0_);
    Assembler::operKernel(stiff, mesh_, fe, R1_);

    //fill P
    Eigen::SparseLU<SpMat> solver;
    solver.compute(R0_);
    auto X2 = solver.solve(R1_);
    P_ = R1_.transpose() * X2;
}


template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillPsiQuad(){
    for(UInt i = 0; i < Integrator::NNODES; ++i)
        PsiQuad_.row(i) = reference_eval_point<EL_NNODES, mydim>(Integrator::NODES[i]);
}

template<UInt ORDER, UInt mydim, UInt ndim>
Real DataProblem<ORDER, mydim, ndim>::FEintegrate_exponential(const VectorXr& g) const{

    using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1> >;

    Real total_sum = 0.;

    for(UInt triangle = 0; triangle < mesh_.num_elements(); ++triangle){

        Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.getElement(triangle);

// (3) -------------------------------------------------
        Eigen::Matrix<Real, EL_NNODES,1> sub_g;
        for (UInt i = 0; i < EL_NNODES; i++){
            sub_g[i] = g[tri_activated[i].getId()];
        }

// (4) -------------------------------------------------
        Eigen::Matrix<Real, Integrator::NNODES,1> expg = (PsiQuad_*sub_g).array().exp();

        total_sum += expg.dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]))*tri_activated.getMeasure();

    }

    return total_sum;
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat
DataProblem<ORDER, mydim, ndim>::computePsi(const std::vector<UInt>& indices) const{

    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    UInt nnodes = mesh_.num_nodes();
    UInt nlocations = indices.size();
    SpMat psi(nlocations, nnodes);

    std::vector<coeff> triplets;
    triplets.reserve(EL_NNODES*nlocations);

    for(auto it = indices.cbegin(); it != indices.cend(); it++)
    {
        //std::cout << *it << std::endl;
        //operator<<(std::cout, deData_.data(*it));
        Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(deData_.data(*it));

        if(tri_activated.getId() == Identifier::NVAL)
        {
            std::cout << "WARNING: the following observation is not in the domain";
            operator<<(std::cout, deData_.data(*it));
        }
        else
        {
            for(UInt node = 0; node < EL_NNODES ; ++node)
            {
                Real evaluator = tri_activated.evaluate_point(deData_.data(*it), Eigen::Matrix<Real,EL_NNODES,1>::Unit(node));
                triplets.emplace_back(it-indices.cbegin(), tri_activated[node].getId(), evaluator);
            }
        }
    }

    psi.setFromTriplets(triplets.begin(),triplets.end());

    psi.prune(tolerance);
    psi.makeCompressed();

    return psi;
}

// #####################################################################################################################
// ################################################ SPACE TIME PROBLEM #################################################
// #####################################################################################################################

template<UInt ORDER, UInt mydim, UInt ndim>
DataProblem_time<ORDER, mydim, ndim>::DataProblem_time(const std::vector<Point<ndim>>& data, const std::vector<Real>& data_time,
                                                       const UInt& order, const VectorXr& fvec, Real heatStep, UInt heatIter,
                                                       const std::vector<Real>& lambda, const std::vector<Real>& lambda_time,
                                                       const UInt& nfolds, const UInt& nsim, const std::vector<Real>& stepProposals,
                                                       Real tol1, Real tol2, bool print, UInt search, const RNumericMatrix& points,
                                                       const RIntegerMatrix& sides, const RIntegerMatrix& elements,
                                                       const RIntegerMatrix& neighbors, const std::vector<Real>& mesh_time,
                                                       bool isTime) :
        DataProblem<ORDER, mydim, ndim>(data, order, fvec, heatStep, heatIter, lambda, nfolds, nsim, stepProposals,
                                        tol1, tol2, print, search, points, sides, elements, neighbors, isTime),
        deData_time_(data_time, lambda_time), mesh_time_(mesh_time), spline(mesh_time) {

    std::vector<Point<ndim>>& data_ = this->deData_.data();
    std::vector<Real>& data_time_ = deData_time_.data();
    //const Real t_min = *std::min_element(mesh_time_.cbegin(), mesh_time_.cend());
    const Real t_min = mesh_time_.front();
    //const Real t_max = *std::max_element(mesh_time_.cbegin(), mesh_time_.cend());
    const Real t_max = mesh_time_.back();

    // REMOVE POINTS NOT IN THE DOMAIN
    for (auto it = data_.begin(); it != data_.end();) {
        //std::cout << it->getId() << ": " << it->coord()[0] << " " << it->coord()[1] << " " << data_time_[it - data_.begin()] << std::endl;
        Element<this->EL_NNODES, mydim, ndim> tri_activated = this->mesh_.findLocation(data_[it - data_.begin()]);
        if (tri_activated.getId() == Identifier::NVAL || (data_time_[it - data_.begin()] < t_min || data_time_[it - data_.begin()] > t_max)) {
            data_time_.erase(data_time_.begin() + (it - data_.begin()));
            it = data_.erase(it);
            //std::cout << "WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n";
        } else {
            it++;
        }
    }

    std::cout << "WARNING: " << data.size() - data_.size() << " observations removed." << std::endl;
    std::cout << "WARNING: " << data_.size() << " observations used in the algorithm." << std::endl;
    std::ofstream ost("../data/space_time/N_obs.txt");
    ost << data_.size();
    ost.close();

    deData_time_.createMap(data_);
    //deData_time_.printMap(std::cout);

    if(this->isFvecEmpty())
        createMap_Heat(data_time_);

    //! Last part to compute the space matrices
    fillGlobalPsi();

    //! Computation of time matrices
    fillGlobalPhi();
    fillTimeMass();
    fillTimeSecondDerivative();

    //! Computation of space penalty matrix
    //makeLumped();
    fillPenaltySpace();

    //! Assembling space-time matrices
    Upsilon_ = computeUpsilon(this->GlobalPsi_, GlobalPhi_);
    //Upsilon_ = computeUpsilon(this->GlobalPsi_, GlobalPhi_, deData_time_.getMap());
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillGlobalPsi() {
    std::vector <UInt> idx(deData_time_.dataSize());
    std::iota(idx.begin(), idx.end(), 0);
/*
    std::vector<Point<ndim>>& data_ = this->deData_.data();
    idx.reserve(deData_time_.getID_noD().size());

    auto it = data_.cbegin();
    for (UInt i : deData_time_.getID_noD()) {
        while (it->id() != i)
            ++it;
        idx.push_back(it - data_.cbegin());
    }

    for (UInt i : idx)
        std::cout << i << " ";
    std::cout << std::endl;
*/
    this->GlobalPsi_ = this->computePsi(idx);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillGlobalPhi(void)
{
    //Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
    UInt M = spline.num_knots()-SPLINE_DEGREE-1;
    UInt m = deData_time_.timeSize_noD();

    GlobalPhi_.resize(m, M);
    Real value;

    for(UInt i = 0; i < m; ++i)
    {
        for(UInt j = 0; j < M; ++j)
        {
            value = spline.BasisFunction(j, this->deData_time_.data_noD()[i]);
            if(value != 0)
            {
                GlobalPhi_.coeffRef(i,j) = value;
            }
        }
    }
/*
    std::cout << "phi" << std::endl;
    std::cout << std::setw(7);
    for (UInt i = 0; i < GlobalPhi_.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < GlobalPhi_.cols(); ++j)
            std::cout << GlobalPhi_.coeff(i, j) << std::setw(7);
        std::cout << std::endl;
    }
*/
    GlobalPhi_.makeCompressed();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillTimeMass(void)
{
    Assembler::operKernel(spline, K0_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillTimeSecondDerivative(void)
{
    //Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
    Assembler::operKernel(spline, Pt_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::makeLumped(void) {
    R0_Lumped = this->R0_;
    VectorXr diag;
    UInt l = R0_Lumped.rows();
    diag.resize(l);
    for (UInt i = 0; i < l; ++i) {
        for (int j = 0; j < l; ++j)
            diag[i]+=R0_Lumped.coeff(i,j);
    }
    std::cout<<"diag[0] "<<diag[0]<<std::endl;
    R0_Lumped = diag.asDiagonal();
    std::cout<<"R0_Lumped: "<<R0_Lumped.coeff(0,0)<<" "<<R0_Lumped.coeff(0,1)<<" "<<R0_Lumped.coeff(1,1)<<std::endl;
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillPenaltySpace(void)
{
    // Update R1_
    SpMat R1_temp;
    UInt N_r = this->R1_.rows(), N_c = this->R1_.cols();
    UInt M_r = K0_.rows(), M_c = K0_.cols();

    SpMat K0temp(K0_);
    K0temp.setIdentity();

    //R1_temp = kroneckerProduct(this->getStiffness(), K0temp);
    R1_temp = kroneckerProduct(K0temp, this->getStiffness());
    R1_temp.makeCompressed();

    // Update R0
    SpMat R0_temp;
    //R0_temp = kroneckerProduct(this->getMass(), K0temp);
    R0_temp = kroneckerProduct(K0temp, this->getMass());
    R0_temp.makeCompressed();

    // Compute Space Penalty
    Ps_.resize(N_r * M_r, N_c * M_c);
    Eigen::SparseLU<SpMat> factorized_R0(R0_temp);
    Ps_ = (R1_temp).transpose()*factorized_R0.solve(R1_temp);     // R == _R1^t*R0^{-1}*R1
    Ps_.makeCompressed();
}

// LUMPED VERSION
/*
template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillPenaltySpace(void)
{
    // Update R1_
    SpMat R1_temp;
    UInt N_r = this->R1_.rows(), N_c = this->R1_.cols();
    UInt M_r = K0_.rows(), M_c = K0_.cols();
    SpMat K0temp;
    K0temp.resize(M_r,M_c);
//    for (UInt i = 0; i < K0temp.rows(); ++i) {
//        for (int j = 0; j < K0temp.cols(); ++j){
//            std::cout<<K0temp.coeff(i,j)<<" ";
//        }
//        std::cout<<std::endl;
//    }
    K0temp.setIdentity();
    R1_temp = kroneckerProduct(this->getStiffness(), K0temp);
    R1_temp.makeCompressed();
    // Update R0
    SpMat R0_temp;
    R0_temp = kroneckerProduct(getLumped(),K0temp);
    R0_temp.makeCompressed();
    // Compute Space Penalty
    Ps_.resize(N_r * M_r, N_c * M_c);
    Eigen::SparseLU<SpMat> factorized_R0(R0_temp);
    Ps_ = (R1_temp).transpose()*factorized_R0.solve(R1_temp);     // R == _R1^t*R0^{-1}*R1
    Ps_.makeCompressed();
}
*/

template<UInt ORDER, UInt mydim, UInt ndim>
MatrixXr DataProblem_time<ORDER, mydim, ndim>::fillPhiQuad(UInt time_node) const {
    MatrixXr phi;
    phi.resize(Integrator_t::NNODES,SPLINE_DEGREE+1);
    Real t_a = mesh_time_[time_node], t_b = mesh_time_[time_node+1];
    std::array<Real,Integrator_t::NNODES> ref_nodes;
    for(UInt k = 0; k < Integrator_t::NNODES; ++k)
        ref_nodes[k]=((t_b-t_a)*Integrator_t::NODES[k]+t_a+t_b)/2;
    for(UInt j = 0; j < phi.cols(); j++){
        for(UInt i = 0; i < phi.rows(); i++)
            phi(i,j) = spline.BasisFunction(time_node+j, ref_nodes[i]);
    }
    return phi;
}

template<UInt ORDER, UInt mydim, UInt ndim>
Real DataProblem_time<ORDER, mydim, ndim>::FEintegrate_exponential(const VectorXr& g) const{
    // Kronecker product of the Gauss quadrature rules weights
    VectorXr weights_kronecker;
    weights_kronecker.resize(Integrator::NNODES*Integrator_t::NNODES);
    UInt k=0;
    for (UInt i = 0;  i < Integrator_t::NNODES; i++) {
        for (UInt j = 0;  j < Integrator::NNODES; j++){
            weights_kronecker[k] = Integrator::WEIGHTS[j]*Integrator_t::WEIGHTS[i];
            ++k;
        }
    }

    // Initialization
    Real int1 = 0.;
    const MatrixXr& PsiQuad = this->getPsiQuad(); //It is always the same
    UInt global_idx = 0; //index that keeps track of the first B-spline basis function active in the current time-interval
    for (int time_step = 0; time_step < getNumNodes_time()-1;  ++time_step) {
        MatrixXr PhiQuad = fillPhiQuad(time_step); //PhiQuad changes at each time interval
        MatrixXr Phi_kronecker_Psi = kroneckerProduct_Matrix(PhiQuad,PsiQuad);
        for(UInt triangle = 0; triangle < this->getNumElements(); triangle++) {
            Element<this->EL_NNODES, mydim, ndim> tri_activated = this->getElement(triangle);
//// (1) -------------------------------------------------
            VectorXr sub_g;
            sub_g.resize(Phi_kronecker_Psi.cols());
            UInt k=0; //index for sub_g
            for (int j = global_idx; j < global_idx+PhiQuad.cols(); ++j) {
                for (UInt i = 0; i < PsiQuad.cols(); ++i){
                    sub_g[k++]=g[tri_activated[i].getId()+this->getNumNodes()*j];
                }
            }

            VectorXr expg = (Phi_kronecker_Psi*sub_g).array().exp();
            int1 += expg.dot(weights_kronecker.transpose()) * tri_activated.getMeasure() * (getMesh_time()[time_step+1]-getMesh_time()[time_step])/2;

        }
        ++global_idx;
    }
    return int1;
}


template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const SpMat &psi, const SpMat &phi) const {

    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    const UInt phi_r = phi.rows();
    const UInt phi_c = phi.cols();
    const UInt psi_c = psi.cols();

    if (deData_time_.getID_noD().size() != deData_time_.dataSize())
        std::cout << "WARNING: " << deData_time_.dataSize() - deData_time_.getID_noD().size() << " spatial duplicates!" << std::endl;

    if (deData_time_.data_noD().size() != deData_time_.dataSize())
        std::cout << "WARNING: " << deData_time_.dataSize() - deData_time_.data_noD().size() << " temporal duplicates!" << std::endl;

/*
    std::cout << "psi" << std::endl;
    for (UInt i = 0; i < psi.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < psi_c; ++j)
            std::cout << psi.coeff(i, j) << std::setw(3);
        std::cout << std::endl;
    }

    std::cout << "phi" << std::endl;
    std::cout << std::setw(7);
    for (UInt i = 0; i < phi.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < phi_c; ++j)
            std::cout << phi.coeff(i, j) << std::setw(7);
        std::cout << std::endl;
    }
*/
    std::vector<coeff> Upsilon_tripletList;
    Upsilon_tripletList.reserve(deData_time_.dataSize() * phi_c * psi_c);

    UInt global_row_counter = 0;
    for (UInt i = 0; i < phi_r; ++i) {
        for (UInt j = 0; j < deData_time_.dataSize(); ++j) {
            if (deData_time_.data()[j] == deData_time_.data_noD()[i]) {
                //std::cout << "data[j]: " << deData_time_.data()[j] << "data_noD[i]: " << deData_time_.data_noD()[i] << std::endl;
                //std::cout << "datum: " << j << std::endl;
                SpMat localKProd_(1,phi_c * psi_c);
                localKProd_ = kroneckerProduct(phi.row(i), psi.row(j));
                for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                    Upsilon_tripletList.emplace_back(global_row_counter,idx,localKProd_.coeff(0,idx));
                ++global_row_counter;
            }
        }
    }

    SpMat upsilon(deData_time_.dataSize(), phi_c * psi_c);
    upsilon.setFromTriplets(Upsilon_tripletList.begin(), Upsilon_tripletList.end());

    upsilon.prune(tolerance);
    upsilon.makeCompressed();
/*
    std::cout << "upsilon" << std::endl;
    for (UInt i = 0; i < upsilon.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < upsilon.cols(); ++j)
            std::cout << upsilon.coeff(i, j) << std::setw(3);
        std::cout << std::endl;
    }
*/
    //Upsilon_tripletList.clear();
    //localPhi_tripletList.clear();

    return upsilon;
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const SpMat& psi, const SpMat& phi,
                                                           const std::map<UInt, std::set<UInt>>& data_noD) const {
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    const UInt psi_r = psi.rows();
    const UInt psi_c = psi.cols();
    const UInt phi_c = phi.cols();
/*
    std::cout << "psi" << std::endl;
    for (UInt i = 0; i < psi_r; ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < psi_c; ++j)
            std::cout << psi.coeff(i, j) << std::setw(3);
        std::cout << std::endl;
    }

    std::cout << "phi" << std::endl;
    for (UInt i = 0; i < phi.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < phi_c; ++j)
            std::cout << phi.coeff(i, j) << std::setw(7);
        std::cout << std::endl;
    }
*/
    std::vector<coeff> Upsilon_tripletList; // needed to build the Upsilon_ matrix in the end
    std::vector<coeff> localPhi_tripletList; // needed to build the "current" localPhi matrix for each spatial location

    UInt global_row_counter = 0; // current number of computed rows of the Upsilon_ matrix

    std::map<UInt, std::set<UInt>>::const_iterator current_it = data_noD.cbegin();

    // Construction of the Upsilon_ matrix (location by location)
    for (UInt i = 0; i < psi_r; ++i) { // row by row (of psi)
        localPhi_tripletList.clear();
        localPhi_tripletList.reserve(current_it->second.size() * phi_c);
        UInt local_row_counter = 0; // current number of rows to compute
        for (UInt j : current_it->second) {
            for (UInt k = 0; k < phi_c; ++k) {
                localPhi_tripletList.emplace_back(local_row_counter, k, phi.coeff(j, k));
            }
            ++local_row_counter;
        }
        SpMat localPhi(local_row_counter, phi_c);
        localPhi.setFromTriplets(localPhi_tripletList.begin(), localPhi_tripletList.end());

        SpMat localKProd(local_row_counter, psi_c * phi_c);
        localKProd = kroneckerProduct(psi.row(i), localPhi);

        Upsilon_tripletList.reserve(global_row_counter * psi_c * phi_c + localKProd.rows() * psi_c * phi_c);
        for(UInt idx = 0; idx < localKProd.outerSize(); ++idx)
            for(SpMat::InnerIterator it(localKProd,idx); it; ++it) {
                Upsilon_tripletList.emplace_back(it.row() + global_row_counter, it.col(), it.value());
            }

        global_row_counter += local_row_counter;

        if(current_it != data_noD.cend())
            ++current_it;
    }

    SpMat upsilon(global_row_counter, psi_c * phi_c);
    upsilon.setFromTriplets(Upsilon_tripletList.begin(), Upsilon_tripletList.end());

    upsilon.prune(tolerance);
    upsilon.makeCompressed();
/*
    std::cout << "upsilon" << std::endl;
    for (UInt i = 0; i < upsilon.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < upsilon.cols(); ++j)
            std::cout << upsilon.coeff(i, j) << std::setw(3);
        std::cout << std::endl;
    }
*/
    //Upsilon_tripletList.clear();
    //localPhi_tripletList.clear();

    return upsilon;
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const std::vector<UInt>& indices) const {

    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    SpMat psi = this->computePsi(indices);

    SpMat phi(indices.size(), GlobalPhi_.cols());
    std::vector<coeff> Phi_tripletList;
    Phi_tripletList.reserve(indices.size() * GlobalPhi_.cols());

    UInt index_GlobalPhi_;
    for (UInt i = 0; i < indices.size(); ++i) {
        index_GlobalPhi_ = 0;
        while(data_time(indices[i])!=deData_time_.data_noD()[index_GlobalPhi_])
            ++index_GlobalPhi_;
        //std::cout << indices[i] << " " << index_GlobalPhi_ << " " << data_time(indices[i]) << " " << deData_time_.data_noD()[index_GlobalPhi_] << std::endl;
        //std::cout << GlobalPhi_.coeff(index_GlobalPhi_, 4) << std::endl;
        for (UInt j = 0; j < GlobalPhi_.cols(); ++j) {
            Phi_tripletList.emplace_back(i, j, GlobalPhi_.coeff(index_GlobalPhi_, j));
        }
    }

    //std::cout << phi.rows() << " " << indices.size() << " " << GlobalPhi_.rows() << std::endl;
    //std::cout << phi.cols() << " " << GlobalPhi_.cols();

    phi.setFromTriplets(Phi_tripletList.begin(), Phi_tripletList.end());

    phi.prune(tolerance);
    phi.makeCompressed();

    std::vector<coeff> Upsilon_tripletList;
    Upsilon_tripletList.reserve(indices.size() * phi.cols() * psi.cols());

    // Construction of the Upsilon_ matrix (location by location)
    for (UInt i = 0; i < indices.size(); ++i) {
        SpMat localKProd_(1, phi.cols() * psi.cols());
        localKProd_ = kroneckerProduct(phi.row(i), psi.row(i));
        for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
            Upsilon_tripletList.emplace_back(i, idx, localKProd_.coeff(0, idx));
    }

    SpMat upsilon(indices.size(), phi.cols() * psi.cols());
    upsilon.setFromTriplets(Upsilon_tripletList.begin(), Upsilon_tripletList.end());

    upsilon.prune(tolerance);
    upsilon.makeCompressed();

    return upsilon;
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::createMap_Heat(const std::vector<Real>& data_time){
    //This function creates a Map useful to compute the Heat_Initialization Process
    //According to the time discretization only a few splines are active for each data_time hence
    //that observation must be taken into consideration in the according part of Heat_Init

    const Real tol = 0;

    const UInt M = spline.num_knots()-SPLINE_DEGREE-1;

    for (int j = 0; j < M; ++j) data_Heat_[j]={};
    for (int i = 0; i < data_time.size(); ++i) {
        for (int j = 0; j < M; ++j) {
            if(std::abs(spline.BasisFunction(j,data_time[i])) >= tol)
                data_Heat_[j].push_back(i);
        }
    }
}

#endif //DEV_FDAPDE_DATA_PROBLEM_IMP_H

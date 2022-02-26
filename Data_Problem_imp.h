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
            //std::cout << it - data_.begin() << std::endl;
        }
    }

    // FILL MATRICES
    fillFEMatrices();
    fillPsiQuad();

    //if(!isTime)
        //fillGlobalPsi();
    if(!isTime) {
        std::vector<UInt> v(deData_.dataSize());
        std::iota(v.begin(), v.end(), 0);
        GlobalPsi_ = computePsi(v);
    }
}
/*
template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillGlobalPsi() {
    std::vector <UInt> v(deData_.dataSize());
    std::iota(v.begin(), v.end(), 0);
    GlobalPsi_ = computePsi(v);
}
*/
template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem<ORDER, mydim, ndim>::fillFEMatrices() {
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
                                                       bool isTime, bool isTimeDiscrete, bool flagMass, bool flagLumped) :
        DataProblem<ORDER, mydim, ndim>(data, order, fvec, heatStep, heatIter, lambda, nfolds, nsim, stepProposals,
                                        tol1, tol2, print, search, points, sides, elements, neighbors, isTime),
        deData_time_(data_time, lambda_time), mesh_time_(mesh_time), spline_(mesh_time), flagMass_(flagMass),
        flagLumped_(flagLumped) {

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

    // FILL SPACE MATRIX
    std::vector<UInt> v(this->deData_.dataSize());
    std::iota(v.begin(), v.end(), 0);
    this->GlobalPsi_ = this->computePsi(v);

    if(isTimeDiscrete) {
        deData_time_.setTimes2Locations();
        //Upsilon_indices_.resize(deData_time_.getNTimes());
        deData_time_.printTimes2Locations(std::cout);
    }

    if(this->isFvecEmpty())
        setDataHeat();

    //! Last part to compute the space matrices
    //fillGlobalPsi();

    //! Computation of time matrices
    fillGlobalPhi();
    fillTimeMass();
    fillTimeSecondDerivative();

    //! Computation of space and time penalty matrices
    fillPenaltySpace();
    fillPenaltyTime();

    //! Assembling space-time matrices
    Upsilon_ = computeUpsilon(GlobalPhi_, this->GlobalPsi_);
    //Upsilon_ = computeUpsilon(this->GlobalPsi_, GlobalPhi_, deData_time_.getMap());
}
/*
template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillGlobalPsi() {
    std::vector <UInt> idx(deData_time_.dataSize());
    std::iota(idx.begin(), idx.end(), 0);

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

    this->GlobalPsi_ = this->computePsi(idx);
}
*/
template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillGlobalPhi(void)
{
    //Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
    //UInt M = spline_.num_knots()-SPLINE_DEGREE-1;
    const UInt M = getSplineNumber();
    const UInt m = deData_time_.getNTimes();

    GlobalPhi_.resize(m, M);
    Real value;

    for(UInt i = 0; i < m; ++i)
    {
        for(UInt j = 0; j < M; ++j)
        {
            value = spline_.BasisFunction(j, this->deData_time_.time(i));
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
    Spline<SPLINE_DEGREE, 0> spline_0(mesh_time_);
    Assembler::operKernel(spline_, K0_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillTimeSecondDerivative(void)
{
    //Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
    Assembler::operKernel(spline_, Pt_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::makeLumped(const SpMat& mass) const {

    VectorXr diag = mass * VectorXr::Ones(mass.cols());
    SpMat lumped_mass(diag.asDiagonal());

    return lumped_mass;
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillPenaltySpace()
{
    // Update R1_
    SpMat R1_temp;
    UInt N_r = this->R1_.rows(), N_c = this->R1_.cols();
    UInt M_r = K0_.rows(), M_c = K0_.cols();

    SpMat K0temp(K0_);
    if(!flagMass_)
        K0temp.setIdentity();

    R1_temp = kroneckerProduct(K0temp, this->getStiffness());
    R1_temp.makeCompressed();

    // Update R0
    SpMat R0_temp;
    R0_temp = kroneckerProduct(K0temp, this->getMass());
    R0_temp.makeCompressed();

    if(flagLumped_)
        R0_temp = makeLumped(R0_temp);

    // Compute Space Penalty
    Ps_.resize(N_r * M_r, N_c * M_c);
    Eigen::SparseLU<SpMat> factorized_R0(R0_temp);
    Ps_ = (R1_temp).transpose()*factorized_R0.solve(R1_temp);     // R == _R1^t*R0^{-1}*R1
    Ps_.makeCompressed();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillPenaltyTime() {
    SpMat mass_temp(this->getMass());
    if(!flagMass_)
        mass_temp.setIdentity();
    Pt_ = kroneckerProduct(getPt(),mass_temp);
    Pt_.makeCompressed();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::setDataHeat() {
    //This function creates a Map useful to compute the Heat_Initialization Process
    //const UInt M = spline_.num_knots()-SPLINE_DEGREE-1;
    const UInt M = getSplineNumber();
    data_Heat_.resize(M);

    //! ### POSSIBLE PARALLEL openMP ###
    for (int i = 0; i < deData_time_.getNTimes(); ++i) {
        for (int j = 0; j < M; ++j) {
            if(spline_.BasisFunction(j,deData_time_.time(i))!= 0) //std::abs(spline.BasisFunction(j,data_time[i]))>=tol)
                data_Heat_[j].push_back(i);
        }
    }
}

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
            phi(i,j) = spline_.BasisFunction(time_node+j, ref_nodes[i]);
    }
    return phi;
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const SpMat &phi, const SpMat &psi) const {

    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    const UInt phi_r = phi.rows();
    const UInt phi_c = phi.cols();
    const UInt psi_c = psi.cols();

    //if (deData_time_.getID_noD().size() != deData_time_.dataSize())
        //std::cout << "WARNING: " << deData_time_.dataSize() - deData_time_.getID_noD().size() << " spatial duplicates!" << std::endl;

    if (deData_time_.getNTimes() != deData_time_.dataSize())
        std::cout << "WARNING: " << deData_time_.dataSize() - deData_time_.getNTimes() << " temporal duplicates!" << std::endl;
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

    if(deData_time_.getNTimes() != deData_time_.dataSize()) // time duplicates: phi_r < psi_r
    {
        //UInt global_row_counter = 0;
        for(UInt i = 0; i < phi_r; ++i) {
            const std::vector<UInt>& v = deData_time_.getTimes2Locations(i);
            for(UInt j : v) {
                //Upsilon_indices_[j] = global_row_counter;
                SpMat localKProd_(1, phi_c * psi_c);
                std::cout << "riga phi: " << i << ", riga psi: " << j << std::endl;
                localKProd_ = kroneckerProduct(phi.row(i), psi.row(j));
                for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                    Upsilon_tripletList.emplace_back(j,idx,localKProd_.coeff(0, idx));
                //++global_row_counter;
            }
        }
    }
    else // NO time duplicates: phi_r = psi_r = #observations
    {
        for(UInt i = 0; i < phi_r; ++i) {
            SpMat localKProd_(1, phi_c * psi_c);
            localKProd_ = kroneckerProduct(phi.row(i), psi.row(i));
            for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                Upsilon_tripletList.emplace_back(i,idx,localKProd_.coeff(0, idx));
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

    return upsilon;
}

/* // ALTERNATIVE VERSION (WITHOUT Upsilon_indices_)
template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const std::vector<UInt>& indices) const
{
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    SpMat psi = this->computePsi(indices);
    SpMat phi = GlobalPhi_;

    const UInt phi_r = phi.rows();
    const UInt phi_c = phi.cols();
    const UInt psi_c = psi.cols();

    std::vector<coeff> Upsilon_tripletList;
    Upsilon_tripletList.reserve(indices.size() * phi_c * psi_c);

    if(deData_time_.getNTimes() != deData_time_.dataSize()) // time duplicates
    {
        UInt global_row_counter = 0;
        for(UInt i = 0; i < phi_r; ++i) {
            for(UInt j = 0; j < deData_time_.getTimes2Locations(i).size(); ++j) {
                auto it = std::find(indices.begin(), indices.end(), deData_time_.getTimes2Locations(i)[j]);
                //auto it = std::lower_bound(indices.begin(), indices.end(), deData_time_.getTimes2Locations(i)[j]);
                if(it != indices.end()) {
                //if(*it == deData_time_.getTimes2Locations(i)[j]) {
                    SpMat localKProd_(1, phi_c * psi_c);
                    localKProd_ = kroneckerProduct(phi.row(i), psi.row(indices[it-indices.begin()]));
                    for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                        Upsilon_tripletList.emplace_back(global_row_counter,idx,localKProd_.coeff(0, idx));
                }
                ++global_row_counter;
            }
        }
    }
    else // NO time duplicates
    {
        for(UInt i = 0; i < indices.size(); ++i) {
            SpMat localKProd_(1, phi_c * psi_c);
            localKProd_ = kroneckerProduct(phi.row(indices[i]), psi.row(i));
            for (UInt idx = 0; idx < localKProd_.outerSize(); ++idx)
                Upsilon_tripletList.emplace_back(i,idx,localKProd_.coeff(0, idx));
        }
    }

    SpMat upsilon(indices.size(), phi_c * psi_c);
    upsilon.setFromTriplets(Upsilon_tripletList.begin(), Upsilon_tripletList.end());

    upsilon.prune(tolerance);
    upsilon.makeCompressed();

    return upsilon;
}
*/

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const std::vector<UInt>& indices) const
{
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    Eigen::SparseMatrix<Real,Eigen::RowMajor> upsilon(indices.size(), Upsilon_.cols());

    for(UInt i = 0; i < indices.size(); ++i) {
        //upsilon.row(i) = Upsilon_.row(Upsilon_indices_[indices[i]]);
        upsilon.row(i) = Upsilon_.row(indices[i]);
    }

    upsilon.prune(tolerance);
    upsilon.makeCompressed();
/*
    std::cout << "upsilon(indices)" << std::endl;
    for (UInt i = 0; i < upsilon.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < upsilon.cols(); ++j)
            std::cout << upsilon.coeff(i, j) << std::setw(3);
        std::cout << std::endl;
    }
*/
    return upsilon;

}



































/*
template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const SpMat& psi, const SpMat& phi,
                                                           const std::map<UInt, std::set<UInt>>& data_noD) const {
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    const UInt psi_r = psi.rows();
    const UInt psi_c = psi.cols();
    const UInt phi_c = phi.cols();

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

    std::cout << "upsilon" << std::endl;
    for (UInt i = 0; i < upsilon.rows(); ++i) {
        std::cout << i << "-th row: " << std::endl;
        std::cout << std::setw(7);
        for (UInt j = 0; j < upsilon.cols(); ++j)
            std::cout << upsilon.coeff(i, j) << std::setw(3);
        std::cout << std::endl;
    }

    //Upsilon_tripletList.clear();
    //localPhi_tripletList.clear();

    return upsilon;
}
*/
/*
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
*/



#endif //DEV_FDAPDE_DATA_PROBLEM_IMP_H

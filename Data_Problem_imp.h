//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_DATA_PROBLEM_IMP_H
#define DEV_FDAPDE_DATA_PROBLEM_IMP_H

#include "Integration.h"

template<UInt ORDER, UInt mydim, UInt ndim>
DataProblem<ORDER, mydim, ndim>::DataProblem(const std::vector<Point<ndim>>& data, const UInt& order,
                                             const VectorXr& fvec, Real heatStep, UInt heatIter,
                                             const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                                             const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print,
                                             UInt search, const RNumericMatrix& points, const RIntegerMatrix& sides,
                                             const RIntegerMatrix& elements, const RIntegerMatrix& neighbors, bool isTime) :
    deData_(data, order, fvec, heatStep, heatIter, lambda, nfolds, nsim, stepProposals, tol1, tol2, print, search),
    mesh_(points, sides, elements, neighbors, search) {
    /*
    // PROJECTION
    if(mydim == 2 && ndim == 3){
        std::cout << "##### DATA PROJECTION #####\n";
        projection<ORDER, mydim, ndim> projection(mesh_, data);
        data = projection.computeProjection();
    }

    // REMOVE POINTS NOT IN THE DOMAIN
    for(auto it = data.begin(); it != data.end(); ){
        Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(data[it - data.begin()]);
        if(tri_activated.getId() == Identifier::NVAL)
        {
            it = data.erase(it);
            std::cout << "WARNING: an observation is not in the domain. It is removed and the algorithm proceeds.\n";
        }
        else {
            it++;
        }
    }
    */
    // FILL MATRICES
    fillFEMatrices();
    fillPsiQuad();
    if(!isTime)   //!***
        fillGlobalPsi(); //!***
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
    P_ = R1_.transpose()* X2;
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

        Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.findLocation(deData_.data(*it));

        if(tri_activated.getId() == Identifier::NVAL)
        {
            std::cout << "WARNING: the following observation is not in the domain";
            //(deData_.getDatum(*it)).print(std::cout);
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
DataProblem_time<ORDER, mydim, ndim>::DataProblem_time(const std::vector<Point<ndim>>& data, std::vector<Real>& data_time,
                                                       const UInt& order, const VectorXr& fvec, Real heatStep, UInt heatIter,
                                                       const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                                                       const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print,
                                                       UInt search, const RNumericMatrix& points, const RIntegerMatrix& sides,
                                                       const RIntegerMatrix& elements, const RIntegerMatrix& neighbors,
                                                       const std::vector<Real>& mesh_time, bool isTime) :
        DataProblem<ORDER, mydim, ndim>(data, order, fvec, heatStep, heatIter, lambda, nfolds, nsim, stepProposals,
                                        tol1, tol2, print, search, points, sides, elements, neighbors, isTime),
        deData_time_(data, data_time), mesh_time_(mesh_time), spline(mesh_time) {
    //! Last part to compute the space matrices
    fillGlobalPsi();
    //! Computation of time matrices
    fillGlobalPhi();
    fillTimeMass();
    fillTimeSecondDerivative();
    fillPhiQuad();

    //! Assembling space-time matrices
    Upsilon_ = computeUpsilon(this->GlobalPsi_, GlobalPhi_, deData_time_.getMap());
}

template<UInt ORDER, UInt mydim, UInt ndim>
void DataProblem_time<ORDER, mydim, ndim>::fillGlobalPsi() {
    std::vector <UInt> v=deData_time_.getID_noDup();
    this->GlobalPsi_ = this->computePsi(v);
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
        for(UInt j = 0; j<M; ++j)
        {
            value = spline.BasisFunction(j, this->deData_time.data_noD()[i]);
            if(value != 0)
            {
                GlobalPhi_.coeffRef(i,j) = value;
            }
        }
    }
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
    Spline<SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
    Assembler::operKernel(spline, Pt_);
}

template<UInt ORDER, UInt mydim, UInt ndim>
SpMat DataProblem_time<ORDER, mydim, ndim>::computeUpsilon(const SpMat& psi, const SpMat& phi,
                                                      const std::map<UInt, std::set<UInt>>& data_noD) const{
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    const UInt psi_r = psi.rows();
    const UInt psi_c = psi.cols();
    const UInt phi_c = phi.cols();

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

    //Upsilon_tripletList.clear();
    //localPhi_tripletList.clear();

    return upsilon;
}

#endif //DEV_FDAPDE_DATA_PROBLEM_IMP_H

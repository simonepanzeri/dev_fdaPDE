//
// Created by simonepanzeri on 24/11/2021.
//

#ifndef DEV_FDAPDE_DE_DATA_H
#define DEV_FDAPDE_DE_DATA_H

#include "FdaPDE.h"
#include "Point.h"

template<UInt ndim>
class DEData {
private:
    // Data = spatial locations.
    std::vector<Point<ndim>> data_;
    // Finite element order.
    UInt order_;
    // Initial coefficients for the density.
    VectorXr fvec_;
    // Time step parameter for the heat diffusion process.
    Real heatStep_;
    // Number of iterations for the heat diffusion process.
    UInt heatIter_;
    // Penalization parameters (in space). The best one is chosen with k fold cross validation.
    std::vector<Real> lambda_;
    // Number of folds for cross validation.
    UInt Nfolds_;
    // Number of simulations of the optimization algorithm.
    UInt nsim_;
    // Optimization parameters for fixed step methods.
    std::vector<Real> stepProposals_;
    // Tolerances for optimization algorithm termination criteria.
    Real tol1_;
    Real tol2_;
    // A boolean: true if the user wants to see the value of the functional printed during the optimization descent.
    bool print_;
    // Integer specifying the search algorithm type (tree or naive search algorithm).
    UInt search_;

public:
    //! Constructors
    DEData() {};

    explicit DEData(const std::vector<Point<ndim>>& data, const UInt& order, const VectorXr& fvec, Real heatStep,
                    UInt heatIter, const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                    const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print, UInt search);

    //! Getters
    //! A method to access the data.
    std::vector<Point<ndim>>& data() {return data_;}
    //! A const method to access the data.
    const std::vector<Point<ndim>>& data() const {return data_;}
    //! A method to access a datum.
    Point<ndim>& data(UInt i) {return data_[i];}
    //! A const method to access a datum.
    const Point<ndim>& data(UInt i) const {return data_[i];}
    //! A method returning the number of observations.
    UInt dataSize() const {return data_.size();}
    //! A method returning the the input order.
    UInt getOrder() const {return order_;}
    //! A method returning the initial coefficients for the density.
    VectorXr getFvec() const {return fvec_;}
    //! A method returning the heat diffusion process alpha parameter.
    Real getHeatStep() const {return heatStep_;}
    //! A method returning the number of iterations for the heat diffusion process.
    UInt getHeatIter() const {return heatIter_;}
    //! A method returning a bool which says if there is a user's initial density.
    bool isFvecEmpty() const {return fvec_.size() == 0;}
    //! A method returning the penalization parameters (in space).
    Real getLambda(UInt i) const {return lambda_[i];}
    //! A method returning the number of lambdas (in space).
    UInt getNlambda()  const {return lambda_.size();}
    //! A method returning the number of folds for CV.
    UInt getNfolds()  const {return Nfolds_;}
    //! A method returning the number of iterations to use in the optimization algorithm.
    UInt getNsimulations() const {return nsim_;}
    //! A method returning the number of parameters for fixed step methods.
    UInt getNstepProposals() const {return stepProposals_.size();}
    //! A method returning a parameter for fixed step methods.
    Real getStepProposals(UInt i) const {return stepProposals_[i];}
    //! A method returning the tolerance for optimization algorithm first termination criteria.
    Real getTol1() const {return tol1_;}
    //! A method returning the tolerance for optimization algorithm second termination criteria.
    Real getTol2() const {return tol2_;}
    //! A method returning the boolean print member.
    bool Print() const {return print_;}
    //! A method returning the integer that specifies the search algorithm type.
    UInt getSearch() const {return search_;}

    //! Print
    //! A method printing data.
    void printData(std::ostream & out) const;
};

template <UInt ndim>
class DEData_time {
private:
    //! A vector containing the time observations from data (it may contain duplicates).
    std::vector<Real> data_time_;
    //! A vector containing the time observations (without duplicates).
    std::vector<Real> data_time_noD_;
    //! A data structure containing in data_noD.first the IDs of the locations with no duplicates
    //! and in data_noD.second the time indices from data_time_noD_ in which each spatial observation is observed.
    std::map<UInt, std::set<UInt>> data_noD_;
    //! Penalization parameters (in time). The best one is chosen with k fold cross validation.
    std::vector<Real> lambda_time_;

    //! A method to insert in data_noD_ the point ID i as key and the position that the time instant t (at which i is
    //! observed) take in data_time_noD_ as value in the set corresponding to the i-th key (considering duplicated IDs
    //! with respect to i, if any).
    void insert(UInt i, Real t) {data_noD_[i].insert(std::find(data_time_noD_.begin(), data_time_noD_.end(), t)-data_time_noD_.begin());}
    void insert(UInt i) {insert(i,i);}
    //! An helper function for the construction of the data_noD_ map that checks whether a certain spatial location with
    //! ID equal to k has already been inserted in the map, by relying on a helper set.
    bool isAlready(UInt k, const std::set<UInt>& set_helper) const {return (set_helper.find(k)!=set_helper.end());}

public:
    //! Constructor
    DEData_time(const std::vector<Real>& data_time, const std::vector<Real>& lambda_time);

    //! A method to clear data_time_noD_ and data_noD_, removing data that are not inside the spatio-temporal domain of interest.
    void createMap (const std::vector<Point<ndim>>& data);

    //! Getters
    //! A method to access the data.
    std::vector<Real>& data() {return data_time_;}
    //! A const method to access the data.
    const std::vector<Real>& data() const {return data_time_;}
    //! A const method to access the data stored in the data_noD_ map.
    const std::map<UInt, std::set<UInt>>& getMap() const {return data_noD_;}
    //! A method to access the data (without duplicates).
    std::vector<Real>& data_noD() {return data_time_noD_;}
    //! A const method to access the time data (without duplicates).
    const std::vector<Real>& data_noD() const {return data_time_noD_;}
    //! A method returning the number of observations.
    UInt dataSize() const {return data_time_.size();}
    //! A method returning the number of data time evaluations.
    UInt timeSize_noD() const {return data_time_noD_.size();}
    //! A method returning the penalization parameters (in time).
    Real getLambda_time(UInt i) const {return lambda_time_[i];}
    //! A method returning the number of lambdas (in time).
    UInt getNlambda_time()  const {return lambda_time_.size();}

    //! The following two methods need to be used coupled:
    //! A method returning the IDs of spatial data points (without duplicates).
    std::vector<UInt> getID_noD() const;
    //! A method returning the time indices in which point id appears.
    const std::set<UInt>& getTimesIndices(UInt id) const;

    //! Print
    //! A method printing data.
    void printData(std::ostream& out) const;
    //! A method printing the map.
    void printMap(std::ostream& out) const;
};

#include "DE_Data_imp.h"

#endif //DEV_FDAPDE_DE_DATA_H

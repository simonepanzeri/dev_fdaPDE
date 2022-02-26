//
// Created by simonepanzeri on 02/12/2021.
//

#ifndef DEV_FDAPDE_DESCENT_DIRECTION_H
#define DEV_FDAPDE_DESCENT_DIRECTION_H

#include "FdaPDE.h"

// This file contains the direction search technique useful for the optimization algorithm of the Density Estimation problem

/*! @brief An abstract class for computing the descent direction. The father is pure virtual; the right direction
is computed inside the children according to a proper method chosen.
*/

template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionBase{
protected:
    // to give generality if you want to add other children: T is a typename for FunctionalProblem (spatial problem, by
    // default) or for FunctionalProblem_time (spatio-temporal problem)
    const T& funcProblem_;

public:
    //! A constructor
    DirectionBase(const T& fp): funcProblem_(fp){};
    //! A destructor.
    virtual ~DirectionBase(){};
    //! A pure virtual clone method.
    virtual std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const = 0;
    //! A pure virtual method to compute the descent direction.
    virtual VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) = 0;
    //! A pure virtual method to reset all the old parameters.
    virtual void resetParameters() = 0;
};


//! @brief A class for computing the gradient descent direction.
template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionGradient : public DirectionBase<ORDER, mydim, ndim, T>{
public:
    //! A delegating constructor.
    DirectionGradient(const T& fp):
            DirectionBase<ORDER, mydim, ndim, T>(fp){};
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const override;
    //! A method to compute the gradient descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters. In the gradient method they aren't.
    void resetParameters() override {};

};


//! @brief A class for computing the BFGS descent direction.
template<UInt ORDER, UInt mydim, UInt ndim, class T = FunctionalProblem<ORDER, mydim, ndim>>
class DirectionBFGS : public DirectionBase<ORDER, mydim, ndim, T>{
private:
    MatrixXr HInit_, HOld_;
    VectorXr gOld_, gradOld_;
    bool updateH_;

public:
    //! A constructor.
    DirectionBFGS(const T& fp, UInt k):
            DirectionBase<ORDER, mydim, ndim, T>(fp), HInit_(MatrixXr::Identity(k, k)), HOld_(MatrixXr::Identity(k, k)), updateH_(false){};
    //! A copy constructor: it just creates a DirectionBFGS object with the same features of rhs but it initializes the matrices.
    DirectionBFGS(const DirectionBFGS<ORDER, mydim, ndim, T>& rhs);
    //! Clone method overridden.
    std::unique_ptr<DirectionBase<ORDER, mydim, ndim, T>> clone() const override;
    //! A method to compute the BFGS descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters.
    void resetParameters() override;

};

//! ####################################################################################################################
//! ######################################## SPACE-TIME PROBLEM ########################################################
//! ####################################################################################################################
/*
template<UInt ORDER, UInt mydim, UInt ndim>
class DirectionBase_time{
protected:
    // to give generality if you want to add other children
    const FunctionalProblem_time<ORDER, mydim, ndim>& funcProblem_;

public:
    //! A constructor
    DirectionBase_time(const FunctionalProblem_time<ORDER, mydim, ndim>& fp): funcProblem_(fp){};
    //! A destructor.
    virtual ~DirectionBase_time(){};
    //! A pure virtual clone method.
    virtual std::unique_ptr<DirectionBase_time<ORDER, mydim, ndim>> clone() const = 0;
    //! A pure virtual method to compute the descent direction.
    virtual VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) = 0;
    //! A pure virtual method to reset all the old parameters.
    virtual void resetParameters() = 0;
};


//! @brief A class for computing the gradient descent direction.
template<UInt ORDER, UInt mydim, UInt ndim>
class DirectionGradient_time : public DirectionBase_time<ORDER, mydim, ndim>{
public:
    //! A delegating constructor.
    DirectionGradient_time(const FunctionalProblem_time<ORDER, mydim, ndim>& fp):
            DirectionBase_time<ORDER, mydim, ndim>(fp){};
    //! Clone method overridden.
    std::unique_ptr<DirectionBase_time<ORDER, mydim, ndim>> clone() const override;
    //! A method to compute the gradient descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters. In the gradient method they aren't.
    void resetParameters() override {};

};

//! @brief A class for computing the BFGS descent direction.
template<UInt ORDER, UInt mydim, UInt ndim>
class DirectionBFGS_time : public DirectionBase_time<ORDER, mydim, ndim>{
private:
    MatrixXr HInit_, HOld_;
    VectorXr gOld_, gradOld_;
    bool updateH_;

public:
    //! A constructor.
    DirectionBFGS_time(const FunctionalProblem_time<ORDER, mydim, ndim>& fp, UInt k):
            DirectionBase_time<ORDER, mydim, ndim>(fp), HInit_(MatrixXr::Identity(k, k)), HOld_(MatrixXr::Identity(k, k)), updateH_(false){};
    //! A copy constructor: it just creates a DirectionBFGS object with the same features of rhs but it initializes the matrices.
    DirectionBFGS_time(const DirectionBFGS_time<ORDER, mydim, ndim>& rhs);
    //! Clone method overridden.
    std::unique_ptr<DirectionBase_time<ORDER, mydim, ndim>> clone() const override;
    //! A method to compute the BFGS descent direction.
    VectorXr computeDirection(const VectorXr& g, const VectorXr& grad) override;
    //! A method to reset all the old parameters.
    void resetParameters() override;

};
*/
#include "Descent_Direction_imp.h"

#endif //DEV_FDAPDE_DESCENT_DIRECTION_H

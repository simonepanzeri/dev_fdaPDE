//
// Created by simonepanzeri on 02/12/2021.
//

#ifndef DEV_FDAPDE_DESCENT_DIRECTION_FACTORY_H
#define DEV_FDAPDE_DESCENT_DIRECTION_FACTORY_H

#include "FdaPDE.h"

//template<typename T, typename... Args>
//std::unique_ptr<T> make_unique(Args&&... args)
//{
//    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
//}

//! brief@ A Factory class: a class for the choice of the step mehod for the optimization algorithm.
template<UInt ORDER, UInt mydim, UInt ndim>
class DescentDirection_factory
{
public:
    //! A method that builds a pointer to the right object for the direction choice, taking as parameters a string and others objects needed for constructor.
    static std::unique_ptr<DirectionBase<ORDER,  mydim,  ndim>>
    createDirectionSolver(const DataProblem<ORDER, mydim, ndim>& dp,
    const FunctionalProblem<ORDER, mydim, ndim>& fp, const std::string& d)
    {
        if (d=="Gradient")
            return make_unique<DirectionGradient<ORDER,mydim,ndim>>(fp);
        else if (d=="BFGS")
            return make_unique<DirectionBFGS<ORDER,mydim,ndim>>(fp, dp.getNumNodes());
        else{

           //  Rprintf("Unknown direction option - using gradient direction");

            return make_unique<DirectionGradient<ORDER,mydim,ndim>>(fp);
        }
    }

};

//! ####################################################################################################################
//! ######################################## SPACE-TIME PROBLEM ########################################################
//! ####################################################################################################################

template<UInt ORDER, UInt mydim, UInt ndim>
class DescentDirection_factory_time{
public:
    //! A method that builds a pointer to the right object for the direction choice, taking as parameters a string and others objects needed for constructor.
    static std::unique_ptr<DirectionBase_time<ORDER,  mydim,  ndim>>
    createDirectionSolver(const DataProblem_time<ORDER, mydim, ndim>& dp,
                          const FunctionalProblem_time<ORDER, mydim, ndim>& fp, const std::string& d)
    {
        if (d=="Gradient")
            return make_unique<DirectionGradient_time<ORDER,mydim,ndim>>(fp);
        //else if (d=="BFGS")
        //    return make_unique<DirectionBFGS<ORDER,mydim,ndim>>(fp, dp.getNumNodes());
        else{

            //  Rprintf("Unknown direction option - using gradient direction");

            return make_unique<DirectionGradient_time<ORDER,mydim,ndim>>(fp);
        }
    }

};


#endif //DEV_FDAPDE_DESCENT_DIRECTION_FACTORY_H

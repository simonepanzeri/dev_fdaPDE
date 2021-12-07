//
// Created by simonepanzeri on 01/12/2021.
//

#ifndef DEV_FDAPDE_DENSITY_INITIALIZATION_FACTORY_H
#define DEV_FDAPDE_DENSITY_INITIALIZATION_FACTORY_H

#include "FdaPDE.h"
#include "Preprocess_Phase.h"
#include "Preprocess_Phase.h"

//template<typename T, typename... Args>
//std::unique_ptr<T> make_unique(Args&&... args)
//{
//    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
//}

template<UInt ORDER, UInt mydim, UInt ndim>
class DensityInitialization_factory
{
public:
    //! A method that builds a pointer to the right object for the initialization choice.
    static std::unique_ptr<DensityInitialization<ORDER,  mydim,  ndim>>
    createInitializationSolver(const DataProblem<ORDER, mydim, ndim>& dp,
                               const FunctionalProblem<ORDER, mydim, ndim>& fp){

        if(!dp.isFvecEmpty())
            return make_unique<UserInitialization<ORDER, mydim, ndim>>(dp);
        else
            return make_unique<HeatProcess<ORDER, mydim, ndim>>(dp, fp);

    }

};

//! ####################################################################################################################
//! ######################################## SPACE-TIME PROBLEM ########################################################
//! ####################################################################################################################

template<UInt ORDER, UInt mydim, UInt ndim>
class DensityInitialization_factory_time{
public:
    //! A method that builds a pointer to the right object for the initialization choice.
    static std::unique_ptr<DensityInitialization_time<ORDER,  mydim,  ndim>>
    createInitializationSolver(const DataProblem_time<ORDER, mydim, ndim>& dp,
                               const FunctionalProblem_time<ORDER, mydim, ndim>& fp){

        //if(!dp.isFvecEmpty())
            return make_unique<UserInitialization_time<ORDER, mydim, ndim>>(dp);
        //else
        //    return make_unique<HeatProcess<ORDER, mydim, ndim>>(dp, fp);

    }

};
#endif //DEV_FDAPDE_DENSITY_INITIALIZATION_FACTORY_H

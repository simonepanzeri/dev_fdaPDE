//
// Created by simonepanzeri on 02/12/2021.
//

#ifndef DEV_FDAPDE_OPTIMIZATION_ALGORITHM_FACTORY_H
#define DEV_FDAPDE_OPTIMIZATION_ALGORITHM_FACTORY_H

#include "FdaPDE.h"

//template<typename T, typename... Args>
//std::shared_ptr<T> make_shared(Args&&... args)
//{
//    return std::make_shared<T>(new T(std::forward<Args>(args)...));
//}

//!brief@ A Factory class: a class for the choice of the step mehod for the optimization algorithm.
template<UInt ORDER, UInt mydim, UInt ndim>
class MinimizationAlgorithm_factory
{
public:
    //! A method that builds a pointer to the right object for the step choice, taking as parameters a string and others objects needed for constructor.
    static std::shared_ptr<MinimizationAlgorithm<ORDER,  mydim,  ndim>>
    createStepSolver(const DataProblem<ORDER, mydim, ndim>& dp,
    const FunctionalProblem<ORDER, mydim, ndim>& fp,
    const std::string& d, const std::string& s)
    {
        if(s == "Fixed_Step") return std::make_shared<FixedStep<ORDER, mydim, ndim>>(dp, fp, d);

        else if(s == "Backtracking_Method") return std::make_shared<BacktrackingMethod<ORDER, mydim, ndim>>(dp, fp, d);

        else if(s == "Wolfe_Method") return std::make_shared<WolfeMethod<ORDER, mydim, ndim>>(dp, fp, d);

        else{

            //Rprintf("Unknown step option - using fixed step\n");

            return std::make_shared<FixedStep<ORDER, mydim, ndim>>(dp, fp,  std::move(d));
        }
    }

};

//! ####################################################################################################################
//! ######################################## SPACE-TIME PROBLEM ########################################################
//! ####################################################################################################################

template<UInt ORDER, UInt mydim, UInt ndim>
class MinimizationAlgorithm_factory_time{
public:
    //! A method that builds a pointer to the right object for the step choice, taking as parameters a string and others objects needed for constructor.
    static std::shared_ptr<MinimizationAlgorithm_time<ORDER,  mydim,  ndim>>
    createStepSolver(const DataProblem_time<ORDER, mydim, ndim>& dp,
                     const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
                     const std::string& d, const std::string& s)
    {
        if(s == "Fixed_Step") {
            std::cout << "Fixed_Step" << std::endl;
            return std::make_shared<FixedStep_time<ORDER, mydim, ndim>>(dp, fp, d);
        }
        else if(s == "Backtracking_Method") {
            std::cout << "Backtracking_Method" << std::endl;
            return std::make_shared<BacktrackingMethod_time<ORDER, mydim, ndim>>(dp, fp, d);
        }
        else if(s == "Wolfe_Method") {
            std::cout << "Wolfe_Method" << std::endl;
            return std::make_shared<WolfeMethod_time<ORDER, mydim, ndim>>(dp, fp, d);
        }
        else{
            std::cout << "Unknown step option - using fixed step" << std::endl;
            //Rprintf("Unknown step option - using fixed step\n");

            return std::make_shared<FixedStep_time<ORDER, mydim, ndim>>(dp, fp,  std::move(d));
        }
    }

};


#endif //DEV_FDAPDE_OPTIMIZATION_ALGORITHM_FACTORY_H

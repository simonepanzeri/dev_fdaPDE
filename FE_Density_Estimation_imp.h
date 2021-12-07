//
// Created by simonepanzeri on 01/12/2021.
//

#ifndef DEV_FDAPDE_FE_DENSITY_ESTIMATION_IMP_H
#define DEV_FDAPDE_FE_DENSITY_ESTIMATION_IMP_H

template<UInt ORDER, UInt mydim, UInt ndim>
FEDE<ORDER, mydim, ndim>::
FEDE(const DataProblem<ORDER, mydim, ndim>& dp,
     const FunctionalProblem<ORDER, mydim, ndim>& fp,
     std::shared_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> ma, const std::string& p):
        dataProblem_(dp), funcProblem_(fp), minAlgo_(ma){

        preprocess_ = Preprocess_factory<ORDER, mydim, ndim>::createPreprocessSolver(dp, fp, ma, p);

}


template<UInt ORDER, UInt mydim, UInt ndim>
void
FEDE<ORDER, mydim, ndim>::apply(){

    // perform the preprocess phase
    //! Rprintf("##### PREPROCESS PHASE #####\n");
    //! preprocess_ -> performPreprocessTask();

    // collect preprocess results
    VectorXr gInit;
    //! std::tie(fInit_, gInit, bestLambda_) = preprocess_ -> getPreprocessParameter();

    //! CV_errors_ = preprocess_ -> getCvError();

    // final minimization descent
    //! Rprintf("##### FINAL STEP #####\n");

    gcoeff_ = minAlgo_->apply_core(dataProblem_.getGlobalPsi(), bestLambda_, gInit);

}

//! ####################################################################################################################
//! ######################################## SPACE-TIME PROBLEM ########################################################
//! ####################################################################################################################

template<UInt ORDER, UInt mydim, UInt ndim>
FEDE_time<ORDER, mydim, ndim>::
FEDE_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
     const FunctionalProblem_time<ORDER, mydim, ndim>& fp,
     std::shared_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> ma, const std::string& p):
        dataProblem_(dp), funcProblem_(fp), minAlgo_(ma){

    preprocess_ = Preprocess_factory_time<ORDER, mydim, ndim>::createPreprocessSolver(dp, fp, ma, p);

}


template<UInt ORDER, UInt mydim, UInt ndim>
void
FEDE_time<ORDER, mydim, ndim>::apply(){
    std::cout<<"==== PREPROCESS PHASE ===="<<std::endl;
    // perform the preprocess phase
    //! Rprintf("##### PREPROCESS PHASE #####\n");
    preprocess_ -> performPreprocessTask();

    std::cout<<"==== COLLECT PREPROCESS RESULT ===="<<std::endl;
    // collect preprocess results
    VectorXr gInit;
    std::tie(fInit_, gInit, bestLambda_S, bestLambda_T) = preprocess_ -> getPreprocessParameter();

    std::cout<<" dim of gInit vector: "<<gInit.size()<<std::endl;

    std::cout<<"==== CV PREPROCESS PHASE ===="<<std::endl;
    CV_errors_ = preprocess_ -> getCvError();

    // final minimization descent
    //! Rprintf("##### FINAL STEP #####\n");
    std::cout<<"==== FINAL STEP ===="<<std::endl;
    gcoeff_ = minAlgo_->apply_core(dataProblem_.getUpsilon(), bestLambda_S, bestLambda_T, gInit);

}

#endif //DEV_FDAPDE_FE_DENSITY_ESTIMATION_IMP_H

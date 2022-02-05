//
// Created by simonepanzeri on 02/12/2021.
//

#ifndef DEV_FDAPDE_OPTIMIZATION_ALGORITHM_IMP_H
#define DEV_FDAPDE_OPTIMIZATION_ALGORITHM_IMP_H

#include "Data_Generation.h"

template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm<ORDER, mydim, ndim>::
MinimizationAlgorithm(const DataProblem<ORDER, mydim, ndim>& dp,
                      const FunctionalProblem<ORDER, mydim, ndim>& fp, const std::string& d):
        dataProblem_(dp), funcProblem_(fp){

    direction_ = DescentDirection_factory<ORDER,  mydim,  ndim>::createDirectionSolver(dp, fp, d);

};


template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm<ORDER, mydim, ndim>::
MinimizationAlgorithm(const MinimizationAlgorithm<ORDER, mydim, ndim>& rhs):
        dataProblem_(rhs.dataProblem_), funcProblem_(rhs.funcProblem_){

    direction_ = rhs.direction_->clone();
};


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>>
FixedStep<ORDER, mydim, ndim>::clone() const {

    return make_unique<FixedStep<ORDER, mydim, ndim>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
FixedStep<ORDER,mydim,ndim>::apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const{

    // termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1 +1, dllik = toll1 +1, dpen = toll1 +1;

    // to save the current point
    VectorXr g_curr;

    // variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen, pen_old;
    UInt i;

    for(UInt e = 0; e < this->dataProblem_.getNstepProposals(); e++){
        // start always with the initial point
        g_curr = g;

        std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
        norm_grad = std::sqrt(grad.dot(grad));

///*************************************************************************
        std::vector<Real> data_loss, data_llik_1, data_llik_2, data_pen, data_norm_grad;
        data_loss.reserve(this->dataProblem_.getNsimulations()+1);
        data_llik_1.reserve(this->dataProblem_.getNsimulations()+1);
        data_llik_2.reserve(this->dataProblem_.getNsimulations()+1);
        data_pen.reserve(this->dataProblem_.getNsimulations()+1);
        data_norm_grad.reserve(this->dataProblem_.getNsimulations()+1);
        data_loss.push_back(loss);
        data_pen.push_back(pen);
        data_norm_grad.push_back(norm_grad);
        data_llik_1.push_back(-(Psi*g_curr).sum());
        data_llik_2.push_back(loss - data_llik_1[0] - lambda*pen);

        writeSolution<VectorXr>(g_curr, "../data/space/g_sol/g_sol_" + std::to_string(0) + ".txt");
///*************************************************************************

        if(this->dataProblem_.Print()){
            std::cout << "Loss: " << loss << ", llik: " << llik << ", pen: " << pen << ", norm_Lp: " << norm_grad << std::endl;
            //Rprintf("loss %f, llik %f, pen %f, norm_Lp %f\n", loss, llik, pen, norm_grad);
        }

        for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen > toll1) && norm_grad > toll2; i++){

            // Termination criteria variables
            loss_old = loss;
            llik_old = llik;
            pen_old = pen;

            // Compute a descent direction
            d = this->direction_->computeDirection(g_curr, grad);

            // Update the point
            g_curr = g_curr + this->dataProblem_.getStepProposals(e)*d;

            // Update termination criteria variables
            std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
            dloss = std::abs((loss - loss_old)/loss_old);
            dllik = std::abs((llik - llik_old)/llik_old);
            dpen = std::abs((pen - pen_old)/pen_old);
            norm_grad = std::sqrt(grad.dot(grad));

///*************************************************************************
            data_loss.push_back(loss);
            data_pen.push_back(pen);
            data_norm_grad.push_back(norm_grad);
            data_llik_1.push_back(-(Psi*g_curr).sum());
            data_llik_2.push_back(loss - data_llik_1.back() - lambda*pen);

            if ((i+1) < 10 || (i+1)%50 == 0)
                writeSolution<VectorXr>(g_curr, "../data/space/g_sol/g_sol_" + std::to_string(i+1) + ".txt");
///*************************************************************************

            // Check if the step is ok
            if((loss_old - loss) < 0){
                if(this->dataProblem_.Print()){
                    std::cout << "Iter: " << i+1 << ", loss: " << loss << ", llik: " << llik << ", pen: " << pen << ", norm_Lp: " << norm_grad << std::endl;
                    // Rprintf("The loss function increases: not good. Try decreasing the optimization parameter.\n", norm_grad);
                }
                break;
            }

            if(this->dataProblem_.Print()){
                std::cout << "Iter: " << i+1 << ", loss: " << loss << ", llik: " << llik << ", pen: " << pen << ", norm_Lp: " << norm_grad << std::endl;
                //Rprintf("Iter %d, loss %f, llik %f, pen %f, norm_Lp %f\n", i+1, loss, llik, pen, norm_grad);
            }
        }

        this->direction_->resetParameters();

        if ((loss_old - loss) < 0){
            std::cout << "Iter: " << i+1 << ", loss: " << loss << ", llik: " << llik << ", pen: " << pen << ", norm_Lp: " << norm_grad << std::endl;
            writeSolution<std::vector<Real>>(data_loss, "../data/space/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen, "../data/space/functional/data_pen.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda}, "../data/space/functional/lambda.txt");
        }
        else if(dloss <= toll1 && dllik <= toll1 && dpen <= toll1){
            writeSolution<std::vector<Real>>(data_loss, "../data/space/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen, "../data/space/functional/data_pen.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda}, "../data/space/functional/lambda.txt");
            if(this->dataProblem_.Print()){
                std::cout << "The algorithm reaches the tolerance in terms of the functional. Norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen: " << dpen << std::endl;
                //Rprintf("The algorithm reaches the tolerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
            }
            return g_curr;
        }
        else if(norm_grad <= toll2){
            writeSolution<std::vector<Real>>(data_loss, "../data/space/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen, "../data/space/functional/data_pen.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda}, "../data/space/functional/lambda.txt");
            if(this->dataProblem_.Print()){
                std::cout << "The algorithm reaches the tolerance in terms of the slope. Norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen: " << dpen << std::endl;
                //Rprintf("The algorithm reaches the tolerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
            }
            return g_curr;
        }
        else if(i == this->dataProblem_.getNsimulations()){
            writeSolution<std::vector<Real>>(data_loss, "../data/space/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen, "../data/space/functional/data_pen.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda}, "../data/space/functional/lambda.txt");
            if(this->dataProblem_.Print()){
                std::cout << "The algorithm reaches the maximum number of iterations. Norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen: " << dpen << std::endl;
                //Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
            }
            return g_curr;
        }
    }

    // If you arrive here you don't have a good gradient parameter
    std::cout << "ERROR: The loss function increases: not good. Try decreasing the optimization parameter" << std::endl;
    // Rprintf("ERROR: The loss function increases: not good. Try decreasing the optimization parameter");
    //std::abort();
    return VectorXr::Constant(g.size(),0);
}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
AdaptiveStep<ORDER,mydim,ndim>::apply_core(const SpMat& Psi, Real lambda, const VectorXr& g) const{

    // termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1 +1, dllik = toll1 +1, dpen = toll1 +1;

    // to save the current point
    VectorXr g_curr = g;

    // variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen, pen_old, step;
    UInt i;

    std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
    norm_grad = std::sqrt(grad.dot(grad));

    if(this->dataProblem_.Print()){
        // Rprintf("loss %f, llik %f, pen %f, norm_Lp %f\n", loss, llik, pen, norm_grad);
    }

    for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen > toll1) && norm_grad > toll2; i++){

        // Termination criteria variables
        loss_old = loss;
        llik_old = llik;
        pen_old = pen;

        // Compute a descent direction
        d = this->direction_->computeDirection(g_curr, grad);

        // Compute a step
        step = computeStep(g_curr, loss, grad, d, lambda, Psi);

        // Update the point
        g_curr = g_curr + step*d;

        // Update termination criteria variables
        std::tie(loss, grad, llik, pen) = this->funcProblem_.computeFunctional_g(g_curr, lambda, Psi);
        dloss = std::abs((loss - loss_old)/loss_old);
        dllik = std::abs((llik - llik_old)/llik_old);
        dpen = std::abs((pen - pen_old)/pen_old);
        norm_grad = std::sqrt(grad.dot(grad));

        if(this->dataProblem_.Print()){
            // Rprintf("Iter %d, loss %f, llik %f, pen %f, norm_Lp %f\n", i+1, loss, llik, pen, norm_grad);
        }

    }

    this->direction_->resetParameters();

    if(dloss <= toll1 && dllik <= toll1 && dpen <= toll1){
        if(this->dataProblem_.Print()){
          //  Rprintf("The algorithm reaches the tollerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
    else if(norm_grad <= toll2){
        if(this->dataProblem_.Print()){
          //  Rprintf("The algorithm reaches the tollerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
    else{
        if(this->dataProblem_.Print()){
         //   Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
}


template<UInt ORDER, UInt mydim, UInt ndim>
Real
BacktrackingMethod<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda, const SpMat& Psi) const{

    Real ro = 0.5, alpha = 1/ro, c = 0.5;

    Real loss_new, llik_new, pen_new, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir =  grad.dot(dir);

    do{
        // update step
        alpha *= ro;

        slope = c*alpha*(grad_dir);

        // Update the point
        new_point = g + alpha*dir;

        // functional in the new point
        std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);

    } while(loss_new > (loss + slope));

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>>
BacktrackingMethod<ORDER, mydim, ndim>::clone() const {

    return make_unique<BacktrackingMethod<ORDER, mydim, ndim>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim>
Real
WolfeMethod<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda, const SpMat& Psi) const{

    Real alpha = 1, alphamax = 0, alphamin = 0, c1 = 1e-4, c2 = 0.9;

    Real loss_new, llik_new, pen_new, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir = grad.dot(dir);
    slope = c1*alpha*grad_dir;


    // Update the point
    new_point = g + alpha*dir;

    // functional in the new point
    std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);

    bool again = true;

    while(again){

        again = false;

        while(loss_new > (loss + slope)){
            // update step
            alphamax = alpha;
            alpha = 0.5*(alphamin + alphamax);

            // try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);
            slope = c1*alpha*grad_dir;
        }

        if(grad_new.dot(dir) < c2*grad_dir && std::abs(grad_dir)>1e-2){

            again = true;

            // update step
            alphamin = alpha;
            alpha = alphamax==0 ? 2*alphamin : 0.5*(alphamin+alphamax);

            // try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new) = this->funcProblem_.computeFunctional_g(new_point, lambda, Psi);
            slope =  alpha*c1*grad_dir;
        }
    }

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>>
WolfeMethod<ORDER, mydim, ndim>::clone() const {

    return make_unique<WolfeMethod<ORDER, mydim, ndim>>(*this);

}

//! ####################################################################################################################
//! ######################################## SPACE-TIME PROBLEM ########################################################
//! ####################################################################################################################

template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm_time<ORDER, mydim, ndim>::
MinimizationAlgorithm_time(const DataProblem_time<ORDER, mydim, ndim>& dp,
                           const FunctionalProblem_time<ORDER, mydim, ndim>& fp, const std::string& d):
        dataProblem_(dp), funcProblem_(fp){

    direction_ = DescentDirection_factory_time<ORDER,  mydim,  ndim>::createDirectionSolver(dp, fp, d);

};


template<UInt ORDER, UInt mydim, UInt ndim>
MinimizationAlgorithm_time<ORDER, mydim, ndim>::
MinimizationAlgorithm_time(const MinimizationAlgorithm_time<ORDER, mydim, ndim>& rhs):
        dataProblem_(rhs.dataProblem_), funcProblem_(rhs.funcProblem_){

    direction_ = rhs.direction_->clone();
};


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>>
FixedStep_time<ORDER, mydim, ndim>::clone() const {

    return make_unique<FixedStep_time<ORDER, mydim, ndim>>(*this);

}


template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
FixedStep_time<ORDER,mydim,ndim>::apply_core(const SpMat& Upsilon, Real lambda_S, Real lambda_T,const VectorXr& g) const{

    //std::cout << "==== INTO MINIMIZATION APPLY CORE ====" << std::endl;
    // termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1+1, dllik = toll1+1, dpen_S = toll1+1, dpen_T = toll1+1;

    // to save the current point
    VectorXr g_curr;

    // variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen_S, pen_T, pen_old_S, pen_old_T;
    UInt i;

    for(UInt e = 0; e < this->dataProblem_.getNstepProposals(); e++){
        // start always with the initial point
        g_curr = g;

        //std::cout << "==== COMPUTATION OF FUNCTIONAL ====" << std::endl;
        std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
        norm_grad = std::sqrt(grad.dot(grad));
        //std::cout << "norm_grad: " << norm_grad<<std::endl;

///*************************************************************************
        std::vector<Real> data_loss, data_llik_1, data_llik_2, data_pen_S, data_pen_T, data_norm_grad;
        data_loss.reserve(this->dataProblem_.getNsimulations()+1);
        data_llik_1.reserve(this->dataProblem_.getNsimulations()+1);
        data_llik_2.reserve(this->dataProblem_.getNsimulations()+1);
        data_pen_S.reserve(this->dataProblem_.getNsimulations()+1);
        data_pen_T.reserve(this->dataProblem_.getNsimulations()+1);
        data_norm_grad.reserve(this->dataProblem_.getNsimulations()+1);
        data_loss.push_back(loss);
        data_pen_S.push_back(pen_S);
        data_pen_T.push_back(pen_T);
        data_norm_grad.push_back(norm_grad);
        data_llik_1.push_back(-(Upsilon*g_curr).sum());
        data_llik_2.push_back(loss - data_llik_1[0] - lambda_S*pen_S - lambda_T*pen_T);

        writeSolution<VectorXr>(g_curr, "../data/space_time/g_sol/g_sol_" + std::to_string(0) + ".txt");
///*************************************************************************

        if(this->dataProblem_.Print()){
            //std::cout << "Loss: " << loss << ", llik: " << llik << ", pen_S: " << pen_S << ", pen_T: " << pen_T << ", norm_Lp: " << norm_grad << std::endl;
            //Rprintf("loss %f, llik %f, pen_S %f, pen_T %f, norm_Lp %f\n", loss, llik, pen_S, pen_T, norm_grad);
        }

        for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen_S > toll1 || dpen_T > toll1) && norm_grad > toll2; i++){
            // Termination criteria variables
            loss_old = loss;
            llik_old = llik;
            pen_old_S = pen_S;
            pen_old_T = pen_T;

            // Compute a descent direction
            //std::cout << "==== COMPUTATION OF THE DIRECTION ====" << std::endl;
            d = this->direction_->computeDirection(g_curr, grad);

            // Update the point
            //std::cout << "==== POINT UPDATE ====" << std::endl;
            g_curr = g_curr + this->dataProblem_.getStepProposals(e)*d;

            // Update termination criteria variables
            //std::cout << "==== TERMINATION CRITERIA UPDATE ====" << std::endl;
            std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
            dloss = std::abs((loss - loss_old)/loss_old);
            dllik = std::abs((llik - llik_old)/llik_old);
            dpen_S = std::abs((pen_S - pen_old_S)/pen_old_S);
            dpen_T = std::abs((pen_T - pen_old_T)/pen_old_T);
            norm_grad = std::sqrt(grad.dot(grad));

///*************************************************************************
            data_loss.push_back(loss);
            data_pen_S.push_back(pen_S);
            data_pen_T.push_back(pen_T);
            data_norm_grad.push_back(norm_grad);
            data_llik_1.push_back(-(Upsilon*g_curr).sum());
            data_llik_2.push_back(loss - data_llik_1.back() - lambda_S*pen_S - lambda_T*pen_T);

            if ((i+1) < 10 || (i+1)%50 == 0)
                writeSolution<VectorXr>(g_curr, "../data/space_time/g_sol/g_sol_" + std::to_string(i+1) + ".txt");
///*************************************************************************

            // Check if the step is ok
            if((loss_old - loss) < 0){
                if(this->dataProblem_.Print()){
                    std::cout << "The loss function increases: not good. Try decreasing the optimization parameter " << std::endl;
                    //Rprintf("The loss function increases: not good. Try decreasing the optimization parameter.\n", norm_grad);
                }
                break;
            }

            if(this->dataProblem_.Print()){
                //std::cout << "Iter: " << i+1 << ", loss: " << loss << ", llik: " << llik << ", pen_S: " << pen_S << ", pen_T: " << pen_T << ", norm_Lp: " << norm_grad << std::endl;
                //Rprintf("Iter %d, loss %f, llik %f, pen_S %f, pen_T %f, norm_Lp %f\n", i+1, loss, llik, pen_S, pen_T, norm_grad);
            }
        }

        this->direction_->resetParameters();

        if ((loss_old - loss) < 0){
            //std::cout << "Iter: " << i+1 << ", loss: " << loss << ", llik: " << llik << ", pen_S: " << pen_S << ", pen_T: " << pen_T << ", norm_Lp: " << norm_grad << std::endl;
            writeSolution<std::vector<Real>>(data_loss, "../data/space_time/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen_S, "../data/space_time/functional/data_pen_S.txt");
            writeSolution<std::vector<Real>>(data_pen_T, "../data/space_time/functional/data_pen_T.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space_time/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space_time/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space_time/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda_S, lambda_T}, "../data/space_time/functional/lambdas.txt");
        }
        else if(dloss <= toll1 && dllik <= toll1 && dpen_S <= toll1 && dpen_T <= toll1){
///*************************************************************************
            writeSolution<std::vector<Real>>(data_loss, "../data/space_time/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen_S, "../data/space_time/functional/data_pen_S.txt");
            writeSolution<std::vector<Real>>(data_pen_T, "../data/space_time/functional/data_pen_T.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space_time/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space_time/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space_time/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda_S, lambda_T}, "../data/space_time/functional/lambdas.txt");
///*************************************************************************
            if(this->dataProblem_.Print()){
                std::cout << "The algorithm reaches the tolerance in terms of the functional." << std::endl << "Iter: " << i+1 << ", norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen_S: " << dpen_S << ", dpen_T: " << dpen_T << std::endl;
                //Rprintf("The algorithm reaches the tolerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
            }
            return g_curr;
        }
        else if(norm_grad <= toll2){
///*************************************************************************
            writeSolution<std::vector<Real>>(data_loss, "../data/space_time/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen_S, "../data/space_time/functional/data_pen_S.txt");
            writeSolution<std::vector<Real>>(data_pen_T, "../data/space_time/functional/data_pen_T.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space_time/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space_time/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space_time/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda_S, lambda_T}, "../data/space_time/functional/lambdas.txt");
///*************************************************************************
            if(this->dataProblem_.Print()){
                std::cout << "The algorithm reaches the tolerance in terms of the slope." << std::endl << "Iter: " << i+1 << ", norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen_S: " << dpen_S << ", dpen_T: " << dpen_T << std::endl;
                //Rprintf("The algorithm reaches the tolerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
            }
            return g_curr;
        }
        else if(i == this->dataProblem_.getNsimulations()){
///*************************************************************************
            writeSolution<std::vector<Real>>(data_loss, "../data/space_time/functional/data_loss.txt");
            writeSolution<std::vector<Real>>(data_pen_S, "../data/space_time/functional/data_pen_S.txt");
            writeSolution<std::vector<Real>>(data_pen_T, "../data/space_time/functional/data_pen_T.txt");
            writeSolution<std::vector<Real>>(data_norm_grad, "../data/space_time/functional/data_norm_grad.txt");
            writeSolution<std::vector<Real>>(data_llik_1, "../data/space_time/functional/data_llik_1.txt");
            writeSolution<std::vector<Real>>(data_llik_2, "../data/space_time/functional/data_llik_2.txt");
            writeSolution<std::vector<Real>>({lambda_S, lambda_T}, "../data/space_time/functional/lambdas.txt");
///*************************************************************************
            if(this->dataProblem_.Print()){
                std::cout << "The algorithm reaches the maximum number of iterations." << std::endl << "Iter: " << i+1 << ", norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen_S: " << dpen_S << ", dpen_T: " << dpen_T << std::endl;
                //Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen_S: %f, dpen_T: %f\n", norm_grad, dloss, dllik, dpen_S, dpen_T);
            }
            return g_curr;
        }
    }

    // If you arrive here you don't have a good gradient parameter
    std::cout<<"ERROR: The loss function increases: not good. Try decreasing the optimization parameter"<<std::endl;
    // Rprintf("ERROR: The loss function increases: not good. Try decreasing the optimization parameter");
    //std::abort();
    return VectorXr::Constant(g.size(),0);

}

template<UInt ORDER, UInt mydim, UInt ndim>
VectorXr
AdaptiveStep_time<ORDER,mydim,ndim>::apply_core(const SpMat& Upsilon, Real lambda_S, Real lambda_T, const VectorXr& g) const{

    //std::cout << "==== INTO MINIMIZATION APPLY CORE ====" << std::endl;
    // termination criteria variables
    const Real toll1 = this->dataProblem_.getTol1(), toll2 = this->dataProblem_.getTol2();
    Real norm_grad, dloss = toll1+1, dllik = toll1+1, dpen_S = toll1+1, dpen_T = toll1+1;

    // to save the current point
    VectorXr g_curr = g;

    // variables
    VectorXr grad, d;
    Real loss, loss_old, llik, llik_old, pen_S, pen_T, pen_old_S, pen_old_T, step;
    UInt i;

    //std::cout << "==== COMPUTATION OF FUNCTIONAL ====" << std::endl;
    std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
    norm_grad = std::sqrt(grad.dot(grad));
    //std::cout << "norm_grad: " << norm_grad << std::endl;

///*************************************************************************
    std::vector<Real> data_loss, data_llik_1, data_llik_2, data_pen_S, data_pen_T, data_norm_grad;
    data_loss.reserve(this->dataProblem_.getNsimulations()+1);
    data_llik_1.reserve(this->dataProblem_.getNsimulations()+1);
    data_llik_2.reserve(this->dataProblem_.getNsimulations()+1);
    data_pen_S.reserve(this->dataProblem_.getNsimulations()+1);
    data_pen_T.reserve(this->dataProblem_.getNsimulations()+1);
    data_norm_grad.reserve(this->dataProblem_.getNsimulations()+1);
    data_loss.push_back(loss);
    data_pen_S.push_back(pen_S);
    data_pen_T.push_back(pen_T);
    data_norm_grad.push_back(norm_grad);
    data_llik_1.push_back(-(Upsilon*g_curr).sum());
    data_llik_2.push_back(loss - data_llik_1[0] - lambda_S*pen_S - lambda_T*pen_T);

    writeSolution<VectorXr>(g_curr, "../data/space_time/g_sol/g_sol_" + std::to_string(0) + ".txt");
///*************************************************************************

    if(this->dataProblem_.Print()){
        //std::cout << "Loss: " << loss << ", llik: " << llik << ", pen_S: " << pen_S << ", pen_T: " << pen_T << ", norm_Lp: " << norm_grad << std::endl;
        // Rprintf("loss %f, llik %f, pen %f, norm_Lp %f\n", loss, llik, pen, norm_grad);
    }

    for(i = 0; i < this->dataProblem_.getNsimulations() && (dloss > toll1 || dllik > toll1 || dpen_S > toll1 || dpen_T > toll1) && norm_grad > toll2; i++){

        // Termination criteria variables
        loss_old = loss;
        llik_old = llik;
        pen_old_S = pen_S;
        pen_old_T = pen_T;

        // Compute a descent direction
        //std::cout << "==== COMPUTATION OF THE DIRECTION ====" << std::endl;
        d = this->direction_->computeDirection(g_curr, grad);

        // Compute a step
        //std::cout << "==== POINT UPDATE ====" << std::endl;
        step = computeStep(g_curr, loss, grad, d, lambda_S, lambda_T, Upsilon);

        // Update the point
        g_curr = g_curr + step*d;

        // Update termination criteria variables
        //std::cout << "==== TERMINATION CRITERIA UPDATE ====" << std::endl;
        std::tie(loss, grad, llik, pen_S, pen_T) = this->funcProblem_.computeFunctional_g(g_curr, lambda_S, lambda_T, Upsilon);
        dloss = std::abs((loss - loss_old)/loss_old);
        dllik = std::abs((llik - llik_old)/llik_old);
        dpen_S = std::abs((pen_S - pen_old_S)/pen_old_S);
        dpen_T = std::abs((pen_T - pen_old_T)/pen_old_T);
        norm_grad = std::sqrt(grad.dot(grad));
///*************************************************************************
        data_loss.push_back(loss);
        data_pen_S.push_back(pen_S);
        data_pen_T.push_back(pen_T);
        data_norm_grad.push_back(norm_grad);
        data_llik_1.push_back(-(Upsilon*g_curr).sum());
        data_llik_2.push_back(loss - data_llik_1.back() - lambda_S*pen_S - lambda_T*pen_T);

        if ((i+1) < 10 || (i+1)%50 == 0)
            writeSolution<VectorXr>(g_curr, "../data/space_time/g_sol/g_sol_" + std::to_string(i+1) + ".txt");
///*************************************************************************

        if(this->dataProblem_.Print()){
            //std::cout << "Iter: " << i+1 << ", loss: " << loss << ", llik: " << llik << ", pen_S: " << pen_S << ", pen_T: " << pen_T << ", norm_Lp: " << norm_grad << std::endl;
            // Rprintf("Iter %d, loss %f, llik %f, pen %f, norm_Lp %f\n", i+1, loss, llik, pen, norm_grad);
        }

    }

    this->direction_->resetParameters();

    if(dloss <= toll1 && dllik <= toll1 && dpen_S <= toll1 && dpen_T <= toll1){
///*************************************************************************
        writeSolution<std::vector<Real>>(data_loss, "../data/space_time/functional/data_loss.txt");
        writeSolution<std::vector<Real>>(data_pen_S, "../data/space_time/functional/data_pen_S.txt");
        writeSolution<std::vector<Real>>(data_pen_T, "../data/space_time/functional/data_pen_T.txt");
        writeSolution<std::vector<Real>>(data_norm_grad, "../data/space_time/functional/data_norm_grad.txt");
        writeSolution<std::vector<Real>>(data_llik_1, "../data/space_time/functional/data_llik_1.txt");
        writeSolution<std::vector<Real>>(data_llik_2, "../data/space_time/functional/data_llik_2.txt");
        writeSolution<std::vector<Real>>({lambda_S, lambda_T}, "../data/space_time/functional/lambdas.txt");
///*************************************************************************
        if(this->dataProblem_.Print()){
            std::cout << "The algorithm reaches the tolerance in terms of the functional." << std::endl << "Iter: " << i+1 << ", norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen_S: " << dpen_S << ", dpen_T: " << dpen_T << std::endl;
            //Rprintf("The algorithm reaches the tolerance in terms of the functional. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
    else if(norm_grad <= toll2){
///*************************************************************************
        writeSolution<std::vector<Real>>(data_loss, "../data/space_time/functional/data_loss.txt");
        writeSolution<std::vector<Real>>(data_pen_S, "../data/space_time/functional/data_pen_S.txt");
        writeSolution<std::vector<Real>>(data_pen_T, "../data/space_time/functional/data_pen_T.txt");
        writeSolution<std::vector<Real>>(data_norm_grad, "../data/space_time/functional/data_norm_grad.txt");
        writeSolution<std::vector<Real>>(data_llik_1, "../data/space_time/functional/data_llik_1.txt");
        writeSolution<std::vector<Real>>(data_llik_2, "../data/space_time/functional/data_llik_2.txt");
        writeSolution<std::vector<Real>>({lambda_S, lambda_T}, "../data/space_time/functional/lambdas.txt");
///*************************************************************************
        if(this->dataProblem_.Print()){
            std::cout << "The algorithm reaches the tolerance in terms of the slope." << std::endl << "Iter: " << i+1 << ", norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen_S: " << dpen_S << ", dpen_T: " << dpen_T << std::endl;
            //Rprintf("The algorithm reaches the tolerance in terms of the slope. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
    else{
///*************************************************************************
        writeSolution<std::vector<Real>>(data_loss, "../data/space_time/functional/data_loss.txt");
        writeSolution<std::vector<Real>>(data_pen_S, "../data/space_time/functional/data_pen_S.txt");
        writeSolution<std::vector<Real>>(data_pen_T, "../data/space_time/functional/data_pen_T.txt");
        writeSolution<std::vector<Real>>(data_norm_grad, "../data/space_time/functional/data_norm_grad.txt");
        writeSolution<std::vector<Real>>(data_llik_1, "../data/space_time/functional/data_llik_1.txt");
        writeSolution<std::vector<Real>>(data_llik_2, "../data/space_time/functional/data_llik_2.txt");
        writeSolution<std::vector<Real>>({lambda_S, lambda_T}, "../data/space_time/functional/lambdas.txt");
///*************************************************************************
        if(this->dataProblem_.Print()){
            std::cout << "The algorithm reaches the maximum number of iterations." << std::endl << "Iter: " << i+1 << ", norm of Lp: " << norm_grad << ", dloss: " << dloss << ", dllik: " << dllik << ", dpen_S: " << dpen_S << ", dpen_T: " << dpen_T << std::endl;
            //Rprintf("The algorithm reaches the maximum number of iterations. Norm of Lp: %f, dloss: %f, dllik: %f, dpen: %f\n", norm_grad, dloss, dllik, dpen);
        }
        return g_curr;
    }
}


template<UInt ORDER, UInt mydim, UInt ndim>
Real
BacktrackingMethod_time<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda_S, Real lambda_T, const SpMat& Upsilon) const{

    Real ro = 0.5, alpha = 1/ro, c = 0.5;

    Real loss_new, llik_new, pen_new_S, pen_new_T, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir = grad.dot(dir);

    do{
        // update step
        alpha *= ro;

        slope = c*alpha*(grad_dir);

        // Update the point
        new_point = g + alpha*dir;

        // functional in the new point
        std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);

    } while(loss_new > (loss + slope));

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>>
BacktrackingMethod_time<ORDER, mydim, ndim>::clone() const {

    return make_unique<BacktrackingMethod_time<ORDER, mydim, ndim>>(*this);

}

template<UInt ORDER, UInt mydim, UInt ndim>
Real
WolfeMethod_time<ORDER,mydim,ndim>::computeStep(const VectorXr& g, Real loss, const VectorXr& grad, const VectorXr& dir, Real lambda_S, Real lambda_T, const SpMat& Upsilon) const{

    Real alpha = 1, alphamax = 0, alphamin = 0, c1 = 1e-4, c2 = 0.9;

    Real loss_new, llik_new, pen_new_S, pen_new_T, slope, grad_dir;
    VectorXr grad_new, new_point;

    grad_dir = grad.dot(dir);
    slope = c1*alpha*grad_dir;


    // Update the point
    new_point = g + alpha*dir;

    // functional in the new point
    std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);

    bool again = true;

    while(again){

        again = false;

        while(loss_new > (loss + slope)){
            // update step
            alphamax = alpha;
            alpha = 0.5*(alphamin + alphamax);

            // try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);
            slope = c1*alpha*grad_dir;
        }

        if(grad_new.dot(dir) < c2*grad_dir && std::abs(grad_dir)>1e-2){

            again = true;

            // update step
            alphamin = alpha;
            alpha = alphamax==0 ? 2*alphamin : 0.5*(alphamin+alphamax);

            // try with the new point
            new_point = g + alpha*dir;
            std::tie(loss_new, grad_new, llik_new, pen_new_S, pen_new_T) = this->funcProblem_.computeFunctional_g(new_point, lambda_S, lambda_T, Upsilon);
            slope =  alpha*c1*grad_dir;
        }
    }

    return alpha;
}


template<UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>>
WolfeMethod_time<ORDER, mydim, ndim>::clone() const {

    return make_unique<WolfeMethod_time<ORDER, mydim, ndim>>(*this);

}

#endif //DEV_FDAPDE_OPTIMIZATION_ALGORITHM_IMP_H

#include "RObjects.h"
#include "Functional_Problem.h"
#include "FE_Density_Estimation.h"
#include "Optimization_Algorithm_Factory.h"

#include "Data_Generation.h"

int main() {
    //! DATA GENERATION
    // Template parameters
    const UInt ORDER = 1;
    const UInt mydim = 2;
    const UInt ndim = 2;

    const UInt spline_degree = 3;

    // Spatio-temporal data
    std::vector<Point<ndim>> data_locations;
    std::vector<Real> data_time;
    std::tie(data_locations, data_time) = createSTLocations<ndim>();

    // Mesh
    std::vector<Real> points_v = readMesh<Real>("../data/space_time/mesh/nodes.txt");
    RNumericMatrix points_m(points_v.data(), points_v.size()/ndim, ndim);

    std::vector<UInt> sides_v = readMesh<UInt>("../data/space_time/mesh/sides.txt");
    for (size_t i = 0; i < sides_v.size(); ++i)
        --sides_v[i];
    RIntegerMatrix sides_m(sides_v.data(), sides_v.size()/2, 2);

    std::vector<UInt> elements_v = readMesh<UInt>("../data/space_time/mesh/triangles.txt");
    for (size_t i = 0; i < elements_v.size(); ++i)
        --elements_v[i];
    RIntegerMatrix elements_m(elements_v.data(), elements_v.size()/3, 3);

    std::vector<UInt> neighbors_v = readMesh<UInt>("../data/space_time/mesh/neighbors.txt");
    for (size_t i = 0; i < neighbors_v.size(); ++i) {
        if (neighbors_v[i] != -1)
            --neighbors_v[i];
    }
    RIntegerMatrix neighbors_m(neighbors_v.data(),neighbors_v.size()/3,3);

    std::vector<Real> mesh_time = readMesh<Real>("../data/space_time/mesh/mesh_time.txt");

    // Parameters
    const std::vector<Real> lambda{1e-1, 1e-2};
    const std::vector<Real> lambda_time{1e-2, 1e-3};
    const std::vector<Real> stepProposals{0.001};

    const Real tol1 = 0.00001, tol2 = 0.;
    const UInt nsim = 5000;

    //VectorXr fvec = VectorXr::Ones(points_v.size()/ndim * (mesh_time.size()+spline_degree-1));
    //fvec = std::exp(1)*fvec;
    VectorXr fvec;

    //std::string direction_method = "Gradient";
    std::string direction_method = "BFGS";

    std::string step_method = "Fixed_Step";
    //std::string step_method = "Backtracking_Method";
    //std::string step_method = "Wolfe_Method";

    //std::string preprocess_method = "NoCrossValidation";
    //std::string preprocess_method = "RightCV";
    std::string preprocess_method = "SimplifiedCV";

    UInt n_folds;
    if (preprocess_method == "SimplifiedCV")
        n_folds = lambda.size() * lambda_time.size();
    else
        n_folds = 7;

    //! SPATIO-TEMPORAL DENSITY ESTIMATION
/*
    // Density estimation data
    DEData<ndim> deData(data_locations, 2, fvec, 0.023, 100, lambda, n_folds, nsim,
                        stepProposals, tol1, tol2, false, 1);

    // Density estimation data time
    DEData_time deData_time(data_time, lambda_time);
    deData_time.setTimes2Locations();
    deData_time.printTimes2Locations(std::cout);
*/
    // Data problem time
    DataProblem_time<ORDER, mydim, ndim> dataProblem_time(data_locations, data_time, 2, fvec, 0.023,
                                                          500, lambda, lambda_time, n_folds, nsim, stepProposals,
                                                          tol1, tol2, true, 1, points_m, sides_m,
                                                          elements_m, neighbors_m, mesh_time,
                                                          1,0,0,0);
    std::cout << "DataProblem_time DONE" << std::endl;
/*
    for (auto it = dataProblem_time.data().cbegin(); it != dataProblem_time.data().cend(); ++it) {
        std::cout << "time: " << dataProblem_time.data_time(it - dataProblem_time.data().cbegin());
        std::cout << '\t' << "id ";
        operator << (std::cout, *it);
    }
*/
    // Functional problem time
    FunctionalProblem_time<ORDER, mydim, ndim> functionalProblem_time(dataProblem_time);
/*
    VectorXr g = VectorXr::Ones(dataProblem_time.getNumNodes() * dataProblem_time.getSplineNumber());
    std::tuple<Real, VectorXr, Real, Real, Real> result = functionalProblem_time.computeFunctional_g(g, 0.01, 0.01, dataProblem_time.getUpsilon());
    std::cout << "Functional problem time (g): " << std::endl;
    std::cout << "Penalized log-likelihood: " << std::get<0>(result) << std::endl;
    std::cout << "Log-likelihood: " << std::get<2>(result) << std::endl;
    std::cout << "Spatial penalization: " << std::get<3>(result) << std::endl;
    std::cout << "Temporal penalization: " << std::get<4>(result) << std::endl;
*/
    std::cout << "FunctionalProblem_time DONE" << std::endl;

    // Minimization algorithm time
    std::shared_ptr<MinimizationAlgorithm_time<ORDER, mydim, ndim>> minimizationAlgo_time =
            MinimizationAlgorithm_factory_time<ORDER, mydim, ndim>::createStepSolver(dataProblem_time, functionalProblem_time,
                                                                                     direction_method, step_method);
    std::cout << "MinimizationAlgorithm_time DONE" << std::endl;

    // Finite element density estimation time
    FEDE_time<ORDER, mydim, ndim> fede_time(dataProblem_time, functionalProblem_time, minimizationAlgo_time,
                                            preprocess_method);
    std::cout << "FEDE_time DONE" << std::endl;

    // Perform the whole task
    fede_time.apply();
    std::cout << "fede_time.apply() DONE" << std::endl;

    // Collect results
    VectorXr g_sol = fede_time.getDensity_g();
    //std::vector<const VectorXr*> f_init = fede_time.getInitialDensity();
    Real lambda_sol_S = fede_time.getBestLambda_S();
    Real lambda_sol_T = fede_time.getBestLambda_T();
    std::vector<Real> CV_errors = fede_time.getCvError();
    //const std::vector<Point<ndim>>& data = dataProblem_time.data();
    //const std::vector<Real>& data_t = dataProblem_time.data_time();

    std::cout << std::endl << std::endl;
    std::cout << "==== SUMMARY ====" << std::endl;
    std::cout << "Best lambda_S: " << lambda_sol_S << std::endl;
    std::cout << "Best lambda_T: " << lambda_sol_T << std::endl;
    if (preprocess_method != "NoCrossValidation") {
        std::cout << "CV_errors: ";
        for (Real i: CV_errors)
            std::cout << i << " ";
        std::cout << std::endl;
    }

    UInt ns = points_v.size()/ndim;
    UInt M = dataProblem_time.getSplineNumber();

    //std::vector<Real> t{0, M_PI/8, M_PI/4, 3*M_PI/8, M_PI/2, 5*M_PI/8, 3*M_PI/4, 7*M_PI/8, M_PI};
    std::vector<Real> t{mesh_time};

    VectorXr result = createSTSolution(ns, M, mesh_time, t, g_sol);

    // Copy results in a .txt file
    writeSolution<VectorXr>(result.array().exp(), "../data/space_time/solution/solution_st.txt");
/*
    //! How the solution varies iteration after iteration
    for (UInt iter = 0; iter <= nsim; ++iter) {
        std::ifstream ist("../data/space_time/g_sol/g_sol_"+std::to_string(iter)+".txt");
        if(ist) {
            VectorXr sol;
            sol.resize(g_sol.size());
            Real val;
            for (UInt i = 0; i < sol.size(); ++i) {
                ist >> val;
                sol[i] = val;
            }
            VectorXr res = createSTSolution(ns, M, mesh_time, t, sol);
            writeSolution<VectorXr>(res.array().exp(), "../data/space_time/solution/solution_st_"+std::to_string(iter)+".txt");
            ist.close();
        }
    }
*/
    return 0;
}

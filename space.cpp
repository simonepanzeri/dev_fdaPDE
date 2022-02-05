#include "RObjects.h"
#include "Functional_Problem.h"
#include "FE_Density_Estimation.h"
#include "Optimization_Algorithm_Factory.h"

#include "Data_Generation.h"

int main() {

    //! DATA GENERATION
    // Template parameters
    const UInt ORDER = 1;
    const UInt mydim = 3;
    const UInt ndim = 3;

    // Spatial data
    std::vector<Point<ndim>> data = createLocations<ndim>();

    // Mesh
    std::vector<Real> points_v = readMesh<Real>("../data/space/mesh/nodes_3D.txt");
    RNumericMatrix points_m(points_v.data(), points_v.size()/ndim, ndim);

    std::vector<UInt> sides_v = readMesh<UInt>("../data/space/mesh/faces.txt");
    for (size_t i = 0; i < sides_v.size(); ++i)
        --sides_v[i];
    RIntegerMatrix sides_m(sides_v.data(), sides_v.size()/3, 3);

    std::vector<UInt> elements_v = readMesh<UInt>("../data/space/mesh/tetrahedrons.txt");
    for (size_t i = 0; i < elements_v.size(); ++i)
        --elements_v[i];
    RIntegerMatrix elements_m(elements_v.data(), elements_v.size()/4, 4);

    std::vector<UInt> neighbors_v = readMesh<UInt>("../data/space/mesh/neighbors_3D.txt");
    for (size_t i = 0; i < neighbors_v.size(); ++i) {
        if (neighbors_v[i] != -1)
            --neighbors_v[i];
    }
    RIntegerMatrix neighbors_m(neighbors_v.data(),neighbors_v.size()/4,4);

    // Parameters
    const std::vector<Real> lambda{0.01};
    const std::vector<Real> stepProposals{0.001};
    const Real tol1 = 0.00001, tol2 = 0.;
    const UInt nsim = 5000;

    VectorXr fvec = VectorXr::Ones(points_v.size()/ndim);
    fvec = std::exp(1)*fvec;

    std::string direction_method = "Gradient";
    std::string step_method = "Fixed_Step";
    std::string preprocess_method = "NoCrossValidation";

    //! DENSITY ESTIMATION
    // Density estimation data
    DEData<ndim> deData(data, 2, fvec, 0.023, 100, lambda, 2, nsim, stepProposals,
                        tol1, tol2, false, 1);

    // Data problem
    DataProblem<ORDER, mydim, ndim> dataProblem(data, 2, fvec, 0.023, 100, lambda, 2,
                                                nsim, stepProposals, tol1, tol2, false, 1,
                                                points_m, sides_m, elements_m, neighbors_m);

    // Functional problem
    FunctionalProblem<ORDER, mydim, ndim> functionalProblem(dataProblem);

    VectorXr g = VectorXr::Ones(dataProblem.getNumNodes());
    std::tuple<Real, VectorXr, Real, Real> result = functionalProblem.computeFunctional_g(g, 0.01, dataProblem.getGlobalPsi());
    std::cout << "Functional problem: " << std::endl;
    std::cout << "Penalized log-likelihood: " << std::get<0>(result) << std::endl;
    std::cout << "Log-likelihood: " << std::get<2>(result) << std::endl;
    std::cout << "Spatial penalization: " << std::get<3>(result) << std::endl;

    // Minimization algorithm
    std::shared_ptr<MinimizationAlgorithm<ORDER, mydim, ndim>> minimizationAlgo =
                MinimizationAlgorithm_factory<ORDER, mydim, ndim>::createStepSolver(dataProblem, functionalProblem,
                                                                                    direction_method, step_method);

    // Finite element density estimation
    FEDE<ORDER, mydim, ndim> fede(dataProblem, functionalProblem, minimizationAlgo, preprocess_method);

    // Perform the whole task
    fede.apply();

    // Collect results
    VectorXr g_sol = fede.getDensity_g();
    //std::vector<const VectorXr*> f_init = fede.getInitialDensity();
    //Real lambda_sol = fede.getBestLambda();
    //std::vector<Real> CV_errors = fede.getCvError();
    //const std::vector<Point<ndim> >& data = dataProblem.data();

    // Copy results in a .txt file
    writeSolution<VectorXr>(g_sol, "../data/space/solution/solution_s.txt");

    return 0;
}
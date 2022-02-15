//
// Created by simonepanzeri on 14/12/2021.
//

#ifndef DEV_FDAPDE_DATA_GENERATION_H
#define DEV_FDAPDE_DATA_GENERATION_H

#include "FdaPDE.h"
#include "Spline.h"

#include <fstream>

// Generate mesh data (nodes, sides, triangles, neighbors, mesh_time) from .txt files
template<typename T>
std::vector<T> readMesh(const std::string filename) {
    std::ifstream ist(filename);
    std::vector<T> result;
    if(!ist)
        std::cerr << "Problem with file " << filename << " opening";
    else {
        UInt dimension;
        ist >> dimension;
        result.reserve(dimension);
        T n;
        while (ist >> n)
            result.push_back(n);
        ist.close();
    }

    return result;
}

// Generate spatial data
template<UInt ndim>
std::vector<Point<ndim>> createLocations() {
    std::vector<Point<ndim>> locations;
    std::ifstream ist("../data/space/data_s.txt");
    if(!ist)
        std::cerr << "Problem with file data_s.txt opening";
    else {
        UInt dimension;
        ist >> dimension;
        locations.reserve(dimension/ndim);
        for (UInt n = 0; n < dimension/ndim; ++n) {
            Real x, y, z;
            ist >> x;
            ist >> y;
            ist >> z;

            locations.emplace_back(Point<ndim>(n, {x, y, z}));
        }
        ist.close();
    }

    return locations;
}

// Generate spatio-temporal data
template<UInt ndim>
std::pair<std::vector<Point<ndim>>, std::vector<Real>> createSTLocations() {
    std::vector<Point<ndim>> locations;
    std::vector<Real> time;
    std::ifstream ist("../data/space_time/data_st.txt");
    if(!ist)
        std::cerr << "Problem with file data_st.txt opening";
    else {
        UInt dimension;
        ist >> dimension;
        locations.reserve(dimension/(ndim+1));
        time.reserve(dimension/(ndim+1));
        for (UInt n = 0; n < dimension/(ndim+1); ++n) {
            Real x, y, t;
            ist >> x;
            ist >> y;
            ist >> t;

            locations.emplace_back(Point<ndim>(n, {x, y}));
            time.push_back(t);
        }
        ist.close();
    }

    return std::make_pair(locations, time);
}

/*
// Generate the solution (evaluating at mesh_time time instants)
VectorXr createSTSolution (const UInt ns, const UInt M, const std::vector<Real>& mesh_time, const std::vector<Real>& t,
                           const VectorXr& g_sol) {
    const UInt nt = mesh_time.size();
    const UInt n = t.size();

    MatrixXr phi(n,M);
    Spline<3,2> spline(mesh_time);

    for (UInt i = 0; i < n; ++i) {
        //std::cout << "first block" << std::endl;
        for (UInt j = 0; j < M; ++j) {
            phi(i,j) = spline.BasisFunction(j, t[i]);
        }
    }

    VectorXr result;
    result.resize(ns*n);
    for(UInt j = 0; j < n; ++j) {
        //std::cout << "second block" << std::endl;
        for(UInt k = 0; k < ns; ++k) {
            result[k+j*ns] = g_sol[k*nt]*phi(j,0);
        }
    }

    for(UInt i = 1; i < M; ++i) {
        //std::cout << "third block" << std::endl;
        for(UInt j = 0; j < n; ++j) {
            if(phi(j,i)!=0) {
                for(UInt k = 0; k < ns; ++k) {
                    result[k+j*ns] = result[k+j*ns] + g_sol[k*nt+i]*phi(j,i);
                }
            }
        }
    }

    return result;
}
*/

// Generate the solution (evaluating at t time instants)
VectorXr createSTSolution (const UInt ns, const UInt M, const std::vector<Real>& mesh_time, const std::vector<Real>& t,
                           const VectorXr& g_sol) {
    const UInt nt = mesh_time.size();
    const UInt n = t.size();

    MatrixXr phi(M, n);
    Spline<3,2> spline(mesh_time);

    for (UInt i = 0; i < n; ++i)
    {
        for (UInt j = 0; j < M; ++j)
        {
            phi(j,i) = spline.BasisFunction(j, t[i]);
        }
    }

    VectorXr result;
    result.resize(ns*n);
    for(UInt j = 0; j < n; ++j) {
        for(UInt k = 0; k < ns; ++k) {
            result[k+j*ns] = g_sol[k]*phi(0,j);
        }
    }
    for(UInt i = 1; i < M; i++) {
        for(UInt j = 0; j < n; ++j) {
            if(phi(i,j) != 0) {
                for(UInt k = 0; k < ns; ++k) {
                    result[k+j*ns] = result[k+j*ns] + g_sol[k+ns*i]*phi(i,j);
                }
            }
        }
    }

    return result;
}

// Export the solution to .txt file
template<typename T>
void writeSolution (const T& sol, const std::string filename) {
    std::ofstream ost(filename);
    for (const Real& elem : sol)
        ost << elem << " ";
    ost.close();
}

#endif //DEV_FDAPDE_DATA_GENERATION_H

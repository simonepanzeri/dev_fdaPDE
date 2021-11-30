#include "FdaPDE.h"

#include "RObjects.h"
#include "Point.h"
#include "DE_Data.h"
#include "Data_Problem.h"
#include "Functional_Problem.h"

std::vector<Real> createNodes();
std::vector<UInt> createSides();
std::vector<UInt> createElements();
std::vector<UInt> createNeighbors();

int main() {
    const UInt N = 2;
    const UInt Npoints = 10;
    std::vector<Point<N>> locations;

    locations.reserve(Npoints);
    locations.emplace_back(Point<N>(0, {.1, .1}));
    locations.emplace_back(Point<N>(1, {.3, .3}));
    locations.emplace_back(Point<N>(2, {.1, .1}));
    locations.emplace_back(Point<N>(3, {.5, .1}));
    locations.emplace_back(Point<N>(4, {.5, .1}));
    locations.emplace_back(Point<N>(5, {.1, .1}));
    locations.emplace_back(Point<N>(6, {.1, .1}));
    locations.emplace_back(Point<N>(7, {.3, .3}));
    locations.emplace_back(Point<N>(8, {.3, .3}));
    locations.emplace_back(Point<N>(9, {.1, .1}));

    std::vector<Real> data_time{0.75, 0.25, 0.35, 0.15, 0.75, 0.65, 0.25, 0.35, 0.75, 0.95};

    std::vector<Real> lambda{0.001, 0.01, 0.1};
    std::vector<Real> stepProposals{0.1, 1., 10., 100.};

    VectorXr fvec{{1.0, 2.0}};

    DEData<N> dedata(locations, 2, fvec, 0.023, 100, lambda, 7, 1000, stepProposals,
                     0.00001, 0.0001, false, 1);

    std::vector<Real> time_mesh{.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.};
    DEData_time<N> dedatatime(locations, data_time);

    dedatatime.printMap(std::cout);

    std::vector<Real> points_v = createNodes();
    RNumericMatrix points_m(points_v.data(), 37, 2);

    std::vector<UInt> sides_v = createSides();
    for (size_t i = 0; i < sides_v.size(); ++i)
        --sides_v[i];
    RIntegerMatrix sides_m(sides_v.data(), 37, 2);

    std::vector<UInt> elements_v = createElements();
    for (size_t i = 0; i < elements_v.size(); ++i)
        --elements_v[i];
    RIntegerMatrix elements_m(elements_v.data(), 36, 3);

    std::vector<UInt> neighbors_v = createNeighbors();
    for (size_t i = 0; i < neighbors_v.size(); ++i) {
        if (neighbors_v[i] != -1)
            --neighbors_v[i];
    }
    RIntegerMatrix neighbors_m(neighbors_v.data(),36,3);

    DataProblem<1, 2, 2> dataProb(locations, 2, fvec, 0.023, 100, lambda, 7, 1000, stepProposals,
    0.00001, 0.0001, false, 1, points_m, sides_m, elements_m, neighbors_m);

    DataProblem_time<1, 2, 2> dataProb_time(locations, data_time, 2, fvec, 0.023, 100, lambda, 7, 1000, stepProposals,
                                  0.00001, 0.0001, false, 1, points_m, sides_m, elements_m, neighbors_m, time_mesh);



    return 0;
}

std::vector<Real> createNodes() {
    return std::vector<Real> { // for each node (column): first row = x; second row = y
            0.0000000, 0.1111111, 0.2222222, 0.3333333, 0.4444444, 0.5555556, 0.6666667, 0.7777778, 0.8888889, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 1.0000000, 0.0000000, 0.1111111, 0.2222222, 0.3333333, 0.4444444, 0.5555556, 0.6666667, 0.7777778, 0.8888889, 1.0000000, 0.5000000,
            0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1111111, 0.1111111, 0.2222222, 0.2222222, 0.3333333, 0.3333333, 0.4444444, 0.4444444, 0.5555556, 0.5555556, 0.6666667, 0.6666667, 0.7777778, 0.7777778, 0.8888889, 0.8888889, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.5000000
    };
}

std::vector<UInt> createSides() {
    return std::vector<UInt> { // for each side (column): first row = first point; second row = second point
        11,  1,  2, 15, 13,  3,  2, 13,  4,  3,  4,  5, 37, 37, 17, 37,  2, 19, 37, 21, 25, 23, 29, 21, 30, 23, 28, 27, 25, 30, 37, 29, 19, 37, 31, 16, 37,  7,  8,  9, 14,  7,  8, 10, 12,  9, 14, 18, 16, 12, 18, 20, 32, 37, 33, 35, 34, 24, 33, 22, 34, 37, 24, 26, 26, 36, 22, 20, 32,  6, 37,  5,
         1,  2, 11, 13,  3, 15, 13, 11, 15,  4,  5, 37,  4, 15, 15, 17,  3, 37, 21, 19, 23, 29, 25, 30, 23, 21, 27, 25, 28, 29, 30, 28, 17, 31, 30, 37,  7, 16,  9, 14,  8,  8, 16, 12,  9, 10, 16, 37, 18, 14, 20, 37, 37, 33, 32, 34, 24, 35, 22, 34, 33, 22, 26, 35, 36, 35, 24, 22, 31,  7,  6,  6
    };
}

std::vector<UInt> createElements() {
    return std::vector<UInt> { // for each element (column): first row = first node; second row = second node; third row = third node
        11, 15, 11,  4,  4, 15, 17,  3, 19, 25, 21, 28, 30, 30, 28, 17, 37, 16,  8,  7, 10, 14, 18, 12, 18, 32, 35, 33, 22, 24, 26, 22, 22, 37,  6,  5,
         1, 13,  2, 15,  5,  4, 15, 13, 37, 23, 30, 27, 29, 21, 25, 37, 31, 37,  9,  8, 12, 16, 37, 14, 20, 37, 34, 22, 33, 26, 36, 24, 37, 32,  7,  6,
         2,  3, 13,  3, 37, 37, 37,  2, 21, 29, 23, 25, 23, 37, 29, 19, 30,  7, 14, 16,  9,  8, 16,  9, 37, 33, 24, 34, 37, 35, 35, 34, 20, 31, 37, 37
    };
}

std::vector<UInt> createNeighbors() {
    return std::vector<UInt> { // for each node (column): neighbors in the rows
        -1,  8,  8,  2, 36,  5,  6,  3, 14, 13, 13, -1, 10,  9, 10,  9, -1, 35, 24, 22, 24, 20, 18, 19, 33, 29, 32, 32, 26, 31, -1, 27, 25, -1, 18, 35,
         3,  4, -1, -1,  6,  7, 16, -1, -1, 15, -1, 15, 11, 17, -1, -1, 14, 20, 22, 18, -1, 19, -1, 21, 23, -1, 30, -1, 33, 27, 30, 28, -1, 17, 36,  5,
        -1, -1,  1,  6, -1,  4, -1,  2, 16, -1, 14, -1, -1, 11, 12,  7, 34, 23, -1, -1, -1, -1, 25, -1, -1, 34, -1, 29, 28, -1, -1, -1, 29, 26, -1, -1
    };
}
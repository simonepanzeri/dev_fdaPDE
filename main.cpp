#include "FdaPDE.h"

#include "RObjects.h"
#include "Point.h"
#include "DE_Data.h"
#include "Data_Problem.h"

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
    RNumericMatrix points_m(points_v.data(), 2, 37);

    std::vector<UInt> sides_v = createSides();
    RIntegerMatrix sides_m(sides_v.data(), 2, 72);

    std::vector<UInt> elements_v = createElements();
    RIntegerMatrix elements_m(elements_v.data(), 3, 36);

    std::vector<UInt> neighbors_v = createNeighbors();
    RIntegerMatrix neighbors_m(neighbors_v.data(),3,36);

    /*DataProblem<1, 2, 2> dataProb(locations, 2, fvec, 0.023, 100, lambda, 7, 1000, stepProposals,
    0.00001, 0.0001, false, 1, points_m, sides_m, elements_m, neighbors_m);*/

    return 0;
}

std::vector<Real> createNodes() {
    return std::vector<Real> {
            0.0000000, 0.0000000,
            0.1111111, 0.0000000,
            0.2222222, 0.0000000,
            0.3333333, 0.0000000,
            0.4444444, 0.0000000,
            0.5555556, 0.0000000,
            0.6666667, 0.0000000,
            0.7777778, 0.0000000,
            0.8888889, 0.0000000,
            1.0000000, 0.0000000,
            0.0000000, 0.1111111,
            1.0000000, 0.1111111,
            0.0000000, 0.2222222,
            1.0000000, 0.2222222,
            0.0000000, 0.3333333,
            1.0000000, 0.3333333,
            0.0000000, 0.4444444,
            1.0000000, 0.4444444,
            0.0000000, 0.5555556,
            1.0000000, 0.5555556,
            0.0000000, 0.6666667,
            1.0000000, 0.6666667,
            0.0000000, 0.7777778,
            1.0000000, 0.7777778,
            0.0000000, 0.8888889,
            1.0000000, 0.8888889,
            0.0000000, 1.0000000,
            0.1111111, 1.0000000,
            0.2222222, 1.0000000,
            0.3333333, 1.0000000,
            0.4444444, 1.0000000,
            0.5555556, 1.0000000,
            0.6666667, 1.0000000,
            0.7777778, 1.0000000,
            0.8888889, 1.0000000,
            1.0000000, 1.0000000,
            0.5000000, 0.5000000
    };
}

std::vector<UInt> createSides() {
    return std::vector<UInt> {
            11,  1,
             1,  2,
             2, 11,
            15, 13,
            13,  3,
             3, 15,
             2, 13,
            13, 11,
             4, 15,
             3,  4,
             4,  5,
             5, 37,
            37,  4,
            37, 15,
            17, 15,
            37, 17,
             2,  3,
            19, 37,
            37, 21,
            21, 19,
            25, 23,
            23, 29,
            29, 25,
            21, 30,
            30, 23,
            23, 21,
            28, 27,
            27, 25,
            25, 28,
            30, 29,
            37, 30,
            29, 28,
            19, 17,
            37, 31,
            31, 30,
            16, 37,
            37,  7,
             7, 16,
             8,  9,
             9, 14,
            14,  8,
             7,  8,
             8, 16,
            10, 12,
            12,  9,
             9, 10,
            14, 16,
            18, 37,
            16, 18,
            12, 14,
            18, 20,
            20, 37,
            32, 37,
            37, 33,
            33, 32,
            35, 34,
            34, 24,
            24, 35,
            33, 22,
            22, 34,
            34, 33,
            37, 22,
            24, 26,
            26, 35,
            26, 36,
            36, 35,
            22, 24,
            20, 22,
            32, 31,
             6,  7,
            37,  6,
             5,  6
    };
}

std::vector<UInt> createElements() {
    return std::vector<UInt> {
            11,  1,  2,
            15, 13,  3,
            11,  2, 13,
             4, 15,  3,
             4,  5, 37,
            15,  4, 37,
            17, 15, 37,
             3, 13,  2,
            19, 37, 21,
            25, 23, 29,
            21, 30, 23,
            28, 27, 25,
            30, 29, 23,
            30, 21, 37,
            28, 25, 29,
            17, 37, 19,
            37, 31, 30,
            16, 37,  7,
             8,  9, 14
    };
}

std::vector<UInt> createNeighbors() {
    return std::vector<UInt> {
        -1,    3,   -1,
         8,    4,   -1,
         8,   -1,    1,
         2,   -1,    6,
        36,    6,   -1,
         5,    7,    4,
         6,   16,   -1,
         3,   -1,    2,
        14,   -1,   16,
        13,   15,   -1,
        13,   -1,   14,
        -1,   15,   -1,
        10,   11,   -1,
         9,   17,   11,
        10,   -1,   12,
         9,   -1,    7,
        -1,   14,   34,
        35,   20,   23,
        24,   22,   -1,
        22,   18,   -1,
        24,   -1,   -1,
        20,   19,   -1,
        18,   -1,   25,
        19,   21,   -1,
        33,   23,   -1,
        29,   -1,   34,
        32,   30,   -1,
        32,   -1,   29,
        26,   33,   28,
        31,   27,   -1,
        -1,   30,   -1,
        27,   28,   -1,
        25,   -1,   29,
        -1,   17,   26,
        18,   36,   -1,
        35,    5,   -1
    };
}
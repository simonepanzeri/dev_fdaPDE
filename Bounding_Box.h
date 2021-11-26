//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_BOUNDING_BOX_H
#define DEV_FDAPDE_BOUNDING_BOX_H

#include "FdaPDE.h"
#include "Mesh_Objects.h"

//! Policy used when boxes have to be stored in the tree (NDIMP is the physical dimension of Box: 2 -> 2D, 3 -> 3D)
template<int NDIMP>
class Box {
protected:
    //! A vector of rectangle corner coordinates.
    // First NDIMP values are the coordinates of the rectangle corner with minimum coordinates,
	// followed by the coordinates of the opposite corner. (2D: xmin, ymin, xmax, ymax)
	std::vector<Real> x_;

public:
    //! Default constructor
    Box();
    //! Another constructor
    Box(std::vector<Real> const& coord);

    template<UInt NNODES, int NDIME, int NDIMPP>
    Box(Element<NNODES, NDIME, NDIMPP> const& element);

    //! A method returning the i-th coordinate value.
    inline Real operator[](int const& i) { return x_[i]; }
    //! A const method returning the i-th coordinate value.
    inline Real operator[](int const& i) const { return x_[i]; }
    //! A method returning the number of physical space dimension.
    inline static constexpr int dp() { return NDIMP; }
    //! A method returning the number of dimensions used for the search (2*NDIMP).
    inline static constexpr int dt() { return 2*NDIMP; }
    //! A method returning the size of coordinate array (composed of min and max).
    inline static constexpr int coordsize() { return 2*NDIMP; }
    //! A method setting coordinate values.
    void set(std::vector<Real> const& data);
    //! A method to get coordinate values.
    std::vector<Real> get() const {return x_; };
    //! A method to print minimum box point and maximum box point.
    void print(std::ostream& out) const;
};

#include "Bounding_Box_imp.h"

#endif //DEV_FDAPDE_BOUNDING_BOX_H

//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_PROJECTION_H
#define DEV_FDAPDE_PROJECTION_H

#include "FdaPDE.h"
#include "Mesh.h"

template <UInt ORDER,UInt mydim, UInt ndim>
class projection{
public:
    projection(const MeshHandler<ORDER,mydim,ndim>&, const std::vector<Point<ndim>>&);
    std::vector<Point<ndim>> computeProjection();
};

template<UInt ORDER>
class projection<ORDER, 2, 3>{
private:
    const MeshHandler<ORDER, 2, 3>& mesh_;
    const std::vector<Point<3>>& deData_; // the points to be projected
    const UInt num_points;

    std::vector<UInt> computeNodePatch(UInt) const;

public:
    projection(const MeshHandler<ORDER, 2, 3>& m, const std::vector<Point<3>>& d) : mesh_(m), deData_(d),
        num_points(d.size()) {};

    std::vector<Point<3>> computeProjection();
};

#include "Projection_imp.h"

#endif //DEV_FDAPDE_PROJECTION_H

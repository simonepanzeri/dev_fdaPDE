//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_AD_TREE_H
#define DEV_FDAPDE_AD_TREE_H

#include "FdaPDE.h"
#include "Mesh_Objects.h"
#include "Domain.h"
#include "Tree_Node.h"
#include "Tree_Header.h"
#include "Bounding_Box.h"
#include "Exception_Handling.h"

/**	\class ADTree
 * 	\brief Alternating binary range searching tree.
 *	\param Shape: template parameter, the original shape
 */
template<class Shape>
class ADTree {
protected:
    /** The header.
     *
     *  It contains general information about the tree.
     */
    TreeHeader<Shape> header_;
    /// Vector of tree nodes.
    std::vector<TreeNode<Shape>> data_;
    /** \brief Adds a point to the tree.
     * 	It throws:
     * 	<ul>
     * 	<li> a TreeDomainError exception if the point is out of domain;
     * 	<li> a TreeAlloc exception if there is no more space in the tree to add the node;
     * 	<li> a LevRuntimeError if you exceed the limit set for the tree levels due to the inclusion of the node.
     * 	</ul>
     */
    int adtrb(Id shapeid, std::vector<Real> const& coords);
    //! A method to handle a TreeDomainError exception.
    int handledomerr(Id shapeid, std::vector<Real> const& coords);
    //! A method to handle a TreeAlloc exception.
    int handletreealloc(Id shapeid, std::vector<Real> const& coords);
    //! A method to handle a LevRuntimeError exception.
    int handleleverr(Id shapeid, std::vector<Real> const& coords);
    /** Searches dimension associated to a given level.
     *
     * 	\param[in] lev The given level.
     * 	\param[in] dim The number of dimensions used for the search.
     */
    inline int searchdim(int const& lev, int const& dim) const {
        return (lev % dim);
    }
    /** Finds delta associated to division at a given level.
     *
     * 	\param[in] lev The given level.
     * 	\param[in] dim The number of dimensions used for the search.
     */
    inline double delta(int const& lev, int const& dim) const {
        return std::pow(0.5, int(lev/dim)+1);
    }

    /** It fills all the locations of the tree. Object's coordinates are stored to perform searching operations.
   *  See mesh_handler to verify what points and triangle must contain.
   */
    //void setTree(SEXP Rmesh);

    void setTree(const RNumericMatrix& points, const RIntegerMatrix& triangle);

public:

    ADTree(const RNumericMatrix& points, const RIntegerMatrix& elements);

    /// Returns a reference to the tree header.
    inline TreeHeader<Shape> gettreeheader() const { return header_; }
    /** A Method to add a node to the tree.
     * 	It calls the handlers of the exceptions that can be thrown by adtrb().
     *
     * 	\param[in] coords Coordinates of the point.
     *
     *	The location of the current node in the tree is returned.
     */
    int addtreenode(Id shapeid, std::vector<Real> const& coords);
    /** A method to get out the information stored at a given node.
     *
     * 	\param[in] loc Location of the searched node.
     * 	\param[out] coord Bounding box coordinates of the object stored in the loc-th location.
     *  \param[out] id Id of searched node.
     */
    inline void gettri(int const& loc, std::vector<Real>& coord, Id& id);
    /** A method to get out the node stored at a given location.
     *
     * 	\param[in] loc Location of the searched node.
     */
    inline TreeNode<Shape> gettreenode(int const& loc) const{return data_[loc];};

    /** A method to find all points or (bounding) boxes that intersect a given box.
     *
     * 	\param[in] region Box where searching described by the representative point obtained through a corner transformation.
     * 	\param[out] found Indices of the found elements.
     *
     * 	This function returns true if it has completed successfully, false otherwise.
     */
    bool search(std::vector<Real> const& region, std::set<int>& found)const;
    /// Deletes a specified location in the tree.
    //void deltreenode(int const& index);
    //! A method to get the j-th coordinate of the bounding box of the p-th object stored in the node.
    inline Real pointcoord(int const& p, int const& j) const { return data_[p].getcoord(j); }
    //! A method to get the the Id of the original object of the p-th treenode.
    inline Id pointId (int const& p) const { return data_[p].getid(); }
    //! A method to output information contained in the tree header.
    template<class S>
    friend std::ostream& operator<<(std::ostream& ostr, ADTree<S> const& myadt);
};

#include "AD_Tree_imp.h"

#endif //DEV_FDAPDE_AD_TREE_H

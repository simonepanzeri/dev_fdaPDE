//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_TREE_HEADER_H
#define DEV_FDAPDE_TREE_HEADER_H

#include "FdaPDE.h"
#include "Mesh_Objects.h"
#include "Bounding_Box.h"
#include "Domain.h"
#include "Tree_Node.h"

//! It contains general information about the tree.
template<class T>
class TreeHeader {
protected:
    // Tree memory locations.
    int tree_loc_;

    // Tree levels.
    int tree_lev_;

    // Number of physical space dimensions (typically 2 or 3).
    int ndimp_;

    // Number of dimensions used for the search (2*ndimp because we use box).
    int ndimt_;

    //	Number of logical locations currently used in the tree. Initialized to 0.
    int nele_;

    //! Tree indices: the use of iava and iend avoids the necessity of initializing the stack of available locations.
    int iava_;
    int iend_;

    // Tree's domain.
    Domain<T> tree_domain_;

    //! A protected constructor.
    TreeHeader(int const& ntree, Domain<T> const& d) :
        tree_loc_(ntree), tree_lev_(0), ndimp_(T::dp()), ndimt_(T::dt()), nele_(0), iava_(1), iend_(1), tree_domain_(d) {}
    //! A method trying to set the number of tree memory locations (throws a LocLengthError exception if nt is out of range).
    void stml(int const & nt) { tree_loc_ = nt; }
public:
    //! Default constructor.
    TreeHeader() : tree_loc_(0), tree_lev_(0), ndimp_(T::dp()), ndimt_(T::dt()), nele_(0), iava_(1), iend_(1),
                   tree_domain_(){}

    // constructor in case there is already tree information
    TreeHeader(int const& tree_loc, int const& tree_lev, int const& ndimp, int const& ndimt,
               int const& nele, int const& iava, int const& iend, Domain<T> const& tree_domain):
        tree_loc_(tree_loc), tree_lev_(tree_lev), ndimp_(ndimp), ndimt_(ndimt),
        nele_(nele), iava_(iava), iend_(iend), tree_domain_(tree_domain) {}

    //! A method to get the number of tree memory locations.
    inline int gettreeloc() const { return tree_loc_; }
    //! A method to set the number of tree memory locations (handles a LocLengthError exception).
    void settreeloc(int const& nt) { stml(nt); }
    //! A method to get the number of tree levels.
    inline int gettreelev() const { return tree_lev_; }
    //! A method to set the number of tree levels.
    inline void settreelev(int const& nl) { tree_lev_ = nl; }
    //! A method to get the number of physical space dimension.
    inline int getndimp() const { return ndimp_; }
    //! A method to get the number of dimensions used for the search.
    inline int getndimt() const { return ndimt_; }
    //! A method to get the number of logical locations currently used in the tree.
    inline int getnele() const { return nele_; }
    //! A method to set the number of logical locations currently used in the tree.
    inline void setnele(int const & ne) { nele_ = ne; }
    //! A method to get the next available location in the stack.
    inline int getiava() const { return iava_; }
    //! A method to set the next available location in the stack.
    inline void setiava(int const & ia) { iava_ = ia; }
    //! A method to get the next available location in the tree free store.
    inline int getiend() const { return iend_; }
    //! A method to set the next available location in the tree free store.
    inline void setiend(int const & ie) { iend_ = ie; }
    //! A method to get the i-th coordinate of the origin of the tree's (domain) bounding box. --> queste falle friend
    inline double domainorig(int const & i) const { return tree_domain_.orig(i); }
    //! A method to get the i-th scaling factor of the tree's bounding box.
    inline double domainscal(int const & i) const { return tree_domain_.scal(i); }
    //! Output operator.
    /* It prints out:
     *	- number of tree memory locations;
     *	- number of tree levels;
     *	- number of physical space dimensions;
     *	- number of pieces of information carried by the tree;
     *	- number of dimensions used for the search;
     *	- number of logical locations currently used in the tree;
     *	- tree domain.
     */
    template<class S>
    friend std::ostream& operator<<(std::ostream& ostr, TreeHeader<S> const& head) {
        ostr << "General informations about the tree" << std::endl;
        ostr << "----------------------------------" << std::endl;
        ostr << "Tree memory locations: " << head.tree_loc_ << std::endl;
        ostr << "Number of tree levels: " << head.tree_lev_ << std::endl;
        ostr << "Number of physical space dimension: " << head.ndimp_ << std::endl;
        ostr << "Number of dimensions used for the search: " << head.ndimt_ << std::endl;
        ostr << "Number of logical locations currently used in the tree: " << head.nele_ << std::endl;
        ostr << head.tree_domain_;
        ostr << "----------------------------------" << std::endl;

        return ostr;
    }

    template<class S>
    friend TreeHeader<S> createtreeheader(int const& nt, Domain<S> const &d) {TreeHeader<T> hd(nt, d); return hd;}
};

#endif //DEV_FDAPDE_TREE_HEADER_H

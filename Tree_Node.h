//
// Created by simonepanzeri on 25/11/2021.
//

#ifndef DEV_FDAPDE_TREE_NODE_H
#define DEV_FDAPDE_TREE_NODE_H

#include "FdaPDE.h"
#include "Bounding_Box.h"
#include "Mesh_Objects.h"

//! Define a tree node. The template parameter T is the original shape; in the treenode are stored the bbox of the
//! original object and an index that identify that object.

template<class T>
class TreeNode {
protected:
    // Position of the father node. (It's used by the algorithm for deleting a tree node.)
    //`int father_;

    //! Bounding Box of the object to be stored
    Box<T::dp()> box_;

    //! Positions of left and right children.
    int children_[2];

    //! The id of the Element that create the Box.
    //Be careful! This is not the id of the Treenode but id of Element
    Id id_;

    // do per scontato che ci sia un oggetto mesh che contenga tutte le forme (salvate in qualche modo) e che le
    // identifichi attraverso un id di tipo Uint!
    // generalizzando si potrebbe usare un puntatore a Shape (parametro template della forma generica, può essere un
    // triangolo, o l'id stesso se si vuole tornare al caso precedente), in teoria poi non devo distruggere la memoria
    // perchè la forma è salvata in una struttra mesh che deve rimanere inalterata
    // Shape * id_;
    // un'altra possibilità è salvare la forma anche nella struttura mesh con un shared_ptr e mettere anche qui uno
    // shared_ptr<Shape>

public:
    //!	Default constructor.
    TreeNode(): box_() { //father_(0),
        children_[0] = 0;
        children_[1] = 0;
        id_ = std::numeric_limits<UInt>::max();
    }
    //!	Another constructor (T is the shape of the id; it is needed for Box constructor; it works with Element or Box).
    TreeNode(Id const id, T shape): box_(shape) { //father_(0),
        children_[0] = 0;
        children_[1] = 0;
        id_ = id;
    }

    TreeNode(Id const id, const Box<T::dp()>& shape): box_(shape) { //father_(0),
        children_[0] = 0;
        children_[1] = 0;
        id_ = id;
    }

    // constructor in case there is already tree information
    TreeNode(Box<T::dp()> const& box, Id const& id, int const& left_child, int const& right_child):
            box_(box), id_(id) {
        children_[0] = left_child;
        children_[1] = right_child;
        //father_(0)
    }
    //! A method setting the father.
    //inline void setfather(int const & ifth) { father_ = ifth; }

    //! A method setting a child.
    inline void setchild(short int const& flag, int const& child) { children_[flag] = child; }
    //! A method returning the father.
    //inline int getfather() const { return father_; }
    //! A method returning a child.
    inline int getchild(short int const& flag) const { return children_[flag]; }
    //! A method setting the coordinates of the bbox stored in the node.
    inline void setcoords(std::vector<Real> const& data) { box_.set(data); }
    //! A method to get the i-th coordinate of the bounding box of the object stored in the node.
    inline Real getcoord(int const& i) const { return box_[i]; }
    //! A method to set the id stored in the node.
    inline void setid(Id id) { id_ = id; }
    //! A method to get the id stored in the node.
    inline Id getid() const { return id_; }
    //! A method returning a reference to box_.
    inline Box<T::dp()>& getbox() { return box_; }
    //! A method printing information about the treenode.
    void print(std::ostream& out) const {
        //Be careful! this id is not the id of the Treenode but id of Triangle
        out << "------------------------------" << std::endl;
        out << "Shape Id: --" << id_ << "--" <<std::endl;
        //out << "Father:  " << father_ << std::endl;
        out << "Left Children:  " << children_[0] << std::endl;
        out << "Right Children:  " << children_[1] << std::endl;
        out << "Box: ";
        box_.print(out);
        out << std::endl;
    }
};

#endif //DEV_FDAPDE_TREE_NODE_H

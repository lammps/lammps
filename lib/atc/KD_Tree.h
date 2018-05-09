#ifndef KD_TREE_H
#define KD_TREE_H

#include "Array2D.h"
#include "MatrixDef.h"
#include "MatrixLibrary.h"
#include <cmath>
#include <vector>
#include <utility>

class Node {
  public:
    Node() {
      // Zero argument constructor just for default initialization.
    }

    Node(int node, double x, double y, double z) 
      : index_(node)
    {
      coords_[0] = x;
      coords_[1] = y;
      coords_[2] = z;
    }

    int index_;
    double coords_[3];

    double distanceSquared(Node query) {
      return pow(coords_[0] - query.coords_[0], 2) 
              + pow(coords_[1] - query.coords_[1], 2) 
              + pow(coords_[2] - query.coords_[2], 2);
    }

    double distanceInDimension(Node query, int dimension) {
      return pow(coords_[dimension] - query.coords_[dimension], 2);
    }

    bool lessThanInDimension(Node query, int dimension) {
      if (dimension == 0) return Node::compareX(*this, query);
      else if (dimension == 1) return Node::compareY(*this, query);
      else if (dimension == 2) return Node::compareZ(*this, query);
      else return false;
    }

    bool operator==(const Node &rhs) {
      return ((*this).coords_[0] == rhs.coords_[0] &&
              (*this).coords_[1] == rhs.coords_[1] &&
              (*this).coords_[2] == rhs.coords_[2]);
    }

    static bool compareX(Node one, Node two) { return one.coords_[0] < two.coords_[0]; }
    static bool compareY(Node one, Node two) { return one.coords_[1] < two.coords_[1]; }
    static bool compareZ(Node one, Node two) { return one.coords_[2] < two.coords_[2]; }

};

typedef std::pair<int,std::vector<Node> > Elem;

class KD_Tree {
  public:
    static KD_Tree* create_KD_tree(const int nNodesPerElement, const int nNodes, 
                                   const DENS_MAT *points, const int nElems, 
                                   const Array2D<int> &conn);

    KD_Tree(std::vector<Node> *points, std::vector<Elem> *elems,
            int dimension=0);

    ~KD_Tree();

    std::vector<int> find_nearest(DENS_VEC query) { 
      // Create a fake Node and find the nearest Node in the tree to it.
      return find_nearest_elements(Node(-1, query(0), query(1), query(2)));
    }
    
    std::vector<int> find_nearest_elements(Node query, int dimension=0);
    
    std::vector<std::vector<int> > getElemIDs(int depth);

  private:
    Node value_;

    std::vector<Node> *sortedPts_;
    std::vector<Elem> *candElems_;
    KD_Tree *leftChild_;
    KD_Tree *rightChild_;

};

#endif

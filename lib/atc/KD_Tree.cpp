#include "KD_Tree.h"
#include <assert.h>

KD_Tree *KD_Tree::create_KD_tree(const int nNodesPerElem, const int nNodes, 
                                 const DENS_MAT *nodalCoords, const int nElems, 
                                 const Array2D<int> &conn) {
  vector<Node> *points = new vector<Node>(); // Initialize an empty list of Nodes
  for (int node = 0; node < nNodes; node++) { // Insert all nodes into list
    points->push_back(Node(node, (*nodalCoords)(0, node),
                                 (*nodalCoords)(1, node),
                                 (*nodalCoords)(2, node)));
  }
  vector<Elem> *elements = new vector<Elem>();
  for (int elem = 0; elem < nElems; elem++) {
    vector<Node> nodes = vector<Node>();
    for (int node = 0; node < nNodesPerElem; node++) {
      nodes.push_back((*points)[conn(node, elem)]);
    }
    elements->push_back(Elem(elem, nodes));
  }
  return new KD_Tree(points, elements);
}

KD_Tree::~KD_Tree() {
  delete sortedPts_;
  delete leftChild_;
  delete rightChild_;
}

KD_Tree::KD_Tree(vector<Node> *points, vector<Elem> *elements, 
                 int dimension)
  : candElems_(elements) {
  // Set up comparison functions
  bool (*compare)(Node, Node);
  if (dimension == 0) {
    compare = Node::compareX;
  } else if (dimension == 1) {
    compare = Node::compareY;
  } else {
    compare = Node::compareZ;
  }

  // Sort points by their coordinate in the current dimension
  sort(points->begin(), points->end(), compare); 
  sortedPts_ = points;

  // Pick the median point as the root of the tree
  size_t nNodes = points->size();
  size_t med = nNodes/2;
  value_ = (*sortedPts_)[med];
  
  // Recursively construct the left sub-tree
  vector<Node> *leftPts   = new vector<Node>;
  vector<Elem> *leftElems = new vector<Elem>;
  // Recursively construct the right sub-tree
  vector<Node> *rightPts   = new vector<Node>;
  vector<Elem> *rightElems = new vector<Elem>;
  for (vector<Elem>::iterator elit = candElems_->begin();
       elit != candElems_->end(); elit++) {
    // Identify elements that should be kept on either side
    bool foundElemLeft  = false;
    bool foundElemRight = false;
    for (vector<Node>::iterator ndit = elit->second.begin();
         ndit != elit->second.end(); ndit++) {
      // Search this node
      if (compare(*ndit, value_)) {
        if (find(leftPts->begin(), leftPts->end(), *ndit) == leftPts->end()) {
          leftPts->push_back(*ndit);
        }
        foundElemLeft = true;
      }
      if (compare(value_, *ndit)) {
        if (find(rightPts->begin(), rightPts->end(), *ndit) == rightPts->end()) {
          rightPts->push_back(*ndit);
        }
        foundElemRight = true;
      }
    }
    if (foundElemLeft) leftElems->push_back(*elit);
    if (foundElemRight) rightElems->push_back(*elit);
  }

  // Create child tree, or NULL if there's nothing to create
  if (candElems_->size() - leftElems->size() < 4 || leftElems->size() == 0) {
    leftChild_ = NULL;
  } else {
    leftChild_ = new KD_Tree(leftPts, leftElems, (dimension+1) % 3);
  }
  // Create child tree, or NULL if there's nothing to create
  if (candElems_->size() - rightElems->size() < 4 || rightElems->size() == 0) {
    rightChild_ = NULL;
  } else {
    rightChild_ = new KD_Tree(rightPts, rightElems, (dimension+1) % 3);
  }
}

vector<int> KD_Tree::find_nearest_elements(Node query, int dimension) {
  // if the root coordinate is less than the query coordinate
  
 
  // If the query point is less that the value (split) point of this
  // tree, either recurse to the left or return this node's elements
  // if there is no left child.
  if (query.lessThanInDimension(value_, dimension)) {
    if (leftChild_ == NULL) {
      vector<int> result = vector<int>();
      for (vector<Elem>::iterator elem = candElems_->begin();
           elem != candElems_->end(); elem++) {
        result.push_back(elem->first);
      }
      return result;
    }
    return leftChild_->find_nearest_elements(query, (dimension+1) % 3);
  } else {
    if (rightChild_ == NULL) {
      vector<int> result = vector<int>();
      for (vector<Elem>::iterator elem = candElems_->begin();
           elem != candElems_->end(); elem++) {
        result.push_back(elem->first);
      }
      return result;
    }
    return rightChild_->find_nearest_elements(query, (dimension+1) % 3);
  }
}

vector<vector<int> > KD_Tree::getElemIDs(int depth) {
  
  vector<vector<int> > result;
  vector<vector<int> > temp;
  
  assert(depth >= 0 );
  if (depth == 0) {
    vector<int> candElemIDs;
    vector<Elem>::iterator it;
    for(it = candElems_->begin(); it != candElems_->end(); ++it) {
      candElemIDs.push_back((*it).first);
    }

    sort(candElemIDs.begin(), candElemIDs.end());
    result.push_back(candElemIDs);

  } else if (leftChild_ == NULL || rightChild_ == NULL) {
    // Insert all nodes at this level once,
    // then insert a bunch of empty vectors.
    temp = this->getElemIDs(0);
    result.insert(result.end(), temp.begin(), temp.end());
    int numRequested = floor(pow(2,depth));
    for (int i = 0; i < numRequested - 1; ++i) {
      vector<int> emptyVec;
      result.push_back(emptyVec);
    }
  } else {
    --depth;
    temp = leftChild_->getElemIDs(depth);
    result.insert(result.end(), temp.begin(), temp.end());
    temp = rightChild_->getElemIDs(depth);
    result.insert(result.end(), temp.begin(), temp.end());
  }
  
  return result;
}

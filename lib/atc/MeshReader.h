#ifndef MESH_READER_H
#define MESH_READER_H

#include "Array2D.h"
#include "MatrixDef.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"
#include "FE_Mesh.h"

#include <sstream>

using namespace std;

namespace ATC {

  class MeshReader {

  public:
    /** constructor, takes a filename */
    MeshReader(string filename, Array<bool> periodicity, double tol=1.e-8);

    /** destructor */
    ~MeshReader();

    /** creates handle to new mesh object */
    FE_Mesh* create_mesh();

  private:
    /** reads .mesh format file */
    void read_mesh_file();

    /** reads .exo format file */
    void read_exo_file();

    /** helper function for parsing mesh type string */
    int int_from_str(string str); 

    /** Data members for storing necessary information */
    string meshfile_;
    Array<bool> periodicity_;
    string elementType_;
    int nNodes_;
    int nElements_;
    int nNodeSets_;
    Array2D<int> * conn_;
    DENS_MAT * nodeCoords_;
    Array<pair<string,set<int> > > * nodeSets_;
    double coordTol_;
  };

}; // end namespace ATC

#endif

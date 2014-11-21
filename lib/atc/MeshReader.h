#ifndef MESH_READER_H
#define MESH_READER_H

#include "Array2D.h"
#include "MatrixDef.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"
#include "FE_Mesh.h"
#include <sstream>
#include <set>
#include <utility>
#include <string>

namespace ATC {

  class MeshReader {

  public:
    /** constructor, takes a filename */
    MeshReader(std::string filename, Array<bool> periodicity, double tol=1.e-8);

    /** destructor */
    ~MeshReader();

    /** creates handle to new mesh object */
    FE_Mesh* create_mesh();

  private:
    int number_of_vertices(std::string str);

    /** reads .mesh format file */
    void read_mesh_file();

    /** reads .exo format file */
    void read_exo_file();

    /** Data members for storing necessary information */
    std::string meshfile_;
    ATC_matrix::Array<bool> periodicity_;
    std::string elementType_;
    int nNodes_;
    int nElements_;
    int nNodeSets_;
    ATC_matrix::Array2D<int> * conn_;
    DENS_MAT * nodeCoords_;
    ATC_matrix::Array<std::pair<std::string,std::set<int> > > * nodeSets_;
  };

}; // end namespace ATC

#endif

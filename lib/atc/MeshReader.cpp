#include "MeshReader.h"
#include "LammpsInterface.h"
#include "Utility.h"
#ifdef HAS_EXODUS
//#include <stdio.h> 
//#include "netcdf.h" 
#include "exodusII.h"
#endif

using ATC_Utility::to_string;
using std::ifstream;
using std::istringstream;
using std::stringstream;
using std::map;
using std::pair;
using std::set;
using std::string;

namespace ATC {
  /** constructor, takes a filename */
  MeshReader::MeshReader(string filename, 
                         Array<bool> periodicity, 
                         double tol)
    : meshfile_(filename),
      periodicity_(periodicity),
      nNodes_(0),
      nElements_(0)
  {
    conn_ = new Array2D<int>();
    nodeCoords_ = new DENS_MAT;
    nodeSets_ = new Array< std::pair< string,set<int> > >();

    size_t idx = filename.rfind('.');
    if (idx == string::npos) {
      throw ATC_Error("Given mesh file is of unknown type.");
    }

    string ext = filename.substr(idx+1);
    
    if      (ext == "mesh"){ read_mesh_file(); }
    else { throw ATC_Error("mesh file is of unknown type."); }
  }

  /** destructor */
  MeshReader::~MeshReader()
  {
    if (conn_) delete conn_;
    if (nodeCoords_) delete nodeCoords_;
    if (nodeSets_) delete nodeSets_;
  }

  /** creates handle to new mesh object */
  FE_Mesh* MeshReader::create_mesh()
  {
    return new FE_3DMesh(elementType_,
                         nNodes_, nElements_,
                         conn_, nodeCoords_,
                         periodicity_,
                         nodeSets_);
  }

  int MeshReader::number_of_vertices(string str)
  {
    string temp;
    int number=0;
    for (unsigned int i=0; i < str.size(); i++) {
      if (isdigit(str[i])) {
        for (unsigned int a=i; a<str.size(); a++) {
          temp += str[a];               
        }
        break;
      }
    }
    istringstream(temp) >> number;
    return number;
  }


  /** reads .mesh format file */
  void MeshReader::read_mesh_file() {
    ifstream in;
    in.open(meshfile_.c_str(), ifstream::in);
    string header;
    while (getline(in,header)) {
      istringstream words(header);
      string section;
      words >> section;
      if (section == "Coordinates") {
        words >> nNodes_;
        nodeCoords_->reset(3, nNodes_, false);
        string line;
        for (int i = 0; i < nNodes_; ++i) {
          getline(in,line);
          istringstream coords(line);
          coords >> (*nodeCoords_)(0, i);
          coords >> (*nodeCoords_)(1, i);
          coords >> (*nodeCoords_)(2, i);
        }
        stringstream ss;
        ss << "read " << nNodes_ << " nodes";
        ATC::LammpsInterface::instance()->print_msg_once(ss.str());
      }
      else if (section == "Elements") {
        words >> nElements_;
        words >> elementType_;
        int nVerts = number_of_vertices(elementType_);
        conn_->reset(nVerts, nElements_);
        string line;
        for (int i = 0; i < nElements_; ++i) {
          getline(in,line);
          istringstream verts(line);
          for (int j = 0; j < nVerts; ++j ) {
            int node;
            verts >> node;
            (*conn_)(j,i) = node-1;
          }
        }
        stringstream ss;
        ss << "read " << nElements_ << " " << elementType_ << " elements";
        ATC::LammpsInterface::instance()->print_msg_once(ss.str());
      }
      else if (section == "Nodesets") {
        words >> nNodeSets_;
        nodeSets_->reset(nNodeSets_);
        string infoline;
        string setName;
        int nNodes;
        string line;
        for (int i = 0; i < nNodeSets_; ++i) {
          getline(in,infoline);
          istringstream info(infoline);
          info >> setName;
          info >> nNodes;
          (*nodeSets_)(i).first = setName;
          getline(in,line);
          istringstream nodes(line);
          for (int j = 0; j < nNodes; ++j ) {
            int node;
            nodes >> node;
            (*nodeSets_)(i).second.insert(node-1);
          }
        }
      }
    }
    in.close();
    if (nodeCoords_->size() == 0) {
      throw ATC_Error("Could not find mesh file, or mesh file empty.");
    }
  }

  /** reads .exo format file */
  void MeshReader::read_exo_file() {
#ifndef HAS_EXODUS
    throw ATC_Error("Reading ExodusII .exo files not supported."); 
#else
    int CPU_word_size=0,IO_word_size=0;
    float version;
    int exoid = ex_open (meshfile_.c_str(), EX_READ, 
      &CPU_word_size, &IO_word_size, &version);
    if (exoid < 0) { throw ATC_Error("couldn't open "+meshfile_); }
    int nsd,nelemblk,nfsets;
    char title[MAX_LINE_LENGTH+1];
    int error = ex_get_init (exoid, title, 
      &nsd, &nNodes_, &nElements_, &nelemblk, &nNodeSets_, &nfsets);
    if (error > 0) { throw ATC_Error("problem with init "+meshfile_+" "+title); }
    // coordinates 
    float x[nNodes_], y[nNodes_], z[nNodes_];
    error = ex_get_coord (exoid, x, y, z);
    if (error > 0) { throw ATC_Error("problem with getting coordinates "+meshfile_); }
    nodeCoords_->reset(nsd,nNodes_);
    DENS_MAT & nodes = *nodeCoords_;
    for (int i = 0; i < nNodes_; ++i) {
      nodes(0,i) = x[i];  // this is a float to double conversion 
      nodes(1,i) = y[i];  // this is a float to double conversion 
      nodes(2,i) = z[i];  // this is a float to double conversion 
    }
    ATC::LammpsInterface::instance()->print_msg_once("read "+to_string(nNodes_)+
" nodes");
    // elements
    int blkIds[nelemblk],nblkelem[nelemblk],nnpe[nelemblk],na[nelemblk];
    error = ex_get_elem_blk_ids(exoid, blkIds);
    char etype[MAX_STR_LENGTH+1];
    string lastType;
    for (int i=0; i<nelemblk; i++){ 
      error = ex_get_elem_block (exoid, blkIds[i], etype,
        &(nblkelem[i]), &(nnpe[i]), &(na[i]));
      elementType_ = etype;
      if (i > 0 && elementType_ != lastType )
        { throw ATC_Error(meshfile_+" is composed of multiple types"); }
      lastType = etype;
    }
    int nVerts = number_of_vertices(elementType_);
    conn_->reset(nVerts, nElements_);
    int n = 0;
    for (int i=0; i<nelemblk; i++) { 
      int bconn[nnpe[i]*nblkelem[i]];
      error = ex_get_elem_conn (exoid, blkIds[i], &bconn); 
      for (int j=0; j<nblkelem[i]; j++) { 
        for (int k=0; k<nnpe[i]; k++) { 
          (*conn_)(k,n) = bconn[k+j*nnpe[i]]-1;
        }
        n++;
      }
      ATC::LammpsInterface::instance()->print_msg_once("read "+to_string(n)+" "+elementType_+" elements, block "+to_string(i+1)+"/"+to_string(nelemblk));
    }
    // nodesets
    int nsetIds[nNodeSets_];
    error = ex_get_node_set_ids (exoid, nsetIds);
    int nnodes,ndist;
    //nodeSets_ = new Array< pair< string,set<int> > >();
    nodeSets_->reset(nNodeSets_);
    for (int i=0; i<nNodeSets_; i++) { 
      (*nodeSets_)(i).first = to_string(nsetIds[i]);
      error = ex_get_node_set_param (exoid, nsetIds[i], &nnodes, &ndist);
      int nodes[nnodes];
      error = ex_get_node_set (exoid, nsetIds[i], nodes);
      for (int j=0; j<nnodes; j++) { 
        (*nodeSets_)(i).second.insert(nodes[j]-1);
      }
    }
    error = ex_close(exoid);
    if (error > 0) { throw ATC_Error("problem with closing "+meshfile_); }
#endif
  }
}; // end namespace ATC

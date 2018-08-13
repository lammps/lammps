// ATC header files
#include "FE_Element.h"
#include "FE_Mesh.h"
#include "LammpsInterface.h"
#include "ATC_Error.h"
#include "OutputManager.h"
#include "Utility.h"
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <functional>

// Other headers
#include <iostream>


using namespace std;
using ATC_Utility::dbl_geq;
using ATC_Utility::parse_min;
using ATC_Utility::parse_max;
using ATC_Utility::parse_minmax;
using ATC_Utility::split_values;
using ATC_Utility::plane_coords;
using ATC_Utility::norm3;
using ATC_Utility::to_string;
using ATC_Utility::tolerance;





namespace ATC {

  // constants
  const static double tangentTolerance = 0.01;
  const static double tol = 1.0e-10;


  // =============================================================
  //   class FE_Mesh
  // =============================================================
  FE_Mesh::FE_Mesh()
    : decomposition_(false), 
    lammpsPartition_(false),
    partitioned_(false),
    nNodes_(0), 
    nNodesUnique_(0),
    feElement_(NULL),
    twoDimensional_(false),
    hasPlanarFaces_(false)

  {
  }

  // -------------------------------------------------------------
  FE_Mesh::~FE_Mesh()
  {
    if (feElement_) delete feElement_;
  }

  // -------------------------------------------------------------
  //  modify
  // -------------------------------------------------------------
  bool FE_Mesh::modify(int narg, char **arg)
  {
    bool match = false;

    if (strcmp(arg[0],"mesh")==0) 
    {
     /*! \page man_mesh_create_faceset_box fix_modify AtC mesh create_faceset box
        \section syntax
         fix_modify AtC mesh create_faceset <id> box
         <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> <in|out> [units]
         - <id> = id to assign to the collection of FE faces
         - <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of
         the bounding box that is coincident with the desired FE faces
         - <in|out> = "in" gives inner faces to the box, 
                      "out" gives the outer faces to the box
         - units = option to specify real as opposed to lattice units
        \section examples
         <TT> fix_modify AtC mesh create_faceset obndy box -4.0 4.0 -12 12 -12 12 out </TT>
        \section description
          Command to assign an id to a set of FE faces.
        \section restrictions
        Only viable for rectangular grids.
        \section related
        \section default
        The default options are units = lattice and the use of outer faces
       */
     /*! \page man_mesh_create_faceset_plane fix_modify AtC mesh create_faceset plane 
        \section syntax
         fix_modify AtC mesh create_faceset <id> plane 
         <x|y|z> <val1> <x|y|z> <lval2> <uval2> [units]
         - <id> = id to assign to the collection of FE faces
         - <x|y|z> = coordinate directions that define plane on which faceset lies 
         - <val1>,<lval2>,<uval2> = plane is specified as the x|y|z=val1 plane bounded by
              the segments x|y|z = [lval2,uval2]  
         - units = option to specify real as opposed to lattice units
        \section examples
         <TT> fix_modify AtC mesh create_faceset xyplane plane y 0 x -4 0 </TT>
        \section description
          Command to assign an id to a set of FE faces. 
        \section restrictions
        Only viable for rectangular grids. 
        \section related
        \section default
        The default option is units = lattice.
       */
      if (strcmp(arg[1],"create_faceset")==0) 
      {
        int argIdx = 2;
        string tag = arg[argIdx++];
        if (strcmp(arg[argIdx],"plane")==0) 
        {
          argIdx++;
          int ndir, idir[3], isgn;
          double xlimits[3][2];
          parse_plane(argIdx, narg, arg, ndir, idir, isgn, xlimits);
          if (xlimits[idir[1]][0] == xlimits[idir[1]][1]) 
            split_values(xlimits[idir[1]][0],xlimits[idir[1]][1]);
          if (xlimits[idir[2]][0] == xlimits[idir[2]][1]) 
            split_values(xlimits[idir[2]][0],xlimits[idir[2]][1]);
          parse_units(argIdx,narg,arg,
            xlimits[0][0],xlimits[0][1],
            xlimits[1][0],xlimits[1][1],
            xlimits[2][0],xlimits[2][1]);
          if (ndir > 1) {
            create_faceset(tag, xlimits[idir[0]][0], idir[0], isgn, 
              idir[1], xlimits[idir[1]][0], xlimits[idir[1]][1], 
              idir[2], xlimits[idir[2]][0], xlimits[idir[2]][1]);
          }
          else {
            create_faceset(tag, xlimits[idir[0]][0], idir[0], isgn);
          }
          match = true;
        }
        // bounding_box
        else 
        {
          if (strcmp(arg[argIdx],"box")==0) argIdx++;
          double xmin = parse_min(arg[argIdx++]);
          double xmax = parse_max(arg[argIdx++]);
          double ymin = parse_min(arg[argIdx++]);
          double ymax = parse_max(arg[argIdx++]);
          double zmin = parse_min(arg[argIdx++]);
          double zmax = parse_max(arg[argIdx++]);
          bool outward = true;
          if (narg > argIdx && (strcmp(arg[argIdx++],"in") == 0)) 
            outward = false;
          parse_units(argIdx,narg,arg, xmin,xmax,ymin,ymax,zmin,zmax);
          create_faceset(tag, xmin, xmax, ymin, ymax, zmin, zmax, outward);
          match = true;
        }
      }
   /*! \page man_mesh_create_nodeset fix_modify AtC mesh create_nodeset
      \section syntax
      fix_modify AtC mesh create_nodeset <id>
      <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>
      - <id> = id to assign to the collection of FE nodes
      - <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of
      the bounding box that contains only the desired nodes
      \section examples
      <TT> fix_modify AtC mesh create_nodeset lbc -12.1 -11.9 -12 12 -12 12 </TT>
      \section description
      Command to assign an id to a set of FE nodes to be used subsequently
      in defining boundary conditions.
      \section restrictions
      \section related
      \section default
      Coordinates are assumed to be in lattice units.
    */
      else if (strcmp(arg[1],"create_nodeset")==0) {
        int argIdx = 2;
        string tag  = arg[argIdx++];
        double xmin, xmax, ymin, ymax, zmin, zmax;
        if (strcmp(arg[argIdx],"plane")==0) {
          argIdx++;
          int ndir, idir[3], isgn;
          double xlimits[3][2];
          parse_plane(argIdx, narg, arg, ndir, idir, isgn, xlimits);
          xmin = xlimits[0][0];
          xmax = xlimits[0][1];
          ymin = xlimits[1][0];
          ymax = xlimits[1][1];
          zmin = xlimits[2][0];
          zmax = xlimits[2][1];
        }
        else {
          xmin = parse_min(arg[argIdx++]);
          xmax = parse_max(arg[argIdx++]);
          ymin = parse_min(arg[argIdx++]);
          ymax = parse_max(arg[argIdx++]);
          zmin = parse_min(arg[argIdx++]);
          zmax = parse_max(arg[argIdx++]);
        }
        if (xmin == xmax) split_values(xmin,xmax);
        if (ymin == ymax) split_values(ymin,ymax);
        if (zmin == zmax) split_values(zmin,zmax);
        parse_units(argIdx,narg,arg, xmin,xmax,ymin,ymax,zmin,zmax);
        create_nodeset(tag, xmin, xmax, ymin, ymax, zmin, zmax);
        match = true;
      }
   /*! \page man_mesh_add_to_nodeset fix_modify AtC mesh add_to_nodeset
      \section syntax
      fix_modify AtC mesh add_to_nodeset <id>
      <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>
      - <id> = id of FE nodeset to be added to 
      - <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of
      the bounding box that contains the desired nodes to be added
      \section examples
      <TT> fix_modify AtC mesh add_to_nodeset lbc -11.9 -11 -12 12 -12 12 </TT>
      \section description
      Command to add nodes to an already existing FE nodeset.
      \section restrictions
      \section related
      \section default
      Coordinates are assumed to be in lattice units.
    */
      else if (strcmp(arg[1],"add_to_nodeset")==0) {
        string tag  = arg[2];
        double xmin = parse_min(arg[3]);
        double xmax = parse_max(arg[4]);
        double ymin = parse_min(arg[5]);
        double ymax = parse_max(arg[6]);
        double zmin = parse_min(arg[7]);
        double zmax = parse_max(arg[8]);
        add_to_nodeset(tag, xmin, xmax, ymin, ymax, zmin, zmax);
        match = true;
      }
   /*! \page man_mesh_create_elementset fix_modify AtC mesh create_elementset
      \section syntax
      fix_modify AtC create_elementset <id>
      <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>
      - <id> = id to assign to the collection of FE element
      - <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of
      the bounding box that contains only the desired elements
      \section examples
      <TT> fix_modify AtC mesh create_elementset middle -4.1 4.1 -100 100 -100 1100 </TT>
      \section description
      Command to assign an id to a set of FE elements to be used subsequently
      in defining material and mesh-based operations.
      \section restrictions
      Only viable for rectangular grids.
      \section related
      \section default
      Coordinates are assumed to be in lattice units.
    */
      else if (strcmp(arg[1],"create_elementset")==0) {
        int argIdx = 2;
        string tag  = arg[argIdx++];
        double xmin = 0;
        double xmax = 0;
        double ymin = 0;
        double ymax = 0;
        double zmin = 0;
        double zmax = 0;
        if (narg > 4) {
          xmin = parse_min(arg[argIdx++]);
          xmax = parse_max(arg[argIdx++]);
          ymin = parse_min(arg[argIdx++]);
          ymax = parse_max(arg[argIdx++]);
          zmin = parse_min(arg[argIdx++]);
          zmax = parse_max(arg[argIdx++]);
          if (xmin == xmax) split_values(xmin,xmax);
          if (ymin == ymax) split_values(ymin,ymax);
          if (zmin == zmax) split_values(zmin,zmax);
          parse_units(argIdx,narg,arg, xmin,xmax,ymin,ymax,zmin,zmax);
        }
        else {
          string regionName = arg[argIdx++];
          LammpsInterface::instance()->
          region_bounds(regionName.c_str(),xmin,xmax,ymin,ymax,zmin,zmax);
        }
        create_elementset(tag, xmin, xmax, ymin, ymax, zmin, zmax);
        match = true;
      }
/*! \page man_mesh_nodeset_to_elementset fix_modify AtC mesh nodeset_to_elementset
      \section syntax
      fix_modify AtC nodeset_to_elementset <nodeset_id> <elementset_id> <max/min>
      - <nodeset_id> = id of desired nodeset from which to create elementset
      - <elementset_id> = id to assign to the collection of FE element
      - <max/min> = flag to choose either the maximal or minimal elementset
      \section examples
      <TT> fix_modify AtC mesh nodeset_to_elementset myNodeset myElementset min </TT>
      \section description
      Command to create an elementset from an existing nodeset.  Either the minimal element set
      of elements with all nodes in the set, or maximal element set with all elements with at
      least one node in the set, can be created
      \section restrictions
      None.
      \section related
      \section default
      Unless specified, the maximal element set is created
    */
      else if (strcmp(arg[1],"nodeset_to_elementset")==0) {
        string nodesetId = arg[2];
        string elementsetId = arg[3];
        set<int> elementSet;
        if (narg < 5) {
          this->nodeset_to_maximal_elementset(nodesetId,elementSet);
        }
        else {
          if (strcmp(arg[4],"max")==0) {
            this->nodeset_to_maximal_elementset(nodesetId,elementSet);
          }
          else if (strcmp(arg[4],"min")==0) {
            this->nodeset_to_minimal_elementset(nodesetId,elementSet);
          }
          else {
            return false;
          }
        }
        elementSetMap_[elementsetId] = elementSet;
        match = true;
      }
   /*! \page man_mesh_output fix_modify AtC mesh output
      \section syntax
      fix_modify AtC mesh output <file_prefix> 
      \section examples
      <TT> fix_modify AtC mesh output meshData </TT> \n
      \section description
      Command to output mesh and associated data: nodesets, facesets, and
      elementsets. This data is only output once upon initialization since
      currently the mesh is static. Creates (binary, "gold" format) Ensight 
      output of mesh data.
      \section restrictions
      none
      \section related
      \section default
      none
    */
      else if (strcmp(arg[1],"output")==0) {
        string outputPrefix  = arg[2];
        output(outputPrefix);
        match = true;
      }
    }
    return match;
  }
  // -------------------------------------------------------------
  void FE_Mesh::parse_units(int & argIdx, int narg, char ** arg,
    double & xmin, double & xmax, double & ymin, double & ymax, 
    double & zmin, double & zmax)
  {
    if (narg > argIdx && (strcmp(arg[argIdx++],"units") == 0)) {}
    else 
    { // scale from lattice units to physical units
      xmin *= xscale_; xmax *= xscale_;
      ymin *= yscale_; ymax *= yscale_;
      zmin *= zscale_; zmax *= zscale_;
    }
  }
  // -------------------------------------------------------------
  //   parse plane
  // -------------------------------------------------------------
  void FE_Mesh::parse_plane(int & argIdx, int narg, char ** arg,  
    int & ndir, int * idir, int & isgn, double xlimits[][2])
  {
    ndir = 0;
    xlimits[0][0] = parse_min("-INF");
    xlimits[0][1] = parse_max("INF");
    xlimits[1][0] = parse_min("-INF");
    xlimits[1][1] = parse_max("INF");
    xlimits[2][0] = parse_min("-INF");
    xlimits[2][1] = parse_max("INF");
    string n(arg[argIdx++]);
    int i1dir = -1, i2dir = -1, i3dir = -1;
    string_to_index(n, i1dir, isgn);
    idir[ndir++] = i1dir;
    idir[ndir]   = (i1dir+1) % 3;
    idir[ndir+1] = (i1dir+2) % 3;
    xlimits[i1dir][0] =  parse_minmax(arg[argIdx++]);
    xlimits[i1dir][1] = xlimits[i1dir][0];
    if (narg > argIdx ) {
      if (string_to_index(arg[argIdx],i2dir)) {
        argIdx++;
        xlimits[i2dir][0] = parse_min(arg[argIdx++]); 
        xlimits[i2dir][1] = parse_max(arg[argIdx++]);
        idir[ndir++] = i2dir;
      }
      if (narg > argIdx ) {
        if (string_to_index(arg[argIdx],i3dir)) {
          argIdx++;
          xlimits[i3dir][0] = parse_min(arg[argIdx++]); 
          xlimits[i3dir][1] = parse_max(arg[argIdx++]);
        }
      }
      else {
        i3dir = 0;
        if      ((i1dir == 0 && i2dir == 1) || (i1dir == 1 && i2dir == 0) ) {
          i3dir = 2;
        }
        else if ((i1dir == 0 && i2dir == 2) || (i1dir == 2 && i2dir == 0) ) {
          i3dir = 1;
        }
      }
      idir[ndir++] = i3dir;
    }
    if ((idir[0]==idir[1]) || (idir[0]==idir[2]) || (idir[1]==idir[2]) ) {
      throw ATC_Error( "inconsistent directions in plane:"+to_string(idir[0]+1)+" "+to_string(idir[1]+1)+" "+to_string(idir[2]+1));  
    }
  }
  // -------------------------------------------------------------
  //   initialize
  // -------------------------------------------------------------
  void FE_Mesh::initialize(void) 
  {

    bool aligned = is_aligned();
    if (!aligned) {
      feElement_->set_projection_guess(CENTROID_LINEARIZED);
      ATC::LammpsInterface::instance()->print_msg_once("WARNING: mesh is not aligned with the coordinate directions atom-to-element mapping will be expensive");
      // if HEX8 -> orient();
    }
    bool twoD = is_two_dimensional(); 
    if (twoD) {
      feElement_->set_projection_guess(TWOD_ANALYTIC);
      if (feElement_->order()< 3) hasPlanarFaces_ = true; 
      ATC::LammpsInterface::instance()->print_msg_once(" mesh is two dimensional");
    }
  }
  //-----------------------------------------------------------------
  
  //-----------------------------------------------------------------
  void FE_Mesh::write_mesh(string meshFile)
  {
    ofstream out;
    out.open(meshFile.c_str());
    DENS_MAT & x = *(coordinates());
    Array2D<int> & conn = *(connectivity());
    int nNodes = x.nCols(); // transpose
    int ndm = x.nRows();
    int nElems = conn.nCols();
    int nodesPerElem = conn.nRows(); // transpose
    out << "Coordinates " << nNodes << "\n";
    for (int n = 0; n < nNodes; n++) {
      for (int i = 0; i < ndm; i++) {
        out << "  " << std::setprecision(16) << x(i,n);
      }
      out << "\n";
    }
    out << "\n";
    string type = element_type();
    out << "Elements " << nElems << " " << type << "\n";
    for (int n = 0; n < nElems; n++) {
      for (int i = 0; i < nodesPerElem; i++) {
        out << 1+conn(i,n) << " ";
      }
      out << "\n";
    }
    out << "\n";
    if (nodeSetMap_.size()) {
      out << "Nodesets " << nodeSetMap_.size() <<"\n";
      NODE_SET_MAP::const_iterator niter;
      map<string,DENS_MAT> nodesets;
      for (niter = nodeSetMap_.begin(); niter != nodeSetMap_.end(); niter++) {
        string name = niter->first;
        const set<int> & nset = niter->second;
        out << name << " " << nset.size() << "\n";
        set<int>::const_iterator iter;
        for (iter = nset.begin(); iter != nset.end(); iter++) {
          out << *iter << "  " ;   
        }
        out << "\n";
      }
    }
  }
  // -------------------------------------------------------------
  //   test whether almost structured
  // -------------------------------------------------------------
  bool FE_Mesh::is_aligned(void) const
  {
    vector<bool> foundBestMatch(nSD_,false);
    vector<DENS_VEC> tangents(nSD_);
    DENS_VEC xi0(nSD_);
    xi0 = 0;
    DENS_MAT eltCoords;
    for (int ielem = 0; ielem < nElts_; ielem++) {
       element_coordinates(ielem,eltCoords);
       feElement_->tangents(eltCoords,xi0,tangents,true); 

      for (unsigned i = 0; i < tangents.size(); i++) {
        // find maximum value for which global axis its closest to
        int maxIndex = 0;
        double maxValue = abs(tangents[i](0));
        for (int j = 1; j < nSD_; j++) {
          if (abs(tangents[i](j)) > maxValue) {
            maxValue = abs(tangents[i](j));
            maxIndex = j;
          }
        }

        // make sure no other tangent is associated with this direction
        if (foundBestMatch[maxIndex]) {
          return false;
        }
        else {
          foundBestMatch[maxIndex] = true;
        }

        // compute deviation from a perfectly aligned vector
        double error = 0.;
        for (int j = 1; j < nSD_; j++) {
          if (j != maxIndex) {
            error += abs(tangents[i](j));
          }
        }
        error /= maxValue;
        if (error > tangentTolerance) {
          return false;
        }
      }
    }
    return true;
  }

  // -------------------------------------------------------------
  //   element_type
  // -------------------------------------------------------------
  string FE_Mesh::element_type(void) const  {
      int npe = feElement_->num_elt_nodes(); 
      if      (npe == 4)  { return "TET4"; }
      else if (npe == 8)  { return "HEX8"; }
      else if (npe == 20) { return "HEX20"; }
      else if (npe == 27) { return "HEX27"; }
      return "UNKNOWN";
  }

  // -------------------------------------------------------------
  //   element_coordinates
  // -------------------------------------------------------------
  void FE_Mesh::element_coordinates(const int eltID,
                                    DENS_MAT & xCoords) const
  {
    const int nne = num_nodes_per_element();
    xCoords.reset(nSD_, nne, false);
    for (int inode=0; inode<nne; inode++) {
      const int id = element_connectivity_global(eltID, inode);
      for (int isd=0; isd<nSD_; isd++) {
        xCoords(isd,inode) = nodalCoords_(isd,id);
      }
    }
  
  }
  // -------------------------------------------------------------
  //   position
  // -------------------------------------------------------------
  void FE_Mesh::position(const int eltID,
                         const VECTOR & xi,
                               DENS_VEC & x) const
  {
    const int nne = num_nodes_per_element();
    DENS_VEC N;
    feElement_->shape_function(xi,N);
    x.reset(nSD_); 
    for (int inode=0; inode<nne; inode++) {
      const int id = element_connectivity_global(eltID, inode);
      for (int isd=0; isd<nSD_; isd++) {
        x(isd) += nodalCoords_(isd,id)*N(inode);
      }
    }
  }

  // -------------------------------------------------------------
  // element size in each direction
  // -------------------------------------------------------------
  void FE_Mesh::bounding_box(const int ielem,
                              DENS_VEC & xmin, DENS_VEC & xmax) 
  {
    xmin.reset(nSD_);
    xmax.reset(nSD_);
    int nne = num_nodes_per_element();
    for (int isd=0; isd<nSD_; isd++) {
      int id = element_connectivity_global(ielem, 0);
      double x = nodalCoords_(isd,id);
      xmin(isd) = x;
      xmax(isd) = x;
      for (int inode=1; inode<nne; inode++) {
        id = element_connectivity_global(ielem, inode);
        x = nodalCoords_(isd,id);
        xmin(isd) = min(xmin(isd), x );
        xmax(isd) = max(xmax(isd), x );
      }
    }
  }

  // -------------------------------------------------------------
  // element size in each direction 
  // -------------------------------------------------------------
  void FE_Mesh::element_size(const int ielem,
                              double & hx, double & hy, double & hz) 
  {
    DENS_VEC xmin(nSD_), xmax(nSD_);
    bounding_box(ielem,xmin,xmax);
    hx = xmax(0)-xmin(0);
    hy = xmax(1)-xmin(1);
    hz = xmax(2)-xmin(2);
  }

  // -------------------------------------------------------------
  //   face_coordinates
  // -------------------------------------------------------------
  void FE_Mesh::face_coordinates(const PAIR face, DENS_MAT & xCoords) const
  {
    const int eltID=face.first, faceID=face.second;
    const int nnf = num_nodes_per_face();
    const Array2D <int> & local_conn = local_face_connectivity();

    xCoords.reset(nSD_, nnf, false);

    for (int inode=0; inode < nnf; inode++) 
    {
      int id = element_connectivity_global(eltID, local_conn(faceID,inode));
      for (int isd=0; isd<nSD_; isd++) 
        xCoords(isd,inode) = nodalCoords_(isd,id);
    }
  }

  // -------------------------------------------------------------
  //   nodal_coordinates
  // -------------------------------------------------------------
  DENS_VEC FE_Mesh::nodal_coordinates(const int nodeID) const
  {
    DENS_VEC xCoords(nSD_, false);
    const int id = uniqueToGlobalMap_(nodeID);
    for (int isd=0; isd<nSD_; isd++)
      xCoords(isd) = nodalCoords_(isd, id);
    return xCoords;
  }
  DENS_VEC FE_Mesh::global_coordinates(const int nodeID) const
  {
    DENS_VEC xCoords(nSD_, false);
    for (int isd=0; isd<nSD_; isd++)
      xCoords(isd) = nodalCoords_(isd, nodeID);
    return xCoords;
  }

  // -------------------------------------------------------------
  //   query_nodeset
  // -------------------------------------------------------------
  bool FE_Mesh::query_nodeset(const string & name) const
  {
    if (name == "all")  return true;
    if (nodeSetMap_.find(name) == nodeSetMap_.end()) return false;
    return true;
  }

  // -------------------------------------------------------------
  //   get_nodeset
  // -------------------------------------------------------------
  const set<int> & FE_Mesh::nodeset(const string & name) const
  {
    NODE_SET_MAP::const_iterator iter = nodeSetMap_.find(name);
    if (name == "all") return nodeSetAll_;
    else if (iter == nodeSetMap_.end()) 
      throw ATC_Error( "No nodeset with name " + name + " found.");  
    else return iter->second;
  }

  // -------------------------------------------------------------
  //   get_elementset
  // -------------------------------------------------------------
  const set<int> & FE_Mesh::elementset(const string & name) const
  {
    NODE_SET_MAP::const_iterator iter = elementSetMap_.find(name);
    if (name == "all") return elementSetAll_;
    else if (iter == elementSetMap_.end()) 
      throw ATC_Error( "No elementset with name " + name + " found.");  
    else return iter->second;
  }

  // -------------------------------------------------------------
  //   nodeset_to_minimal_elementset
  // -------------------------------------------------------------
  void FE_Mesh::nodeset_to_minimal_elementset
    (const string & name, set<int> & elemSet) const
  {
    if (name == "all") {
      for (int ielem = 0; ielem < nElts_; ielem++) {
        elemSet.insert(ielem);
      }
    }
    else {
      NODE_SET_MAP::const_iterator iter = nodeSetMap_.find(name);
      if (iter == nodeSetMap_.end()) 
        throw ATC_Error( "No nodeset with name " + name + " found.");
      nodeset_to_minimal_elementset(iter->second,elemSet);
      if (elemSet.size()==0) {
        throw ATC_Error("No elements found in minimal condensation of nodeset " + name);
      }
    }
  }

  // -------------------------------------------------------------
  //   nodeset_to_minimal_elementset
  // -------------------------------------------------------------
  void FE_Mesh::nodeset_to_minimal_elementset
  (const set<int> & nodeSet, set<int> & elemSet) const
  {
    int npe = num_nodes_per_element();
    for (int ielem=0; ielem < nElts_; ielem++) {
      int inode = 0;
      bool in = true;
      while (in && inode < npe) {
        int node = element_connectivity_unique(ielem, inode);
        set<int>::const_iterator iter = nodeSet.find(node);
        if (iter == nodeSet.end()) { in=false; }
        inode++;
      }
      if (in) elemSet.insert(ielem);
    }
  }

  // -------------------------------------------------------------
  //   nodeset_to_maximal_elementset
  // -------------------------------------------------------------
  void FE_Mesh::nodeset_to_maximal_elementset(const string &name, set<int> &elemSet) const
  {
    if (name == "all") {
      for (int ielem = 0; ielem < nElts_; ielem++) {
        elemSet.insert(ielem);
      }
    }
    else {
      NODE_SET_MAP::const_iterator iter = nodeSetMap_.find(name);
      if (iter == nodeSetMap_.end()) 
        throw ATC_Error( "No nodeset with name " + name + " found.");
      nodeset_to_maximal_elementset(iter->second,elemSet);
      if (elemSet.size()==0) {
        throw ATC_Error("No elements found in maximal condensation of nodeset " + name);
      }
    }
  }

  // -------------------------------------------------------------
  //   nodeset_to_maximal_elementset
  // -------------------------------------------------------------
  void FE_Mesh::nodeset_to_maximal_elementset(const set<int> &nodeSet, set<int> &elemSet) const
  {
    int npe = num_nodes_per_element();
    for (int ielem = 0; ielem < nElts_; ielem++) {
      int inode = 0;
      bool in = false;
      while (!in && inode < npe) {
        int node = element_connectivity_unique(ielem, inode);
        set<int>::const_iterator iter = nodeSet.find(node);
        if (iter != nodeSet.end()) { in = true; }
        inode++;
      }
      if (in) elemSet.insert(ielem);
    }
  }

  // -------------------------------------------------------------
  //   elementset_to_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::elementset_to_nodeset
    (const set<int> & elemSet, set<int>  nodeSet) const
  {
    int npe = num_nodes_per_element();
    set<int>::const_iterator itr;
    for (itr = elemSet.begin(); itr != elemSet.end(); itr++) 
    {
      int ielem = *itr;
      for (int inode=0; inode < npe; inode++) 
      {
        int node = element_connectivity_global(ielem, inode);
        nodeSet.insert(node);
      }
    }
  }

  // -------------------------------------------------------------
  //   elementset_to_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::elementset_to_nodeset
    (const string & name, set<int>  nodeSet) const
  {
    if (name == "all")
      for (int ielem = 0; ielem < nElts_; ielem++) 
        nodeSet.insert(ielem);

    else 
    {
      ELEMENT_SET_MAP::const_iterator iter = elementSetMap_.find(name);
      if (iter == elementSetMap_.end()) 
        throw ATC_Error( "No elementset with name " + name + " found.");

      int npe = num_nodes_per_element();
      const set<int> &elemSet = iter->second;
      set<int>::const_iterator itr;
      for (itr = elemSet.begin(); itr != elemSet.end(); itr++) 
      {
        int ielem = *itr;
        for (int inode=0; inode < npe; inode++) 
        {
          int node = element_connectivity_unique(ielem, inode);
          nodeSet.insert(node);
          inode++;
        }
      }
    }
  }
  set<int>  FE_Mesh::elementset_to_nodeset
    (const string & name) const
  {
    set<int> nset;
    elementset_to_nodeset(name,nset);
    return nset;
  }

  // -------------------------------------------------------------
  //   elementset_to_minimal_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::elementset_to_minimal_nodeset
    (const string & name, set<int> & nodeSet) const
  {
    // return:  set - complement_of_set
    if (name == "all")  { return;}
    else 
    {
      elementset_to_nodeset(name,nodeSet);
      set<int> compElemSet;
      elementset_complement(name,compElemSet);
      int npe = num_nodes_per_element();
      set<int>::const_iterator itr;
      for (itr = compElemSet.begin(); itr != compElemSet.end(); itr++) 
      {
        int ielem = *itr;
        for (int inode=0; inode < npe; inode++) 
        {
          int node = element_connectivity_unique(ielem, inode);
          nodeSet.erase(node);
          inode++;
        }
      }
    }
  }

  // -------------------------------------------------------------
  //   elementset_complement
  // -------------------------------------------------------------
  void FE_Mesh::elementset_complement
    (const string & name, set<int> & cElemSet) const
  {
    // return:  set - complement_of_set
    if (name == "all")  { return;}
    else 
    {
      ELEMENT_SET_MAP::const_iterator iter = elementSetMap_.find(name);
      if (iter == elementSetMap_.end()) 
        throw ATC_Error( "No elementset with name " + name + " found.");

      const set<int> &elemSet = iter->second;
      for (int ielem = 0; ielem < nElts_; ielem++) 
      {
        if(elemSet.find(ielem) == elemSet.end() ) cElemSet.insert(ielem);
      }
    }
  }

  // -------------------------------------------------------------
  //   elementset_complement
  // -------------------------------------------------------------
  void FE_Mesh::elementset_complement
    (const set<int> & elemSet, set<int> & cElemSet) const
  {
    for (int ielem = 0; ielem < nElts_; ielem++) 
    {
      if(elemSet.find(ielem) == elemSet.end() ) cElemSet.insert(ielem);
    }
  }
  // -------------------------------------------------------------
  //   faceset_to_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::faceset_to_nodeset(const string &name, set<int> &nodeSet) const
  {
    if (name == "all") {
      for (int inode = 0; inode < nNodesUnique_; inode++)
        nodeSet.insert(inode);
    }
    else 
    {
      FACE_SET_MAP::const_iterator faceset = faceSetMap_.find(name);
      if (faceset == faceSetMap_.end()) 
        throw ATC_Error( "No faceset with name " + name + " found.");
      const set<PAIR> & faceSet = faceset->second;
      set<PAIR>::const_iterator iter;
      Array <int> conn;
      for (iter = faceSet.begin(); iter != faceSet.end(); iter++) 
      {
        PAIR face = *iter;
        face_connectivity_unique(face,conn);
        for (int i = 0; i < conn.size() ; ++i) {
          int inode = conn(i);
          nodeSet.insert(inode);
        }
      }
    }
  }
  void FE_Mesh::faceset_to_nodeset(const set<PAIR> &faceSet, set<int> &nodeSet) const
  {
    set<PAIR>::const_iterator iter;
    Array <int> conn;
    for (iter = faceSet.begin(); iter != faceSet.end(); iter++) 
    {
      PAIR face = *iter;
      face_connectivity_unique(face,conn);
      for (int i = 0; i < conn.size() ; ++i) {
        int inode = conn(i);
        nodeSet.insert(inode);
      }
    }
  }
  // -------------------------------------------------------------
  //   faceset_to_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::faceset_to_nodeset_global(const string &name, set<int> &nodeSet) const
  {
    if (name == "all") {
      for (int inode = 0; inode < nNodes_; inode++)
        nodeSet.insert(inode);
    }
    else 
    {
      FACE_SET_MAP::const_iterator faceset = faceSetMap_.find(name);
      if (faceset == faceSetMap_.end()) 
        throw ATC_Error( "No faceset with name " + name + " found.");
      const set<PAIR> & faceSet = faceset->second;
      set<PAIR>::const_iterator iter;
      Array <int> conn;
      for (iter = faceSet.begin(); iter != faceSet.end(); iter++) 
      {
        PAIR face = *iter;
        face_connectivity(face,conn);
        for (int i = 0; i < conn.size() ; ++i) {
          int inode = conn(i);
          nodeSet.insert(inode);
        }
      }
    }
  }
  void FE_Mesh::faceset_to_nodeset_global(const set<PAIR> &faceSet, set<int> &nodeSet) const
  {
    set<PAIR>::const_iterator iter;
    Array <int> conn;
    for (iter = faceSet.begin(); iter != faceSet.end(); iter++) 
    {
      PAIR face = *iter;
      face_connectivity(face,conn);
      for (int i = 0; i < conn.size() ; ++i) {
        int inode = conn(i);
        nodeSet.insert(inode);
      }
    }
  }

  // -------------------------------------------------------------
  //   get_faceset
  // -------------------------------------------------------------
  const set<PAIR> &FE_Mesh::faceset(const string & name) const
  {
    FACE_SET_MAP::const_iterator iter = faceSetMap_.find(name);
    if (iter == faceSetMap_.end())
    {
      throw ATC_Error( "No faceset with name " + name + " found.");
    }
    return iter->second;
  }

  // -------------------------------------------------------------
  //   create_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::create_nodeset(const string & name,
                               const set<int> & nodeSet)
  {
    // Make sure we don't already have a nodeset with this name
    NODE_SET_MAP::iterator iter = nodeSetMap_.find(name);
    if (iter != nodeSetMap_.end()) {
      string message("A nodeset with name " + name + " is already defined.");
      throw ATC_Error( message);
    }
    nodeSetMap_[name] = nodeSet;

    if (ATC::LammpsInterface::instance()->rank_zero()) { 
      stringstream ss;
      ss << "created nodeset " << name 
         << " with " << nodeSet.size() << " nodes";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    }
  }
  void FE_Mesh::create_nodeset(const string & name,
                               double xmin, 
                               double xmax,
                               double ymin, 
                               double ymax,
                               double zmin, 
                               double zmax)
  {
    // Make sure we don't already have a nodeset with this name
    NODE_SET_MAP::iterator iter = nodeSetMap_.find(name);
    if (iter != nodeSetMap_.end()) {
      string message("A nodeset with name " + name + " is already defined.");
      throw ATC_Error( message);
    }

    set<int> nodeSet;

    // Loop over nodes and add their unique id's to the set if they're
    // in the correct range
    for (int inode = 0; inode < nNodes_; inode++) {
      double x = nodalCoords_(0,inode);
      double y = nodalCoords_(1,inode);
      double z = nodalCoords_(2,inode);
      if ( (xmin <= x) && (x <= xmax) &&
           (ymin <= y) && (y <= ymax) &&
           (zmin <= z) && (z <= zmax) ) {
        int uid = globalToUniqueMap_(inode);
        nodeSet.insert(uid);
      }
    }
    if (nodeSet.size() == 0) {
      string message("nodeset " + name + " has zero size.");
      throw ATC_Error( message);
    }

    nodeSetMap_[name] = nodeSet;

    if (ATC::LammpsInterface::instance()->rank_zero()) { 
      stringstream ss;
      ss << "created nodeset " << name 
         << " with " << nodeSet.size() << " nodes";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    }
  }

  // -------------------------------------------------------------
  //   add_to_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::add_to_nodeset(const string & name,
                               double xmin, 
                               double xmax,
                               double ymin, 
                               double ymax,
                               double zmin, 
                               double zmax)
  {
    // Make sure we already have a nodeset with this name
    NODE_SET_MAP::iterator iter = nodeSetMap_.find(name);
    if (iter == nodeSetMap_.end()) {
      string message("A nodeset with name " +name + " is not already defined.");
      throw ATC_Error( message);
    }

    set<int> nodeSet;

    // Loop over nodes and add their unique id's to the set if they're
    // in the correct range
    for (int inode = 0; inode < nNodes_; inode++) {
      double x = nodalCoords_(0,inode);
      double y = nodalCoords_(1,inode);
      double z = nodalCoords_(2,inode);
      if ( (xmin <= x) && (x <= xmax) &&
           (ymin <= y) && (y <= ymax) &&
           (zmin <= z) && (z <= zmax) ) {
        int uid = globalToUniqueMap_(inode);
        nodeSet.insert(uid);
      }
    }
    if (nodeSet.size() == 0) {
      string message("nodeset " + name + " has zero size.");
      throw ATC_Error( message);
    }

    nodeSetMap_[name].insert(nodeSet.begin(),nodeSet.end());

    if (ATC::LammpsInterface::instance()->rank_zero()) { 
      stringstream ss;
      ss   << "added " << nodeSet.size() << " nodes to nodeset " << name ;
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    }
  }

  // -------------------------------------------------------------
  //   create_faceset
  // -------------------------------------------------------------
  void FE_Mesh::create_faceset(const string & name,
                               double xmin, 
                               double xmax,
                               double ymin, 
                               double ymax,
                               double zmin, 
                               double zmax, 
                               bool outward)
  {
    // Make sure we don't already have a nodeset with this name
    FACE_SET_MAP::iterator iter = faceSetMap_.find(name);
    if (iter != faceSetMap_.end()) 
      throw ATC_Error( "A faceset with name " + name + " is already defined.");

    set<PAIR> faceSet;
    // Loop over face and add their unique id's to the set if they concide
    // with region
    const int nf = num_faces_per_element();
    const int npf = num_nodes_per_face();
    const Array2D<int> & face_conn = local_face_connectivity();
    for (int ielem = 0; ielem < nElts_; ielem++) 
    {
      for (int iface = 0; iface < nf; iface++) 
      {
        bool in = true;
        bool on_xmin = true, on_xmax = true; 
        bool on_ymin = true, on_ymax = true;
        bool on_zmin = true, on_zmax = true;
        bool x_neg = false, x_pos = false;
        bool y_neg = false, y_pos = false;
        bool z_neg = false, z_pos = false;
        double x,y,z;
        for (int inode = 0; inode < npf; inode++) 
        {
          x = nodalCoords_(0,connectivity_(face_conn(iface,inode),ielem));
          y = nodalCoords_(1,connectivity_(face_conn(iface,inode),ielem));
          z = nodalCoords_(2,connectivity_(face_conn(iface,inode),ielem));
          
          if ( x + tol < xmin) { in = false; break; }
          if ( x - tol > xmax) { in = false; break; }
          if ( y + tol < ymin) { in = false; break; }
          if ( y - tol > ymax) { in = false; break; }
          if ( z + tol < zmin) { in = false; break; }
          if ( z - tol > zmax) { in = false; break; }

          on_xmin = on_xmin &&  fabs(x-xmin) <= tol;
          on_xmax = on_xmax &&  fabs(x-xmax) <= tol;
          on_ymin = on_ymin &&  fabs(y-ymin) <= tol;
          on_ymax = on_ymax &&  fabs(y-ymax) <= tol;
          on_zmin = on_zmin &&  fabs(z-zmin) <= tol;
          on_zmax = on_zmax &&  fabs(z-zmax) <= tol;
        }
        if (in) {
          // note based on structured grid
          if (outward) 
          {
            if (on_xmin && iface==0) { x_neg = true;}
            if (on_xmax && iface==1) { x_pos = true;}
            if (on_ymin && iface==2) { y_neg = true;}
            if (on_ymax && iface==3) { y_pos = true;}
            if (on_zmin && iface==4) { z_neg = true;}
            if (on_zmax && iface==5) { z_pos = true;}
          }
          else 
          {
            if (on_xmin && iface==1) { x_pos = true;}
            if (on_xmax && iface==0) { x_neg = true;}
            if (on_ymin && iface==3) { y_pos = true;}
            if (on_ymax && iface==2) { y_neg = true;}
            if (on_zmin && iface==5) { z_pos = true;}
            if (on_zmax && iface==4) { z_neg = true;}
          }
  
          if (  (x_neg || x_pos) || (y_neg || y_pos) || (z_neg || z_pos) ) {
            PAIR face(ielem,iface);
            faceSet.insert(face);
          }
        }
      }
    }
    if (faceSet.empty()) throw ATC_Error( "faceset "+name+" is empty.");

    faceSetMap_[name] = faceSet;
    if (ATC::LammpsInterface::instance()->comm_rank() == 0) { 
      stringstream ss;
      ss   << "created faceset " << name 
           << " with " << faceSet.size() << " faces";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    }
  }

  void FE_Mesh::create_faceset(const string & name,
                               double xRef, 
                               int nIdx, int nSgn,
                               int nIdx2, double x2lo, double x2hi,
                               int nIdx3, double x3lo, double x3hi)
  {
    double xtol = tolerance(xRef);
    // Make sure we don't already have a faceset with this name
    FACE_SET_MAP::iterator iter = faceSetMap_.find(name);
    if (iter != faceSetMap_.end()) 
      throw ATC_Error( "A faceset with name "+name+" is already defined.");

    bool finite2 = (nIdx2 >= 0);
    bool finite3 = (nIdx3 >= 0);

    set<PAIR> faceSet;
    // Loop over faces i& add unique id's to the set if concide w/ plane
    int nf = num_faces_per_element();
    int npf = num_nodes_per_face();
    const Array2D<int> & face_conn = local_face_connectivity();
    for (int ielem = 0; ielem < nElts_; ielem++) 
    {
      for (int iface = 0; iface < nf; iface++) 
      {
        bool in = true;
        // all nodes must be on the plane
        for (int inode = 0; inode < npf; inode++) {
          int node = connectivity_(face_conn(iface,inode),ielem);
          double x = nodalCoords_(nIdx,node);
          if ( fabs(x-xRef) > xtol){ in = false; break;}
          if (finite2) {
            double y = nodalCoords_(nIdx2,node);
            if ( y < x2lo || y > x2hi){ in = false; break;}
          }
          if (finite3) {
            double y = nodalCoords_(nIdx3,node);
            if ( y < x3lo || y > x3hi){ in = false; break;}
          }
        }
        // check correct orientation
        if (in) 
        {
          if ( (nIdx == 0 && iface==0 && nSgn == -1) 
            || (nIdx == 0 && iface==1 && nSgn ==  1)
            || (nIdx == 1 && iface==2 && nSgn == -1)
            || (nIdx == 1 && iface==3 && nSgn ==  1)
            || (nIdx == 2 && iface==4 && nSgn == -1)
            || (nIdx == 3 && iface==5 && nSgn ==  1) )
          {
            PAIR face(ielem,iface);
            faceSet.insert(face);
          }
        }
      }
    }

    if (faceSet.empty()) 
      throw ATC_Error( "faceset "+name+" is empty.");

    faceSetMap_[name] = faceSet;
    if (ATC::LammpsInterface::instance()->comm_rank() == 0) { 
      stringstream ss;
      ss   << "created faceset " << name 
           << " with " << faceSet.size() << " faces";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    }
  }

  // -------------------------------------------------------------
  //   create_elementset
  // -------------------------------------------------------------
  void FE_Mesh::create_elementset(const string & name,
                               double xmin, 
                               double xmax,
                               double ymin, 
                               double ymax,
                               double zmin, 
                               double zmax)
  {
    // Make sure we don't already have a elementset with this name
    ELEMENT_SET_MAP::iterator iter = elementSetMap_.find(name);
    if (iter != elementSetMap_.end()) {
      string message("An elementset with name "+name+" is already defined.");
      throw ATC_Error( message);
    }

    set<int> nodeSet;

    // Loop over nodes and add their unique id's to the set if they're
    // in the correct range
    for (int inode = 0; inode < nNodes_; inode++) {
      double x = nodalCoords_(0,inode);
      double y = nodalCoords_(1,inode);
      double z = nodalCoords_(2,inode);
      if ( (xmin <= x) && (x <= xmax) &&
           (ymin <= y) && (y <= ymax) &&
           (zmin <= z) && (z <= zmax) ) {
        int uid = globalToUniqueMap_(inode);
        nodeSet.insert(uid);
      }
    }
    if (nodeSet.size() == 0) {
      string message("elementset " + name + " has zero size.");
      throw ATC_Error( message);
    }

    // create a minimal element set from all the nodes included in the region
    set<int> elemSet;
    int npe = num_nodes_per_element();
    for (int ielem=0; ielem < nElts_; ielem++) 
    {
      int inode = 0;
      bool in = true;
      while (in && inode < npe) 
      {
        int node = connectivityUnique_(inode, ielem);
        set<int>::const_iterator iter = nodeSet.find(node);
        if (iter == nodeSet.end()) { in=false; }
        inode++;
      }
      if (in) elemSet.insert(ielem);
    }
    if (elemSet.size() == 0) {
      string message("element set " + name + " has zero size.");
      throw ATC_Error( message);
    }
    elementSetMap_[name] = elemSet;

    if (ATC::LammpsInterface::instance()->comm_rank() == 0) { 
      stringstream ss;
      ss   << "created elementset " << name 
           << " with " << elemSet.size() << " elements";
      ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    }
  }
  // -------------------------------------------------------------
  //   get_numIPsPerElement()
  // -------------------------------------------------------------
  int FE_Mesh::num_ips_per_element() const
  {
    return feElement_->num_ips();
  }
  // -------------------------------------------------------------
  //   get_numNodesPerElement()
  // -------------------------------------------------------------
  int FE_Mesh::num_nodes_per_element() const
  {
    return feElement_->num_elt_nodes();
  }
  // -------------------------------------------------------------
  //   get_numFacesPerElement()
  // -------------------------------------------------------------
  int FE_Mesh::num_faces_per_element() const
  {
    return feElement_->num_faces();
  }
  // -------------------------------------------------------------
  //   get_num_ips_per_face()
  // -------------------------------------------------------------
  int FE_Mesh::num_ips_per_face() const
  {
    return feElement_->num_face_ips();
  }
  // -------------------------------------------------------------
  //   get_num_nodes_per_face()
  // -------------------------------------------------------------
  int FE_Mesh::num_nodes_per_face() const
  {
    return feElement_->num_face_nodes();
  }
  // -------------------------------------------------------------
  //   mappings from element id to associated nodes
  // -------------------------------------------------------------
  
  
  void FE_Mesh::element_connectivity_global(const int eltID,
                                            Array<int> & nodes) const
  {
    const int npe = num_nodes_per_element();
    nodes.reset(npe);

    // use connectivity arrays
    if (decomposition_ && partitioned_) {
      for (int inode = 0; inode < npe; inode++) {
        nodes(inode) = myConnectivity_(inode, map_elem_to_myElem(eltID)); 
      }
    } else {
      for (int inode = 0; inode < npe; inode++) {
        nodes(inode) = connectivity_(inode, eltID); 
      }
    }
  }
  // -------------------------------------------------------------
  //   
  // -------------------------------------------------------------
  void FE_Mesh::element_connectivity_unique(const int eltID,
                                            Array<int> & nodes) const
  {
    const int npe = num_nodes_per_element();
    nodes.reset(npe);

    // use connectivity arrays
    if (decomposition_ && partitioned_) {
      for (int inode = 0; inode < npe; inode++) {
        nodes(inode) = myConnectivityUnique_(inode, map_elem_to_myElem(eltID)); 
      }
    } else {
      for (int inode = 0; inode < npe; inode++) {
        nodes(inode) = connectivityUnique_(inode, eltID); 
      }
    }
  }

  // -------------------------------------------------------------
  //  
  // -------------------------------------------------------------
  int FE_Mesh::element_connectivity_global(const int eltID,
                                            const int inode) const
  {
    if (decomposition_ && partitioned_) {
      return myConnectivity_(inode, map_elem_to_myElem(eltID));
    } else {
      return connectivity_(inode, eltID);
    }
  }
  // -------------------------------------------------------------
  //   
  // -------------------------------------------------------------
  int FE_Mesh::element_connectivity_unique(const int eltID,
                                            const int inode) const
  {
    if (decomposition_ && partitioned_) {
      return myConnectivityUnique_(inode, map_elem_to_myElem(eltID));
    } else {
      return connectivityUnique_(inode, eltID);
    }
  }
  // -------------------------------------------------------------
  //  
  // -------------------------------------------------------------
  AliasArray<int> FE_Mesh::element_connectivity_global(const int eltID) const
  {
    if (decomposition_ && partitioned_) {
      return myConnectivity_.column(map_elem_to_myElem(eltID));
    } else {
      return connectivity_.column(eltID);
    }
  }
  // -------------------------------------------------------------
  //   
  // -------------------------------------------------------------
  AliasArray<int> FE_Mesh::element_connectivity_unique(const int eltID) const
  {
    if (decomposition_ && partitioned_) {
      return myConnectivityUnique_.column(map_elem_to_myElem(eltID));
    } else {
      return connectivityUnique_.column(eltID);
    }
  }
  // -------------------------------------------------------------
  //   local_face_connectivity()
  // -------------------------------------------------------------
  const Array2D<int> &FE_Mesh::local_face_connectivity() const
  {
    return feElement_->local_face_conn();
  }

  // -------------------------------------------------------------
  //   maps to/from partitioned element data
  // -------------------------------------------------------------

  int FE_Mesh::map_elem_to_myElem(int elemID) const
  {
    
    return elemToMyElemMap_.find(elemID)->second;
  }

  int FE_Mesh::map_myElem_to_elem(int myElemID) const
  {
    return myElts_[myElemID]; 
  }
  
  // -------------------------------------------------------------
  //   shape function evaluation
  // -------------------------------------------------------------
  
  // set quadrature scheme pass-through
  void FE_Mesh::set_quadrature(FeIntQuadrature type) 
  { 
    feElement_->set_quadrature(type); 
  }

  // shape function evaluation
  void FE_Mesh::shape_functions(const VECTOR &x,
                                DENS_VEC &N,
                                Array<int> &nodeList) const
  {
    // get element id from global coordinates
    int eltID = map_to_element(x);
    
    // call appropriate function below, with eltID
    shape_functions(x,eltID,N,nodeList);
  }

  void FE_Mesh::shape_functions(const DENS_VEC &x, 
                                DENS_VEC &N,
                                DENS_MAT &dNdx,
                                Array<int> &nodeList) const
  {
    // get element id from global coordinates
    int eltID = map_to_element(x);

    // call appropriate function below, with eltID
    shape_functions(x,eltID,N,dNdx,nodeList);
  }

  void FE_Mesh::shape_functions(const VECTOR &x,
                                const int eltID,
                                DENS_VEC &N,
                                Array<int> &nodeList) const
  {
    // Get element node coordinates from mesh
    DENS_MAT eltCoords;
    element_coordinates(eltID, eltCoords);
    
    // pass through call
    feElement_->shape_function(eltCoords,x,N);
    
    // determine nodes which correspond to shape function indices
    element_connectivity_unique(eltID,nodeList);
  }

  void FE_Mesh::shape_functions(const DENS_VEC &x,
                                const int eltID,
                                DENS_VEC &N,
                                DENS_MAT &dNdx,
                                Array<int> &nodeList) const
  {
    // Get element node coordinates from mesh
    DENS_MAT eltCoords;
    element_coordinates(eltID,eltCoords);

    // pass through call
    feElement_->shape_function(eltCoords,x,N,dNdx);
    
    // determine nodes which correspond to shp function indices
    element_connectivity_unique(eltID,nodeList);
  }

  void FE_Mesh::shape_function_derivatives(const DENS_VEC &x,
                                const int eltID,
                                DENS_MAT &dNdx,
                                Array<int> &nodeList) const
  {
    // Get element node coordinates from mesh
    DENS_MAT eltCoords;
    element_coordinates(eltID,eltCoords);

    // pass through call
    feElement_->shape_function_derivatives(eltCoords,x,dNdx);
    
    // determine nodes which correspond to shp function indices
    element_connectivity_unique(eltID,nodeList);
  }

  void FE_Mesh::shape_function(const int eltID,
                               DENS_MAT &N, 
                               DIAG_MAT &weights) const
  {
    // unused data (but required to calc weights)
    vector<DENS_MAT> dN;

    // call below function with dN to avoid duplicated code
    shape_function(eltID,N,dN,weights);
  }

  void FE_Mesh::shape_function(int eltID,
                               DENS_MAT &N, 
                               vector<DENS_MAT> &dN,
                               DIAG_MAT &weights) const
  {
    // Get element node coordinates from mesh
    DENS_MAT eltCoords;
    element_coordinates(eltID,eltCoords);

    // pass through call
    feElement_->shape_function(eltCoords,N,dN,weights);
  }

  void FE_Mesh::face_shape_function(const PAIR &face,
                                    DENS_MAT &N,
                                    DENS_MAT &n,
                                    DIAG_MAT &weights) const
  {
    int eltID = face.first;
    int faceID = face.second;
    
    // Get element node coordinates from mesh
    DENS_MAT eltCoords;
    element_coordinates(eltID,eltCoords);

    // pass through call
    feElement_->face_shape_function(eltCoords,faceID,N,n,weights);
  }

  void FE_Mesh::face_shape_function(const PAIR &face,
                                    DENS_MAT &N,
                                    vector<DENS_MAT> &dN,
                                    vector<DENS_MAT> &Nn,
                                    DIAG_MAT &weights) const
  {
    int eltID = face.first;
    int faceID = face.second;
    
    // Get element node coordinates from mesh
    DENS_MAT eltCoords;
    element_coordinates(eltID,eltCoords);

    // pass through call
    feElement_->face_shape_function(eltCoords,faceID,N,dN,Nn,weights);
  }

  double FE_Mesh::face_normal(const PAIR &face,
                              int ip,
                              DENS_VEC &normal)  const
  {
    int eltID = face.first;
    int faceID = face.second;
    
    // Get element node coordinates from mesh
    DENS_MAT eltCoords;
    element_coordinates(eltID,eltCoords);

    // pass through call
    double J = feElement_->face_normal(eltCoords,faceID,ip,normal);
    return J;
  }

  //-----------------------------------------------------------------------
  void FE_Mesh::output(string prefix) const
  {
    set<int> otypes;
    otypes.insert(ENSIGHT);
    otypes.insert(FULL_GNUPLOT);
    OutputManager meshSets(prefix,otypes);
    meshSets.write_geometry(&nodalCoords_, & connectivity_);
    OUTPUT_LIST subsetData;
    int size = nNodesUnique_;
    // material
//  DENS_MAT material(nNodes_,1);
//  material = 1;
//  subsetData["material"] = &material;
    //string name = "global_to_unique_map";
    // nodesets
    NODE_SET_MAP::const_iterator niter;
    map<string,DENS_MAT> nodesets;
    for (niter = nodeSetMap_.begin(); niter != nodeSetMap_.end(); niter++) {
      string name = niter->first;
      const set<int> & nset = niter->second;
      string nodeset = "nodeset_"+name;
      nodesets[nodeset].reset(size,1);
      set<int>::const_iterator iter;
      for (iter = nset.begin(); iter != nset.end(); iter++) {
        (nodesets[nodeset])(*iter,0) = 1;  
      }
      subsetData[nodeset] = & nodesets[nodeset];
    }
    // facesets
    FACE_SET_MAP::const_iterator fiter;
    map<string,DENS_MAT> facesets;
    for (fiter = faceSetMap_.begin(); fiter != faceSetMap_.end(); fiter++) {
      string name = fiter->first;
      string faceset = "faceset_"+name;
      facesets[faceset].reset(size,1);
      set<int> nset;
      faceset_to_nodeset(name,nset);
      set<int>::const_iterator iter;
      for (iter = nset.begin(); iter != nset.end(); iter++) {
        (facesets[faceset])(*iter,0) = 1;
      }
      subsetData[faceset] = & facesets[faceset];
    }
    // elementsets
    ELEMENT_SET_MAP::const_iterator eiter;
    map<string,DENS_MAT> elemsets;
    for (eiter = elementSetMap_.begin();eiter != elementSetMap_.end();eiter++) {
      string name = eiter->first;
      const set<int> & eset = eiter->second;
      string elemset = "elementset_"+name;
      elemsets[elemset].reset(size,1);
      set<int>::const_iterator iter;
      for (iter = eset.begin(); iter != eset.end(); iter++) {
        Array<int> nodes;
        element_connectivity_unique(*iter, nodes);
        for(int i = 0; i < nodes.size(); ++i) {
          (elemsets[elemset])(nodes(i),0) = 1; 
        }
      }
      subsetData[elemset] = & elemsets[elemset];
    }
    // output
    const int * map = globalToUniqueMap_.data();
    if (subsetData.size() > 0 )
      meshSets.write_data(0.0,&subsetData,map);
    else
      if (LammpsInterface::instance()->comm_rank() == 0) {
        stringstream ss;
        ss << "Warning mesh output requested without any mesh entities, output suppressed";
        ATC::LammpsInterface::instance()->print_msg(ss.str());
      }
  }

  bool FE_Mesh::is_owned_elt(int elt) const
  {
    return (find(myElts_.begin(), myElts_.end(), elt) != myElts_.end()); 
  }


  // -------------------------------------------------------------
  // -------------------------------------------------------------
  //   class FE_3DMesh
  // -------------------------------------------------------------
  // -------------------------------------------------------------
  FE_3DMesh::FE_3DMesh(const string elementType,
                       const int nNodes, 
                       const int nElements,
                       const Array2D<int> *connectivity, 
                       const DENS_MAT *nodalCoordinates,
                       const Array<bool> periodicity,
                       const Array< pair< string, set<int> > > *nodeSets):
    FE_Mesh(),
    minEltSize_(0),
    tree_(NULL)
  {
    // Pick which element class to make
    if (elementType == "HEX8") {
      feElement_ = new FE_ElementHex(8,4,2);
    } else if (elementType == "HEX20") {
      feElement_ = new FE_ElementHex(20,8,3);
    } else if (elementType == "HEX27") {
      feElement_ = new FE_ElementHex(27,9,3);
    } else if (elementType == "TET4") {
      feElement_ = new FE_ElementTet(4,3,2);
      hasPlanarFaces_ = true;
    } else {
      throw ATC_Error("Unrecognized element type specified.");
    }

    
    nSD_ = 3;
    nNodes_ = nNodes;
    nNodesUnique_ = nNodes;
    nElts_ = nElements;
    xscale_ = yscale_ = zscale_ = 1;
    periodicity_ = periodicity;

    // Connectivity and coordinates describe the mesh geometry.
    connectivity_.reset(connectivity->nRows(),connectivity->nCols());
    connectivity_ = (*connectivity);
    nodalCoords_.reset(nodalCoordinates->nRows(),nodalCoordinates->nCols());
    nodalCoords_ = (*nodalCoordinates);

    // set minimum element size
    minEltSize_ = 1.e20; 
    for (int i=0; i< connectivity_.nCols(); ++i) {
      int n1 = connectivity_(0,i);
      int n2 = connectivity_(1,i);
      double dx[3] = {fabs(nodalCoords_(0,n1)-nodalCoords_(0,n2)),
                      fabs(nodalCoords_(1,n1)-nodalCoords_(1,n2)),
                      fabs(nodalCoords_(2,n1)-nodalCoords_(2,n2))};
      minEltSize_ = min(minEltSize_,norm3(dx));
    }

    // create global-unique maps
    coordTol_ = this->coordinate_tolerance();
    setup_periodicity();

    // Create the "all" elementset, the "all" nodeset, and the read-in nodesets.
    for (int elem = 0; elem < nElts_; elem++) elementSetAll_.insert(elem); 
    for (int node = 0; node < nNodesUnique_; node++) nodeSetAll_.insert(node);
    const Array<pair<string,set<int> > > & sets = *nodeSets;
    for (int nodeSet = 0; nodeSet < sets.size(); ++nodeSet) {
      const set<int> & nset = sets(nodeSet).second;
      set<int> copy; 
      if (compactRemap_.size() > 0) {
        for (set<int>::iterator itr = nset.begin(); itr != nset.end(); itr++) {
          copy.insert(globalToUniqueMap_(compactRemap_(*itr)));
        }
      }
      else {
        for (set<int>::iterator itr = nset.begin(); itr != nset.end(); itr++) {
          copy.insert(globalToUniqueMap_(*itr));
        }
      }
      create_nodeset(sets(nodeSet).first, copy);
    }
    
    // Insert nodes and elements into KD-tree for PIE search.
    if (tree_ == NULL) {
      tree_ = KD_Tree::create_KD_tree(feElement_->num_elt_nodes(), nNodes_, 
        &nodalCoords_, nElts_, connectivity_);
    }

  }

  FE_3DMesh::~FE_3DMesh() {
    if (tree_) delete tree_; 
  }

  // -------------------------------------------------------------
  //  setup_periodicity
  // -------------------------------------------------------------
  void FE_3DMesh::setup_periodicity(double tol)
  {
    // unique <-> global id maps
    globalToUniqueMap_.reset(nNodes_);
    uniqueToGlobalMap_.reset(nNodes_);

    for (int node = 0; node < nNodes_; node++) {
      globalToUniqueMap_(node) = node;
    }

    // recursive fix and setup global to unique map
    if (periodicity_(2)) { fix_periodicity(3); }
    if (periodicity_(1)) { fix_periodicity(2); }
    if (periodicity_(0)) { fix_periodicity(1); }
    

    // renumber to compact unique numbering
    // unique nodes map to the same id with global to unique
    if (nNodesUnique_ < nNodes_) {
      int n = 0;
      int m = nNodesUnique_;
      compactRemap_.reset(nNodes_); // old to new global numbering
      compactRemap_ = -1;
      for (int node = 0; node < nNodes_; node++) {
        if (globalToUniqueMap_(node) == node) { // unique nodes
          compactRemap_(node) = n++;
        }
        else { // periodic nodes
          compactRemap_(node) = m++;
        }
      }
      if (n != nNodesUnique_) throw ATC_Error("didn't compact numbering");
      int npe = num_nodes_per_element();
      // temporary copy
      DENS_MAT  coor = DENS_MAT(nodalCoords_);
      Array2D<int> conn(connectivity_);
      Array<int> oldGlobalToUniqueMap = globalToUniqueMap_;
      for (int node = 0; node < nNodes_; node++) {
        // unique = remap * global2unique  global
        globalToUniqueMap_(compactRemap_(node)) = compactRemap_(oldGlobalToUniqueMap(node));
        if (compactRemap_(node) < 0) throw ATC_Error("mis-map to compact numbering");
        // swap coordinates
        for (int i = 0; i < nSD_; i++) {
          nodalCoords_(i,compactRemap_(node)) = coor(i,node);
        }
      }
      for (int elem = 0; elem < nElts_; elem++) {
        for (int i = 0; i < npe; i++) {
          connectivity_(i,elem) = compactRemap_(conn(i,elem));
        }
      }
    }

    for (int node = 0; node < nNodes_; node++) {
      uniqueToGlobalMap_(globalToUniqueMap_(node)) = node;
    }
    set_unique_connectivity();
    //amended_ = true; // would always be true since setup_always called
  }

  void FE_3DMesh::fix_periodicity(int idir)
  {
    set<int> topNodes,botNodes;
    int ntop = find_boundary_nodes( idir,topNodes);
    int nbot = find_boundary_nodes(-idir,botNodes);
    if (ntop != nbot)
      throw ATC_Error("can't fix periodicity, number of top and bottom nodes are different ");
    bool match = match_nodes(idir,topNodes,botNodes,globalToUniqueMap_);
    if (!match) {
      stringstream ss;
      ss << "can't match periodic nodes with tolerance " << coordTol_;
      throw ATC_Error(ss.str());
    }
  }

  int FE_3DMesh::find_boundary_nodes(int idir, set<int> & nodes)
  {
    nodes.clear();
    double limit = 0;
    int idm = abs(idir)-1;
    DENS_MAT & x = nodalCoords_;
    if (idir > 0) limit = x.row_max(idm);
    else          limit = x.row_min(idm);
    for (int i=0; i < x.nCols(); ++i) {
      double xi = x(idm,i);
      if (fabs(xi-limit) < coordTol_) nodes.insert(i);
    }
//  stringstream ss;
//  ss << "found " << nodes.size() << " nodes at x_" << abs(idir) << "= "<< limit ;
//  ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    return nodes.size();
  }

  bool FE_3DMesh::match_nodes(int idir, set<int> & nodes1, set<int> & nodes2,
     Array<int> & map)
  {
    int i1=0,i2=1;
    plane_coords(idir-1,i1,i2);
    vector<bool> found(nNodes_,false);
    DENS_MAT & x = nodalCoords_;
    for (set<int>::iterator it1=nodes1.begin(); it1 != nodes1.end(); it1++) {
      int n1 = *it1;
      double x1 = x(i1,n1);
      double x2 = x(i2,n1);
      for (set<int>::iterator it2=nodes2.begin(); it2 != nodes2.end(); it2++) {
        int n2 = *it2;
        if (!found[n2]) {
          double y1 = x(i1,n2);
          double y2 = x(i2,n2);
          if (fabs(x1-y1) < coordTol_ && fabs(x2-y2) < coordTol_) {
            map(n1) = n2;
            found[n2] = true;
            // coincidence
            x(i1,n2) = x1;
            x(i2,n2) = x2;
          }
        }
      }
      if (map(n1) == n1) return false;
    }
    nNodesUnique_ -= nodes1.size();
    stringstream ss;
    ss << "condensed " << nodes1.size() << " periodic nodes in the " << abs(idir) << " direction";
    ATC::LammpsInterface::instance()->print_msg_once(ss.str());
    return true;
  }
  
  void FE_3DMesh::set_unique_connectivity(void) 
  {
    int numEltNodes = feElement_->num_elt_nodes();
    connectivityUnique_.reset(numEltNodes, nElts_);
    uniqueNodeToElementMap_.reset(nNodes_);
    for (int node = 0; node < nNodes_; node++) {
      uniqueNodeToElementMap_(node) = vector<int>();
    }
    for (int elem = 0; elem < nElts_; elem++) {
      for (int node = 0; node < numEltNodes; node++) {
        int global_node = connectivity_(node, elem);
        int unique_node = globalToUniqueMap_(global_node);
        connectivityUnique_(node,elem) = unique_node;
        uniqueNodeToElementMap_(unique_node).push_back(elem);
      }
    }
  }

  // orient the local coordinate with the global one
  bool FE_3DMesh::orient(int idm)
  {
    int numEltNodes = feElement_->num_elt_nodes();
    if (numEltNodes != 8) throw ATC_Error("can't currently orient non HEX8 elements");
    DENS_MAT x;
    for (int elem = 0; elem < nElts_; elem++) {
      element_coordinates(elem,x);
      double xmax = x.row_max(idm);
      double xmin = x.row_min(idm);
      set<int> top,bot;
      for (int node = 0; node < numEltNodes; node++) {
        // find top nodes
        if      ((x(idm,node) - xmax) < coordTol_) top.insert(node);
        // find bottom nodes
        else if ((x(idm,node) - xmin) < coordTol_) bot.insert(node);
        else return false;
      }
      // order by rh rule
      // start with one
    }
    throw ATC_Error("not completely implemented function: FE_3DMesh::orient");
    return true;
  }

  // -------------------------------------------------------------
  //  amend mesh for cut at specified faces
  // -------------------------------------------------------------
  void FE_3DMesh::cut_mesh(const set<PAIR> & faceSet, const set<int> & nodeSet)
  {
    int nsd = nSD_;
    // get all nodes to create an old to new number map
    set<int> dupNodes;
    map<int,int> oldToNewMap;
    faceset_to_nodeset_global(faceSet,dupNodes);
    // remove edge nodes
    Array<int> & node_map = globalToUniqueMap_;
    set<int>::const_iterator itr;
    for (itr = dupNodes.begin(); itr != dupNodes.end(); itr++) {
      int gnode = *itr;
      int unode = node_map(gnode);
      if (nodeSet.find(unode) != nodeSet.end()) {
        dupNodes.erase(gnode);
      }
    }
    int nNodesAdded = dupNodes.size();
    // copy coordinates and nodemap
    DENS_MAT &coordinates = nodalCoords_;
    int nNodesNew = coordinates.nCols() + nNodesAdded;
    coordinates.resize(nsd,nNodesNew,true);
    node_map.resize(nNodesNew,true);
    // add duplicate coordinates
    int iNew = nNodes_;
    int iNewUnique = nNodesUnique_;
    for (itr = dupNodes.begin(); itr != dupNodes.end(); 
         itr++,iNew++) {
      int iOld = *itr;
      oldToNewMap[iOld] = iNew;      // global ids
      if (iOld == node_map(iOld)) {  // non-image atom
        node_map(iNew) = iNewUnique++; 
      } else { 
        node_map(iNew) = -1; 
      }
      for(int j = 0; j < nsd; j++) {
        coordinates(j,iNew) = coordinates(j,iOld);
      }
    }
    nNodes_ = iNew;
    nNodesUnique_ = iNewUnique;
    for (itr = dupNodes.begin(); itr != dupNodes.end(); 
         itr++,iNew++) {
      int iOld = *itr;
      iNew = oldToNewMap[iOld];        // global ids
      int iOldImage = node_map(iOld);
      if (iOld != iOldImage) {         // image atom
        int iNewImage = oldToNewMap[iOldImage];
        node_map(iNew) = node_map(iNewImage);
      }
    }
    // update element connectivities
    const int nnf = num_nodes_per_face();
    const Array2D <int> & local_conn = feElement_->local_face_conn();
    set< PAIR >::iterator iter;
    for (iter = faceSet.begin(); iter != faceSet.end(); iter++) 
    {
      PAIR face = *iter;
      int eltID=face.first, faceID=face.second;
      for (int inode=0; inode < nnf; inode++) {
        int lid = local_conn(faceID,inode);
        int id = connectivity_(lid, eltID);
        if (oldToNewMap.find(id) != oldToNewMap.end() ) {
          int new_id = (*oldToNewMap.find(id)).second;
          connectivity_(lid,eltID)        = new_id;
          connectivityUnique_(lid, eltID) = node_map(new_id); 
        }
      }
    }
  }


  // -------------------------------------------------------------
  //  amend mesh for deleted elements
  // -------------------------------------------------------------
  void FE_3DMesh::delete_elements(const set<int> &elementList)
  {
    int nPE = num_nodes_per_element();
    set<int> elementsNew;
    elementset_complement(elementList,elementsNew);
    int nElementsNew = elementsNew.size();
    set<int> newToOld;
    map<int,int> oldToNewMap;
    elementset_to_nodeset(elementsNew,newToOld);
    int nNodesNew = newToOld.size(); 
    set<int>::const_iterator itr;

    // coordinates & node map (from nodes to data)
    const DENS_MAT &coordinates = nodal_coordinates();
    const Array<int> & node_map = global_to_unique_map();

    DENS_MAT *newCoords = new DENS_MAT(nSD_,nNodesNew);
    Array<int> *newNodeMap = new Array<int> (nNodesNew);
    Array2D<int> * newConnectivity = new Array2D<int>(nPE,nElementsNew);
    int k = 0, i = 0;
    for (itr = newToOld.begin(); itr != newToOld.end(); itr++) {
      int node = *itr;
      oldToNewMap[node] = i++;
      (*newNodeMap)(k) = node_map(node);
      for(int j = 0; j < nSD_; j++) {
        (*newCoords)(j,k) = coordinates(j,node);
      }
      k++;
    }
    // nNodes_ = nNodesNew; ???
    // connectivity
    k = 0;
    for (itr = elementsNew.begin(); itr != elementsNew.end(); itr++) {
      int ielem = *itr;
      for(int j = 0; j < nPE; j++) {
        int old_node = connectivity_(j,ielem);
        map<int,int>::iterator map_itr = oldToNewMap.find(old_node);
        if (map_itr == oldToNewMap.end()) { 
          stringstream ss; 
          ss << "map failure " << old_node << "\n"; 
          ATC::LammpsInterface::instance()->print_msg(ss.str());
        }
        int node = map_itr->second;
        (*newConnectivity)(j,k) = node;
      }
      k++;
    }
    delete newCoords;
    delete newNodeMap;
    delete newConnectivity;
  }

  // -------------------------------------------------------------
  //  partition_mesh 
  // -------------------------------------------------------------
  
  void FE_3DMesh::partition_mesh()
  {
    if (lammpsPartition_) {
      lammps_partition_mesh();
    }

    if (partitioned_) return;


    int nProcs, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs); // get the number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // get the current processor's rank

    // use the KD tree for partitioning, getting more blocks than
    // processors
    if (tree_ == NULL) {
      tree_ = KD_Tree::create_KD_tree(feElement_->num_elt_nodes(), 
                                      nNodes_, &nodalCoords_, 
                                      nElts_, connectivity_);
    }

    // Get the divisions of the elements from the KD-tree.
    int depth = ceil(log(nProcs)/log(2));
      // In order to make sure there are enough divisions to evenly
      // divide between all processors, we get the next-highest
      // power of 2.
    vector<vector<int> > procEltLists = tree_->getElemIDs(depth);
    int numEltLists = procEltLists.size();
    
    // Make sure the KD tree is behaving as expected.
    assert(numEltLists >= nProcs);
 
    // If the KD-tree was not able to return enough divisions,
    // duplicate the largest list.
      
      // elements, then the division would be more even. 
    vector<vector<int> >::iterator it;
    if (numNonempty(procEltLists) < nProcs) {
      // Find the list of maximum size and assign it to empty processors
      vector<int> maxSizeList (0);
      for (it = procEltLists.begin(); it != procEltLists.end(); ++it) {
        if (it->size() > maxSizeList.size())
          maxSizeList.assign(it->begin(), it->end());
      }
      for (it = procEltLists.begin(); it != procEltLists.end(); ++it) {
        if (it->empty()) {
          if (numNonempty(procEltLists) >= nProcs) break;
          it->assign(maxSizeList.begin(), maxSizeList.end());
        }
      }
    }
   
    // We will store the owning processor for each element.
    int * eltToOwners = new int[nElts_];
    for (int i = 0; i < nElts_; ++i) {
      eltToOwners[i] = -1; // -1 here means the element is unowned.
    }

    // Prune elements that appear on more than one processor.
    // Simultaneously fill in the ownership array.
    prune_duplicate_elements( procEltLists, eltToOwners);

    // If we have more lists than processors, get rid of the
    // extras and redistribute their elements.
    if (numNonempty(procEltLists) > nProcs) {
      redistribute_extra_proclists(procEltLists, eltToOwners, nProcs);
    }

    // Sort the lists so that the fuller ones are in the front.
    sort(procEltLists.begin(), procEltLists.end(), vectorCompSize);

    // Assign each processor a list of elements.
    myElts_ = procEltLists[myrank];

    //mesh_distribute(eltStartIDs);
    delete[] eltToOwners;

    // We should do nodes later.

    if (decomposition_) distribute_mesh_data();
    partitioned_ = true;
  }

  void FE_3DMesh::departition_mesh()
  {
    if (!partitioned_) return;
    partitioned_ = false;
  }

  void FE_3DMesh::prune_duplicate_elements(vector<vector<int> > & procEltLists, 
                                           int * eltToOwners)
  {
    int procID = 0;
    vector<vector<int> >::iterator topIt;
    vector<int> * conflictingProc;
    vector<int>::iterator it, toErase;
    for (topIt = procEltLists.begin(); topIt != procEltLists.end(); ++topIt) {
      // Simultaneously fill in eltToOwners and prune, if an element
      // appears on multiple processors.
      for (it = topIt->begin(); it != topIt->end(); ++it) {
        // If the element has no corresponding processor in eltToOwners,
        // record it as belonging to processor *it.
        if (eltToOwners[*it] == -1) {
          eltToOwners[*it] = procID;
        } 
        else {
          // If it does have a processor in eltToOwners, then we need
          // to remove it from either processor *topIt or from the processor
          // listed in eltToOwners. We discriminate based on size.
          conflictingProc = &(procEltLists[eltToOwners[*it]]);
          if (conflictingProc->size() <= procEltLists[procID].size()) {
            // Delete element from processor *topIt, if it has more elements.
            it = topIt->erase(it);
            --it;
          } 
          else {
            // Delete element from conflicting processor otherwise.
            toErase = find(conflictingProc->begin(), conflictingProc->end(), *it);
            conflictingProc->erase(toErase);
            eltToOwners[*it] = procID;
          }
        }
      }
      ++procID;
    }
  }

  // -------------------------------------------------------------
  //  lammps_partition_mesh 
  // -------------------------------------------------------------
  
  void FE_3DMesh::lammps_partition_mesh()
  {
    if (LammpsInterface::instance()->domain_triclinic()) {
      LammpsInterface::instance()->print_msg_once("Cannot use lammps partitioning, domain is triclinic");
      return;
    }
    LammpsInterface::instance()->print_msg_once("Warning: Using native lammps partitioning");
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double boxxlo, boxxhi, boxylo, boxyhi, boxzlo, boxzhi;
    LammpsInterface::instance()->sub_bounds(xlo, xhi, ylo, yhi, zlo, zhi);
    LammpsInterface::instance()->box_bounds(boxxlo, boxxhi, boxylo, boxyhi, boxzlo, boxzhi);

    
    myElts_.clear();
    double xCent, yCent, zCent;
    
    // Assign elements to processors based on the centroid of the element.
    int numNodes = num_nodes_per_element();
    for (int i = 0; i < nElts_; ++i)
    {
      xCent = 0.0;
      yCent = 0.0;
      zCent = 0.0;
      for (int j = 0; j < numNodes; ++ j)
      {
        xCent += nodalCoords_(0, connectivity_(j,i));
        yCent += nodalCoords_(1, connectivity_(j,i));
        zCent += nodalCoords_(2, connectivity_(j,i));
      }
      xCent /= numNodes;
      yCent /= numNodes;
      zCent /= numNodes;
      if (xCent < boxxlo) xCent = boxxlo;
      if (xCent < boxxhi) xCent = boxxhi;
      if (yCent < boxylo) yCent = boxylo;
      if (yCent < boxyhi) yCent = boxyhi;
      if (zCent < boxzlo) zCent = boxzlo;
      if (zCent < boxzhi) zCent = boxzhi; 
      if ( dbl_geq(xCent, xlo) &&
          ((xhi == boxxhi) || !dbl_geq(xCent, xhi)) &&
           dbl_geq(yCent, ylo) &&
          ((yhi == boxyhi) || !dbl_geq(yCent, yhi)) &&
           dbl_geq(zCent, zlo) &&
          ((zhi == boxzhi) || !dbl_geq(zCent, zhi))) {
        myElts_.push_back(i);
      }
    }


    // if decomposing add in myAtomElts list based on nodal locations, i.e., all elements with a local node
    // myElts: for FE assembly
    // myAndGhost : for atom ops like restrict
    if (decomposition_) {
      set<int> elms;
      vector<int>::const_iterator itr;
      for (itr=myElts_.begin(); itr!=myElts_.end(); itr++) {elms.insert(*itr); }
      set<int> nodes;
      elementset_to_nodeset(elms,nodes);
      elms.clear();
      nodeset_to_maximal_elementset(nodes,elms);
      myAndGhostElts_.clear();
      set<int>::const_iterator iter;
      for (iter=elms.begin(); iter!=elms.end(); iter++) 
        {myAndGhostElts_.push_back(*iter);}
      distribute_mesh_data();
    }
    partitioned_ = true;
    return;
    
  }

  void FE_3DMesh::redistribute_extra_proclists(vector<vector<int> > &procEltLists,
                                               int *eltToOwners, int nProcs)
  {
      DENS_MAT faceAdjacencies(nElts_, num_faces_per_element());
      faceAdjacencies = -1; // Set all values to -1, indicating uninitialized/uncalculated

      int currentElt, adjacentElt, procID;
      
      // Put all of the hobos onto one master list, allHomelessElts.
      list<int> allHomelessElts;
      vector<int> oneHomelessList;
      vector<vector<int> >::iterator current;
      int nHoboLists = numNonempty(procEltLists) - nProcs;
      for (int i = 0; i < nHoboLists; ++i) {
        current = min_element(procEltLists.begin(), procEltLists.end(), vectorCompSize);
        oneHomelessList = *current;
        allHomelessElts.insert(allHomelessElts.end(), 
                               oneHomelessList.begin(), oneHomelessList.end());
        current->clear();
      }
      
      // Make sure the hobos lose their association with their old processor.
      list<int>::iterator it;
      for (it = allHomelessElts.begin(); it != allHomelessElts.end(); it++){
        eltToOwners[*it] = -1;
      }

      // Figure out which elements the hobos are adjacent to. That way, they
      // will know what processors they can be redistributed to.
      compute_face_adjacencies(allHomelessElts, faceAdjacencies); 

      // Place homeless elements onto lists that correspond to actual processors.
      while (!allHomelessElts.empty()) {
        currentElt = allHomelessElts.back();
        
        // This will store the ID of the processor with the fewest elements
        // so far that has an element adjacent to currentElt.
        PAIR smallestProc(-1, INT_MAX);
        
        // Iterate over the faces, check the processors of adjacent elements,
        // and slate the element to go on the adjacent processor with the fewest 
        // elements.
        for (int localFaceID = 0; localFaceID < num_faces_per_element(); ++localFaceID) {
          adjacentElt = faceAdjacencies(currentElt, localFaceID); 
          
          // This means that there is no adjacency through this face.
          if (adjacentElt >= nElts_) continue;
          
          procID = eltToOwners[adjacentElt];
          // The procID > -1 check makes sure we're not adjacent to another
          // homeless element by this face, in which case it won't have a
          // processor to put currentElt onto yet.
          if (procID > -1 && ((int) procEltLists[procID].size()) < smallestProc.second) {
            smallestProc = PAIR(procID, procEltLists[procID].size());
          } 
        }
        
        allHomelessElts.pop_back();

        // If we couldn't find an adjacent element that had a processor,
        // skip for now and come back to it later.
        if (smallestProc.first == -1) {
          allHomelessElts.push_front(currentElt);
        }
        // Otherwise, put it onto the processor with the fewest elements that 
        // we found.
        else {
          procEltLists[smallestProc.first].push_back(currentElt);
          eltToOwners[currentElt] = smallestProc.first;
        }
      }
  }

  
  void FE_3DMesh::compute_face_adjacencies(const list<int> &elts, 
                                         DENS_MAT &faceAdjacencies)
  {
    list<int>::const_iterator cit;
    for (cit = elts.begin(); cit != elts.end(); ++cit) {
      // For each element, look at ever face, get the nodes on that face,
      // and find out which elements are in the intersection of elements
      // containing that node. Those two elements are adjacent, and we
      // mark it as such in the faceAdjacencies array.
      for (int localFaceID = 0; localFaceID < num_faces_per_element(); ++localFaceID) {
        Array<int> faceNodes;
        face_connectivity(PAIR(*(cit), localFaceID), faceNodes);
        // Put the first node's elements into the accumulator to start.
        vector<int> vIntersect = uniqueNodeToElementMap_(faceNodes(0));
        vector<int> vCurrent;
        // set_intersect requires a vector large enough to contain the 
        // max possible intersect, which cannot be larger than the entirety 
        // of the first vector involved. 
        vector<int> vTemp(vIntersect.size(), -1);
        // Find the intersection of each of the nodes' element vectors.
        for (int ithOnFace = 1; ithOnFace < num_nodes_per_face(); ++ithOnFace) {
          vCurrent = uniqueNodeToElementMap_(faceNodes(ithOnFace)); // Vector of elements for this node
          set_intersection(vCurrent.begin(), vCurrent.end(),
                           vIntersect.begin(), vIntersect.end(), vTemp.begin());
          vIntersect = vTemp;
          // Because we initialized the vector to be larger than necessary, maybe,
          // we remove all of the meaningless values used in initialization.
          while (vIntersect.back() == -1)
            vIntersect.pop_back();
          vTemp.clear();
          vTemp.resize(vIntersect.size(),-1);
        }
        // This means there is an adjacent face, since there
        // are two elements sharing all of the nodes on the face.
        if (vIntersect.size() == 2) {
          // We want to choose the element id of NOT the current
          // element to be listed as the adjacency.
          
          //    well, but that requires more complicated memoization and 
          //    this doesn't take much extra time.
          if (*cit == vIntersect[0]) {
            faceAdjacencies(*cit, localFaceID) = vIntersect[1];
          } 
          else {
            faceAdjacencies(*cit, localFaceID) = vIntersect[0];
          }
        } 
        // This means the element is on the border.
        else if (vIntersect.size() == 1) {
          faceAdjacencies(*cit, localFaceID) = INT_MAX;
        } 
        else {
          // This should never ever happen! The nodes should at least
          // share one element, since they are all on one face!
          // There should also never be more than two elements on the
          // same face... that would defy mathematics and physics in 
          // every way.
        }
      }
    }
  }

  // Sometimes we want to count the number of vectors that actually 
  // have stuff in them. We use this.
  int FE_3DMesh::numNonempty(vector<vector<int> > & v)
  {
    int result = 0;
    vector<vector<int> >::iterator it;
    for (it = v.begin(); it != v.end(); ++it){
      if (!(it->empty())){
        result++;
      }
    }
    return result;
  }

  int FE_3DMesh::map_to_element(const DENS_VEC &query) const
  {
    DENS_MAT eltCoords;
    vector<int> candidates = tree_->find_nearest(query);
    vector<int> matches = vector<int>();

    // Search through each of the nearest elements
    for (vector<int>::iterator elem = candidates.begin(); 
                               elem < candidates.end(); elem++) {
      if (contains_point(*elem, query)) {
        matches.push_back(*elem); // Keep track of the elements 
                                  // which contain it
      }
    }

    // Return the index of the parent element which does contain the point
    if (matches.size() == 1) {  // x is conclusively in a single element
      return matches[0];
    } else if (matches.size() == 0) {  // not so much
      throw ATC_Error("FE_3DMesh::map_to_element could not find an element");
    } else {  // definitely not so much
      throw ATC_Error("FE_3DMesh::map_to_element found multiple elements");
    }
  }

  bool FE_3DMesh::contains_point(const int eltID,
                                 const DENS_VEC &x) const
  {
    DENS_MAT eltCoords;
    element_coordinates(eltID, eltCoords);
    return feElement_->contains_point(eltCoords, x);
  }

  //-----------------------------------------------------------------------
  void FE_3DMesh::distribute_mesh_data() 
  {
    myNElts_ = myElts_.size();

    //create elemToMyElemMap_
    elemToMyElemMap_.clear();
    for (int myID = 0; myID < myNElts_; ++myID) {
      int baseID = myElts_[myID];
      elemToMyElemMap_[baseID] = myID;
    }

    // create myConnectivity_, myConnectivityUnique_
    int numEltNodes = feElement_->num_elt_nodes();
    myConnectivity_.reset(numEltNodes, myNElts_);
    myConnectivityUnique_.reset(numEltNodes, myNElts_);
    for (int elem = 0; elem < myNElts_; ++elem) {
      for (int node = 0; node < numEltNodes; ++node) {
        myConnectivity_(node,elem) = connectivity_(node,map_myElem_to_elem(elem));
        myConnectivityUnique_(node,elem) = connectivityUnique_(node,map_myElem_to_elem(elem));
      }
    }
  }

  // -------------------------------------------------------------
  // -------------------------------------------------------------
  //   class FE_Rectangular3DMesh
  // -------------------------------------------------------------
  // -------------------------------------------------------------

  FE_Rectangular3DMesh::FE_Rectangular3DMesh(
                                     const Array<double> & hx,
                                     const Array<double> & hy,
                                     const Array<double> & hz,
                                     const double xmin, const double xmax,
                                     const double ymin, const double ymax,
                                     const double zmin, const double zmax,
                                     const Array<bool> periodicity,
                                     const double xscale,
                                     const double yscale,
                                     const double zscale)
    : hx_(hx), hy_(hy), hz_(hz)
  {
    tree_ = NULL;
    hasPlanarFaces_ = true;
    xscale_ = xscale;
    yscale_ = yscale;
    zscale_ = zscale;

    borders_[0][0] = xmin;
    borders_[0][1] = ymin;
    borders_[0][2] = zmin;
    borders_[1][0] = xmax;
    borders_[1][1] = ymax;
    borders_[1][2] = zmax;
    L_[0] = xmax-xmin; 
    L_[1] = ymax-ymin;
    L_[2] = zmax-zmin;
    n_[0] = hx_.size(); 
    n_[1] = hy_.size();
    n_[2] = hz_.size();
    // Compute region size and element size
    double Lx = 0;
    for (int i = 0; i < n_[0]; ++i) { Lx += hx_(i); }
    double Ly = 0;
    for (int i = 0; i < n_[1]; ++i) { Ly += hy_(i); }
    double Lz = 0;
    for (int i = 0; i < n_[2]; ++i) { Lz += hz_(i); }
    // rescale to fit box
    double ax = L_[0]/Lx;
    for (int i = 0; i < n_[0]; ++i) { hx_(i) *= ax; }
    double ay = L_[1]/Ly;
    for (int i = 0; i < n_[1]; ++i) { hy_(i) *= ay; }
    double az = L_[2]/Lz;
    for (int i = 0; i < n_[2]; ++i) { hz_(i) *= az; }

    // fill node locations
    nSD_ = 3;
    x_.reserve(nSD_);
    for (int i = 0; i < nSD_; ++i) {x_.push_back(Array<double>(n_[i]+1)); }
    Array<double> & xI = x_[0]; 
    xI(0) = xmin;
    for (int i = 0; i < n_[0]; ++i) { xI(i+1) = xI(i)+hx_(i); }
    Array<double> & yI = x_[1]; 
    yI(0) = ymin;
    for (int i = 0; i < n_[1]; ++i) { yI(i+1) = yI(i)+hy_(i); }
    Array<double> & zI = x_[2]; 
    zI(0) = zmin;
    for (int i = 0; i < n_[2]; ++i) { zI(i+1) = zI(i)+hz_(i); }

    // Member data setup
    nElts_ = n_[0] * n_[1] * n_[2];
    nNodes_ = (n_[0]+1) * (n_[1]+1) * (n_[2]+1);

    periodicity_ = Array<bool>(periodicity);

    feElement_ = new FE_ElementRect();
    nodalCoords_.reset(3, nNodes_);
    connectivity_.reset(feElement_->num_elt_nodes(), nElts_);

    // Fill nodal coordinates
    double x[3] = {xmin,ymin,zmin};
    int inode = 0;
    for (int k = 0; k <= n_[2]; ++k) {
      for (int j = 0; j <= n_[1]; ++j) {
        for (int i = 0; i <= n_[0]; ++i) {
          for (int m = 0; m < 3; ++m) {
            nodalCoords_(m,inode) = x[m];
          }
          ++inode;
          if (i < n_[0]) x[0] += hx_(i);
        }
        if (j < n_[1]) x[1] += hy_(j);
        x[0] = xmin;
      }
      if (k < n_[2]) x[2] += hz_(k);
      x[1] = ymin;
    }

    // Compute element connectivities
    int ielt = 0;
    int noffx = 1;
    int noffy = n_[0] + 1;
    int noffz = (n_[0]+1) * (n_[1]+1);
    for (int k = 0; k < n_[2]; ++k) {
      for (int j = 0; j < n_[1]; ++j) {
        for (int i = 0; i < n_[0]; ++i) {
          int i1 = i + j*noffy + k*noffz;
          connectivity_(0,ielt) = i1;
          connectivity_(1,ielt) = i1 + noffx;
          connectivity_(2,ielt) = i1 + noffx + noffy;
          connectivity_(3,ielt) = i1 + noffy;
          connectivity_(4,ielt) = i1 + noffz;
          connectivity_(5,ielt) = i1 + noffx + noffz;
          connectivity_(6,ielt) = i1 + noffx + noffy + noffz;
          connectivity_(7,ielt) = i1 + noffy + noffz;
          ++ielt;
        }
      }
    }
    setup_periodicity();
  }

  // -------------------------------------------------------------
  //  partition_mesh 
  // -------------------------------------------------------------
  
  void FE_Rectangular3DMesh::partition_mesh()
  {
    if (lammpsPartition_) {
      lammps_partition_mesh();
    }

    if (partitioned_) return;


    // Currently this code has been rather naively copied from 
    // FE_Uniform3DMesh::partition_mesh()

    // Determine dimensions of mesh in order to partition according to largest dimension.
    double xmin = borders_[0][0];
    double xmax = borders_[1][0];
    double ymin = borders_[0][1];
    double ymax = borders_[1][1];
    double zmin = borders_[0][2];
    double zmax = borders_[1][2];

    int numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    int processorRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &processorRank);
      
    // Spatially partition along the largest dimension.
    procs_.clear();
    if (max(max(L_[0], L_[1]), L_[2]) == L_[0]) {
      partitionAxis_ = 0;
      for (int i = 0; i < numProcs; ++i) {
        procs_.push_back(xmin + (L_[0]*i)/numProcs);
      }
      procs_.push_back(xmax);
    }
    else if (max(max(L_[0], L_[1]), L_[2]) == L_[1]) {
      partitionAxis_ = 1;
      for (int i = 0; i < numProcs; ++i) {
        procs_.push_back(ymin + (L_[1]*i)/numProcs);
      }
      procs_.push_back(ymax);
    }
    else {
      partitionAxis_ = 2;
      for (int i = 0; i < numProcs; ++i) {
        procs_.push_back(zmin + (L_[2]*i)/numProcs);
      }
      procs_.push_back(zmax);
    }

    // Distribute each node to the correct processor
    myNodes_.clear();
    for (int i = 0; i < nNodes_; ++i) {
      // Allocate node to this processor if it lies between processor's left and right boundaries.
      if ( dbl_geq(nodalCoords_(partitionAxis_, i), procs_[processorRank]) &&
          !dbl_geq(nodalCoords_(partitionAxis_, i), procs_[processorRank + 1])) {
        myNodes_.push_back(i);
      }
      // Also allocate nodes on the right boundary to the last processor.
      if ((processorRank == numProcs - 1) &&
           dbl_geq(nodalCoords_(partitionAxis_, i), procs_[processorRank + 1])) {
        myNodes_.push_back(i);
      }
    }
   
    // Distribute each element to the correct processor - assign it to the processor
    // which owns its node of lowest index. (this partitioning scheme is unambiguous) 
    myElts_.clear();
    for (int i = 0; i < nElts_; ++i) {
      int min = INT_MAX;
      for (int j = 0; j < connectivity_.nRows(); j++) {
        if (connectivity_(j, i) < min)
          min = connectivity_(j, i);
      }
      if (find(myNodes_.begin(), myNodes_.end(), min) != myNodes_.end()) {
        myElts_.push_back(i);
      }
    }

    /* Commented out because ghost nodes are never used and dx_ is not a member of
       FE_Rectangular3DMesh.

    // Compute the facesets that describes the left and right boundaries 
    // in order to determine ghost nodes.
    int leftMult = 0;
    while ((leftMult+1)*dx_[partitionAxis_] < procs_[processorRank]) {
      ++leftMult;
    }
    int rightMult = 0;
    while ((rightMult)*dx_[partitionAxis_] < procs_[processorRank+1]) {
      ++rightMult;
    }
    // Compute our ghost nodes - nodes that we need that belong to adjacent processors,
    // and our shared nodes - our nodes that are ghosted on the adjacent processors.
    for (int i = 0; i < nNodes_; ++i) {
      if (nodalCoords_(partitionAxis_, i) == leftMult*dx_[partitionAxis_])
        ghostNodesL_.push_back(i);
      else if (nodalCoords_(partitionAxis_, i) == rightMult*dx_[partitionAxis_])
        ghostNodesR_.push_back(i);
      else if (nodalCoords_(partitionAxis_, i) == (leftMult+1)*dx_[partitionAxis_])
        shareNodesL_.push_back(i);
      else if (nodalCoords_(partitionAxis_, i) == (rightMult-1)*dx_[partitionAxis_])
        shareNodesR_.push_back(i);
    }*/
    if (decomposition_) distribute_mesh_data();
    partitioned_ = true;
  }

  void FE_Rectangular3DMesh::departition_mesh()
  {
    if (!partitioned_) return;
    partitioned_ = false;
  }

  // -------------------------------------------------------------
  //   setup_periodicity
  // -------------------------------------------------------------
  void FE_Rectangular3DMesh::setup_periodicity()
  {
    nNodesUnique_ = 1;
    for (int i = 0; i < 3; i++) {
      nNodesUnique_ *= (n_[i] + 1 - periodicity_(i));
    } 
    
    // form maximal nodeset
    for (int i = 0; i < nNodesUnique_; i++) {
      nodeSetAll_.insert(i); 
    }
    
    // Create global-to-unique map: globalToUniqueMap_(ig) = iu
    globalToUniqueMap_.reset(nNodes_);
    uniqueToGlobalMap_.reset(nNodesUnique_);
    for (int k = 0; k <= n_[2]; ++k) {
      int kper = (k == n_[2] && periodicity_(2)) ? 0 : k;
      for (int j = 0; j <= n_[1]; ++j) {
        int jper = (j == n_[1] && periodicity_(1)) ? 0 : j;
        for (int i = 0; i <= n_[0]; ++i) {
          int iper = (i == n_[0] && periodicity_(0)) ? 0 : i;
          int id = i + j*(n_[0]+1) + k*(n_[0]+1)*(n_[1]+1);
          int uid = iper + jper*(n_[0]+1-periodicity_(0)) 
            + kper*(n_[0]+1-periodicity_(0))*(n_[1]+1-periodicity_(1));
          globalToUniqueMap_(id) = uid;
          uniqueToGlobalMap_(uid) = id;
        }
      }
    }
    set_unique_connectivity();

    // form maximal elementset
    for (int i = 0; i < nElts_; i++) {
      elementSetAll_.insert(i); 
    }
  }

  int FE_Rectangular3DMesh::map_to_element(const DENS_VEC &x) const
  {
    int ix[3]; // ix[i] is the element in the ith direction
    for (int i = 0; i < 3; i++) {
      // map to box
      double y = x(i);
      if (periodicity_(i)) {
         double diff = y-borders_[0][i];
         int shift = int(diff/L_[i]);
         if (diff < 0.) shift--;
         y -= shift*L_[i];
      } 
      // project into element
      ix[i] = x_[i].index(y);
      if (fabs(y-borders_[0][i]) < tol) { ix[i] = 0; } // on the lower boundary
      if (ix[i] < 0 || ix[i] >= n_[i]) {
        string msg = "FE_Rectangular3DMesh:: point maps outside of mesh, coordinate "
         + index_to_string(i) + "=" + to_string(x(i)) + " image=" + to_string(y)
         + " not in " + to_string(borders_[0][i]) + ":" + to_string(borders_[1][i]);
        throw ATC_Error(msg);
      }
    }
    int elt = ix[2]*(n_[0]*n_[1]) + ix[1]*n_[0] + ix[0];
    return elt;
  }
  // -------------------------------------------------------------
  // -------------------------------------------------------------
  //   class FE_Uniform3DMesh
  // -------------------------------------------------------------
  // -------------------------------------------------------------

  FE_Uniform3DMesh::FE_Uniform3DMesh(const int nx,
                                     const int ny,
                                     const int nz,
                                     const double xmin, const double xmax,
                                     const double ymin, const double ymax,
                                     const double zmin, const double zmax,
                                     const Array<bool> periodicity,
                                     const double xscale,
                                     const double yscale,
                                     const double zscale)
  {
    hasPlanarFaces_ = true;
    tree_ = NULL;
    xscale_ = xscale;
    yscale_ = yscale;
    zscale_ = zscale;
    n_[0] = nx;
    n_[1] = ny;
    n_[2] = nz;

    borders_[0][0] = xmin;
    borders_[1][0] = xmax;
    borders_[0][1] = ymin;
    borders_[1][1] = ymax;
    borders_[0][2] = zmin;
    borders_[1][2] = zmax;

    periodicity_ = Array<bool>(periodicity);

    // Compute region size and element size
    for (int i = 0; i < 3; i++) {
      L_[i] = borders_[1][i] - borders_[0][i];
      dx_[i] = L_[i]/n_[i];
    }

    // Member data setup
    nSD_ = 3;
    nElts_ = n_[0] * n_[1] * n_[2];
    nNodes_ = (n_[0]+1) * (n_[1]+1) * (n_[2]+1);

    feElement_ = new FE_ElementRect();

    nodalCoords_.reset(3, nNodes_);
    connectivity_.reset(feElement_->num_elt_nodes(), nElts_);

    // Fill nodal coordinates
    double ix[3];
    int inode = 0;
    for (int k = 0; k <= n_[2]; ++k) {
      ix[2] = borders_[0][2] + k*dx_[2];
      for (int j = 0; j <= n_[1]; ++j) {
        ix[1] = borders_[0][1] + j*dx_[1];
        for (int i = 0; i <= n_[0]; ++i) {
          ix[0] = borders_[0][0] + i*dx_[0];
          for (int m = 0; m < 3; ++m) {
            nodalCoords_(m,inode) = ix[m];
          }
          ++inode;
        }
      }
    }

    // Compute element connectivities
    int ielt = 0;
    int noffx = 1;
    int noffy = n_[0] + 1;
    int noffz = (n_[0]+1) * (n_[1]+1);
    for (int k = 0; k < n_[2]; ++k) {
      for (int j = 0; j < n_[1]; ++j) {
        for (int i = 0; i < n_[0]; ++i) {
          int i1 = i + j*noffy + k*noffz;
          connectivity_(0,ielt) = i1;
          connectivity_(1,ielt) = i1 + noffx;
          connectivity_(2,ielt) = i1 + noffx + noffy;
          connectivity_(3,ielt) = i1 + noffy;
          connectivity_(4,ielt) = i1 + noffz;
          connectivity_(5,ielt) = i1 + noffx + noffz;
          connectivity_(6,ielt) = i1 + noffx + noffy + noffz;
          connectivity_(7,ielt) = i1 + noffy + noffz;
          ++ielt;
        }
      }
    }

    setup_periodicity();
  }

  // -------------------------------------------------------------
  //  destructor
  // -------------------------------------------------------------
  FE_Uniform3DMesh::~FE_Uniform3DMesh() 
  {
    // Clean up is currently unimplemented
  }

  // -------------------------------------------------------------
  //  partition_mesh 
  // -------------------------------------------------------------
  
  void FE_Uniform3DMesh::partition_mesh()
  {
    if (lammpsPartition_) {
      lammps_partition_mesh();
    }

    if (partitioned_) return;


    // Determine dimensions of mesh in order to partition according to largest dimension.
    double xmin = borders_[0][0];
    double xmax = borders_[1][0];
    double ymin = borders_[0][1];
    double ymax = borders_[1][1];
    double zmin = borders_[0][2];
    double zmax = borders_[1][2];

    int numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    int processorRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &processorRank);
      
    // Spatially partition along the largest dimension.
    procs_.clear();
    if (max(max(L_[0], L_[1]), L_[2]) == L_[0]) {
      partitionAxis_ = 0;
      for (int i = 0; i < numProcs; ++i) {
        procs_.push_back(xmin + (L_[0]*i)/numProcs);
      }
      procs_.push_back(xmax);
    }
    else if (max(max(L_[0], L_[1]), L_[2]) == L_[1]) {
      partitionAxis_ = 1;
      for (int i = 0; i < numProcs; ++i) {
        procs_.push_back(ymin + (L_[1]*i)/numProcs);
      }
      procs_.push_back(ymax);
    }
    else {
      partitionAxis_ = 2;
      for (int i = 0; i < numProcs; ++i) {
        procs_.push_back(zmin + (L_[2]*i)/numProcs);
      }
      procs_.push_back(zmax);
    }

    // Distribute each node to the correct processor
    myNodes_.clear();
    for (int i = 0; i < nNodes_; ++i) {
      // Allocate node to this processor if it lies between processor's left and right boundaries.
      if ( dbl_geq(nodalCoords_(partitionAxis_, i), procs_[processorRank]) &&
          !dbl_geq(nodalCoords_(partitionAxis_, i), procs_[processorRank + 1])) {
        myNodes_.push_back(i);
      }
      // Also allocate nodes on the right boundary to the last processor.
      if ((processorRank == numProcs - 1) &&
           dbl_geq(nodalCoords_(partitionAxis_, i), procs_[processorRank + 1])) {
        myNodes_.push_back(i);
      }
    }
   
    // Distribute each element to the correct processor - assign it to the processor
    // which owns its node of lowest index. (this partitioning scheme is unambiguous) 
    myElts_.clear();
    for (int i = 0; i < nElts_; ++i) {
      int min = INT_MAX;
      for (int j = 0; j < connectivity_.nRows(); j++) {
        if (connectivity_(j, i) < min)
          min = connectivity_(j, i);
      }
      if (find(myNodes_.begin(), myNodes_.end(), min) != myNodes_.end()) {
        myElts_.push_back(i);
      }
    }

    // Compute the facesets that describes the left and right boundaries 
    // in order to determine ghost nodes.
    int leftMult = 0;
    while ((leftMult+1)*dx_[partitionAxis_] < procs_[processorRank]) {
      ++leftMult;
    }
    int rightMult = 0;
    while ((rightMult)*dx_[partitionAxis_] < procs_[processorRank+1]) {
      ++rightMult;
    }
    // Compute our ghost nodes - nodes that we need that belong to adjacent processors,
    // and our shared nodes - our nodes that are ghosted on the adjacent processors.
    for (int i = 0; i < nNodes_; ++i) {
      if (nodalCoords_(partitionAxis_, i) == leftMult*dx_[partitionAxis_])
        ghostNodesL_.push_back(i);
      else if (nodalCoords_(partitionAxis_, i) == rightMult*dx_[partitionAxis_])
        ghostNodesR_.push_back(i);
      else if (nodalCoords_(partitionAxis_, i) == (leftMult+1)*dx_[partitionAxis_])
        shareNodesL_.push_back(i);
      else if (nodalCoords_(partitionAxis_, i) == (rightMult-1)*dx_[partitionAxis_])
        shareNodesR_.push_back(i);
    }
    if (decomposition_) distribute_mesh_data();
    partitioned_ = true;
  }

  void FE_Uniform3DMesh::departition_mesh()
  {
    if (!partitioned_) return;
    partitioned_ = false;
  }

  // -------------------------------------------------------------
  //   map_to_element
  // -------------------------------------------------------------
  int FE_Uniform3DMesh::map_to_element(const DENS_VEC &x) const
  {
    // countx[i] is the number of the element, where 1 is the 
    // element adjacent to the lower border, in the ith direction
    int countx[3];
    for (int i = 0; i < 3; ++i) {
      // handle points on upper boundary; not sure why this is 
      // hard-coded in, though...
      if (fabs(x(i)-borders_[1][i]) < tol) {
        countx[i] = n_[i] - 1;
      } else {
        // find the x, y, and z bins for the point in the mesh
        countx[i] = (int)floor((x(i)-borders_[0][i])/dx_[i]);
      }
      
      // handle points out-of-range [0:nx-1] w/ periodicity
      if (countx[i] < 0 || countx[i] >= n_[i]) {
        if (periodicity_(i)) {
          countx[i] = countx[i] % n_[i];
          // handles c++ ambiguous mod problems
          if (countx[i] < 0) countx[i] += n_[i]; 
        } else {
          string msg = " point maps outside " 
                       "of mesh, coordinate "                     + 
                       index_to_string(i) + " " + to_string(x(i)) + 
                       " not in " + to_string(borders_[0][i])     + 
                       " : "      + to_string(borders_[1][i]);
          throw ATC_Error(FILELINE,msg);
        }
      }
    }

    int elt = countx[2]*(n_[0]*n_[1]) + 
              countx[1]*n_[0]          + 
              countx[0];
    return elt;
  }


} // namespace ATC

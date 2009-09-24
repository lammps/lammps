// ATC header files
#include "FE_Element.h"
#include "FE_Mesh.h"
#include "LammpsInterface.h"
#include "ATC_Error.h"
#include "OutputManager.h"
#include "StringManip.h"

// Other headers
#include <iostream>

using namespace std;
using namespace ATC_STRING;

namespace ATC {

  const static double tol = 1.0e-10;

  // -------------------------------------------------------------
  // -------------------------------------------------------------
  //   class FE_Mesh
  // -------------------------------------------------------------
  // -------------------------------------------------------------

  FE_Mesh::FE_Mesh():
    nNodesUnique_(0), nNodes_(0)
  {
    feElement_ = NULL;
  }

  FE_Mesh::~FE_Mesh()
  {
    if (feElement_)         delete feElement_;
  }

  // -------------------------------------------------------------
  //  modify
  // -------------------------------------------------------------
  bool FE_Mesh::modify(int narg, char **arg)
  {
    bool match = false;

    if (strcmp(arg[0],"mesh")==0) 
    {
     /*! \page man_mesh_faceset fix_modify AtC mesh create_faceset
        \section syntax
         fix_modify AtC create_faceset <id> 
         <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> <in|out> [units]
         - <id> = id to assign to the collection of FE faces
         - <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of
         the bounding box that is coincident with the desired FE faces
         - <in|out> = "in" gives inner faces to the box, 
                      "out" gives the outer faces to the box
         - units = option to specify real as opposed to lattice units
        \section examples
         <TT> fix_modify AtC mesh create_faceset obndy -4.0 4.0 -12 12 -12 12 out </TT>
        \section description
          Command to assign an id to a set of FE faces to be used subsequently
          in defining flux boundary conditions.
        \section restrictions
        Only viable for rectangular grids. Also "INF" is not currrently handled.
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
          string n(arg[argIdx++]);
          int idir, isgn;
          string_to_index(n, idir, isgn);
          double x =  atof(arg[argIdx++]);
          if (narg > argIdx && (strcmp(arg[argIdx++],"units") == 0)) 
          {} 
          else 
          {
            if      (idir == 0) { x *= xscale_; }
            else if (idir == 1) { x *= yscale_; }
            else if (idir == 2) { x *= zscale_; }
          }
          create_faceset(tag, x, idir, isgn);
          match = true;
        }
        // bounding_box
        else 
        {
          if (strcmp(arg[argIdx],"box")==0) argIdx++;
          double xmin = atof(arg[argIdx++]);
          double xmax = atof(arg[argIdx++]);
          double ymin = atof(arg[argIdx++]);
          double ymax = atof(arg[argIdx++]);
          double zmin = atof(arg[argIdx++]);
          double zmax = atof(arg[argIdx++]);
          bool outward = true;
          if (narg > argIdx && (strcmp(arg[argIdx++],"in") == 0)) 
            outward = false;

          if (narg > argIdx && (strcmp(arg[argIdx++],"units") == 0)) 
          {}
          else 
          { // scale from lattice units to physical units
            xmin *= xscale_;
            xmax *= xscale_;
            ymin *= yscale_;
            ymax *= yscale_;
            zmin *= zscale_;
            zmax *= zscale_;
          }
          
          create_faceset(tag, xmin, xmax, ymin, ymax, zmin, zmax, outward);
          match = true;
        }
      }
   /*! \page man_mesh_nodeset fix_modify AtC mesh create_nodeset
      \section syntax
      fix_modify AtC create_nodeset <id>
      <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>
      - <id> = id to assign to the collection of FE nodes
      - <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> = coordinates of
      the bounding box that contains only the desired nodes
      \section examples
      <TT> fix_modify AtC mesh create_faceset left -4.1 3.9 -100 100 -100 100 </TT>
      \section description
      Command to assign an id to a set of FE nodes to be used subsequently
      in defining boundary conditions.
      \section restrictions
      Only viable for rectangular grids. Also "INF" is not currrently handled.
      \section related
      \section default
      Coordinates are assumed to be in lattice units.
    */
      else if (strcmp(arg[1],"create_nodeset")==0) {
        string tag  = arg[2];
        double xmin = atof(arg[3]);
        double xmax = atof(arg[4]);
        double ymin = atof(arg[5]);
        double ymax = atof(arg[6]);
        double zmin = atof(arg[7]);
        double zmax = atof(arg[8]);
        create_nodeset(tag, xmin, xmax, ymin, ymax, zmin, zmax);
        match = true;
      }
   /*! \page man_mesh_elemset fix_modify AtC mesh create_elementset
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
      Only viable for rectangular grids. Also "INF" is not currrently handled.
      \section related
      \section default
      Coordinates are assumed to be in lattice units.
    */
      else if (strcmp(arg[1],"create_elementset")==0) {
        string tag  = arg[2];
        double xmin = atof(arg[3]);
        double xmax = atof(arg[4]);
        double ymin = atof(arg[5]);
        double ymax = atof(arg[6]);
        double zmin = atof(arg[7]);
        double zmax = atof(arg[8]);
        create_elementset(tag, xmin, xmax, ymin, ymax, zmin, zmax);
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
  //   element_coordinates
  // -------------------------------------------------------------
  void FE_Mesh::element_coordinates(const int eltID,
                                    DENS_MAT & xCoords) const
  {
    const int nne = get_nNodesPerElement();
    xCoords.reset(nSD_, nne, false);
    for (int inode=0; inode<nne; inode++) 
    {
      const int id = connectivity_(inode, eltID);
      for (int isd=0; isd<nSD_; isd++) 
        xCoords(isd,inode) = nodalCoords_(isd,id);
    }
  
  }
  // -------------------------------------------------------------
  //   face_coordinates
  // -------------------------------------------------------------
  void FE_Mesh::face_coordinates(const PAIR face, DENS_MAT & xCoords) const
  {
    const int eltID=face.first, faceID=face.second;
    const int nnf = get_nNodesPerFace();
    const Array2D <int> & local_conn = feElement_->local_face_conn();

    xCoords.reset(nSD_, nnf, false);

    for (int inode=0; inode < nnf; inode++) 
    {
      int id = connectivity_(local_conn(faceID,inode), eltID);
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
  const set<int> & FE_Mesh::get_nodeset(const string & name) const
  {
    NODE_SET_MAP::const_iterator iter = nodeSetMap_.find(name);
    if (name == "all") return nodeSetAll_;
    else if (iter == nodeSetMap_.end()) 
      throw ATC_Error(0, "No nodeset with name " + name + " found.");  
    else return iter->second;
  }

  // -------------------------------------------------------------
  //   get_elementset
  // -------------------------------------------------------------
  const set<int> & FE_Mesh::get_elementset(const string & name) const
  {
    NODE_SET_MAP::const_iterator iter = elementSetMap_.find(name);
    if (name == "all") return elementSetAll_;
    else if (iter == elementSetMap_.end()) 
      throw ATC_Error(0, "No elementset with name " + name + " found.");  
    else return iter->second;
  }

  // -------------------------------------------------------------
  //   nodeset_to_minimal_elementset
  // -------------------------------------------------------------
  void FE_Mesh::nodeset_to_minimal_elementset
    (const string & name, set<int> & elemSet) const
  {
    if (name == "all")
      for (int ielem = 0; ielem < nElts_; ielem++) 
        elemSet.insert(ielem);

    else 
    {
      NODE_SET_MAP::const_iterator iter = nodeSetMap_.find(name);
      if (iter == nodeSetMap_.end()) 
        throw ATC_Error(0, "No nodeset with name " + name + " found.");

      int npe = get_nNodesPerElement();
      const set<int> &nodeSet = iter->second;
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
    }
  }

  // -------------------------------------------------------------
  //   nodeset_to_maximal_elementset
  // -------------------------------------------------------------
  void FE_Mesh::nodeset_to_maximal_elementset(const string &name, set<int> &elemSet) const
  {
    if (name == "all")
      for (int ielem = 0; ielem < nElts_; ielem++)
        elemSet.insert(ielem);

    else 
    {
      NODE_SET_MAP::const_iterator iter = nodeSetMap_.find(name);
      if (iter == nodeSetMap_.end()) 
        throw ATC_Error(0, "No nodeset with name " + name + " found.");

      int npe = get_nNodesPerElement();
      const set<int> & nodeSet = iter->second;
      for (int ielem = 0; ielem < nElts_; ielem++) 
      {
        int inode = 0;
        bool in = false;
        while (!in && inode < npe) 
        {
          int node = connectivityUnique_(inode, ielem);
          set<int>::const_iterator iter = nodeSet.find(node);
          if (iter != nodeSet.end()) { in = true; }
          inode++;
        }
        if (in) elemSet.insert(ielem);
      }
    }
  }

  // -------------------------------------------------------------
  //   elementset_to_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::elementset_to_nodeset
    (const set<int> & elemSet, set<int> & nodeSet) const
  {
    int npe = get_nNodesPerElement();
    set<int>::const_iterator itr;
    for (itr = elemSet.begin(); itr != elemSet.end(); itr++) 
    {
      int ielem = *itr;
      for (int inode=0; inode < npe; inode++) 
      {
        int node = connectivity_(inode, ielem);
        nodeSet.insert(node);
      }
    }
  }

  // -------------------------------------------------------------
  //   elementset_to_nodeset
  // -------------------------------------------------------------
  void FE_Mesh::elementset_to_nodeset
    (const string & name, set<int> & nodeSet) const
  {
    if (name == "all")
      for (int ielem = 0; ielem < nElts_; ielem++) 
        nodeSet.insert(ielem);

    else 
    {
      ELEMENT_SET_MAP::const_iterator iter = elementSetMap_.find(name);
      if (iter == elementSetMap_.end()) 
        throw ATC_Error(0, "No elementset with name " + name + " found.");

      int npe = get_nNodesPerElement();
      const set<int> &elemSet = iter->second;
      set<int>::const_iterator itr;
      for (itr = elemSet.begin(); itr != elemSet.end(); itr++) 
      {
        int ielem = *itr;
        for (int inode=0; inode < npe; inode++) 
        {
          int node = connectivityUnique_(inode, ielem);
          nodeSet.insert(node);
          inode++;
        }
      }
    }
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
      int npe = get_nNodesPerElement();
      set<int>::const_iterator itr;
      for (itr = compElemSet.begin(); itr != compElemSet.end(); itr++) 
      {
        int ielem = *itr;
        for (int inode=0; inode < npe; inode++) 
        {
          int node = connectivityUnique_(inode, ielem);
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
        throw ATC_Error(0, "No elementset with name " + name + " found.");

      int npe = get_nNodesPerElement();
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
    int npe = get_nNodesPerElement();
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
        throw ATC_Error(0, "No faceset with name " + name + " found.");
      const set<PAIR> & faceSet = faceset->second;
      set<PAIR>::const_iterator iter;
      Array <int> conn;
      for (iter = faceSet.begin(); iter != faceSet.end(); iter++) 
      {
        pair<int,int> face = *iter;
        face_connectivity_unique(face,conn);
        for (int i = 0; i < conn.get_length() ; ++i) {
          int inode = conn(i);
          nodeSet.insert(inode);
        }
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
        throw ATC_Error(0, "No faceset with name " + name + " found.");
      const set<PAIR> & faceSet = faceset->second;
      set<PAIR>::const_iterator iter;
      Array <int> conn;
      for (iter = faceSet.begin(); iter != faceSet.end(); iter++) 
      {
        pair<int,int> face = *iter;
        face_connectivity(face,conn);
        for (int i = 0; i < conn.get_length() ; ++i) {
          int inode = conn(i);
          nodeSet.insert(inode);
        }
      }
    }
  }


  // -------------------------------------------------------------
  //   get_faceset
  // -------------------------------------------------------------
  const set<PAIR> &FE_Mesh::get_faceset(const string & name) const
  {
    FACE_SET_MAP::const_iterator iter = faceSetMap_.find(name);
    if (iter == faceSetMap_.end())
    {
      throw ATC_Error(0, "No faceset with name " + name + " found.");
    }
    return iter->second;
  }

  // -------------------------------------------------------------
  //   create_nodeset
  // -------------------------------------------------------------
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
      throw ATC_Error(0, message);
    }

    xmin *= xscale_;
    xmax *= xscale_;
    ymin *= yscale_;
    ymax *= yscale_;
    zmin *= zscale_;
    zmax *= zscale_;

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
      throw ATC_Error(0, message);
    }

    nodeSetMap_[name] = nodeSet;

    if (ATC::LammpsInterface::instance()->comm_rank() == 0) { 
      cout << " ATC:: created nodeset " << name 
           << " with " << nodeSet.size() << " nodes\n";
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
      throw ATC_Error(0, "A faceset with name " + name + " is already defined.");

    set<PAIR> faceSet;
    // Loop over face and add their unique id's to the set if they concide
    // with region
    const int nf = get_nFacesPerElement();
    const int npf = get_nNodesPerFace();
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
          // NOTE this is based on structured grid
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
    if (faceSet.empty()) throw ATC_Error(0, "faceset "+name+" is empty.");

    faceSetMap_[name] = faceSet;
    if (ATC::LammpsInterface::instance()->comm_rank() == 0) { 
      cout << " ATC:: created faceset " << name 
           << " with " << faceSet.size() << " faces\n";
    }
  }

  void FE_Mesh::create_faceset(const string & name,
                               double xRef, 
                               int nIdx, int nSgn)
  {
    // Make sure we don't already have a nodeset with this name
    FACE_SET_MAP::iterator iter = faceSetMap_.find(name);
    // NOTE change this to add to a faceset
    if (iter != faceSetMap_.end()) 
      throw ATC_Error(0, "A faceset with name "+name+" is already defined.");
    if (periodicFlag_[nIdx]==1) 
      throw ATC_Error(0,"Faceset definition via plane in periodic direction.");

    set<PAIR> faceSet;
    // Loop over faces i& add unique id's to the set if concide w/ plane
    int nf = get_nFacesPerElement();
    int npf = get_nNodesPerFace();
    const Array2D<int> & face_conn = local_face_connectivity();
    for (int ielem = 0; ielem < nElts_; ielem++) 
    {
      for (int iface = 0; iface < nf; iface++) 
      {
        bool in = true;
        bool neg = false, pos = false;
        // all nodes must be on the plane
        for (int inode = 0; inode < npf; inode++) {
          double x 
             = nodalCoords_(nIdx,connectivity_(face_conn(iface,inode),ielem));
          if ( fabs(x-xRef) > tol){ in = false;}
        }
        // check correct orientation
        bool add = false;
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
      throw ATC_Error(0, "faceset "+name+" is empty.");

    faceSetMap_[name] = faceSet;
    if (ATC::LammpsInterface::instance()->comm_rank() == 0) { 
      cout << " ATC:: created faceset " << name 
           << " with " << faceSet.size() << " faces\n";
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
      throw ATC_Error(0, message);
    }

    xmin *= xscale_;
    xmax *= xscale_;
    ymin *= yscale_;
    ymax *= yscale_;
    zmin *= zscale_;
    zmax *= zscale_;

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
      throw ATC_Error(0, message);
    }

    // create a minimal element set from all the nodes included in the region
    set<int> elemSet;
    int npe = get_nNodesPerElement();
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
      throw ATC_Error(0, message);
    }
    elementSetMap_[name] = elemSet;

    if (ATC::LammpsInterface::instance()->comm_rank() == 0) { 
      cout << " ATC:: created elementset " << name 
           << " with " << elemSet.size() << " elements\n";
    }
  }
  // -------------------------------------------------------------
  //   get_nIPsPerElement()
  // -------------------------------------------------------------
  int FE_Mesh::get_nIPsPerElement() const
  {
    return feElement_->num_ips();
  }
  // -------------------------------------------------------------
  //   get_nNodesPerElement()
  // -------------------------------------------------------------
  int FE_Mesh::get_nNodesPerElement() const
  {
    return feElement_->num_elt_nodes();
  }
  // -------------------------------------------------------------
  //   get_nFacesPerElement()
  // -------------------------------------------------------------
  int FE_Mesh::get_nFacesPerElement() const
  {
    return feElement_->num_faces();
  }
  // -------------------------------------------------------------
  //   get_nIPsPerFace()
  // -------------------------------------------------------------
  int FE_Mesh::get_nIPsPerFace() const
  {
    return feElement_->num_face_ips();
  }
  // -------------------------------------------------------------
  //   get_nNodesPerFace()
  // -------------------------------------------------------------
  int FE_Mesh::get_nNodesPerFace() const
  {
    return feElement_->num_face_nodes();
  }
  // -------------------------------------------------------------
  //   mappings from element id to associated nodes
  // -------------------------------------------------------------
  void FE_Mesh::element_connectivity_global(const int eltID,
                                            Array<int> & nodes) const
  {
    const int npe = get_nNodesPerElement();
    nodes.reset(npe);

    // use connectivity arrays
    for (int inode = 0; inode < npe; inode++) 
      nodes(inode) = connectivity_(inode, eltID);
  }
  // -------------------------------------------------------------
  //   
  // -------------------------------------------------------------
  void FE_Mesh::element_connectivity_unique(const int eltID,
                                            Array<int> & nodes) const
  {
    const int npe = get_nNodesPerElement();
    nodes.reset(npe);

    // use connectivity arrays
    for (int inode = 0; inode < npe; inode++)
      nodes(inode) = connectivityUnique_(inode, eltID);
  }

  // -------------------------------------------------------------
  //   local_face_connectivity()
  // -------------------------------------------------------------
  const Array2D<int> &  FE_Mesh::local_face_connectivity() const
  {
    return feElement_->local_face_conn();
  }

  // -------------------------------------------------------------
  //   shape function evaluation
  // -------------------------------------------------------------
  void FE_Mesh::shape_functions(const VECTOR &x,
                                DENS_VEC &shp,
                                int & eltID,
                                Array<int> &node_list) const
  {
    // get element id and local coordinates from global coordinates
    DENS_VEC xi;
    eltID = map_to_element(x, xi);

    // call FE_Engine shape function routines
    feElement_->shape_function(eltID,xi,shp);

    // determine nodes which correspond to shp function indices
    element_connectivity_unique(eltID,node_list);

    return;
  }

  void FE_Mesh::shape_functions(const VECTOR &x,
                                DENS_VEC &shp,
                                int & eltID,
                                Array<int> &node_list,
                                const Array<bool> &periodicity) const
  {
    // get element id and local coordinates from global coordinates
    DENS_VEC xi;
    eltID = map_to_element(x, xi, periodicity);

    // call FE_Engine shape function routines
    feElement_->shape_function(eltID,xi,shp);

    // determine nodes which correspond to shp function indices
    element_connectivity_unique(eltID,node_list);

    return;
  }

  void FE_Mesh::shape_functions(const VECTOR &x, 
                                DENS_VEC &shp,
                                DENS_MAT &dshp,
				int & eltID,
                                Array<int>& node_list) const
  {
    // get element id and local coordinates from global coordinates
    DENS_VEC xi;
    eltID = map_to_element(x,xi);

    // call FE_Engine shape function routines
    feElement_->shape_function(eltID, xi, shp, dshp);
    
    // determine nodes which correspond to shp function indices
    element_connectivity_unique(eltID,node_list);
    return;
  }

  void FE_Mesh::shape_function(int eltID, DENS_MAT &N, DIAG_MAT &weights) const
  {
    vector<DENS_MAT> dN; // NOTE dummy/placeholder
    // pass through call
    feElement_->shape_function(eltID, N, dN, weights);
  }

  void FE_Mesh::shape_function(int eltID, DENS_MAT &N, vector<DENS_MAT> &dN,
                               DIAG_MAT &weights) const
  {
    // pass through call
    feElement_->shape_function(eltID,N,dN,weights);
  }

  void FE_Mesh::face_shape_function(const PAIR &face,
                                    DenseMatrix<double> &N,
                                    vector< DenseMatrix<double> > &dN,
                                    vector< DenseMatrix<double> > &Nn,
                                    DiagonalMatrix<double> &weights) const
  {
    // pass through call
    feElement_->face_shape_function(face,N,dN,Nn,weights);
  }

  void FE_Mesh::face_shape_function(const PAIR &face,
                                    DenseMatrix<double> &N,
                                    DenseMatrix<double>  &n,
                                    DiagonalMatrix<double> &weights) const
  {
    // pass through call
    feElement_->face_shape_function(face,N,n,weights);
  }

  void FE_Mesh::set_quadrature(int type) {feElement_->set_quadrature(type);}

  //-----------------------------------------------------------------------
  void FE_Mesh::output(string prefix) const
  {
    OutputManager meshSets(prefix,GNUPLOT);
    meshSets.write_geometry(nodalCoords_, & connectivity_);
    OUTPUT_LIST subsetData;
    int size = nNodesUnique_;
    // material
//  DENS_MAT material(nNodes_,1);
//  material = 1;
//  subsetData["material"] = &material;
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
        for(int i = 0; i < nodes.get_length(); ++i) {
          (elemsets[elemset])(nodes(i),0) = 1; 
        }
      }
      subsetData[elemset] = & elemsets[elemset];
    }
    // output
    const int * map = globalToUniqueMap_.get_data();
    if (subsetData.size() > 0 )
      meshSets.write_data(0.0,&subsetData,map);
    else
      if (LammpsInterface::instance()->comm_rank() == 0)
        cout << "ATC: Warning mesh output requested without any mesh entities, output suppressed";
  }

  // -------------------------------------------------------------
  // -------------------------------------------------------------
  //   class FE_Uniform3DMesh
  // -------------------------------------------------------------
  // -------------------------------------------------------------

  FE_Uniform3DMesh::FE_Uniform3DMesh(const int nx,
                                     const int ny,
                                     const int nz,
                                     double xmin, double xmax,
                                     double ymin, double ymax,
                                     double zmin, double zmax,
                                     double xscale,
                                     double yscale,
                                     double zscale,
                                     int xperiodic,
                                     int yperiodic,
                                     int zperiodic)
  {
    nx_[0] = nx;
    nx_[1] = ny;
    nx_[2] = nz;

    borders_[0][0] = xmin;
    borders_[1][0] = xmax;
    borders_[0][1] = ymin;
    borders_[1][1] = ymax;
    borders_[0][2] = zmin;
    borders_[1][2] = zmax;

    xscale_ = xscale;
    yscale_ = yscale;
    zscale_ = zscale;

    // Compute region size and element size
    for (int i = 0; i < 3; i++) {
      Lx_[i] = borders_[1][i] - borders_[0][i];
      dx_[i] = Lx_[i]/nx_[i];
    }

    // Member data setup
    nSD_ = 3;
    nElts_ = nx_[0] * nx_[1] * nx_[2];
    nNodes_ = (nx_[0]+1) * (nx_[1]+1) * (nx_[2]+1);

    periodicFlag_[0] = xperiodic;
    periodicFlag_[1] = yperiodic;
    periodicFlag_[2] = zperiodic;

    feElement_ = new FE_ElementHex(this);

    nodalCoords_.reset(3, nNodes_);
    connectivity_.reset(feElement_->num_elt_nodes(), nElts_);

    // Fill nodal coordinates
    double ix[3];
    int inode = 0;
    for (int k = 0; k <= nx_[2]; ++k) {
      ix[2] = borders_[0][2] + k*dx_[2];
      for (int j = 0; j <= nx_[1]; ++j) {
        ix[1] = borders_[0][1] + j*dx_[1];
        for (int i = 0; i <= nx_[0]; ++i) {
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
    int noffy = nx_[0] + 1;
    int noffz = (nx_[0]+1) * (nx_[1]+1);
    for (int k = 0; k < nx_[2]; ++k) {
      for (int j = 0; j < nx_[1]; ++j) {
        for (int i = 0; i < nx_[0]; ++i) {
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

    // Setup periodicity
    setup_periodicity();
  }

  FE_Uniform3DMesh::~FE_Uniform3DMesh() 
  {
    // Nothing to do here
  }


  // -------------------------------------------------------------
  //   setup_periodicity
  // -------------------------------------------------------------
  void FE_Uniform3DMesh::setup_periodicity()
  {
    // int to flag whether fem mesh is periodic in each direction
    int perfem[3];
    for (int i = 0; i < 3; i++) {
      perfem[i] = periodicFlag_[i]; //1; 
    }

    // determine number of unique nodes
    nNodesUnique_ = 1;
    for (int i = 0; i < 3; i++) {
      nNodesUnique_ *= (nx_[i] + 1 - perfem[i]);
    } 

    // form maximal nodeset
    for (int i = 0; i < nNodesUnique_; i++) {
      nodeSetAll_.insert(i); 
    }

    // Create global-to-unique map: globalToUniqueMap_(ig) = iu
    globalToUniqueMap_.reset(nNodes_);
    uniqueToGlobalMap_.reset(nNodesUnique_);
    for (int k = 0; k <= nx_[2]; ++k) {
      int kper = (k == nx_[2] && perfem[2]) ? 0 : k;
      for (int j = 0; j <= nx_[1]; ++j) {
        int jper = (j == nx_[1] && perfem[1]) ? 0 : j;
        for (int i = 0; i <= nx_[0]; ++i) {
          int iper = (i == nx_[0] && perfem[0]) ? 0 : i;
          int id = i + j*(nx_[0]+1) + k*(nx_[0]+1)*(nx_[1]+1);
          int uid = iper + jper*(nx_[0]+1-perfem[0]) 
            + kper*(nx_[0]+1-perfem[0])*(nx_[1]+1-perfem[1]);
          globalToUniqueMap_(id) = uid;
          uniqueToGlobalMap_(uid) = id;
        }
      }
    }  

    // Map connectivity list
    int numEltNodes = feElement_->num_elt_nodes();
    connectivityUnique_.reset(numEltNodes, nElts_);
    for (int inode = 0; inode < numEltNodes; inode++) {
      for (int ielt = 0; ielt < nElts_; ielt++) {
        connectivityUnique_(inode,ielt) = 
          globalToUniqueMap_(connectivity_(inode,ielt));
      }
    }

    // form maximal elementset
    for (int i = 0; i < nElts_; i++) {
      elementSetAll_.insert(i); 
    }

  }

  // -------------------------------------------------------------
  //   map_to_element
  // -------------------------------------------------------------

  int FE_Uniform3DMesh::map_to_element(const DENS_VEC & x) const
  {

    // ix[i] is the element in the ith direction
    int ix[3];
    for (int i = 0; i < 3; i++) 
    {
      // handle points on upper boundary
      if (fabs(x(i)-borders_[1][i]) < tol) {
        ix[i] = nx_[i] - 1; 
      }
      // find bin
      else {
        ix[i] = (int)floor((x(i)-borders_[0][i])/dx_[i]); 
      }
      // handle points out-of-range [0:nx-1]  w/ periodicity
      if (ix[i] < 0 || ix[i] > nx_[i]-1) {
        if (periodicFlag_[i]) { 
          ix[i] = (ix[i] +nx_[i]) % nx_[i] ; // assume within a box length
        }
        else {
          string msg = "FE_Uniform3DMesh:: point maps outside of mesh, coordinate "
           + index_to_string(i) + " " + tostring(x(i))
           + " not in " + tostring(borders_[0][i]) + " : " + tostring(borders_[1][i]);
          throw ATC_Error(0, msg);
        }
      }
    }
    int elt = ix[2]*(nx_[0]*nx_[1]) + ix[1]*nx_[0] + ix[0];
    return elt;
  }

  int FE_Uniform3DMesh::map_to_element(const DENS_VEC & x,
                                       DENS_VEC & xi) const
  {

    xi.reset(3);
    // ix[i] is the element in the ith direction
    int ix[3];
    for (int i = 0; i < 3; i++) {
      // handle points on upper boundary
      if (fabs(x(i)-borders_[1][i]) < tol) {
        ix[i] = nx_[i] - 1; }
      // find bin
      else {
        ix[i] = (int)floor((x(i)-borders_[0][i])/dx_[i]); }

      // handle points out-of-range [0:nx-1]  w/ periodicity
      int shift = 0;
      if (ix[i] < 0 || ix[i] > nx_[i]-1) {
        if (periodicFlag_[i]) {
          if (ix[i] < 0) shift =  1;
          else           shift = -1;
          ix[i] = (ix[i]+shift*nx_[i]) % nx_[i] ; // assume within a box length
        } 
        else {
          string msg = "FE_Uniform3DMesh:: point maps outside of mesh, coordinate "
           + index_to_string(i) + " " + tostring(x(i))
           + " not in " + tostring(borders_[0][i]) + " : " + tostring(borders_[1][i]);
          throw ATC_Error(0, msg);
        }
      }
      double x_lower = borders_[0][i] + (ix[i] - shift*nx_[i])*dx_[i];
      xi(i) = 2.0*(x(i) - x_lower)/dx_[i] - 1.0;
    }
    int elt = ix[2]*(nx_[0]*nx_[1]) + ix[1]*nx_[0] + ix[0];
    return elt;
  }

  int FE_Uniform3DMesh::map_to_element(const DENS_VEC & x,
                                       DENS_VEC & xi,
                                       const Array<bool> & periodicFlag) const
  {

    xi.reset(3);
    // ix[i] is the element in the ith direction
    int ix[3];
    for (int i = 0; i < 3; i++) {
      // handle points on upper boundary
      if (fabs(x(i)-borders_[1][i]) < tol) {
        ix[i] = nx_[i] - 1; }
      // find bin
      else {
        ix[i] = (int)floor((x(i)-borders_[0][i])/dx_[i]); }

      // handle points out-of-range [0:nx-1]  w/ periodicity
      int shift = 0;
      if (ix[i] < 0 || ix[i] > nx_[i]-1) {
        if (periodicFlag(i)) {
          if (ix[i] < 0) shift =  1;
          else           shift = -1;
          ix[i] = (ix[i]+shift*nx_[i]) % nx_[i] ; // assume within a box length
        } 
        else {
          string msg = "FE_Uniform3DMesh:: point maps outside of mesh, coordinate "
           + index_to_string(i) + " " + tostring(x(i))
           + " not in " + tostring(borders_[0][i]) + " : " + tostring(borders_[1][i]);
          throw ATC_Error(0, msg);
        }
      }
      double x_lower = borders_[0][i] + (ix[i] - shift*nx_[i])*dx_[i];
      xi(i) = 2.0*(x(i) - x_lower)/dx_[i] - 1.0;
    }
    int elt = ix[2]*(nx_[0]*nx_[1]) + ix[1]*nx_[0] + ix[0];
    return elt;
  }

} // namespace ATC

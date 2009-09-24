#ifndef FE_MESH_H
#define FE_MESH_H

// ATC_Transfer headers
#include "Array.h"
#include "Array2D.h"
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

// Other headers
#include <vector>
#include <map>
#include <set>
#include <utility>

using namespace std;

namespace ATC {
  // Forward declarations
  class FE_Element;

  /**
   *  @class  FE_Mesh
   *  @brief  Base class for a finite element mesh
   */
  class FE_Mesh {

  public:

    /** constructor */
    FE_Mesh();
  
    /** destructor */
    virtual ~FE_Mesh();

    /** parser/modifier */
    bool modify(int narg, char **arg);

    /** evaluate shape function at real coordinates */
    void shape_functions(const VECTOR &x, 
                         DENS_VEC& shp, 
                         int & eltID,
                         Array<int>& node_list) const;

    /** evaluate shape function at real coordinates */
    void shape_functions(const VECTOR &x, 
                         DENS_VEC& shp,
                         int & eltID,
                         Array<int>& node_list,
                         const Array<bool>& periodicity) const;

    /** evaluate shape function at real coordinates */
    void shape_functions(const VECTOR &x, 
                         DENS_VEC & shp,
                         DENS_MAT & dshp,
                         int & eltID,
                         Array<int>& node_list) const;

    /** evaluate shape functions for all ip's on an element */
    // N is numIPsInElement X numNodesInElement
    void shape_function(const int eltID,
                        DENS_MAT &N,
                        DIAG_MAT &weights) const;

    /** evaluate shape functions for all ip's on an element */
    // N is numIPsInElement X numNodesInElement
    void shape_function(const int eltID,
                        DENS_MAT &N,
                        vector<DENS_MAT> &dN,
                        DIAG_MAT &weights) const;

    /** evaluate shape functions for all ip's on an face*/
    // N is numIPsInFace X numNodesInElement
    void face_shape_function(const PAIR &  face,
                             DENS_MAT &N,
                             vector<DENS_MAT> &dN,
                             vector<DENS_MAT> &Nn,
                             DIAG_MAT &weights) const;

    void face_shape_function(const PAIR &  face,
                             DENS_MAT &N,
                             DENS_MAT &n,
                             DIAG_MAT &weights) const;

    /** return connectivity (global ID numbers) for element eltID */
    void element_connectivity_global(const int eltID,
                                     Array<int> & nodes) const;

    void element_connectivity_unique(const int eltID,
                                     Array<int> & nodes) const;

    void face_connectivity(const pair<int,int> & faceID,
                                  Array<int> & nodes) const
    { int nNodesPerFace = get_nNodesPerFace();
      nodes.reset(nNodesPerFace);
      int eltID = faceID.first;
      int localFace = faceID.second;
      const Array2D<int> & localFaceConn = local_face_connectivity();
      for(int i = 0; i < nNodesPerFace; ++i) {
        nodes(i) = connectivity_(localFaceConn(localFace,i),eltID);
      }
    }
    void face_connectivity_unique(const pair<int,int> & faceID,
                                  Array<int> & nodes) const
    { int nNodesPerFace = get_nNodesPerFace();
      nodes.reset(nNodesPerFace);
      int eltID = faceID.first;
      int localFace = faceID.second;
      const Array2D<int> & localFaceConn = local_face_connectivity();
      for(int i = 0; i < nNodesPerFace; ++i) {
        nodes(i) = connectivityUnique_(localFaceConn(localFace,i),eltID);
      }
    }



    /** 
     *  return spatial coordinates for element nodes on eltID,
     *  indexed xCoords(isd,inode)
     */
    void element_coordinates(const int eltID,
                             DENS_MAT & xCoords) const;

    void face_coordinates(const pair <int,int> face,
                          DENS_MAT & xCoords) const;


    /** access to the nodal coordinate values */
    const DENS_MAT & nodal_coordinates(void) {return nodalCoords_  ;}

    /** access to nodal coordinates of a unique node */
    DENS_VEC nodal_coordinates(const int nodeID) const;

    /** access to nodal coordinates of a node */
    DENS_VEC global_coordinates(const int nodeID) const;

    /** access to the element connectivity values */
    const Array2D<int> & connectivity(void) {return connectivity_  ;}

    /** map spatial location to element */
    virtual int map_to_element(const DENS_VEC & x) const = 0;

    /** map spatial location to element and local coordinates */
    virtual int map_to_element(const DENS_VEC & x, 
                               DENS_VEC & xi) const = 0; 

    /** map spatial location to element and local coordinates */
    virtual int map_to_element(const DENS_VEC & x, 
                               DENS_VEC & xi, 
                               const Array<bool>& periodicity) const = 0; 

    /** map global node numbering to unique node numbering */
    int map_global_to_unique(const int global_id) const {return globalToUniqueMap_(global_id);}
    inline const int* global_to_unique_map_pointer(void) {return globalToUniqueMap_.get_data();}
    inline const Array<int>& global_to_unique_map(void) {return globalToUniqueMap_;}

    /** map unique node numbering a global node numbering */
    int map_unique_to_global(const int unique_id) 
    {return uniqueToGlobalMap_(unique_id);}
    inline const int* unique_to_global_map(void) {return uniqueToGlobalMap_.get_data();}

    /** query whether a nodeset with the given name exists */
    bool query_nodeset(const string & name) const;

    /** get node set (unique ID's) from the string name assigned to the set */
    const set<int> & get_nodeset(const string & name) const;

    /** create node set with tag "name" from nodes in given spatial range */
    void create_nodeset(const string & name,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax);

    /** get element set from the string name assigned to the set */
    const set<int> & get_elementset(const string & name) const;

    /** create element set with tag "name" from nodes in given spatial range */
    void create_elementset(const string & name,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax);

    /** get the minimal element set from a nodeset by name */
    void nodeset_to_minimal_elementset
      (const string & name, set<int> & elemSet) const;
    /** get the maximal element set from a nodeset by name */
    void nodeset_to_maximal_elementset
      (const string & name, set<int> & elemSet) const;
    /** get complement of element set by name */
    void elementset_complement(const string & name, set<int> & elemSet) const;
    void elementset_complement(const set<int> & elemSet, set<int> & elemSetComplement) const;
    /** get the node set from an element set by name */
    void elementset_to_minimal_nodeset
      (const string & name, set<int> & nodeSet) const;
    void elementset_to_nodeset(const string & name, set<int> & nodeSet) const;
    void elementset_to_nodeset(const set<int> & elemSet, set<int> & nodeSet) const;
    /** convert faceset to nodeset in  _unique_node numbering */
    void faceset_to_nodeset(const string &name, set<int> &nodeSet) const;
    void faceset_to_nodeset_global(const string &name, set<int> &nodeSet) const;

    /** get face set from the string name assigned to the set */
    const set< pair <int,int> > & get_faceset(const string & name) const;

    /** create face set with tag "name" from faces aligned with box */
    void create_faceset(const string & name,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax, 
                        bool outward);
    void create_faceset(const string & name, double x, int idir, int isgn);

    /** return number of spatial dimensions */
    int get_nSpatialDimensions() const { return nSD_; };

    /** return total number of nodes */
    int get_nNodes() const { return nNodes_; };

    /** return number of unique nodes */
    int get_nNodesUnique() const { return nNodesUnique_; };

    /** return number of elements */
    int get_nElements() const { return nElts_; };

    /** return number of integration points per element */
    int get_nIPsPerElement() const;

    /** return number of nodes per element */
    int get_nNodesPerElement() const;

    /** return number of faces per element */
    int get_nFacesPerElement() const;
  
    /** return number of nodes per face */
    int get_nNodesPerFace() const;
  
    /** return number of integration points per face */
    int get_nIPsPerFace() const;

    /** return scale in x */
    double get_xscale() const {return xscale_;}
    
    /** return scale in y */
    double get_yscale() const {return yscale_;}

    /** return scale in z */
    double get_zscale() const {return zscale_;}

    /** local face connectivity */
    const Array2D<int> & local_face_connectivity() const;
 
    /** element size in each direction */
    virtual void element_size(const int ielem, 
                              double hx, double hy, double hz) = 0; 

    /** element size in each direction */
    virtual double min_element_size(void) const = 0;

    /** set quadrature on element : a pass-through */
    void set_quadrature(int type); 

    /** output mesh subsets */
    void output(string prefix) const; 
  
  
  protected:

    /** number of spatial dimensions */
    int nSD_;

    /** number of elements */
    int nElts_;

    /** number of nodes */
    int nNodes_;
    int nNodesUnique_;

    /** periodicity flags */
    int periodicFlag_[3];

    /** element type for this mesh */
    FE_Element *feElement_;

    /** Nodal coordinates: nodalCoords_(nsd, numnode) */
    DENS_MAT nodalCoords_;

    /** Element connectivity: connectivity_(neltnode, nelt) */
    Array2D<int> connectivity_;
    Array2D<int> connectivityUnique_;

    /** map of global to unique node ID's */
    Array<int> globalToUniqueMap_;

    /** map of unique to global node ID's */
    Array<int> uniqueToGlobalMap_;

    /** map of string names to node sets */
    NODE_SET_MAP nodeSetMap_;

    /** maximal nodeset */
    set<int> nodeSetAll_;

    /** map of string names to node sets */
    FACE_SET_MAP faceSetMap_;

    /** map of string names to element sets */
    ELEMENT_SET_MAP elementSetMap_;

    /** maximal elementset */
    set<int> elementSetAll_;

    /** length scaling used by lammps */
    double xscale_, yscale_, zscale_;

    /** scratch memory for setting face connectivities */
    Array<int> connScratch_;

  };

  /**
   *  @class  FE_Uniform3DMesh
   *  @brief  Derived class for a uniform mesh, a structured mesh with 
   *          fixed element sizes in x, y, and z directions 
   */

  class FE_Uniform3DMesh : public FE_Mesh {
  
  public:

    /** constructor */
    FE_Uniform3DMesh(const int nx,
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
                     int zperiodic);

    /** destructor */
    ~FE_Uniform3DMesh();

    virtual void element_size(const int ielem, 
                              double hx, double hy, double hz) 
    { hx = Lx_[0]/nx_[0]; hy = Lx_[1]/nx_[1]; hz = Lx_[2]/nx_[2]; }

    virtual double min_element_size(void) const
      {return min(Lx_[0]/nx_[0], min(Lx_[1]/nx_[1], Lx_[2]/nx_[2])); }

    /** map spatial location to element */
    virtual int map_to_element(const DENS_VEC & x) const;
  
    /** map spatial location to element and local coordinates */
    virtual int map_to_element(const DENS_VEC & x, 
                               DENS_VEC & xi) const; 

    /** map spatial location to element and local coordinates */
    virtual int map_to_element(const DENS_VEC & x, 
                               DENS_VEC & xi,
                               const Array<bool> & periodicity) const; 


  protected:

    // Number of elements in each spatial direction
    int nx_[3];

    // Bounds of region on which mesh is defined
    double borders_[2][3];

    // Region size in each direction
    double Lx_[3];

    // Element size in each direction
    double dx_[3];

    /** create global-to-unique node mapping */
    void setup_periodicity();

  };

}; // namespace ATC_Transfer

#endif // FE_MESH_H

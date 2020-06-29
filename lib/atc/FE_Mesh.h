#ifndef FE_MESH_H
#define FE_MESH_H

#include "Array.h"
#include "Array2D.h"
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"
#include "KD_Tree.h"
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <cfloat>
#include <string>
#include <vector>
#include "mpi.h"

namespace ATC {

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

    /** initialization */
    void initialize(void);

    /** write an unstructured mesh */
    void write_mesh(std::string meshFile);



// SJL why? will they ever be overridden by derived classes? in what
//     situation would that be required, or desirable? virtual functions
//     are slightly less efficient because they cannot be hard-linked at
//     compile time.

    bool is_partitioned() const { return partitioned_; }
    virtual void partition_mesh() = 0;
    virtual void departition_mesh() = 0;

    int map_elem_to_myElem(int elemID) const;
    int map_myElem_to_elem(int myElemID) const;
    void set_lammps_partition(bool t) {lammpsPartition_ = t;}
    void set_data_decomposition(bool t) {decomposition_ = t;}

    /** set quadrature on element from within interpolate class */
    void set_quadrature(FeIntQuadrature type);

    /** evaluate shape function at real coordinates */
    void position(const int elem,
                  const VECTOR &xi, 
                        DENS_VEC &x) const;

    /** evaluate shape function at real coordinates */
    void shape_functions(const VECTOR &x, 
                         DENS_VEC &N,
                         Array<int> &nodeList) const;

    /** evaluate shape function at real coordinates */
    void shape_functions(const VECTOR &x, 
                                 DENS_VEC &N,
                                 Array<int> &nodeList,
                                 const Array<bool> &) const;

    /** evaluate shape function at real coordinates */
    void shape_functions(const DENS_VEC &x, 
                         DENS_VEC &N,
                         DENS_MAT &dNdx,
                         Array<int> &nodeList) const;

    /** evaluate shape function at real coordinates */
    void shape_functions(const VECTOR &x, 
                         const int eltID,
                         DENS_VEC &N,
                         Array<int> &nodeList) const;

    /** evaluate shape function at real coordinates */
    void shape_functions(const DENS_VEC &x, 
                         const int eltID,
                         DENS_VEC &N,
                         DENS_MAT &dNdx,
                         Array<int> &nodeList) const;

    /** evaluate shape function at real coordinates */
    void shape_function_derivatives(const DENS_VEC &x, 
                         const int eltID,
                         DENS_MAT &dNdx,
                         Array<int> &nodeList) const;

    /** evaluate shape functions for all ip's on an element */
    // N is numIPsInElement X numNodesInElement
    void shape_function(const int eltID,
                        DENS_MAT &N,
                        DIAG_MAT &weights) const;

    /** evaluate shape functions for all ip's on an element */
    // N is numIPsInElement X numNodesInElement
    void shape_function(const int eltID,
                        DENS_MAT &N,
                        std::vector<DENS_MAT> &dN,
                        DIAG_MAT &weights) const;

    /** evaluate shape functions for all ip's on a face */
    // N is numIPsInFace X numNodesInElement
    void face_shape_function(const PAIR &face,
                             DENS_MAT &N,
                             DENS_MAT &n,
                             DIAG_MAT &weights) const;

    void face_shape_function(const PAIR &face,
                             DENS_MAT &N,
                             std::vector<DENS_MAT> &dN,
                             std::vector<DENS_MAT> &Nn,
                             DIAG_MAT &weights) const;

    /** compute normal vector from the specified face */
    double face_normal(const PAIR &face,
                       const int ip,
                       DENS_VEC &normal) const;


    /** return connectivity (global ID numbers) for element eltID */
    void element_connectivity_global(const int eltID,
                                     Array<int> & nodes) const;

    void element_connectivity_unique(const int eltID,
                                     Array<int> & nodes) const;

    int element_connectivity_global(const int eltID,
                                     const int inode) const;

    int element_connectivity_unique(const int eltID,
                                     const int inode) const;

    AliasArray<int> element_connectivity_global(const int eltID) const;

    AliasArray<int> element_connectivity_unique(const int eltID) const;

    void face_connectivity(const PAIR & faceID,
                           Array<int> & nodes) const
    { int nNodesPerFace = num_nodes_per_face();
      nodes.reset(nNodesPerFace);
      int eltID = faceID.first;
      int localFace = faceID.second;
      const Array2D<int> & localFaceConn = local_face_connectivity();
      for(int i = 0; i < nNodesPerFace; ++i) {
        nodes(i) = element_connectivity_global(eltID, localFaceConn(localFace,i));
      }
    }
    void face_connectivity_unique(const PAIR & faceID,
                                  Array<int> & nodes) const
    { int nNodesPerFace = num_nodes_per_face();
      nodes.reset(nNodesPerFace);
      int eltID = faceID.first;
      int localFace = faceID.second;
      const Array2D<int> & localFaceConn = local_face_connectivity();
      for(int i = 0; i < nNodesPerFace; ++i) {
        nodes(i) = element_connectivity_unique(eltID, localFaceConn(localFace,i));
      }
    }

    /** 
     *  return spatial coordinates for element nodes on eltID,
     *  indexed xCoords(isd,inode)
     */
    void element_coordinates(const int eltID,
                             DENS_MAT & xCoords) const;

    void face_coordinates(const PAIR face,
                          DENS_MAT & xCoords) const;

    /** access to the nodal coordinate values */
    const DENS_MAT & nodal_coordinates(void) const {return nodalCoords_  ;}

    /** access to nodal coordinates of a unique node */
    DENS_VEC nodal_coordinates(const int nodeID) const;

    /** access to nodal coordinates of a node */
    DENS_VEC global_coordinates(const int nodeID) const;

    /** map spatial location to element */
    virtual int map_to_element(const DENS_VEC &x) const = 0;

    /** map global node numbering to unique node numbering */
    int map_global_to_unique(const int global_id) const 
    {
      return globalToUniqueMap_(global_id);
    }
    inline const Array<int>& global_to_unique_map(void) const 
    {
      return globalToUniqueMap_;
    }

    /** map unique node numbering a global node numbering */
    int map_unique_to_global(const int unique_id) 
    {
      return uniqueToGlobalMap_(unique_id);
    }
    inline const Array<int>& unique_to_global_map(void) const 
    {
      return uniqueToGlobalMap_;
    }

    /** query whether a nodeset with the given name exists */
    bool query_nodeset(const std::string & name) const;

    /** get node set (unique ID's) from the string name assigned to the set */
    const std::set<int> & nodeset(const std::string & name) const;

    /** create node set with tag "name" from nodes in given spatial range */
    void create_nodeset(const std::string & name, const std::set<int> & nodeset);
    void create_nodeset(const std::string & name,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax);

    /** add to node set with tag "name" from nodes in given spatial range */
    void add_to_nodeset(const std::string & name,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax);

    /** get element set from the string name assigned to the set */
    const std::set<int> & elementset(const std::string & name) const;

    /** create element set with tag "name" from nodes in given spatial range */
    void create_elementset(const std::string & name,
                           double xmin, double xmax,
                           double ymin, double ymax,
                           double zmin, double zmax);


    /** get the minimal element set from a nodeset by name */
    void nodeset_to_minimal_elementset(const std::string &name, 
                                       std::set<int> &elemSet) const;
    /** get the minimal element set from a set of nodes */
    void nodeset_to_minimal_elementset(const std::set<int> &nodeSet, 
                                       std::set<int> &elemSet) const;
    /** get the maximal element set from a nodeset by name */
    void nodeset_to_maximal_elementset(const std::string &name, 
                                       std::set<int> &elemSet) const;
    /** get the maximal element set from a set of nodes */
    void nodeset_to_maximal_elementset(const std::set<int> &nodeSet, 
                                       std::set<int> &elemSet) const;

    /** get complement of element set by name */
    void elementset_complement(const std::string &name, 
                               std::set<int> &elemSet) const;
    void elementset_complement(const std::set<int> &elemSet, 
                               std::set<int> &elemSetComplement) const;

    /** get the node set from an element set by name */
    void elementset_to_minimal_nodeset(const std::string &name, 
                                       std::set<int> &nodeSet) const;

    void elementset_to_nodeset(const std::string &name, 
                               std::set<int> nodeSet) const;
    void elementset_to_nodeset(const std::set<int> &elemSet, 
                               std::set<int> nodeSet) const;
    std::set<int> elementset_to_nodeset(const std::string &name) const;

    /** convert faceset to nodeset in _unique_ node numbering */
    void faceset_to_nodeset(const std::string &name, 
                            std::set<int> &nodeSet) const;
    void faceset_to_nodeset(const std::set<PAIR> &faceSet, 
                            std::set<int> &nodeSet) const;

    void faceset_to_nodeset_global(const std::string &name, 
                                   std::set<int> &nodeSet) const;
    void faceset_to_nodeset_global(const std::set<PAIR> &faceSet, 
                                   std::set<int> &nodeSet) const;

    /** get face set from the string name assigned to the set */
    const std::set< std::pair<int,int> > & faceset(const std::string & name) const;

    /** create face set with tag "name" from faces aligned with box */
    void create_faceset(const std::string & name,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax, 
                        bool outward);
    /** create face set with tag "name" from faces aligned with plane */
    void create_faceset(const std::string & name, double x, int idir, int isgn,
                        int nIdx2=-1, double x2lo=0.0, double x2hi=0.0,
                        int nIdx3=-1, double x3lo=0.0, double x3hi=0.0);

    /** cut mesh */
    virtual void cut_mesh(const std::set<PAIR> & faceSet, const std::set<int> & nodeSet) = 0;
    
    /** delete elements */
    virtual void delete_elements(const std::set<int> & elementList) = 0;

    /** return number of spatial dimensions */
    int num_spatial_dimensions() const { return nSD_; }

    /** return total number of nodes */
    int num_nodes() const { return nNodes_; }

    /** return number of unique nodes */
    int num_nodes_unique() const { return nNodesUnique_; }

    /** return number of elements */
    int num_elements() const { return nElts_; }

    /** return number of elements partitioned to my processor */
    int my_num_elements() const { return myNElts_; }

    /** return number of integration points per element */
    int num_ips_per_element() const;

    /** return number of nodes per element */
    int num_nodes_per_element() const;

    /** return number of faces per element */
    int num_faces_per_element() const;
  
    /** return number of nodes per face */
    int num_nodes_per_face() const;
  
    /** return number of integration points per face */
    int num_ips_per_face() const;

    /** return a pointer to the connectivity. This function will only work
        when mesh is not partitioned. */
    Array2D<int> * connectivity(void) { return &connectivity_; }
    /** return a pointer to the connectivity */
    DENS_MAT * coordinates(void) { return &nodalCoords_;} 
    /** Engine nodeMap stuff  */
    Array<int> *node_map(void) { return &globalToUniqueMap_;}


    /** return scale in x,y,z */
    double xscale() const { return xscale_; }
    double yscale() const { return yscale_; }
    double zscale() const { return zscale_; }

    /** local face connectivity */
    const Array2D<int> & local_face_connectivity() const;
 
    /** element size in each direction */
    virtual void bounding_box(const int ielem, 
                              DENS_VEC & xmin, DENS_VEC & xmax); 

    /** element size in each direction */
    virtual void element_size(const int ielem, 
                              double & hx, double & hy, double & hz); 

    /** element size in each direction */
    virtual double min_element_size(void) const {return 0.0 ;} 

    /** get nodal coordinates for a given element */ 
    void element_field(const int eltIdx, const DENS_MAT f,
                       DENS_MAT &local_field)
    {
      int dof = f.nCols();
      Array<int> nodes;
      element_connectivity_unique(eltIdx,nodes);
      local_field.reset(num_nodes_per_element(), dof);
      for (int i = 0; i < nodes.size(); i++) {
        for (int j = 0; j < dof; j++) local_field(i,j) = f(nodes(i),j);
      }
    }

    /** almost structured */
    bool is_aligned(void) const;

    /** extruded */
    bool is_two_dimensional(void) const {return twoDimensional_;}

   virtual double coordinate_tolerance(void) const {return 1.e-8;}

    /** element type */
    std::string element_type(void) const ; 

    /** output mesh subsets */
    void output(std::string prefix) const; 
  
    /* Parallelization data members */
 
    /** return element vector for this processor */
    const std::vector<int> & owned_elts() const { return myElts_; }
    const std::vector<int> & owned_and_ghost_elts() const { 
      return (decomposition_) ? myAndGhostElts_: myElts_; }
    bool is_owned_elt(int elt) const;
    
  protected:

    void parse_plane(int & argIdx, int narg, char ** arg,  
      int & ndir, int * idir, int & isgn, double xlimits[][2]);

    void parse_units(int & argIdx, int narg, char ** arg,  
      double & xmin, double & xmax, double & ymin, double & ymax, double & zmin, double & zmax);

    /** will this mesh use data decomposition? */
    bool decomposition_;

    /** should the mesh use the native lammps partitioning? */
    bool lammpsPartition_;

    /** is the element/node data currently partitioned among processors? */
    bool partitioned_;

    /** number of spatial dimensions */
    int nSD_;

    /** number of elements */
    int nElts_;
    /** number of elements partitioned to this processor */
    int myNElts_;

    /** number of nodes */
    int nNodes_;
    int nNodesUnique_;

    /** periodicity flags */
    Array<bool> periodicity_;

    /** element type for this mesh */
    FE_Element *feElement_;

    /** Nodal coordinates: nodalCoords_(nsd, numnode) */
    DENS_MAT nodalCoords_;

    /** Element connectivity: connectivity_(neltnode, nelt) */
    Array2D<int> connectivity_;
    Array2D<int> myConnectivity_;
    Array2D<int> connectivityUnique_;
    Array2D<int> myConnectivityUnique_;

    /** map from unique node id to associated elements for periodic meshes */
    /** this data structure is only ever accessed from an unpartitioned context */
    Array<std::vector<int> > uniqueNodeToElementMap_;

    /** map of global to unique node ID's */
    Array<int> globalToUniqueMap_;
    Array<int> compactRemap_; // for condensing unique numbering

    /** map of unique to global node ID's */
    Array<int> uniqueToGlobalMap_;

    /** map of string names to node sets */
    NODE_SET_MAP nodeSetMap_;

    /** maximal nodeset */
    std::set<int> nodeSetAll_;

    /** map of string names to node sets */
    FACE_SET_MAP faceSetMap_;

    /** map of string names to element sets */
    ELEMENT_SET_MAP elementSetMap_;

    /** maximal elementset */
    std::set<int> elementSetAll_;

    /** length scaling used by lammps */
    double xscale_, yscale_, zscale_;

    /** Processor demarcations */
    std::vector<double> procs_;

    /** Dimension (x=0, y=1, or z=2) */
    int partitionAxis_;

    /** List of nodes for this processor */
    std::vector<int> myNodes_;

    /** List of elements for this processor */
    std::vector<int> myElts_;
    std::vector<int> myAndGhostElts_;

    /** maps between my IDs and the total IDs */
    std::map<int,int> elemToMyElemMap_;
   
    /** Lists of ghost nodes/neighbor ghost nodes */
    std::vector<int> ghostNodesL_; 
    std::vector<int> ghostNodesR_;
    std::vector<int> shareNodesL_;
    std::vector<int> shareNodesR_;

    /** extruded */
    bool twoDimensional_;
    bool hasPlanarFaces_;
    double coordTol_;

  };

  /**
   *  @class  FE_3DMesh
   *  @brief  Derived class for an unstructured 3D mesh
   */
  class FE_3DMesh : public FE_Mesh {
  public:
    /** constructor */
    FE_3DMesh(){}; 

    /** constructor for read-in mesh **/
    // can later be extended to take nodesets, elementsets, etc.
    FE_3DMesh(const std::string elementType, 
              const int nNodes,
              const int nElements,
              const Array2D<int> *connectivity,
              const DENS_MAT *nodalCoordinates,
              const Array<bool> periodicity,
              const Array<std::pair<std::string,std::set<int> > > *nodeSets);

    /** destructor */
    virtual ~FE_3DMesh();

    void partition_mesh(void);

    void departition_mesh(void);

    void lammps_partition_mesh(void);

    /** Removes duplicate elements that appear in more than one vector
         within procEltLists. **/
    void prune_duplicate_elements(std::vector<std::vector<int> > &procEltLists, 
                                  int *eltToOwners);
    
    /** Takes procEltLists, and if there are more than nProcs of them
        it takes the extra elements and distributes them to other vectors
        in procEltLists. */
          
          //       processors if during pruning processors end up
          //       elementless. This is slightly complicated because of
          //       ghost nodes.
    void redistribute_extra_proclists(std::vector<std::vector<int> > &procEltLists,
                                      int *eltToOwners, int nProcs);
                                        
    /** This takes in a dense matrix and a list of elements
        and fills in a standard adjacency list (within the matrix)
        for those elements. **/
          
          //       the set intersection, which does redundant computations
          //       right now, and filling in the adjacencies for both elements
          //       simultaneously when two elements share a face.
    void compute_face_adjacencies(const std::list<int> &elts, 
                                  DENS_MAT &faceAdjacencies);

    /** Counts the number of nonempty vectors in a vector of vectors. **/
    int numNonempty(std::vector<std::vector<int> > &v);
      
    /**  In the partitioning, we want to sort vectors of integers by size, 
          and furthermore we want empty vectors to count as the "largest" 
          possible vector because they dont want to count in the minimum. **/
    struct vectorComparer {
        bool operator() (std::vector<int> l, std::vector<int> r) {
          if (l.empty())
            return false;
          if (r.empty())
            return true;
          return (l.size() < r.size());
        }
    } vectorCompSize;

    virtual double min_element_size(void) const {return minEltSize_; }
    virtual double coordinate_tolerance(void) const {
      return 0.25*(this->min_element_size()); // loose
      //return 0.5;
    }
    virtual void cut_mesh(const std::set<PAIR> &faceSet,
                          const std::set<int> &nodeSet);
    
    virtual void delete_elements(const std::set<int> &elementSet);

    /** map spatial location to element */
    virtual int map_to_element(const DENS_VEC &x) const;
    
    /** sends out data to processors during partitioning */
    void distribute_mesh_data();
  protected:
    /** create global-to-unique node mapping */
    virtual void setup_periodicity(double tol);
    virtual void setup_periodicity() { setup_periodicity(1.e-8); }
    void fix_periodicity  (int idim);
    int find_boundary_nodes(int idim, std::set<int> & nodes);
    bool match_nodes(int idim, std::set<int> & top, std::set<int> & bot,
                     Array<int> & map);
    void set_unique_connectivity(void);
    bool orient(int idir);

    double minEltSize_;
    KD_Tree *tree_;

    /** test if a specified element actually contains the given point */
    bool contains_point(const int eltID, const DENS_VEC & x) const;

  private:
    Array<std::vector<int> > nodeToParentElements_;

  };

  /**
   *  @class  FE_Rectangular3DMesh
   *  @brief  Derived class for a structured mesh with 
   *          variable element sizes in x, y, and z directions 
   */
  class FE_Rectangular3DMesh : public FE_3DMesh {
  public:
    /** constructor */
    FE_Rectangular3DMesh(){}; 
    FE_Rectangular3DMesh(
              const Array<double> & hx,
              const Array<double> & hy,
              const Array<double> & hz,
              const double xmin, const double xmax, 
              const double ymin, const double ymax,
              const double zmin, const double zmax,
              const Array<bool> periodicity,
              const double xscale=1,
              const double yscale=1,
              const double zscale=1);

    /** destructor */
    virtual ~FE_Rectangular3DMesh() {};
    
    void partition_mesh(void);

    void departition_mesh(void);

    /** map spatial location to element */
    virtual int map_to_element(const DENS_VEC &x) const; 

  protected:

    /** Number of elements in each spatial direction */
    int n_[3];

    /** Bounds of region on which mesh is defined */
    double borders_[2][3];

    /** Region size in each direction */
    double L_[3];

    

    /** create global-to-unique node mapping */
    virtual void setup_periodicity(); // note no "tol"

  private: // only used by this class
    /** partitions in x,y,z */
    Array<double> hx_, hy_, hz_;

    /** nodal locations */
    std::vector< Array<double> > x_;
  };

  /**
   *  @class  FE_Uniform3DMesh
   *  @brief  Derived class for a uniform structured mesh with 
   *          fixed element sizes in x, y, and z directions 
   */
  class FE_Uniform3DMesh : public FE_Rectangular3DMesh {
  
  public:

    /** constructor */
    FE_Uniform3DMesh(const int nx,
                     const int ny,
                     const int nz,
                     const double xmin, const double xmax,
                     const double ymin, const double ymax,
                     const double zmin, const double zmax,
                     const Array<bool> periodicity,
                     const double xscale=1,
                     const double yscale=1,
                     const double zscale=1);

    /** destructor */
    virtual ~FE_Uniform3DMesh();

    void partition_mesh(void);

    void departition_mesh(void);
    
    virtual void element_size(const int /* ielem */, 
                              double &hx, double &hy, double &hz)
    { hx = L_[0]/n_[0]; hy = L_[1]/n_[1]; hz = L_[2]/n_[2]; }

    virtual double min_element_size(void) const
    { return std::min(L_[0]/n_[0], std::min(L_[1]/n_[1], L_[2]/n_[2])); }

    /** map spatial location to element */
    virtual int map_to_element(const DENS_VEC &x) const; 

  private: // only used by this class
    /** Element size in each direction */
    double dx_[3];

  };

} // namespace ATC_Transfer

#endif // FE_MESH_H

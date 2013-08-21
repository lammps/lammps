#ifndef FE_ENGINE_H
#define FE_ENGINE_H

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <string>
#include "ATC_TypeDefs.h"
#include "Array.h"
#include "Array2D.h"
#include "FE_Mesh.h"
#include "PhysicsModel.h"
#include "OutputManager.h"
#include "MeshReader.h"
#include "mpi.h"

namespace ATC {
  
  class ATC_Method;
  class FE_Element;
  class XT_Function;
  class KernelFunction;
  
  /**
   *  @class  FE_Engine 
   *  @brief  Base class for computing and assembling mass matrix 
   *                                              and rhs vectors
   */  
  
  class FE_Engine{
  public:
    /** constructor/s  */
    FE_Engine(MPI_Comm comm);
  
    /** destructor */
    ~FE_Engine();
  
    /** initialize */
    void initialize();
  
    MPI_Comm communicator() {return communicator_;}
    void partition_mesh();
    void departition_mesh();
    bool is_partitioned() const { return feMesh_->is_partitioned(); }
    int map_elem_to_myElem(int elemID) const 
    { return feMesh_->map_elem_to_myElem(elemID); }
    int map_myElem_to_elem(int myElemID) const 
    { return feMesh_->map_myElem_to_elem(myElemID); }

    // note: it is misleading to declare the following const
    //       because it touches the nIPsPer* data members, which
    //       are now declared mutable. Why? Well, set_quadrature
    //       has to be called from a const function, and all the
    //       matrices dependent on nIPsPer* are declared mutable
    //       as well (and have been). I think this is because a
    //       const engine needs to be able to deal with various
    //       quadratures and update its data members directly, which
    //       are really convenience-copies of data members that
    //       are more pertinent to other classes (FE_Interpolate,
    //       for the most part) that it uses temporarily for space/
    //       time speedups while doing it's computations.
    //       
    //       I approve of this usage of mutable, but the const/
    //       non-const member function declaring in this class is
    //       really all wrong to begin with.

    /** set quadrature scheme, resize matrices if necessary as per
     *  initialize() */
    void set_quadrature(FeIntQuadrature quadType, bool temp=true) const;

    /** parser/modifier */
    bool modify(int narg, char **arg);
  
    /** finish up */
    void finish();

    /** print out the "global connectivity" of all elements */
    void print_mesh() const;

    //----------------------------------------------------------------
    /** \name output */
    //----------------------------------------------------------------
    /*@{*/
    /** these assume the caller is handling the parallel collection */
    void initialize_output(int rank, std::string outputPrefix, std::set<int> otypes);
  
    /** write geometry */
    void write_geometry(void);
  
    /** write data: data is arrayed over _unique_ nodes 
                    and then mapped by the engine */
    void write_data(double time, FIELDS &soln, OUTPUT_LIST *data=NULL);
    void write_data(double time, OUTPUT_LIST *data);

    void write_restart_file(std::string fileName, RESTART_LIST *data)
    { outputManager_.write_restart_file(fileName,data); }

    void read_restart_file(std::string fileName, RESTART_LIST *data)
    { outputManager_.read_restart_file(fileName,data); }

    void delete_elements(const std::set<int> &elementList);
    void cut_mesh(const std::set<PAIR> &cutFaces, const std::set<int> &edgeNodes);

    void add_global(const std::string name, const double value) 
    { outputManager_.add_global(name,value); }

    void add_field_names(const std::string field, const std::vector<std::string> & names)
    { outputManager_.add_field_names(field,names); }

    void reset_globals() { outputManager_.reset_globals(); }

    /** pass through to access output manager */
    OutputManager *output_manager() { return &outputManager_; }
    /*@}*/

    //----------------------------------------------------------------
    /** \name assembled matrices and vectors */
    //----------------------------------------------------------------
    /*@{*/
    DENS_VEC interpolate_field(const DENS_VEC & x, const FIELD & f) const;

    /** interpolate fields */
    void interpolate_fields(const int               ielem,
                            const FIELDS            &fields,
                            AliasArray<int>         &conn,
                            DENS_MAT                &N,
                            DIAG_MAT                &weights,
                            std::map<FieldName,DENS_MAT> &fieldsAtIPs) const;

    /** interpolate fields & gradients */
    void interpolate_fields(const int       ielem,
                            const FIELDS    &fields,
                            AliasArray<int> &conn,
                            DENS_MAT        &N,
                            DENS_MAT_VEC    &dN,
                            DIAG_MAT        &weights,
                            FIELD_MATS      &fieldsAtIPs,
                            GRAD_FIELD_MATS &grad_fieldsAtIPs) const;

    /** compute a dimensionless stiffness matrix */
    void stiffness_matrix(SPAR_MAT &matrix) const;


    /** compute tangent matrix for a pair of fields - native quadrature */
    void compute_tangent_matrix(
                const RHS_MASK          &rhsMask,
                const std::pair<FieldName,FieldName> row_col,
                const FIELDS            &fields,
                const PhysicsModel      *physicsModel,
                const Array<int>        &elementMaterials,
                SPAR_MAT                &tangent,
                const DenseMatrix<bool> *elementMask=NULL) const;

    /** compute tangent matrix for a pair of fields - given quadrature */
    void compute_tangent_matrix(const RHS_MASK &rhsMask,
                const std::pair<FieldName,FieldName> row_col,
                const FIELDS            &fields,
                const PhysicsModel      *physicsModel,
                const Array<std::set<int> >  &pointMaterialGroups,
                const DIAG_MAT          &weights,
                const SPAR_MAT          &N,
                const SPAR_MAT_VEC      &dN,
                SPAR_MAT                &tangent,
                const DenseMatrix<bool> *elementMask=NULL) const;

    /** compute a consistent mass matrix for a field */
    void compute_mass_matrix(
                const Array<FieldName>  &mask,
                const FIELDS            &fields,
                const PhysicsModel      *physicsModel,
                const Array<int>        &elementMaterials,
                CON_MASS_MATS           &mass_matrix,
                const DenseMatrix<bool> *elementMask=NULL) const;

    /** compute a dimensionless mass matrix */
    void compute_mass_matrix(SPAR_MAT &mass_matrix) const;
    
    /** computes a dimensionless mass matrix for the given-quadrature */
    void compute_mass_matrix(const DIAG_MAT &weights,
                             const SPAR_MAT &N,
                             SPAR_MAT       &mass_matrix) const;
    
    /** compute a single dimensionless mass matrix */
    void compute_lumped_mass_matrix(DIAG_MAT &lumped_mass_matrix) const;

    /** compute lumped mass matrix = diag (\int \rho N_I dV) */
    void compute_lumped_mass_matrix(
                const Array<FieldName>  &mask, 
                const FIELDS            &fields,
                const PhysicsModel      *physicsModel,
                const Array<int>        &elementMaterials,
                MASS_MATS               &mass_matrix,
                const DenseMatrix<bool> *elementMask=NULL) const;

    /** compute dimensional lumped mass matrix using given quadrature */
    void compute_lumped_mass_matrix(
                       const Array<FieldName> &mask, 
                       const FIELDS           &fields,
                       const PhysicsModel     *physicsModel,
                       const Array<std::set<int> > &pointMaterialGroups,
                       const DIAG_MAT         &weights,
                       const SPAR_MAT         &N,
                       MASS_MATS              &mass_matrix) const;

    /** compute an approximation to a finite difference gradient from mesh */
    void compute_gradient_matrix(SPAR_MAT_VEC &grad_matrix) const;

    /** compute energy */
    void compute_energy(const Array<FieldName>  &mask,
                        const FIELDS            &fields,
                        const PhysicsModel      *physicsModel,
                        const Array<int>        &elementMaterials,
                        FIELD_MATS              &energy, 
                        const DenseMatrix<bool> *elementMask=NULL,
                        const IntegrationDomainType domain=FULL_DOMAIN) const;

    /** compute residual or RHS of the dynamic weak eqn */
    void compute_rhs_vector(
                    const RHS_MASK          &rhsMask,
                    const FIELDS            &fields,
                    const PhysicsModel      *physicsModel,
                    const Array<int>        &elementMaterials,
                    FIELDS                  &rhs,
                    const DenseMatrix<bool> *elementMask=NULL) const;

    /** compute RHS for given quadrature */
    void compute_rhs_vector(const RHS_MASK         &rhsMask,
                            const FIELDS           &fields,
                            const PhysicsModel     *physicsModel,
                            const Array<std::set<int> > &pointMaterialGroups,
                            const DIAG_MAT         &weights,
                            const SPAR_MAT         &N,
                            const SPAR_MAT_VEC     &dN,
                            FIELDS                 &rhs) const;

    /** compute pointwise source for given quadrature */
    void compute_source(const Array2D<bool>    &rhsMask,
                        const FIELDS           &fields,
                        const PhysicsModel     *physicsModel,
                        const Array<std::set<int> > &pointMaterialGroups,
                        const DIAG_MAT         &weights,
                        const SPAR_MAT         &N,
                        const SPAR_MAT_VEC     &dN,
                        FIELD_MATS             &sources) const;

    /** compute flux in domain i.e. N^T B_integrand */
    void compute_flux(const RHS_MASK          &rhsMask,
                      const FIELDS            &fields,
                      const PhysicsModel      *physicsModel,
                      const Array<int>        &elementMaterials,
                      GRAD_FIELD_MATS         &flux, 
                      const DenseMatrix<bool> *elementMask=NULL) const;

    /** compute the flux on the MD/FE boundary */
    void compute_boundary_flux(const RHS_MASK     &rhsMask, 
                               const FIELDS       &fields,
                               const PhysicsModel *physicsModel,
                               const Array<int>   &elementMaterials,
                               const std::set<PAIR>    &faceSet, 
                               FIELDS             &rhs) const;

    /** compute the flux on using an L2 interpolation of the flux */
    void compute_boundary_flux(const RHS_MASK         &rhsMask,
                               const FIELDS           &fields,
                               const PhysicsModel     *physicsModel,
                               const Array<int>       &elementMaterials,
                               const Array<std::set<int> > &pointMaterialGroups,
                               const DIAG_MAT         &weights,
                               const SPAR_MAT         &N,
                               const SPAR_MAT_VEC     &dN,
                               const DIAG_MAT         &flux_mask,
                               FIELDS                 &rhs,
                               const DenseMatrix<bool> *elementMask=NULL,
                               const std::set<int>          *nodeSet=NULL) const;

    /** compute prescribed flux given an array of functions of x & t */
    void add_fluxes(const Array<bool>    &fieldMask,  
                    const double         time,
                    const SURFACE_SOURCE &sourceFunctions, 
                    FIELDS               &nodalSources) const;
    void compute_fluxes(const Array<bool>    &fieldMask,  
                    const double         time,
                    const SURFACE_SOURCE &sourceFunctions, 
                    FIELDS               &nodalSources) const
    {
      SURFACE_SOURCE::const_iterator src_iter;
      for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++) {
        _fieldName_ = src_iter->first;
        if (!fieldMask((int)_fieldName_)) continue;
        if (nodalSources[_fieldName_].nRows()==0) { 
           nodalSources[_fieldName_].reset(nNodesUnique_,1); 
        }
      }
      add_fluxes(fieldMask,  time, sourceFunctions, nodalSources);
    }

    /** compute prescribed flux given an array of functions of u, x & t */
    void add_robin_fluxes(const Array2D<bool> &rhsMask,
                          const FIELDS        &fields,
                          const double        time,
                          const ROBIN_SURFACE_SOURCE &sourceFunctions,
                          FIELDS              &nodalSources) const;

    void add_robin_tangent(const Array2D<bool> &rhsMask,
                           const FIELDS        &fields,
                           const double        time,
                           const ROBIN_SURFACE_SOURCE &sourceFunctions,
                           SPAR_MAT            &tangent) const;

    /** compute nodal vector of volume based sources */
    void add_sources(const Array<bool>   &fieldMask, 
                     const double        time,
                     const VOLUME_SOURCE &sourceFunctions, 
                     FIELDS              &nodalSources) const;

    /** compute surface flux of a nodal field */
    void field_surface_flux(const DENS_MAT  &field,
                            const std::set<PAIR> &faceSet,
                            DENS_MAT        &values,
                            const bool      contour=false,
                            const int       axis=2) const;

    /** integrate a nodal field over an element set */
    DENS_VEC integrate(const DENS_MAT  &field, const ESET & eset) const;


    /** integrate a nodal field over an face set */
    DENS_VEC integrate(const DENS_MAT  &field, const FSET & fset) const
    { throw ATC_Error(FILELINE,"unimplemented function"); }

    /*@}*/

    //----------------------------------------------------------------
    /** \name shape functions */
    //----------------------------------------------------------------
    /*@{*/

    /** evaluate shape function at a list of points in R^3 */
    void evaluate_shape_functions(const MATRIX &coords, 
                                  SPAR_MAT &N) const;

    /** evaluate shape function & derivatives at a list of points in R^3 */
    void evaluate_shape_functions(const MATRIX &coords, 
                                  SPAR_MAT &N,
                                  SPAR_MAT_VEC &dN) const;

    /** evaluate shape function at a list of points in R^3 */
    void evaluate_shape_functions(const MATRIX &coords, 
                                  const INT_ARRAY &pointToEltMap,
                                  SPAR_MAT &N) const;

    /** evaluate shape function & derivatives at a list of points in R^3 */
    void evaluate_shape_functions(const MATRIX &coords, 
                                  const INT_ARRAY &pointToEltMap,
                                  SPAR_MAT &N,
                                  SPAR_MAT_VEC &dN) const;   

    /** evaluate shape derivatives at a list of points in R^3 */
    void evaluate_shape_function_derivatives(const MATRIX &coords, 
                                             const INT_ARRAY &pointToEltMap,
                                             SPAR_MAT_VEC &dN) const;

    void shape_functions(const VECTOR &x, 
                         DENS_VEC &shp,
                         Array<int> &node_list) const
    { feMesh_->shape_functions(x,shp,node_list); }

    void shape_functions(const VECTOR & x,
                         DENS_VEC& shp, 
                         DENS_MAT& dshp,
                         Array<int> &node_list) const
    { feMesh_->shape_functions(x,shp,dshp,node_list); }

    void shape_functions(const VECTOR &x,
                         const int eltId,
                         DENS_VEC& shp,
                         Array<int> &node_list) const
    { feMesh_->shape_functions(x,eltId,shp,node_list); }

    void shape_functions(const VECTOR &x,
                         DENS_VEC& shp,
                         Array<int> &node_list,
                         int &eltId) const
    { feMesh_->shape_functions(x,shp,node_list,eltId); }

    void shape_functions(const VECTOR &x,
                         const int eltId,
                         DENS_VEC &shp, 
                         DENS_MAT &dshp,
                         Array<int> &node_list) const
    { feMesh_->shape_functions(x,eltId,shp,dshp,node_list); }
    /*@}*/

    //----------------------------------------------------------------
    /** \name kernel functions */
    //----------------------------------------------------------------
    /** evaluate kernel function */
    void evaluate_kernel_functions(const MATRIX &pt_coords,
                                   SPAR_MAT &N) const;

    /**  kernel matrix bandwidth */
    int kernel_matrix_bandwidth(const MATRIX &pt_coords) const;

    //----------------------------------------------------------------
    /** \name nodeset */
    //----------------------------------------------------------------
    /** pass through */
    void create_nodeset(const std::string &name, const std::set<int> &nodeset) 
    { feMesh_->create_nodeset(name,nodeset); }

    //----------------------------------------------------------------
    /** \name accessors */
    //----------------------------------------------------------------
    /*@{*/
    /** even though these are pass-throughs there is a necessary 
     *  translation */
    /** return number of unique nodes */
    int num_nodes() const { return feMesh_->num_nodes_unique(); }
  
    /** return number of total nodes */
    int nNodesTotal() const { return feMesh_->num_nodes(); }
  
    /** return number of elements */
    int num_elements() const { return feMesh_->num_elements(); }
    int my_num_elements() const { return feMesh_->my_num_elements(); }

    /** return number of nodes per element */
    int num_nodes_per_element() const { return feMesh_->num_nodes_per_element(); }

    /** return element connectivity */
    void element_connectivity(const int eltID, 
                              Array<int> & nodes) const
    { feMesh_->element_connectivity_unique(eltID, nodes); }
  
    /** return face connectivity */
    void face_connectivity(const PAIR &faceID, 
                           Array<int> &nodes) const
    { feMesh_->face_connectivity_unique(faceID, nodes); }

    /** in lieu of pass-throughs const accessors ... */
    /** return const ptr to mesh */
    const FE_Mesh* fe_mesh() const { return feMesh_; }
  
    /** return number of spatial dimensions */
    int nsd() const { return feMesh_->num_spatial_dimensions(); }
  
    /** return if the FE mesh has been created */
    int has_mesh() const { return feMesh_!=NULL; }
  
    /** get nodal coordinates for a given element */
    void element_coordinates(const int eltIdx, DENS_MAT &coords)
    { feMesh_->element_coordinates(eltIdx,coords); }

    /** get nodal coordinates for a given element */
    void element_field(const int eltIdx, const DENS_MAT field, 
                       DENS_MAT &local_field)
    { feMesh_->element_field(eltIdx, field, local_field); }

    /** access list of elements to be deleted */
    const std::set<int> &null_elements(void) const
    { return nullElements_; } 

    /** access to the amended nodal coordinate values */
    const DENS_MAT &nodal_coordinates(void) const 
    { return (*feMesh_->coordinates()); }

    /** map global node numbering to unique node numbering for 
     *  amended mesh */
    int map_global_to_unique(const int global_id) const 
    { return (*feMesh_->node_map())(global_id); }
    
    int number_of_global_nodes(void) const { return nNodes_; }

    /*@}*/

    /** set kernel */
    
    void set_kernel(KernelFunction* ptr);
    KernelFunction *kernel(int i=0) { return kernelFunction_; }

  private:
    //----------------------------------------------------------------
    /** mesh setup commands (called from modify) */
    //----------------------------------------------------------------
    /*@{*/

    MPI_Comm communicator_;

    /** finite element mesh */
    FE_Mesh *feMesh_;

    /** auxillary kernel function */
    KernelFunction *kernelFunction_;
 
    /** initialized flag */ 
    bool initialized_;

    /** create a uniform, structured mesh */
    void create_mesh(Array<double> &dx, 
                     Array<double> &dy, 
                     Array<double> &dz, 
                     const char *regionName,
                     Array<bool> periodic);

    void create_mesh(int nx, int ny, int nz, 
                     const char *regionName,
                     Array<bool> periodic);

    /** read an unstructured mesh from a file */
    void read_mesh(std::string meshFile, Array<bool> & periodicity);
    /*@}*/

    /** data that can be used for a subset of original mesh */
    std::set<int> nullElements_;
    
    /** faces upon which nodes are duplicated */
    std::set<PAIR> cutFaces_;
    std::set<int> cutEdge_;
  
    /** workspace */
    int nNodesPerElement_;
    int nSD_;
    int nElems_;
    int nNodes_;       /** number of global nodes */
    int nNodesUnique_; /** number of unique nodes */
    mutable int nIPsPerElement_;
    mutable int nIPsPerFace_;
    mutable FeIntQuadrature quadrature_;
    mutable FIELDS::const_iterator _fieldItr_; 
    mutable FieldName _fieldName_;
    
    /** sized arrays */
    mutable DIAG_MAT _weights_;
    mutable DENS_MAT _N_, _Nw_;
    mutable DENS_MAT_VEC _dN_, _dNw_;
    mutable DIAG_MAT _fweights_;
    mutable DENS_MAT _fN_;
    mutable DENS_MAT_VEC _fdN_, _nN_;

    /** unsized arrays */
    mutable DENS_MAT _Nmat_;
    mutable FIELD_MATS _fieldsAtIPs_;
    mutable GRAD_FIELD_MATS _gradFieldsAtIPs_; 
    mutable DENS_MAT _Nfluxes_;
    mutable AliasArray<int> _conn_;
    mutable DENS_MAT_VEC _Bfluxes_;
  
    /** output object */
    OutputManager outputManager_;
  
  };

}; // end namespace ATC

#endif

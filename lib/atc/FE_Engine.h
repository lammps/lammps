/** fe_engine :
   * computes and assembles mass matrix, rhs vectors
   * initial conditions handled in atc_transfer
*/

/** field structure:
   a field is a dense matrix of numPoints X number of DOF in the field
   a gradient is a std::vector of fields of length numSpatialDimensions
   a set of fields is a map of fieldName->field of length numFields
   a set of gradients is a map of fieldName->gradient of length numFields
   Note:  shape functions follow similar conventions with a shape function being
     a field of numPoints X numNodes 
     and a shape function gradient being a std::vector of shape functions
     of length numSpatialDimensions, although this is modified in calls when
     numPoints = 1
   Note:  the convention between shape function format and field format allows 
     for shape functions to matmat nodal fields, creating matrices of 
     numPoints X numElementsInField for evaluating data at atomic/quadarture points
*/

/** current internal limitations:
   * 3 spatial dimensions
   * 8 node bricks
   * structured mesh
   * no stiffness matrix 
       (i.e. no implicit integration or special treatment of linear problems)
   * lumped mass
*/



/** terminology:
density rate = flux
             = Grad.FLUX(f,D_x f) 
             + SOURCE(f,D_x f) + PRESCRIBED_SOURCE(x,t) 
             + EXTRINSIC_SOURCE(f,D_x f,f_e,D_x f_e)
*/

#ifndef FE_ENGINE_H
#define FE_ENGINE_H

// Other headers
#include <vector>
#include <map>

// ATC_Transfer headers
#include "Array.h"
#include "Array2D.h"
#include "FE_Mesh.h"
#include "PhysicsModel.h"
#include "OutputManager.h"

using namespace std;

namespace ATC {
  
  // Forward declarations
  class ATC_Transfer;
  class FE_Element;
  class XT_Function;
 
  class FE_Engine{
  public:
    /** constructor/s  */
    FE_Engine(ATC_Transfer * atcTransfer);
  
    /** destructor */
    ~FE_Engine();
  
    /** initialize */
    void initialize();
  
    /** parser/modifier */
    bool modify(int narg, char **arg);
  
    /** finish up */
    void finish();

    //----------------------------------------------------------------
    /** \name output */
    //----------------------------------------------------------------
    /*@{*/
    /** these assume the caller is handling the parallel collection */
    void initialize_output(int rank,
      string outputPrefix, OutputType otype = ENSIGHT);
  
    /** write geometry */
    void write_geometry(void);
  
    /** write data: data is arrayed over _unique_ nodes 
                  & then mapped by the engine */
    void write_data(double time, FIELDS &soln, OUTPUT_LIST *data=NULL);
    void write_data(double time, OUTPUT_LIST *data);
    void write_restart_file(string fileName, OUTPUT_LIST *data)
    {outputManager_.write_restart_file(fileName,data);};
    void read_restart_file(string fileName, OUTPUT_LIST *data)
    {outputManager_.read_restart_file(fileName,data);};

    void delete_elements(const set<int> & elementList);

    void add_global(string name, double value) 
      {outputManager_.add_global(name,value);}
    void reset_globals() {outputManager_.reset_globals();}

    /** pass through to access output manager */
    OutputManager* output_manager() {return &outputManager_;}
    /*@}*/

    //----------------------------------------------------------------
    /** \name assembled matrices and vectors */
    //----------------------------------------------------------------
    /*@{*/
    /** compute a dimensionless stiffness matrix */
    void compute_stiffness_matrix(SPAR_MAT &matrix) const;

    /** compute a dimensionless mass matrix */
    void compute_mass_matrix(SPAR_MAT &mass_matrix) const;
    /** computes a dimensionless mass matrix for the given-quadrature */
    void compute_mass_matrix(const DIAG_MAT &weights,
                             const SPAR_MAT &N,
                             SPAR_MAT &mass_matrix) const;
    /** compute a single dimensionless mass matrix */
    void compute_lumped_mass_matrix(DIAG_MAT &lumped_mass_matrix) const;

    /** compute lumped mass matrix = diag (\int \rho N_I dV) */
    void compute_lumped_mass_matrix(const Array<FieldName>   &mask, 
                                    const FIELDS &fields,
                                    const PhysicsModel * physicsModel,
                                    const Array<int>   & elementMaterials,
                                    map<FieldName, DIAG_MAT> &mass_matrix,
                                    const Array<bool> *element_mask=NULL) const;
    /** compute dimensional lumped mass matrix using given quadrature */
    void compute_lumped_mass_matrix(const Array<FieldName>   &mask, 
                                    const FIELDS &fields,
                                    const PhysicsModel * physicsModel,
                                    const Array<set<int> > & pointMaterialGroups,
                                    const DIAG_MAT      &weights,
                                    const SPAR_MAT      &N,
                                    map<FieldName, DIAG_MAT> &mass_matrix) const;

    /** compute an approximation to a finite difference gradient from mesh */
    void compute_gradient_matrix(GRAD_SHPFCN &grad_matrix) const;


    /** compute energy */
    void compute_energy(const Array<FieldName> &mask,
                        const FIELDS &fields,
                        const PhysicsModel * physicsModel,
                        const Array<int>   & elementMaterials,
                        FIELDS &energy, 
                        const Array<bool> *element_mask=NULL) const;


    /** compute residual or RHS of the dynamic weak eqn */
    void compute_rhs_vector(const Array2D<bool> &rhs_mask,
                            const FIELDS &fields,
                            const PhysicsModel * physicsModel,
                            const Array<int>   & elementMaterials,
                            FIELDS &rhs, 
                            const Array<bool> *element_mask=NULL) const;

    /** compute RHS for given quadrature */
    void compute_rhs_vector(const Array2D<bool> &rhs_mask,
                            const FIELDS        &fields,
                            const PhysicsModel * physicsModel,
                            const Array<set<int> > & pointMaterialGroups,
                            const DIAG_MAT      &weights,
                            const SPAR_MAT      &N,
                            const GRAD_SHPFCN   &dN,
                            FIELDS              &rhs) const;

    /** compute flux in domain i.e. N^T B_integrand */
    void compute_flux(const Array2D<bool> & rhs_mask,
                        const FIELDS &fields,
                        const PhysicsModel * physicsModel,
                        const Array<int>   & elementMaterials,
                        GRAD_FIELDS &flux, 
                        const Array<bool> *element_mask=NULL) const;

    /** compute the flux on the MD/FE boundary */
    void compute_boundary_flux(
      const Array2D<bool> & rhs_mask, 
      const FIELDS        & fields,
      const PhysicsModel  * physicsModel,
      const Array<int>    & elementMaterials,
      const set<PAIR>     & faceSet, 
      FIELDS              & rhs) const;

    /** compute the flux on using an L2 interpolation of the flux */
    void compute_boundary_flux( 
      const Array2D<bool>     & rhs_mask,
      const FIELDS            & fields,
      const PhysicsModel      * physicsModel,
      const Array<int>        & elementMaterials,
      const Array<set<int> >  & pointMaterialGroups,
      const DIAG_MAT          & weights,
      const SPAR_MAT          & N,
      const GRAD_SHPFCN       & dN,
      const DIAG_MAT          & flux_mask,
      FIELDS                  & rhs ) const;

    /** compute prescribed flux given an array of functions of x & t */
    void  add_fluxes(const Array<bool> &fieldMask, 
                     const double time,
                     const SURFACE_SOURCE & sourceFunctions, 
                     FIELDS &nodalSources) const;

    /** compute nodal vector of volume based sources */
    void add_sources(const Array<bool> &fieldMask, 
                     const double time,
                     const VOLUME_SOURCE &sourceFunctions, 
                     FIELDS &nodalSources) const;

    /** compute surface flux of a nodal field */
    void field_surface_flux(const DENS_MAT & field,
                            const set<PAIR> &faceSet,
                            DENS_MAT & values,
                            const bool contour = false,
                            const int axis = 2) const;

    /*@}*/

    //----------------------------------------------------------------
    /** \name shape functions */
    //----------------------------------------------------------------
    /*@{*/
    /** evaluate shape function at a list of points in R^3 */
    void evaluate_shape_functions(const MATRIX &coords, 
                                  SPAR_MAT &N,
                                  Array<int> & pointToEltMap) const;
  
    /** evaluate shape function & derivatives at a list of points in R^3 */
    void evaluate_shape_functions( const MATRIX &coords, 
                                   SPAR_MAT &N,
                                   GRAD_SHPFCN &dN,
                                   Array<int> & pointToEltMap) const;
  
    /** evaluate all shape function & derivatives at a specific R^3 location */ 
    void evaluate_shape_functions(const VECTOR & x,
                                  Array<int>& node_index,
                                  DENS_VEC& shp, 
                                  DENS_MAT& dshp,
                                  int & eltID) const;

    /** pass through */
    void shape_functions(const VECTOR &x, 
                         DENS_VEC& shp, 
                         int & eltID,
                         Array<int>& node_list) const
    { feMesh_->shape_functions(x,shp,eltID,node_list); }

    void shape_functions(const VECTOR &x, 
                         DENS_VEC& shp, 
                         int & eltID,
                         Array<int>& node_list,
                         const Array<bool>& periodicity) const
    { feMesh_->shape_functions(x,shp,eltID,node_list, periodicity); }
    /*@}*/

    //----------------------------------------------------------------
    /** \name accessors */
    //----------------------------------------------------------------
    /*@{*/
    /** even though these are pass-throughs there is a necessary translation */
    /** return number of unique nodes */
    int get_nNodes() const { return feMesh_->get_nNodesUnique(); };
  
    /** return number of total nodes */
    int get_nNodesTotal() const { return feMesh_->get_nNodes(); };
  
    /** return number of elements */
    int get_nElements() const { return feMesh_->get_nElements(); };

    /** return element connectivity */
    void element_connectivity(const int eltID, Array<int> & nodes) const
    { feMesh_->element_connectivity_unique(eltID, nodes); }
  
    /** return face connectivity */
    void face_connectivity(const PAIR &faceID, Array<int> &nodes) const
    {  feMesh_->face_connectivity_unique(faceID, nodes); }

    /** in lieu of pass-throughs const accessors ... */
    // return const ptr to mesh
    const FE_Mesh* get_feMesh() const { return feMesh_;}
  
    // return number of spatial dimensions
    int get_nsd() const { return feMesh_->get_nSpatialDimensions(); }
  
    // return if the FE mesh has been created
    int fe_mesh_exist() const { return feMesh_!=NULL; }
  
    // get nodal coordinates for a given element
    void element_coordinates(const int eltIdx, DENS_MAT &coords)
    { feMesh_->element_coordinates(eltIdx,coords); }

    // access list of elements to be deleted
    set<int> & null_elements(void) { return nullElements_; }
    /*@}*/

  private:
    //----------------------------------------------------------------
    /** mesh setup commands (called from modify) */
    //----------------------------------------------------------------
    /*@{*/
 
    /** initialized flag */ 
    bool initialized_;

    /** create a uniform, structured mesh */
    void create_mesh(int nx, int ny, int nz, char * regionName,
                     int xperiodic, int yperiodic, int zperiodic);
    /*@}*/

    /** ATC transfer object */
    ATC_Transfer * atcTransfer_;
  
    /** finite element mesh */
    FE_Mesh * feMesh_;

    /** data that can be used for a subset of original mesh */
    set<int> nullElements_;
    bool amendedMeshData_;
    const Array2D<int> * connectivity_;
    const Array<int>   * nodeMap_;
    const DENS_MAT     * coordinates_;
  
    /** workspace */
    int nNodesPerElement_;
    int nIPsPerElement_;
    int nIPsPerFace_;
    int nSD_;
    int nElems_;
  
    /** output object */
    OutputManager outputManager_;
  
    /** base name for output files */
    string outputPrefix_;
  
    /** output frequency (NOTE will move to "Transfer") */
    int outputFrequency_;
  
    /** list of output timesteps */
    vector<double> outputTimes_;
  };
}; // end namespace ATC

#endif

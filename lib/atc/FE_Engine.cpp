#include "FE_Engine.h"
#include "ATC_Transfer.h"
#include "FE_Element.h"
#include "XT_Function.h"
#include "LammpsInterface.h"
#include "PhysicsModel.h"

#include <stdio.h>
#include <map>

using namespace std;

namespace ATC{

  //-----------------------------------------------------------------
  FE_Engine::FE_Engine(ATC_Transfer * atcTransfer)
    : atcTransfer_(atcTransfer),
      feMesh_(NULL),
      initialized_(false),
      amendedMeshData_(false),
      outputManager_()
  {
  }

  FE_Engine::~FE_Engine()
  {
    if (feMesh_)       delete feMesh_;
    if (connectivity_ && amendedMeshData_) delete connectivity_;
    if (nodeMap_      && amendedMeshData_) delete nodeMap_;
    if (coordinates_  && amendedMeshData_) delete coordinates_;
  }

  //-----------------------------------------------------------------
  void FE_Engine::initialize()
  {
    if (! feMesh_) 
      throw ATC_Error(0,"FE_Engine has no mesh");

    if (! initialized_ ) {
      // mesh data
      nodeMap_      = &(feMesh_->global_to_unique_map());
      connectivity_ = &(feMesh_->connectivity());
      coordinates_  = &(feMesh_->nodal_coordinates());

      // set up work spaces
      nNodesPerElement_ = feMesh_->get_nNodesPerElement();
      nIPsPerElement_ = feMesh_->get_nIPsPerElement();
      nIPsPerFace_ = feMesh_->get_nIPsPerFace();
      nSD_ = feMesh_->get_nSpatialDimensions();
      nElems_ = feMesh_->get_nElements();

      // remove specified elements
      if (nullElements_.size() > 0) delete_elements(nullElements_);
 
      initialized_ = true;
    }
  }

  //-----------------------------------------------------------------
  bool FE_Engine::modify(int narg, char **arg)
  {
    bool match = false;
    /*! \page man_fem_mesh fix_modify AtC fem create mesh
      \section syntax
      fix_modify AtC fem create mesh <nx> <ny> <nz> <region-id> 
      <f|p> <f|p> <f|p>
      - nx ny nz = number of elements in x, y, z
      - region-id = id of region that is to be meshed
      - f p p  = perioidicity flags for x, y, z
      \section examples
      <TT> fix_modify AtC fem  create mesh 10 1  1  feRegion p p p </TT>
      \section description
      Creates a uniform mesh in a rectangular region
      \section restrictions
      creates only uniform rectangular grids in a rectangular region
      \section related
      \section default
      none
    */
    if (strcmp(arg[0],"fem")==0) {
      // create mesh
      if (strcmp(arg[1],"create")==0) {
        if (strcmp(arg[2],"mesh")==0) {
          int nx = atoi(arg[3]);
          int ny = atoi(arg[4]);
          int nz = atoi(arg[5]);
          int xperiodic = 0;
          if (strcmp(arg[7],"p")==0) xperiodic = 1;
          int yperiodic = 0;
          if (strcmp(arg[8],"p")==0) yperiodic = 1;
          int zperiodic = 0;
          if (strcmp(arg[9],"p")==0) zperiodic = 1;
          if (feMesh_) throw ATC_Error(0,"FE Engine already has a mesh");
          create_mesh(nx, ny, nz, arg[6], xperiodic, yperiodic, zperiodic );

          // construct prescribed data manager NOTE move to ATC_Transfer?
          atcTransfer_->initialize_mesh_data();

          match = true;
        }
      }
    }
    /*! \page man_mesh_delete_elements fix_modify AtC mesh delete_elements
      \section syntax
      fix_modify AtC mesh delete_elements <element_set>
      - <element_set> = name of an element set
      \section examples
      <TT> fix_modify AtC delete_elements gap </TT>
      \section description
      Deletes a group of elements from the mesh.
      \section restrictions
      \section related
      \section default
      none
    */
    else if ( (strcmp(arg[0],"mesh")==0)  
           && (strcmp(arg[1],"delete_elements")==0) ) {
      string esetName = arg[2];
      set<int> elemSet = feMesh_->get_elementset(esetName);
      nullElements_.insert(elemSet.begin(), elemSet.end());
      match = true;
    }
    // FE_Mesh
    else {
      match = feMesh_->modify(narg,arg);
    }
    return match;
  }

  //-----------------------------------------------------------------
  void FE_Engine::finish()
  {
    // Nothing to do
  }

  //-----------------------------------------------------------------
  //  write geometry
  //-----------------------------------------------------------------
  void FE_Engine::initialize_output(int rank,
    string outputPrefix, OutputType otype)
  {
    outputManager_.initialize(outputPrefix, otype);
    if (! feMesh_) throw ATC_Error(0,"output needs mesh");
    if (! initialized_) initialize();
    if (rank == 0) outputManager_.write_geometry(*coordinates_, connectivity_);
  }

  //-----------------------------------------------------------------
  //  write geometry
  //-----------------------------------------------------------------
  void FE_Engine::write_geometry(void)
  {
    outputManager_.write_geometry(*coordinates_, connectivity_);
  }

  // -------------------------------------------------------------
  //  write data  
  // -------------------------------------------------------------
  void FE_Engine::write_data(double time, FIELDS &soln, OUTPUT_LIST *data)
  {
    outputManager_.write_data(time, &soln, data, nodeMap_->get_data());
  }

  // -------------------------------------------------------------
  //  write data  
  // -------------------------------------------------------------
  void FE_Engine::write_data(double time, OUTPUT_LIST *data)
  {
    outputManager_.write_data(time, data, nodeMap_->get_data());
  }

  // -------------------------------------------------------------
  //  amend mesh for deleted elements
  // -------------------------------------------------------------
  void FE_Engine::delete_elements(const set<int> & elementList)
  {
    int nsd = feMesh_->get_nSpatialDimensions();
    set<int> elementsNew;
    feMesh_->elementset_complement(elementList,elementsNew);
    int nElementsNew = elementsNew.size();
    set<int> newToOld;
    map<int,int> oldToNewMap;
    feMesh_->elementset_to_nodeset(elementsNew,newToOld);
    int nNodesNew = newToOld.size(); 
    set<int>::const_iterator itr;


    // coordinates & node map (from nodes to data)
    const DENS_MAT &coordinates = feMesh_->nodal_coordinates();
    const Array<int> & node_map = feMesh_->global_to_unique_map();
    if (coordinates_ && amendedMeshData_) delete coordinates_;
    coordinates_ = new DENS_MAT(nsd,nNodesNew);
    DENS_MAT & coor = * (const_cast<DENS_MAT *> ( coordinates_));
    if (nodeMap_ && amendedMeshData_) delete nodeMap_;
    nodeMap_     = new Array<int> (nNodesNew);
    Array<int> & nmap = * (const_cast<Array<int> *> ( nodeMap_));
    int k = 0, i = 0;
    for (itr = newToOld.begin(); itr != newToOld.end(); itr++) {
      int node = *itr;
      oldToNewMap[node]=i++;
      nmap(k) = node_map(node);
      for(int j = 0; j < nsd; j++) {
        coor(j,k) = coordinates(j,node);
      }
      k++;
    }
    // connectivity
    const Array2D<int> & connectivity = feMesh_->connectivity();
    int nNodesPerElement = connectivity.nRows();
    int nElements = connectivity.nCols();
    if (connectivity_ && amendedMeshData_) delete connectivity_;
    connectivity_ = new Array2D<int>(nNodesPerElement,nElementsNew);
    Array2D<int> & conn = * (const_cast<Array2D<int> *> ( connectivity_));
    k = 0;
    for (itr = elementsNew.begin(); itr != elementsNew.end(); itr++) {
      int ielem = *itr;
      for(int j = 0; j < nNodesPerElement; j++) {
        int old_node = connectivity(j,ielem);
        map<int,int>::iterator map_itr = oldToNewMap.find(old_node);
        if (map_itr == oldToNewMap.end()) { 
          cout << "map failure " << old_node << "\n"; 
        }
        int node = map_itr->second;
        conn(j,k) = node;
      }
      k++;
    }
    amendedMeshData_ = true;
  }

  // -------------------------------------------------------------
  // compute dimensionless stiffness matrix using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_stiffness_matrix(SPAR_MAT &matrix) const
  {
    //int nfields = field_mask.get_length(); 
    int nNodes = feMesh_->get_nNodesUnique();
    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nIPsPerElement = feMesh_->get_nIPsPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();

    DENS_MAT N(nIPsPerElement,nNodesPerElement);
    vector<DENS_MAT> dN(nsd);
    dN.assign(nsd, DENS_MAT(nIPsPerElement,nNodesPerElement));
    DiagonalMatrix<double> weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);


    // assemble consistent mass (nnodes X nnodes)
    matrix.reset(nNodes,nNodes);
    DENS_MAT elementMassMatrix(nNodesPerElement,nNodesPerElement);
    for (int ielem = 0; ielem < nelems; ++ielem) 
    {  
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, dN, weights);
      // perform quadrature
      elementMassMatrix = dN[0].transMat(weights*dN[0]); 
      for (int i = 1; i < nsd; ++i) {
        elementMassMatrix += dN[i].transMat(weights*dN[i]); 
      }
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);
      for (int i = 0; i < nNodesPerElement; ++i) 
      {
        int inode = conn(i);
        for (int j = 0; j < nNodesPerElement; ++j) 
        {
          int jnode = conn(j);
          matrix.add(inode, jnode, elementMassMatrix(i,j));
        }
      }
    }
  }
  // -------------------------------------------------------------
  // compute dimensionless consistent mass using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_mass_matrix(SPAR_MAT &mass_matrix) const
  {
    int nNodes = feMesh_->get_nNodesUnique();
    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nIPsPerElement = feMesh_->get_nIPsPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();

    DENS_MAT N(nIPsPerElement,nNodesPerElement);
    DiagonalMatrix<double> weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);

    // assemble consistent mass (nnodes X nnodes)
    mass_matrix.reset(nNodes,nNodes);
    DENS_MAT elementMassMatrix(nNodesPerElement,nNodesPerElement);
    for (int ielem = 0; ielem < nelems; ++ielem) 
    {  
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, weights);
      // perform quadrature
      elementMassMatrix = N.transMat(weights*N); 
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);
      for (int i = 0; i < nNodesPerElement; ++i) 
      {
        int inode = conn(i);
        for (int j = 0; j < nNodesPerElement; ++j) 
        {
          int jnode = conn(j);
          mass_matrix.add(inode, jnode, elementMassMatrix(i,j));
        }
      }
    }
  }

  // -------------------------------------------------------------
  // compute consistent mass using given quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_mass_matrix(const DIAG_MAT &weights,
                                      const SPAR_MAT &N,
                                      SPAR_MAT &mass_matrix) const
  {
    int nn = N.nCols();
    int nips = N.nRows();
    
    DENS_MAT tmp_mass_matrix_local(nn,nn), tmp_mass_matrix(nn,nn);
    if (nips>0) { tmp_mass_matrix_local = N.transMat(weights*N); } 
    // share information between processors
    int count = nn*nn; 
    LammpsInterface::instance()->allsum(
      tmp_mass_matrix_local.get_ptr(),
      tmp_mass_matrix.get_ptr(), count);

    // create sparse from dense
    mass_matrix.reset(tmp_mass_matrix);

  }

  // -------------------------------------------------------------
  // compute dimensionless lumped mass using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_lumped_mass_matrix(DiagonalMatrix<double> & M) const
  {
    // initialize
    int nNodes = feMesh_->get_nNodesUnique();
    M.reset(nNodes,nNodes); 

    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nIPsPerElement = feMesh_->get_nIPsPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();

    DENS_MAT N(nNodesPerElement,nIPsPerElement);
    DIAG_MAT weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);

    // assemble lumped diagonal mass (nnodes X nnodes)
    int inode = -1;
    DENS_VEC weightVec(nIPsPerElement);
    DENS_VEC Nvec(nNodesPerElement);
    for (int ielem = 0; ielem < nelems; ++ielem) 
    {  
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, weights);
      int nips = weights.nRows();  // different than nIPsPerElement?
      weightVec = 1.;
      weightVec = weights*weightVec;
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);
      Nvec = N*weightVec;
      for (int i = 0; i < nNodesPerElement; ++i) 
      {
        inode = conn(i);
        M(inode,inode) += Nvec(i);
      }
    }
  }

  // -------------------------------------------------------------
  // compute physical lumped mass using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_lumped_mass_matrix(
    const Array<FieldName>& field_mask,
    const FIELDS          & fields,
    const PhysicsModel * physicsModel,
    const Array<int> &elementMaterials,
    map<FieldName, DIAG_MAT>& M, // mass matrix
    const Array<bool> *element_mask) const
  {
    int nfields = field_mask.get_length();
    // initialize
    int nNodes = feMesh_->get_nNodesUnique();
    for (int j = 0; j < nfields; ++j) { 
      M[field_mask(j)].reset(nNodes,nNodes); 
    }

    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nIPsPerElement = feMesh_->get_nIPsPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();

    DENS_MAT N(nNodesPerElement,nIPsPerElement);
    DIAG_MAT weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);
    FIELDS fieldsAtIPs, localElementFields;
    DENS_MAT Nmat;

    // assemble lumped diagonal mass (nnodes X nnodes)
    for (int ielem = 0; ielem < nelems; ++ielem)
    {
      // if element is masked, skip it
      if (element_mask && !(*element_mask)(ielem)) continue;
      // material id
      int imat = elementMaterials(ielem);
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, weights);
      int nips = weights.nRows();  
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);

      // interpolate fields at IPs
      FIELDS::const_iterator field;
      for (field = fields.begin(); field != fields.end(); field++) 
      {
        // field values at all nodes
        const DENS_MAT &vI = field->second;
        // field values at element nodes
        DENS_MAT &vIe = localElementFields[field->first];
        // field values at integration points -> to be computed
        DENS_MAT &vP = fieldsAtIPs[field->first];

        int numFieldDOF = vI.nCols();
        vIe.reset(nNodesPerElement, numFieldDOF);
        // gather local field
        for (int i = 0; i < nNodesPerElement; i++) {
          for (int j = 0; j < numFieldDOF; j++) {
            vIe(i,j) = vI(conn(i),j);
          }
        }
          
        // interpolate field at integration points
        vP = N*vIe;
      }

      // compute energy/momentum/mass densities
      FIELDS capacities;
      physicsModel->M_integrand(field_mask, fieldsAtIPs, capacities, imat);

      // integrate & assemble
      for (int j = 0; j < nfields; ++j) {
        FieldName thisFieldName = field_mask(j);
        Nmat = N.transMat(weights*capacities[thisFieldName]);
        for (int i = 0; i < nNodesPerElement; ++i) {  
          int inode = conn(i);  
          M[thisFieldName](inode,inode) += Nmat(i,0);// assume all dof same
        }
      }
    }
  }

  // -------------------------------------------------------------
  // compute physical lumped mass using given quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_lumped_mass_matrix(
    const Array<FieldName>   &field_mask, 
    const FIELDS &fields,
    const PhysicsModel * physicsModel,
    const Array<set<int> > & pointMaterialGroups,
    const DIAG_MAT      &weights, // NOTE use these as a mask
    const SPAR_MAT      &N,
    map<FieldName, DIAG_MAT> &M) const // mass matrices 
  {
    int nips = weights.nCols();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nfields = field_mask.get_length();
    int nNodes = feMesh_->get_nNodesUnique();
    // initialize
    FieldName thisFieldName;
    map<FieldName, DIAG_MAT> M_local;
    for (int j = 0; j < nfields; ++j) {
      thisFieldName = field_mask(j);
      M_local[thisFieldName].reset(nNodes,nNodes);
      M[thisFieldName].reset(nNodes,nNodes);
    }
    
    if (nips>0) {
      // compute fields at given ips
      FIELDS fieldsAtIPs;
      for (int j = 0; j < nfields; ++j) {
        thisFieldName = field_mask(j);
        const FIELD & field = (fields.find(thisFieldName))->second;
        int numFieldDOF = field.nCols();
        fieldsAtIPs[thisFieldName].reset(nips,numFieldDOF);
        fieldsAtIPs[thisFieldName] = N*field;
      }
    
      // treat single material point sets specially
      int nMatls = pointMaterialGroups.size();
      int atomMatls = 0;
      for (int imat = 0; imat < nMatls; imat++) {
        const set<int> & indices = pointMaterialGroups(imat);
        if ( indices.size() > 0) atomMatls++;
      }
      bool singleMaterial = ( atomMatls == 1 );
    
      // setup data structures
      FIELDS capacities;
      // evaluate physics model
      if (singleMaterial) {
        physicsModel->M_integrand(field_mask, fieldsAtIPs, capacities);
      }
      else {
        throw ATC_Error(0,"unimplemented function in FE_Engine::compute_mass_matrix");
      }
    
      // integrate & assemble
      for (int j = 0; j < nfields; ++j) {
        FieldName thisFieldName = field_mask(j);
        M_local[thisFieldName].reset( // assume all columns same
          column(N.transMat(weights*capacities[thisFieldName]),0) );
      }
    }
    
    // Share information between processors
    for (int j = 0; j < nfields; ++j) {
      FieldName thisFieldName = field_mask(j);
      int count = M_local[thisFieldName].size();
      LammpsInterface::instance()->allsum(
        M_local[thisFieldName].get_ptr(),
        M[thisFieldName].get_ptr(), count);
    }
  }

  //-----------------------------------------------------------------
  // compute assembled average gradient evaluated at the nodes
  //-----------------------------------------------------------------
  void FE_Engine::compute_gradient_matrix( 
    GRAD_SHPFCN & grad_matrix) const
  {
    int nNodesUnique = feMesh_->get_nNodesUnique();
    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();
    int nIPsPerElement = feMesh_->get_nIPsPerElement();

    DENS_MAT N(nIPsPerElement,nNodesPerElement);
    vector< DENS_MAT > dN(nsd);
    vector< DENS_MAT > tmp_grad_matrix(nsd);
    for (int i = 0; i < nsd; i++) {
      dN[i].reset(nNodesPerElement,nNodesPerElement);
      tmp_grad_matrix[i].reset(nNodesUnique,nNodesUnique);
    }
  
    DiagonalMatrix<double> weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);
  
    // element wise assembly
    Array<int> count(nNodesUnique); count = 0;
    feMesh_->set_quadrature(FE_Element::NODAL_QUADRATURE);
    for (int ielem = 0; ielem < nelems; ++ielem) {
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, dN, weights);
      int nips = weights.nRows();
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);

      for (int j = 0; j < nIPsPerElement; ++j) {
        int jnode = conn(j); // NOTE assumes ip order is node order
        count(jnode) += 1;
        for (int i = 0; i < nNodesPerElement; ++i) {
          int inode = conn(i);
          for (int k = 0; k < nsd; ++k) {
            tmp_grad_matrix[k](jnode,inode) += dN[k](j,i);
          }
        }
      }
    }
    feMesh_->set_quadrature(FE_Element::GAUSSIAN_QUADRATURE); // reset to default
    for (int inode = 0; inode < nNodesUnique; ++inode) {
      //cout << inode << " count " << count(inode) << "\n";
      for (int jnode = 0; jnode < nNodesUnique; ++jnode) {
        for (int k = 0; k < nsd; ++k) {
          tmp_grad_matrix[k](jnode,inode) /= count(jnode);
        }
      }
    }
    // compact dense matrices
    grad_matrix.resize(nsd);
    for (int k = 0; k < nsd; ++k) {
      grad_matrix[k].reset(tmp_grad_matrix[k]);
    }
  }

  // -------------------------------------------------------------
  // compute energy per node using native quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_energy(const Array<FieldName> &mask,
                                     const FIELDS &fields,
                                     const PhysicsModel * physicsModel,
                                     const Array<int>   & elementMaterials,
                                     FIELDS &energy,
                                     const Array<bool> *element_mask) const
  {
    //NOTE move to initialization and make workspace
    const int nNodesPerElement = feMesh_->get_nNodesPerElement();
    const int nIPsPerElement   = feMesh_->get_nIPsPerElement();
    const int nsd              = feMesh_->get_nSpatialDimensions();
    const int nelems           = feMesh_->get_nElements();
    DENS_MAT Nmat;

    //NOTE move to initialization and make workspace
    DENS_MAT N(nIPsPerElement,nNodesPerElement);
    vector<DENS_MAT> dN(nsd);
    dN.assign(nsd, DENS_MAT(nIPsPerElement,nNodesPerElement));
  
    FIELDS fieldsAtIPs, localElementFields;
    GRAD_FIELDS grad_fieldsAtIPs;
    FIELDS::const_iterator field;
    FieldName thisFieldName;
    int numFieldDOF;

    for (int n = 0; n < mask.get_length(); n++) {
      FieldName thisFieldName = mask(n);
      const FIELD & field = (fields.find(thisFieldName))->second;
      energy[thisFieldName].reset(field.nRows(), field.nCols());
    }

    DIAG_MAT weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);

    // element wise assembly
    int inode = -1;
    for (int ielem=0; ielem < nelems; ielem++) 
    {
      // if element is masked, skip it
      if (element_mask && !(*element_mask)(ielem)) continue;
      // material id
      int imat = elementMaterials(ielem);
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, dN, weights);
      int nips = weights.nRows();
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);
      // compute fields and gradients of fields ips of this element

      for (int n = 0; n < mask.get_length(); n++) {
        FieldName thisFieldName = mask(n);
        FIELDS::const_iterator field = (fields.find(thisFieldName));
        // field values at all nodes
        const DENS_MAT &vI = field->second;
        // field values at element nodes
        DENS_MAT &vIe = localElementFields[field->first];
        // field values at integration points -> to be computed
        DENS_MAT &vP = fieldsAtIPs[field->first];
        // gradients of field at integration points -> to be computed
        vector<DENS_MAT> &dvP = grad_fieldsAtIPs[field->first];

        numFieldDOF = vI.nCols();
        vIe.reset(nNodesPerElement, numFieldDOF);
        // gather local field
        for (int i = 0; i < nNodesPerElement; i++)
          for (int j = 0; j < numFieldDOF; j++)
            vIe(i,j) = vI(conn(i),j);
          
        // interpolate field at integration points
        vP = N*vIe;
        dvP.assign(nsd, DENS_MAT(nIPsPerElement, numFieldDOF));
        for (int j = 0; j < nsd; ++j) dvP[j] = dN[j]*vIe;
      }

      FIELDS energies;
      physicsModel->E_integrand(mask, fieldsAtIPs, grad_fieldsAtIPs, energies, imat);
      // assemble
      for (int n = 0; n < mask.get_length(); n++) {
        FieldName thisFieldName = mask(n);
        FIELDS::const_iterator field = (fields.find(thisFieldName));
        numFieldDOF = (field->second).nCols();
        Nmat = N.transMat(weights*energies[thisFieldName]);
        for (int i = 0; i < nNodesPerElement; ++i) {  
          inode = conn(i);  
          for (int k = 0; k < numFieldDOF; ++k) {
            energy[thisFieldName](inode,k) += Nmat(i,k);
          }
        }
      }
    } // element loop
  }  // function call

  // -------------------------------------------------------------
  // compute rhs using native quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_rhs_vector(const Array2D<bool> &rhs_mask,
                                     const FIELDS &fields,
                                     const PhysicsModel * physicsModel,
                                     const Array<int>   & elementMaterials,
                                     FIELDS &rhs,
                                     const Array<bool> *element_mask) const
  {
    //NOTE move to initialization and make workspace
    const int nNodesPerElement = feMesh_->get_nNodesPerElement();
    const int nIPsPerElement   = feMesh_->get_nIPsPerElement();
    const int nsd              = feMesh_->get_nSpatialDimensions();
    const int nelems           = feMesh_->get_nElements();
    DENS_MAT Nmat;

    //NOTE move to initialization and make workspace
    DENS_MAT N(nIPsPerElement,nNodesPerElement);
    vector<DENS_MAT> dN(nsd);
    dN.assign(nsd, DENS_MAT(nIPsPerElement,nNodesPerElement));
  
    FIELDS fieldsAtIPs, localElementFields;
    GRAD_FIELDS grad_fieldsAtIPs;
    FIELDS::const_iterator field;
    FieldName thisFieldName;
    int numFieldDOF;

    for (field = fields.begin(); field != fields.end(); field++) {
      thisFieldName = field->first;
      if (rhs_mask(thisFieldName,FLUX) || rhs_mask(thisFieldName,SOURCE)) {
        int nrows = (field->second).nRows();
        int ncols = (field->second).nCols();
        rhs[thisFieldName].reset(nrows,ncols);
      }
    }

    DIAG_MAT weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);

    // element wise assembly
    int inode = -1;
    for (int ielem=0; ielem < nelems; ielem++) 
    {
      // if element is masked, skip it
      if (element_mask && !(*element_mask)(ielem)) continue;
      // material id
      int imat = elementMaterials(ielem);
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, dN, weights);
      int nips = weights.nRows();
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);
      // compute fields and gradients of fields ips of this element

      for (field = fields.begin(); field != fields.end(); field++) 
      {
        // field values at all nodes
        const DENS_MAT &vI = field->second;
        // field values at element nodes
        DENS_MAT &vIe = localElementFields[field->first];
        // field values at integration points -> to be computed
        DENS_MAT &vP = fieldsAtIPs[field->first];
        // gradients of field at integration points -> to be computed
        vector<DENS_MAT> &dvP = grad_fieldsAtIPs[field->first];

        numFieldDOF = vI.nCols();
        vIe.reset(nNodesPerElement, numFieldDOF);
        // gather local field
        for (int i = 0; i < nNodesPerElement; i++)
          for (int j = 0; j < numFieldDOF; j++)
            vIe(i,j) = vI(conn(i),j);
          
        // interpolate field at integration points
        vP = N*vIe;
        dvP.assign(nsd, DENS_MAT(nIPsPerElement, numFieldDOF));
        for (int j = 0; j < nsd; ++j) dvP[j] = dN[j]*vIe;
      }

      // NOTE move scaling of N by weights up here

      // evaluate physics model
      // N_fluxes is a set of fields
      if(physicsModel->has_N_integrand()) {
        FIELDS N_fluxes;
        physicsModel->N_integrand(rhs_mask, fieldsAtIPs, grad_fieldsAtIPs, N_fluxes, imat);
        // assemble
        for (field = fields.begin(); field != fields.end(); field++) 
        {
          thisFieldName = field->first;
          // NOTE why is this line here?
          //      SOURCE should refer only to adding extrinsic coupling sources
          if (!rhs_mask(thisFieldName,SOURCE)) continue;
          numFieldDOF = (field->second).nCols();
          Nmat = N.transMat(weights*N_fluxes[thisFieldName]);
          for (int i = 0; i < nNodesPerElement; ++i) {  
            inode = conn(i);  
            for (int k = 0; k < numFieldDOF; ++k) {
              rhs[thisFieldName](inode,k) += Nmat(i,k);
            }
          }
        }
      }
        
      // evaluate Physics model
      // B_fluxes is a set of field gradients
      if (physicsModel->has_B_integrand()) {
        GRAD_FIELDS B_fluxes;
        physicsModel->B_integrand(rhs_mask, fieldsAtIPs, grad_fieldsAtIPs, B_fluxes, imat);
        
        // assemble
        for (field = fields.begin(); field != fields.end(); field++) 
        {
          thisFieldName = field->first;
          numFieldDOF = (field->second).nCols();
          if (!rhs_mask(thisFieldName,FLUX)) continue;
          
          for (int j = 0; j < nsd; j++) 
          {
            Nmat = dN[j].transMat(weights*B_fluxes[thisFieldName][j]);
            for (int i = 0; i < nNodesPerElement; ++i) 
            {  
              inode = conn(i); 
              for (int k = 0; k < numFieldDOF; ++k) 
              {
                rhs[thisFieldName](inode,k) += Nmat(i,k);
              } // scatter k
            } // scatter i
          } // loop over nsd
        } //loop on fields
      } // if B_integrand
    } // element loop
  }  // function call

  // -------------------------------------------------------------
  // compute rhs using given quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_rhs_vector(const Array2D<bool> &rhs_mask,
                                     const FIELDS        &fields,
                                     const PhysicsModel  *physicsModel,
                                     const Array<set<int> >&pointMaterialGroups,
                                     const DIAG_MAT      &weights,
                                     const SPAR_MAT      &N,
                                     const GRAD_SHPFCN   &dN,
                                     FIELDS              &rhs) const
  {
    int nips = weights.nCols();
    int nNodesUnique = feMesh_->get_nNodesUnique();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();

    FieldName thisFieldName;
    FIELDS::const_iterator field;
    FIELDS rhs_local;
    for (field = rhs.begin(); field != rhs.end(); field++) {
      thisFieldName = field->first;
      if (rhs_mask(thisFieldName,FLUX) || rhs_mask(thisFieldName,SOURCE)) {
        int nrows = (field->second).nRows();
        int ncols = (field->second).nCols();
        rhs      [thisFieldName].reset(nrows,ncols);
        rhs_local[thisFieldName].reset(nrows,ncols);
      }
    }

    if (nips>0) {
      // compute fields and gradients of fields at given ips
      GRAD_FIELDS grad_fieldsAtIPs;
      FIELDS fieldsAtIPs;
      int numFieldDOF;
      for (field = fields.begin(); field != fields.end(); field++) {
        thisFieldName = field->first;
        const FIELD & field = (fields.find(thisFieldName))->second;
        numFieldDOF = field.nCols();
        grad_fieldsAtIPs[thisFieldName].assign(nsd,DENS_MAT(nips,numFieldDOF));
        fieldsAtIPs[thisFieldName].reset(nips,numFieldDOF);
        fieldsAtIPs[thisFieldName] = N*field;
        for (int j = 0; j < nsd; ++j) {
          grad_fieldsAtIPs[thisFieldName][j] = dN[j]*field;
        }
      }

      // treat single material point sets specially
      int nMatls = pointMaterialGroups.size();
      int atomMatls = 0;
      for (int imat = 0; imat < nMatls; imat++) {
        const set<int> & indices = pointMaterialGroups(imat);
        if ( indices.size() > 0) atomMatls++;
      }
      bool singleMaterial = ( atomMatls == 1 );
    
      // setup data structures
      FIELDS N_fluxes;
      GRAD_FIELDS B_fluxes;
      if(physicsModel->has_B_integrand()) {
        for (field = fields.begin(); field != fields.end(); field++) {
          thisFieldName = field->first;
          int index = (int) thisFieldName;
          if ( rhs_mask(index,FLUX) ) {
            B_fluxes[thisFieldName].reserve(nsd);
          }
        }
      }
      // evaluate physics model
      if (singleMaterial) 
      {
        if(physicsModel->has_N_integrand()) {
          physicsModel->N_integrand(rhs_mask, fieldsAtIPs, grad_fieldsAtIPs, 
            N_fluxes);
        }
        if(physicsModel->has_B_integrand()) {
          physicsModel->B_integrand(rhs_mask, fieldsAtIPs, grad_fieldsAtIPs, 
            B_fluxes);
        }
      }
      else 
      {
        // NOTE this copy in fields/copy out fluxes is extremely slow 
        // NOTE alternates: 
        // 1) permanent workspace with per-row mapped clones per matl
        //    from caller/atc
        // 2) per point calls to N/B_integrand
        // 3) collect internalToAtom into materials and use mem-cont clones
        //    what about dof for matrices and data storage: clone _matrix_

        // for each material group: 
        // set up storage
        FIELDS      groupN_fluxes, fieldsGroup;
        GRAD_FIELDS groupB_fluxes, grad_fieldsGroup;
        if(physicsModel->has_B_integrand()) {
          for (field = fields.begin(); field != fields.end(); field++) {
            thisFieldName = field->first;
            int index = (int) thisFieldName;
            int ndof = (field->second).nCols();
            grad_fieldsGroup[thisFieldName].assign(nsd,DENS_MAT(nips,ndof));
            N_fluxes[thisFieldName].reset(nips,ndof);
            B_fluxes[thisFieldName].assign(nsd,DENS_MAT(nips,ndof));
            if ( rhs_mask(index,FLUX) ) {
              groupB_fluxes[thisFieldName].reserve(nsd);
            }
          }
        }
        // copy fields
        for ( int imat = 0; imat < pointMaterialGroups.get_length(); imat++) 
        {
          const set<int> & indices = pointMaterialGroups(imat);
          int npts = indices.size();
          int i = 0;
          for (set<int>::const_iterator iter=indices.begin();
                                        iter != indices.end(); iter++, i++ ) {
            for (field = fields.begin(); field != fields.end(); field++) {
              thisFieldName = field->first;
              int ndof = (field->second).nCols();
              fieldsGroup[thisFieldName].reset(npts,ndof);
              for (int j = 0; j < nsd; ++j) {
                (grad_fieldsGroup[thisFieldName][j]).reset(npts,ndof);
              }
              for (int dof = 0; dof < ndof; ++dof) {
                fieldsGroup[thisFieldName](i,dof) 
                  = fieldsAtIPs[thisFieldName](*iter,dof);
                for (int j = 0; j < nsd; ++j) {
                  grad_fieldsGroup[thisFieldName][j](i,dof) 
                    = grad_fieldsAtIPs[thisFieldName][j](*iter,dof);
                }
              }
            }
          }
          // calculate forces, & assemble 
          if(physicsModel->has_N_integrand()) {
            physicsModel->N_integrand(rhs_mask, fieldsGroup, grad_fieldsGroup, 
              groupN_fluxes, imat);
            for (field = fields.begin(); field != fields.end(); field++) {
              thisFieldName = field->first;
              int index = (int) thisFieldName;
              int ndof = (field->second).nCols();
              if ( rhs_mask(index,SOURCE) ) {
                int i = 0;
                for (set<int>::const_iterator iter=indices.begin();
                                          iter != indices.end(); iter++, i++ ) {
                  for (int dof = 0; dof < ndof; ++dof) {
                    N_fluxes[thisFieldName](*iter,dof) 
                      = groupN_fluxes[thisFieldName](i,dof);
                  }
                }
              }
            }
          }
          // calculate forces, & assemble 
          if(physicsModel->has_B_integrand()) {
            physicsModel->B_integrand(rhs_mask, fieldsGroup, grad_fieldsGroup, 
              groupB_fluxes, imat);
            for (field = fields.begin(); field != fields.end(); field++) {
              thisFieldName = field->first;
              int index = (int) thisFieldName;
              int ndof = (field->second).nCols();
              if ( rhs_mask(index,FLUX) ) {
                int i = 0;
                for (set<int>::const_iterator iter=indices.begin();
                                          iter != indices.end(); iter++, i++ ) {
                  for (int dof = 0; dof < ndof; ++dof) {
                    for (int j = 0; j < nsd; ++j) {
                      B_fluxes[thisFieldName][j](*iter,dof) 
                        = groupB_fluxes[thisFieldName][j](i,dof);
                    }
                  }
                }
              }
            }
          }
        }
      }
    
      // assemble
      for (field = fields.begin(); field != fields.end(); field++) 
      {
        thisFieldName = field->first;
        int index = (int) thisFieldName;
        if(physicsModel->has_N_integrand()) {
          if ( rhs_mask(index,SOURCE) ) {
            rhs_local[thisFieldName] 
              += N.transMat(weights*N_fluxes[thisFieldName]);
          }
        }
        if(physicsModel->has_B_integrand()) {
          if ( rhs_mask(index,FLUX) ) {
            for (int j = 0; j < nsd; ++j) {
              rhs_local[thisFieldName] 
                += dN[j].transMat(weights*B_fluxes[thisFieldName][j]);
            }
          }
        }
      }
    }
    
    // Share information between processors
    for (field = rhs.begin(); field != rhs.end(); field++) {
      thisFieldName = field->first;
      if (rhs_mask(thisFieldName,FLUX) || rhs_mask(thisFieldName,SOURCE)) {
        int count = rhs_local[thisFieldName].size();
        //rhs_local[thisFieldName].print("RHS LCL");
        LammpsInterface::instance()->allsum(
          rhs_local[thisFieldName].get_ptr(),
          rhs[thisFieldName].get_ptr(), count);
      }
    }
  }

  // -------------------------------------------------------------
  // compute flux for post processing
  // -------------------------------------------------------------
  void FE_Engine::compute_flux(const Array2D<bool> &rhs_mask,
                                     const FIELDS &fields,
                                     const PhysicsModel * physicsModel,
                                     const Array<int>   & elementMaterials,
                                     GRAD_FIELDS &flux,
                                     const Array<bool> *element_mask) const
  {
    //NOTE move to initialization and make workspace
    const int nNodesPerElement = feMesh_->get_nNodesPerElement();
    const int nIPsPerElement   = feMesh_->get_nIPsPerElement();
    const int nsd              = feMesh_->get_nSpatialDimensions();
    const int nelems           = feMesh_->get_nElements();
    DENS_MAT Nmat;

    //NOTE move to initialization and make workspace
    DENS_MAT N(nIPsPerElement,nNodesPerElement);
    vector<DENS_MAT> dN(nsd);
    dN.assign(nsd, DENS_MAT(nIPsPerElement,nNodesPerElement));
  
    FIELDS fieldsAtIPs, localElementFields;
    GRAD_FIELDS grad_fieldsAtIPs;
    FIELDS::const_iterator field;
    FieldName thisFieldName;
    int numFieldDOF;

    for (field = fields.begin(); field != fields.end(); field++) {
      thisFieldName = field->first;
      if (rhs_mask(thisFieldName,FLUX)) {
        int nrows = (field->second).nRows();
        int ncols = (field->second).nCols();
        flux[thisFieldName].assign(nsd, DENS_MAT(nrows,ncols));
      }
    }

    DIAG_MAT weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);

    // element wise assembly
    int inode = -1;
    for (int ielem=0; ielem < nelems; ielem++) 
    {
      // if element is masked, skip it
      if (element_mask && !(*element_mask)(ielem)) continue;
      // material id
      int imat = elementMaterials(ielem);
      // evaluate shape functions
      feMesh_->shape_function(ielem, N, dN, weights);
      int nips = weights.nRows();
      // get connectivity
      feMesh_->element_connectivity_unique(ielem, conn);
      // compute fields and gradients of fields ips of this element

      for (field = fields.begin(); field != fields.end(); field++) 
      {
        // field values at all nodes
        const DENS_MAT &vI = field->second;
        // field values at element nodes
        DENS_MAT &vIe = localElementFields[field->first];
        // field values at integration points -> to be computed
        DENS_MAT &vP = fieldsAtIPs[field->first];
        // gradients of field at integration points -> to be computed
        vector<DENS_MAT> &dvP = grad_fieldsAtIPs[field->first];

        numFieldDOF = vI.nCols();
        vIe.reset(nNodesPerElement, numFieldDOF);
        // gather local field
        for (int i = 0; i < nNodesPerElement; i++)
          for (int j = 0; j < numFieldDOF; j++)
            vIe(i,j) = vI(conn(i),j);
          
        // interpolate field at integration points
        vP = N*vIe;
        dvP.assign(nsd, DENS_MAT(nIPsPerElement, numFieldDOF));
        for (int j = 0; j < nsd; ++j) dvP[j] = dN[j]*vIe;
      }

      // NOTE move scaling of N by weights up here

      // evaluate Physics model
      // B_fluxes is a set of field gradients
      if (physicsModel->has_B_integrand()) {
        GRAD_FIELDS B_fluxes;
        physicsModel->B_integrand(rhs_mask, fieldsAtIPs, grad_fieldsAtIPs, B_fluxes, imat);
        
        // assemble
        for (field = fields.begin(); field != fields.end(); field++) 
        {
          thisFieldName = field->first;
          numFieldDOF = (field->second).nCols();
          if (!rhs_mask(thisFieldName,FLUX)) continue;
          
          for (int j = 0; j < nsd; j++) 
          {
            Nmat = N.transMat(weights*B_fluxes[thisFieldName][j]);
            for (int i = 0; i < nNodesPerElement; ++i) 
            {  
              inode = conn(i); 
              for (int k = 0; k < numFieldDOF; ++k) 
              {
                flux[thisFieldName][j](inode,k) += Nmat(i,k);
              } // scatter k
            } // scatter i
          } // loop over nsd
        } //loop on fields
      } // if B_integrand
    } // element loop
  }  // function call

  //-----------------------------------------------------------------
  // boundary flux using native quadrature
  //-----------------------------------------------------------------
  void FE_Engine::compute_boundary_flux(
    const Array2D<bool> & rhs_mask,
    const FIELDS        & fields,
    const PhysicsModel  * physicsModel,
    const Array<int>    & elementMaterials,
    const set< pair <int,int> > &  faceSet,
    FIELDS & rhs) const
  {
    //NOTE move to initialization and make workspace
    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nIPsPerFace = feMesh_->get_nIPsPerFace();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();
    DENS_MAT Nmat(nNodesPerElement,3);

    FieldName thisFieldName;
    FIELDS::const_iterator field;
    for (field = fields.begin(); field != fields.end(); field++) {
      thisFieldName = field->first;
      if (rhs_mask(thisFieldName,FLUX)) {
        int nrows = (field->second).nRows();
        int ncols = (field->second).nCols();
        rhs      [thisFieldName].reset(nrows,ncols);
      }
    }

    //NOTE move to initialization and make workspace
    // sizing working arrays
    DENS_MAT N(nIPsPerFace,nNodesPerElement);
    vector< DENS_MAT > dN(nsd), nN(nsd);
    for (int i = 0; i < nsd; i++) {
      dN[i].reset(nIPsPerFace,nNodesPerElement);
      nN[i].reset(nIPsPerFace,nNodesPerElement);
    }
  
    GRAD_FIELDS grad_fieldsAtIPs;
    FIELDS fieldsAtIPs;
    FIELDS localElementFields;
    int numFieldDOF;
    for (field = fields.begin(); field != fields.end(); field++) {
      thisFieldName = field->first;
      numFieldDOF = field->second.nCols();
      grad_fieldsAtIPs[thisFieldName].reserve(nsd);
      for (int i = 0; i < nsd; ++i) {
        grad_fieldsAtIPs[thisFieldName].push_back(DENS_MAT(nIPsPerFace,numFieldDOF));
      }
      fieldsAtIPs[thisFieldName].reset(nIPsPerFace,numFieldDOF);
      localElementFields[thisFieldName].reset(nNodesPerElement,numFieldDOF);
    }
  
    DiagonalMatrix<double> weights(nIPsPerFace,nIPsPerFace);
    Array<int> conn(nNodesPerElement);
  
    // element wise assembly
    int inode = -1;
    set< pair <int,int> >::iterator iter;
    for (iter = faceSet.begin(); iter != faceSet.end(); iter++) {
      // evaluate shape functions at ips
      feMesh_->face_shape_function(*iter, N, dN, nN, weights);
      int nips = weights.nRows();
      // get connectivity
      int ielem = iter->first;
      int imat = elementMaterials(ielem);
      feMesh_->element_connectivity_unique(ielem, conn);

      // interpolate fields and gradients of fields ips of this element
      for (field = fields.begin(); field != fields.end(); field++) {
        thisFieldName = field->first;
        const DenseMatrix<double> * thisField  
          = (const DENS_MAT *) &(field->second);
        numFieldDOF = thisField->nCols();
        //NOTE an hopefully unnecessary copy (should alias)
        for (int i = 0; i < nNodesPerElement; i++) {
          for (int j = 0; j < numFieldDOF; j++) {
            localElementFields[thisFieldName](i,j) = (*thisField)(conn(i),j);
          }
        }
        //  ips X dof    = ips X nodes * nodes X dof
        fieldsAtIPs[thisFieldName] = N*localElementFields[thisFieldName];
        for (int j = 0; j < nsd; ++j) {
          grad_fieldsAtIPs[thisFieldName][j] = dN[j]*localElementFields[thisFieldName];
        }
      }

      // Evaluate- physics model
      // do nothing for N_integrand
      // nN : precomputed and held by ATC_Transfer
      if(physicsModel->has_B_integrand()) {
        map<FieldName, vector< DENS_MAT > > B_fluxes;
        for (field = fields.begin(); field != fields.end(); field++) {
          thisFieldName = field->first;
          if ( rhs_mask(thisFieldName,FLUX) ) {
            B_fluxes[thisFieldName].reserve(nsd);
          }
        }
        physicsModel->B_integrand(rhs_mask, fieldsAtIPs, grad_fieldsAtIPs, B_fluxes,imat);
        // assemble
        for (field = fields.begin(); field != fields.end(); field++) {
          thisFieldName = field->first;
          int index = (int) thisFieldName;
          const DenseMatrix<double> * thisField  
            = (const DENS_MAT *) &(field->second);
          if ( rhs_mask(index,FLUX) ) {
            numFieldDOF = thisField->nCols();
            Nmat.reset(nNodesPerElement,numFieldDOF);
            for (int j = 0; j < nsd; j++) {
              // nodes X dof = nodes X ips * ips X dof  ???
              Nmat = nN[j].transMat(weights*B_fluxes[thisFieldName][j]);
              for (int i = 0; i < nNodesPerElement; ++i) {
                inode = conn(i);
                for (int k = 0; k < numFieldDOF; ++k) {
                  rhs[thisFieldName](inode,k) += Nmat(i,k);
                }
              }
            }
          }
        }
      }

    }
  }

  // -------------------------------------------------------------
  // compute boundary flux using given quadrature and interpolation
  // -------------------------------------------------------------
  void FE_Engine::compute_boundary_flux( 
    const Array2D<bool>                  & rhs_mask,
    const FIELDS                         & fields,
    const PhysicsModel                   * physicsModel,
    const Array< int >                   & elementMaterials,
    const Array< set<int> >              & pointMaterialGroups,
    const DIAG_MAT                       & weights,
    const SPAR_MAT                       & N,
    const GRAD_SHPFCN                    & dN,
    const DIAG_MAT                       & flux_mask,
    FIELDS                               & flux) const
  {
    int nNodesUnique = feMesh_->get_nNodesUnique();
    FieldName thisFieldName;
    
    map<FieldName, DENS_MAT > rhs;
    map<FieldName, DENS_MAT > rhsSubset;
    map<FieldName, DENS_MAT > rhsSubset_local;

    FIELDS::const_iterator field;
    for (field = fields.begin(); field != fields.end(); field++) {
      thisFieldName = field->first;
      if (rhs_mask(thisFieldName,FLUX)) {
        int nrows = (field->second).nRows();
        int ncols = (field->second).nCols();
        rhs      [thisFieldName].reset(nrows,ncols);
        rhsSubset[thisFieldName].reset(nrows,ncols);
      }
    }

    // R_I    = - int_Omega    Delta N_I .q dV
    compute_rhs_vector(rhs_mask, fields, physicsModel, elementMaterials, rhs);
    // R_I^md = - int_Omega^md Delta N_I .q dV
    compute_rhs_vector(rhs_mask, fields, physicsModel, pointMaterialGroups, 
      weights, N, dN, rhsSubset);
    
    // flux_I = 1/Delta V_I V_I^md R_I + R_I^md
    //   where : Delta V_I = int_Omega N_I dV
    for (int n = 0; n < NUM_FIELDS; ++n) {
      FieldName thisFieldName = (FieldName) n;
      if ( rhs_mask(thisFieldName,FLUX) ) { 
        flux[thisFieldName] 
          = rhsSubset[thisFieldName] - flux_mask*rhs[thisFieldName];
      }
    }
  }
  
  //-----------------------------------------------------------------
  // prescribed boundary flux using native quadrature 
  // integrate \int_delV N_I s(x,t).n dA
  //-----------------------------------------------------------------
  void FE_Engine::add_fluxes(const Array<bool> &fieldMask,
                             const double time,
                             const SURFACE_SOURCE & sourceFunctions,
                             FIELDS &nodalSources) const
  {
    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nNodesPerFace = feMesh_->get_nNodesPerFace();
    int nIPsPerFace = feMesh_->get_nIPsPerFace();
    int nsd = feMesh_->get_nSpatialDimensions();
    FieldName thisField;
    int nFieldDOF = 1;
    XT_Function * f = NULL; 

    //NOTE move arrays to initialization and make workspace
    // sizing working arrays
    vector<DENS_MAT> dN, nN;
    dN.assign(nsd, DENS_MAT(nIPsPerFace, nNodesPerElement));
    nN.assign(nsd, DENS_MAT(nIPsPerFace, nNodesPerElement));

    DENS_MAT N(nIPsPerFace,nNodesPerElement);
    DENS_MAT xCoords(nsd,nNodesPerElement);
    DENS_MAT xAtIPs(nsd,nIPsPerFace);
    DENS_MAT Nmat(nNodesPerElement,nsd);
    DIAG_MAT weights(nIPsPerFace,nIPsPerFace);
    Array<int> conn(nNodesPerElement);
    DENS_MAT faceSource(nIPsPerFace,nFieldDOF);

    // element wise assembly
    SURFACE_SOURCE::const_iterator src_iter;
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++)
    {
      FieldName thisFieldName = src_iter->first;
      if (!fieldMask((int)thisFieldName)) continue; 
      
      typedef map<PAIR,Array<XT_Function*> > FSET;
      const FSET *fset = (const FSET *)&(src_iter->second);
      FSET::const_iterator fset_iter;
      for (fset_iter = fset->begin(); fset_iter != fset->end(); fset_iter++)
      {
        const PAIR &face = fset_iter->first;
        const int elem = face.first;
        const Array <XT_Function*> &fs = fset_iter->second;
        feMesh_->element_connectivity_unique(elem, conn);
        // evaluate location at ips
        feMesh_->face_shape_function(face, N, dN, nN, weights);
        feMesh_->element_coordinates(elem, xCoords);
        MultAB(xCoords,N,xAtIPs,0,1); //xAtIPs = xCoords*(N.transpose());

        // interpolate prescribed flux at ips of this element
        // NOTE this seems totally redundant
        FSET::const_iterator face_iter = fset->find(face);
        if (face_iter == fset->end()) 
        {  
          cout << "face not found" << std::endl;
          // NOTE: do something there
        }
        // NOTE why can't I know the size earlier??
        // NOTE only normal flux here
        nFieldDOF = (face_iter->second).size();
        faceSource.reset(nIPsPerFace,nFieldDOF);
        for (int ip = 0; ip < nIPsPerFace; ++ip) 
        {
          for (int idof = 0; idof<nFieldDOF; ++idof) 
          {
            f = fs(idof);
            if (!f) continue;
            faceSource(ip,idof) = f->f(column(xAtIPs,ip).get_ptr(),time);
          }
        }
        // assemble
        Nmat = N.transMat(weights*faceSource);
        for (int i = 0; i < nNodesPerElement; ++i) 
        {
          int inode = conn(i);
          for (int idof = 0; idof < nFieldDOF; ++idof) 
          {
            nodalSources[thisFieldName](inode,idof) += Nmat(i,idof);
          }  // end assemble nFieldDOF
        }  // end assemble nNodesPerElement
      }  // end fset loop
    } // field loop
  }

  //-----------------------------------------------------------------
  // prescribed volume flux using native quadrature 
  // integrate \int_V N_I s(x,t) dV
  //-----------------------------------------------------------------
  void FE_Engine::add_sources(const Array<bool> &fieldMask,
                              const double time,
                              const VOLUME_SOURCE &sourceFunctions,
                              FIELDS &nodalSources) const
  {
    int nNodes = feMesh_->get_nNodes();
    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nIPsPerElement = feMesh_->get_nIPsPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();
    int nFieldDOF = 1;
    XT_Function * f = NULL; 

    //NOTE move arrays to initialization and make workspace
    // sizing working arrays
    DENS_MAT elemSource(nIPsPerElement,nFieldDOF);
    DENS_MAT N(nIPsPerElement,nNodesPerElement);
    DENS_MAT xCoords(nsd,nNodesPerElement);
    DENS_MAT xAtIPs(nsd,nIPsPerElement);
    DENS_MAT Nwg(nNodesPerElement,1);
    DENS_MAT Nmat(nNodesPerElement,nsd);
    DiagonalMatrix<double> weights(nIPsPerElement,nIPsPerElement);
    Array<int> conn(nNodesPerElement);

    // element wise assembly
    for (int ielem = 0; ielem < nelems; ++ielem) {
      feMesh_->element_connectivity_unique(ielem, conn);

      // evaluate location at ips
      feMesh_->shape_function(ielem, N, weights);
      feMesh_->element_coordinates(ielem, xCoords);
      xAtIPs =xCoords*(N.transpose());

      VOLUME_SOURCE ::const_iterator src_iter;
      for (src_iter = sourceFunctions.begin(); 
           src_iter != sourceFunctions.end(); src_iter++) {
        FieldName thisField = src_iter->first;
        int index = (int) thisField;
        if ( fieldMask(index) ) {
          const Array2D<XT_Function *> * thisSource 
            = (const Array2D<XT_Function *> *) &(src_iter->second);
          nFieldDOF = thisSource->nCols();
          elemSource.reset(nIPsPerElement,nFieldDOF);

          // interpolate prescribed flux at ips of this element
          for (int ip = 0; ip < nIPsPerElement; ++ip) {
            for (int idof = 0; idof < nFieldDOF; ++idof) {
              f = (*thisSource)(ielem,idof);
              if (f) {
                elemSource(ip,idof) = f->f(column(xAtIPs,ip).get_ptr(),time);
              }
            }
          }

          // assemble
          Nmat.reset(nNodesPerElement,nFieldDOF);
          Nmat = N.transMat(weights*elemSource);
          for (int i = 0; i < nNodesPerElement; ++i) {
            int inode = conn(i);
            for (int idof = 0; idof < nFieldDOF; ++idof) {
              nodalSources[thisField](inode,idof) += Nmat(i,idof);
            }
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------
  // boundary integral of a nodal field
  //-----------------------------------------------------------------
  void FE_Engine::field_surface_flux(
    const DENS_MAT & field,
    const set< PAIR > &  faceSet,
    DENS_MAT & values,
    const bool contour,
    const int axis) const
  {
    int nNodesPerElement = feMesh_->get_nNodesPerElement();
    int nIPsPerFace = feMesh_->get_nIPsPerFace();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nelems = feMesh_->get_nElements();
    int numFieldDOF = field.nCols();

    double a[3] = {0,0,0};
    a[axis] = 1;

    // sizing working arrays
    DENS_MAT N(nIPsPerFace,nNodesPerElement);
    DENS_MAT n(nsd,nIPsPerFace);
    DENS_MAT fieldsAtIPs(nIPsPerFace,numFieldDOF);
    DENS_MAT localElementFields(nNodesPerElement,numFieldDOF);
    DiagonalMatrix<double> weights(nIPsPerFace,nIPsPerFace);
    Array<int> conn(nNodesPerElement);

    DENS_MAT Nmat(nsd,numFieldDOF); 
    DENS_MAT integrals(numFieldDOF,nsd); 

    // element wise assembly
    set< pair <int,int> >::iterator iter;
    for (iter = faceSet.begin(); iter != faceSet.end(); iter++) {
      // evaluate shape functions at ips
      feMesh_->face_shape_function(*iter, N, n, weights);
      // cross n with axis to get tangent
      if (contour) {
        double t[3];
        for (int i = 0; i < nIPsPerFace; i++) {
          t[0] = a[1]*n(2,i) - a[2]*n(1,i);
          t[1] = a[2]*n(0,i) - a[0]*n(2,i);
          t[2] = a[0]*n(1,i) - a[1]*n(0,i);
          n(0,i) = t[0];
          n(1,i) = t[1];
          n(2,i) = t[2];
        }
      }
      int nips = weights.nRows();
      // get connectivity
      int ielem = iter->first;
      feMesh_->element_connectivity_unique(ielem, conn);

      // interpolate fields and gradients of fields ips of this element
      for (int i = 0; i < nNodesPerElement; i++) {
        for (int j = 0; j < numFieldDOF; j++) {
          localElementFields(i,j) = field(conn(i),j);
        }
      }
      //  ips X dof    = ips X nodes * nodes X dof
      fieldsAtIPs = N*localElementFields;

      // assemble : integral(k,j) = sum_ip n(j,ip) wg(ip,ip) field(ip,k)
      Nmat = n*weights*fieldsAtIPs;
      for (int j = 0; j < nsd; j++) {
        for (int k = 0; k < numFieldDOF; ++k) {
          integrals(k,j) += Nmat(j,k);
        }
      }
    }
    // (S.n)_1 = S_1j n_j = S_11 n_1 + S_12 n_2 + S_13 n_3
    // (S.n)_2 = S_2j n_j = S_21 n_1 + S_22 n_2 + S_23 n_3
    // (S.n)_3 = S_3j n_j = S_31 n_1 + S_32 n_2 + S_33 n_3
    if (numFieldDOF == 9) { // tensor
      values.reset(nsd,1);
      values(0,0) = integrals(0,0)+integrals(1,1)+integrals(2,2);
      values(1,0) = integrals(3,0)+integrals(4,1)+integrals(5,2);
      values(2,0) = integrals(6,0)+integrals(7,1)+integrals(8,2);
    }
    else if (numFieldDOF == 6) { // sym tensor
      values.reset(nsd,1);
      values(0,0) = integrals(0,0)+integrals(3,1)+integrals(4,2);
      values(1,0) = integrals(3,0)+integrals(1,1)+integrals(5,2);
      values(2,0) = integrals(4,1)+integrals(5,1)+integrals(2,2);
    }
    // (v.n) = v_j n_j 
    else if (numFieldDOF == 3) { // vector
      values.reset(1,1);
      values(0,0) = integrals(0,0)+integrals(1,1)+integrals(2,2);
    }
    // s n = s n_j e_j
    else if (numFieldDOF == 1) { // scalar
      values.reset(nsd,1);
      values(0,0) = integrals(0,0);
      values(1,0) = integrals(0,1);
      values(2,0) = integrals(0,2);
    }
    else {
      string msg = "FE_Engine::field_surface_flux unsupported field width: ";
      msg += tostring(numFieldDOF);
      throw ATC_Error(0,msg);
    }
  }


  //-----------------------------------------------------------------
  // evaluate shape functions at given points
  //-----------------------------------------------------------------
  void FE_Engine::evaluate_shape_functions(
    const MATRIX & pt_coords,
    SPAR_MAT & N,
    Array<int> & pointToEltMap) const
  {
    // Get shape function and derivatives at atomic locations
    int npe = feMesh_->get_nNodesPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nnodes = feMesh_->get_nNodesUnique();
    int npts = pt_coords.nCols();

    pointToEltMap.reset(npts);
    
    // loop over point (atom) coordinates
    DENS_VEC x(nsd);
    Array<int> node_index(npe);
    DENS_VEC shp(npe);
    
    N.reset(npts,nnodes);
    for (int i = 0; i < npts; ++i) {
      for (int k = 0; k < nsd; ++k) {
        x(k) = pt_coords(k,i);
      }
      int eltID;
      feMesh_->shape_functions(x,shp,eltID,node_index);
      for (int j = 0; j < npe; ++j) {
        int jnode = node_index(j);
        N.add(i,jnode,shp(j));
      }
      pointToEltMap(i) = eltID;
    }
    N.compress();   
  }

  //-----------------------------------------------------------------
  // evaluate shape functions and their spatial derivatives at given points
  //-----------------------------------------------------------------
  void FE_Engine::evaluate_shape_functions(
    const MATRIX & pt_coords,
    SPAR_MAT & N,
    GRAD_SHPFCN & dN,
    Array<int> & pointToEltMap) const
  {
    // Get shape function and derivatives at atomic locations
    int npe = feMesh_->get_nNodesPerElement();
    int nsd = feMesh_->get_nSpatialDimensions();
    int nnodes = feMesh_->get_nNodesUnique();
    int npts = pt_coords.nCols();
    
    pointToEltMap.reset(npts);

    // loop over point (atom) coordinates
    DENS_VEC x(nsd);
    Array<int> node_index(npe);
    DENS_VEC shp(npe);
    DENS_MAT dshp(nsd,npe);
    
    for (int k = 0; k < nsd; ++k) { 
      dN[k].reset(npts,nnodes); 
    }
          
    N.reset(npts,nnodes); 
    for (int i = 0; i < npts; ++i) {
      for (int k = 0; k < nsd; ++k) {
        x(k) = pt_coords(k,i);
      }
      int eltID;
      feMesh_->shape_functions(x,shp,dshp,eltID,node_index);
      for (int j = 0; j < npe; ++j) {
        int jnode = node_index(j);
        N.add(i,jnode,shp(j));
        for (int k = 0; k < nsd; ++k) {
          dN[k].add(i,jnode,dshp(k,j));
        }
      }
      pointToEltMap(i) = eltID;
    }
    N.compress();   
    for (int k = 0; k < nsd; ++k) { 
      dN[k].compress();   
    }
  }

  //-----------------------------------------------------------------
  // evaluate shape functions at given points
  //-----------------------------------------------------------------
  void FE_Engine::evaluate_shape_functions(const VECTOR & x,
                                           Array<int>& node_index,
                                           DENS_VEC& shp, 
                                           DENS_MAT& dshp,
                                           int & eltID) const
  {
    // Get shape function and derivatives at a specific point
    int nsd = feMesh_->get_nSpatialDimensions();
    int npe = feMesh_->get_nNodesPerElement();

    dshp.reset(nsd,npe);

    feMesh_->shape_functions(x,shp,dshp,eltID,node_index);
  
  }

  //-----------------------------------------------------------------
  // create a uniform rectangular mesh on a rectangular region
  //-----------------------------------------------------------------
  void FE_Engine::create_mesh(int nx, int ny, int nz, char * regionName,
                              int xperiodic, int yperiodic, int zperiodic)
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xscale, yscale, zscale;
    double boxxlo, boxxhi, boxylo, boxyhi, boxzlo, boxzhi;

    // check to see if region exists and is it a box, and if so get the bounds
    bool isBox;
    isBox = atcTransfer_->get_region_bounds(regionName,
                                            xmin, xmax,
                                            ymin, ymax,
                                            zmin, zmax,
                                            xscale,
                                            yscale,
                                            zscale);
    if (!isBox) throw ATC_Error(0,"Region for FE mesh is not a box");
    
    feMesh_ = new FE_Uniform3DMesh(nx,ny,nz,
                                   xmin,xmax,
                                   ymin,ymax,
                                   zmin,zmax,
                                   xscale,
                                   yscale,
                                   zscale,
                                   xperiodic,
                                   yperiodic,
                                   zperiodic);
    if (LammpsInterface::instance()->comm_rank()==0) { 
      printf(" ATC:: created FEM Mesh with %i Global Nodes, %i Unique Nodes, and %i Elements\n",
             feMesh_->get_nNodes(),feMesh_->get_nNodesUnique(),feMesh_->get_nElements()); }
  }

};

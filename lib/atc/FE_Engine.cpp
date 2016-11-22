#include "FE_Engine.h"

#include "ATC_Transfer.h"
#include "FE_Element.h"
#include "Function.h"
#include "PhysicsModel.h"
#include "KernelFunction.h"
#include "Utility.h"
#include "MPI_Wrappers.h"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string>





using namespace std;
using ATC_Utility::is_numeric;
using ATC_Utility::to_string;;
using MPI_Wrappers::allsum;
using MPI_Wrappers::int_allsum;
using MPI_Wrappers::rank_zero;
using MPI_Wrappers::print_msg;
using MPI_Wrappers::print_msg_once;

namespace ATC{

  static const double tol_sparse = 1.e-30;//tolerance for compaction from dense

  //-----------------------------------------------------------------
  FE_Engine::FE_Engine(MPI_Comm comm)
    : communicator_(comm),
      feMesh_(NULL),
      initialized_(false),
      outputManager_()
  {
    // Nothing to do here
  }
  
  //-----------------------------------------------------------------
  FE_Engine::~FE_Engine()
  {
    if (feMesh_) delete feMesh_;
  }
  //-----------------------------------------------------------------
  
  void FE_Engine::initialize()
  {
    if (!feMesh_) throw ATC_Error("FE_Engine has no mesh");

    if (!initialized_) {
      // set up work spaces
      nNodesPerElement_ = feMesh_->num_nodes_per_element();
      nIPsPerElement_ = feMesh_->num_ips_per_element();
      nIPsPerFace_ = feMesh_->num_ips_per_face();
      nSD_ = feMesh_->num_spatial_dimensions();
      nElems_ = feMesh_->num_elements();
      nNodesUnique_ = feMesh_->num_nodes_unique();
      nNodes_ = feMesh_->num_nodes();
    
      // arrays & matrices
      _weights_.reset(nIPsPerElement_,nIPsPerElement_);
      _N_.reset(nIPsPerElement_,nNodesPerElement_);
      _dN_.assign(nSD_, DENS_MAT(nIPsPerElement_,nNodesPerElement_));
      _Nw_.reset(nIPsPerElement_,nNodesPerElement_);
      _dNw_.assign(nSD_, DENS_MAT(nIPsPerElement_,nNodesPerElement_));
      _Bfluxes_.assign(nSD_, DENS_MAT());
      // faces
      _fweights_.reset(nIPsPerElement_,nIPsPerElement_);
      _fN_.reset(nIPsPerFace_,nNodesPerElement_); 
      _fdN_.assign(nSD_, DENS_MAT(nIPsPerFace_, nNodesPerElement_));
      _nN_.assign(nSD_, DENS_MAT(nIPsPerFace_, nNodesPerElement_));

      // remove specified elements
      if (nullElements_.size() > 0) delete_elements(nullElements_);
 
      initialized_ = true;
    }
  }

  void FE_Engine::partition_mesh()
  {
    if (is_partitioned()) return;

    feMesh_->partition_mesh();  
    
    // now do all FE_Engine data structure partitioning

    // partition nullElements_
    /*for (vector<int>::iterator elemsIter = feMesh_->processor_elts().begin();
         elemsIter != feMesh_->processor_elts().end();
         ++elemsIter)
    {
      if (nullElements_.find(*elemsIter) != nullElements_.end()) {
        myNullElements_.insert(map_elem_to_myElem(*elemsIter));
      }
    }*/
  }

  void FE_Engine::departition_mesh()
  {
    if (!is_partitioned()) return;

    feMesh_->departition_mesh();
  }

  void FE_Engine::set_quadrature(FeIntQuadrature type, bool temp) const
  {
    if (!feMesh_) throw ATC_Error("FE_Engine has no mesh");

    feMesh_->set_quadrature(type);
    if (!temp) quadrature_ = type;

    int nIPsPerElement_new = feMesh_->num_ips_per_element();
    int nIPsPerFace_new = feMesh_->num_ips_per_face();
  
    if (nIPsPerElement_ != nIPsPerElement_new) {
      // arrays & matrices
      nIPsPerElement_ = nIPsPerElement_new;
      _weights_.resize(nIPsPerElement_,nIPsPerElement_);
      _N_.resize(nIPsPerElement_,nNodesPerElement_);
      _dN_.assign(nSD_, DENS_MAT(nIPsPerElement_,nNodesPerElement_));
      _Nw_.reset(nIPsPerElement_,nNodesPerElement_);
      _dNw_.assign(nSD_, DENS_MAT(nIPsPerElement_,nNodesPerElement_));
    }
    if (nIPsPerFace_ != nIPsPerFace_new) {
      // faces
      nIPsPerFace_ = nIPsPerFace_new;
      _fweights_.reset(nIPsPerElement_,nIPsPerElement_);
      _fN_.reset(nIPsPerFace_,nNodesPerElement_); 
      _fdN_.assign(nSD_, DENS_MAT(nIPsPerFace_, nNodesPerElement_));
      _nN_.assign(nSD_, DENS_MAT(nIPsPerFace_, nNodesPerElement_));
    }
    
  }

  //-----------------------------------------------------------------
  bool FE_Engine::modify(int narg, char **arg)
  {
    bool match = false;
    /*! \page man_mesh_create fix_modify AtC mesh create
      \section syntax
      fix_modify AtC mesh create <nx> <ny> <nz> <region-id> 
      <f|p> <f|p> <f|p> \n
      - nx ny nz = number of elements in x, y, z
      - region-id = id of region that is to be meshed
      - f p p = periodicity flags for x, y, z
      \section examples
      <TT> fix_modify AtC mesh create 10 1  1  feRegion p p p </TT> \n
      \section description
      Creates a uniform mesh in a rectangular region
      \section restrictions
      Creates only uniform rectangular grids in a rectangular region
      \section related
      \ref man_mesh_quadrature
      \section default
      When created, mesh defaults to gauss2 (2-point Gaussian) quadrature. 
      Use "mesh quadrature" command to change quadrature style.
    */
    int argIdx = 0;
    if (strcmp(arg[argIdx],"mesh")==0) { 
      argIdx++;
      // create mesh
      if (strcmp(arg[argIdx],"create")==0) {
        if (feMesh_) throw ATC_Error("FE Engine already has a mesh");
        argIdx++;
        int nx = atoi(arg[argIdx++]);
        int ny = atoi(arg[argIdx++]);
        int nz = atoi(arg[argIdx++]);
        string box = arg[argIdx++];

        Array<bool> periodicity(3);
        periodicity(0) = (strcmp(arg[argIdx++],"p")==0) ? true : false; 
        periodicity(1) = (strcmp(arg[argIdx++],"p")==0) ? true : false;
        periodicity(2) = (strcmp(arg[argIdx++],"p")==0) ? true : false;

        if (argIdx < narg ) {
          Array<double> dx(nx),dy(ny),dz(nz);

          dx = 0;
          dy = 0;
          dz = 0;
          while (argIdx < narg) {
            if      (strcmp(arg[argIdx],"dx")==0) {
              parse_partitions(++argIdx, narg, arg, 0, dx);
            }
            else if (strcmp(arg[argIdx],"dy")==0) {
              parse_partitions(++argIdx, narg, arg, 1, dy);
            }
            else if (strcmp(arg[argIdx],"dz")==0) {
              parse_partitions(++argIdx, narg, arg, 2, dz);
            }
          }

          create_mesh(dx, dy, dz, box.c_str(), periodicity);
        }
        else  {
          create_mesh(nx, ny, nz, box.c_str(), periodicity);
        }
        quadrature_ = GAUSS2;
        match = true;
      }
    /*! \page man_mesh_quadrature fix_modify AtC mesh quadrature
      \section syntax
      fix_modify AtC mesh quadrature <quad>
      - quad = one of <nodal|gauss1|gauss2|gauss3|face> --- when a mesh is created it defaults to gauss2, use this call to change it after the fact
      \section examples
      <TT> fix_modify AtC mesh quadrature face </TT>
      \section description
      (Re-)assigns the quadrature style for the existing mesh.
      \section restrictions
      \section related
      \ref man_mesh_create
      \section default
      none 
    */
      else if (strcmp(arg[argIdx],"quadrature")==0) {
        argIdx++;
        string quadStr = arg[argIdx];
        FeIntQuadrature quadEnum = string_to_FIQ(quadStr);
        set_quadrature(quadEnum);
        match = true;
      }
    /*! \page man_mesh_read fix_modify AtC mesh read 
      \section syntax
      fix_modify AtC mesh read <filename> <f|p> <f|p> <f|p> 
      - filename = name of file containing mesh to be read 
      - f p p = periodicity flags for x, y, z
      \section examples
      <TT> fix_modify AtC mesh read myComponent.mesh p p p </TT> \n
      <TT> fix_modify AtC mesh read myOtherComponent.exo </TT>
      \section description
      Reads a mesh from a text or exodus file, and assigns periodic
      boundary conditions if needed. 
      \section restrictions
      \section related
      \section default
      periodicity flags are false by default 
    */
      else if (strcmp(arg[argIdx],"read")==0) {
        argIdx++;
        string meshFile = arg[argIdx++];
        
        Array<bool> periodicity(3);
        periodicity = false;
        if (argIdx < narg)  {
          periodicity(0) = (strcmp(arg[argIdx++],"p")==0) ? true : false; 
          periodicity(1) = (strcmp(arg[argIdx++],"p")==0) ? true : false;
          periodicity(2) = (strcmp(arg[argIdx++],"p")==0) ? true : false;
        }
        read_mesh(meshFile,periodicity);

        if (periodicity(0) || periodicity(1) || periodicity(2)) {
          meshFile = "periodic_"+meshFile;
          stringstream ss; 
          ss << "writing periodicity corrected mesh: " << meshFile;
          print_msg(communicator_,ss.str());
          feMesh_->write_mesh(meshFile);
          feMesh_->output(meshFile);
        }
        match = true;
      }
    /*! \page man_mesh_write fix_modify AtC mesh write 
      \section syntax
      fix_modify AtC mesh write <filename> 
      - filename = name of file to write mesh
      \section examples
      <TT> fix_modify AtC mesh write myMesh.mesh </TT> \n
      \section description
      Writes a mesh to a text file. 
      \section restrictions
      \section related
      \section default
    */
      else if (strcmp(arg[argIdx],"write")==0) {
        argIdx++;
        string meshFile = arg[argIdx];
        feMesh_->write_mesh(meshFile);
        match = true;
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
      else if (strcmp(arg[argIdx],"delete_elements")==0) {
        argIdx++;
        string esetName = arg[argIdx];
        set<int> elemSet = feMesh_->elementset(esetName);
        nullElements_.insert(elemSet.begin(), elemSet.end());
        match = true;
      }
      else if (strcmp(arg[argIdx],"cut")==0) {
        argIdx++;
        string fsetName = arg[argIdx++];
        set<PAIR> faceSet = feMesh_->faceset(fsetName);
        cutFaces_.insert(faceSet.begin(), faceSet.end());
        if (narg > argIdx && strcmp(arg[argIdx],"edge")==0) {
          argIdx++;
          string nsetName = arg[argIdx];
          set<int> nodeSet = feMesh_->nodeset(nsetName);
          cutEdge_.insert(nodeSet.begin(), nodeSet.end());
        }
        // cut mesh
        if (cutFaces_.size() > 0) cut_mesh(cutFaces_,cutEdge_);

        match = true;
      }
      else if (strcmp(arg[argIdx],"lammps_partition")==0) {
        feMesh_->set_lammps_partition(true);
        match = true;
      }
      else if (strcmp(arg[argIdx],"data_decomposition")==0) {
        feMesh_->set_data_decomposition(true);
        match = true;
      }
      else {
        if ( ! feMesh_ ) throw ATC_Error("need mesh for parsing");
        match = feMesh_->modify(narg,arg);
      }
    }
    // FE_Mesh
    else {
      if ( ! feMesh_ ) throw ATC_Error("need mesh for parsing");
      match = feMesh_->modify(narg,arg);
    }
    return match;
  }
  //-----------------------------------------------------------------
  void FE_Engine::parse_partitions(int & argIdx, int narg, char ** arg, 
    int idof, Array<double> & dx ) const
  {
    double x[3] = {0,0,0}; 
    int nx = dx.size();
    // parse relative values for each element
    if (is_numeric(arg[argIdx])) {
      for (int i = 0; i < nx; ++i) { 
        if (is_numeric(arg[argIdx])) { dx(i) = atof(arg[argIdx++]); }
        else throw ATC_Error("not enough element partitions");
      }
    }
    // each segment of the piecewise funcion is length-normalized separately
    else if (strcmp(arg[argIdx],"position-number-density")==0) { 
      argIdx++;
      double y[nx],w[nx];
      int n[nx];
      int nn = 0;
      while (argIdx < narg) { 
        if (! is_numeric(arg[argIdx])) break;
        y[nn]   = atof(arg[argIdx++]);
        n[nn]   = atoi(arg[argIdx++]);
        w[nn++] = atof(arg[argIdx++]);
      }
      if (n[nn-1] != nx)  throw ATC_Error("total element partitions do not match");
      int k = 0;
      for (int i = 1; i < nn; ++i) { 
        int dn = n[i]-n[i-1];
        double dy = y[i]-y[i-1];
        double w0 = w[i-1];
        double dw = w[i]-w0;
        double lx = 0;
        double l[dn];
        for (int j = 0; j < dn; ++j) {
          double x = (j+0.5)/dn; 
          double dl = w0+x*dw;
          lx += dl;
          l[j]= dl;
        }
        double scale = dy/lx;
        for (int j = 0; j < dn; ++j) {
          dx(k++) = scale*l[j];
        }
      }
    }
    // construct relative values from a density function
    // evaluate for a domain (0,1)
    else {
      XT_Function * f = XT_Function_Mgr::instance()->function(&(arg[argIdx]),narg-argIdx);
      argIdx++;
      while (argIdx < narg) { 
        if (! is_numeric(arg[argIdx++])) break;
      }
      for (int i = 0; i < nx; ++i) { 
        x[idof] = (i+0.5)/nx; dx(i) = f->f(x,0.); 
      }
    }
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
    string outputPrefix, set<int> otypes)
  {
    outputManager_.initialize(outputPrefix, otypes);
    if (!feMesh_) throw ATC_Error("output needs mesh");
    if (!initialized_) initialize();
    if (!feMesh_->coordinates() || !feMesh_->connectivity()) 
      throw ATC_Error("output mesh not properly initialized");
    if (!feMesh_->coordinates()->nCols() || 
        !feMesh_->connectivity()->nCols()) 
      throw ATC_Error("output mesh is empty");
    if (rank == 0) 
      outputManager_.write_geometry(feMesh_->coordinates(), 
                                    feMesh_->connectivity());
     outputManager_.print_custom_names();
  }

  //-----------------------------------------------------------------
  //  write geometry
  //-----------------------------------------------------------------
  void FE_Engine::write_geometry(void)
  {
    outputManager_.write_geometry(feMesh_->coordinates(), 
                                  feMesh_->connectivity());
  }

  // -------------------------------------------------------------
  //  write data  
  // -------------------------------------------------------------
  void FE_Engine::write_data(double time, 
                             FIELDS &soln, 
                             OUTPUT_LIST *data)
  {
    outputManager_.write_data(
                        time, &soln, data,
                        (feMesh_->node_map())->data());
  }

  // -------------------------------------------------------------
  //  write data  
  // -------------------------------------------------------------
  void FE_Engine::write_data(double time, OUTPUT_LIST *data)
  {
    outputManager_.write_data(
                        time, data, 
                        feMesh_->node_map()->data());
  }

  // -------------------------------------------------------------
  //  amend mesh for deleted elements
  // -------------------------------------------------------------
  void FE_Engine::delete_elements(const set<int> & elementList)
  {
    feMesh_->delete_elements(elementList);
  }

  // -------------------------------------------------------------
  //  amend mesh for cut at specified faces
  // -------------------------------------------------------------
  void FE_Engine::cut_mesh(const set<PAIR> &faceSet, 
                           const set<int> &nodeSet)
  {
    feMesh_->cut_mesh(faceSet,nodeSet);
  }
  // -------------------------------------------------------------
  //  interpolate one value
  // -------------------------------------------------------------
  DENS_VEC  FE_Engine::interpolate_field(const DENS_VEC & x, const FIELD & f) const
  {
    DENS_VEC N;
    Array<int> nodelist;
    feMesh_->shape_functions(x, N, nodelist);
    const DENS_MAT &vI = f.quantity();
    int dof = vI.nCols();
    DENS_MAT vIe(nNodesPerElement_, dof);
    for (int i = 0; i < nNodesPerElement_; i++)
      for (int j = 0; j < dof; j++)
        vIe(i,j) = vI(nodelist(i),j);
    DENS_VEC vP;
    vP = N*vIe;
    return vP;
  }

  // -------------------------------------------------------------
  //  interpolate fields and gradients
  //  Currently, this function will break if called with an unowned ielem.
  //  Currently, this function is only called with owned ielems.
  // -------------------------------------------------------------
  void FE_Engine::interpolate_fields(
    const int ielem,
    const FIELDS &fields,
    AliasArray<int> & conn,
    DENS_MAT & N,
    DENS_MAT_VEC & dN,
    DIAG_MAT & _weights_,
    FIELD_MATS &fieldsAtIPs,
    GRAD_FIELD_MATS &gradFieldsAtIPs) const
  {
    // evaluate shape functions
    feMesh_->shape_function(ielem, N, dN, _weights_);
    // get connectivity
    _conn_ = feMesh_->element_connectivity_unique(ielem);
    // compute fields and gradients of fields ips of this element
    FIELD_MATS localElementFields;
    for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++)  {
      // field values at all nodes
      _fieldName_ = _fieldItr_->first;
      const DENS_MAT &vI = (_fieldItr_->second).quantity();
      int dof = vI.nCols();
      // field values at integration points -> to be computed
      DENS_MAT &vP = fieldsAtIPs[_fieldName_];
      // gradients of field at integration points -> to be computed
      DENS_MAT_VEC &dvP = gradFieldsAtIPs[_fieldName_];
      
      if (_fieldName_ == ELECTRON_WAVEFUNCTION_ENERGIES ) {
        vP = vI;
        continue;
      }
      // field values at element nodes
      DENS_MAT &vIe = localElementFields[_fieldName_];
      // gather local field
      vIe.reset(nNodesPerElement_, dof);
      for (int i = 0; i < nNodesPerElement_; i++)
        for (int j = 0; j < dof; j++)
          vIe(i,j) = vI(conn(i),j);
      // interpolate field at integration points
      vP = N*vIe;
      // gradients
      dvP.assign(nSD_, DENS_MAT(nIPsPerElement_, dof));
      for (int j = 0; j < nSD_; ++j) dvP[j] = dN[j]*vIe;
    }
  }
  // -------------------------------------------------------------
  //  interpolate fields 
  //  Currently, this function will break if called with an unowned ielem.
  //  Currently, this function is only called with owned ielems.
  // -------------------------------------------------------------
  void FE_Engine::interpolate_fields(
    const int ielem,
    const FIELDS &fields,
    AliasArray<int> & conn,
    DENS_MAT & N,
    DIAG_MAT & _weights_,
    FIELD_MATS &fieldsAtIPs) const
  {
    // evaluate shape functions
    feMesh_->shape_function(ielem, N, _weights_);
    // get connectivity
    _conn_ = feMesh_->element_connectivity_unique(ielem);
    // compute fields and gradients of fields ips of this element
    FIELD_MATS localElementFields;
    for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++)  {
      // field values at all nodes
      _fieldName_ = _fieldItr_->first;
      const DENS_MAT &vI = (_fieldItr_->second).quantity();
      int dof = vI.nCols();
      // field values at integration points -> to be computed
      DENS_MAT &vP = fieldsAtIPs[_fieldName_];
      // field values at element nodes
      DENS_MAT &vIe = localElementFields[_fieldName_];
      
      if (_fieldName_ == ELECTRON_WAVEFUNCTION_ENERGIES ) {
        vP = vI;
        continue;
      }

      // gather local field
      vIe.reset(nNodesPerElement_, dof);
      for (int i = 0; i < nNodesPerElement_; i++)
        for (int j = 0; j < dof; j++)
          vIe(i,j) = vI(conn(i),j);
      // interpolate field at integration points
      vP = N*vIe;
    }
  }

  // -------------------------------------------------------------
  // compute dimensionless stiffness matrix using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::stiffness_matrix(SPAR_MAT &matrix) const
  {
    // assemble consistent mass 
    matrix.reset(nNodesUnique_,nNodesUnique_);// zero since partial fill
    DENS_MAT elementMassMatrix(nNodesPerElement_,nNodesPerElement_); 

    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // evaluate shape functions
      feMesh_->shape_function(ielem, _N_, _dN_, _weights_); // _N_ unused
      // perform quadrature
      elementMassMatrix = _dN_[0].transMat(_weights_*_dN_[0]); 
      for (int i = 1; i < nSD_; ++i) {
        elementMassMatrix += _dN_[i].transMat(_weights_*_dN_[i]); 
      }
      // get connectivity
      _conn_ = feMesh_->element_connectivity_unique(ielem);
      for (int i = 0; i < nNodesPerElement_; ++i) 
      {
        int inode = _conn_(i);
        for (int j = 0; j < nNodesPerElement_; ++j) 
        {
          int jnode = _conn_(j);
          matrix.add(inode, jnode, elementMassMatrix(i,j));
        }
      }
    }
#ifdef ISOLATE_FE
    sparse_allsum(communicator_,matrix);
#else
    LammpsInterface::instance()->sparse_allsum(matrix);
#endif
    matrix.compress();
  }

  // -------------------------------------------------------------
  // compute tangent using native quadrature for one (field,field) pair
  // -------------------------------------------------------------
  void FE_Engine::compute_tangent_matrix(
                                     const Array2D<bool> &rhsMask,
                                     const pair<FieldName,FieldName> row_col,
                                     const FIELDS &fields,
                                     const PhysicsModel * physicsModel,
                                     const Array<int>   & elementMaterials,
                                     SPAR_MAT &tangent,
                                     const DenseMatrix<bool> *elementMask) const
  {

    tangent.reset(nNodesUnique_,nNodesUnique_);
    FieldName rowField = row_col.first;
    FieldName colField = row_col.second;
    bool BB = rhsMask(rowField,FLUX);
    bool NN = rhsMask(rowField,SOURCE);
    DENS_MAT elementMatrix(nNodesPerElement_,nNodesPerElement_); 
    DENS_MAT coefsAtIPs;

    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // if element is masked, skip it
      if (elementMask && !(*elementMask)(ielem,0)) continue;
      // material id
      int imat = elementMaterials(ielem);
      const Material * mat = physicsModel->material(imat);
      // interpolate fields and gradients (nonlinear only)
      interpolate_fields(ielem,fields,_conn_,_N_,_dN_,_weights_,
        _fieldsAtIPs_,_gradFieldsAtIPs_);
        
      // evaluate Physics model
      
      if (! (physicsModel->null(rowField,imat)) ) {
        if (BB && physicsModel->weak_equation(rowField)-> 
          has_BB_tangent_coefficients() ) {
          physicsModel->weak_equation(rowField)->
            BB_tangent_coefficients(colField, _fieldsAtIPs_, mat, coefsAtIPs);
          DIAG_MAT D(column(coefsAtIPs,0)); 
          D = _weights_*D;
          elementMatrix = _dN_[0].transMat(D*_dN_[0]);
          for (int i = 1; i < nSD_; i++) {
            elementMatrix += _dN_[i].transMat(D*_dN_[i]);
          } 
        }
        else {
          elementMatrix.reset(nNodesPerElement_,nNodesPerElement_);
        }
        if (NN && physicsModel->weak_equation(rowField)->
          has_NN_tangent_coefficients() ) {
          physicsModel->weak_equation(rowField)->
            NN_tangent_coefficients(colField, _fieldsAtIPs_, mat, coefsAtIPs);
          DIAG_MAT D(column(coefsAtIPs,0)); 
          D = _weights_*D;
          elementMatrix += _N_.transMat(D*_N_);
        }
        // assemble
        for (int i = 0; i < nNodesPerElement_; ++i) 
        {  
          int inode = _conn_(i); 
          for (int j = 0; j < nNodesPerElement_; ++j)
          {
            int jnode = _conn_(j);
            tangent.add(inode, jnode, elementMatrix(i,j));
          }
        } 
      }
    } 
#ifdef ISOLATE_FE
    sparse_allsum(communicator_,tangent);
#else
    LammpsInterface::instance()->sparse_allsum(tangent);
#endif
    tangent.compress(); 
  }  

  // -------------------------------------------------------------
  // compute tangent using given quadrature for one (field,field) pair
  // -------------------------------------------------------------
  void FE_Engine::compute_tangent_matrix(const Array2D<bool> &rhsMask,
                            const pair<FieldName,FieldName> row_col,
                            const FIELDS        &fields,
                            const PhysicsModel * physicsModel,
                            const Array<set<int> > & pointMaterialGroups,
                            const DIAG_MAT      &weights,
                            const SPAR_MAT      &N,
                            const SPAR_MAT_VEC  &dN,
                            SPAR_MAT &tangent,
                            const DenseMatrix<bool> *elementMask ) const
  {
    int nn = nNodesUnique_;
    FieldName rowField = row_col.first;
    FieldName colField = row_col.second;
    bool BB = rhsMask(rowField,FLUX);
    bool NN = rhsMask(rowField,SOURCE);
    DENS_MAT K(nn,nn);
    DENS_MAT coefsAtIPs;
    int nips = weights.nCols();
    if (nips>0) {
      // compute fields and gradients of fields at given ips
      
      GRAD_FIELD_MATS gradFieldsAtIPs;
      FIELD_MATS fieldsAtIPs;
      for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
        _fieldName_          = _fieldItr_->first;
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int dof = field.nCols();
        gradFieldsAtIPs[_fieldName_].assign(nSD_,DENS_MAT(nips,dof));
        fieldsAtIPs[_fieldName_] = N*field;
        for (int j = 0; j < nSD_; ++j) {
          gradFieldsAtIPs[_fieldName_][j] = (*dN[j])*field;
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
      
      if (! singleMaterial ) throw ATC_Error("FE_Engine::compute_tangent_matrix-given quadrature can not handle multiple atom material currently");
      if (singleMaterial)
      {
        int imat = 0; 
        const Material *   mat = physicsModel->material(imat);
        // evaluate Physics model
        if (! (physicsModel->null(rowField,imat)) ) {
          if (BB && physicsModel->weak_equation(rowField)->
            has_BB_tangent_coefficients() ) {
            physicsModel->weak_equation(rowField)->
              BB_tangent_coefficients(colField, fieldsAtIPs, mat, coefsAtIPs);
            DIAG_MAT D(column(coefsAtIPs,0)); 
            D = weights*D;
            K = (*dN[0]).transMat(D*(*dN[0]));
            for (int i = 1; i < nSD_; i++) {
              K += (*dN[i]).transMat(D*(*dN[i]));
            } 
          }
          if (NN && physicsModel->weak_equation(rowField)->
            has_NN_tangent_coefficients() ) {
            physicsModel->weak_equation(rowField)->
              NN_tangent_coefficients(colField, fieldsAtIPs, mat, coefsAtIPs);
            DIAG_MAT D(column(coefsAtIPs,0)); 
            D = weights*D;
            K += N.transMat(D*N);
          }
        }
      }
    }
    // share information between processors
    int count = nn*nn; 
    DENS_MAT K_sum(nn,nn);
    allsum(communicator_, K.ptr(), K_sum.ptr(), count);
    // create sparse from dense
    tangent.reset(K_sum,tol_sparse);
    tangent.compress();
  }  

  // -------------------------------------------------------------
  //  compute a consistent mass matrix for native quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_mass_matrix(
     const Array<FieldName>& field_mask,
     const FIELDS &fields,
     const PhysicsModel * physicsModel,
     const Array<int>   & elementMaterials,
     
     CON_MASS_MATS & massMatrices,
     const DenseMatrix<bool> *elementMask) const
  {
    int nfields = field_mask.size();
    vector<FieldName> usedFields;

    
    DENS_MAT elementMassMatrix(nNodesPerElement_,nNodesPerElement_); 
    
    // (JAT, 04/21/11) FIX THIS
    DENS_MAT capacity;

    // zero, use incoming matrix as template if possible
    for (int j = 0; j < nfields; ++j) {
      _fieldName_ = field_mask(j);
      SPAR_MAT & M = massMatrices[_fieldName_].set_quantity();
      
      // compresses 2May11
      if   (M.has_template()) { M = 0; }
      else                    { M.reset(nNodesUnique_,nNodesUnique_); }
      M.reset(nNodesUnique_,nNodesUnique_);
    }

    // element wise assembly
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin(); 
         elemsIter != myElems.end();
         ++elemsIter) 
    {
      int ielem = *elemsIter;
      // if element is masked, skip 
      if (elementMask && !(*elementMask)(ielem,0)) continue;
      // material id
      int imat = elementMaterials(ielem);
      const Material * mat = physicsModel->material(imat);
      // interpolate fields
      interpolate_fields(ielem,fields,_conn_,_N_,_weights_,_fieldsAtIPs_);
      for (int j = 0; j < nfields; ++j) {
        _fieldName_ = field_mask(j);
        SPAR_MAT & M = massMatrices[_fieldName_].set_quantity();
        // skip null weakEqns by material
        if (! (physicsModel->null(_fieldName_,imat)) ) {
          physicsModel->weak_equation(_fieldName_)->
            M_integrand(_fieldsAtIPs_, mat, capacity);
          DIAG_MAT rho(column(capacity,0)); 
          elementMassMatrix = _N_.transMat(_weights_*rho*_N_);
          // assemble
          for (int i = 0; i < nNodesPerElement_; ++i) {  
            int inode = _conn_(i); 
            for (int j = 0; j < nNodesPerElement_; ++j) {
              int jnode = _conn_(j);
              M.add(inode, jnode, elementMassMatrix(i,j));
            }
          }
        }
      }
    }

    for (int j = 0; j < nfields; ++j) {
      _fieldName_ = field_mask(j);
      SPAR_MAT & M = massMatrices[_fieldName_].set_quantity();
#ifdef ISOLATE_FE
    sparse_allsum(communicator_,M);
#else
    LammpsInterface::instance()->sparse_allsum(M);
#endif
    }
    // fix zero diagonal entries due to null material elements
    for (int j = 0; j < nfields; ++j) {
      _fieldName_ = field_mask(j);
      
      SPAR_MAT & M = massMatrices[_fieldName_].set_quantity(); 
      for (int inode = 0; inode < nNodesUnique_; ++inode) {
        if (! M.has_entry(inode,inode)) {
          M.set(inode,inode,1.0); 
        } 
      }
      M.compress();
    }

  } 

  // -------------------------------------------------------------
  // compute dimensionless consistent mass using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_mass_matrix(SPAR_MAT &massMatrix) const
  {
    // assemble nnodes X nnodes matrix

    massMatrix.reset(nNodesUnique_,nNodesUnique_);// zero since partial fill
    DENS_MAT elementMassMatrix(nNodesPerElement_,nNodesPerElement_); 
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // evaluate shape functions
      feMesh_->shape_function(ielem, _N_, _weights_);
      // perform quadrature
      elementMassMatrix = _N_.transMat(_weights_*_N_); 
      // get connectivity
      _conn_ = feMesh_->element_connectivity_unique(ielem);
      for (int i = 0; i < nNodesPerElement_; ++i) {
        int inode = _conn_(i);
        for (int j = 0; j < nNodesPerElement_; ++j) {
          int jnode = _conn_(j);
          massMatrix.add(inode, jnode, elementMassMatrix(i,j));
        }
      }
    }
    // Assemble partial results
#ifdef ISOLATE_FE
    sparse_allsum(communicator_,massMatrix);
#else
    LammpsInterface::instance()->sparse_allsum(massMatrix);
#endif
  }

  // -------------------------------------------------------------
  // compute dimensionless consistent mass using given quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_mass_matrix(const DIAG_MAT &weights,
                                      const SPAR_MAT &N,
                                      SPAR_MAT &massMatrix) const
  {
    int nn   = N.nCols();
    int nips = N.nRows();
    DENS_MAT tmp_mass_matrix_local(nn,nn), tmp_mass_matrix(nn,nn);
    if (nips>0) { tmp_mass_matrix_local = N.transMat(weights*N); } 
    // share information between processors
    int count = nn*nn; 
    allsum(communicator_,
      tmp_mass_matrix_local.ptr(),
      tmp_mass_matrix.ptr(), count);
    // create sparse from dense
    massMatrix.reset(tmp_mass_matrix,tol_sparse);
  }

  // -------------------------------------------------------------
  // compute dimensionless lumped mass using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_lumped_mass_matrix(DIAG_MAT & M) const
  {
    M.reset(nNodesUnique_,nNodesUnique_); // zero since partial fill
    // assemble lumped diagonal mass 
    DENS_VEC Nvec(nNodesPerElement_);
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // evaluate shape functions
      feMesh_->shape_function(ielem, _N_, _weights_);
      CLON_VEC w(_weights_);
      // get connectivity
      _conn_ = feMesh_->element_connectivity_unique(ielem);
      Nvec = w*_N_;
      for (int i = 0; i < nNodesPerElement_; ++i) {
        int inode = _conn_(i);
        M(inode,inode) += Nvec(i);
      }
    }
    // Assemble partial results
    allsum(communicator_,MPI_IN_PLACE, M.ptr(), M.size());
  }

  // -------------------------------------------------------------
  // compute physical lumped mass using native quadrature 
  // -------------------------------------------------------------
  void FE_Engine::compute_lumped_mass_matrix(
    const Array<FieldName>& field_mask,
    const FIELDS          & fields,
    const PhysicsModel * physicsModel,
    const Array<int> &elementMaterials,
    MASS_MATS & massMatrices, // mass matrix
    const DenseMatrix<bool> *elementMask) const
  {
    int nfields = field_mask.size();
    // zero initialize for assembly
    for (int j = 0; j < nfields; ++j) { 
      DIAG_MAT & M = massMatrices[field_mask(j)].set_quantity();
      M.reset(nNodesUnique_,nNodesUnique_); 
    }
    // assemble diagonal matrix
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // if element is masked, skip it
      if (elementMask && !(*elementMask)(ielem,0)) continue;
      // material id
      int imat = elementMaterials(ielem);
      const Material * mat = physicsModel->material(imat);
      interpolate_fields(ielem,fields,_conn_,_N_,_weights_,_fieldsAtIPs_);
      // compute densities, integrate & assemble
      DENS_MAT capacity;
      for (int j = 0; j < nfields; ++j) {
        _fieldName_ = field_mask(j);
        if (! physicsModel->null(_fieldName_,imat)) { 
          physicsModel->weak_equation(_fieldName_)->
            M_integrand(_fieldsAtIPs_, mat, capacity);
          _Nmat_ = _N_.transMat(_weights_*capacity);
          DIAG_MAT & M = massMatrices[_fieldName_].set_quantity();
          for (int i = 0; i < nNodesPerElement_; ++i) {
            int inode = _conn_(i);  
            M(inode,inode) += _Nmat_(i,0);

          }
        }
      }
    }
    // Assemble partial results
    for (int j = 0; j < nfields; ++j) {
      _fieldName_ = field_mask(j);
      DIAG_MAT & M = massMatrices[_fieldName_].set_quantity();
      allsum(communicator_,MPI_IN_PLACE, M.ptr(), M.size());
    }

    // fix zero diagonal entries due to null material elements
    for (int inode = 0; inode < nNodesUnique_; ++inode) {
      for (int j = 0; j < nfields; ++j) {
        _fieldName_ = field_mask(j);
        DIAG_MAT & M = massMatrices[_fieldName_].set_quantity();
        if (M(inode,inode) == 0.0) { 
          M(inode,inode) = 1.0; 
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
    const DIAG_MAT      &weights, 
    const SPAR_MAT      &N,
    MASS_MATS &M) const // mass matrices 
  {
    int nips = weights.nCols();
    int nfields = field_mask.size();
    // initialize
    map<FieldName, DIAG_MAT> M_local;
    for (int j = 0; j < nfields; ++j) {
      _fieldName_ = field_mask(j);
      M_local[_fieldName_].reset(nNodesUnique_,nNodesUnique_);
      M[_fieldName_].reset(nNodesUnique_,nNodesUnique_);
    }
    
    if (nips>0) {
      // compute fields at given ips
      // compute all fields since we don't the capacities dependencies
      FIELD_MATS fieldsAtIPs;
      for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
        _fieldName_ = _fieldItr_->first;
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int dof = field.nCols();
        fieldsAtIPs[_fieldName_].reset(nips,dof);
        fieldsAtIPs[_fieldName_] = N*field;
      }
    
      // treat single material point sets specially
      int nMatls = pointMaterialGroups.size();
      int atomMatls = 0;
      int iatomMat = 0;
      for (int imat = 0; imat < nMatls; imat++) {
        const set<int> & indices = pointMaterialGroups(imat);
        if ( indices.size() > 0) {
          iatomMat = imat;
          atomMatls++;
        }
      }
      if (atomMatls == 0) {
        throw ATC_Error("no materials in atom region - atoms may have migrated to FE-only region");
      }
      bool singleMaterial = ( atomMatls == 1 );
      if (! singleMaterial ) {
        stringstream ss; ss << " WARNING: multiple materials in atomic region";
       print_msg(communicator_,ss.str());
      }
    
      // setup data structures
      FIELD_MATS capacities;
      // evaluate physics model & compute capacity|density for requested fields
      if (singleMaterial) {
        const Material * mat = physicsModel->material(iatomMat);
        for (int j = 0; j < nfields; ++j) {
          _fieldName_ = field_mask(j);
          // mass matrix needs to be invertible so null matls have cap=1
          if (physicsModel->null(_fieldName_,iatomMat)) {
             throw ATC_Error("null material not supported for atomic region (single material)");
             const FIELD &  f = (fields.find(_fieldName_))->second;
             capacities[_fieldName_].reset(f.nRows(),f.nCols());
             capacities[_fieldName_] = 1.;
             
          } 
          else {
            physicsModel->weak_equation(_fieldName_)->
              M_integrand(fieldsAtIPs, mat, capacities[_fieldName_]);
          }
        }
      }
      else {
        FIELD_MATS groupCapacities, fieldsGroup;
        for (int j = 0; j < nfields; ++j) {
          _fieldName_ = field_mask(j);
          capacities[_fieldName_].reset(nips,1);
        }
        for ( int imat = 0; imat < pointMaterialGroups.size(); imat++) {
          const Material * mat = physicsModel->material(imat);
          const set<int> & indices = pointMaterialGroups(imat);
          int npts = indices.size();
          if (npts > 0) {
            for (int j = 0; j < nfields; ++j) {
              _fieldName_ = field_mask(j);
              groupCapacities[_fieldName_].reset(npts,1);
              const FIELD &  f = (fields.find(_fieldName_))->second;
              int ndof = f.nCols();
              fieldsGroup[_fieldName_].reset(npts,ndof);
              int i = 0;
              for (set<int>::const_iterator iter=indices.begin();
                   iter != indices.end(); iter++, i++ ) {
                for (int dof = 0; dof < ndof; ++dof) {
                  fieldsGroup[_fieldName_](i,dof) 
                    = fieldsAtIPs[_fieldName_](*iter,dof);
                }
              }
            }
            for (int j = 0; j < nfields; ++j) {
              _fieldName_ = field_mask(j);
              if (physicsModel->null(_fieldName_,imat)) {
                throw ATC_Error("null material not supported for atomic region (multiple materials)");
                const FIELD &  f = (fields.find(_fieldName_))->second;
                groupCapacities[_fieldName_].reset(f.nRows(),f.nCols());
                groupCapacities[_fieldName_] = 1.;
              } 
              else {
                physicsModel->weak_equation(_fieldName_)->
                  M_integrand(fieldsGroup, mat, groupCapacities[_fieldName_]);
              }
              int i = 0;
              for (set<int>::const_iterator iter=indices.begin();
                   iter != indices.end(); iter++, i++ ) {
                capacities[_fieldName_](*iter,0) 
                  = groupCapacities[_fieldName_](i,0);
              }
            }
          }
        }
      }
      
      // integrate & assemble
      for (int j = 0; j < nfields; ++j) {
        _fieldName_ = field_mask(j);
        M_local[_fieldName_].reset( // assume all columns same
          column(N.transMat(weights*capacities[_fieldName_]),0) );
      }
    }
    
    // Share information between processors
    for (int j = 0; j < nfields; ++j) {
      _fieldName_ = field_mask(j);
      DIAG_MAT & myMassMat(M[_fieldName_].set_quantity());
      int count = M_local[_fieldName_].size();
      allsum(communicator_, M_local[_fieldName_].ptr(), myMassMat.ptr(), count);
    }
  }
  
  //-----------------------------------------------------------------
  // compute assembled average gradient evaluated at the nodes
  //-----------------------------------------------------------------
  void FE_Engine::compute_gradient_matrix(SPAR_MAT_VEC & grad_matrix) const
  {
    // zero
    DENS_MAT_VEC tmp_grad_matrix(nSD_);  
    for (int i = 0; i < nSD_; i++) {
      tmp_grad_matrix[i].reset(nNodesUnique_,nNodesUnique_); 
    }
    // element wise assembly
    Array<int> count(nNodesUnique_); count = 0;
    set_quadrature(NODAL);
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // evaluate shape functions
      feMesh_->shape_function(ielem, _N_, _dN_, _weights_);
      // get connectivity
      _conn_ = feMesh_->element_connectivity_unique(ielem);
      for (int j = 0; j < nIPsPerElement_; ++j) {
        int jnode = _conn_(j); 
        count(jnode) += 1;
        for (int i = 0; i < nNodesPerElement_; ++i) {
          int inode = _conn_(i);
          for (int k = 0; k < nSD_; ++k) {
            tmp_grad_matrix[k](jnode,inode) += _dN_[k](j,i);

          }
        }
      }
    }
    // Assemble partial results
    for (int k = 0; k < nSD_; ++k) {
      allsum(communicator_,MPI_IN_PLACE, tmp_grad_matrix[k].ptr(), tmp_grad_matrix[k].size());
    }
    int_allsum(communicator_,MPI_IN_PLACE, count.ptr(), count.size());
    set_quadrature(quadrature_); //reset to default
    for (int inode = 0; inode < nNodesUnique_; ++inode) {
      for (int jnode = 0; jnode < nNodesUnique_; ++jnode) {
        for (int k = 0; k < nSD_; ++k) {
          tmp_grad_matrix[k](jnode,inode) /= count(jnode);
        }
      }
    }
    // compact dense matrices
    for (int k = 0; k < nSD_; ++k) {
      grad_matrix[k]->reset(tmp_grad_matrix[k],tol_sparse);
    }
  }

  // -------------------------------------------------------------
  // compute energy per node using native quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_energy(const Array<FieldName> &mask,
    const FIELDS &fields,
    const PhysicsModel * physicsModel,
    const Array<int>   & elementMaterials,
    FIELD_MATS &energies,
    const DenseMatrix<bool> *elementMask,
    const IntegrationDomainType domain) const
  {
    // Zero out all fields
    for (int n = 0; n < mask.size(); n++) {
      _fieldName_ = mask(n);
      _fieldItr_ = fields.find(_fieldName_);
      const DENS_MAT & field = (_fieldItr_->second).quantity();
      energies[_fieldName_].reset(field.nRows(), 1);
    }

    DENS_MAT elementEnergy(nNodesPerElement_,1); // workspace
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // if element is masked, skip it
      if (domain != FULL_DOMAIN && elementMask && !(*elementMask)(ielem,0)) continue;
      // material id
      int imat = elementMaterials(ielem);
      const Material * mat = physicsModel->material(imat);
      interpolate_fields(ielem,fields,_conn_,_N_,_dN_,_weights_,
        _fieldsAtIPs_,_gradFieldsAtIPs_);
      
      // assemble
      for (int n = 0; n < mask.size(); n++) {
        _fieldName_ = mask(n);
        if (physicsModel->null(_fieldName_,imat)) continue;
        if( ! (physicsModel->weak_equation(_fieldName_)-> has_E_integrand()))  continue;
        physicsModel->weak_equation(_fieldName_)->
          E_integrand(_fieldsAtIPs_, _gradFieldsAtIPs_, mat, elementEnergy);
        _fieldItr_ = fields.find(_fieldName_);
        _Nmat_ = _N_.transMat(_weights_*elementEnergy);
        DENS_MAT & energy = energies[_fieldName_];
        for (int i = 0; i < nNodesPerElement_; ++i) {  
          int inode = _conn_(i);  
          energy(inode,0) += _Nmat_(i,0);
        }
      }
    }

    for (int n = 0; n < mask.size(); n++) {
      _fieldName_ = mask(n);
      DENS_MAT& myEnergy(energies[_fieldName_]);
      allsum(communicator_,MPI_IN_PLACE, myEnergy.ptr(), myEnergy.size());
    }
  }

  // -------------------------------------------------------------
  // compute rhs using native quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_rhs_vector(const Array2D<bool> &rhsMask,
    const FIELDS &fields,
    const PhysicsModel * physicsModel,
    const Array<int>   & elementMaterials,
    FIELDS &rhs,
    bool freeOnly,
    const DenseMatrix<bool> *elementMask) const
  {
    vector<FieldName> usedFields;


    // size and zero output
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if (rhsMask(_fieldName_, FLUX) || rhsMask(_fieldName_, SOURCE)) {
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        rhs[_fieldName_].reset(field.nRows(), field.nCols());

        // Save field names for easy lookup later.
        usedFields.push_back(_fieldName_);
      }
    }

    // Iterate over elements partitioned onto current processor.
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin(); 
         elemsIter != myElems.end();
         ++elemsIter) 
    {
      int ielem = *elemsIter;
      // Skip masked elements 
      if (elementMask && !(*elementMask)(ielem,0)) continue;
      int imat = elementMaterials(ielem); // material id
      const Material * mat = physicsModel->material(imat);
      
      // interpolate field values to integration points
      interpolate_fields(ielem,fields,_conn_,_N_,_dN_,_weights_,
        _fieldsAtIPs_,_gradFieldsAtIPs_);
      
      // rescale by _weights_, a diagonal matrix
      _Nw_ = _weights_*_N_;
      for (int j = 0; j < nSD_; ++j) _dNw_[j] = _weights_*_dN_[j];
      
      // evaluate physics model and assemble
      // _Nfluxes is a set of fields
      for (int n = 0; n < rhsMask.nRows(); n++) {
        _fieldName_ = FieldName(n);
        if (!rhsMask(_fieldName_,SOURCE)) continue;
        if (physicsModel->null(_fieldName_,imat)) continue;

        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        DENS_MAT & myRhs(rhs[_fieldName_].set_quantity());
       
        int dof = field.nCols();
        bool has = physicsModel->weak_equation(_fieldName_)->
          N_integrand(_fieldsAtIPs_,_gradFieldsAtIPs_, mat, _Nfluxes_);
        if (!has) continue;
        _Nmat_ = _Nw_.transMat(_Nfluxes_);
        
        for (int i = 0; i < nNodesPerElement_; ++i) {  
          int inode = _conn_(i);  
          for (int k = 0; k < dof; ++k) {
            myRhs(inode,k) += _Nmat_(i,k);
          }
        }
      }
      // _Bfluxes_ is a set of field gradients
      for (int n = 0; n < rhsMask.nRows(); n++) {
        _fieldName_ = FieldName(n);
        if (!rhsMask(_fieldName_,FLUX)) continue;
        if (physicsModel->null(_fieldName_,imat)) continue;
        
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        DENS_MAT & myRhs(rhs[_fieldName_].set_quantity());
         
        int dof = field.nCols();
        physicsModel->weak_equation(_fieldName_)->
          B_integrand(_fieldsAtIPs_, _gradFieldsAtIPs_, mat, _Bfluxes_);

        for (int j = 0; j < nSD_; j++) {
          _Nmat_ = _dNw_[j].transMat(_Bfluxes_[j]);
          for (int i = 0; i < nNodesPerElement_; ++i) {  
            int inode = _conn_(i); 
            for (int k = 0; k < dof; ++k) {
              myRhs(inode,k) += _Nmat_(i,k);
            } 
          }
        }
      } 
    }

    vector<FieldName>::iterator _fieldIter_;
    for (_fieldIter_ = usedFields.begin(); _fieldIter_ != usedFields.end(); 
         ++_fieldIter_) {
      // myRhs is where we put the global result for this field.
      _fieldName_ = *_fieldIter_;
      DENS_MAT & myRhs(rhs[_fieldName_].set_quantity());

      // Sum matrices in-place across all processors into myRhs.
      allsum(communicator_,MPI_IN_PLACE, myRhs.ptr(), myRhs.size());
    }
  }  

  // -------------------------------------------------------------
  // compute rhs using given quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_rhs_vector(const Array2D<bool> &rhsMask,
    const FIELDS        &fields,
    const PhysicsModel  *physicsModel,
    const Array<set<int> >&pointMaterialGroups,
    const DIAG_MAT      &weights,
    const SPAR_MAT      &N,
    const SPAR_MAT_VEC  &dN,
    FIELDS              &rhs) const
  {
    FIELD_MATS rhs_local;
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if (rhsMask(_fieldName_,FLUX) || rhsMask(_fieldName_,SOURCE)) {
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int nrows = field.nRows();
        int ncols = field.nCols();
        rhs      [_fieldName_].reset(nrows,ncols);
        rhs_local[_fieldName_].reset(nrows,ncols);
      }
    }
    
    int nips = weights.nCols();

    if (nips>0) {
      // compute fields and gradients of fields at given ips
      
      GRAD_FIELD_MATS gradFieldsAtIPs;
      FIELD_MATS fieldsAtIPs;
      for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
        _fieldName_          = _fieldItr_->first;
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int dof = field.nCols();
        gradFieldsAtIPs[_fieldName_].assign(nSD_,DENS_MAT(nips,dof));
        fieldsAtIPs[_fieldName_] = N*field;
        for (int j = 0; j < nSD_; ++j) {
          gradFieldsAtIPs[_fieldName_][j] = (*dN[j])*field;
        }
      }

      // setup data structures
      FIELD_MATS Nfluxes; 
      GRAD_FIELD_MATS Bfluxes; 
      for (int n = 0; n < rhsMask.nRows(); n++) {
        _fieldName_ = FieldName(n);
        int index = (int) _fieldName_;
        if ( rhsMask(index,FLUX) ) {
          Bfluxes[_fieldName_].assign(nSD_, DENS_MAT()); 
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
      
      // evaluate physics model
      if (singleMaterial) 
      {
        int imat = 0;
        const Material *   mat = physicsModel->material(imat);
        for (int n = 0; n < rhsMask.nRows(); n++) {
          _fieldName_ = FieldName(n);
          int index = (int) _fieldName_;
          if (! physicsModel->null(_fieldName_,imat)) {
            if ( rhsMask(index,SOURCE) ) {
              physicsModel->weak_equation(_fieldName_)->
                N_integrand(fieldsAtIPs, gradFieldsAtIPs, mat, Nfluxes[_fieldName_]);
            }
            if ( rhsMask(index,FLUX) ) {
              physicsModel->weak_equation(_fieldName_)->
                B_integrand(fieldsAtIPs, gradFieldsAtIPs, mat, Bfluxes[_fieldName_]);
            }
          }
        }
      }
      else 
      {
        
        
        // 1) permanent workspace with per-row mapped clones per matl
        //    from caller/atc
        // 2) per point calls to N/B_integrand
        // 3) collect internalToAtom into materials and use mem-cont clones
        //    what about dof for matrices and data storage: clone _matrix_
        
        // for each material group: 
        // set up storage
        DENS_MAT        group_Nfluxes;
        DENS_MAT_VEC    group_Bfluxes; 
        group_Bfluxes.reserve(nSD_); 
        FIELD_MATS      fieldsGroup;
        GRAD_FIELD_MATS gradFieldsGroup;
        for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
          _fieldName_          = _fieldItr_->first;
          const DENS_MAT & field = (_fieldItr_->second).quantity();
          int ndof = field.nCols();
          gradFieldsGroup[_fieldName_].assign(nSD_,DENS_MAT(nips,ndof));
          
          Nfluxes[_fieldName_].reset(nips,ndof);
          Bfluxes[_fieldName_].assign(nSD_,DENS_MAT(nips,ndof));
          //}
        }

        // copy fields
        for ( int imat = 0; imat < pointMaterialGroups.size(); imat++) 
        {
          const set<int> & indices = pointMaterialGroups(imat);
          const Material *   mat = physicsModel->material(0);
          int npts = indices.size();
          int i = 0;
          for (set<int>::const_iterator iter=indices.begin();
                                        iter != indices.end(); iter++, i++ ) {
            for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
              _fieldName_          = _fieldItr_->first;
              const DENS_MAT & field = (_fieldItr_->second).quantity();
              int ndof = field.nCols();
              fieldsGroup[_fieldName_].reset(npts,ndof);
              for (int j = 0; j < nSD_; ++j) {
                (gradFieldsGroup[_fieldName_][j]).reset(npts,ndof);
              }
              for (int dof = 0; dof < ndof; ++dof) {
                fieldsGroup[_fieldName_](i,dof) 
                  = fieldsAtIPs[_fieldName_](*iter,dof);
                for (int j = 0; j < nSD_; ++j) {
                  gradFieldsGroup[_fieldName_][j](i,dof) 
                    = gradFieldsAtIPs[_fieldName_][j](*iter,dof);
                }
              }
            }
          }
          // calculate forces, & assemble 
          for (int n = 0; n < rhsMask.nRows(); n++) {
            _fieldName_ = FieldName(n);
            _fieldItr_ = fields.find(_fieldName_);
            int index = (int) _fieldName_;
            if (physicsModel->null(_fieldName_,imat)) continue;
            if ( rhsMask(index,SOURCE) ) {
              const DENS_MAT & field = (_fieldItr_->second).quantity();
              int ndof = field.nCols();
              bool has = physicsModel->weak_equation(_fieldName_)->
                N_integrand(fieldsGroup, gradFieldsGroup, mat, group_Nfluxes);
              if (! has) throw ATC_Error("atomic source can not be null currently");
              int i = 0;
              for (set<int>::const_iterator iter=indices.begin();
                                   iter != indices.end(); iter++, i++ ) {
                for (int dof = 0; dof < ndof; ++dof) {
                  Nfluxes[_fieldName_](*iter,dof) += group_Nfluxes(i,dof);
                }
              }
            }
          }
          // calculate forces, & assemble 
          for (int n = 0; n < rhsMask.nRows(); n++) {
            _fieldName_ = FieldName(n);
            if (physicsModel->null(_fieldName_,imat)) continue;
            int index = (int) _fieldName_;
            if ( rhsMask(index,FLUX) ) {
              _fieldItr_ = fields.find(_fieldName_);
              const DENS_MAT & field = (_fieldItr_->second).quantity();   
              int ndof = field.nCols();
              physicsModel->weak_equation(_fieldName_)->
                B_integrand(fieldsGroup, gradFieldsGroup, mat, group_Bfluxes);
              int i = 0;
              for (set<int>::const_iterator iter=indices.begin();
                                   iter != indices.end(); iter++, i++ ) {
                for (int dof = 0; dof < ndof; ++dof) {
                  for (int j = 0; j < nSD_; ++j) {
                    Bfluxes[_fieldName_][j](*iter,dof) += group_Bfluxes[j](i,dof);
                  }
                }
              }
            }
          }
        }
      } // endif multiple materials
      
      // assemble : N/Bfluxes -> rhs_local
      for (int n = 0; n < rhsMask.nRows(); n++) {
        _fieldName_ = FieldName(n);
        int index = (int) _fieldName_;
        if ( rhsMask(index,SOURCE) && Nfluxes[_fieldName_].nCols() > 0 ) {
          rhs_local[_fieldName_] 
             += N.transMat(weights*Nfluxes[_fieldName_]);
        }
        if ( rhsMask(index,FLUX) && (Bfluxes[_fieldName_][0]).nCols() > 0 ) {
          for (int j = 0; j < nSD_; ++j) {
            rhs_local[_fieldName_] 
              += dN[j]->transMat(weights*Bfluxes[_fieldName_][j]);
          }
        }
      }
    } // end nips > 0
    
    // Share information between processors: rhs_local -> rhs
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if (rhsMask(_fieldName_,FLUX) || rhsMask(_fieldName_,SOURCE)) {
        DENS_MAT & myRhs(rhs[_fieldName_].set_quantity());
        int count = rhs_local[_fieldName_].size();
        allsum(communicator_, rhs_local[_fieldName_].ptr(), myRhs.ptr(), count);
      }
    }
  }

  // -------------------------------------------------------------
  // compute sources using given quadrature
  // -------------------------------------------------------------
  void FE_Engine::compute_source(const Array2D<bool> &rhsMask,
    const FIELDS        &fields,
    const PhysicsModel  *physicsModel,
    const Array<set<int> >&pointMaterialGroups,
    const DIAG_MAT      &weights,
    const SPAR_MAT      &N,
    const SPAR_MAT_VEC  &dN,
    FIELD_MATS          &sources) const
  {
    int nips = weights.nCols();


    if (nips>0) {
      FIELD_MATS  Nfluxes;
      // compute fields and gradients of fields at given ips
      
      GRAD_FIELD_MATS gradFieldsAtIPs;
      FIELD_MATS fieldsAtIPs;
      for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
        _fieldName_          = _fieldItr_->first;
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int dof = field.nCols();
        gradFieldsAtIPs[_fieldName_].assign(nSD_,DENS_MAT(nips,dof));
        fieldsAtIPs[_fieldName_] = N*field;
        for (int j = 0; j < nSD_; ++j) {
          gradFieldsAtIPs[_fieldName_][j] = (*dN[j])*field;
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
      
      
      // evaluate physics model
      if (singleMaterial) 
      {
        int imat = 0;
        const Material *   mat = physicsModel->material(imat);
        for (int n = 0; n < rhsMask.nRows(); n++) {
          _fieldName_ = FieldName(n);
          int index = (int) _fieldName_;
          if (! physicsModel->null(_fieldName_,imat)) {
            if ( rhsMask(index,SOURCE) ) {
              bool has = physicsModel->weak_equation(_fieldName_)->
                N_integrand(fieldsAtIPs, gradFieldsAtIPs, mat, Nfluxes[_fieldName_]);
              if (! has) throw ATC_Error("atomic source can not be null currently");
            }
          }
        }
      }
      else 
      {
        throw ATC_Error("compute source does not handle multiple materials currently");
      }
      
      // assemble
      for (int n = 0; n < rhsMask.nRows(); n++) {
        _fieldName_ = FieldName(n);
        int index = (int) _fieldName_;
        if ( rhsMask(index,SOURCE) ) {
          sources[_fieldName_] =weights*Nfluxes[_fieldName_];
        }
      }
    }
    
    // no need to share information between processors
  }

  // -------------------------------------------------------------
  // compute flux for post processing
  // -------------------------------------------------------------
  void FE_Engine::compute_flux(const Array2D<bool> &rhsMask,
    const FIELDS &fields,
    const PhysicsModel * physicsModel,
    const Array<int>   & elementMaterials,
    GRAD_FIELD_MATS &fluxes,
    const DenseMatrix<bool> *elementMask) const
  {
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if (rhsMask(_fieldName_,FLUX)) {
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();

        fluxes[_fieldName_].assign(nSD_, DENS_MAT(field.nRows(),field.nCols()));
      }
    }

    // element wise assembly
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      // if element is masked, skip it
      if (elementMask && !(*elementMask)(ielem,0)) continue;
      // material id
      int imat = elementMaterials(ielem);
      const Material * mat = physicsModel->material(imat);
      interpolate_fields(ielem,fields,_conn_,_N_,_dN_,_weights_,
        _fieldsAtIPs_,_gradFieldsAtIPs_);
      _Nw_ = _weights_*_N_;
      // evaluate Physics model & assemble
      for (int n = 0; n < rhsMask.nRows(); n++) {
        _fieldName_ = FieldName(n);
        if (!rhsMask(_fieldName_,FLUX)) continue;
        if (physicsModel->null(_fieldName_,imat)) continue;
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int dof = field.nCols();
        physicsModel->weak_equation(_fieldName_)->
          B_integrand(_fieldsAtIPs_, _gradFieldsAtIPs_, mat, _Bfluxes_);
        for (int j = 0; j < nSD_; j++) {
          _Nmat_ = _Nw_.transMat(_Bfluxes_[j]);
          for (int i = 0; i < nNodesPerElement_; ++i) {  
            int inode = _conn_(i); 
            for (int k = 0; k < dof; ++k) {
              fluxes[_fieldName_][j](inode,k) += _Nmat_(i,k); 
            } 
          } 
        }
      } 
    } 
    // Assemble partial results
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if (!rhsMask(_fieldName_,FLUX)) continue;
      for (int j = 0; j < nSD_; j++) {
        allsum(communicator_,MPI_IN_PLACE, fluxes[_fieldName_][j].ptr(), fluxes[_fieldName_][j].size());  
      } 
    }
  }  

  //-----------------------------------------------------------------
  // boundary flux using native quadrature
  //-----------------------------------------------------------------
  void FE_Engine::compute_boundary_flux(const Array2D<bool> & rhsMask,
    const FIELDS        & fields,
    const PhysicsModel  * physicsModel,
    const Array<int>    & elementMaterials,
    const set< pair <int,int> > &  faceSet,
    FIELDS & rhs) const
  {
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if (rhsMask(_fieldName_,FLUX)) {
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        rhs[_fieldName_].reset(field.nRows(),field.nCols());
      }
    }

    FIELD_MATS localElementFields; 
    GRAD_FIELD_MATS gradFieldsAtIPs; 
    FIELD_MATS fieldsAtIPs;

    for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
      _fieldName_ = _fieldItr_->first;
      const FIELD & field = _fieldItr_->second;
      int dof = field.nCols();
      gradFieldsAtIPs[_fieldName_].reserve(nSD_);
      for (int i = 0; i < nSD_; ++i) {
        gradFieldsAtIPs[_fieldName_].push_back(DENS_MAT(nIPsPerFace_,dof));
      }
      fieldsAtIPs[_fieldName_].reset(nIPsPerFace_,dof);
      localElementFields[_fieldName_].reset(nNodesPerElement_,dof);
    }
  
    // element wise assembly
    set< PAIR >::iterator iter;
    for (iter = faceSet.begin(); iter != faceSet.end(); iter++) {
      // get connectivity
      int ielem = iter->first;
      // if this is not our element, do not do calculations
      if (!feMesh_->is_owned_elt(ielem)) continue;
      int imat = elementMaterials(ielem);
      const Material * mat = physicsModel->material(imat);
      _conn_ = feMesh_->element_connectivity_unique(ielem);

      // evaluate shape functions at ips
      feMesh_->face_shape_function(*iter, _fN_, _fdN_, _nN_, _fweights_);

      // interpolate fields and gradients of fields ips of this element

      for (_fieldItr_ = fields.begin(); _fieldItr_ != fields.end(); _fieldItr_++) {
        _fieldName_          = _fieldItr_->first;
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int dof = field.nCols();
        
        for (int i = 0; i < nNodesPerElement_; i++) {
          for (int j = 0; j < dof; j++) {
            localElementFields[_fieldName_](i,j) = field(_conn_(i),j);
          }
        }
        //  ips X dof    = ips X nodes * nodes X dof
        fieldsAtIPs[_fieldName_] = _fN_*localElementFields[_fieldName_];
        for (int j = 0; j < nSD_; ++j) {
          gradFieldsAtIPs[_fieldName_][j] = _fdN_[j]*localElementFields[_fieldName_];
        }
      }

      // Evaluate- physics model
      // do nothing for N_integrand
      // nN : precomputed and held by ATC_Transfer
      // assemble
      for (int n = 0; n < rhsMask.nRows(); n++) {
        _fieldName_ = FieldName(n);
        int index = (int) _fieldName_;
        if ( rhsMask(index,FLUX) ) {
          _fieldItr_ = fields.find(_fieldName_);
          const DENS_MAT & field = (_fieldItr_->second).quantity();
          DENS_MAT & myRhs(rhs[_fieldName_].set_quantity());
          physicsModel->weak_equation(_fieldName_)->
            B_integrand(fieldsAtIPs, gradFieldsAtIPs, mat, _Bfluxes_);
          int dof = field.nCols();
          for (int j = 0; j < nSD_; j++) {
            _Nmat_ = _nN_[j].transMat(_fweights_*_Bfluxes_[j]);
            for (int i = 0; i < nNodesPerElement_; ++i) {
              int inode = _conn_(i);
              for (int k = 0; k < dof; ++k) {
                myRhs(inode,k) += _Nmat_(i,k);
              }
            }
          }
        }
      }
    }
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      int index = (int) _fieldName_;
      if ( rhsMask(index,FLUX) ) {
        DENS_MAT & myRhs(rhs[_fieldName_].set_quantity());
        allsum(communicator_,MPI_IN_PLACE, myRhs.ptr(), myRhs.size());
      }
    }
  }

  // -------------------------------------------------------------
  // compute boundary flux using given quadrature and interpolation
  // -------------------------------------------------------------
  void FE_Engine::compute_boundary_flux(const Array2D<bool>      & rhsMask,
    const FIELDS             & fields,
    const PhysicsModel       * physicsModel,
    const Array< int >       & elementMaterials,
    const Array< set<int> >  & pointMaterialGroups,
    const DIAG_MAT           & _weights_,
    const SPAR_MAT           & N,
    const SPAR_MAT_VEC       & dN,
    const DIAG_MAT           & flux_mask,
    FIELDS                   & flux,
    const DenseMatrix<bool>  * elementMask,
    const set<int>           * nodeSet) const
  {
    
    
    FIELDS rhs, rhsSubset;
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if (rhsMask(_fieldName_,FLUX)) {
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        int nrows = field.nRows();
        int ncols = field.nCols();
        rhs      [_fieldName_].reset(nrows,ncols);
        rhsSubset[_fieldName_].reset(nrows,ncols);
      }
    }
    
    // R_I    = - int_Omega    Delta _N_I .q dV
    compute_rhs_vector(rhsMask, fields, physicsModel, elementMaterials, rhs, elementMask);

    // R_I^md = - int_Omega^md Delta _N_I .q dV
    compute_rhs_vector(rhsMask, fields, physicsModel, pointMaterialGroups, 
      _weights_, N, dN, rhsSubset);

    // flux_I = 1/Delta V_I V_I^md R_I + R_I^md
    //   where : Delta V_I = int_Omega _N_I dV
    for (int n = 0; n < rhsMask.nRows(); n++) {
      _fieldName_ = FieldName(n);
      if ( rhsMask(_fieldName_,FLUX) ) { 
        if (nodeSet) {
          _fieldItr_ = fields.find(_fieldName_);
          const DENS_MAT & field = (_fieldItr_->second).quantity();
          int nrows = field.nRows();
          int ncols = field.nCols();
          DENS_MAT & myFlux(flux[_fieldName_].set_quantity());
          const DENS_MAT & myRhsSubset(rhsSubset[_fieldName_].quantity());
          const DENS_MAT & myRhs(rhs[_fieldName_].quantity());
          myFlux.reset(nrows,ncols);
          set<int>::const_iterator iset;
          for (iset = nodeSet->begin(); iset != nodeSet->end(); iset++) {
            for (int j = 0; j < ncols; j++) {
              myFlux(*iset,j) = myRhsSubset(*iset,j) - flux_mask(*iset,*iset)*myRhs(*iset,j);
            }
          }
        }
        else {
          flux[_fieldName_] 
            = rhsSubset[_fieldName_].quantity() - flux_mask*(rhs[_fieldName_].quantity());
        }
      }
    }
  }

  /** integrate a nodal field over an element set */
  DENS_VEC FE_Engine::integrate(const DENS_MAT  &field, const ESET & eset) const
  {
    int dof = field.nCols();
    DENS_MAT eField(nNodesPerElement_, dof);
    int nips = nIPsPerElement_; 
    DENS_MAT ipField(nips, dof);
    DENS_VEC integral(dof); integral = 0;
    for (ESET::const_iterator itr = eset.begin(); itr != eset.end(); ++itr) {
      int ielem = *itr;
      // if this is not our element, do not do calculations
      if (!feMesh_->is_owned_elt(ielem)) continue;
      feMesh_->shape_function(ielem,_N_, _weights_);
      _conn_ = feMesh_->element_connectivity_unique(ielem);
      for (int i = 0; i < nNodesPerElement_; i++) {
        for (int j = 0; j < dof; j++) {
          eField(i,j) = field(_conn_(i),j); }}
      ipField = _N_*eField;
      for (int i = 0; i < nips; ++i) {  
        for (int j = 0; j < dof; ++j) {
          integral(j) += ipField(i,j)*_weights_[i]; 
        }
      }
    } 
    // assemble partial results
    allsum(communicator_,MPI_IN_PLACE, integral.ptr(), integral.size());
    return integral;
  }

  //-----------------------------------------------------------------
  // Robin boundary flux using native quadrature 
  // integrate \int_delV _N_I s(x,t).n dA
  //-----------------------------------------------------------------
  void FE_Engine::add_robin_fluxes(const Array2D<bool> &rhsMask,
    const FIELDS        & fields,
    const double time,
    const ROBIN_SURFACE_SOURCE & sourceFunctions,
    FIELDS &nodalSources) const
  {
    
    // sizing working arrays
    DENS_MAT xCoords(nSD_,nNodesPerElement_);
    DENS_MAT faceSource;
    DENS_MAT localField;
    DENS_MAT xAtIPs(nSD_,nIPsPerFace_);
    DENS_MAT uAtIPs(nIPsPerFace_,1); 
    FIELDS myNodalSources;

    // element wise assembly
    ROBIN_SURFACE_SOURCE::const_iterator src_iter;
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++)
    {
      _fieldName_ = src_iter->first;
      if (!rhsMask((int)_fieldName_,ROBIN_SOURCE)) continue; 
      DENS_MAT & s(nodalSources[_fieldName_].set_quantity());
      myNodalSources[_fieldName_] = DENS_MAT(s.nRows(),s.nCols());
    }
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++)
    {
      _fieldName_ = src_iter->first;
      if (!rhsMask((int)_fieldName_,ROBIN_SOURCE)) continue; 
      
      typedef map<PAIR,Array<UXT_Function*> > FSET;
      const FSET *fset = (const FSET *)&(src_iter->second);
      FSET::const_iterator fset_iter;
      for (fset_iter = fset->begin(); fset_iter != fset->end(); fset_iter++)
      {
        const PAIR &face = fset_iter->first;
        const int elem = face.first;
        // if this is not our element, do not do calculations
        if (!feMesh_->is_owned_elt(elem)) continue;
        const Array <UXT_Function*> &fs = fset_iter->second;
        _conn_ = feMesh_->element_connectivity_unique(elem);
        // evaluate location at ips
        feMesh_->face_shape_function(face, _fN_, _fdN_, _nN_, _fweights_);
        feMesh_->element_coordinates(elem, xCoords);
        xAtIPs = xCoords*(_fN_.transpose()); 
        // collect field
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        feMesh_->element_field(elem, field, localField);
        uAtIPs = _fN_*localField;

        // interpolate prescribed flux at ips of this element
        FSET::const_iterator face_iter = fset->find(face);
        int nFieldDOF = (face_iter->second).size();
        faceSource.reset(nIPsPerFace_,nFieldDOF);
        for (int ip = 0; ip < nIPsPerFace_; ++ip) {
          for (int idof = 0; idof<nFieldDOF; ++idof) {
            UXT_Function * f = fs(idof);
            if (!f) continue;

            faceSource(ip,idof) = f->f(&(uAtIPs(ip,0)),
                                       column(xAtIPs,ip).ptr(),time);

            DENS_MAT coefsAtIPs(nIPsPerFace_,1);
            coefsAtIPs(ip,idof) = f->dfdu(&(uAtIPs(ip,0)),
                                       column(xAtIPs,ip).ptr(),time);
            faceSource(ip,idof) -= coefsAtIPs(ip,0)*uAtIPs(ip,0);
          }
        }
        
        // assemble
        DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
        _Nmat_ = _fN_.transMat(_fweights_*faceSource);
        for (int i = 0; i < nNodesPerElement_; ++i) {
          int inode = _conn_(i);
          for (int idof = 0; idof < nFieldDOF; ++idof) {
            s(inode,idof) += _Nmat_(i,idof);
          }  
        } 
      }  
    } 
    // assemble partial result matrices
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++) {
      _fieldName_ = src_iter->first;
      if (!rhsMask((int) _fieldName_,ROBIN_SOURCE)) continue; 
      DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
      allsum(communicator_,MPI_IN_PLACE, s.ptr(), s.size());  
      DENS_MAT & src(nodalSources[_fieldName_].set_quantity());
      src += s;
    }
  }

  //-----------------------------------------------------------------
  // Robin boundary flux stiffness using native quadrature 
  // integrate \int_delV _N_I ds/du(x,t).n dA
  //-----------------------------------------------------------------
  void FE_Engine::add_robin_tangent(const Array2D<bool> &rhsMask,
    const FIELDS        & fields,
    const double time,
    const ROBIN_SURFACE_SOURCE & sourceFunctions,
    SPAR_MAT &tangent) const
  {
    
    // sizing working arrays
    DENS_MAT xCoords(nSD_,nNodesPerElement_);
    DENS_MAT coefsAtIPs;
    DENS_MAT localField;
    DENS_MAT xAtIPs(nSD_,nIPsPerFace_);
    DENS_MAT uAtIPs(nIPsPerFace_,1);

    // element wise assembly
    ROBIN_SURFACE_SOURCE::const_iterator src_iter;
    SPAR_MAT K(tangent.nRows(), tangent.nCols());
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++)
    {
      _fieldName_ = src_iter->first;
      if (!rhsMask((int)_fieldName_,ROBIN_SOURCE)) continue; 
      
      typedef map<PAIR,Array<UXT_Function*> > FSET;
      const FSET *fset = (const FSET *)&(src_iter->second);
      FSET::const_iterator fset_iter;
      for (fset_iter = fset->begin(); fset_iter != fset->end(); fset_iter++)
      {
        const PAIR &face = fset_iter->first;
        const int elem = face.first;
        // if this is not our element, do not do calculations
        if (!feMesh_->is_owned_elt(elem)) continue;
        const Array <UXT_Function*> &fs = fset_iter->second;
        _conn_ = feMesh_->element_connectivity_unique(elem);
        // evaluate location at ips
        feMesh_->face_shape_function(face, _fN_, _fdN_, _nN_, _fweights_);
        feMesh_->element_coordinates(elem, xCoords);
        xAtIPs = xCoords*(_fN_.transpose()); 
        // collect field
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        feMesh_->element_field(elem, field, localField);
        uAtIPs = _fN_*localField;

        // interpolate prescribed flux at ips of this element
        FSET::const_iterator face_iter = fset->find(face);
        int nFieldDOF = (face_iter->second).size();
        coefsAtIPs.reset(nIPsPerFace_,nFieldDOF);
        for (int ip = 0; ip < nIPsPerFace_; ++ip) {
          for (int idof = 0; idof<nFieldDOF; ++idof) {
            UXT_Function * f = fs(idof);
            if (!f) continue;
            coefsAtIPs(ip,idof) = f->dfdu(&(uAtIPs(ip,0)),
                                       column(xAtIPs,ip).ptr(),time);
          }
        }
        
        // assemble
        DIAG_MAT D(column(coefsAtIPs,0)); 
        D *= -1; 
        D *= _fweights_; 
        _Nmat_ = _fN_.transMat(D*_fN_);
        for (int i = 0; i < nNodesPerElement_; ++i) {
          int inode = _conn_(i);
          for (int j = 0; j < nNodesPerElement_; ++j) {
            int jnode = _conn_(j);
            K.add(inode, jnode, _Nmat_(i,j));
          }  
        } 
      }  
    } 
    // assemble partial result matrices
#ifdef ISOLATE_FE
    sparse_allsum(communicator_,K);
#else
    LammpsInterface::instance()->sparse_allsum(K);
#endif
    tangent += K;
  }

  //-----------------------------------------------------------------
  // Robin boundary flux using native quadrature 
  // integrate \int_delV _N_I s(x,t).n dA
  //-----------------------------------------------------------------
  void FE_Engine::add_open_fluxes(const Array2D<bool> &rhsMask,
                                  const FIELDS        & fields,
                                  const OPEN_SURFACE  & openFaces,
                                  FIELDS &nodalSources,
                                  const FieldName Velocity) const
  {
    
    // sizing working arrays
    DENS_MAT xCoords(nSD_,nNodesPerElement_);
    DENS_MAT faceSource;
    DENS_MAT localField;
    DENS_MAT xAtIPs(nSD_,nIPsPerFace_);
    DENS_MAT uAtIPs; // (nIPsPerFace_,field nCols);
    DENS_MAT aAtIPs(nIPsPerFace_,1);
    FIELDS myNodalSources;

    // throw error if electron velocity is not defined

    // element wise assembly
    OPEN_SURFACE::const_iterator face_iter;
    for (face_iter=openFaces.begin(); face_iter!=openFaces.end(); face_iter++)
    {
      _fieldName_ = face_iter->first;
      if (!rhsMask((int)_fieldName_,OPEN_SOURCE)) continue; 
      DENS_MAT & s(nodalSources[_fieldName_].set_quantity());
      myNodalSources[_fieldName_] = DENS_MAT(s.nRows(),s.nCols());
    }
    for (face_iter=openFaces.begin(); face_iter!=openFaces.end(); face_iter++)
    {
      _fieldName_ = face_iter->first;
      if (!rhsMask((int)_fieldName_,OPEN_SOURCE)) continue; 
      
      typedef set<PAIR> FSET;
      const FSET *fset = (const FSET *)&(face_iter->second);
      FSET::const_iterator fset_iter;
      for (fset_iter = fset->begin(); fset_iter != fset->end(); fset_iter++)
      {
        const PAIR &face = *fset_iter;
        const int elem = face.first;
        // if this is not our element, do not do calculations
        if (!feMesh_->is_owned_elt(elem)) continue;
        // get velocity, multiply with field (vector gives tensor)
        // dot with face normal
        _conn_ = feMesh_->element_connectivity_unique(elem);
        // evaluate location at ips
        feMesh_->face_shape_function(face, _fN_, _fdN_, _nN_, _fweights_);
        // collect field at ips
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        feMesh_->element_field(elem, field, localField);
        int nFieldDOF = field.nCols();
        // "u" is the quantity being advected
        uAtIPs.reset(nIPsPerFace_,nFieldDOF);
        uAtIPs = _fN_*localField;
        // collect velocity at ips for "advection" = v.n
        _fieldItr_ = fields.find(Velocity);
        const DENS_MAT & advection = (_fieldItr_->second).quantity(); // advection velocity
        feMesh_->element_field(elem, advection, localField);
        for (int j = 0; j < nSD_; ++j) { // nSD_ == nDOF for the velocity
          for (int ip = 0; ip < nIPsPerFace_; ++ip) {
            for (int I = 0; I < nNodesPerElement_; ++I) {
              aAtIPs(ip,0) += (_nN_[j])(ip,I)*localField(I,j);
            }
          }
        }
        // construct open flux at ips of this element
        // flux = field \otimes advection_vector \dot n
        faceSource.reset(nIPsPerFace_,nFieldDOF);
        for (int ip = 0; ip < nIPsPerFace_; ++ip) {
          for (int idof = 0; idof<nFieldDOF; ++idof) {
            // multiply field DOF value times advection vector
            faceSource(ip,idof) = aAtIPs(ip,1)*uAtIPs(ip,idof);//(v.n) u
          }
        }
        // assemble contributions for this direction of the face normal
        // _nN_[j](ip,I) nodal shapefunction(I) at ip X face normal(j)
        // _fweights_ diagonal nips X nips ==> facesource is nips X ndofs
        DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
        _Nmat_ = _fN_.transMat(_fweights_*faceSource);
        // s_Ii = \sum_ip N_I (u_i v.n)_ip wg_ip
        for (int i = 0; i < nNodesPerElement_; ++i) {
          int inode = _conn_(i);
          for (int idof = 0; idof < nFieldDOF; ++idof) {
            s(inode,idof) += _Nmat_(i,idof);
          }  
        }
      }  
    } 
    // assemble partial result matrices
    for (face_iter=openFaces.begin(); face_iter!=openFaces.end(); face_iter++) {
      _fieldName_ = face_iter->first;
      if (!rhsMask((int) _fieldName_,OPEN_SOURCE)) continue; 
      DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
      allsum(communicator_,MPI_IN_PLACE, s.ptr(), s.size());  
      DENS_MAT & src(nodalSources[_fieldName_].set_quantity());
      src += s;
    }
  }

  //-----------------------------------------------------------------
  // Open boundary flux stiffness using native quadrature 
  // integrate \int_delV _N_I ds/du(x,t).n dA
  //-----------------------------------------------------------------
  void FE_Engine::add_open_tangent(const Array2D<bool> &rhsMask,
    const FIELDS        & fields,
    const OPEN_SURFACE  & openFaces,
    SPAR_MAT &tangent,
    const FieldName Velocity) const
  {
    
    // sizing working arrays
    DENS_MAT xCoords(nSD_,nNodesPerElement_);
    DENS_MAT faceSource;
    DENS_MAT localField;
    DENS_MAT xAtIPs(nSD_,nIPsPerFace_);
    DENS_MAT uAtIPs; // (nIPsPerFace_,field nCols);
    DENS_MAT aAtIPs(nIPsPerFace_,nSD_);
    SPAR_MAT K(tangent.nRows(), tangent.nCols());
    // element wise assembly
    OPEN_SURFACE::const_iterator face_iter;
    for (face_iter=openFaces.begin(); face_iter!=openFaces.end(); face_iter++)
    {
      _fieldName_ = face_iter->first;
      if (!rhsMask((int)_fieldName_,OPEN_SOURCE)) continue; 
      bool convective = false;
      if (_fieldName_ == Velocity) convective = true;
      
      typedef set<PAIR> FSET;
      const FSET *fset = (const FSET *)&(face_iter->second);
      FSET::const_iterator fset_iter;
      for (fset_iter = fset->begin(); fset_iter != fset->end(); fset_iter++)
      {
        const PAIR &face = *fset_iter;
        const int elem = face.first;
        // if this is not our element, do not do calculations
        if (!feMesh_->is_owned_elt(elem)) continue;
        _conn_ = feMesh_->element_connectivity_unique(elem);
        // evaluate location at ips
        feMesh_->face_shape_function(face, _fN_, _fdN_, _nN_, _fweights_);
        // collect field at ips
        _fieldItr_ = fields.find(_fieldName_);
        const DENS_MAT & field = (_fieldItr_->second).quantity();
        feMesh_->element_field(elem, field, localField);
        int nFieldDOF = field.nCols();
        // "u" is the quantity being advected
        uAtIPs.reset(nIPsPerFace_,nFieldDOF);
        uAtIPs = _fN_*localField;
        // collect velocity at ips for "advection" = v.n
        _fieldItr_ = fields.find(Velocity);
        const DENS_MAT & advection = (_fieldItr_->second).quantity(); // advection velocity
        feMesh_->element_field(elem, advection, localField);
        for (int j = 0; j < nSD_; ++j) { // nSD_ == nDOF for the velocity
          for (int ip = 0; ip < nIPsPerFace_; ++ip) {
            for (int I = 0; I < nNodesPerElement_; ++I) {
              aAtIPs(ip,0) += (_nN_[j])(ip,I)*localField(I,j);
            }
          }
        }
        // K_IJ = \sum_ip N_I ( v.n I + u (x) n ) N_J wg_ip
        DENS_MAT D(nFieldDOF,nFieldDOF);
        DENS_VEC n(nSD_);
        for (int ip = 0; ip < nIPsPerFace_; ++ip) {
          for (int idof = 0; idof<nFieldDOF; ++idof) {
            D(idof,idof) -= aAtIPs(ip,0);
            if (convective) {
              feMesh_->face_normal(face,ip,n);
              for (int jdof = 0; jdof<nFieldDOF; ++jdof) {
                D(idof,jdof) -=  uAtIPs(ip,idof)*n(jdof);
              }
            }
          }
        }
        // assemble
        _Nmat_ = _fN_.transMat(D*_fN_);
        for (int i = 0; i < nNodesPerElement_; ++i) {
          int inode = _conn_(i);
          for (int j = 0; j < nNodesPerElement_; ++j) {
            int jnode = _conn_(j);
            K.add(inode, jnode, _Nmat_(i,j));
          }  
        } 
      }  
    } 
    // assemble partial result matrices
#ifdef ISOLATE_FE
    sparse_allsum(communicator_,K);
#else
    LammpsInterface::instance()->sparse_allsum(K);
#endif
    tangent += K;
  }
  
  //-----------------------------------------------------------------
  // prescribed boundary flux using native quadrature 
  // integrate \int_delV _N_I s(x,t).n dA
  //-----------------------------------------------------------------
  void FE_Engine::add_fluxes(const Array<bool> &fieldMask,
    const double time,
    const SURFACE_SOURCE & sourceFunctions,
    FIELDS &nodalSources) const
  {
    
    // sizing working arrays
    DENS_MAT xCoords(nSD_,nNodesPerElement_);
    DENS_MAT xAtIPs(nSD_,nIPsPerFace_);
    DENS_MAT faceSource;
    FIELDS myNodalSources;

    // element wise assembly
    SURFACE_SOURCE::const_iterator src_iter;
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++)
    {
      _fieldName_ = src_iter->first;
      if (!fieldMask((int)_fieldName_)) continue; 
      DENS_MAT & s(nodalSources[_fieldName_].set_quantity());
      myNodalSources[_fieldName_] = DENS_MAT(s.nRows(),s.nCols());
    }
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++)
    {
      _fieldName_ = src_iter->first;
      if (!fieldMask((int)_fieldName_)) continue; 
      
      typedef map<PAIR,Array<XT_Function*> > FSET;
      const FSET *fset = (const FSET *)&(src_iter->second);
      FSET::const_iterator fset_iter;
      for (fset_iter = fset->begin(); fset_iter != fset->end(); fset_iter++)
      {
        const PAIR &face = fset_iter->first;
        const int elem = face.first;
        // if this is not our element, do not do calculations
        if (!feMesh_->is_owned_elt(elem)) continue;
        const Array <XT_Function*> &fs = fset_iter->second;
        _conn_ = feMesh_->element_connectivity_unique(elem);
        // evaluate location at ips
        feMesh_->face_shape_function(face, _fN_, _fdN_, _nN_, _fweights_);
        feMesh_->element_coordinates(elem, xCoords);

        MultAB(xCoords,_fN_,xAtIPs,0,1); //xAtIPs = xCoords*(N.transpose());

        // interpolate prescribed flux at ips of this element
        
        FSET::const_iterator face_iter = fset->find(face);
        if (face_iter == fset->end()) 
        {  
          stringstream ss;
          ss << "face not found" << std::endl;
          print_msg(communicator_,ss.str());
          
        }
        
        
        int nFieldDOF = (face_iter->second).size();
        faceSource.reset(nIPsPerFace_,nFieldDOF);
        for (int ip = 0; ip < nIPsPerFace_; ++ip) {
          for (int idof = 0; idof<nFieldDOF; ++idof) {
            XT_Function * f = fs(idof);
            if (!f) continue;
            faceSource(ip,idof) = f->f(column(xAtIPs,ip).ptr(),time);
          }
        }
        
        // assemble
        DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
        _Nmat_ = _fN_.transMat(_fweights_*faceSource);
        for (int i = 0; i < nNodesPerElement_; ++i) 
        {
          int inode = _conn_(i);
          for (int idof = 0; idof < nFieldDOF; ++idof) 
          {
            s(inode,idof) += _Nmat_(i,idof);
          }  // end assemble nFieldDOF
        }  // end assemble nNodesPerElement_
      }  // end fset loop
    } // field loop
    // assemble partial result matrices
    for (src_iter=sourceFunctions.begin(); src_iter!=sourceFunctions.end(); src_iter++) {
      _fieldName_ = src_iter->first;
      if (!fieldMask((int)_fieldName_)) continue; 
      DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
      allsum(communicator_,MPI_IN_PLACE, s.ptr(), s.size());  
      DENS_MAT & src(nodalSources[_fieldName_].set_quantity());
      src += s;
    }
  }

  //-----------------------------------------------------------------
  // prescribed volume flux using native quadrature 
  // integrate \int_V _N_I s(x,t) dV
  //-----------------------------------------------------------------
  void FE_Engine::add_sources(const Array<bool> &fieldMask,
    const double time,
    const VOLUME_SOURCE &sourceFunctions,
    FIELDS &nodalSources) const
  {
    
    // sizing working arrays
    DENS_MAT elemSource;
    DENS_MAT xCoords(nSD_,nNodesPerElement_);
    DENS_MAT xAtIPs(nSD_,nIPsPerElement_);
    FIELDS myNodalSources;

    for (VOLUME_SOURCE::const_iterator src_iter = sourceFunctions.begin(); 
         src_iter != sourceFunctions.end(); src_iter++) {
      _fieldName_ = src_iter->first;
      int index = (int) _fieldName_;
      if (! fieldMask(index)) continue;
      DENS_MAT & s(nodalSources[_fieldName_].set_quantity());
      myNodalSources[_fieldName_] = DENS_MAT(s.nRows(),s.nCols());
    }
    vector<int> myElems = feMesh_->owned_elts();
    for (vector<int>::iterator elemsIter = myElems.begin();
         elemsIter != myElems.end();
         ++elemsIter)
    {
      int ielem = *elemsIter;
      _conn_ = feMesh_->element_connectivity_unique(ielem);
      // evaluate location at ips
      feMesh_->shape_function(ielem, _N_, _weights_);
      feMesh_->element_coordinates(ielem, xCoords);
      xAtIPs =xCoords*(_N_.transpose());
      for (VOLUME_SOURCE::const_iterator src_iter = sourceFunctions.begin(); 
           src_iter != sourceFunctions.end(); src_iter++) {
        _fieldName_ = src_iter->first;
        int index = (int) _fieldName_; 
        if ( fieldMask(index) ) {
          const Array2D<XT_Function *> * thisSource 
            = (const Array2D<XT_Function *> *) &(src_iter->second);
          int nFieldDOF = thisSource->nCols();
          elemSource.reset(nIPsPerElement_,nFieldDOF);
          // interpolate prescribed flux at ips of this element
          for (int ip = 0; ip < nIPsPerElement_; ++ip) {
            for (int idof = 0; idof < nFieldDOF; ++idof) {
              XT_Function * f = (*thisSource)(ielem,idof);
              if (f) {
                elemSource(ip,idof) = f->f(column(xAtIPs,ip).ptr(),time);
              }
            }
          }
          // assemble
          _Nmat_ = _N_.transMat(_weights_*elemSource);
          DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
          
          for (int i = 0; i < nNodesPerElement_; ++i) {
            int inode = _conn_(i);
            for (int idof = 0; idof < nFieldDOF; ++idof) {
              s(inode,idof) += _Nmat_(i,idof);
            }
          }
        }
      }
    }
    // Aggregate unmasked nodal sources on all processors.
    for (VOLUME_SOURCE::const_iterator src_iter = sourceFunctions.begin(); 
         src_iter != sourceFunctions.end(); src_iter++) {
      _fieldName_ = src_iter->first;
      int index = (int) _fieldName_;
      if (!fieldMask(index)) continue;
      DENS_MAT & s(myNodalSources[_fieldName_].set_quantity());
      allsum(communicator_,MPI_IN_PLACE, s.ptr(), s.size());  
      DENS_MAT & src(nodalSources[_fieldName_].set_quantity());
      src += s;
    }
  }
  //-----------------------------------------------------------------
  // previously computed nodal sources
  //-----------------------------------------------------------------
  void FE_Engine::add_sources(const Array<bool> &fieldMask,
    const double time,
    const FIELDS &sources,
    FIELDS &nodalSources) const
  {
    FIELDS::const_iterator itr;
    for (itr=sources.begin(); itr!=sources.end(); itr++) {
      _fieldName_ = itr->first;
      if (!fieldMask((int)_fieldName_)) continue; 
      DENS_MAT & src(nodalSources[_fieldName_].set_quantity());
      const DENS_MAT &  s((sources.find(_fieldName_)->second).quantity());
      for (int inode = 0; inode < nNodesUnique_; ++inode) {
        for (int j = 0; j < src.nCols(); ++j) {
          src(inode,j) += s(inode,j);
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
    int dof = field.nCols();

    double a[3] = {0,0,0};
    a[axis] = 1;

    // sizing working arrays
    DENS_MAT n(nSD_,nIPsPerFace_);
    DENS_MAT localElementFields(nNodesPerElement_,dof);
    DENS_MAT integrals(dof,nSD_); 
    DENS_MAT fieldsAtIPs; 
                          // SJL shouldn't this just be _fieldsAtIPs_
                          // the data member?

    // element wise assembly
    set< pair <int,int> >::iterator iter;
    for (iter = faceSet.begin(); iter != faceSet.end(); iter++) {
      int ielem = iter->first;
      // if this is not our element, do not do calculations
      //if (!feMesh_->is_owned_elt(ielem)) continue;

      // evaluate shape functions at ips
      feMesh_->face_shape_function(*iter, _N_, n, _fweights_);
      // cross n with axis to get tangent
      if (contour) {
        double t[3];
        for (int i = 0; i < nIPsPerFace_; i++) {
          t[0] = a[1]*n(2,i) - a[2]*n(1,i);
          t[1] = a[2]*n(0,i) - a[0]*n(2,i);
          t[2] = a[0]*n(1,i) - a[1]*n(0,i);
          n(0,i) = t[0];
          n(1,i) = t[1];
          n(2,i) = t[2];
        }
      }

      // get connectivity
      _conn_ = feMesh_->element_connectivity_unique(ielem);

      // interpolate fields and gradients of fields ips of this element
      for (int i = 0; i < nNodesPerElement_; i++) {
        for (int j = 0; j < dof; j++) {
          localElementFields(i,j) = field(_conn_(i),j);
        }
      }
      //  ips X dof    = ips X nodes * nodes X dof
      fieldsAtIPs = _N_*localElementFields;

      // assemble : integral(k,j) = sum_ip n(j,ip) wg(ip,ip) field(ip,k)
      _Nmat_ = n*_fweights_*fieldsAtIPs;
      for (int j = 0; j < nSD_; j++) {
        for (int k = 0; k < dof; ++k) {
          integrals(k,j) += _Nmat_(j,k);
        }
      }
    }
    // assemble partial results
    //LammpsInterface::instance()->allsum(MPI_IN_PLACE, integrals.ptr(), integrals.size());
    // (S.n)_1 = S_1j n_j = S_11 n_1 + S_12 n_2 + S_13 n_3
    // (S.n)_2 = S_2j n_j = S_21 n_1 + S_22 n_2 + S_23 n_3
    // (S.n)_3 = S_3j n_j = S_31 n_1 + S_32 n_2 + S_33 n_3
    if (dof == 9) { // tensor
      values.reset(nSD_,1);
      values(0,0) = integrals(0,0)+integrals(1,1)+integrals(2,2);
      values(1,0) = integrals(3,0)+integrals(4,1)+integrals(5,2);
      values(2,0) = integrals(6,0)+integrals(7,1)+integrals(8,2);
    }
    else if (dof == 6) { // sym tensor
      values.reset(nSD_,1);
      values(0,0) = integrals(0,0)+integrals(3,1)+integrals(4,2);
      values(1,0) = integrals(3,0)+integrals(1,1)+integrals(5,2);
      values(2,0) = integrals(4,1)+integrals(5,1)+integrals(2,2);
    }
    // (v.n) = v_j n_j 
    else if (dof == 3) { // vector
      values.reset(1,1);
      values(0,0) = integrals(0,0)+integrals(1,1)+integrals(2,2);
    }
    // s n = s n_j e_j
    else if (dof == 1) { // scalar
      values.reset(nSD_,1);
      values(0,0) = integrals(0,0);
      values(1,0) = integrals(0,1);
      values(2,0) = integrals(0,2);
    }
    else {
      string msg = "FE_Engine::field_surface_flux unsupported field width: ";
      msg += to_string(dof);
      throw ATC_Error(msg);
    }
  }

  //-----------------------------------------------------------------
  // evaluate shape functions at given points
  //-----------------------------------------------------------------
  
  void FE_Engine::evaluate_shape_functions(
    const MATRIX & pt_coords,
    SPAR_MAT & N) const
  {
    // Get shape function and derivatives at atomic locations
    int nnodes = feMesh_->num_nodes_unique();
    int npts = pt_coords.nRows();
    
    // loop over point (atom) coordinates
    DENS_VEC x(nSD_);
    Array<int> node_index(nNodesPerElement_);
    DENS_VEC shp(nNodesPerElement_);
    
    N.reset(npts,nnodes);
    for (int i = 0; i < npts; ++i) {
      for (int k = 0; k < nSD_; ++k) { x(k) = pt_coords(i,k); }
      feMesh_->shape_functions(x,shp,node_index);
      for (int j = 0; j < nNodesPerElement_; ++j) {
        int jnode = node_index(j);
        N.add(i,jnode,shp(j));
      }
    }
    N.compress();   
  }

  //-----------------------------------------------------------------
  // evaluate shape functions and their spatial derivatives at given points
  //-----------------------------------------------------------------
  void FE_Engine::evaluate_shape_functions(
    const MATRIX & pt_coords,
    SPAR_MAT & N,
    SPAR_MAT_VEC & dN) const
  {
    // Get shape function and derivatives at atomic locations
    int nnodes = feMesh_->num_nodes_unique();
    int npts = pt_coords.nRows();

    // loop over point (atom) coordinates
    DENS_VEC x(nSD_);
    Array<int> node_index(nNodesPerElement_);
    DENS_VEC shp(nNodesPerElement_);
    DENS_MAT dshp(nSD_,nNodesPerElement_);
    
    for (int k = 0; k < nSD_; ++k) {
      dN[k]->reset(npts,nnodes);    
    }
          
    N.reset(npts,nnodes); 
    for (int i = 0; i < npts; ++i) {
      for (int k = 0; k < nSD_; ++k) { x(k) = pt_coords(i,k); }
      feMesh_->shape_functions(x,shp,dshp,node_index);
      for (int j = 0; j < nNodesPerElement_; ++j) {
        int jnode = node_index(j);
        N.add(i,jnode,shp(j));
        for (int k = 0; k < nSD_; ++k) {
          dN[k]->add(i,jnode,dshp(k,j));
        }
      }
    }
    N.compress();   
    for (int k = 0; k < nSD_; ++k) { 
      dN[k]->compress();   
    }
  }

  //-----------------------------------------------------------------
  // evaluate shape functions at given points
  //-----------------------------------------------------------------
    void FE_Engine::evaluate_shape_functions(
    const MATRIX & pt_coords,
    const INT_ARRAY & pointToEltMap,
    SPAR_MAT & N) const
  {
    // Get shape function and derivatives at atomic locations
    int nnodes = feMesh_->num_nodes_unique();
    int npts = pt_coords.nRows();
    
    // loop over point (atom) coordinates
    DENS_VEC x(nSD_);
    Array<int> node_index(nNodesPerElement_);
    DENS_VEC shp(nNodesPerElement_);
    
    N.reset(npts,nnodes);
    for (int i = 0; i < npts; ++i) {
      for (int k = 0; k < nSD_; ++k) { x(k) = pt_coords(i,k); }
      int eltId = pointToEltMap(i,0);
      feMesh_->shape_functions(x,eltId,shp,node_index);
      for (int j = 0; j < nNodesPerElement_; ++j) {
        int jnode = node_index(j);
        N.add(i,jnode,shp(j));
      }
    }
    N.compress();   
  }

  //-----------------------------------------------------------------
  // evaluate shape functions and their spatial derivatives at given points
  //-----------------------------------------------------------------
  void FE_Engine::evaluate_shape_functions(
    const MATRIX & pt_coords,
    const INT_ARRAY & pointToEltMap,
    SPAR_MAT & N,
    SPAR_MAT_VEC & dN) const
  {
    // Get shape function and derivatives at atomic locations
    int nnodes = feMesh_->num_nodes_unique();
    int npts = pt_coords.nRows();

    // loop over point (atom) coordinates
    DENS_VEC x(nSD_);
    Array<int> node_index(nNodesPerElement_);
    DENS_VEC shp(nNodesPerElement_);
    DENS_MAT dshp(nSD_,nNodesPerElement_);
    for (int k = 0; k < nSD_; ++k) {
      dN[k]->reset(npts,nnodes);    
    }
          
    N.reset(npts,nnodes); 
    for (int i = 0; i < npts; ++i) {
      for (int k = 0; k < nSD_; ++k) { x(k) = pt_coords(i,k); }
      int eltId = pointToEltMap(i,0);
      feMesh_->shape_functions(x,eltId,shp,dshp,node_index);
      for (int j = 0; j < nNodesPerElement_; ++j) {
        int jnode = node_index(j);
        N.add(i,jnode,shp(j));
        for (int k = 0; k < nSD_; ++k) {
          dN[k]->add(i,jnode,dshp(k,j));
        }
      }
    }
    N.compress();   
    for (int k = 0; k < nSD_; ++k) { 
      dN[k]->compress();   
    }
  }
  //-----------------------------------------------------------------
  // evaluate shape functions spatial derivatives at given points
  //-----------------------------------------------------------------
  void FE_Engine::evaluate_shape_function_derivatives(
    const MATRIX & pt_coords,
    const INT_ARRAY & pointToEltMap,
    SPAR_MAT_VEC & dNdx ) const
  {
    int nnodes = feMesh_->num_nodes_unique();
    int npts = pt_coords.nRows();
    for (int k = 0; k < nSD_; ++k) {
      dNdx[k]->reset(npts,nnodes);    
    }

    // loop over point (atom) coordinates
    DENS_VEC x(nSD_);
    Array<int> node_index(nNodesPerElement_);
    DENS_MAT dshp(nSD_,nNodesPerElement_);
    
    for (int i = 0; i < npts; ++i) {
      for (int k = 0; k < nSD_; ++k) { x(k) = pt_coords(i,k); }
      int eltId = pointToEltMap(i,0);
      feMesh_->shape_function_derivatives(x,eltId,dshp,node_index);
      for (int j = 0; j < nNodesPerElement_; ++j) {
        int jnode = node_index(j);
        for (int k = 0; k < nSD_; ++k) {
          dNdx[k]->add(i,jnode,dshp(k,j));
        }
      }
    }
    for (int k = 0; k < nSD_; ++k) { 
      dNdx[k]->compress();   
    }
  }

  //-------------------------------------------------------------------
  // set kernels
  //-------------------------------------------------------------------
  void FE_Engine::set_kernel(KernelFunction* ptr){kernelFunction_=ptr;}

  //-------------------------------------------------------------------
  // kernel_matrix_bandwidth
  //-------------------------------------------------------------------
  int FE_Engine::kernel_matrix_bandwidth(
    const MATRIX & ptCoords) const
  {

    if (! kernelFunction_) throw ATC_Error("no kernel function specified");
    int npts = ptCoords.nRows();
    int bw = 0;
    if (npts>0) {
      DENS_VEC xI(nSD_),x(nSD_),dx(nSD_);
      for (int i = 0; i < nNodes_; i++) {
        xI = feMesh_->global_coordinates(i);
        for (int j = 0; j < npts; j++) {
          for (int k = 0; k < nSD_; ++k) { x(k) = ptCoords(j,k); } 
          dx = x - xI;
          if (kernelFunction_->in_support(dx)) bw = max(bw,abs(j-map_global_to_unique(i)));
        }
      }
    }
    return bw;
  }
  //-------------------------------------------------------------------
  // evaluate kernels
  //-------------------------------------------------------------------
  void FE_Engine::evaluate_kernel_functions(
    const MATRIX & ptCoords,
    SPAR_MAT & N ) const
  {
#ifdef EXTENDED_ERROR_CHECKING

   print_msg_once(communicator_,"kernel matrix bandwidth " +to_string(kernel_matrix_bandwidth(ptCoords)));
#endif

    if (! kernelFunction_) throw ATC_Error("no kernel function specified");
    int npts = ptCoords.nRows();
    if (npts>0) {
      N.reset(npts,nNodesUnique_); // transpose of N_Ia
      DENS_VEC xI(nSD_),x(nSD_),dx(nSD_);
      for (int i = 0; i < nNodes_; i++) {
        xI = feMesh_->global_coordinates(i);
        for (int j = 0; j < npts; j++) {
          for (int k = 0; k < nSD_; ++k) { x(k) = ptCoords(j,k); } 
          dx = x - xI;
          double val = kernelFunction_->value(dx);
          val *= kernelFunction_->dimensionality_factor();
          if (val > 0) N.add(j,map_global_to_unique(i),val);
        }
      }
      N.compress();
    }
  }

  //-----------------------------------------------------------------
  // create a non-uniform rectangular mesh on a rectangular region
  //-----------------------------------------------------------------

  void FE_Engine::create_mesh(Array<double> & dx, 
    Array<double> & dy, 
    Array<double> & dz, 
    const char * regionName,
    Array<bool> periodicity)
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xscale, yscale, zscale;

    // check to see if region exists and is it a box, and if so get the bounds
    bool isBox;
    isBox = ATC::LammpsInterface::instance()->region_bounds(
                                    regionName,
                                    xmin, xmax,
                                    ymin, ymax,
                                    zmin, zmax,
                                    xscale,
                                    yscale,
                                    zscale);
    if (!isBox) throw ATC_Error("Region for FE mesh is not a box");


    if (dx(0) == 0.0) dx = (xmax-xmin)/dx.size();
    if (dy(0) == 0.0) dy = (ymax-ymin)/dy.size();
    if (dz(0) == 0.0) dz = (zmax-zmin)/dz.size();

    
    feMesh_ = new FE_Rectangular3DMesh(dx,dy,dz,
                                       xmin, xmax,
                                       ymin, ymax,
                                       zmin, zmax,
                                       periodicity,
                                       xscale, yscale, zscale);
    stringstream ss;
    ss << "created structured mesh with " << feMesh_->num_nodes() << " nodes, " 
       << feMesh_->num_nodes_unique() << " unique nodes, and " 
       << feMesh_->num_elements() << " elements";
    print_msg_once(communicator_,ss.str());
#ifdef ATC_VERBOSE
    print_partitions(xmin,xmax,dx);
    print_partitions(ymin,ymax,dy);
    print_partitions(zmin,zmax,dz);
#endif
  }
  //-----------------------------------------------------------------
  void FE_Engine::print_partitions(double xmin, double xmax, Array<double> & dx) const {
    stringstream msg;
    msg.precision(3);
    msg << std::fixed;
    msg <<  "\nindex weight fraction location size[A] size[uc]:\n";
    double sum = 0;
    for (int i = 0; i < dx.size(); ++i) { sum += dx(i); }
    double xn= 1.0/sum;
    double xs= xn*(xmax-xmin);
    double xl= xs/(ATC::LammpsInterface::instance()->xlattice());
    double sumn = 0, sums = 0, suml = 0;
    double x = xmin;
    for (int i = 0; i < dx.size(); ++i) { 
       double dxn = dx(i)*xn; sumn += dxn;
       double dxs = dx(i)*xs; sums += dxs;
       double dxl = dx(i)*xl; suml += dxl;
       msg << std::setw(5) << i+1
           << std::setw(7) << dx(i)
           << std::setw(9) << dxn
           << std::setw(9) << x
           << std::setw(8) << dxs
           << std::setw(9) << dxl << "\n";
       x += dxs;
     }
     msg << "sum  " << setw(7) << sum
                    << setw(9) << sumn
                    << setw(9) << x
                    << setw(8) << sums
                    << setw(9) << suml << "\n";
     print_msg_once(communicator_,msg.str());
  }
  //-----------------------------------------------------------------
  // create a uniform rectangular mesh on a rectangular region
  //-----------------------------------------------------------------
  void FE_Engine::create_mesh(int nx, int ny, int nz, const char * regionName,
                              Array<bool> periodicity)
  {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xscale, yscale, zscale;

    // check to see if region exists and is it a box, and if so get the bounds
    bool isBox;
    isBox = ATC::LammpsInterface::instance()->region_bounds(
                                    regionName,
                                    xmin, xmax,
                                    ymin, ymax,
                                    zmin, zmax,
                                    xscale, yscale, zscale);
    if (!isBox) throw ATC_Error("Region for FE mesh is not a box");
    
    feMesh_ = new FE_Uniform3DMesh(nx,ny,nz,
                                   xmin, xmax,
                                   ymin, ymax,
                                   zmin, zmax,
                                   periodicity,
                                   xscale, yscale, zscale);
    stringstream ss;
    ss << "created uniform mesh with " << feMesh_->num_nodes() << " nodes, " 
       << feMesh_->num_nodes_unique() << " unique nodes, and " 
       << feMesh_->num_elements() << " elements";
    print_msg_once(communicator_,ss.str());
  }

  //-----------------------------------------------------------------
  // read a mesh from an input file using MeshReader object
  //-----------------------------------------------------------------
  void FE_Engine::read_mesh(string meshFile, Array<bool> & periodicity)
  {
    if (feMesh_) throw ATC_Error("FE_Engine already has a mesh");
    feMesh_ = MeshReader(meshFile,periodicity).create_mesh();
    feMesh_->initialize();
  }
};

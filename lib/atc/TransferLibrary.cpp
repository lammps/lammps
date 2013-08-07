// ATC headers
#include "TransferLibrary.h"
#include "ATC_Coupling.h"
#include "PrescribedDataManager.h"
#include "LinearSolver.h"
#include "PerAtomQuantityLibrary.h"
#include "KernelFunction.h"
#include "MoleculeSet.h"

//#include <typeinfo>

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class NodalAtomVolume
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  NodalAtomVolume::NodalAtomVolume(ATC_Method * atc,
                                   SPAR_MAN * shapeFunction) :
    atc_(atc),
    shapeFunction_(shapeFunction),
    lammpsInterface_(atc->lammps_interface()),
    feEngine_(atc->fe_engine()),
    tol_(1.e-10) 
  {
    shapeFunction_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void NodalAtomVolume::reset_quantity() const
  {
    // solve equation \sum_a N_Ia \sum_J N_Ja dV_J = \int_Omega N_I dV
    // form left-hand side
    int nNodes = shapeFunction_->nCols();
    SPAR_MAT lhs(nNodes,nNodes);
    atc_->compute_consistent_md_mass_matrix(shapeFunction_->quantity(),lhs);

    // form right-hand side
    _rhsMatrix_.resize(nNodes,nNodes);
    feEngine_->compute_lumped_mass_matrix(_rhsMatrix_);
    _rhs_.resize(nNodes);
    _rhs_.copy(_rhsMatrix_.ptr(),_rhsMatrix_.size(),1);

    // change entries for all entries if no atoms in shape function support
    double totalVolume = _rhs_.sum();
    double averageVolume = averaging_operation(totalVolume);

    _scale_.resize(nNodes);
    for (int i = 0; i < nNodes; i++) {
      if ((abs(lhs(i,i)) > 0.))  
        _scale_(i) = 1.;
      else
        _scale_(i) = 0.;
    }
    lhs.row_scale(_scale_);
    for (int i = 0; i < nNodes; i++) {
      if (_scale_(i) < 0.5) {
        lhs.set(i,i,1.);
        _rhs_(i) = averageVolume;
      }
    }
    lhs.compress();
    
    // solve equation
    LinearSolver solver(lhs, ATC::LinearSolver::ITERATIVE_SOLVE_SYMMETRIC, true);
    solver.set_max_iterations(lhs.nRows());
    solver.set_tolerance(tol_);
    quantity_.reset(nNodes,1);
    CLON_VEC tempQuantity(quantity_,CLONE_COL,0);
    solver.solve(tempQuantity,_rhs_);
  }

  //--------------------------------------------------------
  //  averaging_operation
  //--------------------------------------------------------
  double NodalAtomVolume::averaging_operation(const double totalVolume) const
  {
    int nLocal[1] = {shapeFunction_->nRows()};
    int nGlobal[1] = {0};
    lammpsInterface_->int_allsum(nLocal,nGlobal,1);
    return totalVolume/(double(nGlobal[0]));
  }

  //--------------------------------------------------------
  //  overloading operation to get number of rows
  //--------------------------------------------------------
  int NodalAtomVolume::nRows() const
  {
    return atc_->num_nodes();
  }

  //--------------------------------------------------------
  //  overloading operation to get number of columns
  //--------------------------------------------------------
  int NodalAtomVolume::nCols() const
  {
    return 1;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class NodalVolume
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  averaging_operation
  //--------------------------------------------------------
  double NodalVolume::averaging_operation(const double totalVolume) const
  {
    int nNodes = shapeFunction_->nCols();
    return totalVolume/nNodes;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class NodalAtomVolumeElement
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  NodalAtomVolumeElement::NodalAtomVolumeElement(ATC_Method * atc,
                                                 SPAR_MAN * shapeFunction,
                                                 PerAtomQuantity<int> * atomElement) :
    atc_(atc),
    shapeFunction_(shapeFunction),
    atomElement_(atomElement),
    feEngine_(atc->fe_engine()),
    tol_(1.e-10) 
  {
    shapeFunction_->register_dependence(this);
    if (!atomElement_) {
      InterscaleManager & interscaleManager = atc_->interscale_manager();
      atomElement_ = interscaleManager.per_atom_int_quantity("AtomElement");
    }
    atomElement_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void NodalAtomVolumeElement::reset_quantity() const
  {
    // Using analyses by G. Wagner and J. Templeton, weights ~ phi*M^{-1}*V
    //  where phi are the dimensionless shape/weighting functions, 
    //  M is the "mass" matrix M_IJ,
    //  V is the vector of nodal and element volumes
    //
    //
    // form atom-element shape functions and elemental volumes
    const FE_Mesh * feMesh = feEngine_->fe_mesh();
    int nElts = feMesh->num_elements();
    int nNodes = shapeFunction_->nCols();
    int nLocal = shapeFunction_->nRows();

    // form "mass" matrix M_IJ (I,J = 0, 1, ..., nNodes+nElts-1)
    int neSize = nNodes+nElts;
    const INT_ARRAY & atomElement(atomElement_->quantity());
    SPAR_MAT nodEltShpFcnMatrix(nLocal,neSize);
    const SPAR_MAT & shapeFunction(shapeFunction_->quantity());
    for(int a = 0; a < nLocal; a++) {
      for(int I = 0; I < nNodes; I++) {
        nodEltShpFcnMatrix.set(a,I,shapeFunction(a,I));
      }
      int thisCol = nNodes+atomElement(a,0);
      nodEltShpFcnMatrix.set(a,thisCol,1);
    }
      
    SPAR_MAT neMassMatrix(neSize,neSize);
    atc_->compute_consistent_md_mass_matrix(nodEltShpFcnMatrix,neMassMatrix);
      
    // form vector of nodal and elemental volumes 
    _nodeVolumesMatrix_.resize(nNodes,nNodes);
    feEngine_->compute_lumped_mass_matrix(_nodeVolumesMatrix_);
    _nodeVolumes_.resize(nNodes);
    _nodeVolumes_.copy(_nodeVolumesMatrix_.ptr(),_nodeVolumesMatrix_.size(),1);
    DENS_VEC _nodeElementVolumes_(neSize);
    for(int I = 0; I < nNodes; I++) {
      _nodeElementVolumes_(I) = _nodeVolumes_(I);
    }
    double averageEltVolume = 0.0;
    for(int E = 0; E < nElts; E++) {
      double minx, maxx, miny, maxy, minz, maxz;
      feMesh->element_coordinates(E,_nodalCoords_);
      minx = _nodalCoords_(0,0); maxx = _nodalCoords_(0,0);
      miny = _nodalCoords_(1,0); maxy = _nodalCoords_(1,0);
      minz = _nodalCoords_(2,0); maxz = _nodalCoords_(2,0);
      for (int j = 1; j < _nodalCoords_.nCols(); ++j) {
        if (_nodalCoords_(0,j)<minx) minx = _nodalCoords_(0,j);
        if (_nodalCoords_(0,j)>maxx) maxx = _nodalCoords_(0,j);
        if (_nodalCoords_(1,j)<miny) miny = _nodalCoords_(1,j);
        if (_nodalCoords_(1,j)>maxy) maxy = _nodalCoords_(1,j);
        if (_nodalCoords_(2,j)<minz) minz = _nodalCoords_(2,j);
        if (_nodalCoords_(2,j)>maxz) maxz = _nodalCoords_(2,j);
      }
      _nodeElementVolumes_(nNodes+E) = (maxx-minx)*(maxy-miny)*(maxz-minz);
      averageEltVolume += (maxx-minx)*(maxy-miny)*(maxz-minz);
    }
    averageEltVolume /= nElts;
      
    // correct entries of mass matrix if no atoms in shape function support
    double totalNodalVolume = _nodeVolumes_.sum();
    double averageNodalVolume = totalNodalVolume/nNodes;
    _scale_.resize(neSize);
    for (int i = 0; i < neSize; i++) {
      if ((abs(neMassMatrix(i,i)) > 0.)) { 
        _scale_(i) = 1.;
      } else {
        printf("No atoms are in support of node/element %i\n",i);
        _scale_(i) = 0.;
      }
    }
    neMassMatrix.row_scale(_scale_);
    for (int i = 0; i < neSize; i++) {
      if (_scale_(i) < 0.5) {
        neMassMatrix.set(i,i,1.);
        if (i < nNodes) { 
          _nodeElementVolumes_(i) = averageNodalVolume;
        } else {
          _nodeElementVolumes_(i) = averageEltVolume;
        }
      }
    }
    neMassMatrix.compress();
      
    // solve equation
    LinearSolver solver(neMassMatrix, ATC::LinearSolver::ITERATIVE_SOLVE_SYMMETRIC, true);
    solver.set_max_iterations(neMassMatrix.nRows());
    double myTol = 1.e-10; 
    solver.set_tolerance(myTol);
    quantity_.resize(neSize,0);
    CLON_VEC tempQuantity(quantity_,CLONE_COL,0);
    solver.solve(tempQuantity,_nodeElementVolumes_);
  }

  //--------------------------------------------------------
  //  overloading operation to get number of rows
  //--------------------------------------------------------
  int NodalAtomVolumeElement::nRows() const
  {
    return atc_->num_nodes();
  }

  //--------------------------------------------------------
  //  overloading operation to get number of columns
  //--------------------------------------------------------
  int NodalAtomVolumeElement::nCols() const
  {
    return 1;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomTypeElement
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomTypeElement::AtomTypeElement(ATC_Coupling * atc,
                                   PerAtomQuantity<int> * atomElement) :
    atomElement_(atomElement),
    nElts_((atc->fe_engine())->num_elements())
  {
    if (!atomElement_) {
      atomElement_ = (atc->interscale_manager()).per_atom_int_quantity("AtomElement");
    }
    atomElement_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtomTypeElement::reset_quantity() const
  {
    // determine which elements contain internal atoms
    quantity_.resize(nElts_,1);
    _quantityLocal_.resize(nElts_,1);
    _quantityLocal_ = 0;
    const INT_ARRAY & atomElement(atomElement_->quantity());
    for (int i = 0; i < atomElement_->nRows(); ++i) {
      _quantityLocal_(atomElement(i,0),0) = 1;
    }
    // swap contributions
    lammpsInterface_->logical_or(_quantityLocal_.ptr(),
                                 quantity_.ptr(),nElts_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElementMask
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ElementMask::ElementMask(ATC_Coupling * atc,
                           MatrixDependencyManager<DenseMatrix, int> * hasInternal,
                           MatrixDependencyManager<DenseMatrix, int> * hasGhost) :
    hasInternal_(hasInternal),
    hasGhost_(hasGhost),
    feEngine_(atc->fe_engine())
  {
    if (!hasInternal_) {
      hasInternal_ = (atc->interscale_manager()).dense_matrix_int("ElementHasInternal");
    }
    if (!hasGhost_) {
      hasGhost_ = (atc->interscale_manager()).dense_matrix_int("ElementHasGhost");
    }
    
    hasInternal_->register_dependence(this);
    if (hasGhost_) hasGhost_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void ElementMask::reset_quantity() const
  {
    const INT_ARRAY & hasInternal(hasInternal_->quantity());
    int nElts = hasInternal.size();
    quantity_.resize(nElts,1);

    if (hasGhost_) {
      const INT_ARRAY & hasGhost(hasGhost_->quantity());
      for (int i = 0; i < nElts; ++i) {
        quantity_(i,0) = !hasInternal(i,0) || hasGhost(i,0);
      }
    }
    else {
      for (int i = 0; i < nElts; ++i) {
        quantity_(i,0) = !hasInternal(i,0);
      }
    }
      
    const set<int> & nullElements = feEngine_->null_elements();
    set<int>::const_iterator iset;
    for (iset = nullElements.begin(); iset != nullElements.end(); iset++) {
      int ielem = *iset;
      quantity_(ielem,0) = false;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ElementMaskNodeSet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ElementMaskNodeSet::ElementMaskNodeSet(ATC_Coupling * atc,
                                         RegulatedNodes * nodeSet) :
    nodeSet_(nodeSet),
    feMesh_((atc->fe_engine())->fe_mesh())
  { 
    nodeSet_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void ElementMaskNodeSet::reset_quantity() const
  {
    quantity_.resize(feMesh_->num_elements(),1);
    quantity_ = false;

    // get the maximal element set corresponding to those nodes
    set<int> elementSet;
    feMesh_->nodeset_to_maximal_elementset(nodeSet_->quantity(),elementSet);

    set<int>::const_iterator iset;
    for (iset = elementSet.begin(); iset != elementSet.end(); iset++) {
      quantity_(*iset,0) = true;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class NodalGeometryType
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  NodalGeometryType::NodalGeometryType(ATC_Coupling * atc,
                                       MatrixDependencyManager<DenseMatrix, int> * hasInternal,
                                       MatrixDependencyManager<DenseMatrix, int> * hasGhost) :
    hasInternal_(hasInternal),
    hasGhost_(hasGhost),
    feEngine_(atc->fe_engine()),
    nNodes_(atc->num_nodes()),
    nElts_((atc->fe_engine())->num_elements())
  {
    if (!hasInternal_) {
      hasInternal_ = (atc->interscale_manager()).dense_matrix_int("ElementHasInternal");
    }
    if (!hasGhost_) {
      hasGhost_ = (atc->interscale_manager()).dense_matrix_int("ElementHasGhost");
    }

    hasInternal_->register_dependence(this);
    if (hasGhost_) hasGhost_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void NodalGeometryType::reset_quantity() const
  {
    const INT_ARRAY & hasInternal(hasInternal_->quantity());
    _nodesInternal_.resize(nNodes_);
    _nodesInternal_ = 0;
    _nodesGhost_.reset(nNodes_);
    _nodesGhost_ = 0;
    Array<int> nodes;
    
    
    
    vector<int> myElems = feEngine_->fe_mesh()->owned_elts();
    if (hasGhost_) {
      const INT_ARRAY & hasGhost(hasGhost_->quantity()) ;
      // iterate through all elements owned by this processor
      
      
      for (vector<int>::iterator elemsIter = myElems.begin();
           elemsIter != myElems.end();
           ++elemsIter)
      {
        int ielem = *elemsIter;
        if (hasInternal(ielem,0) || hasGhost(ielem,0)) {
          feEngine_->element_connectivity(ielem,nodes);
          for (int j = 0; j < nodes.size(); j++) {
            if (hasInternal(ielem,0)) {
              _nodesInternal_(nodes(j)) = 1;
            }
            if (hasGhost(ielem,0)) {
              _nodesGhost_(nodes(j)) = 1;
            }
          }
        }
      }
      // sum up partial result arrays
      
      lammpsInterface_->logical_or(MPI_IN_PLACE, _nodesInternal_.ptr(), _nodesInternal_.size());
      lammpsInterface_->logical_or(MPI_IN_PLACE, _nodesGhost_.ptr(), _nodesGhost_.size());
    }
    else {
      // iterate through all elements owned by this processor
      for (vector<int>::iterator elemsIter = myElems.begin();
           elemsIter != myElems.end();
           ++elemsIter)
      {
        int ielem = *elemsIter;
        if (hasInternal(ielem,0)) {
          feEngine_->element_connectivity(ielem,nodes);
          for (int j = 0; j < nodes.size(); j++) {
            _nodesInternal_(nodes(j)) = 1;
          }
        }
      }
      // sum up partial result arrays
      lammpsInterface_->logical_or(MPI_IN_PLACE, _nodesInternal_.ptr(), _nodesInternal_.size());
    }

    quantity_.resize(nNodes_,1);
    for (int i = 0; i < nNodes_; i++) {
      if (_nodesInternal_(i) && _nodesGhost_(i)) {
        quantity_(i,0) = BOUNDARY;
      }
      else if (_nodesInternal_(i)) {
        quantity_(i,0) = MD_ONLY;
      }
      else {
        quantity_(i,0) = FE_ONLY;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class NodeToSubset
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  NodeToSubset::NodeToSubset(ATC_Method * atc,
                             SetDependencyManager<int> * subsetNodes) :
    LargeToSmallMap(),
    atc_(atc),
    subsetNodes_(subsetNodes)
  {
    subsetNodes_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void NodeToSubset::reset_quantity() const
  {
    int nNodes = atc_->num_nodes();
    const set<int> & subsetNodes(subsetNodes_->quantity());
    quantity_.resize(nNodes,1);
    size_ = 0;

    for (int i = 0; i < nNodes; i++) {
      if (subsetNodes.find(i) != subsetNodes.end()) {
        quantity_(i,0) = size_;
        size_++;
      }
      else {
        quantity_(i,0) = -1;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SubsetToNode
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  SubsetToNode::SubsetToNode(NodeToSubset * nodeToSubset) :
    nodeToSubset_(nodeToSubset)
  {
    nodeToSubset_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void SubsetToNode::reset_quantity() const
  {
    const INT_ARRAY & nodeToSubset(nodeToSubset_->quantity());
    int nNodes = nodeToSubset.nRows();
    int count = 0;
    quantity_.resize(nodeToSubset_->size(),1);

    for (int i = 0; i < nNodes; i++) {
      if (nodeToSubset(i,0) > -1) {
        quantity_(count,0) = i;
        count++;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ReducedSparseMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ReducedSparseMatrix::ReducedSparseMatrix(ATC_Method * atc,
                                           SPAR_MAN * source,
                                           LargeToSmallAtomMap * map) :
    SparseMatrixTransfer<double>(),
    source_(source),
    map_(map)
  {
    source_->register_dependence(this);
    map_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  ReducedSparseMatrix::~ReducedSparseMatrix()
  {
    source_->remove_dependence(this);
    map_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void ReducedSparseMatrix::reset_quantity() const
  {
    const SPAR_MAT & source(source_->quantity());
    const INT_ARRAY & map(map_->quantity());
    quantity_.reset(source.nRows(),source.nCols());

    for (int i = 0; i < source.nRows(); i++) {
      int idx = map(i,0);
      if (idx > -1) {
        source.row(i,_row_,_index_);
        for (int j = 0; j < _row_.size(); j++) {
          quantity_.set(i,_index_(j),_row_(j));
        }
      }
    }
    quantity_.compress();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RowMappedSparseMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void RowMappedSparseMatrix::reset_quantity() const
  {
    const SPAR_MAT & source(source_->quantity());
    const INT_ARRAY & map(map_->quantity());
    quantity_.reset(map_->size(),source.nCols());

    for (int i = 0; i < source.nRows(); i++) {
      int idx = map(i,0);
      if (idx > -1) {
        source.row(i,_row_,_index_);
        for (int j = 0; j < _row_.size(); j++) {
          quantity_.set(idx,_index_(j),_row_(j));
        }
      }
    }
    quantity_.compress();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RowMappedSparseMatrixVector
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  RowMappedSparseMatrixVector::RowMappedSparseMatrixVector(ATC_Method * atc,
                                                           VectorDependencyManager<SPAR_MAT * > * source,
                                                           LargeToSmallAtomMap * map) :
    VectorTransfer<SPAR_MAT * >(),
    source_(source),
    map_(map)
  {
    source_->register_dependence(this);
    map_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  RowMappedSparseMatrixVector::~RowMappedSparseMatrixVector()
  {
    source_->remove_dependence(this);
    map_->remove_dependence(this);
    for (unsigned i = 0; i < quantity_.size(); ++i) {
      if (quantity_[i]) delete quantity_[i];
    }
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void RowMappedSparseMatrixVector::reset_quantity() const
  {
    const vector<SPAR_MAT * > & source(source_->quantity());
    const INT_ARRAY & map(map_->quantity());

    for (unsigned i = 0; i < quantity_.size(); ++i) {
      if (quantity_[i]) delete quantity_[i];
    }
    quantity_.resize(source.size(),NULL);
    for (unsigned i = 0; i < source.size(); i++) {
      quantity_[i] = new SPAR_MAT(map_->size(),source[i]->nCols());
    }

    for (unsigned i = 0; i < source.size(); i++) {
      for (int j = 0; j < source[i]->nRows(); j++) {
        int idx = map(j,0);
        if (idx > -1) {
          source[i]->row(j,_row_,_index_);
          for (int k = 0; k < _row_.size(); k++) {
            quantity_[i]->set(idx,_index_(k),_row_(k));
          }
        }
      }
    }
    for (unsigned i = 0; i < source.size(); i++) {
      quantity_[i]->compress();
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MappedDiagonalMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  MappedDiagonalMatrix::MappedDiagonalMatrix(ATC_Method * atc,
                                             DIAG_MAN * source,
                                             LargeToSmallAtomMap * map) :
    DiagonalMatrixTransfer<double>(),
    source_(source),
    map_(map)
  {
    source_->register_dependence(this);
    map_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  MappedDiagonalMatrix::~MappedDiagonalMatrix()
  {
    source_->remove_dependence(this);
    map_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void MappedDiagonalMatrix::reset_quantity() const
  {
    const DIAG_MAT & source(source_->quantity());
    const INT_ARRAY & map(map_->quantity());
    quantity_.resize(map_->size(),map_->size());

    for (int i = 0; i < source.nRows(); i++) {
      int idx = map(i,0);
      if (idx > -1) {
        quantity_(idx,idx) = source(i,i);
      }
    }
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MappedQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  MappedQuantity::MappedQuantity(ATC_Method * atc,
                                 DENS_MAN * source,
                                 LargeToSmallMap * map) :
    DenseMatrixTransfer<double>(),
    source_(source),
    map_(map)
  {
    source_->register_dependence(this);
    map_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void MappedQuantity::reset_quantity() const
  {
    const DENS_MAT & source(source_->quantity());
    const INT_ARRAY & map(map_->quantity());
    quantity_.resize(map_->size(),source.nCols());

    for (int i = 0; i < source.nRows(); i++) {
      int idx = map(i,0);
      if (idx > -1) {
        for (int j = 0; j < source.nCols(); j++) {
          quantity_(idx,j) = source(i,j);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class RegulatedNodes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  RegulatedNodes::RegulatedNodes(ATC_Coupling * atc,
                                 FieldName fieldName,
                                 MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType) :
    SetTransfer<int>(),
    nodalGeometryType_(nodalGeometryType),
    atc_(atc),
    feEngine_(atc->fe_engine()),
    prescribedDataManager_(atc->prescribed_data_manager())
  {
    if (!nodalGeometryType_) {
      nodalGeometryType_ = (atc_->interscale_manager()).dense_matrix_int("NodalGeometryType");
    }
    if (nodalGeometryType_) {
      nodalGeometryType_->register_dependence(this);
    }
    else {
      throw ATC_Error("TransferLibrary::RegulatedNodes - No Nodal Geometry Type provided");
    }

    const map<FieldName,int> & atcFieldSizes(atc_->field_sizes());
    if (fieldName == NUM_TOTAL_FIELDS) {
      fieldSizes_ = atcFieldSizes;
    }
    else {
      map<FieldName,int>::const_iterator fs_iter = atcFieldSizes.find(fieldName);
      fieldSizes_.insert(pair<FieldName,int>(fieldName,fs_iter->second));
    }
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void RegulatedNodes::reset_quantity() const
  {
    const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
    quantity_.clear();
    
    for (int i = 0; i < nodeType.size(); ++i) {
      if (nodeType(i,0) != FE_ONLY) {
        quantity_.insert(i);
      }
    }
  }

  //--------------------------------------------------------
  //  insert_boundary_nodes
  //--------------------------------------------------------
  void RegulatedNodes::insert_boundary_nodes() const
  {
    const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
    
    for (int i = 0; i < nodeType.size(); ++i) {
      if (nodeType(i,0) == BOUNDARY) {
        quantity_.insert(i);
      }
    }
  }

  //--------------------------------------------------------
  //  insert_boundary_faces
  //--------------------------------------------------------
  void RegulatedNodes::insert_boundary_faces() const
  {
    const set<string> & boundaryFaceNames(atc_->boundary_face_names());
    set<string>::const_iterator iter;
    set<int>::const_iterator nodeIter;

    for (iter = boundaryFaceNames.begin(); iter != boundaryFaceNames.end(); iter++) {
      set<int> nodeSet;
      feEngine_->fe_mesh()->faceset_to_nodeset(*iter,nodeSet);
      for (nodeIter = nodeSet.begin(); nodeIter != nodeSet.end(); nodeIter++) {
        quantity_.insert(*nodeIter);
      }
    }
  }

  //--------------------------------------------------------
  //  insert_fixed_nodes
  //--------------------------------------------------------
  void RegulatedNodes::insert_fixed_nodes() const
  {
    
    const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
    map<FieldName,int>::const_iterator fs_iter;

    for (int i = 0; i < nodeType.size(); ++i) {
      bool isFixed = false;
      for (fs_iter = fieldSizes_.begin(); fs_iter != fieldSizes_.end(); fs_iter++) {
        for (int j = 0; j < fs_iter->second; j++) {
          isFixed = prescribedDataManager_->is_fixed(i,fs_iter->first,j);
          if (isFixed) break;
        }
      }
      if (isFixed && ((nodeType(i,0)==MD_ONLY) || (nodeType(i,0)==BOUNDARY))) {
        quantity_.insert(i);
      }
    }
  }

  //--------------------------------------------------------
  //  insert_face_fluxes
  //--------------------------------------------------------
  void RegulatedNodes::insert_face_fluxes() const
  {
    const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
    set<int>::const_iterator inode;
    map<FieldName,int>::const_iterator fs_iter;
    
    for (fs_iter = fieldSizes_.begin(); fs_iter != fieldSizes_.end(); fs_iter++) {
      for (int j = 0; j < fs_iter->second; j++) {
        set<int> faceFluxNodes = prescribedDataManager_->flux_face_nodes(fs_iter->first,j);
        for (inode = faceFluxNodes.begin(); inode != faceFluxNodes.end(); inode++) {
          if (((nodeType(*inode,0)==MD_ONLY) || (nodeType(*inode,0)==BOUNDARY))) {
            quantity_.insert(*inode);
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //  insert_element_fluxes
  //--------------------------------------------------------
  void RegulatedNodes::insert_element_fluxes() const
  {
    const INT_ARRAY & nodeType(nodalGeometryType_->quantity());
    set<int>::const_iterator inode;
    map<FieldName,int>::const_iterator fs_iter;

    for (fs_iter = fieldSizes_.begin(); fs_iter != fieldSizes_.end(); fs_iter++) {
      for (int j = 0; j < fs_iter->second; j++) {
        set<int> elementFluxNodes = prescribedDataManager_->flux_element_nodes(fs_iter->first,j);
        for (inode = elementFluxNodes.begin(); inode != elementFluxNodes.end(); inode++) {
          if (((nodeType(*inode,0)==MD_ONLY) || (nodeType(*inode,0)==BOUNDARY))) {
            quantity_.insert(*inode);
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class BoundaryNodes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void BoundaryNodes::reset_quantity() const
  {
    quantity_.clear();

    // a) they are a boundary node
    RegulatedNodes::insert_boundary_nodes();

    // b) they are in a specified boundary faceset
    RegulatedNodes::insert_boundary_faces();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FluxNodes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void FluxNodes::reset_quantity() const
  {
    quantity_.clear();
    
    // a) they have a fixed face flux
    RegulatedNodes::insert_face_fluxes();

    // b) they have a fixed element flux
    RegulatedNodes::insert_element_fluxes();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FluxBoundaryNodes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void FluxBoundaryNodes::reset_quantity() const
  {
    FluxNodes::reset_quantity();

    // a) they are a boundary node
    RegulatedNodes::insert_boundary_nodes();

    // b) they are in a specified boundary faceset
    RegulatedNodes::insert_boundary_faces();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AllRegulatedNodes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AllRegulatedNodes::reset_quantity() const
  {
    FluxBoundaryNodes::reset_quantity();

    // a) they are a fixed node
    RegulatedNodes::insert_fixed_nodes();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FixedNodes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void FixedNodes::reset_quantity() const
  {
    quantity_.clear();
    
    // a) they are a fixed node
    RegulatedNodes::insert_fixed_nodes();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FixedBoundaryNodes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void FixedBoundaryNodes::reset_quantity() const
  {
    FixedNodes::reset_quantity();

    // a) they are a boundary node
    RegulatedNodes::insert_boundary_nodes();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DenseMatrixQuotient
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  DenseMatrixQuotient::DenseMatrixQuotient(DENS_MAN* matrixNumerator,
                                           DENS_MAN* matrixDenominator):
  matrixNumerator_(matrixNumerator),
  matrixDenominator_(matrixDenominator)
  {
    matrixNumerator_->register_dependence(this);
    matrixDenominator_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Reset_quantity
  //--------------------------------------------------------
  void DenseMatrixQuotient::reset_quantity() const
  {
      quantity_ = matrixNumerator_->quantity();
      quantity_.divide_zero_safe(matrixDenominator_->quantity());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DenseMatrixSum
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  DenseMatrixSum::DenseMatrixSum(DENS_MAN* matrixOne,
                                 DENS_MAN* matrixTwo):
  matrixOne_(matrixOne), matrixTwo_(matrixTwo)
  {
    matrixOne_->register_dependence(this);
    matrixTwo_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Reset_quantity
  //--------------------------------------------------------
  void DenseMatrixSum::reset_quantity() const
  {
      quantity_ = matrixOne_->quantity();
      quantity_ += (matrixTwo_->quantity());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DenseMatrixDelta
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  DenseMatrixDelta::DenseMatrixDelta(DENS_MAN* matrix,
                                     DENS_MAT* reference):
  matrix_(matrix), reference_(reference)
  {
    matrix_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Reset_quantity
  //--------------------------------------------------------
  void DenseMatrixDelta::reset_quantity() const
  {
      quantity_ = matrix_->quantity();
      quantity_ -= *reference_;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class PointToElementMap
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  PointToElementMap::PointToElementMap(ATC_Method * atc,
                                       MatrixDependencyManager<DenseMatrix, double> * pointPositions) :
    DenseMatrixTransfer<int>(),
    pointPositions_(pointPositions),
    feMesh_((atc->fe_engine())->fe_mesh())
  {
    pointPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  PointToElementMap::~PointToElementMap()
  {
    pointPositions_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void PointToElementMap::reset_quantity() const
  {
    const DENS_MAT & position(pointPositions_->quantity());
    int nsd = position.nCols();
    int nRows = position.nRows();
    quantity_.resize(nRows,nsd);

    DENS_VEC coords(nsd);
    for (int i = 0; i < nRows; i++) {
      for (int j = 0; j < nsd; j++) {
        coords(j) = position(i,j);
      }
      quantity_(i,0) = feMesh_->map_to_element(coords);
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class Interpolant
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  Interpolant::Interpolant(ATC_Method * atc,
                           MatrixDependencyManager<DenseMatrix, int>* pointToElementMap,
                           DENS_MAN* pointPositions) :
    SparseMatrixTransfer<double>(),
    pointToElementMap_(pointToElementMap),
    pointPositions_(pointPositions),
    feEngine_(atc->fe_engine())
  {
    pointToElementMap_->register_dependence(this);
    pointPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void Interpolant::reset_quantity() const
  {
    const DENS_MAT & positions(pointPositions_->quantity());
    const INT_ARRAY & elements(pointToElementMap_->quantity());
    if (positions.nRows() > 0) {
      feEngine_->evaluate_shape_functions(positions,
                                          elements,
                                          this->quantity_);
    }
    else {
      quantity_.resize(0,feEngine_->num_nodes());
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class InterpolantGradient
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  InterpolantGradient::InterpolantGradient(ATC_Method * atc,
                                           MatrixDependencyManager<DenseMatrix, int>* pointToElementMap,
                                           DENS_MAN* pointPositions) :
    VectorTransfer<SPAR_MAT * >(),
    pointToElementMap_(pointToElementMap),
    pointPositions_(pointPositions),
    feEngine_(atc->fe_engine())
  {
    pointToElementMap_->register_dependence(this);
    pointPositions_->register_dependence(this);
    quantity_.resize(atc->nsd(),NULL);
    for (int i = 0; i < atc->nsd(); ++i) {
      quantity_[i] = new SPAR_MAT();
    }
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  InterpolantGradient::~InterpolantGradient() {
    pointToElementMap_->remove_dependence(this);
    pointPositions_->remove_dependence(this);
    for (unsigned i = 0; i < quantity_.size(); ++i) {
      if (quantity_[i]) delete quantity_[i];
    }
  }

  //--------------------------------------------------------
  //  reset_quantity()
  //--------------------------------------------------------
  void InterpolantGradient::reset_quantity() const
  {
    const DENS_MAT & positions(pointPositions_->quantity());
    const INT_ARRAY & elements(pointToElementMap_->quantity());
    if (positions.nRows() > 0) {
      feEngine_->evaluate_shape_function_derivatives(positions,
                                                     elements,
                                                     this->quantity_);
    }
    else {
      for (unsigned i = 0; i < quantity_.size(); ++i) {
        (this->quantity_)[i]->resize(0,feEngine_->num_nodes());
      }
    }
  }

  //--------------------------------------------------------
  //  Class PerAtomShapeFunctionGradient
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  PerAtomShapeFunctionGradient::PerAtomShapeFunctionGradient(ATC_Method * atc,
                                                             MatrixDependencyManager<DenseMatrix, int>* atomToElementMap,
                                                             DENS_MAN* atomPositions,
                                                             const string & tag,
                                                             AtomType atomType) :
    VectorTransfer<SPAR_MAT * >(),
    atomToElementMap_(atomToElementMap),
    atomPositions_(atomPositions),
    feEngine_(atc->fe_engine())
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomToElementMap_) {
      atomToElementMap_ = interscaleManager.per_atom_int_quantity("AtomElement");
    }
    if (!atomPositions_) {
      atomPositions_ = interscaleManager.per_atom_quantity("AtomicCoarseGrainingPositions");
    }
    atomToElementMap_->register_dependence(this);
    atomPositions_->register_dependence(this);

    // storage container
    matrices_.resize(atc->nsd(),NULL);
    for (int i = 0; i < atc->nsd(); ++i) {
      matrices_[i] = new AtcAtomSparseMatrix<double>(atc,feEngine_->num_nodes(),
                                                     feEngine_->num_nodes_per_element(),
                                                     atomType);
      stringstream myint;
      myint << i;
      interscaleManager.add_per_atom_sparse_matrix(matrices_[i],
                                                   tag+myint.str());
      matrices_[i]->register_dependence(this);
    }

    quantity_.resize(atc->nsd(),NULL);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  PerAtomShapeFunctionGradient::~PerAtomShapeFunctionGradient() {
    atomToElementMap_->remove_dependence(this);
    atomPositions_->remove_dependence(this);
    for (unsigned i = 0; i < matrices_.size(); ++i) {
      matrices_[i]->remove_dependence(this);
    }
  }

  //--------------------------------------------------------
  //  reset_quantity()
  //--------------------------------------------------------
  void PerAtomShapeFunctionGradient::reset_quantity() const
  {
    const DENS_MAT & positions(atomPositions_->quantity());
    const INT_ARRAY & elements(atomToElementMap_->quantity());
    for (unsigned i = 0; i < quantity_.size(); ++i) {
      (this->quantity_)[i] = & matrices_[i]->set_quantity();
    }

    if (positions.nRows() > 0) {
      feEngine_->evaluate_shape_function_derivatives(positions,
                                                     elements,
                                                     this->quantity_);
    }
    else {
      for (unsigned i = 0; i < quantity_.size(); ++i) {
        (this->quantity_)[i]->resize(0,feEngine_->num_nodes());
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class InterpolantSmallMolecule
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  InterpolantSmallMolecule::InterpolantSmallMolecule(ATC_Method * atc,
                                                     MatrixDependencyManager<DenseMatrix, int>* moleculeToElementMap,
                                                     DENS_MAN* moleculePositions,
                                                     MoleculeSet * moleculeSet) :
    Interpolant(atc,moleculeToElementMap,moleculePositions),
    moleculeSet_(moleculeSet)
  {
    moleculeSet_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  InterpolantSmallMolecule::~InterpolantSmallMolecule()
  {
    moleculeSet_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void InterpolantSmallMolecule::reset_quantity() const
  {
    Interpolant::reset_quantity();

    // scale rows by fraction of molecules on this proc
    _fractions_.resize((this->quantity_).nRows());
    for (int i = 0; i < moleculeSet_->local_molecule_count(); i++)
      _fractions_(i) = moleculeSet_->local_fraction(i);
    (this->quantity_).row_scale(_fractions_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyKernelAccumulation
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyKernelAccumulation::OnTheFlyKernelAccumulation(ATC_Method * atc,
                                                         PerAtomQuantity<double> * source,
                                                         KernelFunction* kernelFunction,
                                                         DENS_MAN* atomCoarseGrainingPositions):
    atc_(atc),
    source_(source),
    kernelFunction_(kernelFunction), 
    atomCoarseGrainingPositions_(atomCoarseGrainingPositions),
    feMesh_((atc_->fe_engine())->fe_mesh())
  {
    source_->register_dependence(this);
    atomCoarseGrainingPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyKernelAccumulation::reset_quantity() const
  {
    const DENS_MAT & positions(atomCoarseGrainingPositions_->quantity());
    const DENS_MAT & source(source_->quantity());
    int nNodes = feMesh_->num_nodes_unique();
    quantity_.resize(nNodes,source.nCols());
    _quantityLocal_.reset(nNodes,source.nCols());
    
    if (source.nRows()>0) {
      DENS_VEC xI(positions.nCols()),xa(positions.nCols()),xaI(positions.nCols());
      double val;
      for (int i = 0; i < nNodes; i++) {
        xI = feMesh_->nodal_coordinates(i);
        for (int j = 0; j < positions.nRows(); j++) {
          for (int k = 0; k < positions.nCols(); ++k) {
            xa(k) = positions(j,k);
          }
          xaI = xa - xI;
          atc_->lammps_interface()->periodicity_correction(xaI.ptr());
          val = kernelFunction_->value(xaI);
          if (val > 0) {
            for (int k=0; k < source.nCols(); k++) {
              _quantityLocal_(i,k) += val*source(j,k);
            }
          }
        }
      }
    }

    // accumulate across processors
    int count = quantity_.nRows()*quantity_.nCols();
    lammpsInterface_->allsum(_quantityLocal_.ptr(),quantity_.ptr(),count);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyKernelAccumulationNormalized
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyKernelAccumulationNormalized::OnTheFlyKernelAccumulationNormalized(ATC_Method * atc,
                   PerAtomQuantity<double> * source,
                   KernelFunction* kernelFunction,
                   DENS_MAN* atomCoarseGrainingPositions,
                   DIAG_MAN * weights):
    OnTheFlyKernelAccumulation(atc,source,kernelFunction,atomCoarseGrainingPositions),
    weights_(weights)
  {
    weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyKernelAccumulationNormalized::reset_quantity() const
  {
    OnTheFlyKernelAccumulation::reset_quantity();
    quantity_ *= weights_->quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyKernelAccumulationNormalizedReferenced
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyKernelAccumulationNormalizedReferenced::OnTheFlyKernelAccumulationNormalizedReferenced(ATC_Method * atc,
                PerAtomQuantity<double> * source,
                KernelFunction* kernelFunction,
                DENS_MAN* atomCoarseGrainingPositions,
                DIAG_MAN* weights,
                const DENS_MAT * reference):
    OnTheFlyKernelAccumulationNormalized(atc,source,kernelFunction,atomCoarseGrainingPositions,weights),
    reference_(reference)
  {
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyKernelAccumulationNormalizedReferenced::reset_quantity() const
  {
    OnTheFlyKernelAccumulationNormalized::reset_quantity();
    quantity_ -= *reference_;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyKernelAccumulationNormalizedScaled
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyKernelAccumulationNormalizedScaled::OnTheFlyKernelAccumulationNormalizedScaled(ATC_Method * atc,
        PerAtomQuantity<double> * source,
        KernelFunction* kernelFunction,
        DENS_MAN* atomCoarseGrainingPositions,
        DIAG_MAN* weights,
        const double scale):
    OnTheFlyKernelAccumulationNormalized(atc,source,kernelFunction,atomCoarseGrainingPositions,weights),
    scale_(scale)
  {
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyKernelAccumulationNormalizedScaled::reset_quantity() const
  {
    OnTheFlyKernelAccumulationNormalized::reset_quantity();
    quantity_ *= scale_;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AccumulantWeights
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AccumulantWeights::AccumulantWeights(SPAR_MAN* accumulant):
    DiagonalMatrixTransfer<double>(),
    accumulant_(accumulant)
  {
    accumulant_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Reset_quantity
  //--------------------------------------------------------
  void AccumulantWeights::reset_quantity() const
  {
    const SPAR_MAT accumulant(accumulant_->quantity());
    int nNodes = accumulant.nCols();

    // get summation of atoms per node
    _localWeights_.reset(nNodes); _weights_.resize(nNodes);
    if (accumulant.nRows()>0) {
      _localWeights_ = (accumulant_->quantity()).col_sum();
    }
    lammpsInterface_->allsum(_localWeights_.ptr(),_weights_.ptr(),nNodes);
    
    // assign weights
    quantity_.resize(nNodes,nNodes);
    for (int i = 0; i < nNodes; i++) {
      if (_weights_(i) > 0.) {
        quantity_(i,i) = 1./_weights_(i);
      }
      else {
        quantity_(i,i) = 0.;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyKernelWeights
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyKernelWeights::OnTheFlyKernelWeights(DENS_MAN* weights):
    DiagonalMatrixTransfer<double>(),
    weights_(weights)
  {
    weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Reset_quantity
  //--------------------------------------------------------
  void OnTheFlyKernelWeights::reset_quantity() const
  {
    const DENS_MAT & weights(weights_->quantity());
    int nNodes = weights.nRows();
    
    // assign weights
    quantity_.resize(nNodes,nNodes);
    for (int i = 0; i < nNodes; i++) {
      if (weights(i,0) > 0.) {
        quantity_(i,i) = 1./weights(i,0);
      }
      else {
        quantity_(i,i) = 0.;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KernelInverseVolumes
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  KernelInverseVolumes::KernelInverseVolumes(ATC_Method * atc,
                                             KernelFunction* kernelFunction):
    kernelFunction_(kernelFunction), 
    feMesh_((atc->fe_engine())->fe_mesh())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  multiplication by transpose
  //--------------------------------------------------------
  void KernelInverseVolumes::reset_quantity() const
  {
    int nNodes = feMesh_->num_nodes_unique();
    quantity_.resize(nNodes,nNodes);
    quantity_ = kernelFunction_->inv_vol();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyMeshAccumulation
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyMeshAccumulation::OnTheFlyMeshAccumulation(ATC_Method * atc,
                                                     PerAtomQuantity<double> * source,
                                                     DENS_MAN* atomCoarseGrainingPositions):
    atc_(atc),
    source_(source),
    atomCoarseGrainingPositions_(atomCoarseGrainingPositions),
    feMesh_((atc_->fe_engine())->fe_mesh())
  {
    source_->register_dependence(this);
    atomCoarseGrainingPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyMeshAccumulation::reset_quantity() const
  {
    const DENS_MAT & positions(atomCoarseGrainingPositions_->quantity());
    const DENS_MAT & source(source_->quantity());
    int nNodes = feMesh_->num_nodes_unique();
    int nodesPerElement = feMesh_->num_nodes_per_element();
    Array<int> node_list(nodesPerElement);
    DENS_VEC shp(nodesPerElement);
    quantity_.resize(nNodes,source.nCols());
    _quantityLocal_.reset(nNodes,source.nCols());
    DENS_VEC xj(atc_->nsd());
    
    if (source.nRows()>0) {
      for (int j = 0; j < source.nRows(); j++) {
        for (int k = 0; k < atc_->nsd(); k++) {
          xj(k) = positions(j,k);
        }
        feMesh_->shape_functions(xj,shp,node_list);
        for (int I = 0; I < nodesPerElement; I++) {
          int inode = node_list(I);
          for (int k = 0; k < source.nCols(); k++) {
            //quantity_(inode,k) += shp(I)*source(j,k);
            _quantityLocal_(inode,k) += shp(I)*source(j,k);
          }
        } 
      }
    }
    // accumulate across processors
    int count = quantity_.nRows()*quantity_.nCols();
    lammpsInterface_->allsum(_quantityLocal_.ptr(),quantity_.ptr(),count);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyMeshAccumulationNormalized
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyMeshAccumulationNormalized::OnTheFlyMeshAccumulationNormalized(ATC_Method * atc,
                   PerAtomQuantity<double> * source,
                   DENS_MAN* atomCoarseGrainingPositions,
                   DIAG_MAN * weights):
    OnTheFlyMeshAccumulation(atc,source,atomCoarseGrainingPositions),
    weights_(weights)
  {
    weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyMeshAccumulationNormalized::reset_quantity() const
  {
    OnTheFlyMeshAccumulation::reset_quantity();
    quantity_ *= weights_->quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyMeshAccumulationNormalizedReferenced
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyMeshAccumulationNormalizedReferenced::OnTheFlyMeshAccumulationNormalizedReferenced(ATC_Method * atc,
                PerAtomQuantity<double> * source,
                DENS_MAN* atomCoarseGrainingPositions,
                DIAG_MAN* weights,
                const DENS_MAT * reference):
    OnTheFlyMeshAccumulationNormalized(atc,source,atomCoarseGrainingPositions,weights),
    reference_(reference)
  {
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyMeshAccumulationNormalizedReferenced::reset_quantity() const
  {
    OnTheFlyMeshAccumulationNormalized::reset_quantity();
    quantity_ -= *reference_;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyMeshAccumulationNormalizedScaled
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyMeshAccumulationNormalizedScaled::OnTheFlyMeshAccumulationNormalizedScaled(ATC_Method * atc,
        PerAtomQuantity<double> * source,
        DENS_MAN* atomCoarseGrainingPositions,
        DIAG_MAN* weights,
        const double scale):
    OnTheFlyMeshAccumulationNormalized(atc,source,atomCoarseGrainingPositions,weights),
    scale_(scale)
  {
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyMeshAccumulationNormalizedScaled::reset_quantity() const
  {
    OnTheFlyMeshAccumulationNormalized::reset_quantity();
    quantity_ *= scale_;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class NativeShapeFunctionGradient
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  NativeShapeFunctionGradient::NativeShapeFunctionGradient(ATC_Method * atc) :
    VectorTransfer<SPAR_MAT * >(),
    feEngine_(atc->fe_engine())
  {
    quantity_.resize(atc->nsd(),NULL);
    for (int i = 0; i < atc->nsd(); ++i) {
      quantity_[i] = new SPAR_MAT();
    }
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  NativeShapeFunctionGradient::~NativeShapeFunctionGradient() {
    for (unsigned i = 0; i < quantity_.size(); ++i) {
      if (quantity_[i]) delete quantity_[i];
    }
  }

  //--------------------------------------------------------
  //  reset_quantity()
  //--------------------------------------------------------
  void NativeShapeFunctionGradient::reset_quantity() const
  {
    feEngine_->compute_gradient_matrix(quantity_); 
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class OnTheFlyShapeFunctionProlongation
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  OnTheFlyShapeFunctionProlongation::OnTheFlyShapeFunctionProlongation(ATC_Method * atc,
                                                   DENS_MAN * source,
                                                   DENS_MAN * atomCoarseGrainingPositions):
    FeToAtomTransfer(atc,source),
    atomCoarseGrainingPositions_(atomCoarseGrainingPositions),
    feMesh_((atc->fe_engine())->fe_mesh())
  {
    atomCoarseGrainingPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  // destructor
  //--------------------------------------------------------
  OnTheFlyShapeFunctionProlongation::~OnTheFlyShapeFunctionProlongation() 
  {
    atomCoarseGrainingPositions_->remove_dependence(this);
  };

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void OnTheFlyShapeFunctionProlongation::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & positions(atomCoarseGrainingPositions_->quantity());
      const DENS_MAT & source(source_->quantity());
      int nodesPerElement = feMesh_->num_nodes_per_element();
      Array<int> node_list(nodesPerElement);
      DENS_VEC shp(nodesPerElement);
      DENS_VEC xj(positions.nCols());
      int nLocal = positions.nRows();
      quantity_ = 0.;
      for (int j = 0; j < nLocal; j++) {
        for (int k = 0; k < source.nCols(); k++) {
          xj(k) = positions(j,k);
        }
        feMesh_->shape_functions(xj,shp,node_list);
        for (int I = 0; I < nodesPerElement; I++) {
          int inode = node_list(I);
          for (int k = 0; k < source.nCols(); k++) {
            quantity_(j,k) += shp(I)*source(inode,k);
          }
        }
      }
    }
  }
}

// a library for the various transfer operators used in ATC

#ifndef TRANSFER_LIBRARY_H
#define TRANSFER_LIBRARY_H

#include <string>
#include <map>
#include <vector>

// ATC_Method headers
#include "TransferOperator.h"

namespace ATC {

  // forward declarations
  class ATC_Method;
  class ATC_Coupling;
  class FE_Engine;
  class FE_Mesh;
  class PrescribedDataManager;
  class LargeToSmallAtomMap;
  class MoleculeSet;

  /**
   *  @class  NodalAtomVolume
   *  @brief  Computes the nodal volumes which coarse grain the volume per atom
   */
  
  class NodalAtomVolume : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    NodalAtomVolume(ATC_Method * atc, SPAR_MAN * shapeFunction);
    
    // destructor
    virtual ~NodalAtomVolume() {shapeFunction_->remove_dependence(this);};

    /** apply transfer operator */
    virtual void reset_quantity() const;

    /** overload function to get the number rows */
    virtual int nRows() const;

    /** overload function to get the number columns */
    virtual int nCols() const;

  protected:

    /** applies the needed averaging operation to get the appropriate average volume */
    virtual double averaging_operation(const double totalVolume) const;

    /** pointer to atc method object */
    const ATC_Method * atc_;

    /** shape function matrix used to base volumes on */
    SPAR_MAN * shapeFunction_;

    /** pointer to lammps interface */
    const LammpsInterface * lammpsInterface_;

    /** pointer to the FE engine */
    const FE_Engine * feEngine_;

    /** tolerance used in the linear solver */
    const double tol_;

    // workspace variables
    mutable DIAG_MAT _rhsMatrix_;
    mutable DENS_VEC _rhs_;
    mutable DENS_VEC _scale_;

  private:

    // do not define
    NodalAtomVolume();

  };

  /**
   *  @class  NodalVolume
   *  @brief  Computes the nodal volumes associated with the support of each nodal shape function
   */
  
  class NodalVolume : public NodalAtomVolume {

  public:
    
    // constructor
    NodalVolume(ATC_Method * atc, SPAR_MAN * shapeFunction) : NodalAtomVolume(atc,shapeFunction) {};
    
    // destructor
    virtual ~NodalVolume() {};

  protected:

    /** applies the needed averaging operation to get the appropriate average volume */
    virtual double averaging_operation(const double totalVolume) const;

  private:

    // do not define
    NodalVolume();

  };

  /**
   *  @class  NodalAtomVolumeElement
   *  @brief  Computes the nodal volumes which coarse grain the volume per atom based on element volumes
   */
  
  class NodalAtomVolumeElement : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    NodalAtomVolumeElement(ATC_Method * atc, SPAR_MAN * shapeFunction,
                           PerAtomQuantity<int> * atomElement=NULL);
    
    // destructor
    virtual ~NodalAtomVolumeElement() {
      shapeFunction_->remove_dependence(this);
      atomElement_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

    /** overload function to get the number rows */
    virtual int nRows() const;

    /** overload function to get the number columns */
    virtual int nCols() const;

  protected:

    /** pointer to atc method object */
    ATC_Method * atc_;

    /** shape function matrix used to base volumes on */
    SPAR_MAN * shapeFunction_;

    /** object containing the elements associated with each atom */
    PerAtomQuantity<int> * atomElement_;

    /** pointer to the FE engine */
    const FE_Engine * feEngine_;

    /** tolerance used in the linear solver */
    const double tol_;

    // workspace variables
    mutable DIAG_MAT _nodeVolumesMatrix_;
    mutable DENS_VEC _nodeVolumes_;
    mutable DENS_VEC _nodeElementVolumes_;
    mutable DENS_MAT _nodalCoords_;
    mutable DENS_VEC _scale_;

  private:

    // do not define
    NodalAtomVolumeElement();

  };

  /**
   *  @class  AtomTypeElement
   *  @brief  determines if the given atom type is in each element
   */
  class AtomTypeElement : public DenseMatrixTransfer<int> {

  public:
    
    // constructor
    AtomTypeElement(ATC_Coupling * atc,
                    PerAtomQuantity<int> * atomElement = NULL);
    
    // destructor
    virtual ~AtomTypeElement() {
      atomElement_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** map from atom to the element in which it resides */
    PerAtomQuantity<int> * atomElement_;

    /** number of element */
    int nElts_;

    // workspace
    mutable INT_ARRAY _quantityLocal_;

  private:

    // do not define
    AtomTypeElement();

  };

  /**
   *  @class  ElementMask
   *  @brief  determines which elements should be used for FE quadrature
   */
  class ElementMask : public DenseMatrixTransfer<bool> {

  public:
    
    // constructor
    ElementMask(ATC_Coupling * atc,
                MatrixDependencyManager<DenseMatrix, int> * hasInternal = NULL,
                MatrixDependencyManager<DenseMatrix, int> * hasGhost = NULL);
    
    // destructor
    virtual ~ElementMask() {
      hasInternal_->remove_dependence(this);
      if (hasGhost_) hasGhost_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** shape function matrix over the internal atoms used to help define the geometry */
    MatrixDependencyManager<DenseMatrix, int> * hasInternal_;

    /** shape function matrix over the ghost atoms used to help define the geometry */
    MatrixDependencyManager<DenseMatrix, int> * hasGhost_;

    /** finite element engine */
    const FE_Engine * feEngine_;

  private:

    // do not define
    ElementMask();

  };

  /**
   *  @class  AtomElementMask
   *  @brief  determines which elements should be used for FE quadrature based on presence of atoms
   */
  class AtomElementMask : public DenseMatrixTransfer<bool> {

  public:
    
    // constructor
    AtomElementMask(ATC_Coupling * atc,
                    MatrixDependencyManager<DenseMatrix, int> * hasAtoms = NULL);
    
    // destructor
    virtual ~AtomElementMask() {
      hasAtoms_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** shape function matrix over the atoms used to help define the geometry */
    MatrixDependencyManager<DenseMatrix, int> * hasAtoms_;

    /** finite element engine */
    const FE_Engine * feEngine_;

  private:

    // do not define
    AtomElementMask();

  };

  /**
   *  @class  NodalGeometryType
   *  @brief  Computes the computational geometry associated with each node with respect to the FE and MD regions
   */
  
  class NodalGeometryType : public DenseMatrixTransfer<int> {

  public:
    
    // constructor
    NodalGeometryType(ATC_Coupling * atc,
                      MatrixDependencyManager<DenseMatrix, int> * hasInternal = NULL,
                      MatrixDependencyManager<DenseMatrix, int> * hasGhost = NULL);
    
    // destructor
    virtual ~NodalGeometryType() {
      hasInternal_->remove_dependence(this);
      if (hasGhost_) hasGhost_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** shape function matrix over the internal atoms used to help define the geometry */
    MatrixDependencyManager<DenseMatrix, int> * hasInternal_;

    /** shape function matrix over the ghost atoms used to help define the geometry */
    MatrixDependencyManager<DenseMatrix, int> * hasGhost_;

    /** finite element engine */
    const FE_Engine * feEngine_;

    /** number of nodes in the mesh */
    int nNodes_;

    /** number of elements in the mesh */
    int nElts_;

    // workspace
    // ints so they can be shared through MPI
    mutable Array<int> _nodesInternal_, _nodesGhost_;

  private:

    // do not define
    NodalGeometryType();

  };

  /**
   *  @class  NodalGeometryTypeElementSet
   *  @brief  Divdes nodes into MD, FE, and boundary based on if they are connected to nodes in internal and not internal
   */
  
  class NodalGeometryTypeElementSet : public DenseMatrixTransfer<int> {

  public:
    
    // constructor
    NodalGeometryTypeElementSet(ATC_Coupling * atc,
                                MatrixDependencyManager<DenseMatrix, int> * hasInternal = NULL);
    
    // destructor
    virtual ~NodalGeometryTypeElementSet() {
      hasInternal_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** shape function matrix over the internal atoms used to help define the geometry */
    MatrixDependencyManager<DenseMatrix, int> * hasInternal_;

    /** finite element engine */
    const FE_Engine * feEngine_;

    /** number of nodes in the mesh */
    int nNodes_;

    /** number of elements in the mesh */
    int nElts_;

    // workspace
    // ints so they can be shared through MPI
    mutable Array<int> _nodesInternal_, _nodesGhost_;

  private:

    // do not define
    NodalGeometryTypeElementSet();

  };

  /**
   *  @class  LargeToSmallMap
   *  @brief  mapping from a larger set to a smaller set
   */

  class LargeToSmallMap : public DenseMatrixTransfer<int> {

  public:
    
    // constructor
    LargeToSmallMap() : size_(0) {};
    
    // destructor
    virtual ~LargeToSmallMap() {};

    /** get the number of elements associated with the map */
    virtual int size() const {this->quantity(); return size_;};

  protected:

    /** number of nodes in the map */
    mutable int size_;

  };

  /**
   *  @class  NodeToSubset
   *  @brief  mapping from all nodes to a subset 
   */

  class NodeToSubset : public LargeToSmallMap {

  public:
    
    // constructor
    NodeToSubset(ATC_Method * atc,
                 SetDependencyManager<int> * subsetNodes);
    
    // destructor
    virtual ~NodeToSubset() {
      subsetNodes_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:
    
    /** pointer to atc to get the number of nodes */
    const ATC_Method * atc_;

    /** map from atom to the element in which it resides */
    SetDependencyManager<int> * subsetNodes_;

  private:

    // do not define
    NodeToSubset();

  };

  /**
   *  @class  SubsetToNode
   *  @brief  mapping from a subset of nodes to all nodes 
   */

  class SubsetToNode : public DenseMatrixTransfer<int> {

  public:
    
    // constructor
    SubsetToNode(NodeToSubset * nodeToSubset);
    
    // destructor
    virtual ~SubsetToNode() {
      nodeToSubset_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** map from atom to the element in which it resides */
    NodeToSubset * nodeToSubset_;

  private:

    // do not define
    SubsetToNode();

  };

  /**
   *  @class  ReducedSparseMatrix
   *  @brief  reduction of a sparse matrix to only those columns consistent with the map
   */
  
  class ReducedSparseMatrix : public SparseMatrixTransfer<double> {

  public:
    
    // constructor
    ReducedSparseMatrix(ATC_Method * atc,
                          SPAR_MAN * source,
                          LargeToSmallAtomMap * map);
    
    // destructor
    virtual ~ReducedSparseMatrix();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** original quantity */
    SPAR_MAN * source_;

    /** mapping transfer */
    LargeToSmallAtomMap * map_;

    // workspace
    mutable DenseVector<double> _row_; // vector storing the row values of a matrix
    mutable DenseVector<INDEX> _index_; // vector storing the column indices of the values

  private:

    // do not define
    ReducedSparseMatrix();

  };

  /**
   *  @class  RowMappedSparseMatrix
   *  @brief  mapping of rows from a sparse matrix to a new sparse matrix
   */
  
  class RowMappedSparseMatrix : public ReducedSparseMatrix {

  public:
    
    // constructor
    RowMappedSparseMatrix(ATC_Method * atc,
                          SPAR_MAN * source,
                          LargeToSmallAtomMap * map) :
               ReducedSparseMatrix(atc,source,map) {};
    
    // destructor
    virtual ~RowMappedSparseMatrix() {};

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

  private:

    // do not define
    RowMappedSparseMatrix();

  };

  /**
   *  @class  RowMappedSparseMatrixVector
   *  @brief  mapping of rows from a vector sparse matrices to a new vector of sparse matrices
   */
  
  class RowMappedSparseMatrixVector : public VectorTransfer<SPAR_MAT * > {

  public:
    
    // constructor
    RowMappedSparseMatrixVector(ATC_Method * atc,
                                VectorDependencyManager<SPAR_MAT * > * source,
                                LargeToSmallAtomMap * map);
    
    // destructor
    virtual ~RowMappedSparseMatrixVector();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** original quantity */
    VectorDependencyManager<SPAR_MAT * > * source_;

    /** mapping transfer */
    LargeToSmallAtomMap * map_;

    // workspace
    mutable DenseVector<double> _row_; // vector storing the row values of a matrix
    mutable DenseVector<INDEX> _index_; // vector storing the column indices of the values

  private:

    // do not define
    RowMappedSparseMatrixVector();

  };

  /**
   *  @class  MappedDiagonalMatrix
   *  @brief  mapping of a diagonal matrix to a new diagronal matrix
   */
  
  class MappedDiagonalMatrix : public DiagonalMatrixTransfer<double> {

  public:
    
    // constructor
    MappedDiagonalMatrix(ATC_Method * atc,
                         DIAG_MAN * source,
                         LargeToSmallAtomMap * map);
    
    // destructor
    virtual ~MappedDiagonalMatrix();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** original quantity */
    DIAG_MAN * source_;

    /** mapping transfer */
    LargeToSmallAtomMap * map_;

  private:

    // do not define
    MappedDiagonalMatrix();

  };

  /**
   *  @class  MappedQuantity
   *  @brief  generic reduced mapping
   */
  
  class MappedQuantity : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    MappedQuantity(ATC_Method * atc,
                   DENS_MAN * source,
                   LargeToSmallMap * map);
    
    // destructor
    virtual ~MappedQuantity() {
      source_->remove_dependence(this);
      map_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** original quantity */
    DENS_MAN * source_;

    /** mapping transfer */
    LargeToSmallMap * map_;

  private:

    // do not define
    MappedQuantity();

  };

  /**
   *  @class  RegulatedNodes
   *  @brief  set of all nodes being controlled by an atomic regulator
   */

  class RegulatedNodes : public SetTransfer<int> {

  public:
    
    // constructor
    RegulatedNodes(ATC_Coupling * atc,
                   FieldName fieldName = NUM_TOTAL_FIELDS,
                   MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType = NULL);
    
    // destructor
    virtual ~RegulatedNodes() {
      nodalGeometryType_->remove_dependence(this);
    };

  protected:

    /** recomputes the set of regulated nodes based on prescribed data and atom locations */
    virtual void reset_quantity() const;

    /** map from fields to their sizes */
    std::map<FieldName,int> fieldSizes_;

    /** map from atom to element in which it resides */
    MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType_;

    /** pointer to atc object */
    ATC_Coupling * atc_;

    /** pointer to the finite element engine */
    const FE_Engine * feEngine_;

    /** pointer to the prescribed data manager */
    const PrescribedDataManager * prescribedDataManager_;

    // methods
    /** adds boundary nodes */
    void insert_boundary_nodes() const;

    /** adds nodes associated with boundary facesets */
    void insert_boundary_faces() const;

    /** adds fixed nodes */
    void insert_fixed_nodes() const;

    /** adds nodes associated with face fluxes */
    void insert_face_fluxes() const;

    /** adds nodes associated with element fluxes */
    void insert_element_fluxes() const;

  private:

    // do not define
    RegulatedNodes();

  };

  /**
   *  @class  FluxNodes
   *  @brief  set of all nodes being controlled by an atomic regulator for fluxes
   */

  class FluxNodes : public RegulatedNodes {

  public:
    
    // constructor
    FluxNodes(ATC_Coupling * atc,
              FieldName fieldName = NUM_TOTAL_FIELDS,
              MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType = NULL) :
                RegulatedNodes(atc,fieldName,nodalGeometryType) {};
    
    // destructor
    virtual ~FluxNodes() {};

  protected:

    /** recomputes the set of regulated nodes based on prescribed data and atom locations */
    virtual void reset_quantity() const;

  private:

    // do not define
    FluxNodes();

  };

  /**
   *  @class  BoundaryNodes
   *  @brief  set of all nodes associated with the FE/MD boundary
   */

  class BoundaryNodes : public RegulatedNodes {

  public:
    
    // constructor
    BoundaryNodes(ATC_Coupling * atc,
                  FieldName fieldName = NUM_TOTAL_FIELDS,
                  MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType = NULL) :
                    RegulatedNodes(atc,fieldName,nodalGeometryType) {};
    
    // destructor
    virtual ~BoundaryNodes() {};

  protected:

    /** recomputes the set of regulated nodes based on prescribed data and atom locations */
    virtual void reset_quantity() const;

  private:

    // do not define
    BoundaryNodes();

  };

  /**
   *  @class  FluxBoundaryNodes
   *  @brief  set of all nodes being controlled by an atomic regulator for fluxes
   */

  class FluxBoundaryNodes : public FluxNodes {

  public:
    
    // constructor
    FluxBoundaryNodes(ATC_Coupling * atc,
                      FieldName fieldName = NUM_TOTAL_FIELDS,
                      MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType = NULL) :
                       FluxNodes(atc,fieldName,nodalGeometryType) {};
    
    // destructor
    virtual ~FluxBoundaryNodes() {};

  protected:

    /** recomputes the set of regulated nodes based on prescribed data and atom locations */
    virtual void reset_quantity() const;

  private:

    // do not define
    FluxBoundaryNodes();

  };

    /**
   *  @class  AllRegulatedNodes
   *  @brief  set of all nodes being controlled by an atomic regulator when localization is used
   */

  class AllRegulatedNodes : public FluxBoundaryNodes {

  public:
    
    // constructor
    AllRegulatedNodes(ATC_Coupling * atc,
                      FieldName fieldName = NUM_TOTAL_FIELDS,
                      MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType = NULL) :
                        FluxBoundaryNodes(atc,fieldName,nodalGeometryType) {};
    
    // destructor
    virtual ~AllRegulatedNodes() {};

  protected:

    /** recomputes the set of regulated nodes based on prescribed data and atom locations */
    virtual void reset_quantity() const;

  private:

    // do not define
    AllRegulatedNodes();

  };

    /**
   *  @class  FixedNodes
   *  @brief  set of all nodes being controlled by an atomic regulator which are fixed
   */

  class FixedNodes : public RegulatedNodes {

  public:
    
    // constructor
    FixedNodes(ATC_Coupling * atc,
               FieldName fieldName = NUM_TOTAL_FIELDS,
               MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType = NULL) :
                 RegulatedNodes(atc,fieldName,nodalGeometryType) {};
    
    // destructor
    virtual ~FixedNodes() {};

  protected:

    /** recomputes the set of regulated nodes based on prescribed data and atom locations */
    virtual void reset_quantity() const;

  private:

    // do not define
    FixedNodes();

  };

  /**
   *  @class  FixedBoundaryNodes
   *  @brief  set of all nodes being controlled by an atomic regulator which are fixed, including for coupling
   */

  class FixedBoundaryNodes : public FixedNodes {

  public:
    
    // constructor
    FixedBoundaryNodes(ATC_Coupling * atc,
                       FieldName fieldName = NUM_TOTAL_FIELDS,
                       MatrixDependencyManager<DenseMatrix, int> * nodalGeometryType = NULL) :
                         FixedNodes(atc,fieldName,nodalGeometryType) {};
    
    // destructor
    virtual ~FixedBoundaryNodes() {};

  protected:

    /** recomputes the set of regulated nodes based on prescribed data and atom locations */
    virtual void reset_quantity() const;

  private:

    // do not define
    FixedBoundaryNodes();

  };

  /**
   *  @class  ElementMaskNodeSet
   *  @brief  determines which elements should be used for FE quadrature based on them possessing certain nodes
   */
  class ElementMaskNodeSet : public DenseMatrixTransfer<bool> {

  public:
    
    // constructor
    ElementMaskNodeSet(ATC_Coupling * atc,
                       SetDependencyManager<int> * nodeSet);

    // destructor
    virtual ~ElementMaskNodeSet() {
      nodeSet_->remove_dependence(this);
    };

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** transfer determining used nodes */
    SetDependencyManager<int> * nodeSet_;

    /** finite element mesh */
    const FE_Mesh * feMesh_;

  private:

    // do not define
    ElementMaskNodeSet();

  };

  /**
   *  @class  PointToElementMap
   *  @brief  Class for identifying the element associated with a given point
   */

  class PointToElementMap : public DenseMatrixTransfer<int> {

  public:

    // constructor
    PointToElementMap(ATC_Method * atc,
                      MatrixDependencyManager<DenseMatrix, double> * pointPositions);

    // destructor
    virtual ~PointToElementMap();

  protected:

    /** resets the data if necessary */
    virtual void reset_quantity() const;

    /** point positions */
    MatrixDependencyManager<DenseMatrix, double> * pointPositions_;

    /** pointer to finite element mesh */
    const FE_Mesh * feMesh_;

  private:

    // do not define
    PointToElementMap();

  };

  /**
   *  @class  Interpolant
   *  @brief  constructs the spatial values of the shape functions
   */
  
  class Interpolant : public SparseMatrixTransfer<double> {

  public:

    // constructor
    Interpolant(ATC_Method * atc,
                MatrixDependencyManager<DenseMatrix, int>* pointToElementMap,
                DENS_MAN* pointPositions); 
    
    // destructor
    virtual ~Interpolant() {
      pointToElementMap_->remove_dependence(this);
      pointPositions_->remove_dependence(this);
    };

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const;

  protected:

    /** map from point coordinates to element */
    MatrixDependencyManager<DenseMatrix, int>* pointToElementMap_;

    /** coordinates used to evaluate shape functions */
    DENS_MAN* pointPositions_;

    /** pointer to the FE engine */
    const FE_Engine * feEngine_;

  private:

    // do not define
    Interpolant();

  };

  /**
   *  @class  InterpolantGradient
   *  @brief  constructs the spatial derivatives of the shape functions
   */
  
  class InterpolantGradient : public VectorTransfer<SPAR_MAT * > {

  public:

    // constructor
    InterpolantGradient(ATC_Method * atc,
                        MatrixDependencyManager<DenseMatrix, int>* pointToElementMap,
                        DENS_MAN* pointPositions); 
    
    // destructor
    virtual ~InterpolantGradient();

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const;

  protected:

    /** map from point coordinates to element */
    MatrixDependencyManager<DenseMatrix, int>* pointToElementMap_;

    /** coordinates used to evaluate shape functions */
    DENS_MAN* pointPositions_;

    /** pointer to the FE engine */
    const FE_Engine * feEngine_;

  private:

    // do not define
    InterpolantGradient();

  };

  /**
   *  @class  PerAtomShapeFunctionGradient
   *  @brief  constructs the spatial derivatives of the shape functions at each atom
   */
  
  class PerAtomShapeFunctionGradient : public VectorTransfer<SPAR_MAT * > {

  public:

    // constructor
    PerAtomShapeFunctionGradient(ATC_Method * atc,
                                 MatrixDependencyManager<DenseMatrix, int>* atomToElementMap = NULL,
                                 DENS_MAN* atomPositions = NULL,
                                 const std::string & tag = "AtomicShapeFunctionGradient",
                                 AtomType atomType = INTERNAL); 
    
    // destructor
    virtual ~PerAtomShapeFunctionGradient();

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const;

  protected:

    /** map from point coordinates to element */
    MatrixDependencyManager<DenseMatrix, int>* atomToElementMap_;

    /** coordinates used to evaluate shape functions */
    DENS_MAN* atomPositions_;

    /** pointer to the FE engine */
    const FE_Engine * feEngine_;

    /** container for dependency managers */
    std::vector<AtcAtomSparseMatrix<double> * > matrices_;

  private:

    // do not define
    PerAtomShapeFunctionGradient();

  };

  /**
   *  @class  InterpolantSmallMolecule
   *  @brief  constructs the spatial values of the shape functions for small molecules
   */
  
  class InterpolantSmallMolecule : public Interpolant {

  public:

    // constructor
    InterpolantSmallMolecule(ATC_Method * atc,
                             MatrixDependencyManager<DenseMatrix, int>* moleculeToElementMap,
                             DENS_MAN* moleculePositions,
                             MoleculeSet * moleculeSet); 
    
    // destructor
    virtual ~InterpolantSmallMolecule();

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const;

  protected:

    /** pointer to atc base class */
    MoleculeSet * moleculeSet_;

    // workspace
    mutable DENS_VEC _fractions_;

  private:

    // do not define
    InterpolantSmallMolecule();

  };



  /**
   *  @class  DenseMatrixQuotient
   *  @brief  quotient of two matrices
   */

  class DenseMatrixQuotient : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    DenseMatrixQuotient(DENS_MAN * matrixNumerator, DENS_MAN * matrixDenominator);
    
    // destructor
    virtual ~DenseMatrixQuotient() {
      matrixNumerator_->remove_dependence(this);
      matrixDenominator_->remove_dependence(this);
    };

    /** apply transfer operator  */
    virtual void reset_quantity() const;

    /** get number of columns */
    virtual int nCols() const {return matrixNumerator_->nCols();};

  protected:
  
    DENS_MAN* matrixNumerator_; 
    DENS_MAN* matrixDenominator_;

  private:

    // do not define
    DenseMatrixQuotient();

  };

  /**
   *  @class  DenseMatrixSum
   *  @brief  sum of two matrices
   */

  class DenseMatrixSum : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    DenseMatrixSum(DENS_MAN * matrixOne, DENS_MAN * matrixTwo);
    
    // destructor
    virtual ~DenseMatrixSum() {
      matrixOne_->remove_dependence(this);
      matrixTwo_->remove_dependence(this);
    };

    /** apply transfer operator  */
    virtual void reset_quantity() const;

    /** get number of columns */
    virtual int nCols() const {return matrixOne_->nCols();};

  protected:
  
    DENS_MAN* matrixOne_; 
    DENS_MAN* matrixTwo_;

  private:

    // do not define
    DenseMatrixSum();

  };

  /**
   *  @class  DenseMatrixDelta
   *  @brief change relative to a reference value
   */

  class DenseMatrixDelta : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    DenseMatrixDelta(DENS_MAN * current, DENS_MAT * reference);
    
    // destructor
    virtual ~DenseMatrixDelta() {
      matrix_->remove_dependence(this);
    };

    /** apply transfer operator  */
    virtual void reset_quantity() const;

  protected:
  
    DENS_MAN* matrix_; 
    DENS_MAT* reference_;

  private:

    // do not define
    DenseMatrixDelta();

  };

  /**
   *  @class  OnTheFlyKernelAccumulation
   *  @brief  implements the accumulant on the fly
   */
  
  class OnTheFlyKernelAccumulation : public DenseMatrixTransfer<double> {

  public:

    // constructor
    OnTheFlyKernelAccumulation(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               KernelFunction* kernelFunction,
                               DENS_MAN* atomCoarseGrainingPositions); 
    
    // destructor
    virtual ~OnTheFlyKernelAccumulation() {
      source_->remove_dependence(this);
      atomCoarseGrainingPositions_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

  protected:

    /** pointer to atc base class */
    const ATC_Method * atc_;

    /** pointer to source data */
    PerAtomQuantity<double> * source_;

    /** map from atom coordinates to element */
    KernelFunction* kernelFunction_;

    /** atomic coordinates used for coarse graining operations */
    DENS_MAN* atomCoarseGrainingPositions_;

    /** access to mesh */
    const FE_Mesh * feMesh_;

    // workspace
    mutable DENS_MAT _quantityLocal_;

  private:

    // do not define
    OnTheFlyKernelAccumulation();

  };

  /**
   *  @class  OnTheFlyKernelAccumulationNormalized
   *  @brief  implements a normalized accumulant on the fly
   */
  
  class OnTheFlyKernelAccumulationNormalized : public OnTheFlyKernelAccumulation {

  public:

    // constructor
    OnTheFlyKernelAccumulationNormalized(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               KernelFunction* kernelFunction,
                               DENS_MAN* atomCoarseGrainingPositions,
                               DIAG_MAN* weights);
    
    // destructor
    virtual ~OnTheFlyKernelAccumulationNormalized() {
      weights_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

  protected:

    DIAG_MAN * weights_;

  private:

    // do not define
    OnTheFlyKernelAccumulationNormalized();

  };

  /**
   *  @class  OnTheFlyKernelAccumulationNormalizedReferenced
   *  @brief  implements a normalized referenced accumulant on the fly
   */
  
  class OnTheFlyKernelAccumulationNormalizedReferenced : public OnTheFlyKernelAccumulationNormalized {

  public:

    // constructor
    OnTheFlyKernelAccumulationNormalizedReferenced(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               KernelFunction* kernelFunction,
                               DENS_MAN* atomCoarseGrainingPositions,
                               DIAG_MAN* weights,
                               DENS_MAN * reference);
    
    // destructor
    virtual ~OnTheFlyKernelAccumulationNormalizedReferenced() {
      reference_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

  protected:

    /** reference value */
    DENS_MAN * reference_;

  private:

    // do not define
    OnTheFlyKernelAccumulationNormalizedReferenced();

  };

  /**
   *  @class  OnTheFlyKernelNormalizedAccumulationScaled
   *  @brief  implements a scaled accumulant on the fly
   */
  
  class OnTheFlyKernelAccumulationNormalizedScaled : public OnTheFlyKernelAccumulationNormalized {

  public:

    // constructor
    OnTheFlyKernelAccumulationNormalizedScaled(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               KernelFunction* kernelFunction,
                               DENS_MAN* atomCoarseGrainingPositions,
                               DIAG_MAN* weights,
                               const double scale);
    
    // destructor
    virtual ~OnTheFlyKernelAccumulationNormalizedScaled() {
      atomCoarseGrainingPositions_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

  protected:

    const double scale_;

  private:

    // do not define
    OnTheFlyKernelAccumulationNormalizedScaled();

  };

  /**
   *  @class  AccumulantWeights
   *  @brief  weights for kernel function accumulants
   */

  class AccumulantWeights : public DiagonalMatrixTransfer<double> {

  public:
    
    // constructor
    AccumulantWeights(SPAR_MAN * accumulant);
    
    // destructor
    virtual ~AccumulantWeights() {
      accumulant_->remove_dependence(this);
    };

    /** apply transfer operator  */
    virtual void reset_quantity() const;

  protected:
  
    SPAR_MAN* accumulant_;

    // workspace
    mutable DENS_VEC _localWeights_;
    mutable DENS_VEC _weights_;

  private:

    // do not define
    AccumulantWeights();

  };

  /**
   *  @class  OnTheFlyKernelWeights
   *  @brief  weights for on-the-fly kernel function
   */

  class OnTheFlyKernelWeights : public DiagonalMatrixTransfer<double> {

  public:
    
    // constructor
    OnTheFlyKernelWeights(DENS_MAN * weights);
    
    // destructor
    virtual ~OnTheFlyKernelWeights() {
      weights_->remove_dependence(this);
    };

    /** apply transfer operator  */
    virtual void reset_quantity() const;

  protected:
  
    DENS_MAN* weights_;

  private:

    // do not define
    OnTheFlyKernelWeights();

  };

  /**
   *  @class  KernelInverseVolumes
   *  @brief  inverse volumes for kernel function accumulants
   */

  class KernelInverseVolumes : public DiagonalMatrixTransfer<double> {

  public:
    
    // constructor
    KernelInverseVolumes(ATC_Method * atc,
                         KernelFunction* kernelFunction);
    
    // destructor
    virtual ~KernelInverseVolumes() {};

    /** apply transfer operator  */
    virtual void reset_quantity() const;

  protected:
  
    /** kernel function being used */
    KernelFunction* kernelFunction_;

    /** access to mesh */
    const FE_Mesh * feMesh_;

  private:

    // do not define
    KernelInverseVolumes();

  };

  /**
   *  @class  OnTheFlyMeshAccumulation
   *  @brief  implements the mesh-based accumulant on the fly
   */
  
  class OnTheFlyMeshAccumulation : public DenseMatrixTransfer<double> {

  public:

    // constructor
    OnTheFlyMeshAccumulation(ATC_Method * atc,
                             PerAtomQuantity<double> * source,
                             DENS_MAN* atomCoarseGrainingPositions); 
    
    // destructor
    virtual ~OnTheFlyMeshAccumulation() {
      source_->remove_dependence(this);
      atomCoarseGrainingPositions_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

    /** access nCols_ */
    virtual int nCols() const { return source_->nCols(); } 

  protected:

    /** pointer to atc base class */
    const ATC_Method * atc_;

    /** pointer to source data */
    PerAtomQuantity<double> * source_;

    /** atomic coordinates used for coarse graining operations */
    DENS_MAN* atomCoarseGrainingPositions_;

    /** access to mesh */
    const FE_Mesh * feMesh_;

    // workspace
    mutable DENS_MAT _quantityLocal_;

  private:

    // do not define
    OnTheFlyMeshAccumulation();

  };

  /**
   *  @class  OnTheFlyMeshAccumulationNormalized
   *  @brief  implements a normalized mesh-based accumulant on the fly
   */
  
  class OnTheFlyMeshAccumulationNormalized : public OnTheFlyMeshAccumulation {

  public:

    // constructor
    OnTheFlyMeshAccumulationNormalized(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               DENS_MAN* atomCoarseGrainingPositions,
                               DIAG_MAN* weights);
    
    // destructor
    virtual ~OnTheFlyMeshAccumulationNormalized() {
      weights_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

  protected:

    DIAG_MAN * weights_;

  private:

    // do not define
    OnTheFlyMeshAccumulationNormalized();

  };

  /**
   *  @class  OnTheFlyMeshAccumulationNormalizedReferenced
   *  @brief  implements a normalized referenced mesh-based accumulant on the fly
   */
  
  class OnTheFlyMeshAccumulationNormalizedReferenced : public OnTheFlyMeshAccumulationNormalized {

  public:

    // constructor
    OnTheFlyMeshAccumulationNormalizedReferenced(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               DENS_MAN* atomCoarseGrainingPositions,
                               DIAG_MAN* weights,
                               DENS_MAN * reference);
    
    // destructor
    virtual ~OnTheFlyMeshAccumulationNormalizedReferenced() {
      reference_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

  protected:

    /** reference value */
    DENS_MAN * reference_;

  private:

    // do not define
    OnTheFlyMeshAccumulationNormalizedReferenced();

  };

  /**
   *  @class  OnTheFlyMeshAccumulationNormalizedScaled
   *  @brief  implements a scaled mesh-based accumulant on the fly
   */
  
  class OnTheFlyMeshAccumulationNormalizedScaled : public OnTheFlyMeshAccumulationNormalized {

  public:

    // constructor
    OnTheFlyMeshAccumulationNormalizedScaled(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               DENS_MAN* atomCoarseGrainingPositions,
                               DIAG_MAN* weights,
                               const double scale);
    
    // destructor
    virtual ~OnTheFlyMeshAccumulationNormalizedScaled() {
      atomCoarseGrainingPositions_->remove_dependence(this);
    };

    /** do nothing */
    virtual void reset_quantity() const;

  protected:

    const double scale_;

  private:

    // do not define
    OnTheFlyMeshAccumulationNormalizedScaled();

  };

  /**
   *  @class  NativeShapeFunctionGradient
   *  @brief  constructs the spatial derivatives of the shape functions
   */
  
  class NativeShapeFunctionGradient : public VectorTransfer<SPAR_MAT * > {

  public:

    // constructor
    NativeShapeFunctionGradient(ATC_Method * atc);
    
    // destructor
    virtual ~NativeShapeFunctionGradient();

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const;

  protected:

    /** pointer to the FE engine */
    const FE_Engine * feEngine_;

  private:

    // do not define
    NativeShapeFunctionGradient();

  };

  /**
   *  @class  OnTheFlyShapeFunctionProlongation
   *  @brief  implements the interpolant on the fly
   */
  
  class OnTheFlyShapeFunctionProlongation : public FeToAtomTransfer {

  public:

    // constructor
    OnTheFlyShapeFunctionProlongation(ATC_Method * atc,
                                      DENS_MAN * source,
                                      DENS_MAN * atomCoarseGrainingPositions); 
    
    // destructor
    virtual ~OnTheFlyShapeFunctionProlongation();

    /** do nothing */
    virtual void reset() const;

  protected:

    /** pointer to atc base class */
    //const ATC_Method * atc_;

    /** pointer to source data */
    //DENS_MAN * source_;

    /** atomic coordinates used for coarse graining operations */
    DENS_MAN* atomCoarseGrainingPositions_;

    /** access to mesh */
    const FE_Mesh * feMesh_;

    // workspace
    //mutable DENS_MAT _quantityLocal_;

  private:

    // do not define
    OnTheFlyShapeFunctionProlongation();

  };

}
#endif

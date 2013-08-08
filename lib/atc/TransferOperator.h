// A class for defining transfer operations

#ifndef TRANSFER_OPERATOR_H
#define TRANSFER_OPERATOR_H

// ATC_Method headers
#include "PerAtomQuantityLibrary.h"
#include <set>
#include <vector>

using namespace std;

namespace ATC {

  // forward declarations
  class ATC_Method;
  class KernelFunction;
  class FE_Mesh;

  /**
   *  @class  DenseMatrixTransfer 
   *  @brief  Class for defining objects that generate dense matrix quantities from other matrix quantities
   */
  template <typename T>
  class DenseMatrixTransfer : public MatrixDependencyManager<DenseMatrix, T> {

  public:
    
    // constructor
    DenseMatrixTransfer() :
    MatrixDependencyManager<DenseMatrix, T>(),
      lammpsInterface_(LammpsInterface::instance()) {};
    
    // destructor
    virtual ~DenseMatrixTransfer() {};

    /** apply transfer operator */
    virtual const DenseMatrix<T> & quantity() const {
      if (this->need_reset()) {
        this->reset_quantity(); 
        MatrixDependencyManager<DenseMatrix, T>::needReset_ = false;
      } 
      return MatrixDependencyManager<DenseMatrix, T>::quantity_;
    };


        /** sets the quantity to a given value */
    virtual void operator=(const DenseMatrix<T> & target)
      {throw ATC_Error("DenseMatrixTransfer::set_quantity - Cannot modify transfer-based matrices");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("DenseMatrixTransfer::operator= - Cannot modify transfer-based matrices");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DenseMatrix<T> & addition)
      {throw ATC_Error("DenseMatrixTransfer::operator+= - Cannot modify transfer-based matrices");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("DenseMatrixTransfer::operator+= - Cannot modify transfer-based matrices");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DenseMatrix<T> & subtraction)
      {throw ATC_Error("DenseMatrixTransfer::operator-= - Cannot modify transfer-based matrices");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("DenseMatrixTransfer::operator-= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DenseMatrix<T> & multiplier)
      {throw ATC_Error("DenseMatrixTransfer::operator*= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("DenseMatrixTransfer::operator*= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DenseMatrix<T> & divisor)
      {throw ATC_Error("DenseMatrixTransfer::operator/= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("DenseMatrixTransfer::operator/= - Cannot modify transfer-based matrices");};

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const = 0;

    /** pointer to LammpsInterface for MPI */
    LammpsInterface * lammpsInterface_;

  };

  /**
   *  @class  SparseMatrixTransfer 
   *  @brief  Class for defining objects that generate dense matrix quantities from other matrix quantities
   */
  template <typename T>
  class SparseMatrixTransfer : public MatrixDependencyManager<SparseMatrix, T> {

  public:
    
    // constructor
    SparseMatrixTransfer() :
    MatrixDependencyManager<SparseMatrix, T>(),
      lammpsInterface_(LammpsInterface::instance()) {};
    
    // destructor
    virtual ~SparseMatrixTransfer() {};

    /** apply transfer operator */
    virtual const SparseMatrix<T> & quantity() const {if (this->need_reset()) {this->reset_quantity(); MatrixDependencyManager<SparseMatrix, T>::needReset_ = false;} return MatrixDependencyManager<SparseMatrix, T>::quantity_;};


        /** sets the quantity to a given value */
    virtual void operator=(const SparseMatrix<T> & target)
      {throw ATC_Error("SparseMatrixTransfer::set_quantity - Cannot modify transfer-based matrices");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("SparseMatrixTransfer::operator= - Cannot modify transfer-based matrices");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const SparseMatrix<T> & addition)
      {throw ATC_Error("SparseMatrixTransfer::operator+= - Cannot modify transfer-based matrices");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("SparseMatrixTransfer::operator+= - Cannot modify transfer-based matrices");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const SparseMatrix<T> & subtraction)
      {throw ATC_Error("SparseMatrixTransfer::operator-= - Cannot modify transfer-based matrices");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("SparseMatrixTransfer::operator-= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const SparseMatrix<T> & multiplier)
      {throw ATC_Error("SparseMatrixTransfer::operator*= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("SparseMatrixTransfer::operator*= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const SparseMatrix<T> & divisor)
      {throw ATC_Error("SparseMatrixTransfer::operator/= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("SparseMatrixTransfer::operator/= - Cannot modify transfer-based matrices");};

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const = 0;

    /** pointer to LammpsInterface for MPI */
    LammpsInterface * lammpsInterface_;

  };

  /**
   *  @class  DiagonalMatrixTransfer 
   *  @brief  Class for defining objects that generate diagonal matrix quantities from other matrix quantities
   */
  template <typename T>
  class DiagonalMatrixTransfer : public MatrixDependencyManager<DiagonalMatrix, T> {

  public:
    
    // constructor
    DiagonalMatrixTransfer() :
    MatrixDependencyManager<DiagonalMatrix, T>(),
      lammpsInterface_(LammpsInterface::instance()) {};
    
    // destructor
    virtual ~DiagonalMatrixTransfer() {};

    /** apply transfer operator */
    virtual const DiagonalMatrix<T> & quantity() const {if (this->need_reset()) {this->reset_quantity(); MatrixDependencyManager<DiagonalMatrix, T>::needReset_ = false;} return MatrixDependencyManager<DiagonalMatrix, T>::quantity_;};


        /** sets the quantity to a given value */
    virtual void operator=(const DiagonalMatrix<T> & target)
      {throw ATC_Error("DiagonalMatrixTransfer::set_quantity - Cannot modify transfer-based matrices");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("DiagonalMatrixTransfer::operator= - Cannot modify transfer-based matrices");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DiagonalMatrix<T> & addition)
      {throw ATC_Error("DiagonalMatrixTransfer::operator+= - Cannot modify transfer-based matrices");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("DiagonalMatrixTransfer::operator+= - Cannot modify transfer-based matrices");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DiagonalMatrix<T> & subtraction)
      {throw ATC_Error("DiagonalMatrixTransfer::operator-= - Cannot modify transfer-based matrices");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("DiagonalMatrixTransfer::operator-= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DiagonalMatrix<T> & multiplier)
      {throw ATC_Error("DiagonalMatrixTransfer::operator*= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("DiagonalMatrixTransfer::operator*= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DiagonalMatrix<T> & divisor)
      {throw ATC_Error("DiagonalMatrixTransfer::operator/= - Cannot modify transfer-based matrices");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("DiagonalMatrixTransfer::operator/= - Cannot modify transfer-based matrices");};

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const = 0;

    /** pointer to LammpsInterface for MPI */
    LammpsInterface * lammpsInterface_;

  };

  /**
   *  @class  SetTransfer 
   *  @brief  Class for defining objects that generate sets using prescribed algorithms
   */
  template <typename T>
  class SetTransfer : public SetDependencyManager<T> {

  public:
    
    // constructor
    SetTransfer() :
      SetDependencyManager<T>() {};
    
    // destructor
    virtual ~SetTransfer() {};

    /** apply transfer operator */
    virtual const set<T> & quantity() const {if (this->need_reset()) {this->reset_quantity(); SetDependencyManager<T>::needReset_ = false;} return SetDependencyManager<T>::quantity_;};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual set<T> & set_quantity()
      {throw ATC_Error("SetTransfer::set_quantity - Cannot modify protected quantities"); return this->quantity_;};

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const = 0;

  };

  /**
   *  @class  VectorTransfer 
   *  @brief  Class for defining objects that generate sets using prescribed algorithms
   */
  template <typename T>
  class VectorTransfer : public VectorDependencyManager<T> {

  public:
    
    // constructor
    VectorTransfer() :
      VectorDependencyManager<T>() {};
    
    // destructor
    virtual ~VectorTransfer() {};

    /** apply transfer operator */
    virtual const vector<T> & quantity() const {if (this->need_reset()) {this->reset_quantity(); VectorDependencyManager<T>::needReset_ = false;} return VectorDependencyManager<T>::quantity_;};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual vector<T> & set_quantity()
      {throw ATC_Error("VectorTransfer::set_quantity - Cannot modify protected quantities"); return this->quantity_;};

  protected:

    /** does the actual computation of the quantity */
    virtual void reset_quantity() const = 0;

  };

  /**
   *  @class  AtomToFeTransfer 
   *  @brief  Class for defining objects to transfer atomistic quantities to FE quantities 
   */

  class AtomToFeTransfer : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    AtomToFeTransfer(ATC_Method * atc,
                     PerAtomQuantity<double> * source);
    
    // destructor
    virtual ~AtomToFeTransfer();

  protected:

    /** pointer to atc */
    ATC_Method * atc_;

    /** pointer to source atomic quantity data */
    PerAtomQuantity<double> * source_;

  private:

    // do not define
    AtomToFeTransfer();

  };

  /**
   *  @class  AtomDiagonalMatrixToFeTransfer 
   *  @brief  Class for defining objects to transfer atomistic quantities to FE quantities 
   */

  class AtomDiagonalMatrixToFeTransfer : public DenseMatrixTransfer<double> {

  public:
    
    // constructor
    AtomDiagonalMatrixToFeTransfer(ATC_Method * atc,
                                   PerAtomDiagonalMatrix<double> * source);
    
    // destructor
    virtual ~AtomDiagonalMatrixToFeTransfer();

  protected:

    /** pointer to atc */
    ATC_Method * atc_;

    /** pointer to source atomic quantity data */
    PerAtomDiagonalMatrix<double> * source_;

  private:

    // do not define
    AtomDiagonalMatrixToFeTransfer();

  };

  /**
   *  @class  FeToAtomTransfer
   *  @brief  Class for defining objects to transfer FE quantities to atomistic quantities
   */

  class FeToAtomTransfer : public ProtectedAtomQuantity<double> {

  public:
    
    // constructor
    FeToAtomTransfer(ATC_Method * atc,
                     DENS_MAN * source);
    
    // destructor
    virtual ~FeToAtomTransfer();

  protected:

    /** pointer to source Fe matrix data */
    DENS_MAN * source_;

  private:

    // do not define
    FeToAtomTransfer();

  };

  /**
   *  @class  FeToAtomDiagonalMatrix
   *  @brief  Class for defining objects to transfer FE quantities to atomistic diagonal matrices
   */

  class FeToAtomDiagonalMatrix : public ProtectedAtomDiagonalMatrix<double> {

  public:
    
    // constructor
    FeToAtomDiagonalMatrix(ATC_Method * atc,
                           DENS_MAN * source);
    
    // destructor
    virtual ~FeToAtomDiagonalMatrix();

  protected:

    /** pointer to source Fe matrix data */
    DENS_MAN * source_;

  private:

    // do not define
    FeToAtomDiagonalMatrix();

  };

  /**
   *  @class  MatToMatTransfer
   *  @brief  Class for defining objects that transfer quantities between materials
   */
  template <typename T>
  class MatToMatTransfer : public DenseMatrixTransfer<T> {

  public:
    
    // constructor
    MatToMatTransfer(MatrixDependencyManager<DenseMatrix, T> * source) :
      source_(source) {source_->register_dependence(this);};
    
    // destructor
    virtual ~MatToMatTransfer() {source_->remove_dependence(this);};

  protected:

    /** pointer to source matrix data */
    MatrixDependencyManager<DenseMatrix, T> * source_;

  private:

    // do not define
    MatToMatTransfer();

  };

  /**
   *  @class  AtfShapeFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic quantities to FE using shape functions
   *          (implements restrict_volumetric_quantity)
   */
  
  class AtfShapeFunctionRestriction : public AtomToFeTransfer {

  public:
    
    // constructor
    AtfShapeFunctionRestriction(ATC_Method * atc,
                                PerAtomQuantity<double> * source,
                                SPAR_MAN * shapeFunction);
    
    // destructor
    virtual ~AtfShapeFunctionRestriction();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** reference to shape function matrix */
    SPAR_MAN * shapeFunction_;

    /** persistant workspace */
    
    
    mutable DENS_MAT _workspace_;

    /** applies restriction operation across all processors */
    virtual void global_restriction() const;

    /** applies restriction operation on this processor */
    virtual void local_restriction(const DENS_MAT & sourceMatrix,
                                   const SPAR_MAT & shapeFunctionMatrix) const;

  private:

    // do not define
    AtfShapeFunctionRestriction();

  };

  /**
   *  @class  AdmtfShapeFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic diagonal matrices to FE using shape functions=
   */

  class AdmtfShapeFunctionRestriction : public AtomDiagonalMatrixToFeTransfer {

  public:
    
    // constructor
    AdmtfShapeFunctionRestriction(ATC_Method * atc,
                                PerAtomDiagonalMatrix<double> * source,
                                SPAR_MAN * shapeFunction);
    
    // destructor
    virtual ~AdmtfShapeFunctionRestriction();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** reference to shape function matrix */
    SPAR_MAN * shapeFunction_;

    /** persistant workspace */
    
    
    mutable DENS_MAT _workspace_;

    /** applies restriction operation across all processors */
    virtual void global_restriction() const;

    /** applies restriction operation on this processor */
    virtual void local_restriction(const DENS_MAT & sourceMatrix,
                                   const SPAR_MAT & shapeFunctionMatrix) const;

  private:

    // do not define
    AdmtfShapeFunctionRestriction();

  };

  /**
   *  @class  AtfProjection
   *  @brief  
   */

  class AtfProjection : public AtomToFeTransfer {

  public:
    
    // constructor
    AtfProjection(ATC_Method * atc,
                  PerAtomQuantity<double> * source,
                  SPAR_MAN * accumulant,
                  DIAG_MAN * weights = NULL);
    
    // destructor
    virtual ~AtfProjection();

    /** apply transfer operator */
    virtual void reset_quantity() const;

    /** get number of columns */
    virtual int nCols() const {return source_->nCols();};

  protected:

    /** reference to shape function matrix */
    SPAR_MAN * accumulant_;
    DIAG_MAN * weights_;
    DENS_MAT * reference_;

    /** persistant workspace */
    
    
    mutable DENS_MAT _workspace_;

    /** applies restriction operation across all processors */
    virtual void global_restriction() const;

    /** applies restriction operation on this processor */
    virtual void local_restriction(const DENS_MAT & sourceMatrix,
                                   const SPAR_MAT & shapeFunctionMatrix) const;

  private:

    // do not define
    AtfProjection();

  };
  class AtfProjectionScaled : public AtfProjection {

  public:
    
    // constructor
    AtfProjectionScaled(ATC_Method * atc,
                  PerAtomQuantity<double> * source,
                  SPAR_MAN * accumulant,
                  const double scale,
                  DIAG_MAN * weights = NULL);
    
    // destructor
    virtual ~AtfProjectionScaled();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:
    /** reference to shape function matrix */
    double     scale_;

  private:

    // do not define
    AtfProjectionScaled();

  };

  /**
   *  @class  AtfProjectionReferenced
   *  @brief  
   */

  class AtfProjectionReferenced : public AtfProjection {

  public:
    
    // constructor
    AtfProjectionReferenced(ATC_Method * atc,
                  PerAtomQuantity<double> * source,
                  SPAR_MAN * accumulant,
                  const DENS_MAT * reference,
                  DIAG_MAN * weights = NULL);
    
    // destructor
    virtual ~AtfProjectionReferenced();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    const DENS_MAT * reference_;

  private:

    // do not define
    AtfProjectionReferenced();

  };


  /**
   *  @class  AtfWeightedShapeFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic quantities to FE using shape functions
   *          including approximate quadrature weights
   *          (implements restrict_unscaled)
   */

  class AtfWeightedShapeFunctionRestriction : public AtfShapeFunctionRestriction {

  public:
    
    // constructor
    AtfWeightedShapeFunctionRestriction(ATC_Method * atc,
                                        PerAtomQuantity<double> * source,
                                        SPAR_MAN * shapeFunction,
                                        DIAG_MAN * weights);
    
    // destructor
    virtual ~AtfWeightedShapeFunctionRestriction() {weights_->remove_dependence(this);};

  protected:

    /** reference to diagonal weighting matrix */
    DIAG_MAN * weights_;

    /** applies restriction operation on this processor */
    virtual void local_restriction(const DENS_MAT & sourceMatrix,
                                   const SPAR_MAT & shapeFunctionMatrix) const;

  private:

    // do not define
    AtfWeightedShapeFunctionRestriction();

  };

  /**
   *  @class  AtfNodeWeightedShapeFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic quantities to FE using shape functions
   *          including weighting at the mesh nodes
   */

  class AtfNodeWeightedShapeFunctionRestriction : public AtfShapeFunctionRestriction {

  public:
    
    // constructor
    AtfNodeWeightedShapeFunctionRestriction(ATC_Method * atc,
                                            PerAtomQuantity<double> * source,
                                            SPAR_MAN * shapeFunction,
                                            DIAG_MAN * weights);
    
    // destructor
    virtual ~AtfNodeWeightedShapeFunctionRestriction() {weights_->remove_dependence(this);};

  protected:

    /** reference to diagonal weighting matrix */
    DIAG_MAN * weights_;

    /** applies restriction operation across all processors */
    virtual void global_restriction() const;

  private:

    // do not define
    AtfNodeWeightedShapeFunctionRestriction();

  };

  /**
   *  @class  AtfShapeFunctionProjection
   *  @brief  Class for defining objects that transfer restricted atomistic quantities to FE using shape functions
   *          (implements project/project_volumetric_quantity assuming restrict_unscaled/
   *           restrict_volumetric_quantity has been applied)
   */

  class AtfShapeFunctionProjection : public MatToMatTransfer<double> {

  public:
    
    // constructor
    AtfShapeFunctionProjection(ATC_Method * atc, 
                               DENS_MAN * source,
                               FieldName thisField);
    
    // destructor
    virtual ~AtfShapeFunctionProjection();

    /** apply transfer operator */
    virtual void reset_quantity() const ;

  protected:

    /** pointer to atc object for mass matrix inversion */
    ATC_Method * atc_;

    /** field name to define mass matrix to use */
    FieldName thisField_;

  private:

    // do not define
    AtfShapeFunctionProjection();

  };

  /**
   *  @class  AtfShapeFunctionMdProjection
   *  @brief  Class for defining objects that transfer restricted atomistic quantities to FE using shape functions for the MD region
   *          (implements project_md/project_md_volumetric_quantity assuming 
   *           restrict_unscaled/restrict_volumetric_quantity has been applied)
   */

  class AtfShapeFunctionMdProjection : public MatToMatTransfer<double> {

  public:
    
    // constructor
    AtfShapeFunctionMdProjection(ATC_Method * atc, 
                                 DENS_MAN * source,
                                 FieldName thisField);
    
    // destructor
    virtual ~AtfShapeFunctionMdProjection();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** pointer to atc object for mass matrix inversion */
    ATC_Method * atc_;

    /** field name to define mass matrix to use */
    FieldName thisField_;

  private:

    // do not define
    AtfShapeFunctionMdProjection();

  };

  /**
   *  @class  AtfShapeFunctionMdProjectionScaled
   *  @brief  Class for defining objects that transfer restricted atomistic quantities to FE using shape functions for the MD region with a scaling factor
   *          (implements project_md/project_md_volumetric_quantity assuming 
   *           restrict_unscaled/restrict_volumetric_quantity has been applied)
   */

  class AtfShapeFunctionMdProjectionScaled : public AtfShapeFunctionMdProjection {

  public:
    
    // constructor
    AtfShapeFunctionMdProjectionScaled(ATC_Method * atc,
                                           DENS_MAN * source,
                                           double scale,
                                           FieldName thisField);
    
    // destructor
    virtual ~AtfShapeFunctionMdProjectionScaled();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** scale multiplier */
   double scale_;

  private:

    // do not define
    AtfShapeFunctionMdProjectionScaled();

  };

  /**
   *  @class  AtfShapeFunctionMdProjectionReferenced
   *  @brief  Class for defining objects that transfer restricted atomistic quantities to FE using shape functions for the MD region with respect to a reference
   *          (implements project_md/project_md_volumetric_quantity assuming 
   *           restrict_unscaled/restrict_volumetric_quantity has been applied)
   */

  class AtfShapeFunctionMdProjectionReferenced : public AtfShapeFunctionMdProjection {

  public:
    
    // constructor
    AtfShapeFunctionMdProjectionReferenced(ATC_Method * atc,
                                           DENS_MAN * source,
                                           const DENS_MAT * reference,
                                           FieldName thisField);
    
    // destructor
    virtual ~AtfShapeFunctionMdProjectionReferenced();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** reference value */
    const DENS_MAT * reference_;

  private:

    // do not define
    AtfShapeFunctionMdProjectionReferenced();

  };

  /**
   *  @class  AtfKernelFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic quantities to FE using kernel functions
   */

  class AtfKernelFunctionRestriction : public AtomToFeTransfer {

  public:
    
    // constructor
    AtfKernelFunctionRestriction(ATC_Method * atc,
                                 PerAtomQuantity<double> * source,
                                 PerAtomQuantity<double> * coarseGrainingPositions,
                                 KernelFunction * kernelFunction);
    
    // destructor
    virtual ~AtfKernelFunctionRestriction();

    /** apply transfer operator */
    virtual void reset_quantity() const;

  protected:

    /** pointer to the kernel function evaluator */
    KernelFunction * kernelFunction_;

    /** reference positions for coarse graining operations */
    PerAtomQuantity<double> * coarseGrainingPositions_;

    /** pointer to the mesh being used */
    const FE_Mesh * feMesh_;

    /** persistant workspace */
    
    
    mutable DENS_MAT _workspace_;
    mutable DENS_VEC _xI_, _xa_, _xaI_;

    /** applies restriction operation across all processors */
    virtual void global_restriction() const;

    /** applies restriction operation on this processor */
    virtual void local_restriction(const DENS_MAT & sourceMatrix,
                                   const DENS_MAT & positions,
                                   const KernelFunction * kernelFunction) const;

  private:

    // do not define
    AtfKernelFunctionRestriction();

  };

  /**
   *  @class  AtfWeightedKernelFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic quantities to mesh using kernel functions
   *          including approximate quadrature weights
   */

  class AtfWeightedKernelFunctionRestriction : public AtfKernelFunctionRestriction {

  public:
    
    // constructor
    AtfWeightedKernelFunctionRestriction(ATC_Method * atc,
                                         PerAtomQuantity<double> * source,
                                         PerAtomQuantity<double> * coarseGrainingPositions,
                                         KernelFunction * kernelFunction,
                                         DIAG_MAN * weights);
    
    // destructor
    virtual ~AtfWeightedKernelFunctionRestriction() {weights_->remove_dependence(this);};

  protected:

    /** reference to diagonal weighting matrix */
    DIAG_MAN * weights_;

    /** applies restriction operation on this processor */
    virtual void local_restriction(const DENS_MAT & sourceMatrix,
                                   const DENS_MAT & positions,
                                   const KernelFunction * kernelFunction) const;

  private:

    // do not define
    AtfWeightedKernelFunctionRestriction();

  };

  /**
   *  @class  AtfNodeWeightedKernelFunctionRestriction
   *  @brief  Class for defining objects that transfer atomistic quantities to a mesh using kernel functions
   *          including weighting at the mesh nodes
   */

  class AtfNodeWeightedKernelFunctionRestriction : public AtfKernelFunctionRestriction {

  public:
    
    // constructor
    AtfNodeWeightedKernelFunctionRestriction(ATC_Method * atc,
                                             PerAtomQuantity<double> * source,
                                             PerAtomQuantity<double> * coarseGrainingPositions,
                                             KernelFunction * kernelFunction,
                                             DIAG_MAN * weights);
    
    // destructor
    virtual ~AtfNodeWeightedKernelFunctionRestriction() {weights_->remove_dependence(this);};

  protected:

    /** reference to diagonal weighting matrix */
    DIAG_MAN * weights_;

    /** applies restriction operation across all processors */
    virtual void global_restriction() const;

  private:

    // do not define
    AtfNodeWeightedKernelFunctionRestriction();

  };

  /**
   *  @class  FtaShapeFunctionProlongation
   *  @brief  Class for defining objects that transfer FE quantities to atoms using shape functions
   *          (implements prolong)
   */

  class FtaShapeFunctionProlongation : public FeToAtomTransfer {

  public:
    
    // constructor
    FtaShapeFunctionProlongation(ATC_Method * atc,
                                 DENS_MAN * source,
                                 SPAR_MAN * shapeFunction);
    
    // destructor
    virtual ~FtaShapeFunctionProlongation();

  protected:

    /** apply transfer operator if needed */
    virtual void reset() const;

    /** reference to shape function matrix */
    SPAR_MAN * shapeFunction_;

  private:

    // do not define
    FtaShapeFunctionProlongation();

  };

  /**
   *  @class  FtaShapeFunctionProlongationDiagonalMatrix
   *  @brief  Class for defining objects that transfer FE quantities to atomic diagonal matrices using shape functions
   *          (implements prolong)
   */

  class FtaShapeFunctionProlongationDiagonalMatrix : public FeToAtomDiagonalMatrix {

  public:
    
    // constructor
    FtaShapeFunctionProlongationDiagonalMatrix(ATC_Method * atc,
                                               DENS_MAN * source,
                                               SPAR_MAN * shapeFunction);
    
    // destructor
    virtual ~FtaShapeFunctionProlongationDiagonalMatrix();

  protected:

    /** apply transfer operator if needed */
    virtual void reset() const;

    /** reference to shape function matrix */
    SPAR_MAN * shapeFunction_;

    // workspace
    mutable DENS_MAT _temp_; // temporary storage for dense matrix

  private:

    // do not define
    FtaShapeFunctionProlongationDiagonalMatrix();

  };

  /**
   *   transfer for a matrix to a gradient mitigated by sparse matrices
   *
   */
  class MatToGradBySparse : public MatToMatTransfer<double> {

  public: 

    //constructor 
    MatToGradBySparse(ATC_Method * atc,
                      DENS_MAN * source,
                      VectorDependencyManager<SPAR_MAT * > * gradientMatrices);
    //destructor
    virtual ~MatToGradBySparse();

    // apply transfer operator  
    virtual void reset_quantity() const;

  protected:
    
    // pointer to sparseMatrix
    VectorDependencyManager<SPAR_MAT * > * gradientMatrices_;

  private: 

    // do not define
    MatToGradBySparse();
  
  };

  /**
   * transfer from dense to dense by diagonal matrix multiplier for anything
  **/ 
  class DiagonalMatrixMultiply : public MatToMatTransfer<double> {

  public: 

    //constructor 
    DiagonalMatrixMultiply(DENS_MAN * source,
                           DIAG_MAN * diagonalMatrix);
    //destructor
    virtual ~DiagonalMatrixMultiply();

    // apply transfer operator  
    virtual void reset_quantity() const;

  protected:
    
    // pointer to sparseMatrix
    DIAG_MAN * diagonalMatrix_;

  private: 

    // do not define
    DiagonalMatrixMultiply();
  
  };

#ifdef ATC_WHO
/**
  // class sparse matrix multiplier for anything
  //delete later
**/ 
  class SparseMatrixMultiply : public MatToMatTransfer<double> {

  public: 

    //constructor 
    SparseMatrixMultiply(ATC_Method * atc,
                         DENS_MAN * source,
                         SPAR_MAN * sparseMatrix);
    //destructor
    virtual ~SparseMatrixMultiply();

  protected:
    
    // pointer to sparseMatrix
    SPAR_MAN * sparseMatrix_;
    
    // apply transfer operator  
    virtual const DENS_MAT & quantity() const;

  private: 

    // do not define
    SparseMatrixMultiply();
  
  };

#endif

}
#endif

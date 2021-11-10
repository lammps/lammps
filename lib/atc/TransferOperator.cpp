// ATC headers
#include "TransferOperator.h"
#include "ATC_Method.h"
#include "ATC_Coupling.h"
#include "KernelFunction.h"
#include "FE_Mesh.h"
//#include "AtomToMoleculeTransfer.h"

//#include <typeinfo>

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomToFeTransfer
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomToFeTransfer::AtomToFeTransfer(ATC_Method * atc,
                                     PerAtomQuantity<double> * source) :
    DenseMatrixTransfer<double>(),
    atc_(atc),
    source_(source)
  {
    source_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtomToFeTransfer::~AtomToFeTransfer()
  {
    source_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomDiagonalMatrixToFeTransfer
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomDiagonalMatrixToFeTransfer::AtomDiagonalMatrixToFeTransfer(ATC_Method * atc,
                                                                 PerAtomDiagonalMatrix<double> * source) :
    DenseMatrixTransfer<double>(),
    atc_(atc),
    source_(source)
  {
    source_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtomDiagonalMatrixToFeTransfer::~AtomDiagonalMatrixToFeTransfer()
  {
    source_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FeToAtomTransfer
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  FeToAtomTransfer::FeToAtomTransfer(ATC_Method * atc,
                                     DENS_MAN * source,
                                     AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,source->nCols(),atomType),
    source_(source)
  {
    source_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  FeToAtomTransfer::~FeToAtomTransfer()
  {
    source_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FeToAtomDiagonalMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  FeToAtomDiagonalMatrix::FeToAtomDiagonalMatrix(ATC_Method * atc,
                                                 DENS_MAN * source) :
    ProtectedAtomDiagonalMatrix<double>(atc),
    source_(source)
  {
    source_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  FeToAtomDiagonalMatrix::~FeToAtomDiagonalMatrix()
  {
    source_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfShapeFunctionRestriction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtfShapeFunctionRestriction::AtfShapeFunctionRestriction(ATC_Method * atc,
                                                           PerAtomQuantity<double> * source,
                                                           SPAR_MAN * shapeFunction) :
      AtomToFeTransfer(atc,source),
      shapeFunction_(shapeFunction)
  {
    shapeFunction_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtfShapeFunctionRestriction::~AtfShapeFunctionRestriction()
  {
    shapeFunction_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfShapeFunctionRestriction::reset_quantity() const
  {
    global_restriction();
  }

  //--------------------------------------------------------
  //  global_restriction
  //--------------------------------------------------------
  void AtfShapeFunctionRestriction::global_restriction() const
  {
    // computes nodeData = N*atomData where N are the shape functions
    const DENS_MAT & sourceMatrix(source_->quantity());
    // reallocate memory only if sizing has changed
    const SPAR_MAT & shapeFunctionMatrix(shapeFunction_->quantity());
    quantity_.resize(shapeFunctionMatrix.nCols(),sourceMatrix.nCols());

    local_restriction(sourceMatrix,shapeFunctionMatrix);

    // communicate for total restriction
    int count = quantity_.nRows()*quantity_.nCols();
    lammpsInterface_->allsum(_workspace_.ptr(),quantity_.ptr(),count);
  }

  //--------------------------------------------------------
  //  local_restriction
  //--------------------------------------------------------
  void AtfShapeFunctionRestriction::local_restriction(const DENS_MAT & sourceMatrix,
                                                      const SPAR_MAT & shapeFunctionMatrix) const
  {
    if (sourceMatrix.nRows() > 0)
      _workspace_ = shapeFunctionMatrix.transMat(sourceMatrix);
    else
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AdmtfShapeFunctionRestriction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AdmtfShapeFunctionRestriction::AdmtfShapeFunctionRestriction(ATC_Method * atc,
                                                               PerAtomDiagonalMatrix<double> * source,
                                                               SPAR_MAN * shapeFunction) :
      AtomDiagonalMatrixToFeTransfer(atc,source),
      shapeFunction_(shapeFunction)
  {
    shapeFunction_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AdmtfShapeFunctionRestriction::~AdmtfShapeFunctionRestriction()
  {
    shapeFunction_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AdmtfShapeFunctionRestriction::reset_quantity() const
  {
    global_restriction();
  }

  //--------------------------------------------------------
  //  global_restriction
  //--------------------------------------------------------
  void AdmtfShapeFunctionRestriction::global_restriction() const
  {
    // computes nodeData = N*atomData where N are the shape functions
    const CLON_VEC & sourceMatrix(source_->quantity());
    // reallocate memory only if sizing has changed
    const SPAR_MAT & shapeFunctionMatrix(shapeFunction_->quantity());
    quantity_.resize(shapeFunctionMatrix.nCols(),sourceMatrix.nCols());

    local_restriction(sourceMatrix,shapeFunctionMatrix);

    // communicate for total restriction
    int count = quantity_.nRows()*quantity_.nCols();
    lammpsInterface_->allsum(_workspace_.ptr(),quantity_.ptr(),count);
  }

  //--------------------------------------------------------
  //  local_restriction
  //--------------------------------------------------------
  void AdmtfShapeFunctionRestriction::local_restriction(const DENS_MAT & sourceMatrix,
                                                        const SPAR_MAT & shapeFunctionMatrix) const
  {
    if (sourceMatrix.nRows() > 0)
      _workspace_ = shapeFunctionMatrix.transMat(sourceMatrix);
    else
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DiagonalMatrixMultiply
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  DiagonalMatrixMultiply::DiagonalMatrixMultiply(DENS_MAN * source,
                                                 DIAG_MAN * diagonalMatrix) :
    MatToMatTransfer<double>(source),
    diagonalMatrix_(diagonalMatrix)
  {
    diagonalMatrix_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  DiagonalMatrixMultiply::~DiagonalMatrixMultiply()
  {
    diagonalMatrix_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void DiagonalMatrixMultiply::reset_quantity() const
  {
    quantity_ = (diagonalMatrix_->quantity())*(source_->quantity());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfProjection
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtfProjection::AtfProjection(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               SPAR_MAN * accumulant,
                               DIAG_MAN * weights) :
      AtomToFeTransfer(atc,source),
      accumulant_(accumulant),
      weights_(weights)
  {
    if(accumulant_) accumulant_->register_dependence(this);
    if(weights_) weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtfProjection::~AtfProjection()
  {
    if(accumulant_) accumulant_->remove_dependence(this);
    if(weights_) weights_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfProjection::reset_quantity() const
  {
    global_restriction();
  }

  //--------------------------------------------------------
  //  global_restriction
  //--------------------------------------------------------
  void AtfProjection::global_restriction() const
  {
    // computes nodeData = N*atomData where N are the shape functions
    const DENS_MAT & sourceMatrix(source_->quantity());
    // reallocate memory only if sizing has changed
    const SPAR_MAT & accumulantMatrix(accumulant_->quantity());
    quantity_.resize(accumulantMatrix.nCols(),sourceMatrix.nCols());

    local_restriction(sourceMatrix,accumulantMatrix);

    // communicate for total restriction
    int count = quantity_.nRows()*quantity_.nCols();
    lammpsInterface_->allsum(_workspace_.ptr(),quantity_.ptr(),count);
    if (weights_) {
      CLON_VEC w(weights_->quantity());
      quantity_ *= w;
    }
  }

  //--------------------------------------------------------
  //  local_restriction
  //--------------------------------------------------------
  void AtfProjection::local_restriction(const DENS_MAT & sourceMatrix,
                                        const SPAR_MAT & accumulantMatrix) const
  {
    if (sourceMatrix.nRows() > 0)
      _workspace_ = accumulantMatrix.transMat(sourceMatrix);
    else
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
  }
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtfProjectionScaled::AtfProjectionScaled(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               SPAR_MAN * accumulant,
                               const double scale,
                               DIAG_MAN * weights) :
      AtfProjection(atc,source,accumulant,weights),
      scale_(scale)
  {
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtfProjectionScaled::~AtfProjectionScaled()
  {
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfProjectionScaled::reset_quantity() const
  {
    global_restriction();
    quantity_ *= scale_;
  }
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtfProjectionReferenced::AtfProjectionReferenced(ATC_Method * atc,
                               PerAtomQuantity<double> * source,
                               SPAR_MAN * accumulant,
                               DENS_MAN * reference,
                               DIAG_MAN * weights) :
      AtfProjection(atc,source,accumulant,weights),
      reference_(reference)
  {
    if(reference_) reference_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtfProjectionReferenced::~AtfProjectionReferenced()
  {
    if(reference_) reference_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfProjectionReferenced::reset_quantity() const
  {
    global_restriction();
    quantity_ -= (reference_->quantity());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfWeightedShapeFunctionRestriction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  AtfWeightedShapeFunctionRestriction::AtfWeightedShapeFunctionRestriction(ATC_Method * atc,
                                                                           PerAtomQuantity<double> * source,
                                                                           SPAR_MAN * shapeFunction,
                                                                           DIAG_MAN * weights) :
    AtfShapeFunctionRestriction(atc,source,shapeFunction),
    weights_(weights)
  {
    weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  local_restriction
  //--------------------------------------------------------
  void AtfWeightedShapeFunctionRestriction::local_restriction(const DENS_MAT & sourceMatrix,
                                                              const SPAR_MAT & shapeFunctionMatrix) const
  {
    if (sourceMatrix.nRows() > 0) {
      const DIAG_MAT & weightsMatrix(weights_->quantity());
      _workspace_ = shapeFunctionMatrix.transMat(weightsMatrix*sourceMatrix);
    }
    else {
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfNodeWeightedShapeFunctionRestriction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  AtfNodeWeightedShapeFunctionRestriction::AtfNodeWeightedShapeFunctionRestriction(ATC_Method * atc,
                                                                                   PerAtomQuantity<double> * source,
                                                                                   SPAR_MAN * shapeFunction,
                                                                                   DIAG_MAN * weights) :
    AtfShapeFunctionRestriction(atc,source,shapeFunction),
    weights_(weights)
  {
    weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  global_restriction
  //--------------------------------------------------------
  void AtfNodeWeightedShapeFunctionRestriction::global_restriction() const
  {
    AtfShapeFunctionRestriction::global_restriction();
    quantity_ *= weights_->quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfShapeFunctionProjection
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  AtfShapeFunctionProjection::AtfShapeFunctionProjection(ATC_Method * atc,
                                                         DENS_MAN * source,
                                                         FieldName thisField) :
    MatToMatTransfer<double>(source),
    atc_(atc),
    thisField_(thisField)
  {
    DIAG_MAN & massMat(atc_->mass_mat(thisField_));
    massMat.register_dependence(this);
  }

  //--------------------------------------------------------
  // Destructor
  //--------------------------------------------------------
  AtfShapeFunctionProjection::~AtfShapeFunctionProjection()
  {
    DIAG_MAN & massMat(atc_->mass_mat(thisField_));
    massMat.remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfShapeFunctionProjection::reset_quantity() const
  {
    const DENS_MAT & sourceMatrix(source_->quantity());
    atc_->apply_inverse_mass_matrix(sourceMatrix,quantity_,thisField_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfKernelFunctionRestriction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtfKernelFunctionRestriction::AtfKernelFunctionRestriction(ATC_Method * atc,
                                                             PerAtomQuantity<double> * source,
                                                             PerAtomQuantity<double> * coarseGrainingPositions,
                                                             KernelFunction * kernelFunction) :
      AtomToFeTransfer(atc,source),
      kernelFunction_(kernelFunction),
      coarseGrainingPositions_(coarseGrainingPositions),
      feMesh_(((atc->fe_engine())->fe_mesh()))
  {
    coarseGrainingPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtfKernelFunctionRestriction::~AtfKernelFunctionRestriction()
  {
    coarseGrainingPositions_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfKernelFunctionRestriction::reset_quantity() const
  {
      global_restriction();
  }

  //--------------------------------------------------------
  //  global_restriction
  //--------------------------------------------------------
  void AtfKernelFunctionRestriction::global_restriction() const
  {
    // computes nodeData = N*atomData where N are the shape functions
    const DENS_MAT & sourceMatrix(source_->quantity());
    const DENS_MAT & positions(coarseGrainingPositions_->quantity());
    // reallocate memory only if sizing has changed
    quantity_.resize(atc_->num_nodes(),sourceMatrix.nCols());

    local_restriction(sourceMatrix,positions,
                      kernelFunction_);

    // communicate for total restriction
    int count = quantity_.nRows()*quantity_.nCols();
    lammpsInterface_->allsum(_workspace_.ptr(),quantity_.ptr(),count);
  }

  //--------------------------------------------------------
  //  local_restriction
  //--------------------------------------------------------
  void AtfKernelFunctionRestriction::local_restriction(const DENS_MAT & sourceMatrix,
                                                       const DENS_MAT & positions,
                                                       const KernelFunction * kernelFunction) const
  {
    if (sourceMatrix.nRows() > 0) {
      _xI_.resize(positions.nCols()); _xaI_.resize(positions.nCols());
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
      double val;
      for (int i = 0; i < quantity_.nRows(); i++) {
        _xI_ = feMesh_->nodal_coordinates(i);
        for (int j = 0; j < sourceMatrix.nRows(); j++) {
          for (int k = 0; k < positions.nCols(); k++)
            _xaI_(k) = _xI_(k) - positions(j,k);
          atc_->lammps_interface()->periodicity_correction(_xaI_.ptr());
          val = kernelFunction->value(_xaI_);
          if (val > 0) {
            for (int k = 0; k < sourceMatrix.nCols(); k++)
              _workspace_(i,k) += val*sourceMatrix(j,k);
          }
        }
      }
    }
    else
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfWeightedKernelFunctionRestriction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  AtfWeightedKernelFunctionRestriction::AtfWeightedKernelFunctionRestriction(ATC_Method * atc,
                                                                             PerAtomQuantity<double> * source,
                                                                             PerAtomQuantity<double> * coarseGrainingPositions,
                                                                             KernelFunction * kernelFunction,
                                                                             DIAG_MAN * weights) :
    AtfKernelFunctionRestriction(atc,source,coarseGrainingPositions,kernelFunction),
    weights_(weights)
  {
    weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  local_restriction
  //--------------------------------------------------------
  void AtfWeightedKernelFunctionRestriction::local_restriction(const DENS_MAT & sourceMatrix,
                                                               const DENS_MAT & positions,
                                                               const KernelFunction * kernelFunction) const
  {
    const DIAG_MAT & weights(weights_->quantity());
    if (sourceMatrix.nRows() > 0) {
      _xI_.resize(positions.nCols()); _xaI_.resize(positions.nCols());
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
      double val;
      for (int i = 0; i < quantity_.nRows(); i++) {
        _xI_ = feMesh_->nodal_coordinates(i);
        for (int j = 0; j < sourceMatrix.nRows(); j++) {
          for (int k = 0; k < positions.nCols(); k++)
            _xaI_(k) = _xI_(k) - positions(j,k);
          atc_->lammps_interface()->periodicity_correction(_xaI_.ptr());
          val = kernelFunction->value(_xaI_);
          if (val > 0) {
            val *= weights(j,j);
            for (int k = 0; k < sourceMatrix.nCols(); k++)
              _workspace_(i,k) += val*sourceMatrix(j,k);
          }
        }
      }
    }
    else {
      _workspace_.reset(quantity_.nRows(),quantity_.nCols());
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfNodeWeightedKernelFunctionRestriction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  AtfNodeWeightedKernelFunctionRestriction::AtfNodeWeightedKernelFunctionRestriction(ATC_Method * atc,
                                                                                     PerAtomQuantity<double> * source,
                                                                                     PerAtomQuantity<double> * coarseGrainingPositions,
                                                                                     KernelFunction * kernelFunction,
                                                                                     DIAG_MAN * weights) :
    AtfKernelFunctionRestriction(atc,source,coarseGrainingPositions,kernelFunction),
    weights_(weights)
  {
    weights_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  global_restriction
  //--------------------------------------------------------
  void AtfNodeWeightedKernelFunctionRestriction::global_restriction() const
  {
    AtfKernelFunctionRestriction::global_restriction();
    quantity_ *= weights_->quantity();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfShapeFunctionMdProjection
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  AtfShapeFunctionMdProjection::AtfShapeFunctionMdProjection(ATC_Method * atc,
                                                             DENS_MAN * source,
                                                             FieldName thisField) :
    MatToMatTransfer<double>(source),
    atc_(atc),
    thisField_(thisField)
  {
    DIAG_MAN & massMat(atc_->set_mass_mat_md(thisField_));
    massMat.register_dependence(this);
  }

  //--------------------------------------------------------
  // Destructor
  //--------------------------------------------------------
  AtfShapeFunctionMdProjection::~AtfShapeFunctionMdProjection()
  {
    DIAG_MAN & massMat(atc_->set_mass_mat_md(thisField_));
    massMat.remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfShapeFunctionMdProjection::reset_quantity() const
  {
    const DENS_MAT & sourceMatrix(source_->quantity());
    atc_->apply_inverse_md_mass_matrix(sourceMatrix,quantity_,thisField_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtfShapeFunctionMdProjectionScaled
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  AtfShapeFunctionMdProjectionScaled::AtfShapeFunctionMdProjectionScaled(ATC_Method * atc,
                                                                         DENS_MAN * source,
                                                                         double scale,
                                                                         FieldName thisField) :
    AtfShapeFunctionMdProjection(atc,source,thisField),
    scale_(scale)
  {
  }

  //--------------------------------------------------------
  // Destructor
  //--------------------------------------------------------
  AtfShapeFunctionMdProjectionScaled::~AtfShapeFunctionMdProjectionScaled()
  {
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfShapeFunctionMdProjectionScaled::reset_quantity() const
  {
    AtfShapeFunctionMdProjection::reset_quantity();
    quantity_ *= scale_;
  }

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtfShapeFunctionMdProjectionReferenced::AtfShapeFunctionMdProjectionReferenced(ATC_Method * atc,
    DENS_MAN * source,
    DENS_MAN * reference,
    FieldName thisField) :
    AtfShapeFunctionMdProjection(atc,source,thisField),
    reference_(reference)
  {
    reference_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtfShapeFunctionMdProjectionReferenced::~AtfShapeFunctionMdProjectionReferenced()
  {
    reference_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void AtfShapeFunctionMdProjectionReferenced::reset_quantity() const
  {
    quantity_ = source_->quantity();
    quantity_ -= reference_->quantity();
    atc_->apply_inverse_md_mass_matrix(quantity_,thisField_);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FtaShapeFunctionProlongation
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  FtaShapeFunctionProlongation::FtaShapeFunctionProlongation(ATC_Method * atc,
                                                             DENS_MAN * source,
                                                             SPAR_MAN * shapeFunction,
                                                             AtomType atomType) :
    FeToAtomTransfer(atc,source,atomType),
    shapeFunction_(shapeFunction)
  {
    shapeFunction_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  FtaShapeFunctionProlongation::~FtaShapeFunctionProlongation()
  {
    shapeFunction_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void FtaShapeFunctionProlongation::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      if (atc_.nlocal() > 0) {
        const DENS_MAT & sourceMatrix(source_->quantity());
        const SPAR_MAT & shapeFunctionMatrix(shapeFunction_->quantity());
        quantity_ = shapeFunctionMatrix*sourceMatrix;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FtaShapeFunctionProlongationDiagonalMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  FtaShapeFunctionProlongationDiagonalMatrix::FtaShapeFunctionProlongationDiagonalMatrix(ATC_Method * atc,
                                                                                         DENS_MAN * source,
                                                                                         SPAR_MAN * shapeFunction) :
    FeToAtomDiagonalMatrix(atc,source),
    shapeFunction_(shapeFunction)
  {
    shapeFunction_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  FtaShapeFunctionProlongationDiagonalMatrix::~FtaShapeFunctionProlongationDiagonalMatrix()
  {
    shapeFunction_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void FtaShapeFunctionProlongationDiagonalMatrix::reset() const
  {
    if (need_reset()) {
      PerAtomDiagonalMatrix<double>::reset();
      if (atc_.nlocal() > 0) {
        const DENS_MAT & sourceMatrix(source_->quantity());
        const SPAR_MAT & shapeFunctionMatrix(shapeFunction_->quantity());
        _temp_ = shapeFunctionMatrix*sourceMatrix;
        for (int i = 0; i < quantity_.size(); ++i) {
          quantity_(i,i) = _temp_(i,0);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MatToGradBySparse
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  // Constructor
  //--------------------------------------------------------
  MatToGradBySparse::MatToGradBySparse(ATC_Method * /* atc */,
                                       DENS_MAN * source,
                                       VectorDependencyManager<SPAR_MAT * > * gradientMatrices) :
    MatToMatTransfer<double>(source),
    gradientMatrices_(gradientMatrices)
  {
    gradientMatrices_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  MatToGradBySparse::~MatToGradBySparse()
  {
    gradientMatrices_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void MatToGradBySparse::reset_quantity() const
  {
    const DENS_MAT & source(source_->quantity());
    const SPAR_MAT_VEC & gradientMatrices(gradientMatrices_->quantity());
    int nsd = gradientMatrices.size();
    int nNodes = source.nRows();
    int nCols = source.nCols();
    quantity_.reset(nNodes,nCols*nsd);

    int index = 0;
    for (int n = 0; n < nCols; n++) {
      for (int m = 0; m < nsd; m++) {
        CLON_VEC inData(source,CLONE_COL,n);
        CLON_VEC outData(quantity_,CLONE_COL,index);
        outData = (*(gradientMatrices[m]))*inData;
        ++index;
      }
    }
  }

}

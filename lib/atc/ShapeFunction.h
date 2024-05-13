// a data structure for managing interpolation/restriction operations

#ifndef SHAPE_FUNCTION_H
#define SHAPE_FUNCTION_H

#include "DependencyManager.h"
#include "PaqAtcUtility.h"
#include <map>
#include <vector>
#include <set>
#include <pair>

namespace ATC {

  // forward declarations
  class FE_Engine;
  class ATC_Method;
  class LammpsInterface;

  // type defs
  typedef pair< vector<double>, vector<int> > POINT_VAL;
  typedef pair< vector<vector<double> >, vector<int> > POINT_GRAD;

  /**
   *  @class FeEngineInterface
   *  @class Base class for defining interfaces to the finite element engine to handle different shape functions
   */

  class FeEngineInterface {

  public:

    // constructor
    FeEngineInterface(FE_Engine * feEngine, DENS_MAN * coordinates) : feEngine_(feEngine), coordinates_(coordinates) {};

    // destructor
    virtual ~FeEngineInterface() {};

    /** evaluate shape function at a set of coordinates */
    virtual void evaluate_shape_function(int index,
                                         POINT_VAL & data) = 0;

  protected:

    /** pointer to the engine */
    FE_Engine * feEngine_;

    /** quantity defining locations */
    DENS_MAN * coordinates_;

  private:

    // do not define
    FeEngineInterface();

  };

  /**
   *  @class FeEngineInterfacePu
   *  @class Interfaces to the finite element engine to handle partition of unity (PU) shape functions
   */

  class FeEngineInterfacePu : public FeEngineInterface {

  public:

    // constructor
    FeEngineInterfacePu(FE_Engine * feEngine) : FeEngineInterface(feEngine) {};

    // destructor
    virtual ~FeEngineInterfacePu() {};

    /** evaluate shape function at a set of coordinates */
    virtual void evaluate_shape_function(int index,
                                         POINT_VAL & data);

  protected:

    /** pointer to the point to element map */
    MatrixDependencyManager<DenseMatrix, int> * pointToElementMap_;

  private:

    // do not define
    FeEngineInterfacePu();

  };

  /**
   *  @class FeEngineInterfaceMls
   *  @class Interfaces to the finite element engine to handle moving least squares (MLS) shape functions
   */

  class FeEngineInterfaceMls : public FeEngineInterface {

  public:

    // constructor
    FeEngineInterfaceMls(FE_Engine * feEngine) : FeEngineInterface(feEngine) {};

    // destructor
    virtual ~FeEngineInterfaceMls() {};

    /** evaluate shape function at a set of coordinates */
    virtual void evaluate_shape_function(int index,
                                         POINT_VAL & data);

  protected:


  private:

    // do not define
    FeEngineInterfaceMls();

  };

  /**
   *  @class ShapeFunctionBase
   *  @class Base class for defining shape functions for restriction and interpolation
   */

  class ShapeFunctionBase : public DependencyManager {

  public:

    // constructor
    // discriminate engine interface based on point to element map presence?
    // could also add a manager and pass the interface in on construction
    ShapeFunctionBase(FE_Engine * feEngine,
                      DENS_MAN * coordinates) : DependencyManager(), feEngineInterface_(feEngine), coordinates_(coordinates)
      {coordinates_->register_dependence(this)};

    // destructor
    virtual ~ShapeFunctionBase() {coordinates_->remove_dependence(this)};

    /** use the shape function as a restriction operator */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output) = 0;

    /** apply the restriction on a subset of points */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const set<int> & points) = 0;

    /** apply the restriction on a subset of points discriminating with booleans */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const Array<bool> & points) = 0;

    /** use the shape function as an interpolation operator */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output) = 0;

    /** apply the interpolation on a subset of points */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const set<int> & points) = 0;

    /** apply the interpolation on a subset of points discriminating with booleans */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const Array<bool> & points) = 0;

  protected:

    /** apply a reset if needed */
    virtual void reset() = 0;

    /** object to interface with the engine */
    FeEngineInterface * feEngineInterface_;

    /** quantity defining locations */
    DENS_MAN * coordinates_;

  private:

    // do not define
    ShapeFunctionBase();

  };

  /**
   *  @class ShapeFunction
   *  @class Defines general shape functions for restriction and interpolation
   */

  class ShapeFunction : public ShapeFunctionBase {

  public:

    // constructor
    // discriminate engine interface based on point to element map presence?
    ShapeFunction(FE_Engine * feEngine,
                  DENS_MAN * coordinates) : ShapeFunctionBase(feEngine,coordinates) {values.reserve(coordinates->nRows();};

    // destructor
    virtual ~ShapeFunction() {};

    /** use the shape function as a restriction operator */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output);

    /** apply the restriction on a subset of points */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const set<int> & points);

    /** apply the restriction on a subset of points discriminating with booleans */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const Array<bool> & points);

    /** use the shape function as an interpolation operator */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output);

    /** apply the interpolation on a subset of points */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const set<int> & points);

    /** apply the interpolation on a subset of points discriminating with booleans */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const Array<bool> & points);

  protected:

    /** apply a reset if needed */
    virtual void reset();

    /** storage for shape function values, indexed by rows of coordinates */
    vector<POINT_VAL> values;

  private:

    // do not define
    ShapeFunction();

  };

  /**
   *  @class ShapeFunctionAtomic
   *  @class Defines shape functions for restriction and interpolation for atomic data
   */

  class ShapeFunctionAtomic : public ShapeFunctionBase {

  public:

    // constructor
    // discriminate engine interface based on point to element map presence?
    ShapeFunctionAtomic(ATC_Method * atc,
                        DENS_MAN * coordinates,
                        AtomType atomType) : ShapeFunctionBase(atc->fe_engine(),coordinates), atc_(atc,atomType), quantityToLammps_(atc_.atc_to_lammps_map()) {};

    // destructor
    virtual ~ShapeFunctionAtomic() {};

    /** use the shape function as a restriction operator */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output);

    /** apply the restriction on a subset of points */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const set<int> & points);

    /** apply the restriction on a subset of points discriminating with booleans */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const Array<bool> & points);

    /** use the shape function as an interpolation operator */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output);

    /** apply the interpolation on a subset of points */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const set<int> & points);

    /** apply the interpolation on a subset of points discriminating with booleans */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const Array<bool> & points);

  protected:

    /** apply a reset if needed */
    virtual void reset();

    /** storage for shape function values, map is based on lammps global index */
    map<int,POINT_VAL> values;

    /** interface to atc functionality */
    PaqAtcUtility atc_;

    /** interface to lammps functionality */
    LammpsInterface * lammpsInterface_;

    /** map from this quantity's AtC indexing to Lammps local indexing for atomic arrays */
    const Array<int> & quantityToLammps_;

  private:

    // do not define
    ShapeFunctionAtomic();

  };

  /**
   *  @class ShapeFunctionAtomicMask
   *  @class Defines shape functions for restriction and interpolation for atomic data, but only on a subset of either ghost or internal
   */

  class ShapeFunctionAtomicMask : public ShapeFunctionAtomic {

  public:

    // constructor
    // discriminate engine interface based on point to element map presence?
    ShapeFunctionAtomicMask(ATC_Method * atc,
                            DENS_MAN * coordinates,
                            AtomType atomType,
                            PerAtomQuantity<bool> * mask) : ShapeFunctionAtomic(atc,coordinates,atomType), mask_(mask) {mask_->register_dependence(this);};

    // destructor
      virtual ~ShapeFunctionAtomicMask() {mask_->remove_dependence(this);};



  protected:

    /** apply a reset if needed */
    virtual void reset();

    /** mask to screen out atoms */
    PerAtomQuantity<bool> * mask_;

  private:

    // do not define
    ShapeFunctionAtomicMask();

  };

    /**
   *  @class ShapeFunctionGrad
   *  @class Defines general shape gradients functions for restriction and interpolation
   */

  class ShapeFunctionGrad : public ShapeFunctionBase {

  public:

    // constructor
    // discriminate engine interface based on point to element map presence?
    ShapeFunctionGrad(FE_Engine * feEngine,
                      DENS_MAN * coordinates) : ShapeFunctionBase(feEngine,coordinates) {values.reserve(coordinates->nRows();};

    // destructor
    virtual ~ShapeFunctionGrad() {};

    /** use the shape function as a restriction operator */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output);

    /** apply the restriction on a subset of points */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const set<int> & points);

    /** apply the restriction on a subset of points discriminating with booleans */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const Array<bool> & points);

    /** use the shape function as an interpolation operator */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output);

    /** apply the interpolation on a subset of points */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const set<int> & points);

    /** apply the interpolation on a subset of points discriminating with booleans */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const Array<bool> & points);

  protected:

    /** apply a reset if needed */
    virtual void reset();

    /** storage for shape function gradients, indexed by rows of coordinates */
    vector<POINT_GRAD> values;

  private:

    // do not define
    ShapeFunctionGrad();

  };

  /**
   *  @class ShapeFunctionGradAtomic
   *  @class Defines shape functions gradients for restriction and interpolation for atomic data
   */

  class ShapeFunctionGradAtomic : public ShapeFunctionBase {

  public:

    // constructor
    // discriminate engine interface based on point to element map presence?
    ShapeFunctionGradAtomic(ATC_Method * atc,
                            DENS_MAN * coordinates,
                            AtomType atomType) : ShapeFunctionBase(atc->fe_engine(),coordinates), atc_(atc,atomType), quantityToLammps_(atc_.atc_to_lammps_map()) {};

    // destructor
    virtual ~ShapeFunctionGradAtomic() {};

    /** use the shape function as a restriction operator */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output);

    /** apply the restriction on a subset of points */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const set<int> & points);

    /** apply the restriction on a subset of points discriminating with booleans */
    virtual void restrict(const DENS_MAT & input,
                          DENS_MAT & output,
                          const Array<bool> & points);

    /** use the shape function as an interpolation operator */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output);

    /** apply the interpolation on a subset of points */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const set<int> & points);

    /** apply the interpolation on a subset of points discriminating with booleans */
    virtual void interpolation(const DENS_MAT & input,
                               DENS_MAT & output,
                               const Array<bool> & points);

  protected:

    /** apply a reset if needed */
    virtual void reset();

    /** storage for shape function values, map is based on lammps global index */
    map<int,POINT_GRAD> values;

    /** interface to atc functionality */
    PaqAtcUtility atc_;

    /** interface to lammps functionality */
    LammpsInterface * lammpsInterface_;

    /** map from this quantity's AtC indexing to Lammps local indexing for atomic arrays */
    const Array<int> & quantityToLammps_;

  private:

    // do not define
    ShapeFunctionGradAtomic();

  };

  /**
   *  @class ShapeFunctionGradAtomicMask
   *  @class Defines shape function gradients s for restriction and interpolation for atomic data, but only on a subset of either ghost or internal
   */

  class ShapeFunctionGradAtomicMask : public ShapeFunctionGradAtomic {

  public:

    // constructor
    // discriminate engine interface based on point to element map presence?
    ShapeFunctionGradAtomicMask(ATC_Method * atc,
                            DENS_MAN * coordinates,
                            AtomType atomType,
                            PerAtomQuantity<bool> * mask) : ShapeFunctionGradAtomic(atc,coordinates,atomType), mask_(mask) {mask_->register_dependence(this);};

    // destructor
      virtual ~ShapeFunctionGradAtomicMask() {mask_->remove_dependence(this);};



  protected:

    /** apply a reset if needed */
    virtual void reset();

    /** mask to screen out atoms */
    PerAtomQuantity<bool> * mask_;

  private:

    // do not define
    ShapeFunctionGradAtomicMask();

  };

};





#endif

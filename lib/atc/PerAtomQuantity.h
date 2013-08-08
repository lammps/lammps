// A class for defining atomic quantities for interscale operations

#ifndef PER_ATOM_QUANTITY_H
#define PER_ATOM_QUANTITY_H

// ATC_Method headers
#include "LammpsInterface.h"
#include "DependencyManager.h"
#include "PaqAtcUtility.h"
#include <set>
#include <vector>

using namespace std;

namespace ATC {

  // forward declarations
  class ATC_Method;
  template <typename TT> class ClonedAtomQuantity;

  /**
   *  @class  PerAtomQuantity 
   *  @brief  Base class for objects that manage atomic quantities and their AtC interface
   */
  template <typename T>
  class PerAtomQuantity : public MatrixDependencyManager<DenseMatrix, T> {

  public:
    
    // constructor
    PerAtomQuantity(ATC_Method * atc, int nCols = 1, AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~PerAtomQuantity();

    /** access to a constant dense matrix of the quantity, indexed by prescribed counts */
    virtual const DenseMatrix<T> & quantity() const
      {reset(); return MatrixDependencyManager<DenseMatrix, T>::quantity();};

    /** access to a non-constant dens matrix of the quantity, indexed by prescribed atom counts */
    virtual DenseMatrix<T> & set_quantity()
      {reset(); return MatrixDependencyManager<DenseMatrix, T>::set_quantity();}

    /** number of columns in quantity */
    INDEX nCols() const {return nCols_;};

    /** sets the Lammps quantity to a given value, input is indexed by AtC atom counts */
    virtual void operator=(const DenseMatrix<T> & target) {PerAtomQuantity<T>::reset(); MatrixDependencyManager<DenseMatrix, T>::operator=(target);};

    /** adds the given data to the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator+=(const DenseMatrix<T> & addition) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator+=(addition);};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator+=(addition);};

    /** subtracts the given data from the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator-=(const DenseMatrix<T> & subtraction) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator-=(subtraction);};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(double subtraction) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator-=(subtraction);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DenseMatrix<T> & multiplier) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator*=(multiplier);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator*=(multiplier);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DenseMatrix<T> & divisor) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator/=(divisor);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor) {reset(); MatrixDependencyManager<DenseMatrix, T>::operator/=(divisor);};

    /** sets quantity to lammps data, if needed, should be called in pre_exchange */
    virtual void prepare_exchange() {this->set_lammps_to_quantity();};

    /** resets AtC local quantity after exchange */
    virtual void post_exchange() {(this->quantity_).resize(atc_.nlocal(),nCols_); this->set_quantity_to_lammps();};
    
    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return nCols_;};

    /** packs up data for parallel transfer when atoms change processors */
    virtual int pack_exchange(int i, double *buffer);

    /** unpacks data after parallel transfer when atoms change processors */
    virtual int unpack_exchange(int i, double *buffer);

    // pack/unpack_comm only valid if the quantity is over all real and processor ghost atoms
    /** packs up data for parallel transfer to ghost atoms on other processors */
    virtual int pack_comm(int index, double *buf,
                          int pbc_flag, int *pbc);

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    virtual int unpack_comm(int index, double *buf);

    /** returns per-atom size of communicated data */
    virtual int size_comm() const {return nCols_;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const string & tag);

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j);

    /** access type of data this quantity is applied to */
    AtomType atom_type() const {return atomType_;};

    /** specialized reset to account for quantities which lammps can change */
    virtual void lammps_force_reset() {};

    /** resets local storage */
    virtual void reset_nlocal() { this->force_reset();}

  protected:
    
    /** utility object to access ATC methods */
    PaqAtcUtility atc_;

    /** pointer to access Lammps data */
    LammpsInterface * lammpsInterface_;

    /** type of atoms this quantity applies to */
    AtomType atomType_;

    /** number of columns of the per atom quantity, must be defined in derived class */
    int nCols_;

    /** map from this quantity's AtC indexing to Lammps indexing for atomic arrays */
    const Array<int> & quantityToLammps_;

    /** resets data, if necessary */
    virtual void reset() const {if (this->needReset_) {(this->quantity_).resize(atc_.nlocal(),nCols_); this->needReset_ = false;}};

    /** sets the quantity based on a lammps pointer */
    virtual void set_lammps_to_quantity() const;

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const;

    
    /** gets appropriate data for lammps pointer */
    virtual T * lammps_scalar() const = 0;

    /** gets appropriate data for lammps pointer */
    virtual T ** lammps_vector() const = 0;

    /** point to lammps-style array for data */
    T * lammpsScalar_;

    /** pointer to lammps-style double array for data */
    T ** lammpsVector_;

  private:
    
    // do not define
    PerAtomQuantity();
    
  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LammpsAtomQuantity
  //    A base class for defining objects that manage
  //    quantities but the lammps data forms the 
  //    absolute definition for the contained data.
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class LammpsAtomQuantity : public PerAtomQuantity<T> {

  public:

    // used as a friend so that clones can use it's methods
    template <typename TT>
    friend class ClonedAtomQuantity;

    // constructor
    LammpsAtomQuantity(ATC_Method * atc,int nCols = 1, AtomType atomType = INTERNAL) :
      PerAtomQuantity<T>(atc,nCols,atomType) {};

    // destructor
    virtual ~LammpsAtomQuantity() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual DenseMatrix<T> & set_quantity()
      {throw ATC_Error("LammpsAtomQuantity::set_quantity - Cannot modify shallow per atom quantities outside of manager class"); return this->quantity_;};

    /** sets the Lammps quantity to a given value, input is indexed by AtC atom counts */
    virtual void operator=(const DenseMatrix<T> & target)
      {PerAtomQuantity<T>::operator=(target); this->set_lammps_to_quantity();};

    /** adds the given data to the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator+=(const DenseMatrix<T> & addition)
      {PerAtomQuantity<T>::operator+=(addition); this->set_lammps_to_quantity();};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {PerAtomQuantity<T>::operator+=(addition); this->set_lammps_to_quantity();};

    /** subtracts the given data from the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator-=(const DenseMatrix<T> & subtraction)
      {PerAtomQuantity<T>::operator-=(subtraction); this->set_lammps_to_quantity();};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {PerAtomQuantity<T>::operator-=(subtraction); this->set_lammps_to_quantity();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DenseMatrix<T> & multiplier)
      {PerAtomQuantity<T>::operator*=(multiplier); this->set_lammps_to_quantity();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {PerAtomQuantity<T>::operator*=(multiplier); this->set_lammps_to_quantity();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DenseMatrix<T> & divisor)
      {PerAtomQuantity<T>::operator/=(divisor); this->set_lammps_to_quantity();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {PerAtomQuantity<T>::operator/=(divisor); this->set_lammps_to_quantity();};

    /** sets quantity to lammps data, but not needed */
    virtual void prepare_exchange() {};

    /** packs data for parallel transfer to ghost atoms on other processors */
    virtual int pack_comm(int index, double *buf,
                          int pbc_flag, int *pbc);

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    virtual int unpack_comm(int index, double *buf);

  protected:

    virtual void reset() const;

    
    /** gets appropriate data for lammps pointer */
    virtual T * lammps_scalar() const = 0;

    /** gets appropriate data for lammps pointer */
    virtual T ** lammps_vector() const = 0;

  private:

    // do not define
    LammpsAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ShallowAtomQuantity
  //    A base class for defining objects that manage
  //    quantities but do not own their own lammps data.
  //    The lammps data forms the absolute definition
  //    for the contained data.
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ShallowAtomQuantity : public LammpsAtomQuantity<T> {

  public:

    // constructor
    ShallowAtomQuantity(ATC_Method * atc,int nCols = 1, AtomType atomType = INTERNAL) :
      LammpsAtomQuantity<T>(atc,nCols,atomType) {};

    // destructor
    virtual ~ShallowAtomQuantity() {};

    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return 0;};

    /** packs up data for parallel transfer when atoms change processors */
    virtual int pack_exchange(int i, double *buffer) {return 0;};

    /** unpacks data after parallel transfer when atoms change processors */
    virtual int unpack_exchange(int i, double *buffer) {return 0;};

    /** packs up data for parallel transfer to ghost atoms on other processors */
    virtual int pack_comm(int index, double *buf, 
                          int pbc_flag, int *pbc) {return 0;};

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    virtual int unpack_comm(int index, double *buf) {return 0;};

    /** returns size of per-atom communication */
    virtual int size_comm() {return 0;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const string & tag) {};

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j) {};

  protected:


  private:

    // do not define
    ShallowAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ClonedLammpsAtomQuantity
  //    A base class for defining objects that manage
  //    quantities defined at atoms based on data in 
  //    a LammpsAtomQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  template <typename T>
  class ClonedAtomQuantity : public ShallowAtomQuantity<T> {

  public:

    // constructor
    ClonedAtomQuantity(ATC_Method * atc,
                       LammpsAtomQuantity<T> * target,
                       AtomType atomType=INTERNAL) :
      ShallowAtomQuantity<T>(atc,target->nCols(),atomType),
      target_(target)
      {
        target_->register_dependence(this);
      };

    // destructor
    virtual ~ClonedAtomQuantity() {};

  protected:

    /** reference to originating per atom quantity */
    LammpsAtomQuantity<T> * target_;

    /** sets quantity based on root lammps data */
    virtual void set_quantity_to_lammps() const
    {
      target_->reset(); // make sure target is all the way up to date
      ShallowAtomQuantity<T>::set_quantity_to_lammps(); // change appropriate data
    }

    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const
      {
        target_->reset(); // make sure target is all the way up to date
        ShallowAtomQuantity<T>::set_lammps_to_quantity(); // change appropriate data
        target_->propagate_reset();
        this->needReset_ = false;
      };

    /** gets appropriate pointer for lammps data */
    virtual T * lammps_scalar() const
      {return target_->lammps_scalar();};

    /** gets appropriate pointer for lammps data */
    virtual T ** lammps_vector() const
      {return target_->lammps_vector();};

  private:

    // do not define
    ClonedAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedClonedAtomQuantity
  //    A base class for defining objects that manage
  //    quantities defined at atoms based on data in 
  //    a pointer managed in the standard lammps way.
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  template <typename T>
  class ProtectedClonedAtomQuantity : public ShallowAtomQuantity<T> {

  public:

    // constructor
    ProtectedClonedAtomQuantity(ATC_Method * atc,
                                int nCols,
                                AtomType atomType=INTERNAL) :
      ShallowAtomQuantity<T>(atc,nCols,atomType) {};

    // destructor
    virtual ~ProtectedClonedAtomQuantity() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual DenseMatrix<T> & set_quantity()
      {throw ATC_Error("ProtectedClonedAtomQuantity::set_quantity - Cannot modify protected per atom quantities"); return this->quantity_;};

    /** sets the quantity to a given value */
    virtual void operator=(const DenseMatrix<T> & target)
      {throw ATC_Error("ProtectedClonedAtomQuantity::set_quantity - Cannot modify protected per atom quantities");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator= - Cannot modify protected per atom quantities");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DenseMatrix<T> & addition)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator+= - Cannot modify protected per atom quantities");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator+= - Cannot modify protected per atom quantities");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DenseMatrix<T> & subtraction)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator-= - Cannot modify protected per atom quantities");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator-= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DenseMatrix<T> & multiplier)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DenseMatrix<T> & divisor)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator/= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("ProtectedClonedAtomQuantity::operator/= - Cannot modify protected per atom quantities");};

  protected:

    /** sets lammps data based on the quantity, not needed by this class */
    virtual void set_lammps_to_quantity() const {};

    /** gets appropriate pointer for lammps data */
    virtual T * lammps_scalar() const {return NULL;};

    /** gets appropriate pointer for lammps data */
    virtual T ** lammps_vector() const {return NULL;};

  private:

    // do not define
    ProtectedClonedAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtcAtomQuantity
  //    A funcational base class for defining objects that 
  //    manage quantities defined at atoms and their AtC 
  //    interface that are defined by AtC classes
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class AtcAtomQuantity : public PerAtomQuantity<T> {

  public:

    // constructor
    AtcAtomQuantity(ATC_Method * atc,
                    int nCols = 1,
                    AtomType atomType = INTERNAL) : PerAtomQuantity<T>(atc,nCols,atomType) {};

    // destructor
    virtual ~AtcAtomQuantity() {};

  protected:

    // these get the appropriate pointer to local data mimicing lammps storage
    /** gets appropriate data for lammps pointer */
    virtual T * lammps_scalar() const {return this->lammpsScalar_;};

    /** gets appropriate data for lammps pointer */
    virtual T ** lammps_vector() const {return this->lammpsVector_;};

  private:

    // do not define
    AtcAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedAtomQuantity
  //    A base class for defining objects that manage
  //    quantities defined at atoms internally and do not
  //    allow for reset externally
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ProtectedAtomQuantity : public AtcAtomQuantity<T> {

  public:

    // constructor
    ProtectedAtomQuantity(ATC_Method * atc,
                          int nCols = 1,
                          AtomType atomType = INTERNAL)
      : AtcAtomQuantity<T>(atc, nCols, atomType) {};

    // destructor
    virtual ~ProtectedAtomQuantity() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual DenseMatrix<T> & set_quantity()
      {throw ATC_Error("ProtectedAtomQuantity::set_quantity - Cannot modify protected per atom quantities"); return this->quantity_;};

    /** sets the quantity to a given value */
    virtual void operator=(const DenseMatrix<T> & target)
      {throw ATC_Error("ProtectedAtomQuantity::set_quantity - Cannot modify protected per atom quantities");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("ProtectedAtomQuantity::operator= - Cannot modify protected per atom quantities");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DenseMatrix<T> & addition)
      {throw ATC_Error("ProtectedAtomQuantity::operator+= - Cannot modify protected per atom quantities");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("ProtectedAtomQuantity::operator+= - Cannot modify protected per atom quantities");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DenseMatrix<T> & subtraction)
      {throw ATC_Error("ProtectedAtomQuantity::operator-= - Cannot modify protected per atom quantities");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("ProtectedAtomQuantity::operator-= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DenseMatrix<T> & multiplier)
      {throw ATC_Error("ProtectedAtomQuantity::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("ProtectedAtomQuantity::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DenseMatrix<T> & divisor)
      {throw ATC_Error("ProtectedAtomQuantity::operator/= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("ProtectedAtomQuantity::operator/= - Cannot modify protected per atom quantities");};

  protected:

    /** resets the data if necessary */
    virtual void reset() const = 0;

    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const {this->reset(); PerAtomQuantity<T>::set_lammps_to_quantity();};

  private:

    // do not define
    ProtectedAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedMappedAtomQuantity
  //    A base class for defining objects that manage
  //    quantities defined at atoms internally and do not
  //    allow for reset externally, but are mapped onto a
  //    subset of internal atoms
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ProtectedMappedAtomQuantity : public ProtectedAtomQuantity<T> {

  public:

    // constructor
    ProtectedMappedAtomQuantity(ATC_Method * atc,
                                MatrixDependencyManager<DenseMatrix, int> * atomMap,
                                int nCols = 1,
                                AtomType atomType = INTERNAL)
      : ProtectedAtomQuantity<T>(atc, nCols, atomType), atomMap_(atomMap)
        {atomMap_->register_dependence(this);};

    // destructor
    virtual ~ProtectedMappedAtomQuantity() {atomMap_->remove_dependence(this);};

  protected:

    // methods
    /** resets data, if necessary */
    virtual void reset() const {if (this->needReset_) {(this->quantity_).resize(atomMap_->size(),this->nCols()); this->needReset_ = false;}};

    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const;

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const;

    // data
    /** map from atoms of atomType to desired subset */
    MatrixDependencyManager<DenseMatrix, int> * atomMap_;

  private:

    // do not define
    ProtectedMappedAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedLammpsAtomQuantity
  //    A base class for defining objects that manage
  //    quantities defined at atoms using lammps storage as
  //    the definition and do not allow for reset externally
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ProtectedLammpsAtomQuantity : public LammpsAtomQuantity<T> {

  public:

    // constructor
    ProtectedLammpsAtomQuantity(ATC_Method * atc,
                                int nCols = 1,
                                AtomType atomType = INTERNAL)
      : LammpsAtomQuantity<T>(atc, nCols, atomType) {};

    // destructor
    virtual ~ProtectedLammpsAtomQuantity() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual DenseMatrix<T> & set_quantity()
      {throw ATC_Error("ProtectedLammpsAtomQuantity::set_quantity - Cannot modify protected per atom quantities"); return this->quantity_;};

        /** sets the quantity to a given value */
    virtual void operator=(const DenseMatrix<T> & target)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::set_quantity - Cannot modify protected per atom quantities");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator= - Cannot modify protected per atom quantities");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DenseMatrix<T> & addition)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator+= - Cannot modify protected per atom quantities");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator+= - Cannot modify protected per atom quantities");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DenseMatrix<T> & subtraction)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator-= - Cannot modify protected per atom quantities");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator-= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DenseMatrix<T> & multiplier)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DenseMatrix<T> & divisor)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator/= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("ProtectedLammpsAtomQuantity::operator/= - Cannot modify protected per atom quantities");};

  protected:

    /** resets the data if necessary */
    virtual void reset() const = 0;

    // these get the appropriate pointer to local data mimicing lammps storage
    /** gets appropriate data for lammps pointer */
    virtual T * lammps_scalar() const {return this->lammpsScalar_;};

    /** gets appropriate data for lammps pointer */
    virtual T ** lammps_vector() const {return this->lammpsVector_;};

  private:

    // do not define
    ProtectedLammpsAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ConstantQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ConstantQuantity : public ProtectedAtomQuantity<T> {

  public:

    // constructor
    ConstantQuantity(ATC_Method * atc,
                     T constant,
                     int nCols = 1,
                     AtomType atomType = INTERNAL) : ProtectedAtomQuantity<T>(atc,nCols,atomType), constant_(constant) {};

    // destructor
    virtual ~ConstantQuantity() {};

    // eliminate MPI, use resets instead
    /** sets quantity to lammps data, if needed, should be called in pre_exchange */
    virtual void prepare_exchange() {};

    /** resets local AtC storage */
    virtual void post_exchange() {this->set_reset();};

    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return 0;};

    /** packs up data for parallel transfer */
    virtual int pack_exchange(int i, double *buffer) {return 0;};

    /** unpacks data after parallel transfer */
    virtual int unpack_exchange(int i, double *buffer) {return 0;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const string & tag) {};

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j) {};

  protected:

    /** resets the data if necessary */
    virtual void reset() const {if (this->need_reset()) {PerAtomQuantity<T>::reset(); this->quantity_ = constant_;}};

    /** constant to set data to */
    T constant_;

  private:

    // do not define
    ConstantQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ConstantQuantityMapped
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ConstantQuantityMapped : public ProtectedMappedAtomQuantity<T> {

  public:

    // constructor
    ConstantQuantityMapped(ATC_Method * atc,
                           T constant,
                           MatrixDependencyManager<DenseMatrix, int> * atomMap,
                           int nCols = 1,
                           AtomType atomType = INTERNAL) : ProtectedMappedAtomQuantity<T>(atc,atomMap,nCols,atomType), constant_(constant) {};

    // destructor
    virtual ~ConstantQuantityMapped() {};

    // eliminate MPI, use resets instead
    /** sets quantity to lammps data, if needed, should be called in pre_exchange */
    virtual void prepare_exchange() {};

    /** resets local AtC storage */
    virtual void post_exchange() {this->set_reset();};

    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return 0;};

    /** packs up data for parallel transfer */
    virtual int pack_exchange(int i, double *buffer) {return 0;};

    /** unpacks data after parallel transfer */
    virtual int unpack_exchange(int i, double *buffer) {return 0;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const string & tag) {};

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j) {};

  protected:

    /** resets the data if necessary */
    virtual void reset() const {if (this->need_reset()) {ProtectedMappedAtomQuantity<T>::reset(); this->quantity_ = constant_;}};

    /** constant to set data to */
    T constant_;

  private:

    // do not define
    ConstantQuantityMapped();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SummedAtomicQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------
  // Had to define all functions in header, not sure why (JAT, 12/14/11)
  template <typename T>
  class SummedAtomicQuantity : public ProtectedAtomQuantity<T> {

  public:

    // constructor
    SummedAtomicQuantity(ATC_Method * atcTransfer,
                         PerAtomQuantity<T> * quantity1,
                         PerAtomQuantity<T> * quantity2) :
      ProtectedAtomQuantity<T>(atcTransfer, quantity1->nCols(), quantity1->atom_type()),
      quantity1_(quantity1),
      quantity2_(quantity2)
      {
        if (quantity1_->nCols() != quantity2_->nCols())
          throw ATC_Error("SummedAtomicQuantity - dependencies do not have same number of columns");
        if (quantity1_->atom_type() != quantity2_->atom_type())
          throw ATC_Error("SummedAtomicQuantity - dependencies do not have same atom type");
        quantity1_->register_dependence(this);
        quantity2_->register_dependence(this);
      };
    
    // destructor
    virtual ~SummedAtomicQuantity()
      {
        quantity1_->remove_dependence(this);
        quantity2_->remove_dependence(this);
      };

  protected:

    /** resets the data if necessary */
    virtual void reset() const
      {
        if (this->need_reset()) {
          PerAtomQuantity<T>::reset();
          const DenseMatrix<T> & quantity1(quantity1_->quantity());
          const DenseMatrix<T> & quantity2(quantity2_->quantity());
          DenseMatrix<T> & myQuantity(this->quantity_);
          myQuantity = quantity1;
          myQuantity += quantity2;
        }
      };

    /** first quantity */
    PerAtomQuantity<T> * quantity1_;

    /** second quantity */
    PerAtomQuantity<T> * quantity2_;

  private:

    // do not define
    SummedAtomicQuantity();

  };

  /**
   *  @class  PerAtomDiagonalMatrix 
   *  @brief  Base class for objects that manage atomic diagonal matrices and their AtC interface
   */
  template <typename T>
  class PerAtomDiagonalMatrix : public MatrixDependencyManager<DiagonalMatrix, T> {

  public:
    
    // constructor
    PerAtomDiagonalMatrix(ATC_Method * atc, AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~PerAtomDiagonalMatrix();

    /** access to a constant dense matrix of the quantity, indexed by prescribed counts */
    virtual const DiagonalMatrix<T> & quantity() const
      {reset(); return MatrixDependencyManager<DiagonalMatrix, T>::quantity();};

    /** access to a non-constant dens matrix of the quantity, indexed by prescribed atom counts */
    virtual DiagonalMatrix<T> & set_quantity()
      {reset(); return MatrixDependencyManager<DiagonalMatrix, T>::set_quantity();}

    /** sets the Lammps quantity to a given value, input is indexed by AtC atom counts */
    virtual void operator=(const DiagonalMatrix<T> & target) {PerAtomDiagonalMatrix<T>::reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator=(target);};

    /** adds the given data to the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator+=(const DiagonalMatrix<T> & addition) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator+=(addition);};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator+=(addition);};

    /** subtracts the given data from the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator-=(const DiagonalMatrix<T> & subtraction) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator-=(subtraction);};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(double subtraction) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator-=(subtraction);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DiagonalMatrix<T> & multiplier) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator*=(multiplier);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator*=(multiplier);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DiagonalMatrix<T> & divisor) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator/=(divisor);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor) {this->reset(); MatrixDependencyManager<DiagonalMatrix, T>::operator/=(divisor);};

    /** sets quantity to lammps data, if needed, should be called in pre_exchange */
    virtual void prepare_exchange() {this->set_lammps_to_quantity();};

    /** resets AtC local quantity after exchange */
    virtual void post_exchange() {(this->quantity_).resize(atc_.nlocal()); this->set_quantity_to_lammps();};
    
    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return 1;};

    /** packs up data for parallel transfer when atoms change processors */
    virtual int pack_exchange(int i, double *buffer);

    /** unpacks data after parallel transfer when atoms change processors */
    virtual int unpack_exchange(int i, double *buffer);

    // pack/unpack_comm only valid if the quantity is over all real and processor ghost atoms
    /** packs up data for parallel transfer to ghost atoms on other processors */
    virtual int pack_comm(int index, double *buf,
                          int pbc_flag, int *pbc);

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    virtual int unpack_comm(int index, double *buf);

    /** returns per-atom size of communicated data */
    virtual int size_comm() const {return 1;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const string & tag);

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j);

    /** access type of data this quantity is applied to */
    AtomType atom_type() const {return atomType_;};

    /** specialized reset to account for quantities which lammps can change */
    virtual void lammps_force_reset() {};

    /** resets local storage */
    virtual void reset_nlocal() { this->force_reset();}

  protected:
    
    /** utility object to access ATC methods */
    PaqAtcUtility atc_;

    /** pointer to access Lammps data */
    LammpsInterface * lammpsInterface_;

    /** type of atoms this quantity applies to */
    AtomType atomType_;

    /** map from this quantity's AtC indexing to Lammps indexing for atomic arrays */
    const Array<int> & quantityToLammps_;

    /** resets data, if necessary */
    virtual void reset() const {if (this->needReset_) {(this->quantity_).resize(atc_.nlocal()); this->needReset_ = false;}};

    /** sets the quantity based on a lammps pointer */
    virtual void set_lammps_to_quantity() const;

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const;

    
    /** gets appropriate data for lammps pointer */
    virtual T * lammps_scalar() const = 0;

    /** point to lammps-style array for data */
    T * lammpsScalar_;

  private:
    
    // do not define
    PerAtomDiagonalMatrix();
    
  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtcAtomDiagonalMatrix
  //    A funcational base class for defining objects that 
  //    manage diagonal matrices defined at atoms and their  
  //    AtC interface that are defined by AtC classes
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class AtcAtomDiagonalMatrix : public PerAtomDiagonalMatrix<T> {

  public:

    // constructor
    AtcAtomDiagonalMatrix(ATC_Method * atc,
                          AtomType atomType = INTERNAL) : PerAtomDiagonalMatrix<T>(atc,atomType) {};

    // destructor
    virtual ~AtcAtomDiagonalMatrix() {};

  protected:

    // these get the appropriate pointer to local data mimicing lammps storage
    /** gets appropriate data for lammps pointer */
    virtual T * lammps_scalar() const {return this->lammpsScalar_;};

  private:

    // do not define
    AtcAtomDiagonalMatrix();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedAtomDiagonalMatrix
  //    A base class for defining objects that manage
  //    diagonal matrixes defined at atoms internally and 
  //    do not allow for reset externally
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ProtectedAtomDiagonalMatrix : public AtcAtomDiagonalMatrix<T> {

  public:

    // constructor
    ProtectedAtomDiagonalMatrix(ATC_Method * atc,
                                AtomType atomType = INTERNAL)
      : AtcAtomDiagonalMatrix<T>(atc, atomType) {};

    // destructor
    virtual ~ProtectedAtomDiagonalMatrix() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual DiagonalMatrix<T> & set_quantity()
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::set_quantity - Cannot modify protected per atom quantities"); return this->quantity_;};

    /** sets the quantity to a given value */
    virtual void operator=(const DiagonalMatrix<T> & target)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::set_quantity - Cannot modify protected per atom quantities");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator= - Cannot modify protected per atom quantities");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DiagonalMatrix<T> & addition)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator+= - Cannot modify protected per atom quantities");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator+= - Cannot modify protected per atom quantities");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DiagonalMatrix<T> & subtraction)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator-= - Cannot modify protected per atom quantities");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator-= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DiagonalMatrix<T> & multiplier)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DiagonalMatrix<T> & divisor)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator/= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("ProtectedAtomDiagonalMatrix::operator/= - Cannot modify protected per atom quantities");};

  protected:

    /** resets the data if necessary */
    virtual void reset() const = 0;

    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const {this->reset(); PerAtomDiagonalMatrix<T>::set_lammps_to_quantity();};

  private:

    // do not define
    ProtectedAtomDiagonalMatrix();

  };

  /**
   *  @class  PerAtomSparseMatrix
   *  @brief  Base class for objects that manage atomic sparse matrices and their AtC interface
   */
  template <typename T>
  class PerAtomSparseMatrix : public MatrixDependencyManager<SparseMatrix, T> {

  public:
    
    // constructor
    PerAtomSparseMatrix(ATC_Method * atc, int nCols = 1, int maxEntriesPerRow = 1, AtomType atomType = INTERNAL);
    
    // destructor
    virtual ~PerAtomSparseMatrix();

    /** access to a constant dense matrix of the quantity, indexed by prescribed counts */
    virtual const SparseMatrix<T> & quantity() const
      {reset(); return MatrixDependencyManager<SparseMatrix, T>::quantity();};

    /** access to a non-constant dens matrix of the quantity, indexed by prescribed atom counts */
    virtual SparseMatrix<T> & set_quantity()
      {reset(); return MatrixDependencyManager<SparseMatrix, T>::set_quantity();}

    /** number of columns in quantity */
    INDEX nCols() const {return nCols_;};

    /** maximum number of entries per row */
    INDEX max_entries_per_row() const {return maxEntriesPerRow_;};

    /** sets the Lammps quantity to a given value, input is indexed by AtC atom counts */
    virtual void operator=(const SparseMatrix<T> & target) {PerAtomSparseMatrix<T>::reset(); MatrixDependencyManager<SparseMatrix, T>::operator=(target);};

    /** adds the given data to the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator+=(const SparseMatrix<T> & addition) {reset(); MatrixDependencyManager<SparseMatrix, T>::operator+=(addition);};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition) {reset(); MatrixDependencyManager<SparseMatrix, T>::operator+=(addition);};

    /** subtracts the given data from the Lammps quantity, input is indexed by AtC atom counts */
    virtual void operator-=(const SparseMatrix<T> & subtraction) {reset(); MatrixDependencyManager<SparseMatrix, T>::operator-=(subtraction);};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(double subtraction) {reset(); MatrixDependencyManager<SparseMatrix, T>::operator-=(subtraction);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier) {reset(); MatrixDependencyManager<SparseMatrix, T>::operator*=(multiplier);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const SparseMatrix<T> & divisor) {reset(); MatrixDependencyManager<SparseMatrix, T>::operator/=(divisor);};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor) {reset(); MatrixDependencyManager<SparseMatrix, T>::operator/=(divisor);};

    /** sets quantity to lammps data, if needed, should be called in pre_exchange */
    virtual void prepare_exchange() {this->set_lammps_to_quantity();};

    /** resets AtC local quantity after exchange */
    virtual void post_exchange() {(this->quantity_).reset(atc_.nlocal(),nCols_); this->set_quantity_to_lammps();};
    
    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return 2*maxEntriesPerRow_;};

    /** packs up data for parallel transfer when atoms change processors */
    virtual int pack_exchange(int i, double *buffer);

    /** unpacks data after parallel transfer when atoms change processors */
    virtual int unpack_exchange(int i, double *buffer);

    // pack/unpack_comm only valid if the quantity is over all real and processor ghost atoms
    /** packs up data for parallel transfer to ghost atoms on other processors */
    virtual int pack_comm(int index, double *buf,
                          int pbc_flag, int *pbc);

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    virtual int unpack_comm(int index, double *buf);

    /** returns per-atom size of communicated data */
    virtual int size_comm() const {return 2*maxEntriesPerRow_;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const string & tag);

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j);

    /** access type of data this quantity is applied to */
    AtomType atom_type() const {return atomType_;};

    /** specialized reset to account for quantities which lammps can change */
    virtual void lammps_force_reset() {};

    /** resets local storage */
    //WIP_JAT revert this to force_reset when reset_nlocal is functioning correctly (see comment by ATC_Method::reset_nlocal - interscaleManager_.reset_nlocal())
    virtual void reset_nlocal() { this->set_reset();}

  protected:
    
    /** utility object to access ATC methods */
    PaqAtcUtility atc_;

    /** pointer to access Lammps data */
    LammpsInterface * lammpsInterface_;

    /** type of atoms this quantity applies to */
    AtomType atomType_;

    /** number of columns of the per atom quantity, must be defined in derived class */
    int nCols_;

    /** maximum number of entries that can be stored in a row */
    int maxEntriesPerRow_;

    /** map from this quantity's AtC indexing to Lammps indexing for atomic arrays */
    const Array<int> & quantityToLammps_;

    /** resets data, if necessary */
    virtual void reset() const {if (this->needReset_) {(this->quantity_).reset(atc_.nlocal(),nCols_); this->needReset_ = false;}};

    /** sets the quantity based on a lammps pointer */
    virtual void set_lammps_to_quantity() const;

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const;

    
    /** gets appropriate data for lammps pointer */
    virtual T ** lammps_vector() const = 0;

    /** gets appropriate data for lammps pointer to column indices */
    virtual int ** lammps_column_indices() const = 0;

    /** pointer to lammps-style double array for data */
    T ** lammpsVector_;

    /** stores column indices in lammps-style array */
    int ** lammpsColIndices_;

    // workspace
    mutable DenseVector<T> _values_;
    mutable DenseVector<INDEX> _colIndices_;

  private:
    
    // do not define
    PerAtomSparseMatrix();
    
  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtcAtomSparseMatrix
  //    A funcational base class for defining objects that 
  //    manage sparse matrices defined at atoms and their  
  //    AtC interface that are defined by AtC classes
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class AtcAtomSparseMatrix : public PerAtomSparseMatrix<T> {

  public:

    // constructor
    AtcAtomSparseMatrix(ATC_Method * atc,
                        int nCols = 1, int maxEntriesPerRow = 1,
                        AtomType atomType = INTERNAL) : 
                  PerAtomSparseMatrix<T>(atc,nCols,maxEntriesPerRow,atomType) {};

    // destructor
    virtual ~AtcAtomSparseMatrix() {};

  protected:

    // these get the appropriate pointer to local data mimicing lammps storage
    /** gets appropriate data for lammps pointer */
    virtual T ** lammps_vector() const {return this->lammpsVector_;};
    /** gets appropriate data for lammps pointer to column indices */
    virtual int ** lammps_column_indices() const {return this->lammpsColIndices_;};

  private:

    // do not define
    AtcAtomSparseMatrix();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedAtomSparseMatrix
  //    A base class for defining objects that manage
  //    sparse matrixes defined at atoms internally and 
  //    do not allow for reset externally
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ProtectedAtomSparseMatrix : public AtcAtomSparseMatrix<T> {

  public:

    // constructor
    ProtectedAtomSparseMatrix(ATC_Method * atc,
                              int nCols = 1, int maxEntriesPerRow = 1,
                              AtomType atomType = INTERNAL)
      : AtcAtomSparseMatrix<T>(atc, nCols, maxEntriesPerRow, atomType) {};

    // destructor
    virtual ~ProtectedAtomSparseMatrix() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual SparseMatrix<T> & set_quantity()
      {throw ATC_Error("ProtectedAtomSparseMatrix::set_quantity - Cannot modify protected per atom quantities"); return this->quantity_;};

    /** sets the quantity to a given value */
    virtual void operator=(const SparseMatrix<T> & target)
      {throw ATC_Error("ProtectedAtomSparseMatrix::set_quantity - Cannot modify protected per atom quantities");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const T & target)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator= - Cannot modify protected per atom quantities");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const SparseMatrix<T> & addition)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator+= - Cannot modify protected per atom quantities");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(T addition)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator+= - Cannot modify protected per atom quantities");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const SparseMatrix<T> & subtraction)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator-= - Cannot modify protected per atom quantities");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(T subtraction)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator-= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const SparseMatrix<T> & multiplier)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(T multiplier)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator*= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const SparseMatrix<T> & divisor)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator/= - Cannot modify protected per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(T divisor)
      {throw ATC_Error("ProtectedAtomSparseMatrix::operator/= - Cannot modify protected per atom quantities");};

  protected:

    /** resets the data if necessary */
    virtual void reset() const = 0;

    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const {this->reset(); PerAtomSparseMatrix<T>::set_lammps_to_quantity();};

  private:

    // do not define
    ProtectedAtomSparseMatrix();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedMappedSparseMatrix
  //    A base class for defining objects that manage
  //    sparse matrices defined at atoms internally and do 
  //    not allow for reset externally, but are mapped in
  //    at least one of their dimensions
  //--------------------------------------------------------
  //--------------------------------------------------------
  template <typename T>
  class ProtectedMappedAtomSparseMatrix : public ProtectedAtomSparseMatrix<T> {

  public:

    // constructor
    ProtectedMappedAtomSparseMatrix(ATC_Method * atc,
                                    int nCols = 1,
                                    int maxEntriesPerRow = 1,
                                    AtomType atomType = INTERNAL)
      : ProtectedAtomSparseMatrix<T>(atc, nCols, maxEntriesPerRow, atomType)
        {};

    // destructor
    virtual ~ProtectedMappedAtomSparseMatrix() {};

    /** sets quantity to lammps data, if needed, should be called in pre_exchange */
    virtual void prepare_exchange() {};

    /** resets AtC local quantity after exchange */
    virtual void post_exchange() {this->set_reset();};
    
    /** returns how much lammps memory is used in this function */
    virtual int memory_usage() const {return 0;};

    /** packs up data for parallel transfer when atoms change processors */
    virtual int pack_exchange(int i, double *buffer) {return 0;};

    /** unpacks data after parallel transfer when atoms change processors */
    virtual int unpack_exchange(int i, double *buffer) {return 0;};

    // pack/unpack_comm only valid if the quantity is over all real and processor ghost atoms
    /** packs up data for parallel transfer to ghost atoms on other processors */
    virtual int pack_comm(int index, double *buf,
                          int pbc_flag, int *pbc) {return 0;};

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    virtual int unpack_comm(int index, double *buf) {return 0;};

    /** returns per-atom size of communicated data */
    virtual int size_comm() const {return 0;};

    /** changes size of temperary lammps storage data if transfer is being used */
    virtual void grow_lammps_array(int nmax, const string & tag) {};

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    virtual void copy_lammps_array(int i, int j) {};

  protected:

    // methods
    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const {};

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const {};

    /** gets appropriate data for lammps pointer */
    virtual T ** lammps_vector() const {return NULL;};

    /** gets appropriate data for lammps pointer to column indices */
    virtual int ** lammps_column_indices() const {return NULL;};

  private:

    // do not define
    ProtectedMappedAtomSparseMatrix();

  };
}

#include "PerAtomQuantity-inl.h"
#endif

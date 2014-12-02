// A class for defining an alternate atomic mass, needed for templating

#ifndef FUNDAMENTAL_ATOM_QUANTITY_H
#define FUNDAMENTAL_ATOM_QUANTITY_H

#include <string>

#include "PerAtomQuantity.h"

namespace ATC {

  // forward declarations
  class ATC_Method;

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FundamentalAtomQuantity
  //    A base class for defining objects that manage
  //    quantities defined at atoms based on lammps data
  //--------------------------------------------------------
  //--------------------------------------------------------

  class FundamentalAtomQuantity : public ShallowAtomQuantity<double> {

  public:

    // constructor
    FundamentalAtomQuantity(ATC_Method * atc,
                            LammpsInterface::FundamentalAtomQuantity atomQuantity,
                            AtomType atomType=INTERNAL);

    // destructor
    virtual ~FundamentalAtomQuantity() {};

    /** specialized reset to account for quantities which lammps can change */
    virtual void lammps_force_reset() {this->force_reset();};

  protected:

    /** enumerates the type of atom quantity being considered */
    LammpsInterface::FundamentalAtomQuantity atomQuantity_;

    /** converts from Lammps units to ATC units */
    double unitsConversion_;

    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const;

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const
      {ShallowAtomQuantity<double>::set_quantity_to_lammps(); if (unitsConversion_!=1.) quantity_ *= unitsConversion_;};

    /** gets appropriate pointer for lammps data */
    virtual double * lammps_scalar() const
      {return lammpsInterface_->atom_scalar(atomQuantity_);};

    /** gets appropriate pointer for lammps data */
    virtual double ** lammps_vector() const
      {return lammpsInterface_->atom_vector(atomQuantity_);};

  private:

    // do not define
    FundamentalAtomQuantity();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomMass
  //--------------------------------------------------------
  //--------------------------------------------------------

  class AtomMass : public FundamentalAtomQuantity {

  public:

    // constructor
    AtomMass(ATC_Method * atc,AtomType atomType = INTERNAL);

    // destructor
    virtual ~AtomMass() {};

    /** sets the quantity to a given value */
    virtual void operator=(const DENS_MAT & target)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** sets the quantity to a given constant value */
    virtual void operator=(const double & target)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DENS_MAT & addition)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(double addition)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DENS_MAT & subtraction)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(double subtracts)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DENS_MAT & multiplier)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(double multiplier)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DENS_MAT & divisor)
      {throw ATC_Error("Cannot modify type-based atom mass");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(double divisor)
      {throw ATC_Error("Cannot modify type-based atom mass");};

  protected:

    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const {};

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const;

    /** gets appropriate pointer for lammps data */
    virtual double * lammps_scalar() const
      {return NULL;};

    /** gets appropriate pointer for lammps data */
    virtual double ** lammps_vector() const
      {return NULL;};

  private:

    // do not define
    AtomMass();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ComputedAtomQuantity
  //    A base class for defining objects that manage
  //    quantities defined at atoms by Lammps computes
  //    The compute associated with the tag must already
  //    be initialized.
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ComputedAtomQuantity : public ShallowAtomQuantity<double> {

  public:

    // constructor
    ComputedAtomQuantity(ATC_Method * atc,
                         const std::string & tag,
                         double unitsConversion = 1.,
                         AtomType atomType=INTERNAL);

    // destructor
    virtual ~ComputedAtomQuantity() {};

    /** resets compute, must be this way to accomodate atom sorting between runs */
    virtual void post_exchange() {this->needReset_ = true;};

    /** specialized reset to account for forcing lammps to perform the compute */
    virtual void force_reset();

    /** specialized reset to account for quantities which lammps can change */
    virtual void lammps_force_reset() {this->force_reset();};

    // remove operations that change the lammps data
    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual DENS_MAT & set_quantity()
      {throw ATC_Error("ComputedAtomQuantity::set_quantity - Cannot modify computed per atom quantities"); return quantity_;};

    /** sets the quantity to a given constant value */
    virtual void operator=(const DENS_MAT & target)
      {throw ATC_Error("ComputedAtomQuantity::operator= - Cannot modify computed per atom quantities");};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const DENS_MAT & addition)
      {throw ATC_Error("ComputedAtomQuantity::operator+= - Cannot modify computed per atom quantities");};

    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(double addition)
      {throw ATC_Error("ComputedAtomQuantity::operator+= - Cannot modify computed per atom quantities");};

    /** subtracts the given data from the Lammps quantity */
    virtual void operator-=(const DENS_MAT & subtraction)
      {throw ATC_Error("ComputedAtomQuantity::operator-= - Cannot modify computed per atom quantities");};

    /** subtracts the scalar data from the Lammps quantity for AtC atoms */
    virtual void operator-=(double subtraction)
      {throw ATC_Error("ComputedAtomQuantity::operator-= - Cannot modify computed per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const DENS_MAT & multiplier)
      {throw ATC_Error("ComputedAtomQuantity::operator*= - Cannot modify computed per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(double multiplier)
      {throw ATC_Error("ComputedAtomQuantity::operator*= - Cannot modify computed per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const DENS_MAT & divisor)
      {throw ATC_Error("ComputedAtomQuantity::operator/= - Cannot modify computed per atom quantities");};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(double divisor)
      {throw ATC_Error("ComputedAtomQuantity::operator/= - Cannot modify computed per atom quantities");};

  protected:

    /** pointer to Lammps compute, meant as rapid indexing only (do not use!) */
    COMPUTE_POINTER computePointer_;  

    /** tag for Lammps compute */
    std::string computeTag_;

    /** units conversion from LAMMPS to ATC units */
    double unitsConversion_;

    /** sets the quantity based on a lammps pointer */
    virtual void set_quantity_to_lammps() const
    {ShallowAtomQuantity<double>::set_quantity_to_lammps(); if (unitsConversion_!=1.) quantity_ *= unitsConversion_;};

    /** gets appropriate data for lammps pointer */
    virtual double * lammps_scalar() const {return lammpsInterface_->compute_vector_peratom(computePointer_);};

    /** gets appropriate data for lammps pointer */
    virtual double ** lammps_vector() const {return lammpsInterface_->compute_array_peratom(computePointer_);};

    // not needed if no MPI
    /** sets lammps data based on the quantity */
    virtual void set_lammps_to_quantity() const
      {throw ATC_Error("ComputedAtomQuantity::set_lammps_to_quantity - Cannot modify a compute's LAMMPS data");};

  private:

    // do not define
    ComputedAtomQuantity();

  };

}
#endif

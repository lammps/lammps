#include "FundamentalAtomicQuantity.h"

using std::string;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FundamentalAtomQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  FundamentalAtomQuantity::FundamentalAtomQuantity(ATC_Method * atc,
                                                   LammpsInterface::FundamentalAtomQuantity atomQuantity,
                                                   AtomType atomType) :
    ShallowAtomQuantity<double>(atc,0,atomType),
    atomQuantity_(atomQuantity),
    unitsConversion_(lammpsInterface_->atom_quantity_conversion(atomQuantity_))
  {
    nCols_ = lammpsInterface_->atom_quantity_ndof(atomQuantity_);
  }

  //--------------------------------------------------------
  //  set_lammps_to_quantity
  //--------------------------------------------------------
  void FundamentalAtomQuantity::set_lammps_to_quantity() const
  {
    if (unitsConversion_==1.) { // is there a  way to avoid equal testing a double?
      PerAtomQuantity<double>::set_lammps_to_quantity();
    }
    else { // perform unit conversion
      if (quantity_.nRows()>0) {
        // full matrix copy
        if (atomType_ == ALL || atomType_ == PROC_GHOST) {
          if (nCols_==1) { // scalar
            double * lammpsQuantity = lammps_scalar();

            for (int i = 0; i < atc_.nlocal_total(); i++)
              lammpsQuantity[i] = quantity_(i,0)/unitsConversion_;
          }
          else{ // vector
            double ** lammpsQuantity = lammps_vector();

            for (int i = 0; i < atc_.nlocal_total(); i++)
              for (int j = 0; j < nCols_; j++)
                lammpsQuantity[i][j] = quantity_(i,j)/unitsConversion_;
          }
        }
        // mapped copy
        else {
          int atomIndex;

          if (nCols_==1) { // scalar
            double * lammpsQuantity = lammps_scalar();
            for (int i = 0; i < quantity_.nRows(); i++) {
              atomIndex = quantityToLammps_(i);
              lammpsQuantity[atomIndex] = quantity_(i,0)/unitsConversion_;
            }
          }
          else{ // vector
            double ** lammpsQuantity = lammps_vector();
            for (int i = 0; i < quantity_.nRows(); i++) {
              atomIndex = quantityToLammps_(i);
              for (int j = 0; j < nCols_; j++) {
                lammpsQuantity[atomIndex][j] = quantity_(i,j)/unitsConversion_;
              }
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomMass
  //    Access-only operations when mass is
  //    defined per type.
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomMass::AtomMass(ATC_Method * atc,AtomType atomType) :
    FundamentalAtomQuantity(atc,LammpsInterface::ATOM_MASS,atomType)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  set_quantity_to_lammps
  //--------------------------------------------------------
  void AtomMass::set_quantity_to_lammps() const
  {
    const int * type = lammpsInterface_->atom_type();
    const double * mass = lammpsInterface_->atom_mass();

    if (atomType_ == ALL || atomType_ == PROC_GHOST) {
      for (int i = 0; i < quantity_.nRows(); i++)
        quantity_(i,0) = mass[type[i]];
    }
    else {
      int atomIndex;
      for (int i = 0; i < quantity_.nRows(); i++) {
        atomIndex = quantityToLammps_(i);
        quantity_(i,0) = mass[type[atomIndex]];
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ComputedAtomQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  ComputedAtomQuantity::ComputedAtomQuantity(ATC_Method * atc,
                                             const string & tag,
                                             double unitsConversion,
                                             AtomType atomType) :
    ShallowAtomQuantity<double>(atc,0,atomType),
    computePointer_(nullptr),
    computeTag_(tag),
    unitsConversion_(unitsConversion)
  {
    // register compute with lammps interface and provide pointer for syncing
    computePointer_ = lammpsInterface_->compute_pointer(computeTag_.c_str());
    nCols_ = lammpsInterface_->compute_ncols_peratom(computePointer_);
  }

  //--------------------------------------------------------
  //  force_reset
  //--------------------------------------------------------
  void ComputedAtomQuantity::force_reset()
  {
    // only reset if the compute needs it this timestep
    if (lammpsInterface_->compute_matchstep(computePointer_,lammpsInterface_->ntimestep())) {
      if (!isFixed_) {
        lammpsInterface_->reset_invoked_flag(computePointer_);
      }
      ShallowAtomQuantity<double>::force_reset();
    }
  }

}

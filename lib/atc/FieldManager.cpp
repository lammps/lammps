#include "FieldManager.h"
#include "ATC_Method.h"
#include "LammpsInterface.h"
#include "PerAtomQuantity.h"
#include "TransferOperator.h"


namespace ATC {

typedef PerAtomQuantity<double> PAQ;

//-----------------------------------------------------------------------------
//*
//-----------------------------------------------------------------------------
  FieldManager::FieldManager(ATC_Method * atc):
    atc_(atc),
    interscaleManager_(atc->interscale_manager())
  {};


//-----------------------------------------------------------------------------
//* restricted_atom_quantity
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::restricted_atom_quantity(FieldName field, string name, PAQ * atomicQuantity) 
  {
    if (name == "default") { name = field_to_restriction_name(field); }
    DENS_MAN * quantity = interscaleManager_.dense_matrix(name);

    if (!quantity){
      if      (field == CHARGE_DENSITY) {
        atomicQuantity = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE);
      }
      else if (field == MASS_DENSITY) {
        atomicQuantity = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);
      }
      else {
        
        if (!atomicQuantity) {
          throw ATC_Error("FieldManager::restricted_atom_quantity - need to supply PAQ if restricted quantity does not already exist");
        }
      }
      quantity = new AtfShapeFunctionRestriction(atc_,atomicQuantity,atc_->accumulant());
      interscaleManager_.add_dense_matrix(quantity,name);
    }
    return quantity;
  }

//-----------------------------------------------------------------------------
//* restricted_atom_quantity
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::projected_atom_quantity(FieldName field,string name, PAQ * atomic, FieldName massMat, DIAG_MAN * normalization)
  {
    if (atc_->use_md_mass_normalization()) {
      if (name == "default") { name = field_to_intrinsic_name(field); }
      DENS_MAN * quantity = interscaleManager_.dense_matrix(name);
      if (!quantity) {
        DENS_MAN * restricted = restricted_atom_quantity(field,field_to_restriction_name(field),atomic);
        quantity = new AtfShapeFunctionMdProjection(atc_,restricted,massMat);
        interscaleManager_.add_dense_matrix(quantity,name);
      }
      return quantity;
    }
    else {
      if (name == "default") { name = field_to_string(field); }
      DENS_MAN * quantity = interscaleManager_.dense_matrix(name);
      if (quantity) return quantity;

      if (atc_->kernel_on_the_fly()) {
        if (atc_->kernel_based()) {
          quantity = new OnTheFlyKernelAccumulationNormalized(atc_, atomic, 
                                                              atc_->kernel_function(), 
                                                              atc_->atom_coarsegraining_positions(), 
                                                              normalization);
        } else {
          quantity = new OnTheFlyMeshAccumulationNormalized(atc_, atomic, 
                                                            atc_->atom_coarsegraining_positions(), 
                                                            normalization);
        }
      } else {
        quantity = new AtfProjection(atc_, atomic, 
                                     atc_->accumulant(), 
                                     normalization);
      }
      interscaleManager_.add_dense_matrix(quantity,name);
      return quantity;
    }
  }

//-----------------------------------------------------------------------------
//* referenced_projected_atom_quantity
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::referenced_projected_atom_quantity(FieldName field,string name, PAQ * atomic, const DENS_MAT * reference, FieldName massMat, DIAG_MAN * normalization)
  {
    if (name == "default") { name = field_to_string(field); }
    DENS_MAN * quantity = interscaleManager_.dense_matrix(name);
    if (quantity) return quantity;

    if (atc_->use_md_mass_normalization()) {
      DENS_MAN * restricted = restricted_atom_quantity(field,field_to_restriction_name(field),atomic);
      quantity = new AtfShapeFunctionMdProjectionReferenced(atc_,restricted,reference,massMat);
    }
    else if (atc_->kernel_on_the_fly()) {
      if (atc_->kernel_based()) {
        quantity = new OnTheFlyKernelAccumulationNormalizedReferenced(atc_, 
          atomic, 
          atc_->kernel_function(), 
          atc_->atom_coarsegraining_positions(), 
          normalization,
          reference);
      } else {
        quantity = new OnTheFlyMeshAccumulationNormalizedReferenced(atc_, 
          atomic, 
          atc_->atom_coarsegraining_positions(), 
          normalization,
          reference);
      }
    } else {
      quantity = new AtfProjectionReferenced(atc_, atomic, 
          atc_->accumulant(),
          reference,
          normalization);
    }
    interscaleManager_.add_dense_matrix(quantity,name);
    return quantity;
  }

//-----------------------------------------------------------------------------
//* scaled_projected_atom_quantity
//-----------------------------------------------------------------------------

  DENS_MAN * FieldManager::scaled_projected_atom_quantity(FieldName field,string name, PAQ * atomic, double scale, FieldName massMat, DIAG_MAN * normalization)
  {
    if (name == "default") { name = field_to_string(field); }
    DENS_MAN * quantity = interscaleManager_.dense_matrix(name);
    if (quantity) return quantity;

    if (atc_->use_md_mass_normalization()) {
      DENS_MAN * restricted = restricted_atom_quantity(field,field_to_restriction_name(field),atomic);
      quantity = new AtfShapeFunctionMdProjectionScaled(atc_,restricted,scale,massMat);
    }
    else if (atc_->kernel_on_the_fly()) {
      if (atc_->kernel_based()) {
        quantity = new OnTheFlyKernelAccumulationNormalizedScaled(atc_, atomic, 
          atc_->kernel_function(), 
          atc_->atom_coarsegraining_positions(), 
          normalization,
          scale);
      } else {
        quantity = new OnTheFlyMeshAccumulationNormalizedScaled(atc_, atomic, 
          atc_->atom_coarsegraining_positions(), 
          normalization,
          scale);
      }
    } else {
      quantity = new AtfProjectionScaled(atc_, atomic,
          atc_->accumulant(),
          scale,
          normalization);
    }
    interscaleManager_.add_dense_matrix(quantity,name);
    return quantity;
  }

//-----------------------------------------------------------------------------
//* CHARGE_DENSITY
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::charge_density(string name)
  {
    FundamentalAtomQuantity * atomic = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE);
    return projected_atom_quantity(CHARGE_DENSITY,name,atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* MASS_DENSITY
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::mass_density(string name)
  {
    FundamentalAtomQuantity * atomic = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);
    return projected_atom_quantity(MASS_DENSITY,name,atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* SPECIES_CONCENTRATION
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::species_concentration(string name)
  {
    DENS_MAN * c = NULL;
#ifdef ATC_VERBOSE
    atc_->print_tracked();
#endif
    PAQ * atomSpecies = atomic_species_vector(); 

    if (atc_->kernel_on_the_fly()) {
      if (atc_->kernel_based()) {
        c = new OnTheFlyKernelAccumulationNormalized(atc_, atomSpecies, 
          atc_->kernel_function(), 
          atc_->atom_coarsegraining_positions(), 
          atc_->accumulant_inverse_volumes());
      } else {
        c = new OnTheFlyMeshAccumulationNormalized(atc_, atomSpecies, 
          atc_->atom_coarsegraining_positions(), 
          atc_->accumulant_inverse_volumes());
      }
    } else {
      c = new AtfProjection(atc_, atomSpecies, 
        atc_->accumulant(), 
        atc_->accumulant_inverse_volumes());
    }
    if (name == "default") { name = field_to_string(SPECIES_CONCENTRATION); }
    interscaleManager_.add_dense_matrix(c,name);
    return c;
  }
//-----------------------------------------------------------------------------
//* NUMBER_DENSITY
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::number_density(string name)
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("atomNumber");
    if (!atomic) {
      atomic = new AtomNumber(atc_);
      interscaleManager_.add_per_atom_quantity(atomic, "atomNumber");
    }
    return projected_atom_quantity(NUMBER_DENSITY,name,atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* MOMENTUM 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::momentum(string name)
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("atomMomentum");
    if (!atomic) {
      atomic = new AtomicMomentum(atc_);
      interscaleManager_.add_per_atom_quantity(atomic, "atomMomentum");
    }
    return projected_atom_quantity(MOMENTUM,name,atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* VELOCITY 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::velocity(string name)
  {
    if (name == "default") { name = field_to_string(VELOCITY); }
    DENS_MAN * v = interscaleManager_.dense_matrix(name);
    if (v) return v;

    if (atc_->use_md_mass_normalization()) {
      PAQ * atomic = interscaleManager_.per_atom_quantity("atomMomentum");
      if (!atomic) {
        atomic = new AtomicMomentum(atc_);
        interscaleManager_.add_per_atom_quantity(atomic, "atomMomentum");
      }
      DENS_MAN * restricted = restricted_atom_quantity(VELOCITY,field_to_restriction_name(VELOCITY),atomic);
      v = new AtfShapeFunctionMdProjection(atc_,restricted,VELOCITY);
    }
    else {
      DENS_MAN* p = interscaleManager_.dense_matrix(field_to_string(MOMENTUM));
      if (!p) p = nodal_atomic_field(MOMENTUM);
      DENS_MAN* m = interscaleManager_.dense_matrix(field_to_string(MASS_DENSITY));
      if (!m) m = nodal_atomic_field(MASS_DENSITY);
      v = new DenseMatrixQuotient(p,m);
    }
    interscaleManager_.add_dense_matrix(v,field_to_string(VELOCITY));
    return v;
  }
//-----------------------------------------------------------------------------
//* PROJECTED_VELOCITY 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::projected_velocity(string name)
  {
    FundamentalAtomQuantity * atomic = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);
    return projected_atom_quantity(PROJECTED_VELOCITY,name,atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* DISPLACEMENT 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::displacement(string name)
  {
    if (name == "default") { name = field_to_string(DISPLACEMENT); }
    DENS_MAN * u = interscaleManager_.dense_matrix(name);
    if (u) return u;

    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicMassWeightedDisplacement");
    if (!atomic) {
      FundamentalAtomQuantity * atomMasses = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);    
      FundamentalAtomQuantity * atomPositions = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_POSITION);    
      atomic = new AtomicMassWeightedDisplacement(atc_,atomPositions, atomMasses, atc_->atom_reference_positions(), INTERNAL);
      interscaleManager_.add_per_atom_quantity(atomic,"AtomicMassWeightedDisplacement");
    }
    if (atc_->use_md_mass_normalization()) {
      DENS_MAN * restricted = restricted_atom_quantity(DISPLACEMENT,field_to_restriction_name(DISPLACEMENT),atomic);
      u = new AtfShapeFunctionMdProjection(atc_,restricted,VELOCITY);
    }
    else {
      DENS_MAN * q = NULL;
      if (atc_->kernel_on_the_fly()) {
        if (atc_->kernel_based()) {
          q = new OnTheFlyKernelAccumulationNormalized(atc_, atomic,
            atc_->kernel_function(),
            atc_->atom_coarsegraining_positions(),
            atc_->accumulant_inverse_volumes());
        } else {
          q = new OnTheFlyMeshAccumulationNormalized(atc_, atomic,
            atc_->atom_coarsegraining_positions(),
            atc_->accumulant_inverse_volumes());
        }
      } else {
        q = new AtfProjection(atc_, atomic,
                              atc_->accumulant(),
                              atc_->accumulant_inverse_volumes());
      }
      interscaleManager_.add_dense_matrix(q,"CoarseGrainedAMWD");
      DENS_MAN* m  = interscaleManager_.dense_matrix(field_to_string(MASS_DENSITY));
      u = new DenseMatrixQuotient(q,m);
    }
    interscaleManager_.add_dense_matrix(u,name);
    return u;
  }
//-----------------------------------------------------------------------------
//* REFERENCE_POTENTIAL_ENERGY 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::reference_potential_energy(string name)
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicReferencePotential");
    if (!atomic) {
      atomic = new AtcAtomQuantity<double>(atc_);
      interscaleManager_.add_per_atom_quantity(atomic, "AtomicReferencePotential");
      atomic->set_memory_type(PERSISTENT);
    }
    return projected_atom_quantity(REFERENCE_POTENTIAL_ENERGY,"NodalAtomicReferencePotential",atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* POTENTIAL_ENERGY 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::potential_energy(string name)
  {
    PerAtomQuantity<double> * atomic = interscaleManager_.per_atom_quantity("AtomicPotentialEnergy");
    const DENS_MAT * reference = atc_->nodal_ref_potential_energy();
    return referenced_projected_atom_quantity(POTENTIAL_ENERGY,name,atomic,reference,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* TEMPERATURE 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::temperature(string name)
 {
    double Tcoef = 1./((atc_->nsd())*(atc_->lammps_interface())->kBoltzmann());
    PAQ * atomic = per_atom_quantity("AtomicTwiceFluctuatingKineticEnergy");
    return scaled_projected_atom_quantity(TEMPERATURE,name,atomic,Tcoef,TEMPERATURE,atc_->accumulant_weights());
  }
//-----------------------------------------------------------------------------
//* KINETIC_TEMPERATURE
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::kinetic_temperature(string name)
  {
    double Tcoef = 1./((atc_->nsd())*(atc_->lammps_interface())->kBoltzmann());
    PAQ * atomic = per_atom_quantity("AtomicTwiceKineticEnergy");
    return scaled_projected_atom_quantity(KINETIC_TEMPERATURE,name,atomic,Tcoef,MASS_DENSITY,atc_->accumulant_weights());  
  }
//-----------------------------------------------------------------------------
//* THERMAL_ENERGY 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::thermal_energy(string name)
  {
    double Ecoef = 0.5*atc_->ke_scale();
    PAQ * atomic = per_atom_quantity("AtomicTwiceFluctuatingKineticEnergy");
    return scaled_projected_atom_quantity(THERMAL_ENERGY,name,atomic,Ecoef,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* KINETIC_ENERGY 
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::kinetic_energy(string name)
  {
    double Ecoef = 0.5*atc_->ke_scale();
    PAQ * atomic = per_atom_quantity("AtomicTwiceKineticEnergy");
    return scaled_projected_atom_quantity(KINETIC_ENERGY,name,atomic,Ecoef,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* CHARGE_FLUX
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::charge_flux(string name)
  {

    PAQ * atomic = per_atom_quantity("AtomicChargeVelocity");
    return projected_atom_quantity(CHARGE_FLUX,name,atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }
//-----------------------------------------------------------------------------
//* SPECIES_FLUX
//-----------------------------------------------------------------------------
  DENS_MAN * FieldManager::species_flux(string name)
  {
    PAQ * atomic = per_atom_quantity("AtomicSpeciesVelocity");
    return projected_atom_quantity(SPECIES_FLUX,name,atomic,MASS_DENSITY,atc_->accumulant_inverse_volumes());
  }



//=============================================================================
//* PER ATOM QUANTITIES
//=============================================================================

//-----------------------------------------------------------------------------
//* 2 KE '
//-----------------------------------------------------------------------------
  PAQ * FieldManager::atomic_twice_fluctuating_kinetic_energy()
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicTwiceFluctuatingKineticEnergy");
    if (!atomic) {
      FundamentalAtomQuantity * atomMass = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);
      FundamentalAtomQuantity * atomVelocity = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);
      PAQ * vbar = per_atom_quantity(field_to_prolongation_name(VELOCITY));
      atomic = new TwiceFluctuatingKineticEnergy(atc_,atomVelocity,atomMass,vbar);
      interscaleManager_.add_per_atom_quantity(atomic, "AtomicTwiceFluctuatingKineticEnergy");
    }
    return atomic;
  }
//-----------------------------------------------------------------------------
//* 2 KE 
//-----------------------------------------------------------------------------
  PAQ * FieldManager::atomic_twice_kinetic_energy()
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicTwiceKineticEnergy");
    if (!atomic) {
      FundamentalAtomQuantity * atomMass = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_MASS);
      FundamentalAtomQuantity * atomVelocity = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);
      atomic = new TwiceKineticEnergy(atc_,atomVelocity,atomMass);
      interscaleManager_.add_per_atom_quantity(atomic, "AtomicTwiceKineticEnergy");
    }
    return atomic;
  }
//-----------------------------------------------------------------------------
//* v'
//-----------------------------------------------------------------------------
  PAQ * FieldManager::atomic_fluctuating_velocity()
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicFluctuatingVelocity");
    if (!atomic) {
      FundamentalAtomQuantity * atomVelocity = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY);
      PAQ * atomMeanVelocity = per_atom_quantity(field_to_prolongation_name(VELOCITY));
      atomic = new FluctuatingVelocity(atc_,atomVelocity,atomMeanVelocity);
      interscaleManager_.add_per_atom_quantity(atomic, "AtomicFluctuatingVelocity");
    }
    return atomic;
  }
//-----------------------------------------------------------------------------
//* q v'
//-----------------------------------------------------------------------------
  PAQ * FieldManager::atomic_charge_velocity()
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicChargeVelocity");
    if (!atomic) {
      PAQ * atomVelocity = atomic_fluctuating_velocity();
      FundamentalAtomQuantity * atomCharge = interscaleManager_.fundamental_atom_quantity(LammpsInterface::ATOM_CHARGE);
      atomic = new ChargeVelocity(atc_,atomVelocity,atomCharge);
      interscaleManager_.add_per_atom_quantity(atomic, "AtomicChargeVelocity");
    }
    return atomic;
  }
//-----------------------------------------------------------------------------
//* m^a v'
//-----------------------------------------------------------------------------
  PAQ * FieldManager::atomic_species_velocity()
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicSpeciesVelocity");
    if (!atomic) {
      PAQ * atomVelocity = atomic_fluctuating_velocity();
      PAQ * atomSpecies = atomic_species_vector(); 
      atomic = new SpeciesVelocity(atc_,atomVelocity,atomSpecies);
      interscaleManager_.add_per_atom_quantity(atomic, "AtomicSpeciesVelocity");
    }
    return atomic;
  }

//-----------------------------------------------------------------------------
//* [0 1 0 0 ] for type 2 atom
//-----------------------------------------------------------------------------
  PAQ * FieldManager::atomic_species_vector()
  {
    PAQ * atomic = interscaleManager_.per_atom_quantity("AtomicSpeciesVector");
    if (!atomic) {
      atomic = new AtomTypeVector(atc_,atc_->type_list(),atc_->group_list());
      interscaleManager_.add_per_atom_quantity(atomic,"AtomicSpeciesVector");
    }
    return atomic;
  }

   PAQ * FieldManager::prolonged_field(FieldName field)
   {
     PAQ * quantity = interscaleManager_.per_atom_quantity(field_to_prolongation_name(field));
     if (!quantity) {
       
      DENS_MAN * coarseQuantity = interscaleManager_.dense_matrix(field_to_string(field));
       if (!coarseQuantity) coarseQuantity = nodal_atomic_field(field);
       if (!coarseQuantity) throw ATC_Error("can not prolong quantity: " + field_to_string(field) + " no field registered");
       if (atc_->kernel_on_the_fly()) {
         quantity = new OnTheFlyShapeFunctionProlongation(atc_,
                   coarseQuantity,atc_->atom_coarsegraining_positions());
       } else {
         quantity = new FtaShapeFunctionProlongation(atc_,
                  coarseQuantity,atc_->interpolant());
       }
       interscaleManager_.add_per_atom_quantity(quantity,
                     field_to_prolongation_name(field));
     }
     return quantity;
   }

}

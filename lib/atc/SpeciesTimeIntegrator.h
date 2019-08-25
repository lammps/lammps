#ifndef SPECIES_TIME_INTEGRATOR_H
#define SPECIES_TIME_INTEGRATOR_H

#include <map>
#include <utility>
#include <string>

// ATC headers
#include "TimeIntegrator.h"

namespace ATC {

  // forward declarations
  class ATC_CouplingMass;
  class SpeciesIntegrationMethod;
  
  /**
   *  @class  SpeciesTimeIntegrator
   *  @brief  Class for various time integrators for species FE quantities in the Eulerian frame
   *          (handles parsing and stores basic data structures)
   */
  
  class SpeciesTimeIntegrator : public TimeIntegrator {
    
  public:
    
    // constructor
    SpeciesTimeIntegrator(ATC_CouplingMass * atc,
                          TimeIntegrationType timeIntegrationType);
    
    // destructor
    virtual ~SpeciesTimeIntegrator(){};
    
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);
    
    /** create objects to implement requested numerical method */
    virtual void construct_methods();
    
    /** pack persistent fields */
    virtual void pack_fields(RESTART_LIST & data);
    
    DENS_MAN & nodal_atomic_species_concentration_filtered() { return  nodalAtomicSpeciesConcentrationFiltered_; }

  protected:
        
    /** sets of molecules tracked */
    const std::map<std::string,std::pair<MolSize,int> > & moleculeIds_;

    /** data */
    DENS_MAN nodalAtomicSpeciesConcentrationFiltered_;
    
  private:
    
    // DO NOT define this
    SpeciesTimeIntegrator();
    
  };

  /**
   *  @class  SpeciesIntegrationMethod
   *  @brief  Class for species time integration methods which update FE quantities in time
   */

  class SpeciesIntegrationMethod : public TimeIntegrationMethod {
  
  public:
  
    // constructor
    SpeciesIntegrationMethod(SpeciesTimeIntegrator * speciesTimeIntegrator,
      const std::map<std::string,std::pair<MolSize,int> > & moleculeIds);
        
    // destructor
    virtual ~SpeciesIntegrationMethod() {nodalAtomicMassDensity_=NULL;};

    /** create and get necessary transfer operators */
    virtual void construct_transfers();
        
  protected:
    /** time filtering object */
    TimeFilter * timeFilter_;

    /** finite element mass density field */
    DENS_MAN & massDensity_;
      
    /** atomic nodal mass density field */
// OBSOLETE?
    DENS_MAN & nodalAtomicMassDensityOut_;
    DENS_MAN * nodalAtomicMassDensity_;
      
    /** finite element mass density field */
    DENS_MAN & speciesConcentration_;
    DENS_MAN * nodalAtomicSpeciesConcentration_;
    DENS_MAN & nodalAtomicSpeciesConcentrationFiltered_;

    /** sets of molecules tracked */
    const std::map<std::string,std::pair<MolSize,int> > & moleculeIds_;

  private:

    // DO NOT define this
    SpeciesIntegrationMethod();
  
  };

  /**
   *  @class  SpeciesTimeIntegratorFractionalStep
   *  @brief  FractionalStep integration for FE species quantities
   *          (Uses 2nd order FractionalStep integration to update
   *           the FE mass density field)
   */

  class SpeciesTimeIntegratorFractionalStep : public SpeciesIntegrationMethod {
  
  public:
  
    // constructor
    SpeciesTimeIntegratorFractionalStep(SpeciesTimeIntegrator * speciesTimeIntegrator,
      const std::map<std::string,std::pair<MolSize,int> > & moleculeIds);
        
    // destructor
    virtual ~SpeciesTimeIntegratorFractionalStep(){};

    /** pre time integration */
    virtual void initialize();
        
    // time step methods, corresponding to ATC_Transfer
    virtual void pre_initial_integrate1(double dt);
    virtual void pre_final_integrate1(double dt);
    virtual void post_final_integrate2(double dt);

    /** post processing step before output */
    virtual void post_process();
        
  private:
  
    // DO NOT define this
    SpeciesTimeIntegratorFractionalStep();
  
  };
  /**
   *  @class  SpeciesTimeIntegratorFractionalStepFiltered
   *  @brief  FractionalStep integration for FE species quantities
   *          (Uses 2nd order FractionalStep integration to update
   *           the FE mass density field)
   */

  class SpeciesTimeIntegratorFractionalStepFiltered : public SpeciesTimeIntegratorFractionalStep {
  
  public:
  
    // constructor
    SpeciesTimeIntegratorFractionalStepFiltered(SpeciesTimeIntegrator * speciesTimeIntegrator,
      const std::map<std::string,std::pair<MolSize,int> > & moleculeIds);
        
    // destructor
    virtual ~SpeciesTimeIntegratorFractionalStepFiltered(){};

    /** pre time integration */
    virtual void initialize() {};
        
    // time step methods, corresponding to ATC_Transfer
    /** first part of pre_final_integrate */
    virtual void pre_final_integrate1(double dt);

    /** post processing step before output */
    virtual void post_process(){};

  private:
  
    // DO NOT define this
    SpeciesTimeIntegratorFractionalStepFiltered();
  
  };
};

#endif

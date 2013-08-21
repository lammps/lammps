#ifndef ATC_COUPLING_MASS_H
#define ATC_COUPLING_MASS_H

/** owned field/s: MASS_DENSITY */

// ATC headers
#include "ATC_Coupling.h"

// Other headers
#include <map>
#include <string>

namespace ATC {

  // Forward declarations
  class FE_Engine;
  class SpeciesTimeIntegrator;
  class ChargeRegulator;
  class ConcentrationRegulator;

  /**
   *  @class ATC_CouplingMass
   *  @brief A class for atom-continuum transfers & control for species transport
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ATC_CouplingMass
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ATC_CouplingMass : public ATC_Coupling {

  public:
  
    // constructor
    ATC_CouplingMass(std::string groupName, 
                     double **& perAtomArray,
                     LAMMPS_NS::Fix * thisFix,
                     std::string matParamFile,
                     ExtrinsicModelType extrinsic = NO_MODEL);
      
    // destructor
    virtual ~ATC_CouplingMass();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);
  
    /** pre time integration */
    virtual void initialize();

    /** prior to force */
    virtual void pre_force();
    /** prior to exchange */
    virtual void pre_exchange();
    virtual void reset_atoms() { resetNlocal_=true;}
#ifdef OBSOLETE

    /** first time, after atomic velocity but before position integration */
    virtual void mid_init_integrate();
    /** first time, after atomic integration */
    virtual void post_init_integrate();
#endif
    /** second time, after atomic integration */
    virtual void post_final_integrate();

    /** compute vector for output */
    virtual double compute_vector(int n);

    /** output routines */
    virtual void output();

  protected:

    // functions
    //---------------------------------------------------------------
    /** initialization routines */
    //---------------------------------------------------------------
    /** constructs all data which is updated with time integration, i.e. fields */
    //virtual void construct_time_integration_data();
    /** create methods, e.g. time integrators, filters */
    virtual void construct_methods();
    /** set up data which is dependency managed */
    virtual void construct_transfers();

    /** sets the position/velocity of the ghost atoms */
    virtual void set_ghost_atoms(){};

    /** compute the mass matrix components coming from MD integration */
    virtual void compute_md_mass_matrix(FieldName thisField,
                                        DIAG_MAT & massMats);

    /** physics specific filter initialization */
    void init_filter();

    /** field mask for velocity integration */
    Array2D<bool> speciesMask_;

    // Add in fields for restarting
    virtual void  read_restart_data(std::string fileName_, RESTART_LIST & data);
    virtual void write_restart_data(std::string fileName_, RESTART_LIST & data);

    // DATA structures for tracking individual species and molecules
    FIELD_POINTERS atomicFields_;

    bool resetNlocal_;

    
    
    //      i.e. we only need the correct shape function matrix for restriction
  };

};

#endif

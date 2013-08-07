#ifndef ATC_COUPLING_ENERGY_H
#define ATC_COUPLING_ENERGY_H

// ATC headers
#include "ATC_Coupling.h"
#include "ThermalTimeIntegrator.h"
//TEMP_JAT - remove when needs_reset is moved to base class
#include "AtomicRegulator.h"

// Other headers
#include <map>

namespace ATC {

  class AtfShapeFunctionRestriction;
  class AtfShapeFunctionMdProjection;

  /**
   *  @class ATC_CouplingEnergy
   *  @brief A class for atom-continuum transfers & control involving heat transport
   *         (owned field/s: TEMPERATURE)
   */

  class ATC_CouplingEnergy : public ATC_Coupling {

  public:
  
    // constructor
    ATC_CouplingEnergy(std::string groupName, 
                       double ** & perAtomArray,
                       LAMMPS_NS::Fix * thisFix,
                       std::string matParamFile,
                       ExtrinsicModelType extrinsic = NO_MODEL);

    // destructor
    virtual ~ATC_CouplingEnergy();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();

    /** flags whether a methods reset is required */
    
    virtual bool reset_methods() const {
      bool resetMethods = ATC_Method::reset_methods() || atomicRegulator_->need_reset();
      _ctiIt_ = timeIntegrators_.find(TEMPERATURE);
      if (_ctiIt_ == timeIntegrators_.end()) return resetMethods;
      return resetMethods || (_ctiIt_->second)->need_reset();
    };

    /** post time integration */
    virtual void finish();

    /** first time substep routines */
    virtual void pre_init_integrate();
    /** first time, after atomic velocity but before position integration */
    virtual void mid_init_integrate();
    /** first time, after atomic integration */
    virtual void post_init_integrate();

    /** second time substep routine */
    virtual void post_final_integrate();

    /** compute vector for output */
    virtual double compute_vector(int n);

    /** output */
    virtual void output();

    
    /** set up atom to material identification */
    virtual void reset_atom_materials();

  protected:

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
    
    /** adds resetting of any thermostat arrays associated with local atom count */
    virtual void reset_nlocal();

    /** compute the mass matrix components coming from MD integration */
    virtual void compute_md_mass_matrix(FieldName thisField,
                                        DIAG_MAT & massMats);

    /** operator to compute mass matrix from MD */
    AtfShapeFunctionRestriction * nodalAtomicHeatCapacity_;

    /** physics specific filter initialization */
    void init_filter();

    /** field mask for velocity integration */
    Array2D<bool> temperatureMask_;

    
    double compute_lambda_power(int gid);

    /** kinetic temperature for post-processing */
    AtfShapeFunctionMdProjection * nodalAtomicKineticTemperature_;

    /** configurational temperature for post-processing */
    AtfShapeFunctionMdProjection * nodalAtomicConfigurationalTemperature_;

    /** workspace matrices for output */
    DENS_MAT _keTemp_, _peTemp_;

    // Add in fields for restarting
    virtual void  read_restart_data(string fileName_, RESTART_LIST & data);
    virtual void write_restart_data(string fileName_, RESTART_LIST & data);
    void pack_thermal_fields(RESTART_LIST & data);

  };
    
};
#endif

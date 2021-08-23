#ifndef ATC_COUPLING_MOMENTUM_ENERGY_H
#define ATC_COUPLING_MOMENTUM_ENERGY_H

// ATC headers
#include "ATC_Coupling.h"
#include "Kinetostat.h"
#include "Thermostat.h"
#include "ElasticTimeIntegrator.h"
#include "ThermalTimeIntegrator.h"

// Other headers
#include <string>

namespace ATC {

  class AtfShapeFunctionRestriction;
  class AtfShapeFunctionMdProjection;

  /**
   *  @class ATC_CouplingMomentumEnergy
   *  @brief A class for atom-continuum transfers & control involving momentum and heat transport
   *         (owned field/s: DISPLACEMENT, VELOCITY, TEMPERATURE)
   */

  class ATC_CouplingMomentumEnergy : public ATC_Coupling {

  public:

    // constructor
    ATC_CouplingMomentumEnergy(std::string groupName,
                               double ** & perAtomArray,
                               LAMMPS_NS::Fix * thisFix,
                               std::string matParamFile,
                               ExtrinsicModelType extrinsic = NO_MODEL);

    // destructor
    virtual ~ATC_CouplingMomentumEnergy();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();

    /** compute scalar for output - added energy */
    virtual double compute_scalar(void);

    /** compute vector for output */
    virtual double compute_vector(int n);
    double kinetic_energy();
    double potential_energy();

    /** output */
    virtual void output();

  protected:

    //---------------------------------------------------------------
    /** initialization routines */
    //---------------------------------------------------------------
    /** constructs all data which is updated with time integration, i.e. fields */
    //virtual void construct_time_integration_data();
    /** set up data which is dependency managed */
    virtual void construct_transfers();

    /** physics specific filter initialization */
    void init_filter();
    /** kinetic temperature for post-processing */
    AtfShapeFunctionMdProjection * nodalAtomicKineticTemperature_;

    /** configurational temperature for post-processing */
    AtfShapeFunctionMdProjection * nodalAtomicConfigurationalTemperature_;

    /** workspace matrices for output */
    DENS_MAT _keTemp_, _peTemp_;


    // data
    double refPE_;
  };

};
#endif

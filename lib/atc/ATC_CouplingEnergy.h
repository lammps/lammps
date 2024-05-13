#ifndef ATC_COUPLING_ENERGY_H
#define ATC_COUPLING_ENERGY_H

// ATC headers
#include "ATC_Coupling.h"
#include "ThermalTimeIntegrator.h"
//TEMP_JAT - remove when needs_reset is moved to base class
#include "AtomicRegulator.h"

// Other headers
#include <map>
#include <string>

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

    /** compute vector for output */
    virtual double compute_vector(int n);

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

    /** sets the position/velocity of the ghost atoms */
    virtual void set_ghost_atoms(){};

    /** physics specific filter initialization */
    void init_filter();


    double compute_lambda_power(int gid);

    /** kinetic temperature for post-processing */
    AtfShapeFunctionMdProjection * nodalAtomicKineticTemperature_;

    /** configurational temperature for post-processing */
    AtfShapeFunctionMdProjection * nodalAtomicConfigurationalTemperature_;

    /** workspace matrices for output */
    DENS_MAT _keTemp_, _peTemp_;

  };

};
#endif

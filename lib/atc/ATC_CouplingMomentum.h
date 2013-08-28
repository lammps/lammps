#ifndef ATC_COUPLING_MOMENTUM_H
#define ATC_COUPLING_MOMENTUM_H

/** owned field/s: DISPLACEMENT, VELOCITY */

// ATC headers
#include "ATC_Coupling.h"
#include "Kinetostat.h"
#include "ElasticTimeIntegrator.h"

// Other headers
#include <string>

namespace ATC {

  // Forward declarations
  class AtfShapeFunctionRestriction;
  class AtfShapeFunctionMdProjection;

  /**
   *  @class ATC_CouplingMomentum
   *  @brief A class for atom-continuum transfers & control involving mechanical motion
    *         (owned field/s: DISPLACEMENT, VELOCITY)
   */

  class ATC_CouplingMomentum : public ATC_Coupling {

  public:
  
    // constructor
    ATC_CouplingMomentum(std::string groupName, 
                         double **& perAtomArray,
                         LAMMPS_NS::Fix * thisFix,
                         std::string matParamFile,
                         PhysicsType intrinsicModel,
                         ExtrinsicModelType extrinsicModel = NO_MODEL);
      
    // destructor
    virtual ~ATC_CouplingMomentum();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);
  
    /** pre time integration */
    virtual void initialize();

    /** pre/post atomic force calculation in minimize */
    virtual void min_pre_force();
    virtual void min_post_force();

    /** compute scalar for output - added energy */
    virtual double compute_scalar(void);

    /** compute vector for output */
    virtual double compute_vector(int n);
    double kinetic_energy(const IntegrationDomainType domain=FULL_DOMAIN); // const;
    double potential_energy(const IntegrationDomainType domain=FULL_DOMAIN) const;

    /** output routines */
    virtual void output(void);

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

    // data
    double refPE_;

    /** mass matrix computed with atomic quadrature for KE output */
    MASS_MATS Ma_;
  };

};

#endif

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

    /** flags whether a methods reset is required */
    virtual bool reset_methods() const {
      bool resetMethods = ATC_Method::reset_methods() ||atomicRegulator_->need_reset();
      _ctiIt_ = timeIntegrators_.find(VELOCITY);
      if (_ctiIt_ == timeIntegrators_.end()) return resetMethods;
      return resetMethods || (_ctiIt_->second)->need_reset();
    };

    /** post time integration */
    virtual void finish();
#ifdef OBSOLETE
    

    /** first time, after atomic velocity but before position integration */
    virtual void mid_init_integrate();
    /** first time, after atomic integration */
    virtual void post_init_integrate();
#endif
    /** second time, before atomic integration */
    virtual void pre_final_integrate();
    /** second time, after atomic integration */
    virtual void post_final_integrate();

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

    /** adds resetting of any kinetostat arrays associated with local atom count */
    virtual void reset_nlocal();

    /** compute the mass matrix components coming from MD integration */
    virtual void compute_md_mass_matrix(FieldName thisField,
                                        DIAG_MAT & massMats);

    /** operator to compute the mass matrix for the momentum equation from MD integration */
    AtfShapeFunctionRestriction * nodalAtomicMass_;

    /** operator to compute the dimensionless mass matrix from MD integration */
    AtfShapeFunctionRestriction * nodalAtomicCount_;

    /** physics specific filter initialization */
    void init_filter();

    /** field mask for velocity integration */
    Array2D<bool> velocityMask_;

    // Add in fields for restarting
    virtual void  read_restart_data(std::string fileName_, RESTART_LIST & data);
    virtual void write_restart_data(std::string fileName_, RESTART_LIST & data);
    void pack_elastic_fields(RESTART_LIST & data);

    // data
    double refPE_;

    /** mass matrix computed with atomic quadrature for KE output */
    MASS_MATS Ma_;
  };

};

#endif

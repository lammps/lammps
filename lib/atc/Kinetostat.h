#ifndef KINETOSTAT_H
#define KINETOSTAT_H

// ATC_Transfer headers
#include "AtomicRegulator.h"

// other headers
#include <map>
#include <set>

namespace ATC {

  /**
   *  @class  Kinetostat
   *  @brief  Manager class for atom-continuum control of momentum and position
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class Kinetostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  class Kinetostat : public AtomicRegulator {
  
  public:

    /** kinetostat types */
    enum KinetostatType {
      NONE=0,
      GLC_DISPLACEMENT,
      GLC_VELOCITY,
      FORCE
    };

    enum KinetostatCouplingType {
      UNCOUPLED=0,
      FLUX,
      FIXED
    };
  
    // constructor
    Kinetostat(ATC_Transfer *atcTransfer);
        
    // destructor
    ~Kinetostat(){};
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();

    // data access, intended for method objects
    /** reset the nodal force to a prescribed value */
    void reset_lambda_force(DENS_MAT & target);
    /** return the nodal force induced by lambda */
    DENS_MAT & get_nodal_atomic_lambda_force() { return nodalAtomicLambdaForce_;}
    /** return value of filtered lambda */
    DENS_MAT & get_lambda_force_filtered() { return lambdaForceFiltered_;}
    /** access to kinetostat type */
    KinetostatType get_kinetostat_type() const
      { return kinetostatType_;};
    KinetostatCouplingType get_coupling_mode() const
      { return couplingMode_;};
        
  protected:

    /** kinetostat type flag */
    KinetostatType kinetostatType_;
    /** kinetostat copuling type flag */
    KinetostatCouplingType couplingMode_;

    // kinetostat data
    /** lambda force applied to atoms */
    DENS_MAT nodalAtomicLambdaForce_;
    /** filtered lambda force */
    DENS_MAT lambdaForceFiltered_;
        
  private:

    // DO NOT define this
    Kinetostat();

  };
  
  /**
   *  @class  GlcKinetostat
   *  @brief  Base class for implementation of kinetostat algorithms based on Gaussian least constraints (GLC)
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GlcKinetostat
  //    base class for all thermostats of general form of a
  //    Gaussian least constraint (GLC)
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class GlcKinetostat : public RegulatorShapeFunction {
  
  public:
  
    GlcKinetostat(Kinetostat *kinetostat);
        
    ~GlcKinetostat(){};

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();

  protected:

    // methods
    /** apply forces to atoms */
    virtual void apply_to_atoms(double ** atomicQuantity,
                                const DENS_MAT & lambdaAtom,
                                double dt=0.);
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights(DIAG_MAT & weights);
    /** apply any required corrections for localized kinetostats */
    virtual void apply_localization_correction(const DENS_MAT & source,
                                               DENS_MAT & nodalField,
                                               double weight = 1.){};

    // member data
    /** pointer to thermostat object for data */
    Kinetostat * kinetostat_;
    /** pointer to a time filtering object */
    TimeFilter * timeFilter_;
    /** stress induced by lambda */
    DENS_MAT & nodalAtomicLambdaForce_;
    /** filtered lambda force */
    DENS_MAT & lambdaForceFiltered_;
    /** atomic force induced by lambda */
    DENS_MAT & lambdaForce_;
    /** MD mass matrix */
    MATRIX & mdMassMatrix_;
    /** mass of ATC internal atoms on this processor */
    DENS_VEC atomicMass_;
    /** reference to ATC map from global nodes to overlap nodes */
    Array<int> & nodeToOverlapMap_;
    /** nodeset corresponding to Hoover coupling */
    set<pair<int,int> > hooverNodes_;

  private:
    
    // DO NOT define this
    GlcKinetostat();

  };
  
  /**
   *  @class  DisplacementGlc
   *  @brief  Enforces GLC on atomic position based on FE displacement
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DisplacementGlc
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class DisplacementGlc : public GlcKinetostat {
  
  public:
  
    DisplacementGlc(Kinetostat * kinetostat);
        
    ~DisplacementGlc(){};
        
    /** applies kinetostat to atoms */
    virtual void apply_post_predictor(double dt);

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);
        
  protected:
        
    // methods
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights(DIAG_MAT & weights);
    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double dt);
    /** compute force induced by lambda */
    virtual void compute_lambda_force(DENS_MAT & lambdaForce, DENS_MAT & lambdaAtom, double dt);
    /** sets up and solves kinetostat equations */
    virtual void compute_kinetostat(double dt);
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double dt);
    /** apply any required corrections for localized kinetostats */
    virtual void apply_localization_correction(const DENS_MAT & source,
                                               DENS_MAT & nodalField,
                                               double weight = 1.);

    // data
    /** clone of FE displacement field */
    DENS_MAT & nodalDisplacements_;
    /** pointer to lammps atomic positions */
    double ** x_;

  private:
    
    // DO NOT define this
    DisplacementGlc();
  
  };

  /**
   *  @class  DisplacementGlcFiltered
   *  @brief  Enforces GLC on time filtered atomic position based on FE displacement
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DisplacementGlcFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class DisplacementGlcFiltered : public DisplacementGlc {
  
  public:
  
    DisplacementGlcFiltered(Kinetostat * kinetostat);
        
    ~DisplacementGlcFiltered(){};

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);
        
  protected:
        
    // methods
    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double dt);
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double dt);

    // data
    /** clone of FE nodal atomic displacement field */
    DENS_MAT & nodalAtomicDisplacements_;

  private:
    
    // DO NOT define this
    DisplacementGlcFiltered();
  
  };

  /**
   *  @class  VelocityGlc
   *  @brief  Enforces GLC on atomic velocity based on FE velocity
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class VelocityGlc
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class VelocityGlc : public GlcKinetostat {
  
  public:
  
    VelocityGlc(Kinetostat * kinetostat);
        
    ~VelocityGlc(){};
        
    /** applies kinetostat to atoms */
    virtual void apply_mid_predictor(double dt);

    /** applies kinetostat to atoms */
    virtual void apply_post_corrector(double dt);

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);
        
  protected:
  
    // methods
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights(DIAG_MAT & weights);
    /** compute force induced by lambda */
    virtual void compute_lambda_force(DENS_MAT & lambdaForce, DENS_MAT & lambdaAtom, double dt);
    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double dt);
    /** sets up and solves kinetostat equations */
    virtual void compute_kinetostat(double dt);
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double dt);
    /** apply any required corrections for localized kinetostats */
    virtual void apply_localization_correction(const DENS_MAT & source,
                                               DENS_MAT & nodalField,
                                               double weight = 1.);

    // data
    /** clone of FE velocity field */
    DENS_MAT & nodalVelocities_;
    /** pointer to lammps atomic velocities */
    double ** v_;

  private:

    // DO NOT define this
    VelocityGlc();
  
  };

  /**
   *  @class  VelocityGlcFiltered
   *  @brief  Enforces GLC on time filtered atomic velocity based on FE velocity
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class VelocityGlcFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class VelocityGlcFiltered : public VelocityGlc {
  
  public:
  
    VelocityGlcFiltered(Kinetostat * kinetostat);
        
    ~VelocityGlcFiltered(){};

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);
        
  protected:

    // methods
    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double dt);
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double dt);
  
    // data
    /** clone of FE nodal atomic velocity field */
    DENS_MAT & nodalAtomicVelocities_;

  private:

    // DO NOT define this
    VelocityGlcFiltered();
        
  };

  /**
   *  @class  StressFlux
   *  @brief  Enforces GLC on atomic forces based on FE stresses or accelerations
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class StressFlux
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class StressFlux : public GlcKinetostat {
  
  public:
  
    StressFlux(Kinetostat * kinetostat);
        
    ~StressFlux();

    /** applies kinetostat to atoms in the mid-predictor phase */
    virtual void apply_mid_predictor(double dt);

    /** applies kinetostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies kinetostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);

    /** sets filtered ghost force to prescribed value */
    void reset_filtered_ghost_force(DENS_MAT & targetForce);
    /** returns reference to filtered ghost force */
    DENS_MAT & get_filtered_ghost_force() {return nodalGhostForceFiltered_;};
        
  protected:

    // data
    /** nodal force */
    DENS_MAT & nodalForce_;
    /** nodal force due to atoms */
    DENS_MAT & nodalAtomicForce_;
    /** nodal ghost force */
    DENS_MAT nodalGhostForce_;
    /** filtered ghost force */
    DENS_MAT nodalGhostForceFiltered_;
    /** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAT & momentumSource_;
    /** pointer to lammps atomic velocities */
    double ** v_;
    /** pointer to lammps atomic forces */
    double ** f_;
#if false
    /** initial lammps atomic forces */
    double ** f0_;
#endif

    // methods
    /** compute force induced by lambda */
    virtual void compute_lambda_force(DENS_MAT & lambdaForce);
    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double dt);
    /** sets up and solves kinetostat equations */
    virtual void compute_kinetostat(double dt);
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double dt);
    /** apply forces to atoms */
    virtual void apply_to_atoms(double ** atomicVelocity, 
                                const DENS_MAT & lambdaForce,
                                double dt);
    /** adds in finite element rhs contributions */
    virtual void add_to_rhs(FIELDS & rhs);
    /** computes restricted force on ghost atoms */
    void compute_ghost_force(DENS_MAT & nodalGhostForce);

  private:

    // DO NOT define this
    StressFlux();

  };

  /**
   *  @class  StressFluxFiltered
   *  @brief  Enforces GLC on time filtered atomic forces based on FE stresses or accelerations
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class StressFluxFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class StressFluxFiltered : public StressFlux {
  
  public:
  
    StressFluxFiltered(Kinetostat * kinetostat);
        
    ~StressFluxFiltered(){};

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);
        
  protected:

    // data
    DENS_MAT & nodalAtomicVelocity_;

    // methods
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** adds in finite element rhs contributions */
    virtual void add_to_rhs(FIELDS & rhs);
    /** apply forces to atoms */
    virtual void apply_to_atoms(double ** atomicVelocity, 
                                const DENS_MAT & lambdaForce,
                                double dt);

  private:

    // DO NOT define this
    StressFluxFiltered();

  };

};

#endif

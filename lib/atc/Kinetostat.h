#ifndef KINETOSTAT_H
#define KINETOSTAT_H

#include "AtomicRegulator.h"
#include "PerAtomQuantityLibrary.h"
#include <map>
#include <set>
#include <utility>
#include <string>

namespace ATC {

  // forward declarations
  class FundamentalAtomQuantity;
  class AtfShapeFunctionRestriction;
  template <typename T>
    class ProtectedAtomQuantity;

  /**
   *  @class  Kinetostat
   *  @brief  Manager class for atom-continuum control of momentum and position
   */

  class Kinetostat : public AtomicRegulator {
  
  public:

    // constructor
    Kinetostat(ATC_Coupling *atc,
               const std::string & regulatorPrefix = "");
        
    // destructor
    virtual ~Kinetostat(){};
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** instantiate up the desired method(s) */
    virtual void construct_methods();

    // data access, intended for method objects
    /** reset the nodal force to a prescribed value */
    virtual void reset_lambda_contribution(const DENS_MAT & target);

  private:

    // DO NOT define this
    Kinetostat();

  };

  /**
   *  @class  KinetostatShapeFunction
   *  @brief  Base class for implementation of kinetostat algorithms based on FE shape functions
   */
  
  class KinetostatShapeFunction : public RegulatorShapeFunction {
  
  public:
  
    KinetostatShapeFunction(AtomicRegulator *kinetostat,
                            const std::string & regulatorPrefix = "");
        
    virtual ~KinetostatShapeFunction(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

  protected:

    // methods
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights();

    // member data
    /** MD mass matrix */
    DIAG_MAN & mdMassMatrix_;
    /** pointer to a time filtering object */
    TimeFilter * timeFilter_;
    /** stress induced by lambda */
    DENS_MAN * nodalAtomicLambdaForce_;
    /** filtered lambda force */
    DENS_MAN * lambdaForceFiltered_;
    /** atomic force induced by lambda */
    ProtectedAtomQuantity<double> * atomKinetostatForce_;
    /** lambda prolonged to the atoms */
    ProtectedAtomQuantity<double> * atomLambda_;

    /** pointer to atom velocities */
    FundamentalAtomQuantity * atomVelocities_;
    /** pointer to atom velocities */
    FundamentalAtomQuantity * atomMasses_;

    // workspace
    DENS_MAT _nodalAtomicLambdaForceOut_; // matrix for output only

  private:
    
    // DO NOT define this
    KinetostatShapeFunction();

  };
  
  /**
   *  @class  GlcKinetostat
   *  @brief  Base class for implementation of kinetostat algorithms based on Gaussian least constraints (GLC)
   */
  
  class GlcKinetostat : public KinetostatShapeFunction {
  
  public:
  
    GlcKinetostat(AtomicRegulator *kinetostat);
        
    virtual ~GlcKinetostat(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

  protected:

    // methods
    /** apply forces to atoms */
    virtual void apply_to_atoms(PerAtomQuantity<double> * quantity,
                                const DENS_MAT & lambdaAtom,
                                double dt=0.);

    /** apply any required corrections for localized kinetostats */
    virtual void apply_localization_correction(const DENS_MAT & /* source */,
                                               DENS_MAT & /* nodalField */,
                                               double /* weight */){};
    virtual void apply_localization_correction(const DENS_MAT & /* source */,
                                               DENS_MAT & /* nodalField */){};

    // member data
    /** nodeset corresponding to Hoover coupling */
    std::set<std::pair<int,int> > hooverNodes_;

    
    /** pointer to atom positions */
    FundamentalAtomQuantity * atomPositions_;

  private:
    
    // DO NOT define this
    GlcKinetostat();

  };
  
  /**
   *  @class  DisplacementGlc
   *  @brief  Enforces GLC on atomic position based on FE displacement
   */
  
  class DisplacementGlc : public GlcKinetostat {
  
  public:
  
    DisplacementGlc(AtomicRegulator * kinetostat);
        
    virtual ~DisplacementGlc(){};
        
    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies kinetostat to atoms */
    virtual void apply_post_predictor(double dt);

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /** determine if local shape function matrices are needed */
    virtual bool use_local_shape_functions() const {return (!atomicRegulator_->use_lumped_lambda_solve()) && atomicRegulator_->use_localized_lambda();};
        
  protected:
        
    // methods
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights();

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
    /** restricted atomic displacements at the nodes */
    DENS_MAN * nodalAtomicMassWeightedDisplacement_;
    /** clone of FE displacement field */
    DENS_MAN & nodalDisplacements_;

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
  
    DisplacementGlcFiltered(AtomicRegulator * kinetostat);
        
    virtual ~DisplacementGlcFiltered(){};

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
        
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
    DENS_MAN & nodalAtomicDisplacements_;

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
  
    VelocityGlc(AtomicRegulator * kinetostat);
        
    virtual ~VelocityGlc(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();
        
    /** applies kinetostat to atoms */
    virtual void apply_mid_predictor(double dt);

    /** applies kinetostat to atoms */
    virtual void apply_post_corrector(double dt);

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
    
    /** determine if local shape function matrices are needed */
    virtual bool use_local_shape_functions() const {return (!atomicRegulator_->use_lumped_lambda_solve()) && atomicRegulator_->use_localized_lambda();};

  protected:
  
    // methods
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights();

    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double dt);
    /** sets up and solves kinetostat equations */
    virtual void compute_kinetostat(double dt);
    /** applies kinetostat correction to atoms */
    virtual void apply_kinetostat(double dt);
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double dt);
    /** apply any required corrections for localized kinetostats */
    virtual void apply_localization_correction(const DENS_MAT & source,
                                               DENS_MAT & nodalField,
                                               double weight = 1.);

    // data
    /** restricted atomic displacements at the nodes */
    DENS_MAN * nodalAtomicMomentum_;
    /** clone of FE velocity field */
    DENS_MAN & nodalVelocities_;

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
  
    VelocityGlcFiltered(AtomicRegulator * kinetostat);
        
    virtual ~VelocityGlcFiltered(){};

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
        
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
    DENS_MAN & nodalAtomicVelocities_;

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
  
    StressFlux(AtomicRegulator * kinetostat);
        
    virtual ~StressFlux();

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** applies kinetostat to atoms in the pre-predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies kinetostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields);

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /** sets filtered ghost force to prescribed value */
    void reset_filtered_ghost_force(DENS_MAT & targetForce);
    /** returns reference to filtered ghost force */
    DENS_MAN & filtered_ghost_force() {return nodalGhostForceFiltered_;};

    /** determine if local shape function matrices are needed */
    virtual bool use_local_shape_functions() const {return ((!atomicRegulator_->use_lumped_lambda_solve()) && atomicRegulator_->use_localized_lambda());};
        
  protected:

    // data
    /** nodal force */
    DENS_MAN & nodalForce_;
    /** nodal force due to atoms */
    DENS_MAN * nodalAtomicForce_;
    /** nodal ghost force */
    AtfShapeFunctionRestriction * nodalGhostForce_;
    /** filtered ghost force */
    DENS_MAN nodalGhostForceFiltered_;
    /** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAN & momentumSource_;

    // methods
    /** does initial filtering operations before main computation */
    virtual void apply_pre_filtering(double dt);
    /** sets up and solves kinetostat equations */
    virtual void compute_kinetostat(double dt);
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** computes the nodal FE force applied by the kinetostat */
    virtual void compute_nodal_lambda_force(double dt);
    /** apply forces to atoms */
    virtual void apply_to_atoms(PerAtomQuantity<double> * atomVelocities,
                                const DENS_MAT & lambdaForce,
                                double dt);
    /** adds in finite element rhs contributions */
    virtual void add_to_rhs(FIELDS & rhs);

    // workspace
    DENS_MAT _deltaVelocity_; // change in velocity during time integration
  private:

    // DO NOT define this
    StressFlux();

  };

  /**
   *  @class  StressFluxGhost
   *  @brief  Enforces GLC on atomic forces based on FE stresses or accelerations, using
   *          the ghost forces to prescribe the FE boundary stress
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class StressFluxGhost
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class StressFluxGhost : public StressFlux {
  
  public:
  
    StressFluxGhost(AtomicRegulator * kinetostat);
        
    virtual ~StressFluxGhost() {};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** compute boundary flux, requires kinetostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields);
        
  protected:

    // methods
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    /** adds in finite element rhs contributions */
    virtual void add_to_rhs(FIELDS & rhs);

  private:

    // DO NOT define this
    StressFluxGhost();

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
  
    StressFluxFiltered(AtomicRegulator * kinetostat);
        
    virtual ~StressFluxFiltered(){};

    /** adds in finite element rhs contributions */
    virtual void add_to_rhs(FIELDS & rhs);

    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);
        
  protected:

    // data
    DENS_MAN & nodalAtomicVelocity_;

    // methods
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs, double dt);
    
    /** apply forces to atoms */
    virtual void apply_to_atoms(PerAtomQuantity<double> * quantity,
                                const DENS_MAT & lambdaAtom,
                                double dt);

  private:

    // DO NOT define this
    StressFluxFiltered();

  };

  /**
   *  @class  KinetostatGlcFs
   *  @brief  Base class for implementation of kinetostat algorithms based on Gaussian least constraints (GLC)
   *          when fractional step time integration is used
   */
  
  class KinetostatGlcFs : public KinetostatShapeFunction {
  
  public:
  
    KinetostatGlcFs(AtomicRegulator *kinetostat,
                    const std::string & regulatorPrefix = "");
        
    virtual ~KinetostatGlcFs(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);
    
    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /* flag for performing the full lambda prediction calculation */
    bool full_prediction();

  protected:

    // methods
    /** determine mapping from all nodes to those to which the kinetostat applies */
    void compute_rhs_map();

    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs,
                                    double dt) = 0;

    /** apply forces to atoms */
    virtual void apply_to_atoms(PerAtomQuantity<double> * atomicVelocity,
                                const DENS_MAN * nodalAtomicEnergy,
                                const DENS_MAT & lambdaForce,
                                DENS_MAT & nodalAtomicLambdaPower,
                                double dt);

    /** add contributions from kinetostat to FE energy */
    virtual void add_to_momentum(const DENS_MAT & nodalLambdaForce,
                                 DENS_MAT & deltaMomemtum,
                                 double dt) = 0;

    /* sets up and solves the linear system for lambda */
    virtual void compute_lambda(double dt);

    // member data
    /** reference to AtC FE velocity */
    DENS_MAN & velocity_;

    /** nodal atomic momentum */
    DENS_MAN * nodalAtomicMomentum_;

    /** hack to determine if first timestep has been passed */
    bool isFirstTimestep_;

    /** local version of velocity used as predicted final veloctiy */
    PerAtomQuantity<double> * atomPredictedVelocities_;

    /** predicted nodal atomic momentum */
    AtfShapeFunctionRestriction * nodalAtomicPredictedMomentum_;

    /** FE momentum change from kinetostat forces */
    DENS_MAT deltaMomentum_;

    /** right-hand side data for thermostat equation */
    DENS_MAT rhs_;

    /** fraction of timestep over which constraint is exactly enforced */
    double dtFactor_;

    // workspace
    DENS_MAT _lambdaForceOutput_; // force applied by lambda in output format
    DENS_MAT _velocityDelta_; // change in velocity when lambda force is applied

  private:
    
    // DO NOT define this
    KinetostatGlcFs();

  };

  /**
   *  @class  KinetostatFlux
   *  @brief  Implementation of kinetostat algorithms based on Gaussian least constraints (GLC)
   *          which apply stresses when fractional step time integration is used
   */
  
  class KinetostatFlux : public KinetostatGlcFs {
  
  public:
  
    KinetostatFlux(AtomicRegulator *kinetostat,
                   const std::string & regulatorPrefix = "");
        
    virtual ~KinetostatFlux(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** enables resetting of filtered ghost force */
    void reset_filtered_ghost_force(DENS_MAT & target);

  protected:

    // methods
    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs,
                                    double dt);

    /** add contributions from kinetostat to FE energy */
    virtual void add_to_momentum(const DENS_MAT & nodalLambdaForce,
                                 DENS_MAT & deltaMomemtum,
                                 double dt);

    /** sets up the transfer which is the set of nodes being regulated */
    virtual void construct_regulated_nodes();

    // member data
    /** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAN & momentumSource_;
    
    /** force from ghost atoms restricted to nodes */
    DENS_MAN * nodalGhostForce_;

    /** filtered nodal ghost force */
    DENS_MAN * nodalGhostForceFiltered_;

  private:
    
    // DO NOT define this
    KinetostatFlux();

  };

  /**
   *  @class  KinetostatFluxGhost
   *  @brief  Implements ghost-atom boundary flux and other loads for fractional-step based kinetostats
   */
  
  class KinetostatFluxGhost : public KinetostatFlux {
  
  public:
  
    KinetostatFluxGhost(AtomicRegulator *kinetostat,
                        const std::string & regulatorPrefix = "");
        
    virtual ~KinetostatFluxGhost(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** compute boundary flux */
    virtual void compute_boundary_flux(FIELDS & fields);

  protected:

    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs,
                                    double dt);

    /** add contributions from kinetostat to FE energy */
    virtual void add_to_momentum(const DENS_MAT & nodalLambdaForce,
                                 DENS_MAT & deltaMomemtum,
                                 double dt);

  private:
    
    // DO NOT define this
    KinetostatFluxGhost();

  };

    /**
   *  @class  KinetostatFixed
   *  @brief  Implementation of kinetostat algorithms based on Gaussian least constraints (GLC)
   *          which perform Hoover coupling when fractional step time integration is used
   */
  
  class KinetostatFixed : public KinetostatGlcFs {
  
  public:
  
    KinetostatFixed(AtomicRegulator *kinetostat,
                    const std::string & regulatorPrefix = "");
        
    virtual ~KinetostatFixed(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & /* fields */)
      {boundaryFlux_[VELOCITY] = 0.;};

    /** determine if local shape function matrices are needed */
    virtual bool use_local_shape_functions() const {return atomicRegulator_->use_localized_lambda();};

  protected:

    // methods
    /** initialize data for tracking the change in nodal atomic velocity */
    virtual void initialize_delta_nodal_atomic_momentum(double dt);

    /** compute the change in nodal atomic velocity */
    virtual void compute_delta_nodal_atomic_momentum(double dt);

    /** sets up appropriate rhs for kinetostat equations */
    virtual void set_kinetostat_rhs(DENS_MAT & rhs,
                                    double dt);

    /** add contributions from kinetostat to FE energy */
    virtual void add_to_momentum(const DENS_MAT & nodalLambdaForce,
                                 DENS_MAT & deltaMomemtum,
                                 double dt);

    /* sets up and solves the linear system for lambda */
    virtual void compute_lambda(double dt);

    /** flag for halving the applied force to mitigate numerical errors */
    bool halve_force();

    /** sets up the transfer which is the set of nodes being regulated */
    virtual void construct_regulated_nodes();

    // member data
    /** change in FE momentum over a timestep */
    DENS_MAT deltaFeMomentum_;

    /** initial FE momentum used to compute change */
    DENS_MAT initialFeMomentum_;

    /** change in restricted atomic FE momentum over a timestep */
    DENS_MAT deltaNodalAtomicMomentum_;

    /** initial restricted atomic FE momentum used to compute change */
    DENS_MAT initialNodalAtomicMomentum_;

    /** filtered nodal atomic momentum */
    DENS_MAN nodalAtomicMomentumFiltered_;

    /** coefficient to account for effect of time filtering on rhs terms */
    double filterCoefficient_;

    // workspace
    DENS_MAT _tempNodalAtomicMomentumFiltered_; // stores filtered momentum change in atoms for persistence during predictor

  private:
    
    // DO NOT define this
    KinetostatFixed();

  };

  /**
   *  @class  KinetostatFluxFixed
   *  @brief  Class for kinetostatting using the velocity matching constraint one one set of nodes and the flux matching constraint on another
   */

  class KinetostatFluxFixed : public RegulatorMethod {

  public:

    KinetostatFluxFixed(AtomicRegulator * kinetostat,
                        bool constructThermostats = true);
        
    virtual ~KinetostatFluxFixed();

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** pre-run initialization of method data */
    virtual void initialize();

    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);
    
    /** get data for output */
    virtual void output(OUTPUT_LIST & outputData);

    /** compute boundary flux, requires kinetostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields)
      {kinetostatBcs_->compute_boundary_flux(fields);};

  protected:

    // data
    /** kinetostat for imposing the fluxes */
    KinetostatFlux * kinetostatFlux_;

    /** kinetostat for imposing fixed nodes */
    KinetostatFixed * kinetostatFixed_;

    /** pointer to whichever kinetostat should compute the flux, based on coupling method */
    KinetostatGlcFs * kinetostatBcs_;

  private:

    // DO NOT define this
    KinetostatFluxFixed();
  };

}

#endif

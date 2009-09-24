/** Thermostat : a class for atom-continuum control of energy/temperature */

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

// ATC_Transfer headers
#include "AtomicRegulator.h"

// other headers
#include <map>
#include <set>

namespace ATC {

  /**
   *  @class  Thermtostat
   *  @brief  Manager class for atom-continuum control of thermal energy
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class Thermostat
  //--------------------------------------------------------
  //--------------------------------------------------------

  class Thermostat : public AtomicRegulator {
  
  public:

    /** thermostat types */
    enum ThermostatType {
      NONE=0,
      RESCALE,
      HOOVER,
      FLUX
    };
  
    // constructor
    Thermostat(ATC_Transfer * atcTransfer);
        
    // destructor
    ~Thermostat(){};
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();
        
    // data access, intended for method objects
    /** reset the nodal power to a prescribed value */
    void reset_lambda_power(DENS_MAT & target);
    /** return the nodal power induced by lambda */
    DENS_VEC & get_nodal_atomic_lambda_power(){return nodalAtomicLambdaPower_;};
    /** return value filtered lambda  */
    DENS_VEC & get_filtered_lambda_power(){return lambdaPowerFiltered_;};
    /** access to thermostat type */
    ThermostatType get_thermostat_type() const {return thermostatType_;};
        
  protected:

    /** thermostat type flag */
    ThermostatType thermostatType_;

    // thermostat data
    /** lambda power applied to atoms */
    DENS_VEC nodalAtomicLambdaPower_;
    /** filtered lambda power */
    DENS_VEC lambdaPowerFiltered_;

  private:
    
    // DO NOT define this
    Thermostat();
        
  };

  /**
   *  @class  ThermostatShapeFunction
   *  @brief  Base class for implementation of thermostat algorithms using the shape function matrices
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatShapeFunction
  //    base class for all thermostats of general form
  //    of N^T w N lambda = rhs
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatShapeFunction : public RegulatorShapeFunction {
  
  public:
  
    ThermostatShapeFunction(Thermostat * thermostat);
        
    ~ThermostatShapeFunction(){};

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();

  protected:

    // methods
    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights(DIAG_MAT & weights);

    // member data
    /** pointer to thermostat object for data */
    Thermostat * thermostat_;
    /** pointer to lammps atomic velocities */
    double ** v_;
    /** MD mass matrix */
    MATRIX & mdMassMatrix_;
    /** mass of ATC internal atoms on this processor */
    DENS_VEC atomicMass_;

  private:

    // DO NOT define this
    ThermostatShapeFunction();

  };

  /**
   *  @class  ThermostatRescale
   *  @brief  Enforces constraint on atomic kinetic energy based on FE temperature
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatRescale
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatRescale : public ThermostatShapeFunction {
  
  public:
  
    ThermostatRescale(Thermostat * thermostat);
        
    ~ThermostatRescale(){};
        
    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);
        
  protected:

    /** apply solution to atomic quantities */
    void apply_to_atoms(double ** atomicVelocity,
                        double dt);
        
    /** FE temperature field */
    DENS_MAT & nodalTemperature_;

  private:

    // DO NOT define this
    ThermostatRescale();
  
  };

  /**
   *  @class  ThermostatGlc
   *  @brief  Base class for implementation of thermostat algorithms based on Gaussian least constraints (GLC)
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatGlc
  //    base clas for all thermostats of general form of a
  //    Gaussian least constraint (GLC)
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatGlc : public ThermostatShapeFunction {
  
  public:
  
    ThermostatGlc(Thermostat * thermostat);
        
    ~ThermostatGlc(){};

  protected:

    // methods
        
    /** compute force induced by lambda */
    virtual void compute_lambda_force(double * const * atomicQuantity,
                                      DENS_MAT & lambdaForce);

    /** apply forces to atoms */
    virtual void apply_to_atoms(double ** atomicVelocity, 
                                const DENS_MAT & lambdaForce,
                                double dt);

    // member data
    /** pointer to a time filtering object */
    TimeFilter * timeFilter_;

    /** power induced by lambda */
    DENS_VEC & nodalAtomicLambdaPower_;

    /** filtered lambda power */
    DENS_VEC & lambdaPowerFiltered_;

    /** atomic force induced by lambda */
    DENS_MAT & lambdaForce_;

    /** bool to determine if node is fixed or not */
    Array2D<bool> & isFixedNode_;

  private:

    // DO NOT define this
    ThermostatGlc();

  };

  /**
   *  @class  ThermostatPower
   *  @brief  Enforces GLC on atomic forces based on FE power
   *          when using fractional step time integration
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPower
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatPower : public ThermostatGlc {
  
  public:
  
    ThermostatPower(Thermostat * thermostat);
        
    ~ThermostatPower();

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();
        
    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);
        
  protected:

    // methods
    /** destroys allocated memory */
    void destroy();

    /** update the nodal quantities after predictor step */
    virtual void update_nodal_quantities_predictor(double dt);

    /** compute the change in nodal atomic temperature */
    virtual void compute_delta_nodal_atomic_temperature(double dt);

    /** update the nodal quantities after pre-corrector step */
    virtual void update_nodal_quantities_pre_corrector(double dt);

    /** undoes the update the nodal quantities after pre-corrector step */
    virtual void undo_update_nodal_quantities_pre_corrector(double dt);

    /** update the nodal quantities after post-corrector step */
    virtual void update_nodal_quantities_post_corrector(double dt);

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs_nodes,
                                    double dt);

    /** apply forces to atoms */
    virtual void apply_to_atoms(double ** atomicVelocity,
                                const DENS_MAT & lambdaForce,
                                double dt);

    /** solves the non-linear equation for lambda iteratively */
    void iterate_lambda(const MATRIX & rhs,
                        double * const * atomicVelocity,
                        double dt);

    /** set weighting factor for in matrix Nhat^T * weights * Nhat */
    virtual void set_weights(DIAG_MAT & weights);

    // data
    /** reference to AtC FE temperature */
    DENS_MAT & nodalTemperature_;

    /** reference to AtC restricted atomic temperature */
    DENS_MAT & nodalAtomicTemperature_;

    /** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAT & heatSource_;

    /** reference to weights for Nhat */
    DIAG_MAT & nodalWeights_;

    /** change in FE temperature over a timestep */
    DENS_MAT deltaTemperature_;

    /** change in restricted atomic FE temperature over a timestep */
    DENS_MAT deltaNodalAtomicTemperature_;

    /** FE temperature rate of change */
    DENS_MAT nodalTemperatureRoc_;

    /** local version of velocity used as predicted final veloctiy */
    double ** myAtomicVelocity_;

    /** pointer to lammps atomic forces */
    double ** f_;

    /** number of total atoms on this processor */
    int nLocalTotal_;

    /** maximum number of iterations used in iterative solve for lambda */
    int lambdaMaxIterations_;

    /** tolerance used in iterative solve for lambda */
    double lambdaTolerance_;

    /** coefficient to account for effect of time filtering on rhs terms */
    double filterCoefficient_;
    /**  reference to ATC unity shape function on ghost atoms */
    SPAR_MAT & Nstar_;

    /** maps ghost atom and LAMMPS atom ids */
    Array<int> & ghostToAtom_;

  private:

    // DO NOT define this
    ThermostatPower();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatHoover
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatHoover : public ThermostatPower {
  
  public:
  
    ThermostatHoover(Thermostat * thermostat);
        
    ~ThermostatHoover();

    /** reset number of local atoms, as well as atomic data */
    virtual void reset_nlocal();
        
    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** applies thermostat to atoms in the post-corrector phase */
    virtual void apply_post_corrector(double dt);
        
  protected:

    // methods
    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields) {boundaryFlux_[TEMPERATURE] = 0.;};

    /** adds Hoover power terms to nodal variables after the predictor phase */
    void add_to_lambda_power_predictor(double dt);

    /** adds Hoover power terms to nodal variables after the pre-corrector phase */
    void add_to_lambda_power_pre_corrector(double dt);

    /** adds Hoover power terms to nodal variables after the post-corrector phase */
    void add_to_lambda_power_post_corrector(double dt);

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_hoover_rhs(DENS_MAT & rhs_nodes,
                                double dt);

    /** add Hoover contributions to lambda power */
    void add_to_lambda_power(DENS_MAT & myLambdaForce);

    // data
    /** force coming from Hoover contribution */
    DENS_MAT lambdaForceHoover_;

    /** reference to ATC map from global nodes to overlap nodes */
    Array<int> & nodeToOverlapMap_;

  private:

    // DO NOT define this
    ThermostatHoover();

  };
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPowerFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatPowerFiltered : public ThermostatPower {
  
  public:
    
    ThermostatPowerFiltered(Thermostat * thermostat);
    
    ~ThermostatPowerFiltered(){};
    
  protected:

    /** update the nodal quantities after predictor step */
    virtual void update_nodal_quantities_predictor(double dt);

    /** compute the change in nodal atomic temperature */
    virtual void compute_delta_nodal_atomic_temperature(double dt);

    /** update the nodal quantities after pre-corrector step */
    virtual void update_nodal_quantities_pre_corrector(double dt);

    /** update the nodal quantities after post-corrector step */
    virtual void update_nodal_quantities_post_corrector(double dt);
    
    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs_nodes,
                                    double dt);
    
  private:
    
    // DO NOT define this
    ThermostatPowerFiltered();
    
  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPowerVerlet
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatPowerVerlet : public ThermostatGlc {
  
  public:

    ThermostatPowerVerlet(Thermostat * thermostat);
        
    ~ThermostatPowerVerlet(){};
        
    /** applies thermostat to atoms in the predictor phase */
    virtual void apply_pre_predictor(double dt);

    /** applies thermostat to atoms in the pre-corrector phase */
    virtual void apply_pre_corrector(double dt);

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);

  protected:

    /** nodal temperature rate of change */
    DENS_MAT & nodalTemperatureRoc_;

    /** reference to ATC sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAT & heatSource_;

    /** pointer to lammps atomic forces */
    double ** f_;

    /** sets up and solves thermostat equations */
    virtual void compute_thermostat(double dt);

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs_nodes);

    /** computes the nodal FE power applied by the thermostat */
    virtual void compute_nodal_lambda_power(double dt);

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);

  private:

    // DO NOT define this
    ThermostatPowerVerlet();

  };
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatHooverVerlet
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatHooverVerlet : public ThermostatPowerVerlet {
  
  public:

    ThermostatHooverVerlet(Thermostat * thermostat);
        
    ~ThermostatHooverVerlet(){};

  protected:

    /** reference to ATC map from global nodes to overlap nodes */
    Array<int> & nodeToOverlapMap_;

    /** sets up and solves thermostat equations */
    virtual void compute_thermostat(double dt);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields) {boundaryFlux_[TEMPERATURE] = 0.;};

    /** sets up Hoover component of the thermostat */
    void set_hoover_rhs(DENS_MAT & rhs_nodes);

    /** add Hoover contributions to lambda power */
    void add_to_lambda_power(DENS_MAT & myLambdaForce,
                             double dt);

  private:

    // DO NOT implement this
    ThermostatHooverVerlet();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatPowerVerletFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatPowerVerletFiltered : public ThermostatPowerVerlet {
  
  public:

    ThermostatPowerVerletFiltered(Thermostat * thermostat);
        
    ~ThermostatPowerVerletFiltered(){};

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields);

    /** get data for output */
    virtual void output(double dt, OUTPUT_LIST & outputData);

  protected:

    /** sets up appropriate rhs for thermostat equations */
    virtual void set_thermostat_rhs(DENS_MAT & rhs_nodes);

    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & rhs);

    /** nodal temperature 2nd rate of change (i.e. second time derivative) */
    DENS_MAT & nodalTemperature2Roc_;

    /** reference to ATC rate of change sources coming from prescribed data, AtC coupling, and extrinsic coupling */
    DENS_MAT heatSourceRoc_;

    /** references to ATC field rates of changing for inverting the filtered heat sources */
    FIELDS & fieldsRoc_;
   
    /** flux rate of changes for inverting filtered fluxes */
    FIELDS fluxRoc_;

    /** time scale for the time filter */
    double filterScale_;

  private:

    // DO NOT define this
    ThermostatPowerVerletFiltered();

  };

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ThermostatHooverVerletFiltered
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  class ThermostatHooverVerletFiltered : public ThermostatPowerVerletFiltered {
  
  public:

    ThermostatHooverVerletFiltered(Thermostat * thermostat);
        
    ~ThermostatHooverVerletFiltered(){};

  protected:

    /** reference to ATC map from global nodes to overlap nodes */
    Array<int> & nodeToOverlapMap_;

    /** sets up and solves thermostat equations */
    virtual void compute_thermostat(double dt);

    /** compute boundary flux, requires thermostat input since it is part of the coupling scheme */
    virtual void compute_boundary_flux(FIELDS & fields) {boundaryFlux_[TEMPERATURE] = 0.;};

    /** sets up Hoover component of the thermostat */
    void set_hoover_rhs(DENS_MAT & rhs_nodes);

    /** add Hoover contributions to lambda power */
    void add_to_lambda_power(DENS_MAT & myLambdaForce,
                             double dt);

  private:

    // DO NOT implement this
    ThermostatHooverVerletFiltered();

  };

};

#endif

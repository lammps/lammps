#ifndef CHARGE_REGULATOR_H
#define CHARGE_REGULATOR_H

#include "AtomicRegulator.h"
#include "PoissonSolver.h"
#include <map>
#include <utility>
#include <string>
#include <vector>

namespace ATC {

  /**
   *  @class  ChargeRegulator
   *  @brief  Manager class for atom-continuum control of charge and potential
   */

  class ChargeRegulatorMethod;
  class ChargeRegulator : public AtomicRegulator {
  
  public:

    /** charge regulator types */
    enum ChargeRegulatorType {
      NONE=0,
      FEEDBACK,
      IMAGE_CHARGE,
      EFFECTIVE_CHARGE
    };
    enum ChargedSurfaceType {
      INSULATOR=0,
      DIELECTRIC,
      CONDUCTOR
    };

    /** parser data */
    struct ChargeRegulatorParameters {
      ChargeRegulatorParameters():
        method(NONE),
        value(0),
        groupBit(0),
        depth(0),
        permittivity(0),
        surfaceType(INSULATOR) { };
      ChargeRegulatorType method;
      double value; 
      int groupBit;
      std::string groupTag;
      double depth;
      double permittivity; 
      ChargedSurfaceType surfaceType;
      FSET faceset; 
    };

    /** constructor */
    ChargeRegulator(ATC_Coupling *atc);
        
    /** destructor */
    ~ChargeRegulator();
        
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);
    virtual void construct_methods();
    virtual void initialize();
    virtual void apply_pre_force(double dt);
    virtual void apply_post_force(double dt);
    virtual void output(OUTPUT_LIST & outputData) const;
    virtual double compute_vector(int /* n */) const {return 0;} // TODO

    void assign_poisson_solver(PoissonSolver * solver) { poissonSolver_ = solver;}  
    PoissonSolver * poisson_solver(void) { return poissonSolver_;} 

  protected:

    /** poisson solver */
    PoissonSolver * poissonSolver_;


    /** registry charge regulators */
    std::map<std::string,ChargeRegulatorMethod *> regulators_;
    std::map<std::string,ChargeRegulatorParameters> parameters_;

  private:
    ChargeRegulator(); // DO NOT define this
  };

  /**
   *  @class  ChargeRegulatorMethod
   *  @brief  Base class for implementation of ChargeRegulator algorithms 
   */

  class ChargeRegulatorMethod : public RegulatorShapeFunction {
  
  public:
    ChargeRegulatorMethod(ChargeRegulator *chargeRegulator,
      ChargeRegulator::ChargeRegulatorParameters & parameters);
    ~ChargeRegulatorMethod(){};
    virtual void initialize(void); 
    void set_greens_functions();
    virtual void apply_pre_force(double /* dt */){}; 
    virtual void apply_post_force(double /* dt */){}; 
    virtual void set_weights() {};
    const DENS_VEC & total_influence() const { return sum_;}
    virtual void output(OUTPUT_LIST & outputData);

  protected:
    int nlocal();
    /** pointer to charge regulator object for data */
    ChargeRegulator * chargeRegulator_;
    /** interscale manager */
    class InterscaleManager * interscaleManager_;
    /** direct interface to lammps */
    class LammpsInterface * lammpsInterface_;
    /** poisson solver */
    PoissonSolver * poissonSolver_;
    /** (subset) of Green's functions */
    std::vector<SparseVector<double> > greensFunctions_;
    /** cutoff radius */
    double rC_,rCsq_;
    /** conversion constants */
    double qV2e_, qqrd2e_;
    /** member data */
    XT_Function * targetValue_; 
    double targetPhi_; 
    // controlling 
    FSET surface_;
    NSET nodes_;
    int atomGroupBit_;
    bool boundary_; // atoms on boundary 
    DENS_VEC point_;
    DENS_VEC normal_;
    double depth_;
    ChargeRegulator::ChargedSurfaceType surfaceType_;
    double permittivity_;
    bool initialized_;
    DENS_VEC sum_;
    /** workspace */
    mutable DENS_MAT _atomElectricalForce_;
  private:
    ChargeRegulatorMethod(); // DO NOT define this
  };

  /**
   *  @class  ChargeRegulatorMethodFeedback
   *  @brief  ChargeRegulator with feedback and ad hoc positions
   *          can be used for conductor (1 charged sheet) or insulator (2 sheets)
   */

  class ChargeRegulatorMethodFeedback : public ChargeRegulatorMethod {
  
  public:
  
    /** constructor */
    ChargeRegulatorMethodFeedback(ChargeRegulator *chargeRegulator,
      ChargeRegulator::ChargeRegulatorParameters & parameters);
        
    /** destructor */
    ~ChargeRegulatorMethodFeedback(){};

    /** instantiate all needed data */
    virtual void construct_transfers();

    /** initialize */
    virtual void initialize(void); 

    /** set influence nodes and atoms */  
    void set_influence();

    void set_influence_matrix(void); 

    /** post first verlet step */
    virtual void apply_pre_force(double dt); 

    void apply_feedback_charges();

  protected:

    int nControlNodes_; 
    NSET & controlNodes_;
    // measurement/controlled
    int influenceGroupBit_;
    int nInfluenceAtoms_;
    NSET influenceAtoms_, influenceAtomsIds_; 
    int nInfluenceNodes_;
    NSET influenceNodes_;
    
    
    DENS_MAT invG_;
    DENS_MAT invNNT_;
    DENS_MAT NT_;
    DENS_MAT NTinvNNTinvG_;

  private:
    ChargeRegulatorMethodFeedback(); // DO NOT define this
  };

  /**
   *  @class  ChargeRegulatorMethodImageCharge
   *  @brief  ChargeRegulator with image charges
   */

  class ChargeRegulatorMethodImageCharge : public ChargeRegulatorMethod {
  
  public:
  
    /** constructor */
    ChargeRegulatorMethodImageCharge(ChargeRegulator *chargeRegulator,
      ChargeRegulator::ChargeRegulatorParameters & parameters);
        
    /** destructor */
    ~ChargeRegulatorMethodImageCharge(){};

    /** initialize */
    virtual void initialize(void); 

    /** post first verlet step */
    virtual void apply_post_force(double dt); 

  protected:
    double reflect(DENS_VEC & x) const { 
      double dn = (x-point_).dot(normal_);
      x -= 2*dn*normal_; 
      return dn;
    }
    // internal functions
    void apply_local_forces();
    void correct_charge_densities();
    void correct_forces();

    double layerDepth_;
    double permittivityRatio_;
    NSET & imageNodes_; 
    DENS_MAT imageTransferOp_; 

  private:
    ChargeRegulatorMethodImageCharge(); // DO NOT define this
  };

  /**
   *  @class  ChargeRegulatorMethodEffectiveCharge
   *  @brief  ChargeRegulator with effective charges at FE quadrature pts
   */

  class ChargeRegulatorMethodEffectiveCharge : public ChargeRegulatorMethod {

  typedef std::map<int,std::pair<DENS_VEC,double> > NODE_TO_XF_MAP;
  
  public:
  
    /** constructor */
    ChargeRegulatorMethodEffectiveCharge(ChargeRegulator *chargeRegulator,
      ChargeRegulator::ChargeRegulatorParameters & parameters);
        
    /** destructor */
    ~ChargeRegulatorMethodEffectiveCharge(){};

    /** initialize */
    virtual void initialize(void); 

    /** post first verlet step */
    virtual void apply_post_force(double dt); 

  protected:
    // internal functions
    void apply_local_forces();

    // member data
    double chargeDensity_;
    std::map<int,double> nodalChargePotential_;

    NODE_TO_XF_MAP nodeXFMap_;  

    bool useSlab_;


  private:
    ChargeRegulatorMethodEffectiveCharge(); // DO NOT define this
  };
  
};

#endif

#ifndef CONCENTRATION_REGULATOR_H
#define CONCENTRATION_REGULATOR_H

#include <map>
#include <string>

#include "AtomicRegulator.h"
#include "LammpsInterface.h"

namespace ATC {

  /**
   *  @class  ConcentrationRegulator
   *  @brief  Manager class for atom-continuum control of charge and potential
   */

  class ConcentrationRegulatorMethod;
  class ConcentrationRegulator : public AtomicRegulator {

  public:
    enum ConcentrationRegulatorType {NONE=0,TRANSITION};
    /** parser data */
    struct ConcentrationRegulatorParameters {
      ConcentrationRegulatorParameters():
        method(NONE),
        type(0),
        groupbit(0),
        transitionType(0),
        value(0),
        frequency(0),
        transitionInterval(0),
        maxEnergy(0),
        maxExchanges(0),
        maxAttempts(0) { };
      ConcentrationRegulatorType method;
      int type;
      int groupbit;
      int transitionType;
      double value;
      int frequency;
      int transitionInterval;
      double maxEnergy;
      int maxExchanges;
      int maxAttempts;
      std::string transitionTag;
      ESET elemset;
    };

    /** constructor */
    ConcentrationRegulator(ATC_Coupling *atc);
    /** destructor */
    ~ConcentrationRegulator();
    /** parser/modifier */
    virtual bool modify(int narg, char **arg);
    /** construct methods */
    virtual void construct_methods();
    /** pre time integration */
    virtual void initialize();

    //WIP_JAT need a nicer way to consistently handle sets of regulators, not sure how yet
    // application steps
    /** apply the regulator in the pre-predictor phase */
    virtual void apply_pre_predictor(double /* dt */, int /* timeStep */){};
    /** apply the regulator in the mid-predictor phase */
    virtual void apply_mid_predictor(double /* dt */, int /* timeStep */){};
    /** apply the regulator in the post-predictor phase */
    virtual void apply_post_predictor(double /* dt */, int /* timeStep */){};
    /** apply the regulator in the pre-correction phase */
    virtual void apply_pre_corrector(double /* dt */, int /* timeStep */){};
    /** apply the regulator in the post-correction phase */
    virtual void apply_post_corrector(double /* dt */, int /* timeStep */){};
    /** compute the thermal boundary flux, must be consistent with regulator */
    virtual void compute_boundary_flux(FIELDS & /* fields */){};
    /** add contributions (if any) to the finite element right-hand side */
    virtual void add_to_rhs(FIELDS & /* rhs */){};

    /** prior to exchanges */
    virtual void pre_force();
    /** prior to exchanges */
    virtual void pre_exchange();
    /** following a run */
    virtual void finish();
    /** ATC output */
    virtual void output(OUTPUT_LIST & outputData) const;
    /** scalar output */
    virtual double compute_vector(int n) const;
    virtual int size_vector(int s) const;

    bool needs_temperature() { return regulators_.size() > 0; }

  protected:
    /** registry charge regulators */
    std::map<std::string,ConcentrationRegulatorMethod *> regulators_;
    std::map<std::string,ConcentrationRegulatorParameters> parameters_;

  private:
    ConcentrationRegulator(); // DO NOT define this
  };

  /**
   *  @class  ConcentrationRegulatorMethod
   *  @brief  Base class for ConcentrationRegulator algorithms
   */

  class ConcentrationRegulatorMethod : public RegulatorShapeFunction {

  public:
    ConcentrationRegulatorMethod(ConcentrationRegulator *c)
      : RegulatorShapeFunction(c) {};
    ~ConcentrationRegulatorMethod() {};
    void initialize(void) {};
    virtual void pre_force() {};
    virtual void pre_exchange() {};
    virtual void finish() {};
    virtual void set_weights() {};
    virtual double compute_vector(int /* n */) const { return 0;}
    virtual void output(OUTPUT_LIST & /* outputData */){};
  private:
    ConcentrationRegulatorMethod(); // DO NOT define this
  };

  /**
   *  @class  ConcentrationRegulatorMethodTransition
   *  @brief  GCMC + thermodynamic integration
   */
  class ConcentrationRegulatorMethodTransition : public ConcentrationRegulatorMethod {

  public:
    /** constructor */
    ConcentrationRegulatorMethodTransition(
      ConcentrationRegulator *concentrationRegulator,
      ConcentrationRegulator::ConcentrationRegulatorParameters & parameters);
    /** destructor */
    ~ConcentrationRegulatorMethodTransition(){
      if (randomNumberGenerator_) delete randomNumberGenerator_;
    }
    /** initialize */
    void initialize(void);
    /** prior to force evaluation  */
    virtual void pre_force();
    /** prior to exchanges */
    virtual void pre_exchange();
    /** "thermo" output */
    virtual double compute_vector(int n) const;
  protected:
    /** set transition state: epsilon and charge */
    int mask(int /* type */, int /* groupbit */) {return 0;}
    int count(void) const;
    int excess(void) const;
    double energy(int id) const;
    double energy(double * x) const;
    double insertion_location(DENS_VEC & x) const;
    double deletion_id(ID_PAIR & id) const ;
    double deletion_id_consistent(ID_PAIR & id) const ;
    double deletion_id_free(ID_PAIR & id) const ;
    double insertion_factor(int c) const // a ramp
    {
      if (c < transitionInterval_) return ((double) c)/transitionInterval_;
      else                         return 1.0;
    }
    void transition();
    bool accept(double energy, double T = 0) const;
    bool delete_atoms(int n);
    bool insert_atoms(int n);
    int pick_element() const;
    void pick_coordinates(const int elem, DENS_VEC & xi, DENS_VEC& x ) const;
    void pick_velocity(DENS_VEC & v, const double T ) const;
    double uniform() const;
    double normal() const;

    /** pointer to conc regulator object for data */
    ConcentrationRegulator * concentrationRegulator_;

    /** interscale manager */
    class InterscaleManager * interscaleManager_;

    /** interscale manager */
    class LammpsInterface * lammpsInterface_;

    /** member data */
    class AtomInElementSet * list_;
    double targetConcentration_;
    double targetCount_;
    ESET elemset_;
    POTENTIAL p_;
    RNG_POINTER randomNumberGenerator_;
    DENS_VEC volumes_;
    double V_;
    double q0_;
    int controlType_;
    int controlIndex_;
    int transitionType_;
    int transitionInterval_;;
    int transitionCounter_;
    int nInTransition_;
    double transitionFactor_;
    DENS_VEC epsilon0_;
    ID_LIST deletionIds_, insertionIds_;
    int controlMask_;
    int frequency_;
    double maxEnergy_;
    int maxExchanges_;
    int maxAttempts_;
    int nexchanges_;
    bool initialized_;
    double sigma_;

    void sync_random_number_generators() const;
    mutable int _rngUniformCounter_; // for parallel consistency
    mutable int _rngNormalCounter_; // for parallel consistency

  private:
    ConcentrationRegulatorMethodTransition(); // DO NOT define this
  };
};
#endif

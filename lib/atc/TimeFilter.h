// The type of filter used should be determined by the 
// integrator since filtering much match the time integration scheme


#ifndef TIME_FILTER_H
#define TIME_FILTER_H

#include <set>

#include "ATC_TypeDefs.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"

namespace ATC {

  // forward declarations
  class  ATC_Method;
  class  TimeFilter;

    /**
   *  @class  TimeFilterManager
   *  @brief  Handles parsing and parameter storage for time filters
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterManager
  //--------------------------------------------------------
  //--------------------------------------------------------

  class TimeFilterManager {
  
  public:

    /** enumeration for the functional form underlying the filter */
    enum TimeFilterType {
      NO_FILTER=0, // default
      EXPONENTIAL_FILTER,
      STEP_FILTER
    };

    /** enumeration for the functional form underlying the filter */
    enum FilterIntegrationType {
      INSTANTANEOUS=0, //default
      CRANK_NICHOLSON,
      IMPLICIT_EXPLICIT,
      EXPLICIT_IMPLICIT,
      EXPLICIT,
      IMPLICIT,
      IMPLICIT_UPDATE
    };
  
    // constructor
    TimeFilterManager(ATC_Method * atc);
    
    // destructor
    ~TimeFilterManager();
        
    /** parser/modifier */
    bool modify(int narg, char **arg);

    /** pre time integration */
    void initialize();

    /** get filter base function */
    TimeFilterType filter_type() const {return filterType_;};
        
    /** return filtering time scale */
    double filter_scale() const {return filterScale_;};
        
    /** check if dynamics should be filtering */
    bool filter_dynamics() const {return useFilter_;};

    /** check if variables should be filtered */
    bool filter_variables() const {return (useFilter_ || equilibrateFilter_);};

    /** flag for if reset is needed */
    bool need_reset() const {return needReset_;};

    /** flag if ending equilibration */
    bool end_equilibrate() const {return endEquilibrate_;};

    /** get pointer to ATC transfer methods */
    ATC_Method * atc() {return atc_;};

    /** construct the appropriate time filter */
    TimeFilter * construct(const FilterIntegrationType type = CRANK_NICHOLSON);
  
  protected:

    TimeFilterManager(){};
  
    /** pointer to access ATC methods */
    ATC_Method * atc_;

    /** description of underlying function form of filter */
    TimeFilterType filterType_;
  
    /** filtering time scale */
    double filterScale_;
        
    /** flag to see if filtering is active */
    bool useFilter_;

    /** flag to see if we are equilibrating the filtered variables */
    bool equilibrateFilter_;
        
    /** flag to reset data */
    bool needReset_;

    /** flag to denote switch from equilibration to integration */
    bool endEquilibrate_;

    /** set to store all time filters for later deletion */
    std::set<TimeFilter * > timeFilterSet_;
  
  };

  /**
   *  @class  TimeFilter
   *  @brief  Base class for various temporal filters of atomistic quantities
   *          default behavior is no filter
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilter
  //--------------------------------------------------------
  //--------------------------------------------------------


  class TimeFilter {
  
  public:
  
    // constructor
    TimeFilter(TimeFilterManager & timeFilterManager);
    
    // destructor
    virtual ~TimeFilter(){};

    /** pre time integration */
    virtual void initialize(){};

    /** pre time integration with a target for an initial condition */
    virtual void initialize(const MATRIX & target){initialize();};
        
    /** Step 1:
        apply first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 const MATRIX & unFilteredQuantity,
                                 double dt)
      { TimeFilter::unFilteredQuantityOld_ = unFilteredQuantity;}
        
    /** Step 2:
         apply second step in a time filter update in pre integration phase */
    virtual void apply_pre_step2(MATRIX & filteredQuantity,
                                 const MATRIX & unFilteredQuantity,
                                 double dt) {};
        
    /** Step 3:
        apply first step in a time filter update in post integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                  const MATRIX & unFilteredQuantity,
                                  double dt) 
      { filteredQuantity = unFilteredQuantity;};
        
    /** Step 4:
        apply second step in a time filter update in post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  const MATRIX & unFilteredQuantity,
                                  double dt)
      { filteredQuantity = unFilteredQuantity;}
                                                  
    /** coefficient multipling unfiltered terms in apply_pre_step1 method */
    virtual double unfiltered_coefficient_pre_s1(double dt){return 0.;};

    /** coefficient multipling old filtered terms in apply_pre_step1 method */
    virtual double filtered_coefficient_pre_s1(double dt){return 0.;};
        
    /** coefficient multipling unfiltered terms in apply_post_step1 method */
    virtual double unfiltered_coefficient_post_s1(double dt){return 0.;};

    /** coefficient multipling old filtered terms in apply_post_step1 method */
    virtual double filtered_coefficient_post_s1(double dt){return 0.;};
        
    /** coefficient multipling unfiltered terms in apply_pre_step2 method */
    virtual double unfiltered_coefficient_pre_s2(double dt){return 0.;};

    /** coefficient multipling old filtered terms in apply_pre_step2 method */
    virtual double filtered_coefficient_pre_s2(double dt){return 0.;};
        
    /** coefficient multipling unfiltered terms in apply_post_step2 method */
    virtual double unfiltered_coefficient_post_s2(double dt){return 0.;};

    /** coefficient multipling old filtered terms in apply_post_step2 method */
    virtual double filtered_coefficient_post_s2(double dt){return 0.;};
  
    /** rate of filtered quantity to be called in post integration phase */
    virtual void rate(MATRIX & rate,
                      const MATRIX & filteredQuantity,
                      const MATRIX & unFilteredQuantity,
                      double dt = 0.0)
      { rate = 1/dt*(unFilteredQuantity - TimeFilter::unFilteredQuantityOld_);}; 

  protected:
  
    TimeFilter(){};

    /** pointer to access ATC methods */
    ATC_Method * atc_;

    /** filtering time scale */
    double filterScale_;

    /** filter type */
    TimeFilterManager::TimeFilterType filterType_;

    /** member data to track old unfiltered values */
    DENS_MAT unFilteredQuantityOld_;
  
  };

  /**
   *  @class  TimeFilterExponential
   *  @brief  Base class for filters using an exponential kernel,
   *          derived classes implement specific integration schemes
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterExponential
  //--------------------------------------------------------
  //--------------------------------------------------------

  class TimeFilterExponential : public TimeFilter {
  
  public:
  
    // constructor
    TimeFilterExponential(TimeFilterManager & timeFilterManager);
    
    // destructor
    virtual ~TimeFilterExponential(){};
    /** apply first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 const MATRIX & unFilteredQuantity,
                                 double dt) {};
        
    /** apply second step in a time filter update in pre integration phase */
    virtual void apply_pre_step2(MATRIX & filteredQuantity,
                                 const MATRIX & unFilteredQuantity,
                                 double dt) {};
        
    /** apply first step in a time filter update in post integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                  const MATRIX & unFilteredQuantity,
                                  double dt) {};
        
    /** apply second step in a time filter update in post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  const MATRIX & unFilteredQuantity,
                                  double dt) {};

    /** time rate of filtered quantity */
    virtual void rate(MATRIX & rate,
        const MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double dt = 0)
    { double tau = TimeFilter::filterScale_;
      rate = 1/tau*(unfilteredQuantity - filteredQuantity); };

  protected:

    TimeFilterExponential(){};
  
    //--------------------------------------------------------
    //--------------------------------------------------------
    //  filter integration functions not associated
    //  with any particular class
    //--------------------------------------------------------
    //--------------------------------------------------------
  
    void update_filter(MATRIX & filteredQuantity,
                       const MATRIX & unfilteredQuantity,
                       MATRIX & unfilteredQuantityOld,
                       double tau,
                       double dt)
    {
      filteredQuantity *= 1./(1./dt+1./(2*tau))*((1./dt-1./(2*tau)));
      filteredQuantity += 1./(1./dt+1./(2*tau))*( 1./(2*tau)*(unfilteredQuantity+unfilteredQuantityOld));
      unfilteredQuantityOld = unfilteredQuantity;
    }

    void add_to_filter(MATRIX & filteredQuantity,
                       const MATRIX & unfilteredQuantity,
                       MATRIX & unfilteredQuantityOld,
                       double tau,
                       double dt)
    {
      filteredQuantity += 1./(1./dt+1./(2.*tau))*( 1./(2.*tau))*unfilteredQuantity;
      unfilteredQuantityOld += unfilteredQuantity;
    };

    double unfiltered_coef(double tau, double dt)
    { return 1./(1./dt+1./(2.*tau))*( 1./(2.*tau)); };

    double filtered_coef(double tau, double dt)
    { return 1./(1./dt+1./(2.*tau))*( 1./dt-1./(2*tau) ); };
  
    void update_filter_implicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    { 
      filteredQuantity /= 1.0 + dt/tau;
      filteredQuantity +=  (dt)/(tau+dt)*unfilteredQuantity;
    }

    void add_to_filter_implicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    { filteredQuantity += (1./(1.+dt/tau))*(dt/tau)*unfilteredQuantity; };
  
    double unfiltered_coef_implicit(double tau, double dt)
    { return (1./(1.+dt/tau))*(dt/tau); };

    double filtered_coef_implicit(double tau, double dt)
    { return (1./(1.+dt/tau)); };
  
    void update_filter_explicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    {
      filteredQuantity *= (1.-(dt/tau));
      filteredQuantity += (dt/tau)*unfilteredQuantity;
    }
  
    void add_to_filter_explicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    { filteredQuantity += (dt/tau)*unfilteredQuantity; };
  
    double unfiltered_coef_explicit(double tau, double dt)
    { return (dt/tau); };

    double filtered_coef_explicit(double tau, double dt)
    { return (1.-(dt/tau)); };

  };

  /**
   *  @class  TimeFilterCrankNicolson
   *  @brief  Time Filter using Crank-Nicolson advancement of filtered quantity ODE's
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterCrankNicolson
  //--------------------------------------------------------
  //--------------------------------------------------------
  class TimeFilterCrankNicolson : public TimeFilterExponential {
  
  public:
  
    // constructor
    TimeFilterCrankNicolson(TimeFilterManager & timeFilterManager);
        
    // destructor
    virtual ~TimeFilterCrankNicolson(){};

    /** pre time integration */
    virtual void initialize(){throw ATC_Error("TimeFilterCrankNicolson::initialize() an initial condition is required for this time filter");};
        
    /** pre time integration with an initial condition */
    virtual void initialize(const MATRIX & target);
        
    /** applies first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { update_filter(filteredQuantity,unFilteredQuantity,unFilteredQuantityOld_,TimeFilter::filterScale_,dt); };

    /** applies first step in a time filter update after the pre integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { add_to_filter(filteredQuantity,unFilteredQuantity,unFilteredQuantityOld_,TimeFilter::filterScale_,dt); };
        
    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { update_filter(filteredQuantity,unFilteredQuantity,unFilteredQuantityOld_,TimeFilter::filterScale_,dt); };
        
    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double unfiltered_coefficient_pre_s1(double dt){return unfiltered_coef(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_pre_step1 method */
    virtual double filtered_coefficient_pre_s1(double dt){return filtered_coef(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double unfiltered_coefficient_post_s2(double dt){return unfiltered_coef(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_post_step2 method */
    virtual double filtered_coefficient_post_s2(double dt){return filtered_coef(TimeFilter::filterScale_,dt);};
        
  protected:
  
    TimeFilterCrankNicolson();
  
  };
  
  /**
   *  @class  TimeFilterExplicit
   *  @brief  Time Filter using explicit advancement of filtered quantity ODE's
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterExplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  class TimeFilterExplicit : public TimeFilterExponential {
  
  public:
  
    // constructor
    TimeFilterExplicit(TimeFilterManager & timeFilterManager);
        
    // destructor
    virtual ~TimeFilterExplicit(){};
        
    /** applies first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    {update_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };

    /** applies first step in a time filter update after the pre integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { add_to_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { update_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double unfiltered_coefficient_pre_s1(double dt){return unfiltered_coef_explicit(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_pre_step1 method */
    virtual double filtered_coefficient_pre_s1(double dt){return filtered_coef_explicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double unfiltered_coefficient_post_s2(double dt){return unfiltered_coef_explicit(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_post_step2 method */
    virtual double filtered_coefficient_post_s2(double dt){return filtered_coef_explicit(TimeFilter::filterScale_,dt);};
        
  protected:
  
    TimeFilterExplicit();
        
  };
  
  /**
   *  @class  TimeFilterImplicit
   *  @brief  Time Filter using implicit advancement of filtered quantity ODE's
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterImplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  class TimeFilterImplicit : public TimeFilterExponential {
  
  public:
  
    // constructor
    TimeFilterImplicit(TimeFilterManager & timeFilterManager);
        
    // destructor
    virtual ~TimeFilterImplicit(){};
        
    /** applies first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { update_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { update_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };

    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double unfiltered_coefficient_pre_s1(double dt){return unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_pre_step1 method */
    virtual double filtered_coefficient_pre_s1(double dt){return filtered_coef_implicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double unfiltered_coefficient_post_s2(double dt){return unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_post_step2 method */
    virtual double filtered_coefficient_post_s2(double dt){return filtered_coef_implicit(TimeFilter::filterScale_,dt);};
        
  protected:
  
    TimeFilterImplicit();
  
  };
  
  /**
   *  @class  TimeFilterImplicitExplicit
   *  @brief  Time Filter using two-step implicit/explicit advancement of filtered quantity ODE's
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterImplicitExplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  class TimeFilterImplicitExplicit : public TimeFilterExponential {
  
  public:
  
    // constructor
    TimeFilterImplicitExplicit(TimeFilterManager & timeFilterManager);
        
    // destructor
    virtual ~TimeFilterImplicitExplicit(){};
        
    /** applies first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { update_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,0.5*dt); };
        
    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { update_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,0.5*dt); };
        
    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double unfiltered_coefficient_pre_s1(double dt){return unfiltered_coef_implicit(TimeFilter::filterScale_,0.5*dt);};

    /** return coefficient multipling old filtered terms in the apply_pre_step1 method */
    virtual double filtered_coefficient_pre_s1(double dt){return filtered_coef_implicit(TimeFilter::filterScale_,0.5*dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double unfiltered_coefficient_post_s2(double dt){return unfiltered_coef_explicit(TimeFilter::filterScale_,0.5*dt);};

    /** return coefficient multipling old filtered terms in the apply_post_step2 method */
    virtual double filtered_coefficient_post_s2(double dt){return filtered_coef_explicit(TimeFilter::filterScale_,0.5*dt);};
        
  protected:
  
    TimeFilterImplicitExplicit();
  
  };

    /**
   *  @class  TimeFilterExplicitImplicit
   *  @brief  Time Filter using two-step explicit/implicit advancement of filtered quantity ODE's
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterExplicitImplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  class TimeFilterExplicitImplicit : public TimeFilterExponential {
  
  public:
  
    // constructor
    TimeFilterExplicitImplicit(TimeFilterManager & timeFilterManager);
        
    // destructor
    virtual ~TimeFilterExplicitImplicit(){};
        
    /** applies first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { update_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,0.5*dt); };

    /** applies first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step2(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { add_to_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,0.5*dt); };

    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { update_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,0.5*dt); };

    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { add_to_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,0.5*dt); };
        
    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double unfiltered_coefficient_pre_s1(double dt){return unfiltered_coef_explicit(TimeFilter::filterScale_,0.5*dt);};

    /** return coefficient multipling old filtered terms in the apply_pre_step1 method */
    virtual double filtered_coefficient_pre_s1(double dt){return filtered_coef_explicit(TimeFilter::filterScale_,0.5*dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double unfiltered_coefficient_post_s1(double dt){return unfiltered_coef_implicit(TimeFilter::filterScale_,0.5*dt);};

    /** return coefficient multipling old filtered terms in the apply_post_step2 method */
    virtual double filtered_coefficient_post_s1(double dt){return filtered_coef_implicit(TimeFilter::filterScale_,0.5*dt);};

  protected:
  
    TimeFilterExplicitImplicit();
  
  };

  /**
   *  @class  TimeFilterImplicitUpdate
   *  @brief  Time Filter using implicit advancement of filtered quantity ODE's but adds on contribution at the end of the second step
   */
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterImplicitUpdate
  //--------------------------------------------------------
  //--------------------------------------------------------
  class TimeFilterImplicitUpdate : public TimeFilterExponential {
  
  public:
  
    // constructor
    TimeFilterImplicitUpdate(TimeFilterManager & timeFilterManager);
        
    // destructor
    virtual ~TimeFilterImplicitUpdate(){};
        
    /** applies first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { update_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { add_to_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double unfiltered_coefficient_pre_s1(double dt){return unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_pre_step1 method */
    virtual double filtered_coefficient_pre_s1(double dt){return filtered_coef_implicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double unfiltered_coefficient_post_s1(double dt){return unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};

    /** return coefficient multipling old filtered terms in the apply_post_step2 method */
    virtual double filtered_coefficient_post_s1(double dt){return filtered_coef_implicit(TimeFilter::filterScale_,dt);};
        
  protected:
  
    TimeFilterImplicitUpdate();
  
  };
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterStep
  //--------------------------------------------------------
  //--------------------------------------------------------

  class TimeFilterStep : public TimeFilter {
  
  public:
  
    // constructor
    TimeFilterStep(TimeFilterManager & timeFilterManager);
    
    // destructor
    virtual ~TimeFilterStep(){};

    /** pre time integration */
    virtual void initialize(const MATRIX & target);

    /** apply first step in a time filter update in the pre integration phase */
    virtual void apply_pre_step1(MATRIX & filteredQuantity,
                                 const MATRIX & unFilteredQuantity,
                                 double dt) {};
        
    /** apply second step in a time filter update in pre integration phase */
    virtual void apply_pre_step2(MATRIX & filteredQuantity,
                                 const MATRIX & unFilteredQuantity,
                                 double dt) {};
        
    /** apply first step in a time filter update in post integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                  const MATRIX & unFilteredQuantity,
                                  double dt) {};
        
    /** apply second step in a time filter update in post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  const MATRIX & unFilteredQuantity,
                                  double dt) 
      { update_filter(filteredQuantity, unFilteredQuantity,
           TimeFilter::unFilteredQuantityOld_, TimeFilter::filterScale_, dt); 
      }

    /** time rate of filtered quantity */
    virtual void rate(MATRIX & rate,
        const MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double dt = 0)
    { rate = 1/elapsedTime_*(unfilteredQuantity - filteredQuantity); }

  protected:

    TimeFilterStep(){};
   
    double elapsedTime_;
  
    void update_filter(MATRIX & filteredQuantity,
                       const MATRIX & unfilteredQuantity,
                       MATRIX & unfilteredQuantitySum,
           double tau,
           double dt)
    {

// this average and the next
      if (elapsedTime_ == 0.0) { // a reset
        elapsedTime_ = dt;
        unfilteredQuantitySum = unfilteredQuantity*dt;
        filteredQuantity = unfilteredQuantity;
      }
      else { // a running average
        elapsedTime_ += dt;
        unfilteredQuantitySum += unfilteredQuantity*dt;
        filteredQuantity = unfilteredQuantitySum; 
        filteredQuantity /= elapsedTime_;
      }
      if (elapsedTime_ >= tau && tau > 0) { 
        elapsedTime_ = 0.0;
      }
    };
  };



};

#endif

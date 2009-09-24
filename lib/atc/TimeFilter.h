// The type of filter used should be determined by the 
// integrator since filtering much match the time integration scheme

#ifndef TIME_FILTER_H
#define TIME_FILTER_H

// ATC_Transfer headers
#include "MatrixLibrary.h"
#include "ATC_Error.h"

using namespace std;
namespace ATC {

  // forward declarations
  class  ATC_Transfer;
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
      CRANK_NICHOLSON=0, // default
      IMPLICIT_EXPLICIT,
      EXPLICIT_IMPLICIT,
      EXPLICIT,
      IMPLICIT,
      IMPLICIT_UPDATE
    };
  
    // constructor
    TimeFilterManager(ATC_Transfer * atcTransfer);
    
    // destructor
    ~TimeFilterManager(){};
        
    /** parser/modifier */
    bool modify(int narg, char **arg);

    /** pre time integration */
    void initialize();

    /** get filter base function */
    TimeFilterType get_filter_type() const {return filterType_;};
        
    /** return filtering time scale */
    double get_filter_scale() const {return filterScale_;};
        
    /** check if dynamics should be filtering */
    bool filter_dynamics() const {return useFilter_;};

    /** check if variables should be filtered */
    bool filter_variables() const {return (useFilter_ || equilibrateFilter_);};

    /** flag for if reset is needed */
    bool need_reset() const {return needReset_;};

    /** flag if ending equilibration */
    bool end_equilibrate() const {return endEquilibrate_;};

    /** get pointer to ATC transfer methods */
    ATC_Transfer * get_atc_transfer() {return atcTransfer_;};

    /** construct the appropriate time filter */
    TimeFilter * construct(const FilterIntegrationType type = CRANK_NICHOLSON);
  
  protected:

    TimeFilterManager(){};
  
    /** pointer to access ATC methods */
    ATC_Transfer * atcTransfer_;

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
    virtual void initialize(const MATRIX & target){};
        
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
                                  double dt) {};
        
    /** Step 4:
        apply second step in a time filter update in post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  const MATRIX & unFilteredQuantity,
                                  double dt)
      { filteredQuantity = unFilteredQuantity;}
                                                  
    /** coefficient multipling unfiltered terms in apply_pre_step1 method */
    virtual double get_unfiltered_coefficient_pre_s1(double dt){return 0.;};
        
    /** coefficient multipling unfiltered terms in apply_post_step1 method */
    virtual double get_unfiltered_coefficient_post_s1(double dt){return 0.;};
        
    /** coefficient multipling unfiltered terms in apply_pre_step2 method */
    virtual double get_unfiltered_coefficient_pre_s2(double dt){return 0.;};
        
    /** coefficient multipling unfiltered terms in apply_post_step2 method */
    virtual double get_unfiltered_coefficient_post_s2(double dt){return 0.;};
  
    /** rate of filtered quantity to be called in post integration phase */
    virtual void rate(MATRIX & rate,
                      const MATRIX & filteredQuantity,
                      const MATRIX & unFilteredQuantity,
                      double dt = 0.0)
      { rate = 1/dt*(unFilteredQuantity - TimeFilter::unFilteredQuantityOld_);}; 

  protected:
  
    TimeFilter(){};

    /** pointer to access ATC methods */
    ATC_Transfer * atcTransfer_;

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
      filteredQuantity = 1./(1./dt+1./(2*tau))*( 1./(2*tau)*
                                                 (unfilteredQuantity+unfilteredQuantityOld) +
                                                 (1./dt-1./(2*tau))*filteredQuantity);
      unfilteredQuantityOld = unfilteredQuantity;
    };

    void add_to_filter(MATRIX & filteredQuantity,
                       const MATRIX & unfilteredQuantity,
                       MATRIX & unfilteredQuantityOld,
                       double tau,
                       double dt)
    {
      filteredQuantity += 1./(1./dt+1./(2.*tau))*( 1./(2.*tau))*unfilteredQuantity;
      unfilteredQuantityOld += unfilteredQuantity;
    };

    double get_unfiltered_coef(double tau, double dt)
    { return 1./(1./dt+1./(2.*tau))*( 1./(2.*tau)); };
  
    void update_filter_implicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    // TODO: replace the rest of these like below:
    { filteredQuantity = (1./(1.+dt/tau))*((dt/tau)*unfilteredQuantity + filteredQuantity); };
    //    { 
    //      filteredQuantity /= 1.0 + dt/tau;
    //      filteredQuantity +=  (dt)/(tau+dt)*unfilteredQuantity;
    //    }

    void add_to_filter_implicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    { filteredQuantity += (1./(1.+dt/tau))*(dt/tau)*unfilteredQuantity; };
  
    double get_unfiltered_coef_implicit(double tau, double dt)
    { return (1./(1.+dt/tau))*(dt/tau); };
  
    void update_filter_explicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    { filteredQuantity = (dt/tau)*unfilteredQuantity + (1.-dt/tau)*filteredQuantity; };
  
    void add_to_filter_explicit(MATRIX & filteredQuantity,
                                const MATRIX & unfilteredQuantity,
                                double tau,
                                double dt)
    { filteredQuantity += (dt/tau)*unfilteredQuantity; };
  
    double get_unfiltered_coef_explicit(double dt, double tau)
    { return (dt/tau); };

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
    virtual double get_unfiltered_coefficient_pre_s1(double dt){return get_unfiltered_coef(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double get_unfiltered_coefficient_post_s2(double dt){return get_unfiltered_coef(TimeFilter::filterScale_,dt);};
        
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
    virtual double get_unfiltered_coefficient_pre_s1(double dt){return get_unfiltered_coef_explicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double get_unfiltered_coefficient_post_s2(double dt){return get_unfiltered_coef_explicit(TimeFilter::filterScale_,dt);};
        
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
    virtual double get_unfiltered_coefficient_pre_s1(double dt){return get_unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double get_unfiltered_coefficient_post_s2(double dt){return get_unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};
        
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
    { update_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { update_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double get_unfiltered_coefficient_pre_s1(double dt){return get_unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double get_unfiltered_coefficient_post_s2(double dt){return get_unfiltered_coef_explicit(TimeFilter::filterScale_,dt);};
        
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
    { update_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };

    /** applies second step in a time filter update in the pre integration phase */
    virtual void apply_pre_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { add_to_filter_explicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };

    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step1(MATRIX & filteredQuantity,
                                 MATRIX const & unFilteredQuantity,
                                 double dt)
    { update_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };

    /** applies second step in a time filter update in the post integration phase */
    virtual void apply_post_step2(MATRIX & filteredQuantity,
                                  MATRIX const & unFilteredQuantity,
                                  double dt)
    { add_to_filter_implicit(filteredQuantity,unFilteredQuantity,TimeFilter::filterScale_,dt); };
        
    /** return coefficient multipling unfiltered terms in the apply_pre_step1 method */
    virtual double get_unfiltered_coefficient_pre_s1(double dt){return get_unfiltered_coef_explicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double get_unfiltered_coefficient_post_s1(double dt){return get_unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};

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
    virtual double get_unfiltered_coefficient_pre_s1(double dt){return get_unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};
        
    /** return coefficient multipling unfiltered terms in the apply_post_step2 method */
    virtual double get_unfiltered_coefficient_post_s1(double dt){return get_unfiltered_coef_implicit(TimeFilter::filterScale_,dt);};
        
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
      elapsedTime_ += dt;
      if (elapsedTime_ > tau) {
        elapsedTime_ = dt;
        unfilteredQuantitySum = unfilteredQuantity*dt;
        filteredQuantity = unfilteredQuantity;
      }
      else {
        unfilteredQuantitySum += unfilteredQuantity*dt;
        filteredQuantity = unfilteredQuantitySum; 
        filteredQuantity /= elapsedTime_;
      }
    };
  };

};

#endif

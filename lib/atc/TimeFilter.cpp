#include "ATC_Transfer.h"
#include "TimeFilter.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterManager
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterManager::TimeFilterManager(ATC_Method * atc) :
    atc_(atc),
    filterType_(NO_FILTER),
    filterScale_(0.),
    useFilter_(false),
    equilibrateFilter_(false),
    needReset_(true),
    endEquilibrate_(false)
  {
    // do nothing
  }
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterManager::~TimeFilterManager()
  {
    set<TimeFilter * >::iterator it;
    for (it = timeFilterSet_.begin(); it != timeFilterSet_.end(); it++)
      if (*it) delete *it;
  }
  
  //--------------------------------------------------------
  //  modify
  //    parses input commands
  //--------------------------------------------------------
  bool TimeFilterManager::modify(int narg, char ** arg)
  {
    bool foundMatch = false;
  
    // filter scale size
    /*! \page man_filter_scale fix_modify AtC filter scale
      \section syntax
      fix_modify AtC filter scale <scale> \n
      - scale (real) = characteristic time scale of the filter  \n

      \section examples
      <TT> fix_modify AtC filter scale 10.0 </TT> \n

      \section description
      Filters the MD dynamics to construct a more appropriate continuous field.  Equilibrating first filters the time derivatives without changing the dynamics to provide a better initial condition to the filtered dynamics
      \section restrictions
      only for be used with specific transfers:
      thermal, two_temperature

      \section related
      \ref man_time_filter
      \ref man_filter_type

      \section default
      0.
    */
    if (strcmp(arg[0],"scale")==0) {
      filterScale_ = atof(arg[1]);
      if (filterScale_<=0.)
        throw ATC_Error("Bad filtering time scale");
      foundMatch = true;
    }

    // time filtering activation switch
    /*! \page man_time_filter fix_modify AtC filter
      \section syntax
      fix_modify AtC filter <on | off | equilibrate> \n
      - on | off (keyword) = turns filter on or off\n
      - equilibrate = runs dynamics without filtering but initializes filtered quantities
      \section examples
      <TT> fix_modify atc transfer filter on </TT> \n

      \section description
      Filters the MD dynamics to construct a more appropriate continuous field.  Equilibrating first filters the time derivatives without changing the dynamics to provide a better initial condition to the filtered dynamics

      \section restrictions
      only for be used with specific transfers:
      thermal, two_temperature

      \section related
      \ref man_filter_scale \n
      \ref man_equilibrium_start

      \section default
      off
    */
    else if (strcmp(arg[0],"on")==0) { 
      if (filterScale_<=0. && filterType_ != STEP_FILTER)
        throw ATC_Error("Filtering time scale not initialized");
      useFilter_ = true;
      if (!equilibrateFilter_) {
        needReset_ = true;
        endEquilibrate_ = false;
      }
      else {
        endEquilibrate_ = true;
      }
      equilibrateFilter_ = false;
      foundMatch = true;
    }
    else if (strcmp(arg[0],"off")==0) { 
      useFilter_ = false;
      equilibrateFilter_ = false;
      endEquilibrate_ = false;
      needReset_ = true;
      foundMatch = true;
    }
    else if (strcmp(arg[0],"equilibrate")==0) {
      if (filterScale_<=0. && filterType_ != STEP_FILTER)
        throw ATC_Error("Filtering time scale not initialized");
      equilibrateFilter_ = true;
      endEquilibrate_ = false;
      useFilter_ = false;
      needReset_ = true;
      foundMatch = true;
    }

    // filter type 
    /*! \page man_filter_type fix_modify AtC filter type 
      \section syntax
      fix_modify AtC filter type <exponential | no_filter> \n

      \section examples
      <TT> fix_modify AtC filter type exponential </TT> \n

      \section description
      Specifies the type of time filter used.
      \section restrictions
      only for be used with specific transfers:
      thermal, two_temperature

      \section related
      \ref man_time_filter
      \ref man_filter_scale

      \section default
      No default. 
    */
    else if (strcmp(arg[0],"type")==0) { 
      if (strcmp(arg[1],"exponential")==0) { 
        filterType_ = EXPONENTIAL_FILTER;
      }
      else if (strcmp(arg[1],"step")==0) { 
        filterType_ = STEP_FILTER;
      }
      else if (strcmp(arg[1],"no_filter")==0) { 
        filterType_ = NO_FILTER;
      }
      else throw ATC_Error("Not a supported time filter type");
      foundMatch = true;
    }


    return foundMatch;
  }

  //--------------------------------------------------------
  //  initialize
  //    filter set up before a run
  //--------------------------------------------------------
  void TimeFilterManager::initialize()
  {
    needReset_ = false;
    endEquilibrate_ = false;
  }

  //--------------------------------------------------------
  //  construct
  //    instantiates the filter
  //--------------------------------------------------------
  
  TimeFilter * TimeFilterManager::construct(const FilterIntegrationType intType)
  {
    TimeFilter * newTimeFilter;
    if (useFilter_ || equilibrateFilter_) {
      if (filterType_ == EXPONENTIAL_FILTER) {
        if (intType == IMPLICIT_EXPLICIT) {
          newTimeFilter = new TimeFilterImplicitExplicit(*this);     
        }
        else if (intType == EXPLICIT_IMPLICIT) {
          newTimeFilter = new TimeFilterExplicitImplicit(*this);     
        }
        else if (intType == EXPLICIT) {
          newTimeFilter = new TimeFilterExplicit(*this);     
        }
        else if (intType == IMPLICIT) {
          newTimeFilter = new TimeFilterImplicit(*this);     
        }
        else if (intType == IMPLICIT_UPDATE) {
          newTimeFilter = new TimeFilterImplicitUpdate(*this);     
        }
        else if (intType == CRANK_NICHOLSON) {
          newTimeFilter = new TimeFilterCrankNicolson(*this);     
        }
        else { // default to return base class
          newTimeFilter = new TimeFilter(*this); 
        }
      }
      else if (filterType_ == STEP_FILTER) {
        newTimeFilter = new TimeFilterStep(*this); 
      }
    }
    else { // default to return base class
      newTimeFilter = new TimeFilter(*this); 
    }
    timeFilterSet_.insert(newTimeFilter);
    return newTimeFilter;
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilter
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilter::TimeFilter(TimeFilterManager & timeFilterManager) :
    atc_(timeFilterManager.atc()),
    filterScale_(timeFilterManager.filter_scale())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterExponential
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterExponential::TimeFilterExponential(TimeFilterManager & timeFilterManager) :
    TimeFilter(timeFilterManager)
  {
    TimeFilter::filterType_ = TimeFilterManager::EXPONENTIAL_FILTER;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterCrankNicoslon
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterCrankNicolson::TimeFilterCrankNicolson(TimeFilterManager & timeFilterManager) :
    TimeFilterExponential(timeFilterManager)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize
  //    resets memory as needed
  //--------------------------------------------------------
  void TimeFilterCrankNicolson::initialize(const MATRIX & target)
  {
    TimeFilterExponential::initialize(target);
    unFilteredQuantityOld_.reset(target.nRows(),target.nCols());
    unFilteredQuantityOld_ = target;
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterExplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterExplicit::TimeFilterExplicit(TimeFilterManager & timeFilterManager) :
    TimeFilterExponential(timeFilterManager)
  {
    // do nothing
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterImplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterImplicit::TimeFilterImplicit(TimeFilterManager & timeFilterManager) :
    TimeFilterExponential(timeFilterManager)
  {
    // do nothing
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterImplicitExplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterImplicitExplicit::TimeFilterImplicitExplicit(TimeFilterManager & timeFilterManager) :
    TimeFilterExponential(timeFilterManager)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterExplicitImplicit
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterExplicitImplicit::TimeFilterExplicitImplicit(TimeFilterManager & timeFilterManager) :
    TimeFilterExponential(timeFilterManager)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterImplicitUpdate
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterImplicitUpdate::TimeFilterImplicitUpdate(TimeFilterManager & timeFilterManager) :
    TimeFilterExponential(timeFilterManager)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TimeFilterStep
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  TimeFilterStep::TimeFilterStep(TimeFilterManager & timeFilterManager) :
    TimeFilter(timeFilterManager),
    elapsedTime_(0.0)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  initialize
  //    resets memory as needed
  //--------------------------------------------------------
  void TimeFilterStep::initialize(const MATRIX & target)
  {
    TimeFilter::initialize(target);
    unFilteredQuantityOld_.reset(target.nRows(),target.nCols());
    unFilteredQuantityOld_ = target;
  }


};

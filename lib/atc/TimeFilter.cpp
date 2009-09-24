// ATC_Transfer headers
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
    TimeFilterManager::TimeFilterManager(ATC_Transfer * atcTransfer) :
        atcTransfer_(atcTransfer),
        filterType_(NO_FILTER),
        filterScale_(0.),
        useFilter_(false),
        equilibrateFilter_(false),
        endEquilibrate_(false),
        needReset_(true)
    {
        // do nothing
    }
  
    //--------------------------------------------------------
    //  modify
    //    parses input commands
    //--------------------------------------------------------
    bool TimeFilterManager::modify(int narg, char ** arg)
    {
        bool foundMatch = false;
  
	// filter scale size
        /*! \page man_filter_scale fix_modify AtC transfer filter scale
          \section syntax
          fix_modify AtC transfer filter scale <scale> \n
          - scale (real) = characteristic time scale of the filter  \n

          \section examples
          <TT> fix_modify AtC transfer filter scale 10.0 </TT> \n

          \section description
	  Filters the MD dynamics to construct a more appropriate continuous field.  Equilibrating first filters the time derivatives without changing the dynamics to provide a better initial condition to the filtered dynamics
          \section restrictions
          only for be used with specific transfers:
	  thermal, two_temperature

          \section related
	  \ref man_time_filter

          \section default
          0.
        */
        if (strcmp(arg[0],"scale")==0) {
            filterScale_ = atof(arg[1]);
            if (filterScale_<=0.)
                throw ATC_Error(0,"Bad filtering time scale");
            foundMatch = true;
        }

        // time filtering activation switch
        /*! \page man_time_filter fix_modify AtC transfer filter
          \section syntax
          fix_modify AtC transfer filter <on | off | equilibrate> \n
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
            if (filterScale_<=0.)
                throw ATC_Error(0,"Filtering time scale not initialized");
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
            if (filterScale_<=0.)
                throw ATC_Error(0,"Filtering time scale not initialized");
            equilibrateFilter_ = true;
            endEquilibrate_ = false;
            useFilter_ = false;
            needReset_ = true;
            foundMatch = true;
        }

        /** Example command for to set time filter scale:
            fix_modify atc transfer filter type step*/
        else if (strcmp(arg[0],"type")==0) { 
            if (strcmp(arg[1],"exponential")==0) { 
                filterType_ = EXPONENTIAL_FILTER;
            }
            else if (strcmp(arg[1],"no_filter")==0) { 
                filterType_ = NO_FILTER;
            }
            else throw ATC_Error(0,"Not a supported time filter type");
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
    // NOTE currently the caller is responsible for deletion
    TimeFilter * TimeFilterManager::construct(const FilterIntegrationType intType)
    {
        if (useFilter_ || equilibrateFilter_) {
            if      (filterType_ == EXPONENTIAL_FILTER) {
                if (intType == IMPLICIT_EXPLICIT) {
                    return new TimeFilterImplicitExplicit(*this);     
                }
                else if (intType == EXPLICIT_IMPLICIT) {
                    return new TimeFilterExplicitImplicit(*this);     
                }
                else if (intType == EXPLICIT) {
                    return new TimeFilterExplicit(*this);     
                }
                else if (intType == IMPLICIT) {
                    return new TimeFilterImplicit(*this);     
                }
                else if (intType == IMPLICIT_UPDATE) {
                    return new TimeFilterImplicitUpdate(*this);     
                }
                else {
                    return new TimeFilterCrankNicolson(*this);     
                }
            }
        }
    
        // default to return base class
        return new TimeFilter(*this);     
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
        atcTransfer_(timeFilterManager.get_atc_transfer()),
        filterScale_(timeFilterManager.get_filter_scale())
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
    //  initialize
    //    resets memory as needed
    //--------------------------------------------------------
    void TimeFilterExponential::initialize(const MATRIX & target)
    {
        TimeFilter::initialize(target);
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
        TimeFilter(timeFilterManager)
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
    }


};

// ATC transfer headers
#include "PaqAtcUtility.h"
#include "ATC_Method.h"
#include "FE_Engine.h"

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class PaqAtcUtility
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  PaqAtcUtility::PaqAtcUtility(ATC_Method * atc,
                               AtomType atomType) :
    atc_(atc),
    atomType_(atomType),
    myNlocal(nullptr)
  {
    switch (atomType_) {
      case ALL:
        myNlocal = &ATC_Method::nlocal_total;
        break;
      case INTERNAL:
        myNlocal = &ATC_Method::nlocal;
        break;
      case GHOST:
        myNlocal = &ATC_Method::nlocal_ghost;
        break;
      case PROC_GHOST:
        myNlocal = &ATC_Method::nproc_ghost;
        break;
      default: // default cases to avoid compiler warnings
        break;
    }
  }

  //--------------------------------------------------------
  //  nlocal
  //--------------------------------------------------------
  int PaqAtcUtility::nlocal() const
  {
    return (atc_->*myNlocal)();
  }

  //--------------------------------------------------------
  //  nlocal_total
  //--------------------------------------------------------
  int PaqAtcUtility::nlocal_total() const
  {
    return atc_->nlocal_total();
  }

  //--------------------------------------------------------
  //  nsd
  //--------------------------------------------------------
  int PaqAtcUtility::nsd() const
  {
    return atc_->nsd();
  }

  //--------------------------------------------------------
  //  dt
  //--------------------------------------------------------
  double PaqAtcUtility::dt() const
  {
    return atc_->dt();
  }

  //--------------------------------------------------------
  //  atc_to_lammps_map
  //--------------------------------------------------------
  const Array<int> & PaqAtcUtility::atc_to_lammps_map() const
  {
    return (atomType_==INTERNAL) ? (atc_->internal_to_atom_map()) : (atc_->ghost_to_atom_map());
  }

  //--------------------------------------------------------
  //  fe_engine
  //--------------------------------------------------------
  const FE_Engine * PaqAtcUtility::fe_engine() const
  {
    return atc_->fe_engine();
  }

};

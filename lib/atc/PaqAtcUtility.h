// A class for defining atomic information used by the per-atom quantities
// this is required so they can be templatized

#ifndef PAQ_ATC_UTILITY_H
#define PAQ_ATC_UTILITY_H

#include "ATC_TypeDefs.h"
#include "Array.h"
#include "ATC_Error.h"

namespace ATC {

  // forward declarations
  class ATC_Method;
  class FE_Engine;

  class PaqAtcUtility {

  public:

    // constructor
    PaqAtcUtility(ATC_Method * atc, AtomType atomType=NO_ATOMS);

    // destructor
    ~PaqAtcUtility() {};

    // get the number of local atoms associated with this type
    int nlocal() const;

    // get the total number of atoms on this processor
    int nlocal_total() const;

    /** get the number of spatial dimensions */
    int nsd() const;

    /** get the timestep size */
    double dt() const;

    /** get the map from atc to lammps indexing for this type */
    const Array<int> & atc_to_lammps_map() const;

    /** access to the engine */
    const FE_Engine * fe_engine() const;

  protected:

    /** pointer to atc object */
    ATC_Method * atc_;

    /** type of atoms this quantity applies to */
    AtomType atomType_;

    /** function pointer to number of local atoms for this quantity */
    int (ATC_Method::*myNlocal)() const;

  private:

    // do not define
    PaqAtcUtility();

  };

};

#endif

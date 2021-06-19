/* ----------------------------------------------------------------------

 Coupling operators.
  
 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_COUPLINGOPERATOR_H
#define SELM_COUPLINGOPERATOR_H

//#include "driver_selm.h"
//#include "lammps.h"

#include "SELM_Lagrangian.h"
#include "SELM_Eulerian.h"

namespace LAMMPS_NS {

class LAMMPS;
class DriverSELM;

class SELM_CouplingOperator {

 public:

  SELM_CouplingOperator();
  SELM_CouplingOperator(int, char **);

  virtual ~SELM_CouplingOperator();

  virtual void setGlobalRefs(LAMMPS *lmps, DriverSELM *fx);

  virtual void computeOperator(const char        *operatorName,
                               SELM_Lagrangian   *SELM_LagrangianData,
                               SELM_Eulerian     *SELM_EulerianData) = 0;

  virtual void parse_ParameterFile(const char *filename) = 0;
  virtual void setup() = 0;

  virtual void writeSimulationDataToDisk(const char *baseFilename, int timeIndex) = 0;

  char    nameStr[1000]; /* name associated with this instance of the object */
  char    typeStr[1000]; /* shows type of class as string */
  int     type;          /* key for type of class */

  DriverSELM *driverSELM;
  LAMMPS     *lammps;

  int     flagWriteSimulationData; /* indicates simulation to be written to disk by this object */
  int     saveSkipSimulationData;  /* Number of time steps to skip between saving the simulation data*/

 protected:

 private:

};

}

#endif


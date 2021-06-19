/* ----------------------------------------------------------------------

 Lagrangian mechanics abstract.
 
 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_LAGRANGIAN_H
#define SELM_LAGRANGIAN_H

#include "SELM_Package.h"

namespace LAMMPS_NS {

class DriverSELM;
class LAMMPS;

class SELM_Lagrangian {

 public:

  /* ============== Function Prototypes ============== */
  SELM_Lagrangian();
  SELM_Lagrangian(int, char **);
  SELM_Lagrangian(class LAMMPS *lmps, class DriverSELM *fx);

  virtual ~SELM_Lagrangian();

  virtual void setGlobalRefs(LAMMPS *lmps, DriverSELM *fx);

  //virtual void parse_ParameterFile(char *filename) = 0;

  virtual void setup() = 0;

  virtual void setControlPtsDataFromLammpsData() = 0;
  virtual void setLammpsDataFromControlPtsData() = 0;

  virtual void writeSimulationDataToDisk(const char *filename, int timeIndex) = 0;

  /* ============== Variables ============== */
  char     nameStr[1000]; /* name associated with this instance of the object */
  char     typeStr[1000]; /* shows type of class as string */
  int      type;          /* key for type of class */

  LAMMPS     *lammps;
  DriverSELM *driverSELM;

  int      flagWriteSimulationData; /* determines if any simulation data is written to disk */
  int      saveSkipSimulationData;  /* Number of time steps to skip between saving the simulation data */

 protected:

 private:

};

}

#endif

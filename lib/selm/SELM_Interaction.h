/* ----------------------------------------------------------------------
 Custom interaction abstract.

 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_INTERACTION_H
#define SELM_INTERACTION_H

//#include "driver_selm.h"
//#include "lammps.h"

namespace LAMMPS_NS {

class DriverSELM;
class LAMMPS;

class SELM_Interaction {

 public:
  SELM_Interaction();
  SELM_Interaction(int, char **);
  virtual ~SELM_Interaction();

  virtual void setGlobalRefs(LAMMPS *lmps, DriverSELM *fx);

  virtual void parse_ParameterFile(const char *filename) = 0;
  virtual void setup() = 0;

  virtual void computeForceAndEnergy() = 0;

  virtual void writeSimulationDataToDisk(const char *baseFilename, int timeIndex) = 0;

  char nameStr[1000]; /* name associated with this instance of the object */
  char typeStr[1000]; /* shows type of class as string */
  int  type;          /* key for type of class */

  LAMMPS  *lammps;    /* reference to lammps */
  DriverSELM *driverSELM;   /* reference to SELM related fix routines */

  int  flagWriteSimulationData; /* indicates if simulation to be written by this object */
  int  saveSkipSimulationData;  /* number of time steps to skip between writing data */

 protected:

 private:

};

}

#endif

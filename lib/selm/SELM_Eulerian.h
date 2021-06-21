/* ----------------------------------------------------------------------

 Eulerian mechanics abstract type.
 
 Paul J. Atzberger
 http://atzberger.org/

------------------------------------------------------------------------- */

#ifndef SELM_EULERIAN_H
#define SELM_EULERIAN_H

namespace LAMMPS_NS {

class LAMMPS;
class DriverSELM;

class SELM_Eulerian {

 public:
  SELM_Eulerian();
  SELM_Eulerian(int, char **);

  virtual ~SELM_Eulerian();

  virtual void setGlobalRefs(LAMMPS *lmps, DriverSELM *fx);

  //virtual void parse_ParameterFile(char *filename) = 0;
  virtual void setup() = 0;

  virtual void writeSimulationDataToDisk(const char *baseFilename, int timeIndex) = 0;


  char nameStr[1000]; /* name associated with this instance of the object */
  char typeStr[1000]; /* shows type of class as string */
  int  type;          /* key for type of class */

  LAMMPS  *lammps;
  DriverSELM *driverSELM;

  int  flagWriteSimulationData; /* indicates if this object writes simulation data */
  int  saveSkipSimulationData;  /* number of time steps to skip between writing data to disk */

 protected:

 private:

};

}

#endif

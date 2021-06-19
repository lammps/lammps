/* ----------------------------------------------------------------------

 Integrator abstract.

 Paul J. Atzberger
 http://atzberger.org/
 
 ------------------------------------------------------------------------- */

#ifndef SELM_INTEGRATOR_H
#define SELM_INTEGRATOR_H

//#include "driver_selm.h"
//#include "lammps.h"

namespace LAMMPS_NS {

class DriverSELM;
class LAMMPS;

class SELM_Integrator {

 public:
  SELM_Integrator();
  SELM_Integrator(int, char **);
  virtual ~SELM_Integrator();

  virtual void setGlobalRefs(LAMMPS *lmps, DriverSELM *fx);

  virtual void parse_ParameterFile(char *filename) = 0;
  virtual void setup() = 0;
  virtual void setup_LAMMPS(int vflag);
  virtual void post_force(int vflag);
  virtual void set_LAMMPS_mask(int *mask_ptr) = 0; /* determine if both initial and final integrate are needed */

  virtual void integrate_initialize() = 0;
  virtual void integrate_initial() = 0;
  virtual void integrate_final() = 0;
  
  virtual void writeSimulationDataToDisk(char *baseFilename, int timeIndex) = 0;
  
  // optional methods (chosen to faciliate select integrator types, such as shear)
  virtual void pre_exchange() {} // default is to do nothing, sub-class can implement if needed
  virtual void end_of_step() {} // default is to do nothing, sub-class can implement if needed
  virtual void init_from_fix() {} // trigger for when fix->init() called, do this for the integrator
  //virtual int setmask() {};  // this is handled by set_LAMMPS_mask() above
      
  // additional variables
  char nameStr[1000]; /* name associated with this instance of the object */
  char typeStr[1000]; /* shows type of class as string */
  int  type;          /* key for type of class */

  LAMMPS  *lammps;    /* reference to lammps */
  //DriverSELM *driver_selm;   /* reference to SELM related fix routines */
  DriverSELM *driver_selm; /* we now use instead driver_selm to handle the data */

  int  flagWriteSimulationData; /* indicates if simulation to be written by this object */
  int  saveSkipSimulationData;  /* number of time steps to skip between writing data */

 protected:

 private:

};

}

#endif

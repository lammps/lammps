/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods (SELMs) Package

  Paul J. Atzberger
  http://atzberger.org/
  
  Please cite the follow paper when referencing this package
  
  "Fluctuating Hydrodynamics Methods for Dynamic Coarse-Grained Implicit-Solvent Simulations in LAMMPS," 
  Wang, Y. and Sigurdsson, J. K. and Atzberger, P. J., SIAM Journal on Scientific Computing, 38(5), 2016.
  
  @article{atz_selm_lammps_fluct_hydro,
    title = {Fluctuating Hydrodynamics Methods for Dynamic 
    Coarse-Grained Implicit-Solvent Simulations in LAMMPS},
    author = {Wang, Y. and Sigurdsson, J. K. and Atzberger, P. J.},
    journal = {SIAM Journal on Scientific Computing},
    volume = {38},
    number = {5},
    pages = {S62-S77},
    year = {2016},
    doi = {10.1137/15M1026390},
    URL = {https://doi.org/10.1137/15M1026390},
  }  
    
  For latest releases, examples, and additional information see 
  http://mango-selm.org/
 
------------------------------------------------------------------------- */

#ifndef LMP_DRIVER_SELM_H
#define LMP_DRIVER_SELM_H

#ifndef FFT_FFTW3 /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

#include "Atz_XML_Package.h"

#include "SELM_Eulerian.h"
#include "SELM_Integrator.h"
#include "SELM_Interaction.h"
#include "SELM_Lagrangian.h"
#include "SELM_CouplingOperator.h"

#include "driver_SELM_XML_Handler.h"

// forward declaration
//class SELM_Eulerian;
//class SELM_Lagrangian;
//class SELM_CouplingOperator;
//class SELM_Interaction;
//class SELM_Integrator;

#ifdef USE_PACKAGE_FFTW3
  #include <fftw3.h>
#endif
#include "fix_selm.h"

//#include "mpi.h"
#include "fix.h"
#include "lammps.h"
#include "random_mars.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstddef>

using namespace std;

namespace LAMMPS_NS {

struct DriverSELM {

 public:

  const char *error_str_code;

  /* ================= constants ================= */
  /* constants (values for non-integers defined in .cpp) */
  const int MAX_STR_LEN;

  /* parameter file types */
  const int    PARAM_FILE_TYPE_NULL;
  const int    PARAM_FILE_TYPE_TXT;
  const int    PARAM_FILE_TYPE_XML;

  /* ================= function prototypes ================= */
  DriverSELM(class FixSELM *, class LAMMPS *, int, char **);

  DriverSELM(); /* used to construct empty object,
                primarily for testing purposes */
  ~DriverSELM();

  // =========================== Fix Related Function Calls ==================== */
  int          setmask();
  virtual void init();
  void         setup(int vflag);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void         reset_dt();
  void         post_force(int vflag);
  
  // additional methods/hooks
  void pre_exchange();
  void end_of_step();

  /* =========================== SELM Function Calls =========================== */
  void init_attributes();
  void init_from_fix();
 
  /* Generates messages when encountering errors. */
  void packageError(int code, void *extras);

  /* SELM parse parameter file */
  void SELM_parse_ParameterFile_TXT(char *filename);

  /* SELM parse parameter file */
  void SELM_parse_ParameterFile_XML(char *filename);

  /* write the SELM simulation data to disk */
  void writeAllSimulationData(int timeIndex);
  void writeInfo();
  void writeFinalInfo();
  const string currentDateTime();

  /* =========================== Variables =========================== */
  LAMMPS                                                    *lammps; /* lammps data */

  friend class FixSELM;    
  FixSELM                                                   *fixSELM; /* reference to the fix directly */  // friend class so can access all members
  
  int                                                        SELM_integrator_mask; /* mask set by SELM codes */

  int                                                        SELM_Version;
  string                                                     SELM_SVN_Version;
  string                                                     SELM_Compile_Date_Time;
  string                                                     SELM_Run_Description;

  char                                                      *SELM_BasePath;
  char                                                      *SELM_dir_sim_data;
  char                                                      *SELM_BaseFilename;

  int                                                        SELM_Seed;
  RanMars                                                   *random; /* random number generator */

  int                                                        SELM_Eulerian_List_N;
  SELM_Eulerian                                            **SELM_Eulerian_List;

  int                                                        SELM_Lagrangian_List_N;
  SELM_Lagrangian                                          **SELM_Lagrangian_List;

  int                                                        SELM_CouplingOperator_List_N;
  SELM_CouplingOperator                                    **SELM_CouplingOperator_List;

  int                                                        SELM_Interaction_List_N;
  SELM_Interaction                                         **SELM_Interaction_List;

  SELM_Integrator                                           *SELM_IntegratorData;

};

}

#endif


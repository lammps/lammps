/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods (SELMs) Package (library version)

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
 
------------------------------------------------------------------------- 
*/

/* SELM_includes */
#include "driver_selm.h"
#include "driver_SELM_XML_Handler.h"

#include "Atz_XML_Package.h"
#include "SELM_Eulerian.h"
#include "SELM_Lagrangian.h"
#include "SELM_CouplingOperator.h"
#include "SELM_Interaction.h"
#include "SELM_Integrator.h"

#include "SELM_Parser1.h"

/* LAMMPS includes */
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "comm.h"
#include "universe.h"
#include "version.h" 
#include "random_mars.h"
#include "citeme.h"

/* include these to help trigger recompile each time */
/*
#include "style_angle.h"
#include "style_atom.h"
#include "style_bond.h"
#include "style_command.h"
#include "style_compute.h"
#include "style_dihedral.h"
#include "style_dump.h"
#include "style_fix.h"
#include "style_improper.h"
#include "style_integrate.h"
#include "style_kspace.h"
#include "style_minimize.h"
#include "style_pair.h"
#include "style_region.h"
*/

/* C/C++ includes */
#include <cstdlib>
#include <cstddef>
#include <cstdio>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* =========================== Class definitions =========================== */
DriverSELM::DriverSELM() : error_str_code("fix_selm.cpp"), MAX_STR_LEN(10000), PARAM_FILE_TYPE_NULL(0), PARAM_FILE_TYPE_TXT(1), PARAM_FILE_TYPE_XML(2) {
  /* WARNING: May need to modify LAMMPS codes so that we have
              Fix(NULL, 0, NULL) acts like empty constructor */

  const char *error_str_func = "DriverSELM()";

  int flagWarningMessage = 1;
  
  /* init constants */
  //error_str_code = "fix_selm.cpp";

  /* ================= constants ================= */
  /* constants (values for non-integers defined in .cpp) */
  //MAX_STR_LEN = 10000;

  /* parameter file types */
  //PARAM_FILE_TYPE_NULL = 0;
  //PARAM_FILE_TYPE_TXT  = 1;
  //PARAM_FILE_TYPE_XML  = 2;
  
  /* initialize select attributes */
  init_attributes();

  /* create empty class */
  lammps = NULL;
  fixSELM = NULL;

  stringstream message;
  message << "Empty DriverSELM created. This should only be used" << endl;
  message << "for testing purposes.  This object does not contain" << endl;
  message << "the needed LAMMPS data structure references." << endl;
  message << "Comment out the error generation if you really want to use this." << endl;
  //SELM_Package::packageWarning(error_str_code, error_str_func, message);
  SELM_Package::packageError(error_str_code, error_str_func, message);

}

/* =========================== Class definitions =========================== */
DriverSELM::DriverSELM(FixSELM *fixSELM, LAMMPS *lmp, int narg, char **arg) : error_str_code("fix_selm.cpp"), MAX_STR_LEN(10000), PARAM_FILE_TYPE_NULL(0), PARAM_FILE_TYPE_TXT(1), PARAM_FILE_TYPE_XML(2) {

  const char *error_str_func = "DriverSELM()";

  char tempStr[1000];
  int  paramFileType;
  int  N;
  
  /* ================= constants ================= */
  /* constants (values for non-integers defined in .cpp) */
  //MAX_STR_LEN = 10000;

  /* parameter file types */
  //PARAM_FILE_TYPE_NULL = 0;
  //PARAM_FILE_TYPE_TXT  = 1;
  //PARAM_FILE_TYPE_XML  = 2;
  
  /* ================= init ================= */

  this->fixSELM = fixSELM;
      
  /* initialize select attributes */
  init_attributes();

  /* == Display some version information */
  cout << endl;
  cout << "USER-SELM Package (SVN Version = " << SELM_SVN_Version;
  cout << ", Compile Date-Time = " << SELM_Compile_Date_Time << ") " << endl;
  cout << "For examples and additional information see http://mango-selm.org/" << endl;
  cout << "Simulation Start Date-Time = " << currentDateTime() << endl;

  /* == Save reference to the lammps data */
  lammps = lmp;
  SELM_Package::setLAMMPS(lammps);

  /* == Setup LAMMPS related flags. */
  fixSELM->time_integrate = 1; /* set to 1 for fix performing integration, 0 if fix does not */

  SELM_integrator_mask = 0; /* initial value for the mask for the time step integrator */

  /* == Parse options for the SELM integrator. */
  if (narg < 4) {
    stringstream message;
    message << "Fix SELM requires filename for parameters." <<  endl;
    SELM_Package::packageError(this->error_str_code, error_str_func, message);
  }

  paramFileType = PARAM_FILE_TYPE_XML;
  switch (paramFileType) {

    case 1: //PARAM_FILE_TYPE_TXT:
      /* === Parser for the parameter files for SELM. */
      SELM_parse_ParameterFile_TXT(arg[3]); /* parse the parameter file */
    break;

    case 2: //PARAM_FILE_TYPE_XML:
      /* === Parser XML File for the parameters for SELM. */
      SELM_parse_ParameterFile_XML(arg[3]); /* parse the parameter file */
    break;

    default:
      stringstream message;
      message << "The specified parameter file type is not recognized." << endl;
      message << "paramFileType = " << paramFileType << endl;
      SELM_Package::packageError(this->error_str_code, error_str_func, message);
      break;

  } /* end of switch */

  /* === Setup SELM data structures (parse additional files).
   */

  /* == Setup the random number generator */
  //random = new RanMars(lammps, SELM_Seed + fixSELM->comm->me);
  random = new RanMars(lammps, SELM_Seed);  // WARNING: Be careful for MPI, since each seed the same (use comm->me) 

  /* == Setup the Lagrangian data structures */
  N = SELM_Lagrangian_List_N;
  for (int I = 0; I < N; I++) {
    SELM_Lagrangian_List[I]->setup();
  }

  /* == Setup the Eulerian data structures */
  N = SELM_Eulerian_List_N;
  for (int I = 0; I < N; I++) {
    SELM_Eulerian_List[I]->setup();
  }

  /* == Setup the CouplingOperators data structures */
  N = SELM_CouplingOperator_List_N;
  for (int I = 0; I < N; I++) {
    SELM_CouplingOperator_List[I]->setup();
  }

  /* == Setup the CouplingOperators data structures */
  N = SELM_Interaction_List_N;
  for (int I = 0; I < N; I++) {
    SELM_Interaction_List[I]->setup();
  }

  /* == Setup the Integrator data structures */
  this->SELM_IntegratorData->setup();

  /* == Write the initial integration data */
  writeAllSimulationData(lammps->update->ntimestep);


  /* == Write some data to XML info file */
  /* record the LAMMPS version and SELM versions
   * among other information */
  writeInfo();

}

/* destructor */
DriverSELM::~DriverSELM() {

  /* write some final information to disk */
  writeFinalInfo();

}

void DriverSELM::setup(int vflag) { /* lammps setup */

  SELM_IntegratorData->setup_LAMMPS(vflag);

}

/* ---------------------------------------------------------------------- */
int DriverSELM::setmask()
{
  /*
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
   */

  SELM_integrator_mask = 0;

  /* pass to integrator to handle */
  SELM_IntegratorData->set_LAMMPS_mask(&SELM_integrator_mask);

  return SELM_integrator_mask; /* syncronize the SELM mask with that returned to LAMMPS */
}

/* ---------------------------------------------------------------------- */
void DriverSELM::pre_exchange()
{

  /* pass to integrator to handle */
  SELM_IntegratorData->pre_exchange();

}

/* ---------------------------------------------------------------------- */
void DriverSELM::end_of_step()
{

  /* pass to integrator to handle */
  SELM_IntegratorData->end_of_step();

}

/* ---------------------------------------------------------------------- */

void DriverSELM::init()
{

  const char *error_str_func = "init()";


  /* == Initialize the SELM integrators. */

  /* update->integrate->step; (time step from LAMMPS) */

  /* == Check the integration style is Verlet, if not report error. */
  if (strcmp(lammps->update->integrate_style, "verlet") != 0) {
    stringstream message;
    message << "SELM requires for now use of the verlet integrate_style." <<  endl;
    SELM_Package::packageError(this->error_str_code, error_str_func, message);
  }

  /* == Initialize data structures for SELM. */
  /* integrator trigger fix_init() for any associated initialization */
  SELM_IntegratorData->init_from_fix();

}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */
//void DriverSELM::integrate_initialize() {

//  SELM_IntegratorData->integrate_initialize();

//}


void DriverSELM::initial_integrate(int vflag)
{

  /* == Get the current information about the box deformation and shear. */


  /* == Get the current information about the particles. */


  /* == Determine the hydrodynamic and particle forces. */

  /* == Perform the update of the fluid and particles */
  SELM_IntegratorData->integrate_initial();

  /* if only initial called then write data to disk here,
   * otherwise wait for the final integration step to be
   * called. */
  if (SELM_integrator_mask & FINAL_INTEGRATE == 0) {
    int timeIndex = lammps->update->ntimestep;
    writeAllSimulationData(timeIndex);
  } else {
    /* wait until integrate_final is called */
  }

}

/* ---------------------------------------------------------------------- */

void DriverSELM::final_integrate()
{

  /* == Perform the update of the Eulerian-Lagrangian system */
  SELM_IntegratorData->integrate_final();

  /* == Write the simulation data to disk */
  if (SELM_integrator_mask & FINAL_INTEGRATE == 0) {
    /* do nothing, should not encounter this case. */
    /* maybe report error. */
  } else {
    int timeIndex = lammps->update->ntimestep;
    writeAllSimulationData(timeIndex);
  }

}

/* ---------------------------------------------------------------------- */

void DriverSELM::reset_dt()
{

  const char *error_str_func = "reset_dt()";

  stringstream message;
  message << "The fix_SELM does not implement reset_dt() yet." <<  endl;
  SELM_Package::packageError(this->error_str_code, error_str_func, message);

}

void DriverSELM::post_force(int vflag) {
  SELM_IntegratorData->post_force(vflag);
}


/*****************************************************************************************/
/* Supporting SELM codes */
/*****************************************************************************************/
void DriverSELM::init_attributes() {

  #ifndef SVN_REV
    SELM_SVN_Version        = "(was not available)"; /* indicates not SVN version info. available */
  #else
    SELM_SVN_Version        = SVN_REV;
  #endif

  #ifndef COMPILE_DATE_TIME
    SELM_Compile_Date_Time  = "(was not available)"; /* indicates date and time info. not available */
  #else
    SELM_Compile_Date_Time  = COMPILE_DATE_TIME;
  #endif

  SELM_Run_Description    = "";

  SELM_BasePath           = NULL;
  SELM_dir_sim_data       = NULL;  
  SELM_BaseFilename       = NULL;

}

void DriverSELM::packageError(int code, void *extras) {
    std::exit(code);
}

/* parses the parameters from a file */
void DriverSELM::SELM_parse_ParameterFile_XML(char *filename) {

  const char *error_str_func = "SELM_parse_ParameterFile_XML()";

  DriverSELM *DriverSELM_Data = NULL; /* for testing */

  Driver_SELM_XML_Handler *driver_SELM_DataHandler
    = new Driver_SELM_XML_Handler(this);

  Atz_XML_SAX_Handler_Multilevel *dataHandler
    = new Atz_XML_SAX_Handler_Multilevel(driver_SELM_DataHandler);

  Atz_XML_Parser::parse(filename, dataHandler);

  DriverSELM_Data = (DriverSELM *) driver_SELM_DataHandler->XML_getData();

}

/* parses the parameters from a file */
void DriverSELM::SELM_parse_ParameterFile_TXT(char *filename) {

  const char *error_str_code = "fix_SELM.c";
  const char *error_str_func = "SELM_parse_ParameterFile_TXT()";

}

void DriverSELM::writeAllSimulationData(int timeIndex) {
  
  char *dir_output = this->SELM_dir_sim_data;
  
  // write the Lagrangian DOF to disk
  for (int k = 0; k < SELM_Lagrangian_List_N; k++) {
    if ( (this->SELM_Lagrangian_List[k]->flagWriteSimulationData) &&
         (timeIndex % this->SELM_Lagrangian_List[k]->saveSkipSimulationData == 0) ) {
      // WARNING: changed from this->SELM_BaseFilename (might need to change downstream codes)   
      this->SELM_Lagrangian_List[k]->writeSimulationDataToDisk(dir_output, 
                                                               timeIndex);
    }
  } // end k loop

  // write the Eulerian DOF to disk
  for (int k = 0; k < SELM_Eulerian_List_N; k++) {
    if ( (this->SELM_Eulerian_List[k]->flagWriteSimulationData) &&
         (timeIndex % this->SELM_Eulerian_List[k]->saveSkipSimulationData == 0) ) {
      // WARNING: changed from this->SELM_BaseFilename (might need to change downstream codes)   
      this->SELM_Eulerian_List[k]->writeSimulationDataToDisk(dir_output, 
                                                             timeIndex);
    }
  } // end k loop

  // write the Integrator data to disk
  if ( (this->SELM_IntegratorData->flagWriteSimulationData) &&
       (timeIndex % this->SELM_IntegratorData->saveSkipSimulationData == 0) ) {
    // WARNING: changed from this->SELM_BaseFilename (might need to change downstream codes)   
    this->SELM_IntegratorData->writeSimulationDataToDisk(dir_output, 
                                                         timeIndex);
  } // end if

}

void DriverSELM::writeInfo() {

  char filename[10000];

  const char *error_str_func = "writeInfo()";

  FILE *fid;

  /* open the file for writing the data */
  sprintf(filename, "%s.SELM_Info", SELM_BaseFilename);
  fid = fopen(filename,"w");

  if (fid == NULL) {
    stringstream message;
    message << "Could not open file to write error occured." << endl;
    message << "  filename = " << filename << endl;
    SELM_Package::packageError(this->error_str_code, error_str_func, message);
  }

  /* write XML format */
  stringstream output;
  output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  output << "<SELM_Info>" << endl;
  //output << "<LAMMPS_Version value=" << "\"" <<  lammps->universe->version << "\"" << "/>" << endl;
  // PJA: note the new version of LAMMPS defines hard-code
  output << "<LAMMPS_Version value=" << "\"" << LAMMPS_VERSION << "\"" << "/>" << endl;  
  output << "<SELM_SVN_Version value=" << "\"" <<  SELM_SVN_Version << "\"" << "/>" << endl;
  output << "<SELM_Compile_Date_Time value=" << "\"" <<  SELM_Compile_Date_Time << "\"" << "/>" << endl;
  output << "<Simulation_Start_Date_Time value=" << "\"" <<  currentDateTime() << "\"" << "/>" << endl;
  output << "</SELM_Info>" << endl;

  /* write the output to disk */
  fprintf(fid,"%s",output.str().c_str());

   /* close the file */
  fclose(fid);

}

void DriverSELM::writeFinalInfo() {

  char filename[10000];

  const char *error_str_func = "writeInfo()";

  FILE *fid;

  /* open the file for writing the data */
  sprintf(filename, "%s.SELM_InfoExtra", SELM_BaseFilename);
  fid = fopen(filename,"w");

  if (fid == NULL) {
    stringstream message;
    message << "Could not open file to write error occured." << endl;
    message << "  filename = " << filename << endl;
    SELM_Package::packageError(this->error_str_code, error_str_func, message);
  }

  /* delete final line of current file */

  /* write some additional information */
  /* write XML format */
  stringstream output;
  output << "<SELM_InfoExtra>" << endl;
  output << "<Simulation_Stop_Date_Time value=" << "\"" <<  currentDateTime() << "\"" << "/>" << endl;
  output << "</SELM_InfoExtra>" << endl;

  /* write the output to disk */
  fprintf(fid,"%s",output.str().c_str());

  /* close the file */
  fclose(fid);

}

// Get current date/time, (format such as YYYY-MM-DD.HH:mm:ss)
const string DriverSELM::currentDateTime() {

  time_t     now = time(0);
  struct tm  tstruct;

  char       buf[80];

  tstruct = *localtime(&now);

  // Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
  // for more information about date/time format
  //strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
  strftime(buf, sizeof(buf), "%m-%d-%Y %X", &tstruct);

  return buf;
}

// pass along the initialization
void DriverSELM::init_from_fix() {
  SELM_IntegratorData->init_from_fix();
}




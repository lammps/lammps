/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#include <cstdlib>
#include "stdio.h"
#include "string.h"
#include "wrapper_selm.h"

#include "SELM_Parser1.h"

#include "SELM_Interaction.h"
#include "SELM_Interaction_CUSTOM1.h"


#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "lammps.h"
#include "domain.h"
#include "wrapper_selm.h"

#include <malloc.h>

#ifdef DYNAMIC_LINK_PLUG_IN /* determine if plug-in features enabled */
  #include <iostream>
  #include <dlfcn.h>
#endif


using namespace LAMMPS_NS;
using std::cout;
using std::cerr;

/* ---------------------------------------------------------------------- */
const int   SELM_Interaction_CUSTOM1::TYPE     = SELM_Interaction_Types::TYPE_CUSTOM1;
const char* SELM_Interaction_CUSTOM1::TYPE_STR = SELM_Interaction_Types::TYPE_STR_CUSTOM1;

const char *SELM_Interaction_CUSTOM1::error_str_code = "SELM_Interaction_CUSTOM1.cpp";

SELM_Interaction_CUSTOM1::SELM_Interaction_CUSTOM1() {

  /* WARNING: may need to call parent constructor */

  /* generic initialization */
  init();

}

SELM_Interaction_CUSTOM1::SELM_Interaction_CUSTOM1(int narg, char **arg) : SELM_Interaction(narg, arg) {

  const char *error_str_func = "SELM_Interaction_CUSTOM1()";

  /* generic initialization */
  init();

}


SELM_Interaction_CUSTOM1::SELM_Interaction_CUSTOM1(class LAMMPS *lmps, class DriverSELM *fx) {

  /* generic initialization */
  init();

  setGlobalRefs(lammps, driverSELM);

}

void SELM_Interaction_CUSTOM1::init() {

  type = SELM_Interaction_CUSTOM1::TYPE;
  strcpy(typeStr, SELM_Interaction_CUSTOM1::TYPE_STR);

  strcpy(nameStr, "No Name");

  setGlobalRefs(NULL,NULL);

}

void SELM_Interaction_CUSTOM1::init_library() {

  const char *error_str_func = "init_library()";

  #ifndef DYNAMIC_LINK_PLUG_IN /* determine if plug-in features enabled */
    stringstream message;
    message << "  Dynamic link plug-in functionality is not enabled. " << endl;
    message << "  To us this class codes must be re-compiled with DYNAMIC_LINK_PLUG_IN flag. " << endl;

    SELM_Package::packageWarning(error_str_code, error_str_func, message);

  #else

    void* libHandle = dlopen(libName.c_str(), RTLD_LAZY);

    if (!libHandle) { /* library_handle is null */
      stringstream message;

      message << "Trouble opening the dynamic link library : " << dlerror() << endl;
      message << "library_name = " << libName << endl;

      SELM_Package::packageError(error_str_code, error_str_func, message);
      return;
    }

    // load the symbol
    typedef void (*libInitType)();  /* ideally pass the parsed parameter data as Map */

    // reset errors
    dlerror();
    libInitType libInit = (libInitType) dlsym(libHandle, "init");
    const char *dlsym_error = dlerror();

    if (dlsym_error) {
      message << "Trouble loading symbol from the dynamic link library : " << dlsym_error << endl;
      message << "library_name = " << libName << endl;

      dlclose(libHandle);

      SELM_Package::packageError(error_str_code, error_str_func, message);
    }

    /* call the library initialization routine */
    libInit();

  #endif


}


void SELM_Interaction_CUSTOM1::parse_ParameterFile(const char *filename) {

}


void SELM_Interaction_CUSTOM1::setup() {

  const char *error_str_func = "setup()";

  /* load and initialize the dynamic link library */
  init_library();

}

void SELM_Interaction_CUSTOM1::computeForceAndEnergy() {

}

void SELM_Interaction_CUSTOM1::writeSimulationDataToDisk(const char *baseFilename, int timeIndex) {

}

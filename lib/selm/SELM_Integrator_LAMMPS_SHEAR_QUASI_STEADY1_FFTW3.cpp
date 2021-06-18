/*
--------------------------------------------------------------------------------
  SELM fluctuating hydrodynamics coupling with shear boundary conditions.
    
  Paul J. Atzberger
  http://atzberger.org/
  
  Please cite the following paper when using these methods
  
  "Incorporating shear into stochastic Eulerian-Lagrangian methods for 
  rheological studies of complex fluids and soft materials," 
  Paul J. Atzberger, Physica D: Nonlinear Phenomena, 265, (2013).
    
  @article{atzberger_selm_shear_2013,
    title={Incorporating shear into stochastic Eulerian-Lagrangian methods 
           for rheological studies of complex fluids and soft materials },
    author={Paul J. Atzberger},
    journal={Physica D: Nonlinear Phenomena},
    year={2013},
    pages={57 - 70},
    volume={265},  
    doi={http://doi.org/10.1016/j.physd.2013.09.002},
    issn={0167-2789},
  }
  
  For examples and additional information see http://mango-selm.org/

--------------------------------------------------------------------------------
*/

#ifndef FFT_FFTW /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

/* standard C/C++ includes */
#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include "string.h"

/* SELM includes */
#include "driver_selm.h"

#include "SELM_Parser1.h"

#include "SELM_Integrator.h"
#include "SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3.h"

#include "SELM_Eulerian.h"
#include "SELM_Eulerian_Types.h"
#include "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3.h"

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE.h"
#include "SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID.h"

#include "SELM_CouplingOperator.h"
#include "SELM_CouplingOperator_Types.h"
#include "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1.h"

/* LAMMPS includes */
//#include "respa.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "irregular.h"
#include "domain.h"
#include "lattice.h"
#include "force.h"
#include "modify.h"
#include "math_const.h"
#include "kspace.h"
#include "input.h"
#include "variable.h"
#include "error.h"

#include "citeme.h"

#ifdef USE_PACKAGE_FFTW3
  #include "fftw3.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

/* ============ constants definition =========== */
double SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::UNIT_pi  = 3.141592653589793;

const int    SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::TYPE     = SELM_Integrator_Types::TYPE_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3;
const char*  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::TYPE_STR = SELM_Integrator_Types::TYPE_STR_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3;

const char*  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_STR_NULL      = "NULL";
const char*  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_STR_RM_SHEAR1 = "RM_SHEAR1";
const char*  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_STR_RM_OSC1   = "RM_OSC1";

const char* SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::error_str_code = "SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3.cpp";

enum{NONE=0,FINAL,DELTA,SCALE,VEL,ERATE,TRATE,VOLUME,WIGGLE,VARIABLE};

static const char cite_selm_shear_str[] =
  "USER-SELM Package: Fluctuating Hydrodynamics \n\n"
  "@article{atzberger_selm_shear_2013,\n"
  "title={Incorporating shear into stochastic Eulerian-Lagrangian methods\n"
  "for rheological studies of complex fluids and soft materials },\n"
  "author={Paul J. Atzberger},\n"
  "journal={Physica D: Nonlinear Phenomena},\n"
  "year={2013},\n"
  "pages={57 - 70},\n"
  "volume={265},\n"  
  "doi={http://doi.org/10.1016/j.physd.2013.09.002},\n"
  "issn={0167-2789},\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */
SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3() : SELM_Integrator()
{

  init();

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras
      = (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType *)
        malloc(sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType));
        
  // set to zero
  memset(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras,
         0,sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType));        
        
  // will need to setup the rest of the way from data at later stage
  
}

SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3(int narg, char **arg) : SELM_Integrator(narg, arg)
{

  const char *error_str_func = "SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3()";

  /* generic initialization */
  init();
  
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras
      = (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType *)
        malloc(sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType));
        
  // init to zero
  memset(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras,
         0,sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType));                
        
  printf("Creating SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3(int narg, char **arg) : for this type of constructor not yet implemented. \n");
  
  //if (lmps->citeme) lmps->citeme->add(cite_selm_shear_str);            
  //SELM_Package::packageError();
  std::exit(1);
  
}

SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3(LAMMPS *lmps, DriverSELM *fx)
{

  /* generic initialization */
  init();
  
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras
      = (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType *)
        malloc(sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType));  
        
  // init to zero
  memset(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras,
         0,sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType));                

  setGlobalRefs(lmps, fx);
  
}


SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::~SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3()
{

  // -- from fix_deform.cpp
  if(set) {
    for (int i = 0; i < 6; i++) {
      delete [] set[i].hstr;
      delete [] set[i].hratestr;
    }
  }
  delete [] set;

  delete irregular;
 
  // reset domain's h_rate = 0.0, since this fix may have made it non-zero
  double *h_rate = lammps->domain->h_rate;
  double *h_ratelo = lammps->domain->h_ratelo;

  h_rate[0] = h_rate[1] = h_rate[2] =
  h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;
  // --
    
}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::setup_internals()
{

  const char *error_str_func = "setup_internals()";
  
  Atom *atom = lammps->atom;
  Domain *domain = lammps->domain;
  Update *update = lammps->update;

  int flagShearMode;
  
  int    shearDir;
  int    shearVelDir;
  double shearRate;
  double shearDist;
  double shearDist_last;

  //double shearOmega;
  //double shearRateAmplitude;
  
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3;
  //SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  ShearData_RM_SHEAR1_Type *shearData_RM_SHEAR1;
  ShearData_RM_OSC1_Type   *shearData_RM_OSC1;

  /* generic initialization */
  //init();
  //setGlobalRefs(lmps, fx);
     
  flagShearMode = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode;

  /* de-reference variables */
  switch (flagShearMode) {

    case SHEAR_MODE_TYPE_RM_SHEAR1:

      shearData_RM_SHEAR1 = (ShearData_RM_SHEAR1_Type *)
          SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

      /* retrieve shear values required to generate the flow */
      shearDir       = shearData_RM_SHEAR1->shearDir;
      shearVelDir    = shearData_RM_SHEAR1->shearVelDir;
      shearRate      = shearData_RM_SHEAR1->shearRate;
      shearDist      = shearData_RM_SHEAR1->shearDist;
      shearDist_last = shearData_RM_SHEAR1->shearDist_last;
      
      SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
        = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0];

      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
        = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

      /* setup the mesh in accordance of the current shear conditions */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir
        = shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir
        = shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate
        = shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
        = shearDist;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
        = shearDist_last;

      break;

    case SHEAR_MODE_TYPE_RM_OSC1:

      shearData_RM_OSC1
        = (ShearData_RM_OSC1_Type *)
        SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;
       
      /* update the shear rate (oscillation) */
      shearDir           = shearData_RM_OSC1->shearDir;
      shearVelDir        = shearData_RM_OSC1->shearVelDir;
      //shearOmega         = shearData_RM_OSC1->shearOmega;
      //shearRateAmplitude = shearData_RM_OSC1->shearRateAmplitude;
      shearRate          = shearData_RM_OSC1->shearRate;
      shearDist          = shearData_RM_OSC1->shearDist;
      shearDist_last     = shearData_RM_OSC1->shearDist_last;
      
      SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
        = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0];

      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
        = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;      

      /* setup the mesh in accordance of the current shear conditions */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir
        = shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir
        = shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate
        = shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
        = shearDist;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
        = shearDist_last;

      break;

    } /* end of switch */
    
  //printf("flagShearMode = %d",flagShearMode);
  
  /* create "set" data structure from "fix_deform.cpp" and initialize for bookkeeping. */  
  // create Set for storing box book-keeping information
  set = new Set[6];
  memset(set,0,6*sizeof(Set));  // PJA: sets all values to zero
  for (int i = 0; i < 6; i++) { // to be extra safe we set styles explicitly to NONE
    set[i].style = NONE;
  }
  
  /* -- start of "fix deform" create() related items */
  /* setup domain settings so signals other LAMMPS classes of modifications
     (we follow approaches used in "fix deform" class).
   */
  driver_selm->fixSELM->no_change_box = 1;
  driver_selm->fixSELM->restart_global = 1;
  driver_selm->fixSELM->pre_exchange_migrate = 1;

  // // read options from end of input line
  // // no x remap effectively moves atoms within box, so set restart_pbc
  // options(narg-iarg,&arg[iarg]);
  // from options()
  flipflag = 1; // flips allowed once box gets too big
  remapflag = Domain::V_REMAP; // remap the velocity of particles when bonds cross the boundary

  if (remapflag != Domain::X_REMAP) driver_selm->fixSELM->restart_pbc = 1;
  driver_selm->fixSELM->restart_pbc = 1;

  // setup dimflags used by other classes to check for volume-change conflicts
  for (int i = 0; i < 6; i++)
    if (set[i].style == NONE) dimflag[i] = 0;
    else dimflag[i] = 1;

  //if (dimflag[0]) box_change |= BOX_CHANGE_X;
  //if (dimflag[1]) box_change |= BOX_CHANGE_Y;
  //if (dimflag[2]) box_change |= BOX_CHANGE_Z;
  //if (dimflag[3]) box_change |= BOX_CHANGE_YZ;
  if (dimflag[4]) driver_selm->fixSELM->box_change |= driver_selm->fixSELM->BOX_CHANGE_XZ;
  //if (dimflag[5]) box_change |= BOX_CHANGE_XY;

  driver_selm->fixSELM->box_change |= driver_selm->fixSELM->BOX_CHANGE_XZ;
  
  // // reneighboring only forced if flips can occur due to shape changes
  //if (flipflag && (set[3].style || set[4].style || set[5].style))
  // force_reneighbor = 1;  
  //next_reneighbor = -1;
  
  driver_selm->fixSELM->force_reneighbor = 1;
  driver_selm->fixSELM->next_reneighbor = -1;
  
  flip = 0;
  
  //if (force_reneighbor) irregular = new Irregular(lmp);
  //else irregular = nullptr;  
  irregular = new Irregular(lammps);
          
  /* further initialization will occur from init_from_fix() call */
  
}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::init_from_fix()
{

  const char *error_str_func = "init_from_fix()";
  
  Atom *atom = lammps->atom;
  Domain *domain = lammps->domain;
  Update *update = lammps->update;

  int flagShearMode;
  
  int index;

  int    shearDir;
  int    shearVelDir;
  double shearRate;
  double shearDist;
  double shearDist_last;

  double shearOmega;
  double shearRateAmplitude;
  
  if (lammps->citeme) lammps->citeme->add(cite_selm_shear_str);
  
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3;
  //SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  ShearData_RM_SHEAR1_Type *shearData_RM_SHEAR1;
  ShearData_RM_OSC1_Type   *shearData_RM_OSC1;

  /* generic initialization */
  flagShearMode = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode;

  /* de-reference variables */
  /* PJA: Note: below is not currently really needed. */
  switch (flagShearMode) {

    case SHEAR_MODE_TYPE_RM_SHEAR1:
    
      SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
        = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0];

      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
        = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;    

      shearData_RM_SHEAR1 = (ShearData_RM_SHEAR1_Type *)
          SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

      /* retrieve shear values required to generate the flow */
      shearDir       = shearData_RM_SHEAR1->shearDir;
      shearVelDir    = shearData_RM_SHEAR1->shearVelDir;
      shearRate      = shearData_RM_SHEAR1->shearRate;
      shearDist      = shearData_RM_SHEAR1->shearDist;
      shearDist_last = shearData_RM_SHEAR1->shearDist_last;
      
      /* setup the mesh in accordance of the current shear conditions */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir
        = shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir
        = shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate
        = shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
        = shearDist;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
        = shearDist_last;

      break;

    case SHEAR_MODE_TYPE_RM_OSC1:
    
      SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
        = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0];

      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
        = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;        

      shearData_RM_OSC1
        = (ShearData_RM_OSC1_Type *)
        SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

      /* update the shear rate (oscillation) */
      shearDir           = shearData_RM_OSC1->shearDir;
      shearVelDir        = shearData_RM_OSC1->shearVelDir;      
      shearOmega         = shearData_RM_OSC1->shearOmega;
      shearRateAmplitude = shearData_RM_OSC1->shearRateAmplitude;
      shearRate          = shearData_RM_OSC1->shearRate;
      shearDist          = shearData_RM_OSC1->shearDist;
      shearDist_last     = shearData_RM_OSC1->shearDist_last;
      
      SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
        = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0];

      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
        = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

      /* setup the mesh in accordance of the current shear conditions */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir
        = shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir
        = shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate
        = shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
        = shearDist;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
        = shearDist_last;

      break;

    } /* end of switch */
    
  // setup values for set and related init  
  triclinic = domain->triclinic;  // required to be one or issue error
  if (triclinic == 0) {
    stringstream message; message << "requires a 'triclinic' domain." << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }
   
  index = 4; // xz shear
  set[index].style = VEL;
  //set[index].vel = utils::numeric(FLERR,arg[iarg+2],false,lmp);
  //set[index].vel = shearRate*(domain->boxhi[2] - domain->boxlo[2]); // atz_set the shear rate here for now (might need conversion later) 
  set[index].vel = shearRate*(domain->boxhi[2] - domain->boxlo[2]); // atz_set the shear rate here for now (might need conversion later) 
    
  //---- init()
  
  // Kspace setting  
  if (lammps->force->kspace) kspace_flag = 1;
  else kspace_flag = 0;
  
  double delt = (update->endstep - update->beginstep)*update->dt;

  // set start/stop values for box size and shape
  // if single run, start is current values
  // if multiple runs enabled via run start/stop settings,
  //   start is value when fix deform was issued
  // if VARIABLE or NONE, no need to set
  for (int i = 0; i < 3; i++) {
  
    if (update->firststep == update->beginstep) {
      set[i].lo_start = domain->boxlo[i];
      set[i].hi_start = domain->boxhi[i];
      set[i].vol_start = domain->xprd * domain->yprd * domain->zprd;
    } else {
      set[i].lo_start = set[i].lo_initial;
      set[i].hi_start = set[i].hi_initial;
      set[i].vol_start = set[i].vol_initial;
    }  
    
    if (set[i].style == VEL) {
      set[i].lo_stop = set[i].lo_start - 0.5*delt*set[i].vel;
      set[i].hi_stop = set[i].hi_start + 0.5*delt*set[i].vel;
    }

  }

  for (int i = 3; i < 6; i++) {
       
    if (update->firststep == update->beginstep) {
      if (i == 5) set[i].tilt_start = domain->xy;
      else if (i == 4) set[i].tilt_start = domain->xz;
      else if (i == 3) set[i].tilt_start = domain->yz;
    } else set[i].tilt_start = set[i].tilt_initial;
    
    if (set[i].style == VEL) {
      set[i].tilt_stop = set[i].tilt_start + delt*set[i].vel;
    }
    
  }

  /* init procedure */  
  // set domain->h_rate values for use by domain and other fixes/computes
  // initialize all rates to 0.0
  // cannot set here for TRATE,VOLUME,WIGGLE,VARIABLE since not constant
  h_rate = domain->h_rate;
  h_ratelo = domain->h_ratelo;
  
  for (int i = 0; i < 3; i++) {
    h_rate[i] = h_ratelo[i] = 0.0;
  }
    
  /* for xz we have */ 
  for (int i = 3; i < 6; i++) {
    h_rate[i] = 0.0;
    if (set[i].style == FINAL || set[i].style == DELTA ||
        set[i].style == VEL || set[i].style == ERATE) {
      if (delt != 0.0)
        h_rate[i] = (set[i].tilt_stop - set[i].tilt_start) / delt;
      else h_rate[i] = 0.0;
    }
  } 

  flip = 0; // init flag to zero for the initialization 

  /* -- end of fix deform init() related items */
  
  /* -- also added this code, since domain explicitly looks for fix_deform in domain.cpp */
  /*
  deform_flag = deform_vremap = deform_groupbit = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      deform_flag = 1;
      if (((FixDeform *) modify->fix[i])->remapflag == Domain::V_REMAP) {
        deform_vremap = 1;
        deform_groupbit = modify->fix[i]->groupbit;
      }
    }
  */

  domain->deform_flag = 1;
  domain->deform_vremap = 1;
  domain->deform_groupbit = driver_selm->fixSELM->groupbit;
  
  // region inits
  // for (int i = 0; i < nregion; i++) regions[i]->init();

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::pre_exchange()
{

  Atom *atom = lammps->atom;
  Domain *domain = lammps->domain;

  /* codes from "fix_deform.cpp" to help handle atom flipping */
  if (flip == 0) return;

  domain->yz = set[3].tilt_target = set[3].tilt_flip;
  domain->xz = set[4].tilt_target = set[4].tilt_flip;
  domain->xy = set[5].tilt_target = set[5].tilt_flip;
  
  // atz -- (added to move the bottom of the box to make shear symmetric, with the flip)
  double lx = domain->boxhi[0] - domain->boxlo[0];
  domain->boxlo[0] = set[0].lo_target = set[0].lo_start - 0.5*domain->xz;
  domain->boxhi[0] = set[0].hi_target = set[0].lo_target + lx;
  //set[0].lo_target = set[0].lo_start + delta*(set[0].lo_stop - set[0].lo_start);
  //set[0].hi_target = set[0].lo_target + lx;
  // --
  
  domain->set_global_box();
  domain->set_local_box();

  domain->image_flip(flipxy,flipxz,flipyz);

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  domain->x2lamda(atom->nlocal);
  irregular->migrate_atoms();
  domain->lamda2x(atom->nlocal);

  flip = 0;

  // PJA: Possibly update later the SELM codes with the flip information, if needed.
  
}


// handled above
//int SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::setmask()
//{
//  /* based on codes from "fix_deform.cpp" to help handle atom flipping */
//  int mask = 0;
//  if (force_reneighbor) mask |= PRE_EXCHANGE;
//  mask |= END_OF_STEP;
//  return mask;
//}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::end_of_step() {
  /* based on codes from "fix_deform.cpp", simplified for the shear xz case. */
  
  const char *error_str_func = "end_of_step()";
  
  Atom *atom = lammps->atom;
  Domain *domain = lammps->domain;  
  Update *update = lammps->update;
    
  ShearData_RM_SHEAR1_Type *shearData_RM_SHEAR1;
  ShearData_RM_OSC1_Type   *shearData_RM_OSC1;
  
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3                                                       *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;
      
  int flagShearMode = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode;
  
  /* de-reference variables */
  SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
    = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0]; // assumes type

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras 
    = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;
  
  if (flagShearMode == SHEAR_MODE_TYPE_RM_SHEAR1) {
      shearData_RM_SHEAR1 = (ShearData_RM_SHEAR1_Type *)
          SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;
  } else if (flagShearMode == SHEAR_MODE_TYPE_RM_OSC1) {
      shearData_RM_OSC1
        = (ShearData_RM_OSC1_Type *)
        SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;
  } else {
    stringstream message; message << "Unknown case of flagShearMode = " << flagShearMode << endl; 
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  // -- related to fix_deform.cpp --
  int i;  
  double delta = update->ntimestep - update->beginstep;
  
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  
  for (i = 0; i < 3; i++) {
    if (set[i].style == NONE) {
      set[i].lo_target = domain->boxlo[i];
      set[i].hi_target = domain->boxhi[i];
    }
  }
  
  // set new box size
  // for triclinic, set new box shape
  double *h = domain->h;
  h_rate = domain->h_rate;
  h_ratelo = domain->h_ratelo;
  
  if (triclinic) {
           
    for (i = 3; i < 6; i++) {

      if (set[i].style == NONE) {
        if (i == 5) set[i].tilt_target = domain->xy;
        else if (i == 4) set[i].tilt_target = domain->xz;
        else if (i == 3) set[i].tilt_target = domain->yz;
      } else if ((set[i].style == VEL) && (flagShearMode == SHEAR_MODE_TYPE_RM_SHEAR1)) { 
        // we cut other cases, this is what gets called for VEL
        set[i].tilt_target = set[i].tilt_start + delta*(set[i].tilt_stop - set[i].tilt_start);  
        // shear rate assumed constant, so nothing to update for SELM codes currently.       
      } else if ((set[i].style == VEL) && (flagShearMode == SHEAR_MODE_TYPE_RM_OSC1)) {

        double dt = update->dt; double t = update->ntimestep*dt;
        double omega = shearData_RM_OSC1->shearOmega;
        double A = shearData_RM_OSC1->shearRateAmplitude;
        double shearRate;

        set[i].tilt_target = set[i].tilt_start + A*sin(2.0*UNIT_pi*omega*t); 
        // note this is different than SELM Shear 2013 paper, 
        // we control now via displacement (allows omega = 0.0)
        
        double tilt_vel = (A*2.0*UNIT_pi*omega)*cos(2.0*UNIT_pi*omega*t);
        double lz = domain->boxhi[2] - domain->boxlo[2];
        
        shearRate = tilt_vel/lz;  // the shearRate = vel/L 
                
        //h_rate[i] = TWOPI/set[i].tperiod * set[i].amplitude *
        //  cos(TWOPI*delt/set[i].tperiod);
        //h_ratelo[i] = -0.5*h_rate[i];
        h_rate[2] = tilt_vel; // shearRate;
        h_ratelo[2] = -0.5*h_rate[2]; // PJA: double-check the rates here...

        // update the target values for SELM codes
        //SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last 
        //  = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;
        //SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
        //  = set[i].tilt_target;
        shearData_RM_OSC1->shearRate = shearRate;
        SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate = shearRate;        
      }

      // tilt_target can be large positive or large negative value
      // add/subtract box lengths until tilt_target is closest to current value
      int idenom = 0;
      
      if (i == 5) idenom = 0;
      else if (i == 4) idenom = 0;
      else if (i == 3) idenom = 1;
      
      double denom = set[idenom].hi_target - set[idenom].lo_target;

      double current = h[i]/h[idenom];

      while (set[i].tilt_target/denom - current > 0.0)
        set[i].tilt_target -= denom;
      while (set[i].tilt_target/denom - current < 0.0)
        set[i].tilt_target += denom;
      if (fabs(set[i].tilt_target/denom - 1.0 - current) <
          fabs(set[i].tilt_target/denom - current))
        set[i].tilt_target -= denom;      
    }

    // atz -- (added to move the bottom of the box to make shear symmetric)
    double lx = domain->boxhi[0] - domain->boxlo[0];
    set[0].lo_target = set[0].lo_start - 0.5*set[4].tilt_target;
    set[0].hi_target = set[0].lo_target + lx;
    //set[0].lo_target = set[0].lo_start + delta*(set[0].lo_stop - set[0].lo_start);
    //set[0].hi_target = set[0].lo_target + lx;
    // --
 
  }

  // if any tilt ratios exceed 0.5, set flip = 1 and compute new tilt values
  // do not flip in x or y if non-periodic (can tilt but not flip)
  //   this is b/c the box length would be changed (dramatically) by flip
  // if xz tilt exceeded, adjust C vector by one A vector
  // flip is performed on next timestep, before reneighboring in pre-exchange()
  if (triclinic && flipflag) {
    double xprd = set[0].hi_target - set[0].lo_target;
    double yprd = set[1].hi_target - set[1].lo_target;
    double xprdinv = 1.0 / xprd;
    double yprdinv = 1.0 / yprd;
    if (set[3].tilt_target*yprdinv < -0.5 ||
                                     set[3].tilt_target*yprdinv > 0.5 ||
        set[4].tilt_target*xprdinv < -0.5 ||
                                     set[4].tilt_target*xprdinv > 0.5 ||
        set[5].tilt_target*xprdinv < -0.5 ||
                                     set[5].tilt_target*xprdinv > 0.5) {
      set[3].tilt_flip = set[3].tilt_target;
      set[4].tilt_flip = set[4].tilt_target;
      set[5].tilt_flip = set[5].tilt_target;

      flipxy = flipxz = flipyz = 0;

      if (domain->yperiodic) {
        if (set[3].tilt_flip*yprdinv < -0.5) {
          set[3].tilt_flip += yprd;
          set[4].tilt_flip += set[5].tilt_flip;
          flipyz = 1;
        } else if (set[3].tilt_flip*yprdinv > 0.5) {
          set[3].tilt_flip -= yprd;
          set[4].tilt_flip -= set[5].tilt_flip;
          flipyz = -1;
        }
      }
      if (domain->xperiodic) {
        if (set[4].tilt_flip*xprdinv < -0.5) {
          set[4].tilt_flip += xprd;
          flipxz = 1;
        }
        if (set[4].tilt_flip*xprdinv > 0.5) {
          set[4].tilt_flip -= xprd;
          flipxz = -1;
        }
        if (set[5].tilt_flip*xprdinv < -0.5) {
          set[5].tilt_flip += xprd;
          flipxy = 1;
        }
        if (set[5].tilt_flip*xprdinv > 0.5) {
          set[5].tilt_flip -= xprd;
          flipxy = -1;
        }
      }

      flip = 0;
      if (flipxy || flipxz || flipyz) flip = 1;
      if (flip) driver_selm->fixSELM->next_reneighbor = update->ntimestep + 1;
    }
  }

  // convert atoms and rigid bodies to lamda coords
  if (remapflag == Domain::X_REMAP) {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & driver_selm->fixSELM->groupbit)
        domain->x2lamda(x[i],x[i]);
  }

  // reset global and local box to new size/shape
  // only if deform fix is controlling the dimension
  if (triclinic) {
    if (set[4].style) {
      domain->xz = set[4].tilt_target;
      domain->boxlo[0] = set[0].lo_target;
      domain->boxhi[0] = set[0].hi_target;
    }
    // PJA: changes the xlo and xhi, so bottom also moves
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert atoms and rigid bodies back to box coords
  if (remapflag == Domain::X_REMAP) {
    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & driver_selm->fixSELM->groupbit)
        domain->lamda2x(x[i],x[i]);
  }
  
  // redo KSpace coeffs since box has changed
  if (kspace_flag) lammps->force->kspace->setup();

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::setGlobalRefs(LAMMPS *lmps, DriverSELM *fx) {

  /* setup the LAMMPS related information */
  lammps  = lmps;
  driver_selm = fx;

  //if (force_reneighbor)
  //  irregular = new Irregular(lammps);
  // (setup in constructor)
  //else
  //  irregular = NULL;

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::init() {

  type = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::TYPE;
  strcpy(typeStr, SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::TYPE_STR);
  
  flagWriteSimulationData = 0;
  saveSkipSimulationData  = -1;

}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::parse_ParameterFile(char *baseFilename) {

  const char *error_str_code = "SELM_Integrator_SHEAR_QUASI_STEADY1_FFTW3.cpp";
  const char *error_str_func = "parse_ParameterFile()";

  /* nothing to do currently, now XML is preferred method to load parameters */
    
}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::setup() {

  /* setup the integrator data structure */
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras
    = (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType *)
      malloc(sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType));

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->deltaT
  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->deltaT;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->maxTimeStepIndex
    = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->maxTimeStepIndex;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->mu
    = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->mu;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->rho
      = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->rho;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->KB
      = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->KB;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->T
        = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->T;

  strcpy(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearModeStr,
         SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->flagShearModeStr);

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode
  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->flagShearMode;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData
  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->shearData;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->maxTimeStepIndex
  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->maxTimeStepIndex;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagStochasticDriving
  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->flagStochasticDriving;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagIncompressibleFluid
  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->flagIncompressibleFluid;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagUpdateControlPts = 1;

  flagWriteSimulationData = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->flagWriteSimulationData;
  saveSkipSimulationData  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params->saveSkipSimulationData;

  free(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params);
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagInitializedNumericalMethod = 1;
  
  // now since initialized, setup the internal init() for the class
  setup_internals();
  
}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::packageError(int code, void *extras) {
    std::exit(code);
}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::set_LAMMPS_mask(int *mask_ptr) {

  /* this determines if LAMMPS calls both the initial and final integrator commands */
  /* this command augments the mask flags to control this behavior */
  // WARNING: Might need to add the above bits to ensure integrates also initial() and final() time-step 

  //(*mask_ptr) = setmask(); // set the mask value (handled now below to make more explicit) 

  (*mask_ptr) = 0;
  if (driver_selm->fixSELM->force_reneighbor) (*mask_ptr) |= PRE_EXCHANGE;
  (*mask_ptr) |= END_OF_STEP;

  (*mask_ptr) |= INITIAL_INTEGRATE; /* tells LAMMPS to call integrate_initial via fix_SELM */
  (*mask_ptr) |= FINAL_INTEGRATE;   /* tells LAMMPS to call integrate_final via fix_SELM */

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::integrate_initialize() {

  const char *error_str_func = "integrate_initialize()";

}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::integrate_initial() {

  const char *error_str_func = "integrate_initial()";

  /* synchronize the lammps domain with the fluid domain for shear deformation */
  syncShearDomainWithLammpsDomain();

  /* perform the calculations before the integrator is called */
  IB_appl1_start_time_step_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3();

  SELM_updateFluidAndStructures_initial();

  /* now done internally during fluid-structure update */
  /* synchronize the lammps domain with the fluid domain for shear deformation */
  /* setLammpsDomainFromShearDomain(); */

}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::integrate_final() {

  const char *error_str_func = "integrate_final()";

  SELM_updateFluidAndStructures_final();

}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_updateParticlesTest_initial()
{

  const char *error_str_func = "SELM_updateParticlesTest_initial()";


  Atom   *atom    = lammps->atom;
  int    igroup   = driver_selm->fixSELM->igroup;
  int    groupbit = driver_selm->fixSELM->groupbit;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;
  double *rmass   = atom->rmass;
  double *mass    = atom->mass;
  int    *type    = atom->type;
  int    *mask    = atom->mask;
  int     nlocal  = atom->nlocal;

  double  dtfm    = 1.0/0.0;
  double  dtv     = 1.0/0.0;
  double  dtf     = 1.0/0.0;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }


}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_updateParticlesTest_final()
{

  /*
  Atom   *atom    = lammps->atom;
  int    igroup   = driver_selm->fixSELM->igroup;
  int    groupbit = driver_selm->fixSELM->groupbit;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;
  double *rmass   = atom->rmass;
  double *mass    = atom->mass;
  int    *type    = atom->type;
  int    *mask    = atom->mask;
  int     nlocal  = atom->nlocal;
  */
  
  // PJA: nothing to do currently, but retain the hook. 

}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_updateFluidAndStructures_initial()
{

  const char *error_str_func = "SELM_updateFluidAndStructures_initial()";

  int d,j,k,ell,I,k1,k2,k3;
  int num_dim = -1;  
  int *numMeshPtsPerDir = NULL;
  double meshDeltaX = 0.0;
  double deltaT = 0.0;
  
  double a_k_j[3] = {0,0,0};
   
  double normalizeFluidDensity_k;

  double mu,rho,KB,T;

  int flagShearMode = -1;

  int    shearDir = -1,shearVelDir = -1;
  double shearRate = 0,shearDist = 0,shearDist_last = 0;
  //double shearOmega,shearRateAmplitude;
  
  double L_shearDir = 0, L_shearVelDir = 0;

  int a,b,k_a,k_b,N_a,N_b;
  int vec_k[3] = {-1,-1,-1};

  double rep_d_2_w_j__d_q_a__d_q_b;

  double delta_a_velDir,delta_b_velDir;
  double delta_a_ell,delta_b_ell;
  double delta_ell_shearDir;
  
  double dot_gamma_t = 0.0;

  int timeIndex = -1;

  double curTime = 0.0;

  int flagUseFiniteDifference;

  double *meshCenterX0 = NULL;
  double shearDistMod = -1;

  int    flagImposeUnitCellRep = 0;
  int    flipToLessShearedUnitCell = 0;
  double XX1[3],periodL[3],vec_L[3];
  
  int flag_debug = 0;
  double *controlPt_X = NULL;

  ShearData_RM_SHEAR1_Type *shearData_RM_SHEAR1 = NULL;
  ShearData_RM_OSC1_Type   *shearData_RM_OSC1 = NULL;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3                                                       *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3 = NULL;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras = NULL;
  
  class Domain *domain = lammps->domain;
  class Atom   *atom   = lammps->atom;  

  // get information from lammps
  timeIndex = lammps->update->ntimestep; /* ### */
  
  if (flag_debug > 0) {
    std::cout << "lammps->update->ntimestep " << lammps->update->ntimestep << std::endl;
  }

#ifdef USE_PACKAGE_FFTW3

  /* == Dereference the data (WARNING: assumes only one Eulerian DOF for now) */
  if (driver_selm->SELM_Eulerian_List[0]->type == SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

    SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
      = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0];

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
      = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

    num_dim          = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;
    meshDeltaX       = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
    numMeshPtsPerDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
    meshCenterX0     = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0;

  } else {
    stringstream message;
    message << "Expecting mesh type: " << SELM_Eulerian_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3 << endl;
    message << "Instead mesh type was: " << driver_selm->SELM_Eulerian_List[0]->typeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  if (this->type == SELM_Integrator_Types::TYPE_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3) {

    /* time stepping parameters */
    deltaT = lammps->update->dt;
    /* LAMMPS now overrides any deltaT given in an XML file */ 
    SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->deltaT = deltaT; 
    
    curTime       = timeIndex*deltaT;
    flagShearMode = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode;

    /* fluid parameters */
    mu  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->mu;
    rho = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->rho;
    KB  = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->KB;
    T   = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->T;

    switch (flagShearMode) {

    case SHEAR_MODE_TYPE_RM_SHEAR1:

      shearData_RM_SHEAR1 = (ShearData_RM_SHEAR1_Type *)
          SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

      /* retrieve shear values required to generate the flow */
      shearDir       = shearData_RM_SHEAR1->shearDir;
      shearVelDir    = shearData_RM_SHEAR1->shearVelDir;
      shearRate      = shearData_RM_SHEAR1->shearRate;
      shearDist      = shearData_RM_SHEAR1->shearDist;
      shearDist_last = shearData_RM_SHEAR1->shearDist_last;

      /* setup the mesh in accordance of the current shear conditions */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir
        = shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir
        = shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate
        = shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
        = shearDist;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
        = shearDist_last;

      break;

    case SHEAR_MODE_TYPE_RM_OSC1:

      shearData_RM_OSC1
        = (ShearData_RM_OSC1_Type *)
        SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

      /* retrieve shear values required to generate the flow */
      shearDir           = shearData_RM_OSC1->shearDir;
      shearVelDir        = shearData_RM_OSC1->shearVelDir;      
      //shearOmega         = shearData_RM_OSC1->shearOmega;
      //shearRateAmplitude = shearData_RM_OSC1->shearRateAmplitude;
      shearRate          = shearData_RM_OSC1->shearRate;
      shearDist          = shearData_RM_OSC1->shearDist;
      shearDist_last     = shearData_RM_OSC1->shearDist_last;

      /* setup the mesh in accordance of the current shear conditions */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir
        = shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir
        = shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate
        = shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
        = shearDist;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
        = shearDist_last;

      break;

    } /* end of switch */

  } else {
    stringstream message;
    message << "Expecting time stepper: " << SELM_Integrator_Types::TYPE_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3 << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }
  
  /* == Update the Control Points from the LAMMPS data */
  /* synchronize the lammps data and control points data to one another */
  for (I = 0; I < driver_selm->SELM_Lagrangian_List_N; I++) {
    driver_selm->SELM_Lagrangian_List[I]->setControlPtsDataFromLammpsData();
  }
  
  /* == Add the Pseudoforce to the force density of the fluid */
  /* compute the pseudo-force term arising from the shear boundary condition */
  addPseudoForceTerm(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                     SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                     SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                     mu,shearDir,shearVelDir,shearRate,
                     SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m);

  /* == Compute the DFT of the force density of the fluid */
  normalizeFluidDensity_k = 1;
  for (d = 0; d < num_dim; d++) {
    fftw_execute(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_DFT_plan[d]);
    normalizeFluidDensity_k = normalizeFluidDensity_k * numMeshPtsPerDir[d];
  }

  /* == Compute the steady-state fluid velocity */
  if (num_dim == 2) {

    stringstream message;
    message << "2D Steady-State Sheared Stokes flow not implemented currently." << endl;
    message << "Fluid domain is assumed 3D periodic with Lees-Edwards boundary conditions." << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

  } /* end num_dim == 2 */

  /* == Compute the steady-state fluid velocity */
  if (num_dim == 3) {
  
    // initialize 
    for (j = 0; j < num_dim; j++) {
      a_k_j[j] = 0.0;
      
      if (flag_debug > 0) {
        std::cout << "a_k_j[" << j << "] = " << a_k_j[j] << std::endl;
      }      

    }
    
    /* compute the steady-state Stokes flow: w = L^{-1}g */
    L_shearDir = numMeshPtsPerDir[shearDir]*meshDeltaX;
    dot_gamma_t = shearDist/L_shearDir;
    vec_L[0] = numMeshPtsPerDir[0]*meshDeltaX; vec_L[1] = numMeshPtsPerDir[1]*meshDeltaX; vec_L[2] = numMeshPtsPerDir[2]*meshDeltaX;
    
    for (k3 = 0; k3 < numMeshPtsPerDir[2]; k3++) {
      for (k2 = 0; k2 < numMeshPtsPerDir[1]; k2++) {
        for (k1 = 0; k1 < numMeshPtsPerDir[0]; k1++) {

          I = (k3*numMeshPtsPerDir[1]*numMeshPtsPerDir[0]) + (k2*numMeshPtsPerDir[0]) + k1;

          vec_k[0] = k1; vec_k[1] = k2; vec_k[2] = k3;

          flagUseFiniteDifference = 1;
          if (flagUseFiniteDifference) {

            if ((k1 != 0) || (k2 != 0) || (k3 != 0)) {

              /*
                 Stokes flow subject to Lees-Edwards boundary
                 conditions (periodic images are sheared
                 relative to base image).
                 
                 Finite difference approximation of the Laplacian
                 operator is used under the change of variable
                 which moves in the coordinate system of a
                 sheared box deforming with the flow.  This
                 coordinate system allows for the flow field
                 to be represented on a periodic domain with
                 a jump boundary condition, see paper.  The Laplacian
                 becomes under this change of variable for
                 shear along the z-axis (j = 2) as follows:

                 TeX:
                 $\mathbf{x} = \phi(\mathbf{q}) = (q_1 + \dot{\gamma}tq_3,q_2,q_3)$,
                 $\mathbf{q} = \phi^{-1}(\mathbf{x}) = (x_1 - \dot{\gamma}tx_3,x_2,x_3)$
                 and Laplacian can be expressed as
                 $$\mathcal{L}_{\mathbf{x}} \mathbf{u}^{(j)} 
                 = \sum_{\ell=1}^3 
                 \frac{\partial^2 \mathbf{u}^{(j)} }
                 {\partial \mathbf{x}^{(\ell),2}}
                 = \mathcal{L}_{\mathbf{q}}(t) \mathbf{w}^{(d)}
                 = \sum_{\ell=1}^3 \sum_{a,b = 1}^3 
                 (\delta_{\ell,a} - \dot{\gamma}t \delta_{\ell,3} \delta_{a,1})
                 \frac{\partial^2 \mathbf{w}^{(j)}}
                 {\partial \mathbf{q}^{(a)} \partial \mathbf{q}^{(b)}}
                 (\delta_{\ell,b} - \dot{\gamma}t \delta_{\ell,3} \delta_{b,1})$$

                 We shall use a finite difference approximation of the
                 partial derivatives.  The boundary conditions are
                 incorporated by adding the appropriate jump condition
                 to the values used in the finite difference stencil.
                 Mathematically, this can be taken into account by
                 introducing a source term in the steady-state Stokes
                 equations (move known values to RHS in the linear
                 equations) or even by a linear flow in quasi-steady 
                 regime.  The matrix is then cylcic and inverted
                 using the FFT.

               */
              for (j = 0; j < num_dim; j++) {

                /* compute the fourier representation
                   of the finite difference operator */
                a_k_j[j] = 0.0;
                for (ell = 0; ell < num_dim; ell++) {
                  for (a = 0; a < num_dim; a++) {
                    for (b = 0; b < num_dim; b++) {

                      k_a = vec_k[a]; k_b = vec_k[b];

                      N_a = numMeshPtsPerDir[a]; N_b = numMeshPtsPerDir[b];
                      
                      if (a == ell) {
                        delta_a_ell = 1;
                      } else {
                        delta_a_ell = 0;
                      }
                      
                     if (b == ell) {
                        delta_b_ell = 1;
                      } else {
                        delta_b_ell = 0;
                      }

                      if (a == shearVelDir) {
                        delta_a_velDir = 1;
                      } else {
                        delta_a_velDir = 0;
                      }

                      if (b == shearVelDir) {
                        delta_b_velDir = 1;
                      } else {
                        delta_b_velDir = 0;
                      }

                      if (ell == shearDir) {
                        delta_ell_shearDir = 1;
                      } else {
                        delta_ell_shearDir = 0;
                      }

                      if (a != b) { 
                        /* uses composition of central
                           differences in q_a and q_b directions */                        
                        rep_d_2_w_j__d_q_a__d_q_b = -(sin(2.0*UNIT_pi*k_a/N_a)*sin(2.0*UNIT_pi*k_b/N_b))
                                                  /(meshDeltaX*meshDeltaX);

                      } else { /* uses second order central differences in q_a direction */
                        rep_d_2_w_j__d_q_a__d_q_b = -2.0*(1.0 - cos(2.0*UNIT_pi*k_a/N_a))/(meshDeltaX*meshDeltaX);
                      }

                      a_k_j[j] += mu*(delta_a_ell - dot_gamma_t*delta_ell_shearDir*delta_a_velDir)
                                    *(delta_b_ell - dot_gamma_t*delta_ell_shearDir*delta_b_velDir)
                                    *rep_d_2_w_j__d_q_a__d_q_b;

                    } /* end of b loop */
                  } /* end of a loop */

                } /* end of ell loop */

              } /* end of j loop */

              for (j = 0; j < num_dim; j++) {

                if (a_k_j[j] > 1e-9) {

                  stringstream message;

                  message << "a_k_j[" << j << "] = " << a_k_j[j] << " > 0 failed to hold." << endl;
                  message << "Rough check on the definiteness of the difference operator." << endl;
                  message << "Need to be careful of round-off errors if persists." << endl;

                  SELM_Package::packageError(error_str_code, error_str_func, message);

                }

                SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[j][I][0]
                  = -SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_k[j][I][0]/(a_k_j[j]*normalizeFluidDensity_k);
                SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[j][I][1]
                  = -SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_k[j][I][1]/(a_k_j[j]*normalizeFluidDensity_k);
              } /* end of j loop */

            } else { /* set the zero mode to zero */

              for (d = 0; d < num_dim; d++) {
                SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d][I][0] = 0.0;
                SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d][I][1] = 0.0;
              }

            }

          } /* end of use finite difference approximation */

        } /* end k1 loop */
      } /* end k2 loop */
    } /* end k3 loop */

  } /* end num_dim == 3 */

  /* == Compute the effective stochastic forcing (time integrated average) for the fluid */
  if (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagStochasticDriving) {

    /* -- Compute the the time averaged stochastic fluctuations of the fluid */
    computeTimeAvgStochFluct(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                             SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                             SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                             mu,KB,T,
                             SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->deltaT,
                             shearRate,shearDir,shearVelDir,shearDist,
                             driver_selm->random,
                             SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_k);

    /* -- Add the time average stochastic fluctuations to the fluid velocity field */
    for (k3 = 0; k3 < numMeshPtsPerDir[2]; k3++) {
      for (k2 = 0; k2 < numMeshPtsPerDir[1]; k2++) {
        for (k1 = 0; k1 < numMeshPtsPerDir[0]; k1++) {

          I = (k3*numMeshPtsPerDir[1]*numMeshPtsPerDir[0]) + (k2*numMeshPtsPerDir[0]) + k1;

          for (d = 0; d < num_dim; d++) {
            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d][I][0]
                                                                      += SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_k[d][I][0];
            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d][I][1]
                                                                      += SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_k[d][I][1];
          } /* end of d loop */

        } /* end k1 loop */
      } /* end k2 loop */
    } /* end k3 loop */

  } /* end of check for flagStochasticDriving */

  /* == For incompressible fluids apply projection operator */
  if (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagIncompressibleFluid) {

     projectField(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                  shearDir,shearVelDir,shearDist,
                  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k);

    if (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure) {

      if (timeIndex <= 5) {
        printf("WARNING: %s : %s \n", error_str_code, error_str_func);
        printf("Note that the explicit pressure field p(x) is mainly for visualization and not used in other calculations. \n");
        // The explicit pressure field below for visualization may not yet be debugged fully (double-check)
        // The simulation calculation uses instead projection operator, so does not get used there.
        // This is maninly for visualization purposes.
      }

      /* Compute pressure associated with the force density.  This only includes
         pressure for the 'deterministic part' of the force field. */
      computePressure(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                      shearDir,shearVelDir,shearDist,
                      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_k,
                      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_k);

    } /*end of flagComputePressureField */

  } /* end of flagIncompressibleFluid */

  /* == Compute the IDFT of the velocity field of the fluid */
  for (d = 0; d < num_dim; d++) {
    fftw_execute(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_IDFT_plan[d]);
  }

  /* == Compute the IDFT of the pressure field of the fluid */
  if ((SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagIncompressibleFluid)
      && (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure)) {
    fftw_execute(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_IDFT_plan);
  }

  /* == Update the control point locations (if flag set) */
  if (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagUpdateControlPts) {

      /* The average velocity is computed from the sheared four point region about the
         structure mapped to the sheared coordinate plane, which requires positions of
         the structures assumed to be described in physical space to be adjusted accordingly
         to get the corresponding structure positions in the sheared coordinate frame.
         The velocity field on the periodic mesh also needs to be handled carefully
         with periodic images being adjusted with the appropriate Lees-Edwards jumps
         when crossing a periodic boundary.              
      */

      /************************************************************************************************/
      /* Loop over the coupling operators to determine the Lagrangian velocities via Lambda operator.
      /************************************************************************************************/
      /* loop over the coupling operators DOF */
      for (int couplingOpI = 0;couplingOpI < driver_selm->SELM_CouplingOperator_List_N; couplingOpI++) {

        SELM_CouplingOperator *couplingOp = driver_selm->SELM_CouplingOperator_List[couplingOpI];

        /* maybe make below more generic */
        if (couplingOp->type == SELM_CouplingOperator_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1) {

          SELM_Lagrangian *lagrangian = NULL;
          SELM_Eulerian   *eulerian   = NULL;

          SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 *op
            = (SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 *)couplingOp;

          for (int k = 0; k < op->numCoupleList; k++) {
            lagrangian = op->lagrangianList[k];
            eulerian   = op->eulerianList[k];

            computeControlPtsVel_SHEAR_FFTW3(lagrangian,eulerian,couplingOp);
          }

          /* update the lagrangian DOF */
          if (lagrangian->type == SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_ANGLE_STYLE) {

            SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE
              = (SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *)lagrangian;

            for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numControlPts; k++) {
              for (d = 0; d < num_dim; d++) {
                SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX[(k * num_dim) + d]
                  += SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->pt_Vel[(k* num_dim) + d]*deltaT;
              }
            }

            /* Impose that the control point locations are always described within the
               rectangular unit cell with Lees-Edwards boundary conditions.  Note the
               Stokes flow is considered in the sheared unit cell.  We note that this is
               also important since stress calculations must be handled carefully when
               bonds cross a periodic boundary.
             */
            flagImposeUnitCellRep = 0;  
            if (flagImposeUnitCellRep == 1) {
              /* For situations where we need to impose 
                 the unit cell representation explicitly.
                 Currently our modified fix_deform codes 
                 handle this between steps. */

              for (d = 0; d < num_dim; d++) {
                periodL[d] = numMeshPtsPerDir[d] * meshDeltaX;
              }

              k = floor(fabs(shearDist / periodL[shearVelDir]));
              if (shearDist > 0) {
                shearDistMod = shearDist - k*periodL[shearVelDir];
              } else {
                shearDistMod = shearDist + k*periodL[shearVelDir];
              }

              for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numControlPts; k++) {
                controlPt_X = (SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX + (k * num_dim));

                bri1_unitCellRectImageShearPeriodic((double *) periodL,
                                                    (double *) meshCenterX0, shearDir,
                                                    shearVelDir, shearDistMod, num_dim,
                                                    controlPt_X, (double *) XX1);

                for (d = 0; d < num_dim; d++) {
                  SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX[num_dim * k + d] = XX1[d];
                }

              } /* end of k loop */

            } /* end flagImposeUnitCellRep */

            /* Adjustments may be made later to the structures to avoid excessive drifting in
               stress calculations for linear shearing.
             */
             
          } else if (lagrangian->type == SELM_Lagrangian_Types::TYPE_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE) {

            SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE
              = (SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *)lagrangian;

            for (k = 0; k < SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numControlPts; k++) {
              for (d = 0; d < num_dim; d++) {
                SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX[(k * num_dim) + d]
                  += SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->pt_Vel[(k* num_dim) + d] * deltaT;
              }
            }

            /* Impose that the control point locations are always described within the
               rectangular unit cell with Lees-Edwards boundary conditions.  Note the
               Stokes flow is considered in the sheared unit cell.  We note that this is
               also important since stress calculations must be handled carefully when
               bonds cross a periodic boundary.
             */
            flagImposeUnitCellRep = 0;  
            if (flagImposeUnitCellRep == 1) {
              /* For situations where we need to impose 
                 the unit cell representation explicitly.
                 Currently our modified fix_deform codes 
                 handle this between steps. */

              for (d = 0; d < num_dim; d++) {
                periodL[d] = numMeshPtsPerDir[d] * meshDeltaX;
              }

              k = floor(fabs(shearDist / periodL[shearVelDir]));
              if (shearDist > 0) {
                shearDistMod = shearDist - k*periodL[shearVelDir];
              } else {
                shearDistMod = shearDist + k*periodL[shearVelDir];
              }

              for (k = 0; k < SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numControlPts; k++) {
                controlPt_X = (SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX + (k * num_dim));

                bri1_unitCellRectImageShearPeriodic((double *) periodL,
                                                    (double *) meshCenterX0, shearDir,
                                                    shearVelDir, shearDistMod, num_dim,
                                                    controlPt_X, (double *) XX1);

                for (d = 0; d < num_dim; d++) {
                  SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX[num_dim * k + d] = XX1[d];
                }

              } /* end of k loop */

            } /* end flagImposeUnitCellRep */

            /* Adjustments may be made later to the structures to avoid excessive drifting in
             * stress calculations for linear shearing.
             */

          } else {
            stringstream message;
            message << "Lagrangian Type is not recognized" << endl;
            SELM_Package::packageError(error_str_code, error_str_func, message);
          }

        } else { /* end of TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 */
          stringstream message;
          message << "Integrator does not know how to handle coupling operator type. " << endl;
          message << "  couplingOp->typeStr = " << couplingOp->typeStr << endl;
          SELM_Package::packageError(error_str_code, error_str_func, message);
        }

      } /* end of couplingOpI loop */
      
      /* synchronize the lammps data and control points data to one another */
      /* @optimize: only update the points that changed from above, for now we update all points */
      for (int I = 0; I < driver_selm->SELM_Lagrangian_List_N; I++) {
        driver_selm->SELM_Lagrangian_List[I]->setLammpsDataFromControlPtsData();
      }

   } /* end flagUpdateControlPts */

  /**********************************************************************************/

  /* == Update the shearing of the coordinate system for use in next time step */
  shearDist_last = shearDist;
  L_shearDir = numMeshPtsPerDir[shearDir]*meshDeltaX;
  shearDist = lammps->domain->xz; // we sync with LAMMPS using the xz shear */
  
  /* By symmetry of the system we can always choose the coordinate system for
     the propogator of the system which is sheared at most by one mesh-width spacing.
     The basic idea for the quasi-steady-state Stokes flow is to treat the
     fluid mesh always for a box which is sheared at most one period in the
     shear velocity direction.  In other words, for any arbitrary shear displacement
     the unit cell can always be mapped to one sheared by at most one period in
     the shear velocity direction by using the periodicity.

     WARNING: It may be important not to modify the shearDist variable as this could
     impact how the structure deformations are interpreted over the mesh.
  */
  flagShearMode = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode;
  switch (flagShearMode) {

  case SHEAR_MODE_TYPE_RM_SHEAR1:

    shearData_RM_SHEAR1
      = (ShearData_RM_SHEAR1_Type *)
    SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

    /* we now use LAMMPS methods and the modified fix_deform codes to handle this between steps. */
    flipToLessShearedUnitCell = 0;
    
    /* record the change in shearDist */
    shearData_RM_SHEAR1->shearDist_last = shearDist_last;
    shearData_RM_SHEAR1->shearDist      = shearDist;

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
    = shearData_RM_SHEAR1->shearDist_last;

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
    = shearData_RM_SHEAR1->shearDist;

    break;

  case SHEAR_MODE_TYPE_RM_OSC1:

    shearData_RM_OSC1
      = (ShearData_RM_OSC1_Type *)
      SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;
      
    /* we now use LAMMPS methods and the modified fix_deform codes to handle this between steps. */
    flipToLessShearedUnitCell = 0;

    /* record the change in shearDist */
    shearData_RM_OSC1->shearDist_last = shearDist_last;
    shearData_RM_OSC1->shearDist      = shearDist;

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last
      = shearData_RM_OSC1->shearDist_last;

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist
      = shearData_RM_OSC1->shearDist;

    break;

  } /* end of switch */

#else

  stringstream message;
  message << "Compiler flag USE_PACKAGE_FFTW3 was not set." << endl;
  message << "None of the FFTW3 routines were compiled in the codes." << endl;
  SELM_Package::packageError(error_str_code, error_str_func, message);

#endif

}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_updateFluidAndStructures_final()
{

  const char *error_str_func = "SELM_updateFluidAndStructures_final()";
  
  // nothing to do currently, but retain the hook

}

/*====================================== Private ===================================== */
void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::addPseudoForceTerm(int num_dim,
                                                                          double meshDeltaX,
                                                                          int *numMeshPtsPerDir,
                                                                          double mu,
                                                                          int shearDir,
                                                                          int shearVelDir,
                                                                          double shearRate,
                                                                          fftw_complex **fluidForceDensity_m) {

  const char *error_str_func = "addPseudoForceTerm()";

  int a,b,j1,j2,m1,m2,m3,I_bot,I_top;
  int vec_m_top[3],vec_m_bot[3];
  double L_shearDir,contr_top_d_2_e_d_x_a_d_x_b,contr_bot_d_2_e_d_x_a_d_x_b;

  if (num_dim == 3) {

    /* compute the contribution from the finite difference stencil along the top and
       bottom shear plane of the domain
     */
    for (j1 = 0; j1 < numMeshPtsPerDir[(shearDir + 1) % num_dim]; j1++) {
      for (j2 = 0; j2 < numMeshPtsPerDir[(shearDir + 2) % num_dim]; j2++) {

        vec_m_top[(shearDir + 1) % num_dim] = j1;
        vec_m_top[(shearDir + 2) % num_dim] = j2;
        vec_m_top[(shearDir + 0) % num_dim] = (numMeshPtsPerDir[shearDir] - 1);

        vec_m_bot[(shearDir + 1) % num_dim] = j1;
        vec_m_bot[(shearDir + 2) % num_dim] = j2;
        vec_m_bot[(shearDir + 0) % num_dim] = 0;

        m1 = vec_m_top[0];
        m2 = vec_m_top[1];
        m3 = vec_m_top[2];
        I_top = (m3 * numMeshPtsPerDir[1] * numMeshPtsPerDir[0]) + (m2
            * numMeshPtsPerDir[0]) + m1;

        m1 = vec_m_bot[0];
        m2 = vec_m_bot[1];
        m3 = vec_m_bot[2];
        I_bot = (m3 * numMeshPtsPerDir[1] * numMeshPtsPerDir[0]) + (m2
            * numMeshPtsPerDir[0]) + m1;

        L_shearDir = numMeshPtsPerDir[shearDir] * meshDeltaX;

        contr_top_d_2_e_d_x_a_d_x_b = 0.0;
        contr_bot_d_2_e_d_x_a_d_x_b = 0.0;
        for (a = 0; a < num_dim; a++) {
          for (b = 0; b < num_dim; b++) {

            if (a != b) {
              contr_top_d_2_e_d_x_a_d_x_b += 0.0;
              contr_bot_d_2_e_d_x_a_d_x_b += 0.0;
            } else { /* a == b */

              /* only contributes in the shear direction */
              if (a == shearDir) {
                contr_top_d_2_e_d_x_a_d_x_b += mu * shearRate * L_shearDir
                    / (meshDeltaX * meshDeltaX);
                contr_bot_d_2_e_d_x_a_d_x_b += -mu * shearRate * L_shearDir
                    / (meshDeltaX * meshDeltaX);
              }

            } /* end else a == b */

          } /* end of b loop */
        } /* end of a loop */

        fluidForceDensity_m[shearVelDir][I_top][0] += contr_top_d_2_e_d_x_a_d_x_b;
        fluidForceDensity_m[shearVelDir][I_top][1] += 0.0;
        fluidForceDensity_m[shearVelDir][I_bot][0] += contr_bot_d_2_e_d_x_a_d_x_b;
        fluidForceDensity_m[shearVelDir][I_bot][1] += 0.0;

      } /* end of j2 loop */
    } /* end of j1 loop */

  } else {/* end num_dim == 3 */

    stringstream message;

    message << "num_dim = " << num_dim << endl;
    message << "Pseudo force computation at shear boundary not implemented in this case." << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

  }

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::computeTimeAvgStochFluct(int num_dim,
                                                                         double meshDeltaX,
                                                                         int *numMeshPtsPerDir,
                                                                         double mu,
                                                                         double KB,
                                                                         double T,
                                                                         double deltaT,
                                                                         double shearRate,
                                                                         int shearDir,
                                                                         int shearVelDir,
                                                                         double shearDist,
                                                                         class RanMars *random,
                                                                         fftw_complex **fluidStochForceDensity_k) {

  const char *error_str_func = "computeTimeAvgStochFluct()";

  int d,j,ell,I,I1,I2,k1,k2,k3,kk1,kk2,kk3;

  double G_real_kk,G_imag_kk,sqrt_G_real_kk,sqrt_G_imag_kk;

  double eta_real,eta_imag;

  double meshVolume,cellVolume;

  int a,b,k_a,k_b,N_a,N_b;

  int vec_k[3];

  double rep_d_2_w_j__d_q_a__d_q_b;
  double a_k_j[3];

  double delta_a_velDir,delta_b_velDir,delta_a_ell,delta_b_ell,delta_ell_shearDir;
  
  double L_shearDir,dot_gamma_t;
  double vec_L[3];

  int prodNumMeshPts;

  /* == Compute the DFT of the force density of the fluid */
  /* normalizeFluidDensity_k = 1;
     for (d = 0; d < num_dim; d++) {
     normalizeFluidDensity_k
     = normalizeFluidDensity_k*numMeshPtsPerDir[d];
     }
   */
  if (num_dim == 3) {

    prodNumMeshPts = 1;
    for (d = 0; d < num_dim; d++) {
      prodNumMeshPts = prodNumMeshPts * numMeshPtsPerDir[d];
    }

    meshVolume = 1.0;
    for (d = 0; d < num_dim; d++) {
      meshVolume = meshVolume * numMeshPtsPerDir[d] * meshDeltaX;
    }

    cellVolume = 1.0;
    for (d = 0; d < num_dim; d++) {
      cellVolume = cellVolume * meshDeltaX;
    }

    L_shearDir = numMeshPtsPerDir[shearDir] * meshDeltaX;
    dot_gamma_t = shearDist/L_shearDir;
    vec_L[0] = numMeshPtsPerDir[0]*meshDeltaX; vec_L[1] = numMeshPtsPerDir[1]*meshDeltaX; vec_L[2] = numMeshPtsPerDir[2]*meshDeltaX;

    for (k3 = 0; k3 < numMeshPtsPerDir[2]; k3++) {
      for (k2 = 0; k2 < numMeshPtsPerDir[1]; k2++) {
        for (k1 = 0; k1 < numMeshPtsPerDir[0]; k1++) {

          I = (k3 * numMeshPtsPerDir[1] * numMeshPtsPerDir[0]) + (k2
              * numMeshPtsPerDir[0]) + k1;

          kk1 = (numMeshPtsPerDir[0] - k1) % numMeshPtsPerDir[0];
          kk2 = (numMeshPtsPerDir[1] - k2) % numMeshPtsPerDir[1];
          kk3 = (numMeshPtsPerDir[2] - k3) % numMeshPtsPerDir[2];

          I2 = (kk3 * numMeshPtsPerDir[1] * numMeshPtsPerDir[0]) + (kk2
              * numMeshPtsPerDir[0]) + kk1;

          vec_k[0] = k1; vec_k[1] = k2; vec_k[2] = k3;

          /* Compute the effective fluid fluctuations to construct
             appropriate thermal forcing of the particle equations
             as prescribed by the principle of detailed balance.

             Hydrodynamic coupling:

             H = R_a \wp L^{-1} R_a^T/ meshDeltaX^d ,

             where R_a is delta function averaging procedure.
             \wp is the projection operator.

             Detailed balance requires:

             QQ' = 2*KB*T*H.

             This requires for the fluid generate random variable
             with covariance:

             Q = sqrt(2*KB*T/meshDeltaX^d)*R_a \wp L^{-1/2}.

             where we used that \wp L^{-1} = \wp L^{-1} \wp^T.

             WARNING: Need to set the conjugate modes to be equal
             and self-conjugate to be real.
             
             See the paper for details.

           */
          if ((k1 != 0) || (k2 != 0) || (k3 != 0)) {
            /*
                 Stokes flow subject to Lees-Edwards boundary
                 conditions (periodic images are sheared
                 relative to base image).
                 Finite difference approximation of the Laplacian
                 operator is used under the change of variable
                 which moves in the coordinate system of a
                 sheared box deforming with the flow.  This
                 coordinate system allows for the flow field
                 to be represented on a periodic domain with
                 a jump boundary condition, see paper for details.
                 The Laplacian becomes under this change of variable for
                 shear along the z-axis (j = 2) as follows:

                 TeX:
                 $\mathbf{x} = \phi(\mathbf{q}) = (q_1 + \dot{\gamma}tq_3,q_2,q_3)$,
                 $\mathbf{q} = \phi^{-1}(\mathbf{x}) = (x_1 - \dot{\gamma}tx_3,x_2,x_3)$
                 and Laplacian can be expressed as
                 $$\mathcal{L}_{\mathbf{x}} \mathbf{u}^{(j)} 
                 = \sum_{\ell=1}^3 
                 \frac{\partial^2 \mathbf{u}^{(j)} }
                 {\partial \mathbf{x}^{(\ell),2}}
                 = \mathcal{L}_{\mathbf{q}}(t) \mathbf{w}^{(d)}
                 = \sum_{\ell=1}^3 \sum_{a,b = 1}^3 
                 (\delta_{\ell,a} - \dot{\gamma}t \delta_{\ell,3} \delta_{a,1})
                 \frac{\partial^2 \mathbf{w}^{(j)}}
                 {\partial \mathbf{q}^{(a)} \partial \mathbf{q}^{(b)}}
                 (\delta_{\ell,b} - \dot{\gamma}t \delta_{\ell,3} \delta_{b,1})$$

                 We shall use a finite difference approximation of the
                 partial derivatives.  The boundary conditions are
                 incorporated by adding the appropriate jump condition
                 to the values used in the finite difference stencil.
                 Mathematically, this can be taken into account by
                 introducing a source term in the steady-state Stokes
                 equations (move known values to RHS in the linear
                 equations) or even by a linear flow in quasi-steady 
                 regime.  The matrix is then cylcic and inverted
                 using the FFT.
            */
            L_shearDir = numMeshPtsPerDir[shearDir] * meshDeltaX;
            for (j = 0; j < num_dim; j++) {

              /* compute the fourier representation
              of the finite difference operator */
              a_k_j[j] = 0.0;
              for (ell = 0; ell < num_dim; ell++) {
              
                for (a = 0; a < num_dim; a++) {
                  for (b = 0; b < num_dim; b++) {

                    k_a = vec_k[a]; k_b = vec_k[b];

                    N_a = numMeshPtsPerDir[a]; N_b = numMeshPtsPerDir[b];
                    
                    if (a == ell) {
                      delta_a_ell = 1;
                    } else {
                      delta_a_ell = 0;
                    }
                    
                    if (b == ell) {
                      delta_b_ell = 1;
                    } else {
                      delta_b_ell = 0;
                    }

                    if (a == shearVelDir) {
                      delta_a_velDir = 1;
                    } else {
                      delta_a_velDir = 0;
                    }

                    if (b == shearVelDir) {
                      delta_b_velDir = 1;
                    } else {
                      delta_b_velDir = 0;
                    }

                    if (ell == shearDir) {
                      delta_ell_shearDir = 1;
                    } else {
                      delta_ell_shearDir = 0;
                    }

                    if (a != b) { /* uses composition of central
                                     differences in x_a and x_b directions */
                      rep_d_2_w_j__d_q_a__d_q_b = -(sin(2*UNIT_pi*k_a/N_a)*sin(2*UNIT_pi*k_b/N_b))
                                                  /(meshDeltaX*meshDeltaX);

                    } else { /* uses second order central differences in x_a direction */
                      rep_d_2_w_j__d_q_a__d_q_b = -2.0*(1 - cos(2*UNIT_pi*k_a/N_a))/(meshDeltaX*meshDeltaX);
                    }

                    a_k_j[j] += mu*(delta_a_ell - dot_gamma_t*delta_ell_shearDir*delta_a_velDir)
                                  *(delta_b_ell - dot_gamma_t*delta_ell_shearDir*delta_b_velDir)
                                  *rep_d_2_w_j__d_q_a__d_q_b;

                  } /* end of b loop */
                } /* end of a loop */
                
              } /* end of ell loop */

            } /* end of j loop */

            for (j = 0; j < num_dim; j++) {

              if (a_k_j[j] > 1e-9) {
                stringstream message;

                message << "a_k_j[" << j << "] = " << a_k_j[j] << " > 0 failed to hold." << endl;
                message << "Rough check on the definiteness of the difference operator." << endl;
                message << "need to be careful of round-off errors here." << endl;
                SELM_Package::packageError(error_str_code, error_str_func, message);

              }

            } /* end of the j loop */

            /*
               a         = 2.0*mu/(meshDeltaX*meshDeltaX);
               c_k1      = 1 - cos(2*UNIT_pi*k1/numMeshPtsPerDir[0]);
               c_k2      = 1 - cos(2*UNIT_pi*k2/numMeshPtsPerDir[1]);
               c_k3      = 1 - cos(2*UNIT_pi*k3/numMeshPtsPerDir[2]);

               a_k       = a*(c_k1 + c_k2 + c_k3);
             */

            /* (1.0/a_k) is -L^{-1} *//* WARNING: check this carefully,
               since it has wavenumber dependence */

            for (j = 0; j < num_dim; j++) {

              /* Note that since the discrete Fourier transform we use is not unitary
               an extra factor related to the number of lattice sites in each direction
               arises to obtain the correct covariance structure in
               physical space.  This can be readily deduced by consider the
               covariance structure of a random field having Kronecker delta correlation
               in space represented by applying our discrete Fourier transform.

               The self-conjugate modes are treated differently, since they 
               give different results than the modes with distinct conjugates.

              */
              if (I == I2) { /* indicates mode is self-conjugate */

                /* variance for mode (components) */

                /* sigma_real^2 = deltaT*(2*KB*T/dX^d)*(L^{-1})/(N^d)
                   sigma_imag^2 = 0
                */
                G_real_kk = deltaT * (2.0 * KB * T / cellVolume) * (1.0
                    / (-a_k_j[j] * prodNumMeshPts));
                G_imag_kk = 0.0;

              } else { /* mode is not self-conjugate */

                /* variance for mode (components) */

                /* sigma_real^2 = deltaT*(2*KB*T/dX^d)*(L^{-1})/(2N^d)
                   sigma_imag^2 = deltaT*(2*KB*T/dX^d)*(L^{-1})/(2N^d)
                */
                G_real_kk = deltaT * (2.0 * KB * T / cellVolume) * (1.0
                    / (-a_k_j[j] * 2.0 * prodNumMeshPts));
                G_imag_kk = deltaT * (2.0 * KB * T / cellVolume) * (1.0
                    / (-a_k_j[j] * 2.0 * prodNumMeshPts));

              }

              sqrt_G_real_kk = sqrt(G_real_kk);
              sqrt_G_imag_kk = sqrt(G_imag_kk);

              eta_real = (double) driver_selm->random->gaussian();
              eta_imag = (double) driver_selm->random->gaussian();

              fluidStochForceDensity_k[j][I][0] = sqrt_G_real_kk * eta_real
                  / deltaT;

              fluidStochForceDensity_k[j][I][1] = sqrt_G_imag_kk * eta_imag
                  / deltaT;

            } /* end of j loop */

          } else { /* zero mode set to zero */

            for (d = 0; d < num_dim; d++) {
              fluidStochForceDensity_k[d][I][0] = 0.0;

              fluidStochForceDensity_k[d][I][1] = 0.0;
            }

          } /* end check for zero mode */

        } /* end k1 loop */
      } /* end k2 loop */
    } /* end k3 loop */

    /* -- Set the conjugate modes to be equal.
       Sets each mode equal to its conjugate mode
       this enforces the constraint that the velocity
       field is real valued.

       (@optimize) Algorithm is redundant below, sets each
       mode twice.

       (@optimize) The part should be unnecessary
       in the case thermal effects are not used since
       the force field modes come from a real field.
       Check flag and skip in the future.

    */
    for (k3 = 0; k3 < numMeshPtsPerDir[2]; k3++) {
      for (k2 = 0; k2 < numMeshPtsPerDir[1]; k2++) {
        for (k1 = 0; k1 < numMeshPtsPerDir[0]; k1++) {

          I1 = (k3 * numMeshPtsPerDir[1] * numMeshPtsPerDir[0]) + (k2
              * numMeshPtsPerDir[0]) + k1;

          kk1 = (numMeshPtsPerDir[0] - k1) % numMeshPtsPerDir[0];
          kk2 = (numMeshPtsPerDir[1] - k2) % numMeshPtsPerDir[1];
          kk3 = (numMeshPtsPerDir[2] - k3) % numMeshPtsPerDir[2];

          I2 = (kk3 * numMeshPtsPerDir[1] * numMeshPtsPerDir[0]) + (kk2
              * numMeshPtsPerDir[0]) + kk1;

          for (d = 0; d < num_dim; d++) {

            /* set the real parts equal */
            fluidStochForceDensity_k[d][I1][0]
                                            = fluidStochForceDensity_k[d][I2][0];

            /* set the imaginary parts (conj to each other) */
            /* if the condition requires u_k to be conj to
               itself then set the imaginary part is zero,
               otherwise set the imaginary part conj. */
            if (I1 == I2) { /* self-conjugate mode */

              /* imaginary part set to zero */
              fluidStochForceDensity_k[d][I1][1] = 0;

            } else { /* non-self-conjugate mode */

              /* imaginary part set to conjugate */
              fluidStochForceDensity_k[d][I1][1]
                                              = -fluidStochForceDensity_k[d][I2][1];

            }

          } /* of d loop */

        } /* end k1 loop */
      } /* end k2 loop */
    } /* end k3 loop */

  /* end num_dim  == 3 */
  } else {
    printf("WARNING: %s : %s", error_str_code, error_str_func);
    printf("num_dim = %d \n", num_dim);
    printf("Stochastic field generation computation not implemented currently when num_dim is not 3D.\n");
  }


}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::projectField(int num_dim, double meshDeltaX,
                                                                    int *numMeshPtsPerDir,
                                                                    int shearDir, int shearVelDir,
                                                                    double shearDist,
                                                                    fftw_complex **field_u_k) {

  const char *error_str_func = "projectField()";

  int d,I,k1,k2,k3;
  
  double g_k_sq,g_k_real[3],g_k_imag[3],q_k_real[3],q_k_imag[3],q_k_sq;
  double vec_proj_real[3],vec_proj_imag[3];
  
  double L_shearDir,dot_gamma_t; 
  double normalizeFluidDensity_k;
  
  double conj_q_k_T_dot_w_k__real,conj_q_k_T_dot_w_k__imag,a_real,a_imag;

  if (num_dim == 3) {
    /* PJA: Check elsewhere */
    /*
     if (shearDir != 2) {
       printf("WARNING: %s : %s \n",error_str_code, error_str_func);
       printf("Only shearDir == 2 is currently debugged. \n");
     }
    */
    L_shearDir = numMeshPtsPerDir[shearDir]*meshDeltaX;
    dot_gamma_t = shearDist/L_shearDir;

    for (k3 = 0; k3 < numMeshPtsPerDir[2]; k3++) {
      for (k2 = 0; k2 < numMeshPtsPerDir[1]; k2++) {
        for (k1 = 0; k1 < numMeshPtsPerDir[0]; k1++) {

          I = (k3*numMeshPtsPerDir[1]*numMeshPtsPerDir[0]) 
            + (k2*numMeshPtsPerDir[0]) + k1;

          if ((k1 != 0) || (k2 != 0) || (k3 != 0)) {
            /* modify the projection method */

            /* note that g_k_real = 0 */
            g_k_imag[0] = sin(2.0*UNIT_pi*k1/numMeshPtsPerDir[0])/(meshDeltaX);
            g_k_imag[1] = sin(2.0*UNIT_pi*k2/numMeshPtsPerDir[1])/(meshDeltaX);
            g_k_imag[2] = sin(2.0*UNIT_pi*k3/numMeshPtsPerDir[2])/(meshDeltaX);

            g_k_sq = (g_k_imag[0]*g_k_imag[0]) + (g_k_imag[1]*g_k_imag[1]) + (g_k_imag[2]*g_k_imag[2]);

            /* construct the shear coordinate diverence operator D(t) Fourier representation 
               as D[w_k] = conj(q_k)^T*w_k.  
               
               This will allow us a natural interpretation in terms of projection, since 
               conj(q_k)^Tw_k = <q_k,w_k> dot-product for complex vectors.
               
               Note: the sign for q_k = -g_k, since g_k is purely imaginary and we want 
               divergence representated as a complex-valued dot-product.
            
            */
            q_k_imag[0] = -g_k_imag[0]; q_k_imag[1] = -g_k_imag[1]; q_k_imag[2] = -g_k_imag[2];

            q_k_imag[shearDir] += dot_gamma_t*g_k_imag[shearVelDir];  // note sign + from the conjuage

            q_k_sq = (q_k_imag[0]*q_k_imag[0]) + (q_k_imag[1]*q_k_imag[1]) + (q_k_imag[2]*q_k_imag[2]);

            /* use projection to enforce the incompressibility constraint on u_det_k_np1.

               Note that projection requires we compute

               v_k = (I - (q_k*conj(q_k)^T/\|q_k\|^2))*w_k. 

               This removes the compressible model q_k.
               
               Let q_k = \alpha + i*\beta = i*\beta, since \alpha = 0.
               w_k = \gamma + i*\lambda
               
               conj(q_k)^T*w_k = -i*\beta^T*\gamma + \beta^T*\lambda.
               
               Then use complex multiplication to combine into the projection
               expression to obtain the final v_k.
               
               Note:
               \alpha = 0, \beta = q_k_imag[]; \gamma = field_u_k[d][I][0], \lambda = field_u_k[d][I][1].

            */
            conj_q_k_T_dot_w_k__real = 0.0; conj_q_k_T_dot_w_k__imag = 0.0;            
            for (d = 0; d < num_dim; d++) {
              conj_q_k_T_dot_w_k__real += q_k_imag[d]*field_u_k[d][I][1];
              conj_q_k_T_dot_w_k__imag += -q_k_imag[d]*field_u_k[d][I][0];
            }
            
            // a = conj(q_k)^T*w_k/\|q_k\|^2.
            a_real = conj_q_k_T_dot_w_k__real/q_k_sq; 
            a_imag = conj_q_k_T_dot_w_k__imag/q_k_sq; 
                                                                                      
            for (d = 0; d < num_dim; d++) {
              vec_proj_real[d] = -q_k_imag[d]*a_imag; // q_k*conj(q_k)^T*w_k/\|q_k\|^2.
              field_u_k[d][I][0] -= vec_proj_real[d]; // (I - q_k*conj(q_k)^T*w_k/\|q_k\|^2)*w_k

              vec_proj_imag[d] = q_k_imag[d]*a_real;  // q_k*conj(q_k)^T*w_k/\|q_k\|^2.
              field_u_k[d][I][1] -= vec_proj_imag[d]; // (I - q_k*conj(q_k)^T*w_k/\|q_k\|^2)*w_k

            } /* end of d loop */

          } /* end of check not zero mode */

        } /* end k1 loop */
      } /* end k2 loop */
    } /* end k3 loop */

    /* end num_dim  == 3 */
  } else {
    printf("WARNING: %s : %s", error_str_code, error_str_func);
    printf("num_dim = %d \n", num_dim);
    printf("Incompressible computation not implemented currently for num_dim not 3D. \n");
  }

}


/* Compute the Fourier coefficients of the pressure field.
 This is obtain from:

 p_k = -(g_k^T/(|g_k|^2))*f_k

 NOTE: The pressure reported projects the force density,
 which corresponds to the pressure in the hydrodynamic
 equations.  This is NOT the same as the "pressure" obtained
 from the Hodge decomposition and imposition of incompressibility
 for the fluid velocity field, which is applied to the steady-state
 solution.  For the uniform mesh case we may perform the projection
 at either stage to enforce incompressibility by the fact that the
 projection operator and discrete Laplacian commute.

 */
void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::computePressure(int num_dim,
                                                                double meshDeltaX,
                                                                int *numMeshPtsPerDir,
                                                                int shearDir,
                                                                int shearVelDir,
                                                                double shearDist,
                                                                fftw_complex **forceDensity_k,
                                                                fftw_complex *pressure_k) {

  const char *error_str_func = "computePressure()";

  int d,I,k1,k2,k3;

  double q_k_sq;

  double q_k_dot_u_real_over_q_k_sq,q_k_dot_u_imag_over_q_k_sq;
  double q_k_dot_f_real_over_q_k_sq,q_k_dot_f_imag_over_q_k_sq;

  double normalizeFluidDensity_k;

  double L_shearDir;
  
  double g_k_real[3],g_k_imag[3],q_k_real[3],q_k_imag[3];
  double g_k_sq;
  double q_k_dot_u_real,q_k_dot_u_imag,q_k_dot_f_real,q_k_dot_f_imag;

  normalizeFluidDensity_k = 1;
  for (d = 0; d < num_dim; d++) {
    normalizeFluidDensity_k = normalizeFluidDensity_k * numMeshPtsPerDir[d];
  }

  if (num_dim == 3) {
    /* PJA: check elsewhere */
    /*
    if (shearDir != 2) {
      printf("WARNING: %s : %s \n",error_str_code, error_str_func);
      printf("Only shearDir == 2 is currently debugged. \n");
    }
    */
    L_shearDir = numMeshPtsPerDir[shearDir] * meshDeltaX;

    for (k3 = 0; k3 < numMeshPtsPerDir[2]; k3++) {
      for (k2 = 0; k2 < numMeshPtsPerDir[1]; k2++) {
        for (k1 = 0; k1 < numMeshPtsPerDir[0]; k1++) {

          I = (k3 * numMeshPtsPerDir[1] * numMeshPtsPerDir[0]) + (k2
              * numMeshPtsPerDir[0]) + k1;

          if ((k1 != 0) || (k2 != 0) || (k3 != 0)) {

            /* note that g_k_real = 0 */
            g_k_imag[0] = sin(2.0 * UNIT_pi * k1 / numMeshPtsPerDir[0])/ meshDeltaX;
            g_k_imag[1] = sin(2.0 * UNIT_pi * k2 / numMeshPtsPerDir[1])/ meshDeltaX;
            g_k_imag[2] = sin(2.0 * UNIT_pi * k3 / numMeshPtsPerDir[2])/ meshDeltaX;

            g_k_sq = (g_k_imag[0] * g_k_imag[0]) + (g_k_imag[1] * g_k_imag[1])
                   + (g_k_imag[2] * g_k_imag[2]);

            /* construct the shear coordinate divergence
               operator D(t) Fourier representation */
            q_k_imag[0] = g_k_imag[0];
            q_k_imag[1] = g_k_imag[1];
            q_k_imag[2] = g_k_imag[2];

            q_k_imag[shearDir] = g_k_imag[shearDir] - (shearDist / L_shearDir)
                               *g_k_imag[shearVelDir];

            q_k_sq = (q_k_imag[0] * q_k_imag[0]) + (q_k_imag[1] * q_k_imag[1])
                   + (q_k_imag[2] * q_k_imag[2]);

            q_k_dot_f_real = 0; /* real part of g^T f */
            for (d = 0; d < num_dim; d++) {
              q_k_dot_f_real += -q_k_imag[d] * forceDensity_k[d][I][1]
                                                                    / normalizeFluidDensity_k;
            }
            q_k_dot_f_real_over_q_k_sq = q_k_dot_f_real / q_k_sq;

            q_k_dot_f_imag = 0;
            for (d = 0; d < num_dim; d++) {
              q_k_dot_f_imag += q_k_imag[d] * forceDensity_k[d][I][0]
                                                                   / normalizeFluidDensity_k;
            }
            q_k_dot_f_imag_over_q_k_sq = q_k_dot_f_imag / q_k_sq;

            /* perform projection */
            /*
            for (d = 0; d < num_dim; d++) {
            vec_proj_real[d]            = -g_k_imag[d]*g_k_dot_f_imag_over_g_k_sq;
            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d][I][0]
            += vec_proj_real[d];

            vec_proj_imag[d]            = g_k_imag[d]*g_k_dot_u_real_over_g_k_sq;
            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d][I][1]
            += vec_proj_imag[d];

            }
            */

            pressure_k[I][0] = -q_k_dot_f_real_over_q_k_sq;

            pressure_k[I][1] = -q_k_dot_f_imag_over_q_k_sq;

            /* end check for not zero mode */
          } else { /* for zero mode do the following */
            pressure_k[I][0] = 0.0;
            pressure_k[I][1] = 0.0;
          }

        } /* end k1 loop */
      } /* end k2 loop */
    } /* end k3 loop */

    /* end num_dim  == 3 */
  } else {
    printf("WARNING: %s : %s", error_str_code, error_str_func);
    printf("num_dim = %d \n", num_dim);
    printf(
        "Incompressible computation not implemented currently for num_dim not 3D. \n");
  }

}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::bri1_unitCellRectImageShearPeriodic(double *periodL, double *meshCenterX0,
                                                                                    int shearDir, int shearVelDir, double shearDist,
                                                                                    int num_dim, double *X_orig, double *X_unitCell) {

  const char *error_str_func = "bri1_unitCellRectImageShearPeriodic()";

  IB_appl1_unitCellRectImageShearPeriodic(periodL, meshCenterX0,
                                          shearDir, shearVelDir, shearDist,
                                          num_dim, X_orig, X_unitCell);

}



void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::IB_appl1_unitCellRectImageShearPeriodic(double *periodL, double *meshCenterX0,
                                                                                        int shearDir, int shearVelDir,
                                                                                        double shearDist, int num_dim,
                                                                                        double *X_orig, double *X_unitCell) {

  const char *error_str_func = "IB_appl1_unitCellRectImageShearPeriodic()";

  int d;
  int k_shearDir, k_shearVelDir, k_shearRemDir;

  double XX1[3];

  double L_shearDir;
  double L_shearVelDir;
  double L_shearRemDir;

  double diffShearDir;
  double diffShearVelDir;
  double diffShearRemDir;

  int shearRemDir; /* remaining direction */

  /* We shall use perspective of the system described by
   * the rectangular unit cell which is used
   * to tile space with shifts in the shearVelDir in the
   * shearDir in the tiling.  See Shear paper discussion and figures.
   */

  if (shearDir == shearVelDir) {
    stringstream message;

    message << "The shear direction and shear velocity direction" << endl;
    message << "are not allowed to be in the same direction." << endl;
    message << endl;
    message << "Note: One way this error can occur is when one does" << endl;
    message << "not want to use the shear features in a simulation" << endl;
    message << "and simply sets the shear displacement and directions" << endl;
    message << "all set to zero.  Instead, make a valid choice for the" << endl;
    message << "shear and simply be sure to set the shearDist to zero." << endl;
    message << "For example, shearDir = 2, shearVelDir = 0, shearDist = 0." << endl;
    message << endl;
    message << "The values used in the calling routine were:" << endl;
    message << "  shearDir      = " << shearDir << endl;
    message << "  shearVelDir   = " << shearVelDir << endl;
    message << "  shearDist     = " << shearDist << endl;
    /* message << "  L_shearVelDir = " << L_shearVelDir << endl; */

    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  for (d = 0; d < num_dim; d++) {
    if ((d != shearDir) && (d != shearVelDir)) {
      shearRemDir = d;
    }
  }

  L_shearDir = periodL[shearDir];
  L_shearVelDir = periodL[shearVelDir];
  L_shearRemDir = periodL[shearRemDir];

  if (shearDist > L_shearVelDir) {
    stringstream message;
    message << "We assume that shearDist never exceeds domain length." << endl;
    message << "We require shearDist < L_shearVelDir." << endl;
    message << "  shearDist     = " << shearDist << endl;
    message << "  L_shearVelDir = " << L_shearVelDir << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  /* determine the unit cell image of X1 and call it XX1 */
  for (d = 0; d < num_dim; d++) {
    XX1[d] = X_orig[d];
  }

  diffShearDir = XX1[shearDir] - meshCenterX0[shearDir];
  k_shearDir = floor(fabs(diffShearDir / L_shearDir) + 0.5); /* round */
  if (diffShearDir < 0) {
    k_shearDir = -1 * k_shearDir;
  }

  X_unitCell[shearVelDir] = XX1[shearVelDir] - k_shearDir * shearDist; /* adjust cell alignments */
  X_unitCell[shearDir] = XX1[shearDir] - k_shearDir * L_shearDir; /* move point down cell rows */

  diffShearVelDir = X_unitCell[shearVelDir] - meshCenterX0[shearVelDir];
  k_shearVelDir = floor(fabs(diffShearVelDir / L_shearVelDir) + 0.5); /* round */
  if (diffShearVelDir < 0) {
    k_shearVelDir = -1 * k_shearVelDir;
  }

  X_unitCell[shearVelDir] = X_unitCell[shearVelDir] - k_shearVelDir
      * L_shearVelDir; /* move point to unit cell in shearVelDir */

  diffShearRemDir = XX1[shearRemDir] - meshCenterX0[shearRemDir];
  k_shearRemDir = floor(fabs(diffShearRemDir / L_shearRemDir) + 0.5); /* round */
  if (diffShearRemDir < 0) {
    k_shearRemDir = -1 * k_shearRemDir;
  }

  X_unitCell[shearRemDir] = XX1[shearRemDir] - k_shearRemDir * L_shearRemDir; /* move point to unit cell in shearVelDir */

  /* X_unitCell == XX1 */

}


/* == Compute the fluid velocity of the control points for FFTW3 uniform mesh
   with Lees-Edwards boundary conditions. */
void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::computeControlPtsVel_SHEAR_FFTW3(SELM_Lagrangian       *SELM_LagrangianData,
                                                                                        SELM_Eulerian         *SELM_EulerianData,
                                                                                        SELM_CouplingOperator *SELM_CouplingOperatorsData) {

  const char *error_str_func = "computeControlPtsVel_SHEAR_FFTW3()";

  int k;
  

  /* call the coupling operators */
  SELM_CouplingOperatorsData->computeOperator("GAMMA",
                                              (SELM_Lagrangian *) SELM_LagrangianData,
                                              (SELM_Eulerian *)   SELM_EulerianData);

  /* process the operator data according to the specific type of eulerian and lagrangian data */
  if (SELM_LagrangianData->type == SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_ANGLE_STYLE) {

    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE
      = (SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *)SELM_LagrangianData;

    /* copy the results of this operation to the control point velocities */
    for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numEntriesOpGammaVel; k++) {
      SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->pt_Vel[k] = SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->opGammaVel[k];
    }
    
  } else if (SELM_LagrangianData->type == SELM_Lagrangian_Types::TYPE_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE) {    
  
    SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE
      = (SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *)SELM_LagrangianData;

    /* copy the results of this operation to the control point velocities */
    for (k = 0; k < SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numEntriesOpGammaVel; k++) {
      SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->pt_Vel[k] = SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->opGammaVel[k];
    }
  
  } else if (SELM_LagrangianData->type == SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_STYLE_ELLIPSOID) {    
  
    SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID
      = (SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID *)SELM_LagrangianData;

    /* copy the results of this operation to the control point velocities */
    for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numEntriesOpGammaVel; k++) {
      SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->pt_Vel[k] = SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->opGammaVel[k];
    }

  } else {
    stringstream message;
    message << "No methods implemented to handle this type of Lagrangian data yet." << endl;
    message << "  typeStr = " << SELM_LagrangianData->typeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

}


/* setup for the time step (apply force to the mesh, etc...*/
void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::IB_appl1_start_time_step_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3() {

  const char *error_str_func = "IB_appl1_start_time_step_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3";

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3                       *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3;
  //SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE                         *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  //SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType *SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras;

  int I,j,d,N;
  int num_dim;
  int flagShearMode = -1;
  //double shearRate,ShearDist;
  //int shearDir,shearVelDir;

  #ifdef USE_PACKAGE_FFTW3

  /* == Dereference the data */
  if (driver_selm->SELM_Eulerian_List[0]->type == SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

    SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3    = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *) driver_selm->SELM_Eulerian_List[0];
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

    /* get shear data from appropriate data structure */
    num_dim                               = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

    /* @optimization: much of the below can be simplified here.  Just done for now for safety to sync data structures. */
    flagShearMode = SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode;
    if (flagShearMode == SHEAR_MODE_TYPE_RM_SHEAR1) {

      ShearData_RM_SHEAR1_Type *shearData_RM_SHEAR1 = (ShearData_RM_SHEAR1_Type *)SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

      /* set the shear information for the Eulerian DOF */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate = shearData_RM_SHEAR1->shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir = shearData_RM_SHEAR1->shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir = shearData_RM_SHEAR1->shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist = shearData_RM_SHEAR1->shearDist;

      /* use the shear information consistent with the Eulerian DOF */
      //shearRate = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate;
      //shearDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir;
      //shearVelDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir;
      //shearDist = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;
      
    } else if (flagShearMode == SHEAR_MODE_TYPE_RM_OSC1) {

      ShearData_RM_OSC1_Type *shearData_RM_OSC1 = (ShearData_RM_OSC1_Type *)SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->shearData;

      /* set the shear information for the Eulerian DOF */
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate = shearData_RM_OSC1->shearRate;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir = shearData_RM_OSC1->shearDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir = shearData_RM_OSC1->shearVelDir;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist = shearData_RM_OSC1->shearDist;

      /* use the shear information consistent with the Eulerian DOF */
      //shearRate = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate;
      //shearDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir;
      //shearVelDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir;
      //shearDist = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;      

    } else {
      stringstream message;
      message << "Expecting shear mode of type: " << SHEAR_MODE_TYPE_STR_RM_SHEAR1 << endl;
      message << "Instead shear mode was of type: " << SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearModeStr << endl;
      SELM_Package::packageError(error_str_code, error_str_func, message);
    }

  } else {
    stringstream message;
    message << "Expecting mesh of type: " << SELM_Eulerian_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3 << endl;
    message << "Instead mesh was used of type: " << driver_selm->SELM_Eulerian_List[0]->typeStr <<  endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  /* == Initialize the fluid force density to zero (do this only once) */
  N = 1;
  for (d = 0; d < SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
    N = N * SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d];
  }
  for (d = 0; d < SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
    for (j = 0; j < N; j++) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m[d][j][0] = 0.0;
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m[d][j][1] = 0.0;
    }
  }

  for (int couplingOpI = 0; couplingOpI < driver_selm->SELM_CouplingOperator_List_N; couplingOpI++) {

    /* == Assumes the control point forces were computed by LAMMPS already. */
    SELM_CouplingOperator *couplingOp = driver_selm->SELM_CouplingOperator_List[couplingOpI];

    /* maybe make below more generic */
    if (couplingOp->type == SELM_CouplingOperator_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1) {

      SELM_Lagrangian *lagrangian = NULL;
      SELM_Eulerian   *eulerian   = NULL;

      SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 *op
        = (SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 *)couplingOp;

      /* loop over the coupled entities */
      for (int k = 0; k < op->numCoupleList; k++) {
        lagrangian = op->lagrangianList[k];
        eulerian = op->eulerianList[k];

        /* == Compute the force density arising from the control points. */
        op->computeOperator("LAMBDA",
                            lagrangian,
                            eulerian);

        /* == Dereference the data */
        if (lagrangian->type == SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_ANGLE_STYLE) {

          SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE
            = (SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *) driver_selm->SELM_Lagrangian_List[k];

          /* get data from the time stepping data structures */
          //deltaT = this->SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->deltaT;
          
        } else if (lagrangian->type == SELM_Lagrangian_Types::TYPE_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE) {

          SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE
            = (SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *) driver_selm->SELM_Lagrangian_List[k];

          /* get data from the time stepping data structures */
          //deltaT = this->SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->deltaT;
          
        } else if (lagrangian->type == SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_STYLE_ELLIPSOID) {

          SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID
            = (SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID *) driver_selm->SELM_Lagrangian_List[k];

          /* get data from the time stepping data structures */
          //deltaT = this->SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->deltaT;          

        } else {
          stringstream message;
          message << "Expecting control points of type: "
                  << SELM_Lagrangian_Types::TYPE_STR_LAMMPS_ATOM_ANGLE_STYLE << endl;
          message << "Instead mesh was used of type: "
                  << lagrangian->typeStr << endl;
          SELM_Package::packageError(error_str_code, error_str_func, message);
        } /* end else */

      } /* end of k loop */

    } /* end if LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 */

  } /* end of loop over coupling operators */

#else

  stringstream message;
  message << "Compiler flag USE_PACKAGE_FFTW3 was not set." << endl;
  message << "None of the FFTW3 routines were not compiled in the codes." << endl;
  SELM_Package::packageError(error_str_code, error_str_func, message);

#endif


}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::writeSimulationDataToDisk(char *baseFilename, int timeIndex) {

  const char *error_str_func = "writeSimulationDataToDisk()";

  FILE *fid;
  char  filename[10000];

  /* open the file for writing the data */
  sprintf(filename,"%s_%.9d.SELM_Integrator_%s", baseFilename, timeIndex, typeStr);
  fid = fopen(filename,"w");

  if (fid == NULL) {

    stringstream message;
    message << "Could not open file, error occured." << endl;
    message << "  filename = " << filename << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

  }

  fprintf(fid, "-- SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3 : Simulation Data -- \n");
  fprintf(fid, "\n");

  fprintf(fid,"flagShearMode %d \n",
          SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras->flagShearMode);

  /* close the file */
  fclose(fid);
}


void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::syncShearDomainWithLammpsDomain() {

  const char *error_str_func = "syncShearDomainWithLammpsDomain()";
  
  Domain *domain = lammps->domain;
  
  double  shearRate, shearDist, meshDeltaX;
  int     num_dim;
  int    *numMeshPtsPerDir;
  double *meshCenterX0;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3                                                       *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  /* get the triclinic flag : 0 - orthogonal, 1 - triclinic box */
  
  /* get data from the integrator */
  /* == Dereference the data */
  if (driver_selm->SELM_Eulerian_List[0]->type == SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

    SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3
      = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *)driver_selm->SELM_Eulerian_List[0];

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
      = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

    num_dim          = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;
    meshDeltaX       = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
    numMeshPtsPerDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
    meshCenterX0     = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0;

    shearRate        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate;
    shearDist        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;
           
  } else {
    stringstream message;
    message << "Expecting mesh type: %s." << endl;
    message << "  SELM_Eulerian_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3" << endl;
    message << "Instead mesh type was: " << endl;
    message << "  " << driver_selm->SELM_Eulerian_List[0]->typeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  /* only perform synchronization if the shearRate is non-zero */
  if (shearRate != 0) {

    if (domain->triclinic == 0) {
      stringstream message;
      message << "For simulations with a deforming box domain (Lees-Edwards conditions)" << endl;
      message << "There was a non-zero shear rate and shear distance used in integrator" << endl;
      message << "while the LAMMPS codes had an orthogonal box specified." << endl;
      message << "To use this feature LAMMPS must us a triclinic box." << endl;
      message << "lammps->domain->triclinic = " << domain->triclinic << endl;
      SELM_Package::packageError(error_str_code, error_str_func, message);
    }

  } else {

    /* check the box is undeformed in LAMMPS, if undeformed in SELM */
    if (shearDist == 0) {

      if (domain->triclinic == 0) {
        if ((domain->xy != 0.0) && (domain->xz != 0.0) && (domain->yz != 0.0)) {
          stringstream message;
          message << "For simulations with a deforming box domain (Lees-Edwards conditions)" << endl;
          message << "The SELM and LAMMPS deformation must be setup to be the same." << endl;
          message << "This requires shearDist, shearVelDir, shearDir, be consistent with" << endl;
          message << "the xy, xz, yz in LAMMPS" << endl;
          message << "If this feature is not to be used, then setup an orthogonal box in LAMMPS" << endl;
          SELM_Package::packageError(error_str_code, error_str_func, message);
        }
      }  /* end triclinic == 0 */

    } /* end shearDist == 0 */

  } /* end of else */

}





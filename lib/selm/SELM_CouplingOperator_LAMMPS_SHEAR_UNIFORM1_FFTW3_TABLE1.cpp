/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */

/* standard includes */
#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include "string.h"
#include <malloc.h>

/* SELM includes */
#include "wrapper_selm.h"
#include "SELM_Parser1.h"

#include "SELM_CouplingOperator.h"
#include "SELM_CouplingOperator_Types.h"
#include "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1.h"

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_Types.h"
#include "SELM_Lagrangian_CONTROLPTS_BASIC1.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE.h"
#include "SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID.h"

#include "SELM_Eulerian_Types.h"
#include "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.h"

/* LAMMPS includes */
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace std;

/* ============ constants definition =========== */
const int   SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::TYPE     = SELM_CouplingOperator_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1;
const char *SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::TYPE_STR = SELM_CouplingOperator_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1;


/* coupling operator type */
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::COUPLING_OP_TYPE_STR_NULL     = "NULL";
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::COUPLING_OP_TYPE_STR_GAMMA    = "GAMMA";
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::COUPLING_OP_TYPE_STR_LAMBDA   = "LAMBDA";
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::COUPLING_OP_TYPE_STR_UPSILON  = "UPSILON";

/* operator type */
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_NULL        = "NULL";
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_T_KERNEL_1  = "T_KERNEL_1";
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_T_FAXEN_1   = "T_FAXEN_1";
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_TR_FAXEN_1  = "TR_FAXEN_1";

double SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::UNIT_pi                      = 3.141592653589793;

/* error information */
const char* SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::error_str_code                = "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1.cpp";

/* ---------------------------------------------------------------------- */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::initConstants() {

  /* set the constants for the operator types */
  /*
  const char OPERATOR_TYPE_STR_NULL[]       = "NULL";
  const int  OPERATOR_TYPE_NULL             = 0;

  const char OPERATOR_TYPE_STR_T_KERNEL_1[] = "T_KERNEL_1";
  const int  OPERATOR_TYPE_T_KERNEL_1       = 1;

  const char OPERATOR_TYPE_STR_T_FAXEN_1[]  = "T_FAXEN_1";
  const int  OPERATOR_TYPE_T_FAXEN_1        = 2;

  const char OPERATOR_TYPE_STR_TR_FAXEN_1[] = "TR_FAXEN_1";
  const int  OPERATOR_TYPE_TR_FAXEN_1       = 3;
  */

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::init() {

  /* initialize the constants for this class */
  initConstants();

  /* initialize the variables for this class */
  type                       = SELM_CouplingOperator_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1;
  strcpy(typeStr, SELM_CouplingOperator_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1);

  numEntriesOpGammaResults   = 0;
  opGammaResults             = NULL;

  numEntriesOpLambdaResults  = 0;
  opLambdaResults            = NULL;

  numEntriesOpUpsilonResults = 0;
  opUpsilonResults           = NULL;

  operatorType               = OPERATOR_TYPE_NULL;
  strcpy(operatorTypeStr, SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_NULL);
  operatorData               = NULL;

  flagWriteSimulationData    = 0;
  saveSkipSimulationData     = 0;

}


SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1() : SELM_CouplingOperator()
{

  /* generic initializations for this class */
  init();

}

SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1(int narg, char **arg) : SELM_CouplingOperator(narg, arg)
{

  const char *error_str_func = "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1()";

  /* generic initializations for this class */
  init();

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::parse_ParameterFile(const char *baseFilename) {

  const char *error_str_func = "parse_ParameterFile()";

  char filename[10000];
  char weightTableFilename[10000];

  FILE *fid;

  char tempStr[10000];
  char c;
  int  k;

  /* open the file */
  sprintf(filename, "%s.SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1", baseFilename);
  fid = fopen(filename, "r");

  if (fid == 0) {
    stringstream message;
	message << "Could not open file, error occured." << endl;
	message << "  filename = %s" << filename << endl;
	SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  /* read the first line of the text file. */
  /* this contains comments (# symbol in given example file). */
  c = 0;
  while (c != '\n') {
    fscanf(fid, "%c", &c);
  }

  /* read the second line of the text file. */
  /* this contains comments (# symbol in given example file). */
  c = 0;
  while (c != '\n') {
    fscanf(fid, "%c", &c);
  }

  /* read the operator type */
  fscanf(fid,"%s",  (char *)tempStr);
  fscanf(fid,"%s",  operatorTypeStr);

  operatorType = getOperatorTypeFromStr(operatorTypeStr);

  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1:

    /* read the weight filename */
    fscanf(fid,"%s",  (char *)tempStr);
    fscanf(fid,"%s",  weightTableFilename);

    /* @@@ out of date use of data structures */
    //SELM_weightTable = NULL;
    //readWeightTable(weightTableFilename,
    //                &SELM_weightTable);
    break;

  default:
    stringstream message;
    message << "Invalid operator type was specified." << endl;
	message << "operatorTypeStr = " << operatorTypeStr << endl;
	message << "May not be implemented yet" << endl;
	SELM_Package::packageError(error_str_code, error_str_func, message);

  } /* end of switch */



  /* close the file */
  fclose(fid);

}

int SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::getOperatorTypeFromStr(const char *operatorTypeStr) {

  int operatorType_local;

  if (strcmp(operatorTypeStr, SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_T_KERNEL_1) == 0) {
    operatorType_local = SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_T_KERNEL_1;
  }

  if (strcmp(operatorTypeStr, SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_T_FAXEN_1) == 0) {
    operatorType_local = SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_T_FAXEN_1;
  }

  if (strcmp(operatorTypeStr, SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_STR_TR_FAXEN_1) == 0) {
    operatorType_local = SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_TR_FAXEN_1;
  }

  return operatorType_local;
}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::setup() {

  const char *error_str_func = "setup()";

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::packageError(int code, void *extras) {
    exit(code);
}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperator(const char      *couplingOpTypeStr,
                                                                               SELM_Lagrangian *SELM_LagrangianData,
                                                                               SELM_Eulerian   *SELM_EulerianData) {

  const char *error_str_func = "computeOperator(const char *, SELM_Lagrangian*, SELM_Eulerian*)";

  int flagNotFoundYet;
  int couplingOpType;
  int flagHandledOp = 0;

  /* ---- determine the operator type ---- */
  couplingOpType    = COUPLING_OP_TYPE_NULL;
  flagNotFoundYet   = 1;

  if ((flagNotFoundYet) && (strcmp(couplingOpTypeStr, COUPLING_OP_TYPE_STR_GAMMA)) == 0) {
    couplingOpType          = COUPLING_OP_TYPE_GAMMA;
    flagNotFoundYet = 0;
  }

  if ((flagNotFoundYet) && (strcmp(couplingOpTypeStr, COUPLING_OP_TYPE_STR_LAMBDA)) == 0) {
    couplingOpType          = COUPLING_OP_TYPE_LAMBDA;
    flagNotFoundYet = 0;
  }

  if ((flagNotFoundYet) && (strcmp(couplingOpTypeStr, COUPLING_OP_TYPE_STR_UPSILON)) == 0) {
    couplingOpType          = COUPLING_OP_TYPE_UPSILON;
    flagNotFoundYet = 0;
  }

  /* ------- branch based on the operator type to compute ---- */
  switch (couplingOpType) {

  case COUPLING_OP_TYPE_GAMMA:
    computeOperatorGamma(SELM_LagrangianData,
                         SELM_EulerianData);
    flagHandledOp = 1;
    break;

  case COUPLING_OP_TYPE_LAMBDA:
    computeOperatorLambda(SELM_LagrangianData,
                          SELM_EulerianData);
    flagHandledOp = 1;
    break;

    /*
  case COUPLING_OP_TYPE_UPSILON:
    computeOperatorUpsilon(SELM_LagrangianData,
                           SELM_EulerianData);
    flagHandledOp = 1;
    break;
    */

  } /* ---------- end of switch ----------- */

  /* if pair is not compatible */
  if (flagHandledOp != 1) {
    stringstream message;
    message << "Coupling operator not supported by this class." << endl;
    message << "couplingOpTypeStr = " << couplingOpTypeStr << endl;
    message << "LagrangianTypeStr = " << SELM_LagrangianData->typeStr << endl;
    message << "EulerianTypeStr   = " << SELM_EulerianData->typeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

}




void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorGamma(SELM_Lagrangian *SELM_LagrangianData, SELM_Eulerian *SELM_EulerianData) {

  const char *error_str_func = "computeOperatorGamma(SELM_Lagrangian*, SELM_Eulerian*)";
  int flagHandledPair = 0;

  /* -------------- determine the control points type -------------- */
  switch (SELM_LagrangianData->type) {

  case SELM_Lagrangian_Types::TYPE_CONTROLPTS_BASIC1:

    /* -------------- determine the mesh type -------------- */
    switch (SELM_EulerianData->type) {

    case SELM_Eulerian_Types::TYPE_FLUID_SHEAR_UNIFORM1_FFTW3:
      computeOperatorGamma((SELM_Lagrangian_CONTROLPTS_BASIC1 *) SELM_LagrangianData,
                           (SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3 *) SELM_EulerianData);
      flagHandledPair = 1;
      break;

    } /* ---------- end of switch ----------- */

    break;

    case SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_ANGLE_STYLE:

      /* -------------- determine the mesh type -------------- */
      switch (SELM_EulerianData->type) {

      case SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3:
        computeOperatorGamma((SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *) SELM_LagrangianData,
                             (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *) SELM_EulerianData);
        flagHandledPair = 1;
        break;

      } /* ---------- end of switch ----------- */

      break;
      
    case SELM_Lagrangian_Types::TYPE_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE:

      /* -------------- determine the mesh type -------------- */
      switch (SELM_EulerianData->type) {

      case SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3:
        computeOperatorGamma((SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *) SELM_LagrangianData,
                             (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *) SELM_EulerianData);
        flagHandledPair = 1;
        break;

      } /* ---------- end of switch ----------- */

      break;      

    case SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_STYLE_ELLIPSOID:

      /* -------------- determine the mesh type -------------- */
      switch (SELM_EulerianData->type) {

      case SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3:
        computeOperatorGamma((SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID *) SELM_LagrangianData,
                             (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *) SELM_EulerianData);
        flagHandledPair = 1;
        break;

      } /* ---------- end of switch ----------- */

      break;

  } /* ---------- end of switch ----------- */

  /* report error if this pair combination is unhandled */
  if (flagHandledPair != 1) {
    stringstream message;
    message << "Coupling operator not supported yet for this." << endl;
    message << "combination of Lagrangian and Eulerian degrees of freedom." << endl;
    message << "LagrangianTypeStr = " << SELM_LagrangianData->typeStr << endl;
    message << "EulerianTypeStr   = " << SELM_EulerianData->typeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorGamma(SELM_Lagrangian_CONTROLPTS_BASIC1         *SELM_LagrangianData_CONTROLPTS_BASIC1,
                                                                                    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3) {

  const char *error_str_func = "computeOperatorGamma(SELM_Lagrangian_CONTROLPTS_BASIC1*, SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3*)";

  int     num_dim;
  int     d;
  int     N;

  int     numPts;
  double *X_list;

  stringstream message;

  operatorDataType_T_KERNEL_1 *opData;

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ExtrasType  *SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras;

  IB_appl1_computeSmoothedVelFieldExtrasType *IB_appl1_computeSmoothedVelFieldExtras;

  SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras
    = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  /* setup the smoother data */
  IB_appl1_computeSmoothedVelFieldExtras
  = (IB_appl1_computeSmoothedVelFieldExtrasType *)
  malloc(sizeof(IB_appl1_computeSmoothedVelFieldExtrasType));

  for (d = 0; d < num_dim; d++) {
    IB_appl1_computeSmoothedVelFieldExtras->meshCenterX0[d]     = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0[d];
    IB_appl1_computeSmoothedVelFieldExtras->numMeshPtsPerDir[d] = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d];
  }
  IB_appl1_computeSmoothedVelFieldExtras->meshDeltaX            = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  IB_appl1_computeSmoothedVelFieldExtras->num_dim               = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  IB_appl1_computeSmoothedVelFieldExtras->operatorType
  = operatorType;

  IB_appl1_computeSmoothedVelFieldExtras->u_m
  = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m; /* yes, verified vel is correct */

  numPts = SELM_LagrangianData_CONTROLPTS_BASIC1->numControlPts;
  X_list = SELM_LagrangianData_CONTROLPTS_BASIC1->pt_X;

  /* Computes the different table based
   * operators Gamma.  This includes Gamma for
   * translational and rotational degrees of freedom.
   * The behavior depends on the initialization of
   * this class. */
  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1: /* Basic point-like kernel for control points (translation only) */

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    IB_appl1_computeSmoothedVelFieldExtras->SELM_weightTable
      = opData->weightTable;

    if (SELM_LagrangianData_CONTROLPTS_BASIC1->opGammaVel == NULL) {

      N                                                           = num_dim*numPts;
      SELM_LagrangianData_CONTROLPTS_BASIC1->numEntriesOpGammaVel = N;
      SELM_LagrangianData_CONTROLPTS_BASIC1->opGammaVel           = (double *)malloc(sizeof(double)*N);

    } else { /* check array of the correct size */

      /* if size mis-match re-allocate the storage array for the operator results */
      if (SELM_LagrangianData_CONTROLPTS_BASIC1->numEntriesOpGammaVel != num_dim*numPts) {
        free(SELM_LagrangianData_CONTROLPTS_BASIC1->opGammaVel);
        N                                                           = num_dim*numPts;
        SELM_LagrangianData_CONTROLPTS_BASIC1->numEntriesOpGammaVel = N;
        SELM_LagrangianData_CONTROLPTS_BASIC1->opGammaVel           = (double *)malloc(sizeof(double)*N);
      }

    } /* end else */

    /* Compute the kernel times the velocity field at each point. */
    IB_appl1_computeSmoothedVelField(num_dim, numPts, X_list,
                                     IB_appl1_computeSmoothedVelFieldExtras,
                                     &SELM_LagrangianData_CONTROLPTS_BASIC1->opGammaVel);

    break;

  case OPERATOR_TYPE_T_FAXEN_1:

    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation only)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

    break;

  case OPERATOR_TYPE_TR_FAXEN_1:
    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation and rotation)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  default:
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  } /* end of switch for operator type */

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorGamma(SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE,
                                                        SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

  const char *error_str_func = "computeOperatorGamma(SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE*, SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3*)";

  int     num_dim;
  int     d;
  int     N;

  int     numPts;
  double *X_list;

  stringstream message;

  operatorDataType_T_KERNEL_1 *opData;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  IB_appl1_computeSmoothedVelFieldExtrasType *IB_appl1_computeSmoothedVelFieldExtras;

  SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
    = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  /* setup the smoother data */
  IB_appl1_computeSmoothedVelFieldExtras
  = (IB_appl1_computeSmoothedVelFieldExtrasType *)
  malloc(sizeof(IB_appl1_computeSmoothedVelFieldExtrasType));

  for (d = 0; d < num_dim; d++) {
    IB_appl1_computeSmoothedVelFieldExtras->meshCenterX0[d]     = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0[d];
    IB_appl1_computeSmoothedVelFieldExtras->numMeshPtsPerDir[d] = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d];
  }
  IB_appl1_computeSmoothedVelFieldExtras->meshDeltaX            = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  IB_appl1_computeSmoothedVelFieldExtras->num_dim               = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  IB_appl1_computeSmoothedVelFieldExtras->operatorType
  = operatorType;

  IB_appl1_computeSmoothedVelFieldExtras->u_m
  = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m; /* yes, verified vel is correct */

  numPts = SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numControlPts;
  X_list = SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX;

  /* Computes the different table based
   * operators Gamma.  This includes Gamma for
   * translational and rotational degrees of freedom.
   * The behavior depends on the initialization of
   * this class. */
  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1: /* Basic point-like kernel for control points (translation only) */

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    IB_appl1_computeSmoothedVelFieldExtras->SELM_weightTable
      = opData->weightTable;

    if (SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->opGammaVel == NULL) {

      N                                                                 = num_dim*numPts;
      SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numEntriesOpGammaVel = N;
      SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->opGammaVel           = (double *)malloc(sizeof(double)*N);

    } else { /* check array of the correct size */

      /* if size mis-match re-allocate the storage array for the operator results */
      if (SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numEntriesOpGammaVel != num_dim*numPts) {
        free(SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->opGammaVel);
        N                                                                 = num_dim*numPts;
        SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numEntriesOpGammaVel = N;
        SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->opGammaVel           = (double *)malloc(sizeof(double)*N);
      }

    } /* end else */

    /* Compute the kernel times the velocity field at each point. */
    IB_appl1_computeSmoothedVelField(num_dim, numPts, X_list,
                                     IB_appl1_computeSmoothedVelFieldExtras,
                                     &SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->opGammaVel);

    break;

  case OPERATOR_TYPE_T_FAXEN_1:

    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation only)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

    break;

  case OPERATOR_TYPE_TR_FAXEN_1:
    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation and rotation)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  default:
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  } /* end of switch for operator type */

  /* free the dynamically allocated data structure */
  free(IB_appl1_computeSmoothedVelFieldExtras);

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorGamma(SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE         *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE,
                                                        SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

  const char *error_str_func = "computeOperatorGamma(SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE*, SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3*)";

  int     num_dim;
  int     d;
  int     N;

  int     numPts;
  double *X_list;

  stringstream message;

  operatorDataType_T_KERNEL_1 *opData;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  IB_appl1_computeSmoothedVelFieldExtrasType *IB_appl1_computeSmoothedVelFieldExtras;

  SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
    = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  /* setup the smoother data */
  IB_appl1_computeSmoothedVelFieldExtras
  = (IB_appl1_computeSmoothedVelFieldExtrasType *)
  malloc(sizeof(IB_appl1_computeSmoothedVelFieldExtrasType));

  for (d = 0; d < num_dim; d++) {
    IB_appl1_computeSmoothedVelFieldExtras->meshCenterX0[d]     = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0[d];
    IB_appl1_computeSmoothedVelFieldExtras->numMeshPtsPerDir[d] = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d];
  }
  IB_appl1_computeSmoothedVelFieldExtras->meshDeltaX            = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  IB_appl1_computeSmoothedVelFieldExtras->num_dim               = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  IB_appl1_computeSmoothedVelFieldExtras->operatorType
  = operatorType;

  IB_appl1_computeSmoothedVelFieldExtras->u_m
  = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m; /* yes, verified vel is correct */

  numPts = SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numControlPts;
  X_list = SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX;

  /* Computes the different table based
   * operators Gamma.  This includes Gamma for
   * translational and rotational degrees of freedom.
   * The behavior depends on the initialization of
   * this class. */
  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1: /* Basic point-like kernel for control points (translation only) */

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    IB_appl1_computeSmoothedVelFieldExtras->SELM_weightTable
      = opData->weightTable;

    if (SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->opGammaVel == NULL) {

      N                                                                 = num_dim*numPts;
      SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numEntriesOpGammaVel = N;
      SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->opGammaVel           = (double *)malloc(sizeof(double)*N);

    } else { /* check array of the correct size */

      /* if size mis-match re-allocate the storage array for the operator results */
      if (SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numEntriesOpGammaVel != num_dim*numPts) {
        free(SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->opGammaVel);
        N                                                                 = num_dim*numPts;
        SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numEntriesOpGammaVel = N;
        SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->opGammaVel           = (double *)malloc(sizeof(double)*N);
      }

    } /* end else */

    /* Compute the kernel times the velocity field at each point. */
    IB_appl1_computeSmoothedVelField(num_dim, numPts, X_list,
                                     IB_appl1_computeSmoothedVelFieldExtras,
                                     &SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->opGammaVel);

    break;

  case OPERATOR_TYPE_T_FAXEN_1:

    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation only)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

    break;

  case OPERATOR_TYPE_TR_FAXEN_1:
    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation and rotation)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  default:
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  } /* end of switch for operator type */

  /* free the dynamically allocated data structure */
  free(IB_appl1_computeSmoothedVelFieldExtras);

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorGamma(SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID         *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID,
                                                        SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

  const char *error_str_func = "computeOperatorGamma(SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID*, SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3*)";

  int     num_dim;
  int     d;
  int     N;

  int     numPts;
  double *X_list;

  stringstream message;

  operatorDataType_T_KERNEL_1 *opData;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  IB_appl1_computeSmoothedVelFieldExtrasType *IB_appl1_computeSmoothedVelFieldExtras;

  SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
    = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  /* setup the smoother data */
  IB_appl1_computeSmoothedVelFieldExtras
  = (IB_appl1_computeSmoothedVelFieldExtrasType *)
  malloc(sizeof(IB_appl1_computeSmoothedVelFieldExtrasType));

  for (d = 0; d < num_dim; d++) {
    IB_appl1_computeSmoothedVelFieldExtras->meshCenterX0[d]     = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0[d];
    IB_appl1_computeSmoothedVelFieldExtras->numMeshPtsPerDir[d] = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d];
  }
  IB_appl1_computeSmoothedVelFieldExtras->meshDeltaX            = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  IB_appl1_computeSmoothedVelFieldExtras->num_dim               = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  IB_appl1_computeSmoothedVelFieldExtras->operatorType
  = operatorType;

  IB_appl1_computeSmoothedVelFieldExtras->u_m
  = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m; /* yes, verified vel is correct */

  numPts = SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numControlPts;
  X_list = SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->ptsX;

  /* Computes the different table based
   * operators Gamma.  This includes Gamma for
   * translational and rotational degrees of freedom.
   * The behavior depends on the initialization of
   * this class. */
  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1: /* Basic point-like kernel for control points (translation only) */

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    IB_appl1_computeSmoothedVelFieldExtras->SELM_weightTable
      = opData->weightTable;

    if (SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->opGammaVel == NULL) {

      N                                                                 = num_dim*numPts;
      SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numEntriesOpGammaVel = N;
      SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->opGammaVel           = (double *)malloc(sizeof(double)*N);

    } else { /* check array of the correct size */

      /* if size mis-match re-allocate the storage array for the operator results */
      if (SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numEntriesOpGammaVel != num_dim*numPts) {
        free(SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->opGammaVel);
        N                                                                 = num_dim*numPts;
        SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numEntriesOpGammaVel = N;
        SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->opGammaVel           = (double *)malloc(sizeof(double)*N);
      }

    } /* end else */

    /* Compute the kernel times the velocity field at each point. */
    IB_appl1_computeSmoothedVelField(num_dim, numPts, X_list,
                                     IB_appl1_computeSmoothedVelFieldExtras,
                                     &SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->opGammaVel);

    break;

  case OPERATOR_TYPE_T_FAXEN_1:

    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation only)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

    break;

  case OPERATOR_TYPE_TR_FAXEN_1:
    /* Faxel-like kernel on sphere based on Lebedev quadratures (translation and rotation)*/
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  default:
    message << "Invalid operator type specified." << endl;
    message << "operatorTypeStr = " << operatorTypeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    break;

  } /* end of switch for operator type */

  /* free the dynamically allocated data structure */
  free(IB_appl1_computeSmoothedVelFieldExtras);

}




void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorLambda(SELM_Lagrangian *SELM_LagrangianData,
                                                         SELM_Eulerian   *SELM_EulerianData) {

  const char *error_str_func = "computeOperatorLambda(SELM_Lagrangian*, SELM_Eulerian*)";

  int          flagHandledPair = 0;

  /* -------------- determine the lagrangian degrees of freedom type -------------- */
  switch (SELM_LagrangianData->type) {

  case SELM_Lagrangian_Types::TYPE_CONTROLPTS_BASIC1:

    /* -------------- determine the eulerian degrees of freedom type -------------- */
    switch (SELM_EulerianData->type) {

    case SELM_Eulerian_Types::TYPE_FLUID_SHEAR_UNIFORM1_FFTW3:
      computeOperatorLambda((SELM_Lagrangian_CONTROLPTS_BASIC1 *)        SELM_LagrangianData,
                            (SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3 *) SELM_EulerianData);
      flagHandledPair = 1;
      break;

    } /* ---------- end of switch ----------- */

    break;

    case SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_ANGLE_STYLE:

      /* -------------- determine the eulerian degrees of freedom type -------------- */
      switch (SELM_EulerianData->type) {

      case SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3:
        computeOperatorLambda((SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *)  SELM_LagrangianData,
                              (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *) SELM_EulerianData);
        flagHandledPair = 1;
        break;

      } /* ---------- end of switch ----------- */

      break;

  } /* ---------- end of switch ----------- */

  /* if pair is not compatible */
  if (flagHandledPair != 1) {
    stringstream message;
    message << "Coupling operator is not supported yet for this" << endl;
    message << "combination of Lagrangian and Eulerian degrees of freedom." << endl;
    message << "LagrangianTypeStr = " << SELM_LagrangianData->typeStr << endl;
    message << "EulerianTypeStr   = " << SELM_EulerianData->typeStr << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorLambda(SELM_Lagrangian_CONTROLPTS_BASIC1         *SELM_LagrangianData_CONTROLPTS_BASIC1,
                                                         SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3) {

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras;

  int I;
  int i, j, k, d;
  int N;
  int flagDebugTests = 0;
  int flagDebugStillCompIBForces = 0;

  double *meshCenterX0;
  int *numMeshPtsPerDir;
  double L_shearDir;

  int shearDir;
  int shearVelDir;
  double shearRate;
  double shearDist;

  int num_dim;

  double deltaT;

  double meshDeltaX;

  /* shift the control points into the sheared coordinate system */
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras = SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim          = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  meshDeltaX       = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  numMeshPtsPerDir = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  meshCenterX0     = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0;

  shearRate        = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearRate;
  shearDir         = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearDir;
  shearVelDir      = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir;
  shearDist        = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;

  L_shearDir       = numMeshPtsPerDir[shearDir] * meshDeltaX;



  for (k = 0; k < SELM_LagrangianData_CONTROLPTS_BASIC1->numControlPts; k++) {
    SELM_LagrangianData_CONTROLPTS_BASIC1->pt_X[(k * num_dim) + shearVelDir] += -(shearDist
                                                                   / L_shearDir) * (SELM_LagrangianData_CONTROLPTS_BASIC1->pt_X[(k * num_dim) + shearDir]
                                                                   - meshCenterX0[shearDir]);
  } /* end of k loop */

  IB_appl1_applyControlPtsForceToMesh_FFTW3(SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                            SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                            SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                                            SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                            SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m,
                                            SELM_LagrangianData_CONTROLPTS_BASIC1);

  /* shift the control points back into physical coordinates */
  for (k = 0; k < SELM_LagrangianData_CONTROLPTS_BASIC1->numControlPts; k++) {
    SELM_LagrangianData_CONTROLPTS_BASIC1->pt_X[(k * num_dim) + shearVelDir] += (shearDist
                                                                   / L_shearDir) * (SELM_LagrangianData_CONTROLPTS_BASIC1->pt_X[(k * num_dim) + shearDir]
                                                                   - meshCenterX0[shearDir]);
  } /* end of k loop */

  /* == Compute any additional force density contributions acting on the fluid. */
  /*  (*IB_appl1Extras->userAppl_computeFluidForcesFunc)(IB_appl1Extras->userAppl_computeFluidForcesExtras); */

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorLambda(SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE   *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE,
                                                         SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  int I;
  int i, j, k, d;
  int N;
  int flagDebugTests = 0;
  int flagDebugStillCompIBForces = 0;

  double *meshCenterX0;
  int    *numMeshPtsPerDir;
  double  L_shearDir;

  int     shearDir;
  int     shearVelDir;
  double  shearRate;
  double  shearDist;

  int num_dim;

  double deltaT;

  double meshDeltaX;

  /* shift the control points into the sheared coordinate system */
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
    = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim          = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  meshDeltaX       = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  numMeshPtsPerDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  meshCenterX0     = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0;

  shearRate        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate;
  shearDir         = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir;
  shearVelDir      = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir;
  shearDist        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;

  L_shearDir       = numMeshPtsPerDir[shearDir] * meshDeltaX;

  for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numControlPts; k++) {
    SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX[(k * num_dim) + shearVelDir] += -(shearDist
                                                    / L_shearDir) * (SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX[(k * num_dim) + shearDir]
                                                    - meshCenterX0[shearDir]);
  } /* end of k loop */

  IB_appl1_applyControlPtsForceToMesh_FFTW3(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m,
                                            SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE);

  /* shift the control points back into physical coordinates */
  for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numControlPts; k++) {
    SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX[(k * num_dim) + shearVelDir] += (shearDist
                                                                   / L_shearDir) * (SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX[(k * num_dim) + shearDir]
                                                                   - meshCenterX0[shearDir]);
  } /* end of k loop */

  /* == Compute any additional force density contributions acting on the fluid. */
  /*  (*IB_appl1Extras->userAppl_computeFluidForcesFunc)(IB_appl1Extras->userAppl_computeFluidForcesExtras); */

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorLambda(SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE  *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE,
                                                         SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  int I;
  int i, j, k, d;
  int N;
  int flagDebugTests = 0;
  int flagDebugStillCompIBForces = 0;

  double *meshCenterX0;
  int    *numMeshPtsPerDir;
  double  L_shearDir;

  int     shearDir;
  int     shearVelDir;
  double  shearRate;
  double  shearDist;

  int num_dim;

  double deltaT;

  double meshDeltaX;

  /* shift the control points into the sheared coordinate system */
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
    = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim          = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  meshDeltaX       = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  numMeshPtsPerDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  meshCenterX0     = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0;

  shearRate        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate;
  shearDir         = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir;
  shearVelDir      = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir;
  shearDist        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;

  L_shearDir       = numMeshPtsPerDir[shearDir] * meshDeltaX;

  for (k = 0; k < SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numControlPts; k++) {
    SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX[(k * num_dim) + shearVelDir] += -(shearDist
                                                    / L_shearDir) * (SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX[(k * num_dim) + shearDir]
                                                    - meshCenterX0[shearDir]);
  } /* end of k loop */

  IB_appl1_applyControlPtsForceToMesh_FFTW3(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m,
                                            SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE);

  /* shift the control points back into physical coordinates */
  for (k = 0; k < SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numControlPts; k++) {
    SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX[(k * num_dim) + shearVelDir] += (shearDist
                                                                   / L_shearDir) * (SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX[(k * num_dim) + shearDir]
                                                                   - meshCenterX0[shearDir]);
  } /* end of k loop */

  /* == Compute any additional force density contributions acting on the fluid. */
  /*  (*IB_appl1Extras->userAppl_computeFluidForcesFunc)(IB_appl1Extras->userAppl_computeFluidForcesExtras); */

}



void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorLambda(SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID  *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID,
                                                         SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3) {

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  int I;
  int i, j, k, d;
  int N;
  int flagDebugTests = 0;
  int flagDebugStillCompIBForces = 0;

  double *meshCenterX0;
  int    *numMeshPtsPerDir;
  double  L_shearDir;

  int     shearDir;
  int     shearVelDir;
  double  shearRate;
  double  shearDist;

  int num_dim;

  double deltaT;

  double meshDeltaX;

  /* shift the control points into the sheared coordinate system */
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
    = SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras;

  num_dim          = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim;

  meshDeltaX       = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
  numMeshPtsPerDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  meshCenterX0     = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0;

  shearRate        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate;
  shearDir         = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir;
  shearVelDir      = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir;
  shearDist        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;

  L_shearDir       = numMeshPtsPerDir[shearDir] * meshDeltaX;

  for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numControlPts; k++) {
    SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->ptsX[(k * num_dim) + shearVelDir] += -(shearDist
                                                    / L_shearDir) * (SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->ptsX[(k * num_dim) + shearDir]
                                                    - meshCenterX0[shearDir]);
  } /* end of k loop */

  IB_appl1_applyControlPtsForceToMesh_FFTW3(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m,
                                            SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID);

  /* shift the control points back into physical coordinates */
  for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numControlPts; k++) {
    SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->ptsX[(k * num_dim) + shearVelDir] += (shearDist
                                                                   / L_shearDir) * (SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->ptsX[(k * num_dim) + shearDir]
                                                                   - meshCenterX0[shearDir]);
  } /* end of k loop */

  /* == Compute any additional force density contributions acting on the fluid. */
  /*  (*IB_appl1Extras->userAppl_computeFluidForcesFunc)(IB_appl1Extras->userAppl_computeFluidForcesExtras); */

}


void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::computeOperatorUpsilon(SELM_Lagrangian *SELM_LagrangianData) {

}



/*
==============================================================================
 void IB_appl1_compute_SELM_WEIGHT_FUNC1(vecX, extras)
------------------------------------------------------------------------------

 Compute the SELM weight function from table for the given mesh shift.
 This is useful in computing nearly invariant kernal functions from
 tabulated data.

==============================================================================
*/
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_compute_SELM_WEIGHT_FUNC1(int                               num_dim,
                                                                       int                               numPts,
                                                                       double                           *X_list,
                                                                       double                            deltaX,
                                                                       controlPts_SELM_weightTableType  *SELM_weightTable,
                                                                       int                              *numEval,
                                                                       double                          **eval_ptr) {

  const char *error_str_func = "IB_appl1_compute_SELM_WEIGHT_FUNC1()";

  int    I, j, k, d;

  int    index;
  double r[3];

  int    N;

  double val[3];
  double X[3];
  double X_cm[3];

  double meshDeltaX;

  double *eval;

  /* check eval array is large enough */
  if ((*eval_ptr) == NULL) {
    (*numEval)  = numPts;
    (*eval_ptr) = (double *)malloc(sizeof(double)*numPts);
  }

  if ((*numEval) != numPts) {

    stringstream message;
    message << "  Evaluation array allocated is not large enough." << endl;
    message << "  Could indicate evaluation array was not allocated."  << endl;
    message << "  If the eval. array is set to NULL this routine" << endl;
    message << "  will allocate memory for it." << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

  }


  /* compute the delta function at the specific collection of points */
  eval = (*eval_ptr);
  for (I = 0; I < numPts; I++) {

    /* get the displacement from the center of mass */
    for (d = 0; d < num_dim; d++) {
      X[d]    = X_list[I*num_dim + d];
      X_cm[d] = 0.0; /* we assume the X_list already has the X = X_m - X_cm computed */
    } /* end of d loop */

    N          = 1;
    meshDeltaX = deltaX;
    weightFromTable(num_dim, N, X,
                                    X_cm, meshDeltaX,
                                    SELM_weightTable,
                                    eval_ptr);

  } /* end of I loop */

}




void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::weightFromTable(int num_dim, int numPts, double *ptsX,
                                                                    double *X_cm, double meshDeltaX,
                                                                    controlPts_SELM_weightTableType *SELM_weightTable,
                                                                    double **weightFuncVals_ptr) {

  /* Note: we assume the table has ordering with indexing
   *
   * index = j1 + j2*N + j3*N*N.
   *
   */
  const char *error_str_func = "weightFromTable()";

  double X[3];
  double X0[3];
  double X1[3];
  double alpha[3];
  double *dX;
  int    *vecN;

  double  R_0;

  double *X_cm_list;

  int     e1,e2,e3;
  int     vec_e[3];

  int    vecJ[3];
  int    vecJ0[3];

  double sumVal;

  int    i, j, k, d;
  int    I;
  int    index;

  double val;
  double w[3];

  double weight;

  double norm_X_sq;

  double *weightFuncVals;

  /* dereference values */
  dX        = SELM_weightTable->dX;
  vecN      = SELM_weightTable->vecN;
  R_0       = SELM_weightTable->R_0;
  X_cm_list = SELM_weightTable->X_list;

  if ((num_dim != SELM_weightTable->num_dim) || (num_dim != 3)) {
    stringstream message;
    message << "Only implemented for three dimensions." << endl;
    message << "SELM_weightTable->num_dim = " << SELM_weightTable->num_dim << endl;
    message << "num_dim = " << num_dim << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  if (*weightFuncVals_ptr == NULL) {
    weightFuncVals = (double *)malloc(sizeof(double)*numPts);
  } else {
    weightFuncVals = *weightFuncVals_ptr;
  }

  /* -- loop over the points and compute the weights from table */
  for (I = 0; I < numPts; I++) {

    /* compute the displacement */
    for (d = 0; d < num_dim; d++) {
      /* non-dimensionalize and use lattice symmetries to put in first octant. */
      X[d] = fabs(ptsX[I*num_dim + d] - X_cm[d])/meshDeltaX;
    }

    /* compute the norm of X */
    norm_X_sq = 0;
    for (d = 0; d < num_dim; d++) {
      norm_X_sq += X[d]*X[d];
    }

    if (norm_X_sq >= R_0*R_0) {

      weight = 0; /* assume compact support of weight function */

    } else {

      /* -- Perform interpolation from the weight table
       * We use the lattice symmetry and always consider
       * X in the first octant.
       */

      /* find X0 and weights */
      for (d = 0; d < num_dim; d++) {
        vecJ0[d] = floor(X[d]/dX[d]);
      }

      index      = vecJ0[0] + vecJ0[1]*vecN[0] + vecJ0[2]*vecN[0]*vecN[1];

      for (d = 0; d < num_dim; d++) {
        X0[d]    = X_cm_list[index*num_dim + d];
        alpha[d] = (X[d] - X0[d])/dX[d];
      }

      sumVal = 0.0;
      for (e1 = 0; e1 <= 1; e1++) {
        for (e2 = 0; e2 <= 1; e2++) {
          for (e3 = 0; e3 <= 1; e3++) {

            vec_e[0] = e1; vec_e[1] = e2; vec_e[2] = e3;

            for (d = 0; d < num_dim; d++) {
              vecJ[d] = vecJ0[d] + vec_e[d];
            }

            index    = vecJ[0] + vecJ[1]*vecN[0] + vecJ[2]*vecN[0]*vecN[1];
            val      = SELM_weightTable->weight_X_list[index];

            for (d = 0; d < num_dim; d++) {
              X1[d] = SELM_weightTable->X_list[index*num_dim + d];
              w[d]  = vec_e[d]*alpha[d] + (1 - vec_e[d])*(1 - alpha[d]);
            }

            sumVal += val*w[0]*w[1]*w[2];

          } /* e1 loop */
        } /* e2 loop */
      } /* e3 loop */

      weight = sumVal; /* give weights */

    } /* end of norm check */

    weightFuncVals[I] = weight;

  } /* end of I loop */

  /* -- return the values */
  (*weightFuncVals_ptr) = weightFuncVals;

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::readWeightTable(const char *filename,
                                                                               controlPts_SELM_weightTableType **SELM_weightTable_ptr) {

  const char *error_str_func = "readSELM_weightTable()";

  FILE                            *fid;
  controlPts_SELM_weightTableType *SELM_weightTable;

  char                             c;

  int                   i,j,k,d;
  int                   num_dim;
  int                   numTableVals;

  char                  tempStr[10000];

  if (*SELM_weightTable_ptr == NULL) {
    SELM_weightTable = (controlPts_SELM_weightTableType *)malloc(sizeof(controlPts_SELM_weightTableType));
  } else {
    SELM_weightTable = (*SELM_weightTable_ptr);
  }

  /* open the file for parsing */
  fid = fopen(filename,"r");

  if (fid == 0) {
    stringstream message;
    message << "Could not open file, error occured." << endl;
    message << "  filename = " << filename <<  endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
    packageError(1, 0);
  }

  /* parse comments here */
  /* read the first line of the text file. */
  /* this contains comments (# symbol in given example file). */
  c = 0;
  while (c != '\n') {
    fscanf(fid, "%c", &c);
  }

  /* read the second line of the text file. */
  /* this contains comments (# symbol in given example file). */
  c = 0;
  while (c != '\n') {
    fscanf(fid, "%c", &c);
  }

  fscanf(fid, "%s", (char *) tempStr);
  fscanf(fid, "%s", (char *) tempStr); 
  strcpy(SELM_weightTable->name,tempStr); // PJA: trying to avoid char* issue --> map to static length string

  fscanf(fid, "%s", (char *) tempStr);
  fscanf(fid, "%d", &SELM_weightTable->num_dim);
  num_dim = SELM_weightTable->num_dim;

  fscanf(fid, "%s", (char *) tempStr);
  fscanf(fid, "%lf", &SELM_weightTable->R_0);

  fscanf(fid, "%s", (char *) tempStr);
  for (d = 0; d < num_dim; d++) {
    fscanf(fid, "%d", &SELM_weightTable->vecN[d]);
  }

  fscanf(fid, "%s", (char *) tempStr);
  for (d = 0; d < num_dim; d++) {
    fscanf(fid, "%lf", &SELM_weightTable->dX[d]);
  }

  fscanf(fid, "%s", (char *) tempStr);
  fscanf(fid, "%d", &SELM_weightTable->numTableVals);
  numTableVals                    = SELM_weightTable->numTableVals;

  SELM_weightTable->X_list        = (double *)malloc(sizeof(double)*num_dim*numTableVals);
  fscanf(fid, "%s", (char *) tempStr);
  for (k = 0; k < numTableVals*num_dim; k++) {
    fscanf(fid, "%lf", &SELM_weightTable->X_list[k]);
  } /* end of k loop */

  SELM_weightTable->weight_X_list = (double *)malloc(sizeof(double)*num_dim*numTableVals);
  fscanf(fid, "%s", (char *) tempStr);
  for (k = 0; k < numTableVals*num_dim; k++) {
    fscanf(fid, "%lf", &SELM_weightTable->weight_X_list[k]);
  } /* end of k loop */

  /* return the results */
  (*SELM_weightTable_ptr) = SELM_weightTable;
}




/* Computes the function required inside the spherical average
 * for the translational motion of the sphere.
 */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_userFunc_TR_PARTICLE_Vel_sphFunc_Xcm(int num_dim, int numPts,
                                                                                 double *X_list, void *userExtras,
                                                                                 int *funcVal_num_dim,
                                                                                 double **funcVal) {

  const char *error_str_func = "IB_appl1_userFunc_TR_PARTICLE_Vel_sphFunc_Xcm()";

  int k, d;

  int deltaFuncType;

  double X[3];
  double Xcm[3];

  IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_XcmType *extras =
      (IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_XcmType *) userExtras;

  IB_appl1_computeSmoothedVelFieldExtrasType *smoothedVelExtras =
      extras->smoothedVelExtras;

  /* Compute the kernel times the velocity field at each point. */
  IB_appl1_computeSmoothedVelField(num_dim, numPts, X_list,
                                   smoothedVelExtras, funcVal);

  (*funcVal_num_dim) = 3;

}

/* Force density associated with the translation motion of the sphere.
 */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_userFunc_TR_PARTICLE_Force_sphFunc_Xcm(int      num_dim,
                                                                                    int      numPts,
                                                                                    double  *X_list,
                                                                                    void    *userExtras,
                                                                                    int     *funcVal_num_dim,
                                                                                    double **funcVal) {

  const char *error_str_func = "IB_appl1_userFunc_TR_PARTICLE_Force_sphFunc_Xcm()";

  int    I, j, k, d;

  int    deltaFuncType;

  int    numX = 1;
  double X[3];
  double X_m[3];
  double controlPtForce[3];

  double deltaX;
  double meshDeltaX;
  int     flagKernalType;

  int     numEval = 0;

  double  meshFactor;

  double  eval;
  double *eval_ptr;

  /* controlPts_SELM_weightTableType *SELM_weightTable; */

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType *extras =
    (IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType *) userExtras;

  IB_appl1_computeSmoothedForceFieldExtrasType *smoothedForceExtras =
    extras->smoothedForceExtras;

  for (d = 0; d < num_dim; d++) {
    X_m[d]            = extras->X_m[d];
    controlPtForce[d] = extras->controlPtForce[d];
  }

  deltaX = extras->deltaX;

  /* Compute the kernel times the velocity field at each point. */
  /*deltaFuncType = CONTROL_PT_TYPE_IBPK_DELTA_FUNC; */

  if ((*funcVal) == NULL) {
    (*funcVal_num_dim) = num_dim;
    (*funcVal)         = (double *)malloc(sizeof(double)*numPts*num_dim);
  }

  if ((*funcVal_num_dim) != num_dim) {
    stringstream message;
    message << "  funcVal_num != num_dim" << endl;
    message << "Could indicate array for funcVal not allocated." << endl;
    message << "If funcVal== NULL we allocate memory for it here." << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  /* loop over the points and compute the function values to be averaged
   * over the surface of the sphere.
   */
  for (I = 0; I < numPts; I++) {

    numX = 1;
    for (d = 0; d < num_dim; d++) {
      X[d] = X_m[d] - X_list[I*num_dim + d];
    }

    numEval        = 1;
    eval_ptr       = &eval;
    IB_appl1_compute_SELM_WEIGHT_FUNC1(num_dim,
                                       numX,
                                       X,
                                       deltaX,
                                       smoothedForceExtras->SELM_weightTable,
                                       &numEval,
                                       &eval_ptr);

    /* weights include meshDeltaX term implicitly
     * w(r)*dX^3 = eval(r)
     */
    meshDeltaX = deltaX;
    meshFactor = 1.0/(meshDeltaX*meshDeltaX*meshDeltaX);

    for (d = 0; d < num_dim; d++) {
      (*funcVal)[I*num_dim + d] = eval*meshFactor*controlPtForce[d];
    }

  } /* end of I loop */

}



/* Computes the function required inside the spherical average
 * for the translational motion of the sphere.
 */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_userFunc_TR_PARTICLE_Vel_sphFunc_Theta(int      num_dim, int   numPts,
                                                                                    double  *X_list,  void *userExtras,
                                                                                    int     *funcVal_num_dim,
                                                                                    double **funcVal_ptr) {

  const char *error_str_func = "IB_appl1_userFunc_TR_PARTICLE_Vel_sphFunc_Theta()";

  int     k, d;

  double  X[3];
  double  V[3];
  double  *X_cm;

  double *vel_w;
  double *funcVal;

  double  sphereR;
  double  preFactor;

  IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_ThetaType *extras =
      (IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_ThetaType *) userExtras;

  IB_appl1_computeSmoothedVelFieldExtrasType *smoothedVelExtras =
      extras->smoothedVelExtras;

  /* get the SELM point data */
  sphereR = extras->sphereR;
  X_cm    = extras->X_cm;

  if ((*funcVal_ptr) == NULL) {
    (*funcVal_ptr)     = (double *)malloc(sizeof(double)*numPts*num_dim);
    (*funcVal_num_dim) = num_dim;
  }

  funcVal = (*funcVal_ptr);

  /* Compute the kernel times the velocity field at each point. */
  vel_w         = NULL;
  IB_appl1_computeSmoothedVelField(num_dim, numPts, X_list,
                                   smoothedVelExtras, &vel_w);

  /* now compute z x w which gives the integrand in the average */
  preFactor = 3.0/(2.0*sphereR*sphereR);
  for (k = 0; k < numPts; k++) {

    for (d = 0; d < num_dim; d++) {
      X[d] = X_list[k*num_dim + d] - X_cm[d];
      V[d] = vel_w[k*num_dim + d];
    }

    /* f = z x w */
    funcVal[k*num_dim + 0] = X[1]*V[2] - X[2]*V[1];
    funcVal[k*num_dim + 1] = X[2]*V[0] - X[0]*V[2];
    funcVal[k*num_dim + 2] = X[0]*V[1] - X[1]*V[0];

    /* f = 3/(2*R^2)*<...>*/
    for (d = 0; d < num_dim; d++) {
      funcVal[k*num_dim + d] = preFactor*funcVal[k*num_dim + d];
    }

  } /* end of k loop */

  free(vel_w);

  /* return the results */
  (*funcVal_num_dim) = 3;
  (*funcVal_ptr)     = funcVal;

}


/* Force density associated with the rotational motion of the sphere.
 */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_userFunc_TR_PARTICLE_Force_sphFunc_Theta(int      num_dim,
                                                                                      int      numPts,
                                                                                      double  *X_list,
                                                                                      void    *userExtras,
                                                                                      int     *funcVal_num_dim,
                                                                                      double **funcVal_ptr) {

  const char *error_str_func = "IB_appl1_userFunc_TR_PARTICLE_Force_sphFunc_Theta()";

  int    I, j, k, d;

  int    deltaFuncType;

  int    numX = 1;
  double X[3];
  double X_m[3];

  double X_cm[3];
  double Z[3];

  double controlPtTorque[3];

  double deltaX;
  double sphereR;

  double preFactor;

  int     numEval = 0;

  int     flagKernalType;

  double  *funcVal;

  double  eval;
  double *eval_ptr;

  double  meshFactor;
  double  meshDeltaX;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType *extras =
    (IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType *) userExtras;

  IB_appl1_computeSmoothedForceFieldExtrasType *smoothedForceExtras =
    extras->smoothedForceExtras;

  sphereR = extras->sphereR;
  for (d = 0; d < num_dim; d++) {
    X_cm[d]            = extras->X_cm[d];
    X_m[d]             = extras->X_m[d];
    controlPtTorque[d] = extras->controlPtTorque[d];
  }

  preFactor = -3.0/(2.0*sphereR*sphereR);

  deltaX = extras->deltaX;

  /* Compute the kernel times the velocity field at each point. */
  /*deltaFuncType = CONTROL_PT_TYPE_IBPK_DELTA_FUNC; */
  funcVal = (*funcVal_ptr);
  if (funcVal == NULL) {
    (*funcVal_num_dim) = num_dim;
    funcVal            = (double *)malloc(sizeof(double)*numPts*num_dim);
    (*funcVal_ptr)     = funcVal;
  }

  if ((*funcVal_num_dim) != num_dim) {
    stringstream message;
    message << "  funcVal_num != num_dim" << endl;
    message << "Could indicate array for funcVal not allocated." << endl;
    message << "If funcVal== NULL we allocate memory for it here." << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  /* loop over the points and compute the function values to be averaged
   * over the surface of the sphere.
   */
  for (I = 0; I < numPts; I++) {

    numX = 1;
    for (d = 0; d < num_dim; d++) {
      X[d] = X_m[d] - X_list[I*num_dim + d];
      Z[d] = X_list[I*num_dim + d] - X_cm[d];
    }

    numEval        = 1;
    eval_ptr       = &eval;

    IB_appl1_compute_SELM_WEIGHT_FUNC1(num_dim,
                                       numX,
                                       X,
                                       deltaX,
                                       smoothedForceExtras->SELM_weightTable,
                                       &numEval,
                                       &eval_ptr);

    /* factor needed here since, w(x)dX^3 = eval(x) (weight includes volume element) */
    meshDeltaX = deltaX;
    meshFactor = 1.0/(meshDeltaX*meshDeltaX*meshDeltaX);

    /* compute f = eta * (z x Torque) */
    funcVal[I*num_dim + 0] = Z[1]*controlPtTorque[2] - Z[2]*controlPtTorque[1];
    funcVal[I*num_dim + 1] = Z[2]*controlPtTorque[0] - Z[0]*controlPtTorque[2];
    funcVal[I*num_dim + 2] = Z[0]*controlPtTorque[1] - Z[1]*controlPtTorque[0];

    for (d = 0; d < num_dim; d++) {
      funcVal[I*num_dim + d] = meshFactor*preFactor*eval*funcVal[I*num_dim + d];
    }

  } /* end of I loop */

  /* return the data */
  (*funcVal_num_dim) = num_dim;
  (*funcVal_ptr)     = funcVal;

}



/* Computes values of a smoothed velocity field */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_computeSmoothedVelField(int   num_dim, int numPts, double *X_list,
                                                                                                void *userExtras, double **Vel_list) {

  const char *error_str_func = "IB_appl1_computeSmoothedVelField()";

  int     i, j, k, d;

  double  meshDeltaX;
  double  L;
  double  meshBaseX0[3];
  double  meshCenterX0[3];
  double  meshX0[3];
  double *meshSELMX;

  int     localMeshSize;
  int     numMeshPtsPerDir[3];

  int     numMeshIndices;
  int     flagDebugSELM;

  double *meshVel;
  double *meshX;

  int     meshJ[3];
  int     meshI0[3];

  double  R_0;


  double  weightFuncSum;
  double *weightFuncVals;

  int     N;
  int     j1, j2, j3;
  int     jj1, jj2, jj3;

  int     I, I1, I2;
  int     index;

  double  ptX[3];

  int     pt_size = 1;

  int     numQuadNodes;
  int     ptI;
  double  integralEstimate[3];
  double  ptVel[3];

  int     meshI;
  int    *meshSELMIndices;

  int     numDir;
  int    *meshSELMJ;
  int     numMeshSELMIndices;

  fftw_complex **u_m;

  IB_appl1_computeSmoothedVelFieldExtrasType *extras =
      (IB_appl1_computeSmoothedVelFieldExtrasType *)userExtras;

  controlPts_SELM_weightTableType *SELM_weightTable;


  /* == initialize variables */
  num_dim    = extras->num_dim;
  meshDeltaX = extras->meshDeltaX;

  for (d = 0; d < num_dim; d++) {
    numMeshPtsPerDir[d] = extras->numMeshPtsPerDir[d];
    meshCenterX0[d]     = extras->meshCenterX0[d];
  }

  u_m = extras->u_m;

  if ((*Vel_list) == NULL) {
    (*Vel_list) = (double *) malloc(sizeof(double)*numPts*num_dim);
  }

  /* loop over the points and compute the velocity value */
  for (ptI = 0; ptI < numPts; ptI++) {

    for (d = 0; d < num_dim; d++) {
      ptX[d] = X_list[(ptI * num_dim) + d];
    }

    /* == compute the value of the smoothed velocity field
     */
    switch (extras->operatorType) {

    case OPERATOR_TYPE_T_KERNEL_1:
    /* case BRI1_SELM_TR_PARTICLE_WEIGHT_FUNC1_KERNAL_TYPE: */

      /* -- compute the weights on a patch around the control point */
      SELM_weightTable
      = extras->SELM_weightTable;

      /* == Compute the patch around the center of mass. */

      /* -- determine the mesh indices */

      /* perform index arithmetic to find a mesh point X0 nearest to X */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
        meshI0[d]      = floor(((ptX[d] - meshBaseX0[d])/meshDeltaX) + 0.5);
        meshX0[d]      = meshBaseX0[d] + meshI0[d]*meshDeltaX;
      } /* end of d loop */

      /* now determine the indices for all mesh points with square
       * patch of size R_0
       */
      R_0                  = SELM_weightTable->R_0;
      numDir               = 2*ceil(R_0) + 1;
      numMeshSELMIndices   = numDir*numDir*numDir;
      meshSELMIndices      = (int *)malloc(sizeof(int)*numMeshSELMIndices);
      meshSELMJ            = (int *)malloc(sizeof(int)*numMeshSELMIndices*num_dim);

      meshSELMX            = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);

      /* compute base value for the mesh points */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
      }

      for (j1 = 0; j1 < numDir; j1++) {
        for (j2 = 0; j2 < numDir; j2++) {
          for (j3 = 0; j3 < numDir; j3++) {

            index    = j1 + j2*numDir + j3*numDir*numDir;

            meshJ[0] = meshI0[0] + j1 - ((numDir - 1)/2);
            meshJ[1] = meshI0[1] + j2 - ((numDir - 1)/2);
            meshJ[2] = meshI0[2] + j3 - ((numDir - 1)/2);

            /* compute the index modulo numMeshPtsPerDir */
            for (d = 0; d < num_dim; d++) {

              /* compute the mesh point locations before the adjustments are made to make it mod L */
              meshSELMX[index*num_dim + d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;

              /* make the indices mod N */
              meshSELMJ[index*num_dim + d] = (meshJ[d] % numMeshPtsPerDir[d]);
              if (meshSELMJ[index*num_dim + d] < 0) {
                meshSELMJ[index*num_dim + d] = meshSELMJ[index*num_dim + d] + numMeshPtsPerDir[d];
              }

            } /* end of d loop */

            /* compute the index on the mesh */
            meshI                    = meshSELMJ[index*num_dim + 0]
                                     + meshSELMJ[index*num_dim + 1]*numMeshPtsPerDir[0]
                                     + meshSELMJ[index*num_dim + 2]*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];
            meshSELMIndices[index]   = meshI;

          } /* end of j3 loop */
        } /* end of j2 loop */
      } /* end of j1 loop */

      /* -- compute the weights on the patch of mesh points */
      weightFuncVals = NULL;
      weightFromTable(num_dim, numMeshSELMIndices, meshSELMX,
                                      ptX, meshDeltaX,
                                      SELM_weightTable,
                                      &weightFuncVals);

      weightFuncSum = 0;
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncSum += weightFuncVals[I];
      }

      /* renormalize the weight values, since they should sum to one */
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncVals[I] = weightFuncVals[I]/weightFuncSum;
      }

      /* -- compute the average velocity by looping over the mesh */
      for (d = 0; d < num_dim; d++) {
        ptVel[d] = 0.0;
      }

      for (I = 0; I < numMeshSELMIndices; I++) {

        meshI = meshSELMIndices[I];

        for (d = 0; d < num_dim; d++) {
          ptVel[d] += u_m[d][meshI][0]*weightFuncVals[I];
        }

      } /* end of I loop */

      flagDebugSELM = 0;
      if (flagDebugSELM) {
        printf("WARNING: %s : %s \n", error_str_code, error_str_func);
        printf("Debug information being printed. \n");
        printf(" \n");
        fflush(stdout);

        printf("ptX_raw = [");
        for (I = 0; I < num_dim; I++) {
          printf("%0.8g", ptX[I]);
          if (I < num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

        printf("meshSELMX_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMX[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

        printf("weightFuncVals_raw = [");
        for (I = 0; I < numMeshSELMIndices; I++) {
          printf("%0.8g", weightFuncVals[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

      } /* end of flagDebugSELM */

      /* free dynamic memory */
      free(weightFuncVals);
      free(meshSELMIndices);
      free(meshSELMJ);
      free(meshSELMX);

      break;


    default:

      stringstream message;
      message << "  funcVal_num != num_dim" << endl;
      message << "  The smoothing kernal type was unrecognized." << endl;
      message << "  operatorTypeStr = " << operatorTypeStr << endl;
      message << "  operatorType    = " << operatorType << endl;

      SELM_Package::packageError(error_str_code, error_str_func, message);

    } /* end of switch */

    /* store the computed velocity value */
    for (d = 0; d < num_dim; d++) {
      (*Vel_list)[(ptI*num_dim) + d] = ptVel[d];
    }

  } /* end of ptI loop */

}


/*
 ==============================================================================
 void computeLebedevSphereAvg()
 ------------------------------------------------------------------------------

 Computes the average of a function over the surface of a sphere using the
 Lebedev quadratures.

 V.I. Lebedev, and D.N. Laikov, "A quadrature formula for the sphere of
 the 131st algebraic order of accuracy", Doklady Mathematics, Vol. 59,
 No. 3, 1999, pp. 477-481.

 ==============================================================================
 */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_computeLebedevSphereAvg(int numQuadNodes, double sphereR, int num_dim,
                                                                    double *Xcm,
                                                                    void   (*userFunc)(int num_dim, int numQuadNodes,
                                                                    double *nodeX, void *userExtras,
                                                                    int *funcVal_num_dim,
                                                                    double **funcVal),
                                                                    void *userExtras,
                                                                    int  *integralVal_num_dim, double **integralEstimate_ptr) {

  const char *error_str_func = "IB_appl1_computeLebedevSphereAvg()";

  int     I, i, j, k, d;

  double  integralEstimate = 0.0;

  double *nodeX;
  double *X;

  double *weight;
  double *funcVal;
  int     funcVal_num_dim;

  int     flagDEBUG;

  /* == Get the Lebedev quadrature nodes and weights */
  nodeX  = NULL;
  weight = NULL;
  IB_appl1_getLebedevSphereQuad(numQuadNodes, &nodeX, &weight);

  X = (double *)malloc(sizeof(double)*numQuadNodes*num_dim);

  /* Compute the quadrature nodes on the surface of the sphere
   * of radius R centered at the point X_cm.
   */
  for (j = 0; j < numQuadNodes; j++) {
    for (d = 0; d < num_dim; d++) {
      X[j*num_dim + d] = Xcm[d] + (nodeX[j*num_dim + d]*sphereR);
    }
  }

  flagDEBUG = 0;
  if (flagDEBUG) {
    printf("WARNING: %s : %s \n", error_str_code, error_str_func);
    printf("Debug information being printed. \n");
    printf(" \n");
    fflush(stdout);

    printf("weight_raw = [");
    for (I = 0; I < numQuadNodes; I++) {
      printf("%0.8g", weight[I]);
      if (I < numQuadNodes - 1) {
        printf(", ");
      }
    }
    printf("];\n");
    fflush(stdout);

    printf("nodeX_raw = [");
    for (I = 0; I < numQuadNodes*num_dim; I++) {
      printf("%0.8g", nodeX[I]);
      if (I < numQuadNodes*num_dim - 1) {
        printf(", ");
      }
    }
    printf("];\n");
    fflush(stdout);

    printf("nodeXadj_raw = [");
    for (I = 0; I < numQuadNodes*num_dim; I++) {
      printf("%0.8g", X[I]);
      if (I < numQuadNodes*num_dim - 1) {
        printf(", ");
      }
    }
    printf("];\n");
    fflush(stdout);

  } /* flagDEBUG */

  /* == Compute the user defined function at the adjusted quadrature nodes X*/
  funcVal = NULL;
  (*userFunc)(num_dim, numQuadNodes, X, userExtras, &funcVal_num_dim,
              &funcVal);

  if ((*integralEstimate_ptr) == NULL) {
    (*integralVal_num_dim)  = funcVal_num_dim;
    (*integralEstimate_ptr) = (double *) malloc(sizeof(double)
                                                * (*integralVal_num_dim));
  }

  if ((*integralVal_num_dim) != funcVal_num_dim) {
    stringstream message;
    message << "  funcVal_num != num_dim" << endl;
    message << "  The expected number of dimensions for function" << endl;
    message << "  evaluation and those actually computed differ." << endl;
    message << "  integralVal_num_dim = " << integralVal_num_dim  << endl;
    message << "  funcVal_num_dim     = " << funcVal_num_dim      << endl;

    SELM_Package::packageError(error_str_code, error_str_func, message);

  }

  /* == Compute the quadrature */

  /* compute an estimate of the integral over the sphere */
  for (d = 0; d < funcVal_num_dim; d++) {

    integralEstimate = 0.0;
    for (k = 0; k < numQuadNodes; k++) {
      integralEstimate += weight[k]*funcVal[(k*funcVal_num_dim) + d];
    } /* end of k loop */

    /* == Return the estimate of the integral */
    (*integralEstimate_ptr)[d] = integralEstimate/(4*UNIT_pi);

  } /* end of d loop */

  /* compute the average of the integral over the sphere */
  /* integralEstimate/(4*UNITS1_pi*sphereR*sphereR); */

  /* free dynamic memory */
  free(X);
  free(nodeX);
  free(weight);

  free(funcVal);

}

/* Computes the Lebedev quadrature */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_getLebedevSphereQuad(int numQuadNodes, double **nodeX_ptr,
                                                                 double **weight_ptr) {

  const char *error_str_func = "IB_appl1_getLebedevSphereQuad()";

  int numQuadNodes_table;

  int num_dim = 3;

  double *nodeX;
  double *weight;

  if ((*nodeX_ptr) == NULL) {
    /* allocate memory */
    (*nodeX_ptr) = (double *) malloc(sizeof(double) * numQuadNodes * num_dim);
  }

  if ((*weight_ptr) == NULL) {
    /* allocate memory */
    (*weight_ptr) = (double *) malloc(sizeof(double) * numQuadNodes);
  }

  nodeX = (*nodeX_ptr);
  weight = (*weight_ptr);

  /*****************************************************/
  /* Table of quadrature weights for given numQuadNodes
   /*****************************************************/

  switch (numQuadNodes) {

  /*---------------------------------------------------*/
  case 14:
    numQuadNodes_table = 14;

    weight[0] = 0.837758041;
    nodeX[0] = 1;
    nodeX[1] = 0;
    nodeX[2] = 0;

    weight[1] = 0.837758041;
    nodeX[3] = -1;
    nodeX[4] = 0;
    nodeX[5] = 0;

    weight[2] = 0.837758041;
    nodeX[6] = 0;
    nodeX[7] = 1;
    nodeX[8] = 0;

    weight[3] = 0.837758041;
    nodeX[9] = 0;
    nodeX[10] = -1;
    nodeX[11] = 0;

    weight[4] = 0.837758041;
    nodeX[12] = 0;
    nodeX[13] = 0;
    nodeX[14] = 1;

    weight[5] = 0.837758041;
    nodeX[15] = 0;
    nodeX[16] = 0;
    nodeX[17] = -1;

    weight[6] = 0.9424777961;
    nodeX[18] = 0.5773502692;
    nodeX[19] = 0.5773502692;
    nodeX[20] = 0.5773502692;

    weight[7] = 0.9424777961;
    nodeX[21] = -0.5773502692;
    nodeX[22] = 0.5773502692;
    nodeX[23] = 0.5773502692;

    weight[8] = 0.9424777961;
    nodeX[24] = 0.5773502692;
    nodeX[25] = -0.5773502692;
    nodeX[26] = 0.5773502692;

    weight[9] = 0.9424777961;
    nodeX[27] = 0.5773502692;
    nodeX[28] = 0.5773502692;
    nodeX[29] = -0.5773502692;

    weight[10] = 0.9424777961;
    nodeX[30] = -0.5773502692;
    nodeX[31] = -0.5773502692;
    nodeX[32] = 0.5773502692;

    weight[11] = 0.9424777961;
    nodeX[33] = 0.5773502692;
    nodeX[34] = -0.5773502692;
    nodeX[35] = -0.5773502692;

    weight[12] = 0.9424777961;
    nodeX[36] = -0.5773502692;
    nodeX[37] = 0.5773502692;
    nodeX[38] = -0.5773502692;

    weight[13] = 0.9424777961;
    nodeX[39] = -0.5773502692;
    nodeX[40] = -0.5773502692;
    nodeX[41] = -0.5773502692;

    break;

    /*---------------------------------------------------*/
  case 110:
    numQuadNodes_table = 110;

    weight[0] = 0.04810746585;
    nodeX[0] = 1;
    nodeX[1] = 0;
    nodeX[2] = 0;

    weight[1] = 0.04810746585;
    nodeX[3] = -1;
    nodeX[4] = 0;
    nodeX[5] = 0;

    weight[2] = 0.04810746585;
    nodeX[6] = 0;
    nodeX[7] = 1;
    nodeX[8] = 0;

    weight[3] = 0.04810746585;
    nodeX[9] = 0;
    nodeX[10] = -1;
    nodeX[11] = 0;

    weight[4] = 0.04810746585;
    nodeX[12] = 0;
    nodeX[13] = 0;
    nodeX[14] = 1;

    weight[5] = 0.04810746585;
    nodeX[15] = 0;
    nodeX[16] = 0;
    nodeX[17] = -1;

    weight[6] = 0.1230717353;
    nodeX[18] = 0.5773502692;
    nodeX[19] = 0.5773502692;
    nodeX[20] = 0.5773502692;

    weight[7] = 0.1230717353;
    nodeX[21] = -0.5773502692;
    nodeX[22] = 0.5773502692;
    nodeX[23] = 0.5773502692;

    weight[8] = 0.1230717353;
    nodeX[24] = 0.5773502692;
    nodeX[25] = -0.5773502692;
    nodeX[26] = 0.5773502692;

    weight[9] = 0.1230717353;
    nodeX[27] = 0.5773502692;
    nodeX[28] = 0.5773502692;
    nodeX[29] = -0.5773502692;

    weight[10] = 0.1230717353;
    nodeX[30] = -0.5773502692;
    nodeX[31] = -0.5773502692;
    nodeX[32] = 0.5773502692;

    weight[11] = 0.1230717353;
    nodeX[33] = 0.5773502692;
    nodeX[34] = -0.5773502692;
    nodeX[35] = -0.5773502692;

    weight[12] = 0.1230717353;
    nodeX[36] = -0.5773502692;
    nodeX[37] = 0.5773502692;
    nodeX[38] = -0.5773502692;

    weight[13] = 0.1230717353;
    nodeX[39] = -0.5773502692;
    nodeX[40] = -0.5773502692;
    nodeX[41] = -0.5773502692;

    weight[14] = 0.1031917341;
    nodeX[42] = 0.1851156353;
    nodeX[43] = 0.1851156353;
    nodeX[44] = 0.9651240351;

    weight[15] = 0.1031917341;
    nodeX[45] = -0.1851156353;
    nodeX[46] = 0.1851156353;
    nodeX[47] = 0.9651240351;

    weight[16] = 0.1031917341;
    nodeX[48] = 0.1851156353;
    nodeX[49] = -0.1851156353;
    nodeX[50] = 0.9651240351;

    weight[17] = 0.1031917341;
    nodeX[51] = 0.1851156353;
    nodeX[52] = 0.1851156353;
    nodeX[53] = -0.9651240351;

    weight[18] = 0.1031917341;
    nodeX[54] = -0.1851156353;
    nodeX[55] = -0.1851156353;
    nodeX[56] = 0.9651240351;

    weight[19] = 0.1031917341;
    nodeX[57] = -0.1851156353;
    nodeX[58] = 0.1851156353;
    nodeX[59] = -0.9651240351;

    weight[20] = 0.1031917341;
    nodeX[60] = 0.1851156353;
    nodeX[61] = -0.1851156353;
    nodeX[62] = -0.9651240351;

    weight[21] = 0.1031917341;
    nodeX[63] = -0.1851156353;
    nodeX[64] = -0.1851156353;
    nodeX[65] = -0.9651240351;

    weight[22] = 0.1031917341;
    nodeX[66] = -0.1851156353;
    nodeX[67] = 0.9651240351;
    nodeX[68] = 0.1851156353;

    weight[23] = 0.1031917341;
    nodeX[69] = 0.1851156353;
    nodeX[70] = -0.9651240351;
    nodeX[71] = 0.1851156353;

    weight[24] = 0.1031917341;
    nodeX[72] = 0.1851156353;
    nodeX[73] = 0.9651240351;
    nodeX[74] = -0.1851156353;

    weight[25] = 0.1031917341;
    nodeX[75] = -0.1851156353;
    nodeX[76] = -0.9651240351;
    nodeX[77] = 0.1851156353;

    weight[26] = 0.1031917341;
    nodeX[78] = -0.1851156353;
    nodeX[79] = 0.9651240351;
    nodeX[80] = -0.1851156353;

    weight[27] = 0.1031917341;
    nodeX[81] = 0.1851156353;
    nodeX[82] = -0.9651240351;
    nodeX[83] = -0.1851156353;

    weight[28] = 0.1031917341;
    nodeX[84] = -0.1851156353;
    nodeX[85] = -0.9651240351;
    nodeX[86] = -0.1851156353;

    weight[29] = 0.1031917341;
    nodeX[87] = 0.1851156353;
    nodeX[88] = 0.9651240351;
    nodeX[89] = 0.1851156353;

    weight[30] = 0.1031917341;
    nodeX[90] = 0.9651240351;
    nodeX[91] = 0.1851156353;
    nodeX[92] = 0.1851156353;

    weight[31] = 0.1031917341;
    nodeX[93] = -0.9651240351;
    nodeX[94] = 0.1851156353;
    nodeX[95] = 0.1851156353;

    weight[32] = 0.1031917341;
    nodeX[96] = 0.9651240351;
    nodeX[97] = -0.1851156353;
    nodeX[98] = 0.1851156353;

    weight[33] = 0.1031917341;
    nodeX[99] = 0.9651240351;
    nodeX[100] = 0.1851156353;
    nodeX[101] = -0.1851156353;

    weight[34] = 0.1031917341;
    nodeX[102] = -0.9651240351;
    nodeX[103] = -0.1851156353;
    nodeX[104] = 0.1851156353;

    weight[35] = 0.1031917341;
    nodeX[105] = -0.9651240351;
    nodeX[106] = 0.1851156353;
    nodeX[107] = -0.1851156353;

    weight[36] = 0.1031917341;
    nodeX[108] = 0.9651240351;
    nodeX[109] = -0.1851156353;
    nodeX[110] = -0.1851156353;

    weight[37] = 0.1031917341;
    nodeX[111] = -0.9651240351;
    nodeX[112] = -0.1851156353;
    nodeX[113] = -0.1851156353;

    weight[38] = 0.1249450969;
    nodeX[114] = 0.6904210484;
    nodeX[115] = 0.6904210484;
    nodeX[116] = 0.2159572918;

    weight[39] = 0.1249450969;
    nodeX[117] = -0.6904210484;
    nodeX[118] = 0.6904210484;
    nodeX[119] = 0.2159572918;

    weight[40] = 0.1249450969;
    nodeX[120] = 0.6904210484;
    nodeX[121] = -0.6904210484;
    nodeX[122] = 0.2159572918;

    weight[41] = 0.1249450969;
    nodeX[123] = 0.6904210484;
    nodeX[124] = 0.6904210484;
    nodeX[125] = -0.2159572918;

    weight[42] = 0.1249450969;
    nodeX[126] = -0.6904210484;
    nodeX[127] = -0.6904210484;
    nodeX[128] = 0.2159572918;

    weight[43] = 0.1249450969;
    nodeX[129] = -0.6904210484;
    nodeX[130] = 0.6904210484;
    nodeX[131] = -0.2159572918;

    weight[44] = 0.1249450969;
    nodeX[132] = 0.6904210484;
    nodeX[133] = -0.6904210484;
    nodeX[134] = -0.2159572918;

    weight[45] = 0.1249450969;
    nodeX[135] = -0.6904210484;
    nodeX[136] = -0.6904210484;
    nodeX[137] = -0.2159572918;

    weight[46] = 0.1249450969;
    nodeX[138] = -0.6904210484;
    nodeX[139] = 0.2159572918;
    nodeX[140] = 0.6904210484;

    weight[47] = 0.1249450969;
    nodeX[141] = 0.6904210484;
    nodeX[142] = -0.2159572918;
    nodeX[143] = 0.6904210484;

    weight[48] = 0.1249450969;
    nodeX[144] = 0.6904210484;
    nodeX[145] = 0.2159572918;
    nodeX[146] = -0.6904210484;

    weight[49] = 0.1249450969;
    nodeX[147] = -0.6904210484;
    nodeX[148] = -0.2159572918;
    nodeX[149] = 0.6904210484;

    weight[50] = 0.1249450969;
    nodeX[150] = -0.6904210484;
    nodeX[151] = 0.2159572918;
    nodeX[152] = -0.6904210484;

    weight[51] = 0.1249450969;
    nodeX[153] = 0.6904210484;
    nodeX[154] = -0.2159572918;
    nodeX[155] = -0.6904210484;

    weight[52] = 0.1249450969;
    nodeX[156] = -0.6904210484;
    nodeX[157] = -0.2159572918;
    nodeX[158] = -0.6904210484;

    weight[53] = 0.1249450969;
    nodeX[159] = 0.6904210484;
    nodeX[160] = 0.2159572918;
    nodeX[161] = 0.6904210484;

    weight[54] = 0.1249450969;
    nodeX[162] = 0.2159572918;
    nodeX[163] = 0.6904210484;
    nodeX[164] = 0.6904210484;

    weight[55] = 0.1249450969;
    nodeX[165] = -0.2159572918;
    nodeX[166] = 0.6904210484;
    nodeX[167] = 0.6904210484;

    weight[56] = 0.1249450969;
    nodeX[168] = 0.2159572918;
    nodeX[169] = -0.6904210484;
    nodeX[170] = 0.6904210484;

    weight[57] = 0.1249450969;
    nodeX[171] = 0.2159572918;
    nodeX[172] = 0.6904210484;
    nodeX[173] = -0.6904210484;

    weight[58] = 0.1249450969;
    nodeX[174] = -0.2159572918;
    nodeX[175] = -0.6904210484;
    nodeX[176] = 0.6904210484;

    weight[59] = 0.1249450969;
    nodeX[177] = -0.2159572918;
    nodeX[178] = 0.6904210484;
    nodeX[179] = -0.6904210484;

    weight[60] = 0.1249450969;
    nodeX[180] = 0.2159572918;
    nodeX[181] = -0.6904210484;
    nodeX[182] = -0.6904210484;

    weight[61] = 0.1249450969;
    nodeX[183] = -0.2159572918;
    nodeX[184] = -0.6904210484;
    nodeX[185] = -0.6904210484;

    weight[62] = 0.120580249;
    nodeX[186] = 0.3956894731;
    nodeX[187] = 0.3956894731;
    nodeX[188] = 0.8287699813;

    weight[63] = 0.120580249;
    nodeX[189] = -0.3956894731;
    nodeX[190] = 0.3956894731;
    nodeX[191] = 0.8287699813;

    weight[64] = 0.120580249;
    nodeX[192] = 0.3956894731;
    nodeX[193] = -0.3956894731;
    nodeX[194] = 0.8287699813;

    weight[65] = 0.120580249;
    nodeX[195] = 0.3956894731;
    nodeX[196] = 0.3956894731;
    nodeX[197] = -0.8287699813;

    weight[66] = 0.120580249;
    nodeX[198] = -0.3956894731;
    nodeX[199] = -0.3956894731;
    nodeX[200] = 0.8287699813;

    weight[67] = 0.120580249;
    nodeX[201] = -0.3956894731;
    nodeX[202] = 0.3956894731;
    nodeX[203] = -0.8287699813;

    weight[68] = 0.120580249;
    nodeX[204] = 0.3956894731;
    nodeX[205] = -0.3956894731;
    nodeX[206] = -0.8287699813;

    weight[69] = 0.120580249;
    nodeX[207] = -0.3956894731;
    nodeX[208] = -0.3956894731;
    nodeX[209] = -0.8287699813;

    weight[70] = 0.120580249;
    nodeX[210] = -0.3956894731;
    nodeX[211] = 0.8287699813;
    nodeX[212] = 0.3956894731;

    weight[71] = 0.120580249;
    nodeX[213] = 0.3956894731;
    nodeX[214] = -0.8287699813;
    nodeX[215] = 0.3956894731;

    weight[72] = 0.120580249;
    nodeX[216] = 0.3956894731;
    nodeX[217] = 0.8287699813;
    nodeX[218] = -0.3956894731;

    weight[73] = 0.120580249;
    nodeX[219] = -0.3956894731;
    nodeX[220] = -0.8287699813;
    nodeX[221] = 0.3956894731;

    weight[74] = 0.120580249;
    nodeX[222] = -0.3956894731;
    nodeX[223] = 0.8287699813;
    nodeX[224] = -0.3956894731;

    weight[75] = 0.120580249;
    nodeX[225] = 0.3956894731;
    nodeX[226] = -0.8287699813;
    nodeX[227] = -0.3956894731;

    weight[76] = 0.120580249;
    nodeX[228] = -0.3956894731;
    nodeX[229] = -0.8287699813;
    nodeX[230] = -0.3956894731;

    weight[77] = 0.120580249;
    nodeX[231] = 0.3956894731;
    nodeX[232] = 0.8287699813;
    nodeX[233] = 0.3956894731;

    weight[78] = 0.120580249;
    nodeX[234] = 0.8287699813;
    nodeX[235] = 0.3956894731;
    nodeX[236] = 0.3956894731;

    weight[79] = 0.120580249;
    nodeX[237] = -0.8287699813;
    nodeX[238] = 0.3956894731;
    nodeX[239] = 0.3956894731;

    weight[80] = 0.120580249;
    nodeX[240] = 0.8287699813;
    nodeX[241] = -0.3956894731;
    nodeX[242] = 0.3956894731;

    weight[81] = 0.120580249;
    nodeX[243] = 0.8287699813;
    nodeX[244] = 0.3956894731;
    nodeX[245] = -0.3956894731;

    weight[82] = 0.120580249;
    nodeX[246] = -0.8287699813;
    nodeX[247] = -0.3956894731;
    nodeX[248] = 0.3956894731;

    weight[83] = 0.120580249;
    nodeX[249] = -0.8287699813;
    nodeX[250] = 0.3956894731;
    nodeX[251] = -0.3956894731;

    weight[84] = 0.120580249;
    nodeX[252] = 0.8287699813;
    nodeX[253] = -0.3956894731;
    nodeX[254] = -0.3956894731;

    weight[85] = 0.120580249;
    nodeX[255] = -0.8287699813;
    nodeX[256] = -0.3956894731;
    nodeX[257] = -0.3956894731;

    weight[86] = 0.1218309174;
    nodeX[258] = 0.4783690288;
    nodeX[259] = 0.8781589106;
    nodeX[260] = 0;

    weight[87] = 0.1218309174;
    nodeX[261] = -0.4783690288;
    nodeX[262] = 0.8781589106;
    nodeX[263] = 0;

    weight[88] = 0.1218309174;
    nodeX[264] = 0.4783690288;
    nodeX[265] = -0.8781589106;
    nodeX[266] = 0;

    weight[89] = 0.1218309174;
    nodeX[267] = -0.4783690288;
    nodeX[268] = -0.8781589106;
    nodeX[269] = 0;

    weight[90] = 0.1218309174;
    nodeX[270] = 0.8781589106;
    nodeX[271] = 0.4783690288;
    nodeX[272] = 0;

    weight[91] = 0.1218309174;
    nodeX[273] = -0.8781589106;
    nodeX[274] = 0.4783690288;
    nodeX[275] = 0;

    weight[92] = 0.1218309174;
    nodeX[276] = 0.8781589106;
    nodeX[277] = -0.4783690288;
    nodeX[278] = 0;

    weight[93] = 0.1218309174;
    nodeX[279] = -0.8781589106;
    nodeX[280] = -0.4783690288;
    nodeX[281] = 0;

    weight[94] = 0.1218309174;
    nodeX[282] = 0.4783690288;
    nodeX[283] = 0;
    nodeX[284] = 0.8781589106;

    weight[95] = 0.1218309174;
    nodeX[285] = -0.4783690288;
    nodeX[286] = 0;
    nodeX[287] = 0.8781589106;

    weight[96] = 0.1218309174;
    nodeX[288] = 0.4783690288;
    nodeX[289] = 0;
    nodeX[290] = -0.8781589106;

    weight[97] = 0.1218309174;
    nodeX[291] = -0.4783690288;
    nodeX[292] = 0;
    nodeX[293] = -0.8781589106;

    weight[98] = 0.1218309174;
    nodeX[294] = 0.8781589106;
    nodeX[295] = 0;
    nodeX[296] = 0.4783690288;

    weight[99] = 0.1218309174;
    nodeX[297] = -0.8781589106;
    nodeX[298] = 0;
    nodeX[299] = 0.4783690288;

    weight[100] = 0.1218309174;
    nodeX[300] = 0.8781589106;
    nodeX[301] = 0;
    nodeX[302] = -0.4783690288;

    weight[101] = 0.1218309174;
    nodeX[303] = -0.8781589106;
    nodeX[304] = 0;
    nodeX[305] = -0.4783690288;

    weight[102] = 0.1218309174;
    nodeX[306] = 0;
    nodeX[307] = 0.4783690288;
    nodeX[308] = 0.8781589106;

    weight[103] = 0.1218309174;
    nodeX[309] = 0;
    nodeX[310] = -0.4783690288;
    nodeX[311] = 0.8781589106;

    weight[104] = 0.1218309174;
    nodeX[312] = 0;
    nodeX[313] = 0.4783690288;
    nodeX[314] = -0.8781589106;

    weight[105] = 0.1218309174;
    nodeX[315] = 0;
    nodeX[316] = -0.4783690288;
    nodeX[317] = -0.8781589106;

    weight[106] = 0.1218309174;
    nodeX[318] = 0;
    nodeX[319] = 0.8781589106;
    nodeX[320] = 0.4783690288;

    weight[107] = 0.1218309174;
    nodeX[321] = 0;
    nodeX[322] = -0.8781589106;
    nodeX[323] = 0.4783690288;

    weight[108] = 0.1218309174;
    nodeX[324] = 0;
    nodeX[325] = 0.8781589106;
    nodeX[326] = -0.4783690288;

    weight[109] = 0.1218309174;
    nodeX[327] = 0;
    nodeX[328] = -0.8781589106;
    nodeX[329] = -0.4783690288;

    break;

    /*---------------------------------------------------*/
  case 1202:
    numQuadNodes_table = 1202;

    weight[0] = 0.00138882175;
    nodeX[0] = 1;
    nodeX[1] = 0;
    nodeX[2] = 0;

    weight[1] = 0.00138882175;
    nodeX[3] = -1;
    nodeX[4] = 0;
    nodeX[5] = 0;

    weight[2] = 0.00138882175;
    nodeX[6] = 0;
    nodeX[7] = 1;
    nodeX[8] = 0;

    weight[3] = 0.00138882175;
    nodeX[9] = 0;
    nodeX[10] = -1;
    nodeX[11] = 0;

    weight[4] = 0.00138882175;
    nodeX[12] = 0;
    nodeX[13] = 0;
    nodeX[14] = 1;

    weight[5] = 0.00138882175;
    nodeX[15] = 0;
    nodeX[16] = 0;
    nodeX[17] = -1;

    weight[6] = 0.01156763662;
    nodeX[18] = 0;
    nodeX[19] = 0.7071067812;
    nodeX[20] = 0.7071067812;

    weight[7] = 0.01156763662;
    nodeX[21] = 0;
    nodeX[22] = -0.7071067812;
    nodeX[23] = 0.7071067812;

    weight[8] = 0.01156763662;
    nodeX[24] = 0;
    nodeX[25] = 0.7071067812;
    nodeX[26] = -0.7071067812;

    weight[9] = 0.01156763662;
    nodeX[27] = 0;
    nodeX[28] = -0.7071067812;
    nodeX[29] = -0.7071067812;

    weight[10] = 0.01156763662;
    nodeX[30] = 0.7071067812;
    nodeX[31] = 0;
    nodeX[32] = 0.7071067812;

    weight[11] = 0.01156763662;
    nodeX[33] = 0.7071067812;
    nodeX[34] = 0;
    nodeX[35] = -0.7071067812;

    weight[12] = 0.01156763662;
    nodeX[36] = -0.7071067812;
    nodeX[37] = 0;
    nodeX[38] = 0.7071067812;

    weight[13] = 0.01156763662;
    nodeX[39] = -0.7071067812;
    nodeX[40] = 0;
    nodeX[41] = -0.7071067812;

    weight[14] = 0.01156763662;
    nodeX[42] = 0.7071067812;
    nodeX[43] = 0.7071067812;
    nodeX[44] = 0;

    weight[15] = 0.01156763662;
    nodeX[45] = -0.7071067812;
    nodeX[46] = 0.7071067812;
    nodeX[47] = 0;

    weight[16] = 0.01156763662;
    nodeX[48] = 0.7071067812;
    nodeX[49] = -0.7071067812;
    nodeX[50] = 0;

    weight[17] = 0.01156763662;
    nodeX[51] = -0.7071067812;
    nodeX[52] = -0.7071067812;
    nodeX[53] = 0;

    weight[18] = 0.01147706708;
    nodeX[54] = 0.5773502692;
    nodeX[55] = 0.5773502692;
    nodeX[56] = 0.5773502692;

    weight[19] = 0.01147706708;
    nodeX[57] = -0.5773502692;
    nodeX[58] = 0.5773502692;
    nodeX[59] = 0.5773502692;

    weight[20] = 0.01147706708;
    nodeX[60] = 0.5773502692;
    nodeX[61] = -0.5773502692;
    nodeX[62] = 0.5773502692;

    weight[21] = 0.01147706708;
    nodeX[63] = 0.5773502692;
    nodeX[64] = 0.5773502692;
    nodeX[65] = -0.5773502692;

    weight[22] = 0.01147706708;
    nodeX[66] = -0.5773502692;
    nodeX[67] = -0.5773502692;
    nodeX[68] = 0.5773502692;

    weight[23] = 0.01147706708;
    nodeX[69] = 0.5773502692;
    nodeX[70] = -0.5773502692;
    nodeX[71] = -0.5773502692;

    weight[24] = 0.01147706708;
    nodeX[72] = -0.5773502692;
    nodeX[73] = 0.5773502692;
    nodeX[74] = -0.5773502692;

    weight[25] = 0.01147706708;
    nodeX[75] = -0.5773502692;
    nodeX[76] = -0.5773502692;
    nodeX[77] = -0.5773502692;

    weight[26] = 0.004637520929;
    nodeX[78] = 0.0371263645;
    nodeX[79] = 0.0371263645;
    nodeX[80] = 0.9986206818;

    weight[27] = 0.004637520929;
    nodeX[81] = -0.0371263645;
    nodeX[82] = 0.0371263645;
    nodeX[83] = 0.9986206818;

    weight[28] = 0.004637520929;
    nodeX[84] = 0.0371263645;
    nodeX[85] = -0.0371263645;
    nodeX[86] = 0.9986206818;

    weight[29] = 0.004637520929;
    nodeX[87] = 0.0371263645;
    nodeX[88] = 0.0371263645;
    nodeX[89] = -0.9986206818;

    weight[30] = 0.004637520929;
    nodeX[90] = -0.0371263645;
    nodeX[91] = -0.0371263645;
    nodeX[92] = 0.9986206818;

    weight[31] = 0.004637520929;
    nodeX[93] = -0.0371263645;
    nodeX[94] = 0.0371263645;
    nodeX[95] = -0.9986206818;

    weight[32] = 0.004637520929;
    nodeX[96] = 0.0371263645;
    nodeX[97] = -0.0371263645;
    nodeX[98] = -0.9986206818;

    weight[33] = 0.004637520929;
    nodeX[99] = -0.0371263645;
    nodeX[100] = -0.0371263645;
    nodeX[101] = -0.9986206818;

    weight[34] = 0.004637520929;
    nodeX[102] = -0.0371263645;
    nodeX[103] = 0.9986206818;
    nodeX[104] = 0.0371263645;

    weight[35] = 0.004637520929;
    nodeX[105] = 0.0371263645;
    nodeX[106] = -0.9986206818;
    nodeX[107] = 0.0371263645;

    weight[36] = 0.004637520929;
    nodeX[108] = 0.0371263645;
    nodeX[109] = 0.9986206818;
    nodeX[110] = -0.0371263645;

    weight[37] = 0.004637520929;
    nodeX[111] = -0.0371263645;
    nodeX[112] = -0.9986206818;
    nodeX[113] = 0.0371263645;

    weight[38] = 0.004637520929;
    nodeX[114] = -0.0371263645;
    nodeX[115] = 0.9986206818;
    nodeX[116] = -0.0371263645;

    weight[39] = 0.004637520929;
    nodeX[117] = 0.0371263645;
    nodeX[118] = -0.9986206818;
    nodeX[119] = -0.0371263645;

    weight[40] = 0.004637520929;
    nodeX[120] = -0.0371263645;
    nodeX[121] = -0.9986206818;
    nodeX[122] = -0.0371263645;

    weight[41] = 0.004637520929;
    nodeX[123] = 0.0371263645;
    nodeX[124] = 0.9986206818;
    nodeX[125] = 0.0371263645;

    weight[42] = 0.004637520929;
    nodeX[126] = 0.9986206818;
    nodeX[127] = 0.0371263645;
    nodeX[128] = 0.0371263645;

    weight[43] = 0.004637520929;
    nodeX[129] = -0.9986206818;
    nodeX[130] = 0.0371263645;
    nodeX[131] = 0.0371263645;

    weight[44] = 0.004637520929;
    nodeX[132] = 0.9986206818;
    nodeX[133] = -0.0371263645;
    nodeX[134] = 0.0371263645;

    weight[45] = 0.004637520929;
    nodeX[135] = 0.9986206818;
    nodeX[136] = 0.0371263645;
    nodeX[137] = -0.0371263645;

    weight[46] = 0.004637520929;
    nodeX[138] = -0.9986206818;
    nodeX[139] = -0.0371263645;
    nodeX[140] = 0.0371263645;

    weight[47] = 0.004637520929;
    nodeX[141] = -0.9986206818;
    nodeX[142] = 0.0371263645;
    nodeX[143] = -0.0371263645;

    weight[48] = 0.004637520929;
    nodeX[144] = 0.9986206818;
    nodeX[145] = -0.0371263645;
    nodeX[146] = -0.0371263645;

    weight[49] = 0.004637520929;
    nodeX[147] = -0.9986206818;
    nodeX[148] = -0.0371263645;
    nodeX[149] = -0.0371263645;

    weight[50] = 0.007042182693;
    nodeX[150] = 0.09140060412;
    nodeX[151] = 0.09140060412;
    nodeX[152] = 0.9916107397;

    weight[51] = 0.007042182693;
    nodeX[153] = -0.09140060412;
    nodeX[154] = 0.09140060412;
    nodeX[155] = 0.9916107397;

    weight[52] = 0.007042182693;
    nodeX[156] = 0.09140060412;
    nodeX[157] = -0.09140060412;
    nodeX[158] = 0.9916107397;

    weight[53] = 0.007042182693;
    nodeX[159] = 0.09140060412;
    nodeX[160] = 0.09140060412;
    nodeX[161] = -0.9916107397;

    weight[54] = 0.007042182693;
    nodeX[162] = -0.09140060412;
    nodeX[163] = -0.09140060412;
    nodeX[164] = 0.9916107397;

    weight[55] = 0.007042182693;
    nodeX[165] = -0.09140060412;
    nodeX[166] = 0.09140060412;
    nodeX[167] = -0.9916107397;

    weight[56] = 0.007042182693;
    nodeX[168] = 0.09140060412;
    nodeX[169] = -0.09140060412;
    nodeX[170] = -0.9916107397;

    weight[57] = 0.007042182693;
    nodeX[171] = -0.09140060412;
    nodeX[172] = -0.09140060412;
    nodeX[173] = -0.9916107397;

    weight[58] = 0.007042182693;
    nodeX[174] = -0.09140060412;
    nodeX[175] = 0.9916107397;
    nodeX[176] = 0.09140060412;

    weight[59] = 0.007042182693;
    nodeX[177] = 0.09140060412;
    nodeX[178] = -0.9916107397;
    nodeX[179] = 0.09140060412;

    weight[60] = 0.007042182693;
    nodeX[180] = 0.09140060412;
    nodeX[181] = 0.9916107397;
    nodeX[182] = -0.09140060412;

    weight[61] = 0.007042182693;
    nodeX[183] = -0.09140060412;
    nodeX[184] = -0.9916107397;
    nodeX[185] = 0.09140060412;

    weight[62] = 0.007042182693;
    nodeX[186] = -0.09140060412;
    nodeX[187] = 0.9916107397;
    nodeX[188] = -0.09140060412;

    weight[63] = 0.007042182693;
    nodeX[189] = 0.09140060412;
    nodeX[190] = -0.9916107397;
    nodeX[191] = -0.09140060412;

    weight[64] = 0.007042182693;
    nodeX[192] = -0.09140060412;
    nodeX[193] = -0.9916107397;
    nodeX[194] = -0.09140060412;

    weight[65] = 0.007042182693;
    nodeX[195] = 0.09140060412;
    nodeX[196] = 0.9916107397;
    nodeX[197] = 0.09140060412;

    weight[66] = 0.007042182693;
    nodeX[198] = 0.9916107397;
    nodeX[199] = 0.09140060412;
    nodeX[200] = 0.09140060412;

    weight[67] = 0.007042182693;
    nodeX[201] = -0.9916107397;
    nodeX[202] = 0.09140060412;
    nodeX[203] = 0.09140060412;

    weight[68] = 0.007042182693;
    nodeX[204] = 0.9916107397;
    nodeX[205] = -0.09140060412;
    nodeX[206] = 0.09140060412;

    weight[69] = 0.007042182693;
    nodeX[207] = 0.9916107397;
    nodeX[208] = 0.09140060412;
    nodeX[209] = -0.09140060412;

    weight[70] = 0.007042182693;
    nodeX[210] = -0.9916107397;
    nodeX[211] = -0.09140060412;
    nodeX[212] = 0.09140060412;

    weight[71] = 0.007042182693;
    nodeX[213] = -0.9916107397;
    nodeX[214] = 0.09140060412;
    nodeX[215] = -0.09140060412;

    weight[72] = 0.007042182693;
    nodeX[216] = 0.9916107397;
    nodeX[217] = -0.09140060412;
    nodeX[218] = -0.09140060412;

    weight[73] = 0.007042182693;
    nodeX[219] = -0.9916107397;
    nodeX[220] = -0.09140060412;
    nodeX[221] = -0.09140060412;

    weight[74] = 0.008627187439;
    nodeX[222] = 0.1531077852;
    nodeX[223] = 0.1531077852;
    nodeX[224] = 0.9762766064;

    weight[75] = 0.008627187439;
    nodeX[225] = -0.1531077852;
    nodeX[226] = 0.1531077852;
    nodeX[227] = 0.9762766064;

    weight[76] = 0.008627187439;
    nodeX[228] = 0.1531077852;
    nodeX[229] = -0.1531077852;
    nodeX[230] = 0.9762766064;

    weight[77] = 0.008627187439;
    nodeX[231] = 0.1531077852;
    nodeX[232] = 0.1531077852;
    nodeX[233] = -0.9762766064;

    weight[78] = 0.008627187439;
    nodeX[234] = -0.1531077852;
    nodeX[235] = -0.1531077852;
    nodeX[236] = 0.9762766064;

    weight[79] = 0.008627187439;
    nodeX[237] = -0.1531077852;
    nodeX[238] = 0.1531077852;
    nodeX[239] = -0.9762766064;

    weight[80] = 0.008627187439;
    nodeX[240] = 0.1531077852;
    nodeX[241] = -0.1531077852;
    nodeX[242] = -0.9762766064;

    weight[81] = 0.008627187439;
    nodeX[243] = -0.1531077852;
    nodeX[244] = -0.1531077852;
    nodeX[245] = -0.9762766064;

    weight[82] = 0.008627187439;
    nodeX[246] = -0.1531077852;
    nodeX[247] = 0.9762766064;
    nodeX[248] = 0.1531077852;

    weight[83] = 0.008627187439;
    nodeX[249] = 0.1531077852;
    nodeX[250] = -0.9762766064;
    nodeX[251] = 0.1531077852;

    weight[84] = 0.008627187439;
    nodeX[252] = 0.1531077852;
    nodeX[253] = 0.9762766064;
    nodeX[254] = -0.1531077852;

    weight[85] = 0.008627187439;
    nodeX[255] = -0.1531077852;
    nodeX[256] = -0.9762766064;
    nodeX[257] = 0.1531077852;

    weight[86] = 0.008627187439;
    nodeX[258] = -0.1531077852;
    nodeX[259] = 0.9762766064;
    nodeX[260] = -0.1531077852;

    weight[87] = 0.008627187439;
    nodeX[261] = 0.1531077852;
    nodeX[262] = -0.9762766064;
    nodeX[263] = -0.1531077852;

    weight[88] = 0.008627187439;
    nodeX[264] = -0.1531077852;
    nodeX[265] = -0.9762766064;
    nodeX[266] = -0.1531077852;

    weight[89] = 0.008627187439;
    nodeX[267] = 0.1531077852;
    nodeX[268] = 0.9762766064;
    nodeX[269] = 0.1531077852;

    weight[90] = 0.008627187439;
    nodeX[270] = 0.9762766064;
    nodeX[271] = 0.1531077852;
    nodeX[272] = 0.1531077852;

    weight[91] = 0.008627187439;
    nodeX[273] = -0.9762766064;
    nodeX[274] = 0.1531077852;
    nodeX[275] = 0.1531077852;

    weight[92] = 0.008627187439;
    nodeX[276] = 0.9762766064;
    nodeX[277] = -0.1531077852;
    nodeX[278] = 0.1531077852;

    weight[93] = 0.008627187439;
    nodeX[279] = 0.9762766064;
    nodeX[280] = 0.1531077852;
    nodeX[281] = -0.1531077852;

    weight[94] = 0.008627187439;
    nodeX[282] = -0.9762766064;
    nodeX[283] = -0.1531077852;
    nodeX[284] = 0.1531077852;

    weight[95] = 0.008627187439;
    nodeX[285] = -0.9762766064;
    nodeX[286] = 0.1531077852;
    nodeX[287] = -0.1531077852;

    weight[96] = 0.008627187439;
    nodeX[288] = 0.9762766064;
    nodeX[289] = -0.1531077852;
    nodeX[290] = -0.1531077852;

    weight[97] = 0.008627187439;
    nodeX[291] = -0.9762766064;
    nodeX[292] = -0.1531077852;
    nodeX[293] = -0.1531077852;

    weight[98] = 0.00970166355;
    nodeX[294] = 0.2180928892;
    nodeX[295] = 0.2180928892;
    nodeX[296] = 0.9512470675;

    weight[99] = 0.00970166355;
    nodeX[297] = -0.2180928892;
    nodeX[298] = 0.2180928892;
    nodeX[299] = 0.9512470675;

    weight[100] = 0.00970166355;
    nodeX[300] = 0.2180928892;
    nodeX[301] = -0.2180928892;
    nodeX[302] = 0.9512470675;

    weight[101] = 0.00970166355;
    nodeX[303] = 0.2180928892;
    nodeX[304] = 0.2180928892;
    nodeX[305] = -0.9512470675;

    weight[102] = 0.00970166355;
    nodeX[306] = -0.2180928892;
    nodeX[307] = -0.2180928892;
    nodeX[308] = 0.9512470675;

    weight[103] = 0.00970166355;
    nodeX[309] = -0.2180928892;
    nodeX[310] = 0.2180928892;
    nodeX[311] = -0.9512470675;

    weight[104] = 0.00970166355;
    nodeX[312] = 0.2180928892;
    nodeX[313] = -0.2180928892;
    nodeX[314] = -0.9512470675;

    weight[105] = 0.00970166355;
    nodeX[315] = -0.2180928892;
    nodeX[316] = -0.2180928892;
    nodeX[317] = -0.9512470675;

    weight[106] = 0.00970166355;
    nodeX[318] = -0.2180928892;
    nodeX[319] = 0.9512470675;
    nodeX[320] = 0.2180928892;

    weight[107] = 0.00970166355;
    nodeX[321] = 0.2180928892;
    nodeX[322] = -0.9512470675;
    nodeX[323] = 0.2180928892;

    weight[108] = 0.00970166355;
    nodeX[324] = 0.2180928892;
    nodeX[325] = 0.9512470675;
    nodeX[326] = -0.2180928892;

    weight[109] = 0.00970166355;
    nodeX[327] = -0.2180928892;
    nodeX[328] = -0.9512470675;
    nodeX[329] = 0.2180928892;

    weight[110] = 0.00970166355;
    nodeX[330] = -0.2180928892;
    nodeX[331] = 0.9512470675;
    nodeX[332] = -0.2180928892;

    weight[111] = 0.00970166355;
    nodeX[333] = 0.2180928892;
    nodeX[334] = -0.9512470675;
    nodeX[335] = -0.2180928892;

    weight[112] = 0.00970166355;
    nodeX[336] = -0.2180928892;
    nodeX[337] = -0.9512470675;
    nodeX[338] = -0.2180928892;

    weight[113] = 0.00970166355;
    nodeX[339] = 0.2180928892;
    nodeX[340] = 0.9512470675;
    nodeX[341] = 0.2180928892;

    weight[114] = 0.00970166355;
    nodeX[342] = 0.9512470675;
    nodeX[343] = 0.2180928892;
    nodeX[344] = 0.2180928892;

    weight[115] = 0.00970166355;
    nodeX[345] = -0.9512470675;
    nodeX[346] = 0.2180928892;
    nodeX[347] = 0.2180928892;

    weight[116] = 0.00970166355;
    nodeX[348] = 0.9512470675;
    nodeX[349] = -0.2180928892;
    nodeX[350] = 0.2180928892;

    weight[117] = 0.00970166355;
    nodeX[351] = 0.9512470675;
    nodeX[352] = 0.2180928892;
    nodeX[353] = -0.2180928892;

    weight[118] = 0.00970166355;
    nodeX[354] = -0.9512470675;
    nodeX[355] = -0.2180928892;
    nodeX[356] = 0.2180928892;

    weight[119] = 0.00970166355;
    nodeX[357] = -0.9512470675;
    nodeX[358] = 0.2180928892;
    nodeX[359] = -0.2180928892;

    weight[120] = 0.00970166355;
    nodeX[360] = 0.9512470675;
    nodeX[361] = -0.2180928892;
    nodeX[362] = -0.2180928892;

    weight[121] = 0.00970166355;
    nodeX[363] = -0.9512470675;
    nodeX[364] = -0.2180928892;
    nodeX[365] = -0.2180928892;

    weight[122] = 0.01043203032;
    nodeX[366] = 0.2839874532;
    nodeX[367] = 0.2839874532;
    nodeX[368] = 0.9158068862;

    weight[123] = 0.01043203032;
    nodeX[369] = -0.2839874532;
    nodeX[370] = 0.2839874532;
    nodeX[371] = 0.9158068862;

    weight[124] = 0.01043203032;
    nodeX[372] = 0.2839874532;
    nodeX[373] = -0.2839874532;
    nodeX[374] = 0.9158068862;

    weight[125] = 0.01043203032;
    nodeX[375] = 0.2839874532;
    nodeX[376] = 0.2839874532;
    nodeX[377] = -0.9158068862;

    weight[126] = 0.01043203032;
    nodeX[378] = -0.2839874532;
    nodeX[379] = -0.2839874532;
    nodeX[380] = 0.9158068862;

    weight[127] = 0.01043203032;
    nodeX[381] = -0.2839874532;
    nodeX[382] = 0.2839874532;
    nodeX[383] = -0.9158068862;

    weight[128] = 0.01043203032;
    nodeX[384] = 0.2839874532;
    nodeX[385] = -0.2839874532;
    nodeX[386] = -0.9158068862;

    weight[129] = 0.01043203032;
    nodeX[387] = -0.2839874532;
    nodeX[388] = -0.2839874532;
    nodeX[389] = -0.9158068862;

    weight[130] = 0.01043203032;
    nodeX[390] = -0.2839874532;
    nodeX[391] = 0.9158068862;
    nodeX[392] = 0.2839874532;

    weight[131] = 0.01043203032;
    nodeX[393] = 0.2839874532;
    nodeX[394] = -0.9158068862;
    nodeX[395] = 0.2839874532;

    weight[132] = 0.01043203032;
    nodeX[396] = 0.2839874532;
    nodeX[397] = 0.9158068862;
    nodeX[398] = -0.2839874532;

    weight[133] = 0.01043203032;
    nodeX[399] = -0.2839874532;
    nodeX[400] = -0.9158068862;
    nodeX[401] = 0.2839874532;

    weight[134] = 0.01043203032;
    nodeX[402] = -0.2839874532;
    nodeX[403] = 0.9158068862;
    nodeX[404] = -0.2839874532;

    weight[135] = 0.01043203032;
    nodeX[405] = 0.2839874532;
    nodeX[406] = -0.9158068862;
    nodeX[407] = -0.2839874532;

    weight[136] = 0.01043203032;
    nodeX[408] = -0.2839874532;
    nodeX[409] = -0.9158068862;
    nodeX[410] = -0.2839874532;

    weight[137] = 0.01043203032;
    nodeX[411] = 0.2839874532;
    nodeX[412] = 0.9158068862;
    nodeX[413] = 0.2839874532;

    weight[138] = 0.01043203032;
    nodeX[414] = 0.9158068862;
    nodeX[415] = 0.2839874532;
    nodeX[416] = 0.2839874532;

    weight[139] = 0.01043203032;
    nodeX[417] = -0.9158068862;
    nodeX[418] = 0.2839874532;
    nodeX[419] = 0.2839874532;

    weight[140] = 0.01043203032;
    nodeX[420] = 0.9158068862;
    nodeX[421] = -0.2839874532;
    nodeX[422] = 0.2839874532;

    weight[141] = 0.01043203032;
    nodeX[423] = 0.9158068862;
    nodeX[424] = 0.2839874532;
    nodeX[425] = -0.2839874532;

    weight[142] = 0.01043203032;
    nodeX[426] = -0.9158068862;
    nodeX[427] = -0.2839874532;
    nodeX[428] = 0.2839874532;

    weight[143] = 0.01043203032;
    nodeX[429] = -0.9158068862;
    nodeX[430] = 0.2839874532;
    nodeX[431] = -0.2839874532;

    weight[144] = 0.01043203032;
    nodeX[432] = 0.9158068862;
    nodeX[433] = -0.2839874532;
    nodeX[434] = -0.2839874532;

    weight[145] = 0.01043203032;
    nodeX[435] = -0.9158068862;
    nodeX[436] = -0.2839874532;
    nodeX[437] = -0.2839874532;

    weight[146] = 0.0109160198;
    nodeX[438] = 0.3491177601;
    nodeX[439] = 0.3491177601;
    nodeX[440] = 0.8696169152;

    weight[147] = 0.0109160198;
    nodeX[441] = -0.3491177601;
    nodeX[442] = 0.3491177601;
    nodeX[443] = 0.8696169152;

    weight[148] = 0.0109160198;
    nodeX[444] = 0.3491177601;
    nodeX[445] = -0.3491177601;
    nodeX[446] = 0.8696169152;

    weight[149] = 0.0109160198;
    nodeX[447] = 0.3491177601;
    nodeX[448] = 0.3491177601;
    nodeX[449] = -0.8696169152;

    weight[150] = 0.0109160198;
    nodeX[450] = -0.3491177601;
    nodeX[451] = -0.3491177601;
    nodeX[452] = 0.8696169152;

    weight[151] = 0.0109160198;
    nodeX[453] = -0.3491177601;
    nodeX[454] = 0.3491177601;
    nodeX[455] = -0.8696169152;

    weight[152] = 0.0109160198;
    nodeX[456] = 0.3491177601;
    nodeX[457] = -0.3491177601;
    nodeX[458] = -0.8696169152;

    weight[153] = 0.0109160198;
    nodeX[459] = -0.3491177601;
    nodeX[460] = -0.3491177601;
    nodeX[461] = -0.8696169152;

    weight[154] = 0.0109160198;
    nodeX[462] = -0.3491177601;
    nodeX[463] = 0.8696169152;
    nodeX[464] = 0.3491177601;

    weight[155] = 0.0109160198;
    nodeX[465] = 0.3491177601;
    nodeX[466] = -0.8696169152;
    nodeX[467] = 0.3491177601;

    weight[156] = 0.0109160198;
    nodeX[468] = 0.3491177601;
    nodeX[469] = 0.8696169152;
    nodeX[470] = -0.3491177601;

    weight[157] = 0.0109160198;
    nodeX[471] = -0.3491177601;
    nodeX[472] = -0.8696169152;
    nodeX[473] = 0.3491177601;

    weight[158] = 0.0109160198;
    nodeX[474] = -0.3491177601;
    nodeX[475] = 0.8696169152;
    nodeX[476] = -0.3491177601;

    weight[159] = 0.0109160198;
    nodeX[477] = 0.3491177601;
    nodeX[478] = -0.8696169152;
    nodeX[479] = -0.3491177601;

    weight[160] = 0.0109160198;
    nodeX[480] = -0.3491177601;
    nodeX[481] = -0.8696169152;
    nodeX[482] = -0.3491177601;

    weight[161] = 0.0109160198;
    nodeX[483] = 0.3491177601;
    nodeX[484] = 0.8696169152;
    nodeX[485] = 0.3491177601;

    weight[162] = 0.0109160198;
    nodeX[486] = 0.8696169152;
    nodeX[487] = 0.3491177601;
    nodeX[488] = 0.3491177601;

    weight[163] = 0.0109160198;
    nodeX[489] = -0.8696169152;
    nodeX[490] = 0.3491177601;
    nodeX[491] = 0.3491177601;

    weight[164] = 0.0109160198;
    nodeX[492] = 0.8696169152;
    nodeX[493] = -0.3491177601;
    nodeX[494] = 0.3491177601;

    weight[165] = 0.0109160198;
    nodeX[495] = 0.8696169152;
    nodeX[496] = 0.3491177601;
    nodeX[497] = -0.3491177601;

    weight[166] = 0.0109160198;
    nodeX[498] = -0.8696169152;
    nodeX[499] = -0.3491177601;
    nodeX[500] = 0.3491177601;

    weight[167] = 0.0109160198;
    nodeX[501] = -0.8696169152;
    nodeX[502] = 0.3491177601;
    nodeX[503] = -0.3491177601;

    weight[168] = 0.0109160198;
    nodeX[504] = 0.8696169152;
    nodeX[505] = -0.3491177601;
    nodeX[506] = -0.3491177601;

    weight[169] = 0.0109160198;
    nodeX[507] = -0.8696169152;
    nodeX[508] = -0.3491177601;
    nodeX[509] = -0.3491177601;

    weight[170] = 0.01121809491;
    nodeX[510] = 0.4121431461;
    nodeX[511] = 0.4121431461;
    nodeX[512] = 0.8125737223;

    weight[171] = 0.01121809491;
    nodeX[513] = -0.4121431461;
    nodeX[514] = 0.4121431461;
    nodeX[515] = 0.8125737223;

    weight[172] = 0.01121809491;
    nodeX[516] = 0.4121431461;
    nodeX[517] = -0.4121431461;
    nodeX[518] = 0.8125737223;

    weight[173] = 0.01121809491;
    nodeX[519] = 0.4121431461;
    nodeX[520] = 0.4121431461;
    nodeX[521] = -0.8125737223;

    weight[174] = 0.01121809491;
    nodeX[522] = -0.4121431461;
    nodeX[523] = -0.4121431461;
    nodeX[524] = 0.8125737223;

    weight[175] = 0.01121809491;
    nodeX[525] = -0.4121431461;
    nodeX[526] = 0.4121431461;
    nodeX[527] = -0.8125737223;

    weight[176] = 0.01121809491;
    nodeX[528] = 0.4121431461;
    nodeX[529] = -0.4121431461;
    nodeX[530] = -0.8125737223;

    weight[177] = 0.01121809491;
    nodeX[531] = -0.4121431461;
    nodeX[532] = -0.4121431461;
    nodeX[533] = -0.8125737223;

    weight[178] = 0.01121809491;
    nodeX[534] = -0.4121431461;
    nodeX[535] = 0.8125737223;
    nodeX[536] = 0.4121431461;

    weight[179] = 0.01121809491;
    nodeX[537] = 0.4121431461;
    nodeX[538] = -0.8125737223;
    nodeX[539] = 0.4121431461;

    weight[180] = 0.01121809491;
    nodeX[540] = 0.4121431461;
    nodeX[541] = 0.8125737223;
    nodeX[542] = -0.4121431461;

    weight[181] = 0.01121809491;
    nodeX[543] = -0.4121431461;
    nodeX[544] = -0.8125737223;
    nodeX[545] = 0.4121431461;

    weight[182] = 0.01121809491;
    nodeX[546] = -0.4121431461;
    nodeX[547] = 0.8125737223;
    nodeX[548] = -0.4121431461;

    weight[183] = 0.01121809491;
    nodeX[549] = 0.4121431461;
    nodeX[550] = -0.8125737223;
    nodeX[551] = -0.4121431461;

    weight[184] = 0.01121809491;
    nodeX[552] = -0.4121431461;
    nodeX[553] = -0.8125737223;
    nodeX[554] = -0.4121431461;

    weight[185] = 0.01121809491;
    nodeX[555] = 0.4121431461;
    nodeX[556] = 0.8125737223;
    nodeX[557] = 0.4121431461;

    weight[186] = 0.01121809491;
    nodeX[558] = 0.8125737223;
    nodeX[559] = 0.4121431461;
    nodeX[560] = 0.4121431461;

    weight[187] = 0.01121809491;
    nodeX[561] = -0.8125737223;
    nodeX[562] = 0.4121431461;
    nodeX[563] = 0.4121431461;

    weight[188] = 0.01121809491;
    nodeX[564] = 0.8125737223;
    nodeX[565] = -0.4121431461;
    nodeX[566] = 0.4121431461;

    weight[189] = 0.01121809491;
    nodeX[567] = 0.8125737223;
    nodeX[568] = 0.4121431461;
    nodeX[569] = -0.4121431461;

    weight[190] = 0.01121809491;
    nodeX[570] = -0.8125737223;
    nodeX[571] = -0.4121431461;
    nodeX[572] = 0.4121431461;

    weight[191] = 0.01121809491;
    nodeX[573] = -0.8125737223;
    nodeX[574] = 0.4121431461;
    nodeX[575] = -0.4121431461;

    weight[192] = 0.01121809491;
    nodeX[576] = 0.8125737223;
    nodeX[577] = -0.4121431461;
    nodeX[578] = -0.4121431461;

    weight[193] = 0.01121809491;
    nodeX[579] = -0.8125737223;
    nodeX[580] = -0.4121431461;
    nodeX[581] = -0.4121431461;

    weight[194] = 0.01138616252;
    nodeX[582] = 0.4718993627;
    nodeX[583] = 0.4718993627;
    nodeX[584] = 0.7447294696;

    weight[195] = 0.01138616252;
    nodeX[585] = -0.4718993627;
    nodeX[586] = 0.4718993627;
    nodeX[587] = 0.7447294696;

    weight[196] = 0.01138616252;
    nodeX[588] = 0.4718993627;
    nodeX[589] = -0.4718993627;
    nodeX[590] = 0.7447294696;

    weight[197] = 0.01138616252;
    nodeX[591] = 0.4718993627;
    nodeX[592] = 0.4718993627;
    nodeX[593] = -0.7447294696;

    weight[198] = 0.01138616252;
    nodeX[594] = -0.4718993627;
    nodeX[595] = -0.4718993627;
    nodeX[596] = 0.7447294696;

    weight[199] = 0.01138616252;
    nodeX[597] = -0.4718993627;
    nodeX[598] = 0.4718993627;
    nodeX[599] = -0.7447294696;

    weight[200] = 0.01138616252;
    nodeX[600] = 0.4718993627;
    nodeX[601] = -0.4718993627;
    nodeX[602] = -0.7447294696;

    weight[201] = 0.01138616252;
    nodeX[603] = -0.4718993627;
    nodeX[604] = -0.4718993627;
    nodeX[605] = -0.7447294696;

    weight[202] = 0.01138616252;
    nodeX[606] = -0.4718993627;
    nodeX[607] = 0.7447294696;
    nodeX[608] = 0.4718993627;

    weight[203] = 0.01138616252;
    nodeX[609] = 0.4718993627;
    nodeX[610] = -0.7447294696;
    nodeX[611] = 0.4718993627;

    weight[204] = 0.01138616252;
    nodeX[612] = 0.4718993627;
    nodeX[613] = 0.7447294696;
    nodeX[614] = -0.4718993627;

    weight[205] = 0.01138616252;
    nodeX[615] = -0.4718993627;
    nodeX[616] = -0.7447294696;
    nodeX[617] = 0.4718993627;

    weight[206] = 0.01138616252;
    nodeX[618] = -0.4718993627;
    nodeX[619] = 0.7447294696;
    nodeX[620] = -0.4718993627;

    weight[207] = 0.01138616252;
    nodeX[621] = 0.4718993627;
    nodeX[622] = -0.7447294696;
    nodeX[623] = -0.4718993627;

    weight[208] = 0.01138616252;
    nodeX[624] = -0.4718993627;
    nodeX[625] = -0.7447294696;
    nodeX[626] = -0.4718993627;

    weight[209] = 0.01138616252;
    nodeX[627] = 0.4718993627;
    nodeX[628] = 0.7447294696;
    nodeX[629] = 0.4718993627;

    weight[210] = 0.01138616252;
    nodeX[630] = 0.7447294696;
    nodeX[631] = 0.4718993627;
    nodeX[632] = 0.4718993627;

    weight[211] = 0.01138616252;
    nodeX[633] = -0.7447294696;
    nodeX[634] = 0.4718993627;
    nodeX[635] = 0.4718993627;

    weight[212] = 0.01138616252;
    nodeX[636] = 0.7447294696;
    nodeX[637] = -0.4718993627;
    nodeX[638] = 0.4718993627;

    weight[213] = 0.01138616252;
    nodeX[639] = 0.7447294696;
    nodeX[640] = 0.4718993627;
    nodeX[641] = -0.4718993627;

    weight[214] = 0.01138616252;
    nodeX[642] = -0.7447294696;
    nodeX[643] = -0.4718993627;
    nodeX[644] = 0.4718993627;

    weight[215] = 0.01138616252;
    nodeX[645] = -0.7447294696;
    nodeX[646] = 0.4718993627;
    nodeX[647] = -0.4718993627;

    weight[216] = 0.01138616252;
    nodeX[648] = 0.7447294696;
    nodeX[649] = -0.4718993627;
    nodeX[650] = -0.4718993627;

    weight[217] = 0.01138616252;
    nodeX[651] = -0.7447294696;
    nodeX[652] = -0.4718993627;
    nodeX[653] = -0.4718993627;

    weight[218] = 0.01146025009;
    nodeX[654] = 0.5273145453;
    nodeX[655] = 0.5273145453;
    nodeX[656] = 0.6662422537;

    weight[219] = 0.01146025009;
    nodeX[657] = -0.5273145453;
    nodeX[658] = 0.5273145453;
    nodeX[659] = 0.6662422537;

    weight[220] = 0.01146025009;
    nodeX[660] = 0.5273145453;
    nodeX[661] = -0.5273145453;
    nodeX[662] = 0.6662422537;

    weight[221] = 0.01146025009;
    nodeX[663] = 0.5273145453;
    nodeX[664] = 0.5273145453;
    nodeX[665] = -0.6662422537;

    weight[222] = 0.01146025009;
    nodeX[666] = -0.5273145453;
    nodeX[667] = -0.5273145453;
    nodeX[668] = 0.6662422537;

    weight[223] = 0.01146025009;
    nodeX[669] = -0.5273145453;
    nodeX[670] = 0.5273145453;
    nodeX[671] = -0.6662422537;

    weight[224] = 0.01146025009;
    nodeX[672] = 0.5273145453;
    nodeX[673] = -0.5273145453;
    nodeX[674] = -0.6662422537;

    weight[225] = 0.01146025009;
    nodeX[675] = -0.5273145453;
    nodeX[676] = -0.5273145453;
    nodeX[677] = -0.6662422537;

    weight[226] = 0.01146025009;
    nodeX[678] = -0.5273145453;
    nodeX[679] = 0.6662422537;
    nodeX[680] = 0.5273145453;

    weight[227] = 0.01146025009;
    nodeX[681] = 0.5273145453;
    nodeX[682] = -0.6662422537;
    nodeX[683] = 0.5273145453;

    weight[228] = 0.01146025009;
    nodeX[684] = 0.5273145453;
    nodeX[685] = 0.6662422537;
    nodeX[686] = -0.5273145453;

    weight[229] = 0.01146025009;
    nodeX[687] = -0.5273145453;
    nodeX[688] = -0.6662422537;
    nodeX[689] = 0.5273145453;

    weight[230] = 0.01146025009;
    nodeX[690] = -0.5273145453;
    nodeX[691] = 0.6662422537;
    nodeX[692] = -0.5273145453;

    weight[231] = 0.01146025009;
    nodeX[693] = 0.5273145453;
    nodeX[694] = -0.6662422537;
    nodeX[695] = -0.5273145453;

    weight[232] = 0.01146025009;
    nodeX[696] = -0.5273145453;
    nodeX[697] = -0.6662422537;
    nodeX[698] = -0.5273145453;

    weight[233] = 0.01146025009;
    nodeX[699] = 0.5273145453;
    nodeX[700] = 0.6662422537;
    nodeX[701] = 0.5273145453;

    weight[234] = 0.01146025009;
    nodeX[702] = 0.6662422537;
    nodeX[703] = 0.5273145453;
    nodeX[704] = 0.5273145453;

    weight[235] = 0.01146025009;
    nodeX[705] = -0.6662422537;
    nodeX[706] = 0.5273145453;
    nodeX[707] = 0.5273145453;

    weight[236] = 0.01146025009;
    nodeX[708] = 0.6662422537;
    nodeX[709] = -0.5273145453;
    nodeX[710] = 0.5273145453;

    weight[237] = 0.01146025009;
    nodeX[711] = 0.6662422537;
    nodeX[712] = 0.5273145453;
    nodeX[713] = -0.5273145453;

    weight[238] = 0.01146025009;
    nodeX[714] = -0.6662422537;
    nodeX[715] = -0.5273145453;
    nodeX[716] = 0.5273145453;

    weight[239] = 0.01146025009;
    nodeX[717] = -0.6662422537;
    nodeX[718] = 0.5273145453;
    nodeX[719] = -0.5273145453;

    weight[240] = 0.01146025009;
    nodeX[720] = 0.6662422537;
    nodeX[721] = -0.5273145453;
    nodeX[722] = -0.5273145453;

    weight[241] = 0.01146025009;
    nodeX[723] = -0.6662422537;
    nodeX[724] = -0.5273145453;
    nodeX[725] = -0.5273145453;

    weight[242] = 0.01147148805;
    nodeX[726] = 0.6209475332;
    nodeX[727] = 0.6209475332;
    nodeX[728] = 0.4783809381;

    weight[243] = 0.01147148805;
    nodeX[729] = -0.6209475332;
    nodeX[730] = 0.6209475332;
    nodeX[731] = 0.4783809381;

    weight[244] = 0.01147148805;
    nodeX[732] = 0.6209475332;
    nodeX[733] = -0.6209475332;
    nodeX[734] = 0.4783809381;

    weight[245] = 0.01147148805;
    nodeX[735] = 0.6209475332;
    nodeX[736] = 0.6209475332;
    nodeX[737] = -0.4783809381;

    weight[246] = 0.01147148805;
    nodeX[738] = -0.6209475332;
    nodeX[739] = -0.6209475332;
    nodeX[740] = 0.4783809381;

    weight[247] = 0.01147148805;
    nodeX[741] = -0.6209475332;
    nodeX[742] = 0.6209475332;
    nodeX[743] = -0.4783809381;

    weight[248] = 0.01147148805;
    nodeX[744] = 0.6209475332;
    nodeX[745] = -0.6209475332;
    nodeX[746] = -0.4783809381;

    weight[249] = 0.01147148805;
    nodeX[747] = -0.6209475332;
    nodeX[748] = -0.6209475332;
    nodeX[749] = -0.4783809381;

    weight[250] = 0.01147148805;
    nodeX[750] = -0.6209475332;
    nodeX[751] = 0.4783809381;
    nodeX[752] = 0.6209475332;

    weight[251] = 0.01147148805;
    nodeX[753] = 0.6209475332;
    nodeX[754] = -0.4783809381;
    nodeX[755] = 0.6209475332;

    weight[252] = 0.01147148805;
    nodeX[756] = 0.6209475332;
    nodeX[757] = 0.4783809381;
    nodeX[758] = -0.6209475332;

    weight[253] = 0.01147148805;
    nodeX[759] = -0.6209475332;
    nodeX[760] = -0.4783809381;
    nodeX[761] = 0.6209475332;

    weight[254] = 0.01147148805;
    nodeX[762] = -0.6209475332;
    nodeX[763] = 0.4783809381;
    nodeX[764] = -0.6209475332;

    weight[255] = 0.01147148805;
    nodeX[765] = 0.6209475332;
    nodeX[766] = -0.4783809381;
    nodeX[767] = -0.6209475332;

    weight[256] = 0.01147148805;
    nodeX[768] = -0.6209475332;
    nodeX[769] = -0.4783809381;
    nodeX[770] = -0.6209475332;

    weight[257] = 0.01147148805;
    nodeX[771] = 0.6209475332;
    nodeX[772] = 0.4783809381;
    nodeX[773] = 0.6209475332;

    weight[258] = 0.01147148805;
    nodeX[774] = 0.4783809381;
    nodeX[775] = 0.6209475332;
    nodeX[776] = 0.6209475332;

    weight[259] = 0.01147148805;
    nodeX[777] = -0.4783809381;
    nodeX[778] = 0.6209475332;
    nodeX[779] = 0.6209475332;

    weight[260] = 0.01147148805;
    nodeX[780] = 0.4783809381;
    nodeX[781] = -0.6209475332;
    nodeX[782] = 0.6209475332;

    weight[261] = 0.01147148805;
    nodeX[783] = 0.4783809381;
    nodeX[784] = 0.6209475332;
    nodeX[785] = -0.6209475332;

    weight[262] = 0.01147148805;
    nodeX[786] = -0.4783809381;
    nodeX[787] = -0.6209475332;
    nodeX[788] = 0.6209475332;

    weight[263] = 0.01147148805;
    nodeX[789] = -0.4783809381;
    nodeX[790] = 0.6209475332;
    nodeX[791] = -0.6209475332;

    weight[264] = 0.01147148805;
    nodeX[792] = 0.4783809381;
    nodeX[793] = -0.6209475332;
    nodeX[794] = -0.6209475332;

    weight[265] = 0.01147148805;
    nodeX[795] = -0.4783809381;
    nodeX[796] = -0.6209475332;
    nodeX[797] = -0.6209475332;

    weight[266] = 0.01147399479;
    nodeX[798] = 0.6569722712;
    nodeX[799] = 0.6569722712;
    nodeX[800] = 0.3698308665;

    weight[267] = 0.01147399479;
    nodeX[801] = -0.6569722712;
    nodeX[802] = 0.6569722712;
    nodeX[803] = 0.3698308665;

    weight[268] = 0.01147399479;
    nodeX[804] = 0.6569722712;
    nodeX[805] = -0.6569722712;
    nodeX[806] = 0.3698308665;

    weight[269] = 0.01147399479;
    nodeX[807] = 0.6569722712;
    nodeX[808] = 0.6569722712;
    nodeX[809] = -0.3698308665;

    weight[270] = 0.01147399479;
    nodeX[810] = -0.6569722712;
    nodeX[811] = -0.6569722712;
    nodeX[812] = 0.3698308665;

    weight[271] = 0.01147399479;
    nodeX[813] = -0.6569722712;
    nodeX[814] = 0.6569722712;
    nodeX[815] = -0.3698308665;

    weight[272] = 0.01147399479;
    nodeX[816] = 0.6569722712;
    nodeX[817] = -0.6569722712;
    nodeX[818] = -0.3698308665;

    weight[273] = 0.01147399479;
    nodeX[819] = -0.6569722712;
    nodeX[820] = -0.6569722712;
    nodeX[821] = -0.3698308665;

    weight[274] = 0.01147399479;
    nodeX[822] = -0.6569722712;
    nodeX[823] = 0.3698308665;
    nodeX[824] = 0.6569722712;

    weight[275] = 0.01147399479;
    nodeX[825] = 0.6569722712;
    nodeX[826] = -0.3698308665;
    nodeX[827] = 0.6569722712;

    weight[276] = 0.01147399479;
    nodeX[828] = 0.6569722712;
    nodeX[829] = 0.3698308665;
    nodeX[830] = -0.6569722712;

    weight[277] = 0.01147399479;
    nodeX[831] = -0.6569722712;
    nodeX[832] = -0.3698308665;
    nodeX[833] = 0.6569722712;

    weight[278] = 0.01147399479;
    nodeX[834] = -0.6569722712;
    nodeX[835] = 0.3698308665;
    nodeX[836] = -0.6569722712;

    weight[279] = 0.01147399479;
    nodeX[837] = 0.6569722712;
    nodeX[838] = -0.3698308665;
    nodeX[839] = -0.6569722712;

    weight[280] = 0.01147399479;
    nodeX[840] = -0.6569722712;
    nodeX[841] = -0.3698308665;
    nodeX[842] = -0.6569722712;

    weight[281] = 0.01147399479;
    nodeX[843] = 0.6569722712;
    nodeX[844] = 0.3698308665;
    nodeX[845] = 0.6569722712;

    weight[282] = 0.01147399479;
    nodeX[846] = 0.3698308665;
    nodeX[847] = 0.6569722712;
    nodeX[848] = 0.6569722712;

    weight[283] = 0.01147399479;
    nodeX[849] = -0.3698308665;
    nodeX[850] = 0.6569722712;
    nodeX[851] = 0.6569722712;

    weight[284] = 0.01147399479;
    nodeX[852] = 0.3698308665;
    nodeX[853] = -0.6569722712;
    nodeX[854] = 0.6569722712;

    weight[285] = 0.01147399479;
    nodeX[855] = 0.3698308665;
    nodeX[856] = 0.6569722712;
    nodeX[857] = -0.6569722712;

    weight[286] = 0.01147399479;
    nodeX[858] = -0.3698308665;
    nodeX[859] = -0.6569722712;
    nodeX[860] = 0.6569722712;

    weight[287] = 0.01147399479;
    nodeX[861] = -0.3698308665;
    nodeX[862] = 0.6569722712;
    nodeX[863] = -0.6569722712;

    weight[288] = 0.01147399479;
    nodeX[864] = 0.3698308665;
    nodeX[865] = -0.6569722712;
    nodeX[866] = -0.6569722712;

    weight[289] = 0.01147399479;
    nodeX[867] = -0.3698308665;
    nodeX[868] = -0.6569722712;
    nodeX[869] = -0.6569722712;

    weight[290] = 0.01150184042;
    nodeX[870] = 0.6841788309;
    nodeX[871] = 0.6841788309;
    nodeX[872] = 0.2525839557;

    weight[291] = 0.01150184042;
    nodeX[873] = -0.6841788309;
    nodeX[874] = 0.6841788309;
    nodeX[875] = 0.2525839557;

    weight[292] = 0.01150184042;
    nodeX[876] = 0.6841788309;
    nodeX[877] = -0.6841788309;
    nodeX[878] = 0.2525839557;

    weight[293] = 0.01150184042;
    nodeX[879] = 0.6841788309;
    nodeX[880] = 0.6841788309;
    nodeX[881] = -0.2525839557;

    weight[294] = 0.01150184042;
    nodeX[882] = -0.6841788309;
    nodeX[883] = -0.6841788309;
    nodeX[884] = 0.2525839557;

    weight[295] = 0.01150184042;
    nodeX[885] = -0.6841788309;
    nodeX[886] = 0.6841788309;
    nodeX[887] = -0.2525839557;

    weight[296] = 0.01150184042;
    nodeX[888] = 0.6841788309;
    nodeX[889] = -0.6841788309;
    nodeX[890] = -0.2525839557;

    weight[297] = 0.01150184042;
    nodeX[891] = -0.6841788309;
    nodeX[892] = -0.6841788309;
    nodeX[893] = -0.2525839557;

    weight[298] = 0.01150184042;
    nodeX[894] = -0.6841788309;
    nodeX[895] = 0.2525839557;
    nodeX[896] = 0.6841788309;

    weight[299] = 0.01150184042;
    nodeX[897] = 0.6841788309;
    nodeX[898] = -0.2525839557;
    nodeX[899] = 0.6841788309;

    weight[300] = 0.01150184042;
    nodeX[900] = 0.6841788309;
    nodeX[901] = 0.2525839557;
    nodeX[902] = -0.6841788309;

    weight[301] = 0.01150184042;
    nodeX[903] = -0.6841788309;
    nodeX[904] = -0.2525839557;
    nodeX[905] = 0.6841788309;

    weight[302] = 0.01150184042;
    nodeX[906] = -0.6841788309;
    nodeX[907] = 0.2525839557;
    nodeX[908] = -0.6841788309;

    weight[303] = 0.01150184042;
    nodeX[909] = 0.6841788309;
    nodeX[910] = -0.2525839557;
    nodeX[911] = -0.6841788309;

    weight[304] = 0.01150184042;
    nodeX[912] = -0.6841788309;
    nodeX[913] = -0.2525839557;
    nodeX[914] = -0.6841788309;

    weight[305] = 0.01150184042;
    nodeX[915] = 0.6841788309;
    nodeX[916] = 0.2525839557;
    nodeX[917] = 0.6841788309;

    weight[306] = 0.01150184042;
    nodeX[918] = 0.2525839557;
    nodeX[919] = 0.6841788309;
    nodeX[920] = 0.6841788309;

    weight[307] = 0.01150184042;
    nodeX[921] = -0.2525839557;
    nodeX[922] = 0.6841788309;
    nodeX[923] = 0.6841788309;

    weight[308] = 0.01150184042;
    nodeX[924] = 0.2525839557;
    nodeX[925] = -0.6841788309;
    nodeX[926] = 0.6841788309;

    weight[309] = 0.01150184042;
    nodeX[927] = 0.2525839557;
    nodeX[928] = 0.6841788309;
    nodeX[929] = -0.6841788309;

    weight[310] = 0.01150184042;
    nodeX[930] = -0.2525839557;
    nodeX[931] = -0.6841788309;
    nodeX[932] = 0.6841788309;

    weight[311] = 0.01150184042;
    nodeX[933] = -0.2525839557;
    nodeX[934] = 0.6841788309;
    nodeX[935] = -0.6841788309;

    weight[312] = 0.01150184042;
    nodeX[936] = 0.2525839557;
    nodeX[937] = -0.6841788309;
    nodeX[938] = -0.6841788309;

    weight[313] = 0.01150184042;
    nodeX[939] = -0.2525839557;
    nodeX[940] = -0.6841788309;
    nodeX[941] = -0.6841788309;

    weight[314] = 0.01154527292;
    nodeX[942] = 0.701260433;
    nodeX[943] = 0.701260433;
    nodeX[944] = 0.1283261867;

    weight[315] = 0.01154527292;
    nodeX[945] = -0.701260433;
    nodeX[946] = 0.701260433;
    nodeX[947] = 0.1283261867;

    weight[316] = 0.01154527292;
    nodeX[948] = 0.701260433;
    nodeX[949] = -0.701260433;
    nodeX[950] = 0.1283261867;

    weight[317] = 0.01154527292;
    nodeX[951] = 0.701260433;
    nodeX[952] = 0.701260433;
    nodeX[953] = -0.1283261867;

    weight[318] = 0.01154527292;
    nodeX[954] = -0.701260433;
    nodeX[955] = -0.701260433;
    nodeX[956] = 0.1283261867;

    weight[319] = 0.01154527292;
    nodeX[957] = -0.701260433;
    nodeX[958] = 0.701260433;
    nodeX[959] = -0.1283261867;

    weight[320] = 0.01154527292;
    nodeX[960] = 0.701260433;
    nodeX[961] = -0.701260433;
    nodeX[962] = -0.1283261867;

    weight[321] = 0.01154527292;
    nodeX[963] = -0.701260433;
    nodeX[964] = -0.701260433;
    nodeX[965] = -0.1283261867;

    weight[322] = 0.01154527292;
    nodeX[966] = -0.701260433;
    nodeX[967] = 0.1283261867;
    nodeX[968] = 0.701260433;

    weight[323] = 0.01154527292;
    nodeX[969] = 0.701260433;
    nodeX[970] = -0.1283261867;
    nodeX[971] = 0.701260433;

    weight[324] = 0.01154527292;
    nodeX[972] = 0.701260433;
    nodeX[973] = 0.1283261867;
    nodeX[974] = -0.701260433;

    weight[325] = 0.01154527292;
    nodeX[975] = -0.701260433;
    nodeX[976] = -0.1283261867;
    nodeX[977] = 0.701260433;

    weight[326] = 0.01154527292;
    nodeX[978] = -0.701260433;
    nodeX[979] = 0.1283261867;
    nodeX[980] = -0.701260433;

    weight[327] = 0.01154527292;
    nodeX[981] = 0.701260433;
    nodeX[982] = -0.1283261867;
    nodeX[983] = -0.701260433;

    weight[328] = 0.01154527292;
    nodeX[984] = -0.701260433;
    nodeX[985] = -0.1283261867;
    nodeX[986] = -0.701260433;

    weight[329] = 0.01154527292;
    nodeX[987] = 0.701260433;
    nodeX[988] = 0.1283261867;
    nodeX[989] = 0.701260433;

    weight[330] = 0.01154527292;
    nodeX[990] = 0.1283261867;
    nodeX[991] = 0.701260433;
    nodeX[992] = 0.701260433;

    weight[331] = 0.01154527292;
    nodeX[993] = -0.1283261867;
    nodeX[994] = 0.701260433;
    nodeX[995] = 0.701260433;

    weight[332] = 0.01154527292;
    nodeX[996] = 0.1283261867;
    nodeX[997] = -0.701260433;
    nodeX[998] = 0.701260433;

    weight[333] = 0.01154527292;
    nodeX[999] = 0.1283261867;
    nodeX[1000] = 0.701260433;
    nodeX[1001] = -0.701260433;

    weight[334] = 0.01154527292;
    nodeX[1002] = -0.1283261867;
    nodeX[1003] = -0.701260433;
    nodeX[1004] = 0.701260433;

    weight[335] = 0.01154527292;
    nodeX[1005] = -0.1283261867;
    nodeX[1006] = 0.701260433;
    nodeX[1007] = -0.701260433;

    weight[336] = 0.01154527292;
    nodeX[1008] = 0.1283261867;
    nodeX[1009] = -0.701260433;
    nodeX[1010] = -0.701260433;

    weight[337] = 0.01154527292;
    nodeX[1011] = -0.1283261867;
    nodeX[1012] = -0.701260433;
    nodeX[1013] = -0.701260433;

    weight[338] = 0.006505581558;
    nodeX[1014] = 0.1072382215;
    nodeX[1015] = 0.9942333548;
    nodeX[1016] = 0;

    weight[339] = 0.006505581558;
    nodeX[1017] = -0.1072382215;
    nodeX[1018] = 0.9942333548;
    nodeX[1019] = 0;

    weight[340] = 0.006505581558;
    nodeX[1020] = 0.1072382215;
    nodeX[1021] = -0.9942333548;
    nodeX[1022] = 0;

    weight[341] = 0.006505581558;
    nodeX[1023] = -0.1072382215;
    nodeX[1024] = -0.9942333548;
    nodeX[1025] = 0;

    weight[342] = 0.006505581558;
    nodeX[1026] = 0.9942333548;
    nodeX[1027] = 0.1072382215;
    nodeX[1028] = 0;

    weight[343] = 0.006505581558;
    nodeX[1029] = -0.9942333548;
    nodeX[1030] = 0.1072382215;
    nodeX[1031] = 0;

    weight[344] = 0.006505581558;
    nodeX[1032] = 0.9942333548;
    nodeX[1033] = -0.1072382215;
    nodeX[1034] = 0;

    weight[345] = 0.006505581558;
    nodeX[1035] = -0.9942333548;
    nodeX[1036] = -0.1072382215;
    nodeX[1037] = 0;

    weight[346] = 0.006505581558;
    nodeX[1038] = 0.1072382215;
    nodeX[1039] = 0;
    nodeX[1040] = 0.9942333548;

    weight[347] = 0.006505581558;
    nodeX[1041] = -0.1072382215;
    nodeX[1042] = 0;
    nodeX[1043] = 0.9942333548;

    weight[348] = 0.006505581558;
    nodeX[1044] = 0.1072382215;
    nodeX[1045] = 0;
    nodeX[1046] = -0.9942333548;

    weight[349] = 0.006505581558;
    nodeX[1047] = -0.1072382215;
    nodeX[1048] = 0;
    nodeX[1049] = -0.9942333548;

    weight[350] = 0.006505581558;
    nodeX[1050] = 0.9942333548;
    nodeX[1051] = 0;
    nodeX[1052] = 0.1072382215;

    weight[351] = 0.006505581558;
    nodeX[1053] = -0.9942333548;
    nodeX[1054] = 0;
    nodeX[1055] = 0.1072382215;

    weight[352] = 0.006505581558;
    nodeX[1056] = 0.9942333548;
    nodeX[1057] = 0;
    nodeX[1058] = -0.1072382215;

    weight[353] = 0.006505581558;
    nodeX[1059] = -0.9942333548;
    nodeX[1060] = 0;
    nodeX[1061] = -0.1072382215;

    weight[354] = 0.006505581558;
    nodeX[1062] = 0;
    nodeX[1063] = 0.1072382215;
    nodeX[1064] = 0.9942333548;

    weight[355] = 0.006505581558;
    nodeX[1065] = 0;
    nodeX[1066] = -0.1072382215;
    nodeX[1067] = 0.9942333548;

    weight[356] = 0.006505581558;
    nodeX[1068] = 0;
    nodeX[1069] = 0.1072382215;
    nodeX[1070] = -0.9942333548;

    weight[357] = 0.006505581558;
    nodeX[1071] = 0;
    nodeX[1072] = -0.1072382215;
    nodeX[1073] = -0.9942333548;

    weight[358] = 0.006505581558;
    nodeX[1074] = 0;
    nodeX[1075] = 0.9942333548;
    nodeX[1076] = 0.1072382215;

    weight[359] = 0.006505581558;
    nodeX[1077] = 0;
    nodeX[1078] = -0.9942333548;
    nodeX[1079] = 0.1072382215;

    weight[360] = 0.006505581558;
    nodeX[1080] = 0;
    nodeX[1081] = 0.9942333548;
    nodeX[1082] = -0.1072382215;

    weight[361] = 0.006505581558;
    nodeX[1083] = 0;
    nodeX[1084] = -0.9942333548;
    nodeX[1085] = -0.1072382215;

    weight[362] = 0.009212586854;
    nodeX[1086] = 0.2582068959;
    nodeX[1087] = 0.9660896433;
    nodeX[1088] = 0;

    weight[363] = 0.009212586854;
    nodeX[1089] = -0.2582068959;
    nodeX[1090] = 0.9660896433;
    nodeX[1091] = 0;

    weight[364] = 0.009212586854;
    nodeX[1092] = 0.2582068959;
    nodeX[1093] = -0.9660896433;
    nodeX[1094] = 0;

    weight[365] = 0.009212586854;
    nodeX[1095] = -0.2582068959;
    nodeX[1096] = -0.9660896433;
    nodeX[1097] = 0;

    weight[366] = 0.009212586854;
    nodeX[1098] = 0.9660896433;
    nodeX[1099] = 0.2582068959;
    nodeX[1100] = 0;

    weight[367] = 0.009212586854;
    nodeX[1101] = -0.9660896433;
    nodeX[1102] = 0.2582068959;
    nodeX[1103] = 0;

    weight[368] = 0.009212586854;
    nodeX[1104] = 0.9660896433;
    nodeX[1105] = -0.2582068959;
    nodeX[1106] = 0;

    weight[369] = 0.009212586854;
    nodeX[1107] = -0.9660896433;
    nodeX[1108] = -0.2582068959;
    nodeX[1109] = 0;

    weight[370] = 0.009212586854;
    nodeX[1110] = 0.2582068959;
    nodeX[1111] = 0;
    nodeX[1112] = 0.9660896433;

    weight[371] = 0.009212586854;
    nodeX[1113] = -0.2582068959;
    nodeX[1114] = 0;
    nodeX[1115] = 0.9660896433;

    weight[372] = 0.009212586854;
    nodeX[1116] = 0.2582068959;
    nodeX[1117] = 0;
    nodeX[1118] = -0.9660896433;

    weight[373] = 0.009212586854;
    nodeX[1119] = -0.2582068959;
    nodeX[1120] = 0;
    nodeX[1121] = -0.9660896433;

    weight[374] = 0.009212586854;
    nodeX[1122] = 0.9660896433;
    nodeX[1123] = 0;
    nodeX[1124] = 0.2582068959;

    weight[375] = 0.009212586854;
    nodeX[1125] = -0.9660896433;
    nodeX[1126] = 0;
    nodeX[1127] = 0.2582068959;

    weight[376] = 0.009212586854;
    nodeX[1128] = 0.9660896433;
    nodeX[1129] = 0;
    nodeX[1130] = -0.2582068959;

    weight[377] = 0.009212586854;
    nodeX[1131] = -0.9660896433;
    nodeX[1132] = 0;
    nodeX[1133] = -0.2582068959;

    weight[378] = 0.009212586854;
    nodeX[1134] = 0;
    nodeX[1135] = 0.2582068959;
    nodeX[1136] = 0.9660896433;

    weight[379] = 0.009212586854;
    nodeX[1137] = 0;
    nodeX[1138] = -0.2582068959;
    nodeX[1139] = 0.9660896433;

    weight[380] = 0.009212586854;
    nodeX[1140] = 0;
    nodeX[1141] = 0.2582068959;
    nodeX[1142] = -0.9660896433;

    weight[381] = 0.009212586854;
    nodeX[1143] = 0;
    nodeX[1144] = -0.2582068959;
    nodeX[1145] = -0.9660896433;

    weight[382] = 0.009212586854;
    nodeX[1146] = 0;
    nodeX[1147] = 0.9660896433;
    nodeX[1148] = 0.2582068959;

    weight[383] = 0.009212586854;
    nodeX[1149] = 0;
    nodeX[1150] = -0.9660896433;
    nodeX[1151] = 0.2582068959;

    weight[384] = 0.009212586854;
    nodeX[1152] = 0;
    nodeX[1153] = 0.9660896433;
    nodeX[1154] = -0.2582068959;

    weight[385] = 0.009212586854;
    nodeX[1155] = 0;
    nodeX[1156] = -0.9660896433;
    nodeX[1157] = -0.2582068959;

    weight[386] = 0.01063521204;
    nodeX[1158] = 0.4172752955;
    nodeX[1159] = 0.9087801317;
    nodeX[1160] = 0;

    weight[387] = 0.01063521204;
    nodeX[1161] = -0.4172752955;
    nodeX[1162] = 0.9087801317;
    nodeX[1163] = 0;

    weight[388] = 0.01063521204;
    nodeX[1164] = 0.4172752955;
    nodeX[1165] = -0.9087801317;
    nodeX[1166] = 0;

    weight[389] = 0.01063521204;
    nodeX[1167] = -0.4172752955;
    nodeX[1168] = -0.9087801317;
    nodeX[1169] = 0;

    weight[390] = 0.01063521204;
    nodeX[1170] = 0.9087801317;
    nodeX[1171] = 0.4172752955;
    nodeX[1172] = 0;

    weight[391] = 0.01063521204;
    nodeX[1173] = -0.9087801317;
    nodeX[1174] = 0.4172752955;
    nodeX[1175] = 0;

    weight[392] = 0.01063521204;
    nodeX[1176] = 0.9087801317;
    nodeX[1177] = -0.4172752955;
    nodeX[1178] = 0;

    weight[393] = 0.01063521204;
    nodeX[1179] = -0.9087801317;
    nodeX[1180] = -0.4172752955;
    nodeX[1181] = 0;

    weight[394] = 0.01063521204;
    nodeX[1182] = 0.4172752955;
    nodeX[1183] = 0;
    nodeX[1184] = 0.9087801317;

    weight[395] = 0.01063521204;
    nodeX[1185] = -0.4172752955;
    nodeX[1186] = 0;
    nodeX[1187] = 0.9087801317;

    weight[396] = 0.01063521204;
    nodeX[1188] = 0.4172752955;
    nodeX[1189] = 0;
    nodeX[1190] = -0.9087801317;

    weight[397] = 0.01063521204;
    nodeX[1191] = -0.4172752955;
    nodeX[1192] = 0;
    nodeX[1193] = -0.9087801317;

    weight[398] = 0.01063521204;
    nodeX[1194] = 0.9087801317;
    nodeX[1195] = 0;
    nodeX[1196] = 0.4172752955;

    weight[399] = 0.01063521204;
    nodeX[1197] = -0.9087801317;
    nodeX[1198] = 0;
    nodeX[1199] = 0.4172752955;

    weight[400] = 0.01063521204;
    nodeX[1200] = 0.9087801317;
    nodeX[1201] = 0;
    nodeX[1202] = -0.4172752955;

    weight[401] = 0.01063521204;
    nodeX[1203] = -0.9087801317;
    nodeX[1204] = 0;
    nodeX[1205] = -0.4172752955;

    weight[402] = 0.01063521204;
    nodeX[1206] = 0;
    nodeX[1207] = 0.4172752955;
    nodeX[1208] = 0.9087801317;

    weight[403] = 0.01063521204;
    nodeX[1209] = 0;
    nodeX[1210] = -0.4172752955;
    nodeX[1211] = 0.9087801317;

    weight[404] = 0.01063521204;
    nodeX[1212] = 0;
    nodeX[1213] = 0.4172752955;
    nodeX[1214] = -0.9087801317;

    weight[405] = 0.01063521204;
    nodeX[1215] = 0;
    nodeX[1216] = -0.4172752955;
    nodeX[1217] = -0.9087801317;

    weight[406] = 0.01063521204;
    nodeX[1218] = 0;
    nodeX[1219] = 0.9087801317;
    nodeX[1220] = 0.4172752955;

    weight[407] = 0.01063521204;
    nodeX[1221] = 0;
    nodeX[1222] = -0.9087801317;
    nodeX[1223] = 0.4172752955;

    weight[408] = 0.01063521204;
    nodeX[1224] = 0;
    nodeX[1225] = 0.9087801317;
    nodeX[1226] = -0.4172752955;

    weight[409] = 0.01063521204;
    nodeX[1227] = 0;
    nodeX[1228] = -0.9087801317;
    nodeX[1229] = -0.4172752955;

    weight[410] = 0.01134884348;
    nodeX[1230] = 0.5700366912;
    nodeX[1231] = 0.8216192371;
    nodeX[1232] = 0;

    weight[411] = 0.01134884348;
    nodeX[1233] = -0.5700366912;
    nodeX[1234] = 0.8216192371;
    nodeX[1235] = 0;

    weight[412] = 0.01134884348;
    nodeX[1236] = 0.5700366912;
    nodeX[1237] = -0.8216192371;
    nodeX[1238] = 0;

    weight[413] = 0.01134884348;
    nodeX[1239] = -0.5700366912;
    nodeX[1240] = -0.8216192371;
    nodeX[1241] = 0;

    weight[414] = 0.01134884348;
    nodeX[1242] = 0.8216192371;
    nodeX[1243] = 0.5700366912;
    nodeX[1244] = 0;

    weight[415] = 0.01134884348;
    nodeX[1245] = -0.8216192371;
    nodeX[1246] = 0.5700366912;
    nodeX[1247] = 0;

    weight[416] = 0.01134884348;
    nodeX[1248] = 0.8216192371;
    nodeX[1249] = -0.5700366912;
    nodeX[1250] = 0;

    weight[417] = 0.01134884348;
    nodeX[1251] = -0.8216192371;
    nodeX[1252] = -0.5700366912;
    nodeX[1253] = 0;

    weight[418] = 0.01134884348;
    nodeX[1254] = 0.5700366912;
    nodeX[1255] = 0;
    nodeX[1256] = 0.8216192371;

    weight[419] = 0.01134884348;
    nodeX[1257] = -0.5700366912;
    nodeX[1258] = 0;
    nodeX[1259] = 0.8216192371;

    weight[420] = 0.01134884348;
    nodeX[1260] = 0.5700366912;
    nodeX[1261] = 0;
    nodeX[1262] = -0.8216192371;

    weight[421] = 0.01134884348;
    nodeX[1263] = -0.5700366912;
    nodeX[1264] = 0;
    nodeX[1265] = -0.8216192371;

    weight[422] = 0.01134884348;
    nodeX[1266] = 0.8216192371;
    nodeX[1267] = 0;
    nodeX[1268] = 0.5700366912;

    weight[423] = 0.01134884348;
    nodeX[1269] = -0.8216192371;
    nodeX[1270] = 0;
    nodeX[1271] = 0.5700366912;

    weight[424] = 0.01134884348;
    nodeX[1272] = 0.8216192371;
    nodeX[1273] = 0;
    nodeX[1274] = -0.5700366912;

    weight[425] = 0.01134884348;
    nodeX[1275] = -0.8216192371;
    nodeX[1276] = 0;
    nodeX[1277] = -0.5700366912;

    weight[426] = 0.01134884348;
    nodeX[1278] = 0;
    nodeX[1279] = 0.5700366912;
    nodeX[1280] = 0.8216192371;

    weight[427] = 0.01134884348;
    nodeX[1281] = 0;
    nodeX[1282] = -0.5700366912;
    nodeX[1283] = 0.8216192371;

    weight[428] = 0.01134884348;
    nodeX[1284] = 0;
    nodeX[1285] = 0.5700366912;
    nodeX[1286] = -0.8216192371;

    weight[429] = 0.01134884348;
    nodeX[1287] = 0;
    nodeX[1288] = -0.5700366912;
    nodeX[1289] = -0.8216192371;

    weight[430] = 0.01134884348;
    nodeX[1290] = 0;
    nodeX[1291] = 0.8216192371;
    nodeX[1292] = 0.5700366912;

    weight[431] = 0.01134884348;
    nodeX[1293] = 0;
    nodeX[1294] = -0.8216192371;
    nodeX[1295] = 0.5700366912;

    weight[432] = 0.01134884348;
    nodeX[1296] = 0;
    nodeX[1297] = 0.8216192371;
    nodeX[1298] = -0.5700366912;

    weight[433] = 0.01134884348;
    nodeX[1299] = 0;
    nodeX[1300] = -0.8216192371;
    nodeX[1301] = -0.5700366912;

    weight[434] = 0.008150269577;
    nodeX[1302] = 0.9827986018;
    nodeX[1303] = 0.1771774023;
    nodeX[1304] = 0.05210639477;

    weight[435] = 0.008150269577;
    nodeX[1305] = -0.9827986018;
    nodeX[1306] = 0.1771774023;
    nodeX[1307] = 0.05210639477;

    weight[436] = 0.008150269577;
    nodeX[1308] = 0.9827986018;
    nodeX[1309] = -0.1771774023;
    nodeX[1310] = 0.05210639477;

    weight[437] = 0.008150269577;
    nodeX[1311] = 0.9827986018;
    nodeX[1312] = 0.1771774023;
    nodeX[1313] = -0.05210639477;

    weight[438] = 0.008150269577;
    nodeX[1314] = -0.9827986018;
    nodeX[1315] = -0.1771774023;
    nodeX[1316] = 0.05210639477;

    weight[439] = 0.008150269577;
    nodeX[1317] = 0.9827986018;
    nodeX[1318] = -0.1771774023;
    nodeX[1319] = -0.05210639477;

    weight[440] = 0.008150269577;
    nodeX[1320] = -0.9827986018;
    nodeX[1321] = 0.1771774023;
    nodeX[1322] = -0.05210639477;

    weight[441] = 0.008150269577;
    nodeX[1323] = -0.9827986018;
    nodeX[1324] = -0.1771774023;
    nodeX[1325] = -0.05210639477;

    weight[442] = 0.008150269577;
    nodeX[1326] = 0.1771774023;
    nodeX[1327] = 0.9827986018;
    nodeX[1328] = 0.05210639477;

    weight[443] = 0.008150269577;
    nodeX[1329] = -0.1771774023;
    nodeX[1330] = 0.9827986018;
    nodeX[1331] = 0.05210639477;

    weight[444] = 0.008150269577;
    nodeX[1332] = 0.1771774023;
    nodeX[1333] = -0.9827986018;
    nodeX[1334] = 0.05210639477;

    weight[445] = 0.008150269577;
    nodeX[1335] = 0.1771774023;
    nodeX[1336] = 0.9827986018;
    nodeX[1337] = -0.05210639477;

    weight[446] = 0.008150269577;
    nodeX[1338] = -0.1771774023;
    nodeX[1339] = -0.9827986018;
    nodeX[1340] = 0.05210639477;

    weight[447] = 0.008150269577;
    nodeX[1341] = 0.1771774023;
    nodeX[1342] = -0.9827986018;
    nodeX[1343] = -0.05210639477;

    weight[448] = 0.008150269577;
    nodeX[1344] = -0.1771774023;
    nodeX[1345] = 0.9827986018;
    nodeX[1346] = -0.05210639477;

    weight[449] = 0.008150269577;
    nodeX[1347] = -0.1771774023;
    nodeX[1348] = -0.9827986018;
    nodeX[1349] = -0.05210639477;

    weight[450] = 0.008150269577;
    nodeX[1350] = 0.05210639477;
    nodeX[1351] = 0.9827986018;
    nodeX[1352] = 0.1771774023;

    weight[451] = 0.008150269577;
    nodeX[1353] = -0.05210639477;
    nodeX[1354] = 0.9827986018;
    nodeX[1355] = 0.1771774023;

    weight[452] = 0.008150269577;
    nodeX[1356] = 0.05210639477;
    nodeX[1357] = -0.9827986018;
    nodeX[1358] = 0.1771774023;

    weight[453] = 0.008150269577;
    nodeX[1359] = 0.05210639477;
    nodeX[1360] = 0.9827986018;
    nodeX[1361] = -0.1771774023;

    weight[454] = 0.008150269577;
    nodeX[1362] = -0.05210639477;
    nodeX[1363] = -0.9827986018;
    nodeX[1364] = 0.1771774023;

    weight[455] = 0.008150269577;
    nodeX[1365] = 0.05210639477;
    nodeX[1366] = -0.9827986018;
    nodeX[1367] = -0.1771774023;

    weight[456] = 0.008150269577;
    nodeX[1368] = -0.05210639477;
    nodeX[1369] = 0.9827986018;
    nodeX[1370] = -0.1771774023;

    weight[457] = 0.008150269577;
    nodeX[1371] = -0.05210639477;
    nodeX[1372] = -0.9827986018;
    nodeX[1373] = -0.1771774023;

    weight[458] = 0.008150269577;
    nodeX[1374] = 0.05210639477;
    nodeX[1375] = 0.1771774023;
    nodeX[1376] = 0.9827986018;

    weight[459] = 0.008150269577;
    nodeX[1377] = -0.05210639477;
    nodeX[1378] = 0.1771774023;
    nodeX[1379] = 0.9827986018;

    weight[460] = 0.008150269577;
    nodeX[1380] = 0.05210639477;
    nodeX[1381] = -0.1771774023;
    nodeX[1382] = 0.9827986018;

    weight[461] = 0.008150269577;
    nodeX[1383] = 0.05210639477;
    nodeX[1384] = 0.1771774023;
    nodeX[1385] = -0.9827986018;

    weight[462] = 0.008150269577;
    nodeX[1386] = -0.05210639477;
    nodeX[1387] = -0.1771774023;
    nodeX[1388] = 0.9827986018;

    weight[463] = 0.008150269577;
    nodeX[1389] = 0.05210639477;
    nodeX[1390] = -0.1771774023;
    nodeX[1391] = -0.9827986018;

    weight[464] = 0.008150269577;
    nodeX[1392] = -0.05210639477;
    nodeX[1393] = 0.1771774023;
    nodeX[1394] = -0.9827986018;

    weight[465] = 0.008150269577;
    nodeX[1395] = -0.05210639477;
    nodeX[1396] = -0.1771774023;
    nodeX[1397] = -0.9827986018;

    weight[466] = 0.008150269577;
    nodeX[1398] = 0.9827986018;
    nodeX[1399] = 0.05210639477;
    nodeX[1400] = 0.1771774023;

    weight[467] = 0.008150269577;
    nodeX[1401] = -0.9827986018;
    nodeX[1402] = 0.05210639477;
    nodeX[1403] = 0.1771774023;

    weight[468] = 0.008150269577;
    nodeX[1404] = 0.9827986018;
    nodeX[1405] = -0.05210639477;
    nodeX[1406] = 0.1771774023;

    weight[469] = 0.008150269577;
    nodeX[1407] = 0.9827986018;
    nodeX[1408] = 0.05210639477;
    nodeX[1409] = -0.1771774023;

    weight[470] = 0.008150269577;
    nodeX[1410] = -0.9827986018;
    nodeX[1411] = -0.05210639477;
    nodeX[1412] = 0.1771774023;

    weight[471] = 0.008150269577;
    nodeX[1413] = 0.9827986018;
    nodeX[1414] = -0.05210639477;
    nodeX[1415] = -0.1771774023;

    weight[472] = 0.008150269577;
    nodeX[1416] = -0.9827986018;
    nodeX[1417] = 0.05210639477;
    nodeX[1418] = -0.1771774023;

    weight[473] = 0.008150269577;
    nodeX[1419] = -0.9827986018;
    nodeX[1420] = -0.05210639477;
    nodeX[1421] = -0.1771774023;

    weight[474] = 0.008150269577;
    nodeX[1422] = 0.1771774023;
    nodeX[1423] = 0.05210639477;
    nodeX[1424] = 0.9827986018;

    weight[475] = 0.008150269577;
    nodeX[1425] = -0.1771774023;
    nodeX[1426] = 0.05210639477;
    nodeX[1427] = 0.9827986018;

    weight[476] = 0.008150269577;
    nodeX[1428] = 0.1771774023;
    nodeX[1429] = -0.05210639477;
    nodeX[1430] = 0.9827986018;

    weight[477] = 0.008150269577;
    nodeX[1431] = 0.1771774023;
    nodeX[1432] = 0.05210639477;
    nodeX[1433] = -0.9827986018;

    weight[478] = 0.008150269577;
    nodeX[1434] = -0.1771774023;
    nodeX[1435] = -0.05210639477;
    nodeX[1436] = 0.9827986018;

    weight[479] = 0.008150269577;
    nodeX[1437] = 0.1771774023;
    nodeX[1438] = -0.05210639477;
    nodeX[1439] = -0.9827986018;

    weight[480] = 0.008150269577;
    nodeX[1440] = -0.1771774023;
    nodeX[1441] = 0.05210639477;
    nodeX[1442] = -0.9827986018;

    weight[481] = 0.008150269577;
    nodeX[1443] = -0.1771774023;
    nodeX[1444] = -0.05210639477;
    nodeX[1445] = -0.9827986018;

    weight[482] = 0.009343135396;
    nodeX[1446] = 0.962424923;
    nodeX[1447] = 0.2475716463;
    nodeX[1448] = 0.1115640957;

    weight[483] = 0.009343135396;
    nodeX[1449] = -0.962424923;
    nodeX[1450] = 0.2475716463;
    nodeX[1451] = 0.1115640957;

    weight[484] = 0.009343135396;
    nodeX[1452] = 0.962424923;
    nodeX[1453] = -0.2475716463;
    nodeX[1454] = 0.1115640957;

    weight[485] = 0.009343135396;
    nodeX[1455] = 0.962424923;
    nodeX[1456] = 0.2475716463;
    nodeX[1457] = -0.1115640957;

    weight[486] = 0.009343135396;
    nodeX[1458] = -0.962424923;
    nodeX[1459] = -0.2475716463;
    nodeX[1460] = 0.1115640957;

    weight[487] = 0.009343135396;
    nodeX[1461] = 0.962424923;
    nodeX[1462] = -0.2475716463;
    nodeX[1463] = -0.1115640957;

    weight[488] = 0.009343135396;
    nodeX[1464] = -0.962424923;
    nodeX[1465] = 0.2475716463;
    nodeX[1466] = -0.1115640957;

    weight[489] = 0.009343135396;
    nodeX[1467] = -0.962424923;
    nodeX[1468] = -0.2475716463;
    nodeX[1469] = -0.1115640957;

    weight[490] = 0.009343135396;
    nodeX[1470] = 0.2475716463;
    nodeX[1471] = 0.962424923;
    nodeX[1472] = 0.1115640957;

    weight[491] = 0.009343135396;
    nodeX[1473] = -0.2475716463;
    nodeX[1474] = 0.962424923;
    nodeX[1475] = 0.1115640957;

    weight[492] = 0.009343135396;
    nodeX[1476] = 0.2475716463;
    nodeX[1477] = -0.962424923;
    nodeX[1478] = 0.1115640957;

    weight[493] = 0.009343135396;
    nodeX[1479] = 0.2475716463;
    nodeX[1480] = 0.962424923;
    nodeX[1481] = -0.1115640957;

    weight[494] = 0.009343135396;
    nodeX[1482] = -0.2475716463;
    nodeX[1483] = -0.962424923;
    nodeX[1484] = 0.1115640957;

    weight[495] = 0.009343135396;
    nodeX[1485] = 0.2475716463;
    nodeX[1486] = -0.962424923;
    nodeX[1487] = -0.1115640957;

    weight[496] = 0.009343135396;
    nodeX[1488] = -0.2475716463;
    nodeX[1489] = 0.962424923;
    nodeX[1490] = -0.1115640957;

    weight[497] = 0.009343135396;
    nodeX[1491] = -0.2475716463;
    nodeX[1492] = -0.962424923;
    nodeX[1493] = -0.1115640957;

    weight[498] = 0.009343135396;
    nodeX[1494] = 0.1115640957;
    nodeX[1495] = 0.962424923;
    nodeX[1496] = 0.2475716463;

    weight[499] = 0.009343135396;
    nodeX[1497] = -0.1115640957;
    nodeX[1498] = 0.962424923;
    nodeX[1499] = 0.2475716463;

    weight[500] = 0.009343135396;
    nodeX[1500] = 0.1115640957;
    nodeX[1501] = -0.962424923;
    nodeX[1502] = 0.2475716463;

    weight[501] = 0.009343135396;
    nodeX[1503] = 0.1115640957;
    nodeX[1504] = 0.962424923;
    nodeX[1505] = -0.2475716463;

    weight[502] = 0.009343135396;
    nodeX[1506] = -0.1115640957;
    nodeX[1507] = -0.962424923;
    nodeX[1508] = 0.2475716463;

    weight[503] = 0.009343135396;
    nodeX[1509] = 0.1115640957;
    nodeX[1510] = -0.962424923;
    nodeX[1511] = -0.2475716463;

    weight[504] = 0.009343135396;
    nodeX[1512] = -0.1115640957;
    nodeX[1513] = 0.962424923;
    nodeX[1514] = -0.2475716463;

    weight[505] = 0.009343135396;
    nodeX[1515] = -0.1115640957;
    nodeX[1516] = -0.962424923;
    nodeX[1517] = -0.2475716463;

    weight[506] = 0.009343135396;
    nodeX[1518] = 0.1115640957;
    nodeX[1519] = 0.2475716463;
    nodeX[1520] = 0.962424923;

    weight[507] = 0.009343135396;
    nodeX[1521] = -0.1115640957;
    nodeX[1522] = 0.2475716463;
    nodeX[1523] = 0.962424923;

    weight[508] = 0.009343135396;
    nodeX[1524] = 0.1115640957;
    nodeX[1525] = -0.2475716463;
    nodeX[1526] = 0.962424923;

    weight[509] = 0.009343135396;
    nodeX[1527] = 0.1115640957;
    nodeX[1528] = 0.2475716463;
    nodeX[1529] = -0.962424923;

    weight[510] = 0.009343135396;
    nodeX[1530] = -0.1115640957;
    nodeX[1531] = -0.2475716463;
    nodeX[1532] = 0.962424923;

    weight[511] = 0.009343135396;
    nodeX[1533] = 0.1115640957;
    nodeX[1534] = -0.2475716463;
    nodeX[1535] = -0.962424923;

    weight[512] = 0.009343135396;
    nodeX[1536] = -0.1115640957;
    nodeX[1537] = 0.2475716463;
    nodeX[1538] = -0.962424923;

    weight[513] = 0.009343135396;
    nodeX[1539] = -0.1115640957;
    nodeX[1540] = -0.2475716463;
    nodeX[1541] = -0.962424923;

    weight[514] = 0.009343135396;
    nodeX[1542] = 0.962424923;
    nodeX[1543] = 0.1115640957;
    nodeX[1544] = 0.2475716463;

    weight[515] = 0.009343135396;
    nodeX[1545] = -0.962424923;
    nodeX[1546] = 0.1115640957;
    nodeX[1547] = 0.2475716463;

    weight[516] = 0.009343135396;
    nodeX[1548] = 0.962424923;
    nodeX[1549] = -0.1115640957;
    nodeX[1550] = 0.2475716463;

    weight[517] = 0.009343135396;
    nodeX[1551] = 0.962424923;
    nodeX[1552] = 0.1115640957;
    nodeX[1553] = -0.2475716463;

    weight[518] = 0.009343135396;
    nodeX[1554] = -0.962424923;
    nodeX[1555] = -0.1115640957;
    nodeX[1556] = 0.2475716463;

    weight[519] = 0.009343135396;
    nodeX[1557] = 0.962424923;
    nodeX[1558] = -0.1115640957;
    nodeX[1559] = -0.2475716463;

    weight[520] = 0.009343135396;
    nodeX[1560] = -0.962424923;
    nodeX[1561] = 0.1115640957;
    nodeX[1562] = -0.2475716463;

    weight[521] = 0.009343135396;
    nodeX[1563] = -0.962424923;
    nodeX[1564] = -0.1115640957;
    nodeX[1565] = -0.2475716463;

    weight[522] = 0.009343135396;
    nodeX[1566] = 0.2475716463;
    nodeX[1567] = 0.1115640957;
    nodeX[1568] = 0.962424923;

    weight[523] = 0.009343135396;
    nodeX[1569] = -0.2475716463;
    nodeX[1570] = 0.1115640957;
    nodeX[1571] = 0.962424923;

    weight[524] = 0.009343135396;
    nodeX[1572] = 0.2475716463;
    nodeX[1573] = -0.1115640957;
    nodeX[1574] = 0.962424923;

    weight[525] = 0.009343135396;
    nodeX[1575] = 0.2475716463;
    nodeX[1576] = 0.1115640957;
    nodeX[1577] = -0.962424923;

    weight[526] = 0.009343135396;
    nodeX[1578] = -0.2475716463;
    nodeX[1579] = -0.1115640957;
    nodeX[1580] = 0.962424923;

    weight[527] = 0.009343135396;
    nodeX[1581] = 0.2475716463;
    nodeX[1582] = -0.1115640957;
    nodeX[1583] = -0.962424923;

    weight[528] = 0.009343135396;
    nodeX[1584] = -0.2475716463;
    nodeX[1585] = 0.1115640957;
    nodeX[1586] = -0.962424923;

    weight[529] = 0.009343135396;
    nodeX[1587] = -0.2475716463;
    nodeX[1588] = -0.1115640957;
    nodeX[1589] = -0.962424923;

    weight[530] = 0.01005124659;
    nodeX[1590] = 0.9402007994;
    nodeX[1591] = 0.3354616289;
    nodeX[1592] = 0.05905888853;

    weight[531] = 0.01005124659;
    nodeX[1593] = -0.9402007994;
    nodeX[1594] = 0.3354616289;
    nodeX[1595] = 0.05905888853;

    weight[532] = 0.01005124659;
    nodeX[1596] = 0.9402007994;
    nodeX[1597] = -0.3354616289;
    nodeX[1598] = 0.05905888853;

    weight[533] = 0.01005124659;
    nodeX[1599] = 0.9402007994;
    nodeX[1600] = 0.3354616289;
    nodeX[1601] = -0.05905888853;

    weight[534] = 0.01005124659;
    nodeX[1602] = -0.9402007994;
    nodeX[1603] = -0.3354616289;
    nodeX[1604] = 0.05905888853;

    weight[535] = 0.01005124659;
    nodeX[1605] = 0.9402007994;
    nodeX[1606] = -0.3354616289;
    nodeX[1607] = -0.05905888853;

    weight[536] = 0.01005124659;
    nodeX[1608] = -0.9402007994;
    nodeX[1609] = 0.3354616289;
    nodeX[1610] = -0.05905888853;

    weight[537] = 0.01005124659;
    nodeX[1611] = -0.9402007994;
    nodeX[1612] = -0.3354616289;
    nodeX[1613] = -0.05905888853;

    weight[538] = 0.01005124659;
    nodeX[1614] = 0.3354616289;
    nodeX[1615] = 0.9402007994;
    nodeX[1616] = 0.05905888853;

    weight[539] = 0.01005124659;
    nodeX[1617] = -0.3354616289;
    nodeX[1618] = 0.9402007994;
    nodeX[1619] = 0.05905888853;

    weight[540] = 0.01005124659;
    nodeX[1620] = 0.3354616289;
    nodeX[1621] = -0.9402007994;
    nodeX[1622] = 0.05905888853;

    weight[541] = 0.01005124659;
    nodeX[1623] = 0.3354616289;
    nodeX[1624] = 0.9402007994;
    nodeX[1625] = -0.05905888853;

    weight[542] = 0.01005124659;
    nodeX[1626] = -0.3354616289;
    nodeX[1627] = -0.9402007994;
    nodeX[1628] = 0.05905888853;

    weight[543] = 0.01005124659;
    nodeX[1629] = 0.3354616289;
    nodeX[1630] = -0.9402007994;
    nodeX[1631] = -0.05905888853;

    weight[544] = 0.01005124659;
    nodeX[1632] = -0.3354616289;
    nodeX[1633] = 0.9402007994;
    nodeX[1634] = -0.05905888853;

    weight[545] = 0.01005124659;
    nodeX[1635] = -0.3354616289;
    nodeX[1636] = -0.9402007994;
    nodeX[1637] = -0.05905888853;

    weight[546] = 0.01005124659;
    nodeX[1638] = 0.05905888853;
    nodeX[1639] = 0.9402007994;
    nodeX[1640] = 0.3354616289;

    weight[547] = 0.01005124659;
    nodeX[1641] = -0.05905888853;
    nodeX[1642] = 0.9402007994;
    nodeX[1643] = 0.3354616289;

    weight[548] = 0.01005124659;
    nodeX[1644] = 0.05905888853;
    nodeX[1645] = -0.9402007994;
    nodeX[1646] = 0.3354616289;

    weight[549] = 0.01005124659;
    nodeX[1647] = 0.05905888853;
    nodeX[1648] = 0.9402007994;
    nodeX[1649] = -0.3354616289;

    weight[550] = 0.01005124659;
    nodeX[1650] = -0.05905888853;
    nodeX[1651] = -0.9402007994;
    nodeX[1652] = 0.3354616289;

    weight[551] = 0.01005124659;
    nodeX[1653] = 0.05905888853;
    nodeX[1654] = -0.9402007994;
    nodeX[1655] = -0.3354616289;

    weight[552] = 0.01005124659;
    nodeX[1656] = -0.05905888853;
    nodeX[1657] = 0.9402007994;
    nodeX[1658] = -0.3354616289;

    weight[553] = 0.01005124659;
    nodeX[1659] = -0.05905888853;
    nodeX[1660] = -0.9402007994;
    nodeX[1661] = -0.3354616289;

    weight[554] = 0.01005124659;
    nodeX[1662] = 0.05905888853;
    nodeX[1663] = 0.3354616289;
    nodeX[1664] = 0.9402007994;

    weight[555] = 0.01005124659;
    nodeX[1665] = -0.05905888853;
    nodeX[1666] = 0.3354616289;
    nodeX[1667] = 0.9402007994;

    weight[556] = 0.01005124659;
    nodeX[1668] = 0.05905888853;
    nodeX[1669] = -0.3354616289;
    nodeX[1670] = 0.9402007994;

    weight[557] = 0.01005124659;
    nodeX[1671] = 0.05905888853;
    nodeX[1672] = 0.3354616289;
    nodeX[1673] = -0.9402007994;

    weight[558] = 0.01005124659;
    nodeX[1674] = -0.05905888853;
    nodeX[1675] = -0.3354616289;
    nodeX[1676] = 0.9402007994;

    weight[559] = 0.01005124659;
    nodeX[1677] = 0.05905888853;
    nodeX[1678] = -0.3354616289;
    nodeX[1679] = -0.9402007994;

    weight[560] = 0.01005124659;
    nodeX[1680] = -0.05905888853;
    nodeX[1681] = 0.3354616289;
    nodeX[1682] = -0.9402007994;

    weight[561] = 0.01005124659;
    nodeX[1683] = -0.05905888853;
    nodeX[1684] = -0.3354616289;
    nodeX[1685] = -0.9402007994;

    weight[562] = 0.01005124659;
    nodeX[1686] = 0.9402007994;
    nodeX[1687] = 0.05905888853;
    nodeX[1688] = 0.3354616289;

    weight[563] = 0.01005124659;
    nodeX[1689] = -0.9402007994;
    nodeX[1690] = 0.05905888853;
    nodeX[1691] = 0.3354616289;

    weight[564] = 0.01005124659;
    nodeX[1692] = 0.9402007994;
    nodeX[1693] = -0.05905888853;
    nodeX[1694] = 0.3354616289;

    weight[565] = 0.01005124659;
    nodeX[1695] = 0.9402007994;
    nodeX[1696] = 0.05905888853;
    nodeX[1697] = -0.3354616289;

    weight[566] = 0.01005124659;
    nodeX[1698] = -0.9402007994;
    nodeX[1699] = -0.05905888853;
    nodeX[1700] = 0.3354616289;

    weight[567] = 0.01005124659;
    nodeX[1701] = 0.9402007994;
    nodeX[1702] = -0.05905888853;
    nodeX[1703] = -0.3354616289;

    weight[568] = 0.01005124659;
    nodeX[1704] = -0.9402007994;
    nodeX[1705] = 0.05905888853;
    nodeX[1706] = -0.3354616289;

    weight[569] = 0.01005124659;
    nodeX[1707] = -0.9402007994;
    nodeX[1708] = -0.05905888853;
    nodeX[1709] = -0.3354616289;

    weight[570] = 0.01005124659;
    nodeX[1710] = 0.3354616289;
    nodeX[1711] = 0.05905888853;
    nodeX[1712] = 0.9402007994;

    weight[571] = 0.01005124659;
    nodeX[1713] = -0.3354616289;
    nodeX[1714] = 0.05905888853;
    nodeX[1715] = 0.9402007994;

    weight[572] = 0.01005124659;
    nodeX[1716] = 0.3354616289;
    nodeX[1717] = -0.05905888853;
    nodeX[1718] = 0.9402007994;

    weight[573] = 0.01005124659;
    nodeX[1719] = 0.3354616289;
    nodeX[1720] = 0.05905888853;
    nodeX[1721] = -0.9402007994;

    weight[574] = 0.01005124659;
    nodeX[1722] = -0.3354616289;
    nodeX[1723] = -0.05905888853;
    nodeX[1724] = 0.9402007994;

    weight[575] = 0.01005124659;
    nodeX[1725] = 0.3354616289;
    nodeX[1726] = -0.05905888853;
    nodeX[1727] = -0.9402007994;

    weight[576] = 0.01005124659;
    nodeX[1728] = -0.3354616289;
    nodeX[1729] = 0.05905888853;
    nodeX[1730] = -0.9402007994;

    weight[577] = 0.01005124659;
    nodeX[1731] = -0.3354616289;
    nodeX[1732] = -0.05905888853;
    nodeX[1733] = -0.9402007994;

    weight[578] = 0.01018093606;
    nodeX[1734] = 0.932082204;
    nodeX[1735] = 0.3173615247;
    nodeX[1736] = 0.1746551678;

    weight[579] = 0.01018093606;
    nodeX[1737] = -0.932082204;
    nodeX[1738] = 0.3173615247;
    nodeX[1739] = 0.1746551678;

    weight[580] = 0.01018093606;
    nodeX[1740] = 0.932082204;
    nodeX[1741] = -0.3173615247;
    nodeX[1742] = 0.1746551678;

    weight[581] = 0.01018093606;
    nodeX[1743] = 0.932082204;
    nodeX[1744] = 0.3173615247;
    nodeX[1745] = -0.1746551678;

    weight[582] = 0.01018093606;
    nodeX[1746] = -0.932082204;
    nodeX[1747] = -0.3173615247;
    nodeX[1748] = 0.1746551678;

    weight[583] = 0.01018093606;
    nodeX[1749] = 0.932082204;
    nodeX[1750] = -0.3173615247;
    nodeX[1751] = -0.1746551678;

    weight[584] = 0.01018093606;
    nodeX[1752] = -0.932082204;
    nodeX[1753] = 0.3173615247;
    nodeX[1754] = -0.1746551678;

    weight[585] = 0.01018093606;
    nodeX[1755] = -0.932082204;
    nodeX[1756] = -0.3173615247;
    nodeX[1757] = -0.1746551678;

    weight[586] = 0.01018093606;
    nodeX[1758] = 0.3173615247;
    nodeX[1759] = 0.932082204;
    nodeX[1760] = 0.1746551678;

    weight[587] = 0.01018093606;
    nodeX[1761] = -0.3173615247;
    nodeX[1762] = 0.932082204;
    nodeX[1763] = 0.1746551678;

    weight[588] = 0.01018093606;
    nodeX[1764] = 0.3173615247;
    nodeX[1765] = -0.932082204;
    nodeX[1766] = 0.1746551678;

    weight[589] = 0.01018093606;
    nodeX[1767] = 0.3173615247;
    nodeX[1768] = 0.932082204;
    nodeX[1769] = -0.1746551678;

    weight[590] = 0.01018093606;
    nodeX[1770] = -0.3173615247;
    nodeX[1771] = -0.932082204;
    nodeX[1772] = 0.1746551678;

    weight[591] = 0.01018093606;
    nodeX[1773] = 0.3173615247;
    nodeX[1774] = -0.932082204;
    nodeX[1775] = -0.1746551678;

    weight[592] = 0.01018093606;
    nodeX[1776] = -0.3173615247;
    nodeX[1777] = 0.932082204;
    nodeX[1778] = -0.1746551678;

    weight[593] = 0.01018093606;
    nodeX[1779] = -0.3173615247;
    nodeX[1780] = -0.932082204;
    nodeX[1781] = -0.1746551678;

    weight[594] = 0.01018093606;
    nodeX[1782] = 0.1746551678;
    nodeX[1783] = 0.932082204;
    nodeX[1784] = 0.3173615247;

    weight[595] = 0.01018093606;
    nodeX[1785] = -0.1746551678;
    nodeX[1786] = 0.932082204;
    nodeX[1787] = 0.3173615247;

    weight[596] = 0.01018093606;
    nodeX[1788] = 0.1746551678;
    nodeX[1789] = -0.932082204;
    nodeX[1790] = 0.3173615247;

    weight[597] = 0.01018093606;
    nodeX[1791] = 0.1746551678;
    nodeX[1792] = 0.932082204;
    nodeX[1793] = -0.3173615247;

    weight[598] = 0.01018093606;
    nodeX[1794] = -0.1746551678;
    nodeX[1795] = -0.932082204;
    nodeX[1796] = 0.3173615247;

    weight[599] = 0.01018093606;
    nodeX[1797] = 0.1746551678;
    nodeX[1798] = -0.932082204;
    nodeX[1799] = -0.3173615247;

    weight[600] = 0.01018093606;
    nodeX[1800] = -0.1746551678;
    nodeX[1801] = 0.932082204;
    nodeX[1802] = -0.3173615247;

    weight[601] = 0.01018093606;
    nodeX[1803] = -0.1746551678;
    nodeX[1804] = -0.932082204;
    nodeX[1805] = -0.3173615247;

    weight[602] = 0.01018093606;
    nodeX[1806] = 0.1746551678;
    nodeX[1807] = 0.3173615247;
    nodeX[1808] = 0.932082204;

    weight[603] = 0.01018093606;
    nodeX[1809] = -0.1746551678;
    nodeX[1810] = 0.3173615247;
    nodeX[1811] = 0.932082204;

    weight[604] = 0.01018093606;
    nodeX[1812] = 0.1746551678;
    nodeX[1813] = -0.3173615247;
    nodeX[1814] = 0.932082204;

    weight[605] = 0.01018093606;
    nodeX[1815] = 0.1746551678;
    nodeX[1816] = 0.3173615247;
    nodeX[1817] = -0.932082204;

    weight[606] = 0.01018093606;
    nodeX[1818] = -0.1746551678;
    nodeX[1819] = -0.3173615247;
    nodeX[1820] = 0.932082204;

    weight[607] = 0.01018093606;
    nodeX[1821] = 0.1746551678;
    nodeX[1822] = -0.3173615247;
    nodeX[1823] = -0.932082204;

    weight[608] = 0.01018093606;
    nodeX[1824] = -0.1746551678;
    nodeX[1825] = 0.3173615247;
    nodeX[1826] = -0.932082204;

    weight[609] = 0.01018093606;
    nodeX[1827] = -0.1746551678;
    nodeX[1828] = -0.3173615247;
    nodeX[1829] = -0.932082204;

    weight[610] = 0.01018093606;
    nodeX[1830] = 0.932082204;
    nodeX[1831] = 0.1746551678;
    nodeX[1832] = 0.3173615247;

    weight[611] = 0.01018093606;
    nodeX[1833] = -0.932082204;
    nodeX[1834] = 0.1746551678;
    nodeX[1835] = 0.3173615247;

    weight[612] = 0.01018093606;
    nodeX[1836] = 0.932082204;
    nodeX[1837] = -0.1746551678;
    nodeX[1838] = 0.3173615247;

    weight[613] = 0.01018093606;
    nodeX[1839] = 0.932082204;
    nodeX[1840] = 0.1746551678;
    nodeX[1841] = -0.3173615247;

    weight[614] = 0.01018093606;
    nodeX[1842] = -0.932082204;
    nodeX[1843] = -0.1746551678;
    nodeX[1844] = 0.3173615247;

    weight[615] = 0.01018093606;
    nodeX[1845] = 0.932082204;
    nodeX[1846] = -0.1746551678;
    nodeX[1847] = -0.3173615247;

    weight[616] = 0.01018093606;
    nodeX[1848] = -0.932082204;
    nodeX[1849] = 0.1746551678;
    nodeX[1850] = -0.3173615247;

    weight[617] = 0.01018093606;
    nodeX[1851] = -0.932082204;
    nodeX[1852] = -0.1746551678;
    nodeX[1853] = -0.3173615247;

    weight[618] = 0.01018093606;
    nodeX[1854] = 0.3173615247;
    nodeX[1855] = 0.1746551678;
    nodeX[1856] = 0.932082204;

    weight[619] = 0.01018093606;
    nodeX[1857] = -0.3173615247;
    nodeX[1858] = 0.1746551678;
    nodeX[1859] = 0.932082204;

    weight[620] = 0.01018093606;
    nodeX[1860] = 0.3173615247;
    nodeX[1861] = -0.1746551678;
    nodeX[1862] = 0.932082204;

    weight[621] = 0.01018093606;
    nodeX[1863] = 0.3173615247;
    nodeX[1864] = 0.1746551678;
    nodeX[1865] = -0.932082204;

    weight[622] = 0.01018093606;
    nodeX[1866] = -0.3173615247;
    nodeX[1867] = -0.1746551678;
    nodeX[1868] = 0.932082204;

    weight[623] = 0.01018093606;
    nodeX[1869] = 0.3173615247;
    nodeX[1870] = -0.1746551678;
    nodeX[1871] = -0.932082204;

    weight[624] = 0.01018093606;
    nodeX[1872] = -0.3173615247;
    nodeX[1873] = 0.1746551678;
    nodeX[1874] = -0.932082204;

    weight[625] = 0.01018093606;
    nodeX[1875] = -0.3173615247;
    nodeX[1876] = -0.1746551678;
    nodeX[1877] = -0.932082204;

    weight[626] = 0.01066054175;
    nodeX[1878] = 0.9043674199;
    nodeX[1879] = 0.4090268427;
    nodeX[1880] = 0.1217235051;

    weight[627] = 0.01066054175;
    nodeX[1881] = -0.9043674199;
    nodeX[1882] = 0.4090268427;
    nodeX[1883] = 0.1217235051;

    weight[628] = 0.01066054175;
    nodeX[1884] = 0.9043674199;
    nodeX[1885] = -0.4090268427;
    nodeX[1886] = 0.1217235051;

    weight[629] = 0.01066054175;
    nodeX[1887] = 0.9043674199;
    nodeX[1888] = 0.4090268427;
    nodeX[1889] = -0.1217235051;

    weight[630] = 0.01066054175;
    nodeX[1890] = -0.9043674199;
    nodeX[1891] = -0.4090268427;
    nodeX[1892] = 0.1217235051;

    weight[631] = 0.01066054175;
    nodeX[1893] = 0.9043674199;
    nodeX[1894] = -0.4090268427;
    nodeX[1895] = -0.1217235051;

    weight[632] = 0.01066054175;
    nodeX[1896] = -0.9043674199;
    nodeX[1897] = 0.4090268427;
    nodeX[1898] = -0.1217235051;

    weight[633] = 0.01066054175;
    nodeX[1899] = -0.9043674199;
    nodeX[1900] = -0.4090268427;
    nodeX[1901] = -0.1217235051;

    weight[634] = 0.01066054175;
    nodeX[1902] = 0.4090268427;
    nodeX[1903] = 0.9043674199;
    nodeX[1904] = 0.1217235051;

    weight[635] = 0.01066054175;
    nodeX[1905] = -0.4090268427;
    nodeX[1906] = 0.9043674199;
    nodeX[1907] = 0.1217235051;

    weight[636] = 0.01066054175;
    nodeX[1908] = 0.4090268427;
    nodeX[1909] = -0.9043674199;
    nodeX[1910] = 0.1217235051;

    weight[637] = 0.01066054175;
    nodeX[1911] = 0.4090268427;
    nodeX[1912] = 0.9043674199;
    nodeX[1913] = -0.1217235051;

    weight[638] = 0.01066054175;
    nodeX[1914] = -0.4090268427;
    nodeX[1915] = -0.9043674199;
    nodeX[1916] = 0.1217235051;

    weight[639] = 0.01066054175;
    nodeX[1917] = 0.4090268427;
    nodeX[1918] = -0.9043674199;
    nodeX[1919] = -0.1217235051;

    weight[640] = 0.01066054175;
    nodeX[1920] = -0.4090268427;
    nodeX[1921] = 0.9043674199;
    nodeX[1922] = -0.1217235051;

    weight[641] = 0.01066054175;
    nodeX[1923] = -0.4090268427;
    nodeX[1924] = -0.9043674199;
    nodeX[1925] = -0.1217235051;

    weight[642] = 0.01066054175;
    nodeX[1926] = 0.1217235051;
    nodeX[1927] = 0.9043674199;
    nodeX[1928] = 0.4090268427;

    weight[643] = 0.01066054175;
    nodeX[1929] = -0.1217235051;
    nodeX[1930] = 0.9043674199;
    nodeX[1931] = 0.4090268427;

    weight[644] = 0.01066054175;
    nodeX[1932] = 0.1217235051;
    nodeX[1933] = -0.9043674199;
    nodeX[1934] = 0.4090268427;

    weight[645] = 0.01066054175;
    nodeX[1935] = 0.1217235051;
    nodeX[1936] = 0.9043674199;
    nodeX[1937] = -0.4090268427;

    weight[646] = 0.01066054175;
    nodeX[1938] = -0.1217235051;
    nodeX[1939] = -0.9043674199;
    nodeX[1940] = 0.4090268427;

    weight[647] = 0.01066054175;
    nodeX[1941] = 0.1217235051;
    nodeX[1942] = -0.9043674199;
    nodeX[1943] = -0.4090268427;

    weight[648] = 0.01066054175;
    nodeX[1944] = -0.1217235051;
    nodeX[1945] = 0.9043674199;
    nodeX[1946] = -0.4090268427;

    weight[649] = 0.01066054175;
    nodeX[1947] = -0.1217235051;
    nodeX[1948] = -0.9043674199;
    nodeX[1949] = -0.4090268427;

    weight[650] = 0.01066054175;
    nodeX[1950] = 0.1217235051;
    nodeX[1951] = 0.4090268427;
    nodeX[1952] = 0.9043674199;

    weight[651] = 0.01066054175;
    nodeX[1953] = -0.1217235051;
    nodeX[1954] = 0.4090268427;
    nodeX[1955] = 0.9043674199;

    weight[652] = 0.01066054175;
    nodeX[1956] = 0.1217235051;
    nodeX[1957] = -0.4090268427;
    nodeX[1958] = 0.9043674199;

    weight[653] = 0.01066054175;
    nodeX[1959] = 0.1217235051;
    nodeX[1960] = 0.4090268427;
    nodeX[1961] = -0.9043674199;

    weight[654] = 0.01066054175;
    nodeX[1962] = -0.1217235051;
    nodeX[1963] = -0.4090268427;
    nodeX[1964] = 0.9043674199;

    weight[655] = 0.01066054175;
    nodeX[1965] = 0.1217235051;
    nodeX[1966] = -0.4090268427;
    nodeX[1967] = -0.9043674199;

    weight[656] = 0.01066054175;
    nodeX[1968] = -0.1217235051;
    nodeX[1969] = 0.4090268427;
    nodeX[1970] = -0.9043674199;

    weight[657] = 0.01066054175;
    nodeX[1971] = -0.1217235051;
    nodeX[1972] = -0.4090268427;
    nodeX[1973] = -0.9043674199;

    weight[658] = 0.01066054175;
    nodeX[1974] = 0.9043674199;
    nodeX[1975] = 0.1217235051;
    nodeX[1976] = 0.4090268427;

    weight[659] = 0.01066054175;
    nodeX[1977] = -0.9043674199;
    nodeX[1978] = 0.1217235051;
    nodeX[1979] = 0.4090268427;

    weight[660] = 0.01066054175;
    nodeX[1980] = 0.9043674199;
    nodeX[1981] = -0.1217235051;
    nodeX[1982] = 0.4090268427;

    weight[661] = 0.01066054175;
    nodeX[1983] = 0.9043674199;
    nodeX[1984] = 0.1217235051;
    nodeX[1985] = -0.4090268427;

    weight[662] = 0.01066054175;
    nodeX[1986] = -0.9043674199;
    nodeX[1987] = -0.1217235051;
    nodeX[1988] = 0.4090268427;

    weight[663] = 0.01066054175;
    nodeX[1989] = 0.9043674199;
    nodeX[1990] = -0.1217235051;
    nodeX[1991] = -0.4090268427;

    weight[664] = 0.01066054175;
    nodeX[1992] = -0.9043674199;
    nodeX[1993] = 0.1217235051;
    nodeX[1994] = -0.4090268427;

    weight[665] = 0.01066054175;
    nodeX[1995] = -0.9043674199;
    nodeX[1996] = -0.1217235051;
    nodeX[1997] = -0.4090268427;

    weight[666] = 0.01066054175;
    nodeX[1998] = 0.4090268427;
    nodeX[1999] = 0.1217235051;
    nodeX[2000] = 0.9043674199;

    weight[667] = 0.01066054175;
    nodeX[2001] = -0.4090268427;
    nodeX[2002] = 0.1217235051;
    nodeX[2003] = 0.9043674199;

    weight[668] = 0.01066054175;
    nodeX[2004] = 0.4090268427;
    nodeX[2005] = -0.1217235051;
    nodeX[2006] = 0.9043674199;

    weight[669] = 0.01066054175;
    nodeX[2007] = 0.4090268427;
    nodeX[2008] = 0.1217235051;
    nodeX[2009] = -0.9043674199;

    weight[670] = 0.01066054175;
    nodeX[2010] = -0.4090268427;
    nodeX[2011] = -0.1217235051;
    nodeX[2012] = 0.9043674199;

    weight[671] = 0.01066054175;
    nodeX[2013] = 0.4090268427;
    nodeX[2014] = -0.1217235051;
    nodeX[2015] = -0.9043674199;

    weight[672] = 0.01066054175;
    nodeX[2016] = -0.4090268427;
    nodeX[2017] = 0.1217235051;
    nodeX[2018] = -0.9043674199;

    weight[673] = 0.01066054175;
    nodeX[2019] = -0.4090268427;
    nodeX[2020] = -0.1217235051;
    nodeX[2021] = -0.9043674199;

    weight[674] = 0.01075216276;
    nodeX[2022] = 0.891240756;
    nodeX[2023] = 0.3854291151;
    nodeX[2024] = 0.2390278479;

    weight[675] = 0.01075216276;
    nodeX[2025] = -0.891240756;
    nodeX[2026] = 0.3854291151;
    nodeX[2027] = 0.2390278479;

    weight[676] = 0.01075216276;
    nodeX[2028] = 0.891240756;
    nodeX[2029] = -0.3854291151;
    nodeX[2030] = 0.2390278479;

    weight[677] = 0.01075216276;
    nodeX[2031] = 0.891240756;
    nodeX[2032] = 0.3854291151;
    nodeX[2033] = -0.2390278479;

    weight[678] = 0.01075216276;
    nodeX[2034] = -0.891240756;
    nodeX[2035] = -0.3854291151;
    nodeX[2036] = 0.2390278479;

    weight[679] = 0.01075216276;
    nodeX[2037] = 0.891240756;
    nodeX[2038] = -0.3854291151;
    nodeX[2039] = -0.2390278479;

    weight[680] = 0.01075216276;
    nodeX[2040] = -0.891240756;
    nodeX[2041] = 0.3854291151;
    nodeX[2042] = -0.2390278479;

    weight[681] = 0.01075216276;
    nodeX[2043] = -0.891240756;
    nodeX[2044] = -0.3854291151;
    nodeX[2045] = -0.2390278479;

    weight[682] = 0.01075216276;
    nodeX[2046] = 0.3854291151;
    nodeX[2047] = 0.891240756;
    nodeX[2048] = 0.2390278479;

    weight[683] = 0.01075216276;
    nodeX[2049] = -0.3854291151;
    nodeX[2050] = 0.891240756;
    nodeX[2051] = 0.2390278479;

    weight[684] = 0.01075216276;
    nodeX[2052] = 0.3854291151;
    nodeX[2053] = -0.891240756;
    nodeX[2054] = 0.2390278479;

    weight[685] = 0.01075216276;
    nodeX[2055] = 0.3854291151;
    nodeX[2056] = 0.891240756;
    nodeX[2057] = -0.2390278479;

    weight[686] = 0.01075216276;
    nodeX[2058] = -0.3854291151;
    nodeX[2059] = -0.891240756;
    nodeX[2060] = 0.2390278479;

    weight[687] = 0.01075216276;
    nodeX[2061] = 0.3854291151;
    nodeX[2062] = -0.891240756;
    nodeX[2063] = -0.2390278479;

    weight[688] = 0.01075216276;
    nodeX[2064] = -0.3854291151;
    nodeX[2065] = 0.891240756;
    nodeX[2066] = -0.2390278479;

    weight[689] = 0.01075216276;
    nodeX[2067] = -0.3854291151;
    nodeX[2068] = -0.891240756;
    nodeX[2069] = -0.2390278479;

    weight[690] = 0.01075216276;
    nodeX[2070] = 0.2390278479;
    nodeX[2071] = 0.891240756;
    nodeX[2072] = 0.3854291151;

    weight[691] = 0.01075216276;
    nodeX[2073] = -0.2390278479;
    nodeX[2074] = 0.891240756;
    nodeX[2075] = 0.3854291151;

    weight[692] = 0.01075216276;
    nodeX[2076] = 0.2390278479;
    nodeX[2077] = -0.891240756;
    nodeX[2078] = 0.3854291151;

    weight[693] = 0.01075216276;
    nodeX[2079] = 0.2390278479;
    nodeX[2080] = 0.891240756;
    nodeX[2081] = -0.3854291151;

    weight[694] = 0.01075216276;
    nodeX[2082] = -0.2390278479;
    nodeX[2083] = -0.891240756;
    nodeX[2084] = 0.3854291151;

    weight[695] = 0.01075216276;
    nodeX[2085] = 0.2390278479;
    nodeX[2086] = -0.891240756;
    nodeX[2087] = -0.3854291151;

    weight[696] = 0.01075216276;
    nodeX[2088] = -0.2390278479;
    nodeX[2089] = 0.891240756;
    nodeX[2090] = -0.3854291151;

    weight[697] = 0.01075216276;
    nodeX[2091] = -0.2390278479;
    nodeX[2092] = -0.891240756;
    nodeX[2093] = -0.3854291151;

    weight[698] = 0.01075216276;
    nodeX[2094] = 0.2390278479;
    nodeX[2095] = 0.3854291151;
    nodeX[2096] = 0.891240756;

    weight[699] = 0.01075216276;
    nodeX[2097] = -0.2390278479;
    nodeX[2098] = 0.3854291151;
    nodeX[2099] = 0.891240756;

    weight[700] = 0.01075216276;
    nodeX[2100] = 0.2390278479;
    nodeX[2101] = -0.3854291151;
    nodeX[2102] = 0.891240756;

    weight[701] = 0.01075216276;
    nodeX[2103] = 0.2390278479;
    nodeX[2104] = 0.3854291151;
    nodeX[2105] = -0.891240756;

    weight[702] = 0.01075216276;
    nodeX[2106] = -0.2390278479;
    nodeX[2107] = -0.3854291151;
    nodeX[2108] = 0.891240756;

    weight[703] = 0.01075216276;
    nodeX[2109] = 0.2390278479;
    nodeX[2110] = -0.3854291151;
    nodeX[2111] = -0.891240756;

    weight[704] = 0.01075216276;
    nodeX[2112] = -0.2390278479;
    nodeX[2113] = 0.3854291151;
    nodeX[2114] = -0.891240756;

    weight[705] = 0.01075216276;
    nodeX[2115] = -0.2390278479;
    nodeX[2116] = -0.3854291151;
    nodeX[2117] = -0.891240756;

    weight[706] = 0.01075216276;
    nodeX[2118] = 0.891240756;
    nodeX[2119] = 0.2390278479;
    nodeX[2120] = 0.3854291151;

    weight[707] = 0.01075216276;
    nodeX[2121] = -0.891240756;
    nodeX[2122] = 0.2390278479;
    nodeX[2123] = 0.3854291151;

    weight[708] = 0.01075216276;
    nodeX[2124] = 0.891240756;
    nodeX[2125] = -0.2390278479;
    nodeX[2126] = 0.3854291151;

    weight[709] = 0.01075216276;
    nodeX[2127] = 0.891240756;
    nodeX[2128] = 0.2390278479;
    nodeX[2129] = -0.3854291151;

    weight[710] = 0.01075216276;
    nodeX[2130] = -0.891240756;
    nodeX[2131] = -0.2390278479;
    nodeX[2132] = 0.3854291151;

    weight[711] = 0.01075216276;
    nodeX[2133] = 0.891240756;
    nodeX[2134] = -0.2390278479;
    nodeX[2135] = -0.3854291151;

    weight[712] = 0.01075216276;
    nodeX[2136] = -0.891240756;
    nodeX[2137] = 0.2390278479;
    nodeX[2138] = -0.3854291151;

    weight[713] = 0.01075216276;
    nodeX[2139] = -0.891240756;
    nodeX[2140] = -0.2390278479;
    nodeX[2141] = -0.3854291151;

    weight[714] = 0.01075216276;
    nodeX[2142] = 0.3854291151;
    nodeX[2143] = 0.2390278479;
    nodeX[2144] = 0.891240756;

    weight[715] = 0.01075216276;
    nodeX[2145] = -0.3854291151;
    nodeX[2146] = 0.2390278479;
    nodeX[2147] = 0.891240756;

    weight[716] = 0.01075216276;
    nodeX[2148] = 0.3854291151;
    nodeX[2149] = -0.2390278479;
    nodeX[2150] = 0.891240756;

    weight[717] = 0.01075216276;
    nodeX[2151] = 0.3854291151;
    nodeX[2152] = 0.2390278479;
    nodeX[2153] = -0.891240756;

    weight[718] = 0.01075216276;
    nodeX[2154] = -0.3854291151;
    nodeX[2155] = -0.2390278479;
    nodeX[2156] = 0.891240756;

    weight[719] = 0.01075216276;
    nodeX[2157] = 0.3854291151;
    nodeX[2158] = -0.2390278479;
    nodeX[2159] = -0.891240756;

    weight[720] = 0.01075216276;
    nodeX[2160] = -0.3854291151;
    nodeX[2161] = 0.2390278479;
    nodeX[2162] = -0.891240756;

    weight[721] = 0.01075216276;
    nodeX[2163] = -0.3854291151;
    nodeX[2164] = -0.2390278479;
    nodeX[2165] = -0.891240756;

    weight[722] = 0.01106243829;
    nodeX[2166] = 0.8676435628;
    nodeX[2167] = 0.4932221185;
    nodeX[2168] = 0.06266250624;

    weight[723] = 0.01106243829;
    nodeX[2169] = -0.8676435628;
    nodeX[2170] = 0.4932221185;
    nodeX[2171] = 0.06266250624;

    weight[724] = 0.01106243829;
    nodeX[2172] = 0.8676435628;
    nodeX[2173] = -0.4932221185;
    nodeX[2174] = 0.06266250624;

    weight[725] = 0.01106243829;
    nodeX[2175] = 0.8676435628;
    nodeX[2176] = 0.4932221185;
    nodeX[2177] = -0.06266250624;

    weight[726] = 0.01106243829;
    nodeX[2178] = -0.8676435628;
    nodeX[2179] = -0.4932221185;
    nodeX[2180] = 0.06266250624;

    weight[727] = 0.01106243829;
    nodeX[2181] = 0.8676435628;
    nodeX[2182] = -0.4932221185;
    nodeX[2183] = -0.06266250624;

    weight[728] = 0.01106243829;
    nodeX[2184] = -0.8676435628;
    nodeX[2185] = 0.4932221185;
    nodeX[2186] = -0.06266250624;

    weight[729] = 0.01106243829;
    nodeX[2187] = -0.8676435628;
    nodeX[2188] = -0.4932221185;
    nodeX[2189] = -0.06266250624;

    weight[730] = 0.01106243829;
    nodeX[2190] = 0.4932221185;
    nodeX[2191] = 0.8676435628;
    nodeX[2192] = 0.06266250624;

    weight[731] = 0.01106243829;
    nodeX[2193] = -0.4932221185;
    nodeX[2194] = 0.8676435628;
    nodeX[2195] = 0.06266250624;

    weight[732] = 0.01106243829;
    nodeX[2196] = 0.4932221185;
    nodeX[2197] = -0.8676435628;
    nodeX[2198] = 0.06266250624;

    weight[733] = 0.01106243829;
    nodeX[2199] = 0.4932221185;
    nodeX[2200] = 0.8676435628;
    nodeX[2201] = -0.06266250624;

    weight[734] = 0.01106243829;
    nodeX[2202] = -0.4932221185;
    nodeX[2203] = -0.8676435628;
    nodeX[2204] = 0.06266250624;

    weight[735] = 0.01106243829;
    nodeX[2205] = 0.4932221185;
    nodeX[2206] = -0.8676435628;
    nodeX[2207] = -0.06266250624;

    weight[736] = 0.01106243829;
    nodeX[2208] = -0.4932221185;
    nodeX[2209] = 0.8676435628;
    nodeX[2210] = -0.06266250624;

    weight[737] = 0.01106243829;
    nodeX[2211] = -0.4932221185;
    nodeX[2212] = -0.8676435628;
    nodeX[2213] = -0.06266250624;

    weight[738] = 0.01106243829;
    nodeX[2214] = 0.06266250624;
    nodeX[2215] = 0.8676435628;
    nodeX[2216] = 0.4932221185;

    weight[739] = 0.01106243829;
    nodeX[2217] = -0.06266250624;
    nodeX[2218] = 0.8676435628;
    nodeX[2219] = 0.4932221185;

    weight[740] = 0.01106243829;
    nodeX[2220] = 0.06266250624;
    nodeX[2221] = -0.8676435628;
    nodeX[2222] = 0.4932221185;

    weight[741] = 0.01106243829;
    nodeX[2223] = 0.06266250624;
    nodeX[2224] = 0.8676435628;
    nodeX[2225] = -0.4932221185;

    weight[742] = 0.01106243829;
    nodeX[2226] = -0.06266250624;
    nodeX[2227] = -0.8676435628;
    nodeX[2228] = 0.4932221185;

    weight[743] = 0.01106243829;
    nodeX[2229] = 0.06266250624;
    nodeX[2230] = -0.8676435628;
    nodeX[2231] = -0.4932221185;

    weight[744] = 0.01106243829;
    nodeX[2232] = -0.06266250624;
    nodeX[2233] = 0.8676435628;
    nodeX[2234] = -0.4932221185;

    weight[745] = 0.01106243829;
    nodeX[2235] = -0.06266250624;
    nodeX[2236] = -0.8676435628;
    nodeX[2237] = -0.4932221185;

    weight[746] = 0.01106243829;
    nodeX[2238] = 0.06266250624;
    nodeX[2239] = 0.4932221185;
    nodeX[2240] = 0.8676435628;

    weight[747] = 0.01106243829;
    nodeX[2241] = -0.06266250624;
    nodeX[2242] = 0.4932221185;
    nodeX[2243] = 0.8676435628;

    weight[748] = 0.01106243829;
    nodeX[2244] = 0.06266250624;
    nodeX[2245] = -0.4932221185;
    nodeX[2246] = 0.8676435628;

    weight[749] = 0.01106243829;
    nodeX[2247] = 0.06266250624;
    nodeX[2248] = 0.4932221185;
    nodeX[2249] = -0.8676435628;

    weight[750] = 0.01106243829;
    nodeX[2250] = -0.06266250624;
    nodeX[2251] = -0.4932221185;
    nodeX[2252] = 0.8676435628;

    weight[751] = 0.01106243829;
    nodeX[2253] = 0.06266250624;
    nodeX[2254] = -0.4932221185;
    nodeX[2255] = -0.8676435628;

    weight[752] = 0.01106243829;
    nodeX[2256] = -0.06266250624;
    nodeX[2257] = 0.4932221185;
    nodeX[2258] = -0.8676435628;

    weight[753] = 0.01106243829;
    nodeX[2259] = -0.06266250624;
    nodeX[2260] = -0.4932221185;
    nodeX[2261] = -0.8676435628;

    weight[754] = 0.01106243829;
    nodeX[2262] = 0.8676435628;
    nodeX[2263] = 0.06266250624;
    nodeX[2264] = 0.4932221185;

    weight[755] = 0.01106243829;
    nodeX[2265] = -0.8676435628;
    nodeX[2266] = 0.06266250624;
    nodeX[2267] = 0.4932221185;

    weight[756] = 0.01106243829;
    nodeX[2268] = 0.8676435628;
    nodeX[2269] = -0.06266250624;
    nodeX[2270] = 0.4932221185;

    weight[757] = 0.01106243829;
    nodeX[2271] = 0.8676435628;
    nodeX[2272] = 0.06266250624;
    nodeX[2273] = -0.4932221185;

    weight[758] = 0.01106243829;
    nodeX[2274] = -0.8676435628;
    nodeX[2275] = -0.06266250624;
    nodeX[2276] = 0.4932221185;

    weight[759] = 0.01106243829;
    nodeX[2277] = 0.8676435628;
    nodeX[2278] = -0.06266250624;
    nodeX[2279] = -0.4932221185;

    weight[760] = 0.01106243829;
    nodeX[2280] = -0.8676435628;
    nodeX[2281] = 0.06266250624;
    nodeX[2282] = -0.4932221185;

    weight[761] = 0.01106243829;
    nodeX[2283] = -0.8676435628;
    nodeX[2284] = -0.06266250624;
    nodeX[2285] = -0.4932221185;

    weight[762] = 0.01106243829;
    nodeX[2286] = 0.4932221185;
    nodeX[2287] = 0.06266250624;
    nodeX[2288] = 0.8676435628;

    weight[763] = 0.01106243829;
    nodeX[2289] = -0.4932221185;
    nodeX[2290] = 0.06266250624;
    nodeX[2291] = 0.8676435628;

    weight[764] = 0.01106243829;
    nodeX[2292] = 0.4932221185;
    nodeX[2293] = -0.06266250624;
    nodeX[2294] = 0.8676435628;

    weight[765] = 0.01106243829;
    nodeX[2295] = 0.4932221185;
    nodeX[2296] = 0.06266250624;
    nodeX[2297] = -0.8676435628;

    weight[766] = 0.01106243829;
    nodeX[2298] = -0.4932221185;
    nodeX[2299] = -0.06266250624;
    nodeX[2300] = 0.8676435628;

    weight[767] = 0.01106243829;
    nodeX[2301] = 0.4932221185;
    nodeX[2302] = -0.06266250624;
    nodeX[2303] = -0.8676435628;

    weight[768] = 0.01106243829;
    nodeX[2304] = -0.4932221185;
    nodeX[2305] = 0.06266250624;
    nodeX[2306] = -0.8676435628;

    weight[769] = 0.01106243829;
    nodeX[2307] = -0.4932221185;
    nodeX[2308] = -0.06266250624;
    nodeX[2309] = -0.8676435628;

    weight[770] = 0.0110722897;
    nodeX[2310] = 0.8581979986;
    nodeX[2311] = 0.4785320676;
    nodeX[2312] = 0.1857505195;

    weight[771] = 0.0110722897;
    nodeX[2313] = -0.8581979986;
    nodeX[2314] = 0.4785320676;
    nodeX[2315] = 0.1857505195;

    weight[772] = 0.0110722897;
    nodeX[2316] = 0.8581979986;
    nodeX[2317] = -0.4785320676;
    nodeX[2318] = 0.1857505195;

    weight[773] = 0.0110722897;
    nodeX[2319] = 0.8581979986;
    nodeX[2320] = 0.4785320676;
    nodeX[2321] = -0.1857505195;

    weight[774] = 0.0110722897;
    nodeX[2322] = -0.8581979986;
    nodeX[2323] = -0.4785320676;
    nodeX[2324] = 0.1857505195;

    weight[775] = 0.0110722897;
    nodeX[2325] = 0.8581979986;
    nodeX[2326] = -0.4785320676;
    nodeX[2327] = -0.1857505195;

    weight[776] = 0.0110722897;
    nodeX[2328] = -0.8581979986;
    nodeX[2329] = 0.4785320676;
    nodeX[2330] = -0.1857505195;

    weight[777] = 0.0110722897;
    nodeX[2331] = -0.8581979986;
    nodeX[2332] = -0.4785320676;
    nodeX[2333] = -0.1857505195;

    weight[778] = 0.0110722897;
    nodeX[2334] = 0.4785320676;
    nodeX[2335] = 0.8581979986;
    nodeX[2336] = 0.1857505195;

    weight[779] = 0.0110722897;
    nodeX[2337] = -0.4785320676;
    nodeX[2338] = 0.8581979986;
    nodeX[2339] = 0.1857505195;

    weight[780] = 0.0110722897;
    nodeX[2340] = 0.4785320676;
    nodeX[2341] = -0.8581979986;
    nodeX[2342] = 0.1857505195;

    weight[781] = 0.0110722897;
    nodeX[2343] = 0.4785320676;
    nodeX[2344] = 0.8581979986;
    nodeX[2345] = -0.1857505195;

    weight[782] = 0.0110722897;
    nodeX[2346] = -0.4785320676;
    nodeX[2347] = -0.8581979986;
    nodeX[2348] = 0.1857505195;

    weight[783] = 0.0110722897;
    nodeX[2349] = 0.4785320676;
    nodeX[2350] = -0.8581979986;
    nodeX[2351] = -0.1857505195;

    weight[784] = 0.0110722897;
    nodeX[2352] = -0.4785320676;
    nodeX[2353] = 0.8581979986;
    nodeX[2354] = -0.1857505195;

    weight[785] = 0.0110722897;
    nodeX[2355] = -0.4785320676;
    nodeX[2356] = -0.8581979986;
    nodeX[2357] = -0.1857505195;

    weight[786] = 0.0110722897;
    nodeX[2358] = 0.1857505195;
    nodeX[2359] = 0.8581979986;
    nodeX[2360] = 0.4785320676;

    weight[787] = 0.0110722897;
    nodeX[2361] = -0.1857505195;
    nodeX[2362] = 0.8581979986;
    nodeX[2363] = 0.4785320676;

    weight[788] = 0.0110722897;
    nodeX[2364] = 0.1857505195;
    nodeX[2365] = -0.8581979986;
    nodeX[2366] = 0.4785320676;

    weight[789] = 0.0110722897;
    nodeX[2367] = 0.1857505195;
    nodeX[2368] = 0.8581979986;
    nodeX[2369] = -0.4785320676;

    weight[790] = 0.0110722897;
    nodeX[2370] = -0.1857505195;
    nodeX[2371] = -0.8581979986;
    nodeX[2372] = 0.4785320676;

    weight[791] = 0.0110722897;
    nodeX[2373] = 0.1857505195;
    nodeX[2374] = -0.8581979986;
    nodeX[2375] = -0.4785320676;

    weight[792] = 0.0110722897;
    nodeX[2376] = -0.1857505195;
    nodeX[2377] = 0.8581979986;
    nodeX[2378] = -0.4785320676;

    weight[793] = 0.0110722897;
    nodeX[2379] = -0.1857505195;
    nodeX[2380] = -0.8581979986;
    nodeX[2381] = -0.4785320676;

    weight[794] = 0.0110722897;
    nodeX[2382] = 0.1857505195;
    nodeX[2383] = 0.4785320676;
    nodeX[2384] = 0.8581979986;

    weight[795] = 0.0110722897;
    nodeX[2385] = -0.1857505195;
    nodeX[2386] = 0.4785320676;
    nodeX[2387] = 0.8581979986;

    weight[796] = 0.0110722897;
    nodeX[2388] = 0.1857505195;
    nodeX[2389] = -0.4785320676;
    nodeX[2390] = 0.8581979986;

    weight[797] = 0.0110722897;
    nodeX[2391] = 0.1857505195;
    nodeX[2392] = 0.4785320676;
    nodeX[2393] = -0.8581979986;

    weight[798] = 0.0110722897;
    nodeX[2394] = -0.1857505195;
    nodeX[2395] = -0.4785320676;
    nodeX[2396] = 0.8581979986;

    weight[799] = 0.0110722897;
    nodeX[2397] = 0.1857505195;
    nodeX[2398] = -0.4785320676;
    nodeX[2399] = -0.8581979986;

    weight[800] = 0.0110722897;
    nodeX[2400] = -0.1857505195;
    nodeX[2401] = 0.4785320676;
    nodeX[2402] = -0.8581979986;

    weight[801] = 0.0110722897;
    nodeX[2403] = -0.1857505195;
    nodeX[2404] = -0.4785320676;
    nodeX[2405] = -0.8581979986;

    weight[802] = 0.0110722897;
    nodeX[2406] = 0.8581979986;
    nodeX[2407] = 0.1857505195;
    nodeX[2408] = 0.4785320676;

    weight[803] = 0.0110722897;
    nodeX[2409] = -0.8581979986;
    nodeX[2410] = 0.1857505195;
    nodeX[2411] = 0.4785320676;

    weight[804] = 0.0110722897;
    nodeX[2412] = 0.8581979986;
    nodeX[2413] = -0.1857505195;
    nodeX[2414] = 0.4785320676;

    weight[805] = 0.0110722897;
    nodeX[2415] = 0.8581979986;
    nodeX[2416] = 0.1857505195;
    nodeX[2417] = -0.4785320676;

    weight[806] = 0.0110722897;
    nodeX[2418] = -0.8581979986;
    nodeX[2419] = -0.1857505195;
    nodeX[2420] = 0.4785320676;

    weight[807] = 0.0110722897;
    nodeX[2421] = 0.8581979986;
    nodeX[2422] = -0.1857505195;
    nodeX[2423] = -0.4785320676;

    weight[808] = 0.0110722897;
    nodeX[2424] = -0.8581979986;
    nodeX[2425] = 0.1857505195;
    nodeX[2426] = -0.4785320676;

    weight[809] = 0.0110722897;
    nodeX[2427] = -0.8581979986;
    nodeX[2428] = -0.1857505195;
    nodeX[2429] = -0.4785320676;

    weight[810] = 0.0110722897;
    nodeX[2430] = 0.4785320676;
    nodeX[2431] = 0.1857505195;
    nodeX[2432] = 0.8581979986;

    weight[811] = 0.0110722897;
    nodeX[2433] = -0.4785320676;
    nodeX[2434] = 0.1857505195;
    nodeX[2435] = 0.8581979986;

    weight[812] = 0.0110722897;
    nodeX[2436] = 0.4785320676;
    nodeX[2437] = -0.1857505195;
    nodeX[2438] = 0.8581979986;

    weight[813] = 0.0110722897;
    nodeX[2439] = 0.4785320676;
    nodeX[2440] = 0.1857505195;
    nodeX[2441] = -0.8581979986;

    weight[814] = 0.0110722897;
    nodeX[2442] = -0.4785320676;
    nodeX[2443] = -0.1857505195;
    nodeX[2444] = 0.8581979986;

    weight[815] = 0.0110722897;
    nodeX[2445] = 0.4785320676;
    nodeX[2446] = -0.1857505195;
    nodeX[2447] = -0.8581979986;

    weight[816] = 0.0110722897;
    nodeX[2448] = -0.4785320676;
    nodeX[2449] = 0.1857505195;
    nodeX[2450] = -0.8581979986;

    weight[817] = 0.0110722897;
    nodeX[2451] = -0.4785320676;
    nodeX[2452] = -0.1857505195;
    nodeX[2453] = -0.8581979986;

    weight[818] = 0.01112159279;
    nodeX[2454] = 0.8396753624;
    nodeX[2455] = 0.4507422593;
    nodeX[2456] = 0.3029466974;

    weight[819] = 0.01112159279;
    nodeX[2457] = -0.8396753624;
    nodeX[2458] = 0.4507422593;
    nodeX[2459] = 0.3029466974;

    weight[820] = 0.01112159279;
    nodeX[2460] = 0.8396753624;
    nodeX[2461] = -0.4507422593;
    nodeX[2462] = 0.3029466974;

    weight[821] = 0.01112159279;
    nodeX[2463] = 0.8396753624;
    nodeX[2464] = 0.4507422593;
    nodeX[2465] = -0.3029466974;

    weight[822] = 0.01112159279;
    nodeX[2466] = -0.8396753624;
    nodeX[2467] = -0.4507422593;
    nodeX[2468] = 0.3029466974;

    weight[823] = 0.01112159279;
    nodeX[2469] = 0.8396753624;
    nodeX[2470] = -0.4507422593;
    nodeX[2471] = -0.3029466974;

    weight[824] = 0.01112159279;
    nodeX[2472] = -0.8396753624;
    nodeX[2473] = 0.4507422593;
    nodeX[2474] = -0.3029466974;

    weight[825] = 0.01112159279;
    nodeX[2475] = -0.8396753624;
    nodeX[2476] = -0.4507422593;
    nodeX[2477] = -0.3029466974;

    weight[826] = 0.01112159279;
    nodeX[2478] = 0.4507422593;
    nodeX[2479] = 0.8396753624;
    nodeX[2480] = 0.3029466974;

    weight[827] = 0.01112159279;
    nodeX[2481] = -0.4507422593;
    nodeX[2482] = 0.8396753624;
    nodeX[2483] = 0.3029466974;

    weight[828] = 0.01112159279;
    nodeX[2484] = 0.4507422593;
    nodeX[2485] = -0.8396753624;
    nodeX[2486] = 0.3029466974;

    weight[829] = 0.01112159279;
    nodeX[2487] = 0.4507422593;
    nodeX[2488] = 0.8396753624;
    nodeX[2489] = -0.3029466974;

    weight[830] = 0.01112159279;
    nodeX[2490] = -0.4507422593;
    nodeX[2491] = -0.8396753624;
    nodeX[2492] = 0.3029466974;

    weight[831] = 0.01112159279;
    nodeX[2493] = 0.4507422593;
    nodeX[2494] = -0.8396753624;
    nodeX[2495] = -0.3029466974;

    weight[832] = 0.01112159279;
    nodeX[2496] = -0.4507422593;
    nodeX[2497] = 0.8396753624;
    nodeX[2498] = -0.3029466974;

    weight[833] = 0.01112159279;
    nodeX[2499] = -0.4507422593;
    nodeX[2500] = -0.8396753624;
    nodeX[2501] = -0.3029466974;

    weight[834] = 0.01112159279;
    nodeX[2502] = 0.3029466974;
    nodeX[2503] = 0.8396753624;
    nodeX[2504] = 0.4507422593;

    weight[835] = 0.01112159279;
    nodeX[2505] = -0.3029466974;
    nodeX[2506] = 0.8396753624;
    nodeX[2507] = 0.4507422593;

    weight[836] = 0.01112159279;
    nodeX[2508] = 0.3029466974;
    nodeX[2509] = -0.8396753624;
    nodeX[2510] = 0.4507422593;

    weight[837] = 0.01112159279;
    nodeX[2511] = 0.3029466974;
    nodeX[2512] = 0.8396753624;
    nodeX[2513] = -0.4507422593;

    weight[838] = 0.01112159279;
    nodeX[2514] = -0.3029466974;
    nodeX[2515] = -0.8396753624;
    nodeX[2516] = 0.4507422593;

    weight[839] = 0.01112159279;
    nodeX[2517] = 0.3029466974;
    nodeX[2518] = -0.8396753624;
    nodeX[2519] = -0.4507422593;

    weight[840] = 0.01112159279;
    nodeX[2520] = -0.3029466974;
    nodeX[2521] = 0.8396753624;
    nodeX[2522] = -0.4507422593;

    weight[841] = 0.01112159279;
    nodeX[2523] = -0.3029466974;
    nodeX[2524] = -0.8396753624;
    nodeX[2525] = -0.4507422593;

    weight[842] = 0.01112159279;
    nodeX[2526] = 0.3029466974;
    nodeX[2527] = 0.4507422593;
    nodeX[2528] = 0.8396753624;

    weight[843] = 0.01112159279;
    nodeX[2529] = -0.3029466974;
    nodeX[2530] = 0.4507422593;
    nodeX[2531] = 0.8396753624;

    weight[844] = 0.01112159279;
    nodeX[2532] = 0.3029466974;
    nodeX[2533] = -0.4507422593;
    nodeX[2534] = 0.8396753624;

    weight[845] = 0.01112159279;
    nodeX[2535] = 0.3029466974;
    nodeX[2536] = 0.4507422593;
    nodeX[2537] = -0.8396753624;

    weight[846] = 0.01112159279;
    nodeX[2538] = -0.3029466974;
    nodeX[2539] = -0.4507422593;
    nodeX[2540] = 0.8396753624;

    weight[847] = 0.01112159279;
    nodeX[2541] = 0.3029466974;
    nodeX[2542] = -0.4507422593;
    nodeX[2543] = -0.8396753624;

    weight[848] = 0.01112159279;
    nodeX[2544] = -0.3029466974;
    nodeX[2545] = 0.4507422593;
    nodeX[2546] = -0.8396753624;

    weight[849] = 0.01112159279;
    nodeX[2547] = -0.3029466974;
    nodeX[2548] = -0.4507422593;
    nodeX[2549] = -0.8396753624;

    weight[850] = 0.01112159279;
    nodeX[2550] = 0.8396753624;
    nodeX[2551] = 0.3029466974;
    nodeX[2552] = 0.4507422593;

    weight[851] = 0.01112159279;
    nodeX[2553] = -0.8396753624;
    nodeX[2554] = 0.3029466974;
    nodeX[2555] = 0.4507422593;

    weight[852] = 0.01112159279;
    nodeX[2556] = 0.8396753624;
    nodeX[2557] = -0.3029466974;
    nodeX[2558] = 0.4507422593;

    weight[853] = 0.01112159279;
    nodeX[2559] = 0.8396753624;
    nodeX[2560] = 0.3029466974;
    nodeX[2561] = -0.4507422593;

    weight[854] = 0.01112159279;
    nodeX[2562] = -0.8396753624;
    nodeX[2563] = -0.3029466974;
    nodeX[2564] = 0.4507422593;

    weight[855] = 0.01112159279;
    nodeX[2565] = 0.8396753624;
    nodeX[2566] = -0.3029466974;
    nodeX[2567] = -0.4507422593;

    weight[856] = 0.01112159279;
    nodeX[2568] = -0.8396753624;
    nodeX[2569] = 0.3029466974;
    nodeX[2570] = -0.4507422593;

    weight[857] = 0.01112159279;
    nodeX[2571] = -0.8396753624;
    nodeX[2572] = -0.3029466974;
    nodeX[2573] = -0.4507422593;

    weight[858] = 0.01112159279;
    nodeX[2574] = 0.4507422593;
    nodeX[2575] = 0.3029466974;
    nodeX[2576] = 0.8396753624;

    weight[859] = 0.01112159279;
    nodeX[2577] = -0.4507422593;
    nodeX[2578] = 0.3029466974;
    nodeX[2579] = 0.8396753624;

    weight[860] = 0.01112159279;
    nodeX[2580] = 0.4507422593;
    nodeX[2581] = -0.3029466974;
    nodeX[2582] = 0.8396753624;

    weight[861] = 0.01112159279;
    nodeX[2583] = 0.4507422593;
    nodeX[2584] = 0.3029466974;
    nodeX[2585] = -0.8396753624;

    weight[862] = 0.01112159279;
    nodeX[2586] = -0.4507422593;
    nodeX[2587] = -0.3029466974;
    nodeX[2588] = 0.8396753624;

    weight[863] = 0.01112159279;
    nodeX[2589] = 0.4507422593;
    nodeX[2590] = -0.3029466974;
    nodeX[2591] = -0.8396753624;

    weight[864] = 0.01112159279;
    nodeX[2592] = -0.4507422593;
    nodeX[2593] = 0.3029466974;
    nodeX[2594] = -0.8396753624;

    weight[865] = 0.01112159279;
    nodeX[2595] = -0.4507422593;
    nodeX[2596] = -0.3029466974;
    nodeX[2597] = -0.8396753624;

    weight[866] = 0.01133655308;
    nodeX[2598] = 0.8165288564;
    nodeX[2599] = 0.5632123021;
    nodeX[2600] = 0.1267774801;

    weight[867] = 0.01133655308;
    nodeX[2601] = -0.8165288564;
    nodeX[2602] = 0.5632123021;
    nodeX[2603] = 0.1267774801;

    weight[868] = 0.01133655308;
    nodeX[2604] = 0.8165288564;
    nodeX[2605] = -0.5632123021;
    nodeX[2606] = 0.1267774801;

    weight[869] = 0.01133655308;
    nodeX[2607] = 0.8165288564;
    nodeX[2608] = 0.5632123021;
    nodeX[2609] = -0.1267774801;

    weight[870] = 0.01133655308;
    nodeX[2610] = -0.8165288564;
    nodeX[2611] = -0.5632123021;
    nodeX[2612] = 0.1267774801;

    weight[871] = 0.01133655308;
    nodeX[2613] = 0.8165288564;
    nodeX[2614] = -0.5632123021;
    nodeX[2615] = -0.1267774801;

    weight[872] = 0.01133655308;
    nodeX[2616] = -0.8165288564;
    nodeX[2617] = 0.5632123021;
    nodeX[2618] = -0.1267774801;

    weight[873] = 0.01133655308;
    nodeX[2619] = -0.8165288564;
    nodeX[2620] = -0.5632123021;
    nodeX[2621] = -0.1267774801;

    weight[874] = 0.01133655308;
    nodeX[2622] = 0.5632123021;
    nodeX[2623] = 0.8165288564;
    nodeX[2624] = 0.1267774801;

    weight[875] = 0.01133655308;
    nodeX[2625] = -0.5632123021;
    nodeX[2626] = 0.8165288564;
    nodeX[2627] = 0.1267774801;

    weight[876] = 0.01133655308;
    nodeX[2628] = 0.5632123021;
    nodeX[2629] = -0.8165288564;
    nodeX[2630] = 0.1267774801;

    weight[877] = 0.01133655308;
    nodeX[2631] = 0.5632123021;
    nodeX[2632] = 0.8165288564;
    nodeX[2633] = -0.1267774801;

    weight[878] = 0.01133655308;
    nodeX[2634] = -0.5632123021;
    nodeX[2635] = -0.8165288564;
    nodeX[2636] = 0.1267774801;

    weight[879] = 0.01133655308;
    nodeX[2637] = 0.5632123021;
    nodeX[2638] = -0.8165288564;
    nodeX[2639] = -0.1267774801;

    weight[880] = 0.01133655308;
    nodeX[2640] = -0.5632123021;
    nodeX[2641] = 0.8165288564;
    nodeX[2642] = -0.1267774801;

    weight[881] = 0.01133655308;
    nodeX[2643] = -0.5632123021;
    nodeX[2644] = -0.8165288564;
    nodeX[2645] = -0.1267774801;

    weight[882] = 0.01133655308;
    nodeX[2646] = 0.1267774801;
    nodeX[2647] = 0.8165288564;
    nodeX[2648] = 0.5632123021;

    weight[883] = 0.01133655308;
    nodeX[2649] = -0.1267774801;
    nodeX[2650] = 0.8165288564;
    nodeX[2651] = 0.5632123021;

    weight[884] = 0.01133655308;
    nodeX[2652] = 0.1267774801;
    nodeX[2653] = -0.8165288564;
    nodeX[2654] = 0.5632123021;

    weight[885] = 0.01133655308;
    nodeX[2655] = 0.1267774801;
    nodeX[2656] = 0.8165288564;
    nodeX[2657] = -0.5632123021;

    weight[886] = 0.01133655308;
    nodeX[2658] = -0.1267774801;
    nodeX[2659] = -0.8165288564;
    nodeX[2660] = 0.5632123021;

    weight[887] = 0.01133655308;
    nodeX[2661] = 0.1267774801;
    nodeX[2662] = -0.8165288564;
    nodeX[2663] = -0.5632123021;

    weight[888] = 0.01133655308;
    nodeX[2664] = -0.1267774801;
    nodeX[2665] = 0.8165288564;
    nodeX[2666] = -0.5632123021;

    weight[889] = 0.01133655308;
    nodeX[2667] = -0.1267774801;
    nodeX[2668] = -0.8165288564;
    nodeX[2669] = -0.5632123021;

    weight[890] = 0.01133655308;
    nodeX[2670] = 0.1267774801;
    nodeX[2671] = 0.5632123021;
    nodeX[2672] = 0.8165288564;

    weight[891] = 0.01133655308;
    nodeX[2673] = -0.1267774801;
    nodeX[2674] = 0.5632123021;
    nodeX[2675] = 0.8165288564;

    weight[892] = 0.01133655308;
    nodeX[2676] = 0.1267774801;
    nodeX[2677] = -0.5632123021;
    nodeX[2678] = 0.8165288564;

    weight[893] = 0.01133655308;
    nodeX[2679] = 0.1267774801;
    nodeX[2680] = 0.5632123021;
    nodeX[2681] = -0.8165288564;

    weight[894] = 0.01133655308;
    nodeX[2682] = -0.1267774801;
    nodeX[2683] = -0.5632123021;
    nodeX[2684] = 0.8165288564;

    weight[895] = 0.01133655308;
    nodeX[2685] = 0.1267774801;
    nodeX[2686] = -0.5632123021;
    nodeX[2687] = -0.8165288564;

    weight[896] = 0.01133655308;
    nodeX[2688] = -0.1267774801;
    nodeX[2689] = 0.5632123021;
    nodeX[2690] = -0.8165288564;

    weight[897] = 0.01133655308;
    nodeX[2691] = -0.1267774801;
    nodeX[2692] = -0.5632123021;
    nodeX[2693] = -0.8165288564;

    weight[898] = 0.01133655308;
    nodeX[2694] = 0.8165288564;
    nodeX[2695] = 0.1267774801;
    nodeX[2696] = 0.5632123021;

    weight[899] = 0.01133655308;
    nodeX[2697] = -0.8165288564;
    nodeX[2698] = 0.1267774801;
    nodeX[2699] = 0.5632123021;

    weight[900] = 0.01133655308;
    nodeX[2700] = 0.8165288564;
    nodeX[2701] = -0.1267774801;
    nodeX[2702] = 0.5632123021;

    weight[901] = 0.01133655308;
    nodeX[2703] = 0.8165288564;
    nodeX[2704] = 0.1267774801;
    nodeX[2705] = -0.5632123021;

    weight[902] = 0.01133655308;
    nodeX[2706] = -0.8165288564;
    nodeX[2707] = -0.1267774801;
    nodeX[2708] = 0.5632123021;

    weight[903] = 0.01133655308;
    nodeX[2709] = 0.8165288564;
    nodeX[2710] = -0.1267774801;
    nodeX[2711] = -0.5632123021;

    weight[904] = 0.01133655308;
    nodeX[2712] = -0.8165288564;
    nodeX[2713] = 0.1267774801;
    nodeX[2714] = -0.5632123021;

    weight[905] = 0.01133655308;
    nodeX[2715] = -0.8165288564;
    nodeX[2716] = -0.1267774801;
    nodeX[2717] = -0.5632123021;

    weight[906] = 0.01133655308;
    nodeX[2718] = 0.5632123021;
    nodeX[2719] = 0.1267774801;
    nodeX[2720] = 0.8165288564;

    weight[907] = 0.01133655308;
    nodeX[2721] = -0.5632123021;
    nodeX[2722] = 0.1267774801;
    nodeX[2723] = 0.8165288564;

    weight[908] = 0.01133655308;
    nodeX[2724] = 0.5632123021;
    nodeX[2725] = -0.1267774801;
    nodeX[2726] = 0.8165288564;

    weight[909] = 0.01133655308;
    nodeX[2727] = 0.5632123021;
    nodeX[2728] = 0.1267774801;
    nodeX[2729] = -0.8165288564;

    weight[910] = 0.01133655308;
    nodeX[2730] = -0.5632123021;
    nodeX[2731] = -0.1267774801;
    nodeX[2732] = 0.8165288564;

    weight[911] = 0.01133655308;
    nodeX[2733] = 0.5632123021;
    nodeX[2734] = -0.1267774801;
    nodeX[2735] = -0.8165288564;

    weight[912] = 0.01133655308;
    nodeX[2736] = -0.5632123021;
    nodeX[2737] = 0.1267774801;
    nodeX[2738] = -0.8165288564;

    weight[913] = 0.01133655308;
    nodeX[2739] = -0.5632123021;
    nodeX[2740] = -0.1267774801;
    nodeX[2741] = -0.8165288564;

    weight[914] = 0.01132241513;
    nodeX[2742] = 0.8015469371;
    nodeX[2743] = 0.543430357;
    nodeX[2744] = 0.2494112162;

    weight[915] = 0.01132241513;
    nodeX[2745] = -0.8015469371;
    nodeX[2746] = 0.543430357;
    nodeX[2747] = 0.2494112162;

    weight[916] = 0.01132241513;
    nodeX[2748] = 0.8015469371;
    nodeX[2749] = -0.543430357;
    nodeX[2750] = 0.2494112162;

    weight[917] = 0.01132241513;
    nodeX[2751] = 0.8015469371;
    nodeX[2752] = 0.543430357;
    nodeX[2753] = -0.2494112162;

    weight[918] = 0.01132241513;
    nodeX[2754] = -0.8015469371;
    nodeX[2755] = -0.543430357;
    nodeX[2756] = 0.2494112162;

    weight[919] = 0.01132241513;
    nodeX[2757] = 0.8015469371;
    nodeX[2758] = -0.543430357;
    nodeX[2759] = -0.2494112162;

    weight[920] = 0.01132241513;
    nodeX[2760] = -0.8015469371;
    nodeX[2761] = 0.543430357;
    nodeX[2762] = -0.2494112162;

    weight[921] = 0.01132241513;
    nodeX[2763] = -0.8015469371;
    nodeX[2764] = -0.543430357;
    nodeX[2765] = -0.2494112162;

    weight[922] = 0.01132241513;
    nodeX[2766] = 0.543430357;
    nodeX[2767] = 0.8015469371;
    nodeX[2768] = 0.2494112162;

    weight[923] = 0.01132241513;
    nodeX[2769] = -0.543430357;
    nodeX[2770] = 0.8015469371;
    nodeX[2771] = 0.2494112162;

    weight[924] = 0.01132241513;
    nodeX[2772] = 0.543430357;
    nodeX[2773] = -0.8015469371;
    nodeX[2774] = 0.2494112162;

    weight[925] = 0.01132241513;
    nodeX[2775] = 0.543430357;
    nodeX[2776] = 0.8015469371;
    nodeX[2777] = -0.2494112162;

    weight[926] = 0.01132241513;
    nodeX[2778] = -0.543430357;
    nodeX[2779] = -0.8015469371;
    nodeX[2780] = 0.2494112162;

    weight[927] = 0.01132241513;
    nodeX[2781] = 0.543430357;
    nodeX[2782] = -0.8015469371;
    nodeX[2783] = -0.2494112162;

    weight[928] = 0.01132241513;
    nodeX[2784] = -0.543430357;
    nodeX[2785] = 0.8015469371;
    nodeX[2786] = -0.2494112162;

    weight[929] = 0.01132241513;
    nodeX[2787] = -0.543430357;
    nodeX[2788] = -0.8015469371;
    nodeX[2789] = -0.2494112162;

    weight[930] = 0.01132241513;
    nodeX[2790] = 0.2494112162;
    nodeX[2791] = 0.8015469371;
    nodeX[2792] = 0.543430357;

    weight[931] = 0.01132241513;
    nodeX[2793] = -0.2494112162;
    nodeX[2794] = 0.8015469371;
    nodeX[2795] = 0.543430357;

    weight[932] = 0.01132241513;
    nodeX[2796] = 0.2494112162;
    nodeX[2797] = -0.8015469371;
    nodeX[2798] = 0.543430357;

    weight[933] = 0.01132241513;
    nodeX[2799] = 0.2494112162;
    nodeX[2800] = 0.8015469371;
    nodeX[2801] = -0.543430357;

    weight[934] = 0.01132241513;
    nodeX[2802] = -0.2494112162;
    nodeX[2803] = -0.8015469371;
    nodeX[2804] = 0.543430357;

    weight[935] = 0.01132241513;
    nodeX[2805] = 0.2494112162;
    nodeX[2806] = -0.8015469371;
    nodeX[2807] = -0.543430357;

    weight[936] = 0.01132241513;
    nodeX[2808] = -0.2494112162;
    nodeX[2809] = 0.8015469371;
    nodeX[2810] = -0.543430357;

    weight[937] = 0.01132241513;
    nodeX[2811] = -0.2494112162;
    nodeX[2812] = -0.8015469371;
    nodeX[2813] = -0.543430357;

    weight[938] = 0.01132241513;
    nodeX[2814] = 0.2494112162;
    nodeX[2815] = 0.543430357;
    nodeX[2816] = 0.8015469371;

    weight[939] = 0.01132241513;
    nodeX[2817] = -0.2494112162;
    nodeX[2818] = 0.543430357;
    nodeX[2819] = 0.8015469371;

    weight[940] = 0.01132241513;
    nodeX[2820] = 0.2494112162;
    nodeX[2821] = -0.543430357;
    nodeX[2822] = 0.8015469371;

    weight[941] = 0.01132241513;
    nodeX[2823] = 0.2494112162;
    nodeX[2824] = 0.543430357;
    nodeX[2825] = -0.8015469371;

    weight[942] = 0.01132241513;
    nodeX[2826] = -0.2494112162;
    nodeX[2827] = -0.543430357;
    nodeX[2828] = 0.8015469371;

    weight[943] = 0.01132241513;
    nodeX[2829] = 0.2494112162;
    nodeX[2830] = -0.543430357;
    nodeX[2831] = -0.8015469371;

    weight[944] = 0.01132241513;
    nodeX[2832] = -0.2494112162;
    nodeX[2833] = 0.543430357;
    nodeX[2834] = -0.8015469371;

    weight[945] = 0.01132241513;
    nodeX[2835] = -0.2494112162;
    nodeX[2836] = -0.543430357;
    nodeX[2837] = -0.8015469371;

    weight[946] = 0.01132241513;
    nodeX[2838] = 0.8015469371;
    nodeX[2839] = 0.2494112162;
    nodeX[2840] = 0.543430357;

    weight[947] = 0.01132241513;
    nodeX[2841] = -0.8015469371;
    nodeX[2842] = 0.2494112162;
    nodeX[2843] = 0.543430357;

    weight[948] = 0.01132241513;
    nodeX[2844] = 0.8015469371;
    nodeX[2845] = -0.2494112162;
    nodeX[2846] = 0.543430357;

    weight[949] = 0.01132241513;
    nodeX[2847] = 0.8015469371;
    nodeX[2848] = 0.2494112162;
    nodeX[2849] = -0.543430357;

    weight[950] = 0.01132241513;
    nodeX[2850] = -0.8015469371;
    nodeX[2851] = -0.2494112162;
    nodeX[2852] = 0.543430357;

    weight[951] = 0.01132241513;
    nodeX[2853] = 0.8015469371;
    nodeX[2854] = -0.2494112162;
    nodeX[2855] = -0.543430357;

    weight[952] = 0.01132241513;
    nodeX[2856] = -0.8015469371;
    nodeX[2857] = 0.2494112162;
    nodeX[2858] = -0.543430357;

    weight[953] = 0.01132241513;
    nodeX[2859] = -0.8015469371;
    nodeX[2860] = -0.2494112162;
    nodeX[2861] = -0.543430357;

    weight[954] = 0.01132241513;
    nodeX[2862] = 0.543430357;
    nodeX[2863] = 0.2494112162;
    nodeX[2864] = 0.8015469371;

    weight[955] = 0.01132241513;
    nodeX[2865] = -0.543430357;
    nodeX[2866] = 0.2494112162;
    nodeX[2867] = 0.8015469371;

    weight[956] = 0.01132241513;
    nodeX[2868] = 0.543430357;
    nodeX[2869] = -0.2494112162;
    nodeX[2870] = 0.8015469371;

    weight[957] = 0.01132241513;
    nodeX[2871] = 0.543430357;
    nodeX[2872] = 0.2494112162;
    nodeX[2873] = -0.8015469371;

    weight[958] = 0.01132241513;
    nodeX[2874] = -0.543430357;
    nodeX[2875] = -0.2494112162;
    nodeX[2876] = 0.8015469371;

    weight[959] = 0.01132241513;
    nodeX[2877] = 0.543430357;
    nodeX[2878] = -0.2494112162;
    nodeX[2879] = -0.8015469371;

    weight[960] = 0.01132241513;
    nodeX[2880] = -0.543430357;
    nodeX[2881] = 0.2494112162;
    nodeX[2882] = -0.8015469371;

    weight[961] = 0.01132241513;
    nodeX[2883] = -0.543430357;
    nodeX[2884] = -0.2494112162;
    nodeX[2885] = -0.8015469371;

    weight[962] = 0.01133825034;
    nodeX[2886] = 0.7773563069;
    nodeX[2887] = 0.5123518486;
    nodeX[2888] = 0.3649832261;

    weight[963] = 0.01133825034;
    nodeX[2889] = -0.7773563069;
    nodeX[2890] = 0.5123518486;
    nodeX[2891] = 0.3649832261;

    weight[964] = 0.01133825034;
    nodeX[2892] = 0.7773563069;
    nodeX[2893] = -0.5123518486;
    nodeX[2894] = 0.3649832261;

    weight[965] = 0.01133825034;
    nodeX[2895] = 0.7773563069;
    nodeX[2896] = 0.5123518486;
    nodeX[2897] = -0.3649832261;

    weight[966] = 0.01133825034;
    nodeX[2898] = -0.7773563069;
    nodeX[2899] = -0.5123518486;
    nodeX[2900] = 0.3649832261;

    weight[967] = 0.01133825034;
    nodeX[2901] = 0.7773563069;
    nodeX[2902] = -0.5123518486;
    nodeX[2903] = -0.3649832261;

    weight[968] = 0.01133825034;
    nodeX[2904] = -0.7773563069;
    nodeX[2905] = 0.5123518486;
    nodeX[2906] = -0.3649832261;

    weight[969] = 0.01133825034;
    nodeX[2907] = -0.7773563069;
    nodeX[2908] = -0.5123518486;
    nodeX[2909] = -0.3649832261;

    weight[970] = 0.01133825034;
    nodeX[2910] = 0.5123518486;
    nodeX[2911] = 0.7773563069;
    nodeX[2912] = 0.3649832261;

    weight[971] = 0.01133825034;
    nodeX[2913] = -0.5123518486;
    nodeX[2914] = 0.7773563069;
    nodeX[2915] = 0.3649832261;

    weight[972] = 0.01133825034;
    nodeX[2916] = 0.5123518486;
    nodeX[2917] = -0.7773563069;
    nodeX[2918] = 0.3649832261;

    weight[973] = 0.01133825034;
    nodeX[2919] = 0.5123518486;
    nodeX[2920] = 0.7773563069;
    nodeX[2921] = -0.3649832261;

    weight[974] = 0.01133825034;
    nodeX[2922] = -0.5123518486;
    nodeX[2923] = -0.7773563069;
    nodeX[2924] = 0.3649832261;

    weight[975] = 0.01133825034;
    nodeX[2925] = 0.5123518486;
    nodeX[2926] = -0.7773563069;
    nodeX[2927] = -0.3649832261;

    weight[976] = 0.01133825034;
    nodeX[2928] = -0.5123518486;
    nodeX[2929] = 0.7773563069;
    nodeX[2930] = -0.3649832261;

    weight[977] = 0.01133825034;
    nodeX[2931] = -0.5123518486;
    nodeX[2932] = -0.7773563069;
    nodeX[2933] = -0.3649832261;

    weight[978] = 0.01133825034;
    nodeX[2934] = 0.3649832261;
    nodeX[2935] = 0.7773563069;
    nodeX[2936] = 0.5123518486;

    weight[979] = 0.01133825034;
    nodeX[2937] = -0.3649832261;
    nodeX[2938] = 0.7773563069;
    nodeX[2939] = 0.5123518486;

    weight[980] = 0.01133825034;
    nodeX[2940] = 0.3649832261;
    nodeX[2941] = -0.7773563069;
    nodeX[2942] = 0.5123518486;

    weight[981] = 0.01133825034;
    nodeX[2943] = 0.3649832261;
    nodeX[2944] = 0.7773563069;
    nodeX[2945] = -0.5123518486;

    weight[982] = 0.01133825034;
    nodeX[2946] = -0.3649832261;
    nodeX[2947] = -0.7773563069;
    nodeX[2948] = 0.5123518486;

    weight[983] = 0.01133825034;
    nodeX[2949] = 0.3649832261;
    nodeX[2950] = -0.7773563069;
    nodeX[2951] = -0.5123518486;

    weight[984] = 0.01133825034;
    nodeX[2952] = -0.3649832261;
    nodeX[2953] = 0.7773563069;
    nodeX[2954] = -0.5123518486;

    weight[985] = 0.01133825034;
    nodeX[2955] = -0.3649832261;
    nodeX[2956] = -0.7773563069;
    nodeX[2957] = -0.5123518486;

    weight[986] = 0.01133825034;
    nodeX[2958] = 0.3649832261;
    nodeX[2959] = 0.5123518486;
    nodeX[2960] = 0.7773563069;

    weight[987] = 0.01133825034;
    nodeX[2961] = -0.3649832261;
    nodeX[2962] = 0.5123518486;
    nodeX[2963] = 0.7773563069;

    weight[988] = 0.01133825034;
    nodeX[2964] = 0.3649832261;
    nodeX[2965] = -0.5123518486;
    nodeX[2966] = 0.7773563069;

    weight[989] = 0.01133825034;
    nodeX[2967] = 0.3649832261;
    nodeX[2968] = 0.5123518486;
    nodeX[2969] = -0.7773563069;

    weight[990] = 0.01133825034;
    nodeX[2970] = -0.3649832261;
    nodeX[2971] = -0.5123518486;
    nodeX[2972] = 0.7773563069;

    weight[991] = 0.01133825034;
    nodeX[2973] = 0.3649832261;
    nodeX[2974] = -0.5123518486;
    nodeX[2975] = -0.7773563069;

    weight[992] = 0.01133825034;
    nodeX[2976] = -0.3649832261;
    nodeX[2977] = 0.5123518486;
    nodeX[2978] = -0.7773563069;

    weight[993] = 0.01133825034;
    nodeX[2979] = -0.3649832261;
    nodeX[2980] = -0.5123518486;
    nodeX[2981] = -0.7773563069;

    weight[994] = 0.01133825034;
    nodeX[2982] = 0.7773563069;
    nodeX[2983] = 0.3649832261;
    nodeX[2984] = 0.5123518486;

    weight[995] = 0.01133825034;
    nodeX[2985] = -0.7773563069;
    nodeX[2986] = 0.3649832261;
    nodeX[2987] = 0.5123518486;

    weight[996] = 0.01133825034;
    nodeX[2988] = 0.7773563069;
    nodeX[2989] = -0.3649832261;
    nodeX[2990] = 0.5123518486;

    weight[997] = 0.01133825034;
    nodeX[2991] = 0.7773563069;
    nodeX[2992] = 0.3649832261;
    nodeX[2993] = -0.5123518486;

    weight[998] = 0.01133825034;
    nodeX[2994] = -0.7773563069;
    nodeX[2995] = -0.3649832261;
    nodeX[2996] = 0.5123518486;

    weight[999] = 0.01133825034;
    nodeX[2997] = 0.7773563069;
    nodeX[2998] = -0.3649832261;
    nodeX[2999] = -0.5123518486;

    weight[1000] = 0.01133825034;
    nodeX[3000] = -0.7773563069;
    nodeX[3001] = 0.3649832261;
    nodeX[3002] = -0.5123518486;

    weight[1001] = 0.01133825034;
    nodeX[3003] = -0.7773563069;
    nodeX[3004] = -0.3649832261;
    nodeX[3005] = -0.5123518486;

    weight[1002] = 0.01133825034;
    nodeX[3006] = 0.5123518486;
    nodeX[3007] = 0.3649832261;
    nodeX[3008] = 0.7773563069;

    weight[1003] = 0.01133825034;
    nodeX[3009] = -0.5123518486;
    nodeX[3010] = 0.3649832261;
    nodeX[3011] = 0.7773563069;

    weight[1004] = 0.01133825034;
    nodeX[3012] = 0.5123518486;
    nodeX[3013] = -0.3649832261;
    nodeX[3014] = 0.7773563069;

    weight[1005] = 0.01133825034;
    nodeX[3015] = 0.5123518486;
    nodeX[3016] = 0.3649832261;
    nodeX[3017] = -0.7773563069;

    weight[1006] = 0.01133825034;
    nodeX[3018] = -0.5123518486;
    nodeX[3019] = -0.3649832261;
    nodeX[3020] = 0.7773563069;

    weight[1007] = 0.01133825034;
    nodeX[3021] = 0.5123518486;
    nodeX[3022] = -0.3649832261;
    nodeX[3023] = -0.7773563069;

    weight[1008] = 0.01133825034;
    nodeX[3024] = -0.5123518486;
    nodeX[3025] = 0.3649832261;
    nodeX[3026] = -0.7773563069;

    weight[1009] = 0.01133825034;
    nodeX[3027] = -0.5123518486;
    nodeX[3028] = -0.3649832261;
    nodeX[3029] = -0.7773563069;

    weight[1010] = 0.01150830253;
    nodeX[3030] = 0.7661621214;
    nodeX[3031] = 0.6394279635;
    nodeX[3032] = 0.06424549224;

    weight[1011] = 0.01150830253;
    nodeX[3033] = -0.7661621214;
    nodeX[3034] = 0.6394279635;
    nodeX[3035] = 0.06424549224;

    weight[1012] = 0.01150830253;
    nodeX[3036] = 0.7661621214;
    nodeX[3037] = -0.6394279635;
    nodeX[3038] = 0.06424549224;

    weight[1013] = 0.01150830253;
    nodeX[3039] = 0.7661621214;
    nodeX[3040] = 0.6394279635;
    nodeX[3041] = -0.06424549224;

    weight[1014] = 0.01150830253;
    nodeX[3042] = -0.7661621214;
    nodeX[3043] = -0.6394279635;
    nodeX[3044] = 0.06424549224;

    weight[1015] = 0.01150830253;
    nodeX[3045] = 0.7661621214;
    nodeX[3046] = -0.6394279635;
    nodeX[3047] = -0.06424549224;

    weight[1016] = 0.01150830253;
    nodeX[3048] = -0.7661621214;
    nodeX[3049] = 0.6394279635;
    nodeX[3050] = -0.06424549224;

    weight[1017] = 0.01150830253;
    nodeX[3051] = -0.7661621214;
    nodeX[3052] = -0.6394279635;
    nodeX[3053] = -0.06424549224;

    weight[1018] = 0.01150830253;
    nodeX[3054] = 0.6394279635;
    nodeX[3055] = 0.7661621214;
    nodeX[3056] = 0.06424549224;

    weight[1019] = 0.01150830253;
    nodeX[3057] = -0.6394279635;
    nodeX[3058] = 0.7661621214;
    nodeX[3059] = 0.06424549224;

    weight[1020] = 0.01150830253;
    nodeX[3060] = 0.6394279635;
    nodeX[3061] = -0.7661621214;
    nodeX[3062] = 0.06424549224;

    weight[1021] = 0.01150830253;
    nodeX[3063] = 0.6394279635;
    nodeX[3064] = 0.7661621214;
    nodeX[3065] = -0.06424549224;

    weight[1022] = 0.01150830253;
    nodeX[3066] = -0.6394279635;
    nodeX[3067] = -0.7661621214;
    nodeX[3068] = 0.06424549224;

    weight[1023] = 0.01150830253;
    nodeX[3069] = 0.6394279635;
    nodeX[3070] = -0.7661621214;
    nodeX[3071] = -0.06424549224;

    weight[1024] = 0.01150830253;
    nodeX[3072] = -0.6394279635;
    nodeX[3073] = 0.7661621214;
    nodeX[3074] = -0.06424549224;

    weight[1025] = 0.01150830253;
    nodeX[3075] = -0.6394279635;
    nodeX[3076] = -0.7661621214;
    nodeX[3077] = -0.06424549224;

    weight[1026] = 0.01150830253;
    nodeX[3078] = 0.06424549224;
    nodeX[3079] = 0.7661621214;
    nodeX[3080] = 0.6394279635;

    weight[1027] = 0.01150830253;
    nodeX[3081] = -0.06424549224;
    nodeX[3082] = 0.7661621214;
    nodeX[3083] = 0.6394279635;

    weight[1028] = 0.01150830253;
    nodeX[3084] = 0.06424549224;
    nodeX[3085] = -0.7661621214;
    nodeX[3086] = 0.6394279635;

    weight[1029] = 0.01150830253;
    nodeX[3087] = 0.06424549224;
    nodeX[3088] = 0.7661621214;
    nodeX[3089] = -0.6394279635;

    weight[1030] = 0.01150830253;
    nodeX[3090] = -0.06424549224;
    nodeX[3091] = -0.7661621214;
    nodeX[3092] = 0.6394279635;

    weight[1031] = 0.01150830253;
    nodeX[3093] = 0.06424549224;
    nodeX[3094] = -0.7661621214;
    nodeX[3095] = -0.6394279635;

    weight[1032] = 0.01150830253;
    nodeX[3096] = -0.06424549224;
    nodeX[3097] = 0.7661621214;
    nodeX[3098] = -0.6394279635;

    weight[1033] = 0.01150830253;
    nodeX[3099] = -0.06424549224;
    nodeX[3100] = -0.7661621214;
    nodeX[3101] = -0.6394279635;

    weight[1034] = 0.01150830253;
    nodeX[3102] = 0.06424549224;
    nodeX[3103] = 0.6394279635;
    nodeX[3104] = 0.7661621214;

    weight[1035] = 0.01150830253;
    nodeX[3105] = -0.06424549224;
    nodeX[3106] = 0.6394279635;
    nodeX[3107] = 0.7661621214;

    weight[1036] = 0.01150830253;
    nodeX[3108] = 0.06424549224;
    nodeX[3109] = -0.6394279635;
    nodeX[3110] = 0.7661621214;

    weight[1037] = 0.01150830253;
    nodeX[3111] = 0.06424549224;
    nodeX[3112] = 0.6394279635;
    nodeX[3113] = -0.7661621214;

    weight[1038] = 0.01150830253;
    nodeX[3114] = -0.06424549224;
    nodeX[3115] = -0.6394279635;
    nodeX[3116] = 0.7661621214;

    weight[1039] = 0.01150830253;
    nodeX[3117] = 0.06424549224;
    nodeX[3118] = -0.6394279635;
    nodeX[3119] = -0.7661621214;

    weight[1040] = 0.01150830253;
    nodeX[3120] = -0.06424549224;
    nodeX[3121] = 0.6394279635;
    nodeX[3122] = -0.7661621214;

    weight[1041] = 0.01150830253;
    nodeX[3123] = -0.06424549224;
    nodeX[3124] = -0.6394279635;
    nodeX[3125] = -0.7661621214;

    weight[1042] = 0.01150830253;
    nodeX[3126] = 0.7661621214;
    nodeX[3127] = 0.06424549224;
    nodeX[3128] = 0.6394279635;

    weight[1043] = 0.01150830253;
    nodeX[3129] = -0.7661621214;
    nodeX[3130] = 0.06424549224;
    nodeX[3131] = 0.6394279635;

    weight[1044] = 0.01150830253;
    nodeX[3132] = 0.7661621214;
    nodeX[3133] = -0.06424549224;
    nodeX[3134] = 0.6394279635;

    weight[1045] = 0.01150830253;
    nodeX[3135] = 0.7661621214;
    nodeX[3136] = 0.06424549224;
    nodeX[3137] = -0.6394279635;

    weight[1046] = 0.01150830253;
    nodeX[3138] = -0.7661621214;
    nodeX[3139] = -0.06424549224;
    nodeX[3140] = 0.6394279635;

    weight[1047] = 0.01150830253;
    nodeX[3141] = 0.7661621214;
    nodeX[3142] = -0.06424549224;
    nodeX[3143] = -0.6394279635;

    weight[1048] = 0.01150830253;
    nodeX[3144] = -0.7661621214;
    nodeX[3145] = 0.06424549224;
    nodeX[3146] = -0.6394279635;

    weight[1049] = 0.01150830253;
    nodeX[3147] = -0.7661621214;
    nodeX[3148] = -0.06424549224;
    nodeX[3149] = -0.6394279635;

    weight[1050] = 0.01150830253;
    nodeX[3150] = 0.6394279635;
    nodeX[3151] = 0.06424549224;
    nodeX[3152] = 0.7661621214;

    weight[1051] = 0.01150830253;
    nodeX[3153] = -0.6394279635;
    nodeX[3154] = 0.06424549224;
    nodeX[3155] = 0.7661621214;

    weight[1052] = 0.01150830253;
    nodeX[3156] = 0.6394279635;
    nodeX[3157] = -0.06424549224;
    nodeX[3158] = 0.7661621214;

    weight[1053] = 0.01150830253;
    nodeX[3159] = 0.6394279635;
    nodeX[3160] = 0.06424549224;
    nodeX[3161] = -0.7661621214;

    weight[1054] = 0.01150830253;
    nodeX[3162] = -0.6394279635;
    nodeX[3163] = -0.06424549224;
    nodeX[3164] = 0.7661621214;

    weight[1055] = 0.01150830253;
    nodeX[3165] = 0.6394279635;
    nodeX[3166] = -0.06424549224;
    nodeX[3167] = -0.7661621214;

    weight[1056] = 0.01150830253;
    nodeX[3168] = -0.6394279635;
    nodeX[3169] = 0.06424549224;
    nodeX[3170] = -0.7661621214;

    weight[1057] = 0.01150830253;
    nodeX[3171] = -0.6394279635;
    nodeX[3172] = -0.06424549224;
    nodeX[3173] = -0.7661621214;

    weight[1058] = 0.01147507935;
    nodeX[3174] = 0.7553584144;
    nodeX[3175] = 0.6269805509;
    nodeX[3176] = 0.1906018223;

    weight[1059] = 0.01147507935;
    nodeX[3177] = -0.7553584144;
    nodeX[3178] = 0.6269805509;
    nodeX[3179] = 0.1906018223;

    weight[1060] = 0.01147507935;
    nodeX[3180] = 0.7553584144;
    nodeX[3181] = -0.6269805509;
    nodeX[3182] = 0.1906018223;

    weight[1061] = 0.01147507935;
    nodeX[3183] = 0.7553584144;
    nodeX[3184] = 0.6269805509;
    nodeX[3185] = -0.1906018223;

    weight[1062] = 0.01147507935;
    nodeX[3186] = -0.7553584144;
    nodeX[3187] = -0.6269805509;
    nodeX[3188] = 0.1906018223;

    weight[1063] = 0.01147507935;
    nodeX[3189] = 0.7553584144;
    nodeX[3190] = -0.6269805509;
    nodeX[3191] = -0.1906018223;

    weight[1064] = 0.01147507935;
    nodeX[3192] = -0.7553584144;
    nodeX[3193] = 0.6269805509;
    nodeX[3194] = -0.1906018223;

    weight[1065] = 0.01147507935;
    nodeX[3195] = -0.7553584144;
    nodeX[3196] = -0.6269805509;
    nodeX[3197] = -0.1906018223;

    weight[1066] = 0.01147507935;
    nodeX[3198] = 0.6269805509;
    nodeX[3199] = 0.7553584144;
    nodeX[3200] = 0.1906018223;

    weight[1067] = 0.01147507935;
    nodeX[3201] = -0.6269805509;
    nodeX[3202] = 0.7553584144;
    nodeX[3203] = 0.1906018223;

    weight[1068] = 0.01147507935;
    nodeX[3204] = 0.6269805509;
    nodeX[3205] = -0.7553584144;
    nodeX[3206] = 0.1906018223;

    weight[1069] = 0.01147507935;
    nodeX[3207] = 0.6269805509;
    nodeX[3208] = 0.7553584144;
    nodeX[3209] = -0.1906018223;

    weight[1070] = 0.01147507935;
    nodeX[3210] = -0.6269805509;
    nodeX[3211] = -0.7553584144;
    nodeX[3212] = 0.1906018223;

    weight[1071] = 0.01147507935;
    nodeX[3213] = 0.6269805509;
    nodeX[3214] = -0.7553584144;
    nodeX[3215] = -0.1906018223;

    weight[1072] = 0.01147507935;
    nodeX[3216] = -0.6269805509;
    nodeX[3217] = 0.7553584144;
    nodeX[3218] = -0.1906018223;

    weight[1073] = 0.01147507935;
    nodeX[3219] = -0.6269805509;
    nodeX[3220] = -0.7553584144;
    nodeX[3221] = -0.1906018223;

    weight[1074] = 0.01147507935;
    nodeX[3222] = 0.1906018223;
    nodeX[3223] = 0.7553584144;
    nodeX[3224] = 0.6269805509;

    weight[1075] = 0.01147507935;
    nodeX[3225] = -0.1906018223;
    nodeX[3226] = 0.7553584144;
    nodeX[3227] = 0.6269805509;

    weight[1076] = 0.01147507935;
    nodeX[3228] = 0.1906018223;
    nodeX[3229] = -0.7553584144;
    nodeX[3230] = 0.6269805509;

    weight[1077] = 0.01147507935;
    nodeX[3231] = 0.1906018223;
    nodeX[3232] = 0.7553584144;
    nodeX[3233] = -0.6269805509;

    weight[1078] = 0.01147507935;
    nodeX[3234] = -0.1906018223;
    nodeX[3235] = -0.7553584144;
    nodeX[3236] = 0.6269805509;

    weight[1079] = 0.01147507935;
    nodeX[3237] = 0.1906018223;
    nodeX[3238] = -0.7553584144;
    nodeX[3239] = -0.6269805509;

    weight[1080] = 0.01147507935;
    nodeX[3240] = -0.1906018223;
    nodeX[3241] = 0.7553584144;
    nodeX[3242] = -0.6269805509;

    weight[1081] = 0.01147507935;
    nodeX[3243] = -0.1906018223;
    nodeX[3244] = -0.7553584144;
    nodeX[3245] = -0.6269805509;

    weight[1082] = 0.01147507935;
    nodeX[3246] = 0.1906018223;
    nodeX[3247] = 0.6269805509;
    nodeX[3248] = 0.7553584144;

    weight[1083] = 0.01147507935;
    nodeX[3249] = -0.1906018223;
    nodeX[3250] = 0.6269805509;
    nodeX[3251] = 0.7553584144;

    weight[1084] = 0.01147507935;
    nodeX[3252] = 0.1906018223;
    nodeX[3253] = -0.6269805509;
    nodeX[3254] = 0.7553584144;

    weight[1085] = 0.01147507935;
    nodeX[3255] = 0.1906018223;
    nodeX[3256] = 0.6269805509;
    nodeX[3257] = -0.7553584144;

    weight[1086] = 0.01147507935;
    nodeX[3258] = -0.1906018223;
    nodeX[3259] = -0.6269805509;
    nodeX[3260] = 0.7553584144;

    weight[1087] = 0.01147507935;
    nodeX[3261] = 0.1906018223;
    nodeX[3262] = -0.6269805509;
    nodeX[3263] = -0.7553584144;

    weight[1088] = 0.01147507935;
    nodeX[3264] = -0.1906018223;
    nodeX[3265] = 0.6269805509;
    nodeX[3266] = -0.7553584144;

    weight[1089] = 0.01147507935;
    nodeX[3267] = -0.1906018223;
    nodeX[3268] = -0.6269805509;
    nodeX[3269] = -0.7553584144;

    weight[1090] = 0.01147507935;
    nodeX[3270] = 0.7553584144;
    nodeX[3271] = 0.1906018223;
    nodeX[3272] = 0.6269805509;

    weight[1091] = 0.01147507935;
    nodeX[3273] = -0.7553584144;
    nodeX[3274] = 0.1906018223;
    nodeX[3275] = 0.6269805509;

    weight[1092] = 0.01147507935;
    nodeX[3276] = 0.7553584144;
    nodeX[3277] = -0.1906018223;
    nodeX[3278] = 0.6269805509;

    weight[1093] = 0.01147507935;
    nodeX[3279] = 0.7553584144;
    nodeX[3280] = 0.1906018223;
    nodeX[3281] = -0.6269805509;

    weight[1094] = 0.01147507935;
    nodeX[3282] = -0.7553584144;
    nodeX[3283] = -0.1906018223;
    nodeX[3284] = 0.6269805509;

    weight[1095] = 0.01147507935;
    nodeX[3285] = 0.7553584144;
    nodeX[3286] = -0.1906018223;
    nodeX[3287] = -0.6269805509;

    weight[1096] = 0.01147507935;
    nodeX[3288] = -0.7553584144;
    nodeX[3289] = 0.1906018223;
    nodeX[3290] = -0.6269805509;

    weight[1097] = 0.01147507935;
    nodeX[3291] = -0.7553584144;
    nodeX[3292] = -0.1906018223;
    nodeX[3293] = -0.6269805509;

    weight[1098] = 0.01147507935;
    nodeX[3294] = 0.6269805509;
    nodeX[3295] = 0.1906018223;
    nodeX[3296] = 0.7553584144;

    weight[1099] = 0.01147507935;
    nodeX[3297] = -0.6269805509;
    nodeX[3298] = 0.1906018223;
    nodeX[3299] = 0.7553584144;

    weight[1100] = 0.01147507935;
    nodeX[3300] = 0.6269805509;
    nodeX[3301] = -0.1906018223;
    nodeX[3302] = 0.7553584144;

    weight[1101] = 0.01147507935;
    nodeX[3303] = 0.6269805509;
    nodeX[3304] = 0.1906018223;
    nodeX[3305] = -0.7553584144;

    weight[1102] = 0.01147507935;
    nodeX[3306] = -0.6269805509;
    nodeX[3307] = -0.1906018223;
    nodeX[3308] = 0.7553584144;

    weight[1103] = 0.01147507935;
    nodeX[3309] = 0.6269805509;
    nodeX[3310] = -0.1906018223;
    nodeX[3311] = -0.7553584144;

    weight[1104] = 0.01147507935;
    nodeX[3312] = -0.6269805509;
    nodeX[3313] = 0.1906018223;
    nodeX[3314] = -0.7553584144;

    weight[1105] = 0.01147507935;
    nodeX[3315] = -0.6269805509;
    nodeX[3316] = -0.1906018223;
    nodeX[3317] = -0.7553584144;

    weight[1106] = 0.01144521609;
    nodeX[3318] = 0.7344305758;
    nodeX[3319] = 0.6031161693;
    nodeX[3320] = 0.3112275947;

    weight[1107] = 0.01144521609;
    nodeX[3321] = -0.7344305758;
    nodeX[3322] = 0.6031161693;
    nodeX[3323] = 0.3112275947;

    weight[1108] = 0.01144521609;
    nodeX[3324] = 0.7344305758;
    nodeX[3325] = -0.6031161693;
    nodeX[3326] = 0.3112275947;

    weight[1109] = 0.01144521609;
    nodeX[3327] = 0.7344305758;
    nodeX[3328] = 0.6031161693;
    nodeX[3329] = -0.3112275947;

    weight[1110] = 0.01144521609;
    nodeX[3330] = -0.7344305758;
    nodeX[3331] = -0.6031161693;
    nodeX[3332] = 0.3112275947;

    weight[1111] = 0.01144521609;
    nodeX[3333] = 0.7344305758;
    nodeX[3334] = -0.6031161693;
    nodeX[3335] = -0.3112275947;

    weight[1112] = 0.01144521609;
    nodeX[3336] = -0.7344305758;
    nodeX[3337] = 0.6031161693;
    nodeX[3338] = -0.3112275947;

    weight[1113] = 0.01144521609;
    nodeX[3339] = -0.7344305758;
    nodeX[3340] = -0.6031161693;
    nodeX[3341] = -0.3112275947;

    weight[1114] = 0.01144521609;
    nodeX[3342] = 0.6031161693;
    nodeX[3343] = 0.7344305758;
    nodeX[3344] = 0.3112275947;

    weight[1115] = 0.01144521609;
    nodeX[3345] = -0.6031161693;
    nodeX[3346] = 0.7344305758;
    nodeX[3347] = 0.3112275947;

    weight[1116] = 0.01144521609;
    nodeX[3348] = 0.6031161693;
    nodeX[3349] = -0.7344305758;
    nodeX[3350] = 0.3112275947;

    weight[1117] = 0.01144521609;
    nodeX[3351] = 0.6031161693;
    nodeX[3352] = 0.7344305758;
    nodeX[3353] = -0.3112275947;

    weight[1118] = 0.01144521609;
    nodeX[3354] = -0.6031161693;
    nodeX[3355] = -0.7344305758;
    nodeX[3356] = 0.3112275947;

    weight[1119] = 0.01144521609;
    nodeX[3357] = 0.6031161693;
    nodeX[3358] = -0.7344305758;
    nodeX[3359] = -0.3112275947;

    weight[1120] = 0.01144521609;
    nodeX[3360] = -0.6031161693;
    nodeX[3361] = 0.7344305758;
    nodeX[3362] = -0.3112275947;

    weight[1121] = 0.01144521609;
    nodeX[3363] = -0.6031161693;
    nodeX[3364] = -0.7344305758;
    nodeX[3365] = -0.3112275947;

    weight[1122] = 0.01144521609;
    nodeX[3366] = 0.3112275947;
    nodeX[3367] = 0.7344305758;
    nodeX[3368] = 0.6031161693;

    weight[1123] = 0.01144521609;
    nodeX[3369] = -0.3112275947;
    nodeX[3370] = 0.7344305758;
    nodeX[3371] = 0.6031161693;

    weight[1124] = 0.01144521609;
    nodeX[3372] = 0.3112275947;
    nodeX[3373] = -0.7344305758;
    nodeX[3374] = 0.6031161693;

    weight[1125] = 0.01144521609;
    nodeX[3375] = 0.3112275947;
    nodeX[3376] = 0.7344305758;
    nodeX[3377] = -0.6031161693;

    weight[1126] = 0.01144521609;
    nodeX[3378] = -0.3112275947;
    nodeX[3379] = -0.7344305758;
    nodeX[3380] = 0.6031161693;

    weight[1127] = 0.01144521609;
    nodeX[3381] = 0.3112275947;
    nodeX[3382] = -0.7344305758;
    nodeX[3383] = -0.6031161693;

    weight[1128] = 0.01144521609;
    nodeX[3384] = -0.3112275947;
    nodeX[3385] = 0.7344305758;
    nodeX[3386] = -0.6031161693;

    weight[1129] = 0.01144521609;
    nodeX[3387] = -0.3112275947;
    nodeX[3388] = -0.7344305758;
    nodeX[3389] = -0.6031161693;

    weight[1130] = 0.01144521609;
    nodeX[3390] = 0.3112275947;
    nodeX[3391] = 0.6031161693;
    nodeX[3392] = 0.7344305758;

    weight[1131] = 0.01144521609;
    nodeX[3393] = -0.3112275947;
    nodeX[3394] = 0.6031161693;
    nodeX[3395] = 0.7344305758;

    weight[1132] = 0.01144521609;
    nodeX[3396] = 0.3112275947;
    nodeX[3397] = -0.6031161693;
    nodeX[3398] = 0.7344305758;

    weight[1133] = 0.01144521609;
    nodeX[3399] = 0.3112275947;
    nodeX[3400] = 0.6031161693;
    nodeX[3401] = -0.7344305758;

    weight[1134] = 0.01144521609;
    nodeX[3402] = -0.3112275947;
    nodeX[3403] = -0.6031161693;
    nodeX[3404] = 0.7344305758;

    weight[1135] = 0.01144521609;
    nodeX[3405] = 0.3112275947;
    nodeX[3406] = -0.6031161693;
    nodeX[3407] = -0.7344305758;

    weight[1136] = 0.01144521609;
    nodeX[3408] = -0.3112275947;
    nodeX[3409] = 0.6031161693;
    nodeX[3410] = -0.7344305758;

    weight[1137] = 0.01144521609;
    nodeX[3411] = -0.3112275947;
    nodeX[3412] = -0.6031161693;
    nodeX[3413] = -0.7344305758;

    weight[1138] = 0.01144521609;
    nodeX[3414] = 0.7344305758;
    nodeX[3415] = 0.3112275947;
    nodeX[3416] = 0.6031161693;

    weight[1139] = 0.01144521609;
    nodeX[3417] = -0.7344305758;
    nodeX[3418] = 0.3112275947;
    nodeX[3419] = 0.6031161693;

    weight[1140] = 0.01144521609;
    nodeX[3420] = 0.7344305758;
    nodeX[3421] = -0.3112275947;
    nodeX[3422] = 0.6031161693;

    weight[1141] = 0.01144521609;
    nodeX[3423] = 0.7344305758;
    nodeX[3424] = 0.3112275947;
    nodeX[3425] = -0.6031161693;

    weight[1142] = 0.01144521609;
    nodeX[3426] = -0.7344305758;
    nodeX[3427] = -0.3112275947;
    nodeX[3428] = 0.6031161693;

    weight[1143] = 0.01144521609;
    nodeX[3429] = 0.7344305758;
    nodeX[3430] = -0.3112275947;
    nodeX[3431] = -0.6031161693;

    weight[1144] = 0.01144521609;
    nodeX[3432] = -0.7344305758;
    nodeX[3433] = 0.3112275947;
    nodeX[3434] = -0.6031161693;

    weight[1145] = 0.01144521609;
    nodeX[3435] = -0.7344305758;
    nodeX[3436] = -0.3112275947;
    nodeX[3437] = -0.6031161693;

    weight[1146] = 0.01144521609;
    nodeX[3438] = 0.6031161693;
    nodeX[3439] = 0.3112275947;
    nodeX[3440] = 0.7344305758;

    weight[1147] = 0.01144521609;
    nodeX[3441] = -0.6031161693;
    nodeX[3442] = 0.3112275947;
    nodeX[3443] = 0.7344305758;

    weight[1148] = 0.01144521609;
    nodeX[3444] = 0.6031161693;
    nodeX[3445] = -0.3112275947;
    nodeX[3446] = 0.7344305758;

    weight[1149] = 0.01144521609;
    nodeX[3447] = 0.6031161693;
    nodeX[3448] = 0.3112275947;
    nodeX[3449] = -0.7344305758;

    weight[1150] = 0.01144521609;
    nodeX[3450] = -0.6031161693;
    nodeX[3451] = -0.3112275947;
    nodeX[3452] = 0.7344305758;

    weight[1151] = 0.01144521609;
    nodeX[3453] = 0.6031161693;
    nodeX[3454] = -0.3112275947;
    nodeX[3455] = -0.7344305758;

    weight[1152] = 0.01144521609;
    nodeX[3456] = -0.6031161693;
    nodeX[3457] = 0.3112275947;
    nodeX[3458] = -0.7344305758;

    weight[1153] = 0.01144521609;
    nodeX[3459] = -0.6031161693;
    nodeX[3460] = -0.3112275947;
    nodeX[3461] = -0.7344305758;

    weight[1154] = 0.01144263581;
    nodeX[3462] = 0.7043837184;
    nodeX[3463] = 0.5693702498;
    nodeX[3464] = 0.4238644782;

    weight[1155] = 0.01144263581;
    nodeX[3465] = -0.7043837184;
    nodeX[3466] = 0.5693702498;
    nodeX[3467] = 0.4238644782;

    weight[1156] = 0.01144263581;
    nodeX[3468] = 0.7043837184;
    nodeX[3469] = -0.5693702498;
    nodeX[3470] = 0.4238644782;

    weight[1157] = 0.01144263581;
    nodeX[3471] = 0.7043837184;
    nodeX[3472] = 0.5693702498;
    nodeX[3473] = -0.4238644782;

    weight[1158] = 0.01144263581;
    nodeX[3474] = -0.7043837184;
    nodeX[3475] = -0.5693702498;
    nodeX[3476] = 0.4238644782;

    weight[1159] = 0.01144263581;
    nodeX[3477] = 0.7043837184;
    nodeX[3478] = -0.5693702498;
    nodeX[3479] = -0.4238644782;

    weight[1160] = 0.01144263581;
    nodeX[3480] = -0.7043837184;
    nodeX[3481] = 0.5693702498;
    nodeX[3482] = -0.4238644782;

    weight[1161] = 0.01144263581;
    nodeX[3483] = -0.7043837184;
    nodeX[3484] = -0.5693702498;
    nodeX[3485] = -0.4238644782;

    weight[1162] = 0.01144263581;
    nodeX[3486] = 0.5693702498;
    nodeX[3487] = 0.7043837184;
    nodeX[3488] = 0.4238644782;

    weight[1163] = 0.01144263581;
    nodeX[3489] = -0.5693702498;
    nodeX[3490] = 0.7043837184;
    nodeX[3491] = 0.4238644782;

    weight[1164] = 0.01144263581;
    nodeX[3492] = 0.5693702498;
    nodeX[3493] = -0.7043837184;
    nodeX[3494] = 0.4238644782;

    weight[1165] = 0.01144263581;
    nodeX[3495] = 0.5693702498;
    nodeX[3496] = 0.7043837184;
    nodeX[3497] = -0.4238644782;

    weight[1166] = 0.01144263581;
    nodeX[3498] = -0.5693702498;
    nodeX[3499] = -0.7043837184;
    nodeX[3500] = 0.4238644782;

    weight[1167] = 0.01144263581;
    nodeX[3501] = 0.5693702498;
    nodeX[3502] = -0.7043837184;
    nodeX[3503] = -0.4238644782;

    weight[1168] = 0.01144263581;
    nodeX[3504] = -0.5693702498;
    nodeX[3505] = 0.7043837184;
    nodeX[3506] = -0.4238644782;

    weight[1169] = 0.01144263581;
    nodeX[3507] = -0.5693702498;
    nodeX[3508] = -0.7043837184;
    nodeX[3509] = -0.4238644782;

    weight[1170] = 0.01144263581;
    nodeX[3510] = 0.4238644782;
    nodeX[3511] = 0.7043837184;
    nodeX[3512] = 0.5693702498;

    weight[1171] = 0.01144263581;
    nodeX[3513] = -0.4238644782;
    nodeX[3514] = 0.7043837184;
    nodeX[3515] = 0.5693702498;

    weight[1172] = 0.01144263581;
    nodeX[3516] = 0.4238644782;
    nodeX[3517] = -0.7043837184;
    nodeX[3518] = 0.5693702498;

    weight[1173] = 0.01144263581;
    nodeX[3519] = 0.4238644782;
    nodeX[3520] = 0.7043837184;
    nodeX[3521] = -0.5693702498;

    weight[1174] = 0.01144263581;
    nodeX[3522] = -0.4238644782;
    nodeX[3523] = -0.7043837184;
    nodeX[3524] = 0.5693702498;

    weight[1175] = 0.01144263581;
    nodeX[3525] = 0.4238644782;
    nodeX[3526] = -0.7043837184;
    nodeX[3527] = -0.5693702498;

    weight[1176] = 0.01144263581;
    nodeX[3528] = -0.4238644782;
    nodeX[3529] = 0.7043837184;
    nodeX[3530] = -0.5693702498;

    weight[1177] = 0.01144263581;
    nodeX[3531] = -0.4238644782;
    nodeX[3532] = -0.7043837184;
    nodeX[3533] = -0.5693702498;

    weight[1178] = 0.01144263581;
    nodeX[3534] = 0.4238644782;
    nodeX[3535] = 0.5693702498;
    nodeX[3536] = 0.7043837184;

    weight[1179] = 0.01144263581;
    nodeX[3537] = -0.4238644782;
    nodeX[3538] = 0.5693702498;
    nodeX[3539] = 0.7043837184;

    weight[1180] = 0.01144263581;
    nodeX[3540] = 0.4238644782;
    nodeX[3541] = -0.5693702498;
    nodeX[3542] = 0.7043837184;

    weight[1181] = 0.01144263581;
    nodeX[3543] = 0.4238644782;
    nodeX[3544] = 0.5693702498;
    nodeX[3545] = -0.7043837184;

    weight[1182] = 0.01144263581;
    nodeX[3546] = -0.4238644782;
    nodeX[3547] = -0.5693702498;
    nodeX[3548] = 0.7043837184;

    weight[1183] = 0.01144263581;
    nodeX[3549] = 0.4238644782;
    nodeX[3550] = -0.5693702498;
    nodeX[3551] = -0.7043837184;

    weight[1184] = 0.01144263581;
    nodeX[3552] = -0.4238644782;
    nodeX[3553] = 0.5693702498;
    nodeX[3554] = -0.7043837184;

    weight[1185] = 0.01144263581;
    nodeX[3555] = -0.4238644782;
    nodeX[3556] = -0.5693702498;
    nodeX[3557] = -0.7043837184;

    weight[1186] = 0.01144263581;
    nodeX[3558] = 0.7043837184;
    nodeX[3559] = 0.4238644782;
    nodeX[3560] = 0.5693702498;

    weight[1187] = 0.01144263581;
    nodeX[3561] = -0.7043837184;
    nodeX[3562] = 0.4238644782;
    nodeX[3563] = 0.5693702498;

    weight[1188] = 0.01144263581;
    nodeX[3564] = 0.7043837184;
    nodeX[3565] = -0.4238644782;
    nodeX[3566] = 0.5693702498;

    weight[1189] = 0.01144263581;
    nodeX[3567] = 0.7043837184;
    nodeX[3568] = 0.4238644782;
    nodeX[3569] = -0.5693702498;

    weight[1190] = 0.01144263581;
    nodeX[3570] = -0.7043837184;
    nodeX[3571] = -0.4238644782;
    nodeX[3572] = 0.5693702498;

    weight[1191] = 0.01144263581;
    nodeX[3573] = 0.7043837184;
    nodeX[3574] = -0.4238644782;
    nodeX[3575] = -0.5693702498;

    weight[1192] = 0.01144263581;
    nodeX[3576] = -0.7043837184;
    nodeX[3577] = 0.4238644782;
    nodeX[3578] = -0.5693702498;

    weight[1193] = 0.01144263581;
    nodeX[3579] = -0.7043837184;
    nodeX[3580] = -0.4238644782;
    nodeX[3581] = -0.5693702498;

    weight[1194] = 0.01144263581;
    nodeX[3582] = 0.5693702498;
    nodeX[3583] = 0.4238644782;
    nodeX[3584] = 0.7043837184;

    weight[1195] = 0.01144263581;
    nodeX[3585] = -0.5693702498;
    nodeX[3586] = 0.4238644782;
    nodeX[3587] = 0.7043837184;

    weight[1196] = 0.01144263581;
    nodeX[3588] = 0.5693702498;
    nodeX[3589] = -0.4238644782;
    nodeX[3590] = 0.7043837184;

    weight[1197] = 0.01144263581;
    nodeX[3591] = 0.5693702498;
    nodeX[3592] = 0.4238644782;
    nodeX[3593] = -0.7043837184;

    weight[1198] = 0.01144263581;
    nodeX[3594] = -0.5693702498;
    nodeX[3595] = -0.4238644782;
    nodeX[3596] = 0.7043837184;

    weight[1199] = 0.01144263581;
    nodeX[3597] = 0.5693702498;
    nodeX[3598] = -0.4238644782;
    nodeX[3599] = -0.7043837184;

    weight[1200] = 0.01144263581;
    nodeX[3600] = -0.5693702498;
    nodeX[3601] = 0.4238644782;
    nodeX[3602] = -0.7043837184;

    weight[1201] = 0.01144263581;
    nodeX[3603] = -0.5693702498;
    nodeX[3604] = -0.4238644782;
    nodeX[3605] = -0.7043837184;

    break;

  default:

    stringstream message;
    message << "The number of quadrature nodes requested" << endl;
    message << "is not yet supported."                    << endl;
    message << "numQuadNodes = " << numQuadNodes          << endl;
    message << "Try values 14, 110, 1202."                << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);

  } /* end of switch */

  /*****************************************************/

}


/* == Compute the force density arising from the control points for FFTW3 uniform mesh. */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_applyControlPtsForceToMesh_FFTW3(int                        num_dim,
                                                                             int                        numMeshPtsPerDir[3],
                                                                             double                     meshDeltaX,
                                                                             double                    *meshCenterX0,
                                                                             fftw_complex              *f_m[3],
                                                                             SELM_Lagrangian_CONTROLPTS_BASIC1    *SELM_LagrangianData_CONTROLPTS_BASIC1) {

  const char *error_str_func = "IB_appl1_applyControlPtsForceToMesh_FFTW3()";

  int i, j, k, d;

  int i1, i2, i3;
  int j1, j2, j3;

  int ii1, ii2, ii3;
  int jj1, jj2, jj3;

  double L;
  int c;

  int N;
  int I, I1, I2;

  int meshI;

  int meshJ[3];

  double meshBaseX0[3];
  int meshI0[3];

  int xMinusIndex;
  int xPlusIndex;

  int yMinusIndex;
  int yPlusIndex;

  int zMinusIndex;
  int zPlusIndex;

  int tensor_cell_order;
  int tensor_cell_size;

  int tensor_face_order;
  int tensor_face_size;

  int localMeshSize;

  double *meshX;

  double  X[3];

  double meshX0[3];
  double *meshForce;
  double controlPtX[3];
  double controlPtF[3];

  double periodBaseX0[3];
  double periodDiffX0[3];
  double baseX0[3];

  int      tensorDataIndex;
  double  *userTensorData;

  int      flagFoundItem;
  int      periodI;

  double **ptForce_ptr;
  double   ptForce[3];

  int      numQuadNodes;
  double   sphereR;

  int      kernalMeshwidths;


  int      flagDebugSELM;

  int      numMeshSELMIndices;
  int     *meshSELMIndices;
  int     *meshSELMJ;

  double  *meshSELMX;

  double   meshFactor;

  int      index;
  int      numDir;

  int      meshSELMForceDensity_num_dim;
  double  *meshSELMForceDensity;
  double  *meshSELMForceDensity_ptr;

  int      flagSELM_Translational_forceDensity;
  int      flagSELM_Rotational_forceDensity;

  int      flagDegFreedomType = -1;
  double   X_cm[3];
  double   F_cm[3];

  double   Theta[3];
  double   T_ang[3];

  double   R_0;
  double  *weightFuncVals;
  double   weightFuncSum;

  operatorDataType_T_KERNEL_1 *opData;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmExtras;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaExtras;

  IB_appl1_computeSmoothedForceFieldExtrasType
    *IB_appl1_computeSmoothedForceFieldExtras;

  if (num_dim > 3) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("  Only dimensions <= 3 are implemented \n");
    printf("  num_dim = %d \n", num_dim);
    packageError(1, 0);
  }

  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1:

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    /* -- Loop over the control points and apply the forces to the mesh
     *    via the approximate Diract delta functions.
     */
    flagDegFreedomType = -1;
    for (k = 0; k < SELM_LagrangianData_CONTROLPTS_BASIC1->numControlPts; k++) {

      /* -- Get control point information */
      for (d = 0; d < num_dim; d++) {
        X_cm[d] = SELM_LagrangianData_CONTROLPTS_BASIC1->pt_X[(k * num_dim) + d];
        F_cm[d] = SELM_LagrangianData_CONTROLPTS_BASIC1->pt_Force[(k * num_dim) + d];
      }

      /* == Compute the force density. */

      /* -- Compute the force density associated with the surface
       *    of a sphere and the interpolations used. */

      /* -- determine the mesh indices */

      /* perform index arithmetic to find a mesh point X0 near X */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
        meshI0[d]      = floor(((X_cm[d] - meshBaseX0[d])/meshDeltaX) + 0.5);
        meshX0[d]      = meshBaseX0[d] + meshI0[d]*meshDeltaX;
      } /* end of d loop */

      /* now determine the indices for all mesh points with square
       * patch of size sphereR + kernalMeshwidths*meshDeltaX (depends on smoother used)
       */
      R_0                  = opData->weightTable->R_0;
      numDir               = 2*ceil(R_0) + 1;
      numMeshSELMIndices   = numDir*numDir*numDir;
      meshSELMIndices      = (int *)malloc(sizeof(int)*numMeshSELMIndices);
      meshSELMJ            = (int *)malloc(sizeof(int)*numMeshSELMIndices*num_dim);

      meshSELMX            = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);
      meshSELMForceDensity = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);

      /* compute base value for the mesh points */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
      }

      for (j1 = 0; j1 < numDir; j1++) {
        for (j2 = 0; j2 < numDir; j2++) {
          for (j3 = 0; j3 < numDir; j3++) {

            index    = j1 + j2*numDir + j3*numDir*numDir;

            meshJ[0] = meshI0[0] + j1 - ((numDir - 1)/2);
            meshJ[1] = meshI0[1] + j2 - ((numDir - 1)/2);
            meshJ[2] = meshI0[2] + j3 - ((numDir - 1)/2);

            /* compute the index modulo numMeshPtsPerDir */
            for (d = 0; d < num_dim; d++) {

              /* compute the mesh point locations before the adjustments are made to make it mod L */
              meshSELMX[index*num_dim + d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;

              /* make the indices mod N */
              meshSELMJ[index*num_dim + d] = (meshJ[d] % numMeshPtsPerDir[d]);
              if (meshSELMJ[index*num_dim + d] < 0) {
                meshSELMJ[index*num_dim + d] = meshSELMJ[index*num_dim + d] + numMeshPtsPerDir[d];
              }

            } /* end of d loop */

            /* compute the index on the mesh */
            meshI                    = meshSELMJ[index*num_dim + 0]
                                                 + meshSELMJ[index*num_dim + 1]*numMeshPtsPerDir[0]
                                                                                                 + meshSELMJ[index*num_dim + 2]*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];
            meshSELMIndices[index]   = meshI;

          } /* end of j3 loop */
        } /* end of j2 loop */
      } /* end of j1 loop */

      /* -- compute the weights on the patch of mesh points */
      weightFuncVals = NULL;
      weightFromTable(num_dim, numMeshSELMIndices, meshSELMX,
                                      X_cm, meshDeltaX,
                                      opData->weightTable,
                                      &weightFuncVals);

      weightFuncSum = 0;
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncSum += weightFuncVals[I];
      }

      /* renormalize the weight values, since they should sum to one */
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncVals[I] = weightFuncVals[I]/weightFuncSum;
      }

      /* -- construct the force density by looping over the mesh points */
      meshFactor = 1.0/(meshDeltaX*meshDeltaX*meshDeltaX); /* weights include meshDeltaX term implicitly */
      for (I = 0; I < numMeshSELMIndices; I++) {

        meshI = meshSELMIndices[I];

        for (d = 0; d < num_dim; d++) {
          meshSELMForceDensity[I*num_dim + d] = F_cm[d]*weightFuncVals[I]*meshFactor;
        }

      } /* end of I loop */

      /* -- apply the force to the mesh */
      for (I = 0; I < numMeshSELMIndices; I++) {
        meshI = meshSELMIndices[I];
        for (d = 0; d < num_dim; d++) {
          f_m[d][meshI][0] += meshSELMForceDensity[I*num_dim + d];
          f_m[d][meshI][1] += 0;
        } /* end of d loop */
      } /* end of I loop */

      /* -- debug test for the force density */
      flagDebugSELM = 0;
      if (flagDebugSELM) {
        for (j1 = 0; j1 < numMeshPtsPerDir[0]; j1++) {
          for (j2 = 0; j2 < numMeshPtsPerDir[1]; j2++) {
            for (j3 = 0; j3 < numMeshPtsPerDir[2]; j3++) {

              I = j1 + j2*numMeshPtsPerDir[0]
                                           + j3*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];

              meshJ[0] = j1;
              meshJ[1] = j2;
              meshJ[2] = j3;

              for (d = 0; d < num_dim; d++) {
                X[d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;
              } /* end of d */

              for (d = 0; d < num_dim; d++) {
                if (d == 0) {
                  f_m[d][I][0] = X[d]; /* meshSELMForceDensity[I*num_dim + d]; */
                  f_m[d][I][1] = 0;
                } else {
                  f_m[d][I][0] += 0;
                  f_m[d][I][1] += 0;
                }
              } /* end of d */

            } /* end of j3 */
          } /* end of j2 */
        } /* end of j1 */
      }

      flagDebugSELM = 0;
      if (flagDebugSELM) {
        printf("WARNING: %s : %s \n", error_str_code, error_str_func);
        printf("Debug information being printed. \n");
        printf(" \n");
        fflush(stdout);

        printf("meshSELMX_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMX[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

        printf("meshSELMForceDensity_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMForceDensity[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

      } /* end of flagDebugSELM */

      /* free dynamic memory */
      free(meshSELMX);
      free(meshSELMForceDensity);

      free(meshSELMIndices);
      free(meshSELMJ);

      free(weightFuncVals);

    } /* end of k loop */

    break;

  default:
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("This operator type is not yet supported in the implementation. \n");
    printf("  operatorTypeStr = %s \n", operatorTypeStr);
    packageError(1, 0);

  } /* end of switch */


}



/* == Compute the force density arising from the control points for FFTW3 uniform mesh. */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_applyControlPtsForceToMesh_FFTW3(int                        num_dim,
                                                                             int                        numMeshPtsPerDir[3],
                                                                             double                     meshDeltaX,
                                                                             double                    *meshCenterX0,
                                                                             fftw_complex              *f_m[3],
                                                                             SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE    *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE) {

  const char *error_str_func = "IB_appl1_applyControlPtsForceToMesh_FFTW3()";

  int i, j, k, d;

  int i1, i2, i3;
  int j1, j2, j3;

  int ii1, ii2, ii3;
  int jj1, jj2, jj3;

  double L;
  int c;

  int N;
  int I, I1, I2;

  int meshI;

  int meshJ[3];

  double meshBaseX0[3];
  int meshI0[3];

  int xMinusIndex;
  int xPlusIndex;

  int yMinusIndex;
  int yPlusIndex;

  int zMinusIndex;
  int zPlusIndex;

  int tensor_cell_order;
  int tensor_cell_size;

  int tensor_face_order;
  int tensor_face_size;

  int localMeshSize;

  double *meshX;

  double  X[3];

  double meshX0[3];
  double *meshForce;
  double controlPtX[3];
  double controlPtF[3];

  double periodBaseX0[3];
  double periodDiffX0[3];
  double baseX0[3];

  int      tensorDataIndex;
  double  *userTensorData;

  int      flagFoundItem;
  int      periodI;

  double **ptForce_ptr;
  double   ptForce[3];

  int      numQuadNodes;
  double   sphereR;

  int      kernalMeshwidths;


  int      flagDebugSELM;

  int      numMeshSELMIndices;
  int     *meshSELMIndices;
  int     *meshSELMJ;

  double  *meshSELMX;

  double   meshFactor;

  int      index;
  int      numDir;

  int      meshSELMForceDensity_num_dim;
  double  *meshSELMForceDensity;
  double  *meshSELMForceDensity_ptr;

  int      flagSELM_Translational_forceDensity;
  int      flagSELM_Rotational_forceDensity;

  int      flagDegFreedomType = -1;
  double   X_cm[3];
  double   F_cm[3];

  double   Theta[3];
  double   T_ang[3];

  double   R_0;
  double  *weightFuncVals;
  double   weightFuncSum;

  operatorDataType_T_KERNEL_1 *opData;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmExtras;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaExtras;

  IB_appl1_computeSmoothedForceFieldExtrasType
    *IB_appl1_computeSmoothedForceFieldExtras;

  if (num_dim > 3) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("  Only dimensions <= 3 are implemented \n");
    printf("  num_dim = %d \n", num_dim);
    packageError(1, 0);
  }

  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1:

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    /* -- Loop over the control points and apply the forces to the mesh
     *    via the approximate Dirac-delta-functions.
     */
    flagDegFreedomType = -1;
    for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->numControlPts; k++) {

      /* -- Get control point information */
      for (d = 0; d < num_dim; d++) {
        X_cm[d] = SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->ptsX[(k * num_dim) + d];
        F_cm[d] = SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE->pt_Force[(k * num_dim) + d];
      }

      /* == Compute the force density. */

      /* -- Compute the force density associated with the surface
       *    of a sphere and the interpolations used. */

      /* -- determine the mesh indices */

      /* perform index arithmetic to find a mesh point X0 near X */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
        meshI0[d]      = floor(((X_cm[d] - meshBaseX0[d])/meshDeltaX) + 0.5);
        meshX0[d]      = meshBaseX0[d] + meshI0[d]*meshDeltaX;
      } /* end of d loop */

      /* now determine the indices for all mesh points with square
       * patch of size sphereR + kernalMeshwidths*meshDeltaX (depends on smoother used)
       */
      R_0                  = opData->weightTable->R_0;
      numDir               = 2*ceil(R_0) + 1;
      numMeshSELMIndices   = numDir*numDir*numDir;
      meshSELMIndices      = (int *)malloc(sizeof(int)*numMeshSELMIndices);
      meshSELMJ            = (int *)malloc(sizeof(int)*numMeshSELMIndices*num_dim);

      meshSELMX            = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);
      meshSELMForceDensity = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);

      /* compute base value for the mesh points */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
      }

      for (j1 = 0; j1 < numDir; j1++) {
        for (j2 = 0; j2 < numDir; j2++) {
          for (j3 = 0; j3 < numDir; j3++) {

            index    = j1 + j2*numDir + j3*numDir*numDir;

            meshJ[0] = meshI0[0] + j1 - ((numDir - 1)/2);
            meshJ[1] = meshI0[1] + j2 - ((numDir - 1)/2);
            meshJ[2] = meshI0[2] + j3 - ((numDir - 1)/2);

            /* compute the index modulo numMeshPtsPerDir */
            for (d = 0; d < num_dim; d++) {

              /* compute the mesh point locations before the adjustments are made to make it mod L */
              meshSELMX[index*num_dim + d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;

              /* make the indices mod N */
              meshSELMJ[index*num_dim + d] = (meshJ[d] % numMeshPtsPerDir[d]);
              if (meshSELMJ[index*num_dim + d] < 0) {
                meshSELMJ[index*num_dim + d] = meshSELMJ[index*num_dim + d] + numMeshPtsPerDir[d];
              }

            } /* end of d loop */

            /* compute the index on the mesh */
            meshI                    = meshSELMJ[index*num_dim + 0]
                                     + meshSELMJ[index*num_dim + 1]*numMeshPtsPerDir[0]
                                     + meshSELMJ[index*num_dim + 2]*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];
            meshSELMIndices[index]   = meshI;

          } /* end of j3 loop */
        } /* end of j2 loop */
      } /* end of j1 loop */

      /* -- compute the weights on the patch of mesh points */
      weightFuncVals = NULL;
      weightFromTable(num_dim, numMeshSELMIndices, meshSELMX,
                      X_cm, meshDeltaX,
                      opData->weightTable,
                      &weightFuncVals);

      weightFuncSum = 0;
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncSum += weightFuncVals[I];
      }

      /* re-normalize the weight values, since they should sum to one */
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncVals[I] = weightFuncVals[I]/weightFuncSum;
      }

      /* -- construct the force density by looping over the mesh points */
      meshFactor = 1.0/(meshDeltaX*meshDeltaX*meshDeltaX); /* weights include meshDeltaX term implicitly */
      for (I = 0; I < numMeshSELMIndices; I++) {

        meshI = meshSELMIndices[I];

        for (d = 0; d < num_dim; d++) {
          meshSELMForceDensity[I*num_dim + d] = F_cm[d]*weightFuncVals[I]*meshFactor;
        }

      } /* end of I loop */

      /* -- apply the force to the mesh */
      for (I = 0; I < numMeshSELMIndices; I++) {
        meshI = meshSELMIndices[I];
        for (d = 0; d < num_dim; d++) {
          f_m[d][meshI][0] += meshSELMForceDensity[I*num_dim + d];
          f_m[d][meshI][1] += 0;
        } /* end of d loop */
      } /* end of I loop */

      /* -- debug test for the force density */
      flagDebugSELM = 0;
      if (flagDebugSELM) {
        for (j1 = 0; j1 < numMeshPtsPerDir[0]; j1++) {
          for (j2 = 0; j2 < numMeshPtsPerDir[1]; j2++) {
            for (j3 = 0; j3 < numMeshPtsPerDir[2]; j3++) {

              I = j1 + j2*numMeshPtsPerDir[0]
                                           + j3*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];

              meshJ[0] = j1;
              meshJ[1] = j2;
              meshJ[2] = j3;

              for (d = 0; d < num_dim; d++) {
                X[d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;
              } /* end of d */

              for (d = 0; d < num_dim; d++) {
                if (d == 0) {
                  f_m[d][I][0] = X[d]; /* meshSELMForceDensity[I*num_dim + d]; */
                  f_m[d][I][1] = 0;
                } else {
                  f_m[d][I][0] += 0;
                  f_m[d][I][1] += 0;
                }
              } /* end of d */

            } /* end of j3 */
          } /* end of j2 */
        } /* end of j1 */
      }

      flagDebugSELM = 0;
      if (flagDebugSELM) {
        printf("WARNING: %s : %s \n", error_str_code, error_str_func);
        printf("Debug information being printed. \n");
        printf(" \n");
        fflush(stdout);

        printf("meshSELMX_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMX[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

        printf("meshSELMForceDensity_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMForceDensity[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

      } /* end of flagDebugSELM */

      /* free dynamic memory */
      free(meshSELMX);
      free(meshSELMForceDensity);

      free(meshSELMIndices);
      free(meshSELMJ);

      free(weightFuncVals);

    } /* end of k loop */

    break;

  default:
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("This operator type is not yet supported in the implementation. \n");
    printf("  operatorTypeStr = %s \n", operatorTypeStr);
    packageError(1, 0);

  } /* end of switch */


}



/* == Compute the force density arising from the control points for FFTW3 uniform mesh. */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_applyControlPtsForceToMesh_FFTW3(int                        num_dim,
                                                                             int                        numMeshPtsPerDir[3],
                                                                             double                     meshDeltaX,
                                                                             double                    *meshCenterX0,
                                                                             fftw_complex              *f_m[3],
                                                                             SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE    *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE) {

  const char *error_str_func = "IB_appl1_applyControlPtsForceToMesh_FFTW3(*,LAMMPS_HYBRID_CHARGE_ANGLE_STYLE)";

  int i, j, k, d;

  int i1, i2, i3;
  int j1, j2, j3;

  int ii1, ii2, ii3;
  int jj1, jj2, jj3;

  double L;
  int c;

  int N;
  int I, I1, I2;

  int meshI;

  int meshJ[3];

  double meshBaseX0[3];
  int meshI0[3];

  int xMinusIndex;
  int xPlusIndex;

  int yMinusIndex;
  int yPlusIndex;

  int zMinusIndex;
  int zPlusIndex;

  int tensor_cell_order;
  int tensor_cell_size;

  int tensor_face_order;
  int tensor_face_size;

  int localMeshSize;

  double *meshX;

  double  X[3];

  double meshX0[3];
  double *meshForce;
  double controlPtX[3];
  double controlPtF[3];

  double periodBaseX0[3];
  double periodDiffX0[3];
  double baseX0[3];

  int      tensorDataIndex;
  double  *userTensorData;

  int      flagFoundItem;
  int      periodI;

  double **ptForce_ptr;
  double   ptForce[3];

  int      numQuadNodes;
  double   sphereR;

  int      kernalMeshwidths;


  int      flagDebugSELM;

  int      numMeshSELMIndices;
  int     *meshSELMIndices;
  int     *meshSELMJ;

  double  *meshSELMX;

  double   meshFactor;

  int      index;
  int      numDir;

  int      meshSELMForceDensity_num_dim;
  double  *meshSELMForceDensity;
  double  *meshSELMForceDensity_ptr;

  int      flagSELM_Translational_forceDensity;
  int      flagSELM_Rotational_forceDensity;

  int      flagDegFreedomType = -1;
  double   X_cm[3];
  double   F_cm[3];

  double   Theta[3];
  double   T_ang[3];

  double   R_0;
  double  *weightFuncVals;
  double   weightFuncSum;

  operatorDataType_T_KERNEL_1 *opData;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmExtras;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaExtras;

  IB_appl1_computeSmoothedForceFieldExtrasType
    *IB_appl1_computeSmoothedForceFieldExtras;

  if (num_dim > 3) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("  Only dimensions <= 3 are implemented \n");
    printf("  num_dim = %d \n", num_dim);
    packageError(1, 0);
  }

  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1:

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    /* -- Loop over the control points and apply the forces to the mesh
     *    via the approximate Dirac-delta-functions.
     */
    flagDegFreedomType = -1;
    for (k = 0; k < SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->numControlPts; k++) {

      /* -- Get control point information */
      for (d = 0; d < num_dim; d++) {
        X_cm[d] = SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->ptsX[(k * num_dim) + d];
        F_cm[d] = SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE->pt_Force[(k * num_dim) + d];
      }

      /* == Compute the force density. */

      /* -- Compute the force density associated with the surface
       *    of a sphere and the interpolations used. */

      /* -- determine the mesh indices */

      /* perform index arithmetic to find a mesh point X0 near X */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
        meshI0[d]      = floor(((X_cm[d] - meshBaseX0[d])/meshDeltaX) + 0.5);
        meshX0[d]      = meshBaseX0[d] + meshI0[d]*meshDeltaX;
      } /* end of d loop */

      /* now determine the indices for all mesh points with square
       * patch of size sphereR + kernalMeshwidths*meshDeltaX (depends on smoother used)
       */
      R_0                  = opData->weightTable->R_0;
      numDir               = 2*ceil(R_0) + 1;
      numMeshSELMIndices   = numDir*numDir*numDir;
      meshSELMIndices      = (int *)malloc(sizeof(int)*numMeshSELMIndices);
      meshSELMJ            = (int *)malloc(sizeof(int)*numMeshSELMIndices*num_dim);

      meshSELMX            = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);
      meshSELMForceDensity = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);

      /* compute base value for the mesh points */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
      }

      for (j1 = 0; j1 < numDir; j1++) {
        for (j2 = 0; j2 < numDir; j2++) {
          for (j3 = 0; j3 < numDir; j3++) {

            index    = j1 + j2*numDir + j3*numDir*numDir;

            meshJ[0] = meshI0[0] + j1 - ((numDir - 1)/2);
            meshJ[1] = meshI0[1] + j2 - ((numDir - 1)/2);
            meshJ[2] = meshI0[2] + j3 - ((numDir - 1)/2);

            /* compute the index modulo numMeshPtsPerDir */
            for (d = 0; d < num_dim; d++) {

              /* compute the mesh point locations before the adjustments are made to make it mod L */
              meshSELMX[index*num_dim + d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;

              /* make the indices mod N */
              meshSELMJ[index*num_dim + d] = (meshJ[d] % numMeshPtsPerDir[d]);
              if (meshSELMJ[index*num_dim + d] < 0) {
                meshSELMJ[index*num_dim + d] = meshSELMJ[index*num_dim + d] + numMeshPtsPerDir[d];
              }

            } /* end of d loop */

            /* compute the index on the mesh */
            meshI                    = meshSELMJ[index*num_dim + 0]
                                     + meshSELMJ[index*num_dim + 1]*numMeshPtsPerDir[0]
                                     + meshSELMJ[index*num_dim + 2]*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];
            meshSELMIndices[index]   = meshI;

          } /* end of j3 loop */
        } /* end of j2 loop */
      } /* end of j1 loop */

      /* -- compute the weights on the patch of mesh points */
      weightFuncVals = NULL;
      weightFromTable(num_dim, numMeshSELMIndices, meshSELMX,
                      X_cm, meshDeltaX,
                      opData->weightTable,
                      &weightFuncVals);

      weightFuncSum = 0;
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncSum += weightFuncVals[I];
      }

      /* re-normalize the weight values, since they should sum to one */
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncVals[I] = weightFuncVals[I]/weightFuncSum;
      }

      /* -- construct the force density by looping over the mesh points */
      meshFactor = 1.0/(meshDeltaX*meshDeltaX*meshDeltaX); /* weights include meshDeltaX term implicitly */
      for (I = 0; I < numMeshSELMIndices; I++) {

        meshI = meshSELMIndices[I];

        for (d = 0; d < num_dim; d++) {
          meshSELMForceDensity[I*num_dim + d] = F_cm[d]*weightFuncVals[I]*meshFactor;
        }

      } /* end of I loop */

      /* -- apply the force to the mesh */
      for (I = 0; I < numMeshSELMIndices; I++) {
        meshI = meshSELMIndices[I];
        for (d = 0; d < num_dim; d++) {
          f_m[d][meshI][0] += meshSELMForceDensity[I*num_dim + d];
          f_m[d][meshI][1] += 0;
        } /* end of d loop */
      } /* end of I loop */

      /* -- debug test for the force density */
      flagDebugSELM = 0;
      if (flagDebugSELM) {
        for (j1 = 0; j1 < numMeshPtsPerDir[0]; j1++) {
          for (j2 = 0; j2 < numMeshPtsPerDir[1]; j2++) {
            for (j3 = 0; j3 < numMeshPtsPerDir[2]; j3++) {

              I = j1 + j2*numMeshPtsPerDir[0]
                                           + j3*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];

              meshJ[0] = j1;
              meshJ[1] = j2;
              meshJ[2] = j3;

              for (d = 0; d < num_dim; d++) {
                X[d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;
              } /* end of d */

              for (d = 0; d < num_dim; d++) {
                if (d == 0) {
                  f_m[d][I][0] = X[d]; /* meshSELMForceDensity[I*num_dim + d]; */
                  f_m[d][I][1] = 0;
                } else {
                  f_m[d][I][0] += 0;
                  f_m[d][I][1] += 0;
                }
              } /* end of d */

            } /* end of j3 */
          } /* end of j2 */
        } /* end of j1 */
      }

      flagDebugSELM = 0;
      if (flagDebugSELM) {
        printf("WARNING: %s : %s \n", error_str_code, error_str_func);
        printf("Debug information being printed. \n");
        printf(" \n");
        fflush(stdout);

        printf("meshSELMX_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMX[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

        printf("meshSELMForceDensity_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMForceDensity[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

      } /* end of flagDebugSELM */

      /* free dynamic memory */
      free(meshSELMX);
      free(meshSELMForceDensity);

      free(meshSELMIndices);
      free(meshSELMJ);

      free(weightFuncVals);

    } /* end of k loop */

    break;

  default:
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("This operator type is not yet supported in the implementation. \n");
    printf("  operatorTypeStr = %s \n", operatorTypeStr);
    packageError(1, 0);

  } /* end of switch */


}



/* == Compute the force density arising from the control points for FFTW3 uniform mesh. */
void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::IB_appl1_applyControlPtsForceToMesh_FFTW3(int                        num_dim,
                                                                             int                        numMeshPtsPerDir[3],
                                                                             double                     meshDeltaX,
                                                                             double                    *meshCenterX0,
                                                                             fftw_complex              *f_m[3],
                                                                             SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID    *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID) {

  const char *error_str_func = "IB_appl1_applyControlPtsForceToMesh_FFTW3(*,LAMMPS_ATOM_STYLE_ELLIPSOID)";

  int i, j, k, d;

  int i1, i2, i3;
  int j1, j2, j3;

  int ii1, ii2, ii3;
  int jj1, jj2, jj3;

  double L;
  int c;

  int N;
  int I, I1, I2;

  int meshI;

  int meshJ[3];

  double meshBaseX0[3];
  int meshI0[3];

  int xMinusIndex;
  int xPlusIndex;

  int yMinusIndex;
  int yPlusIndex;

  int zMinusIndex;
  int zPlusIndex;

  int tensor_cell_order;
  int tensor_cell_size;

  int tensor_face_order;
  int tensor_face_size;

  int localMeshSize;

  double *meshX;

  double  X[3];

  double meshX0[3];
  double *meshForce;
  double controlPtX[3];
  double controlPtF[3];

  double periodBaseX0[3];
  double periodDiffX0[3];
  double baseX0[3];

  int      tensorDataIndex;
  double  *userTensorData;

  int      flagFoundItem;
  int      periodI;

  double **ptForce_ptr;
  double   ptForce[3];

  int      numQuadNodes;
  double   sphereR;

  int      kernalMeshwidths;


  int      flagDebugSELM;

  int      numMeshSELMIndices;
  int     *meshSELMIndices;
  int     *meshSELMJ;

  double  *meshSELMX;

  double   meshFactor;

  int      index;
  int      numDir;

  int      meshSELMForceDensity_num_dim;
  double  *meshSELMForceDensity;
  double  *meshSELMForceDensity_ptr;

  int      flagSELM_Translational_forceDensity;
  int      flagSELM_Rotational_forceDensity;

  int      flagDegFreedomType = -1;
  double   X_cm[3];
  double   F_cm[3];

  double   Theta[3];
  double   T_ang[3];

  double   R_0;
  double  *weightFuncVals;
  double   weightFuncSum;

  operatorDataType_T_KERNEL_1 *opData;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmExtras;

  IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType
    *IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaExtras;

  IB_appl1_computeSmoothedForceFieldExtrasType
    *IB_appl1_computeSmoothedForceFieldExtras;

  if (num_dim > 3) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("  Only dimensions <= 3 are implemented \n");
    printf("  num_dim = %d \n", num_dim);
    packageError(1, 0);
  }

  switch (operatorType) {

  case OPERATOR_TYPE_T_KERNEL_1:

    opData = (operatorDataType_T_KERNEL_1 *) this->operatorData;

    /* -- Loop over the control points and apply the forces to the mesh
     *    via the approximate Dirac-delta-functions.
     */
    flagDegFreedomType = -1;
    for (k = 0; k < SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->numControlPts; k++) {

      /* -- Get control point information */
      for (d = 0; d < num_dim; d++) {
        X_cm[d] = SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->ptsX[(k * num_dim) + d];
        F_cm[d] = SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID->pt_Force[(k * num_dim) + d];
      }

      /* == Compute the force density. */

      /* -- Compute the force density associated with the surface
       *    of a sphere and the interpolations used. */

      /* -- determine the mesh indices */

      /* perform index arithmetic to find a mesh point X0 near X */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
        meshI0[d]      = floor(((X_cm[d] - meshBaseX0[d])/meshDeltaX) + 0.5);
        meshX0[d]      = meshBaseX0[d] + meshI0[d]*meshDeltaX;
      } /* end of d loop */

      /* now determine the indices for all mesh points with square
       * patch of size sphereR + kernalMeshwidths*meshDeltaX (depends on smoother used)
       */
      R_0                  = opData->weightTable->R_0;
      numDir               = 2*ceil(R_0) + 1;
      numMeshSELMIndices   = numDir*numDir*numDir;
      meshSELMIndices      = (int *)malloc(sizeof(int)*numMeshSELMIndices);
      meshSELMJ            = (int *)malloc(sizeof(int)*numMeshSELMIndices*num_dim);

      meshSELMX            = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);
      meshSELMForceDensity = (double *)malloc(sizeof(double)*numMeshSELMIndices*num_dim);

      /* compute base value for the mesh points */
      for (d = 0; d < num_dim; d++) {
        L              = meshDeltaX*numMeshPtsPerDir[d];
        meshBaseX0[d]  = meshCenterX0[d] - (L/2.0);
      }

      for (j1 = 0; j1 < numDir; j1++) {
        for (j2 = 0; j2 < numDir; j2++) {
          for (j3 = 0; j3 < numDir; j3++) {

            index    = j1 + j2*numDir + j3*numDir*numDir;

            meshJ[0] = meshI0[0] + j1 - ((numDir - 1)/2);
            meshJ[1] = meshI0[1] + j2 - ((numDir - 1)/2);
            meshJ[2] = meshI0[2] + j3 - ((numDir - 1)/2);

            /* compute the index modulo numMeshPtsPerDir */
            for (d = 0; d < num_dim; d++) {

              /* compute the mesh point locations before the adjustments are made to make it mod L */
              meshSELMX[index*num_dim + d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;

              /* make the indices mod N */
              meshSELMJ[index*num_dim + d] = (meshJ[d] % numMeshPtsPerDir[d]);
              if (meshSELMJ[index*num_dim + d] < 0) {
                meshSELMJ[index*num_dim + d] = meshSELMJ[index*num_dim + d] + numMeshPtsPerDir[d];
              }

            } /* end of d loop */

            /* compute the index on the mesh */
            meshI                    = meshSELMJ[index*num_dim + 0]
                                     + meshSELMJ[index*num_dim + 1]*numMeshPtsPerDir[0]
                                     + meshSELMJ[index*num_dim + 2]*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];
            meshSELMIndices[index]   = meshI;

          } /* end of j3 loop */
        } /* end of j2 loop */
      } /* end of j1 loop */

      /* -- compute the weights on the patch of mesh points */
      weightFuncVals = NULL;
      weightFromTable(num_dim, numMeshSELMIndices, meshSELMX,
                      X_cm, meshDeltaX,
                      opData->weightTable,
                      &weightFuncVals);

      weightFuncSum = 0;
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncSum += weightFuncVals[I];
      }

      /* re-normalize the weight values, since they should sum to one */
      for (I = 0; I < numMeshSELMIndices; I++) {
        weightFuncVals[I] = weightFuncVals[I]/weightFuncSum;
      }

      /* -- construct the force density by looping over the mesh points */
      meshFactor = 1.0/(meshDeltaX*meshDeltaX*meshDeltaX); /* weights include meshDeltaX term implicitly */
      for (I = 0; I < numMeshSELMIndices; I++) {

        meshI = meshSELMIndices[I];

        for (d = 0; d < num_dim; d++) {
          meshSELMForceDensity[I*num_dim + d] = F_cm[d]*weightFuncVals[I]*meshFactor;
        }

      } /* end of I loop */

      /* -- apply the force to the mesh */
      for (I = 0; I < numMeshSELMIndices; I++) {
        meshI = meshSELMIndices[I];
        for (d = 0; d < num_dim; d++) {
          f_m[d][meshI][0] += meshSELMForceDensity[I*num_dim + d];
          f_m[d][meshI][1] += 0;
        } /* end of d loop */
      } /* end of I loop */

      /* -- debug test for the force density */
      flagDebugSELM = 0;
      if (flagDebugSELM) {
        for (j1 = 0; j1 < numMeshPtsPerDir[0]; j1++) {
          for (j2 = 0; j2 < numMeshPtsPerDir[1]; j2++) {
            for (j3 = 0; j3 < numMeshPtsPerDir[2]; j3++) {

              I = j1 + j2*numMeshPtsPerDir[0]
                                           + j3*numMeshPtsPerDir[0]*numMeshPtsPerDir[1];

              meshJ[0] = j1;
              meshJ[1] = j2;
              meshJ[2] = j3;

              for (d = 0; d < num_dim; d++) {
                X[d] = meshBaseX0[d] + meshJ[d]*meshDeltaX;
              } /* end of d */

              for (d = 0; d < num_dim; d++) {
                if (d == 0) {
                  f_m[d][I][0] = X[d]; /* meshSELMForceDensity[I*num_dim + d]; */
                  f_m[d][I][1] = 0;
                } else {
                  f_m[d][I][0] += 0;
                  f_m[d][I][1] += 0;
                }
              } /* end of d */

            } /* end of j3 */
          } /* end of j2 */
        } /* end of j1 */
      }

      flagDebugSELM = 0;
      if (flagDebugSELM) {
        printf("WARNING: %s : %s \n", error_str_code, error_str_func);
        printf("Debug information being printed. \n");
        printf(" \n");
        fflush(stdout);

        printf("meshSELMX_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMX[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

        printf("meshSELMForceDensity_raw = [");
        for (I = 0; I < numMeshSELMIndices*num_dim; I++) {
          printf("%0.8g", meshSELMForceDensity[I]);
          if (I < numMeshSELMIndices*num_dim - 1) {
            printf(", ");
          }
        }
        printf("];\n");
        fflush(stdout);

      } /* end of flagDebugSELM */

      /* free dynamic memory */
      free(meshSELMX);
      free(meshSELMForceDensity);

      free(meshSELMIndices);
      free(meshSELMJ);

      free(weightFuncVals);

    } /* end of k loop */

    break;

  default:
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("This operator type is not yet supported in the implementation. \n");
    printf("  operatorTypeStr = %s \n", operatorTypeStr);
    packageError(1, 0);

  } /* end of switch */


}



void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::writeSimulationDataToDisk(const char *baseFilename, int timeIndex) {

  const char *error_str_func = "writeSimulationDataToDisk()";

  FILE *fid;
  char  filename[10000];

  /* open the file for writing the data */
  sprintf(filename, "%s_%.9d.SELM_CouplingOperator_%s", baseFilename, timeIndex, typeStr);
  fid = fopen(filename,"w");

  if (fid == NULL) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("Could not open file, error occured. \n");
    printf("  filename = %s \n", filename);
    packageError(1, 0);
  }

  fprintf(fid, "-- SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 : Simulation Data -- \n");
  fprintf(fid, "\n");

  /* close the file */
  fclose(fid);
}



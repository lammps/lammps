/* -------------------------------------------------------------------------
  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#include <cstdlib>
#include "stdio.h"
#include "string.h"
#include "SELM_Parser1.h"
#include "SELM_Eulerian.h"
#include "SELM_Eulerian_Types.h"
#include "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include <malloc.h>

using namespace LAMMPS_NS;

/* ============ constants definition =========== */
const double SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::UNIT_pi  = 3.141592653589793;

const int    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::TYPE     = SELM_Eulerian_Types::TYPE_FLUID_SHEAR_UNIFORM1_FFTW3;
const char*  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::TYPE_STR = SELM_Eulerian_Types::TYPE_STR_FLUID_SHEAR_UNIFORM1_FFTW3;

/* ---------------------------------------------------------------------- */
void SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::init() {

  /* identify the type of the mesh */
  type = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::TYPE;
  strcpy(typeStr, SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::TYPE_STR);

}

SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3() : SELM_Eulerian()
{

  /* generic initialization */
  init();

}

SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3(int narg, char **arg) : SELM_Eulerian(narg, arg)
{

  const char *error_str_code = "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.cpp";
  const char *error_str_func = "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3()";

  /* generic initialization */
  init();

  /*
  printf("Creating SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3 : not yet implemented. \n");
  exit(1);
  */

}

void SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::parse_ParameterFile(const char *baseFilename) {

  const char *error_str_code = "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.cpp";
  const char *error_str_func = "parse_ParameterFile()";

  char filename[10000];

  int j, k, N = 0;

  int paramIndex;
  int areAllSetFlag;

  const int maxNumPatchRefinements = 100;

  SELM_Parser1 *parser;

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ParamsType *params;

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params
  = (SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ParamsType *)
  malloc(sizeof(SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ParamsType));

  params = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params; /* notation */

  sprintf(filename, "%s.SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3", baseFilename);

  parser = new SELM_Parser1();

  /* Setup the parsing data structure */
  SELM_Parser1::paramSpecificationType *paramSpecification =
      (SELM_Parser1::paramSpecificationType *) malloc(sizeof(SELM_Parser1::paramSpecificationType));

  paramSpecification->paramDescrList
  = (SELM_Parser1::paramDescrType *) malloc(sizeof(SELM_Parser1::paramDescrType) * SELM_Parser1::MAX_NUM_PARAMS);

  paramIndex = 0;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "num_dim");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->num_dim;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName,
         "numMeshPtsPerDir");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT_LIST;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar
  = &params->numMeshPtsPerDir;
  paramSpecification->paramDescrList[paramIndex].paramExtras = &N;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "meshDeltaX");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_DOUBLE;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->meshDeltaX;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName,
         "meshCenterX0");
  paramSpecification->paramDescrList[paramIndex].paramType
  = SELM_Parser1::PARAMTYPE_DOUBLE_LIST;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar
  = params->meshCenterX0;
  paramSpecification->paramDescrList[paramIndex].paramExtras = &N;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName,
         "flagUseFluidPressure");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar
  = &params->flagUseFluidPressure;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "flagWriteSimulationData");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->flagWriteSimulationData;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "saveSkipSimulationData");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->saveSkipSimulationData;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "flagWriteFluidVel_VTK");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->flagWriteFluidVel_VTK;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "flagWriteFluidForce_VTK");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->flagWriteFluidForce_VTK;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "flagWriteFluidPressure_VTK");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->flagWriteFluidPressure_VTK;
  paramSpecification->paramDescrList[paramIndex].paramExtras = NULL;
  paramIndex++;

  paramSpecification->numParams = paramIndex;

  /* Parse the parameters */
  parser->parseParameters(filename, paramSpecification);

  /* Check whethor all parameters were set */
  areAllSetFlag = parser->areAllParametersSet(paramSpecification);

  if (areAllSetFlag == 0) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("The following Parameters were not set: \n");
    /* Show parameters parsed */
    parser->printUnsetParameters(paramSpecification);
    packageError(1, 0);
  }

  /* Show parameters parsed */
  printf("Parameters for SELM_Eulerian_%s: \n",typeStr);
  parser->printParameters(paramSpecification);

  /* Clean up the parsing data structures */
  free(paramSpecification->paramDescrList);
  free(paramSpecification);

  /* delete the parser */

  delete parser;

}

void SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::setup() {

  const char *error_str_code = "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.cpp";
  const char *error_str_func = "setup()";

  int I;
  int i, j, d;
  int *numDir;
  int N;
  int num_dim;

  /* Convert typeStr to type for later use.... */
  //driverSELM->SELM_Eulerian_List[0]->typeStr,

#ifdef USE_PACKAGE_FFTW3
  fftw_complex *in;
  fftw_complex *out;
#endif

#ifdef USE_PACKAGE_FFTW3

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras
    = (SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ExtrasType *)
       malloc(sizeof(SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ExtrasType));

  /* setup the mesh using the parameters data structure*/
  num_dim                                                      = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->num_dim;

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim               = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->num_dim;

  /* These parameter values should be set and overrridden by the integrator later, not from outside parsing
   * we set values to try to indicate this is uninitialized. */
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearRate             = 0;
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearDir              = -1;
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir           = -1;
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearDist             = 0;
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last        = 0;

  for (d = 0; d < num_dim; d++) {
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d] = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->numMeshPtsPerDir[d];
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0[d]     = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->meshCenterX0[d];
  }
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX            = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->meshDeltaX;

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure  = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->flagUseFluidPressure;

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->flagComputeStress     = 0;
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->controlPtsStress      = NULL;

  flagWriteSimulationData    = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->flagWriteSimulationData;
  saveSkipSimulationData     = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->saveSkipSimulationData;

  flagWriteFluidVel_VTK      = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->flagWriteFluidVel_VTK;
  flagWriteFluidForce_VTK    = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->flagWriteFluidForce_VTK;
  flagWriteFluidPressure_VTK = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params->flagWriteFluidPressure_VTK;




  /* == Allocate the FFTW3 arrays */
  N = 1;
  for (d = 0; d < SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
    N = N * SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d];
  }
  for (d = 0; d < SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m[d]
                                                        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d]
                                                        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m[d]
                                                            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_k[d]
                                                            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_m[d]
                                                                 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_k[d]
                                                                 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  } /* end of d loop */

  if (SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure) {
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m
    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_k
    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  } else {
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m = NULL;
    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_k = NULL;
  }


  /* == create the FFTW3 "plans" for the Fourier transforms to be applied. */
  /* NOTE: the FFTW3 codes use "row major" ordering as opposed to "column major"
   which is the convention we used.  To handle this without having to recode
   our indexing everywhere, we instead reverse our interpretation of the
   FFTW3 indices, so we allocate plans with sizes reverse what might be
   expected: N3, N2, N1 sized arrays instead of N1, N2, N3.  Once this is
   done the difference of convention in how entries are stored should not
   effect how the transforms actually perform.
   */

  /* force density transforms */
  numDir = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  for (d = 0; d < num_dim; d++) {
    in  = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m[d];
    out = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d];

    if (num_dim == 2) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_DFT_plan[d]
                                                              = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_IDFT_plan[d]
                                                               = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_DFT_plan[d]
                                                              = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                                                                                 FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_IDFT_plan[d]
                                                               = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                                                                                  FFTW_ESTIMATE);
    }

  } /* end of d loop */

  /* fluid drift velocity transforms */
  numDir = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  for (d = 0; d < num_dim; d++) {
    in = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m[d];
    out = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_k[d];

    if (num_dim == 2) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_DFT_plan[d]
                                                                  = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_IDFT_plan[d]
                                                                   = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_DFT_plan[d]
                                                                  = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                                                                                     FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_IDFT_plan[d]
                                                                   = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                                                                                      FFTW_ESTIMATE);
    }

  } /* end of d loop */

  /* fluid stochastic force density transforms */
  numDir = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  for (d = 0; d < num_dim; d++) {
    in = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_m[d];
    out = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_k[d];

    if (num_dim == 2) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_DFT_plan[d]
                                                                       = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_IDFT_plan[d]
                                                                        = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_DFT_plan[d]
                                                                       = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                                                                                          FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_IDFT_plan[d]
                                                                        = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                                                                                           FFTW_ESTIMATE);
    }

  } /* end of d loop */

  if (SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure) {

    /* pressure related transforms */
    numDir = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
    in = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m;
    out = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_k;

    if (num_dim == 2) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_DFT_plan
      = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_IDFT_plan
      = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_DFT_plan
      = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                         FFTW_ESTIMATE);
      SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_IDFT_plan
      = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                         FFTW_ESTIMATE);
    }

  } /* end flagUseFluidPressure */

#endif


  /* free the parameters data structure */
  free(SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params);

}

void SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::packageError(int code, void *extras) {
    exit(code);
}



void SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::userAppl_writeFFTW3VecFieldVTKFile(const char          *filename,
                                                                        int            num_dim,
                                                                        int           *numMeshPtsPerDir,
                                                                        double        *meshCenterX0,
                                                                        double        *meshLengths,
                                                                        int            numIndices,
                                                                        int           *indices,
                                                                        const char    *vec_name,
                                                                        fftw_complex **vec_array) {
  FILE  *fid;
  int    k;
  int    d;
  int    J[3];
  double X[3];
  double deltaX;
  int    index;

  fid = fopen(filename,"w");

  fprintf(fid,"# vtk DataFile Version 1.0\n");
  fprintf(fid,"FFTW3 Vector Field %s. \n",vec_name);
  fprintf(fid,"ASCII\n");
  fprintf(fid,"\n");

  /* legacy formats can be buggy, new XML formats preferred */
  if (numIndices >= 0) {
    fprintf(fid,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid,"POINTS %d float\n",numIndices);
    for (k = 0; k < numIndices; k++) {
      index = indices[k];
      J[2]    = index/(numMeshPtsPerDir[1]*numMeshPtsPerDir[0]);
      J[1]    = (index - J[2]*(numMeshPtsPerDir[1]*numMeshPtsPerDir[0]))
              / numMeshPtsPerDir[0];
      J[0]    = (index - J[2]*(numMeshPtsPerDir[1]*numMeshPtsPerDir[0])
          - J[1]*numMeshPtsPerDir[0]);

      for (d = 0; d < num_dim; d++) {
        deltaX = meshLengths[d]/numMeshPtsPerDir[d];
        X[d]   = meshCenterX0[d] - (meshLengths[d]/2.0) + J[d]*deltaX;
      }

      fprintf(fid,"%g %g %g", X[0], X[1], X[2]);
      fprintf(fid,"\n");
    }
    fprintf(fid,"\n");

    fprintf(fid,"CELLS %d %d\n",numIndices,2*numIndices);
    for (k = 0; k < numIndices; k++) {
      fprintf(fid,"1 %d \n",k);
    }
    fprintf(fid,"\n");

    fprintf(fid,"CELL_TYPES %d \n",numIndices);
    for (k = 0; k < numIndices; k++) {
      fprintf(fid,"1 \n");
    }
    fprintf(fid,"\n");

    fprintf(fid,"POINT_DATA %d \n",numIndices);
    fprintf(fid,"\n");

    fprintf(fid,"VECTORS %s float \n",vec_name);
    for (k = 0; k < numIndices; k++) {
      index = indices[k];
      for (d = 0; d < num_dim; d++) {
        fprintf(fid,"%g ", vec_array[d][index][0]);
      }
      fprintf(fid,"\n");
    }
    fprintf(fid,"\n");

  } /* numIndices >= 0 */

  /* legacy formats can be buggy, new XML formats preferred */
  if (numIndices == -1) {
    numIndices = numMeshPtsPerDir[0]*numMeshPtsPerDir[1]*numMeshPtsPerDir[2];
    fprintf(fid,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid,"POINTS %d float\n",numIndices);
    for (k = 0; k < numIndices; k++) {
      index = k;
      J[2]    = index/(numMeshPtsPerDir[1]*numMeshPtsPerDir[0]);
      J[1]    = (index - J[2]*(numMeshPtsPerDir[1]*numMeshPtsPerDir[0]))
              / numMeshPtsPerDir[0];
      J[0]    = (index - J[2]*(numMeshPtsPerDir[1]*numMeshPtsPerDir[0])
          - J[1]*numMeshPtsPerDir[0]);

      for (d = 0; d < num_dim; d++) {
        deltaX = meshLengths[d]/numMeshPtsPerDir[d];
        X[d]   = meshCenterX0[d] - (meshLengths[d]/2.0) + J[d]*deltaX;
      }

      fprintf(fid,"%g %g %g",X[0],X[1],X[2]);
      fprintf(fid,"\n");
    }
    fprintf(fid,"\n");

    fprintf(fid,"CELLS %d %d\n",numIndices,2*numIndices);
    for (k = 0; k < numIndices; k++) {
      fprintf(fid,"1 %d \n",k);
    }
    fprintf(fid,"\n");

    fprintf(fid,"CELL_TYPES %d \n",numIndices);
    for (k = 0; k < numIndices; k++) {
      fprintf(fid,"1 \n");
    }
    fprintf(fid,"\n");

    fprintf(fid,"POINT_DATA %d \n",numIndices);
    fprintf(fid,"\n");

    fprintf(fid,"VECTORS %s float \n",vec_name);
    for (k = 0; k < numIndices; k++) {
      index = k;
      for (d = 0; d < num_dim; d++) {
        fprintf(fid,"%g ", vec_array[d][index][0]);
      }
      fprintf(fid,"\n");
    }
    fprintf(fid,"\n");

  } /* numIndices == -1 */

  /* legacy formats can be buggy, new XML formats preferred */
  if (numIndices == -2) {
    fprintf(fid,"DATASET RECTILINEAR_GRID\n");
    fprintf(fid,"DIMENSIONS %d %d %d\n",
            numMeshPtsPerDir[0],
            numMeshPtsPerDir[1],
            numMeshPtsPerDir[2]);

    for (d = 0; d < num_dim; d++) {

      if (d == 0) {
        fprintf(fid,"X_COORDINATES\n");
      }
      if (d == 1) {
        fprintf(fid,"Y_COORDINATES\n");
      }
      if (d == 2) {
        fprintf(fid,"Z_COORDINATES\n");
      }

      deltaX = meshLengths[d]/numMeshPtsPerDir[d];
      for (k = 0; k < numMeshPtsPerDir[d]; k++) {
        X[d] = meshCenterX0[d] - (meshLengths[d]/2.0) + k*deltaX;
        fprintf(fid,"%g ",X[d]);
      }
      fprintf(fid,"\n");
    } /* end of d loop */
    fprintf(fid,"\n");


  } /* if numIndices == -1 */

  fclose(fid);

}



void SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3::writeSimulationDataToDisk(const char *baseFilename, int timeIndex) {

  const char *error_str_code = "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.cpp";
  const char *error_str_func = "writeSimulationDataToDisk()";

  FILE  *fid;
  int    d;
  double meshLengths[3];

  char  filename[10000];

  /* open the file for writing the data */
  sprintf(filename, "%s_%.9d.SELM_Eulerian_%s", baseFilename, timeIndex, typeStr);

  fid = fopen(filename,"w");

  if (fid == NULL) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("Could not open file, error occured. \n");
    printf("  filename = %s \n", filename);
    packageError(1, 0);
  }

  fprintf(fid, "-- SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3 : Simulation Data -- \n");
  fprintf(fid, "\n");

  fprintf(fid,"shearDir %d \n",
          SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearDir);

  fprintf(fid,"shearVelDir %d \n",
          SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir);

  fprintf(fid,"shearRate %lf \n",
          SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearRate);

  fprintf(fid,"shearDist %lf \n",
          SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last);

  /* close the file */
  fclose(fid);

  /* write additional data to disk */
  if (flagWriteFluidVel_VTK) {

    sprintf(filename,"%s_SELM_Eulerian_%s_FluidVel_%.9d.vtk", baseFilename, typeStr, timeIndex);

    for (d = 0; d < SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
      meshLengths[d] = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d]
                     * SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
    }

    userAppl_writeFFTW3VecFieldVTKFile(filename,
                                       SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                       (int *)    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                       (double *) SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                       (double *) meshLengths,
                                       -1, NULL,
                                       "fluid_velocity",
                                       SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m);

  } /* end of flagWriteFluidVel_VTK */

  if (flagWriteFluidForce_VTK) {

    sprintf(filename,"%s_SELM_Eulerian_%s_FluidForce_%.9d.vtk", baseFilename, typeStr, timeIndex);

    for (d = 0; d < SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
      meshLengths[d] = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d]
                     * SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
    }

    userAppl_writeFFTW3VecFieldVTKFile(filename,
                                       SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                       (int *)SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                       (double *)SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                       (double *)meshLengths,
                                       -1, NULL,
                                       "fluid_forceDensity",
                                       SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m);

  } /* end of flagWriteFluidForce_VTK */

  if (flagWriteFluidPressure_VTK) {

    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("Writing pressure to VTK format is not yet implemented \n");
    packageError(1,0);

    /*
    sprintf("%s_SELM_Eulerian_%s_FluidPressure.vtk", baseFilename, typeStr);

    for (d = 0; d < SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
      meshLengths[d] = SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d]
                     * SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
    }

    userAppl_writeFFTW3VecFieldVTKFile((char *)filename,
                                       SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                       (int *)    SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                       (double *) SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                       (double *) meshLengths,
                                       -1, NULL,
                                       "fluid_pressure",
                                       SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m);

    */
    
  } /* end of flagWriteFluidPressure_VTK */

}


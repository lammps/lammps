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

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_CONTROLPTS_BASIC1.h"

#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "lammps.h"
#include "domain.h"
#include "wrapper_selm.h"

#include <malloc.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
const int   SELM_Lagrangian_CONTROLPTS_BASIC1::TYPE     = SELM_Lagrangian_Types::TYPE_CONTROLPTS_BASIC1;
const char* SELM_Lagrangian_CONTROLPTS_BASIC1::TYPE_STR = SELM_Lagrangian_Types::TYPE_STR_CONTROLPTS_BASIC1;

SELM_Lagrangian_CONTROLPTS_BASIC1::SELM_Lagrangian_CONTROLPTS_BASIC1() : SELM_Lagrangian() {

  /* generic initialization */
  init();

}

SELM_Lagrangian_CONTROLPTS_BASIC1::SELM_Lagrangian_CONTROLPTS_BASIC1(int narg, char **arg) : SELM_Lagrangian(narg, arg) {

  const char *error_str_code = "SELM_Lagrangian_CONTROLPTS_BASIC1.cpp";
  const char *error_str_func = "SELM_Lagrangian_CONTROLPTS_BASIC1()";

  /* generic initialization */
  init();

}


SELM_Lagrangian_CONTROLPTS_BASIC1::SELM_Lagrangian_CONTROLPTS_BASIC1(class LAMMPS *lmps, class DriverSELM *fx) {

  /* generic initialization */
  init();

  setGlobalRefs(lammps, driverSELM);

}

void SELM_Lagrangian_CONTROLPTS_BASIC1::init() {

  type = SELM_Lagrangian_CONTROLPTS_BASIC1::TYPE;
  strcpy(typeStr, SELM_Lagrangian_CONTROLPTS_BASIC1::TYPE_STR);

  strcpy(nameStr, "No Name");

  num_dim             = 3;

  numControlPts       = 0;
  numControlPts_alloc = 0;

  pt_X                = NULL;
  pt_Vel              = NULL;

  pt_Energy           = 0;
  pt_Force            = NULL;

  pt_type             = NULL;
  pt_type_extras      = NULL;

  numEntriesOpGammaVel = 0;
  opGammaVel           = NULL;

  setGlobalRefs(NULL,NULL);

}


void SELM_Lagrangian_CONTROLPTS_BASIC1::parse_ParameterFile(const char *baseFilename) {

  const char *error_str_code = "SELM_Lagrangian_CONTROLPTS_BASIC1.cpp";
  const char *error_str_func = "parse_ParameterFile()";

  char filename[10000];

  int j, k, N = 0;

  int paramIndex;
  int areAllSetFlag;

  const int maxNumPatchRefinements = 100;

  SELM_Parser1 *parser;

  SELM_Lagrangian_CONTROLPTS_BASIC1_ParamsType *params;

  SELM_Lagrangian_CONTROLPTS_BASIC1_Params
  = (SELM_Lagrangian_CONTROLPTS_BASIC1_ParamsType *)
  malloc(sizeof(SELM_Lagrangian_CONTROLPTS_BASIC1_ParamsType));

  params = SELM_Lagrangian_CONTROLPTS_BASIC1_Params; /* notation */

  sprintf(filename, "%s.SELM_Lagrangian_CONTROLPTS_BASIC1", baseFilename);

  parser = new SELM_Parser1();

  /* Setup the parsing data structure */
  SELM_Parser1::paramSpecificationType *paramSpecification =
      (SELM_Parser1::paramSpecificationType *) malloc(sizeof(SELM_Parser1::paramSpecificationType));

  paramSpecification->paramDescrList
  = (SELM_Parser1::paramDescrType *) malloc(sizeof(SELM_Parser1::paramDescrType) * SELM_Parser1::MAX_NUM_PARAMS);

  paramIndex = 0;

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

  strcpy(paramSpecification->paramDescrList[paramIndex].paramName, "flagWriteControlPts_VTK");
  paramSpecification->paramDescrList[paramIndex].paramType = SELM_Parser1::PARAMTYPE_INT;
  paramSpecification->paramDescrList[paramIndex].paramSetFlag = 0;
  paramSpecification->paramDescrList[paramIndex].paramVar = &params->flagWriteControlPts_VTK;
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
  printf("Parameters for SELM_Lagrangian_%s: \n", typeStr);
  parser->printParameters(paramSpecification);

  /* Clean up the parsing data structures */
  free(paramSpecification->paramDescrList);
  free(paramSpecification);

  /* delete the parser */

  delete parser;

}


void SELM_Lagrangian_CONTROLPTS_BASIC1::setup() {

  const char *error_str_code = "SELM_Lagrangian_CONTROLPTS_BASIC1.cpp";
  const char *error_str_func = "setup()";

  /* setup the object from the parsed params */
  if (SELM_Lagrangian_CONTROLPTS_BASIC1_Params != NULL) {
    flagWriteSimulationData = SELM_Lagrangian_CONTROLPTS_BASIC1_Params->flagWriteSimulationData;
    saveSkipSimulationData  = SELM_Lagrangian_CONTROLPTS_BASIC1_Params->saveSkipSimulationData;
    flagWriteControlPts_VTK = SELM_Lagrangian_CONTROLPTS_BASIC1_Params->flagWriteControlPts_VTK;
  } else {
    printf("WARNING: %s : %s \n", error_str_code, error_str_func);
    printf("  SELM_Lagrangian_CONTROLPTS_BASIC1_Params == NULL \n");
    printf("  No initialization performed. \n");
  }

}



void SELM_Lagrangian_CONTROLPTS_BASIC1::setControlPtsDataFromLammpsData() {

  const char *error_str_code = "SELM_Lagrangian_CONTROLPTS_BASIC1.cpp";
  const char *error_str_func = "setControlPtsDataFromLammps()";

  int    j, k, d, I;
  int    N;

  Atom   *atom    = lammps->atom;

  int     nlocal  = atom->nlocal;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;

  /*
  int    igroup   = fix->igroup;
  int    groupbit = fix->groupbit;
  double *rmass   = atom->rmass;
  double *mass    = atom->mass;
  int    *type    = atom->type;
  int    *mask    = atom->mask;
  */

  /* set the values of the control points
   * using the LAMMPS data
   */
  num_dim         = lammps->domain->dimension;

  numControlPts   = lammps->atom->nlocal;

  /* allocate memory for the data, if needed */
  if (numControlPts_alloc < numControlPts) {

    N         = num_dim*numControlPts;

    if (pt_X != NULL) {
      free(pt_X);
    }
    pt_X      = (double *)malloc(sizeof(double)*N);

    if (pt_Vel != NULL) {
      free(pt_Vel);
    }
    pt_Vel    = (double *)malloc(sizeof(double)*N);

    /*
    if (opGammaVel != NULL) {
      free(opGammaVel);
    }

    numEntriesOpGammaVel = N;
    opGammaVel           = (double *)malloc(sizeof(double)*N);
    */

    pt_Energy = 0;

    if (pt_Force != NULL) {
      free(pt_Force);
    }
    pt_Force  = (double *)malloc(sizeof(double)*N);

    numControlPts_alloc = numControlPts;

  } else {
    /* nothing to do */
  }

  /* copy the data to the control points data structure */
  for (k = 0; k < numControlPts; k++) {
    for (d = 0; d < num_dim; d++) {

      I             = k*num_dim + d;

      pt_X[I]       = x[k][d];
      pt_Vel[I]     = v[k][d];
      pt_Force[I]   = f[k][d];

      /* data structures for SELM operators */
      /*
      opGammaVel[I] = 0;
      */

    } /* end of d loop */
  } /* end of k loop */

}


void SELM_Lagrangian_CONTROLPTS_BASIC1::setLammpsDataFromControlPtsData() {

  const char *error_str_code = "SELM_Lagrangian_CONTROLPTS_BASIC1.cpp";
  const char *error_str_func = "setLammpsDataFromControlPts()";

  int    j, k, d, I;
  int    N;

  Atom   *atom    = lammps->atom;

  int     nlocal  = atom->nlocal;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;

  /*
  int    igroup   = fix->igroup;
  int    groupbit = fix->groupbit;
  double *rmass   = atom->rmass;
  double *mass    = atom->mass;
  int    *type    = atom->type;
  int    *mask    = atom->mask;
  */

  /* set the values of the control points
   * using the LAMMPS data
   */
  num_dim         = lammps->domain->dimension;

  if (numControlPts != lammps->atom->nlocal) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("The control points data and LAMMPS are not synced. \n");
    printf("A different number of control points and local LAMMPS \n");
    printf("points was detected. \n");
    printf("numControlPts = %d \n", numControlPts);
    printf("lammps->atom->nlocal = %d \n", lammps->atom->nlocal);
    packageError(1, 0);
  }

  /* copy the control points data to the LAMMPS data structures */
  for (k = 0; k < numControlPts; k++) {
    for (d = 0; d < num_dim; d++) {

      I           = k*num_dim + d;

      x[k][d]     = pt_X[I];
      v[k][d]     = pt_Vel[I];
      f[k][d]     = pt_Force[I];

    } /* end of d loop */
  } /* end of k loop */

}


void SELM_Lagrangian_CONTROLPTS_BASIC1::packageError(int code, void *extras) {
    exit(code);
}



void SELM_Lagrangian_CONTROLPTS_BASIC1::userAppl_writePtsVTKFile(const char          *filename,
                              int            num_dim,
                              int            numPtsX,
                              const char    *ptsX_name,
                              double        *ptsX,
                              int            numScalarLists,
                              char         **scalarNames,
                              int           *numScalars,
                              double       **scalarLists,
                              int            numVecLists,
                              char         **vecNames,
                              int           *numVecs,
                              double       **vecLists) {

  FILE  *fid;
  int    k;
  int    d;
  int    I;

  int    numEntries;
  char   *scalar_name;
  double *scalar_array;
  char   *vec_name;
  double *vec_array;

  fid = fopen(filename,"w");

  fprintf(fid,"# vtk DataFile Version 1.0\n");
  fprintf(fid,"Pts Data %s (along with data at points). \n", ptsX_name);
  fprintf(fid,"ASCII\n");
  fprintf(fid,"\n");
  fprintf(fid,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fid,"POINTS %d float\n",numPtsX);
  for (k = 0; k < numPtsX; k++) {
    for (d = 0; d < num_dim; d++) {
      fprintf(fid,"%g ",ptsX[k*num_dim + d]);
      fprintf(fid,"\n");
    }
  }
  fprintf(fid,"\n");

  fprintf(fid,"CELLS %d %d\n",numPtsX,2*numPtsX);
  for (k = 0; k < numPtsX; k++) {
    fprintf(fid,"1 %d \n",k);
  }
  fprintf(fid,"\n");

  fprintf(fid,"CELL_TYPES %d\n",numPtsX);
  for (k = 0; k < numPtsX; k++) {
    fprintf(fid,"1 \n");
  }
  fprintf(fid,"\n");

  fprintf(fid,"POINT_DATA %d \n",numPtsX);
  fprintf(fid,"\n");

  /* loop over the scalar arrays */
  for (I = 0; I < numScalarLists; I++) {
    scalar_name   = scalarNames[I];
    scalar_array  = scalarLists[I];
    numEntries    = numScalars[I];

    fprintf(fid,"SCALARS %s float \n",scalar_name);
    fprintf(fid,"LOOKUP_TABLE default \n");
    for (k = 0; k < numEntries; k++) {
      fprintf(fid,"%g \n", scalar_array[k]);
    }
    fprintf(fid,"\n");
  } /* end of I loop */

  /* loop over the vector arrays */
  for (I = 0; I < numVecLists; I++) {
    vec_name   = vecNames[I];
    vec_array  = vecLists[I];
    numEntries = numVecs[I];

    fprintf(fid,"VECTORS %s float \n", vec_name);

    for (k = 0; k < numEntries; k++) {
      for (d = 0; d < num_dim; d++) {
        fprintf(fid,"%g ", vec_array[k*num_dim + d]);
      }
      fprintf(fid,"\n");
    }
    fprintf(fid,"\n");
  } /* end of I loop */

  fclose(fid);

}





void SELM_Lagrangian_CONTROLPTS_BASIC1::writeSimulationDataToDisk(const char *baseFilename, int timeIndex) {

  const char *error_str_code = "SELM_Lagrangian_CONTROLPTS_BASIC1.cpp";
  const char *error_str_func = "writeSimulationDataToDisk()";

  FILE *fid;
  int   i,j,k,d;

  int            numScalarLists;
  char         **scalarNames;
  int           *numScalars;
  double       **scalarLists;

  int            numVectorLists;
  char         **vectorNames;
  int           *numVectors;
  double       **vectorLists;

  char filename[10000];

  /* ===========  write control points data to disk =========== */

  /* open the file for writing the data */
  sprintf(filename,"%s_%.9d.SELM_Lagrangian_%s", baseFilename, timeIndex, typeStr);
  fid = fopen(filename,"w");

  if (fid == NULL) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("Could not open file, error occured. \n");
    printf("  filename = %s \n", filename);
    packageError(1, 0);
  }

  fprintf(fid, "-- SELM_Lagrangian_CONTROLPTS_BASIC1 : Simulation Data -- \n");
  fprintf(fid, "\n");
  fprintf(fid,"numControlPts %d \n",numControlPts);
  fprintf(fid,"num_dim %d \n", num_dim);

  fprintf(fid,"pt_X \n");
  for (j = 0; j < num_dim*numControlPts; j++) {
    fprintf(fid,"%.16g ", pt_X[j]);
  }
  fprintf(fid,"\n");

  fprintf(fid,"pt_Vel \n");
  for (j = 0; j < num_dim*numControlPts; j++) {
    fprintf(fid,"%.16g ", pt_Vel[j]);
  }
  fprintf(fid,"\n");

  fprintf(fid,"pt_Energy %.16g \n", pt_Energy);

  fprintf(fid,"pt_Force \n");
  for (j = 0; j < num_dim*numControlPts; j++) {
    fprintf(fid,"%.16g ", pt_Force[j]);
  }
  fprintf(fid,"\n");

  fprintf(fid,"pt_type \n");
  if (pt_type != NULL) {
    for (j = 0; j < numControlPts; j++) {
      fprintf(fid,"%d ", pt_type[j]);
    }
  } else {
    for (j = 0; j < numControlPts; j++) {
      fprintf(fid,"%d ", 0);
    }
  }
  fprintf(fid,"\n");

  fprintf(fid,"pt_type_extras \n");
  /* nothing implemented at this time */

  /* close the file */
  fclose(fid);


  /* ===========  write any additional data =========== */
  if (flagWriteControlPts_VTK) {

    /* all control points assumed to be lipids */
    sprintf(filename,"%s_SELM_Lagrangian_%s_%.9d.vtk",
            baseFilename, typeStr, timeIndex);

    numScalarLists    = 1;
    numScalars        = (int *)malloc(sizeof(int)*numScalarLists);
    scalarNames       = (char **)malloc(sizeof(char *)*numScalarLists);
    scalarLists       = (double **)malloc(sizeof(double *)*numScalarLists);

    numScalars[0]     = numControlPts;
    scalarNames[0]    = (char *)malloc(sizeof(char)*1000); // @optimzation: PJA: WARNING fixed size for names
    strcpy(scalarNames[0],"control_pts_index");
    scalarLists[0]    = (double *)malloc(sizeof(double)*numScalars[0]);
    for (k = 0; k < numScalars[0]; k++) {

      scalarLists[0][k] = k; /* no special grouping information conveyed, just index */

      /* could be used to group the control points */
      /*
      if ((k % 5) == 0) {
        scalarLists[0][k] = 0.0;
      }

      if ((k % 5) == 1) {
        scalarLists[0][k] = 1.0;
      }

      if ((k % 5) >= 2) {
        scalarLists[0][k] = 2.0;
      }
      */

    } /* end of k loop */

    numVectorLists    = 1;
    numVectors        = (int *)     malloc(sizeof(int)      * numVectorLists);
    vectorNames       = (char **)   malloc(sizeof(char *)   * numVectorLists);
    vectorLists       = (double **) malloc(sizeof(double *) * numVectorLists);

    vectorNames[0]    = (char *)    malloc(sizeof(char)*1000); // @optimzation: PJA: WARNING fixed size for names
    strcpy(vectorNames[0],"control_pts_velocity");
    num_dim           = num_dim;
    numVectors[0]     = numControlPts;
    vectorLists[0]    = (double *)malloc(sizeof(double)*numVectors[0]*num_dim);
    for (k = 0; k < numVectors[0]; k++) {
      for (d = 0; d < num_dim; d++) {
        vectorLists[0][k*num_dim + d] = pt_Vel[k*num_dim + d];
      }
    }

    userAppl_writePtsVTKFile(filename,
                             num_dim,
                             numControlPts,
                             "Control_Pts",
                             pt_X,
                             numScalarLists,
                             scalarNames,
                             numScalars,
                             scalarLists,
                             numVectorLists,
                             vectorNames,
                             numVectors,
                             vectorLists);

    /* free the scalar lists data */
    for (k = 0; k < numScalarLists; k++) {
      free(scalarNames[k]);
      free(scalarLists[k]);
    }
    free(numScalars);
    free(scalarNames);
    free(scalarLists);

    /* free the vector lists data */
    for (k = 0; k < numVectorLists; k++) {
      free(vectorNames[k]);
      free(vectorLists[k]);
    }
    free(numVectors);
    free(vectorNames);
    free(vectorLists);

  } /* end of flagWriteControlPts_VTK */

}


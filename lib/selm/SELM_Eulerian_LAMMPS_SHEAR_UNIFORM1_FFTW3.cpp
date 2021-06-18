/* -------------------------------------------------------------------------
  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------ 
*/

#include <cstdlib>
#include "stdio.h"
#include "string.h"
#include "wrapper_selm.h"
#include "SELM_Package.h"
#include "SELM_Parser1.h"
#include "SELM_Eulerian.h"
#include "SELM_Eulerian_Types.h"
#include "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include <malloc.h>

using namespace LAMMPS_NS;

/* ============ constants definition =========== */
const double SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::UNIT_pi  = 3.141592653589793;

const int    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::TYPE     = SELM_Eulerian_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3;
const char*  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::TYPE_STR = SELM_Eulerian_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3;

/* ---------------------------------------------------------------------- */
void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::init() {

  /* identify the type of the mesh */
  type = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::TYPE;
  strcpy(typeStr, SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::TYPE_STR);

  flagWriteSimulationData    = 0;
  saveSkipSimulationData     = 1;

  flagWriteFluidVel_VTK      = 0;
  flagWriteFluidForce_VTK    = 0;
  flagWriteFluidPressure_VTK = 0;

}

SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3() : SELM_Eulerian()
{

  /* generic initialization */
  init();

}

SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3(int narg, char **arg) : SELM_Eulerian(narg, arg)
{

  const char *error_str_code = "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3.cpp";
  const char *error_str_func = "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3()";

  /* generic initialization */
  init();

  /*
  printf("Creating SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 : not yet implemented. \n");
  exit(1);
  */

}


void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::setup() {

  const char *error_str_code = "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3.cpp";
  const char *error_str_func = "setup()";

  int I;
  int i, j, d;
  int *numDir;
  int N;
  int num_dim;

  /* @@@ Convert typeStr to type for later use.... */
  //driverSELM->SELM_Eulerian_List[0]->typeStr,

#ifdef USE_PACKAGE_FFTW3
  fftw_complex *in;
  fftw_complex *out;
#endif

#ifdef USE_PACKAGE_FFTW3

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras
    = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType *)
       malloc(sizeof(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType));
       
  // initialize to zero
  //memset(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras,
  //       0,sizeof(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ExtrasType));

  /* setup the mesh using the parameters data structure*/
  num_dim                                                                 = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->num_dim;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim               = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->num_dim;

  /* These parameter values may be overrridden by the integrator later, not from outside parsing */
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate             = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->shearRate;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir              = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->shearDir;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir           = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->shearVelDir;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist             = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->shearDist;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last        = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist;

  for (d = 0; d < num_dim; d++) {
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d] = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->numMeshPtsPerDir[d];
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0[d]     = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->meshCenterX0[d];
  }
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX            = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->meshDeltaX;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure  = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->flagUseFluidPressure;

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->flagComputeStress     = 0;
  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->controlPtsStress      = NULL;

  /*
  flagWriteSimulationData    = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->flagWriteSimulationData;
  saveSkipSimulationData     = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->saveSkipSimulationData;

  flagWriteFluidVel_VTK      = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->flagWriteFluidVel_VTK;
  flagWriteFluidForce_VTK    = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->flagWriteFluidForce_VTK;
  flagWriteFluidPressure_VTK = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params->flagWriteFluidPressure_VTK;
  */

  /* == Allocate the FFTW3 arrays */
  N = 1;
  for (d = 0; d < SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
    N = N * SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d];
  }
  for (d = 0; d < SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m[d]
                                                        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d]
                                                        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m[d]
                                                            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_k[d]
                                                            = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_m[d]
                                                                 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_k[d]
                                                                 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  } /* end of d loop */

  if (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure) {
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m
    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_k
    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  } else {
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m = NULL;
    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_k = NULL;
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
  numDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  for (d = 0; d < num_dim; d++) {
    in  = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m[d];
    out = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_k[d];

    if (num_dim == 2) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_DFT_plan[d]
                                                              = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_IDFT_plan[d]
                                                               = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_DFT_plan[d]
                                                              = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                                                                                 FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_IDFT_plan[d]
                                                               = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                                                                                  FFTW_ESTIMATE);
    }

  } /* end of d loop */

  /* fluid drift velocity transforms */
  numDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  for (d = 0; d < num_dim; d++) {
    in = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m[d];
    out = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_k[d];

    if (num_dim == 2) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_DFT_plan[d]
                                                                  = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_IDFT_plan[d]
                                                                   = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_DFT_plan[d]
                                                                  = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                                                                                     FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_IDFT_plan[d]
                                                                   = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                                                                                      FFTW_ESTIMATE);
    }

  } /* end of d loop */

  /* fluid stochastic force density transforms */
  numDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
  for (d = 0; d < num_dim; d++) {
    in = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_m[d];
    out = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_k[d];

    if (num_dim == 2) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_DFT_plan[d]
                                                                       = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_IDFT_plan[d]
                                                                        = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_DFT_plan[d]
                                                                       = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                                                                                          FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidStochForceDensity_IDFT_plan[d]
                                                                        = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                                                                                           FFTW_ESTIMATE);
    }

  } /* end of d loop */

  if (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->flagUseFluidPressure) {

    /* pressure related transforms */
    numDir = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir;
    in = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m;
    out = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_k;

    if (num_dim == 2) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_DFT_plan
      = fftw_plan_dft_2d(numDir[1], numDir[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_IDFT_plan
      = fftw_plan_dft_2d(numDir[1], numDir[0], out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (num_dim == 3) {
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_DFT_plan
      = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], in, out, FFTW_FORWARD,
                         FFTW_ESTIMATE);
      SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_IDFT_plan
      = fftw_plan_dft_3d(numDir[2], numDir[1], numDir[0], out, in, FFTW_BACKWARD,
                         FFTW_ESTIMATE);
    }

  } /* end flagUseFluidPressure */

#endif


  /* free the parameters data structure */
  free(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params);

}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::packageError(int code, void *extras) {
    exit(code);
}


void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::userAppl_writeFFTW3VecFieldVTKFile(const char   *filename,
                                                                                   int     num_dim,
                                                                                   int    *numMeshPtsPerDir,
                                                                                   double *meshCenterX0,
                                                                                   double *meshLengths,
                                                                                   int     numIndices,
                                                                                   int    *indices,
                                                                                   char   *vec_name,
                                                                                   fftw_complex **vec_array) {
  /* setup no shear conditions */
  int shearVelDir  = 0;
  int shearDir     = 2;
  double shearDist = 0;

  /* call the writing routine */
  userAppl_writeFFTW3VecFieldVTKFile(filename,num_dim,
                                     numMeshPtsPerDir,meshCenterX0,meshLengths,
                                     shearVelDir,shearDir,shearDist,
                                     numIndices,indices,vec_name,vec_array);
}


void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::userAppl_writeFFTW3VecFieldVTKFile(const char   *filename,
                                                                                   int     num_dim,
                                                                                   int    *numMeshPtsPerDir,
                                                                                   double *meshCenterX0,
                                                                                   double *meshLengths,
                                                                                   int     shearVelDir,
                                                                                   int     shearDir,
                                                                                   double  shearDist,
                                                                                   int     numIndices,
                                                                                   int    *indices,
                                                                                   char   *vec_name,
                                                                                   fftw_complex **vec_array) {
  FILE  *fid;
  int    k;
  int    d;
  int    J[3];
  double X[3];
  double deltaX;
  int    index;

  double L_shearDir;
  double shearAdj;

  fid = fopen(filename,"w");

  fprintf(fid,"# vtk DataFile Version 1.0\n");
  fprintf(fid,"FFTW3 Vector Field %s. \n",vec_name);
  fprintf(fid,"ASCII\n");
  fprintf(fid,"\n");

  /* not yet tested and debugged */
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

      /* Adjust to obtain the mesh site for the sheared coordinate system used. */
      L_shearDir      = numMeshPtsPerDir[shearDir] * deltaX;
      shearAdj        = (shearDist / L_shearDir) * (X[shearDir] - meshCenterX0[shearDir]);

      X[shearVelDir] += shearAdj;

      /* write the point location */
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

  /* not yet tested and debugged */
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

      /* Adjust to obtain the mesh site for the sheared coordinate system used. */
      L_shearDir      = numMeshPtsPerDir[shearDir] * deltaX;
      shearAdj        = (shearDist / L_shearDir) * (X[shearDir] - meshCenterX0[shearDir]);

      X[shearVelDir] += shearAdj;

      /* write the point location */
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

  /* not yet tested and debugged */
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



void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::writeSimulationDataToDisk(const char *basePath, int timeIndex) {  // @@@ ideally change this to base_dir 

  const char *error_str_code = "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3.cpp";
  const char *error_str_func = "writeSimulationDataToDisk()";

  FILE  *fid;
  int    I,d;
  double meshLengths[3];
  
  int num_vec_fields;
  char *vec_names[10]; // hard-coded (so keep less than 10 arrays, or adjust here)
  fftw_complex **vec_arrays[10];

  char  filename[10000];
  
  /* open the file for writing the data */
  sprintf(filename, "%s/SELM_Eulerian_%s_%.9d.xml", basePath, this->nameStr, timeIndex);

  fid = fopen(filename,"w");

  if (fid == NULL) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("Could not open file, error occured. \n");
    printf("  filename = %s \n", filename);
    packageError(1, 0);
  }
  
  fprintf(fid,"<?xml version=\"1.0\"?> \n");

  fprintf(fid, "<data>\n");

  fprintf(fid,"<shearDir value =\"%d\"/> \n",
          SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir);

  fprintf(fid,"<shearVelDir value=\"%d\"/> \n",
          SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir);

  fprintf(fid,"<shearRate value=\"%lf\"/> \n",
          SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearRate);

//  fprintf(fid,"shearDist %lf \n",
//          SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist_last);

  fprintf(fid,"<shearDist value=\"%lf\"/> \n",
          SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist);
          
  fprintf(fid, "</data>\n");

  /* close the file */
  fclose(fid);

  /* write additional data to disk */
  /*
    sprintf("%s.SELM_Eulerian_%s_EulerianRegionData_FluidVel", baseFilename, typeStr);
    writeFluidVelMeshRegionsData_FFTW3(filename,
                                       runDataExtras_SHEAR_POLYMER_FLUID1->fluidVelMeshRegionList,
                                       IB_appl1Mesh_UNIFORM1_FFTW3_Extras);
  */

  if (flagWriteFluidVel_VTK || flagWriteFluidForce_VTK || flagWriteFluidPressure_VTK) {

    //sprintf(filename,"%s_%s_SELM_Eulerian_%s_FluidVel_%.9d.vtk", baseFilename, this->nameStr, typeStr, timeIndex);
    sprintf(filename,"%s/SELM_Eulerian_%s_mesh_%.9d.vtu", basePath, this->nameStr, timeIndex);

    for (d = 0; d < SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
      meshLengths[d] = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d]
                     * SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
    }

    I = 0;
    if (flagWriteFluidVel_VTK) {
      const char *name = "fluid_velocity"; int nn;
      nn = strlen(name); vec_names[I] = (char *)malloc(sizeof(char)*(nn + 1)); strcpy(vec_names[I],name);
      vec_arrays[I] = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m;    
      I++;
    }
    
    if (flagWriteFluidForce_VTK) {
      const char *name = "fluid_force_density"; int nn;
      nn = strlen(name); vec_names[I] = (char *)malloc(sizeof(char)*(nn + 1)); strcpy(vec_names[I],name); 
      vec_arrays[I] = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m;
      I++;
    } 
    
    if (flagWriteFluidPressure_VTK) {
      //vec_name[I] = "fluid_force_density";      
      //vec_arrays[I] = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidForceDensity_m
      //I++;
      if (timeIndex < 5) {
        stringstream message;
        message << "Writing pressure to VTK format is not yet implemented." << endl;
        SELM_Package::packageError(error_str_code, error_str_func, message);
      }
    }
    
    num_vec_fields = I;
    
    userAppl_write3VecFieldVTUFile((char *)filename,
                                   SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                   (int *)    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                   (double *) SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                   (double *) meshLengths,
                                   SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir,
                                   SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir,
                                   SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist,
                                   num_vec_fields,
                                   vec_names,
                                   vec_arrays);

    // free the vector names array
    for (I = 0; I < num_vec_fields; I++) {
      free(vec_names[I]);
    }

    /*
    userAppl_writeFFTW3VecFieldVTKFile((char *)filename,
                                       SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                       (int *)    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                       (double *) SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                       (double *) meshLengths,
                                       SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearVelDir,
                                       SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDir,
                                       SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->shearDist,
                                       -1, NULL,
                                       "fluid_velocity",
                                       SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidDriftVel_m);
    */

  } /* end of flagWriteFluidVel_VTK / flagWriteFluidForce_VTK / flagWriteFluidPressure_VTK */

  //if (flagWriteFluidPressure_VTK) {

    /*
    sprintf("%s_%s_SELM_Eulerian_%s_FluidPressure.vtk", baseFilename, this->nameStr, typeStr);

    for (d = 0; d < SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim; d++) {
      meshLengths[d] = SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d]
                     * SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshDeltaX;
    }

    userAppl_writeFFTW3VecFieldVTKFile((char *)filename,
                                       SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->num_dim,
                                       (int *)    SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                       (double *) SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                       (double *) meshLengths,
                                       shearVelDir,
                                       shearDir,
                                       shearDist,
                                       -1, NULL,
                                       "fluid_pressure",
                                       SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Extras->fluidPressure_m);

    */
  //} /* end of flagWriteFluidPressure_VTK */


  /*
  if (runDataExtras_SHEAR_POLYMER_FLUID1->flagWriteFluidForceDensityMeshRegions) {
    sprintf(filename,"%s_%.6d.meshRegionsData_fluidForceDensity",
            aibControls->baseFilename,
            aibControls->timeIndex);

    writeFluidForceDensityMeshRegionsData_FFTW3(filename,
                                                runDataExtras_SHEAR_POLYMER_FLUID1->fluidForceDensityMeshRegionList,
                                                IB_appl1Mesh_UNIFORM1_FFTW3_Extras);

    flagWriteVTK = 1;
    if (flagWriteVTK) {
      sprintf(filename,"%s_%.6d_fluidForceDensity.vtk",
              aibControls->baseFilename,
              aibControls->timeIndex);

      for (d = 0; d < num_dim; d++) {
        meshLengths[d] = IB_appl1Mesh_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir[d]
                       * IB_appl1Mesh_UNIFORM1_FFTW3_Extras->meshDeltaX;
      }
      userAppl_writeFFTW3VecFieldVTKFile(filename,
                                         IB_appl1Mesh_UNIFORM1_FFTW3_Extras->num_dim,
                                         IB_appl1Mesh_UNIFORM1_FFTW3_Extras->numMeshPtsPerDir,
                                         IB_appl1Mesh_UNIFORM1_FFTW3_Extras->meshCenterX0,
                                         meshLengths,
                                         -1, NULL,
                                         "fluid_forceDensity",
                                         IB_appl1Mesh_UNIFORM1_FFTW3_Extras->fluidForceDensity_m);
    }

  }

  if (runDataExtras_SHEAR_POLYMER_FLUID1->flagWriteFluidPressureMeshRegions) {
    sprintf(filename,"%s_%.6d.meshRegionsData_fluidPressure",
            aibControls->baseFilename,
            aibControls->timeIndex);

    writeFluidPressureMeshRegionsData_FFTW3(filename,
                                            runDataExtras_SHEAR_POLYMER_FLUID1->fluidPressureMeshRegionList,
                                            IB_appl1Mesh_UNIFORM1_FFTW3_Extras);
  }

  */

}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::userAppl_write3VecFieldVTUFile(const char   *filename,
                                                                               int     num_dim,
                                                                               int    *numMeshPtsPerDir,
                                                                               double *meshCenterX0,
                                                                               double *meshLengths,
                                                                               int     shearVelDir,
                                                                               int     shearDir,
                                                                               double  shearDist,
                                                                               int     num_vec_fields,
                                                                               char   **vec_names,
                                                                               fftw_complex ***vec_arrays) {
  FILE  *fid;
  int    j,k,d,index;
  int    j1,j2,j3;
  int    n1,n2,n3;
  int    J[3],I[8];
  double X[3];
  double deltaX,L_shearDir,shearAdj;
  int    numPoints,numCells;
  int    num_cell_shape;
  
  numPoints = (numMeshPtsPerDir[0] + 1)*(numMeshPtsPerDir[1] + 1)*(numMeshPtsPerDir[2] + 1);
  numCells = numMeshPtsPerDir[0]*numMeshPtsPerDir[1]*numMeshPtsPerDir[2];

  /* write VTK XML file */
  fid = fopen(filename,"w");

  fprintf(fid,"<?xml version=\"1.0\"?> \n");
  fprintf(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"> \n");
  fprintf(fid,"<UnstructuredGrid> \n");
  fprintf(fid,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\"> \n",numPoints,numCells);
  fprintf(fid," \n");
  
  fprintf(fid,"<PointData Scalars=\"point_id\" Vectors=\"\"> \n");
  fprintf(fid," \n");
  fprintf(fid,"<DataArray type=\"Float32\" Name=\"point_id\" NumberOfComponents=\"1\" format=\"ascii\">  \n");
  for (k = 0; k < numPoints; k++) {
    fprintf(fid,"%d ",k);
    if ((k % 1000) == 0) { // put line break every 1000 points
      fprintf(fid," \n");
    }
  }
  fprintf(fid,"</DataArray> \n");
  fprintf(fid," \n");
  fprintf(fid,"</PointData> \n");
  fprintf(fid," \n");
  
  fprintf(fid,"<CellData> \n");
  fprintf(fid," \n");  
  fprintf(fid,"<DataArray type=\"Float32\" Name=\"cell_id\" NumberOfComponents=\"1\" format=\"ascii\"> \n");
  for (k = 0; k < numCells; k++) {
    fprintf(fid,"%d ",k);
    if ((k % 1000) == 0) { // put line break every 1000 points
      fprintf(fid," \n");
    }    
  }
  fprintf(fid," \n");
  fprintf(fid,"</DataArray> \n");
  fprintf(fid," \n");

  /* write the vector fields */
  for (j = 0; j < num_vec_fields; j++) {
    fprintf(fid,"<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"> \n",vec_names[j]);
    for (k = 0; k < numCells; k++) {
      for (d = 0; d < num_dim; d++) {
        fprintf(fid,"%g ",vec_arrays[j][d][k][0]); /* only save the real part */
      }
      fprintf(fid,"\n");
    }
    fprintf(fid,"\n");
    fprintf(fid,"</DataArray> \n");
    fprintf(fid,"\n");    
  }

  fprintf(fid,"<DataArray type=\"Float32\" Name=\"cell_center_pt\" NumberOfComponents=\"3\" format=\"ascii\"> \n");
  for (j3 = 0; j3 < numMeshPtsPerDir[2]; j3++) {
    for (j2 = 0; j2 < numMeshPtsPerDir[1]; j2++) {
      for (j1 = 0; j1 < numMeshPtsPerDir[0]; j1++) {
        J[0] = j1; J[1] = j2; J[2] = j3;
        
        for (d = 0; d < num_dim; d++) {
          deltaX = meshLengths[d]/numMeshPtsPerDir[d];
          X[d]   = meshCenterX0[d] - (meshLengths[d]/2.0) + J[d]*deltaX + 0.5*deltaX; // WARNING: Should double-check check vis. agrees with numerical conventions
        }

        /* Adjust to obtain the mesh site for the sheared coordinate system used. */
        L_shearDir      = numMeshPtsPerDir[shearDir] * deltaX;
        shearAdj        = (shearDist / L_shearDir) * (X[shearDir] - meshCenterX0[shearDir]);

        X[shearVelDir] += shearAdj;

        /* write the point location */
        fprintf(fid,"%g %g %g \n",X[0],X[1],X[2]);
      }
    }
  }
  fprintf(fid," \n");
  fprintf(fid,"</DataArray> \n");
  fprintf(fid," \n");
  fprintf(fid,"</CellData> \n");
  fprintf(fid," \n");
  
  fprintf(fid,"<Points> \n");
  fprintf(fid,"<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"-1000\" RangeMax=\"1000\"> \n");
  for (j3 = 0; j3 < numMeshPtsPerDir[2] + 1; j3++) {
    for (j2 = 0; j2 < numMeshPtsPerDir[1] + 1; j2++) {
      for (j1 = 0; j1 < numMeshPtsPerDir[0] +1; j1++) {
        J[0] = j1; J[1] = j2; J[2] = j3;
        //index = indices[k];
        //J[2]    = index/(numMeshPtsPerDir[1]*numMeshPtsPerDir[0]);
        //J[1]    = (index - J[2]*(numMeshPtsPerDir[1]*numMeshPtsPerDir[0]))
        //        / numMeshPtsPerDir[0];
        //J[0]    = (index - J[2]*(numMeshPtsPerDir[1]*numMeshPtsPerDir[0]) - J[1]*numMeshPtsPerDir[0]);

        for (d = 0; d < num_dim; d++) {
          deltaX = meshLengths[d]/numMeshPtsPerDir[d];
          X[d]   = meshCenterX0[d] - (meshLengths[d]/2.0) + J[d]*deltaX;
        }

        /* Adjust to obtain the mesh site for the sheared coordinate system used. */
        L_shearDir      = numMeshPtsPerDir[shearDir] * deltaX;
        shearAdj        = (shearDist / L_shearDir) * (X[shearDir] - meshCenterX0[shearDir]);

        X[shearVelDir] += shearAdj;

        /* write the point location */
        fprintf(fid,"%g %g %g \n",X[0],X[1],X[2]);
      }
    }
  }
  fprintf(fid," \n");  
  fprintf(fid,"</DataArray> \n");
  fprintf(fid,"</Points> \n");
  fprintf(fid," \n");
  
  fprintf(fid,"<Cells> \n");
  fprintf(fid,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"1e+299\" RangeMax=\"-1e+299\"> \n");
  /* Hexagedrons create */
  /*
  I = np.zeros((x.shape[0],8),dtype=int);
  I[:,0] = i3*n1*n2 + i2*n1 + i1; # base point index
  I[:,1] = I[:,0] + 1;
  I[:,2] = I[:,0] + 1 + n1;
  I[:,3] = I[:,0] + n1;

  I[:,4] = I[:,0] + n1*n2; # base point index
  I[:,5] = I[:,0] + 1 + n1*n2;
  I[:,6] = I[:,0] + 1 + n1 + n1*n2;
  I[:,7] = I[:,0] + n1 + n1*n2;
  */
  num_cell_shape = 8;
  n1 = numMeshPtsPerDir[0] + 1;n2 = numMeshPtsPerDir[1] + 1;n3 = numMeshPtsPerDir[2] + 1;
  for (j3 = 0; j3 < numMeshPtsPerDir[2]; j3++) {
    for (j2 = 0; j2 < numMeshPtsPerDir[1]; j2++) {
      for (j1 = 0; j1 < numMeshPtsPerDir[0]; j1++) {
        J[0] = j1; J[1] = j2; J[2] = j3;

        I[0] = j3*n1*n2 + j2*n1 + j1;
        I[1] = I[0] + 1;
        I[2] = I[0] + 1 + n1;
        I[3] = I[0] + n1;

        I[4] = I[0] + n1*n2; 
        I[5] = I[0] + 1 + n1*n2;
        I[6] = I[0] + 1 + n1 + n1*n2;
        I[7] = I[0] + n1 + n1*n2;

        /* write the point location */
        for (d = 0; d < num_cell_shape; d++) {
          fprintf(fid,"%d ",I[d]);
        }

        fprintf(fid,"\n");
      }
    }
  }
  fprintf(fid,"\n");
  fprintf(fid,"</DataArray> \n");
  fprintf(fid,"\n");
  fprintf(fid,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"-1e+299\" RangeMax=\"1e+299\"> \n");
  for (k = 0; k < numCells; k++) {
    fprintf(fid,"%d ",(k + 1)*8);
    if ((k % 1000) == 0) { // put line break every 1000 points
      fprintf(fid," \n");  
    }    
  }
  fprintf(fid,"\n");  
  fprintf(fid,"</DataArray> \n");
  fprintf(fid,"\n");
  fprintf(fid,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"-1e+299\" RangeMax=\"1e+299\"> \n");
  for (k = 0; k < numCells; k++) {
    fprintf(fid,"12 ");
    if ((k % 1000) == 0) { // put line break every 1000 points
      fprintf(fid," \n");  
    }        
  }
  fprintf(fid,"\n");
  fprintf(fid,"</DataArray> \n");
  fprintf(fid,"\n");
  fprintf(fid,"</Cells> \n");
  
  fprintf(fid,"\n");
  fprintf(fid,"</Piece> \n");
  fprintf(fid,"\n");
  fprintf(fid,"</UnstructuredGrid> \n");
  fprintf(fid,"\n");
  fprintf(fid,"</VTKFile> \n");
        
  fclose(fid);
  
}



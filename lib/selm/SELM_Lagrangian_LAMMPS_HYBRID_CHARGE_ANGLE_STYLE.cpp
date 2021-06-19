/*
 * SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE.cpp
 *
 * http://atzberger.org/
 *
 */
#include <cstdlib>
#include "stdio.h"
#include "string.h"

#include "SELM_Parser1.h"

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE.h"

#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "lammps.h"
#include "domain.h"

#include <malloc.h>

namespace LAMMPS_NS {

const int   SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::TYPE     = SELM_Lagrangian_Types::TYPE_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE;
const char* SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::TYPE_STR = SELM_Lagrangian_Types::TYPE_STR_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE;

const char *SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::error_str_code = "SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE.cpp";

SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE() : SELM_Lagrangian() {

	  /* generic initialization */
	  init();

	}

SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE(int narg, char **arg) : SELM_Lagrangian(narg, arg) {

  const char *error_str_func = "SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE()";

  /* generic initialization */
  init();

}


SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE(class LAMMPS *lmps, class DriverSELM *fx) {

  /* generic initialization */
  init();

  setGlobalRefs(lammps, driverSELM);

}

void SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::init() {

  type = SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::TYPE;
  strcpy(typeStr, SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::TYPE_STR);

  strcpy(nameStr, "No Name");

  num_dim             = 3;

  numControlPts       = 0;
  numControlPts_alloc = 0;

  ptsX                = NULL;
  pt_Vel              = NULL;

  atomID              = NULL;
  moleculeID          = NULL;
  typeID              = NULL;
  atomMass            = NULL;
  atomCharge          = NULL;

  pt_Energy           = 0;
  pt_Force            = NULL;

  pt_type             = NULL;
  pt_type_extras      = NULL;

  numEntriesOpGammaVel = 0;
  opGammaVel           = NULL;

  setGlobalRefs(NULL,NULL);

}


void SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::setup() {

  const char *error_str_func = "setup()";
  printf("WARNING: %s : %s \n", error_str_code, error_str_func);
  printf("  set() is no longer implemented  \n");

}



void SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::setControlPtsDataFromLammpsData() {

  const char *error_str_func = "setControlPtsDataFromLammps()";

  int    j, k, d, I;
  int    N;

  Atom   *atom    = lammps->atom;

  int     nlocal  = atom->nlocal;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;

  int    *type    = atom->type; /* not each Lagrangian data structure corresponds to one LAMMPS type */

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

  int num_LAMMPS_atoms      = lammps->atom->nlocal;
  int countNumPtsOfSameType = 0;

  for (k = 0; k < num_LAMMPS_atoms; k++) {
    if (type[k] == this->typeID[0]) { /* assumed only one type for the entire Lagrangian DOF */
      countNumPtsOfSameType++;
    }
  }

  numControlPts = countNumPtsOfSameType; /* number of points correspond to the Lagrangian DOF */

  /* allocate memory for the data, if needed */
  if (numControlPts_alloc < numControlPts) {

    N         = num_dim*numControlPts;

    if (ptsX != NULL) {
      free(ptsX);
    }
    ptsX      = (double *)malloc(sizeof(double)*N);

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
  int curI = 0;
  for (k = 0; k < num_LAMMPS_atoms; k++) {

    if (type[k] == this->typeID[0]) { /* only set using atoms of the same type as the Lagrangian DOF */

      for (d = 0; d < num_dim; d++) {

        I           = curI * num_dim + d;

        ptsX[I]     = x[k][d];
        pt_Vel[I]   = v[k][d];
        pt_Force[I] = f[k][d];

        /* data structures for SELM operators */
        /*
         opGammaVel[I] = 0;
         */

      } /* end of d loop */

      curI++; /* update the current type index to points */

    } /* end of same type check */

  } /* end of k loop */

}


void SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::setLammpsDataFromControlPtsData() {

  const char *error_str_func = "setLammpsDataFromControlPts()";

  int    j, k, d, I;
  int    N;

  Atom   *atom    = lammps->atom;

  int     nlocal  = atom->nlocal;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;

  int    *type    = atom->type;

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
  num_dim              = lammps->domain->dimension;

  int num_LAMMPS_atoms = lammps->atom->nlocal;

  int num_LAMMPS_set = 0;

  /* copy the control points data to the LAMMPS data structures */
  int curI = 0;
  for (k = 0; k < num_LAMMPS_atoms; k++) {

    if (type[k] == this->typeID[0]) { /* assume only one type per Lagrangian DOF */

      for (d = 0; d < num_dim; d++) {

        I = curI * num_dim + d;

        x[k][d] = ptsX[I];
        v[k][d] = pt_Vel[I];
        f[k][d] = pt_Force[I];

      } /* end of d loop */

      curI++;
      num_LAMMPS_set++;

    }

  } /* end of k loop */


  if (numControlPts != num_LAMMPS_set) {
      stringstream message;
      message << "The control points data and LAMMPS are not synced." << endl;
      message << "A different number of control points and set LAMMPS" << endl;
      message << "points was detected." << endl;
      message << "numControlPts = " << numControlPts << endl;
      message << "num_LAMMPS_set = " << num_LAMMPS_set << endl;
      SELM_Package::packageError(error_str_code, error_str_func, message);
    }


}


void SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::packageError(int code, void *extras) {
  const char *error_str_func = "packageError()";

    stringstream message;
    message << "code =" << code << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
}



void SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::userAppl_writePtsVTKFile(const char          *filename,
                              int            num_dim,
                              int            numPtsX,
                              const char          *ptsX_name,
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





void SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE::writeSimulationDataToDisk(const char *baseFilename, int timeIndex) {

  const char *error_str_func = "writeSimulationDataToDisk()";

  FILE *fid;
  int   i,j,k,d;

  int I;

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
  sprintf(filename,"%s_%s_%.9d.SELM_Lagrangian_%s", baseFilename, this->nameStr, timeIndex, typeStr);
  fid = fopen(filename,"w");

  if (fid == NULL) {
    printf("ERROR: %s : %s \n", error_str_code, error_str_func);
    printf("Could not open file, error occured. \n");
    printf("  filename = %s \n", filename);
    packageError(1, 0);

  }

  fprintf(fid, "-- SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE : Simulation Data -- \n");
  fprintf(fid, "\n");
  fprintf(fid,"numControlPts %d \n",numControlPts);
  fprintf(fid,"num_dim %d \n", num_dim);

  fprintf(fid,"pt_X \n");
  for (j = 0; j < num_dim*numControlPts; j++) {
    fprintf(fid,"%.16g ", ptsX[j]);
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
  if (flagWriteVTK) {

    /* all control points assumed to be lipids */
    sprintf(filename,"%s_%s_SELM_Lagrangian_%s_%.9d.vtk",
            baseFilename,  this->nameStr, typeStr, timeIndex);

    I=0;
    numScalarLists    = 1;
    numScalars        = (int *)malloc(sizeof(int)*numScalarLists);
    scalarNames       = (char **)malloc(sizeof(char *)*numScalarLists);
    scalarLists       = (double **)malloc(sizeof(double *)*numScalarLists);

    numScalars[I]     = numControlPts;
    scalarNames[I]    = (char *)malloc(sizeof(char)*100);
    strcpy(scalarNames[I],"control_pts_index");
    scalarLists[I]    = (double *)malloc(sizeof(double)*numScalars[0]);
    for (k = 0; k < numScalars[I]; k++) {
      //group the control points by its charges
      scalarLists[I][k]=atomCharge[k];

    } /* end of k loop */


    /* setup the vector information associated with the control points */
    I                 = 0;
    numVectorLists    = 2;
    numVectors        = (int *)     malloc(sizeof(int)      * numVectorLists);
    vectorNames       = (char **)   malloc(sizeof(char *)   * numVectorLists);
    vectorLists       = (double **) malloc(sizeof(double *) * numVectorLists);

    /* velocity */
    vectorNames[I]    = (char *)    malloc(sizeof(char)*100);
    strcpy(vectorNames[I],"control_pts_velocity");
    //num_dim           = num_dim;
    numVectors[I]     = numControlPts;
    vectorLists[I]    = (double *)malloc(sizeof(double)*numVectors[I]*num_dim);
    for (k = 0; k < numVectors[I]; k++) {
      for (d = 0; d < num_dim; d++) {
        vectorLists[I][k*num_dim + d] = pt_Vel[k*num_dim + d];
      }
    }
    I++;

    /* force */
    vectorNames[I]    = (char *)    malloc(sizeof(char)*100);
    strcpy(vectorNames[I],"control_pts_force");
    //num_dim           = num_dim;
    numVectors[I]     = numControlPts;
    vectorLists[I]    = (double *)malloc(sizeof(double)*numVectors[I]*num_dim);
    for (k = 0; k < numVectors[I]; k++) {
      for (d = 0; d < num_dim; d++) {
        vectorLists[I][k*num_dim + d] = pt_Force[k*num_dim + d];
      }
    }
    I++;

    /* write the VTK file */
    userAppl_writePtsVTKFile(filename,
                             num_dim,
                             numControlPts,
                             "Control_Pts",
                             ptsX,
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

  } /* end of flagWriteVTK */

}


}

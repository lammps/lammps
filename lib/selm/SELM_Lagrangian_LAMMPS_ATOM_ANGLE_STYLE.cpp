/*-----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#include <cstdlib>
#include "stdio.h"
#include <string.h>

//#include "SELM_Parser1.h"
#include "Atz_XML_Package.h"
#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE.h"

#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "lammps.h"
#include "domain.h"

#include <malloc.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
int   SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::TYPE     = SELM_Lagrangian_Types::TYPE_LAMMPS_ATOM_ANGLE_STYLE;
const char* SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::TYPE_STR = SELM_Lagrangian_Types::TYPE_STR_LAMMPS_ATOM_ANGLE_STYLE;

const char *SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::error_str_code = "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE.cpp";
// could use instead the pre-processor constant to get file automatically!
// char *SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::error_str_code = __FILE__;
// can also give line numbers by using __LINE__  (int type)  [Compiler substitutes]

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPE_NULL                = 0;
const char*  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPESTR_NULL             = "NULL";

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPE_SELM                = 1;
const char*  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPESTR_SELM             = "SELM";

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPE_VTK                 = 2;
const char*  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPESTR_VTK              = "vtk";

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPE_VTK_LEGACY          = 3;
const char*  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPESTR_VTK_LEGACY       = "VTK_legacy";

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPE_ATOM_ID             = 4;
const char*  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPESTR_ATOM_ID          = "atomID";

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPE_VELOCITY            = 5;
const char*  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPESTR_VELOCITY         = "velocity";

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPE_FORCE               = 6;
const char*  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TYPESTR_FORCE            = "force";

const int    SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::OUTPUTFLAG_TOTAL_NUM                = 7;

SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE() : SELM_Lagrangian() {

  /* generic initialization */
  init();

}

SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE(int narg, char **arg) : SELM_Lagrangian(narg, arg) {

  const char *error_str_func = "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE()";

  /* generic initialization */
  init();

}


SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE(class LAMMPS *lmps, class DriverSELM *fx) {

  /* generic initialization */
  init();

  setGlobalRefs(lammps, driverSELM);

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::init() {

  type = SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::TYPE;
  strcpy(typeStr, SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::TYPE_STR);

  strcpy(nameStr, "No Name");

  num_dim             = 3;

  numControlPts       = 0;
  numControlPts_alloc = 0;

  ptsX                = NULL;
  pt_Vel              = NULL;

  atomID              = NULL;
  atomLammpsIndex     = NULL;
  moleculeID          = NULL;
  typeID              = NULL;
  atomMass            = NULL;

  pt_Energy           = 0;
  pt_Force            = NULL;

  pt_type             = NULL;
  pt_type_extras      = NULL;

  numEntriesOpGammaVel = 0;
  opGammaVel           = NULL;

  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE_Params = NULL;

  setGlobalRefs(NULL,NULL);

  // setup output flags
  for (int i = 0; i < OUTPUTFLAG_TOTAL_NUM; i++) {
    outputFlagsStr[i][0] = 0;
    outputFlags[i]       = 0;
  }
  strcpy(outputFlagsStr[OUTPUTFLAG_TYPE_NULL],             OUTPUTFLAG_TYPESTR_NULL);
  strcpy(outputFlagsStr[OUTPUTFLAG_TYPE_SELM],             OUTPUTFLAG_TYPESTR_SELM);
  strcpy(outputFlagsStr[OUTPUTFLAG_TYPE_VTK],              OUTPUTFLAG_TYPESTR_VTK);
  strcpy(outputFlagsStr[OUTPUTFLAG_TYPE_VTK_LEGACY],       OUTPUTFLAG_TYPESTR_VTK_LEGACY);
  strcpy(outputFlagsStr[OUTPUTFLAG_TYPE_ATOM_ID],          OUTPUTFLAG_TYPESTR_ATOM_ID);
  strcpy(outputFlagsStr[OUTPUTFLAG_TYPE_VELOCITY],         OUTPUTFLAG_TYPESTR_VELOCITY);
  strcpy(outputFlagsStr[OUTPUTFLAG_TYPE_FORCE],            OUTPUTFLAG_TYPESTR_FORCE);

}


void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::setup() {

  const char *error_str_func = "setup()";
  
  setControlPtsDataFromLammpsData(); // initialize from lammps to start

}


void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::setControlPtsDataFromLammpsData() {

  const char *error_str_func = "setControlPtsDataFromLammps()";

  int    j, k, d, I;
  int    N;

  Atom   *atom    = lammps->atom;

  int     nlocal  = atom->nlocal;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;
  double *mass    = atom->mass;
  int    *tag     = atom->tag;

  int    *type    = atom->type; /* not each Lagrangian data structure corresponds to one LAMMPS type */
  int    *molecule= atom->molecule; /* not each Lagrangian data structure corresponds to one LAMMPS type */

  /*
  int    igroup   = fix->igroup;
  int    groupbit = fix->groupbit;
  double *rmass   = atom->rmass;

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
    countNumPtsOfSameType++;
  }

  numControlPts = countNumPtsOfSameType; /* number of points correspond to the Lagrangian DOF */

  /* allocate memory for the data, if needed */
  /* @optimization */
  /* Looks like the atom->x[0] is the point to a whole continuous block
   * of memory for the atom->x array but that the atom->x[0][d] indexing
   * is the type defined for this block for convenience and to avoid needing
   * to do excessive index arithmetic.  This means we can likely just
   * map the pointers to ours instead of performing expensive copies as
   * currently implemented.   Need to check with the LAMMPS developers before
   * making this assumption.
   */
  if (numControlPts_alloc < numControlPts) {

    N         = num_dim*numControlPts;
    if (ptsX != NULL) {
      free(ptsX);
    }
    ptsX      = (double *)malloc(sizeof(double)*N);

    N         = num_dim*numControlPts;
    if (pt_Vel != NULL) {
      free(pt_Vel);
    }
    pt_Vel    = (double *)malloc(sizeof(double)*N);

    N         = numControlPts;
    if (atomMass != NULL) {
      free(atomMass);
    }
    atomMass = (double *)malloc(sizeof(double)*N);
    for (k = 0; k < N; k++) { // assume does not change
      atomMass[k] = mass[atom->type[k]]; // mass array based on atom type
    }

    N         = numControlPts;
    if (atomID != NULL) {
      free(atomID);
    }
    atomID    = (int *)malloc(sizeof(int)*N);

    N         = numControlPts;
    if (moleculeID != NULL) {
      free(moleculeID);
    }

    if (atom->molecule != NULL) {  // check that LAMMPS has molecule data 
      moleculeID = (int *)malloc(sizeof(int)*N);
      for (k = 0; k < N; k++) { // assume does not change
	      moleculeID[k] = molecule[k];
      }
    } else {
      // Issue warning that molecule type not used for the atoms      
      stringstream message;
      message<<"Molecule data for atoms is NULL in LAMMPS" << endl;
      message<<"This indicates that no molecule data is available for the atoms specified." <<endl;
      message<<"atom->molecule   = NULL "<<endl;
      message<<"The SELM codes set NULL array for moleculeID." <<endl;
      message<<"SELM: moleculeID = NULL "<<endl;
      SELM_Package::packageWarning(error_str_code, error_str_func, message);
    }

    N         = numControlPts;
    if (typeID != NULL) {
      free(typeID);
    }
    typeID = (int *)malloc(sizeof(int)*N);
    for (k = 0; k < num_LAMMPS_atoms; k++) { // assume does not change
      typeID[k] = type[k];
    }

    N = numControlPts;
    if (atomLammpsIndex != NULL) {
      free(atomLammpsIndex);
    }
    atomLammpsIndex = (int *)malloc(sizeof(int)*N);
    for (k = 0; k < N; k++) { // assume does not change
      atomLammpsIndex[k] = k;
    }

    /*
    if (opGammaVel != NULL) {
      free(opGammaVel);
    }

    numEntriesOpGammaVel = N;
    opGammaVel           = (double *)malloc(sizeof(double)*N);
    */

    pt_Energy = 0;

    N         = num_dim*numControlPts;
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

    // WARNING: Ignore type check set all assuming only one type of Lagrangian throughout (see above)
    //if (type[k] == this->typeID[0]) { /* only set using atoms of the same type as the Lagrangian DOF */

    atomID[curI]          = tag[k];   /* WARNING: need to check the semantics of this assignment (tags -> id) */
    for (d = 0; d < num_dim; d++) {

      I                     = curI * num_dim + d;

      ptsX[I]               = x[k][d];
      pt_Vel[I]             = v[k][d];
      pt_Force[I]           = f[k][d];

      /* data structures for SELM operators */
      /*
       opGammaVel[I] = 0;
       */

    } /* end of d loop */

    curI++; /* update the current type index to points */

  //} /* end of same type check */

  } /* end of k loop */


  int flagComputeTotalForce = 0;
  if (flagComputeTotalForce) {

    double maxAbsForce;
    maxAbsForce = 0.0;
    double avgPtForce[3];
    avgPtForce[0] = 0.0; avgPtForce[1] = 0.0; avgPtForce[2] = 0.0;
    double absF = 0.0; double F;
    for (int I = 0; I < numControlPts; I++) {
      for (int d = 0; d < num_dim; d++) {
        F = pt_Force[I * num_dim + d];
        avgPtForce[d] += F;
        absF += F*F;
      }
      absF = sqrt(absF);
      if (absF > maxAbsForce) {
        maxAbsForce = absF;
      }
    }
    for (int d = 0; d < num_dim; d++) {
      avgPtForce[d] = avgPtForce[d]/numControlPts;
    }

    stringstream message;
    message << "FROM: " << error_str_code << ":" << error_str_func << endl;
    message   << "  timestep = " << lammps->update->ntimestep << endl;
    message   << "  lagrangian->nameStr = " << this->nameStr << endl;
    message   << "  maxAbsForce = " << maxAbsForce << endl;
    for (int d = 0; d < num_dim; d++) {
      message << "  avgPtForce[" << d << "] = " << avgPtForce[d] << endl;
    }
    cout << message.str();

  } /* flagComputeTotalForce */


}


void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::setLammpsDataFromControlPtsData() {

  const char *error_str_func = "setLammpsDataFromControlPts()";

  int    j, k, d, I;
  int    N;

  Atom   *atom    = lammps->atom;

  int     nlocal  = atom->nlocal;

  double **x      = atom->x;
  double **v      = atom->v;
  double **f      = atom->f;

  //int    *type     = atom->type;
  //int    *molecule = atom->molecule;

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

    // WARNING: changed codes to ignore types by assuming codes use just one Lagrangian type for now.
    // see above modifications
    //if (type[k] == this->typeID[0]) { /* assume only one type per Lagrangian DOF */

      for (d = 0; d < num_dim; d++) {

        I       = curI * num_dim + d;

        // should ideally check atomLammpsIndex[curI] == k;

        x[k][d] = ptsX[I];
        v[k][d] = pt_Vel[I];
        f[k][d] = pt_Force[I];

        // tags[k] = atomID[I];

      } /* end of d loop */

      curI++;
      num_LAMMPS_set++;

    //} /* end type check

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


void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::packageError(int code, void *extras) {
    exit(code);
}



void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::userAppl_writePtsVTKFile(const char *filename,
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





void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::writeSimulationDataToDisk(const char *baseFilename, int timeIndex) {

  const char *error_str_func = "writeSimulationDataToDisk()";

  FILE *fid;
  int   i,j,k,d;

  int I, vecCount;

  int            numScalarLists;
  char         **scalarNames;
  int           *numScalars;
  double       **scalarLists;

  int            numVectorLists;
  char         **vectorNames;
  int           *numVectors;
  double       **vectorLists;

  char filename[10000];

  // make sure the controlPts data is up to date
  setControlPtsDataFromLammpsData();

  /* ===========  write control points data to disk in SELM data file=========== */
  if (outputFlags[OUTPUTFLAG_TYPE_SELM]) {
    writeSELM(baseFilename, timeIndex);
  } /* end of SELM */

  /* ===========  VTK Legacy file =========== */
  if (outputFlags[OUTPUTFLAG_TYPE_VTK_LEGACY]) {
    writeVTKLegacy(baseFilename, timeIndex);
  } // end of VTK_LEGACY


  /* ===========  VTK XML file (.vtp) =========== */
  if (outputFlags[OUTPUTFLAG_TYPE_VTK]) {
    // write the .vtp file for particle configuration
    writeVTK(baseFilename, timeIndex);
  } // end of VTK XML format (.vtp)

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::setSimulationOutputFlags(const char *outputFlagStr) {

  char parseStr[1000];
  char tokenStr[100][100]; // assumes not more than 100 tokens!

  if (strcmp(outputFlagStr,"all") == 0) {

    // set all output flags to one
    for (int i = 0; i < OUTPUTFLAG_TOTAL_NUM; i++) {
      outputFlags[i] = 1;
    }

  } else { // loop over tokens and selectively set output
    // loop over the tokens of the string and set the flags
    strcpy(parseStr,outputFlagStr); // avoids destroying original string // @optimization PJA: uses fixed size string
    char *p        = strtok(parseStr, " ");
    int  k         = 0;
    int  numTokens = 0;
    while (p) {
      //printf ("Token: %s\n", p);
      strcpy(tokenStr[k],p);
      p = strtok(NULL, " ");
      k++;
    }
    numTokens = k;

    // reset all flags to default values
    resetSimulationOutputFlags();

    for (k = 0; k < numTokens; k++) {
      for (int i = 0; i < OUTPUTFLAG_TOTAL_NUM; i++) {
        if (strcmp(tokenStr[k],outputFlagsStr[i]) == 0) {
          outputFlags[i] = 1;
        }
      }
    }

  } /* end else */

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::resetSimulationOutputFlags() {

  // set all to zero
  for (int i = 0; i < OUTPUTFLAG_TOTAL_NUM; i++) {
    outputFlags[i] = 0;
  }

  // set default values that are on
  outputFlags[OUTPUTFLAG_TYPE_VTK]      = 1;
  outputFlags[OUTPUTFLAG_TYPE_ATOM_ID]  = 1;
  outputFlags[OUTPUTFLAG_TYPE_FORCE]    = 1;
  outputFlags[OUTPUTFLAG_TYPE_VELOCITY] = 1;

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::writeVTKLegacy(const char *baseFilename, int timeIndex) {

  char filename[10000];

  /* all control points assumed to be lipids */
  sprintf(filename,"%s_%s_SELM_Lagrangian_%s_%.9d.vtk",
          baseFilename,  this->nameStr, typeStr, timeIndex);

  writeVTKLegacy(filename);

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::writeVTKLegacy(const char *filename) {

  const char *error_str_func = "writeVTKLegacy()";

  FILE *fid;
  int   i,j,k,d;

  int I, vecCount;

  int            numScalarLists;
  char         **scalarNames;
  int           *numScalars;
  double       **scalarLists;

  int            numVectorLists;
  char         **vectorNames;
  int           *numVectors;
  double       **vectorLists;

  numScalarLists    = 1;
  numScalars        = (int *)malloc(sizeof(int)*numScalarLists);
  scalarNames       = (char **)malloc(sizeof(char *)*numScalarLists);
  scalarLists       = (double **)malloc(sizeof(double *)*numScalarLists);

  numScalars[0]     = numControlPts;
  scalarNames[0]    = (char *)malloc(sizeof(char)*100);
  strcpy(scalarNames[0],"control_pts_index");
  scalarLists[0]    = (double *)malloc(sizeof(double)*numScalars[0]);
  for (k = 0; k < numScalars[0]; k++) {
    scalarLists[0][k] = k; /* no special grouping information conveyed, just index */
  } /* end of k loop */

  /* setup the vector information associated with the control points */
  I                 = 0;
  vecCount          = 0;
  numVectorLists    = 2;
  numVectors        = (int *)     malloc(sizeof(int)      * numVectorLists);
  vectorNames       = (char **)   malloc(sizeof(char *)   * numVectorLists);
  vectorLists       = (double **) malloc(sizeof(double *) * numVectorLists);

  /* velocity */
  if (outputFlags[OUTPUTFLAG_TYPE_VELOCITY]) {
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
  }

  /* force */
  if (outputFlags[OUTPUTFLAG_TYPE_FORCE]) {
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
  }

  vecCount = I;

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
                           vecCount, //numVectorLists,
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

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::writeSELM(const char *baseFilename, int timeIndex) {
  char filename[10000];

  sprintf(filename,"%s_%s_%.9d.SELM_Lagrangian_%s", baseFilename, this->nameStr, timeIndex, typeStr);
  writeSELM(filename);

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::writeSELM(const char *filename) {

  const char *error_str_func = "writeSELM()";

  FILE *fid;
  int   i,j,k,d;

  int I, vecCount;

  /* open the file for writing the data */
      fid = fopen(filename,"w");

      if (fid == NULL) {
        printf("ERROR: %s : %s \n", error_str_code, error_str_func);
        printf("Could not open file, error occured. \n");
        printf("  filename = %s \n", filename);
        packageError(1, 0);

      }

      fprintf(fid, "-- SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE : Simulation Data -- \n");
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
}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::writeVTK(const char *baseFilename, int timeIndex) {
  char filename[10000];

  sprintf(filename,"%s_%s_%.9d.vtp", baseFilename, this->nameStr, timeIndex);
  writeVTK(filename);

}

void SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE::writeVTK(const char *filename) {

  const char *error_str_func = "writeVTK(filename)";

  int a = 1;

  stringstream extrasStr;

  ofstream fid;
  fid.open(filename);

  if (!fid.is_open()) {  // if file is not open
    stringstream message;
    message << "Could not open file to write error occured." << endl;
    message << "  filename = " << filename << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  Atz_XML_Package::writeXMLHeader(fid);

  extrasStr.str("");
  extrasStr << "type=\""       << "PolyData" << "\" ";
  extrasStr << "version=\""    << "0.1" << "\" ";
  extrasStr << "byte_order=\"" << "LittleEndian" << "\"";
  Atz_XML_Package::writeTagStart(fid,
                                 "VTKFile",
                                 extrasStr.str().c_str());

  Atz_XML_Package::writeTagStart(fid,"PolyData");

  // Piece
  extrasStr.str("");
  extrasStr << "NumberOfPoints=\"" << numControlPts << "\" ";
  extrasStr << "NumberOfVerts=\"0\" ";
  extrasStr << "NumberOfLines=\"0\" ";
  extrasStr << "NumberOfStrips=\"0\" ";
  extrasStr << "NumberOfPolys=\"0\"";

  Atz_XML_Package::writeTagStart(fid,"Piece",
                                 extrasStr.str().c_str());

  // PointData
  if (outputFlags[OUTPUTFLAG_TYPE_VELOCITY]) {
    extrasStr.str("");
    extrasStr << "Vectors=\"" << "Velocity" << "\"";
    Atz_XML_Package::writeTagStart(fid,"PointData",
                                   extrasStr.str().c_str());

  } else if (outputFlags[OUTPUTFLAG_TYPE_FORCE]) {
    extrasStr.str("");
    extrasStr << "Vectors=\"" << "Force" << "\"";
    Atz_XML_Package::writeTagStart(fid,"PointData",
                                   extrasStr.str().c_str());
  } else {
    Atz_XML_Package::writeTagStart(fid,"PointData");
  }

  // Data Arrays
  if ((outputFlags[OUTPUTFLAG_TYPE_ATOM_ID]) && (numControlPts != 0)) {

    // DataArray
    extrasStr.str("");
    extrasStr << "type=\"" << "Int32" << "\" ";
    extrasStr << "Name=\"" << "atomID" << "\" ";
    extrasStr << "NumberOfComponents=\"" << 1 << "\" ";
    extrasStr << "format=\"" << "ascii" << "\"";
    Atz_XML_Package::writeTagStart(fid,"DataArray",
                                   extrasStr.str().c_str());

    for (int I = 0; I < numControlPts; I++) {
      fid << this->atomID[I] << " ";
    }
    fid << endl;

    Atz_XML_Package::writeTagEnd(fid,"DataArray");
  }

  // Data Arrays
  if ((outputFlags[OUTPUTFLAG_TYPE_VELOCITY]) && (numControlPts != 0)) {

    // DataArray
    extrasStr.str("");
    extrasStr << "type=\"" << "Float32" << "\" ";
    extrasStr << "Name=\"" << "Velocity" << "\" ";
    extrasStr << "NumberOfComponents=\"" << num_dim << "\" ";
    extrasStr << "format=\"" << "ascii" << "\"";
    Atz_XML_Package::writeTagStart(fid,"DataArray",
                                   extrasStr.str().c_str());

    for (int I = 0; I < num_dim*numControlPts; I++) {
      fid << this->pt_Vel[I] << " ";
    }
    fid << endl;

    Atz_XML_Package::writeTagEnd(fid,"DataArray");
  }

  if ((outputFlags[OUTPUTFLAG_TYPE_FORCE]) && (numControlPts != 0)) {

    // DataArray
    extrasStr.str("");
    extrasStr << "type=\"" << "Float32" << "\" ";
    extrasStr << "Name=\"" << "Force" << "\" ";
    extrasStr << "NumberOfComponents=\"" << num_dim << "\" ";
    extrasStr << "format=\"" << "ascii" << "\"";
    Atz_XML_Package::writeTagStart(fid,"DataArray",
                                   extrasStr.str().c_str());

    for (int I = 0; I < num_dim*numControlPts; I++) {
      fid << this->pt_Force[I] << " ";
    }
    fid << endl;

    Atz_XML_Package::writeTagEnd(fid,"DataArray");
  }

  Atz_XML_Package::writeTagEnd(fid,"PointData");

  Atz_XML_Package::writeTagStart(fid,"cellData");
  Atz_XML_Package::writeTagEnd(fid,"cellData");

  Atz_XML_Package::writeTagStart(fid,"Points");

  if (numControlPts != 0) {
    // DataArray
    extrasStr.str("");
    extrasStr << "type=\"" << "Float32" << "\" ";
    extrasStr << "NumberOfComponents=\"" << num_dim << "\" ";
    extrasStr << "format=\"" << "ascii" << "\"";
    Atz_XML_Package::writeTagStart(fid,"DataArray",
                                   extrasStr.str().c_str());

    for (int I = 0; I < num_dim*numControlPts; I++) {
      fid << this->ptsX[I] << " ";
    }
    fid << endl;

    Atz_XML_Package::writeTagEnd(fid,"DataArray");
  }

  Atz_XML_Package::writeTagEnd(fid,"Points");

  Atz_XML_Package::writeTagStart(fid,"Verts");
  Atz_XML_Package::writeTagEnd(fid,"Verts");

  Atz_XML_Package::writeTagStart(fid,"Lines");
  Atz_XML_Package::writeTagEnd(fid,"Lines");

  Atz_XML_Package::writeTagStart(fid,"Strips");
  Atz_XML_Package::writeTagEnd(fid,"Strips");

  Atz_XML_Package::writeTagStart(fid,"Polys");
  Atz_XML_Package::writeTagEnd(fid,"Polys");

  Atz_XML_Package::writeTagEnd(fid,"Piece");

  Atz_XML_Package::writeTagEnd(fid,"PolyData");

  Atz_XML_Package::writeTagEnd(fid,"VTKFile");

  fid.close();

}



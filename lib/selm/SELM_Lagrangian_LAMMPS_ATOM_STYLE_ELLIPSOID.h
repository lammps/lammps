/* ----------------------------------------------------------------------
 Lagrangian mechanics.

 Paul J. Atzberger
 http://atzberger.org/
------------------------------------------------------------------------- */

#ifndef SELM_LAGRANGIAN_LAMMPS_ATOM_STYLE_ELLIPSOID_H
#define SELM_LAGRANGIAN_LAMMPS_ATOM_STYLE_ELLIPSOID_H

#ifndef FFT_FFTW3 /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_Types.h"

namespace LAMMPS_NS {

class SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID : public SELM_Lagrangian {

 public:

  /* ======================== Constants ======================= */
  static const int   TYPE;
  static const char* TYPE_STR;

  static const char *error_str_code;

  /* ======================== Function prototypes ======================= */
  SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID();
  SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID(int, char **);
  SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID(class LAMMPS *lmps, class DriverSELM *fx);

  // ~SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID();

  void init();

  //void parse_ParameterFile(char *filename);

  void setup();

  void setControlPtsDataFromLammpsData();
  void setLammpsDataFromControlPtsData();

  void writeSimulationDataToDisk(const char *baseFilename, int timeIndex);

  void setSimulationOutputFlags(const char *outputFlagStr);
  void resetSimulationOutputFlags();

  void writeSELM(const char *baseFilename, int timeIndex);
  void writeSELM(const char *filename);

  void writeVTKLegacy(const char *baseFilename, int timeIndex);
  void writeVTKLegacy(const char *filename);

  void writeVTK(const char *baseFilename, int timeIndex);
  void writeVTK(const char *filename);

  void userAppl_writePtsVTKFile(const char    *filename,
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
                                double       **vecLists);

  void packageError(int errorCode, void *extras);



  /* ======================== Data structure type definitions ======================== */
  typedef struct SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID_ParamsType {

    int                         flagWriteSimulationData;
    int                         saveSkipSimulationData;

    //int                         flagWriteControlPts_VTK;

  } SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID_ParamsType;


  /* ======================== Variables ======================= */

  // output flags
  const static int       OUTPUTFLAG_TYPE_NULL;
  const static char     *OUTPUTFLAG_TYPESTR_NULL;

  const static int       OUTPUTFLAG_TYPE_SELM;
  const static char     *OUTPUTFLAG_TYPESTR_SELM;

  const static int       OUTPUTFLAG_TYPE_VTK;
  const static char     *OUTPUTFLAG_TYPESTR_VTK;

  const static int       OUTPUTFLAG_TYPE_VTK_LEGACY;
  const static char     *OUTPUTFLAG_TYPESTR_VTK_LEGACY;

  const static int       OUTPUTFLAG_TYPE_ATOM_ID;
  const static char     *OUTPUTFLAG_TYPESTR_ATOM_ID;

  const static int       OUTPUTFLAG_TYPE_VELOCITY;
  const static char     *OUTPUTFLAG_TYPESTR_VELOCITY;

  const static int       OUTPUTFLAG_TYPE_FORCE;
  const static char     *OUTPUTFLAG_TYPESTR_FORCE;

  const static int       OUTPUTFLAG_TOTAL_NUM;

  SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID_ParamsType *SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID_Params;

  int            num_dim;             /* Number of spatial dimension of mesh. */

  int            numControlPts;       /* Number of control points. */
  int            numControlPts_alloc; /* Amount of memory allocated for control points (assumes 3D). */

  double        *ptsX;                /* Location of the control points. */
  double        *pt_Vel;              /* Velocity of the control points. */

  int           *atomID;
  int           *atomLammpsIndex;
  int           *moleculeID;
  int           *typeID;
  double        *atomMass;

  double         pt_Energy;           /* Energy associated with the control points. */
  double        *pt_Force;            /* Force associated with the control points. */

  int           *pt_type;             /* Type of the control points. */
  void         **pt_type_extras;      /* Extra data associated with the type. */

  //int            flagWriteVTK;

  char           outputFlagsStr[100][100];  // names of output flags
  int            outputFlags[100];          // value of output flag

  /* store data related to SELM operators */
  int            numEntriesOpGammaVel;
  double        *opGammaVel;          /* Velocity data resulting from the Gamma operator. */

  };

  /* =============================== Variables =========================== */

}

#endif

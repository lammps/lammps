/* ----------------------------------------------------------------------
 Lagrangian mechanics.

 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_LAGRANGIAN_CONTROLPTS_BASIC1_H
#define SELM_LAGRANGIAN_CONTROLPTS_BASIC1_H

#ifndef FFT_FFTW /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_Types.h"

namespace LAMMPS_NS {

class SELM_Lagrangian_CONTROLPTS_BASIC1 : public SELM_Lagrangian {

 public:

  /* ======================== Constants ======================= */
  static const int   TYPE;
  static const char* TYPE_STR;

  /* ======================== Function prototypes ======================= */
  SELM_Lagrangian_CONTROLPTS_BASIC1();
  SELM_Lagrangian_CONTROLPTS_BASIC1(int, char **);
  SELM_Lagrangian_CONTROLPTS_BASIC1(class LAMMPS *lmps, class DriverSELM *fx);

  // ~SELM_Lagrangian_CONTROLPTS_BASIC1();

  void init();

  void parse_ParameterFile(const char *filename);

  void setup();

  void setControlPtsDataFromLammpsData();
  void setLammpsDataFromControlPtsData();

  void writeSimulationDataToDisk(const char *baseFilename, int timeIndex);

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
  typedef struct SELM_Lagrangian_CONTROLPTS_BASIC1_ParamsType {

    int                         flagWriteSimulationData;
    int                         saveSkipSimulationData;

    int                         flagWriteControlPts_VTK;

  } SELM_Lagrangian_CONTROLPTS_BASIC1_ParamsType;


  /* ======================== Variables ======================= */

  SELM_Lagrangian_CONTROLPTS_BASIC1_ParamsType *SELM_Lagrangian_CONTROLPTS_BASIC1_Params;

  int            num_dim;             /* Number of spatial dimension of mesh. */

  int            numControlPts;       /* Number of control points. */
  int            numControlPts_alloc; /* Amount of memory allocated for control points (assumes 3D). */

  double        *pt_X;                /* Location of the control points. */
  double        *pt_Vel;              /* Velocity of the control points. */

  double         pt_Energy;           /* Energy associated with the control points. */
  double        *pt_Force;            /* Force associated with the control points. */

  int           *pt_type;             /* Type of the control points. */
  void         **pt_type_extras;      /* Extra data associated with the type. */

  int            flagWriteControlPts_VTK;

  /* store data related to SELM operators */
  int            numEntriesOpGammaVel;
  double        *opGammaVel;          /* Velocity data resulting from the Gamma operator. */

};

}

#endif

/*

 SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE.h
 
 http://atzberger.org/
  
*/

#ifndef SELM_LAGRANGIAN_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE_H_
#define SELM_LAGRANGIAN_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE_H_

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_Types.h"

namespace LAMMPS_NS {

class SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE: public LAMMPS_NS::SELM_Lagrangian {
public:
    //constant:
    static const int   TYPE;
    static const char* TYPE_STR;

    static const char *error_str_code;

	SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE();
	SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE(int, char **);
	SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE(class LAMMPS *lmps, class DriverSELM *fx);
//	virtual ~SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE();

    void init();

    void parse_ParameterFile(const char *filename);

    void setup();

    void setControlPtsDataFromLammpsData();
    void setLammpsDataFromControlPtsData();

    void writeSimulationDataToDisk(const char *baseFilename, int timeIndex);

    void userAppl_writePtsVTKFile(const char *filename,
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

    int                         flagWriteSimulationData;
    int                         saveSkipSimulationData;

    /* ======================== Variables ======================= */

    int            num_dim;             /* Number of spatial dimension of mesh. */

    int            numControlPts;       /* Number of control points. */
    int            numControlPts_alloc; /* Amount of memory allocated for control points (assumes 3D). */

    double        *ptsX;                /* Location of the control points. */
    double        *pt_Vel;              /* Velocity of the control points. */

    int           *atomID;
    int           *moleculeID;
    int           *typeID;
    double        *atomMass;
    double        *atomCharge;

    double         pt_Energy;           /* Energy associated with the control points. */
    double        *pt_Force;            /* Force associated with the control points. */

    int           *pt_type;             /* Type of the control points. */
    void         **pt_type_extras;      /* Extra data associated with the type. */

    int            flagWriteVTK;

    int            flagMobile;          /* 0: fixed atom, 1: mobile atom

    /* store data related to SELM operators */
    int            numEntriesOpGammaVel;
    double        *opGammaVel;          /* Velocity data resulting from the Gamma operator. */

};

}

#endif /* SELM_LAGRANGIAN_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE_H_ */

/* ----------------------------------------------------------------------
 * Coupling operators.
 * 
 * Paul J. Atzberger
 * http://atzberger.org/
 *
------------------------------------------------------------------------- */

#ifndef SELM_COUPLINGOPERATOR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_H
#define SELM_COUPLINGOPERATOR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_H

#include "SELM_CouplingOperator.h"

#ifndef FFT_FFTW /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

#include "SELM_CouplingOperator.h"
#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_CONTROLPTS_BASIC1.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE.h"
#include "SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID.h"
#include "SELM_Eulerian.h"
#include "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.h"
#include "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3.h"

namespace LAMMPS_NS {

class SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 : public SELM_CouplingOperator {

public:

  /* ================= constants ================= */
  /* constants (values for non-integers defined in .cpp) */

  /* specific operators to consider */
  static const int    TYPE;
  static const char*  TYPE_STR;

  static const int    COUPLING_OP_TYPE_NULL        = 0;
  static const char*  COUPLING_OP_TYPE_STR_NULL;

  static const int    COUPLING_OP_TYPE_GAMMA       = 1;
  static const char*  COUPLING_OP_TYPE_STR_GAMMA;

  static const int    COUPLING_OP_TYPE_LAMBDA      = 2;
  static const char*  COUPLING_OP_TYPE_STR_LAMBDA;

  static const int    COUPLING_OP_TYPE_UPSILON     = 3;
  static const char*  COUPLING_OP_TYPE_STR_UPSILON;

  /* specific operator types to use */
  static const char*  OPERATOR_TYPE_STR_NULL;
  static const int    OPERATOR_TYPE_NULL = 0;

  static const char*  OPERATOR_TYPE_STR_T_KERNEL_1;
  static const int    OPERATOR_TYPE_T_KERNEL_1 = 1;

  static const char*  OPERATOR_TYPE_STR_T_FAXEN_1;
  static const int    OPERATOR_TYPE_T_FAXEN_1 = 2;

  static const char*  OPERATOR_TYPE_STR_TR_FAXEN_1;
  static const int    OPERATOR_TYPE_TR_FAXEN_1 = 3;

  /* misc values */
  static double       UNIT_pi;

  /* error information */
  static const char* error_str_code;

  /* ================= data structure type definitions ================= */
  typedef struct controlPts_SELM_weightTableType {
    char    name[10000];
    int     num_dim;
    double  R_0;
    int     vecN[3];
    double  dX[3];
    int     numTableVals;
    double *X_list;
    double *weight_X_list;
  } controlPts_SELM_weightTableType;

  typedef struct operatorDataType_T_KERNEL_1 {
    char                             weightTableFilename[10000];
    controlPts_SELM_weightTableType *weightTable;
  } operatorDataType_T_KERNEL_1;


  typedef struct IB_appl1_computeSmoothedForceFieldExtrasType {

    int    num_dim;
    double meshDeltaX;
    int    numMeshPtsPerDir[3];
    double meshCenterX0[3];

    controlPts_SELM_weightTableType *SELM_weightTable;

  } IB_appl1_computeSmoothedForceFieldExtrasType;


  typedef struct IB_appl1_computeSmoothedVelFieldExtrasType {

    int    num_dim;
    double meshDeltaX;
    int    numMeshPtsPerDir[3];
    double meshCenterX0[3];

    int                              operatorType;
    controlPts_SELM_weightTableType *SELM_weightTable;

#ifdef USE_PACKAGE_FFTW3
    fftw_complex **u_m;
#endif

  } IB_appl1_computeSmoothedVelFieldExtrasType;


  typedef struct IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_XcmType {
    IB_appl1_computeSmoothedVelFieldExtrasType *smoothedVelExtras;
  } IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_XcmType;

  typedef struct IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_ThetaType {
    IB_appl1_computeSmoothedVelFieldExtrasType *smoothedVelExtras;
    double                                      sphereR;
    double                                      X_cm[3];
  } IB_appl1_userFuncExtras_TR_PARTICLE_Vel_sphFunc_ThetaType;

  typedef struct IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType {
    IB_appl1_computeSmoothedForceFieldExtrasType *smoothedForceExtras;

    double                                        X_m[3];
    double                                        controlPtForce[3];

    double                                        deltaX;

  } IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_XcmType;

  typedef struct IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType {
    IB_appl1_computeSmoothedForceFieldExtrasType *smoothedForceExtras;

    double                                        X_m[3];
    double                                        controlPtTorque[3];

    double                                        deltaX;

    double                                        sphereR;
    double                                        X_cm[3];

  } IB_appl1_userFuncExtras_TR_PARTICLE_Force_sphFunc_ThetaType;


  typedef struct IB_appl1_userFuncExtras_WEIGHT_FUNC1_PARTICLE_Vel_sphFunc_XcmType {
    int a;
    /*IB_appl1_computeSmoothedVelFieldExtrasType *smoothedVelExtras; */
  } IB_appl1_userFuncExtras_WEIGHT_FUNC1_PARTICLE_Vel_sphFunc_XcmType;

  typedef struct IB_appl1_userFuncExtras_WEIGHT_FUNC1_PARTICLE_Force_sphFunc_XcmType {
    int a;
    /*
      IB_appl1_computeSmoothedForceFieldExtrasType *smoothedForceExtras;

      double                                        X_m[3];
      double                                        controlPtForce[3];

      double                                        deltaX;
     */

  } IB_appl1_userFuncExtras_WEIGHT_FUNC1_PARTICLE_Force_sphFunc_XcmType;

  /* ================= function prototypes ================= */
  SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1();
  SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1(int, char **);

  void parse_ParameterFile(const char *filename);

  int getOperatorTypeFromStr(const char *operatorTypeStr);

  void setup();

  void packageError(int errorCode, void *extras);

  void writeSimulationDataToDisk(const char *filename, int timeIndex);

  void initConstants();
  void init();

  void computeOperator(const char      *couplingOpTypeStr,
                       SELM_Lagrangian *SELM_LagrangianData,
                       SELM_Eulerian   *SELM_EulerianData);

  /* Gamma operators */
  void computeOperatorGamma(SELM_Lagrangian *SELM_LagrangianData,
                            SELM_Eulerian       *SELM_EulerianData);

  void computeOperatorGamma(SELM_Lagrangian_CONTROLPTS_BASIC1        *SELM_LagrangianData_BASIC,
                            SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3);

  void computeOperatorGamma(SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE  *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE,
                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3);

  void computeOperatorGamma(SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE,
                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3);
                            
  void computeOperatorGamma(SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID,
                            SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3);                            

  /* Lambda operators */
  void computeOperatorLambda(SELM_Lagrangian *SELM_LagrangianData,
                             SELM_Eulerian   *SELM_EulerianData);

  void computeOperatorLambda(SELM_Lagrangian_CONTROLPTS_BASIC1         *SELM_LagrangianData_BASIC,
                             SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_FLUID_SHEAR_UNIFORM1_FFTW3);

  void computeOperatorLambda(SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE   *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE,
                             SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3 *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3);
                            
  void computeOperatorLambda(SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE,
                             SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3);
                             
  void computeOperatorLambda(SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID,
                             SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3  *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3);                             


  /* Upsilon operators */
  void computeOperatorUpsilon(SELM_Lagrangian *SELM_LagrangianData);

  /* generic functions used */
  void IB_appl1_compute_SELM_WEIGHT_FUNC1(int                               num_dim,
                                          int                               numPts,
                                          double                           *X_list,
                                          double                            deltaX,
                                          SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::controlPts_SELM_weightTableType  *SELM_weightTable,
                                          int                              *numEval,
                                          double                          **eval_ptr);

  void readWeightTable(char const *filename,
                       controlPts_SELM_weightTableType **SELM_weightTable_ptr);

  void weightFromTable(int num_dim, int numPts, double *ptsX,
                       double *X_cm, double meshDeltaX,
                       controlPts_SELM_weightTableType *SELM_weightTable,
                       double **weightFuncVals_ptr);

  void IB_appl1_computeSmoothedVelField(int      num_dim, int numPts, double *X_list,
                                        void *userExtras, double **Vel_list);

  void IB_appl1_computeLebedevSphereAvg(int numQuadNodes, double sphereR, int num_dim,
                                        double *Xcm,
                                        void   (*userFunc)(int num_dim, int numQuadNodes,
                                                           double *nodeX, void *userExtras,
                                                           int *funcVal_num_dim,
                                                           double **funcVal),
                                        void *userExtras,
                                        int  *integralVal_num_dim, double **integralEstimate_ptr);



  void IB_appl1_getLebedevSphereQuad(int numQuadNodes, double **nodeX_ptr,
                                     double **weight_ptr);

  void IB_appl1_userFunc_TR_PARTICLE_Force_sphFunc_Theta(int      num_dim,
                                                         int      numPts,
                                                         double  *X_list,
                                                         void    *userExtras,
                                                         int     *funcVal_num_dim,
                                                         double **funcVal_ptr);

  void IB_appl1_userFunc_TR_PARTICLE_Force_sphFunc_Xcm(int      num_dim,
                                                       int      numPts,
                                                       double  *X_list,
                                                       void    *userExtras,
                                                       int     *funcVal_num_dim,
                                                       double **funcVal);

  void IB_appl1_userFunc_TR_PARTICLE_Vel_sphFunc_Theta(int      num_dim, int   numPts,
                                                       double  *X_list,  void *userExtras,
                                                       int     *funcVal_num_dim,
                                                       double **funcVal_ptr);

  void IB_appl1_userFunc_TR_PARTICLE_Vel_sphFunc_Xcm(int num_dim, int numPts,
                                                     double *X_list, void *userExtras,
                                                     int *funcVal_num_dim,
                                                     double **funcVal);

  void IB_appl1_applyControlPtsForceToMesh_FFTW3(int                       num_dim,
                                                 int                       numMeshPtsPerDir[3],
                                                 double                    meshDeltaX,
                                                 double                   *meshCenterX0,
                                                 fftw_complex             *f_m[3],
                                                 SELM_Lagrangian_CONTROLPTS_BASIC1   *SELM_LagrangianData_CONTROLPTS_BASIC1);

  void IB_appl1_applyControlPtsForceToMesh_FFTW3(int                       num_dim,
                                                   int                       numMeshPtsPerDir[3],
                                                   double                    meshDeltaX,
                                                   double                   *meshCenterX0,
                                                   fftw_complex             *f_m[3],
                                                   SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE   *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE);

  void IB_appl1_applyControlPtsForceToMesh_FFTW3(int                       num_dim,
                                                   int                       numMeshPtsPerDir[3],
                                                   double                    meshDeltaX,
                                                   double                   *meshCenterX0,
                                                   fftw_complex             *f_m[3],
                                                   SELM_Lagrangian_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE   *SELM_LagrangianData_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE);
                                                   
  void IB_appl1_applyControlPtsForceToMesh_FFTW3(int                       num_dim,
                                                   int                       numMeshPtsPerDir[3],
                                                   double                    meshDeltaX,
                                                   double                   *meshCenterX0,
                                                   fftw_complex             *f_m[3],
                                                   SELM_Lagrangian_LAMMPS_ATOM_STYLE_ELLIPSOID   *SELM_LagrangianData_LAMMPS_ATOM_STYLE_ELLIPSOID);                                                   
                                                   

  /* ================= variables ================= */
  int              operatorType;
  char             operatorTypeStr[10000];

  void            *operatorData;


  int              numCoupleList;   /* number of lagrangian and eulerian DOF coupled */
  SELM_Lagrangian **lagrangianList; /* lagrangian info */
  SELM_Eulerian   **eulerianList;   /* eulerian info */

  /* operator data */
  int     numEntriesOpGammaResults;
  double *opGammaResults;

  int     numEntriesOpLambdaResults;
  double *opLambdaResults;

  int     numEntriesOpUpsilonResults;
  double *opUpsilonResults;

};

}


#endif

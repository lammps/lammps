/*
--------------------------------------------------------------------------------
  SELM fluctuating hydrodynamics coupling with shear boundary conditions.
    
  Paul J. Atzberger
  http://atzberger.org/
  
  Please cite the following paper when using these methods
  
  "Incorporating shear into stochastic Eulerian-Lagrangian methods for 
  rheological studies of complex fluids and soft materials," 
  Paul J. Atzberger, Physica D: Nonlinear Phenomena, 265, (2013).
    
  @article{atzberger_selm_shear_2013,
    title={Incorporating shear into stochastic Eulerian-Lagrangian methods 
           for rheological studies of complex fluids and soft materials},
    author={Paul J. Atzberger},
    journal={Physica D: Nonlinear Phenomena},
    year={2013},
    pages={57-70},
    volume={265},  
    doi={http://doi.org/10.1016/j.physd.2013.09.002},
    issn={0167-2789},
  }
  
  For examples and additional information see http://mango-selm.org/

--------------------------------------------------------------------------------
*/

#ifndef SELM_INTEGRATOR_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_H
#define SELM_INTEGRATOR_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_H

#ifndef FFT_FFTW3 /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

#include "driver_selm.h"

#include "SELM_Integrator.h"
#include "SELM_Integrator_Types.h"

#include "SELM_Lagrangian.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE.h"

#include "SELM_Eulerian.h"
#include "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3.h"

#include "SELM_CouplingOperator.h"

#include "lammps.h"
#include "fix.h"
#include "domain.h"
#include "irregular.h"

//#include "fftw.h"
#ifdef USE_PACKAGE_FFTW3
  #include "fftw3.h"
#endif

// forward declaration
/*
class SELM_Eulerian;
class SELM_Lagrangian;
class SELM_CouplingOperator;
class SELM_Interaction;
class SELM_Integrator;
*/

namespace LAMMPS_NS {

class SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3 : public SELM_Integrator {

 public:

  /* ================= constants ================= */
  static const int   TYPE;
  static const char* TYPE_STR;

  /* constants (values for non-integers defined in .cpp) */

  static const int   SHEAR_MODE_TYPE_NULL = 0;
  static const char* SHEAR_MODE_TYPE_STR_NULL;


  /*!\brief Key which controls the mode used for the shearing flow.  This case corresponds
   * to a flow with a constant shear rate.
   */
  static const int   SHEAR_MODE_TYPE_RM_SHEAR1 = 1;
  static const char* SHEAR_MODE_TYPE_STR_RM_SHEAR1;

  /*!\brief Key which controls the mode used for the shearing flow.  This case corresponds
   * to a flow with an oscillating shear rate.
   */
  static const int   SHEAR_MODE_TYPE_RM_OSC1   = 2;
  static const char* SHEAR_MODE_TYPE_STR_RM_OSC1;

  /*!\brief Key which controls the mode used.  In this case no shear applied to the system.
   */
  static const int SELM_INTEGRATOR_TYPE_OVERDAMPED1_NO_SHEAR                = 0;

  /*!\brief Key which controls the mode used for the shearing flow.  This case corresponds
   * to a flow with a constant shear rate.
   */
  static const int SELM_INTEGRATOR_TYPE_OVERDAMPED1_RM_SHEAR1               = 1;

  /*!\brief Key which controls the mode used for the shearing flow.  This case corresponds
   * to a flow with an oscillating shear rate.
   */
  static const int SELM_TYPE_OVERDAMPED1_RM_OSC1                           = 2;

  static double    UNIT_pi;

  static const char *error_str_code;
  
  // from fix_deform.cpp ---
  int remapflag;                   // whether x,v are remapped across PBC
  int dimflag[6];                  // which dims are deformed
  
  //protected:  
  int triclinic,scaleflag,flipflag;
  int flip,flipxy,flipxz,flipyz;   // flags for flipping atoms 
  double *h_rate,*h_ratelo;
  int varflag;                     // 1 if VARIABLE option is used, 0 if not
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes
  class Irregular *irregular;      // for migrating atoms after box flips
  // ---
 
  struct Set {
    int style,substyle;
    double flo,fhi,ftilt;
    double dlo,dhi,dtilt;
    double scale,vel,rate;
    double amplitude,tperiod;
    double lo_initial,hi_initial;
    double lo_start,hi_start,lo_stop,hi_stop,lo_target,hi_target;
    double tilt_initial,tilt_start,tilt_stop,tilt_target,tilt_flip;
    double tilt_min,tilt_max;
    double vol_initial,vol_start;
    int fixed,dynamic1,dynamic2;
    char *hstr,*hratestr;
    int hvar,hratevar;
  };
  Set *set; // for book-keeping box size and target (from fix_deform.h)
  
  /* ================= function prototypes ================= */
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3();
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3(LAMMPS *lmps, DriverSELM *fx);
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3(int, char **);
  ~SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3();

  void init();
  void setGlobalRefs(LAMMPS *lmps, DriverSELM *fx);
  
  void setup_internals();
  
  // from fix_deform.cpp --
  //int setmask(); // handled by another routine we call from fix for this
  virtual void pre_exchange(); // related to hooks for the fix
  virtual void end_of_step(); // related to hooks for the fix
  virtual void init_from_fix(); // related to hooks for the fix
  // --

  void parse_ParameterFile(char *filename);
  void setup();

  void set_LAMMPS_mask(int *mask_ptr);

  void integrate_initialize();
  void integrate_initial();
  void integrate_final();

  void SELM_updateFluidAndStructures_initial();
  void SELM_updateFluidAndStructures_final();

  void SELM_updateParticlesTest_initial();
  void SELM_updateParticlesTest_final();

  void syncShearDomainWithLammpsDomain();

  void packageError(int errorCode, void *extras);

  void writeSimulationDataToDisk(char *baseFilename, int timeIndex);

  /* ================= data structure types ================= */
  typedef struct ShearData_RM_SHEAR1_Type {
    double shearRate;                   /*!< Constant rate of the mesh shearing motion */
    int    shearDir;
    int    shearVelDir;
    double shearDist;
    double shearDist_last;
  } ShearData_RM_SHEAR1_Type;

  typedef struct ShearData_RM_OSC1_Type {
    double    shearOmega;               /*!< Frequency for use with oscillating shear flow. */
    double    shearRateAmplitude;       /*!< Amplitude of the shear rate in the oscillating flow case. */

    double    shearRate;
    int       shearDir;
    int       shearVelDir;
    double    shearDist;
    double    shearDist_last;
  } ShearData_RM_OSC1_Type;

  typedef struct SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType {

#ifdef USE_PACKAGE_FFTW3
    int       maxTimeStepIndex;
    double    deltaT;
    double    mu; /* fluid viscosity */
    double    rho;
    double    KB;
    double    T;

    char      flagShearModeStr[1000];   /*!< Indicates the type of shear flow to consider (steady, oscillating, etc...). */
    int       flagShearMode;            /*!< Indicates the type of shear flow to consider (steady, oscillating, etc...). */
    void     *shearData;                /*!< Data associated with the shearing behavior */

    int       flagStochasticDriving;
    int       flagIncompressibleFluid;

    int       flagWriteSimulationData;
    int       saveSkipSimulationData;

#endif

  } SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType;

  typedef struct SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType {

  #ifdef USE_PACKAGE_FFTW3
    int       flagInitializedNumericalMethod;  /*!< Indicates if numerical scheme initialized. */

    int       maxTimeStepIndex;         /*!< Maximum number of time steps to iterate. */
    double    deltaT;                   /*!< Time step for the numerical method. */
    double    mu;
    double    rho;
    double    KB;
    double    T;

    char      flagShearModeStr[1000];   /*!< Indicates the type of shear flow to consider (steady, oscillating, etc...). */
    int       flagShearMode;            /*!< Indicates the type of shear flow to consider (steady, oscillating, etc...). */
    void     *shearData;                /*!< Data associated with the shearing behavior */

    int       flagStochasticDriving;           /*!< Indicates if stochastic forcing is to be used. */
    int       flagIncompressibleFluid;         /*!< Indicates if projection enforcing incompressibility used. */

    int       flagUpdateControlPts;            /*!< Indicates if the numerical method should determine new control
                                                      point locations or if this will be handled elsewhere in the codes. */
  #endif

  } SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType;

  /* ================= variables ================= */
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType *SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params;
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ExtrasType *SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Extras;

 private:

  /* ================= function prototypes ================= */
  void addPseudoForceTerm(int num_dim,
                          double meshDeltaX,
                          int *numMeshPtsPerDir,
                          double mu,
                          int shearDir,
                          int shearVelDir,
                          double shearRate,
                          fftw_complex **fluidForceDensity_m);

  void computeTimeAvgStochFluct(int num_dim,
                                double meshDeltaX,
                                int *numMeshPtsPerDir,
                                double mu,
                                double KB,
                                double T,
                                double deltaT,
                                double shearRate,
                                int shearDir,
                                int shearVelDir,
                                double shearDist,
                                /*gaussDevStateType *gaussDevState, */
                                class RanMars *random,
                                fftw_complex **fluidStochForceDensity_k);

  void projectField(int num_dim, double meshDeltaX,
                    int *numMeshPtsPerDir,
                    int shearDir, int shearVelDir,
                    double shearDist,
                    fftw_complex **field_u_k);

  void computePressure(int num_dim,
                       double meshDeltaX,
                       int *numMeshPtsPerDir,
                       int shearDir,
                       int shearVelDir,
                       double shearDist,
                       fftw_complex **forceDensity_k,
                       fftw_complex *pressure_k);

  void bri1_unitCellRectImageShearPeriodic(double *periodL, double *meshCenterX0,
                                           int shearDir, int shearVelDir, double shearDist,
                                           int num_dim, double *X_orig, double *X_unitCell);

  void IB_appl1_unitCellRectImageShearPeriodic(double *periodL, double *meshCenterX0,
                                               int shearDir, int shearVelDir,
                                               double shearDist, int num_dim,
                                               double *X_orig, double *X_unitCell);

  void computeControlPtsVel_SHEAR_FFTW3(SELM_Lagrangian       *SELM_LagrangianData_LAMMPS_ATOM_ANGLE_STYLE,
                                        SELM_Eulerian         *SELM_EulerianData_LAMMPS_SHEAR_UNIFORM1_FFTW3,
                                        SELM_CouplingOperator *SELM_CouplingOperatorsData);

  void IB_appl1_start_time_step_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3();


};

}



#endif

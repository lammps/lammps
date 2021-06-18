/* ----------------------------------------------------------------------
 Eulerian sheared mesh.

 Paul J. Atzberger
 http://atzberger.org/

------------------------------------------------------------------------- */

#ifndef SELM_EULERIAN_FLUID_SHEAR_UNIFORM1_FFTW3_H
#define SELM_EULERIAN_FLUID_SHEAR_UNIFORM1_FFTW3_H

#ifndef FFT_FFTW /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

/* == includes */
#include "SELM_Eulerian.h"

#ifdef USE_PACKAGE_FFTW3
  #include "fftw3.h"
#endif

/* == define class */
namespace LAMMPS_NS {

class SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3 : public SELM_Eulerian {

 public:

  /* ======================== Constants ======================== */
  static const double    UNIT_pi;
  static const int       TYPE;
  static const char*     TYPE_STR;

  /* ======================== Function prototypes ======================== */
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3();
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3(int, char **);

  void init();
  void parse_ParameterFile(const char *filename);
  void setup();

  void packageError(int errorCode, void *extras);

  void writeSimulationDataToDisk(const char *filename, int timeIndex);

#ifdef USE_PACKAGE_FFTW3

  void userAppl_writeFFTW3VecFieldVTKFile(const char    *filename,
                                          int            num_dim,
                                          int           *numMeshPtsPerDir,
                                          double        *meshCenterX0,
                                          double        *meshLengths,
                                          int            numIndices,
                                          int           *indices,
                                          const char    *vec_name,
                                          fftw_complex **vec_array);

#endif

  // ~SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3();

  /* ======================== Data structure type definitions ======================== */
  typedef struct SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ParamsType {

  #ifdef USE_PACKAGE_FFTW3
    int                         num_dim;
    int                         numMeshPtsPerDir[3];
    double                      meshDeltaX;
    double                      meshCenterX0[3];

    int                         flagUseFluidPressure;

    int                         flagWriteSimulationData;
    int                         saveSkipSimulationData;

    int                         flagWriteFluidVel_VTK;
    int                         flagWriteFluidForce_VTK;
    int                         flagWriteFluidPressure_VTK;

  #endif

  } SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ParamsType;

  typedef struct SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ExtrasType {

  #ifdef USE_PACKAGE_FFTW3
    int                         num_dim;
    int                         numMeshPtsPerDir[3];
    double                      meshDeltaX;
    double                      meshCenterX0[3];

    /* values should be set by the integrator */
    double                      shearRate;                /*!< Rate of shear in the specified direction. */
    int                         shearDir;                 /*!< The direction along which the velocity field shears. */
    int                         shearVelDir;              /*!< The direction of the velocity field. */
    double                      shearDist_last;           /*!< Shear distance before update. */
    double                      shearDist;                /*!< Distance the reference frame has been sheared. */

    fftw_complex               *fluidDriftVel_m[3];
    fftw_complex               *fluidDriftVel_k[3];
    fftw_plan                   fluidDriftVel_DFT_plan[3];
    fftw_plan                   fluidDriftVel_IDFT_plan[3];

    fftw_complex               *fluidForceDensity_m[3];
    fftw_complex               *fluidForceDensity_k[3];
    fftw_plan                   fluidForceDensity_DFT_plan[3];
    fftw_plan                   fluidForceDensity_IDFT_plan[3];

    fftw_complex               *fluidStochForceDensity_m[3];
    fftw_complex               *fluidStochForceDensity_k[3];
    fftw_plan                   fluidStochForceDensity_DFT_plan[3];
    fftw_plan                   fluidStochForceDensity_IDFT_plan[3];

    fftw_complex               *fluidPressure_m;
    fftw_complex               *fluidPressure_k;
    fftw_plan                   fluidPressure_DFT_plan;
    fftw_plan                   fluidPressure_IDFT_plan;

    int                         flagUseFluidPressure;

    int                         flagComputeStress;
    double                     *controlPtsStress;

  #endif

  } SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ExtrasType;

  /* ======================== Variables ======================== */
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ParamsType *SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Params;
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_ExtrasType *SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_Extras;

  int                         flagWriteFluidVel_VTK;
  int                         flagWriteFluidForce_VTK;
  int                         flagWriteFluidPressure_VTK;


};

}

#endif

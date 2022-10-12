#ifndef EXTRINSIC_MODEL_TWO_TEMPERATURE
#define EXTRINSIC_MODEL_TWO_TEMPERATURE

#include <string>

#include "ExtrinsicModel.h"
#include "FieldEulerIntegrator.h"

namespace ATC {
  class ATC_Coupling;
  class PrescribedDataManager;
  class ExtrinsicModel;
  class PhysicsModel;

  /**
   *  @class  ExtrinsicModelTwoTemperature
   *  @brief  add electron temperature physics to phonon physics
   *          owned field: ELECTRON_TEMPERATURE
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelTwoTemperature
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ExtrinsicModelTwoTemperature : public ExtrinsicModel {

  public:

    // constructor
    ExtrinsicModelTwoTemperature(ExtrinsicModelManager * modelManager,
                   ExtrinsicModelType modelType,
                   std::string matFileName);

    // destructor
    virtual ~ExtrinsicModelTwoTemperature();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg);

    /** pre time integration */
    virtual void initialize();

    /** Predictor phase, executed before Verlet */
    virtual void pre_init_integrate();

    /** Set sources to AtC equation */
    virtual void set_sources(FIELDS & fields, FIELDS & sources);

    /** Add model-specific output data */
    virtual void output(OUTPUT_LIST & outputData);

    /** set up LAMMPS display variables */
    virtual int size_vector(int externalSize);

    /** get LAMMPS display variables */
    virtual bool compute_vector(int n, double & value);

  protected:
    /** electron time integration flag */
    TimeIntegrator::TimeIntegrationType electronTimeIntegration_;

    /** electron time integration flag */
    FieldEulerIntegrator * temperatureIntegrator_;

    /** number of electron time integration subscycles */
    int nsubcycle_;

    /** flag for turning off exchange during equilibration */
    bool exchangeFlag_;

    /** offset/size for LAMMPS display output */
    int baseSize_;

  };

};
#endif

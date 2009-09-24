#ifndef EXTRINSIC_MODEL
#define EXTRINSIC_MODEL

// ATC_Transfer headers
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

using namespace std;
namespace ATC {

  // forward declarations
  class ATC_Transfer;
  class ExtrinsicModel;
  class PhysicsModel;

  /** enumeration for the model types avaiable */
  enum ExtrinsicModelType {
    NO_MODEL=0,
    TWO_TEMPERATURE,
    NUM_MODELS
  };

  /**
   *  @class  ExtrinsicModelManager
   *  @brief  Handles parsing and parameter storage extrinsic models
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModelManager
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ExtrinsicModelManager {
  
  public:

    // constructor
    ExtrinsicModelManager(ATC_Transfer * atcTransfer);

    // destructor
    ~ExtrinsicModelManager();

    /** parser/modifier */
    bool modify(int narg, char **arg);

    /** create_model */
    void create_model(ExtrinsicModelType modelType, string matFileName);

    /** pre time integration */
    void initialize();

    /** set up LAMMPS display variables */
    int size_vector(int intrinsicSize);

    /** get LAMMPS display variables */
    double compute_vector(int n);

    /** post integration run */
    // is this called at end of run or simulation
    void finish();

    // calls during LAMMPS Velocity-Verlet integration
    /** Predictor phase, executed before Verlet */
    void pre_init_integrate(ExtrinsicModelType modelType = NUM_MODELS);
    /** Predictor phase, executed between velocity and position Verlet */
    void mid_init_integrate(ExtrinsicModelType modelType = NUM_MODELS);
    /** Predictor phase, executed after Verlet */
    void post_init_integrate(ExtrinsicModelType modelType = NUM_MODELS);
    /** Corrector phase, executed before Verlet */
    void pre_final_integrate(ExtrinsicModelType modelType = NUM_MODELS);
    /** Corrector phase, executed after Verlet*/
    void post_final_integrate(ExtrinsicModelType modelType = NUM_MODELS);

    /** get source terms for AtC equations */
    void set_sources(FIELDS & fields, FIELDS & sources,
                     ExtrinsicModelType modelType = NUM_MODELS);

    /** return output data to main AtC */
    void output(double dt,OUTPUT_LIST & outputData);

    /** model name enum to string */
    static bool model_to_string(const ExtrinsicModelType index, string & name) 
    {
      switch (index) {
        case NO_MODEL:
          name = "no_model";
          break;
        case TWO_TEMPERATURE:
          name = "two_temperature";
          break;
        default:
          return false;
          break;
      }
      return true; 
    };

    /** string to model enum */
    static bool string_to_model(const string & name, ExtrinsicModelType & index) 
    {
      if      (name=="no_model")
        index = NO_MODEL;
      else if (name=="two_temperature")
        index = TWO_TEMPERATURE;
      else
        return false;
      
      return true;
    };

    /** access to ATC transfer object */
    ATC_Transfer * get_atc_transfer() {return atcTransfer_;};

  protected:


    /** associated ATC_Transfer object */
    ATC_Transfer * atcTransfer_;

    /** equation hanlder */
    vector<ExtrinsicModel *> extrinsicModels_;

  private:
    ExtrinsicModelManager(); // DO NOT define this, only use constructor above

  };

  /**
   *  @class  ExtrinsicModel
   *  @brief  base class for functionality of all extrinsic models
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ExtrinsicModel
  //--------------------------------------------------------
  //--------------------------------------------------------

  class ExtrinsicModel {
  
  public:

    // constructor
    ExtrinsicModel(ExtrinsicModelManager * modelManager,
                   ExtrinsicModelType modelType,
                   string matFileName);

    // destructor
    virtual ~ExtrinsicModel();

    /** parser/modifier */
    virtual bool modify(int narg, char **arg) {return false;};

    /** pre time integration */
    virtual void initialize(){};

    /** set up LAMMPS display variables */
    virtual int size_vector(int externalSize) {return 0;};

    /** get LAMMPS display variables */
    virtual bool compute_vector(int n, double & value) {return false;};

    /** post integration run */
    // is this called at end of run or simulation
    virtual void finish(){};

    /** Predictor phase, executed before Verlet */
    virtual void pre_init_integrate(){};

    /** Predictor phase, executed between velocity and position Verlet */
    virtual void mid_init_integrate(){};

    /** Predictor phase, executed after Verlet */
    virtual void post_init_integrate(){};

    /** Corrector phase, executed before Verlet */
    virtual void pre_final_integrate(){};

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(){};

    /** Corrector phase, executed after Verlet*/
    virtual void post_final_integrate(){};

    /** Set sources to AtC equation */
    virtual void set_sources(FIELDS & fields, FIELDS & sources){};

    /** Add model-specific output data */
    virtual void output(double dt, OUTPUT_LIST & outputData){};

    /** get the fields and their sizes */
    void get_num_fields(map<FieldName,int> & fieldSizes);

    /** return the type of model being used */
    ExtrinsicModelType get_model_type() {return modelType_;};

  protected:

    ExtrinsicModel(){};

    /** ATC transfer object */
    ATC_Transfer * atc_;

    /** model manager object */
    ExtrinsicModelManager * modelManager_;

    /** tag for model type */
    ExtrinsicModelType modelType_;

    /** list of model fields in this model */
    std::map<FieldName, int> fieldSizes_;

    /** definition for physics used in this model */
    PhysicsModel * physicsModel_;

    /** rhs */
    FIELDS rhs_;

    /** rhs mask for coupling with MD */
    Array2D<bool> rhsMaskIntrinsic_;

    GRAD_FIELDS fluxes_;

  };

};
#endif

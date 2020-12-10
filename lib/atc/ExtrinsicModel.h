#ifndef EXTRINSIC_MODEL
#define EXTRINSIC_MODEL

#include <vector>
#include <map>
#include <string>

#include "ATC_TypeDefs.h"
#include "MatrixLibrary.h"

namespace ATC {
  class ATC_Coupling;
  class ExtrinsicModel;
  class PhysicsModel;

  /** enumeration for the model types available */
  enum ExtrinsicModelType {
    NO_MODEL=0,
    TWO_TEMPERATURE,
    DRIFT_DIFFUSION,
    DRIFT_DIFFUSION_EQUILIBRIUM,
    DRIFT_DIFFUSION_SCHRODINGER,
    DRIFT_DIFFUSION_SCHRODINGER_SLICE,
    CONVECTIVE_DRIFT_DIFFUSION,
    CONVECTIVE_DRIFT_DIFFUSION_EQUILIBRIUM,
    CONVECTIVE_DRIFT_DIFFUSION_SCHRODINGER,
    ELECTROSTATIC,
    ELECTROSTATIC_EQUILIBRIUM,
    FEM_EFIELD,
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
    ExtrinsicModelManager(ATC_Coupling * atcTransfer);

    // destructor
    ~ExtrinsicModelManager();

    /** parser/modifier */
    bool modify(int narg, char **arg);

    /** create_model */
    void create_model(ExtrinsicModelType modelType, std::string matFileName);

    /** construct the transfers needed by the model */
    void construct_transfers();

    /** pre time integration */
    void initialize();

    /** set up LAMMPS display variables */
    int size_vector(int intrinsicSize);

    /** get LAMMPS display variables */
    double compute_scalar(void);
    double compute_vector(int n);

    /** post integration run */
    // is this called at end of run or simulation
    void finish();

    // calls during LAMMPS Velocity-Verlet integration
    /** Predictor phase, executed before Verlet */
    void pre_init_integrate(ExtrinsicModelType modelType = NUM_MODELS);

    /** Predictor phase, executed after Verlet */
    void post_init_integrate(ExtrinsicModelType modelType = NUM_MODELS);

    /** Make changes to the forces lammps calculates */
    void post_force(ExtrinsicModelType modelType = NUM_MODELS);

    /** Corrector phase, executed before Verlet */
    void pre_final_integrate(ExtrinsicModelType modelType = NUM_MODELS);
    /** Corrector phase, executed after Verlet*/
    void post_final_integrate(ExtrinsicModelType modelType = NUM_MODELS);

    /** get source terms for AtC equations */
    void set_sources(FIELDS & fields, FIELDS & sources,
                     ExtrinsicModelType modelType = NUM_MODELS);

    /** return output data to main AtC */
    void output(OUTPUT_LIST & outputData);


    /** model name enum to string */
    static bool model_to_string(const ExtrinsicModelType index, std::string & name) 
    {
      switch (index) {
        case NO_MODEL:
          name = "no_model";
          break;
        case TWO_TEMPERATURE:
          name = "two_temperature";
          break;
        case DRIFT_DIFFUSION:
          name = "drift_diffusion";
          break;
        case DRIFT_DIFFUSION_EQUILIBRIUM:
          name = "drift_diffusion-equilibrium";
          break;
        case DRIFT_DIFFUSION_SCHRODINGER:
          name = "drift_diffusion-schrodinger";
          break;
        case DRIFT_DIFFUSION_SCHRODINGER_SLICE:
          name = "drift_diffusion-schrodinger-slice";
          break;
        case CONVECTIVE_DRIFT_DIFFUSION:
          name = "convective_drift_diffusion";
          break;
        case CONVECTIVE_DRIFT_DIFFUSION_EQUILIBRIUM:
          name = "convective_drift_diffusion-equilibrium";
          break;
        case CONVECTIVE_DRIFT_DIFFUSION_SCHRODINGER:
          name = "convective_drift_diffusion-schrodinger";
          break;
        case ELECTROSTATIC:
          name = "electrostatic";
          break;
        case ELECTROSTATIC_EQUILIBRIUM:
          name = "electrostatic-equilibrium";
          break;
        case FEM_EFIELD:
          name = "fem_efield";
          break;
        default:
          return false;
          break;
      }
      return true; 
    };

    /** string to model enum */
    static bool string_to_model(const std::string & name, ExtrinsicModelType & index) 
    {
      if      (name=="no_model")
        index = NO_MODEL;
      else if (name=="two_temperature")
        index = TWO_TEMPERATURE;
      else if (name=="drift_diffusion")
        index = DRIFT_DIFFUSION;
      else if (name=="drift_diffusion-equilibrium")
        index = DRIFT_DIFFUSION_EQUILIBRIUM;
      else if (name=="drift_diffusion-schrodinger")
        index = DRIFT_DIFFUSION_SCHRODINGER;
      else if (name=="drift_diffusion-schrodinger-slice")
        index = DRIFT_DIFFUSION_SCHRODINGER_SLICE;
      else if (name=="convective_drift_diffusion")
        index = CONVECTIVE_DRIFT_DIFFUSION;
      else if (name=="convective_drift_diffusion-equilibrium")
        index = CONVECTIVE_DRIFT_DIFFUSION_EQUILIBRIUM;
      else if (name=="convective_drift_diffusion-schrodinger")
        index = CONVECTIVE_DRIFT_DIFFUSION_SCHRODINGER;
      else if (name=="electrostatic")
        index = ELECTROSTATIC;
      else if (name=="electrostatic-equilibrium")
        index = ELECTROSTATIC_EQUILIBRIUM;
      else if (name=="fem_efield")
        index = FEM_EFIELD;

      else
        return false;
      
      return true;
    };

    /** access to ATC transfer object */
    ATC_Coupling * atc() {return atc_;};

    /** access to model of a specific type */
    const ExtrinsicModel * model(const ExtrinsicModelType type) const;

  protected:

    /** associated ATC_Coupling object */
    ATC_Coupling * atc_;

    /** equation handler */
    std::vector<ExtrinsicModel *> extrinsicModels_;

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
                   std::string matFileName);

    // destructor
    virtual ~ExtrinsicModel();

    /** parser/modifier */
    virtual bool modify(int /* narg */, char ** /* arg */) {return false;};

    /** construct transfers needed by the model */
    virtual void construct_transfers(){};

    /** pre time integration */
    virtual void initialize();

    /** set up LAMMPS display variables */
    virtual int size_vector(int /* externalSize */) {return 0;};

    /** get LAMMPS display variables */
    virtual double compute_scalar(void) { return 0.0; }
    virtual bool compute_vector(int /* n */, double & /* value */) {return false;};

    /** post integration run */
    // is this called at end of run or simulation
    virtual void finish(){};

    /** Predictor phase, executed before Verlet */
    virtual void pre_init_integrate(){};

    /** Predictor phase, executed after Verlet */
    virtual void post_init_integrate(){};

    /** changes to lammps forces */
    virtual void post_force(){};

    /** Corrector phase, executed before Verlet */
    virtual void pre_final_integrate(){};

    /** Corrector phase, Verlet second step for velocity */
    virtual void final_integrate(){};

    /** Corrector phase, executed after Verlet*/
    virtual void post_final_integrate(){};

    /** Set sources to AtC equation */
    virtual void set_sources(FIELDS & /* fields */, FIELDS & /* sources */){};

    /** Add model-specific output data */
    virtual void output(OUTPUT_LIST & /* outputData */){};

    /** get the fields and their sizes */
    void num_fields(std::map<FieldName,int> & fieldSizes);

    /** return the type of model being used */
    ExtrinsicModelType model_type() const {return modelType_;};

  protected:

    ExtrinsicModel(){};

    /** ATC transfer object */
    ATC_Coupling * atc_;

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

    
    
    GRAD_FIELD_MATS fluxes_;

    /** number of nodes */
    int nNodes_;

    /** number of spatial dimensions */
    int nsd_;

  };

}
#endif

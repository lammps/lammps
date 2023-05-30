// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvarscript.h"

#if (__cplusplus >= 201103L)
std::map<std::string, std::function<colvar::cvc* (const std::string& subcv_conf)>> colvar::global_cvc_map = std::map<std::string, std::function<colvar::cvc* (const std::string& subcv_conf)>>();
#endif

colvar::colvar()
{
  prev_timestep = -1L;
  after_restart = false;
  kinetic_energy = 0.0;
  potential_energy = 0.0;

#ifdef LEPTON
  dev_null = 0.0;
#endif

  expand_boundaries = false;

  description = "uninitialized colvar";
  colvar::init_dependencies();
}


/// Compare two cvcs using their names
/// Used to sort CVC array in scripted coordinates
bool colvar::compare_cvc(const colvar::cvc* const i, const colvar::cvc* const j)
{
  return i->name < j->name;
}


int colvar::init(std::string const &conf)
{
  cvm::log("Initializing a new collective variable.\n");
  colvarparse::set_string(conf);

  int error_code = COLVARS_OK;

  colvarmodule *cv = cvm::main();

  get_keyval(conf, "name", this->name,
             (std::string("colvar")+cvm::to_str(cv->variables()->size())));

  if ((cvm::colvar_by_name(this->name) != NULL) &&
      (cvm::colvar_by_name(this->name) != this)) {
    cvm::error("Error: this colvar cannot have the same name, \""+this->name+
                      "\", as another colvar.\n",
               COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  // Initialize dependency members
  // Could be a function defined in a different source file, for space?

  this->description = "colvar " + this->name;

  error_code |= init_components(conf);
  if (error_code != COLVARS_OK) {
    return cvm::get_error();
  }

  size_t i;

#ifdef LEPTON
  error_code |= init_custom_function(conf);
  if (error_code != COLVARS_OK) {
    return cvm::get_error();
  }
#endif

  // Setup colvar as scripted function of components
  if (get_keyval(conf, "scriptedFunction", scripted_function,
    "", colvarparse::parse_silent)) {

    enable(f_cv_scripted);
    cvm::log("This colvar uses scripted function \"" + scripted_function + "\".\n");
    cvm::main()->cite_feature("Scripted functions (Tcl)");

    std::string type_str;
    get_keyval(conf, "scriptedFunctionType", type_str, "scalar");

    x.type(colvarvalue::type_notset);
    int t;
    for (t = 0; t < colvarvalue::type_all; t++) {
      if (type_str == colvarvalue::type_keyword(colvarvalue::Type(t))) {
        x.type(colvarvalue::Type(t));
        break;
      }
    }
    if (x.type() == colvarvalue::type_notset) {
      cvm::error("Could not parse scripted colvar type.", COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }

    cvm::log(std::string("Expecting colvar value of type ")
      + colvarvalue::type_desc(x.type()));

    if (x.type() == colvarvalue::type_vector) {
      int size;
      if (!get_keyval(conf, "scriptedFunctionVectorSize", size)) {
        cvm::error("Error: no size specified for vector scripted function.",
                   COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
      x.vector1d_value.resize(size);
    }

    x_reported.type(x);

    // Sort array of cvcs based on their names
    // Note: default CVC names are in input order for same type of CVC
    std::sort(cvcs.begin(), cvcs.end(), colvar::compare_cvc);

    if(cvcs.size() > 1) {
      cvm::log("Sorted list of components for this scripted colvar:\n");
      for (i = 0; i < cvcs.size(); i++) {
        cvm::log(cvm::to_str(i+1) + " " + cvcs[i]->name);
      }
    }

    // Build ordered list of component values that will be
    // passed to the script
    for (i = 0; i < cvcs.size(); i++) {
      sorted_cvc_values.push_back(&(cvcs[i]->value()));
    }
  }

  if (!(is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function))) {
    colvarvalue const &cvc_value = (cvcs[0])->value();
    if (cvm::debug())
      cvm::log ("This collective variable is a "+
                colvarvalue::type_desc(cvc_value.type())+
                ((cvc_value.size() > 1) ? " with "+
                 cvm::to_str(cvc_value.size())+" individual components.\n" :
                 ".\n"));
    x.type(cvc_value);
    x_reported.type(cvc_value);
  }

  set_enabled(f_cv_scalar, (value().type() == colvarvalue::type_scalar));

  // If using scripted biases, any colvar may receive bias forces
  // and will need its gradient
  if (cvm::scripted_forces()) {
    enable(f_cv_gradient);
  }

  // check for linear combinations
  {
    bool lin = !(is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function));
    for (i = 0; i < cvcs.size(); i++) {

  //     FIXME this is a reverse dependency, ie. cv feature depends on cvc flag
  //     need to clarify this case
  //     if ((cvcs[i])->b_debug_gradients)
  //       enable(task_gradients);

      if ((cvcs[i])->sup_np != 1) {
        if (cvm::debug() && lin)
          cvm::log("Warning: You are using a non-linear polynomial "
                    "combination to define this collective variable, "
                    "some biasing methods may be unavailable.\n");
        lin = false;

        if ((cvcs[i])->sup_np < 0) {
          cvm::log("Warning: you chose a negative exponent in the combination; "
                    "if you apply forces, the simulation may become unstable "
                    "when the component \""+
                    (cvcs[i])->function_type+"\" approaches zero.\n");
        }
      }
    }
    set_enabled(f_cv_linear, lin);
  }

  // Colvar is homogeneous if:
  // - it is linear (hence not scripted)
  // - all cvcs have coefficient 1 or -1
  // i.e. sum or difference of cvcs
  {
    bool homogeneous = is_enabled(f_cv_linear);
    for (i = 0; i < cvcs.size(); i++) {
      if (cvm::fabs(cvm::fabs(cvcs[i]->sup_coeff) - 1.0) > 1.0e-10) {
        homogeneous = false;
      }
    }
    set_enabled(f_cv_homogeneous, homogeneous);
  }

  // A single-component variable almost concides with its CVC object
  if ((cvcs.size() == 1) && is_enabled(f_cv_homogeneous)) {
    if ( !is_enabled(f_cv_scripted) && !is_enabled(f_cv_custom_function) &&
         (cvm::fabs(cvcs[0]->sup_coeff - 1.0) < 1.0e-10) &&
         (cvcs[0]->sup_np == 1) ) {
      enable(f_cv_single_cvc);
    }
  }

  // Colvar is deemed periodic if:
  // - it is homogeneous
  // - all cvcs are periodic
  // - all cvcs have the same period
  if (is_enabled(f_cv_homogeneous) && cvcs[0]->is_enabled(f_cvc_periodic)) {
    bool b_periodic = true;
    period = cvcs[0]->period;
    wrap_center = cvcs[0]->wrap_center;
    for (i = 1; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled(f_cvc_periodic) || cvcs[i]->period != period) {
        b_periodic = false;
        period = 0.0;
        cvm::log("Warning: although one component is periodic, this colvar will "
                 "not be treated as periodic, either because the exponent is not "
                 "1, or because components of different periodicity are defined.  "
                 "Make sure that you know what you are doing!");
      }
    }
    set_enabled(f_cv_periodic, b_periodic);
  }

  // Allow scripted/custom functions to be defined as periodic
  if ( (is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function)) && is_enabled(f_cv_scalar) ) {
    if (get_keyval(conf, "period", period, 0.)) {
      enable(f_cv_periodic);
      get_keyval(conf, "wrapAround", wrap_center, 0.);
    }
  }

  // check that cvcs are compatible

  for (i = 0; i < cvcs.size(); i++) {

    // components may have different types only for scripted functions
    if (!(is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function)) && (colvarvalue::check_types(cvcs[i]->value(),
                                                                cvcs[0]->value())) ) {
      cvm::error("ERROR: you are defining this collective variable "
                 "by using components of different types. "
                 "You must use the same type in order to "
                 "sum them together.\n", COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
  }

  active_cvc_square_norm = 0.;
  for (i = 0; i < cvcs.size(); i++) {
    active_cvc_square_norm += cvcs[i]->sup_coeff * cvcs[i]->sup_coeff;
  }

  // at this point, the colvar's type is defined
  f.type(value());

  x_old.type(value());
  v_fdiff.type(value());
  v_reported.type(value());
  fj.type(value());
  ft.type(value());
  ft_reported.type(value());
  f_old.type(value());
  f_old.reset();

  x_restart.type(value());

  reset_bias_force();

  get_keyval(conf, "timeStepFactor", time_step_factor, 1);
  if (time_step_factor < 0) {
    cvm::error("Error: timeStepFactor must be positive.\n");
    return COLVARS_ERROR;
  }
  if (time_step_factor != 1) {
    enable(f_cv_multiple_ts);
  }

  error_code |= init_grid_parameters(conf);

  // Detect if we have a single component that is an alchemical lambda
  if (is_enabled(f_cv_single_cvc) && cvcs[0]->function_type == "alchLambda") {
    enable(f_cv_external);
  }

  error_code |= init_extended_Lagrangian(conf);
  error_code |= init_output_flags(conf);

  // Now that the children are defined we can solve dependencies
  enable(f_cv_active);

  error_code |= parse_analysis(conf);

  if (cvm::debug())
    cvm::log("Done initializing collective variable \""+this->name+"\".\n");

  return error_code;
}


#ifdef LEPTON
int colvar::init_custom_function(std::string const &conf)
{
  std::string expr, expr_in; // expr_in is a buffer to remember expr after unsuccessful parsing
  std::vector<Lepton::ParsedExpression> pexprs;
  Lepton::ParsedExpression pexpr;
  size_t pos = 0; // current position in config string
  double *ref;

  if (!key_lookup(conf, "customFunction", &expr_in, &pos)) {
    return COLVARS_OK;
  }

  cvm::main()->cite_feature("Custom functions (Lepton)");

  enable(f_cv_custom_function);
  cvm::log("This colvar uses a custom function.\n");

  do {
    expr = expr_in;
    if (cvm::debug())
      cvm::log("Parsing expression \"" + expr + "\".\n");
    try {
      pexpr = Lepton::Parser::parse(expr);
      pexprs.push_back(pexpr);
    }
    catch (...) {
      cvm::error("Error parsing expression \"" + expr + "\".\n", COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }

    try {
      value_evaluators.push_back(
          new Lepton::CompiledExpression(pexpr.createCompiledExpression()));
      // Define variables for cvc values
      // Stored in order: expr1, cvc1, cvc2, expr2, cvc1...
      for (size_t i = 0; i < cvcs.size(); i++) {
        for (size_t j = 0; j < cvcs[i]->value().size(); j++) {
          std::string vn = cvcs[i]->name +
              (cvcs[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
          try {
            ref =&value_evaluators.back()->getVariableReference(vn);
          }
          catch (...) { // Variable is absent from expression
            // To keep the same workflow, we use a pointer to a double here
            // that will receive CVC values - even though none was allocated by Lepton
            ref = &dev_null;
            cvm::log("Warning: Variable " + vn + " is absent from expression \"" + expr + "\".\n");
          }
          value_eval_var_refs.push_back(ref);
        }
      }
    }
    catch (...) {
      cvm::error("Error compiling expression \"" + expr + "\".\n", COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
  } while (key_lookup(conf, "customFunction", &expr_in, &pos));


  // Now define derivative with respect to each scalar sub-component
  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t j = 0; j < cvcs[i]->value().size(); j++) {
      std::string vn = cvcs[i]->name +
          (cvcs[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
      // Element ordering: we want the
      // gradient vector of derivatives of all elements of the colvar
      // wrt to a given element of a cvc ([i][j])
      for (size_t c = 0; c < pexprs.size(); c++) {
        gradient_evaluators.push_back(
            new Lepton::CompiledExpression(pexprs[c].differentiate(vn).createCompiledExpression()));
        // and record the refs to each variable in those expressions
        for (size_t k = 0; k < cvcs.size(); k++) {
          for (size_t l = 0; l < cvcs[k]->value().size(); l++) {
            std::string vvn = cvcs[k]->name +
                (cvcs[k]->value().size() > 1 ? cvm::to_str(l+1) : "");
            try {
              ref = &gradient_evaluators.back()->getVariableReference(vvn);
            }
            catch (...) { // Variable is absent from derivative
              // To keep the same workflow, we use a pointer to a double here
              // that will receive CVC values - even though none was allocated by Lepton
              if (cvm::debug()) {
                cvm::log("Warning: Variable " + vvn + " is absent from derivative of \"" + expr + "\" wrt " + vn + ".\n");
              }
              ref = &dev_null;
            }
            grad_eval_var_refs.push_back(ref);
          }
        }
      }
    }
  }


  if (value_evaluators.size() == 0) {
    cvm::error("Error: no custom function defined.\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  std::string type_str;
  bool b_type_specified = get_keyval(conf, "customFunctionType",
                                     type_str, "scalar", parse_silent);
  x.type(colvarvalue::type_notset);
  int t;
  for (t = 0; t < colvarvalue::type_all; t++) {
    if (type_str == colvarvalue::type_keyword(colvarvalue::Type(t))) {
      x.type(colvarvalue::Type(t));
      break;
    }
  }
  if (x.type() == colvarvalue::type_notset) {
    cvm::error("Could not parse custom colvar type.", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  // Guess type based on number of expressions
  if (!b_type_specified) {
    if (value_evaluators.size() == 1) {
      x.type(colvarvalue::type_scalar);
    } else {
      x.type(colvarvalue::type_vector);
    }
  }

  if (x.type() == colvarvalue::type_vector) {
    x.vector1d_value.resize(value_evaluators.size());
  }

  x_reported.type(x);
  cvm::log(std::string("Expecting colvar value of type ")
    + colvarvalue::type_desc(x.type())
    + (x.type()==colvarvalue::type_vector ? " of size " + cvm::to_str(x.size()) : "")
    + ".\n");

  if (x.size() != value_evaluators.size()) {
    cvm::error("Error: based on custom function type, expected "
               + cvm::to_str(x.size()) + " scalar expressions, but "
               + cvm::to_str(value_evaluators.size()) + " were found.\n");
    return COLVARS_INPUT_ERROR;
  }

  return COLVARS_OK;
}

#else

int colvar::init_custom_function(std::string const &conf)
{

  std::string expr;
  size_t pos = 0;
  if (key_lookup(conf, "customFunction", &expr, &pos)) {
    std::string msg("Error: customFunction requires the Lepton library.");
#if (__cplusplus < 201103L)
    // NOTE: this is not ideal; testing for the Lepton library's version would
    // be more accurate, but also less portable
    msg +=
      std::string("  Note also that recent versions of Lepton require C++11: "
                  "please see https://colvars.github.io/README-c++11.html.");
#endif
    return cvm::error(msg, COLVARS_NOT_IMPLEMENTED);
  }

  return COLVARS_OK;
}

#endif // #ifdef LEPTON


int colvar::init_grid_parameters(std::string const &conf)
{
  int error_code = COLVARS_OK;

  colvarmodule *cv = cvm::main();

  cvm::real default_width = width;

  if (!key_already_set("width")) {
    // The first time, check if the CVC has a width to provide
    default_width = 1.0;
    if (is_enabled(f_cv_single_cvc) && cvcs[0]->is_enabled(f_cvc_width)) {
      cvm::real const cvc_width = cvcs[0]->get_param("width");
      default_width = cvc_width;
    }
  }

  get_keyval(conf, "width", width, default_width);

  if (width <= 0.0) {
    cvm::error("Error: \"width\" must be positive.\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  lower_boundary.type(value());
  upper_boundary.type(value());
  lower_boundary.real_value = 0.0;
  upper_boundary.real_value = width; // Default to 1-wide grids

  if (is_enabled(f_cv_scalar)) {

    if (is_enabled(f_cv_single_cvc)) {
      // Get the default boundaries from the component
      if (cvcs[0]->is_enabled(f_cvc_lower_boundary)) {
        enable(f_cv_lower_boundary);
        enable(f_cv_hard_lower_boundary);
        lower_boundary =
          *(reinterpret_cast<colvarvalue const *>(cvcs[0]->get_param_ptr("lowerBoundary")));
      }
      if (cvcs[0]->is_enabled(f_cvc_upper_boundary)) {
        enable(f_cv_upper_boundary);
        enable(f_cv_hard_upper_boundary);
        upper_boundary =
          *(reinterpret_cast<colvarvalue const *>(cvcs[0]->get_param_ptr("upperBoundary")));
      }
    }

    if (get_keyval(conf, "lowerBoundary", lower_boundary, lower_boundary)) {
      enable(f_cv_lower_boundary);
      // Because this is the user's choice, we cannot assume it is a true
      // physical boundary
      disable(f_cv_hard_lower_boundary);
    }

    if (get_keyval(conf, "upperBoundary", upper_boundary, upper_boundary)) {
      enable(f_cv_upper_boundary);
      disable(f_cv_hard_upper_boundary);
    }

    // Parse legacy wall options and set up a harmonicWalls bias if needed
    cvm::real lower_wall_k = 0.0, upper_wall_k = 0.0;
    cvm::real lower_wall = 0.0, upper_wall = 0.0;
    std::string lw_conf, uw_conf;

    if (get_keyval(conf, "lowerWallConstant", lower_wall_k, 0.0,
                   parse_silent)) {
      cvm::log("Reading legacy options lowerWall and lowerWallConstant: "
               "consider using a harmonicWalls restraint (caution: force constant would then be scaled by width^2).\n");
      if (!get_keyval(conf, "lowerWall", lower_wall)) {
        error_code |= cvm::error("Error: the value of lowerWall must be set "
                                 "explicitly.\n", COLVARS_INPUT_ERROR);
      }
      lw_conf = std::string("\n\
    lowerWallConstant "+cvm::to_str(lower_wall_k*width*width)+"\n\
    lowerWalls "+cvm::to_str(lower_wall)+"\n");
    }

    if (get_keyval(conf, "upperWallConstant", upper_wall_k, 0.0,
                   parse_silent)) {
      cvm::log("Reading legacy options upperWall and upperWallConstant: "
               "consider using a harmonicWalls restraint (caution: force constant would then be scaled by width^2).\n");
      if (!get_keyval(conf, "upperWall", upper_wall)) {
        error_code |= cvm::error("Error: the value of upperWall must be set "
                                 "explicitly.\n", COLVARS_INPUT_ERROR);
      }
      uw_conf = std::string("\n\
    upperWallConstant "+cvm::to_str(upper_wall_k*width*width)+"\n\
    upperWalls "+cvm::to_str(upper_wall)+"\n");
    }

    if (lw_conf.size() && uw_conf.size()) {
      if (lower_wall >= upper_wall) {
        error_code |= cvm::error("Error: the upper wall, "+
                                 cvm::to_str(upper_wall)+
                                 ", is not higher than the lower wall, "+
                                 cvm::to_str(lower_wall)+".\n",
                                 COLVARS_INPUT_ERROR);
      }
    }

    if (lw_conf.size() || uw_conf.size()) {
      cvm::log("Generating a new harmonicWalls bias for compatibility purposes.\n");
      std::string const walls_conf("\n\
harmonicWalls {\n\
    name "+this->name+"w\n\
    colvars "+this->name+"\n"+lw_conf+uw_conf+"\
    timeStepFactor "+cvm::to_str(time_step_factor)+"\n"+
                             "}\n");
      error_code |= cv->append_new_config(walls_conf);
    }
  }

  get_keyval_feature(this, conf, "hardLowerBoundary", f_cv_hard_lower_boundary,
                     is_enabled(f_cv_hard_lower_boundary));

  get_keyval_feature(this, conf, "hardUpperBoundary", f_cv_hard_upper_boundary,
                     is_enabled(f_cv_hard_upper_boundary));

  // consistency checks for boundaries and walls
  if (is_enabled(f_cv_lower_boundary) && is_enabled(f_cv_upper_boundary)) {
    if (lower_boundary >= upper_boundary) {
      error_code |= cvm::error("Error: the upper boundary, "+
                               cvm::to_str(upper_boundary)+
                               ", is not higher than the lower boundary, "+
                               cvm::to_str(lower_boundary)+".\n",
                               COLVARS_INPUT_ERROR);
    }
  }

  get_keyval(conf, "expandBoundaries", expand_boundaries, expand_boundaries);
  if (expand_boundaries && periodic_boundaries()) {
    error_code |= cvm::error("Error: trying to expand boundaries that already "
                             "cover a whole period of a periodic colvar.\n",
                             COLVARS_INPUT_ERROR);
  }

  if (expand_boundaries && is_enabled(f_cv_hard_lower_boundary) &&
      is_enabled(f_cv_hard_upper_boundary)) {
    error_code |= cvm::error("Error: inconsistent configuration "
                             "(trying to expand boundaries, but both "
                             "hardLowerBoundary and hardUpperBoundary "
                             "are enabled).\n", COLVARS_INPUT_ERROR);
  }

  return error_code;
}


int colvar::init_extended_Lagrangian(std::string const &conf)
{
  colvarproxy *proxy = cvm::main()->proxy;
  get_keyval_feature(this, conf, "extendedLagrangian", f_cv_extended_Lagrangian, false);

  if (is_enabled(f_cv_extended_Lagrangian)) {
    cvm::real temp, tolerance, extended_period;

    cvm::log("Enabling the extended Lagrangian term for colvar \""+
             this->name+"\".\n");

    // Mark x_ext as uninitialized so we can initialize it to the colvar value when updating
    x_ext.type(colvarvalue::type_notset);
    v_ext.type(value());
    fr.type(value());
    const bool temp_provided = get_keyval(conf, "extendedTemp", temp,
                                          proxy->target_temperature());
    if (is_enabled(f_cv_external)) {
      // In the case of an "external" coordinate, there is no coupling potential:
      // only the fictitious mass is meaningful
      get_keyval(conf, "extendedMass", ext_mass);
      // Ensure that the computed restraint energy term is zero
      ext_force_k = 0.0;
    } else {
      // Standard case of coupling to a geometric colvar
      if (temp <= 0.0) { // Then a finite temperature is required
        if (temp_provided)
          cvm::error("Error: \"extendedTemp\" must be positive.\n", COLVARS_INPUT_ERROR);
        else
          cvm::error("Error: a positive temperature must be provided, either "
                    "by enabling a thermostat, or through \"extendedTemp\".\n",
                    COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
      get_keyval(conf, "extendedFluctuation", tolerance);
      if (tolerance <= 0.0) {
        cvm::error("Error: \"extendedFluctuation\" must be positive.\n", COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
      ext_force_k = proxy->boltzmann() * temp / (tolerance * tolerance);
      cvm::log("Computed extended system force constant: " + cvm::to_str(ext_force_k) + " [E]/U^2\n");

      get_keyval(conf, "extendedTimeConstant", extended_period, 200.0);
      if (extended_period <= 0.0) {
        cvm::error("Error: \"extendedTimeConstant\" must be positive.\n", COLVARS_INPUT_ERROR);
      }
      ext_mass = (proxy->boltzmann() * temp * extended_period * extended_period)
        / (4.0 * PI * PI * tolerance * tolerance);
      cvm::log("Computed fictitious mass: " + cvm::to_str(ext_mass) + " [E]/(U/fs)^2   (U: colvar unit)\n");
    }
    {
      bool b_output_energy;
      get_keyval(conf, "outputEnergy", b_output_energy, false);
      if (b_output_energy) {
        enable(f_cv_output_energy);
      }
    }

    get_keyval(conf, "extendedLangevinDamping", ext_gamma, 1.0);
    if (ext_gamma < 0.0) {
      cvm::error("Error: \"extendedLangevinDamping\" may not be negative.\n", COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
    if (ext_gamma != 0.0) {
      enable(f_cv_Langevin);
      ext_gamma *= 1.0e-3; // correct as long as input is required in ps-1 and cvm::dt() is in fs
      // Adjust Langevin sigma for slow time step if time_step_factor != 1
      ext_sigma = cvm::sqrt(2.0 * proxy->boltzmann() * temp * ext_gamma * ext_mass / (cvm::dt() * cvm::real(time_step_factor)));
    }

    get_keyval_feature(this, conf, "reflectingLowerBoundary", f_cv_reflecting_lower_boundary, false);
    get_keyval_feature(this, conf, "reflectingUpperBoundary", f_cv_reflecting_upper_boundary, false);
  }

  return COLVARS_OK;
}


int colvar::init_output_flags(std::string const &conf)
{
  {
    bool b_output_value;
    get_keyval(conf, "outputValue", b_output_value, true);
    if (b_output_value) {
      enable(f_cv_output_value);
    }
  }

  {
    bool b_output_velocity;
    get_keyval(conf, "outputVelocity", b_output_velocity, false);
    if (b_output_velocity) {
      enable(f_cv_output_velocity);
    }
  }

  {
    bool temp;
    if (get_keyval(conf, "outputSystemForce", temp, false, colvarparse::parse_silent)) {
      cvm::error("Option outputSystemForce is deprecated: only outputTotalForce is supported instead.\n"
                 "The two are NOT identical: see https://colvars.github.io/totalforce.html.\n", COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
  }

  get_keyval_feature(this, conf, "outputTotalForce", f_cv_output_total_force, false);
  get_keyval_feature(this, conf, "outputAppliedForce", f_cv_output_applied_force, false);
  get_keyval_feature(this, conf, "subtractAppliedForce", f_cv_subtract_applied_force, false);

  return COLVARS_OK;
}

#if (__cplusplus >= 201103L)
// C++11
template<typename def_class_name> int colvar::init_components_type(std::string const &,
                                                                   char const * /* def_desc */,
                                                                   char const *def_config_key) {
  // global_cvc_map is only supported in the C++11 case
  global_cvc_map[def_config_key] = [](const std::string& cvc_conf){return new def_class_name(cvc_conf);};
  // TODO: maybe it is better to do more check to avoid duplication in the map?
  return COLVARS_OK;
}

int colvar::init_components_type_from_global_map(const std::string& conf,
                                                 const char* def_config_key) {
#else
template<typename def_class_name> int colvar::init_components_type(std::string const & conf,
                                                                   char const * /* def_desc */,
                                                                   char const *def_config_key) {
#endif
  size_t def_count = 0;
  std::string def_conf = "";
  size_t pos = 0;
  while ( this->key_lookup(conf,
                           def_config_key,
                           &def_conf,
                           &pos) ) {
    if (!def_conf.size()) continue;
    cvm::log("Initializing "
             "a new \""+std::string(def_config_key)+"\" component"+
             (cvm::debug() ? ", with configuration:\n"+def_conf
              : ".\n"));
    cvm::increase_depth();
    // only the following line is different from init_components_type
    // in the non-C++11 case
#if (__cplusplus >= 201103L)
    cvc *cvcp = global_cvc_map.at(def_config_key)(def_conf);
#else
    cvc *cvcp = new def_class_name(def_conf);
#endif
    if (cvcp != NULL) {
      cvcs.push_back(cvcp);
      cvcp->check_keywords(def_conf, def_config_key);
      cvcp->set_function_type(def_config_key);
      if (cvm::get_error()) {
        cvm::error("Error: in setting up component \""+
                   std::string(def_config_key)+"\".\n", COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
      cvm::decrease_depth();
    } else {
      cvm::decrease_depth();
      cvm::error("Error: in allocating component \""+
                   std::string(def_config_key)+"\".\n",
                 COLVARS_MEMORY_ERROR);
      return COLVARS_MEMORY_ERROR;
    }

    if ( (cvcp->period != 0.0) || (cvcp->wrap_center != 0.0) ) {
      if (! cvcp->is_enabled(f_cvc_periodic)) {
        cvm::error("Error: invalid use of period and/or "
                   "wrapAround in a \""+
                   std::string(def_config_key)+
                   "\" component.\n"+
                   "Period: "+cvm::to_str(cvcp->period) +
                   " wrapAround: "+cvm::to_str(cvcp->wrap_center),
                   COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
    }

    if ( ! cvcs.back()->name.size()) {
      std::ostringstream s;
      s << def_config_key << std::setfill('0') << std::setw(4) << ++def_count;
      cvcs.back()->name = s.str();
      /* pad cvc number for correct ordering when sorting by name */
    }

    cvcs.back()->setup();
    if (cvm::debug()) {
      cvm::log("Done initializing a \""+
               std::string(def_config_key)+
               "\" component"+
               (cvm::debug() ?
                ", named \""+cvcs.back()->name+"\""
                : "")+".\n");
    }
    def_conf = "";
    if (cvm::debug()) {
      cvm::log("Parsed "+cvm::to_str(cvcs.size())+
               " components at this time.\n");
    }
  }

  return COLVARS_OK;
}

int colvar::init_components(std::string const &conf)
{
  int error_code = COLVARS_OK;
  size_t i = 0, j = 0;

  // in the non-C++11 case, the components are initialized directly by init_components_type;
  // in the C++11 case, the components are stored in the global_cvc_map at first
  // by init_components_type, and then the map is iterated to initialize all components.
  error_code |= init_components_type<distance>(conf, "distance", "distance");
  error_code |= init_components_type<distance_vec>(conf, "distance vector", "distanceVec");
  error_code |= init_components_type<cartesian>(conf, "Cartesian coordinates", "cartesian");
  error_code |= init_components_type<distance_dir>(conf, "distance vector "
    "direction", "distanceDir");
  error_code |= init_components_type<distance_z>(conf, "distance projection "
    "on an axis", "distanceZ");
  error_code |= init_components_type<distance_xy>(conf, "distance projection "
    "on a plane", "distanceXY");
  error_code |= init_components_type<polar_theta>(conf, "spherical polar angle theta",
    "polarTheta");
  error_code |= init_components_type<polar_phi>(conf, "spherical azimuthal angle phi",
    "polarPhi");
  error_code |= init_components_type<distance_inv>(conf, "average distance "
    "weighted by inverse power", "distanceInv");
  error_code |= init_components_type<distance_pairs>(conf, "N1xN2-long vector "
    "of pairwise distances", "distancePairs");
  error_code |= init_components_type<dipole_magnitude>(conf, "dipole magnitude",
    "dipoleMagnitude");
  error_code |= init_components_type<coordnum>(conf, "coordination "
    "number", "coordNum");
  error_code |= init_components_type<selfcoordnum>(conf, "self-coordination "
    "number", "selfCoordNum");
  error_code |= init_components_type<groupcoordnum>(conf, "group-coordination "
    "number", "groupCoord");
  error_code |= init_components_type<angle>(conf, "angle", "angle");
  error_code |= init_components_type<dipole_angle>(conf, "dipole angle", "dipoleAngle");
  error_code |= init_components_type<dihedral>(conf, "dihedral", "dihedral");
  error_code |= init_components_type<h_bond>(conf, "hydrogen bond", "hBond");
  error_code |= init_components_type<alpha_angles>(conf, "alpha helix", "alpha");
  error_code |= init_components_type<dihedPC>(conf, "dihedral "
    "principal component", "dihedralPC");
  error_code |= init_components_type<orientation>(conf, "orientation", "orientation");
  error_code |= init_components_type<orientation_angle>(conf, "orientation "
    "angle", "orientationAngle");
  error_code |= init_components_type<orientation_proj>(conf, "orientation "
    "projection", "orientationProj");
  error_code |= init_components_type<tilt>(conf, "tilt", "tilt");
  error_code |= init_components_type<spin_angle>(conf, "spin angle", "spinAngle");
  error_code |= init_components_type<rmsd>(conf, "RMSD", "rmsd");
  error_code |= init_components_type<gyration>(conf, "radius of "
    "gyration", "gyration");
  error_code |= init_components_type<inertia>(conf, "moment of "
    "inertia", "inertia");
  error_code |= init_components_type<inertia_z>(conf, "moment of inertia around an axis", "inertiaZ");
  error_code |= init_components_type<eigenvector>(conf, "eigenvector", "eigenvector");
  error_code |= init_components_type<alch_lambda>(conf, "alchemical coupling parameter", "alchLambda");
  error_code |= init_components_type<alch_Flambda>(conf, "force on alchemical coupling parameter", "alchFLambda");
  error_code |= init_components_type<gspath>(conf, "geometrical path collective variables (s)", "gspath");
  error_code |= init_components_type<gzpath>(conf, "geometrical path collective variables (z)", "gzpath");
  error_code |= init_components_type<linearCombination>(conf, "linear combination of other collective variables", "linearCombination");
  error_code |= init_components_type<gspathCV>(conf, "geometrical path collective variables (s) for other CVs", "gspathCV");
  error_code |= init_components_type<gzpathCV>(conf, "geometrical path collective variables (z) for other CVs", "gzpathCV");
  error_code |= init_components_type<aspathCV>(conf, "arithmetic path collective variables (s) for other CVs", "aspathCV");
  error_code |= init_components_type<azpathCV>(conf, "arithmetic path collective variables (s) for other CVs", "azpathCV");
  error_code |= init_components_type<euler_phi>(conf, "euler phi angle of the optimal orientation", "eulerPhi");
  error_code |= init_components_type<euler_psi>(conf, "euler psi angle of the optimal orientation", "eulerPsi");
  error_code |= init_components_type<euler_theta>(conf, "euler theta angle of the optimal orientation", "eulerTheta");
#ifdef LEPTON
  error_code |= init_components_type<customColvar>(conf, "CV with support of the lepton custom function", "customColvar");
#endif
  error_code |= init_components_type<neuralNetwork>(conf, "neural network CV for other CVs", "NeuralNetwork");

  error_code |= init_components_type<map_total>(conf, "total value of atomic map", "mapTotal");
#if (__cplusplus >= 201103L)
  // iterate over all available CVC in the map
  for (auto it = global_cvc_map.begin(); it != global_cvc_map.end(); ++it) {
    error_code |= init_components_type_from_global_map(conf, it->first.c_str());
    // TODO: is it better to check the error code here?
    if (error_code != COLVARS_OK) {
      cvm::log("Failed to initialize " + it->first + " with the following configuration:\n");
      cvm::log(conf);
      // TODO: should it stop here?
    }
  }
#endif
  if (!cvcs.size() || (error_code != COLVARS_OK)) {
    cvm::error("Error: no valid components were provided "
               "for this collective variable.\n",
               COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  // Check for uniqueness of CVC names (esp. if user-provided)
  for (i = 0; i < cvcs.size(); i++) {
    for (j = i+1; j < cvcs.size(); j++) {
      if (cvcs[i]->name == cvcs[j]->name) {
        cvm::error("Components " + cvm::to_str(i) + " and " + cvm::to_str(j) +\
          " cannot have the same name \"" +  cvcs[i]->name+ "\".\n", COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
    }
  }

  n_active_cvcs = cvcs.size();

  // Store list of children cvcs for dependency checking purposes
  for (i = 0; i < cvcs.size(); i++) {
    add_child(cvcs[i]);
  }

  cvm::log("All components initialized.\n");

  return COLVARS_OK;
}


void colvar::do_feature_side_effects(int id)
{
  switch (id) {
    case f_cv_total_force_calc:
      cvm::request_total_force();
      break;
    case f_cv_collect_atom_ids:
      // Needed for getting gradients vias collect_gradients
      // or via atomic forces e.g. in Colvars Dashboard in VMD
      if (atom_ids.size() == 0) {
        build_atom_list();
      }
      break;
  }
}


void colvar::build_atom_list(void)
{
  // If atomic gradients are requested, build full list of atom ids from all cvcs
  std::list<int> temp_id_list;

  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t j = 0; j < cvcs[i]->atom_groups.size(); j++) {
      cvm::atom_group const &ag = *(cvcs[i]->atom_groups[j]);
      for (size_t k = 0; k < ag.size(); k++) {
        temp_id_list.push_back(ag[k].id);
      }
      if (ag.is_enabled(f_ag_fitting_group) && ag.is_enabled(f_ag_fit_gradients)) {
        cvm::atom_group const &fg = *(ag.fitting_group);
        for (size_t k = 0; k < fg.size(); k++) {
          temp_id_list.push_back(fg[k].id);
        }
      }
    }
  }

  temp_id_list.sort();
  temp_id_list.unique();

  std::list<int>::iterator li;
  for (li = temp_id_list.begin(); li != temp_id_list.end(); ++li) {
    atom_ids.push_back(*li);
  }

  temp_id_list.clear();

  atomic_gradients.resize(atom_ids.size());
  if (atom_ids.size()) {
    if (cvm::debug())
      cvm::log("Colvar: created atom list with " + cvm::to_str(atom_ids.size()) + " atoms.\n");
  } else {
    cvm::log("Warning: colvar components communicated no atom IDs.\n");
  }
}


int colvar::parse_analysis(std::string const &conf)
{

  //   if (cvm::debug())
  //     cvm::log ("Parsing analysis flags for collective variable \""+
  //               this->name+"\".\n");

  runave_length = 0;
  bool b_runave = false;
  if (get_keyval(conf, "runAve", b_runave) && b_runave) {

    enable(f_cv_runave);

    get_keyval(conf, "runAveLength", runave_length, 1000);
    get_keyval(conf, "runAveStride", runave_stride, 1);

    if ((cvm::restart_out_freq % runave_stride) != 0) {
      cvm::error("Error: runAveStride must be commensurate with the restart frequency.\n", COLVARS_INPUT_ERROR);
    }

    get_keyval(conf, "runAveOutputFile", runave_outfile, runave_outfile);
  }

  acf_length = 0;
  bool b_acf = false;
  if (get_keyval(conf, "corrFunc", b_acf) && b_acf) {

    enable(f_cv_corrfunc);

    get_keyval(conf, "corrFuncWithColvar", acf_colvar_name, this->name);
    if (acf_colvar_name == this->name) {
      cvm::log("Calculating auto-correlation function.\n");
    } else {
      cvm::log("Calculating correlation function with \""+
                this->name+"\".\n");
    }

    std::string acf_type_str;
    get_keyval(conf, "corrFuncType", acf_type_str, to_lower_cppstr(std::string("velocity")));
    if (acf_type_str == to_lower_cppstr(std::string("coordinate"))) {
      acf_type = acf_coor;
    } else if (acf_type_str == to_lower_cppstr(std::string("velocity"))) {
      acf_type = acf_vel;
      enable(f_cv_fdiff_velocity);
      colvar *cv2 = cvm::colvar_by_name(acf_colvar_name);
      if (cv2 == NULL) {
        return cvm::error("Error: collective variable \""+acf_colvar_name+
                          "\" is not defined at this time.\n", COLVARS_INPUT_ERROR);
      }
      cv2->enable(f_cv_fdiff_velocity); // Manual dependency to object of same type
    } else if (acf_type_str == to_lower_cppstr(std::string("coordinate_p2"))) {
      acf_type = acf_p2coor;
    } else {
      cvm::log("Unknown type of correlation function, \""+
                        acf_type_str+"\".\n");
      cvm::set_error_bits(COLVARS_INPUT_ERROR);
    }

    get_keyval(conf, "corrFuncOffset", acf_offset, 0);
    get_keyval(conf, "corrFuncLength", acf_length, 1000);
    get_keyval(conf, "corrFuncStride", acf_stride, 1);

    if ((cvm::restart_out_freq % acf_stride) != 0) {
      cvm::error("Error: corrFuncStride must be commensurate with the restart frequency.\n", COLVARS_INPUT_ERROR);
    }

    get_keyval(conf, "corrFuncNormalize", acf_normalize, true);
    get_keyval(conf, "corrFuncOutputFile", acf_outfile, acf_outfile);
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvar::init_dependencies() {
  size_t i;
  if (features().size() == 0) {
    for (i = 0; i < f_cv_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_cv_active, "active", f_type_dynamic);
    // Do not require f_cvc_active in children, as some components may be disabled
    // Colvars must be either a linear combination, or scalar (and polynomial) or scripted/custom
    require_feature_alt(f_cv_active, f_cv_scalar, f_cv_linear, f_cv_scripted, f_cv_custom_function);

    init_feature(f_cv_awake, "awake", f_type_static);
    require_feature_self(f_cv_awake, f_cv_active);

    init_feature(f_cv_gradient, "gradient", f_type_dynamic);
    require_feature_children(f_cv_gradient, f_cvc_gradient);

    init_feature(f_cv_collect_gradient, "collect_gradient", f_type_dynamic);
    require_feature_self(f_cv_collect_gradient, f_cv_gradient);
    require_feature_self(f_cv_collect_gradient, f_cv_scalar);
    require_feature_self(f_cv_collect_gradient, f_cv_collect_atom_ids);
    // The following exclusions could be lifted by implementing the feature
    exclude_feature_self(f_cv_collect_gradient, f_cv_scripted);
    exclude_feature_self(f_cv_collect_gradient, f_cv_custom_function);
    require_feature_children(f_cv_collect_gradient, f_cvc_explicit_gradient);

    init_feature(f_cv_collect_atom_ids, "collect_atom_ids", f_type_dynamic);
    require_feature_children(f_cv_collect_atom_ids, f_cvc_collect_atom_ids);

    init_feature(f_cv_fdiff_velocity, "velocity_from_finite_differences", f_type_dynamic);

    // System force: either trivial (spring force); through extended Lagrangian, or calculated explicitly
    init_feature(f_cv_total_force, "total_force", f_type_dynamic);
    require_feature_alt(f_cv_total_force, f_cv_extended_Lagrangian, f_cv_total_force_calc);

    // Deps for explicit total force calculation
    init_feature(f_cv_total_force_calc, "total_force_calculation", f_type_dynamic);
    require_feature_self(f_cv_total_force_calc, f_cv_scalar);
    require_feature_self(f_cv_total_force_calc, f_cv_linear);
    require_feature_children(f_cv_total_force_calc, f_cvc_inv_gradient);
    require_feature_self(f_cv_total_force_calc, f_cv_Jacobian);

    init_feature(f_cv_Jacobian, "Jacobian_derivative", f_type_dynamic);
    require_feature_self(f_cv_Jacobian, f_cv_scalar);
    require_feature_self(f_cv_Jacobian, f_cv_linear);
    require_feature_children(f_cv_Jacobian, f_cvc_Jacobian);

    init_feature(f_cv_hide_Jacobian, "hide_Jacobian_force", f_type_user);
    require_feature_self(f_cv_hide_Jacobian, f_cv_Jacobian); // can only hide if calculated
    exclude_feature_self(f_cv_hide_Jacobian, f_cv_extended_Lagrangian);

    init_feature(f_cv_extended_Lagrangian, "extended_Lagrangian", f_type_user);
    require_feature_self(f_cv_extended_Lagrangian, f_cv_scalar);
    require_feature_self(f_cv_extended_Lagrangian, f_cv_gradient);

    init_feature(f_cv_Langevin, "Langevin_dynamics", f_type_user);
    require_feature_self(f_cv_Langevin, f_cv_extended_Lagrangian);

    init_feature(f_cv_external, "external", f_type_user);
    require_feature_self(f_cv_external, f_cv_single_cvc);

    init_feature(f_cv_single_cvc, "single_component", f_type_static);

    init_feature(f_cv_linear, "linear", f_type_static);

    init_feature(f_cv_scalar, "scalar", f_type_static);

    init_feature(f_cv_output_energy, "output_energy", f_type_user);

    init_feature(f_cv_output_value, "output_value", f_type_user);

    init_feature(f_cv_output_velocity, "output_velocity", f_type_user);
    require_feature_self(f_cv_output_velocity, f_cv_fdiff_velocity);

    init_feature(f_cv_output_applied_force, "output_applied_force", f_type_user);

    init_feature(f_cv_output_total_force, "output_total_force", f_type_user);
    require_feature_self(f_cv_output_total_force, f_cv_total_force);

    init_feature(f_cv_subtract_applied_force, "subtract_applied_force_from_total_force", f_type_user);
    require_feature_self(f_cv_subtract_applied_force, f_cv_total_force);

    init_feature(f_cv_lower_boundary, "lower_boundary", f_type_user);
    require_feature_self(f_cv_lower_boundary, f_cv_scalar);

    init_feature(f_cv_upper_boundary, "upper_boundary", f_type_user);
    require_feature_self(f_cv_upper_boundary, f_cv_scalar);

    init_feature(f_cv_hard_lower_boundary, "hard_lower_boundary", f_type_user);
    require_feature_self(f_cv_hard_lower_boundary, f_cv_lower_boundary);

    init_feature(f_cv_hard_upper_boundary, "hard_upper_boundary", f_type_user);
    require_feature_self(f_cv_hard_upper_boundary, f_cv_upper_boundary);

    init_feature(f_cv_reflecting_lower_boundary, "reflecting_lower_boundary", f_type_user);
    require_feature_self(f_cv_reflecting_lower_boundary, f_cv_lower_boundary);
    require_feature_self(f_cv_reflecting_lower_boundary, f_cv_extended_Lagrangian);

    init_feature(f_cv_reflecting_upper_boundary, "reflecting_upper_boundary", f_type_user);
    require_feature_self(f_cv_reflecting_upper_boundary, f_cv_upper_boundary);
    require_feature_self(f_cv_reflecting_upper_boundary, f_cv_extended_Lagrangian);

    init_feature(f_cv_grid, "grid", f_type_dynamic);
    require_feature_self(f_cv_grid, f_cv_scalar);

    init_feature(f_cv_runave, "running_average", f_type_user);

    init_feature(f_cv_corrfunc, "correlation_function", f_type_user);

    init_feature(f_cv_scripted, "scripted", f_type_user);

    init_feature(f_cv_custom_function, "custom_function", f_type_user);
    exclude_feature_self(f_cv_custom_function, f_cv_scripted);

    init_feature(f_cv_periodic, "periodic", f_type_static);
    require_feature_self(f_cv_periodic, f_cv_scalar);
    init_feature(f_cv_scalar, "scalar", f_type_static);
    init_feature(f_cv_linear, "linear", f_type_static);
    init_feature(f_cv_homogeneous, "homogeneous", f_type_static);

    // because total forces are obtained from the previous time step,
    // we cannot (currently) have colvar values and total forces for the same timestep
    init_feature(f_cv_multiple_ts, "multiple_timestep", f_type_static);
    exclude_feature_self(f_cv_multiple_ts, f_cv_total_force_calc);

    // check that everything is initialized
    for (i = 0; i < colvardeps::f_cv_ntot; i++) {
      if (is_not_set(i)) {
        cvm::error("Uninitialized feature " + cvm::to_str(i) + " in " + description);
      }
    }
  }

  // Initialize feature_states for each instance
  feature_states.reserve(f_cv_ntot);
  for (i = 0; i < f_cv_ntot; i++) {
    feature_states.push_back(feature_state(true, false));
    // Most features are available, so we set them so
    // and list exceptions below
   }

  feature_states[f_cv_fdiff_velocity].available =
    cvm::main()->proxy->simulation_running();

  return COLVARS_OK;
}


void colvar::setup()
{
  // loop over all components to update masses and charges of all groups
  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
      cvm::atom_group *atoms = cvcs[i]->atom_groups[ig];
      atoms->setup();
      atoms->print_properties(name, i, ig);
      atoms->read_positions();
    }
  }
}


std::vector<std::vector<int> > colvar::get_atom_lists()
{
  std::vector<std::vector<int> > lists;
  for (size_t i = 0; i < cvcs.size(); i++) {
    std::vector<std::vector<int> > li = cvcs[i]->get_atom_lists();
    lists.insert(lists.end(), li.begin(), li.end());
  }
  return lists;
}


std::vector<int> const &colvar::get_volmap_ids()
{
  volmap_ids_.resize(cvcs.size());
  for (size_t i = 0; i < cvcs.size(); i++) {
    if (cvcs[i]->param_exists("mapID") == COLVARS_OK) {
      volmap_ids_[i] =
        *(reinterpret_cast<int const *>(cvcs[i]->get_param_ptr("mapID")));
    } else {
      volmap_ids_[i] = -1;
    }
  }
  return volmap_ids_;
}


colvar::~colvar()
{
  // There is no need to call free_children_deps() here
  // because the children are cvcs and will be deleted
  // just below

  // Clear references to this colvar's cvcs as children
  // for dependency purposes
  remove_all_children();

  for (std::vector<cvc *>::reverse_iterator ci = cvcs.rbegin();
      ci != cvcs.rend();
      ++ci) {
    // clear all children of this cvc (i.e. its atom groups)
    // because the cvc base class destructor can't do it early enough
    // and we don't want to have each cvc derived class do it separately
    (*ci)->remove_all_children();
    delete *ci;
  }
  cvcs.clear();

  while (biases.size() > 0) {
    size_t const i = biases.size()-1;
    cvm::log("Warning: before deleting colvar " + name
             + ", deleting related bias " + biases[i]->name);
    delete biases[i];
  }
  biases.clear();

  // remove reference to this colvar from the module
  colvarmodule *cv = cvm::main();
  for (std::vector<colvar *>::iterator cvi = cv->variables()->begin();
       cvi != cv->variables()->end();
       ++cvi) {
    if ( *cvi == this) {
      cv->variables()->erase(cvi);
      break;
    }
  }

  cv->config_changed();

#ifdef LEPTON
  for (std::vector<Lepton::CompiledExpression *>::iterator cei = value_evaluators.begin();
       cei != value_evaluators.end();
       ++cei) {
    if (*cei != NULL) delete (*cei);
  }
  value_evaluators.clear();

  for (std::vector<Lepton::CompiledExpression *>::iterator gei = gradient_evaluators.begin();
       gei != gradient_evaluators.end();
       ++gei) {
    if (*gei != NULL) delete (*gei);
  }
  gradient_evaluators.clear();
#endif
}



// ******************** CALC FUNCTIONS ********************


// Default schedule (everything is serialized)
int colvar::calc()
{
  // Note: if anything is added here, it should be added also in the SMP block of calc_colvars()
  int error_code = COLVARS_OK;
  if (is_enabled(f_cv_active)) {
    error_code |= update_cvc_flags();
    if (error_code != COLVARS_OK) return error_code;
    error_code |= calc_cvcs();
    if (error_code != COLVARS_OK) return error_code;
    error_code |= collect_cvc_data();
  }
  return error_code;
}


int colvar::calc_cvcs(int first_cvc, size_t num_cvcs)
{
  if (cvm::debug())
    cvm::log("Calculating colvar \""+this->name+"\", components "+
             cvm::to_str(first_cvc)+" through "+cvm::to_str(first_cvc+num_cvcs)+".\n");

  colvarproxy *proxy = cvm::main()->proxy;
  int error_code = COLVARS_OK;

  error_code |= check_cvc_range(first_cvc, num_cvcs);
  if (error_code != COLVARS_OK) {
    return error_code;
  }

  if ((cvm::step_relative() > 0) && (!proxy->total_forces_same_step())){
    // Use Jacobian derivative from previous timestep
    error_code |= calc_cvc_total_force(first_cvc, num_cvcs);
  }
  // atom coordinates are updated by the next line
  error_code |= calc_cvc_values(first_cvc, num_cvcs);
  error_code |= calc_cvc_gradients(first_cvc, num_cvcs);
  error_code |= calc_cvc_Jacobians(first_cvc, num_cvcs);
  if (proxy->total_forces_same_step()){
    // Use Jacobian derivative from this timestep
    error_code |= calc_cvc_total_force(first_cvc, num_cvcs);
  }

  if (cvm::debug())
    cvm::log("Done calculating colvar \""+this->name+"\".\n");

  return error_code;
}


int colvar::collect_cvc_data()
{
  if (cvm::debug())
    cvm::log("Calculating colvar \""+this->name+"\"'s properties.\n");

  colvarproxy *proxy = cvm::main()->proxy;
  int error_code = COLVARS_OK;

  if ((cvm::step_relative() > 0) && (!proxy->total_forces_same_step())){
    // Total force depends on Jacobian derivative from previous timestep
    // collect_cvc_total_forces() uses the previous value of jd
    error_code |= collect_cvc_total_forces();
  }
  error_code |= collect_cvc_values();
  error_code |= collect_cvc_gradients();
  error_code |= collect_cvc_Jacobians();
  if (proxy->total_forces_same_step()){
    // Use Jacobian derivative from this timestep
    error_code |= collect_cvc_total_forces();
  }
  error_code |= calc_colvar_properties();

  if (cvm::debug())
    cvm::log("Done calculating colvar \""+this->name+"\"'s properties.\n");

  return error_code;
}


int colvar::check_cvc_range(int first_cvc, size_t /* num_cvcs */)
{
  if ((first_cvc < 0) || (first_cvc >= ((int) cvcs.size()))) {
    cvm::error("Error: trying to address a component outside the "
               "range defined for colvar \""+name+"\".\n", COLVARS_BUG_ERROR);
    return COLVARS_BUG_ERROR;
  }
  return COLVARS_OK;
}


int colvar::calc_cvc_values(int first_cvc, size_t num_cvcs)
{
  size_t const cvc_max_count = num_cvcs ? num_cvcs : num_active_cvcs();
  size_t i, cvc_count;

  // calculate the value of the colvar

  if (cvm::debug())
    cvm::log("Calculating colvar components.\n");

  // First, calculate component values
  cvm::increase_depth();
  for (i = first_cvc, cvc_count = 0;
       (i < cvcs.size()) && (cvc_count < cvc_max_count);
       i++) {
    if (!cvcs[i]->is_enabled()) continue;
    cvc_count++;
    (cvcs[i])->read_data();
    (cvcs[i])->calc_value();
    if (cvm::debug())
      cvm::log("Colvar component no. "+cvm::to_str(i+1)+
                " within colvar \""+this->name+"\" has value "+
                cvm::to_str((cvcs[i])->value(),
                cvm::cv_width, cvm::cv_prec)+".\n");
  }
  cvm::decrease_depth();

  return COLVARS_OK;
}


int colvar::collect_cvc_values()
{
  x.reset();

  // combine them appropriately, using either a scripted function or a polynomial
  if (is_enabled(f_cv_scripted)) {
    // cvcs combined by user script
    int res = cvm::proxy->run_colvar_callback(scripted_function, sorted_cvc_values, x);
    if (res == COLVARS_NOT_IMPLEMENTED) {
      cvm::error("Scripted colvars are not implemented.");
      return COLVARS_NOT_IMPLEMENTED;
    }
    if (res != COLVARS_OK) {
      cvm::error("Error running scripted colvar");
      return COLVARS_OK;
    }

#ifdef LEPTON
  } else if (is_enabled(f_cv_custom_function)) {

    size_t l = 0; // index in the vector of variable references

    for (size_t i = 0; i < x.size(); i++) {
      // Fill Lepton evaluator variables with CVC values, serialized into scalars
      for (size_t j = 0; j < cvcs.size(); j++) {
        for (size_t k = 0; k < cvcs[j]->value().size(); k++) {
          *(value_eval_var_refs[l++]) = cvcs[j]->value()[k];
        }
      }
      x[i] = value_evaluators[i]->evaluate();
    }
#endif

  } else if (x.type() == colvarvalue::type_scalar) {
    // polynomial combination allowed
    for (size_t i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      x += (cvcs[i])->sup_coeff *
      ( ((cvcs[i])->sup_np != 1) ?
        cvm::integer_power((cvcs[i])->value().real_value, (cvcs[i])->sup_np) :
        (cvcs[i])->value().real_value );
    }
  } else {
    for (size_t i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      x += (cvcs[i])->sup_coeff * (cvcs[i])->value();
    }
  }

  if (cvm::debug())
    cvm::log("Colvar \""+this->name+"\" has value "+
              cvm::to_str(x, cvm::cv_width, cvm::cv_prec)+".\n");

  if (after_restart) {
    if (cvm::proxy->simulation_running()) {
      cvm::real const jump2 = dist2(x, x_restart) / (width*width);
      if (jump2 > 0.25) {
        cvm::error("Error: the calculated value of colvar \""+name+
                   "\":\n"+cvm::to_str(x)+"\n differs greatly from the value "
                   "last read from the state file:\n"+cvm::to_str(x_restart)+
                   "\nPossible causes are changes in configuration, "
                   "wrong state file, or how PBC wrapping is handled.\n",
                   COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
    }
  }

  return COLVARS_OK;
}


int colvar::calc_cvc_gradients(int first_cvc, size_t num_cvcs)
{
  size_t const cvc_max_count = num_cvcs ? num_cvcs : num_active_cvcs();
  size_t i, cvc_count;

  if (cvm::debug())
    cvm::log("Calculating gradients of colvar \""+this->name+"\".\n");

  // calculate the gradients of each component
  cvm::increase_depth();
  for (i = first_cvc, cvc_count = 0;
      (i < cvcs.size()) && (cvc_count < cvc_max_count);
      i++) {
    if (!cvcs[i]->is_enabled()) continue;
    cvc_count++;

    if ((cvcs[i])->is_enabled(f_cvc_gradient)) {
      (cvcs[i])->calc_gradients();
      // if requested, propagate (via chain rule) the gradients above
      // to the atoms used to define the roto-translation
     (cvcs[i])->calc_fit_gradients();
      if ((cvcs[i])->is_enabled(f_cvc_debug_gradient))
        (cvcs[i])->debug_gradients();
    }

    cvm::decrease_depth();

    if (cvm::debug())
      cvm::log("Done calculating gradients of colvar \""+this->name+"\".\n");
  }

  return COLVARS_OK;
}


int colvar::collect_cvc_gradients()
{
  size_t i;
  if (is_enabled(f_cv_collect_gradient)) {
    // Collect the atomic gradients inside colvar object
    for (unsigned int a = 0; a < atomic_gradients.size(); a++) {
      atomic_gradients[a].reset();
    }
    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      cvcs[i]->collect_gradients(atom_ids, atomic_gradients);
    }
  }
  return COLVARS_OK;
}


int colvar::calc_cvc_total_force(int first_cvc, size_t num_cvcs)
{
  size_t const cvc_max_count = num_cvcs ? num_cvcs : num_active_cvcs();
  size_t i, cvc_count;

  if (is_enabled(f_cv_total_force_calc)) {
    if (cvm::debug())
      cvm::log("Calculating total force of colvar \""+this->name+"\".\n");

    cvm::increase_depth();

    for (i = first_cvc, cvc_count = 0;
        (i < cvcs.size()) && (cvc_count < cvc_max_count);
        i++) {
      if (!cvcs[i]->is_enabled()) continue;
      cvc_count++;
      (cvcs[i])->calc_force_invgrads();
    }
    cvm::decrease_depth();


    if (cvm::debug())
      cvm::log("Done calculating total force of colvar \""+this->name+"\".\n");
  }

  return COLVARS_OK;
}


int colvar::collect_cvc_total_forces()
{
  if (is_enabled(f_cv_total_force_calc)) {
    ft.reset();

    if (cvm::step_relative() > 0) {
      // get from the cvcs the total forces from the PREVIOUS step
      for (size_t i = 0; i < cvcs.size();  i++) {
        if (!cvcs[i]->is_enabled()) continue;
            if (cvm::debug())
            cvm::log("Colvar component no. "+cvm::to_str(i+1)+
                " within colvar \""+this->name+"\" has total force "+
                cvm::to_str((cvcs[i])->total_force(),
                cvm::cv_width, cvm::cv_prec)+".\n");
        // linear combination is assumed
        ft += (cvcs[i])->total_force() * (cvcs[i])->sup_coeff / active_cvc_square_norm;
      }
    }

    if (!(is_enabled(f_cv_hide_Jacobian) && is_enabled(f_cv_subtract_applied_force))) {
      // add the Jacobian force to the total force, and don't apply any silent
      // correction internally: biases such as colvarbias_abf will handle it
      // If f_cv_hide_Jacobian is enabled, a force of -fj is present in ft due to the
      // Jacobian-compensating force
      ft += fj;
    }
  }

  return COLVARS_OK;
}


int colvar::calc_cvc_Jacobians(int first_cvc, size_t num_cvcs)
{
  size_t const cvc_max_count = num_cvcs ? num_cvcs : num_active_cvcs();

  if (is_enabled(f_cv_Jacobian)) {
    cvm::increase_depth();
    size_t i, cvc_count;
    for (i = first_cvc, cvc_count = 0;
         (i < cvcs.size()) && (cvc_count < cvc_max_count);
         i++) {
      if (!cvcs[i]->is_enabled()) continue;
      cvc_count++;
      (cvcs[i])->calc_Jacobian_derivative();
    }
    cvm::decrease_depth();
  }

  return COLVARS_OK;
}


int colvar::collect_cvc_Jacobians()
{
  colvarproxy *proxy = cvm::main()->proxy;
  if (is_enabled(f_cv_Jacobian)) {
    fj.reset();
    for (size_t i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
        if (cvm::debug())
          cvm::log("Colvar component no. "+cvm::to_str(i+1)+
            " within colvar \""+this->name+"\" has Jacobian derivative"+
            cvm::to_str((cvcs[i])->Jacobian_derivative(),
            cvm::cv_width, cvm::cv_prec)+".\n");
      // linear combination is assumed
      fj += (cvcs[i])->Jacobian_derivative() * (cvcs[i])->sup_coeff / active_cvc_square_norm;
    }
    fj *= proxy->boltzmann() * proxy->target_temperature();
  }

  return COLVARS_OK;
}


int colvar::calc_colvar_properties()
{
  if (is_enabled(f_cv_fdiff_velocity)) {
    // calculate the velocity by finite differences
    if (cvm::step_relative() == 0) {
      x_old = x;
      v_fdiff.reset(); // Do not pretend we know anything about the actual velocity
      // eg. upon restarting. That would require saving v_fdiff or x_old to the state file
    } else {
      v_fdiff = fdiff_velocity(x_old, x);
      v_reported = v_fdiff;
    }
  }

  if (is_enabled(f_cv_extended_Lagrangian)) {
    // initialize the restraint center in the first step to the value
    // just calculated from the cvcs
    // Do the same if no simulation is running (eg. VMD postprocessing)
    if ((cvm::step_relative() == 0 && !after_restart) || x_ext.type() == colvarvalue::type_notset || !cvm::proxy->simulation_running()) {
      x_ext = x;
      if (is_enabled(f_cv_reflecting_lower_boundary) && x_ext < lower_boundary) {
        cvm::log("Warning: initializing extended coordinate to reflective lower boundary, as colvar value is below.");
        x_ext = lower_boundary;
      }
      if (is_enabled(f_cv_reflecting_upper_boundary) && x_ext > upper_boundary) {
        cvm::log("Warning: initializing extended coordinate to reflective upper boundary, as colvar value is above.");
        x_ext = upper_boundary;
      }

      v_ext.reset(); // (already 0; added for clarity)
    }

    // Special case of a repeated timestep (eg. multiple NAMD "run" statements)
    // revert values of the extended coordinate and velocity prior to latest integration
    if (cvm::proxy->simulation_running() && cvm::step_relative() == prev_timestep) {
      x_ext = prev_x_ext;
      v_ext = prev_v_ext;
    }
    // report the restraint center as "value"
    // These position and velocities come from integration at the _previous timestep_ in update_forces_energy()
    // But we report values at the beginning of the timestep (value at t=0 on the first timestep)
    x_reported = x_ext;
    v_reported = v_ext;
    // the "total force" with the extended Lagrangian is
    // calculated in update_forces_energy() below

  } else {

    if (is_enabled(f_cv_subtract_applied_force)) {
      // correct the total force only if it has been measured
      // TODO add a specific test instead of relying on sq norm
      if (ft.norm2() > 0.0) {
        ft -= f_old;
      }
    }

    x_reported = x;
    ft_reported = ft;
  }

  // At the end of the first update after a restart, we can reset the flag
  after_restart = false;
  return COLVARS_OK;
}


cvm::real colvar::update_forces_energy()
{
  if (cvm::debug())
    cvm::log("Updating colvar \""+this->name+"\".\n");

  // set to zero the applied force
  f.type(value());
  f.reset();
  fr.reset();

  // If we are not active at this timestep, that's all we have to do
  // return with energy == zero
  if (!is_enabled(f_cv_active)) return 0.;

  // add the biases' force, which at this point should already have
  // been summed over each bias using this colvar
  // fb is already multiplied by the relevant time step factor for each bias
  f += fb;

  if (is_enabled(f_cv_Jacobian)) {
    // the instantaneous Jacobian force was not included in the reported total force;
    // instead, it is subtracted from the applied force (silent Jacobian correction)
    // This requires the Jacobian term for the *current* timestep
    // Need to scale it for impulse MTS
    if (is_enabled(f_cv_hide_Jacobian))
      f -= fj * cvm::real(time_step_factor);
  }

  // At this point f is the force f from external biases that will be applied to the
  // extended variable if there is one
  if (is_enabled(f_cv_extended_Lagrangian) && cvm::proxy->simulation_running()) {
    update_extended_Lagrangian();
  }

  if (!is_enabled(f_cv_external)) {
    // Now adding the force on the actual colvar (for those biases that
    // bypass the extended Lagrangian mass)
    f += fb_actual;
  }

  if (cvm::debug())
    cvm::log("Done updating colvar \""+this->name+"\".\n");
  return (potential_energy + kinetic_energy);
}


void colvar::update_extended_Lagrangian()
{
  if (cvm::debug()) {
    cvm::log("Updating extended-Lagrangian degree of freedom.\n");
  }

  if (prev_timestep > -1L) {
    // Keep track of slow timestep to integrate MTS colvars
    // the colvar checks the interval after waking up twice
    cvm::step_number n_timesteps = cvm::step_relative() - prev_timestep;
    if (n_timesteps != 0 && n_timesteps != time_step_factor) {
      cvm::error("Error: extended-Lagrangian " + description + " has timeStepFactor " +
        cvm::to_str(time_step_factor) + ", but was activated after " + cvm::to_str(n_timesteps) +
        " steps at timestep " + cvm::to_str(cvm::step_absolute()) + " (relative step: " +
        cvm::to_str(cvm::step_relative()) + ").\n" +
        "Make sure that this colvar is requested by biases at multiples of timeStepFactor.\n");
      return;
    }
  }

  // Integrate with slow timestep (if time_step_factor != 1)
  cvm::real dt = cvm::dt() * cvm::real(time_step_factor);

  colvarvalue f_ext(fr.type()); // force acting on the extended variable
  f_ext.reset();

  if (is_enabled(f_cv_external)) {
    // There are no forces on the "actual colvar" bc there is no gradient wrt atomic coordinates
    // So we apply this to the extended DOF
    f += fb_actual;
  }

  fr    = f;
  // External force has been scaled for a 1-timestep impulse, scale it back because we will
  // integrate it with the colvar's own timestep factor
  f_ext = f / cvm::real(time_step_factor);

  colvarvalue f_system(fr.type()); // force exterted by the system on the extended DOF

  if (is_enabled(f_cv_external)) {
    // Add "alchemical" force from external variable
    f_system = cvcs[0]->total_force();
    // f is now irrelevant because we are not applying atomic forces in the simulation
    // just driving the external variable lambda
  } else {
    // the total force is applied to the fictitious mass, while the
    // atoms only feel the harmonic force + wall force
    // fr: bias force on extended variable (without harmonic spring), for output in trajectory
    // f_ext: total force on extended variable (including harmonic spring)
    // f: - initially, external biasing force
    //    - after this code block, colvar force to be applied to atomic coordinates
    //      ie. spring force (fb_actual will be added just below)
    f_system = (-0.5 * ext_force_k) * this->dist2_lgrad(x_ext, x);
    f        = -1.0 * f_system;
    // Coupling force is a slow force, to be applied to atomic coords impulse-style
    // over a single MD timestep
    f *= cvm::real(time_step_factor);
  }
  f_ext += f_system;

  if (is_enabled(f_cv_subtract_applied_force)) {
    // Report a "system" force without the biases on this colvar
    // that is, just the spring force (or alchemical force)
    ft_reported = f_system;
  } else {
    // The total force acting on the extended variable is f_ext
    // This will be used in the next timestep
    ft_reported = f_ext;
  }

  // backup in case we need to revert this integration timestep
  // if the same MD timestep is re-run
  prev_x_ext = x_ext;
  prev_v_ext = v_ext;

  // leapfrog: starting from x_i, f_i, v_(i-1/2)
  v_ext  += (0.5 * dt) * f_ext / ext_mass;
  // Because of leapfrog, kinetic energy at time i is approximate
  kinetic_energy = 0.5 * ext_mass * v_ext * v_ext;
  potential_energy = 0.5 * ext_force_k * this->dist2(x_ext, x);
  // leap to v_(i+1/2)
  if (is_enabled(f_cv_Langevin)) {
    v_ext -= dt * ext_gamma * v_ext;
    colvarvalue rnd(x);
    rnd.set_random();
    v_ext += dt * ext_sigma * rnd / ext_mass;
  }
  v_ext  += (0.5 * dt) * f_ext / ext_mass;
  x_ext  += dt * v_ext;

  cvm::real delta = 0; // Length of overshoot past either reflecting boundary
  if ((is_enabled(f_cv_reflecting_lower_boundary) && (delta = x_ext - lower_boundary) < 0) ||
      (is_enabled(f_cv_reflecting_upper_boundary) && (delta = x_ext - upper_boundary) > 0)) {
    x_ext -= 2.0 * delta;
    v_ext *= -1.0;
    if ((is_enabled(f_cv_reflecting_lower_boundary) && (delta = x_ext - lower_boundary) < 0) ||
        (is_enabled(f_cv_reflecting_upper_boundary) && (delta = x_ext - upper_boundary) > 0)) {
      cvm::error("Error: extended coordinate value " + cvm::to_str(x_ext) + " is still outside boundaries after reflection.\n");
    }
  }

  x_ext.apply_constraints();
  this->wrap(x_ext);
  if (is_enabled(f_cv_external)) {
    // Colvar value is constrained to the extended value
    x = x_ext;
    cvcs[0]->set_value(x_ext);
  }
}


int colvar::end_of_step()
{
  if (cvm::debug())
    cvm::log("End of step for colvar \""+this->name+"\".\n");

  if (is_enabled(f_cv_fdiff_velocity)) {
    x_old = x;
  }

  if (is_enabled(f_cv_subtract_applied_force)) {
    f_old = f;
  }

  prev_timestep = cvm::step_relative();

  return COLVARS_OK;
}


void colvar::communicate_forces()
{
  size_t i;
  if (cvm::debug()) {
    cvm::log("Communicating forces from colvar \""+this->name+"\".\n");
    cvm::log("Force to be applied: " + cvm::to_str(f) + "\n");
  }

  if (is_enabled(f_cv_scripted)) {
    std::vector<cvm::matrix2d<cvm::real> > func_grads;
    func_grads.reserve(cvcs.size());
    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      func_grads.push_back(cvm::matrix2d<cvm::real> (x.size(),
                                                     cvcs[i]->value().size()));
    }
    int res = cvm::proxy->run_colvar_gradient_callback(scripted_function, sorted_cvc_values, func_grads);

    if (res != COLVARS_OK) {
      if (res == COLVARS_NOT_IMPLEMENTED) {
        cvm::error("Colvar gradient scripts are not implemented.", COLVARS_NOT_IMPLEMENTED);
      } else {
        cvm::error("Error running colvar gradient script");
      }
      return;
    }

    int grad_index = 0; // index in the scripted gradients, to account for some components being disabled
    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      // cvc force is colvar force times colvar/cvc Jacobian
      // (vector-matrix product)
      (cvcs[i])->apply_force(colvarvalue(f.as_vector() * func_grads[grad_index++],
                             cvcs[i]->value().type()));
    }

#ifdef LEPTON
  } else if (is_enabled(f_cv_custom_function)) {

    size_t r = 0; // index in the vector of variable references
    size_t e = 0; // index of the gradient evaluator

    for (i = 0; i < cvcs.size(); i++) {  // gradient with respect to cvc i
      cvm::matrix2d<cvm::real> jacobian (x.size(), cvcs[i]->value().size());
      for (size_t j = 0; j < cvcs[i]->value().size(); j++) { // j-th element
        for (size_t c = 0; c < x.size(); c++) { // derivative of scalar element c of the colvarvalue

          // Feed cvc values to the evaluator
          for (size_t k = 0; k < cvcs.size(); k++) { //
            for (size_t l = 0; l < cvcs[k]->value().size(); l++) {
              *(grad_eval_var_refs[r++]) = cvcs[k]->value()[l];
            }
          }
          jacobian[c][j] = gradient_evaluators[e++]->evaluate();
        }
      }
      // cvc force is colvar force times colvar/cvc Jacobian
      // (vector-matrix product)
      (cvcs[i])->apply_force(colvarvalue(f.as_vector() * jacobian,
                             cvcs[i]->value().type()));
    }
#endif

  } else if (x.type() == colvarvalue::type_scalar) {

    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      (cvcs[i])->apply_force(f * (cvcs[i])->sup_coeff *
                             cvm::real((cvcs[i])->sup_np) *
                             (cvm::integer_power((cvcs[i])->value().real_value,
                                                 (cvcs[i])->sup_np-1)) );
    }

  } else {

    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      (cvcs[i])->apply_force(f * (cvcs[i])->sup_coeff);
    }
  }

  if (cvm::debug())
    cvm::log("Done communicating forces from colvar \""+this->name+"\".\n");
}


int colvar::set_cvc_flags(std::vector<bool> const &flags)
{
  if (flags.size() != cvcs.size()) {
    cvm::error("ERROR: Wrong number of CVC flags provided.");
    return COLVARS_ERROR;
  }
  // We cannot enable or disable cvcs in the middle of a timestep or colvar evaluation sequence
  // so we store the flags that will be enforced at the next call to calc()
  cvc_flags = flags;
  return COLVARS_OK;
}


void colvar::update_active_cvc_square_norm()
{
  active_cvc_square_norm = 0.0;
  for (size_t i = 0; i < cvcs.size(); i++) {
    if (cvcs[i]->is_enabled()) {
      active_cvc_square_norm += cvcs[i]->sup_coeff * cvcs[i]->sup_coeff;
    }
  }
}


int colvar::update_cvc_flags()
{
  // Update the enabled/disabled status of cvcs if necessary
  if (cvc_flags.size()) {
    n_active_cvcs = 0;
    for (size_t i = 0; i < cvcs.size(); i++) {
      cvcs[i]->set_enabled(f_cvc_active, cvc_flags[i]);
      if (cvcs[i]->is_enabled()) {
        n_active_cvcs++;
      }
    }
    if (!n_active_cvcs) {
      cvm::error("ERROR: All CVCs are disabled for colvar " + this->name +"\n");
      return COLVARS_ERROR;
    }
    cvc_flags.clear();

    update_active_cvc_square_norm();
  }

  return COLVARS_OK;
}


int colvar::update_cvc_config(std::vector<std::string> const &confs)
{
  cvm::log("Updating configuration for colvar \""+name+"\"\n");

  if (confs.size() != cvcs.size()) {
    return cvm::error("Error: Wrong number of CVC config strings.  "
                      "For those CVCs that are not being changed, try passing "
                      "an empty string.", COLVARS_INPUT_ERROR);
  }

  int error_code = COLVARS_OK;
  int num_changes = 0;
  for (size_t i = 0; i < cvcs.size(); i++) {
    if (confs[i].size()) {
      std::string conf(confs[i]);
      cvm::increase_depth();
      error_code |= cvcs[i]->colvar::cvc::init(conf);
      error_code |= cvcs[i]->check_keywords(conf,
                                            cvcs[i]->config_key.c_str());
      cvm::decrease_depth();
      num_changes++;
    }
  }

  if (num_changes == 0) {
    cvm::log("Warning: no changes were applied through modifycvcs; "
             "please check that its argument is a list of strings.\n");
  }

  update_active_cvc_square_norm();

  return error_code;
}


int colvar::cvc_param_exists(std::string const &param_name)
{
  if (is_enabled(f_cv_single_cvc)) {
    return cvcs[0]->param_exists(param_name);
  }
  return cvm::error("Error: calling colvar::cvc_param_exists() for a variable "
                    "with more than one component.\n", COLVARS_NOT_IMPLEMENTED);
}


cvm::real colvar::get_cvc_param(std::string const &param_name)
{
  if (is_enabled(f_cv_single_cvc)) {
    return cvcs[0]->get_param(param_name);
  }
  cvm::error("Error: calling colvar::get_cvc_param() for a variable "
             "with more than one component.\n", COLVARS_NOT_IMPLEMENTED);
  return 0.0;
}


void const *colvar::get_cvc_param_ptr(std::string const &param_name)
{
  if (is_enabled(f_cv_single_cvc)) {
    return cvcs[0]->get_param_ptr(param_name);
  }
  cvm::error("Error: calling colvar::get_cvc_param() for a variable "
             "with more than one component.\n", COLVARS_NOT_IMPLEMENTED);
  return NULL;
}


colvarvalue const *colvar::get_cvc_param_grad(std::string const &param_name)
{
  if (is_enabled(f_cv_single_cvc)) {
    return cvcs[0]->get_param_grad(param_name);
  }
  cvm::error("Error: calling colvar::get_cvc_param_grad() for a variable "
             "with more than one component.\n", COLVARS_NOT_IMPLEMENTED);
  return NULL;
}


int colvar::set_cvc_param(std::string const &param_name, void const *new_value)
{
  if (is_enabled(f_cv_single_cvc)) {
    return cvcs[0]->set_param(param_name, new_value);
  }
  return cvm::error("Error: calling colvar::set_cvc_param() for a variable "
                    "with more than one component.\n", COLVARS_NOT_IMPLEMENTED);
}


// ******************** METRIC FUNCTIONS ********************
// Use the metrics defined by \link colvar::cvc \endlink objects


bool colvar::periodic_boundaries(colvarvalue const &lb, colvarvalue const &ub) const
{
  if (period > 0.0) {
    if ( ((cvm::sqrt(this->dist2(lb, ub))) / this->width)
         < 1.0E-10 ) {
      return true;
    }
  }

  return false;
}

bool colvar::periodic_boundaries() const
{
  if ( (!is_enabled(f_cv_lower_boundary)) || (!is_enabled(f_cv_upper_boundary)) ) {
    // Return false if answer is unknown at this time
    return false;
  }

  return periodic_boundaries(lower_boundary, upper_boundary);
}


cvm::real colvar::dist2(colvarvalue const &x1,
                         colvarvalue const &x2) const
{
  if ( is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function) ) {
    if (is_enabled(f_cv_periodic) && is_enabled(f_cv_scalar)) {
      cvm::real diff = x1.real_value - x2.real_value;
      const cvm::real period_lower_boundary = wrap_center - period / 2.0;
      const cvm::real period_upper_boundary = wrap_center + period / 2.0;
      diff = (diff < period_lower_boundary ? diff + period : ( diff > period_upper_boundary ? diff - period : diff));
      return diff * diff;
    }
  }
  if (is_enabled(f_cv_homogeneous)) {
    return (cvcs[0])->dist2(x1, x2);
  } else {
    return x1.dist2(x2);
  }
}

colvarvalue colvar::dist2_lgrad(colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  if ( is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function) ) {
    if (is_enabled(f_cv_periodic) && is_enabled(f_cv_scalar)) {
      cvm::real diff = x1.real_value - x2.real_value;
      const cvm::real period_lower_boundary = wrap_center - period / 2.0;
      const cvm::real period_upper_boundary = wrap_center + period / 2.0;
      diff = (diff < period_lower_boundary ? diff + period : ( diff > period_upper_boundary ? diff - period : diff));
      return 2.0 * diff;
    }
  }
  if (is_enabled(f_cv_homogeneous)) {
    return (cvcs[0])->dist2_lgrad(x1, x2);
  } else {
    return x1.dist2_grad(x2);
  }
}

colvarvalue colvar::dist2_rgrad(colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  if ( is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function) ) {
    if (is_enabled(f_cv_periodic) && is_enabled(f_cv_scalar)) {
      cvm::real diff = x1.real_value - x2.real_value;
      const cvm::real period_lower_boundary = wrap_center - period / 2.0;
      const cvm::real period_upper_boundary = wrap_center + period / 2.0;
      diff = (diff < period_lower_boundary ? diff + period : ( diff > period_upper_boundary ? diff - period : diff));
      return (-2.0) * diff;
    }
  }
  if (is_enabled(f_cv_homogeneous)) {
    return (cvcs[0])->dist2_rgrad(x1, x2);
  } else {
    return x2.dist2_grad(x1);
  }
}


void colvar::wrap(colvarvalue &x_unwrapped) const
{
  if (!is_enabled(f_cv_periodic)) {
    return;
  }

  if ( is_enabled(f_cv_scripted) || is_enabled(f_cv_custom_function) ) {
    // Scripted functions do their own wrapping, as cvcs might not be periodic
    cvm::real shift = cvm::floor((x_unwrapped.real_value - wrap_center) /
                                 period + 0.5);
    x_unwrapped.real_value -= shift * period;
  } else {
    cvcs[0]->wrap(x_unwrapped);
  }
}


// ******************** INPUT FUNCTIONS ********************

std::istream & colvar::read_state(std::istream &is)
{
  std::streampos const start_pos = is.tellg();

  std::string conf;
  if ( !(is >> colvarparse::read_block("colvar", &conf)) ) {
    // this is not a colvar block
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  {
    std::string check_name = "";
    get_keyval(conf, "name", check_name,
               std::string(""), colvarparse::parse_silent);
    if (check_name.size() == 0) {
      cvm::error("Error: Collective variable in the "
                 "restart file without any identifier.\n", COLVARS_INPUT_ERROR);
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    }

    if (check_name != name)  {
      if (cvm::debug()) {
        cvm::log("Ignoring state of colvar \""+check_name+
                 "\": this colvar is named \""+name+"\".\n");
      }
      is.seekg(start_pos, std::ios::beg);
      return is;
    }
  }

  if ( !(get_keyval(conf, "x", x, x, colvarparse::parse_silent)) ) {
    cvm::log("Error: restart file does not contain "
             "the value of the colvar \""+
             name+"\" .\n");
  } else {
    cvm::log("Restarting collective variable \""+name+"\" from value: "+
             cvm::to_str(x)+"\n");
    x_restart = x;
    after_restart = true;
  }

  if (is_enabled(f_cv_extended_Lagrangian)) {
    if ( !(get_keyval(conf, "extended_x", x_ext,
                      colvarvalue(x.type()), colvarparse::parse_silent)) ||
         !(get_keyval(conf, "extended_v", v_ext,
                      colvarvalue(x.type()), colvarparse::parse_silent)) ) {
      cvm::log("Error: restart file does not contain "
               "\"extended_x\" or \"extended_v\" for the colvar \""+
               name+"\", but you requested \"extendedLagrangian\".\n");
    }
    x_reported = x_ext;
  } else {
    x_reported = x;
  }

  if (is_enabled(f_cv_output_velocity)) {

    if ( !(get_keyval(conf, "v", v_fdiff,
                      colvarvalue(x.type()), colvarparse::parse_silent)) ) {
      cvm::log("Error: restart file does not contain "
               "the velocity for the colvar \""+
               name+"\", but you requested \"outputVelocity\".\n");
    }

    if (is_enabled(f_cv_extended_Lagrangian)) {
      v_reported = v_ext;
    } else {
      v_reported = v_fdiff;
    }
  }

  return is;
}


std::istream & colvar::read_traj(std::istream &is)
{
  std::streampos const start_pos = is.tellg();

  if (is_enabled(f_cv_output_value)) {

    if (!(is >> x)) {
      cvm::log("Error: in reading the value of colvar \""+
                this->name+"\" from trajectory.\n");
      is.clear();
      is.seekg(start_pos, std::ios::beg);
      is.setstate(std::ios::failbit);
      return is;
    }

    if (is_enabled(f_cv_extended_Lagrangian)) {
      is >> x_ext;
      x_reported = x_ext;
    } else {
      x_reported = x;
    }
  }

  if (is_enabled(f_cv_output_velocity)) {

    is >> v_fdiff;

    if (is_enabled(f_cv_extended_Lagrangian)) {
      is >> v_ext;
      v_reported = v_ext;
    } else {
      v_reported = v_fdiff;
    }
  }

  if (is_enabled(f_cv_output_total_force)) {
    is >> ft;
    ft_reported = ft;
  }

  if (is_enabled(f_cv_output_applied_force)) {
    is >> f;
  }

  return is;
}


// ******************** OUTPUT FUNCTIONS ********************

std::ostream & colvar::write_state(std::ostream &os) {

  os << "colvar {\n"
     << "  name " << name << "\n"
     << "  x "
     << std::setprecision(cvm::cv_prec)
     << std::setw(cvm::cv_width)
     << x << "\n";

  if (is_enabled(f_cv_output_velocity)) {
    os << "  v "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << v_reported << "\n";
  }

  if (is_enabled(f_cv_extended_Lagrangian)) {
    os << "  extended_x "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << x_reported << "\n"
       << "  extended_v "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << v_reported << "\n";
  }

  os << "}\n\n";

  if (runave_outfile.size() > 0) {
    cvm::main()->proxy->flush_output_stream(runave_outfile);
  }

  return os;
}


std::ostream & colvar::write_traj_label(std::ostream & os)
{
  size_t const this_cv_width = x.output_width(cvm::cv_width);

  os << " ";

  if (is_enabled(f_cv_output_value)) {

    os << " "
       << cvm::wrap_string(this->name, this_cv_width);

    if (is_enabled(f_cv_extended_Lagrangian) && !is_enabled(f_cv_external)) {
      // extended DOF
      os << " r_"
         << cvm::wrap_string(this->name, this_cv_width-2);
    }
  }

  if (is_enabled(f_cv_output_velocity)) {

    os << " v_"
       << cvm::wrap_string(this->name, this_cv_width-2);

    if (is_enabled(f_cv_extended_Lagrangian) && !is_enabled(f_cv_external)) {
      // extended DOF
      os << " vr_"
         << cvm::wrap_string(this->name, this_cv_width-3);
    }
  }

  if (is_enabled(f_cv_output_energy)) {
    os << " Ep_"
       << cvm::wrap_string(this->name, this_cv_width-3)
       << " Ek_"
       << cvm::wrap_string(this->name, this_cv_width-3);
  }

  if (is_enabled(f_cv_output_total_force)) {
    os << " ft_"
       << cvm::wrap_string(this->name, this_cv_width-3);
  }

  if (is_enabled(f_cv_output_applied_force)) {
    os << " fa_"
       << cvm::wrap_string(this->name, this_cv_width-3);
  }

  return os;
}


std::ostream & colvar::write_traj(std::ostream &os)
{
  os << " ";
  if (is_enabled(f_cv_output_value)) {

    if (is_enabled(f_cv_extended_Lagrangian) && !is_enabled(f_cv_external)) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << x;
    }

    os << " "
       << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
       << x_reported;
  }

  if (is_enabled(f_cv_output_velocity)) {

    if (is_enabled(f_cv_extended_Lagrangian) && !is_enabled(f_cv_external)) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << v_fdiff;
    }

    os << " "
       << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
       << v_reported;
  }

  if (is_enabled(f_cv_output_energy)) {
    os << " "
       << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
       << potential_energy
       << " "
       << kinetic_energy;
  }


  if (is_enabled(f_cv_output_total_force)) {
    os << " "
       << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
       << ft_reported;
  }

  if (is_enabled(f_cv_output_applied_force)) {
    os << " "
       << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
       << applied_force();
  }

  return os;
}


int colvar::write_output_files()
{
  int error_code = COLVARS_OK;

  if (is_enabled(f_cv_corrfunc)) {
    if (acf.size()) {
      if (acf_outfile.size() == 0) {
        acf_outfile = std::string(cvm::output_prefix()+"."+this->name+
                                  ".corrfunc.dat");
      }
      cvm::log("Writing correlation function to file \""+acf_outfile+"\".\n");
      cvm::backup_file(acf_outfile.c_str());
      std::ostream &acf_os = cvm::proxy->output_stream(acf_outfile,
                                                       "colvar ACF file");
      if (!acf_os) {
        error_code |= COLVARS_FILE_ERROR;
      } else {
        error_code |= write_acf(acf_os);
      }
      cvm::proxy->close_output_stream(acf_outfile);
    }
  }

  return error_code;
}



// ******************** ANALYSIS FUNCTIONS ********************

int colvar::analyze()
{
  int error_code = COLVARS_OK;

  if (is_enabled(f_cv_runave)) {
    error_code |= calc_runave();
  }

  if (is_enabled(f_cv_corrfunc)) {
    error_code |= calc_acf();
  }

  return error_code;
}


inline void history_add_value(size_t const           &history_length,
                              std::list<colvarvalue> &history,
                              colvarvalue const      &new_value)
{
  history.push_front(new_value);
  if (history.size() > history_length)
    history.pop_back();
}


inline void history_incr(std::list< std::list<colvarvalue> >           &history,
                         std::list< std::list<colvarvalue> >::iterator &history_p)
{
  if ((++history_p) == history.end())
    history_p = history.begin();
}


int colvar::calc_acf()
{
  // using here an acf_stride-long list of vectors for either
  // coordinates (acf_x_history) or velocities (acf_v_history); each vector can
  // contain up to acf_length values, which are contiguous in memory
  // representation but separated by acf_stride in the time series;
  // the pointer to each vector is changed at every step

  colvar const *cfcv = cvm::colvar_by_name(acf_colvar_name);
  if (cfcv == NULL) {
    return cvm::error("Error: collective variable \""+acf_colvar_name+
                      "\" is not defined at this time.\n", COLVARS_INPUT_ERROR);
  }

  if (acf_x_history.empty() && acf_v_history.empty()) {

    // first-step operations

    if (colvarvalue::check_types(cfcv->value(), value())) {
      cvm::error("Error: correlation function between \""+cfcv->name+
                 "\" and \""+this->name+"\" cannot be calculated, "
                 "because their value types are different.\n",
                 COLVARS_INPUT_ERROR);
    }
    acf_nframes = 0;

    cvm::log("Colvar \""+this->name+"\": initializing correlation function "
             "calculation.\n");

    if (acf.size() < acf_length+1)
      acf.resize(acf_length+1, 0.0);

    size_t i;
    switch (acf_type) {

    case acf_vel:
      // allocate space for the velocities history
      for (i = 0; i < acf_stride; i++) {
        acf_v_history.push_back(std::list<colvarvalue>());
      }
      acf_v_history_p = acf_v_history.begin();
      break;

    case acf_coor:
    case acf_p2coor:
      // allocate space for the coordinates history
      for (i = 0; i < acf_stride; i++) {
        acf_x_history.push_back(std::list<colvarvalue>());
      }
      acf_x_history_p = acf_x_history.begin();
      break;

    case acf_notset:
    default:
      break;
    }

  } else if (cvm::step_relative() > prev_timestep) {

    switch (acf_type) {

    case acf_vel:

      calc_vel_acf((*acf_v_history_p), cfcv->velocity());
      history_add_value(acf_length+acf_offset, *acf_v_history_p,
                        cfcv->velocity());
      history_incr(acf_v_history, acf_v_history_p);
      break;

    case acf_coor:

      calc_coor_acf((*acf_x_history_p), cfcv->value());
      history_add_value(acf_length+acf_offset, *acf_x_history_p,
                        cfcv->value());
      history_incr(acf_x_history, acf_x_history_p);
      break;

    case acf_p2coor:

      calc_p2coor_acf((*acf_x_history_p), cfcv->value());
      history_add_value(acf_length+acf_offset, *acf_x_history_p,
                        cfcv->value());
      history_incr(acf_x_history, acf_x_history_p);
      break;

    case acf_notset:
    default:
      break;
    }
  }

  return COLVARS_OK;
}


void colvar::calc_vel_acf(std::list<colvarvalue> &v_list,
                          colvarvalue const      &v)
{
  // loop over stored velocities and add to the ACF, but only if the
  // length is sufficient to hold an entire row of ACF values
  if (v_list.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  vs_i = v_list.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      ++vs_i;

    // current vel with itself
    *(acf_i) += v.norm2();
    ++acf_i;

    // inner products of previous velocities with current (acf_i and
    // vs_i are updated)
    colvarvalue::inner_opt(v, vs_i, v_list.end(), acf_i);

    acf_nframes++;
  }
}


void colvar::calc_coor_acf(std::list<colvarvalue> &x_list,
                           colvarvalue const &x_now)
{
  // same as above but for coordinates
  if (x_list.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  xs_i = x_list.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      ++xs_i;

    *(acf_i++) += x.norm2();

    colvarvalue::inner_opt(x_now, xs_i, x_list.end(), acf_i);

    acf_nframes++;
  }
}


void colvar::calc_p2coor_acf(std::list<colvarvalue> &x_list,
                             colvarvalue const &x_now)
{
  // same as above but with second order Legendre polynomial instead
  // of just the scalar product
  if (x_list.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  xs_i = x_list.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      ++xs_i;

    // value of P2(0) = 1
    *(acf_i++) += 1.0;

    colvarvalue::p2leg_opt(x_now, xs_i, x_list.end(), acf_i);

    acf_nframes++;
  }
}


int colvar::write_acf(std::ostream &os)
{
  if (!acf_nframes) {
    return COLVARS_OK;
  }

  os.setf(std::ios::scientific, std::ios::floatfield);
  os << "# ";
  switch (acf_type) {
  case acf_vel:
    os << "Velocity";
    break;
  case acf_coor:
    os << "Coordinate";
    break;
  case acf_p2coor:
    os << "Coordinate (2nd Legendre poly)";
    break;
  case acf_notset:
  default:
    break;
  }

  if (acf_colvar_name == name) {
    os << " autocorrelation function for variable \""
       << this->name << "\"\n";
  } else {
    os << " correlation function between variables \"" //
       << this->name << "\" and \"" << acf_colvar_name << "\"\n";
  }

  os << "# Number of samples = ";
  if (acf_normalize) {
    os << (acf_nframes-1) << " (one DoF is used for normalization)\n";
  } else {
    os << acf_nframes << "\n";
  }

  os << "# " << cvm::wrap_string("step", cvm::it_width-2) << " "
     << cvm::wrap_string("corrfunc(step)", cvm::cv_width) << "\n";

  cvm::real const acf_norm = acf.front() / cvm::real(acf_nframes);

  std::vector<cvm::real>::iterator acf_i;
  size_t it = acf_offset;
  for (acf_i = acf.begin(); acf_i != acf.end(); ++acf_i) {
    os << std::setw(cvm::it_width) << acf_stride * (it++) << " "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << ( acf_normalize ?
            (*acf_i)/(acf_norm * cvm::real(acf_nframes)) :
            (*acf_i)/(cvm::real(acf_nframes)) ) << "\n";
  }

  return os.good() ? COLVARS_OK : COLVARS_FILE_ERROR;
}


int colvar::calc_runave()
{
  int error_code = COLVARS_OK;
  colvarproxy *proxy = cvm::main()->proxy;

  if (x_history.empty()) {

    runave.type(value().type());
    runave.reset();

    // first-step operationsf

    if (cvm::debug())
      cvm::log("Colvar \""+this->name+
                "\": initializing running average calculation.\n");

    acf_nframes = 0;

    x_history.push_back(std::list<colvarvalue>());
    x_history_p = x_history.begin();

  } else {

    if ( (cvm::step_relative() % runave_stride) == 0 &&
         (cvm::step_relative() > prev_timestep) ) {

      if ((*x_history_p).size() >= runave_length-1) {

        if (runave_outfile.size() == 0) {
          runave_outfile = std::string(cvm::output_prefix()+"."+
                                       this->name+".runave.traj");
        }

        if (! proxy->output_stream_exists(runave_outfile)) {
          size_t const this_cv_width = x.output_width(cvm::cv_width);
          std::ostream &runave_os = proxy->output_stream(runave_outfile,
                                                         "colvar running average");
          runave_os.setf(std::ios::scientific, std::ios::floatfield);
          runave_os << "# " << cvm::wrap_string("step", cvm::it_width-2)
                    << "   "
                    << cvm::wrap_string("running average", this_cv_width)
                    << " "
                    << cvm::wrap_string("running stddev", this_cv_width)
                    << "\n";
        }

        runave = x;
        std::list<colvarvalue>::iterator xs_i;
        for (xs_i = (*x_history_p).begin();
             xs_i != (*x_history_p).end(); ++xs_i) {
          runave += (*xs_i);
        }
        runave *= 1.0 / cvm::real(runave_length);
        runave.apply_constraints();

        runave_variance = 0.0;
        runave_variance += this->dist2(x, runave);
        for (xs_i = (*x_history_p).begin();
             xs_i != (*x_history_p).end(); ++xs_i) {
          runave_variance += this->dist2(x, (*xs_i));
        }
        runave_variance *= 1.0 / cvm::real(runave_length-1);

        if (runave_outfile.size() > 0) {
          std::ostream &runave_os = proxy->output_stream(runave_outfile);
          runave_os << std::setw(cvm::it_width) << cvm::step_relative()
                    << "   "
                    << std::setprecision(cvm::cv_prec)
                    << std::setw(cvm::cv_width)
                    << runave << " "
                    << std::setprecision(cvm::cv_prec)
                    << std::setw(cvm::cv_width)
                    << cvm::sqrt(runave_variance) << "\n";
        }
      }

      history_add_value(runave_length, *x_history_p, x);
    }
  }

  return error_code;
}

// Static members

std::vector<colvardeps::feature *> colvar::cv_features;

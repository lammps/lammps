/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvarscript.h"
#include <algorithm>



colvar::colvar (std::string const &conf)
{
  size_t i, j;
  cvm::log ("Initializing a new collective variable.\n");

  get_keyval (conf, "name", this->name,
              (std::string ("colvar")+cvm::to_str (cvm::colvars.size()+1)));

  if (cvm::colvar_by_name (this->name) != NULL) {
     cvm::error ("Error: this colvar cannot have the same name, \""+this->name+
                      "\", as another colvar.\n",
               INPUT_ERROR);
  }

  // all tasks disabled by default
  for (i = 0; i < task_ntot; i++) {
    tasks[i] = false;
  }

  kinetic_energy = 0.0;
  potential_energy = 0.0;

  // read the configuration and set up corresponding instances, for
  // each type of component implemented
#define initialize_components(def_desc,def_config_key,def_class_name)   \
  {                                                                     \
    size_t def_count = 0;                                               \
    std::string def_conf = "";                                          \
    size_t pos = 0;                                                     \
    while ( this->key_lookup (conf,                                     \
                              def_config_key,                           \
                              def_conf,                                 \
                              pos) ) {                                  \
      if (!def_conf.size()) continue;                                   \
      cvm::log ("Initializing "                                         \
                "a new \""+std::string (def_config_key)+"\" component"+ \
                (cvm::debug() ? ", with configuration:\n"+def_conf      \
                 : ".\n"));                                             \
      cvm::increase_depth();                                            \
      cvc *cvcp = new colvar::def_class_name (def_conf);                \
      if (cvcp != NULL) {                                               \
        cvcs.push_back (cvcp);                                          \
        cvcp->check_keywords (def_conf, def_config_key);                \
        cvm::decrease_depth();                                          \
      } else {                                                          \
        cvm::error ("Error: in allocating component \""                 \
                          def_config_key"\".\n",                        \
                          MEMORY_ERROR);                                \
      }                                                                 \
      if ( (cvcp->period != 0.0) || (cvcp->wrap_center != 0.0) ) {      \
        if ( (cvcp->function_type != std::string ("distance_z")) &&     \
             (cvcp->function_type != std::string ("dihedral")) &&       \
             (cvcp->function_type != std::string ("spin_angle")) ) {    \
          cvm::error ("Error: invalid use of period and/or "            \
                            "wrapAround in a \""+                       \
                            std::string (def_config_key)+               \
                            "\" component.\n"+                          \
                            "Period: "+cvm::to_str(cvcp->period) +      \
                        " wrapAround: "+cvm::to_str(cvcp->wrap_center), \
                        INPUT_ERROR);                                   \
        }                                                               \
      }                                                                 \
      if ( ! cvcs.back()->name.size())                                  \
        cvcs.back()->name = std::string (def_config_key)+               \
          (cvm::to_str (++def_count));                                  \
      if (cvm::debug())                                                 \
        cvm::log ("Done initializing a \""+                             \
                  std::string (def_config_key)+                         \
                  "\" component"+                                       \
                  (cvm::debug() ?                                       \
                   ", named \""+cvcs.back()->name+"\""                  \
                   : "")+".\n");                                        \
      def_conf = "";                                                    \
      if (cvm::debug())                                                 \
        cvm::log ("Parsed "+cvm::to_str (cvcs.size())+                  \
                  " components at this time.\n");                       \
    }                                                                   \
  }


  initialize_components ("distance",         "distance",       distance);
  initialize_components ("distance vector",  "distanceVec",    distance_vec);
  initialize_components ("distance vector "
                         "direction",        "distanceDir",    distance_dir);
  initialize_components ("distance projection "
                         "on an axis",       "distanceZ",      distance_z);
  initialize_components ("distance projection "
                         "on a plane",       "distanceXY",     distance_xy);
  initialize_components ("average distance weighted by inverse power",
                         "distanceInv", distance_inv);

  initialize_components ("coordination "
                         "number",           "coordNum",       coordnum);
  initialize_components ("self-coordination "
                         "number",           "selfCoordNum",   selfcoordnum);

  initialize_components ("angle",            "angle",          angle);
  initialize_components ("dihedral",         "dihedral",       dihedral);

  initialize_components ("hydrogen bond",    "hBond",          h_bond);

  //  initialize_components ("alpha helix",      "alphaDihedrals", alpha_dihedrals);
  initialize_components ("alpha helix",      "alpha",          alpha_angles);

  initialize_components ("dihedral principal "
                         "component",        "dihedralPC",     dihedPC);

  initialize_components ("orientation",      "orientation",    orientation);
  initialize_components ("orientation "
                         "angle",            "orientationAngle",orientation_angle);
  initialize_components ("orientation "
                         "projection",       "orientationProj",orientation_proj);
  initialize_components ("tilt",             "tilt",           tilt);
  initialize_components ("spin angle",       "spinAngle",      spin_angle);

  initialize_components ("RMSD",             "rmsd",           rmsd);

  //  initialize_components ("logarithm of MSD", "logmsd",         logmsd);

  initialize_components ("radius of "
                         "gyration",         "gyration",       gyration);
  initialize_components ("moment of "
                         "inertia",          "inertia",        inertia);
  initialize_components ("moment of inertia around an axis",
                                             "inertiaZ",       inertia_z);
  initialize_components ("eigenvector",      "eigenvector",    eigenvector);

  if (!cvcs.size()) {
    cvm::error ("Error: no valid components were provided "
                      "for this collective variable.\n",
              INPUT_ERROR);
  }

  cvm::log ("All components initialized.\n");

  // Setup colvar as scripted function of components
  if (get_keyval (conf, "scriptedFunction", scripted_function,
    "", colvarparse::parse_silent)) {

    enable(task_scripted);
    cvm::log("This colvar is a scripted function.");

    std::string type_str;
    get_keyval (conf, "scriptedFunctionType", type_str, "scalar");

    x.type(colvarvalue::type_notset);
    for (i = 0; i < colvarvalue::type_all; i++) {
      if (type_str == colvarvalue::type_keyword[i]) {
        x.type(colvarvalue::Type(i));
        break;
      }
    }
    if (x.type() == colvarvalue::type_notset) {
      cvm::error("Could not parse scripted colvar type.");
      return;
    }
    x_reported.type (x.type());
    cvm::log(std::string("Expecting colvar value of type ")
      + colvarvalue::type_desc[x.type()]);

    // Build ordered list of component values that will be
    // passed to the script
    for (i = 1; i <= cvcs.size(); i++) {
      for (j = 0; j < cvcs.size(); j++) {
        if (cvcs[j]->sup_np == int(i)) {
          sorted_cvc_values.push_back(cvcs[j]->p_value());
          break;
        }
      }
    }
    if (sorted_cvc_values.size() != cvcs.size()) {
      cvm::error("Could not find order numbers for all components"
                  "in componentExp values.");
      return;
    }
  }

  // this is set false if any of the components has an exponent
  // different from 1 in the polynomial
  b_linear = true;

  // these will be set to false if any of the cvcs has them false
  b_inverse_gradients = true;
  b_Jacobian_force    = true;

  // Decide whether the colvar is periodic
  // Used to wrap extended DOF if extendedLagrangian is on
  if (cvcs.size() == 1 && (cvcs[0])->b_periodic && (cvcs[0])->sup_np == 1
                                                && (cvcs[0])->sup_coeff == 1.0 ) {
    this->b_periodic = true;
    this->period = (cvcs[0])->period;
    // TODO write explicit wrap() function for colvars to allow for
    // sup_coeff different from 1
    // this->period = (cvcs[0])->period * (cvcs[0])->sup_coeff;
  } else {
    this->b_periodic = false;
    this->period = 0.0;
  }

  // check the available features of each cvc
  for (i = 0; i < cvcs.size(); i++) {

    if ((cvcs[i])->b_debug_gradients)
      enable (task_gradients);

    if ((cvcs[i])->sup_np != 1) {
      if (cvm::debug() && b_linear)
        cvm::log ("Warning: You are using a non-linear polynomial "
                  "combination to define this collective variable, "
                  "some biasing methods may be unavailable.\n");
      b_linear = false;

      if ((cvcs[i])->sup_np < 0) {
        cvm::log ("Warning: you chose a negative exponent in the combination; "
                  "if you apply forces, the simulation may become unstable "
                  "when the component \""+
                  (cvcs[i])->function_type+"\" approaches zero.\n");
      }
    }

    if ((cvcs[i])->b_periodic && !b_periodic) {
        cvm::log ("Warning: although this component is periodic, the colvar will "
                  "not be treated as periodic, either because the exponent is not "
                  "1, or because multiple components are present. Make sure that "
                  "you know what you are doing!");
    }

    if (! (cvcs[i])->b_inverse_gradients)
      b_inverse_gradients = false;

    if (! (cvcs[i])->b_Jacobian_derivative)
      b_Jacobian_force = false;

    if (!tasks[task_scripted]) {
      // If the combination of components is a scripted function,
      // the components may have different types
      for (size_t j = i; j < cvcs.size(); j++) {
        if ( (cvcs[i])->type() != (cvcs[j])->type() ) {
          cvm::log ("ERROR: you are definining this collective variable "
                            "by using components of different types, \""+
                            colvarvalue::type_desc[(cvcs[i])->type()]+
                            "\" and \""+
                            colvarvalue::type_desc[(cvcs[j])->type()]+
                            "\". "
                            "You must use the same type in order to "
                            " sum them together.\n");
          cvm::set_error_bits(INPUT_ERROR);
        }
      }
    }
  }

  if (!tasks[task_scripted]) {
    colvarvalue::Type const value_type = (cvcs[0])->type();
    if (cvm::debug())
      cvm::log ("This collective variable is a "+
                colvarvalue::type_desc[value_type]+", corresponding to "+
                cvm::to_str (colvarvalue::dof_num[value_type])+
                " internal degrees of freedom.\n");
    x.type (value_type);
    x_reported.type (value_type);
  }

  get_keyval (conf, "width", width, 1.0);
  if (width <= 0.0) {
    cvm::error("Error: \"width\" must be positive.\n", INPUT_ERROR);
  }

  lower_boundary.type (this->type());
  lower_wall.type     (this->type());

  upper_boundary.type (this->type());
  upper_wall.type     (this->type());

  if (this->type() == colvarvalue::type_scalar) {

    if (get_keyval (conf, "lowerBoundary", lower_boundary, lower_boundary)) {
      enable (task_lower_boundary);
    }

    get_keyval (conf, "lowerWallConstant", lower_wall_k, 0.0);
    if (lower_wall_k > 0.0) {
      get_keyval (conf, "lowerWall", lower_wall, lower_boundary);
      enable (task_lower_wall);
    }

    if (get_keyval (conf, "upperBoundary", upper_boundary, upper_boundary)) {
      enable (task_upper_boundary);
    }

    get_keyval (conf, "upperWallConstant", upper_wall_k, 0.0);
    if (upper_wall_k > 0.0) {
      get_keyval (conf, "upperWall", upper_wall, upper_boundary);
      enable (task_upper_wall);
    }
  }

  if (tasks[task_lower_boundary]) {
    get_keyval (conf, "hardLowerBoundary", hard_lower_boundary, false);
  }
  if (tasks[task_upper_boundary]) {
    get_keyval (conf, "hardUpperBoundary", hard_upper_boundary, false);
  }

  // consistency checks for boundaries and walls
  if (tasks[task_lower_boundary] && tasks[task_upper_boundary]) {
    if (lower_boundary >= upper_boundary) {
      cvm::error ("Error: the upper boundary, "+
                        cvm::to_str (upper_boundary)+
                        ", is not higher than the lower boundary, "+
                        cvm::to_str (lower_boundary)+".\n",
                INPUT_ERROR);
    }
  }

  if (tasks[task_lower_wall] && tasks[task_upper_wall]) {
    if (lower_wall >= upper_wall) {
      cvm::error ("Error: the upper wall, "+
                        cvm::to_str (upper_wall)+
                        ", is not higher than the lower wall, "+
                        cvm::to_str (lower_wall)+".\n",
                INPUT_ERROR);
    }

    if (dist2 (lower_wall, upper_wall) < 1.0E-12) {
      cvm::log ("Lower wall and upper wall are equal "
                "in the periodic domain of the colvar: disabling walls.\n");
      disable (task_lower_wall);
      disable (task_upper_wall);
    }
  }

  get_keyval (conf, "expandBoundaries", expand_boundaries, false);
  if (expand_boundaries && periodic_boundaries()) {
    cvm::error ("Error: trying to expand boundaries that already "
                      "cover a whole period of a periodic colvar.\n",
              INPUT_ERROR);
  }
  if (expand_boundaries && hard_lower_boundary && hard_upper_boundary) {
    cvm::error ("Error: inconsistent configuration "
                      "(trying to expand boundaries with both "
                      "hardLowerBoundary and hardUpperBoundary enabled).\n",
              INPUT_ERROR);
}

  {
    bool b_extended_lagrangian;
    get_keyval (conf, "extendedLagrangian", b_extended_lagrangian, false);

    if (b_extended_lagrangian) {
      cvm::real temp, tolerance, period;

      cvm::log ("Enabling the extended Lagrangian term for colvar \""+
                this->name+"\".\n");

      enable (task_extended_lagrangian);

      xr.type (this->type());
      vr.type (this->type());
      fr.type (this->type());

      const bool found = get_keyval (conf, "extendedTemp", temp, cvm::temperature());
      if (temp <= 0.0) {
        if (found)
          cvm::log ("Error: \"extendedTemp\" must be positive.\n");
        else
          cvm::error ("Error: a positive temperature must be provided, either "
                            "by enabling a thermostat, or through \"extendedTemp\".\n",
                      INPUT_ERROR);
      }

      get_keyval (conf, "extendedFluctuation", tolerance);
      if (tolerance <= 0.0) {
        cvm::error("Error: \"extendedFluctuation\" must be positive.\n", INPUT_ERROR);
      }
      ext_force_k = cvm::boltzmann() * temp / (tolerance * tolerance);
      cvm::log ("Computed extended system force constant: " + cvm::to_str(ext_force_k) + " kcal/mol/U^2");

      get_keyval (conf, "extendedTimeConstant", period, 200.0);
      if (period <= 0.0) {
        cvm::error("Error: \"extendedTimeConstant\" must be positive.\n", INPUT_ERROR);
      }
      ext_mass = (cvm::boltzmann() * temp * period * period)
                 / (4.0 * PI * PI * tolerance * tolerance);
      cvm::log ("Computed fictitious mass: " + cvm::to_str(ext_mass) + " kcal/mol/(U/fs)^2   (U: colvar unit)");

      {
        bool b_output_energy;
        get_keyval (conf, "outputEnergy", b_output_energy, false);
        if (b_output_energy) {
          enable (task_output_energy);
        }
      }

      get_keyval (conf, "extendedLangevinDamping", ext_gamma, 1.0);
      if (ext_gamma < 0.0) {
        cvm::error("Error: \"extendedLangevinDamping\" may not be negative.\n", INPUT_ERROR);
      }
      if (ext_gamma != 0.0) {
        enable (task_langevin);
        ext_gamma *= 1.0e-3; // convert from ps-1 to fs-1
        ext_sigma = std::sqrt(2.0 * cvm::boltzmann() * temp * ext_gamma * ext_mass / cvm::dt());
      }
    }
  }

  {
    bool b_output_value;
    get_keyval (conf, "outputValue", b_output_value, true);
    if (b_output_value) {
      enable (task_output_value);
    }
  }

  {
    bool b_output_velocity;
    get_keyval (conf, "outputVelocity", b_output_velocity, false);
    if (b_output_velocity) {
      enable (task_output_velocity);
    }
  }

  {
    bool b_output_system_force;
    get_keyval (conf, "outputSystemForce", b_output_system_force, false);
    if (b_output_system_force) {
      enable (task_output_system_force);
    }
  }

  {
    bool b_output_applied_force;
    get_keyval (conf, "outputAppliedForce", b_output_applied_force, false);
    if (b_output_applied_force) {
      enable (task_output_applied_force);
    }
  }

  if (cvm::b_analysis)
    parse_analysis (conf);

  if (cvm::debug())
    cvm::log ("Done initializing collective variable \""+this->name+"\".\n");
}



void colvar::build_atom_list (void)
{
  // If atomic gradients are requested, build full list of atom ids from all cvcs
  std::list<int> temp_id_list;

  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t j = 0; j < cvcs[i]->atom_groups.size(); j++) {
      for (size_t k = 0; k < cvcs[i]->atom_groups[j]->size(); k++) {
        temp_id_list.push_back (cvcs[i]->atom_groups[j]->at(k).id);
      }
    }
  }

  temp_id_list.sort();
  temp_id_list.unique();

  // atom_ids = std::vector<int> (temp_id_list.begin(), temp_id_list.end());
  unsigned int id_i = 0;
  std::list<int>::iterator li;
  for (li = temp_id_list.begin(); li != temp_id_list.end(); ++li) {
    atom_ids[id_i] = *li;
    id_i++;
  }

  temp_id_list.clear();

  atomic_gradients.resize (atom_ids.size());
  if (atom_ids.size()) {
    if (cvm::debug())
      cvm::log ("Colvar: created atom list with " + cvm::to_str(atom_ids.size()) + " atoms.\n");
  } else {
    cvm::log ("Warning: colvar components communicated no atom IDs.\n");
  }
}


int colvar::parse_analysis (std::string const &conf)
{

  //   if (cvm::debug())
  //     cvm::log ("Parsing analysis flags for collective variable \""+
  //               this->name+"\".\n");

  runave_length = 0;
  bool b_runave = false;
  if (get_keyval (conf, "runAve", b_runave) && b_runave) {

    enable (task_runave);

    get_keyval (conf, "runAveLength", runave_length, 1000);
    get_keyval (conf, "runAveStride", runave_stride, 1);

    if ((cvm::restart_out_freq % runave_stride) != 0) {
      cvm::error("Error: runAveStride must be commensurate with the restart frequency.\n", INPUT_ERROR);
    }

    std::string runave_outfile;
    get_keyval (conf, "runAveOutputFile", runave_outfile,
                std::string (cvm::output_prefix+"."+
                             this->name+".runave.traj"));

    size_t const this_cv_width = x.output_width (cvm::cv_width);
    runave_os.open (runave_outfile.c_str());
    runave_os << "# " << cvm::wrap_string ("step", cvm::it_width-2)
              << "  "
              << cvm::wrap_string ("running average", this_cv_width)
              << " "
              << cvm::wrap_string ("running stddev", this_cv_width)
              << "\n";
  }

  acf_length = 0;
  bool b_acf = false;
  if (get_keyval (conf, "corrFunc", b_acf) && b_acf) {

    enable (task_corrfunc);

    std::string acf_colvar_name;
    get_keyval (conf, "corrFuncWithColvar", acf_colvar_name, this->name);
    if (acf_colvar_name == this->name) {
      cvm::log ("Calculating auto-correlation function.\n");
    } else {
      cvm::log ("Calculating correlation function with \""+
                this->name+"\".\n");
    }

    std::string acf_type_str;
    get_keyval (conf, "corrFuncType", acf_type_str, to_lower_cppstr (std::string ("velocity")));
    if (acf_type_str == to_lower_cppstr (std::string ("coordinate"))) {
      acf_type = acf_coor;
    } else if (acf_type_str == to_lower_cppstr (std::string ("velocity"))) {
      acf_type = acf_vel;
      enable (task_fdiff_velocity);
      if (acf_colvar_name.size())
        (cvm::colvar_by_name (acf_colvar_name))->enable (task_fdiff_velocity);
    } else if (acf_type_str == to_lower_cppstr (std::string ("coordinate_p2"))) {
      acf_type = acf_p2coor;
    } else {
      cvm::log ("Unknown type of correlation function, \""+
                        acf_type_str+"\".\n");
      cvm::set_error_bits(INPUT_ERROR);
    }

    get_keyval (conf, "corrFuncOffset", acf_offset, 0);
    get_keyval (conf, "corrFuncLength", acf_length, 1000);
    get_keyval (conf, "corrFuncStride", acf_stride, 1);

    if ((cvm::restart_out_freq % acf_stride) != 0) {
      cvm::error("Error: corrFuncStride must be commensurate with the restart frequency.\n", INPUT_ERROR);
    }

    get_keyval (conf, "corrFuncNormalize", acf_normalize, true);
    get_keyval (conf, "corrFuncOutputFile", acf_outfile,
                std::string (cvm::output_prefix+"."+this->name+
                             ".corrfunc.dat"));
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvar::enable (colvar::task const &t)
{
  switch (t) {

  case task_output_system_force:
    enable (task_system_force);
    break;

  case task_report_Jacobian_force:
    enable (task_Jacobian_force);
    enable (task_system_force);
    if (cvm::debug())
      cvm::log ("Adding the Jacobian force to the system force, "
                "rather than applying its opposite silently.\n");
    break;

  case task_Jacobian_force:
    // checks below do not apply to extended-system colvars
    if ( !tasks[task_extended_lagrangian] ) {
      enable (task_gradients);

      if (!b_Jacobian_force) {
        cvm::error ("Error: colvar \""+this->name+
                          "\" does not have Jacobian forces implemented.\n",
                  INPUT_ERROR);
      }
      if (!b_linear) {
        cvm::error ("Error: colvar \""+this->name+
                          "\" must be defined as a linear combination "
                          "to calculate the Jacobian force.\n",
                  INPUT_ERROR);
      }
      if (cvm::debug())
        cvm::log ("Enabling calculation of the Jacobian force "
                  "on this colvar.\n");
    }
    fj.type (this->type());
    break;

  case task_system_force:
    if (!tasks[task_extended_lagrangian]) {
      if (!b_inverse_gradients) {
        cvm::error ("Error: one or more of the components of "
                          "colvar \""+this->name+
                          "\" does not support system force calculation.\n",
                  INPUT_ERROR);
      }
      cvm::request_system_force();
    }
    ft.type (this->type());
    ft_reported.type (this->type());
    break;

  case task_output_applied_force:
  case task_lower_wall:
  case task_upper_wall:
    // all of the above require gradients
    enable (task_gradients);
    break;

  case task_fdiff_velocity:
    x_old.type (this->type());
    v_fdiff.type (this->type());
    v_reported.type (this->type());
    break;

  case task_output_velocity:
    enable (task_fdiff_velocity);
    break;

  case task_grid:
    if (this->type() != colvarvalue::type_scalar) {
      cvm::error ("Cannot calculate a grid for collective variable, \""+
                        this->name+"\", because its value is not a scalar number.\n",
                  INPUT_ERROR);
    }
    break;

  case task_extended_lagrangian:
    enable (task_gradients);
    v_reported.type (this->type());
    break;

  case task_lower_boundary:
  case task_upper_boundary:
    if (this->type() != colvarvalue::type_scalar) {
      cvm::error ("Error: this colvar is not a scalar value "
                        "and cannot produce a grid.\n",
                INPUT_ERROR);
    }
    break;

  case task_output_value:
  case task_runave:
  case task_corrfunc:
  case task_ntot:
  case task_langevin:
  case task_output_energy:
  case task_scripted:
    break;

  case task_gradients:
    f.type  (this->type());
    fb.type (this->type());
    break;

  case task_collect_gradients:
    if (this->type() != colvarvalue::type_scalar) {
      cvm::error ("Collecting atomic gradients for non-scalar collective variable \""+
                        this->name+"\" is not yet implemented.\n",
                  INPUT_ERROR);
    }

    enable (task_gradients);
    if (atom_ids.size() == 0) {
      build_atom_list();
    }
    break;
  }

  tasks[t] = true;
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


void colvar::disable (colvar::task const &t)
{
  // check dependencies
  switch (t) {
  case task_gradients:
    disable (task_upper_wall);
    disable (task_lower_wall);
    disable (task_output_applied_force);
    disable (task_system_force);
    disable (task_Jacobian_force);
    break;

  case task_system_force:
    disable (task_output_system_force);
    break;

  case task_Jacobian_force:
    disable (task_report_Jacobian_force);
    break;

  case task_fdiff_velocity:
    disable (task_output_velocity);
    break;

  case task_lower_boundary:
  case task_upper_boundary:
    disable (task_grid);
    break;

  case task_extended_lagrangian:
  case task_report_Jacobian_force:
  case task_output_value:
  case task_output_velocity:
  case task_output_applied_force:
  case task_output_system_force:
  case task_runave:
  case task_corrfunc:
  case task_grid:
  case task_lower_wall:
  case task_upper_wall:
  case task_ntot:
  case task_langevin:
  case task_output_energy:
  case task_collect_gradients:
  case task_scripted:
    break;
  }

  tasks[t] = false;
}


void colvar::setup() {
  // loop over all components to reset masses of all groups
  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
      cvm::atom_group &atoms = *(cvcs[i]->atom_groups[ig]);
      atoms.read_positions();
      atoms.reset_mass(name,i,ig);
    }
  }
}


colvar::~colvar()
{
  for (size_t i = 0; i < cvcs.size(); i++) {
    delete cvcs[i];
  }

  // remove reference to this colvar from the CVM
  for (std::vector<colvar *>::iterator cvi = cvm::colvars.begin();
       cvi != cvm::colvars.end();
       ++cvi) {
    if ( *cvi == this) {
      cvm::colvars.erase (cvi);
      break;
    }
  }
}



// ******************** CALC FUNCTIONS ********************


void colvar::calc()
{
  size_t i, ig;
  if (cvm::debug())
    cvm::log ("Calculating colvar \""+this->name+"\".\n");

  // prepare atom groups for calculation
  if (cvm::debug())
    cvm::log ("Collecting data from atom groups.\n");
  for (i = 0; i < cvcs.size(); i++) {
    for (ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
      cvm::atom_group &atoms = *(cvcs[i]->atom_groups[ig]);
      atoms.reset_atoms_data();
      atoms.read_positions();
      if (atoms.b_center || atoms.b_rotate) {
        atoms.calc_apply_roto_translation();
      }
      // each atom group will take care of its own ref_pos_group, if defined
    }
  }
  if (tasks[task_output_velocity]) {
    for (i = 0; i < cvcs.size(); i++) {
      for (ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
        cvcs[i]->atom_groups[ig]->read_velocities();
      }
    }
  }
  if (tasks[task_system_force]) {
    for (i = 0; i < cvcs.size(); i++) {
      for (ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
        cvcs[i]->atom_groups[ig]->read_system_forces();
      }
    }
  }

  // calculate the value of the colvar

  if (cvm::debug())
    cvm::log ("Calculating colvar components.\n");
  x.reset();

  // First, update component values
  for (i = 0; i < cvcs.size(); i++) {
    cvm::increase_depth();
    (cvcs[i])->calc_value();
    cvm::decrease_depth();
    if (cvm::debug())
      cvm::log ("Colvar component no. "+cvm::to_str (i+1)+
                " within colvar \""+this->name+"\" has value "+
                cvm::to_str ((cvcs[i])->value(),
                cvm::cv_width, cvm::cv_prec)+".\n");
  }

  // Then combine them appropriately
  if (tasks[task_scripted]) {
    // cvcs combined by user script
    int res = cvm::proxy->run_colvar_callback(scripted_function, sorted_cvc_values, x);
    if (res == COLVARS_NOT_IMPLEMENTED) {
      cvm::error("Scripted colvars are not implemented.");
      return;
    }
    if (res != COLVARS_OK) {
      cvm::error("Error running scripted colvar");
      return;
    }
  } else if (x.type() == colvarvalue::type_scalar) {
    // polynomial combination allowed
    for (i = 0; i < cvcs.size(); i++) {
      x += (cvcs[i])->sup_coeff *
      ( ((cvcs[i])->sup_np != 1) ?
        std::pow ((cvcs[i])->value().real_value, (cvcs[i])->sup_np) :
        (cvcs[i])->value().real_value );
    }
  } else {
    // only linear combination allowed
    for (i = 0; i < cvcs.size(); i++) {
      x += (cvcs[i])->sup_coeff * (cvcs[i])->value();
    }
  }

  if (cvm::debug())
    cvm::log ("Colvar \""+this->name+"\" has value "+
              cvm::to_str (x, cvm::cv_width, cvm::cv_prec)+".\n");

  if (tasks[task_gradients]) {

    if (cvm::debug())
      cvm::log ("Calculating gradients of colvar \""+this->name+"\".\n");

    for (i = 0; i < cvcs.size(); i++) {
      // calculate the gradients of each component
      cvm::increase_depth();

      (cvcs[i])->calc_gradients();

      // if requested, propagate (via chain rule) the gradients above
      // to the atoms used to define the roto-translation
      for (ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
        if (cvcs[i]->atom_groups[ig]->b_fit_gradients)
          cvcs[i]->atom_groups[ig]->calc_fit_gradients();
      }

      cvm::decrease_depth();
    }

    if (cvm::debug())
      cvm::log ("Done calculating gradients of colvar \""+this->name+"\".\n");

    if (tasks[task_collect_gradients]) {

      if (tasks[task_scripted]) {
        cvm::error("Collecting atomic gradients is not implemented for "
          "scripted colvars.");
        return;
      }

      // Collect the atomic gradients inside colvar object
      for (unsigned int a = 0; a < atomic_gradients.size(); a++) {
        atomic_gradients[a].reset();
      }
      for (i = 0; i < cvcs.size(); i++) {
        // Coefficient: d(a * x^n) = a * n * x^(n-1) * dx
        cvm::real coeff = (cvcs[i])->sup_coeff * cvm::real ((cvcs[i])->sup_np) *
          std::pow ((cvcs[i])->value().real_value, (cvcs[i])->sup_np-1);

        for (size_t j = 0; j < cvcs[i]->atom_groups.size(); j++) {

          // If necessary, apply inverse rotation to get atomic
          // gradient in the laboratory frame
          if (cvcs[i]->atom_groups[j]->b_rotate) {
            cvm::rotation const rot_inv = cvcs[i]->atom_groups[j]->rot.inverse();

            for (size_t k = 0; k < cvcs[i]->atom_groups[j]->size(); k++) {
              int a = std::lower_bound (atom_ids.begin(), atom_ids.end(),
                  cvcs[i]->atom_groups[j]->at(k).id) - atom_ids.begin();
              atomic_gradients[a] += coeff *
                rot_inv.rotate (cvcs[i]->atom_groups[j]->at(k).grad);
            }

          } else {

            for (size_t k = 0; k < cvcs[i]->atom_groups[j]->size(); k++) {
              int a = std::lower_bound (atom_ids.begin(), atom_ids.end(),
                  cvcs[i]->atom_groups[j]->at(k).id) - atom_ids.begin();
              atomic_gradients[a] += coeff * cvcs[i]->atom_groups[j]->at(k).grad;
            }
          }
        }
      }
    }
  }

  if (tasks[task_system_force]) {

    if (tasks[task_scripted]) {
      // TODO see if this could reasonably be done in a generic way
      // from generic inverse gradients
      cvm::error("System force is not implemented for "
        "scripted colvars.");
      return;
    }
    if (cvm::debug())
      cvm::log ("Calculating system force of colvar \""+this->name+"\".\n");

    ft.reset();

    // if(!tasks[task_extended_lagrangian] && (cvm::step_relative() > 0)) {
    // Disabled check to allow for explicit system force calculation
    // even with extended Lagrangian

    if(cvm::step_relative() > 0) {
      // get from the cvcs the system forces from the PREVIOUS step
      for (i = 0; i < cvcs.size(); i++) {
        (cvcs[i])->calc_force_invgrads();
        // linear combination is assumed
        cvm::increase_depth();
        ft += (cvcs[i])->system_force() / ((cvcs[i])->sup_coeff * cvm::real (cvcs.size()));
        cvm::decrease_depth();
      }
    }

    if (tasks[task_report_Jacobian_force]) {
      // add the Jacobian force to the system force, and don't apply any silent
      // correction internally: biases such as colvarbias_abf will handle it
      ft += fj;
    }

    if (cvm::debug())
      cvm::log ("Done calculating system force of colvar \""+this->name+"\".\n");
  }

  if (tasks[task_fdiff_velocity]) {
    // calculate the velocity by finite differences
    if (cvm::step_relative() == 0)
      x_old = x;
    else {
      v_fdiff = fdiff_velocity (x_old, x);
      v_reported = v_fdiff;
    }
  }

  if (tasks[task_extended_lagrangian]) {

    // initialize the restraint center in the first step to the value
    // just calculated from the cvcs
    // TODO: put it in the restart information
    if (cvm::step_relative() == 0) {
      xr = x;
      vr = 0.0; // (already 0; added for clarity)
    }

    // report the restraint center as "value"
    x_reported = xr;
    v_reported = vr;
    // the "system force" with the extended Lagrangian is just the
    // harmonic term acting on the extended coordinate
    // Note: this is the force for current timestep
    ft_reported = (-0.5 * ext_force_k) * this->dist2_lgrad (xr, x);

  } else {

    x_reported = x;
    ft_reported = ft;
  }

  if (cvm::debug())
    cvm::log ("Done calculating colvar \""+this->name+"\".\n");
}


cvm::real colvar::update()
{
  if (cvm::debug())
    cvm::log ("Updating colvar \""+this->name+"\".\n");


  // set to zero the applied force
  f.reset();

  // add the biases' force, which at this point should already have
  // been summed over each bias using this colvar
  f += fb;


  if (tasks[task_lower_wall] || tasks[task_upper_wall]) {

    // wall force
    colvarvalue fw (this->type());

    // if the two walls are applied concurrently, decide which is the
    // closer one (on a periodic colvar, both walls may be applicable
    // at the same time)
    if ( (!tasks[task_upper_wall]) ||
         (this->dist2 (x, lower_wall) < this->dist2 (x, upper_wall)) ) {

      cvm::real const grad = this->dist2_lgrad (x, lower_wall);
      if (grad < 0.0) {
        fw = -0.5 * lower_wall_k * grad;
        if (cvm::debug())
          cvm::log ("Applying a lower wall force ("+
                    cvm::to_str (fw)+") to \""+this->name+"\".\n");
        f += fw;

      }

    } else {

      cvm::real const grad = this->dist2_lgrad (x, upper_wall);
      if (grad > 0.0) {
        fw = -0.5 * upper_wall_k * grad;
        if (cvm::debug())
          cvm::log ("Applying an upper wall force ("+
                    cvm::to_str (fw)+") to \""+this->name+"\".\n");
        f += fw;
      }
    }
  }

  if (tasks[task_Jacobian_force]) {
    size_t i;
    cvm::increase_depth();
    for (i = 0; i < cvcs.size(); i++) {
      (cvcs[i])->calc_Jacobian_derivative();
    }
    cvm::decrease_depth();

    fj.reset();
    for (i = 0; i < cvcs.size(); i++) {
      // linear combination is assumed
      fj += 1.0 / ( cvm::real (cvcs.size()) *  cvm::real ((cvcs[i])->sup_coeff) ) *
        (cvcs[i])->Jacobian_derivative();
    }
    fj *= cvm::boltzmann() * cvm::temperature();

    // the instantaneous Jacobian force was not included in the reported system force;
    // instead, it is subtracted from the applied force (silent Jacobian correction)
    if (! tasks[task_report_Jacobian_force])
      f -= fj;
  }


  if (tasks[task_extended_lagrangian]) {

    cvm::real dt = cvm::dt();
    cvm::real f_ext;

    // the total force is applied to the fictitious mass, while the
    // atoms only feel the harmonic force
    // fr: extended coordinate force (without harmonic spring), for output in trajectory
    // f_ext: total force on extended coordinate (including harmonic spring)
    // f: - initially, external biasing force
    //    -  after this code block, colvar force to be applied to atomic coordinates, ie. spring force
    fr    = f;
    f_ext = f + (-0.5 * ext_force_k) * this->dist2_lgrad (xr, x);
    f     =     (-0.5 * ext_force_k) * this->dist2_rgrad (xr, x);

    // leapfrog: starting from x_i, f_i, v_(i-1/2)
    vr  += (0.5 * dt) * f_ext / ext_mass;
    // Because of leapfrog, kinetic energy at time i is approximate
    kinetic_energy = 0.5 * ext_mass * vr * vr;
    potential_energy = 0.5 * ext_force_k * this->dist2(xr, x);
    // leap to v_(i+1/2)
    if (tasks[task_langevin]) {
      vr -= dt * ext_gamma * vr.real_value;
      vr += dt * ext_sigma * cvm::rand_gaussian() / ext_mass;
    }
    vr  += (0.5 * dt) * f_ext / ext_mass;
    xr  += dt * vr;
    xr.apply_constraints();
    if (this->b_periodic) this->wrap (xr);
  }


  if (tasks[task_fdiff_velocity]) {
    // set it for the next step
    x_old = x;
  }

  if (cvm::debug())
    cvm::log ("Done updating colvar \""+this->name+"\".\n");
  return (potential_energy + kinetic_energy);
}


void colvar::communicate_forces()
{
  size_t i;
  if (cvm::debug())
    cvm::log ("Communicating forces from colvar \""+this->name+"\".\n");

  if (tasks[task_scripted]) {
    std::vector<colvarvalue> func_grads(cvcs.size());
    int res = cvm::proxy->run_colvar_gradient_callback(scripted_function, sorted_cvc_values, func_grads);

    if (res == COLVARS_NOT_IMPLEMENTED) {
      cvm::error("Colvar gradient scripts are not implemented.");
      return;
    }
    if (res != COLVARS_OK) {
      cvm::error("Error running colvar gradient script");
      return;
    }

    for (i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      // Note: we need a dot product here
      (cvcs[i])->apply_force (f * func_grads[i]);
      cvm::decrease_depth();
    }
  } else if (x.type() == colvarvalue::type_scalar) {

    for (i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->apply_force (f * (cvcs[i])->sup_coeff *
                              cvm::real ((cvcs[i])->sup_np) *
                              (std::pow ((cvcs[i])->value().real_value,
                                      (cvcs[i])->sup_np-1)) );
      cvm::decrease_depth();
    }

  } else {

    for (i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->apply_force (f * (cvcs[i])->sup_coeff);
      cvm::decrease_depth();
    }
  }

  if (cvm::debug())
    cvm::log ("Done communicating forces from colvar \""+this->name+"\".\n");
}



// ******************** METRIC FUNCTIONS ********************
// Use the metrics defined by \link cvc \endlink objects


bool colvar::periodic_boundaries (colvarvalue const &lb, colvarvalue const &ub) const
{
  if ( (!tasks[task_lower_boundary]) || (!tasks[task_upper_boundary]) ) {
    cvm::log ("Error: requesting to histogram the "
                      "collective variable \""+this->name+"\", but a "
                      "pair of lower and upper boundaries must be "
                      "defined.\n");
    cvm::set_error_bits(INPUT_ERROR);
  }

  if (period > 0.0) {
    if ( ((std::sqrt (this->dist2 (lb, ub))) / this->width)
         < 1.0E-10 ) {
      return true;
    }
  }

  return false;
}

bool colvar::periodic_boundaries() const
{
  if ( (!tasks[task_lower_boundary]) || (!tasks[task_upper_boundary]) ) {
    cvm::error ("Error: requesting to histogram the "
                      "collective variable \""+this->name+"\", but a "
                      "pair of lower and upper boundaries must be "
                      "defined.\n");
  }

  return periodic_boundaries (lower_boundary, upper_boundary);
}


cvm::real colvar::dist2 (colvarvalue const &x1,
                         colvarvalue const &x2) const
{
  return (cvcs[0])->dist2 (x1, x2);
}

colvarvalue colvar::dist2_lgrad (colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  return (cvcs[0])->dist2_lgrad (x1, x2);
}

colvarvalue colvar::dist2_rgrad (colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  return (cvcs[0])->dist2_rgrad (x1, x2);
}

cvm::real colvar::compare (colvarvalue const &x1,
                           colvarvalue const &x2) const
{
  return (cvcs[0])->compare (x1, x2);
}

void colvar::wrap (colvarvalue &x) const
{
  (cvcs[0])->wrap (x);
  return;
}


// ******************** INPUT FUNCTIONS ********************

std::istream & colvar::read_restart (std::istream &is)
{
  size_t const start_pos = is.tellg();

  std::string conf;
  if ( !(is >> colvarparse::read_block ("colvar", conf)) ) {
    // this is not a colvar block
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  {
    std::string check_name = "";
    if ( (get_keyval (conf, "name", check_name,
                      std::string (""), colvarparse::parse_silent)) &&
         (check_name != name) )  {
      cvm::error ("Error: the state file does not match the "
                        "configuration file, at colvar \""+name+"\".\n");
    }
    if (check_name.size() == 0) {
      cvm::error ("Error: Collective variable in the "
                        "restart file without any identifier.\n");
    }
  }

  if ( !(get_keyval (conf, "x", x,
                     colvarvalue (x.type()), colvarparse::parse_silent)) ) {
    cvm::log ("Error: restart file does not contain "
              "the value of the colvar \""+
              name+"\" .\n");
  } else {
    cvm::log ("Restarting collective variable \""+name+"\" from value: "+
              cvm::to_str (x)+"\n");
  }

  if (tasks[colvar::task_extended_lagrangian]) {

    if ( !(get_keyval (conf, "extended_x", xr,
                       colvarvalue (x.type()), colvarparse::parse_silent)) &&
         !(get_keyval (conf, "extended_v", vr,
                       colvarvalue (x.type()), colvarparse::parse_silent)) ) {
      cvm::log ("Error: restart file does not contain "
                "\"extended_x\" or \"extended_v\" for the colvar \""+
                name+"\", but you requested \"extendedLagrangian\".\n");
    }
  }

  if (tasks[task_extended_lagrangian]) {
    x_reported = xr;
  } else {
    x_reported = x;
  }

  if (tasks[task_output_velocity]) {

    if ( !(get_keyval (conf, "v", v_fdiff,
                       colvarvalue (x.type()), colvarparse::parse_silent)) ) {
      cvm::log ("Error: restart file does not contain "
                "the velocity for the colvar \""+
                name+"\", but you requested \"outputVelocity\".\n");
    }

    if (tasks[task_extended_lagrangian]) {
      v_reported = vr;
    } else {
      v_reported = v_fdiff;
    }
  }

  return is;
}


std::istream & colvar::read_traj (std::istream &is)
{
  size_t const start_pos = is.tellg();

  if (tasks[task_output_value]) {

    if (!(is >> x)) {
      cvm::log ("Error: in reading the value of colvar \""+
                this->name+"\" from trajectory.\n");
      is.clear();
      is.seekg (start_pos, std::ios::beg);
      is.setstate (std::ios::failbit);
      return is;
    }

    if (tasks[task_extended_lagrangian]) {
      is >> xr;
      x_reported = xr;
    } else {
      x_reported = x;
    }
  }

  if (tasks[task_output_velocity]) {

    is >> v_fdiff;

    if (tasks[task_extended_lagrangian]) {
      is >> vr;
      v_reported = vr;
    } else {
      v_reported = v_fdiff;
    }
  }

  if (tasks[task_output_system_force]) {
    is >> ft;
    ft_reported = ft;
  }

  if (tasks[task_output_applied_force]) {
    is >> f;
  }

  return is;
}


// ******************** OUTPUT FUNCTIONS ********************

std::ostream & colvar::write_restart (std::ostream &os) {

  os << "colvar {\n"
     << "  name " << name << "\n"
     << "  x "
     << std::setprecision (cvm::cv_prec)
     << std::setw (cvm::cv_width)
     << x << "\n";

  if (tasks[task_output_velocity]) {
    os << "  v "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << v_reported << "\n";
  }

  if (tasks[task_extended_lagrangian]) {
    os << "  extended_x "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << xr << "\n"
       << "  extended_v "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << vr << "\n";
  }

  os << "}\n\n";

  return os;
}


std::ostream & colvar::write_traj_label (std::ostream & os)
{
  size_t const this_cv_width = x.output_width (cvm::cv_width);

  os << " ";

  if (tasks[task_output_value]) {

    os << " "
       << cvm::wrap_string (this->name, this_cv_width);

    if (tasks[task_extended_lagrangian]) {
      // extended DOF
      os << " r_"
         << cvm::wrap_string (this->name, this_cv_width-2);
    }
  }

  if (tasks[task_output_velocity]) {

    os << " v_"
       << cvm::wrap_string (this->name, this_cv_width-2);

    if (tasks[task_extended_lagrangian]) {
      // extended DOF
      os << " vr_"
         << cvm::wrap_string (this->name, this_cv_width-3);
    }
  }

  if (tasks[task_output_energy]) {
    os << " Ep_"
       << cvm::wrap_string (this->name, this_cv_width-3)
       << " Ek_"
       << cvm::wrap_string (this->name, this_cv_width-3);
  }

  if (tasks[task_output_system_force]) {
    os << " fs_"
       << cvm::wrap_string (this->name, this_cv_width-3);
  }

  if (tasks[task_output_applied_force]) {
    os << " fa_"
       << cvm::wrap_string (this->name, this_cv_width-3);
  }

  return os;
}


std::ostream & colvar::write_traj (std::ostream &os)
{
  os << " ";

  if (tasks[task_output_value]) {

    if (tasks[task_extended_lagrangian]) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << x;
    }

    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << x_reported;
  }

  if (tasks[task_output_velocity]) {

    if (tasks[task_extended_lagrangian]) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << v_fdiff;
    }

    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << v_reported;
  }

  if (tasks[task_output_energy]) {
    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << potential_energy
       << " "
       << kinetic_energy;
  }


  if (tasks[task_output_system_force]) {
    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << ft_reported;
  }

  if (tasks[task_output_applied_force]) {
    if (tasks[task_extended_lagrangian]) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << fr;
    } else {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << f;
    }
  }

  return os;
}

int colvar::write_output_files()
{
  if (cvm::b_analysis) {

    if (acf.size()) {
      cvm::log ("Writing acf to file \""+acf_outfile+"\".\n");

      std::ofstream acf_os (acf_outfile.c_str());
      if (! acf_os.good()) {
        cvm::error("Cannot open file \""+acf_outfile+"\".\n", FILE_ERROR);
      }
      write_acf (acf_os);
      acf_os.close();
    }

    if (runave_os.good()) {
      runave_os.close();
    }
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}



// ******************** ANALYSIS FUNCTIONS ********************

void colvar::analyse()
{
  if (tasks[task_runave]) {
    calc_runave();
  }

  if (tasks[task_corrfunc]) {
    calc_acf();
  }
}


inline void history_add_value (size_t const           &history_length,
                               std::list<colvarvalue> &history,
                               colvarvalue const      &new_value)
{
  history.push_front (new_value);
  if (history.size() > history_length)
    history.pop_back();
}

inline void history_incr (std::list< std::list<colvarvalue> >           &history,
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

  if (acf_x_history.empty() && acf_v_history.empty()) {

    // first-step operations

    colvar *cfcv = (acf_colvar_name.size() ?
                    cvm::colvar_by_name (acf_colvar_name) :
                    this);
    if (cfcv->type() != this->type()) {
      cvm::error ("Error: correlation function between \""+cfcv->name+
                        "\" and \""+this->name+"\" cannot be calculated, "
                        "because their value types are different.\n",
                        INPUT_ERROR);
    }
    acf_nframes = 0;

    cvm::log ("Colvar \""+this->name+"\": initializing ACF calculation.\n");

    if (acf.size() < acf_length+1)
      acf.resize (acf_length+1, 0.0);

    size_t i;
    switch (acf_type) {

    case acf_vel:
      // allocate space for the velocities history
      for (i = 0; i < acf_stride; i++) {
        acf_v_history.push_back (std::list<colvarvalue>());
      }
      acf_v_history_p = acf_v_history.begin();
      break;

    case acf_coor:
    case acf_p2coor:
      // allocate space for the coordinates history
      for (i = 0; i < acf_stride; i++) {
        acf_x_history.push_back (std::list<colvarvalue>());
      }
      acf_x_history_p = acf_x_history.begin();
      break;

    default:
      break;
    }

  } else {

    colvar *cfcv = (acf_colvar_name.size() ?
                    cvm::colvar_by_name (acf_colvar_name) :
                    this);

    switch (acf_type) {

    case acf_vel:

      if (tasks[task_fdiff_velocity]) {
        // calc() should do this already, but this only happens in a
        // simulation; better do it again in case a trajectory is
        // being read
        v_reported = v_fdiff = fdiff_velocity (x_old, cfcv->value());
      }

      calc_vel_acf ((*acf_v_history_p), cfcv->velocity());
      // store this value in the history
      history_add_value (acf_length+acf_offset, *acf_v_history_p, cfcv->velocity());
      // if stride is larger than one, cycle among different histories
      history_incr (acf_v_history, acf_v_history_p);
      break;

    case acf_coor:

      calc_coor_acf ((*acf_x_history_p), cfcv->value());
      history_add_value (acf_length+acf_offset, *acf_x_history_p, cfcv->value());
      history_incr (acf_x_history, acf_x_history_p);
      break;

    case acf_p2coor:

      calc_p2coor_acf ((*acf_x_history_p), cfcv->value());
      history_add_value (acf_length+acf_offset, *acf_x_history_p, cfcv->value());
      history_incr (acf_x_history, acf_x_history_p);
      break;

    default:
      break;
    }
  }

  if (tasks[task_fdiff_velocity]) {
    // set it for the next step
    x_old = x;
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvar::calc_vel_acf (std::list<colvarvalue> &v_list,
                           colvarvalue const      &v)
{
  // loop over stored velocities and add to the ACF, but only the
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
    colvarvalue::inner_opt (v, vs_i, v_list.end(), acf_i);

    acf_nframes++;
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


void colvar::calc_coor_acf (std::list<colvarvalue> &x_list,
                            colvarvalue const      &x)
{
  // same as above but for coordinates
  if (x_list.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  xs_i = x_list.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      ++xs_i;

    *(acf_i++) += x.norm2();

    colvarvalue::inner_opt (x, xs_i, x_list.end(), acf_i);

    acf_nframes++;
  }
}


void colvar::calc_p2coor_acf (std::list<colvarvalue> &x_list,
                              colvarvalue const      &x)
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

    colvarvalue::p2leg_opt (x, xs_i, x_list.end(), acf_i);

    acf_nframes++;
  }
}


void colvar::write_acf (std::ostream &os)
{
  if (!acf_nframes)
    cvm::log ("Warning: ACF was not calculated (insufficient frames).\n");
  os.setf (std::ios::scientific, std::ios::floatfield);
  os << "# Autocorrelation function for collective variable \""
     << this->name << "\"\n";
  // one frame is used for normalization, the statistical sample is
  // hence decreased
  os << "# nframes = " << (acf_normalize ?
                           acf_nframes - 1 :
                           acf_nframes) << "\n";

  cvm::real const acf_norm = acf.front() / cvm::real (acf_nframes);
  std::vector<cvm::real>::iterator acf_i;
  size_t it = acf_offset;
  for (acf_i = acf.begin(); acf_i != acf.end(); ++acf_i) {
    os << std::setw (cvm::it_width) << acf_stride * (it++) << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << ( acf_normalize ?
            (*acf_i)/(acf_norm * cvm::real (acf_nframes)) :
            (*acf_i)/(cvm::real (acf_nframes)) ) << "\n";
  }
}


void colvar::calc_runave()
{
  if (x_history.empty()) {

    runave.type (x.type());
    runave.reset();

    // first-step operations

    if (cvm::debug())
      cvm::log ("Colvar \""+this->name+
                "\": initializing running average calculation.\n");

    acf_nframes = 0;

    x_history.push_back (std::list<colvarvalue>());
    x_history_p = x_history.begin();

  } else {

    if ( (cvm::step_relative() % runave_stride) == 0) {

      if ((*x_history_p).size() >= runave_length-1) {

        runave = x;
        std::list<colvarvalue>::iterator xs_i;
        for (xs_i = (*x_history_p).begin();
             xs_i != (*x_history_p).end(); ++xs_i) {
          runave += (*xs_i);
        }
        runave *= 1.0 / cvm::real (runave_length);
        runave.apply_constraints();

        runave_variance = 0.0;
        runave_variance += this->dist2 (x, runave);
        for (xs_i = (*x_history_p).begin();
             xs_i != (*x_history_p).end(); ++xs_i) {
          runave_variance += this->dist2 (x, (*xs_i));
        }
        runave_variance *= 1.0 / cvm::real (runave_length-1);

        runave_os << std::setw (cvm::it_width) << cvm::step_relative()
                  << "  "
                  << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
                  << runave << " "
                  << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
                  << std::sqrt (runave_variance) << "\n";
      }

      history_add_value (runave_length, *x_history_p, x);
    }
  }

}



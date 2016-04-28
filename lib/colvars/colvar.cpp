/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvarscript.h"
#include <algorithm>


/// Compare two cvcs using their names
/// Used to sort CVC array in scripted coordinates
bool compare(colvar::cvc *i, colvar::cvc *j) {
  return i->name < j->name;
}


colvar::colvar(std::string const &conf)
  : colvarparse(conf)
{
  size_t i;
  cvm::log("Initializing a new collective variable.\n");

  get_keyval(conf, "name", this->name,
              (std::string("colvar")+cvm::to_str(cvm::colvars.size()+1)));

  if (cvm::colvar_by_name(this->name) != NULL) {
    cvm::error("Error: this colvar cannot have the same name, \""+this->name+
                      "\", as another colvar.\n",
               INPUT_ERROR);
    return;
  }

  // Initialize dependency members
  // Could be a function defined in a different source file, for space?

  this->description = "colvar " + this->name;

  // Initialize static array once and for all
  init_cv_requires();

  kinetic_energy = 0.0;
  potential_energy = 0.0;

  // read the configuration and set up corresponding instances, for
  // each type of component implemented
#define initialize_components(def_desc,def_config_key,def_class_name)   \
  {                                                                     \
    size_t def_count = 0;                                               \
    std::string def_conf = "";                                          \
    size_t pos = 0;                                                     \
    while ( this->key_lookup(conf,                                     \
                              def_config_key,                           \
                              def_conf,                                 \
                              pos) ) {                                  \
      if (!def_conf.size()) continue;                                   \
      cvm::log("Initializing "                                         \
                "a new \""+std::string(def_config_key)+"\" component"+ \
                (cvm::debug() ? ", with configuration:\n"+def_conf      \
                 : ".\n"));                                             \
      cvm::increase_depth();                                            \
      cvc *cvcp = new colvar::def_class_name(def_conf);                \
      if (cvcp != NULL) {                                               \
        cvcs.push_back(cvcp);                                          \
        cvcp->check_keywords(def_conf, def_config_key);                \
        if (cvm::get_error()) {                                         \
          cvm::error("Error: in setting up component \""                \
                      def_config_key"\".\n");                           \
          return;                                                       \
        }                                                               \
        cvm::decrease_depth();                                          \
      } else {                                                          \
        cvm::error("Error: in allocating component \""                  \
                          def_config_key"\".\n",                        \
                          MEMORY_ERROR);                                \
        return;                                                         \
      }                                                                 \
      if ( (cvcp->period != 0.0) || (cvcp->wrap_center != 0.0) ) {      \
        if ( (cvcp->function_type != std::string("distance_z")) &&     \
             (cvcp->function_type != std::string("dihedral")) &&       \
             (cvcp->function_type != std::string("spin_angle")) ) {    \
          cvm::error("Error: invalid use of period and/or "            \
                            "wrapAround in a \""+                       \
                            std::string(def_config_key)+               \
                            "\" component.\n"+                          \
                            "Period: "+cvm::to_str(cvcp->period) +      \
                        " wrapAround: "+cvm::to_str(cvcp->wrap_center), \
                        INPUT_ERROR);                                   \
          return;                                                       \
        }                                                               \
      }                                                                 \
      if ( ! cvcs.back()->name.size()){                                 \
        std::ostringstream s;                                           \
        s << def_config_key << std::setfill('0') << std::setw(4) << ++def_count;\
        cvcs.back()->name = s.str();                                    \
          /* pad cvc number for correct ordering when sorting by name */\
      }                                                                 \
      cvcs.back()->setup();                                             \
      if (cvm::debug())                                                 \
        cvm::log("Done initializing a \""+                             \
                  std::string(def_config_key)+                         \
                  "\" component"+                                       \
                  (cvm::debug() ?                                       \
                   ", named \""+cvcs.back()->name+"\""                  \
                   : "")+".\n");                                        \
      def_conf = "";                                                    \
      if (cvm::debug())                                                 \
        cvm::log("Parsed "+cvm::to_str(cvcs.size())+                  \
                  " components at this time.\n");                       \
    }                                                                   \
  }


  initialize_components("distance",         "distance",       distance);
  initialize_components("distance vector",  "distanceVec",    distance_vec);
  initialize_components("Cartesian coordinates", "cartesian",  cartesian);
  initialize_components("distance vector "
                        "direction",        "distanceDir",    distance_dir);
  initialize_components("distance projection "
                        "on an axis",       "distanceZ",      distance_z);
  initialize_components("distance projection "
                        "on a plane",       "distanceXY",     distance_xy);
  initialize_components("average distance weighted by inverse power",
                        "distanceInv", distance_inv);
  initialize_components("N1xN2-long vector of pairwise distances",
                        "distancePairs", distance_pairs);

  initialize_components("coordination "
                        "number",           "coordNum",       coordnum);
  initialize_components("self-coordination "
                        "number",           "selfCoordNum",   selfcoordnum);

  initialize_components("angle",            "angle",          angle);
  initialize_components("dipole angle",     "dipoleAngle",    dipole_angle);
  initialize_components("dihedral",         "dihedral",       dihedral);

  initialize_components("hydrogen bond",    "hBond",          h_bond);

  //  initialize_components ("alpha helix",      "alphaDihedrals", alpha_dihedrals);
  initialize_components("alpha helix",      "alpha",          alpha_angles);

  initialize_components("dihedral principal "
                        "component",        "dihedralPC",     dihedPC);

  initialize_components("orientation",      "orientation",    orientation);
  initialize_components("orientation "
                        "angle",            "orientationAngle",orientation_angle);
  initialize_components("orientation "
                        "projection",       "orientationProj",orientation_proj);
  initialize_components("tilt",             "tilt",           tilt);
  initialize_components("spin angle",       "spinAngle",      spin_angle);

  initialize_components("RMSD",             "rmsd",           rmsd);

  //  initialize_components ("logarithm of MSD", "logmsd",         logmsd);

  initialize_components("radius of "
                        "gyration",         "gyration",       gyration);
  initialize_components("moment of "
                        "inertia",          "inertia",        inertia);
  initialize_components("moment of inertia around an axis",
                        "inertiaZ",       inertia_z);
  initialize_components("eigenvector",      "eigenvector",    eigenvector);

  if (!cvcs.size()) {
    cvm::error("Error: no valid components were provided "
               "for this collective variable.\n",
               INPUT_ERROR);
    return;
  }

  n_active_cvcs = cvcs.size();

  cvm::log("All components initialized.\n");

  // Store list of children cvcs for dependency checking purposes
  for (i = 0; i < cvcs.size(); i++) {
    add_child(cvcs[i]);
  }

  // Setup colvar as scripted function of components
  if (get_keyval(conf, "scriptedFunction", scripted_function,
    "", colvarparse::parse_silent)) {

    // Make feature available only on user request
    provide(f_cv_scripted);
    enable(f_cv_scripted);
    cvm::log("This colvar uses scripted function \"" + scripted_function + "\".");

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
      cvm::error("Could not parse scripted colvar type.");
      return;
    }
    x_reported.type(x.type());
    cvm::log(std::string("Expecting colvar value of type ")
      + colvarvalue::type_desc(x.type()));

    if (x.type() == colvarvalue::type_vector) {
      int size;
      if (!get_keyval(conf, "scriptedFunctionVectorSize", size)) {
        cvm::error("Error: no size specified for vector scripted function.");
        return;
      }
      x.vector1d_value.resize(size);
    }

    // Sort array of cvcs based on their names
    // Note: default CVC names are in input order for same type of CVC
    std::sort(cvcs.begin(), cvcs.end(), compare);

    if(cvcs.size() > 1) {
      cvm::log("Sorted list of components for this scripted colvar:");
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

  if (!is_enabled(f_cv_scripted)) {
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
  // If using scripted biases, any colvar may receive bias forces
  // and will need its gradient
  if (cvm::scripted_forces()) {
    enable(f_cv_gradient);
  }

  // check for linear combinations
  {
    bool lin = !is_enabled(f_cv_scripted);
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
    feature_states[f_cv_linear]->enabled = lin;
  }

  // Colvar is homogeneous iff:
  // - it is linear (hence not scripted)
  // - all cvcs have coefficient 1 or -1
  // i.e. sum or difference of cvcs
  {
    bool homogeneous = is_enabled(f_cv_linear);
    for (i = 0; i < cvcs.size(); i++) {
      if ((std::fabs(cvcs[i]->sup_coeff) - 1.0) > 1.0e-10) {
        homogeneous = false;
      }
    }
    feature_states[f_cv_homogeneous]->enabled = homogeneous;
  }
  // Colvar is deemed periodic iff:
  // - it is homogeneous
  // - all cvcs are periodic
  // - all cvcs have the same period

  b_periodic = cvcs[0]->b_periodic && is_enabled(f_cv_homogeneous);
  period = cvcs[0]->period;
  for (i = 1; i < cvcs.size(); i++) {
    if (!cvcs[i]->b_periodic || cvcs[i]->period != period) {
      b_periodic = false;
      period = 0.0;
    }
  }
  feature_states[f_cv_periodic]->enabled = b_periodic;

  // check that cvcs are compatible

  for (i = 0; i < cvcs.size(); i++) {
    if ((cvcs[i])->b_periodic && !b_periodic) {
        cvm::log("Warning: although this component is periodic, the colvar will "
                  "not be treated as periodic, either because the exponent is not "
                  "1, or because multiple components are present. Make sure that "
                  "you know what you are doing!");
    }

    // components may have different types only for scripted functions
    if (!is_enabled(f_cv_scripted) && (colvarvalue::check_types(cvcs[i]->value(),
                                                           cvcs[0]->value())) ) {
      cvm::error("ERROR: you are definining this collective variable "
                 "by using components of different types. "
                 "You must use the same type in order to "
                 " sum them together.\n", INPUT_ERROR);
      return;
    }
  }

  active_cvc_square_norm = 0.;
  for (i = 0; i < cvcs.size(); i++) {
    active_cvc_square_norm += cvcs[i]->sup_coeff * cvcs[i]->sup_coeff;
  }

  // at this point, the colvar's type is defined
  f.type(value());
  fb.type(value());

  get_keyval(conf, "width", width, 1.0);
  if (width <= 0.0) {
    cvm::error("Error: \"width\" must be positive.\n", INPUT_ERROR);
  }

  // NOTE: not porting wall stuff to new deps, as this will change to a separate bias
  // the grid functions will wait a little as well

  lower_boundary.type(value());
  lower_wall.type(value());

  upper_boundary.type(value());
  upper_wall.type(value());

  feature_states[f_cv_scalar]->enabled = (value().type() == colvarvalue::type_scalar);

  if (is_enabled(f_cv_scalar)) {

    if (get_keyval(conf, "lowerBoundary", lower_boundary, lower_boundary)) {
      provide(f_cv_lower_boundary);
      enable(f_cv_lower_boundary);
    }

    get_keyval(conf, "lowerWallConstant", lower_wall_k, 0.0);
    if (lower_wall_k > 0.0) {
      get_keyval(conf, "lowerWall", lower_wall, lower_boundary);
      enable(f_cv_lower_wall);
    }

    if (get_keyval(conf, "upperBoundary", upper_boundary, upper_boundary)) {
      provide(f_cv_upper_boundary);
      enable(f_cv_upper_boundary);
    }

    get_keyval(conf, "upperWallConstant", upper_wall_k, 0.0);
    if (upper_wall_k > 0.0) {
      get_keyval(conf, "upperWall", upper_wall, upper_boundary);
      enable(f_cv_upper_wall);
    }
  }

  if (is_enabled(f_cv_lower_boundary)) {
    get_keyval(conf, "hardLowerBoundary", hard_lower_boundary, false);
  }
  if (is_enabled(f_cv_upper_boundary)) {
    get_keyval(conf, "hardUpperBoundary", hard_upper_boundary, false);
  }

  // consistency checks for boundaries and walls
  if (is_enabled(f_cv_lower_boundary) && is_enabled(f_cv_upper_boundary)) {
    if (lower_boundary >= upper_boundary) {
      cvm::error("Error: the upper boundary, "+
                        cvm::to_str(upper_boundary)+
                        ", is not higher than the lower boundary, "+
                        cvm::to_str(lower_boundary)+".\n",
                INPUT_ERROR);
    }
  }

  if (is_enabled(f_cv_lower_wall) && is_enabled(f_cv_upper_wall)) {
    if (lower_wall >= upper_wall) {
      cvm::error("Error: the upper wall, "+
                 cvm::to_str(upper_wall)+
                 ", is not higher than the lower wall, "+
                 cvm::to_str(lower_wall)+".\n",
                 INPUT_ERROR);
    }
  }

  get_keyval(conf, "expandBoundaries", expand_boundaries, false);
  if (expand_boundaries && periodic_boundaries()) {
    cvm::error("Error: trying to expand boundaries that already "
               "cover a whole period of a periodic colvar.\n",
               INPUT_ERROR);
  }
  if (expand_boundaries && hard_lower_boundary && hard_upper_boundary) {
    cvm::error("Error: inconsistent configuration "
               "(trying to expand boundaries with both "
               "hardLowerBoundary and hardUpperBoundary enabled).\n",
               INPUT_ERROR);
  }

  {
    bool b_extended_Lagrangian;
    get_keyval(conf, "extendedLagrangian", b_extended_Lagrangian, false);

    if (b_extended_Lagrangian) {
      cvm::real temp, tolerance, period;

      cvm::log("Enabling the extended Lagrangian term for colvar \""+
                this->name+"\".\n");

      // Make feature available only on user request
      provide(f_cv_extended_Lagrangian);
      enable(f_cv_extended_Lagrangian);

      xr.type(value());
      vr.type(value());
      fr.type(value());

      const bool found = get_keyval(conf, "extendedTemp", temp, cvm::temperature());
      if (temp <= 0.0) {
        if (found)
          cvm::error("Error: \"extendedTemp\" must be positive.\n", INPUT_ERROR);
        else
          cvm::error("Error: a positive temperature must be provided, either "
                     "by enabling a thermostat, or through \"extendedTemp\".\n",
                     INPUT_ERROR);
      }

      get_keyval(conf, "extendedFluctuation", tolerance);
      if (tolerance <= 0.0) {
        cvm::error("Error: \"extendedFluctuation\" must be positive.\n", INPUT_ERROR);
      }
      ext_force_k = cvm::boltzmann() * temp / (tolerance * tolerance);
      cvm::log("Computed extended system force constant: " + cvm::to_str(ext_force_k) + " kcal/mol/U^2");

      get_keyval(conf, "extendedTimeConstant", period, 200.0);
      if (period <= 0.0) {
        cvm::error("Error: \"extendedTimeConstant\" must be positive.\n", INPUT_ERROR);
      }
      ext_mass = (cvm::boltzmann() * temp * period * period)
                 / (4.0 * PI * PI * tolerance * tolerance);
      cvm::log("Computed fictitious mass: " + cvm::to_str(ext_mass) + " kcal/mol/(U/fs)^2   (U: colvar unit)");

      {
        bool b_output_energy;
        get_keyval(conf, "outputEnergy", b_output_energy, false);
        if (b_output_energy) {
          enable(f_cv_output_energy);
        }
      }

      get_keyval(conf, "extendedLangevinDamping", ext_gamma, 1.0);
      if (ext_gamma < 0.0) {
        cvm::error("Error: \"extendedLangevinDamping\" may not be negative.\n", INPUT_ERROR);
      }
      if (ext_gamma != 0.0) {
        enable(f_cv_Langevin);
        ext_gamma *= 1.0e-3; // convert from ps-1 to fs-1
        ext_sigma = std::sqrt(2.0 * cvm::boltzmann() * temp * ext_gamma * ext_mass / cvm::dt());
      }
    }
  }

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
    bool b_output_system_force;
    get_keyval(conf, "outputSystemForce", b_output_system_force, false);
    if (b_output_system_force) {
      enable(f_cv_output_system_force);
    }
  }

  {
    bool b_output_applied_force;
    get_keyval(conf, "outputAppliedForce", b_output_applied_force, false);
    if (b_output_applied_force) {
      enable(f_cv_output_applied_force);
    }
  }

  // Start in active state by default
  enable(f_cv_active);
  // Make sure dependency side-effects are correct
  refresh_deps();

  x_old.type(value());
  v_fdiff.type(value());
  v_reported.type(value());
  fj.type(value());
  ft.type(value());
  ft_reported.type(value());

  if (cvm::b_analysis)
    parse_analysis(conf);

  if (cvm::debug())
    cvm::log("Done initializing collective variable \""+this->name+"\".\n");
}


int colvar::refresh_deps()
{
  // If enabled features are changed upstream, the features below should be refreshed
  if (is_enabled(f_cv_system_force_calc)) {
    cvm::request_system_force();
  }
  if (is_enabled(f_cv_collect_gradient) && atom_ids.size() == 0) {
    build_atom_list();
  }
  return COLVARS_OK;
}


void colvar::build_atom_list(void)
{
  // If atomic gradients are requested, build full list of atom ids from all cvcs
  std::list<int> temp_id_list;

  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t j = 0; j < cvcs[i]->atom_groups.size(); j++) {
      cvm::atom_group &ag = *(cvcs[i]->atom_groups[j]);
      for (size_t k = 0; k < ag.size(); k++) {
        temp_id_list.push_back(ag[k].id);
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
      cvm::error("Error: runAveStride must be commensurate with the restart frequency.\n", INPUT_ERROR);
    }

    std::string runave_outfile;
    get_keyval(conf, "runAveOutputFile", runave_outfile,
                std::string(cvm::output_prefix+"."+
                             this->name+".runave.traj"));

    size_t const this_cv_width = x.output_width(cvm::cv_width);
    cvm::backup_file(runave_outfile.c_str());
    runave_os.open(runave_outfile.c_str());
    runave_os << "# " << cvm::wrap_string("step", cvm::it_width-2)
              << "  "
              << cvm::wrap_string("running average", this_cv_width)
              << " "
              << cvm::wrap_string("running stddev", this_cv_width)
              << "\n";
  }

  acf_length = 0;
  bool b_acf = false;
  if (get_keyval(conf, "corrFunc", b_acf) && b_acf) {

    enable(f_cv_corrfunc);

    std::string acf_colvar_name;
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
      if (acf_colvar_name.size())
        (cvm::colvar_by_name(acf_colvar_name))->enable(f_cv_fdiff_velocity);
    } else if (acf_type_str == to_lower_cppstr(std::string("coordinate_p2"))) {
      acf_type = acf_p2coor;
    } else {
      cvm::log("Unknown type of correlation function, \""+
                        acf_type_str+"\".\n");
      cvm::set_error_bit(INPUT_ERROR);
    }

    get_keyval(conf, "corrFuncOffset", acf_offset, 0);
    get_keyval(conf, "corrFuncLength", acf_length, 1000);
    get_keyval(conf, "corrFuncStride", acf_stride, 1);

    if ((cvm::restart_out_freq % acf_stride) != 0) {
      cvm::error("Error: corrFuncStride must be commensurate with the restart frequency.\n", INPUT_ERROR);
    }

    get_keyval(conf, "corrFuncNormalize", acf_normalize, true);
    get_keyval(conf, "corrFuncOutputFile", acf_outfile,
                std::string(cvm::output_prefix+"."+this->name+
                             ".corrfunc.dat"));
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


void colvar::setup() {
  // loop over all components to reset masses of all groups
  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
      cvm::atom_group &atoms = *(cvcs[i]->atom_groups[ig]);
      atoms.setup();
      atoms.reset_mass(name,i,ig);
      atoms.read_positions();
    }
  }
}


colvar::~colvar()
{
//   Clear references to this colvar's cvcs as children
//   for dependency purposes
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

  // remove reference to this colvar from the CVM
  for (std::vector<colvar *>::iterator cvi = cvm::colvars.begin();
       cvi != cvm::colvars.end();
       ++cvi) {
    if ( *cvi == this) {
      cvm::colvars.erase(cvi);
      break;
    }
  }
}



// ******************** CALC FUNCTIONS ********************


// Default schedule (everything is serialized)
int colvar::calc()
{
  // Note: if anything is added here, it should be added also in the SMP block of calc_colvars()
  int error_code = COLVARS_OK;
  if (is_enabled(f_cv_active)) {
    cvm::combine_errors(error_code, update_cvc_flags());
    cvm::combine_errors(error_code, calc_cvcs());
    cvm::combine_errors(error_code, collect_cvc_data());
  }
  return error_code;
}


int colvar::calc_cvcs(int first_cvc, size_t num_cvcs)
{
  int error_code = COLVARS_OK;
  if (cvm::debug())
    cvm::log("Calculating colvar \""+this->name+"\", components "+
             cvm::to_str(first_cvc)+" through "+cvm::to_str(first_cvc+num_cvcs)+".\n");

  cvm::combine_errors(error_code, check_cvc_range(first_cvc, num_cvcs));
  if (error_code != COLVARS_OK) {
    return error_code;
  }

  cvm::combine_errors(error_code, calc_cvc_values(first_cvc, num_cvcs));
  cvm::combine_errors(error_code, calc_cvc_gradients(first_cvc, num_cvcs));
  cvm::combine_errors(error_code, calc_cvc_sys_forces(first_cvc, num_cvcs));
  cvm::combine_errors(error_code, calc_cvc_Jacobians(first_cvc, num_cvcs));

  if (cvm::debug())
    cvm::log("Done calculating colvar \""+this->name+"\".\n");

  return error_code;
}


int colvar::collect_cvc_data()
{
  if (cvm::debug())
    cvm::log("Calculating colvar \""+this->name+"\"'s properties.\n");

  int error_code = COLVARS_OK;

  cvm::combine_errors(error_code, collect_cvc_values());
  cvm::combine_errors(error_code, collect_cvc_gradients());
  cvm::combine_errors(error_code, collect_cvc_sys_forces());
  cvm::combine_errors(error_code, collect_cvc_Jacobians());

  cvm::combine_errors(error_code, calc_colvar_properties());

  if (cvm::debug())
    cvm::log("Done calculating colvar \""+this->name+"\"'s properties.\n");

  return error_code;
}


int colvar::check_cvc_range(int first_cvc, size_t num_cvcs)
{
  if ((first_cvc < 0) || (first_cvc >= ((int) cvcs.size()))) {
    cvm::error("Error: trying to address a component outside the "
               "range defined for colvar \""+name+"\".\n", BUG_ERROR);
    return BUG_ERROR;
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
  size_t i;

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
  } else if (x.type() == colvarvalue::type_scalar) {
    // polynomial combination allowed
    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      x += (cvcs[i])->sup_coeff *
      ( ((cvcs[i])->sup_np != 1) ?
        std::pow((cvcs[i])->value().real_value, (cvcs[i])->sup_np) :
        (cvcs[i])->value().real_value );
    }
  } else {
    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      x += (cvcs[i])->sup_coeff * (cvcs[i])->value();
    }
  }

  if (cvm::debug())
    cvm::log("Colvar \""+this->name+"\" has value "+
              cvm::to_str(x, cvm::cv_width, cvm::cv_prec)+".\n");

  return COLVARS_OK;
}


int colvar::calc_cvc_gradients(int first_cvc, size_t num_cvcs)
{
  size_t const cvc_max_count = num_cvcs ? num_cvcs : num_active_cvcs();
  size_t i, cvc_count;

  if (is_enabled(f_cv_gradient)) {

    if (cvm::debug())
      cvm::log("Calculating gradients of colvar \""+this->name+"\".\n");

    // calculate the gradients of each component
    cvm::increase_depth();
    for (i = first_cvc, cvc_count = 0;
        (i < cvcs.size()) && (cvc_count < cvc_max_count);
        i++) {
      if (!cvcs[i]->is_enabled()) continue;
      cvc_count++;
      (cvcs[i])->calc_gradients();
      // if requested, propagate (via chain rule) the gradients above
      // to the atoms used to define the roto-translation
      for (size_t ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
        if (cvcs[i]->atom_groups[ig]->b_fit_gradients)
          cvcs[i]->atom_groups[ig]->calc_fit_gradients();

        if (cvcs[i]->is_enabled(f_cvc_debug_gradient)) {
          cvm::log("Debugging gradients for " + cvcs[i]->description);
          cvcs[i]->debug_gradients(cvcs[i]->atom_groups[ig]);
        }
      }
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

    if (is_enabled(f_cv_scripted)) {
      cvm::error("Collecting atomic gradients is not implemented for "
                 "scripted colvars.", COLVARS_NOT_IMPLEMENTED);
      return COLVARS_NOT_IMPLEMENTED;
    }

    // Collect the atomic gradients inside colvar object
    for (unsigned int a = 0; a < atomic_gradients.size(); a++) {
      atomic_gradients[a].reset();
    }
    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      // Coefficient: d(a * x^n) = a * n * x^(n-1) * dx
      cvm::real coeff = (cvcs[i])->sup_coeff * cvm::real((cvcs[i])->sup_np) *
        std::pow((cvcs[i])->value().real_value, (cvcs[i])->sup_np-1);

      for (size_t j = 0; j < cvcs[i]->atom_groups.size(); j++) {

        cvm::atom_group &ag = *(cvcs[i]->atom_groups[j]);

        // If necessary, apply inverse rotation to get atomic
        // gradient in the laboratory frame
        if (ag.b_rotate) {
          cvm::rotation const rot_inv = ag.rot.inverse();

          for (size_t k = 0; k < ag.size(); k++) {
            size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                        ag[k].id) - atom_ids.begin();
            atomic_gradients[a] += coeff * rot_inv.rotate(ag[k].grad);
          }

        } else {

          for (size_t k = 0; k < ag.size(); k++) {
            size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                        ag[k].id) - atom_ids.begin();
            atomic_gradients[a] += coeff * ag[k].grad;
          }
        }
      }
    }
  }
  return COLVARS_OK;
}


int colvar::calc_cvc_sys_forces(int first_cvc, size_t num_cvcs)
{
  size_t const cvc_max_count = num_cvcs ? num_cvcs : num_active_cvcs();
  size_t i, cvc_count;

  if (is_enabled(f_cv_system_force) && !is_enabled(f_cv_extended_Lagrangian)) {
    // If extended Lagrangian is enabled, system force calculation is trivial
    // and done together with integration of the extended coordinate.

    if (is_enabled(f_cv_scripted)) {
      // TODO see if this could reasonably be done in a generic way
      // from generic inverse gradients
      cvm::error("System force is not implemented for "
                 "scripted colvars.", COLVARS_NOT_IMPLEMENTED);
      return COLVARS_NOT_IMPLEMENTED;
    }

    if (cvm::debug())
      cvm::log("Calculating system force of colvar \""+this->name+"\".\n");

    // if (!tasks[task_extended_lagrangian] && (cvm::step_relative() > 0)) {
   // Disabled check to allow for explicit system force calculation
    // even with extended Lagrangian

    if (cvm::step_relative() > 0) {
      cvm::increase_depth();
      // get from the cvcs the system forces from the PREVIOUS step
      for (i = first_cvc, cvc_count = 0;
          (i < cvcs.size()) && (cvc_count < cvc_max_count);
          i++) {
        if (!cvcs[i]->is_enabled()) continue;
        cvc_count++;
        (cvcs[i])->calc_force_invgrads();
      }
      cvm::decrease_depth();
    }

    if (cvm::debug())
      cvm::log("Done calculating system force of colvar \""+this->name+"\".\n");
  }

  return COLVARS_OK;
}


int colvar::collect_cvc_sys_forces()
{
  if (is_enabled(f_cv_system_force) && !is_enabled(f_cv_extended_Lagrangian)) {
    // If extended Lagrangian is enabled, system force calculation is trivial
    // and done together with integration of the extended coordinate.
    ft.reset();

    if (cvm::step_relative() > 0) {
      // get from the cvcs the system forces from the PREVIOUS step
      for (size_t i = 0; i < cvcs.size();  i++) {
        if (!cvcs[i]->is_enabled()) continue;
        // linear combination is assumed
        ft += (cvcs[i])->system_force() * (cvcs[i])->sup_coeff / active_cvc_square_norm;
      }
    }

    if (!is_enabled(f_cv_hide_Jacobian)) {
      // add the Jacobian force to the system force, and don't apply any silent
      // correction internally: biases such as colvarbias_abf will handle it
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
  if (is_enabled(f_cv_Jacobian)) {
    fj.reset();
    for (size_t i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      // linear combination is assumed
      fj += (cvcs[i])->Jacobian_derivative() * (cvcs[i])->sup_coeff / active_cvc_square_norm;
    }
    fj *= cvm::boltzmann() * cvm::temperature();
  }

  return COLVARS_OK;
}


int colvar::calc_colvar_properties()
{
  if (is_enabled(f_cv_fdiff_velocity)) {
    // calculate the velocity by finite differences
    if (cvm::step_relative() == 0)
      x_old = x;
    else {
      v_fdiff = fdiff_velocity(x_old, x);
      v_reported = v_fdiff;
    }
  }

  if (is_enabled(f_cv_extended_Lagrangian)) {

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
    ft_reported = (-0.5 * ext_force_k) * this->dist2_lgrad(xr, x);

  } else {

    x_reported = x;
    ft_reported = ft;
  }

  return COLVARS_OK;
}


cvm::real colvar::update_forces_energy()
{
  if (cvm::debug())
    cvm::log("Updating colvar \""+this->name+"\".\n");

  // set to zero the applied force
  f.type(value());
  f.reset();

  // add the biases' force, which at this point should already have
  // been summed over each bias using this colvar
  f += fb;

  if (is_enabled(f_cv_Jacobian)) {
    // the instantaneous Jacobian force was not included in the reported system force;
    // instead, it is subtracted from the applied force (silent Jacobian correction)
    if (is_enabled(f_cv_hide_Jacobian))
      f -= fj;
  }

  if (is_enabled(f_cv_extended_Lagrangian)) {

    cvm::real dt = cvm::dt();
    cvm::real f_ext;

    // the total force is applied to the fictitious mass, while the
    // atoms only feel the harmonic force
    // fr: bias force on extended coordinate (without harmonic spring), for output in trajectory
    // f_ext: total force on extended coordinate (including harmonic spring)
    // f: - initially, external biasing force
    //    - after this code block, colvar force to be applied to atomic coordinates, ie. spring force
    //      (note: wall potential is added to f after this block)
    fr    = f;
    f_ext = f + (-0.5 * ext_force_k) * this->dist2_lgrad(xr, x);
    f     =     (-0.5 * ext_force_k) * this->dist2_rgrad(xr, x);

    // leapfrog: starting from x_i, f_i, v_(i-1/2)
    vr  += (0.5 * dt) * f_ext / ext_mass;
    // Because of leapfrog, kinetic energy at time i is approximate
    kinetic_energy = 0.5 * ext_mass * vr * vr;
    potential_energy = 0.5 * ext_force_k * this->dist2(xr, x);
    // leap to v_(i+1/2)
    if (is_enabled(f_cv_Langevin)) {
      vr -= dt * ext_gamma * vr.real_value;
      vr += dt * ext_sigma * cvm::rand_gaussian() / ext_mass;
    }
    vr  += (0.5 * dt) * f_ext / ext_mass;
    xr  += dt * vr;
    xr.apply_constraints();
    if (this->b_periodic) this->wrap(xr);
  }


  // Adding wall potential to "true" colvar force, whether or not an extended coordinate is in use
  if (is_enabled(f_cv_lower_wall) || is_enabled(f_cv_upper_wall)) {

    // Wall force
    colvarvalue fw(x);
    fw.reset();

    if (cvm::debug())
      cvm::log("Calculating wall forces for colvar \""+this->name+"\".\n");

    // For a periodic colvar, both walls may be applicable at the same time
    // in which case we pick the closer one
    if ( (!is_enabled(f_cv_upper_wall)) ||
         (this->dist2(x, lower_wall) < this->dist2(x, upper_wall)) ) {

      cvm::real const grad = this->dist2_lgrad(x, lower_wall);
      if (grad < 0.0) {
        fw = -0.5 * lower_wall_k * grad;
        f += fw;
        if (cvm::debug())
          cvm::log("Applying a lower wall force("+
                    cvm::to_str(fw)+") to \""+this->name+"\".\n");
      }

    } else {

      cvm::real const grad = this->dist2_lgrad(x, upper_wall);
      if (grad > 0.0) {
        fw = -0.5 * upper_wall_k * grad;
        f += fw;
        if (cvm::debug())
          cvm::log("Applying an upper wall force("+
                    cvm::to_str(fw)+") to \""+this->name+"\".\n");
      }
    }
  }


  if (is_enabled(f_cv_fdiff_velocity)) {
    // set it for the next step
    x_old = x;
  }

  if (cvm::debug())
    cvm::log("Done updating colvar \""+this->name+"\".\n");
  return (potential_energy + kinetic_energy);
}


void colvar::communicate_forces()
{
  size_t i;
  if (cvm::debug())
    cvm::log("Communicating forces from colvar \""+this->name+"\".\n");

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
  } else if (x.type() == colvarvalue::type_scalar) {

    for (i = 0; i < cvcs.size(); i++) {
      if (!cvcs[i]->is_enabled()) continue;
      (cvcs[i])->apply_force(f * (cvcs[i])->sup_coeff *
                              cvm::real((cvcs[i])->sup_np) *
                              (std::pow((cvcs[i])->value().real_value,
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


int colvar::update_cvc_flags()
{
  // Update the enabled/disabled status of cvcs if necessary
  if (cvc_flags.size()) {
    n_active_cvcs = 0;
    active_cvc_square_norm = 0.;

    for (size_t i = 0; i < cvcs.size(); i++) {
      cvcs[i]->feature_states[f_cvc_active]->enabled = cvc_flags[i];
      if (cvcs[i]->is_enabled()) {
        n_active_cvcs++;
        active_cvc_square_norm += cvcs[i]->sup_coeff * cvcs[i]->sup_coeff;
      }
    }
    if (!n_active_cvcs) {
      cvm::error("ERROR: All CVCs are disabled for colvar " + this->name +"\n");
      return COLVARS_ERROR;
    }
    cvc_flags.resize(0);
  }

  return COLVARS_OK;
}


// ******************** METRIC FUNCTIONS ********************
// Use the metrics defined by \link cvc \endlink objects


bool colvar::periodic_boundaries(colvarvalue const &lb, colvarvalue const &ub) const
{
  if ( (!is_enabled(f_cv_lower_boundary)) || (!is_enabled(f_cv_upper_boundary)) ) {
    cvm::log("Error: checking periodicity for collective variable \""+this->name+"\" "
                    "requires lower and upper boundaries to be defined.\n");
    cvm::set_error_bit(INPUT_ERROR);
  }

  if (period > 0.0) {
    if ( ((std::sqrt(this->dist2(lb, ub))) / this->width)
         < 1.0E-10 ) {
      return true;
    }
  }

  return false;
}

bool colvar::periodic_boundaries() const
{
  if ( (!is_enabled(f_cv_lower_boundary)) || (!is_enabled(f_cv_upper_boundary)) ) {
    cvm::log("Error: checking periodicity for collective variable \""+this->name+"\" "
                    "requires lower and upper boundaries to be defined.\n");
  }

  return periodic_boundaries(lower_boundary, upper_boundary);
}


cvm::real colvar::dist2(colvarvalue const &x1,
                         colvarvalue const &x2) const
{
  if (is_enabled(f_cv_homogeneous)) {
    return (cvcs[0])->dist2(x1, x2);
  } else {
    return x1.dist2(x2);
  }
}

colvarvalue colvar::dist2_lgrad(colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  if (is_enabled(f_cv_homogeneous)) {
    return (cvcs[0])->dist2_lgrad(x1, x2);
  } else {
    return x1.dist2_grad(x2);
  }
}

colvarvalue colvar::dist2_rgrad(colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  if (is_enabled(f_cv_homogeneous)) {
    return (cvcs[0])->dist2_rgrad(x1, x2);
  } else {
    return x2.dist2_grad(x1);
  }
}

void colvar::wrap(colvarvalue &x) const
{
  if (is_enabled(f_cv_homogeneous)) {
    (cvcs[0])->wrap(x);
  }
  return;
}


// ******************** INPUT FUNCTIONS ********************

std::istream & colvar::read_restart(std::istream &is)
{
  size_t const start_pos = is.tellg();

  std::string conf;
  if ( !(is >> colvarparse::read_block("colvar", conf)) ) {
    // this is not a colvar block
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  {
    std::string check_name = "";
    if ( (get_keyval(conf, "name", check_name,
                     std::string(""), colvarparse::parse_silent)) &&
         (check_name != name) )  {
      cvm::error("Error: the state file does not match the "
                 "configuration file, at colvar \""+name+"\".\n");
    }
    if (check_name.size() == 0) {
      cvm::error("Error: Collective variable in the "
                 "restart file without any identifier.\n");
    }
  }

  if ( !(get_keyval(conf, "x", x,
                    colvarvalue(x.type()), colvarparse::parse_silent)) ) {
    cvm::log("Error: restart file does not contain "
             "the value of the colvar \""+
             name+"\" .\n");
  } else {
    cvm::log("Restarting collective variable \""+name+"\" from value: "+
             cvm::to_str(x)+"\n");
  }

  if (is_enabled(f_cv_extended_Lagrangian)) {

    if ( !(get_keyval(conf, "extended_x", xr,
                      colvarvalue(x.type()), colvarparse::parse_silent)) &&
         !(get_keyval(conf, "extended_v", vr,
                      colvarvalue(x.type()), colvarparse::parse_silent)) ) {
      cvm::log("Error: restart file does not contain "
               "\"extended_x\" or \"extended_v\" for the colvar \""+
               name+"\", but you requested \"extendedLagrangian\".\n");
    }
    x_reported = xr;
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
      v_reported = vr;
    } else {
      v_reported = v_fdiff;
    }
  }

  return is;
}


std::istream & colvar::read_traj(std::istream &is)
{
  size_t const start_pos = is.tellg();

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
      is >> xr;
      x_reported = xr;
    } else {
      x_reported = x;
    }
  }

  if (is_enabled(f_cv_output_velocity)) {

    is >> v_fdiff;

    if (is_enabled(f_cv_extended_Lagrangian)) {
      is >> vr;
      v_reported = vr;
    } else {
      v_reported = v_fdiff;
    }
  }

  if (is_enabled(f_cv_output_system_force)) {
    is >> ft;
    ft_reported = ft;
  }

  if (is_enabled(f_cv_output_applied_force)) {
    is >> f;
  }

  return is;
}


// ******************** OUTPUT FUNCTIONS ********************

std::ostream & colvar::write_restart(std::ostream &os) {

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
       << xr << "\n"
       << "  extended_v "
       << std::setprecision(cvm::cv_prec)
       << std::setw(cvm::cv_width)
       << vr << "\n";
  }

  os << "}\n\n";

  return os;
}


std::ostream & colvar::write_traj_label(std::ostream & os)
{
  size_t const this_cv_width = x.output_width(cvm::cv_width);

  os << " ";

  if (is_enabled(f_cv_output_value)) {

    os << " "
       << cvm::wrap_string(this->name, this_cv_width);

    if (is_enabled(f_cv_extended_Lagrangian)) {
      // extended DOF
      os << " r_"
         << cvm::wrap_string(this->name, this_cv_width-2);
    }
  }

  if (is_enabled(f_cv_output_velocity)) {

    os << " v_"
       << cvm::wrap_string(this->name, this_cv_width-2);

    if (is_enabled(f_cv_extended_Lagrangian)) {
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

  if (is_enabled(f_cv_output_system_force)) {
    os << " fs_"
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

    if (is_enabled(f_cv_extended_Lagrangian)) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << x;
    }

    os << " "
       << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
       << x_reported;
  }

  if (is_enabled(f_cv_output_velocity)) {

    if (is_enabled(f_cv_extended_Lagrangian)) {
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


  if (is_enabled(f_cv_output_system_force)) {
    os << " "
       << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
       << ft_reported;
  }

  if (is_enabled(f_cv_output_applied_force)) {
    if (is_enabled(f_cv_extended_Lagrangian)) {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << fr;
    } else {
      os << " "
         << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
         << f;
    }
  }

  return os;
}

int colvar::write_output_files()
{
  if (cvm::b_analysis) {

    if (acf.size()) {
      cvm::log("Writing acf to file \""+acf_outfile+"\".\n");

      cvm::backup_file(acf_outfile.c_str());
      cvm::ofstream acf_os(acf_outfile.c_str());
      if (! acf_os.is_open()) {
        cvm::error("Cannot open file \""+acf_outfile+"\".\n", FILE_ERROR);
      }
      write_acf(acf_os);
      acf_os.close();
    }

    if (runave_os.is_open()) {
      runave_os.close();
    }
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}



// ******************** ANALYSIS FUNCTIONS ********************

void colvar::analyze()
{
  if (is_enabled(f_cv_runave)) {
    calc_runave();
  }

  if (is_enabled(f_cv_corrfunc)) {
    calc_acf();
  }
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
  // coordinates(acf_x_history) or velocities (acf_v_history); each vector can
  // contain up to acf_length values, which are contiguous in memory
  // representation but separated by acf_stride in the time series;
  // the pointer to each vector is changed at every step

  if (acf_x_history.empty() && acf_v_history.empty()) {

    // first-step operations

    colvar *cfcv = (acf_colvar_name.size() ?
                    cvm::colvar_by_name(acf_colvar_name) :
                    this);
    if (colvarvalue::check_types(cfcv->value(), value())) {
      cvm::error("Error: correlation function between \""+cfcv->name+
                 "\" and \""+this->name+"\" cannot be calculated, "
                 "because their value types are different.\n",
                 INPUT_ERROR);
    }
    acf_nframes = 0;

    cvm::log("Colvar \""+this->name+"\": initializing ACF calculation.\n");

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

    default:
      break;
    }

  } else {

    colvar *cfcv = (acf_colvar_name.size() ?
                    cvm::colvar_by_name(acf_colvar_name) :
                    this);

    switch (acf_type) {

    case acf_vel:

      if (is_enabled(f_cv_fdiff_velocity)) {
        // calc() should do this already, but this only happens in a
        // simulation; better do it again in case a trajectory is
        // being read
        v_reported = v_fdiff = fdiff_velocity(x_old, cfcv->value());
      }

      calc_vel_acf((*acf_v_history_p), cfcv->velocity());
      // store this value in the history
      history_add_value(acf_length+acf_offset, *acf_v_history_p, cfcv->velocity());
      // if stride is larger than one, cycle among different histories
      history_incr(acf_v_history, acf_v_history_p);
      break;

    case acf_coor:

      calc_coor_acf((*acf_x_history_p), cfcv->value());
      history_add_value(acf_length+acf_offset, *acf_x_history_p, cfcv->value());
      history_incr(acf_x_history, acf_x_history_p);
      break;

    case acf_p2coor:

      calc_p2coor_acf((*acf_x_history_p), cfcv->value());
      history_add_value(acf_length+acf_offset, *acf_x_history_p, cfcv->value());
      history_incr(acf_x_history, acf_x_history_p);
      break;

    default:
      break;
    }
  }

  if (is_enabled(f_cv_fdiff_velocity)) {
    // set it for the next step
    x_old = x;
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvar::calc_vel_acf(std::list<colvarvalue> &v_list,
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
    colvarvalue::inner_opt(v, vs_i, v_list.end(), acf_i);

    acf_nframes++;
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


void colvar::calc_coor_acf(std::list<colvarvalue> &x_list,
                            colvarvalue const      &x)
{
  // same as above but for coordinates
  if (x_list.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  xs_i = x_list.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      ++xs_i;

    *(acf_i++) += x.norm2();

    colvarvalue::inner_opt(x, xs_i, x_list.end(), acf_i);

    acf_nframes++;
  }
}


void colvar::calc_p2coor_acf(std::list<colvarvalue> &x_list,
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

    colvarvalue::p2leg_opt(x, xs_i, x_list.end(), acf_i);

    acf_nframes++;
  }
}


void colvar::write_acf(std::ostream &os)
{
  if (!acf_nframes)
    cvm::log("Warning: ACF was not calculated (insufficient frames).\n");
  os.setf(std::ios::scientific, std::ios::floatfield);
  os << "# Autocorrelation function for collective variable \""
     << this->name << "\"\n";
  // one frame is used for normalization, the statistical sample is
  // hence decreased
  os << "# nframes = " << (acf_normalize ?
                           acf_nframes - 1 :
                           acf_nframes) << "\n";

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
}


void colvar::calc_runave()
{
  if (x_history.empty()) {

    runave.type(value().type());
    runave.reset();

    // first-step operations

    if (cvm::debug())
      cvm::log("Colvar \""+this->name+
                "\": initializing running average calculation.\n");

    acf_nframes = 0;

    x_history.push_back(std::list<colvarvalue>());
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
        runave *= 1.0 / cvm::real(runave_length);
        runave.apply_constraints();

        runave_variance = 0.0;
        runave_variance += this->dist2(x, runave);
        for (xs_i = (*x_history_p).begin();
             xs_i != (*x_history_p).end(); ++xs_i) {
          runave_variance += this->dist2(x, (*xs_i));
        }
        runave_variance *= 1.0 / cvm::real(runave_length-1);

        runave_os << std::setw(cvm::it_width) << cvm::step_relative()
                  << "  "
                  << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
                  << runave << " "
                  << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
                  << std::sqrt(runave_variance) << "\n";
      }

      history_add_value(runave_length, *x_history_p, x);
    }
  }

}

// Static members

std::vector<cvm::deps::feature *> colvar::cv_features;

// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvardeps.h"


colvardeps::colvardeps()
{
  time_step_factor = 1;
}


colvardeps::~colvardeps() {
  size_t i;

  // Protest if we are deleting an object while a parent object may still depend on it
  if (parents.size()) {
    cvm::log("Warning: destroying \"" + description + "\" before its parents objects:");
    for (i=0; i<parents.size(); i++) {
      cvm::log(parents[i]->description);
    }
  }

  // Do not delete features if it's a static object
  // may change in the future though
//     for (i=0; i<features.size(); i++) {
//       if (features[i] != NULL) delete features[i];
//     }

  remove_all_children();
}


void colvardeps::free_children_deps() {
  // Dereference children requirements of all enabled features
  // Useful when object is destroyed or set inactive
  // CAUTION: when setting the parent object inactive, disable "active" first
  // then call this, to avoid double-dereferencing the deps of "active"

  // Cannot be in the base class destructor because it needs the derived class features()
  size_t i,j,fid;

  if (cvm::debug()) cvm::log("DEPS: freeing children deps for " + description);

  cvm::increase_depth();
  for (fid = 0; fid < feature_states.size(); fid++) {
    if (is_enabled(fid)) {
      for (i=0; i<features()[fid]->requires_children.size(); i++) {
        int g = features()[fid]->requires_children[i];
        for (j=0; j<children.size(); j++) {
          if (cvm::debug()) cvm::log("DEPS: dereferencing children's "
            + children[j]->features()[g]->description);
          children[j]->decr_ref_count(g);
        }
      }
    }
  }
  cvm::decrease_depth();
}


// re-enable children features (and increase ref count accordingly)
// So free_children_deps() can be called whenever an object becomes inactive
void colvardeps::restore_children_deps() {
  size_t i,j,fid;

  cvm::increase_depth();
  for (fid = 0; fid < feature_states.size(); fid++) {
    if (is_enabled(fid)) {
      for (i=0; i<features()[fid]->requires_children.size(); i++) {
        int g = features()[fid]->requires_children[i];
        for (j=0; j<children.size(); j++) {
          if (cvm::debug()) cvm::log("DEPS: re-enabling children's "
            + children[j]->features()[g]->description);
          children[j]->enable(g, false, false);
        }
      }
    }
  }
  cvm::decrease_depth();
}


void colvardeps::provide(int feature_id, bool truefalse) {
  feature_states[feature_id].available = truefalse;
}


void colvardeps::set_enabled(int feature_id, bool truefalse) {
  if (truefalse) {
    enable(feature_id);
  } else {
    disable(feature_id);
  }
}


bool colvardeps::get_keyval_feature(colvarparse *cvp,
                                    std::string const &conf, char const *key,
                                    int feature_id, bool const &def_value,
                                    colvarparse::Parse_Mode const parse_mode)
{
  if (!is_user(feature_id)) {
    cvm::error("Cannot set feature \"" + features()[feature_id]->description + "\" from user input in \"" + description + "\".\n");
    return false;
  }
  bool value;
  bool const found = cvp->get_keyval(conf, key, value, def_value, parse_mode);
  if (value) enable(feature_id);
  return found;
}


int colvardeps::enable(int feature_id,
                       bool dry_run /* default: false */,
                       // dry_run: fail silently, do not enable if available
                       // flag is passed recursively to deps of this feature
                       bool toplevel /* default: true */)
// toplevel: false if this is called as part of a chain of dependency resolution
// this is used to diagnose failed dependencies by displaying the full stack
// only the toplevel dependency will throw a fatal error
{
  int res;
  size_t i, j;
  bool ok;
  feature *f = features()[feature_id];
  feature_state *fs = &feature_states[feature_id];

  if (cvm::debug()) {
    cvm::log("DEPS: " + description +
      (dry_run ? " testing " : " enabling ") +
      "\"" + f->description +"\"");
  }

  if (fs->enabled) {
    if (!(dry_run || toplevel)) {
      // This is a dependency
      // Prevent disabling this feature as long
      // as requirement is enabled
      fs->ref_count++;
      if (cvm::debug())
        cvm::log("DEPS: bumping ref_count to " + cvm::to_str(fs->ref_count));
    }
    // Do not try to further resolve deps
    return COLVARS_OK;
  }

  std::string feature_type_descr = is_static(feature_id) ? "Static" :
    (is_dynamic(feature_id) ? "Dynamic" : "User-controlled");

  if (!fs->available) {
    if (!dry_run) {
      if (toplevel) {
        cvm::error("Error: " + feature_type_descr + " feature unavailable: \""
          + f->description + "\" in " + description + ".");
      } else {
        cvm::log(feature_type_descr + " feature unavailable: \""
          + f->description + "\" in " + description + ".");
      }
    }
    return COLVARS_ERROR;
  }

  if (!toplevel && !is_dynamic(feature_id)) {
    if (!dry_run) {
      cvm::log(feature_type_descr + " feature \"" + f->description
        + "\" may not be enabled as a dependency in " + description + ".\n");
    }
    return COLVARS_ERROR;
  }

  // 1) enforce exclusions
  // reminder: exclusions must be mutual for this to work
  for (i=0; i<f->requires_exclude.size(); i++) {
    feature *g = features()[f->requires_exclude[i]];
    if (cvm::debug())
      cvm::log(f->description + " requires exclude " + g->description);
    if (is_enabled(f->requires_exclude[i])) {
      if (!dry_run) {
        cvm::log("Feature \"" + f->description + "\" is incompatible with \""
        + g->description + "\" in " + description + ".");
        if (toplevel) {
          cvm::error("Error: Failed dependency in " + description + ".");
        }
      }
      return COLVARS_ERROR;
    }
  }

  // 2) solve internal deps (self)
  for (i=0; i<f->requires_self.size(); i++) {
    if (cvm::debug())
      cvm::log(f->description + " requires self " + features()[f->requires_self[i]]->description);
    res = enable(f->requires_self[i], dry_run, false);
    if (res != COLVARS_OK) {
      if (!dry_run) {
        cvm::log("...required by \"" + f->description + "\" in " + description);
        if (toplevel) {
          cvm::error("Error: Failed dependency in " + description + ".");
        }
      }
      return res;
    }
  }

  // 3) solve internal alternate deps
  for (i=0; i<f->requires_alt.size(); i++) {

    // test if one is available; if yes, enable and exit w/ success
    ok = false;
    for (j=0; j<f->requires_alt[i].size(); j++) {
      int g = f->requires_alt[i][j];
      if (cvm::debug())
        cvm::log(f->description + " requires alt " + features()[g]->description);
      res = enable(g, true, false);  // see if available
      if (res == COLVARS_OK) {
        ok = true;
        if (!dry_run) {
          enable(g, false, false); // Require again, for real
          fs->alternate_refs.push_back(g); // We remember we enabled this
          // so we can free it if this feature gets disabled
        }
        break;
      }
    }
    if (!ok) {
      if (!dry_run) {
        cvm::log("\"" + f->description + "\" in " + description
          + " requires one of the following features, none of which can be enabled:\n");
        cvm::log("-----------------------------------------\n");
        cvm::increase_depth();
        for (j=0; j<f->requires_alt[i].size(); j++) {
          int g = f->requires_alt[i][j];
          cvm::log(cvm::to_str(j+1) + ". " + features()[g]->description);
          enable(g, false, false); // Just for printing error output
        }
        cvm::decrease_depth();
        cvm::log("-----------------------------------------");
        if (toplevel) {
          cvm::error("Error: Failed dependency in " + description + ".");
        }
      }
      return COLVARS_ERROR;
    }
  }

  // 4) solve deps in children
  // if the object is inactive, we solve but do not enable: will be enabled
  // when the object becomes active
  cvm::increase_depth();
  for (i=0; i<f->requires_children.size(); i++) {
    int g = f->requires_children[i];
    for (j=0; j<children.size(); j++) {
      res = children[j]->enable(g, dry_run || !is_enabled(), false);
      if (res != COLVARS_OK) {
        if (!dry_run) {
          cvm::log("...required by \"" + f->description + "\" in " + description);
          if (toplevel) {
            cvm::error("Error: Failed dependency in " + description + ".");
          }
        }
        return res;
      }
    }
  }
  cvm::decrease_depth();

  // Actually enable feature only once everything checks out
  if (!dry_run) {
    fs->enabled = true;
    // This should be the only reference
    if (!toplevel) fs->ref_count = 1;
    if (feature_id == 0) {
      // Waking up this object, enable all deps in children
      restore_children_deps();
    }
    do_feature_side_effects(feature_id);
    if (cvm::debug())
      cvm::log("DEPS: feature \"" + f->description + "\" in "
        + description + " enabled, ref_count = 1.");
  }
  return COLVARS_OK;
}


int colvardeps::disable(int feature_id) {
  size_t i, j;
  feature *f = features()[feature_id];
  feature_state *fs = &feature_states[feature_id];

  if (cvm::debug()) cvm::log("DEPS: disabling feature \""
      + f->description + "\" in " + description);

  if (fs->enabled == false) {
    return COLVARS_OK;
  }

  if (fs->ref_count > 1) {
    cvm::error("Error: cannot disable feature \"" + f->description
     + "\" in " + description + " because of " + cvm::to_str(fs->ref_count-1)
     + " remaining references.\n" );
    return COLVARS_ERROR;
  }

  // internal deps (self)
  for (i=0; i<f->requires_self.size(); i++) {
    if (cvm::debug()) cvm::log("DEPS: dereferencing self "
      + features()[f->requires_self[i]]->description);
    decr_ref_count(f->requires_self[i]);
  }

  // alternates
  for (i=0; i<fs->alternate_refs.size(); i++) {
    if (cvm::debug()) cvm::log("DEPS: dereferencing alt "
      + features()[fs->alternate_refs[i]]->description);
    decr_ref_count(fs->alternate_refs[i]);
  }
  // Forget these, now that they are dereferenced
  fs->alternate_refs.clear();

  // deps in children
  // except if the object is inactive, then children dependencies
  // have already been dereferenced by this function
  // (or never referenced if feature was enabled while the object
  // was inactive)
  if (is_enabled()) {
    cvm::increase_depth();
    for (i=0; i<f->requires_children.size(); i++) {
      int g = f->requires_children[i];
      for (j=0; j<children.size(); j++) {
        if (cvm::debug()) cvm::log("DEPS: dereferencing children's "
          + children[j]->features()[g]->description);
        children[j]->decr_ref_count(g);
      }
    }
    cvm::decrease_depth();
  }

  fs->enabled = false;
  fs->ref_count = 0;
  if (feature_id == 0) {
    // Putting this object to sleep
    free_children_deps();
  }
  return COLVARS_OK;
}

int colvardeps::decr_ref_count(int feature_id) {
  int &rc = feature_states[feature_id].ref_count;
  feature *f = features()[feature_id];

  if (cvm::debug())
      cvm::log("DEPS: decreasing reference count of \"" + f->description
        + "\" in " + description + ".\n");

  if (rc <= 0) {
    cvm::error("Error: cannot decrease reference count of feature \"" + f->description
      +  "\" in " + description + ", which is " + cvm::to_str(rc) + ".\n");
    return COLVARS_ERROR;
  }

  rc--;
  if (rc == 0 && f->is_dynamic()) {
    // we can auto-disable this feature
    if (cvm::debug())
      cvm::log("DEPS will now auto-disable dynamic feature \"" + f->description
     + "\" in " + description + ".\n");
    disable(feature_id);
  }
  return COLVARS_OK;
}

void colvardeps::init_feature(int feature_id, const char *description, feature_type type) {
  modify_features()[feature_id]->description = description;
  modify_features()[feature_id]->type = type;
}

// Shorthand macros for describing dependencies
#define f_req_self(f, g) features()[f]->requires_self.push_back(g)
// This macro ensures that exclusions are symmetric
#define f_req_exclude(f, g) features()[f]->requires_exclude.push_back(g); \
                            features()[g]->requires_exclude.push_back(f)
#define f_req_children(f, g) features()[f]->requires_children.push_back(g)
#define f_req_alt2(f, g, h) features()[f]->requires_alt.push_back(std::vector<int>(2));\
  features()[f]->requires_alt.back()[0] = g;                                           \
  features()[f]->requires_alt.back()[1] = h
#define f_req_alt3(f, g, h, i) features()[f]->requires_alt.push_back(std::vector<int>(3));\
  features()[f]->requires_alt.back()[0] = g;                                           \
  features()[f]->requires_alt.back()[1] = h;                                           \
  features()[f]->requires_alt.back()[2] = i
#define f_req_alt4(f, g, h, i, j) features()[f]->requires_alt.push_back(std::vector<int>(4));\
  features()[f]->requires_alt.back()[0] = g;                                           \
  features()[f]->requires_alt.back()[1] = h;                                           \
  features()[f]->requires_alt.back()[2] = i;                                           \
  features()[f]->requires_alt.back()[3] = j

void colvardeps::init_cvb_requires() {
  int i;
  if (features().size() == 0) {
    for (i = 0; i < f_cvb_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_cvb_active, "active", f_type_dynamic);
    f_req_children(f_cvb_active, f_cv_active);

    init_feature(f_cvb_awake, "awake", f_type_static);
    f_req_self(f_cvb_awake, f_cvb_active);

    init_feature(f_cvb_apply_force, "apply force", f_type_user);
    f_req_children(f_cvb_apply_force, f_cv_gradient);

    init_feature(f_cvb_get_total_force, "obtain total force", f_type_dynamic);
    f_req_children(f_cvb_get_total_force, f_cv_total_force);

    init_feature(f_cvb_output_acc_work, "output accumulated work", f_type_user);
    f_req_self(f_cvb_output_acc_work, f_cvb_apply_force);

    init_feature(f_cvb_history_dependent, "history-dependent", f_type_static);

    init_feature(f_cvb_time_dependent, "time-dependent", f_type_static);

    init_feature(f_cvb_scalar_variables, "require scalar variables", f_type_static);
    f_req_children(f_cvb_scalar_variables, f_cv_scalar);

    init_feature(f_cvb_calc_pmf, "calculate a PMF", f_type_static);

    init_feature(f_cvb_calc_ti_samples, "calculate TI samples", f_type_dynamic);
    f_req_self(f_cvb_calc_ti_samples, f_cvb_get_total_force);
    f_req_children(f_cvb_calc_ti_samples, f_cv_grid);

    init_feature(f_cvb_write_ti_samples, "write TI samples ", f_type_user);
    f_req_self(f_cvb_write_ti_samples, f_cvb_calc_ti_samples);

    init_feature(f_cvb_write_ti_pmf, "write TI PMF", f_type_user);
    f_req_self(f_cvb_write_ti_pmf, f_cvb_calc_ti_samples);
  }

  // Initialize feature_states for each instance
  feature_states.reserve(f_cvb_ntot);
  for (i = 0; i < f_cvb_ntot; i++) {
    feature_states.push_back(feature_state(true, false));
    // Most features are available, so we set them so
    // and list exceptions below
  }

  // only compute TI samples when deriving from colvarbias_ti
  feature_states[f_cvb_calc_ti_samples].available = false;
}


void colvardeps::init_cv_requires() {
  size_t i;
  if (features().size() == 0) {
    for (i = 0; i < f_cv_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_cv_active, "active", f_type_dynamic);
    // Do not require f_cvc_active in children, as some components may be disabled
    // Colvars must be either a linear combination, or scalar (and polynomial) or scripted/custom
    f_req_alt4(f_cv_active, f_cv_scalar, f_cv_linear, f_cv_scripted, f_cv_custom_function);

    init_feature(f_cv_awake, "awake", f_type_static);
    f_req_self(f_cv_awake, f_cv_active);

    init_feature(f_cv_gradient, "gradient", f_type_dynamic);
    f_req_children(f_cv_gradient, f_cvc_gradient);

    init_feature(f_cv_collect_gradient, "collect gradient", f_type_dynamic);
    f_req_self(f_cv_collect_gradient, f_cv_gradient);
    f_req_self(f_cv_collect_gradient, f_cv_scalar);
    // The following exlusion could be lifted by implementing the feature
    f_req_exclude(f_cv_collect_gradient, f_cv_scripted);

    init_feature(f_cv_fdiff_velocity, "velocity from finite differences", f_type_dynamic);

    // System force: either trivial (spring force); through extended Lagrangian, or calculated explicitly
    init_feature(f_cv_total_force, "total force", f_type_dynamic);
    f_req_alt2(f_cv_total_force, f_cv_extended_Lagrangian, f_cv_total_force_calc);

    // Deps for explicit total force calculation
    init_feature(f_cv_total_force_calc, "total force calculation", f_type_dynamic);
    f_req_self(f_cv_total_force_calc, f_cv_scalar);
    f_req_self(f_cv_total_force_calc, f_cv_linear);
    f_req_children(f_cv_total_force_calc, f_cvc_inv_gradient);
    f_req_self(f_cv_total_force_calc, f_cv_Jacobian);

    init_feature(f_cv_Jacobian, "Jacobian derivative", f_type_dynamic);
    f_req_self(f_cv_Jacobian, f_cv_scalar);
    f_req_self(f_cv_Jacobian, f_cv_linear);
    f_req_children(f_cv_Jacobian, f_cvc_Jacobian);

    init_feature(f_cv_hide_Jacobian, "hide Jacobian force", f_type_user);
    f_req_self(f_cv_hide_Jacobian, f_cv_Jacobian); // can only hide if calculated

    init_feature(f_cv_extended_Lagrangian, "extended Lagrangian", f_type_user);
    f_req_self(f_cv_extended_Lagrangian, f_cv_scalar);
    f_req_self(f_cv_extended_Lagrangian, f_cv_gradient);

    init_feature(f_cv_Langevin, "Langevin dynamics", f_type_user);
    f_req_self(f_cv_Langevin, f_cv_extended_Lagrangian);

    init_feature(f_cv_linear, "linear", f_type_static);

    init_feature(f_cv_scalar, "scalar", f_type_static);

    init_feature(f_cv_output_energy, "output energy", f_type_user);

    init_feature(f_cv_output_value, "output value", f_type_user);

    init_feature(f_cv_output_velocity, "output velocity", f_type_user);
    f_req_self(f_cv_output_velocity, f_cv_fdiff_velocity);

    init_feature(f_cv_output_applied_force, "output applied force", f_type_user);

    init_feature(f_cv_output_total_force, "output total force", f_type_user);
    f_req_self(f_cv_output_total_force, f_cv_total_force);

    init_feature(f_cv_subtract_applied_force, "subtract applied force from total force", f_type_user);
    f_req_self(f_cv_subtract_applied_force, f_cv_total_force);

    init_feature(f_cv_lower_boundary, "lower boundary", f_type_user);
    f_req_self(f_cv_lower_boundary, f_cv_scalar);

    init_feature(f_cv_upper_boundary, "upper boundary", f_type_user);
    f_req_self(f_cv_upper_boundary, f_cv_scalar);

    init_feature(f_cv_grid, "grid", f_type_dynamic);
    f_req_self(f_cv_grid, f_cv_lower_boundary);
    f_req_self(f_cv_grid, f_cv_upper_boundary);

    init_feature(f_cv_runave, "running average", f_type_user);

    init_feature(f_cv_corrfunc, "correlation function", f_type_user);

    init_feature(f_cv_scripted, "scripted", f_type_user);

    init_feature(f_cv_custom_function, "custom function", f_type_user);
    f_req_exclude(f_cv_custom_function, f_cv_scripted);

    init_feature(f_cv_periodic, "periodic", f_type_static);
    f_req_self(f_cv_periodic, f_cv_scalar);
    init_feature(f_cv_scalar, "scalar", f_type_static);
    init_feature(f_cv_linear, "linear", f_type_static);
    init_feature(f_cv_homogeneous, "homogeneous", f_type_static);

    // because total forces are obtained from the previous time step,
    // we cannot (currently) have colvar values and total forces for the same timestep
    init_feature(f_cv_multiple_ts, "multiple timestep colvar");
    f_req_exclude(f_cv_multiple_ts, f_cv_total_force_calc);
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
}


void colvardeps::init_cvc_requires() {
  size_t i;
  // Initialize static array once and for all
  if (features().size() == 0) {
    for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_cvc_active, "active", f_type_dynamic);
//     The dependency below may become useful if we use dynamic atom groups
//     f_req_children(f_cvc_active, f_ag_active);

    init_feature(f_cvc_scalar, "scalar", f_type_static);

    init_feature(f_cvc_gradient, "gradient", f_type_dynamic);

    init_feature(f_cvc_implicit_gradient, "implicit gradient", f_type_static);
    f_req_children(f_cvc_implicit_gradient, f_ag_implicit_gradient);

    init_feature(f_cvc_inv_gradient, "inverse gradient", f_type_dynamic);
    f_req_self(f_cvc_inv_gradient, f_cvc_gradient);

    init_feature(f_cvc_debug_gradient, "debug gradient", f_type_user);
    f_req_self(f_cvc_debug_gradient, f_cvc_gradient);
    f_req_exclude(f_cvc_debug_gradient, f_cvc_implicit_gradient);

    init_feature(f_cvc_Jacobian, "Jacobian derivative", f_type_dynamic);
    f_req_self(f_cvc_Jacobian, f_cvc_inv_gradient);

    init_feature(f_cvc_com_based, "depends on group centers of mass", f_type_static);

    // init_feature(f_cvc_pbc_minimum_image, "use minimum-image distances with PBCs", f_type_user);

    // Compute total force on first site only to avoid unwanted
    // coupling to other colvars (see e.g. Ciccotti et al., 2005)
    init_feature(f_cvc_one_site_total_force, "compute total force from one group", f_type_user);
    f_req_self(f_cvc_one_site_total_force, f_cvc_com_based);

    init_feature(f_cvc_scalable, "scalable calculation", f_type_static);
    f_req_self(f_cvc_scalable, f_cvc_scalable_com);

    init_feature(f_cvc_scalable_com, "scalable calculation of centers of mass", f_type_static);
    f_req_self(f_cvc_scalable_com, f_cvc_com_based);


    // TODO only enable this when f_ag_scalable can be turned on for a pre-initialized group
    // f_req_children(f_cvc_scalable, f_ag_scalable);
    // f_req_children(f_cvc_scalable_com, f_ag_scalable_com);
  }

  // Initialize feature_states for each instance
  // default as available, not enabled
  // except dynamic features which default as unavailable
  feature_states.reserve(f_cvc_ntot);
  for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
    bool avail = is_dynamic(i) ? false : true;
    feature_states.push_back(feature_state(avail, false));
  }

  // CVCs are enabled from the start - get disabled based on flags
  feature_states[f_cvc_active].enabled = true;

  // Features that are implemented by all cvcs by default
  // Each cvc specifies what other features are available
  feature_states[f_cvc_active].available = true;
  feature_states[f_cvc_gradient].available = true;

  // Use minimum-image distances by default
  feature_states[f_cvc_pbc_minimum_image].enabled = true;

  // Features that are implemented by default if their requirements are
  feature_states[f_cvc_one_site_total_force].available = true;

  // Features That are implemented only for certain simulation engine configurations
  feature_states[f_cvc_scalable_com].available = (cvm::proxy->scalable_group_coms() == COLVARS_OK);
  feature_states[f_cvc_scalable].available = feature_states[f_cvc_scalable_com].available;
}


void colvardeps::init_ag_requires() {
  size_t i;
  // Initialize static array once and for all
  if (features().size() == 0) {
    for (i = 0; i < f_ag_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_ag_active, "active", f_type_dynamic);
    init_feature(f_ag_center, "translational fit", f_type_static);
    init_feature(f_ag_rotate, "rotational fit", f_type_static);
    init_feature(f_ag_fitting_group, "reference positions group", f_type_static);
    init_feature(f_ag_implicit_gradient, "implicit atom gradient", f_type_dynamic);
    init_feature(f_ag_fit_gradients, "fit gradients", f_type_user);
    f_req_exclude(f_ag_fit_gradients, f_ag_implicit_gradient);

    init_feature(f_ag_atom_forces, "atomic forces", f_type_dynamic);

    // parallel calculation implies that we have at least a scalable center of mass,
    // but f_ag_scalable is kept as a separate feature to deal with future dependencies
    init_feature(f_ag_scalable, "scalable group calculation", f_type_static);
    init_feature(f_ag_scalable_com, "scalable group center of mass calculation", f_type_static);
    f_req_self(f_ag_scalable, f_ag_scalable_com);

//     init_feature(f_ag_min_msd_fit, "minimum MSD fit")
//     f_req_self(f_ag_min_msd_fit, f_ag_center)
//     f_req_self(f_ag_min_msd_fit, f_ag_rotate)
//     f_req_exclude(f_ag_min_msd_fit, f_ag_fitting_group)
  }

  // Initialize feature_states for each instance
  // default as unavailable, not enabled
  feature_states.reserve(f_ag_ntot);
  for (i = 0; i < colvardeps::f_ag_ntot; i++) {
    feature_states.push_back(feature_state(false, false));
  }

  // Features that are implemented (or not) by all atom groups
  feature_states[f_ag_active].available = true;
  // f_ag_scalable_com is provided by the CVC iff it is COM-based
  feature_states[f_ag_scalable_com].available = false;
  // TODO make f_ag_scalable depend on f_ag_scalable_com (or something else)
  feature_states[f_ag_scalable].available = true;
  feature_states[f_ag_fit_gradients].available = true;
  feature_states[f_ag_implicit_gradient].available = true;
}


void colvardeps::print_state() {
  size_t i;
  cvm::log("Enabled features of \"" + description + "\" (with reference count)");
  for (i = 0; i < feature_states.size(); i++) {
    if (is_enabled(i))
      cvm::log("- " + features()[i]->description + " ("
        + cvm::to_str(feature_states[i].ref_count) + ")");
  }
  cvm::increase_depth();
  for (i=0; i<children.size(); i++) {
    cvm::log("* child " + cvm::to_str(i+1));
    children[i]->print_state();
  }
  cvm::decrease_depth();
}


void colvardeps::add_child(colvardeps *child) {

  children.push_back(child);
  child->parents.push_back((colvardeps *)this);

  // Solve dependencies of already enabled parent features
  // in the new child

  size_t i, fid;
  cvm::increase_depth();
  for (fid = 0; fid < feature_states.size(); fid++) {
    if (is_enabled(fid)) {
      for (i=0; i<features()[fid]->requires_children.size(); i++) {
        int g = features()[fid]->requires_children[i];
        if (cvm::debug()) cvm::log("DEPS: re-enabling children's "
          + child->features()[g]->description);
        child->enable(g, false, false);
      }
    }
  }
  cvm::decrease_depth();
}


void colvardeps::remove_child(colvardeps *child) {
  int i;
  bool found = false;

  for (i = children.size()-1; i>=0; --i) {
    if (children[i] == child) {
      children.erase(children.begin() + i);
      found = true;
      break;
    }
  }
  if (!found) {
    cvm::error("Trying to remove missing child reference from " + description + "\n");
  }
  found = false;
  for (i = child->parents.size()-1; i>=0; --i) {
    if (child->parents[i] == this) {
      child->parents.erase(child->parents.begin() + i);
      found = true;
      break;
    }
  }
  if (!found) {
    cvm::error("Trying to remove missing parent reference from " + child->description + "\n");
  }
}


void colvardeps::remove_all_children() {
  size_t i;
  int j;
  bool found;

  for (i = 0; i < children.size(); ++i) {
    found = false;
    for (j = children[i]->parents.size()-1; j>=0; --j) {
      if (children[i]->parents[j] == this) {
        children[i]->parents.erase(children[i]->parents.begin() + j);
        found = true;
        break;
      }
    }
    if (!found) {
      cvm::error("Trying to remove missing parent reference from " + children[i]->description + "\n");
    }
  }
  children.clear();
}

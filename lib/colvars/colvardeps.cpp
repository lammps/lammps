#include "colvarmodule.h"
#include "colvardeps.h"


cvm::deps::~deps() {
  size_t i;

  for (i=0; i<feature_states.size(); i++) {
    if (feature_states[i] != NULL) delete feature_states[i];
  }
      // Do not delete features if it's static
//     for (i=0; i<features.size(); i++) {
//       if (features[i] != NULL) delete features[i];
//     }
  remove_all_children();

  // Protest if we are deleting an object while a parent object may still depend on it
  // Another possible strategy is to have the child unlist itself from the parent's children
  if (parents.size()) {
    cvm::log("Warning: destroying " + description + " before its parents objects:");
    for (i=0; i<parents.size(); i++) {
      cvm::log(parents[i]->description);
    }
  }
}


void cvm::deps::provide(int feature_id) {
  feature_states[feature_id]->available = true;
}


int cvm::deps::enable(int feature_id,
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
  feature_state *fs = feature_states[feature_id];

  if (cvm::debug()) {
    cvm::log("DEPS: " + description +
      (dry_run ? " testing " : " requiring ") +
      "\"" + f->description);
  }

  if (fs->enabled) {
    // Do not try to solve deps if already enabled
    return COLVARS_OK;
  }

  if (!fs->available) {
    if (!dry_run) {
      if (toplevel) {
        cvm::error("Error: Feature unavailable: \"" + f->description + "\" in " + description + ".");
      } else {
        cvm::log("Feature unavailable: \"" + f->description + "\" in " + description);
      }
    }
    return COLVARS_ERROR;
  }

  // 1) enforce exclusions
  for (i=0; i<f->requires_exclude.size(); i++) {
    feature *g = features()[f->requires_exclude[i]];
    if (cvm::debug())
      cvm::log(f->description + " requires exclude " + g->description);
    if (is_enabled(f->requires_exclude[i])) {
      if (!dry_run) {
        cvm::log("Features \"" + f->description + "\" is incompatible with \""
        + g->description + "\" in " + description);
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
        if (!dry_run) enable(g, false, false); // Require again, for real
        break;
      }
    }
    if (!ok) {
      if (!dry_run) {
        cvm::log("No dependency satisfied among alternates:");
        cvm::log("-----------------------------------------");
        for (j=0; j<f->requires_alt[i].size(); j++) {
          int g = f->requires_alt[i][j];
          cvm::log(cvm::to_str(j+1) + ". " + features()[g]->description);
          cvm::increase_depth();
          enable(g, false, false); // Just for printing error output
          cvm::decrease_depth();
        }
        cvm::log("-----------------------------------------");
        cvm::log("for \"" + f->description + "\" in " + description);
        if (toplevel) {
          cvm::error("Error: Failed dependency in " + description + ".");
        }
      }
      return COLVARS_ERROR;
    }
  }

  // 4) solve deps in children
  for (i=0; i<f->requires_children.size(); i++) {
    int g = f->requires_children[i];
    if (cvm::debug())
      cvm::log("requires children " + features()[g]->description);
//     cvm::log("children " + cvm::to_str(g));
    for (j=0; j<children.size(); j++) {
//       cvm::log("child " +  children[j]->description);
      cvm::increase_depth();
      res = children[j]->enable(g, dry_run, false);
      cvm::decrease_depth();
      if (res != COLVARS_OK) {
        if (!dry_run) {
          cvm::log("...required by \"" + f->description + "\" in " + description);
        }
        if (toplevel) {
          cvm::error("Error: Failed dependency in " + description + ".");
        }
        return res;
      }
    }
    // If we've just touched the features of child objects, refresh them
    if (!dry_run && f->requires_children.size() != 0) {
      for (j=0; j<children.size(); j++) {
        children[j]->refresh_deps();
      }
    }
  }

  // Actually enable feature only once everything checks out
  if (!dry_run) fs->enabled = true;
  return COLVARS_OK;
}


//     disable() {
//
//       // we need refs to parents to walk up the deps tree!
//       // or refresh
//     }

   // Shorthand macros for describing dependencies
#define f_description(f, d) features()[f]->description = d
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

void cvm::deps::init_cvb_requires() {
  int i;
  if (features().size() == 0) {
    for (i = 0; i < f_cvb_ntot; i++) {
      features().push_back(new feature);
    }
  }

  f_description(f_cvb_active, "active");
  f_req_children(f_cvb_active, f_cv_active);

  f_description(f_cvb_apply_force, "apply force");
  f_req_children(f_cvb_apply_force, f_cv_gradient);

  f_description(f_cvb_get_system_force, "obtain system force");
  f_req_children(f_cvb_get_system_force, f_cv_system_force);

  // Initialize feature_states for each instance
  feature_states.reserve(f_cvb_ntot);
  for (i = 0; i < f_cvb_ntot; i++) {
    feature_states.push_back(new feature_state(true, false));
    // Most features are available, so we set them so
    // and list exceptions below
  }
}


void cvm::deps::init_cv_requires() {
  size_t i;
  if (features().size() == 0) {
    for (i = 0; i < f_cv_ntot; i++) {
      features().push_back(new feature);
    }

    f_description(f_cv_active, "active");
    f_req_children(f_cv_active, f_cvc_active);
    // Colvars must be either a linear combination, or scalar (and polynomial) or scripted
    f_req_alt3(f_cv_active, f_cv_scalar, f_cv_linear, f_cv_scripted);

    f_description(f_cv_gradient, "gradient");
    f_req_children(f_cv_gradient, f_cvc_gradient);

    f_description(f_cv_collect_gradient, "collect gradient");
    f_req_self(f_cv_collect_gradient, f_cv_gradient);
    f_req_self(f_cv_collect_gradient, f_cv_scalar);

    f_description(f_cv_fdiff_velocity, "fdiff_velocity");

    // System force: either trivial (spring force); through extended Lagrangian, or calculated explicitly
    f_description(f_cv_system_force, "system force");
    f_req_alt2(f_cv_system_force, f_cv_extended_Lagrangian, f_cv_system_force_calc);

    // Deps for explicit system force calculation
    f_description(f_cv_system_force_calc, "system force calculation");
    f_req_self(f_cv_system_force_calc, f_cv_scalar);
    f_req_self(f_cv_system_force_calc, f_cv_linear);
    f_req_children(f_cv_system_force_calc, f_cvc_inv_gradient);
    f_req_self(f_cv_system_force_calc, f_cv_Jacobian);

    f_description(f_cv_Jacobian, "Jacobian derivative");
    f_req_self(f_cv_Jacobian, f_cv_scalar);
    f_req_self(f_cv_Jacobian, f_cv_linear);
    f_req_children(f_cv_Jacobian, f_cvc_Jacobian);

    f_description(f_cv_hide_Jacobian, "hide Jacobian force");
    f_req_self(f_cv_hide_Jacobian, f_cv_Jacobian); // can only hide if calculated

    f_description(f_cv_extended_Lagrangian, "extended Lagrangian");

    f_description(f_cv_Langevin, "Langevin dynamics");
    f_req_self(f_cv_Langevin, f_cv_extended_Lagrangian);

    f_description(f_cv_linear, "linear");

    f_description(f_cv_scalar, "scalar");

    f_description(f_cv_output_energy, "output energy");

    f_description(f_cv_output_value, "output value");

    f_description(f_cv_output_velocity, "output velocity");
    f_req_self(f_cv_output_velocity, f_cv_fdiff_velocity);

    f_description(f_cv_output_applied_force, "output applied force");

    f_description(f_cv_output_system_force, "output system force");
    f_req_self(f_cv_output_system_force, f_cv_system_force);

    f_description(f_cv_lower_boundary, "lower boundary");
    f_req_self(f_cv_lower_boundary, f_cv_scalar);

    f_description(f_cv_upper_boundary, "upper boundary");
    f_req_self(f_cv_upper_boundary, f_cv_scalar);

    f_description(f_cv_grid, "grid");
    f_req_self(f_cv_grid, f_cv_lower_boundary);
    f_req_self(f_cv_grid, f_cv_upper_boundary);

    f_description(f_cv_lower_wall, "lower wall");
    f_req_self(f_cv_lower_wall, f_cv_lower_boundary);
    f_req_self(f_cv_lower_wall, f_cv_gradient);

    f_description(f_cv_upper_wall, "upper wall");
    f_req_self(f_cv_upper_wall, f_cv_upper_boundary);
    f_req_self(f_cv_upper_wall, f_cv_gradient);

    f_description(f_cv_runave, "running average");

    f_description(f_cv_corrfunc, "correlation function");

    // The features below are set programmatically
    f_description(f_cv_scripted, "scripted");
    f_description(f_cv_periodic, "periodic");
    f_description(f_cv_scalar, "scalar");
    f_description(f_cv_linear, "linear");
    f_description(f_cv_homogeneous, "homogeneous");
  }

  // Initialize feature_states for each instance
  feature_states.reserve(f_cv_ntot);
  for (i = 0; i < f_cv_ntot; i++) {
    feature_states.push_back(new feature_state(true, false));
    // Most features are available, so we set them so
    // and list exceptions below
   }

  // properties that may NOT be enabled as a dependency
  int unavailable_deps[] = {
    f_cv_lower_boundary,
    f_cv_upper_boundary,
    f_cv_extended_Lagrangian,
    f_cv_Langevin,
    f_cv_scripted,
    f_cv_periodic,
    f_cv_scalar,
    f_cv_linear,
    f_cv_homogeneous
  };
  for (i = 0; i < sizeof(unavailable_deps) / sizeof(unavailable_deps[0]); i++) {
    feature_states[unavailable_deps[i]]->available = false;
  }
}


void cvm::deps::init_cvc_requires() {
  size_t i;
  // Initialize static array once and for all
  if (features().size() == 0) {
    for (i = 0; i < cvm::deps::f_cvc_ntot; i++) {
      features().push_back(new feature);
    }

    f_description(f_cvc_active, "active");
//     The dependency below may become useful if we use dynamic atom groups
//     f_req_children(f_cvc_active, f_ag_active);

    f_description(f_cvc_scalar, "scalar");

    f_description(f_cvc_gradient, "gradient");

    f_description(f_cvc_inv_gradient, "inverse gradient");
    f_req_self(f_cvc_inv_gradient, f_cvc_gradient);

    f_description(f_cvc_debug_gradient, "debug gradient");
    f_req_self(f_cvc_debug_gradient, f_cvc_gradient);

    f_description(f_cvc_Jacobian, "Jacobian");
    f_req_self(f_cvc_Jacobian, f_cvc_inv_gradient);

    f_description(f_cvc_com_based, "depends on group centers of mass");

    f_description(f_cvc_scalable, "scalable calculation");
    f_req_self(f_cvc_scalable, f_cvc_scalable_com);

    f_description(f_cvc_scalable_com, "scalable calculation of centers of mass");
    f_req_self(f_cvc_scalable_com, f_cvc_com_based);


    // TODO only enable this when f_ag_scalable can be turned on for a pre-initialized group
    // f_req_children(f_cvc_scalable, f_ag_scalable);
    // f_req_children(f_cvc_scalable_com, f_ag_scalable_com);
  }

  // Initialize feature_states for each instance
  // default as unavailable, not enabled
  feature_states.reserve(f_cvc_ntot);
  for (i = 0; i < cvm::deps::f_cvc_ntot; i++) {
    feature_states.push_back(new feature_state(false, false));
  }

  // Features that are implemented by all cvcs by default
  // Each cvc specifies what other features are available
  feature_states[f_cvc_active]->available = true;
  feature_states[f_cvc_gradient]->available = true;
  feature_states[f_cvc_scalable_com]->available = (proxy->scalable_group_coms() == COLVARS_OK);
  feature_states[f_cvc_scalable]->available = feature_states[f_cvc_scalable_com]->available;
}


void cvm::deps::init_ag_requires() {
  size_t i;
  // Initialize static array once and for all
  if (features().size() == 0) {
    for (i = 0; i < f_ag_ntot; i++) {
      features().push_back(new feature);
    }

    f_description(f_ag_active, "active");
    f_description(f_ag_center, "translational fit");
    f_description(f_ag_rotate, "rotational fit");
    f_description(f_ag_fitting_group, "reference positions group");
    f_description(f_ag_fit_gradient_group, "fit gradient for main group");
    f_description(f_ag_fit_gradient_ref, "fit gradient for reference group");
    f_description(f_ag_atom_forces, "atomic forces");

    // parallel calculation implies that we have at least a scalable center of mass,
    // but f_ag_scalable is kept as a separate feature to deal with future dependencies
    f_description(f_ag_scalable, "scalable group calculation");
    f_description(f_ag_scalable_com, "scalable group center of mass calculation");
    f_req_self(f_ag_scalable, f_ag_scalable_com);

//     f_description(f_ag_min_msd_fit, "minimum MSD fit")
//     f_req_self(f_ag_min_msd_fit, f_ag_center)
//     f_req_self(f_ag_min_msd_fit, f_ag_rotate)
//     f_req_exclude(f_ag_min_msd_fit, f_ag_fitting_group)
  }

  // Initialize feature_states for each instance
  // default as unavailable, not enabled
  feature_states.reserve(f_ag_ntot);
  for (i = 0; i < cvm::deps::f_ag_ntot; i++) {
    feature_states.push_back(new feature_state(false, false));
  }

  // Features that are implemented (or not) by all atom groups
  feature_states[f_ag_active]->available = true;
  feature_states[f_ag_scalable_com]->available = (proxy->scalable_group_coms() == COLVARS_OK);
  feature_states[f_ag_scalable]->available = feature_states[f_ag_scalable_com]->available;
}


void cvm::deps::print_state() {
  size_t i;
  cvm::log("Enabled features of " + description);
  for (i = 0; i < feature_states.size(); i++) {
    if (feature_states[i]->enabled)
      cvm::log("- " + features()[i]->description);
  }
  for (i=0; i<children.size(); i++) {
    cvm::log("* child " + cvm::to_str(i+1));
    cvm::increase_depth();
    children[i]->print_state();
    cvm::decrease_depth();
  }
}



void cvm::deps::add_child(deps *child) {
  children.push_back(child);
  child->parents.push_back((deps *)this);
}


void cvm::deps::remove_child(deps *child) {
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


void cvm::deps::remove_all_children() {
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

// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
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
                       bool toplevel /* default: true */)
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
      // This is a dependency: prevent disabling this feature as long
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
        + "\" cannot be enabled automatically in " + description + ".");
      if (is_user(feature_id)) {
        cvm::log("Try setting it manually.\n");
      }
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


void colvardeps::init_feature(int feature_id, const char *description_in, feature_type type) {
  modify_features()[feature_id]->description = description_in;
  modify_features()[feature_id]->type = type;
}


// Shorthand functions for describing dependencies
void colvardeps::require_feature_self(int f, int g) {
  features()[f]->requires_self.push_back(g);
}


// Ensure that exclusions are symmetric
void colvardeps::exclude_feature_self(int f, int g) {
  features()[f]->requires_exclude.push_back(g);
  features()[g]->requires_exclude.push_back(f);
}


void colvardeps::require_feature_children(int f, int g) {
  features()[f]->requires_children.push_back(g);
}


void colvardeps::require_feature_alt(int f, int g, int h) {
  features()[f]->requires_alt.push_back(std::vector<int>(2));
  features()[f]->requires_alt.back()[0] = g;
  features()[f]->requires_alt.back()[1] = h;
}


void colvardeps::require_feature_alt(int f, int g, int h, int i) {
  features()[f]->requires_alt.push_back(std::vector<int>(3));
  features()[f]->requires_alt.back()[0] = g;
  features()[f]->requires_alt.back()[1] = h;
  features()[f]->requires_alt.back()[2] = i;
}


void colvardeps::require_feature_alt(int f, int g, int h, int i, int j) {
  features()[f]->requires_alt.push_back(std::vector<int>(4));
  features()[f]->requires_alt.back()[0] = g;
  features()[f]->requires_alt.back()[1] = h;
  features()[f]->requires_alt.back()[2] = i;
  features()[f]->requires_alt.back()[3] = j;
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
  child->parents.push_back(this);

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

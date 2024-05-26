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
#include <sstream>
#include <iomanip>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarparse.h"
#include "colvaratoms.h"


cvm::atom::atom()
{
  index = -1;
  id = -1;
  mass = 1.0;
  charge = 0.0;
  reset_data();
}


cvm::atom::atom(int atom_number)
{
  colvarproxy *p = cvm::main()->proxy;
  index = p->init_atom(atom_number);
  id = p->get_atom_id(index);
  update_mass();
  update_charge();
  reset_data();
}


cvm::atom::atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id)
{
  colvarproxy *p = cvm::main()->proxy;
  index = p->init_atom(residue, atom_name, segment_id);
  id = p->get_atom_id(index);
  update_mass();
  update_charge();
  reset_data();
}


cvm::atom::atom(atom const &a)
  : index(a.index)
{
  colvarproxy *p = cvm::main()->proxy;
  id = p->get_atom_id(index);
  p->increase_refcount(index);
  update_mass();
  update_charge();
  reset_data();
}


cvm::atom::~atom()
{
  if (index >= 0) {
    (cvm::main()->proxy)->clear_atom(index);
  }
}


cvm::atom & cvm::atom::operator = (cvm::atom const &a)
{
  index = a.index;
  id = (cvm::main()->proxy)->get_atom_id(index);
  update_mass();
  update_charge();
  reset_data();
  return *this;
}



cvm::atom_group::atom_group()
{
  init();
}


cvm::atom_group::atom_group(char const *key_in)
{
  key = key_in;
  init();
}


cvm::atom_group::atom_group(std::vector<cvm::atom> const &atoms_in)
{
  init();
  atoms = atoms_in;
  setup();
}


cvm::atom_group::~atom_group()
{
  if (is_enabled(f_ag_scalable) && !b_dummy) {
    (cvm::main()->proxy)->clear_atom_group(index);
    index = -1;
  }

  if (fitting_group) {
    delete fitting_group;
    fitting_group = NULL;
  }

  cvm::main()->unregister_named_atom_group(this);
}


int cvm::atom_group::add_atom(cvm::atom const &a)
{
  if (a.id < 0) {
    return COLVARS_ERROR;
  }

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == a.id) {
      if (cvm::debug())
        cvm::log("Discarding doubly counted atom with number "+
                 cvm::to_str(a.id+1)+".\n");
      return COLVARS_OK;
    }
  }

  // for consistency with add_atom_id(), we update the list as well
  atoms_ids.push_back(a.id);
  atoms.push_back(a);
  total_mass += a.mass;
  total_charge += a.charge;

  return COLVARS_OK;
}


int cvm::atom_group::add_atom_id(int aid)
{
  if (aid < 0) {
    return COLVARS_ERROR;
  }

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      if (cvm::debug())
        cvm::log("Discarding doubly counted atom with number "+
                 cvm::to_str(aid+1)+".\n");
      return COLVARS_OK;
    }
  }

  atoms_ids.push_back(aid);
  return COLVARS_OK;
}


int cvm::atom_group::remove_atom(cvm::atom_iter ai)
{
  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: cannot remove atoms from a scalable group.\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  }

  if (!this->size()) {
    cvm::error("Error: trying to remove an atom from an empty group.\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  } else {
    total_mass -= ai->mass;
    total_charge -= ai->charge;
    atoms_ids.erase(atoms_ids.begin() + (ai - atoms.begin()));
    atoms.erase(ai);
  }

  return COLVARS_OK;
}


int cvm::atom_group::set_dummy()
{
  if (atoms_ids.size() > 0) {
    return cvm::error("Error: setting group with keyword \""+key+
                      "\" and name \""+name+"\" as dummy, but it already "
                      "contains atoms.\n", COLVARS_INPUT_ERROR);
  }
  b_dummy = true;
  return COLVARS_OK;
}


int cvm::atom_group::set_dummy_pos(cvm::atom_pos const &pos)
{
  if (b_dummy) {
    dummy_atom_pos = pos;
  } else {
    return cvm::error("Error: setting dummy position for group with keyword \""+
                      key+"\" and name \""+name+
                      "\", but it is not dummy.\n", COLVARS_INPUT_ERROR);
  }
  return COLVARS_OK;
}


int cvm::atom_group::init()
{
  if (!key.size()) key = "unnamed";
  description = "atom group " + key;
  // These may be overwritten by parse(), if a name is provided

  atoms.clear();
  atom_group::init_dependencies();
  index = -1;

  b_dummy = false;
  b_user_defined_fit = false;
  fitting_group = NULL;

  noforce = false;

  total_mass = 0.0;
  total_charge = 0.0;

  cog.reset();
  com.reset();

  return COLVARS_OK;
}


int cvm::atom_group::init_dependencies() {
  size_t i;
  // Initialize static array once and for all
  if (features().size() == 0) {
    for (i = 0; i < f_ag_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_ag_active, "active", f_type_dynamic);
    init_feature(f_ag_center, "center_to_reference", f_type_user);
    init_feature(f_ag_center_origin, "center_to_origin", f_type_user);
    init_feature(f_ag_rotate, "rotate_to_origin", f_type_user);
    init_feature(f_ag_fitting_group, "fitting_group", f_type_static);
    init_feature(f_ag_explicit_gradient, "explicit_atom_gradient", f_type_dynamic);
    init_feature(f_ag_fit_gradients, "fit_gradients", f_type_user);
    require_feature_self(f_ag_fit_gradients, f_ag_explicit_gradient);

    init_feature(f_ag_atom_forces, "atomic_forces", f_type_dynamic);

    // parallel calculation implies that we have at least a scalable center of mass,
    // but f_ag_scalable is kept as a separate feature to deal with future dependencies
    init_feature(f_ag_scalable, "scalable_group", f_type_dynamic);
    init_feature(f_ag_scalable_com, "scalable_group_center_of_mass", f_type_static);
    require_feature_self(f_ag_scalable_com, f_ag_scalable);

    init_feature(f_ag_collect_atom_ids, "collect_atom_ids", f_type_dynamic);
    exclude_feature_self(f_ag_collect_atom_ids, f_ag_scalable);

    // check that everything is initialized
    for (i = 0; i < colvardeps::f_ag_ntot; i++) {
      if (is_not_set(i)) {
        cvm::error("Uninitialized feature " + cvm::to_str(i) + " in " + description);
      }
    }
  }

  // Initialize feature_states for each instance
  // default as unavailable, not enabled
  feature_states.reserve(f_ag_ntot);
  for (i = 0; i < colvardeps::f_ag_ntot; i++) {
    feature_states.push_back(feature_state(false, false));
  }

  // Features that are implemented (or not) by all atom groups
  feature_states[f_ag_active].available = true;
  feature_states[f_ag_center].available = true;
  feature_states[f_ag_center_origin].available = true;
  feature_states[f_ag_rotate].available = true;

  // f_ag_scalable_com is provided by the CVC iff it is COM-based
  feature_states[f_ag_scalable_com].available = false;
  feature_states[f_ag_scalable].available = true;
  feature_states[f_ag_fit_gradients].available = true;
  feature_states[f_ag_fitting_group].available = true;
  feature_states[f_ag_explicit_gradient].available = true;
  feature_states[f_ag_collect_atom_ids].available = true;

  return COLVARS_OK;
}


int cvm::atom_group::setup()
{
  if (atoms_ids.size() == 0) {
    atoms_ids.reserve(atoms.size());
    for (cvm::atom_iter ai = atoms.begin(); ai != atoms.end(); ai++) {
      atoms_ids.push_back(ai->id);
    }
  }
  for (cvm::atom_iter ai = atoms.begin(); ai != atoms.end(); ai++) {
    ai->update_mass();
    ai->update_charge();
  }
  update_total_mass();
  update_total_charge();
  return COLVARS_OK;
}


void cvm::atom_group::update_total_mass()
{
  if (b_dummy) {
    total_mass = 1.0;
    return;
  }

  if (is_enabled(f_ag_scalable)) {
    total_mass = (cvm::main()->proxy)->get_atom_group_mass(index);
  } else {
    total_mass = 0.0;
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      total_mass += ai->mass;
    }
  }
  if (total_mass < 1e-15) {
    cvm::error("ERROR: " + description + " has zero total mass.\n");
  }
}


void cvm::atom_group::update_total_charge()
{
  if (b_dummy) {
    total_charge = 0.0;
    return;
  }

  if (is_enabled(f_ag_scalable)) {
    total_charge = (cvm::main()->proxy)->get_atom_group_charge(index);
  } else {
    total_charge = 0.0;
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      total_charge += ai->charge;
    }
  }
}


void cvm::atom_group::print_properties(std::string const &colvar_name,
                                       int i, int j)
{
  if (cvm::proxy->updated_masses() && cvm::proxy->updated_charges()) {
    cvm::log("Re-initialized atom group for variable \""+colvar_name+"\":"+
             cvm::to_str(i)+"/"+
             cvm::to_str(j)+". "+ cvm::to_str(atoms_ids.size())+
             " atoms: total mass = "+cvm::to_str(total_mass)+
             ", total charge = "+cvm::to_str(total_charge)+".\n");
  }
}


int cvm::atom_group::parse(std::string const &group_conf)
{
  cvm::log("Initializing atom group \""+key+"\".\n");

  // whether or not to include messages in the log
  // colvarparse::Parse_Mode mode = parse_silent;
  // {
  //   bool b_verbose;
  //   get_keyval (group_conf, "verboseOutput", b_verbose, false, parse_silent);
  //   if (b_verbose) mode = parse_normal;
  // }
  // colvarparse::Parse_Mode mode = parse_normal;

  int parse_error = COLVARS_OK;

  // Optional group name will let other groups reuse atom definition
  if (get_keyval(group_conf, "name", name)) {
    if ((cvm::atom_group_by_name(this->name) != NULL) &&
        (cvm::atom_group_by_name(this->name) != this)) {
      cvm::error("Error: this atom group cannot have the same name, \""+this->name+
                        "\", as another atom group.\n",
                COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
    cvm::main()->register_named_atom_group(this);
    description = "atom group " + name;
  }

  // We need to know about fitting to decide whether the group is scalable
  // and we need to know about scalability before adding atoms
  bool b_defined_center = get_keyval_feature(this, group_conf, "centerToOrigin", f_ag_center_origin, false);
  // Legacy alias
  b_defined_center |= get_keyval_feature(this, group_conf, "centerReference", f_ag_center, is_enabled(f_ag_center_origin), parse_deprecated);
  b_defined_center |= get_keyval_feature(this, group_conf, "centerToReference", f_ag_center, is_enabled(f_ag_center));

  if (is_enabled(f_ag_center_origin) && ! is_enabled(f_ag_center)) {
    return cvm::error("centerToReference may not be disabled if centerToOrigin"
                      "is enabled.\n", COLVARS_INPUT_ERROR);
  }
  // Legacy alias
  bool b_defined_rotate = get_keyval_feature(this, group_conf, "rotateReference", f_ag_rotate, false, parse_deprecated);
  b_defined_rotate |= get_keyval_feature(this, group_conf, "rotateToReference", f_ag_rotate, is_enabled(f_ag_rotate));

  if (is_enabled(f_ag_rotate) || is_enabled(f_ag_center) ||
      is_enabled(f_ag_center_origin)) {
    cvm::main()->cite_feature("Moving frame of reference");
  }

  // is the user setting explicit options?
  b_user_defined_fit = b_defined_center || b_defined_rotate;

  if (is_available(f_ag_scalable_com) && !is_enabled(f_ag_rotate) && !is_enabled(f_ag_center)) {
    enable(f_ag_scalable_com);
  }

  {
    std::string atoms_of = "";
    if (get_keyval(group_conf, "atomsOfGroup", atoms_of)) {
      atom_group * ag = atom_group_by_name(atoms_of);
      if (ag == NULL) {
        cvm::error("Error: cannot find atom group with name " + atoms_of + ".\n");
        return COLVARS_ERROR;
      }
      parse_error |= add_atoms_of_group(ag);
    }
  }

//   if (get_keyval(group_conf, "copyOfGroup", source)) {
//     // Goal: Initialize this as a full copy
//     // for this we'll need an atom_group copy constructor
//     return COLVARS_OK;
//   }

  {
    std::string numbers_conf = "";
    size_t pos = 0;
    while (key_lookup(group_conf, "atomNumbers", &numbers_conf, &pos)) {
      parse_error |= add_atom_numbers(numbers_conf);
      numbers_conf = "";
    }
  }

  {
    std::string index_group_name;
    if (get_keyval(group_conf, "indexGroup", index_group_name)) {
      // use an index group from the index file read globally
      parse_error |= add_index_group(index_group_name);
    }
  }

  {
    std::string range_conf = "";
    size_t pos = 0;
    while (key_lookup(group_conf, "atomNumbersRange",
                      &range_conf, &pos)) {
      parse_error |= add_atom_numbers_range(range_conf);
      range_conf = "";
    }
  }

  {
    std::vector<std::string> psf_segids;
    get_keyval(group_conf, "psfSegID", psf_segids, std::vector<std::string>());
    std::vector<std::string>::iterator psii;
    for (psii = psf_segids.begin(); psii < psf_segids.end(); ++psii) {
      if ( (psii->size() == 0) || (psii->size() > 4) ) {
        cvm::error("Error: invalid PSF segment identifier provided, \""+
                   (*psii)+"\".\n", COLVARS_INPUT_ERROR);
      }
    }

    std::string range_conf = "";
    size_t pos = 0;
    size_t range_count = 0;
    psii = psf_segids.begin();
    while (key_lookup(group_conf, "atomNameResidueRange",
                      &range_conf, &pos)) {
      range_count++;
      if (psf_segids.size() && (range_count > psf_segids.size())) {
        cvm::error("Error: more instances of \"atomNameResidueRange\" than "
                   "values of \"psfSegID\".\n", COLVARS_INPUT_ERROR);
      } else {
        parse_error |= add_atom_name_residue_range(psf_segids.size() ?
          *psii : std::string(""), range_conf);
        if (psf_segids.size()) psii++;
      }
      range_conf = "";
    }
  }

  {
    // read the atoms from a file
    std::string atoms_file_name;
    if (get_keyval(group_conf, "atomsFile", atoms_file_name, std::string(""))) {

      std::string atoms_col;
      if (!get_keyval(group_conf, "atomsCol", atoms_col, std::string(""))) {
        cvm::error("Error: parameter atomsCol is required if atomsFile is set.\n",
                   COLVARS_INPUT_ERROR);
      }

      double atoms_col_value;
      bool const atoms_col_value_defined = get_keyval(group_conf, "atomsColValue", atoms_col_value, 0.0);
      if (atoms_col_value_defined && (!atoms_col_value)) {
        cvm::error("Error: atomsColValue, if provided, must be non-zero.\n", COLVARS_INPUT_ERROR);
      }

      // NOTE: calls to add_atom() and/or add_atom_id() are in the proxy-implemented function
      parse_error |= cvm::load_atoms(atoms_file_name.c_str(), *this, atoms_col, atoms_col_value);
    }
  }

  // Catch any errors from all the initialization steps above
  if (parse_error || cvm::get_error()) return (parse_error || cvm::get_error());

  // checks of doubly-counted atoms have been handled by add_atom() already

  if (get_keyval(group_conf, "dummyAtom", dummy_atom_pos, cvm::atom_pos())) {

    parse_error |= set_dummy();
    parse_error |= set_dummy_pos(dummy_atom_pos);

  } else {

    if (!(atoms_ids.size())) {
      parse_error |= cvm::error("Error: no atoms defined for atom group \""+
                                key+"\".\n", COLVARS_INPUT_ERROR);
    }

    // whether these atoms will ever receive forces or not
    bool enable_forces = true;
    get_keyval(group_conf, "enableForces", enable_forces, true, colvarparse::parse_silent);
    noforce = !enable_forces;
  }

  // Now that atoms are defined we can parse the detailed fitting options
  parse_error |= parse_fitting_options(group_conf);

  if (is_enabled(f_ag_scalable) && !b_dummy) {
    cvm::log("Enabling scalable calculation for group \""+this->key+"\".\n");
    index = (cvm::proxy)->init_atom_group(atoms_ids);
  }

  bool b_print_atom_ids = false;
  get_keyval(group_conf, "printAtomIDs", b_print_atom_ids, false);

  // Calculate all required properties (such as total mass)
  setup();

  if (cvm::debug())
    cvm::log("Done initializing atom group \""+key+"\".\n");

  {
    std::string init_msg;
    init_msg.append("Atom group \""+key+"\" defined with "+
                    cvm::to_str(atoms_ids.size())+" atoms requested");
    if ((cvm::proxy)->updated_masses()) {
      init_msg.append(": total mass = "+
                      cvm::to_str(total_mass));
      if ((cvm::proxy)->updated_charges()) {
        init_msg.append(", total charge = "+
                        cvm::to_str(total_charge));
      }
    }
    init_msg.append(".\n");
    cvm::log(init_msg);
  }

  if (b_print_atom_ids) {
    cvm::log("Internal definition of the atom group:\n");
    cvm::log(print_atom_ids());
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int cvm::atom_group::add_atoms_of_group(atom_group const *ag)
{
  std::vector<int> const &source_ids = ag->atoms_ids;

  if (source_ids.size()) {
    atoms_ids.reserve(atoms_ids.size()+source_ids.size());

    if (is_enabled(f_ag_scalable)) {
      for (size_t i = 0; i < source_ids.size(); i++) {
        add_atom_id(source_ids[i]);
      }
    } else {
      atoms.reserve(atoms.size()+source_ids.size());
      for (size_t i = 0; i < source_ids.size(); i++) {
        // We could use the atom copy constructor, but only if the source
        // group is not scalable - whereas this works in both cases
        // atom constructor expects 1-based atom number
        add_atom(cvm::atom(source_ids[i] + 1));
      }
    }

    if (cvm::get_error()) return COLVARS_ERROR;
  } else {
    cvm::error("Error: source atom group contains no atoms\".\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  }

  return COLVARS_OK;
}


int cvm::atom_group::add_atom_numbers(std::string const &numbers_conf)
{
  std::vector<int> atom_indexes;

  if (numbers_conf.size()) {
    std::istringstream is(numbers_conf);
    int ia;
    while (is >> ia) {
      atom_indexes.push_back(ia);
    }
  }

  if (atom_indexes.size()) {
    atoms_ids.reserve(atoms_ids.size()+atom_indexes.size());

    if (is_enabled(f_ag_scalable)) {
      for (size_t i = 0; i < atom_indexes.size(); i++) {
        add_atom_id((cvm::proxy)->check_atom_id(atom_indexes[i]));
      }
    } else {
      // if we are handling the group on rank 0, better allocate the vector in one shot
      atoms.reserve(atoms.size()+atom_indexes.size());
      for (size_t i = 0; i < atom_indexes.size(); i++) {
        add_atom(cvm::atom(atom_indexes[i]));
      }
    }

    if (cvm::get_error()) return COLVARS_ERROR;
  } else {
    cvm::error("Error: no numbers provided for \""
               "atomNumbers\".\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  }

  return COLVARS_OK;
}


int cvm::atom_group::add_index_group(std::string const &index_group_name)
{
  std::vector<std::string> const &index_group_names =
    cvm::main()->index_group_names;
  std::vector<std::vector<int> *> const &index_groups =
    cvm::main()->index_groups;

  size_t i_group = 0;
  for ( ; i_group < index_groups.size(); i_group++) {
    if (index_group_names[i_group] == index_group_name)
      break;
  }

  if (i_group >= index_group_names.size()) {
    return cvm::error("Error: could not find index group "+
                      index_group_name+" among those already provided.\n",
                      COLVARS_INPUT_ERROR);
  }

  int error_code = COLVARS_OK;

  std::vector<int> const &index_group = *(index_groups[i_group]);

  atoms_ids.reserve(atoms_ids.size()+index_group.size());

  if (is_enabled(f_ag_scalable)) {
    for (size_t i = 0; i < index_group.size(); i++) {
      error_code |=
        add_atom_id((cvm::proxy)->check_atom_id(index_group[i]));
    }
  } else {
    atoms.reserve(atoms.size()+index_group.size());
    for (size_t i = 0; i < index_group.size(); i++) {
      error_code |= add_atom(cvm::atom(index_group[i]));
    }
  }

  return error_code;
}


int cvm::atom_group::add_atom_numbers_range(std::string const &range_conf)
{
  if (range_conf.size()) {
    std::istringstream is(range_conf);
    int initial, final;
    char dash;
    if ( (is >> initial) && (initial > 0) &&
         (is >> dash) && (dash == '-') &&
         (is >> final) && (final > 0) ) {

      atoms_ids.reserve(atoms_ids.size() + (final - initial + 1));

      if (is_enabled(f_ag_scalable)) {
        for (int anum = initial; anum <= final; anum++) {
          add_atom_id((cvm::proxy)->check_atom_id(anum));
        }
      } else {
        atoms.reserve(atoms.size() + (final - initial + 1));
        for (int anum = initial; anum <= final; anum++) {
          add_atom(cvm::atom(anum));
        }
      }

    }
    if (cvm::get_error()) return COLVARS_ERROR;
  } else {
    cvm::error("Error: no valid definition for \"atomNumbersRange\", \""+
               range_conf+"\".\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  }

  return COLVARS_OK;
}


int cvm::atom_group::add_atom_name_residue_range(std::string const &psf_segid,
                                                 std::string const &range_conf)
{
  if (range_conf.size()) {
    std::istringstream is(range_conf);
    std::string atom_name;
    int initial, final;
    char dash;
    if ( (is >> atom_name) && (atom_name.size())  &&
         (is >> initial) && (initial > 0) &&
         (is >> dash) && (dash == '-') &&
         (is >> final) && (final > 0) ) {

      atoms_ids.reserve(atoms_ids.size() + (final - initial + 1));

      if (is_enabled(f_ag_scalable)) {
        for (int resid = initial; resid <= final; resid++) {
          add_atom_id((cvm::proxy)->check_atom_id(resid, atom_name, psf_segid));
        }
      } else {
        atoms.reserve(atoms.size() + (final - initial + 1));
        for (int resid = initial; resid <= final; resid++) {
          add_atom(cvm::atom(resid, atom_name, psf_segid));
        }
      }

      if (cvm::get_error()) return COLVARS_ERROR;
    } else {
      cvm::error("Error: cannot parse definition for \""
                 "atomNameResidueRange\", \""+
                 range_conf+"\".\n");
      return COLVARS_ERROR;
    }
  } else {
    cvm::error("Error: atomNameResidueRange with empty definition.\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}


std::string const cvm::atom_group::print_atom_ids() const
{
  size_t line_count = 0;
  std::ostringstream os("");
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    os << " " << std::setw(9) << atoms_ids[i];
    if (++line_count == 7) {
      os << "\n";
      line_count = 0;
    }
  }
  return os.str();
}


int cvm::atom_group::parse_fitting_options(std::string const &group_conf)
{
  if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {

    if (b_dummy)
      cvm::error("Error: centerToReference or rotateToReference "
                 "cannot be defined for a dummy atom.\n");

    bool b_ref_pos_group = false;
    std::string fitting_group_conf;
    if (key_lookup(group_conf, "refPositionsGroup", &fitting_group_conf)) {
      b_ref_pos_group = true;
      cvm::log("Warning: keyword \"refPositionsGroup\" is deprecated: please use \"fittingGroup\" instead.\n");
    }

    if (b_ref_pos_group || key_lookup(group_conf, "fittingGroup", &fitting_group_conf)) {
      // instead of this group, define another group to compute the fit
      if (fitting_group) {
        cvm::error("Error: the atom group \""+
                   key+"\" has already a reference group "
                   "for the rototranslational fit, which was communicated by the "
                   "colvar component.  You should not use fittingGroup "
                   "in this case.\n", COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
      }
      cvm::log("Within atom group \""+key+"\":\n");
      fitting_group = new atom_group("fittingGroup");
      if (fitting_group->parse(fitting_group_conf) == COLVARS_OK) {
        fitting_group->check_keywords(fitting_group_conf, "fittingGroup");
        if (cvm::get_error()) {
          cvm::error("Error setting up atom group \"fittingGroup\".", COLVARS_INPUT_ERROR);
          return COLVARS_INPUT_ERROR;
        }
      }
      enable(f_ag_fitting_group);
    }

    atom_group *group_for_fit = fitting_group ? fitting_group : this;

    get_keyval(group_conf, "refPositions", ref_pos, ref_pos);

    std::string ref_pos_file;
    if (get_keyval(group_conf, "refPositionsFile", ref_pos_file, std::string(""))) {

      if (ref_pos.size()) {
        cvm::error("Error: cannot specify \"refPositionsFile\" and "
                   "\"refPositions\" at the same time.\n");
      }

      std::string ref_pos_col;
      double ref_pos_col_value=0.0;

      if (get_keyval(group_conf, "refPositionsCol", ref_pos_col, std::string(""))) {
        // if provided, use PDB column to select coordinates
        bool found = get_keyval(group_conf, "refPositionsColValue", ref_pos_col_value, 0.0);
        if (found && ref_pos_col_value == 0.0) {
          cvm::error("Error: refPositionsColValue, "
                     "if provided, must be non-zero.\n", COLVARS_INPUT_ERROR);
          return COLVARS_ERROR;
        }
      }

      ref_pos.resize(group_for_fit->size());
      cvm::load_coords(ref_pos_file.c_str(), &ref_pos, group_for_fit,
                       ref_pos_col, ref_pos_col_value);
    }

    if (ref_pos.size()) {

      if (is_enabled(f_ag_rotate)) {
        if (ref_pos.size() != group_for_fit->size())
          cvm::error("Error: the number of reference positions provided("+
                     cvm::to_str(ref_pos.size())+
                     ") does not match the number of atoms within \""+
                     key+
                     "\" ("+cvm::to_str(group_for_fit->size())+
                     "): to perform a rotational fit, "+
                     "these numbers should be equal.\n", COLVARS_INPUT_ERROR);
      }

      // save the center of geometry of ref_pos and subtract it
      center_ref_pos();

    } else {
      cvm::error("Error: no reference positions provided.\n", COLVARS_INPUT_ERROR);
      return COLVARS_ERROR;
    }

    if (is_enabled(f_ag_rotate) && !noforce) {
      cvm::log("Warning: atom group \""+key+
               "\" will be aligned to a fixed orientation given by the reference positions provided.  "
               "If the internal structure of the group changes too much (i.e. its RMSD is comparable "
               "to its radius of gyration), the optimal rotation and its gradients may become discontinuous.  "
               "If that happens, use fittingGroup (or a different definition for it if already defined) "
               "to align the coordinates.\n");
      // initialize rot member data
      rot.request_group1_gradients(group_for_fit->size());
    }
  }

  // Enable fit gradient calculation only if necessary, and not disabled by the user
  // This must happen after fitting group is defined so that side-effects are performed
  // properly (ie. allocating fitting group gradients)
  {
    bool b_fit_gradients;
    get_keyval(group_conf, "enableFitGradients", b_fit_gradients, true);

    if (b_fit_gradients && (is_enabled(f_ag_center) || is_enabled(f_ag_rotate))) {
      enable(f_ag_fit_gradients);
    }
  }

  return COLVARS_OK;
}


void cvm::atom_group::do_feature_side_effects(int id)
{
  // If enabled features are changed upstream, the features below should be refreshed
  switch (id) {
    case f_ag_fit_gradients:
      if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {
        atom_group *group_for_fit = fitting_group ? fitting_group : this;
        group_for_fit->fit_gradients.assign(group_for_fit->size(), cvm::atom_pos(0.0, 0.0, 0.0));
        rot.request_group1_gradients(group_for_fit->size());
      }
      break;
  }
}


int cvm::atom_group::create_sorted_ids()
{
  // Only do the work if the vector is not yet populated
  if (sorted_atoms_ids.size())
    return COLVARS_OK;

  // Sort the internal IDs
  std::list<int> sorted_atoms_ids_list;
  for (size_t i = 0; i < atoms_ids.size(); i++) {
    sorted_atoms_ids_list.push_back(atoms_ids[i]);
  }
  sorted_atoms_ids_list.sort();
  sorted_atoms_ids_list.unique();
  if (sorted_atoms_ids_list.size() != atoms_ids.size()) {
    return cvm::error("Error: duplicate atom IDs in atom group? (found " +
                      cvm::to_str(sorted_atoms_ids_list.size()) +
                      " unique atom IDs instead of " +
                      cvm::to_str(atoms_ids.size()) + ").\n", COLVARS_BUG_ERROR);
  }

  // Compute map between sorted and unsorted elements
  sorted_atoms_ids.resize(atoms_ids.size());
  sorted_atoms_ids_map.resize(atoms_ids.size());
  std::list<int>::iterator lsii = sorted_atoms_ids_list.begin();
  size_t ii = 0;
  for ( ; ii < atoms_ids.size(); lsii++, ii++) {
    sorted_atoms_ids[ii] = *lsii;
    size_t const pos = std::find(atoms_ids.begin(), atoms_ids.end(), *lsii) -
      atoms_ids.begin();
    sorted_atoms_ids_map[ii] = pos;
  }

  return COLVARS_OK;
}


int cvm::atom_group::overlap(const atom_group &g1, const atom_group &g2){
  for (cvm::atom_const_iter ai1 = g1.begin(); ai1 != g1.end(); ai1++) {
    for (cvm::atom_const_iter ai2 = g2.begin(); ai2 != g2.end(); ai2++) {
      if (ai1->id == ai2->id) {
        return (ai1->id + 1); // 1-based index to allow boolean usage
      }
    }
  }
  return 0;
}


void cvm::atom_group::center_ref_pos()
{
  ref_pos_cog = cvm::atom_pos(0.0, 0.0, 0.0);
  std::vector<cvm::atom_pos>::iterator pi;
  for (pi = ref_pos.begin(); pi != ref_pos.end(); ++pi) {
    ref_pos_cog += *pi;
  }
  ref_pos_cog /= (cvm::real) ref_pos.size();
  for (pi = ref_pos.begin(); pi != ref_pos.end(); ++pi) {
    (*pi) -= ref_pos_cog;
  }
}


void cvm::atom_group::read_positions()
{
  if (b_dummy) return;

  for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    ai->read_position();
  }

  if (fitting_group)
    fitting_group->read_positions();
}


int cvm::atom_group::calc_required_properties()
{
  // TODO check if the com is needed?
  calc_center_of_mass();
  calc_center_of_geometry();

  if (!is_enabled(f_ag_scalable)) {
    if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {
      if (fitting_group) {
        fitting_group->calc_center_of_geometry();
      }

      calc_apply_roto_translation();

      // update COM and COG after fitting
      calc_center_of_geometry();
      calc_center_of_mass();
      if (fitting_group) {
        fitting_group->calc_center_of_geometry();
      }
    }
  }

  // TODO calculate elements of scalable cvc's here before reduction

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


void cvm::atom_group::calc_apply_roto_translation()
{
  // store the laborarory-frame COGs for when they are needed later
  cog_orig = this->center_of_geometry();
  if (fitting_group) {
    fitting_group->cog_orig = fitting_group->center_of_geometry();
  }

  if (is_enabled(f_ag_center)) {
    // center on the origin first
    cvm::atom_pos const rpg_cog = fitting_group ?
      fitting_group->center_of_geometry() : this->center_of_geometry();
    apply_translation(-1.0 * rpg_cog);
    if (fitting_group) {
      fitting_group->apply_translation(-1.0 * rpg_cog);
    }
  }

  if (is_enabled(f_ag_rotate)) {
    // rotate the group (around the center of geometry if f_ag_center is
    // enabled, around the origin otherwise)
    rot.calc_optimal_rotation(fitting_group ?
                              fitting_group->positions() :
                              this->positions(),
                              ref_pos);

    cvm::atom_iter ai;
    for (ai = this->begin(); ai != this->end(); ai++) {
      ai->pos = rot.rotate(ai->pos);
    }
    if (fitting_group) {
      for (ai = fitting_group->begin(); ai != fitting_group->end(); ai++) {
        ai->pos = rot.rotate(ai->pos);
      }
    }
  }

  if (is_enabled(f_ag_center) && !is_enabled(f_ag_center_origin)) {
    // align with the center of geometry of ref_pos
    apply_translation(ref_pos_cog);
    if (fitting_group) {
      fitting_group->apply_translation(ref_pos_cog);
    }
  }
  // update of COM and COG is done from the calling routine
}


void cvm::atom_group::apply_translation(cvm::rvector const &t)
{
  if (b_dummy) {
    cvm::error("Error: cannot translate the coordinates of a dummy atom group.\n", COLVARS_INPUT_ERROR);
    return;
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: cannot translate the coordinates of a scalable atom group.\n", COLVARS_INPUT_ERROR);
    return;
  }

  for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    ai->pos += t;
  }
}


void cvm::atom_group::read_velocities()
{
  if (b_dummy) return;

  if (is_enabled(f_ag_rotate)) {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->read_velocity();
      ai->vel = rot.rotate(ai->vel);
    }

  } else {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->read_velocity();
    }
  }
}


// TODO make this a calc function
void cvm::atom_group::read_total_forces()
{
  if (b_dummy) return;

  if (is_enabled(f_ag_rotate)) {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->read_total_force();
      ai->total_force = rot.rotate(ai->total_force);
    }

  } else {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->read_total_force();
    }
  }
}


int cvm::atom_group::calc_center_of_geometry()
{
  if (b_dummy) {
    cog = dummy_atom_pos;
  } else {
    cog.reset();
    for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
      cog += ai->pos;
    }
    cog /= cvm::real(this->size());
  }
  return COLVARS_OK;
}


int cvm::atom_group::calc_center_of_mass()
{
  if (b_dummy) {
    com = dummy_atom_pos;
    if (cvm::debug()) {
      cvm::log("Dummy atom center of mass = "+cvm::to_str(com)+"\n");
    }
  } else if (is_enabled(f_ag_scalable)) {
    com = (cvm::proxy)->get_atom_group_com(index);
  } else {
    com.reset();
    for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
      com += ai->mass * ai->pos;
    }
    com /= total_mass;
  }
  return COLVARS_OK;
}


int cvm::atom_group::calc_dipole(cvm::atom_pos const &dipole_center)
{
  if (b_dummy) {
    return cvm::error("Error: trying to compute the dipole "
                      "of a dummy group.\n", COLVARS_INPUT_ERROR);
  }
  dip.reset();
  for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
    dip += ai->charge * (ai->pos - dipole_center);
  }
  return COLVARS_OK;
}


void cvm::atom_group::set_weighted_gradient(cvm::rvector const &grad)
{
  if (b_dummy) return;

  scalar_com_gradient = grad;

  if (!is_enabled(f_ag_scalable)) {
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->grad = (ai->mass/total_mass) * grad;
    }
  }
}


void cvm::atom_group::calc_fit_gradients()
{
  if (b_dummy || ! is_enabled(f_ag_fit_gradients)) return;

  if (cvm::debug())
    cvm::log("Calculating fit gradients.\n");

  cvm::atom_group *group_for_fit = fitting_group ? fitting_group : this;

  if (is_enabled(f_ag_center)) {
    // add the center of geometry contribution to the gradients
    cvm::rvector atom_grad;

    for (size_t i = 0; i < this->size(); i++) {
      atom_grad += atoms[i].grad;
    }
    if (is_enabled(f_ag_rotate)) atom_grad = (rot.inverse()).rotate(atom_grad);
    atom_grad *= (-1.0)/(cvm::real(group_for_fit->size()));

    for (size_t j = 0; j < group_for_fit->size(); j++) {
      group_for_fit->fit_gradients[j] = atom_grad;
    }
  }

  if (is_enabled(f_ag_rotate)) {

    // add the rotation matrix contribution to the gradients
    cvm::rotation const rot_inv = rot.inverse();

    for (size_t i = 0; i < this->size(); i++) {

      // compute centered, unrotated position
      cvm::atom_pos const pos_orig =
        rot_inv.rotate((is_enabled(f_ag_center) ? (atoms[i].pos - ref_pos_cog) : (atoms[i].pos)));

      // calculate \partial(R(q) \vec{x}_i)/\partial q) \cdot \partial\xi/\partial\vec{x}_i
      cvm::quaternion const dxdq =
        rot.q.position_derivative_inner(pos_orig, atoms[i].grad);

      for (size_t j = 0; j < group_for_fit->size(); j++) {
        // multiply by {\partial q}/\partial\vec{x}_j and add it to the fit gradients
        for (size_t iq = 0; iq < 4; iq++) {
          group_for_fit->fit_gradients[j] += dxdq[iq] * rot.dQ0_1[j][iq];
        }
      }
    }
  }

  if (cvm::debug())
    cvm::log("Done calculating fit gradients.\n");
}


std::vector<cvm::atom_pos> cvm::atom_group::positions() const
{
  if (b_dummy) {
    cvm::error("Error: positions are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic positions are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  std::vector<cvm::atom_pos> x(this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator xi = x.begin();
  for ( ; ai != this->end(); ++xi, ++ai) {
    *xi = ai->pos;
  }
  return x;
}

std::vector<cvm::atom_pos> cvm::atom_group::positions_shifted(cvm::rvector const &shift) const
{
  if (b_dummy) {
    cvm::error("Error: positions are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic positions are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  std::vector<cvm::atom_pos> x(this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator xi = x.begin();
  for ( ; ai != this->end(); ++xi, ++ai) {
    *xi = (ai->pos + shift);
  }
  return x;
}

std::vector<cvm::rvector> cvm::atom_group::velocities() const
{
  if (b_dummy) {
    cvm::error("Error: velocities are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic velocities are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  std::vector<cvm::rvector> v(this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator vi = v.begin();
  for ( ; ai != this->end(); vi++, ai++) {
    *vi = ai->vel;
  }
  return v;
}

std::vector<cvm::rvector> cvm::atom_group::total_forces() const
{
  if (b_dummy) {
    cvm::error("Error: total forces are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic total forces are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  std::vector<cvm::rvector> f(this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator fi = f.begin();
  for ( ; ai != this->end(); ++fi, ++ai) {
    *fi = ai->total_force;
  }
  return f;
}


// TODO make this an accessor
cvm::rvector cvm::atom_group::total_force() const
{
  if (b_dummy) {
    cvm::error("Error: total total forces are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    return (cvm::proxy)->get_atom_group_total_force(index);
  }

  cvm::rvector f(0.0);
  for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
    f += ai->total_force;
  }
  return f;
}


void cvm::atom_group::apply_colvar_force(cvm::real const &force)
{
  if (cvm::debug()) {
    log("Communicating a colvar force from atom group to the MD engine.\n");
  }

  if (b_dummy) return;

  if (noforce) {
    cvm::error("Error: sending a force to a group that has "
               "\"enableForces\" set to off.\n");
    return;
  }

  if (is_enabled(f_ag_scalable)) {
    (cvm::proxy)->apply_atom_group_force(index, force * scalar_com_gradient);
    return;
  }

  if (is_enabled(f_ag_rotate)) {

    // rotate forces back to the original frame
    cvm::rotation const rot_inv = rot.inverse();
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->apply_force(rot_inv.rotate(force * ai->grad));
    }

  } else {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->apply_force(force * ai->grad);
    }
  }

  if ((is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) && is_enabled(f_ag_fit_gradients)) {

    atom_group *group_for_fit = fitting_group ? fitting_group : this;

    // Fit gradients are already calculated in "laboratory" frame
    for (size_t j = 0; j < group_for_fit->size(); j++) {
      (*group_for_fit)[j].apply_force(force * group_for_fit->fit_gradients[j]);
    }
  }
}


void cvm::atom_group::apply_force(cvm::rvector const &force)
{
  if (cvm::debug()) {
    log("Communicating a colvar force from atom group to the MD engine.\n");
  }

  if (b_dummy) return;

  if (noforce) {
    cvm::error("Error: sending a force to a group that has "
               "\"enableForces\" set to off.\n");
    return;
  }

  if (is_enabled(f_ag_scalable)) {
    (cvm::proxy)->apply_atom_group_force(index, force);
    return;
  }

  if (is_enabled(f_ag_rotate)) {

    cvm::rotation const rot_inv = rot.inverse();
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->apply_force(rot_inv.rotate((ai->mass/total_mass) * force));
    }

  } else {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->apply_force((ai->mass/total_mass) * force);
    }
  }
}


// Static members

std::vector<colvardeps::feature *> cvm::atom_group::ag_features;



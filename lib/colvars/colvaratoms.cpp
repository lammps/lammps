#include "colvaratoms.h"
#include "colvar_rotation_derivative.h"
#include "colvardeps.h"
#include "colvarproxy.h"
#include "colvarmodule.h"

#include <numeric>
#include <algorithm>
#include <list>
#include <iomanip>


std::vector<cvm::real> cvm::atom_group::pos_aos_to_soa(const std::vector<cvm::atom_pos>& aos_in) {
  static_assert(sizeof(cvm::atom_pos) == 3 * sizeof(cvm::real), "The size of cvm::atom_pos requires to be 3 * the size of cvm::real.");
  std::vector<cvm::real> pos_soa(aos_in.size() * 3);
  const size_t offset_x = 0;
  const size_t offset_y = offset_x + aos_in.size();
  const size_t offset_z = offset_y + aos_in.size();
  for (size_t i = 0; i < aos_in.size(); ++i) {
    pos_soa[i + offset_x] = aos_in[i].x;
    pos_soa[i + offset_y] = aos_in[i].y;
    pos_soa[i + offset_z] = aos_in[i].z;
  }
  return pos_soa;
}

cvm::atom_group::simple_atom cvm::atom_group::init_atom_from_proxy(
  colvarproxy* const    p,
  cvm::residue_id const &residue,
  std::string const     &atom_name,
  std::string const     &segment_id) {
  int atom_proxy_index = p->init_atom(residue, atom_name, segment_id);
  if (cvm::get_error()) {
    // Error condition already reported by p->init_atom()
    // return bogus atom, counting on caller to check for errors
    return cvm::atom_group::simple_atom{
      /*.proxy_index = */atom_proxy_index,
      /*.id = */0,
      /*.mass = */0.0,
      /*.charge = */0.0,
      /*.pos = */{0, 0, 0},
      /*.vel = */{0, 0, 0},
      /*.total_force = */{0, 0, 0},
      /*.grad = */{0, 0, 0}};
  }
  const int atom_id = p->get_atom_id(atom_proxy_index);
  const cvm::real atom_mass = p->get_atom_mass(atom_proxy_index);
  const cvm::real atom_charge = p->get_atom_charge(atom_proxy_index);
  return cvm::atom_group::simple_atom{
    /*.proxy_index = */atom_proxy_index,
    /*.id = */atom_id,
    /*.mass = */atom_mass,
    /*.charge = */atom_charge,
    /*.pos = */{0, 0, 0},
    /*.vel = */{0, 0, 0},
    /*.total_force = */{0, 0, 0},
    /*.grad = */{0, 0, 0}};
}

cvm::atom_group::simple_atom cvm::atom_group::init_atom_from_proxy(
  colvarproxy* const    p,
  int atom_number) {
  const int atom_proxy_index = p->init_atom(atom_number);
  const int atom_id = p->get_atom_id(atom_proxy_index);
  const cvm::real atom_mass = p->get_atom_mass(atom_proxy_index);
  const cvm::real atom_charge = p->get_atom_charge(atom_proxy_index);
  return cvm::atom_group::simple_atom{
    /*.proxy_index = */atom_proxy_index,
    /*.id = */atom_id,
    /*.mass = */atom_mass,
    /*.charge = */atom_charge,
    /*.pos = */{0, 0, 0},
    /*.vel = */{0, 0, 0},
    /*.total_force = */{0, 0, 0},
    /*.grad = */{0, 0, 0}};
}

cvm::atom_group::simple_atom cvm::atom_group::init_atom_from_proxy(
  colvarproxy* const    p,
  const simple_atom& atom_in) {
  const int atom_proxy_index = atom_in.proxy_index;
  const int atom_id = p->get_atom_id(atom_proxy_index);
  p->increase_refcount(atom_proxy_index);
  const cvm::real atom_mass = p->get_atom_mass(atom_proxy_index);
  const cvm::real atom_charge = p->get_atom_charge(atom_proxy_index);
  return cvm::atom_group::simple_atom{
    /*.proxy_index = */atom_proxy_index,
    /*.id = */atom_id,
    /*.mass = */atom_mass,
    /*.charge = */atom_charge,
    /*.pos = */{0, 0, 0},
    /*.vel = */{0, 0, 0},
    /*.total_force = */{0, 0, 0},
    /*.grad = */{0, 0, 0}};
}

cvm::atom_group::atom_group():
  b_dummy(false),
  fitting_group(nullptr),
  noforce(false), b_user_defined_fit(false),
  rot_deriv(nullptr), num_atoms(0), index(-1),
  num_ref_pos(0)
{
  key = "unnamed";
  init();
}

cvm::atom_group::atom_group(char const *key_in): atom_group() {
  key = std::string(key_in);
  init();
}

cvm::atom_group::~atom_group() {
  clear_soa();
  if (is_enabled(f_ag_scalable) && !b_dummy) {
    (cvm::main()->proxy)->clear_atom_group(index);
    index = -1;
  }

  if (fitting_group) {
    delete fitting_group;
    fitting_group = NULL;
  }

  if (rot_deriv != nullptr) {
    delete rot_deriv;
    rot_deriv = nullptr;
  }
}

cvm::atom_group::atom_modifier::atom_modifier(cvm::atom_group* ag): m_ag(ag)
{
  if (m_ag->modify_lock.try_lock()) {
    update_from_soa();
  } else {
    throw cvm::error(
      "You are trying to modify the SOA atom group more than once at the same time",
      COLVARS_BUG_ERROR);
  }
}

cvm::atom_group::atom_modifier::~atom_modifier() {
  m_ag->modify_lock.unlock();
  if (cvm::get_error() == COLVARS_OK) {
    sync_to_soa();
  }
}

int cvm::atom_group::init()
{
  if (!key.size()) key = "unnamed";
  description = "atom group " + key;
  // These may be overwritten by parse(), if a name is provided

  // atoms.clear();
  this->clear_soa();
  atom_group::init_dependencies();
  index = -1;

  b_dummy = false;
  b_user_defined_fit = false;
  fitting_group = NULL;
  rot_deriv = nullptr;

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
  for (i = feature_states.size(); i < colvardeps::f_ag_ntot; i++) {
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

int cvm::atom_group::setup() {
  // TODO: What should I do?
  // if (atoms_ids.size() == 0) {
  //   cvm::error("Bug: atoms_ids is empty!", COLVARS_BUG_ERROR);
  // }
  // Update masses and charges
  const colvarproxy *p = cvm::proxy;
  if (!b_dummy && !is_enabled(f_ag_scalable)) {
    atoms_charge.resize(size());
    atoms_mass.resize(size());
    for (size_t i = 0; i < num_atoms; ++i) {
      atoms_charge[i] = p->get_atom_charge(atoms_index[i]);
      atoms_mass[i] = p->get_atom_mass(atoms_index[i]);
    }
  }
  update_total_charge();
  update_total_mass();
  return COLVARS_OK;
}

void cvm::atom_group::atom_modifier::update_from_soa() {
  m_atoms_ids = m_ag->atoms_ids;
  if (!m_ag->is_enabled(f_ag_scalable)) {
    m_atoms.resize(m_ag->size());
    for (size_t i = 0; i < m_atoms.size(); ++i) {
      m_atoms[i] = (*m_ag)[i];
      m_atoms_ids_count[m_ag->atoms_ids[i]] = 1;
    }
  }
  m_total_mass = m_ag->total_mass;
  m_total_charge = m_ag->total_charge;
}

void cvm::atom_group::atom_modifier::sync_to_soa() const {
  m_ag->num_atoms = m_atoms.size();
  m_ag->atoms_ids = m_atoms_ids;
  if (!m_ag->is_enabled(f_ag_scalable)) {
    m_ag->atoms_index.resize(m_ag->num_atoms);
    m_ag->atoms_pos.resize(3 * m_ag->num_atoms);
    m_ag->atoms_charge.resize(m_ag->num_atoms);
    m_ag->atoms_vel.resize(3 * m_ag->num_atoms);
    m_ag->atoms_mass.resize(m_ag->num_atoms);
    m_ag->atoms_grad.resize(3 * m_ag->num_atoms);
    m_ag->atoms_total_force.resize(3 * m_ag->num_atoms);
    m_ag->atoms_weight.resize(m_ag->num_atoms);
    // colvarproxy *p = cvm::main()->proxy;
    if (m_ag->atoms_ids.size() != m_atoms.size()) {
      cvm::error("Number of atom IDs does not match the number of atoms!\n",
                COLVARS_BUG_ERROR);
    }
    // for (size_t i = 0; i < m_ag->atoms_ids.size(); ++i) {
    //   m_atoms[i].proxy_index = p->init_atom(m_ag->atoms_ids[i]);
    // }
    for (size_t i = 0; i < m_ag->num_atoms; ++i) {
      m_ag->atoms_index[i] = m_atoms[i].proxy_index;
      m_ag->pos_x(i) = m_atoms[i].pos.x;
      m_ag->pos_y(i) = m_atoms[i].pos.y;
      m_ag->pos_z(i) = m_atoms[i].pos.z;
      m_ag->vel_x(i) = m_atoms[i].vel.x;
      m_ag->vel_y(i) = m_atoms[i].vel.y;
      m_ag->vel_z(i) = m_atoms[i].vel.z;
      m_ag->grad_x(i) = m_atoms[i].grad.x;
      m_ag->grad_y(i) = m_atoms[i].grad.y;
      m_ag->grad_z(i) = m_atoms[i].grad.z;
      m_ag->total_force_x(i) = m_atoms[i].total_force.x;
      m_ag->total_force_y(i) = m_atoms[i].total_force.y;
      m_ag->total_force_z(i) = m_atoms[i].total_force.z;
      m_ag->atoms_mass[i] = m_atoms[i].mass;
      m_ag->atoms_charge[i] = m_atoms[i].charge;
    }
  }
  // m_ag->total_charge = m_total_charge;
  // m_ag->total_mass = m_total_mass;
  // m_ag->update_total_charge();
  // m_ag->update_total_mass();
}

int cvm::atom_group::atom_modifier::add_atom(const simple_atom& a) {
  if (a.id < 0) {
    return cvm::error("Error: invalid atom number " + cvm::to_str(a.id), COLVARS_INPUT_ERROR);
  }
  // for (size_t i = 0; i < m_atoms_ids.size(); i++) {
  //   if (m_atoms_ids[i] == a.id) {
  //     if (cvm::debug())
  //       cvm::log("Discarding doubly counted atom with number "+
  //                cvm::to_str(a.id+1)+".\n");
  //     return COLVARS_OK;
  //   }
  // }
  auto map_it = m_atoms_ids_count.find(a.id);
  if (map_it != m_atoms_ids_count.end()) {
    if (cvm::debug()) {
      cvm::log("Discarding doubly counted atom with number "+
              cvm::to_str(a.id+1)+".\n");
    }
    map_it->second++;
    return COLVARS_OK;
  } else {
    m_atoms_ids_count[a.id] = 1;
  }
  // for consistency with add_atom_id(), we update the list as well
  m_atoms_ids.push_back(a.id);
  m_atoms.push_back(a);
  m_total_mass += a.mass;
  m_total_charge += a.charge;

  return COLVARS_OK;
}

int cvm::atom_group::atom_modifier::remove_atom(atom_modifier::atom_iter ai)
{
  if (m_ag->is_enabled(f_ag_scalable)) {
    cvm::error("Error: cannot remove atoms from a scalable group.\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  }

  if (!this->size()) {
    cvm::error("Error: trying to remove an atom from an empty group.\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  } else {
    if (ai->proxy_index > 0) {
      colvarproxy *p = cvm::main()->proxy;
      p->clear_atom(ai->proxy_index);
    }
    m_total_mass -= ai->mass;
    m_total_charge -= ai->charge;
    const int id_to_remove = *(m_atoms_ids.begin() + (ai - m_atoms.begin()));
    auto map_it = m_atoms_ids_count.find(id_to_remove);
    if (map_it != m_atoms_ids_count.end()) {
      m_atoms_ids_count.erase(map_it);
    }
    m_atoms_ids.erase(m_atoms_ids.begin() + (ai - m_atoms.begin()));
    m_atoms.erase(ai);
  }

  return COLVARS_OK;
}

int cvm::atom_group::atom_modifier::add_atom_numbers(std::string const &numbers_conf)
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
    m_atoms_ids.reserve(m_atoms_ids.size()+atom_indexes.size());

    if (m_ag->is_enabled(f_ag_scalable)) {
      for (size_t i = 0; i < atom_indexes.size(); i++) {
        add_atom_id((cvm::proxy)->check_atom_id(atom_indexes[i]));
      }
    } else {
      // if we are handling the group on rank 0, better allocate the vector in one shot
      m_atoms.reserve(m_atoms.size()+atom_indexes.size());
      colvarproxy* const p = cvm::main()->proxy;
      for (size_t i = 0; i < atom_indexes.size(); i++) {
        add_atom(init_atom_from_proxy(p, atom_indexes[i]));
      }
    }
    if (cvm::get_error()) return COLVARS_ERROR;
  } else {
    return cvm::error("Error: no numbers provided for \""
               "atomNumbers\".\n", COLVARS_INPUT_ERROR);
  }

  return COLVARS_OK;
}

int cvm::atom_group::atom_modifier::add_atoms_of_group(atom_group const *ag) {
  // TODO: should I check if ag == this??
  std::vector<int> const &source_ids = ag->atoms_ids;
  if (source_ids.size()) {
    m_atoms_ids.reserve(m_atoms_ids.size()+source_ids.size());

    if (m_ag->is_enabled(f_ag_scalable)) {
      for (size_t i = 0; i < source_ids.size(); i++) {
        add_atom_id(source_ids[i]);
      }
    } else {
      m_atoms.reserve(m_atoms.size()+source_ids.size());
      colvarproxy* const p = cvm::main()->proxy;
      for (size_t i = 0; i < source_ids.size(); i++) {
        // We could use the atom copy constructor, but only if the source
        // group is not scalable - whereas this works in both cases
        // atom constructor expects 1-based atom number
        add_atom(init_atom_from_proxy(p, source_ids[i] + 1));
      }
    }

    if (cvm::get_error()) return COLVARS_ERROR;
  } else {
    cvm::error("Error: source atom group contains no atoms\".\n", COLVARS_INPUT_ERROR);
    return COLVARS_ERROR;
  }

  return COLVARS_OK;
}

int cvm::atom_group::atom_modifier::add_index_group(std::string const &index_group_name, bool silent) {
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
    if (silent)
      return COLVARS_INPUT_ERROR;
    else
      return cvm::error("Error: could not find index group "+
                      index_group_name+" among those already provided.\n",
                      COLVARS_INPUT_ERROR);
  }

  int error_code = COLVARS_OK;

  std::vector<int> const &index_group = *(index_groups[i_group]);

  m_atoms_ids.reserve(m_atoms_ids.size()+index_group.size());

  if (m_ag->is_enabled(f_ag_scalable)) {
    for (size_t i = 0; i < index_group.size(); i++) {
      error_code |= add_atom_id((cvm::proxy)->check_atom_id(index_group[i]));
    }
  } else {
    m_atoms.reserve(m_atoms.size()+index_group.size());
    colvarproxy* const p = cvm::main()->proxy;
    for (size_t i = 0; i < index_group.size(); i++) {
      error_code |= add_atom(init_atom_from_proxy(p, index_group[i]));
    }
  }

  return error_code;
}

int cvm::atom_group::atom_modifier::add_atom_numbers_range(std::string const &range_conf) {
  if (range_conf.size()) {
    std::istringstream is(range_conf);
    int initial, final;
    char dash;
    if ( (is >> initial) && (initial > 0) &&
         (is >> dash) && (dash == '-') &&
         (is >> final) && (final > 0) ) {

      m_atoms_ids.reserve(m_atoms_ids.size() + (final - initial + 1));

      if (m_ag->is_enabled(f_ag_scalable)) {
        for (int anum = initial; anum <= final; anum++) {
          add_atom_id((cvm::proxy)->check_atom_id(anum));
        }
      } else {
        m_atoms.reserve(m_atoms.size() + (final - initial + 1));
        colvarproxy* const p = cvm::main()->proxy;
        for (int anum = initial; anum <= final; anum++) {
          add_atom(init_atom_from_proxy(p, anum));
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

int cvm::atom_group::atom_modifier::add_atom_name_residue_range(
  std::string const &psf_segid,
  std::string const &range_conf) {
  if (range_conf.size()) {
    std::istringstream is(range_conf);
    std::string atom_name;
    int initial, final;
    char dash;
    if ( (is >> atom_name) && (atom_name.size())  &&
         (is >> initial) && (initial > 0) &&
         (is >> dash) && (dash == '-') &&
         (is >> final) && (final > 0) ) {

      m_atoms_ids.reserve(m_atoms_ids.size() + (final - initial + 1));

      if (m_ag->is_enabled(f_ag_scalable)) {
        for (int resid = initial; resid <= final; resid++) {
          add_atom_id((cvm::proxy)->check_atom_id(resid, atom_name, psf_segid));
        }
      } else {
        m_atoms.reserve(m_atoms.size() + (final - initial + 1));
        colvarproxy* const p = cvm::main()->proxy;
        for (int resid = initial; resid <= final; resid++) {
          add_atom(init_atom_from_proxy(p, resid, atom_name, psf_segid));
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

void cvm::atom_group::clear_soa() {
  for (size_t i = 0; i < atoms_index.size(); ++i) {
    (cvm::main()->proxy)->clear_atom(atoms_index[i]);
  }
  // Then clear all arrays in SOA
  std::fill(atoms_index.begin(), atoms_index.end(), 0);
  std::fill(atoms_pos.begin(), atoms_pos.end(), 0);
  std::fill(atoms_charge.begin(), atoms_charge.end(), 0);
  std::fill(atoms_vel.begin(), atoms_vel.end(), 0);
  std::fill(atoms_mass.begin(), atoms_mass.end(), 0);
  std::fill(atoms_grad.begin(), atoms_grad.end(), 0);
  std::fill(atoms_total_force.begin(), atoms_total_force.end(), 0);
  std::fill(atoms_weight.begin(), atoms_weight.end(), 0);
  // Reset the number of atoms
  num_atoms = 0;
  // Other fields should be untouched
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

  int error_code = COLVARS_OK;

  // Optional group name will let other groups reuse atom definition
  if (get_keyval(group_conf, "name", name)) {
    if ((cvm::atom_group_soa_by_name(this->name) != NULL) &&
        (cvm::atom_group_soa_by_name(this->name) != this)) {
      cvm::error("Error: this atom group cannot have the same name, \""+this->name+
                        "\", as another atom group.\n",
                COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
    cvm::main()->register_named_atom_group_soa(this);
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
      atom_group * ag = atom_group_soa_by_name(atoms_of);
      if (ag == NULL) {
        cvm::error("Error: cannot find atom group with name " + atoms_of + ".\n");
        return COLVARS_ERROR;
      }
      auto modify_atoms = get_atom_modifier();
      error_code |= modify_atoms.add_atoms_of_group(ag);
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
      auto modify_atoms = get_atom_modifier();
      error_code |= modify_atoms.add_atom_numbers(numbers_conf);
      numbers_conf = "";
    }
  }

  {
    std::string index_group_name;
    if (get_keyval(group_conf, "indexGroup", index_group_name)) {
      // use an index group from the index file read globally
      auto modify_atoms = get_atom_modifier();
      error_code |= modify_atoms.add_index_group(index_group_name);
    }
  }

  {
    std::string range_conf = "";
    size_t pos = 0;
    while (key_lookup(group_conf, "atomNumbersRange",
                      &range_conf, &pos)) {
      auto modify_atoms = get_atom_modifier();
      error_code |= modify_atoms.add_atom_numbers_range(range_conf);
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
        auto modify_atoms = get_atom_modifier();
        error_code |= modify_atoms.add_atom_name_residue_range(psf_segids.size() ?
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

      error_code |= cvm::main()->proxy->load_atoms_pdb(atoms_file_name.c_str(), *this, atoms_col,
                                                       atoms_col_value);
    }
  }

  // Catch any errors from all the initialization steps above
  if (error_code || cvm::get_error()) return (error_code || cvm::get_error());

  // checks of doubly-counted atoms have been handled by add_atom() already

  if (get_keyval(group_conf, "dummyAtom", dummy_atom_pos, cvm::atom_pos())) {

    error_code |= set_dummy();
    error_code |= set_dummy_pos(dummy_atom_pos);

  } else {

    if (!(atoms_ids.size())) {
      error_code |= cvm::error("Error: no atoms defined for atom group \"" + key + "\".\n",
                               COLVARS_INPUT_ERROR);
    }

    // whether these atoms will ever receive forces or not
    bool enable_forces = true;
    get_keyval(group_conf, "enableForces", enable_forces, true, colvarparse::parse_silent);
    noforce = !enable_forces;
  }

  // Now that atoms are defined we can parse the detailed fitting options
  error_code |= parse_fitting_options(group_conf);

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

  if (is_enabled(f_ag_rotate)) setup_rotation_derivative();

  return error_code;
}

int cvm::atom_group::parse_fitting_options(std::string const &group_conf) {
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

    std::vector<cvm::atom_pos> ref_pos_aos;
    get_keyval(group_conf, "refPositions", ref_pos_aos, ref_pos_aos);

    std::string ref_pos_file;
    if (get_keyval(group_conf, "refPositionsFile", ref_pos_file, std::string(""))) {

      if (ref_pos_aos.size()) {
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

      ref_pos_aos.resize(group_for_fit->size());
      cvm::load_coords(ref_pos_file.c_str(), &ref_pos_aos, group_for_fit,
                       ref_pos_col, ref_pos_col_value);
    }

    if (ref_pos_aos.size()) {

      if (is_enabled(f_ag_rotate)) {
        if (ref_pos_aos.size() != group_for_fit->size())
          cvm::error("Error: the number of reference positions provided("+
                     cvm::to_str(ref_pos_aos.size())+
                     ") does not match the number of atoms within \""+
                     key+
                     "\" ("+cvm::to_str(group_for_fit->size())+
                     "): to perform a rotational fit, "+
                     "these numbers should be equal.\n", COLVARS_INPUT_ERROR);
      }

      // Convert ref_pos_aos to SOA
      num_ref_pos = ref_pos_aos.size();
      ref_pos.resize(3 * num_ref_pos);
      for (size_t i = 0; i < num_ref_pos; ++i) {
        ref_pos_x(i) = ref_pos_aos[i].x;
        ref_pos_y(i) = ref_pos_aos[i].y;
        ref_pos_z(i) = ref_pos_aos[i].z;
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

int cvm::atom_group::atom_modifier::add_atom_id(int aid) {
  if (aid < 0) {
    return COLVARS_ERROR;
  }

  // for (size_t i = 0; i < m_atoms_ids.size(); i++) {
  //   if (m_atoms_ids[i] == aid) {
  //     if (cvm::debug())
  //       cvm::log("Discarding doubly counted atom with number "+
  //                cvm::to_str(aid+1)+".\n");
  //     return COLVARS_OK;
  //   }
  // }
  auto map_it = m_atoms_ids_count.find(aid);
  if (map_it != m_atoms_ids_count.end()) {
    if (cvm::debug()) {
      cvm::log("Discarding doubly counted atom with number "+
              cvm::to_str(aid+1)+".\n");
    }
    map_it->second++;
    return COLVARS_OK;
  } else {
    m_atoms_ids_count[aid] = 1;
  }

  m_atoms_ids.push_back(aid);
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

int cvm::atom_group::set_dummy_pos(cvm::atom_pos const &pos) {
  if (b_dummy) {
    dummy_atom_pos = pos;
  } else {
    return cvm::error("Error: setting dummy position for group with keyword \""+
                      key+"\" and name \""+name+
                      "\", but it is not dummy.\n", COLVARS_INPUT_ERROR);
  }
  return COLVARS_OK;
}

void cvm::atom_group::update_total_mass() {
  if (b_dummy) {
    total_mass = 1.0;
    return;
  }

  if (is_enabled(f_ag_scalable)) {
    total_mass = (cvm::main()->proxy)->get_atom_group_mass(index);
  } else {
    // total_mass = 0.0
    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   total_mass += ai->mass;
    // }
    total_mass = std::accumulate(atoms_mass.begin(), atoms_mass.end(), 0.0);
    const double t_m = total_mass;
    std::transform(atoms_mass.begin(), atoms_mass.end(),
                   atoms_weight.begin(), [t_m](cvm::real x){return x/t_m;});
  }
  if (total_mass < 1e-15) {
    cvm::error("Error: " + description + " has zero total mass.\n");
  }
}

void cvm::atom_group::update_total_charge() {
  if (b_dummy) {
    total_charge = 0.0;
    return;
  }

  if (is_enabled(f_ag_scalable)) {
    total_charge = (cvm::main()->proxy)->get_atom_group_charge(index);
  } else {
    // total_charge = 0.0;
    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   total_charge += ai->charge;
    // }
    total_charge = std::accumulate(atoms_charge.begin(), atoms_charge.end(), 0.0);
  }
}

void cvm::atom_group::print_properties(std::string const &colvar_name, int i, int j) {
  if (cvm::proxy->updated_masses() && cvm::proxy->updated_charges()) {
    cvm::log("Re-initialized atom group for variable \""+colvar_name+"\":"+
             cvm::to_str(i)+"/"+
             cvm::to_str(j)+". "+ cvm::to_str(atoms_ids.size())+
             " atoms: total mass = "+cvm::to_str(total_mass)+
             ", total charge = "+cvm::to_str(total_charge)+".\n");
  }
}

std::string const cvm::atom_group::print_atom_ids() const {
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

int cvm::atom_group::create_sorted_ids()
{
  // Only do the work if the vector is not yet populated
  if (sorted_atoms_ids.size())
    return COLVARS_OK;

  // Sort the internal IDs
  std::list<int> sorted_atoms_ids_list(atoms_ids.begin(), atoms_ids.end());
  // for (size_t i = 0; i < atoms_ids.size(); i++) {
  //   sorted_atoms_ids_list.push_back(atoms_ids[i]);
  // }
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

int cvm::atom_group::overlap(const atom_group &g1, const atom_group &g2) {
  // for (cvm::atom_const_iter ai1 = g1.begin(); ai1 != g1.end(); ai1++) {
  //   for (cvm::atom_const_iter ai2 = g2.begin(); ai2 != g2.end(); ai2++) {
  //     if (ai1->id == ai2->id) {
  //       return (ai1->id + 1); // 1-based index to allow boolean usage
  //     }
  //   }
  // }
  for (auto ai1 = g1.atoms_ids.cbegin(); ai1 != g1.atoms_ids.cend(); ++ai1) {
    for (auto ai2 = g2.atoms_ids.cbegin(); ai2 != g2.atoms_ids.cend(); ++ai2) {
      if ((*ai1) == (*ai2)) {
        return (*ai1) + 1; // 1-based index to allow boolean usage
      }
    }
  }
  return 0;
}

void cvm::atom_group::read_positions()
{
  if (b_dummy) return;

  // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
  //   ai->read_position();
  // }
  // TODO: Should I check if this group is scalable?
  if (!is_enabled(f_ag_scalable)) {
    colvarproxy *p = cvm::main()->proxy;
    for (size_t i = 0; i < num_atoms; ++i) {
      const int proxy_index = atoms_index[i];
      const cvm::rvector pos = p->get_atom_position(proxy_index);
      pos_x(i) = pos.x;
      pos_y(i) = pos.y;
      pos_z(i) = pos.z;
    }
  }

  if (fitting_group)
    fitting_group->read_positions();
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

  if (is_enabled(f_ag_fit_gradients) && !b_dummy) {
    // Save the unrotated frame for fit gradients
    // pos_unrotated.resize(size());
    // for (size_t i = 0; i < size(); ++i) {
    //   pos_unrotated[i] = atoms[i].pos;
    // }
    atoms_pos_unrotated = atoms_pos;
  }

  if (is_enabled(f_ag_rotate)) {
    // rotate the group (around the center of geometry if f_ag_center is
    // enabled, around the origin otherwise)
    auto* group_for_fit = fitting_group ? fitting_group : this;
    rot.calc_optimal_rotation_soa(
      group_for_fit->atoms_pos, ref_pos,
      group_for_fit->size(), num_ref_pos);
    const auto rot_mat = rot.matrix();

    // cvm::atom_iter ai;
    // for (ai = this->begin(); ai != this->end(); ai++) {
    //   ai->pos = rot_mat * ai->pos;
    // }
    // if (fitting_group) {
    //   for (ai = fitting_group->begin(); ai != fitting_group->end(); ai++) {
    //     ai->pos = rot_mat * ai->pos;
    //   }
    // }
    this->rotate(rot_mat);
    if (fitting_group/* && (fitting_group != this)*/) {
      // TODO: Should I check if fitting_group != this
      fitting_group->rotate(rot_mat);
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

void cvm::atom_group::setup_rotation_derivative() {
  if (rot_deriv != nullptr) delete rot_deriv;
  auto* group_for_fit = fitting_group ? fitting_group : this;
  rot_deriv = new rotation_derivative(
    rot, group_for_fit->atoms_pos, ref_pos,
    group_for_fit->size(), num_ref_pos
  );
}

void cvm::atom_group::center_ref_pos()
{
  ref_pos_cog = cvm::atom_pos(0.0, 0.0, 0.0);
  // std::vector<cvm::atom_pos>::iterator pi;
  // for (pi = ref_pos.begin(); pi != ref_pos.end(); ++pi) {
  //   ref_pos_cog += *pi;
  // }
  // ref_pos_cog /= (cvm::real) ref_pos.size();
  // for (pi = ref_pos.begin(); pi != ref_pos.end(); ++pi) {
  //   (*pi) -= ref_pos_cog;
  // }
  for (size_t i = 0; i < num_ref_pos; ++i) {
    ref_pos_cog.x += ref_pos_x(i);
    ref_pos_cog.y += ref_pos_y(i);
    ref_pos_cog.z += ref_pos_z(i);
  }
  ref_pos_cog /= num_ref_pos;
  for (size_t i = 0; i < num_ref_pos; ++i) {
    ref_pos_x(i) -= ref_pos_cog.x;
    ref_pos_y(i) -= ref_pos_cog.y;
    ref_pos_z(i) -= ref_pos_cog.z;
  }
}

void cvm::atom_group::rotate(const cvm::rmatrix& rot_mat) {
  for (size_t i = 0; i < num_atoms; ++i) {
    const cvm::real new_x = rot_mat.xx * pos_x(i) + rot_mat.xy * pos_y(i) + rot_mat.xz * pos_z(i);
    const cvm::real new_y = rot_mat.yx * pos_x(i) + rot_mat.yy * pos_y(i) + rot_mat.yz * pos_z(i);
    const cvm::real new_z = rot_mat.zx * pos_x(i) + rot_mat.zy * pos_y(i) + rot_mat.zz * pos_z(i);
    pos_x(i) = new_x;
    pos_y(i) = new_y;
    pos_z(i) = new_z;
  }
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

  // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
  //   ai->pos += t;
  // }
  for (size_t i = 0; i < num_atoms; ++i) {
    pos_x(i) += t.x;
    pos_y(i) += t.y;
    pos_z(i) += t.z;
  }
}

void cvm::atom_group::read_velocities()
{
  if (b_dummy) return;

  if (is_enabled(f_ag_rotate)) {

    const auto rot_mat = rot.matrix();
    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   ai->read_velocity();
    //   ai->vel = rot_mat * ai->vel;
    // }
    colvarproxy *p = cvm::main()->proxy;
    for (size_t i = 0; i < num_atoms; ++i) {
      const int proxy_index = atoms_index[i];
      const cvm::rvector vel = p->get_atom_velocity(proxy_index);
      vel_x(i) = rot_mat.xx * vel.x + rot_mat.xy * vel.y + rot_mat.xz * vel.z;
      vel_y(i) = rot_mat.yx * vel.x + rot_mat.yy * vel.y + rot_mat.yz * vel.z;
      vel_z(i) = rot_mat.zx * vel.x + rot_mat.zy * vel.y + rot_mat.zz * vel.z;
    }

  } else {
    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   ai->read_velocity();
    // }
    colvarproxy *p = cvm::main()->proxy;
    for (size_t i = 0; i < num_atoms; ++i) {
      const int proxy_index = atoms_index[i];
      const cvm::rvector vel = p->get_atom_velocity(proxy_index);
      vel_x(i) = vel.x;
      vel_y(i) = vel.y;
      vel_z(i) = vel.z;
    }
  }
}

void cvm::atom_group::read_total_forces()
{
  if (b_dummy) return;

  if (is_enabled(f_ag_rotate)) {

    const auto rot_mat = rot.matrix();
    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   ai->read_total_force();
    //   ai->total_force = rot_mat * ai->total_force;
    // }
    colvarproxy *p = cvm::main()->proxy;
    for (size_t i = 0; i < num_atoms; ++i) {
      const int proxy_index = atoms_index[i];
      const cvm::rvector total_force = p->get_atom_total_force(proxy_index);
      total_force_x(i) = rot_mat.xx * total_force.x +
                         rot_mat.xy * total_force.y +
                         rot_mat.xz * total_force.z;
      total_force_y(i) = rot_mat.yx * total_force.x +
                         rot_mat.yy * total_force.y +
                         rot_mat.yz * total_force.z;
      total_force_z(i) = rot_mat.zx * total_force.x +
                         rot_mat.zy * total_force.y +
                         rot_mat.zz * total_force.z;
    }

  } else {

    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   ai->read_total_force();
    // }
    colvarproxy *p = cvm::main()->proxy;
    for (size_t i = 0; i < num_atoms; ++i) {
      const int proxy_index = atoms_index[i];
      const cvm::rvector total_force = p->get_atom_total_force(proxy_index);
      total_force_x(i) = total_force.x;
      total_force_y(i) = total_force.y;
      total_force_z(i) = total_force.z;
    }
  }
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

int cvm::atom_group::calc_center_of_geometry()
{
  if (b_dummy) {
    cog = dummy_atom_pos;
  } else {
    cog.reset();
    // for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
    //   cog += ai->pos;
    // }
    for (size_t i = 0; i < num_atoms; ++i) {
      cog.x += pos_x(i);
      cog.y += pos_y(i);
      cog.z += pos_z(i);
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
    // for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
    //   com += ai->mass * ai->pos;
    // }
    for (size_t i = 0; i < num_atoms; ++i) {
      com.x += atoms_mass[i] * pos_x(i);
      com.y += atoms_mass[i] * pos_y(i);
      com.z += atoms_mass[i] * pos_z(i);
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
  // for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
  //   dip += ai->charge * (ai->pos - dipole_center);
  // }
  for (size_t i = 0; i < num_atoms; ++i) {
    dip.x += atoms_charge[i] * (pos_x(i) - dipole_center.x);
    dip.y += atoms_charge[i] * (pos_y(i) - dipole_center.y);
    dip.z += atoms_charge[i] * (pos_z(i) - dipole_center.z);
  }
  return COLVARS_OK;
}

void cvm::atom_group::set_weighted_gradient(cvm::rvector const &grad)
{
  if (b_dummy) return;

  scalar_com_gradient = grad;

  if (!is_enabled(f_ag_scalable)) {
    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   ai->grad = (ai->mass/total_mass) * grad;
    // }
    for (size_t i = 0; i < num_atoms; ++i) {
      grad_x(i) = atoms_weight[i] * grad.x;
      grad_y(i) = atoms_weight[i] * grad.y;
      grad_z(i) = atoms_weight[i] * grad.z;
    }
  }
}

void cvm::atom_group::calc_fit_gradients()
{
  if (b_dummy || ! is_enabled(f_ag_fit_gradients)) return;

  if (cvm::debug())
    cvm::log("Calculating fit gradients.\n");

  cvm::atom_group *group_for_fit = fitting_group ? fitting_group : this;

  auto accessor_main = [this](size_t i){
    // return atoms[i].grad;
    return cvm::rvector{grad_x(i), grad_y(i), grad_z(i)};
  };
  auto accessor_fitting = [&group_for_fit](size_t j, const cvm::rvector& grad){
    // group_for_fit->fit_gradients[j] = grad;
    group_for_fit->fit_gradients_x(j) = grad.x;
    group_for_fit->fit_gradients_y(j) = grad.y;
    group_for_fit->fit_gradients_z(j) = grad.z;
  };
  if (is_enabled(f_ag_center) && is_enabled(f_ag_rotate))
    calc_fit_forces_impl<true, true>(accessor_main, accessor_fitting);
  if (is_enabled(f_ag_center) && !is_enabled(f_ag_rotate))
    calc_fit_forces_impl<true, false>(accessor_main, accessor_fitting);
  if (!is_enabled(f_ag_center) && is_enabled(f_ag_rotate))
    calc_fit_forces_impl<false, true>(accessor_main, accessor_fitting);
  if (!is_enabled(f_ag_center) && !is_enabled(f_ag_rotate))
    calc_fit_forces_impl<false, false>(accessor_main, accessor_fitting);

  if (cvm::debug())
    cvm::log("Done calculating fit gradients.\n");
}

template <bool B_ag_center, bool B_ag_rotate,
          typename main_force_accessor_T, typename fitting_force_accessor_T>
void cvm::atom_group::calc_fit_forces_impl(
  main_force_accessor_T accessor_main,
  fitting_force_accessor_T accessor_fitting) const {
  const cvm::atom_group *group_for_fit = fitting_group ? fitting_group : this;
  // the center of geometry contribution to the gradients
  cvm::rvector atom_grad;
  // the rotation matrix contribution to the gradients
  const auto rot_inv = rot.inverse().matrix();
  // temporary variables for computing and summing derivatives
  cvm::real sum_dxdq[4] = {0, 0, 0, 0};
  cvm::vector1d<cvm::rvector> dq0_1(4);
  // loop 1: iterate over the current atom group
  for (size_t i = 0; i < size(); i++) {
    if (B_ag_center) {
      atom_grad += accessor_main(i);
    }
    if (B_ag_rotate) {
      // calculate \partial(R(q) \vec{x}_i)/\partial q) \cdot \partial\xi/\partial\vec{x}_i
      cvm::quaternion const dxdq =
        rot.q.position_derivative_inner(
          cvm::rvector{pos_unrotated_x(i),
                       pos_unrotated_y(i),
                       pos_unrotated_z(i)},
          accessor_main(i));
      sum_dxdq[0] += dxdq[0];
      sum_dxdq[1] += dxdq[1];
      sum_dxdq[2] += dxdq[2];
      sum_dxdq[3] += dxdq[3];
    }
  }
  if (B_ag_center) {
    if (B_ag_rotate) atom_grad = rot_inv * atom_grad;
    atom_grad *= (-1.0)/(cvm::real(group_for_fit->size()));
  }
  // loop 2: iterate over the fitting group
  if (B_ag_rotate) rot_deriv->prepare_derivative(rotation_derivative_dldq::use_dq);
  for (size_t j = 0; j < group_for_fit->size(); j++) {
    cvm::rvector fitting_force_grad{0, 0, 0};
    if (B_ag_center) {
      fitting_force_grad += atom_grad;
    }
    if (B_ag_rotate) {
      rot_deriv->calc_derivative_wrt_group1<false, true, false>(j, nullptr, &dq0_1);
      // multiply by {\partial q}/\partial\vec{x}_j and add it to the fit gradients
      fitting_force_grad += sum_dxdq[0] * dq0_1[0] +
                            sum_dxdq[1] * dq0_1[1] +
                            sum_dxdq[2] * dq0_1[2] +
                            sum_dxdq[3] * dq0_1[3];
    }
    if (cvm::debug()) {
      cvm::log(cvm::to_str(fitting_force_grad));
    }
    accessor_fitting(j, fitting_force_grad);
  }
}

template <typename main_force_accessor_T, typename fitting_force_accessor_T>
void cvm::atom_group::calc_fit_forces(
  main_force_accessor_T accessor_main,
  fitting_force_accessor_T accessor_fitting) const {
  if (is_enabled(f_ag_center) && is_enabled(f_ag_rotate))
    calc_fit_forces_impl<true, true, main_force_accessor_T, fitting_force_accessor_T>(accessor_main, accessor_fitting);
  if (is_enabled(f_ag_center) && !is_enabled(f_ag_rotate))
    calc_fit_forces_impl<true, false, main_force_accessor_T, fitting_force_accessor_T>(accessor_main, accessor_fitting);
  if (!is_enabled(f_ag_center) && is_enabled(f_ag_rotate))
    calc_fit_forces_impl<false, true, main_force_accessor_T, fitting_force_accessor_T>(accessor_main, accessor_fitting);
  if (!is_enabled(f_ag_center) && !is_enabled(f_ag_rotate))
    calc_fit_forces_impl<false, false, main_force_accessor_T, fitting_force_accessor_T>(accessor_main, accessor_fitting);
}

std::vector<cvm::real> cvm::atom_group::positions() const
{
  if (b_dummy) {
    cvm::error("Error: positions are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic positions are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  // std::vector<cvm::atom_pos> x(this->size(), 0.0);
  // cvm::atom_const_iter ai = this->begin();
  // std::vector<cvm::atom_pos>::iterator xi = x.begin();
  // for ( ; ai != this->end(); ++xi, ++ai) {
  //   *xi = ai->pos;
  // }
  return atoms_pos;
}

std::vector<cvm::real> cvm::atom_group::positions_shifted(cvm::rvector const &shift) const
{
  if (b_dummy) {
    cvm::error("Error: positions are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic positions are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  // std::vector<cvm::atom_pos> x(this->size(), 0.0);
  // cvm::atom_const_iter ai = this->begin();
  // std::vector<cvm::atom_pos>::iterator xi = x.begin();
  // for ( ; ai != this->end(); ++xi, ++ai) {
  //   *xi = (ai->pos + shift);
  // }
  // return x;
  std::vector<cvm::real> shifted = atoms_pos;
  for (size_t i = 0; i < num_atoms; ++i) {
    shifted[i]             += shift.x;
    shifted[i+num_atoms]   += shift.y;
    shifted[i+2*num_atoms] += shift.z;
  }
  return shifted;
}

std::vector<cvm::real> cvm::atom_group::velocities() const {
  if (b_dummy) {
    cvm::error("Error: velocities are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic velocities are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  return atoms_vel;
}

std::vector<cvm::real> cvm::atom_group::total_forces() const {
  if (b_dummy) {
    cvm::error("Error: velocities are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    cvm::error("Error: atomic velocities are not available "
               "from a scalable atom group.\n", COLVARS_INPUT_ERROR);
  }

  return atoms_total_force;
}

cvm::rvector cvm::atom_group::total_force() const {
  if (b_dummy) {
    cvm::error("Error: total total forces are not available "
               "from a dummy atom group.\n", COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_ag_scalable)) {
    return (cvm::proxy)->get_atom_group_total_force(index);
  }
  cvm::rvector f(0, 0, 0);
  for (size_t i = 0; i < num_atoms; ++i) {
    f.x += total_force_x(i);
    f.y += total_force_y(i);
    f.z += total_force_z(i);
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
    const auto rot_inv = rot.inverse().matrix();
    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   ai->apply_force(rot_inv * (force * ai->grad));
    // }
    colvarproxy* const p = cvm::main()->proxy;
    for (size_t i = 0; i < num_atoms; ++i) {
      const int proxy_index = atoms_index[i];
      // The calculation is reordered since A(λ*f) = λ*(Af)
      // where λ is a scalar, f is a vector and A is a matrix.
      const cvm::rvector f{
        force * (rot_inv.xx * grad_x(i) + rot_inv.xy * grad_y(i) + rot_inv.xz * grad_z(i)),
        force * (rot_inv.yx * grad_x(i) + rot_inv.yy * grad_y(i) + rot_inv.yz * grad_z(i)),
        force * (rot_inv.zx * grad_x(i) + rot_inv.zy * grad_y(i) + rot_inv.zz * grad_z(i))
      };
      p->apply_atom_force(proxy_index, f);
    }

  } else {

    // for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    //   ai->apply_force(force * ai->grad);
    // }
    colvarproxy* const p = cvm::main()->proxy;
    for (size_t i = 0; i < num_atoms; ++i) {
      const int proxy_index = atoms_index[i];
      const cvm::rvector f{
        force * grad_x(i),
        force * grad_y(i),
        force * grad_z(i)
      };
      p->apply_atom_force(proxy_index, f);
    }
  }

  if ((is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) && is_enabled(f_ag_fit_gradients)) {

    atom_group *group_for_fit = fitting_group ? fitting_group : this;
    colvarproxy* const p = cvm::main()->proxy;
    // Fit gradients are already calculated in "laboratory" frame
    for (size_t j = 0; j < group_for_fit->size(); j++) {
      // (*group_for_fit)[j].apply_force(force * group_for_fit->fit_gradients[j]);
      const int proxy_index = group_for_fit->atoms_index[j];
      const cvm::rvector f{
        force * group_for_fit->fit_gradients_x(j),
        force * group_for_fit->fit_gradients_y(j),
        force * group_for_fit->fit_gradients_z(j)
      };
      p->apply_atom_force(proxy_index, f);
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

  auto ag_force = get_group_force_object();
  for (size_t i = 0; i < size(); ++i) {
    ag_force.add_atom_force(i, atoms_weight[i] * force);
  }
}

void cvm::atom_group::reset_atoms_data() {
  std::fill(atoms_pos.begin(), atoms_pos.end(), 0);
  std::fill(atoms_vel.begin(), atoms_vel.end(), 0);
  std::fill(atoms_grad.begin(), atoms_grad.end(), 0);
  std::fill(atoms_total_force.begin(), atoms_total_force.end(), 0);
}

void cvm::atom_group::do_feature_side_effects(int id)
{
  // If enabled features are changed upstream, the features below should be refreshed
  switch (id) {
    case f_ag_fit_gradients:
      if (is_enabled(f_ag_center) || is_enabled(f_ag_rotate)) {
        atom_group *group_for_fit = fitting_group ? fitting_group : this;
        // group_for_fit->fit_gradients.assign(group_for_fit->size(), cvm::atom_pos(0.0, 0.0, 0.0));
        // num_fit_gradients = group_for_fit->size();
        group_for_fit->fit_gradients.assign(3 * group_for_fit->size(), 0);
      }
      break;
  }
}

void cvm::atom_group::set_ref_pos_from_aos(const std::vector<cvm::atom_pos>& pos_aos) {
  num_ref_pos = pos_aos.size();
  ref_pos = cvm::atom_group::pos_aos_to_soa(pos_aos);
}

cvm::atom_group::group_force_object cvm::atom_group::get_group_force_object() {
  return cvm::atom_group::group_force_object(this);
}

cvm::atom_group::group_force_object::group_force_object(cvm::atom_group* ag):
m_ag(ag), m_group_for_fit(m_ag->fitting_group ? m_ag->fitting_group : m_ag),
m_has_fitting_force(m_ag->is_enabled(f_ag_center) || m_ag->is_enabled(f_ag_rotate)) {
  if (m_has_fitting_force) {
    if (m_ag->group_forces.size() != 3 * m_ag->size()) {
      m_ag->group_forces.assign(3 * m_ag->size(), 0);
    } else {
      std::fill(m_ag->group_forces.begin(),
                m_ag->group_forces.end(), 0);
    }
  }
}

cvm::atom_group::group_force_object::~group_force_object() {
  if (m_has_fitting_force) {
    apply_force_with_fitting_group();
  }
}

void cvm::atom_group::group_force_object::add_atom_force(
  size_t i, const cvm::rvector& force) {
  if (m_has_fitting_force) {
    // m_ag->group_forces[i] += force;
    m_ag->group_forces_x(i) += force.x;
    m_ag->group_forces_y(i) += force.y;
    m_ag->group_forces_z(i) += force.z;
  } else {
    // Apply the force directly if we don't use fitting
    // (*m_ag)[i].apply_force(force);
    colvarproxy* const p = cvm::main()->proxy;
    const int proxy_index = m_ag->atoms_index[i];
    p->apply_atom_force(proxy_index, force);
  }
}

void cvm::atom_group::group_force_object::apply_force_with_fitting_group() {
  const cvm::rmatrix rot_inv = m_ag->rot.inverse().matrix();
  if (cvm::debug()) {
    cvm::log("Applying force on main group " + m_ag->name + ":\n");
  }
  if (m_ag->b_dummy) return;
  colvarproxy* const p = cvm::main()->proxy;
  if (m_ag->is_enabled(f_ag_rotate)) {
    for (size_t ia = 0; ia < m_ag->size(); ++ia) {
      // const cvm::rvector f_ia = rot_inv * m_ag->group_forces[ia];
      // (*m_ag)[ia].apply_force(f_ia);
      const int proxy_index = m_ag->atoms_index[ia];
      const cvm::rvector f_ia{
        rot_inv.xx * m_ag->group_forces_x(ia) +
        rot_inv.xy * m_ag->group_forces_y(ia) +
        rot_inv.xz * m_ag->group_forces_z(ia),
        rot_inv.yx * m_ag->group_forces_x(ia) +
        rot_inv.yy * m_ag->group_forces_y(ia) +
        rot_inv.yz * m_ag->group_forces_z(ia),
        rot_inv.zx * m_ag->group_forces_x(ia) +
        rot_inv.zy * m_ag->group_forces_y(ia) +
        rot_inv.zz * m_ag->group_forces_z(ia),
      };
      p->apply_atom_force(proxy_index, f_ia);
      if (cvm::debug()) {
        cvm::log(cvm::to_str(f_ia));
      }
    }
  } else {
    for (size_t ia = 0; ia < m_ag->size(); ++ia) {
      const int proxy_index = m_ag->atoms_index[ia];
      const cvm::rvector f_ia{
        m_ag->group_forces_x(ia),
        m_ag->group_forces_y(ia),
        m_ag->group_forces_z(ia)};
      p->apply_atom_force(proxy_index, f_ia);
      if (cvm::debug()) {
        cvm::log(cvm::to_str(f_ia));
      }
    }
  }
  // Gradients are only available with scalar components, so for a scalar component,
  // if f_ag_fit_gradients is disabled, then the forces on the fitting group is not
  // computed. For a vector component, we can only know the forces on the fitting
  // group, but checking this flag can mimic results that the users expect (if
  // "enableFitGradients no" then there is no force on the fitting group).
  if (m_ag->is_enabled(f_ag_fit_gradients)) {
    colvarproxy* const p = cvm::main()->proxy;
    auto accessor_main = [this](size_t i){
      return cvm::rvector(m_ag->group_forces_x(i),
                          m_ag->group_forces_y(i),
                          m_ag->group_forces_z(i));};
    auto accessor_fitting = [this, p](size_t j, const cvm::rvector& fitting_force){
      // (*(m_group_for_fit))[j].apply_force(fitting_force);
      const int proxy_index = m_group_for_fit->atoms_index[j];
      p->apply_atom_force(proxy_index, fitting_force);
    };
    if (cvm::debug()) {
      cvm::log("Applying force on the fitting group of main group" + m_ag->name + ":\n");
    }
    m_ag->calc_fit_forces(accessor_main, accessor_fitting);
    if (cvm::debug()) {
      cvm::log("Done applying force on the fitting group of main group" + m_ag->name + ":\n");
    }
  }
}

std::vector<colvardeps::feature *> cvm::atom_group::ag_features;

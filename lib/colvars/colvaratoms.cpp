/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvaratoms.h"


cvm::atom::atom()
{
  index = -1;
  id = -1;
  reset_data();
}


cvm::atom::atom(int atom_number)
{
  colvarproxy *p = cvm::proxy;
  index = p->init_atom(atom_number);
  if (cvm::debug()) {
    cvm::log("The index of this atom in the colvarproxy arrays is "+
             cvm::to_str(index)+".\n");
  }
  id = p->get_atom_id(index);
  update_mass();
  reset_data();
}


cvm::atom::atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id)
{
  colvarproxy *p = cvm::proxy;
  index = p->init_atom(residue, atom_name, segment_id);
  if (cvm::debug()) {
    cvm::log("The index of this atom in the colvarproxy_namd arrays is "+
             cvm::to_str(index)+".\n");
  }
  id = p->get_atom_id(index);
  update_mass();
  reset_data();
}


cvm::atom::atom(atom const &a)
  : index(a.index)
{
  id = (cvm::proxy)->get_atom_id(index);
  update_mass();
  reset_data();
}


cvm::atom::~atom()
{
  if (index >= 0) {
    (cvm::proxy)->clear_atom(index);
  }
}



// TODO change this arrangement
// Note: "conf" is the configuration of the cvc who is using this atom group;
// "key" is the name of the atom group (e.g. "atoms", "group1", "group2", ...)
cvm::atom_group::atom_group(std::string const &conf,
                            char const        *key_in)
{
  key = key_in;
  cvm::log("Defining atom group \""+
           std::string(key)+"\".\n");
  init();
  // real work is done by parse
  parse(conf);
  setup();
}


cvm::atom_group::atom_group(std::vector<cvm::atom> const &atoms_in)
{
  init();
  atoms = atoms_in;
  setup();
}


cvm::atom_group::atom_group()
{
  init();
}


cvm::atom_group::~atom_group()
{
  if (index >= 0) {
    (cvm::proxy)->clear_atom_group(index);
  }

  if (ref_pos_group) {
    delete ref_pos_group;
    ref_pos_group = NULL;
  }
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
  if (b_scalable) {
    cvm::error("Error: cannot remove atoms from a scalable group.\n", INPUT_ERROR);
    return COLVARS_ERROR;
  }

  if (!this->size()) {
    cvm::error("Error: trying to remove an atom from an empty group.\n", INPUT_ERROR);
    return COLVARS_ERROR;
  } else {
    total_mass -= ai->mass;
    total_charge -= ai->charge;
    atoms_ids.erase(atoms_ids.begin() + (ai - atoms.begin()));
    atoms.erase(ai);
  }

  return COLVARS_OK;
}


int cvm::atom_group::init()
{
  if (!key.size()) key = "atoms";

  atoms.clear();

  b_scalable = false;
  index = -1;

  b_center = false;
  b_rotate = false;
  b_user_defined_fit = false;
  b_fit_gradients = false;
  ref_pos_group = NULL;

  noforce = false;

  total_mass = 0.0;
  total_charge = 0.0;

  cog.reset();
  com.reset();

  return COLVARS_OK;
}


int cvm::atom_group::setup()
{
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

  if (b_scalable) {
    total_mass = (cvm::proxy)->get_atom_group_mass(index);
  } else {
    total_mass = 0.0;
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      total_mass += ai->mass;
    }
  }
}


void cvm::atom_group::reset_mass(std::string &name, int i, int j)
{
  update_total_mass();
  cvm::log("Re-initialized atom group "+name+":"+cvm::to_str(i)+"/"+
           cvm::to_str(j)+". "+ cvm::to_str(atoms_ids.size())+
           " atoms: total mass = "+cvm::to_str(total_mass)+".\n");
}


void cvm::atom_group::update_total_charge()
{
  if (b_dummy) {
    total_charge = 0.0;
    return;
  }

  if (b_scalable) {
    total_charge = (cvm::proxy)->get_atom_group_charge(index);
  } else {
    total_charge = 0.0;
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      total_charge += ai->charge;
    }
  }
}


int cvm::atom_group::parse(std::string const &conf)
{
  std::string group_conf;

  // TODO move this to the cvc class constructor/init

  // save_delimiters is set to false for this call, because "conf" is
  // not the config string of this group, but of its parent object
  // (which has already taken care of the delimiters)
  save_delimiters = false;
  key_lookup(conf, key.c_str(), group_conf, dummy_pos);
  // restoring the normal value, because we do want keywords checked
  // inside "group_conf"
  save_delimiters = true;

  if (group_conf.size() == 0) {
    cvm::error("Error: atom group \""+key+
               "\" is set, but has no definition.\n",
               INPUT_ERROR);
    return COLVARS_ERROR;
  }

  cvm::increase_depth();

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

  // if the cvc allows it, the flag has been set to true by default
  get_keyval(group_conf, "scalable", b_scalable, b_scalable);

  {
    std::string numbers_conf = "";
    size_t pos = 0;
    while (key_lookup(group_conf, "atomNumbers", numbers_conf, pos)) {
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
                      range_conf, pos)) {
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
                   (*psii)+"\".\n", INPUT_ERROR);
      }
    }

    std::string range_conf = "";
    size_t pos = 0;
    size_t range_count = 0;
    psii = psf_segids.begin();
    while (key_lookup(group_conf, "atomNameResidueRange",
                      range_conf, pos)) {
      range_count++;
      if (psf_segids.size() && (range_count > psf_segids.size())) {
        cvm::error("Error: more instances of \"atomNameResidueRange\" than "
                   "values of \"psfSegID\".\n", INPUT_ERROR);
      } else {
        parse_error |= add_atom_name_residue_range(psf_segids.size() ? *psii : std::string(""),
                                                   range_conf);
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
                   INPUT_ERROR);
      }

      double atoms_col_value;
      bool const atoms_col_value_defined = get_keyval(group_conf, "atomsColValue", atoms_col_value, 0.0);
      if (atoms_col_value_defined && (!atoms_col_value)) {
        cvm::error("Error: atomsColValue, if provided, must be non-zero.\n", INPUT_ERROR);
      }

      cvm::load_atoms(atoms_file_name.c_str(), *this, atoms_col, atoms_col_value);
    }
  }

  // Catch any errors from all the initialization steps above
  if (parse_error || cvm::get_error()) return COLVARS_ERROR;

  if (b_scalable) {
    index = (cvm::proxy)->init_atom_group(atoms_ids);
  }

  // checks of doubly-counted atoms have been handled by add_atom() already

  if (get_keyval(group_conf, "dummyAtom", dummy_atom_pos, cvm::atom_pos())) {
    b_dummy = true;
    // note: atoms_ids.size() is used here in lieu of atoms.size(),
    // which can be empty for scalable groups
    if (atoms_ids.size()) {
      cvm::error("Error: cannot set up group \""+
                 key+"\" as a dummy atom "
                 "and provide it with atom definitions.\n", INPUT_ERROR);
    }
  } else {
    b_dummy = false;

    if (!(atoms_ids.size())) {
      cvm::error("Error: no atoms defined for atom group \""+
                 key+"\".\n", INPUT_ERROR);
    }

    // whether these atoms will ever receive forces or not
    bool enable_forces = true;
    // disableForces is deprecated
    if (get_keyval(group_conf, "enableForces", enable_forces, true)) {
      noforce = !enable_forces;
    } else {
      get_keyval(group_conf, "disableForces", noforce, false, colvarparse::parse_silent);
    }
  }

  parse_error |= parse_fitting_options(group_conf);

  // TODO move this to colvarparse object
  check_keywords(group_conf, key.c_str());
  if (cvm::get_error()) {
    cvm::error("Error setting up atom group \""+key+"\".");
    return COLVARS_ERROR;
  }

  // Calculate all required properties (such as total mass)
  setup();

  if (cvm::debug())
    cvm::log("Done initializing atom group \""+key+"\".\n");

  cvm::log("Atom group \""+key+"\" defined, "+
            cvm::to_str(atoms_ids.size())+" atoms initialized: total mass = "+
	    cvm::to_str(total_mass)+", total charge = "+
            cvm::to_str(total_charge)+".\n");

  cvm::decrease_depth();

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
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

    if (b_scalable) {
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
               "atomNumbers\".\n", INPUT_ERROR);
    return COLVARS_ERROR;
  }

  return COLVARS_OK;
}


int cvm::atom_group::add_index_group(std::string const &index_group_name)
{
  std::list<std::string>::iterator names_i = cvm::index_group_names.begin();
  std::list<std::vector<int> >::iterator index_groups_i = cvm::index_groups.begin();
  for ( ; names_i != cvm::index_group_names.end() ; ++names_i, ++index_groups_i) {
    if (*names_i == index_group_name)
      break;
  }

  if (names_i == cvm::index_group_names.end()) {
    cvm::error("Error: could not find index group "+
               index_group_name+" among those provided by the index file.\n",
               INPUT_ERROR);
    return COLVARS_ERROR;
  }

  atoms_ids.reserve(atoms_ids.size()+index_groups_i->size());

  if (b_scalable) {
    for (size_t i = 0; i < index_groups_i->size(); i++) {
      add_atom_id((cvm::proxy)->check_atom_id((*index_groups_i)[i]));
    }
  } else {
    atoms.reserve(atoms.size()+index_groups_i->size());
    for (size_t i = 0; i < index_groups_i->size(); i++) {
      add_atom(cvm::atom((*index_groups_i)[i]));
    }
  }

  if (cvm::get_error())
    return COLVARS_ERROR;

  return COLVARS_OK;
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

      if (b_scalable) {
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
               range_conf+"\".\n", INPUT_ERROR);
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

      if (b_scalable) {
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


int cvm::atom_group::parse_fitting_options(std::string const &group_conf)
{
  bool b_defined_center = get_keyval(group_conf, "centerReference", b_center, false);
  bool b_defined_rotate = get_keyval(group_conf, "rotateReference", b_rotate, false);
  // is the user setting explicit options?
  b_user_defined_fit = b_defined_center || b_defined_rotate;

  get_keyval(group_conf, "enableFitGradients", b_fit_gradients, true);

  if (b_center || b_rotate) {

    if (b_dummy)
      cvm::error("Error: centerReference or rotateReference "
                 "cannot be defined for a dummy atom.\n");

    if (key_lookup(group_conf, "refPositionsGroup")) {
      // instead of this group, define another group to compute the fit
      if (ref_pos_group) {
        cvm::error("Error: the atom group \""+
                   key+"\" has already a reference group "
                   "for the rototranslational fit, which was communicated by the "
                   "colvar component.  You should not use refPositionsGroup "
                   "in this case.\n");
      }
      cvm::log("Within atom group \""+key+"\":\n");
      ref_pos_group = new atom_group(group_conf, "refPositionsGroup");

      // regardless of the configuration, fit gradients must be calculated by refPositionsGroup
      ref_pos_group->b_fit_gradients = this->b_fit_gradients;
      this->b_fit_gradients = false;
    }

    atom_group *group_for_fit = ref_pos_group ? ref_pos_group : this;

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
        if (found && ref_pos_col_value == 0.0)
          cvm::error("Error: refPositionsColValue, "
                     "if provided, must be non-zero.\n", INPUT_ERROR);
      } else {
        // if not, rely on existing atom indices for the group
        group_for_fit->create_sorted_ids();
        ref_pos.resize(group_for_fit->size());
      }

      cvm::load_coords(ref_pos_file.c_str(), ref_pos, group_for_fit->sorted_ids,
                       ref_pos_col, ref_pos_col_value);
    }

    if (ref_pos.size()) {

      if (b_rotate) {
        if (ref_pos.size() != group_for_fit->size())
          cvm::error("Error: the number of reference positions provided("+
                     cvm::to_str(ref_pos.size())+
                     ") does not match the number of atoms within \""+
                     key+
                     "\" ("+cvm::to_str(group_for_fit->size())+
                     "): to perform a rotational fit, "+
                     "these numbers should be equal.\n", INPUT_ERROR);
      }

      // save the center of geometry of ref_pos and subtract it
      center_ref_pos();

    } else {
      cvm::error("Error: no reference positions provided.\n", INPUT_ERROR);
      return COLVARS_ERROR;
    }

    if (b_fit_gradients) {
      group_for_fit->fit_gradients.assign(group_for_fit->size(), cvm::atom_pos(0.0, 0.0, 0.0));
      rot.request_group1_gradients(group_for_fit->size());
    }

    if (b_rotate && !noforce) {
      cvm::log("Warning: atom group \""+key+
               "\" will be aligned to a fixed orientation given by the reference positions provided.  "
               "If the internal structure of the group changes too much (i.e. its RMSD is comparable "
               "to its radius of gyration), the optimal rotation and its gradients may become discontinuous.  "
               "If that happens, use refPositionsGroup (or a different definition for it if already defined) "
               "to align the coordinates.\n");
      // initialize rot member data
      rot.request_group1_gradients(this->size());
    }
  }

  return COLVARS_OK;
}


int cvm::atom_group::create_sorted_ids(void)
{
  // Only do the work if the vector is not yet populated
  if (sorted_ids.size())
    return COLVARS_OK;

  std::list<int> temp_id_list;
  for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    temp_id_list.push_back(ai->id);
  }
  temp_id_list.sort();
  temp_id_list.unique();
  if (temp_id_list.size() != this->size()) {
    cvm::error("Error: duplicate atom IDs in atom group? (found " +
               cvm::to_str(temp_id_list.size()) +
               " unique atom IDs instead of" +
               cvm::to_str(this->size()) + ").\n");
    return COLVARS_ERROR;
  }
  sorted_ids = std::vector<int> (temp_id_list.size());
  unsigned int id_i = 0;
  std::list<int>::iterator li;
  for (li = temp_id_list.begin(); li != temp_id_list.end(); ++li) {
    sorted_ids[id_i] = *li;
    id_i++;
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
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

  if (ref_pos_group)
    ref_pos_group->read_positions();
}


int cvm::atom_group::calc_required_properties()
{
  if (b_dummy) return COLVARS_OK;

  // TODO check if the com is needed?
  calc_center_of_mass();
  if (!b_scalable) {
    // TODO check if calc_center_of_geometry() is needed without a fit?
    calc_center_of_geometry();
    if (b_center || b_rotate) {
      calc_apply_roto_translation();
    }
  }

  // TODO calculate elements of scalable cvc's here before reduction

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

void cvm::atom_group::calc_apply_roto_translation()
{
  atom_group *fit_group = ref_pos_group ? ref_pos_group : this;

  if (b_center) {
    // center on the origin first
    cvm::atom_pos const cog = fit_group->center_of_geometry();
    apply_translation(-1.0 * cog);
  }

  if (b_rotate) {
    // rotate the group (around the center of geometry if b_center is
    // true, around the origin otherwise)
    rot.calc_optimal_rotation(fit_group->positions(), ref_pos);

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->pos = rot.rotate(ai->pos);
    }
  }

  if (b_center) {
    // align with the center of geometry of ref_pos
    apply_translation(ref_pos_cog);
    // update the center of geometry for external use
    cog = ref_pos_cog;
  }

  // recalculate the center of mass
  calc_center_of_mass();
}

void cvm::atom_group::apply_translation(cvm::rvector const &t)
{
  if (b_dummy) {
    cvm::error("Error: cannot translate the coordinates of a dummy atom group.\n", INPUT_ERROR);
    return;
  }

  if (b_scalable) {
    cvm::error("Error: cannot translate the coordinates of a scalable atom group.\n", INPUT_ERROR);
    return;
  }

  for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    ai->pos += t;
  }
}

void cvm::atom_group::apply_rotation(cvm::rotation const &rot)
{
  if (b_dummy) {
    cvm::error("Error: cannot rotate the coordinates of a dummy atom group.\n", INPUT_ERROR);
    return;
  }

  if (b_scalable) {
    cvm::error("Error: cannot rotate the coordinates of a scalable atom group.\n", INPUT_ERROR);
    return;
  }

  for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    ai->pos = rot.rotate(ai->pos - center_of_geometry()) + center_of_geometry();
  }
}


void cvm::atom_group::read_velocities()
{
  if (b_dummy) return;

  if (b_rotate) {

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

void cvm::atom_group::read_system_forces()
{
  if (b_dummy) return;

  if (b_rotate) {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->read_system_force();
      ai->system_force = rot.rotate(ai->system_force);
    }

  } else {

    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
      ai->read_system_force();
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
    cog /= this->size();
  }
  return COLVARS_OK;
}


int cvm::atom_group::calc_center_of_mass()
{
  if (b_dummy) {
    com = dummy_atom_pos;
  } else if (b_scalable) {
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


int cvm::atom_group::calc_dipole(cvm::atom_pos const &com)
{
  if (b_dummy) {
    cvm::error("Error: trying to compute the dipole of an empty group.\n", INPUT_ERROR);
    return COLVARS_ERROR;
  }
  dip.reset();
  for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
    dip += ai->charge * (ai->pos - com);
  }
  return COLVARS_OK;
}


void cvm::atom_group::set_weighted_gradient(cvm::rvector const &grad)
{
  if (b_dummy) return;

  if (b_scalable) {
    scalar_com_gradient = grad;
    return;
  }

  for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    ai->grad = (ai->mass/total_mass) * grad;
  }
}


void cvm::atom_group::calc_fit_gradients()
{
  if (b_dummy) return;

  if ((!b_center) && (!b_rotate)) return; // no fit

  if (cvm::debug())
    cvm::log("Calculating fit gradients.\n");

  atom_group *group_for_fit = ref_pos_group ? ref_pos_group : this;
  group_for_fit->fit_gradients.assign(group_for_fit->size(), cvm::rvector(0.0, 0.0, 0.0));

  if (b_center) {
    // add the center of geometry contribution to the gradients
    for (size_t i = 0; i < this->size(); i++) {
      // need to bring the gradients in original frame first
      cvm::rvector const atom_grad = b_rotate ?
        (rot.inverse()).rotate(atoms[i].grad) :
        atoms[i].grad;
      for (size_t j = 0; j < group_for_fit->size(); j++) {
        group_for_fit->fit_gradients[j] +=
          (-1.0)/(cvm::real(group_for_fit->size())) *
          atom_grad;
      }
    }
  }

  if (b_rotate) {

    // add the rotation matrix contribution to the gradients
    cvm::rotation const rot_inv = rot.inverse();
    cvm::atom_pos const cog = center_of_geometry();

    for (size_t i = 0; i < this->size(); i++) {

      cvm::atom_pos const pos_orig = rot_inv.rotate((b_center ? (atoms[i].pos - cog) : (atoms[i].pos)));

      for (size_t j = 0; j < group_for_fit->size(); j++) {
        // calculate \partial(R(q) \vec{x}_i)/\partial q) \cdot \partial\xi/\partial\vec{x}_i
        cvm::quaternion const dxdq =
          rot.q.position_derivative_inner(pos_orig, atoms[i].grad);
        // multiply by \cdot {\partial q}/\partial\vec{x}_j and add it to the fit gradients
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
               "from a dummy atom group.\n", INPUT_ERROR);
  }

  if (b_scalable) {
    cvm::error("Error: atomic positions are not available "
               "from a scalable atom group.\n", INPUT_ERROR);
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
               "from a dummy atom group.\n", INPUT_ERROR);
  }

  if (b_scalable) {
    cvm::error("Error: atomic positions are not available "
               "from a scalable atom group.\n", INPUT_ERROR);
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
               "from a dummy atom group.\n", INPUT_ERROR);
  }

  if (b_scalable) {
    cvm::error("Error: atomic velocities are not available "
               "from a scalable atom group.\n", INPUT_ERROR);
  }

  std::vector<cvm::rvector> v(this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator vi = v.begin();
  for ( ; ai != this->end(); vi++, ai++) {
    *vi = ai->vel;
  }
  return v;
}

std::vector<cvm::rvector> cvm::atom_group::system_forces() const
{
  if (b_dummy) {
    cvm::error("Error: system forces are not available "
               "from a dummy atom group.\n", INPUT_ERROR);
  }

  if (b_scalable) {
    cvm::error("Error: atomic system forces are not available "
               "from a scalable atom group.\n", INPUT_ERROR);
  }

  std::vector<cvm::rvector> f(this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator fi = f.begin();
  for ( ; ai != this->end(); ++fi, ++ai) {
    *fi = ai->system_force;
  }
  return f;
}

cvm::rvector cvm::atom_group::system_force() const
{
  if (b_dummy) {
    cvm::error("Error: total system forces are not available "
               "from a dummy atom group.\n", INPUT_ERROR);
  }

  if (b_scalable) {
    return (cvm::proxy)->get_atom_group_system_force(index);
  }

  cvm::rvector f(0.0);
  for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
    f += ai->system_force;
  }
  return f;
}


void cvm::atom_group::apply_colvar_force(cvm::real const &force)
{
  if (cvm::debug()) {
    log("Communicating a colvar force from atom group to the MD engine.\n");
  }

  if (b_dummy)
    return;

  if (noforce) {
    cvm::error("Error: sending a force to a group that has "
               "\"enableForces\" set to off.\n");
    return;
  }

  if (b_scalable) {
    (cvm::proxy)->apply_atom_group_force(index, force * scalar_com_gradient);
    return;
  }

  if (b_rotate) {

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

  if ((b_center || b_rotate) && b_fit_gradients) {

    atom_group *group_for_fit = ref_pos_group ? ref_pos_group : this;

    // add the contribution from the roto-translational fit to the gradients
    if (b_rotate) {
      // rotate forces back to the original frame
      cvm::rotation const rot_inv = rot.inverse();
      for (size_t j = 0; j < group_for_fit->size(); j++) {
        (*group_for_fit)[j].apply_force(rot_inv.rotate(force * group_for_fit->fit_gradients[j]));
      }
    } else {
      for (size_t j = 0; j < group_for_fit->size(); j++) {
        (*group_for_fit)[j].apply_force(force * group_for_fit->fit_gradients[j]);
      }
    }
  }

}


void cvm::atom_group::apply_force(cvm::rvector const &force)
{
  if (b_dummy)
    return;

  if (noforce) {
    cvm::error("Error: sending a force to a group that has "
               "\"disableForces\" defined.\n");
    return;
  }

  if (b_scalable) {
    (cvm::proxy)->apply_atom_group_force(index, force);
  }

  if (b_rotate) {

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


void cvm::atom_group::apply_forces(std::vector<cvm::rvector> const &forces)
{
  if (b_dummy)
    return;

  if (noforce)
    cvm::error("Error: sending a force to a group that has "
               "\"disableForces\" defined.\n");

  if (forces.size() != this->size()) {
    cvm::error("Error: trying to apply an array of forces to an atom "
               "group which does not have the same length.\n");
  }

  if (b_rotate) {

    cvm::rotation const rot_inv = rot.inverse();
    cvm::atom_iter ai = this->begin();
    std::vector<cvm::rvector>::const_iterator fi = forces.begin();
    for ( ; ai != this->end(); ++fi, ++ai) {
      ai->apply_force(rot_inv.rotate(*fi));
    }

  } else {

    cvm::atom_iter ai = this->begin();
    std::vector<cvm::rvector>::const_iterator fi = forces.begin();
    for ( ; ai != this->end(); ++fi, ++ai) {
      ai->apply_force(*fi);
    }
  }
}


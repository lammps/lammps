#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvaratoms.h"


// member functions of the "atom" class depend tightly on the MD interface, and are
// thus defined in colvarproxy_xxx.cpp

// in this file only atom_group functions are defined


// Note: "conf" is the configuration of the cvc who is using this atom group;
// "key" is the name of the atom group (e.g. "atoms", "group1", "group2", ...)
cvm::atom_group::atom_group (std::string const &conf,
                             char const        *key)
  : b_center (false), b_rotate (false), b_user_defined_fit (false),
    b_fit_gradients (false),
    ref_pos_group (NULL),
    noforce (false)
{
  cvm::log ("Defining atom group \""+
            std::string (key)+"\".\n");
  // real work is done by parse
  parse (conf, key);
}


cvm::atom_group::atom_group (std::vector<cvm::atom> const &atoms)
  : b_dummy (false), b_center (false), b_rotate (false),
    b_fit_gradients (false), ref_pos_group (NULL),
    noforce (false)
{
  this->reserve (atoms.size());
  for (size_t i = 0; i < atoms.size(); i++) {
    this->push_back (atoms[i]);
  }
  total_mass = 0.0;
  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    total_mass += ai->mass;
  }
}


cvm::atom_group::atom_group()
  : b_dummy (false), b_center (false), b_rotate (false),
    b_fit_gradients (false), ref_pos_group (NULL),
    noforce (false)
{
  total_mass = 0.0;
}


cvm::atom_group::~atom_group()
{
  if (ref_pos_group) {
    delete ref_pos_group;
    ref_pos_group = NULL;
  }
}


void cvm::atom_group::add_atom (cvm::atom const &a)
{
  if (b_dummy) {
    cvm::fatal_error ("Error: cannot add atoms to a dummy group.\n");
  } else {
    this->push_back (a);
    total_mass += a.mass;
  }
}


void cvm::atom_group::reset_mass(std::string &name, int i, int j)
{
  total_mass = 0.0;
  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    total_mass += ai->mass;
  }
  cvm::log ("Re-initialized atom group "+name+":"+cvm::to_str (i)+"/"+
            cvm::to_str (j)+". "+ cvm::to_str (this->size())+
            " atoms: total mass = "+cvm::to_str (this->total_mass)+".\n");
}

void cvm::atom_group::parse (std::string const &conf,
                             char const        *key)
{
  std::string group_conf;

  // save_delimiters is set to false for this call, because "conf" is
  // not the config string of this group, but of its parent object
  // (which has already taken care of the delimiters)
  save_delimiters = false;
  key_lookup (conf, key, group_conf, dummy_pos);
  // restoring the normal value, because we do want keywords checked
  // inside "group_conf"
  save_delimiters = true;

  if (group_conf.size() == 0) {
    cvm::fatal_error ("Error: atom group \""+
                      std::string (key)+"\" is set, but "
                      "has no definition.\n");
  }

  cvm::increase_depth();

  cvm::log ("Initializing atom group \""+std::string (key)+"\".\n");

  // whether or not to include messages in the log
  // colvarparse::Parse_Mode mode = parse_silent;
  // {
  //   bool b_verbose;
  //   get_keyval (group_conf, "verboseOutput", b_verbose, false, parse_silent);
  //   if (b_verbose) mode = parse_normal;
  // }
  colvarparse::Parse_Mode mode = parse_normal;

  {
    //    std::vector<int> atom_indexes;
    std::string numbers_conf = "";
    size_t pos = 0;
    std::vector<int> atom_indexes;
    while (key_lookup (group_conf, "atomNumbers", numbers_conf, pos)) {
      if (numbers_conf.size()) {
        std::istringstream is (numbers_conf);
        int ia;
        while (is >> ia) {
          atom_indexes.push_back (ia);
        }
      }

      if (atom_indexes.size()) {
        this->reserve (this->size()+atom_indexes.size());
        for (size_t i = 0; i < atom_indexes.size(); i++) {
          this->push_back (cvm::atom (atom_indexes[i]));
        }
      } else
        cvm::fatal_error ("Error: no numbers provided for \""
                          "atomNumbers\".\n");

      atom_indexes.clear();
    }

    std::string index_group_name;
    if (get_keyval (group_conf, "indexGroup", index_group_name)) {
      // use an index group from the index file read globally
      std::list<std::string>::iterator names_i = cvm::index_group_names.begin();
      std::list<std::vector<int> >::iterator index_groups_i = cvm::index_groups.begin();
      for ( ; names_i != cvm::index_group_names.end() ; names_i++, index_groups_i++) {
        if (*names_i == index_group_name)
          break;
      }
      if (names_i == cvm::index_group_names.end()) {
        cvm::fatal_error ("Error: could not find index group "+
                          index_group_name+" among those provided by the index file.\n");
      }
      this->reserve (index_groups_i->size());
      for (size_t i = 0; i < index_groups_i->size(); i++) {
        this->push_back (cvm::atom ((*index_groups_i)[i]));
      }
    }
  }

  {
    std::string range_conf = "";
    size_t pos = 0;
    while (key_lookup (group_conf, "atomNumbersRange",
                       range_conf, pos)) {

      if (range_conf.size()) {
        std::istringstream is (range_conf);
        int initial, final;
        char dash;
        if ( (is >> initial) && (initial > 0) &&
             (is >> dash) && (dash == '-') &&
             (is >> final) && (final > 0) ) {
          for (int anum = initial; anum <= final; anum++) {
            this->push_back (cvm::atom (anum));
          }
          range_conf = "";
          continue;
        }
      }

      cvm::fatal_error ("Error: no valid definition for \""
                        "atomNumbersRange\", \""+
                        range_conf+"\".\n");
    }
  }

  std::vector<std::string> psf_segids;
  get_keyval (group_conf, "psfSegID", psf_segids, std::vector<std::string> (), mode);
  for (std::vector<std::string>::iterator psii = psf_segids.begin();
       psii < psf_segids.end(); psii++) {

    if ( (psii->size() == 0) || (psii->size() > 4) ) {
      cvm::fatal_error ("Error: invalid segmend identifier provided, \""+
                        (*psii)+"\".\n");
    }
  }

  {
    std::string range_conf = "";
    size_t pos = 0;
    size_t range_count = 0;
    std::vector<std::string>::iterator psii = psf_segids.begin();
    while (key_lookup (group_conf, "atomNameResidueRange",
                       range_conf, pos)) {
      range_count++;

      if (range_count > psf_segids.size()) {
        cvm::fatal_error ("Error: more instances of \"atomNameResidueRange\" than "
                          "values of \"psfSegID\".\n");
      }

      std::string const &psf_segid = psf_segids.size() ? *psii : std::string ("");

      if (range_conf.size()) {

        std::istringstream is (range_conf);
        std::string atom_name;
        int initial, final;
        char dash;
        if ( (is >> atom_name) && (atom_name.size())  &&
             (is >> initial) && (initial > 0) &&
             (is >> dash) && (dash == '-') &&
             (is >> final) && (final > 0) ) {
          for (int resid = initial; resid <= final; resid++) {
            this->push_back (cvm::atom (resid, atom_name, psf_segid));
          }
          range_conf = "";
        } else {
          cvm::fatal_error ("Error: cannot parse definition for \""
                            "atomNameResidueRange\", \""+
                            range_conf+"\".\n");
        }

      } else {
        cvm::fatal_error ("Error: atomNameResidueRange with empty definition.\n");
      }

      if (psf_segid.size())
        psii++;
    }
  }

  {
    // read the atoms from a file
    std::string atoms_file_name;
    if (get_keyval (group_conf, "atomsFile", atoms_file_name, std::string (""), mode)) {

      std::string atoms_col;
      if (!get_keyval (group_conf, "atomsCol", atoms_col, std::string (""), mode)) {
        cvm::fatal_error ("Error: parameter atomsCol is required if atomsFile is set.\n");
      }

      double atoms_col_value;
      bool const atoms_col_value_defined = get_keyval (group_conf, "atomsColValue", atoms_col_value, 0.0, mode);
      if (atoms_col_value_defined && (!atoms_col_value))
        cvm::fatal_error ("Error: atomsColValue, "
                          "if provided, must be non-zero.\n");

      cvm::load_atoms (atoms_file_name.c_str(), *this, atoms_col, atoms_col_value);
    }
  }

  for (std::vector<cvm::atom>::iterator a1 = this->begin();
       a1 != this->end(); a1++) {
    std::vector<cvm::atom>::iterator a2 = a1;
    ++a2;
    for ( ; a2 != this->end(); a2++) {
      if (a1->id == a2->id) {
        if (cvm::debug())
          cvm::log ("Discarding doubly counted atom with number "+
                    cvm::to_str (a1->id+1)+".\n");
        a2 = this->erase (a2);
        if (a2 == this->end())
          break;
      }
    }
  }

  if (get_keyval (group_conf, "dummyAtom", dummy_atom_pos, cvm::atom_pos(), mode)) {
    b_dummy = true;
    this->total_mass = 1.0;
  } else
    b_dummy = false;

  if (b_dummy && (this->size()))
    cvm::fatal_error ("Error: cannot set up group \""+
                      std::string (key)+"\" as a dummy atom "
                      "and provide it with atom definitions.\n");

#if (! defined (COLVARS_STANDALONE))
  if ( (!b_dummy) && (!cvm::b_analysis) && (!(this->size())) ) {
    cvm::fatal_error ("Error: no atoms defined for atom group \""+
                      std::string (key)+"\".\n");
  }
#endif

  if (!b_dummy) {
    this->total_mass = 0.0;
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      this->total_mass += ai->mass;
    }
  }

  if (!b_dummy) {
    bool enable_forces = true;
    // disableForces is deprecated
    if (get_keyval (group_conf, "enableForces", enable_forces, true, mode)) {
      noforce = !enable_forces;
    } else {
      get_keyval (group_conf, "disableForces", noforce, false, mode);
    }
  }

  // FITTING OPTIONS

  bool b_defined_center = get_keyval (group_conf, "centerReference", b_center, false, mode);
  bool b_defined_rotate = get_keyval (group_conf, "rotateReference", b_rotate, false, mode);
  // is the user setting explicit options?
  b_user_defined_fit = b_defined_center || b_defined_rotate;

  get_keyval (group_conf, "enableFitGradients", b_fit_gradients, true, mode);

  if (b_center || b_rotate) {

    if (b_dummy)
      cvm::fatal_error ("Error: centerReference or rotateReference "
                        "cannot be defined for a dummy atom.\n");

    if (key_lookup (group_conf, "refPositionsGroup")) {
      // instead of this group, define another group to compute the fit
      if (ref_pos_group) {
        cvm::fatal_error ("Error: the atom group \""+
                          std::string (key)+"\" has already a reference group "
                          "for the rototranslational fit, which was communicated by the "
                          "colvar component.  You should not use refPositionsGroup "
                          "in this case.\n");
      }
      cvm::log ("Within atom group \""+std::string (key)+"\":\n");
      ref_pos_group = new atom_group (group_conf, "refPositionsGroup");

      // regardless of the configuration, fit gradients must be calculated by refPositionsGroup
      ref_pos_group->b_fit_gradients = this->b_fit_gradients;
      this->b_fit_gradients = false;
    }

    atom_group *group_for_fit = ref_pos_group ? ref_pos_group : this;

    get_keyval (group_conf, "refPositions", ref_pos, ref_pos, mode);

    std::string ref_pos_file;
    if (get_keyval (group_conf, "refPositionsFile", ref_pos_file, std::string (""), mode)) {

      if (ref_pos.size()) {
        cvm::fatal_error ("Error: cannot specify \"refPositionsFile\" and "
                          "\"refPositions\" at the same time.\n");
      }

      std::string ref_pos_col;
      double ref_pos_col_value;

      if (get_keyval (group_conf, "refPositionsCol", ref_pos_col, std::string (""), mode)) {
        // if provided, use PDB column to select coordinates
        bool found = get_keyval (group_conf, "refPositionsColValue", ref_pos_col_value, 0.0, mode);
        if (found && !ref_pos_col_value)
          cvm::fatal_error ("Error: refPositionsColValue, "
                            "if provided, must be non-zero.\n");
      } else {
        // if not, rely on existing atom indices for the group
        group_for_fit->create_sorted_ids();
        ref_pos.resize (group_for_fit->size());
      }

      cvm::load_coords (ref_pos_file.c_str(), ref_pos, group_for_fit->sorted_ids,
                        ref_pos_col, ref_pos_col_value);
    }

    if (ref_pos.size()) {

      if (b_rotate) {
        if (ref_pos.size() != group_for_fit->size())
          cvm::fatal_error ("Error: the number of reference positions provided ("+
                            cvm::to_str (ref_pos.size())+
                            ") does not match the number of atoms within \""+
                            std::string (key)+
                            "\" ("+cvm::to_str (group_for_fit->size())+
                            "): to perform a rotational fit, "+
                            "these numbers should be equal.\n");
      }

      // save the center of geometry of ref_pos and subtract it
      center_ref_pos();

    } else {
#if (! defined (COLVARS_STANDALONE))
      cvm::fatal_error ("Error: no reference positions provided.\n");
#endif
    }

    if (b_fit_gradients) {
      group_for_fit->fit_gradients.assign (group_for_fit->size(), cvm::atom_pos (0.0, 0.0, 0.0));
      rot.request_group1_gradients (group_for_fit->size());
    }

    if (b_rotate && !noforce) {
      cvm::log ("Warning: atom group \""+std::string (key)+
                "\" will be aligned to a fixed orientation given by the reference positions provided.  "
                "If the internal structure of the group changes too much (i.e. its RMSD is comparable "
                "to its radius of gyration), the optimal rotation and its gradients may become discontinuous.  "
                "If that happens, use refPositionsGroup (or a different definition for it if already defined) "
                "to align the coordinates.\n");
      // initialize rot member data
      rot.request_group1_gradients (this->size());
    }

  }

  if (cvm::debug())
    cvm::log ("Done initializing atom group with name \""+
              std::string (key)+"\".\n");

  this->check_keywords (group_conf, key);

  cvm::log ("Atom group \""+std::string (key)+"\" defined, "+
            cvm::to_str (this->size())+" atoms initialized: total mass = "+
            cvm::to_str (this->total_mass)+".\n");

  cvm::decrease_depth();
}


void cvm::atom_group::create_sorted_ids (void)
{
  // Only do the work if the vector is not yet populated
  if (sorted_ids.size())
    return;

  std::list<int> temp_id_list;
  for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++) {
    temp_id_list.push_back (ai->id);
  }
  temp_id_list.sort();
  temp_id_list.unique();
  if (temp_id_list.size() != this->size()) {
    cvm::fatal_error ("Error: duplicate atom IDs in atom group? (found " +
                      cvm::to_str(temp_id_list.size()) +
                      " unique atom IDs instead of" +
                      cvm::to_str(this->size()) + ").\n");
  }
  sorted_ids = std::vector<int> (temp_id_list.begin(), temp_id_list.end());
  return;
}

void cvm::atom_group::center_ref_pos()
{
  // save the center of geometry of ref_pos and then subtract it from
  // them; in this way it will be possible to use ref_pos also for
  // the rotational fit
  // This is called either by atom_group::parse or by CVCs that set
  // reference positions (eg. RMSD, eigenvector)
  ref_pos_cog = cvm::atom_pos (0.0, 0.0, 0.0);
  std::vector<cvm::atom_pos>::iterator pi = ref_pos.begin();
  for ( ; pi != ref_pos.end(); pi++) {
    ref_pos_cog += *pi;
  }
  ref_pos_cog /= (cvm::real) ref_pos.size();
  for (std::vector<cvm::atom_pos>::iterator pi = ref_pos.begin();
       pi != ref_pos.end(); pi++) {
    (*pi) -= ref_pos_cog;
  }
}

void cvm::atom_group::read_positions()
{
  if (b_dummy) return;

  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    ai->read_position();
  }

  if (ref_pos_group)
    ref_pos_group->read_positions();
}

void cvm::atom_group::calc_apply_roto_translation()
{
  atom_group *fit_group = ref_pos_group ? ref_pos_group : this;

  if (b_center) {
    // center on the origin first
    cvm::atom_pos const cog = fit_group->center_of_geometry();
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->pos -= cog;
    }
  }

  if (b_rotate) {
    // rotate the group (around the center of geometry if b_center is
    // true, around the origin otherwise)
    rot.calc_optimal_rotation (fit_group->positions(), ref_pos);

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->pos = rot.rotate (ai->pos);
    }
  }

  if (b_center) {
    // align with the center of geometry of ref_pos
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->pos += ref_pos_cog;
    }
  }
}

void cvm::atom_group::apply_translation (cvm::rvector const &t)
{
  if (b_dummy) return;

  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    ai->pos += t;
  }
}

void cvm::atom_group::apply_rotation (cvm::rotation const &rot)
{
  if (b_dummy) return;

  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    ai->pos = rot.rotate (ai->pos);
  }
}


void cvm::atom_group::read_velocities()
{
  if (b_dummy) return;

  if (b_rotate) {

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->read_velocity();
      ai->vel = rot.rotate (ai->vel);
    }

  } else {

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->read_velocity();
    }
  }
}

void cvm::atom_group::read_system_forces()
{
  if (b_dummy) return;

  if (b_rotate) {

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->read_system_force();
      ai->system_force = rot.rotate (ai->system_force);
    }

  } else {

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->read_system_force();
    }
  }
}

cvm::atom_pos cvm::atom_group::center_of_geometry() const
{
  if (b_dummy)
    return dummy_atom_pos;

  cvm::atom_pos cog (0.0, 0.0, 0.0);
  for (cvm::atom_const_iter ai = this->begin();
       ai != this->end(); ai++) {
    cog += ai->pos;
  }
  cog /= this->size();
  return cog;
}

cvm::atom_pos cvm::atom_group::center_of_mass() const
{
  if (b_dummy)
    return dummy_atom_pos;

  cvm::atom_pos com (0.0, 0.0, 0.0);
  for (cvm::atom_const_iter ai = this->begin();
       ai != this->end(); ai++) {
    com += ai->mass * ai->pos;
  }
  com /= this->total_mass;
  return com;
}


void cvm::atom_group::set_weighted_gradient (cvm::rvector const &grad)
{
  if (b_dummy) return;

  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    ai->grad = (ai->mass/this->total_mass) * grad;
  }
}


void cvm::atom_group::calc_fit_gradients()
{
  if (b_dummy) return;

  if ((!b_center) && (!b_rotate)) return; // no fit

  if (cvm::debug())
    cvm::log ("Calculating fit gradients.\n");

  atom_group *group_for_fit = ref_pos_group ? ref_pos_group : this;
  group_for_fit->fit_gradients.assign (group_for_fit->size(), cvm::rvector (0.0, 0.0, 0.0));

  if (b_center) {
    // add the center of geometry contribution to the gradients
    for (size_t i = 0; i < this->size(); i++) {
      // need to bring the gradients in original frame first
      cvm::rvector const atom_grad = b_rotate ?
        (rot.inverse()).rotate ((*this)[i].grad) :
        (*this)[i].grad;
      for (size_t j = 0; j < group_for_fit->size(); j++) {
        group_for_fit->fit_gradients[j] +=
          (-1.0)/(cvm::real (group_for_fit->size())) *
          atom_grad;
      }
    }
  }

  if (b_rotate) {

    // add the rotation matrix contribution to the gradients
    cvm::rotation const rot_inv = rot.inverse();
    cvm::atom_pos const cog = this->center_of_geometry();

    for (size_t i = 0; i < this->size(); i++) {

      cvm::atom_pos const pos_orig = rot_inv.rotate ((b_center ? ((*this)[i].pos - cog) : ((*this)[i].pos)));

      for (size_t j = 0; j < group_for_fit->size(); j++) {
        // calculate \partial(R(q) \vec{x}_i)/\partial q) \cdot \partial\xi/\partial\vec{x}_i
        cvm::quaternion const dxdq =
          rot.q.position_derivative_inner (pos_orig, (*this)[i].grad);
        // multiply by \cdot {\partial q}/\partial\vec{x}_j and add it to the fit gradients
        for (size_t iq = 0; iq < 4; iq++) {
          group_for_fit->fit_gradients[j] += dxdq[iq] * rot.dQ0_1[j][iq];
        }
      }
    }
  }
  if (cvm::debug())
    cvm::log ("Done calculating fit gradients.\n");
}


std::vector<cvm::atom_pos> cvm::atom_group::positions() const
{
  if (b_dummy)
    cvm::fatal_error ("Error: positions are not available "
                      "from a dummy atom group.\n");

  std::vector<cvm::atom_pos> x (this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator xi = x.begin();
  for ( ; ai != this->end(); xi++, ai++) {
    *xi = ai->pos;
  }
  return x;
}

std::vector<cvm::atom_pos> cvm::atom_group::positions_shifted (cvm::rvector const &shift) const
{
  if (b_dummy)
    cvm::fatal_error ("Error: positions are not available "
                      "from a dummy atom group.\n");

  std::vector<cvm::atom_pos> x (this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator xi = x.begin();
  for ( ; ai != this->end(); xi++, ai++) {
    *xi = (ai->pos + shift);
  }
  return x;
}

std::vector<cvm::rvector> cvm::atom_group::velocities() const
{
  if (b_dummy)
    cvm::fatal_error ("Error: velocities are not available "
                      "from a dummy atom group.\n");

  std::vector<cvm::rvector> v (this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator vi = v.begin();
  for ( ; ai != this->end(); vi++, ai++) {
    *vi = ai->vel;
  }
  return v;
}

std::vector<cvm::rvector> cvm::atom_group::system_forces() const
{
  if (b_dummy)
    cvm::fatal_error ("Error: system forces are not available "
                      "from a dummy atom group.\n");

  std::vector<cvm::rvector> f (this->size(), 0.0);
  cvm::atom_const_iter ai = this->begin();
  std::vector<cvm::atom_pos>::iterator fi = f.begin();
  for ( ; ai != this->end(); fi++, ai++) {
    *fi = ai->system_force;
  }
  return f;
}

cvm::rvector cvm::atom_group::system_force() const
{
  if (b_dummy)
    cvm::fatal_error ("Error: system forces are not available "
                      "from a dummy atom group.\n");

  cvm::rvector f (0.0);
  for (cvm::atom_const_iter ai = this->begin(); ai != this->end(); ai++) {
    f += ai->system_force;
  }
  return f;
}


void cvm::atom_group::apply_colvar_force (cvm::real const &force)
{
  if (b_dummy)
    return;

  if (noforce)
    cvm::fatal_error ("Error: sending a force to a group that has "
                      "\"enableForces\" set to off.\n");

  if (b_rotate) {

    // rotate forces back to the original frame
    cvm::rotation const rot_inv = rot.inverse();
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->apply_force (rot_inv.rotate (force * ai->grad));
    }

  } else {

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->apply_force (force * ai->grad);
    }
  }

  if ((b_center || b_rotate) && b_fit_gradients) {

    atom_group *group_for_fit = ref_pos_group ? ref_pos_group : this;

    // add the contribution from the roto-translational fit to the gradients
    if (b_rotate) {
      // rotate forces back to the original frame
      cvm::rotation const rot_inv = rot.inverse();
      for (size_t j = 0; j < group_for_fit->size(); j++) {
        (*group_for_fit)[j].apply_force (rot_inv.rotate (force * group_for_fit->fit_gradients[j]));
      }
    } else {
      for (size_t j = 0; j < group_for_fit->size(); j++) {
        (*group_for_fit)[j].apply_force (force * group_for_fit->fit_gradients[j]);
      }
    }
  }

}


void cvm::atom_group::apply_force (cvm::rvector const &force)
{
  if (b_dummy)
    return;

  if (noforce)
    cvm::fatal_error ("Error: sending a force to a group that has "
                      "\"disableForces\" defined.\n");

  if (b_rotate) {

    cvm::rotation const rot_inv = rot.inverse();
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->apply_force (rot_inv.rotate ((ai->mass/this->total_mass) * force));
    }

  } else {

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->apply_force ((ai->mass/this->total_mass) * force);
    }
  }
}


void cvm::atom_group::apply_forces (std::vector<cvm::rvector> const &forces)
{
  if (b_dummy)
    return;

  if (noforce)
    cvm::fatal_error ("Error: sending a force to a group that has "
                      "\"disableForces\" defined.\n");

  if (forces.size() != this->size()) {
    cvm::fatal_error ("Error: trying to apply an array of forces to an atom "
                      "group which does not have the same length.\n");
  }

  if (b_rotate) {

    cvm::rotation const rot_inv = rot.inverse();
    cvm::atom_iter ai = this->begin();
    std::vector<cvm::rvector>::const_iterator fi = forces.begin();
    for ( ; ai != this->end(); fi++, ai++) {
      ai->apply_force (rot_inv.rotate (*fi));
    }

  } else {

    cvm::atom_iter ai = this->begin();
    std::vector<cvm::rvector>::const_iterator fi = forces.begin();
    for ( ; ai != this->end(); fi++, ai++) {
      ai->apply_force (*fi);
    }
  }
}


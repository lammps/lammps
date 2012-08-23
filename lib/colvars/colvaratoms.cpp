#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvaratoms.h"


// atom member functions depend tightly on the MD interface, and are
// thus defined in colvarproxy_xxx.cpp


// atom_group member functions

cvm::atom_group::atom_group (std::string const &conf,
                             char const        *key,
                             atom_group        *ref_pos_group_in)
  : b_center (false), b_rotate (false),
    ref_pos_group (NULL), // this is always set within parse(),
                          // regardless of ref_pos_group_in
    noforce (false)
{
  cvm::log ("Defining atom group \""+
            std::string (key)+"\".\n");
  parse (conf, key, ref_pos_group_in);
}


void cvm::atom_group::parse (std::string const &conf,
                             char const        *key,
                             atom_group        *ref_pos_group_in)
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
  colvarparse::Parse_Mode mode = parse_silent;
  {
    bool b_verbose;
    get_keyval (group_conf, "verboseOutput", b_verbose, false, parse_silent);
    if (b_verbose) mode = parse_normal;
  }

  {
    // get the atoms by numbers
    std::vector<int> atom_indexes;
    if (get_keyval (group_conf, "atomNumbers", atom_indexes, atom_indexes, mode)) {
      if (atom_indexes.size()) {
        this->reserve (this->size()+atom_indexes.size());
        for (size_t i = 0; i < atom_indexes.size(); i++) {
          this->push_back (cvm::atom (atom_indexes[i]));
        }
      } else
        cvm::fatal_error ("Error: no numbers provided for \""
                          "atomNumbers\".\n");
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

  get_keyval (group_conf, "disableForces",   noforce,  false, mode);

  get_keyval (group_conf, "centerReference", b_center, false, mode);
  get_keyval (group_conf, "rotateReference", b_rotate, false, mode);

  if (b_center || b_rotate) {

    if (b_dummy)
      cvm::fatal_error ("Error: cannot set \"centerReference\" or "
                        "\"rotateReference\" with \"dummyAtom\".\n");

    // use refPositionsGroup instead of this group as the one which is
    // used to fit the coordinates
    if (key_lookup (group_conf, "refPositionsGroup")) {
      if (ref_pos_group) {
        cvm::fatal_error ("Error: the atom group \""+
                          std::string (key)+"\" has already a reference group "
                          "for the rototranslational fit, which was communicated by the "
                          "colvar component.  You should not use refPositionsGroup "
                          "in this case.\n");
      }
      cvm::log ("Within atom group \""+std::string (key)+"\":\n");
      ref_pos_group = new atom_group (group_conf, "refPositionsGroup");
    }

    atom_group *ag = ref_pos_group ? ref_pos_group : this;

    if (get_keyval (group_conf, "refPositions", ref_pos, ref_pos, mode)) {
      cvm::log ("Using reference positions from input file.\n");
      if (ref_pos.size() != ag->size()) {
        cvm::fatal_error ("Error: the number of reference positions provided ("+
                          cvm::to_str (ref_pos.size())+
                          ") does not match the number of atoms within \""+
                          std::string (key)+
                          "\" ("+cvm::to_str (ag->size())+").\n");
      }
    }

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
        ag->create_sorted_ids();
      }
      cvm::load_coords (ref_pos_file.c_str(), ref_pos, ag->sorted_ids,
                        ref_pos_col, ref_pos_col_value);
    }

    if (ref_pos.size()) {

      if (b_rotate) {
        if (ref_pos.size() != ag->size())
          cvm::fatal_error ("Error: the number of reference positions provided ("+
                            cvm::to_str (ref_pos.size())+
                            ") does not match the number of atoms within \""+
                            std::string (key)+
                            "\" ("+cvm::to_str (ag->size())+").\n");
      }

      // save the center of mass of ref_pos and then subtract it from
      // them; in this way it is possible to use the coordinates for
      // the rotational fit, if needed
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
    } else {
#if (! defined (COLVARS_STANDALONE))
      if (!cvm::b_analysis)
        cvm::fatal_error ("Error: no reference positions provided.\n");
#endif
    }

    if (b_rotate && !noforce) {
      cvm::log ("Warning: atom group \""+std::string (key)+
                "\" is set to be rotated to a reference orientation: "
                "a torque different than zero on this group "
                "could make the simulation unstable.  "
                "If this happens, set \"disableForces\" to yes "
                "for this group.\n");
    }

  }


  if (cvm::debug())
    cvm::log ("Done initializing atom group with name \""+
              std::string (key)+"\".\n");

  this->check_keywords (group_conf, key);

  cvm::log ("Atom group \""+std::string (key)+"\" defined, "+
            cvm::to_str (this->size())+" initialized: total mass = "+
            cvm::to_str (this->total_mass)+".\n");

  cvm::decrease_depth();
}


cvm::atom_group::atom_group (std::vector<cvm::atom> const &atoms)
  : b_dummy (false), b_center (false), b_rotate (false), 
    ref_pos_group (NULL), noforce (false)
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
    ref_pos_group (NULL), noforce (false)
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


void cvm::atom_group::read_positions()
{
  if (b_dummy) return;

#if (! defined (COLVARS_STANDALONE))
  if (!this->size())
    cvm::fatal_error ("Error: no atoms defined in the requested group.\n");
#endif

  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    ai->read_position();
  }

  if (ref_pos_group)
    ref_pos_group->read_positions();

  atom_group *fit_group = ref_pos_group ? ref_pos_group : this;

  if (b_center) {
    // store aside the current center of geometry (all positions will be
    // set to the closest images to the first one) and then center on
    // the origin
    cvm::atom_pos const cog = fit_group->center_of_geometry();
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->pos -= cog;
    }
  }

  if (b_rotate) {
    // rotate the group (around the center of geometry if b_center is
    // true, around the origin otherwise); store the rotation, in
    // order to bring back the forces to the original frame before
    // applying them
    rot.calc_optimal_rotation (fit_group->positions(), ref_pos);

    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->pos = rot.rotate (ai->pos);
    }
  }

  if (b_center) {
    // use the center of geometry of ref_pos
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

/* This is deprecated.
cvm::atom_pos cvm::atom_group::center_of_geometry (cvm::atom_pos const &ref_pos)
{
  if (b_dummy) {
    cvm::select_closest_image (dummy_atom_pos, ref_pos);
    return dummy_atom_pos;
  }

  cvm::atom_pos cog (0.0, 0.0, 0.0);
  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    cvm::select_closest_image (ai->pos, ref_pos);
    cog += ai->pos;
  }
  cog /= this->size();
  return cog;
} */

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

/* This is deprecated.
cvm::atom_pos cvm::atom_group::center_of_mass (cvm::atom_pos const &ref_pos)
{
  if (b_dummy) {
    cvm::select_closest_image (dummy_atom_pos, ref_pos);
    return dummy_atom_pos;
  }

  cvm::atom_pos com (0.0, 0.0, 0.0);
  for (cvm::atom_iter ai = this->begin();
       ai != this->end(); ai++) {
    cvm::select_closest_image (ai->pos, ref_pos);
    com += ai->mass * ai->pos;
  }
  com /= this->total_mass;
  return com;
}*/

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
                      "\"disableForces\" defined.\n");

  if (b_rotate) {

    // get the forces back from the rotated frame
    cvm::rotation const rot_inv = rot.inverse();
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->apply_force (rot_inv.rotate (force * ai->grad));
    }

  } else {

    // no need to manipulate gradients, they are still in the original
    // frame
    for (cvm::atom_iter ai = this->begin();
         ai != this->end(); ai++) {
      ai->apply_force (force * ai->grad);
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




#include "mpi.h"
#include "lammps.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "output.h"
#include "random_park.h"

#include "fix_colvars.h"

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_lammps.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <sstream>
#include <string>

#define HASH_FAIL  -1

////////////////////////////////////////////////////////////////////////
// local helper functions

// safely move filename to filename.extension
static void my_backup_file(const char *filename, const char *extension)
{
  struct stat sbuf;
  if (stat(filename, &sbuf) == 0) {
    if ( ! extension ) extension = ".BAK";
    char *backup = new char[strlen(filename)+strlen(extension)+1];
    strcpy(backup, filename);
    strcat(backup, extension);
#if defined(_WIN32) && !defined(__CYGWIN__)
    remove(backup);
#endif
    if ( rename(filename,backup) ) {
      char *sys_err_msg = strerror(errno);
      if ( !sys_err_msg ) sys_err_msg = "(unknown error)";
      fprintf(stderr,"Error renaming file %s to %s: %s\n",
	      filename, backup, sys_err_msg);
    }
    delete [] backup;
  }
}

////////////////////////////////////////////////////////////////////////

colvarproxy_lammps::colvarproxy_lammps(LAMMPS_NS::LAMMPS *lmp,
				       const char *conf_file,
				       const char *inp_name,
				       const char *out_name,
				       const int seed,
				       const double tt,
				       struct commdata *coords,
				       struct commdata *forces,
				       struct commdata *oforce,
				       int (*i)(void *,int), void *m)
  : _lmp(lmp),_coords(coords),_forces(forces),_oforce(oforce),
    _idlookup(i),_idmap(m)
{
  if (cvm::debug())
    log("Info: initializing the colvars proxy object.\n");

  _random = new LAMMPS_NS::RanPark(lmp,seed);

  first_timestep=true;
  system_force_requested=false;
  previous_step=-1;
  t_target = tt;

  // set input restart name and strip the extension, if present
  input_prefix_str = std::string(inp_name ? inp_name : "");
  if (input_prefix_str.rfind(".colvars.state") != std::string::npos)
    input_prefix_str.erase(input_prefix_str.rfind(".colvars.state"),
                            std::string(".colvars.state").size());

  // output prefix is always given
  output_prefix_str = std::string(out_name);
  restart_prefix_str = std::string("rest");

  if (_lmp->output->restart_every > 0) {
    restart_prefix_str = std::string(_lmp->output->restart1);
    
    if (restart_prefix_str.rfind(".*") != std::string::npos)
      restart_prefix_str.erase(restart_prefix_str.rfind(".*"),2);
  }

  // initiate the colvarmodule
  colvars = new colvarmodule(conf_file,this);

  if (cvm::debug()) {
    log("colvars_atoms = "+cvm::to_str(colvars_atoms)+"\n");
    log("colvars_atoms_ncopies = "+cvm::to_str(colvars_atoms_ncopies)+"\n");
    log("positions = "+cvm::to_str(positions)+"\n");
    log("total_forces = "+cvm::to_str(total_forces)+"\n");
    log("applied_forces = "+cvm::to_str(applied_forces)+"\n");
    log(cvm::line_marker);
    log("Info: done initializing the colvars proxy object.\n");
  }
}

colvarproxy_lammps::~colvarproxy_lammps()
{
  delete _random;

  if (colvars != NULL) {
    colvars->write_output_files();
    delete colvars;
    colvars = NULL;
  }
}

// trigger colvars computation 
void colvarproxy_lammps::compute()
{
  if (first_timestep) {
    first_timestep = false;
  } else {
    // Use the time step number inherited from GlobalMaster
    if ( _lmp->update->ntimestep - previous_step == 1 ) {
      colvars->it++;
    }
    // Other cases could mean:
    // - run 0
    // - beginning of a new run statement
    // then the internal counter should not be incremented
  }
  previous_step = _lmp->update->ntimestep;

  if (cvm::debug()) {
    cvm::log (cvm::line_marker+
              "colvarproxy_lammps, step no. "+cvm::to_str(colvars->it)+"\n"+
              "Updating internal data.\n");
  }

  // transfer coordinates and clear forces
  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    int j = (*_idlookup)(_idmap,colvars_atoms[i]);
    if (j >=0)
      positions[i] = cvm::rvector(_coords[j].x,_coords[j].y,_coords[j].z);
  }

  if (system_force_requested && cvm::step_relative() > 0) {
    cvm::log (cvm::line_marker+
              "colvarproxy_lammps, system_force request unsupported.\n");

    // sort the array of total forces from the previous step (but only
    // do it if there *is* a previous step!)
    struct commdata *o = _oforce;
    for (size_t i = 0; i < colvars_atoms.size(); i++, o++)
      total_forces[i] = cvm::rvector (o->x, o->y, o->z);
    
#if 0  /* add test for running under minimization */
    if (!found_total_force)
      cvm::fatal_error ("Error: system forces were requested,"
			" but total force on atom "+
			cvm::to_str (colvars_atoms[i]+1) + " was not\n"
			"found. The most probable cause is combination "
			"of energy minimization with a\n"
			"biasing method that requires MD (e.g. ABF). "
			"Always run minimization\n"
			"and ABF separately.");
#endif
    /* XXX: compute "applied force" and figure out this forces mess here */
  }
  
  // call the collective variable module
  colvars->calc();

#if 0
  // /* add restraint energy to total energy */
  reduction->submit();
#endif
}

cvm::rvector colvarproxy_lammps::position_distance(cvm::atom_pos const &pos1,
						   cvm::atom_pos const &pos2)
{
  double xtmp = pos2.x - pos1.x;
  double ytmp = pos2.y - pos1.y;
  double ztmp = pos2.z - pos1.z;
  _lmp->domain->minimum_image(xtmp,ytmp,ztmp);
  return cvm::rvector (xtmp, ytmp, ztmp);
}

cvm::real colvarproxy_lammps::position_dist2(cvm::atom_pos const &pos1,
					     cvm::atom_pos const &pos2)
{
  double xtmp = pos2.x - pos1.x;
  double ytmp = pos2.y - pos1.y;
  double ztmp = pos2.z - pos1.z;
  _lmp->domain->minimum_image(xtmp,ytmp,ztmp);
  return cvm::real(xtmp*xtmp + ytmp*ytmp + ztmp*ztmp);
}


inline void colvarproxy_lammps::select_closest_image(cvm::atom_pos &pos,
						     cvm::atom_pos const &ref)
{
  double xtmp = pos.x - ref.x;
  double ytmp = pos.y - ref.y;
  double ztmp = pos.z - ref.z;
  _lmp->domain->minimum_image(xtmp,ytmp,ztmp);
  pos.x = ref.x + xtmp;
  pos.y = ref.y + ytmp;
  pos.z = ref.z + ztmp;
}

void colvarproxy_lammps::log(std::string const &message)
{
  std::istringstream is (message);
  std::string line;
  while (std::getline (is, line)) {
    if (_lmp->screen)
      fprintf(_lmp->screen,"colvars: %s\n",line.c_str());
    if (_lmp->logfile)
      fprintf(_lmp->logfile,"colvars: %s\n",line.c_str());
  }
}

void colvarproxy_lammps::fatal_error(std::string const &message)
{
  log(message);
  if (!cvm::debug())
    log ("If this error message is unclear, "
	 "try recompiling the colvars library with -DCOLVARS_DEBUG.\n");
  _lmp->error->one(FLERR,
		   "Fatal error in the collective variables module.\n");
}

void colvarproxy_lammps::exit(std::string const &message)
{
  log(message);
  log("Request to exit the simulation made.\n");
}

enum e_pdb_field {
  e_pdb_none,
  e_pdb_occ,
  e_pdb_beta,
  e_pdb_x,
  e_pdb_y,
  e_pdb_z,
  e_pdb_ntot
};


e_pdb_field pdb_field_str2enum (std::string const &pdb_field_str)
{
  e_pdb_field pdb_field = e_pdb_none;

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("O")) {
    pdb_field = e_pdb_occ;
  }

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("B")) {
    pdb_field = e_pdb_beta;
  }

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("X")) {
    pdb_field = e_pdb_x;
  }
  
  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("Y")) {
    pdb_field = e_pdb_y;
  }

  if (colvarparse::to_lower_cppstr (pdb_field_str) ==
      colvarparse::to_lower_cppstr ("Z")) {
    pdb_field = e_pdb_z;
  }

  if (pdb_field == e_pdb_none) {
    cvm::fatal_error ("Error: unsupported PDB field, \""+
                      pdb_field_str+"\".\n");
  }

  return pdb_field;
}

void colvarproxy_lammps::load_coords (char const *pdb_filename,
                                    std::vector<cvm::atom_pos> &pos,
                                    const std::vector<int> &indices,
                                    std::string const pdb_field_str,
                                    double const pdb_field_value)
{

  cvm::fatal_error ("Reading collective variable coordinates "
		    "from a PDB file is currently not supported.\n");

#if 0
  e_pdb_field pdb_field_index;
  bool const use_pdb_field = (pdb_field_str.size() > 0);
  if (use_pdb_field) {
    pdb_field_index = pdb_field_str2enum (pdb_field_str);
  }

  // next index to be looked up in PDB file (if list is supplied)
  std::vector<int>::const_iterator current_index = indices.begin();

  PDB *pdb = new PDB (pdb_filename);
  size_t const pdb_natoms = pdb->num_atoms();
  
  if (pos.size() != pdb_natoms) {

    bool const pos_allocated = (pos.size() > 0);

    size_t ipos = 0, ipdb = 0;
    for ( ; ipdb < pdb_natoms; ipdb++) {

      if (use_pdb_field) {
        // PDB field mode: skip atoms with wrong value in PDB field
        double atom_pdb_field_value = 0.0;

        switch (pdb_field_index) {
        case e_pdb_occ:
          atom_pdb_field_value = (pdb->atom (ipdb))->occupancy();
          break;
        case e_pdb_beta:
          atom_pdb_field_value = (pdb->atom (ipdb))->temperaturefactor();
          break;
        case e_pdb_x:
          atom_pdb_field_value = (pdb->atom (ipdb))->xcoor();
          break;
        case e_pdb_y:
          atom_pdb_field_value = (pdb->atom (ipdb))->ycoor();
          break;
        case e_pdb_z:
          atom_pdb_field_value = (pdb->atom (ipdb))->zcoor();
          break;
        default:
          break;
        }

        if ( (pdb_field_value) &&
             (atom_pdb_field_value != pdb_field_value) ) {
          continue;
        } else if (atom_pdb_field_value == 0.0) {
          continue;
        }

      } else {
        // Atom ID mode: use predefined atom IDs from the atom group
        if (ipdb != *current_index) {
          // Skip atoms not in the list
          continue;
        } else {
          current_index++;
        }
      }
      
      if (!pos_allocated) {
        pos.push_back (cvm::atom_pos (0.0, 0.0, 0.0));
      } else if (ipos >= pos.size()) {
        cvm::fatal_error ("Error: the PDB file \""+
                          std::string (pdb_filename)+
                          "\" contains coordinates for "
                          "more atoms than needed.\n");
      }

      pos[ipos] = cvm::atom_pos ((pdb->atom (ipdb))->xcoor(),
                                 (pdb->atom (ipdb))->ycoor(),
                                 (pdb->atom (ipdb))->zcoor());
      ipos++;
      if (!use_pdb_field && current_index == indices.end())
        break;
    }

    if (ipos < pos.size())
      cvm::fatal_error ("Error: the PDB file \""+
                        std::string (pdb_filename)+
                        "\" contains coordinates for only "+
                        cvm::to_str (ipos)+
                        " atoms, but "+cvm::to_str (pos.size())+
                        " are needed.\n");

    if (current_index != indices.end())
      cvm::fatal_error ("Error: not all atoms found in PDB file.\n");

  } else {

    // when the PDB contains exactly the number of atoms of the array,
    // ignore the fields and just read coordinates
    for (size_t ia = 0; ia < pos.size(); ia++) {
      pos[ia] = cvm::atom_pos ((pdb->atom (ia))->xcoor(),
                               (pdb->atom (ia))->ycoor(),
                               (pdb->atom (ia))->zcoor());
    }
  }

  delete pdb;
#endif
}

void colvarproxy_lammps::load_atoms (char const *pdb_filename,
                                   std::vector<cvm::atom> &atoms,
                                   std::string const pdb_field_str,
                                   double const pdb_field_value)
{
  cvm::fatal_error ("Selecting collective variable atoms "
		    "from a PDB file is currently not supported.\n");

#if 0
  if (pdb_field_str.size() == 0)
    cvm::fatal_error ("Error: must define which PDB field to use "
                      "in order to define atoms from a PDB file.\n");

  PDB *pdb = new PDB (pdb_filename);
  size_t const pdb_natoms = pdb->num_atoms();

  e_pdb_field pdb_field_index = pdb_field_str2enum (pdb_field_str);

  for (size_t ipdb = 0; ipdb < pdb_natoms; ipdb++) {

    double atom_pdb_field_value = 0.0;

    switch (pdb_field_index) {
    case e_pdb_occ:
      atom_pdb_field_value = (pdb->atom (ipdb))->occupancy();
      break;
    case e_pdb_beta:
      atom_pdb_field_value = (pdb->atom (ipdb))->temperaturefactor();
      break;
    case e_pdb_x:
      atom_pdb_field_value = (pdb->atom (ipdb))->xcoor();
      break;
    case e_pdb_y:
      atom_pdb_field_value = (pdb->atom (ipdb))->ycoor();
      break;
    case e_pdb_z:
      atom_pdb_field_value = (pdb->atom (ipdb))->zcoor();
      break;
    default:
      break;
    }

    if ( (pdb_field_value) &&
         (atom_pdb_field_value != pdb_field_value) ) {
      continue;
    } else if (atom_pdb_field_value == 0.0) {
      continue;
    }
     
    atoms.push_back (cvm::atom (ipdb+1));
  }

  delete pdb;
#endif
}


void colvarproxy_lammps::backup_file (char const *filename)
{
  if (std::string (filename).rfind (std::string (".colvars.state"))
      != std::string::npos) {
    my_backup_file (filename, ".old");
  } else {
    my_backup_file (filename, ".BAK");
  }
}


int colvarproxy_lammps::init_lammps_atom (const int &aid, cvm::atom *atom)
{
  int idx = (*_idlookup)(_idmap,aid);
  if (idx < 0)
    return -1;

  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    if (colvars_atoms[i] == aid) {
      // this atom id was already recorded
      colvars_atoms_ncopies[i] += 1;
      return i;
    }
  }

  struct commdata *c = _coords + idx;
  struct commdata *f = _forces + idx;
  struct commdata *o = _oforce + idx;

  // allocate a new slot for this atom
  colvars_atoms_ncopies.push_back (1);
  colvars_atoms.push_back (aid);
  positions.push_back(cvm::rvector(c->x,c->y,c->z));
  total_forces.push_back (cvm::rvector(o->x,o->y,o->z));
  applied_forces.push_back (cvm::rvector(f->x,f->y,f->z));

  atom->id = c->tag;
  atom->mass = _lmp->atom->mass[c->type];

  return colvars_atoms.size()-1;
}

// atom member functions, LAMMPS specific implementations

cvm::atom::atom (const int &id)
{

  if (cvm::debug())
    cvm::log ("Adding atom "+cvm::to_str(id)+
              " for collective variables calculation.\n");

  if (id < 0)
    cvm::fatal_error ("Error: invalid atom ID specified, "+
                      cvm::to_str (id)+"\n");

  int idx = ((colvarproxy_lammps *) cvm::proxy)->init_lammps_atom(id,this);
  if (idx < 0)
    cvm::fatal_error ("Error: atom ID not in fix colvar group, "+
		      cvm::to_str (id)+"\n");

  this->index = idx;
  if (cvm::debug())
    cvm::log ("The index of this atom in the colvarproxy_lammps arrays is "+
              cvm::to_str (this->index)+".\n");

  this->reset_data();
}


/// For AMBER topologies, the segment id is automatically set to
/// "MAIN" (the segment id assigned by NAMD's AMBER topology parser),
/// and is therefore optional when an AMBER topology is used
cvm::atom::atom (cvm::residue_id const &residue,
                 std::string const     &atom_name,
                 std::string const     &segment_id)
{
  cvm::fatal_error ("Creating collective variable atoms "
		    "from a PDB file is currently not supported.\n");
}


// copy constructor
cvm::atom::atom (cvm::atom const &a)
  : index (a.index), id (a.id), mass (a.mass)
{
  // init_lammps_atom() has already been called by a's constructor, no
  // need to call it again

  // need to increment the counter anyway
  colvarproxy_lammps *cp = (colvarproxy_lammps *) cvm::proxy;
  cp->colvars_atoms_ncopies[this->index] += 1;
}


cvm::atom::~atom() 
{
  colvarproxy_lammps *cp = (colvarproxy_lammps *) cvm::proxy;
  if (cp->colvars_atoms_ncopies[this->index] > 0)
    cp->colvars_atoms_ncopies[this->index] -= 1;
}


void cvm::atom::read_position()
{
  colvarproxy_lammps const * const cp = (colvarproxy_lammps *) cvm::proxy;
  this->pos = cp->positions[this->index];
}


void cvm::atom::read_velocity()
{
  cvm::fatal_error("Error: read_velocity is not yet implemented.\n");
}


void cvm::atom::read_system_force()
{
  colvarproxy_lammps const * const cp = (colvarproxy_lammps *) cvm::proxy;
  this->system_force = cp->total_forces[this->index]
    - cp->applied_forces[this->index];
}


void cvm::atom::apply_force (cvm::rvector const &new_force)
{
  colvarproxy_lammps *cp = (colvarproxy_lammps *) cvm::proxy;
  cp->applied_forces[this->index] = cvm::rvector(new_force.x,
						 new_force.y,
						 new_force.z);
}


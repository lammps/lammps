/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_template.h"
#include "atom.h"
#include "molecule.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecTemplate::AtomVecTemplate(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 2;
  mass_type = 1;

  atom->molecule_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) "molecule molindex molatom";
  fields_copy = (char *) "molecule molindex molatom";
  fields_comm = NULL;
  fields_comm_vel = NULL;
  fields_reverse = NULL;
  fields_border = (char *) "molecule molindex molatom";
  fields_border_vel = (char *) "molecule molindex molatom";
  fields_exchange = (char *) "molecule molindex molatom";
  fields_restart = (char *) "molecule molindex molatom";
  fields_create = (char *) "molecule molindex molatom";
  fields_data_atom = (char *) "id molecule molindex molatom type x";
  fields_data_vel = NULL;

  setup_fields();
}

/* ----------------------------------------------------------------------
   process additional arg = molecule template ID
------------------------------------------------------------------------- */

void AtomVecTemplate::process_args(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal atom_style template command");

  int imol = atom->find_molecule(arg[0]);
  if (imol == -1) error->all(FLERR,"Molecule template ID for "
                             "atom_style template does not exist");

  onemols = &atom->molecules[imol];
  nset = atom->molecules[imol]->nset;

  // error check on molecule template fields

  for (int i = 0; i < nset; i++)
    if (onemols[i]->typeflag == 0)
      error->all(FLERR,"Atom style template molecule must have atom types");

  // set bonds_allow,angles_allow,etc based on the molecules in template set
  // similar to how atom_style bond,angle,full set it

  for (int i = 0; i < nset; i++) {
    if (onemols[i]->bondflag) bonds_allow = 1;
    if (onemols[i]->angleflag) angles_allow = 1;
    if (onemols[i]->dihedralflag) dihedrals_allow = 1;
    if (onemols[i]->improperflag) impropers_allow = 1;
  }

  // set nbondtypes,nangletypes,etc based on the molecules in template set
  // do this here b/c data file will typically not contain these settings

  for (int i = 0; i < nset; i++) {
    atom->nbondtypes = MAX(atom->nbondtypes,onemols[i]->nbondtypes);
    atom->nangletypes = MAX(atom->nangletypes,onemols[i]->nangletypes);
    atom->ndihedraltypes = MAX(atom->ndihedraltypes,onemols[i]->ndihedraltypes);
    atom->nimpropertypes = MAX(atom->nimpropertypes,onemols[i]->nimpropertypes);
  }
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecTemplate::create_atom_post(int ilocal)
{
  atom->molindex[ilocal] = -1;
  atom->molatom[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecTemplate::data_atom_post(int ilocal)
{
  int molindex = atom->molindex[ilocal];
  int molatom = atom->molatom[ilocal];

  if (molindex < 0 || molindex >= nset)
    error->one(FLERR,"Invalid template index in Atoms section of data file");
  if (molatom < 0 || molatom >= onemols[molindex]->natoms)
    error->one(FLERR,"Invalid template atom in Atoms section of data file");
}

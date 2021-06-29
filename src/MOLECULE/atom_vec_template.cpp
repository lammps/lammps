// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_template.h"

#include "atom.h"
#include "error.h"
#include "molecule.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecTemplate::AtomVecTemplate(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::TEMPLATE;
  mass_type = PER_TYPE;

  atom->molecule_flag = 1;
  atom->molindex_flag = 1;
  atom->molatom_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) "molecule molindex molatom";
  fields_copy = (char *) "molecule molindex molatom";
  fields_comm = (char *) "";
  fields_comm_vel = (char *) "";
  fields_reverse = (char *) "";
  fields_border = (char *) "molecule molindex molatom";
  fields_border_vel = (char *) "molecule molindex molatom";
  fields_exchange = (char *) "molecule molindex molatom";
  fields_restart = (char *) "molecule molindex molatom";
  fields_create = (char *) "molecule molindex molatom";
  fields_data_atom = (char *) "id molecule molindex molatom type x";
  fields_data_vel = (char *) "id v";

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
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecTemplate::grow_pointers()
{
  molindex = atom->molindex;
  molatom = atom->molatom;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecTemplate::create_atom_post(int ilocal)
{
  molindex[ilocal] = -1;
  molatom[ilocal] = -1;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecTemplate::pack_data_pre(int ilocal)
{
  molindex[ilocal]++;
  molatom[ilocal]++;
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecTemplate::pack_data_post(int ilocal)
{
  molindex[ilocal]--;
  molatom[ilocal]--;
}


/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecTemplate::data_atom_post(int ilocal)
{
  int molindex_one = --molindex[ilocal];
  int molatom_one = --molatom[ilocal];

  if ((molindex_one < -1) || (molindex_one >= nset))
    error->one(FLERR,"Invalid template index in Atoms section of data file");
  if ((molatom_one < -1) || ((molindex_one >= 0) && (molatom_one >= onemols[molindex_one]->natoms)))
    error->one(FLERR,"Invalid template atom in Atoms section of data file");
}

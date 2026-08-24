/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing author: AUTHOR_NAME_TBD (AFFILIATION_TBD)

   Atom style for inertial spin dynamics ("tspin").  It extends atom style
   spin with a spin velocity v_s, a spin mass s_mass, and a spin force f_spin,
   so that the spin modulus sp[3] becomes a dynamical degree of freedom rather
   than a constant of the motion.  f_spin holds the part of the force on the
   spin vector that is a true energy gradient in energy units, as opposed to
   the precession field fm of the other SPIN styles, which is an angular
   frequency in rad.THz.  The spin itself stays in the (direction, modulus)
   representation of atom style spin, so data files, restart files and the
   set command are unchanged.
------------------------------------------------------------------------- */

#include "atom_vec_tspin.h"

#include "atom.h"
#include "domain.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecTSpin::AtomVecTSpin(LAMMPS *lmp) : AtomVec(lmp), AtomVecSpin(lmp)
{
  atom->tsp_flag = 1;

  // extend the field lists inherited from atom style spin
  // v_s and s_mass are not part of the data file format, so this style reads
  //   and writes the same Atoms and Velocities sections as atom style spin
  // both are needed on owned atoms only, so no border communication
  // f_spin is a force, so it is reverse communicated like fm

  fields_grow.insert(fields_grow.end(), {"v_s", "s_mass", "f_spin"});
  fields_copy.insert(fields_copy.end(), {"v_s", "s_mass"});
  fields_exchange.insert(fields_exchange.end(), {"v_s", "s_mass"});
  fields_restart.insert(fields_restart.end(), {"v_s", "s_mass"});
  fields_create.insert(fields_create.end(), {"v_s", "s_mass"});
  fields_reverse.insert(fields_reverse.end(), {"f_spin"});

  // re-process the field lists; AtomVec::setup_fields() is idempotent

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
------------------------------------------------------------------------- */

void AtomVecTSpin::grow_pointers()
{
  AtomVecSpin::grow_pointers();

  v_s = atom->v_s;
  s_mass = atom->s_mass;
  f_spin = atom->f_spin;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecTSpin::force_clear(int n, size_t nbytes)
{
  AtomVecSpin::force_clear(n, nbytes);

  memset(&f_spin[n][0], 0, 3 * nbytes);
}

/* ----------------------------------------------------------------------
   convert read_data file info from general to restricted triclinic
   parent class operates on x, the Velocities section and the spin direction
   child class operates on the spin velocity v_s
------------------------------------------------------------------------- */

void AtomVecTSpin::read_data_general_to_restricted(int nlocal_previous, int nlocal)
{
  AtomVecSpin::read_data_general_to_restricted(nlocal_previous, nlocal);

  for (int i = nlocal_previous; i < nlocal; i++)
    domain->general_to_restricted_vector(v_s[i]);
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecTSpin::property_atom(const std::string &name)
{
  if (name == "smass") return 0;
  if (name == "vsx") return 1;
  if (name == "vsy") return 2;
  if (name == "vsz") return 3;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecTSpin::pack_property_atom(int index, double *buf, int nvalues, int groupbit)
{
  const int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = s_mass[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else {
    const int j = index - 1;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = v_s[i][j];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  }
}

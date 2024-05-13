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

#include "atom_vec_atomic.h"

#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecAtomic::AtomVecAtomic(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

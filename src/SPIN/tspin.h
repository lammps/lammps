/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com

   Shared state of the inertial spin dynamics (tspin) styles.

   These styles run on the ordinary atom_style spin.  The spin itself is the
   per-atom array sp, holding a direction sp[0..2] and a modulus sp[3] in Bohr
   magnetons.  Magnetic interactions continue to use the per-atom fm array in
   rad.THz.  For use with TSPIN, a force producer must encode the complete
   unconstrained-spin force as fm = |S|/hbar * (-dU/dS).

   Inertial spin dynamics adds two quantities that the fixed-modulus styles do
   not need: a spin velocity and a spin mass.  They are kept in the two custom
   per-atom properties below, created on demand through fix property/atom, so
   that they migrate with the atoms and are written to restart files without
   any change to the atom style.  They can also be output directly, for
   instance with

     compute v all property/atom d2_tspin_vs[1] d2_tspin_vs[2] d2_tspin_vs[3]
------------------------------------------------------------------------- */

#ifndef LMP_TSPIN_H
#define LMP_TSPIN_H

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "modify.h"
#include "utils.h"

#include <string>

namespace LAMMPS_NS {

// names of the two custom per-atom properties and the fix that owns them

static constexpr char TSPIN_STATE_ID[] = "TSPIN_STATE";
static constexpr char TSPIN_VS[] = "tspin_vs";
static constexpr char TSPIN_SMASS[] = "tspin_smass";

// smallest initial spin modulus that counts as a magnetic atom when no spin
// mass has been assigned yet.  Once an atom has a positive spin mass it is a
// dynamical spin and remains active even if its modulus passes through zero.

static constexpr double TSPIN_EPS = 1.0e-8;

/* ----------------------------------------------------------------------
   create the two custom per-atom properties unless they already exist
------------------------------------------------------------------------- */

static inline void tspin_create_state(Modify *modify, Error *error)
{
  Fix *state = modify->get_fix_by_id(TSPIN_STATE_ID);
  if (state) {
    if (!utils::strmatch(state->style, "^property/atom"))
      error->all(FLERR, "Fix ID {} is reserved for the internal state of the tspin styles",
                 TSPIN_STATE_ID);
    return;
  }

  // writedata no keeps the internal state out of write_data output; it is
  // carried across runs by restart files instead

  modify->add_fix(std::string(TSPIN_STATE_ID) + " all property/atom d2_" + TSPIN_VS + " 3 d_" +
                  TSPIN_SMASS + " ghost no writedata no");
}

/* ----------------------------------------------------------------------
   index of the spin velocity array and of the spin mass vector
   the indices are stable, the pointers behind them are not, so look the
   pointer up again through atom->darray / atom->dvector on every use
------------------------------------------------------------------------- */

static inline void tspin_bind_state(Atom *atom, Error *error, const std::string &caller,
                                    int &index_vs, int &index_sm)
{
  int flag, cols;
  index_vs = atom->find_custom(TSPIN_VS, flag, cols);
  const bool vs_ok = (index_vs >= 0) && (flag == 1) && (cols == 3);
  index_sm = atom->find_custom(TSPIN_SMASS, flag, cols);
  const bool sm_ok = (index_sm >= 0) && (flag == 1) && (cols == 0);

  if (!vs_ok || !sm_ok)
    error->all(FLERR,
               "{} needs the internal per-atom state of the tspin styles; define one of the "
               "tspin integrators first and do not remove fix {}",
               caller, TSPIN_STATE_ID);
}

}    // namespace LAMMPS_NS

#endif

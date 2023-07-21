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

#include "reset_atoms_image.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ResetAtomsImage::ResetAtomsImage(LAMMPS *lmp) : Command(lmp) {}

/* ----------------------------------------------------------------------
   called as reset_atoms image command in input script
------------------------------------------------------------------------- */

void ResetAtomsImage::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Reset_atoms image command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR, "Cannot use reset_atoms image unless atoms have IDs");
  if (atom->avec->bonds_allow == 0)
    error->all(FLERR, "Cannot use reset_atoms image used when bonds are not allowed");

  // process args

  if (narg < 1) utils::missing_cmd_args(FLERR, "reset_atoms image", error);
  if (narg > 1) error->all(FLERR, "Unknown reset_atoms image keyword: {}", arg[1]);
  int igroup = group->find(arg[0]);
  if (igroup < 0) error->all(FLERR, "Could not find reset_atoms image group {}", arg[0]);
  int groupbit = group->bitmask[igroup];

  if (comm->me == 0) utils::logmesg(lmp, "Resetting image flags ...\n");

  // record wall time for resetting molecule IDs

  MPI_Barrier(world);
  double time1 = platform::walltime();

  // create computes and variables
  // must come before lmp->init so the computes are properly initialized

  auto frags = modify->add_compute("frags_r_i_f all fragment/atom single yes");
  auto chunk = modify->add_compute("chunk_r_i_f all chunk/atom c_frags_r_i_f compress yes");
  auto flags = modify->add_compute("flags_r_i_f all property/atom ix iy iz");
  input->variable->set("ix_r_i_f atom c_flags_r_i_f[1]");
  input->variable->set("iy_r_i_f atom c_flags_r_i_f[2]");
  input->variable->set("iz_r_i_f atom c_flags_r_i_f[3]");
  auto ifmin = modify->add_compute("ifmin_r_i_f all reduce/chunk chunk_r_i_f min "
                                   "v_ix_r_i_f v_iy_r_i_f v_iz_r_i_f");
  auto ifmax = modify->add_compute("ifmax_r_i_f all reduce/chunk chunk_r_i_f max "
                                   "v_ix_r_i_f v_iy_r_i_f v_iz_r_i_f");
  auto cdist = modify->add_compute("cdist_r_i_f all chunk/spread/atom chunk_r_i_f "
                                   "c_ifmax_r_i_f[*] c_ifmin_r_i_f[*]");

  // initialize system since comm->borders() will be invoked

  lmp->init();

  // setup domain, communication
  // exchange will clear map, borders will reset
  // this is the map needed to lookup current global IDs for bond topology

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);

  // trigger computes

  frags->compute_peratom();
  chunk->compute_peratom();
  flags->compute_peratom();
  ifmin->compute_array();
  ifmax->compute_array();
  cdist->compute_peratom();

  // reset image flags for atoms in group

  const int *const mask = atom->mask;
  const int nlocal = atom->nlocal;
  imageint *image = atom->image;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int ix = (image[i] & IMGMASK) - IMGMAX;
      int iy = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int iz = (image[i] >> IMG2BITS) - IMGMAX;
      ix -= (int) floor(0.5 * (cdist->array_atom[i][0] + cdist->array_atom[i][3]));
      iy -= (int) floor(0.5 * (cdist->array_atom[i][1] + cdist->array_atom[i][4]));
      iz -= (int) floor(0.5 * (cdist->array_atom[i][2] + cdist->array_atom[i][5]));
      image[i] = ((imageint) (ix + IMGMAX) & IMGMASK) |
          (((imageint) (iy + IMGMAX) & IMGMASK) << IMGBITS) |
          (((imageint) (iz + IMGMAX) & IMGMASK) << IMG2BITS);
    }
  }

  // cleanup

  modify->delete_compute("cdist_r_i_f");
  modify->delete_compute("ifmax_r_i_f");
  modify->delete_compute("ifmin_r_i_f");
  modify->delete_compute("flags_r_i_f");
  modify->delete_compute("chunk_r_i_f");
  modify->delete_compute("frags_r_i_f");
  input->variable->set("ix_r_i_f delete");
  input->variable->set("iy_r_i_f delete");
  input->variable->set("iz_r_i_f delete");

  // total time

  MPI_Barrier(world);
  if (comm->me == 0)
    utils::logmesg(lmp, "  reset_atoms image CPU = {:.3f} seconds\n", platform::walltime() - time1);
}

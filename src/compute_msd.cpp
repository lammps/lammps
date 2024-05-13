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

#include "compute_msd.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "group.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMSD::ComputeMSD(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg), id_fix(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal compute msd command");

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  create_attribute = 1;
  dynamic_group_allow = 0;

  // optional args

  comflag = 0;
  avflag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "com") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "compute msd com", error);
      comflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "average") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "compute msd average", error);
      avflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown compute msd keyword: {}", arg[iarg]);
  }

  if (group->dynamic[igroup])
    error->all(FLERR, "Compute {} is not compatible with dynamic groups", style);

  // create a new fix STORE style for reference positions
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  id_fix = utils::strdup(id + std::string("_COMPUTE_STORE"));
  fix = dynamic_cast<FixStoreAtom *>(
      modify->add_fix(fmt::format("{} {} STORE/ATOM 3 0 0 1", id_fix, group->names[igroup])));

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset)
    fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;

    double **x = atom->x;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        domain->unmap(x[i], image[i], xoriginal[i]);
      else
        xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;

    // adjust for COM if requested

    if (comflag) {
      double cm[3];
      masstotal = group->mass(igroup);
      group->xcm(igroup, masstotal, cm);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          xoriginal[i][0] -= cm[0];
          xoriginal[i][1] -= cm[1];
          xoriginal[i][2] -= cm[2];
        }
    }

    // initialize counter for average positions if requested

    naverage = 0;
  }

  // displacement vector

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeMSD::~ComputeMSD()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete[] id_fix;
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeMSD::init()
{
  // set fix which stores reference atom coords

  fix = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(id_fix));
  if (!fix) error->all(FLERR, "Could not find compute msd fix with ID {}", id_fix);

  // nmsd = # of atoms in group

  nmsd = group->count(igroup);
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeMSD::compute_vector()
{
  invoked_vector = update->ntimestep;

  // cm = current center of mass

  double cm[3];
  if (comflag)
    group->xcm(igroup, masstotal, cm);
  else
    cm[0] = cm[1] = cm[2] = 0.0;

  // dx,dy,dz = displacement of atom from reference position
  // reference unwrapped position is stored by fix
  // relative to center of mass if comflag is set
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->astore;

  double **x = atom->x;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double dx, dy, dz;
  int xbox, ybox, zbox;

  double msd[4];
  msd[0] = msd[1] = msd[2] = msd[3] = 0.0;

  double xtmp, ytmp, ztmp;

  // update number of averages if requested

  double navfac;
  if (avflag) {
    naverage++;
    navfac = 1.0 / (naverage + 1);
  }

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        xtmp = x[i][0] + xbox * xprd - cm[0];
        ytmp = x[i][1] + ybox * yprd - cm[1];
        ztmp = x[i][2] + zbox * zprd - cm[2];

        // use running average position for reference if requested

        if (avflag) {
          xoriginal[i][0] = (xoriginal[i][0] * naverage + xtmp) * navfac;
          xoriginal[i][1] = (xoriginal[i][1] * naverage + ytmp) * navfac;
          xoriginal[i][2] = (xoriginal[i][2] * naverage + ztmp) * navfac;
        }

        dx = xtmp - xoriginal[i][0];
        dy = ytmp - xoriginal[i][1];
        dz = ztmp - xoriginal[i][2];
        msd[0] += dx * dx;
        msd[1] += dy * dy;
        msd[2] += dz * dz;
        msd[3] += dx * dx + dy * dy + dz * dz;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        xtmp = x[i][0] + h[0] * xbox + h[5] * ybox + h[4] * zbox - cm[0];
        ytmp = x[i][1] + h[1] * ybox + h[3] * zbox - cm[1];
        ztmp = x[i][2] + h[2] * zbox - cm[2];

        // use running average position for reference if requested

        if (avflag) {
          xoriginal[i][0] = (xoriginal[i][0] * naverage + xtmp) * navfac;
          xoriginal[i][1] = (xoriginal[i][0] * naverage + xtmp) * navfac;
          xoriginal[i][2] = (xoriginal[i][0] * naverage + xtmp) * navfac;
        }

        dx = xtmp - xoriginal[i][0];
        dy = ytmp - xoriginal[i][1];
        dz = ztmp - xoriginal[i][2];
        msd[0] += dx * dx;
        msd[1] += dy * dy;
        msd[2] += dz * dz;
        msd[3] += dx * dx + dy * dy + dz * dz;
      }
  }

  MPI_Allreduce(msd, vector, 4, MPI_DOUBLE, MPI_SUM, world);
  if (nmsd) {
    vector[0] /= nmsd;
    vector[1] /= nmsd;
    vector[2] /= nmsd;
    vector[3] /= nmsd;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputeMSD::set_arrays(int i)
{
  double **xoriginal = fix->astore;
  double **x = atom->x;
  xoriginal[i][0] = x[i][0];
  xoriginal[i][1] = x[i][1];
  xoriginal[i][2] = x[i][2];
}

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

#include "change_box.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "irregular.h"
#include "lattice.h"
#include "modify.h"
#include "output.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum{XYZ=0,TILT,BOUNDARY,ORTHO,TRICLINIC,SET,REMAP};
enum{FINAL=0,DELTA,SCALE};
enum{X=0,Y,Z,YZ,XZ,XY};

/* ---------------------------------------------------------------------- */

ChangeBox::ChangeBox(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void ChangeBox::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0)
    error->all(FLERR,"Change_box command before simulation box is defined");
  if (narg < 2) utils::missing_cmd_args(FLERR, "change_box", error);

  if (comm->me == 0) utils::logmesg(lmp,"Changing box ...\n");

  // group

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find change_box group ID {}", arg[0]);
  int groupbit = group->bitmask[igroup];

  // parse operation arguments
  // allocate ops to max possible length
  // volume option does not increment nops

  int dimension = domain->dimension;

  ops = new Operation[narg-1];
  memset(ops,0,(narg-1)*sizeof(Operation));
  nops = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0 || strcmp(arg[iarg],"y") == 0 ||
        strcmp(arg[iarg],"z") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal change_box {} command: missing argument(s)", arg[iarg]);
      ops[nops].style = XYZ;
      if (strcmp(arg[iarg],"x") == 0) ops[nops].dim = X;
      else if (strcmp(arg[iarg],"y") == 0) ops[nops].dim = Y;
      else if (strcmp(arg[iarg],"z") == 0) ops[nops].dim = Z;

      if (dimension == 2 && ops[nops].dim == Z)
        error->all(FLERR,"Cannot change_box in z dimension for 2d simulation");

      if (strcmp(arg[iarg+1],"final") == 0) {
        if (iarg+4 > narg) error->all(FLERR, "Illegal change_box {} final command: missing argument(s)", arg[iarg]);
        ops[nops].flavor = FINAL;
        ops[nops].flo = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        ops[nops].fhi = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        ops[nops].vdim1 = ops[nops].vdim2 = -1;
        nops++;
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"delta") == 0) {
        if (iarg+4 > narg) error->all(FLERR, "Illegal change_box {} delta command: missing argument(s)", arg[iarg]);
        ops[nops].flavor = DELTA;
        ops[nops].dlo = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        ops[nops].dhi = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        ops[nops].vdim1 = ops[nops].vdim2 = -1;
        nops++;
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"scale") == 0) {
        if (iarg+3 > narg) error->all(FLERR, "Illegal change_box {} scale command: missing argument(s)", arg[iarg]);
        ops[nops].flavor = SCALE;
        ops[nops].scale = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        ops[nops].vdim1 = ops[nops].vdim2 = -1;
        nops++;
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"volume") == 0) {
        if (nops == 0 || ops[nops-1].style != XYZ ||
            ops[nops].dim == ops[nops-1].dim)
          error->all(FLERR,"Change_box volume used incorrectly");
        if (ops[nops-1].vdim2 >= 0)
          error->all(FLERR,"Change_box volume used incorrectly");
        else if (ops[nops-1].vdim1 >= 0) ops[nops-1].vdim2 = ops[nops].dim;
        else ops[nops-1].vdim1 = ops[nops].dim;
        iarg += 2;

      } else error->all(FLERR, "Unknown change_box {} argument: {}", arg[iarg], arg[iarg+1]);

    } else if (strcmp(arg[iarg],"xy") == 0 || strcmp(arg[iarg],"xz") == 0 ||
        strcmp(arg[iarg],"yz") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal change_box {} command: missing argument(s)", arg[iarg]);
      ops[nops].style = TILT;
      if (strcmp(arg[iarg],"xy") == 0) ops[nops].dim = XY;
      else if (strcmp(arg[iarg],"xz") == 0) ops[nops].dim = XZ;
      else if (strcmp(arg[iarg],"yz") == 0) ops[nops].dim = YZ;

      if (dimension == 2 && (ops[nops].dim == XZ || ops[nops].dim == YZ))
        error->all(FLERR,"Cannot change_box in xz or yz for 2d simulation");

      if (strcmp(arg[iarg+1],"final") == 0) {
        if (iarg+3 > narg) error->all(FLERR, "Illegal change_box {} final command: missing argument(s)", arg[iarg]);
        ops[nops].flavor = FINAL;
        ops[nops].ftilt = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        nops++;
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"delta") == 0) {
        if (iarg+3 > narg) error->all(FLERR, "Illegal change_box {} delta command: missing argument(s)", arg[iarg]);
        ops[nops].flavor = DELTA;
        ops[nops].dtilt = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        nops++;
        iarg += 3;
      } else error->all(FLERR, "Unknown change_box {} argument: {}", arg[iarg], arg[iarg+1]);

    } else if (strcmp(arg[iarg],"boundary") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "change_box boundary", error);
      ops[nops].style = BOUNDARY;
      ops[nops].boundindex = iarg+1;
      nops++;
      iarg += 4;

    } else if (strcmp(arg[iarg],"ortho") == 0) {
      ops[nops].style = ORTHO;
      nops++;
      iarg += 1;

    } else if (strcmp(arg[iarg],"triclinic") == 0) {
      ops[nops].style = TRICLINIC;
      nops++;
      iarg += 1;

    } else if (strcmp(arg[iarg],"set") == 0) {
      ops[nops].style = SET;
      nops++;
      iarg += 1;

    } else if (strcmp(arg[iarg],"remap") == 0) {
      ops[nops].style = REMAP;
      nops++;
      iarg += 1;

    } else break;
  }

  if (nops == 0) error->all(FLERR, "Unknown change_box keyword: {}", arg[iarg]);

  // move_atoms = 1 if need to move atoms to new procs after box changes
  // anything other than ORTHO or TRICLINIC may cause atom movement

  int move_atoms = 0;
  for (int m = 0; m < nops; m++) {
    if (ops[m].style != ORTHO && ops[m].style != TRICLINIC) move_atoms = 1;
  }

  // error if moving atoms and there is stored per-atom restart state
  // disallowed b/c restart per-atom fix info will not move with atoms

  if (move_atoms && modify->nfix_restart_peratom)
    error->all(FLERR,"Change_box parameter not allowed after "
               "reading restart file with per-atom info");

  // read options from end of input line

  options(narg-iarg,&arg[iarg]);

  // compute scale factors if FINAL,DELTA used since they have distance units

  int flag = 0;
  for (i = 0; i < nops; i++)
    if (ops[i].style == FINAL || ops[i].style == DELTA) flag = 1;

  if (flag && scaleflag) {
    scale[0] = domain->lattice->xlattice;
    scale[1] = domain->lattice->ylattice;
    scale[2] = domain->lattice->zlattice;
  }
  else scale[0] = scale[1] = scale[2] = 1.0;

  // perform sequence of operations
  // first insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs
  // save current box state so can remap atoms from it, if requested

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);
  save_box_state();

  for (int m = 0; m < nops; m++) {
    if (ops[m].style == XYZ) {
      double volume;
      if (domain->dimension == 2) volume = domain->xprd * domain->yprd;
      else volume = domain->xprd * domain->yprd * domain->zprd;

      if (ops[m].flavor == FINAL) {
        domain->boxlo[ops[m].dim] = scale[ops[m].dim]*ops[m].flo;
        domain->boxhi[ops[m].dim] = scale[ops[m].dim]*ops[m].fhi;
        if (ops[m].vdim1 >= 0)
          volume_preserve(ops[m].vdim1,ops[m].vdim2,volume);
        domain->set_initial_box();
        domain->set_global_box();
        domain->set_local_box();
        domain->print_box("  ");

      } else if (ops[m].flavor == DELTA) {
        domain->boxlo[ops[m].dim] += scale[ops[m].dim]*ops[m].dlo;
        domain->boxhi[ops[m].dim] += scale[ops[m].dim]*ops[m].dhi;
        if (ops[m].vdim1 >= 0)
          volume_preserve(ops[m].vdim1,ops[m].vdim2,volume);
        domain->set_initial_box();
        domain->set_global_box();
        domain->set_local_box();
        domain->print_box("  ");

      } else if (ops[m].flavor == SCALE) {
        double mid = 0.5 *
          (domain->boxlo[ops[m].dim] + domain->boxhi[ops[m].dim]);
        double delta = domain->boxlo[ops[m].dim] - mid;
        domain->boxlo[ops[m].dim] = mid + ops[m].scale*delta;
        delta = domain->boxhi[ops[m].dim] - mid;
        domain->boxhi[ops[m].dim] = mid + ops[m].scale*delta;
        if (ops[m].vdim1 >= 0)
          volume_preserve(ops[m].vdim1,ops[m].vdim2,volume);
        domain->set_initial_box();
        domain->set_global_box();
        domain->set_local_box();
        domain->print_box("  ");
      }

    } else if (ops[m].style == TILT) {
      if (domain->triclinic == 0)
        error->all(FLERR,"Cannot change box tilt factors for orthogonal box");

      if (ops[m].flavor == FINAL) {
        if (ops[m].dim == XY) domain->xy = scale[X]*ops[m].ftilt;
        else if (ops[m].dim == XZ) domain->xz = scale[X]*ops[m].ftilt;
        else if (ops[m].dim == YZ) domain->yz = scale[Y]*ops[m].ftilt;
        domain->set_initial_box();
        domain->set_global_box();
        domain->set_local_box();
        domain->print_box("  ");

      } else if (ops[m].flavor == DELTA) {
        if (ops[m].dim == XY) domain->xy += scale[X]*ops[m].dtilt;
        else if (ops[m].dim == XZ) domain->xz += scale[X]*ops[m].dtilt;
        else if (ops[m].dim == YZ) domain->yz += scale[Y]*ops[m].dtilt;
        domain->set_initial_box();
        domain->set_global_box();
        domain->set_local_box();
        domain->print_box("  ");
      }

    } else if (ops[m].style == BOUNDARY) {
      domain->set_boundary(3,&arg[ops[m].boundindex],1);
      if (domain->dimension == 2 && domain->zperiodic == 0)
        error->all(FLERR,
                   "Cannot change box z boundary to "
                   "non-periodic for a 2d simulation");
      domain->set_initial_box();
      domain->set_global_box();
      domain->set_local_box();

    } else if (ops[m].style == ORTHO) {
      if (domain->xy != 0.0 || domain->yz != 0.0 || domain->xz != 0.0)
        error->all(FLERR,"Cannot change box to orthogonal when tilt is non-zero");
      if (output->ndump)
        error->all(FLERR,"Cannot change box ortho/triclinic with dumps defined");
      for (const auto &fix : modify->get_fix_list())
        if (fix->no_change_box)
          error->all(FLERR,"Cannot change box ortho/triclinic with certain fixes defined");
      domain->triclinic = 0;
      domain->set_initial_box();
      domain->set_global_box();
      domain->set_local_box();
      domain->print_box("  ");

    } else if (ops[m].style == TRICLINIC) {
      if (output->ndump)
        error->all(FLERR,"Cannot change box ortho/triclinic with dumps defined");
      for (const auto &fix : modify->get_fix_list())
        if (fix->no_change_box)
          error->all(FLERR,"Cannot change box ortho/triclinic with certain fixes defined");
      domain->triclinic = 1;
      domain->set_lamda_box();
      domain->set_initial_box();
      domain->set_global_box();
      domain->set_local_box();
      domain->print_box("  ");

    } else if (ops[m].style == SET) {
      save_box_state();

    } else if (ops[m].style == REMAP) {

      if (modify->check_rigid_group_overlap(groupbit))
        error->warning(FLERR,"Attempting to remap atoms in rigid bodies");

      // convert atoms to lamda coords, using last box state
      // convert atoms back to box coords, using current box state
      // save current box state

      double **x = atom->x;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          domain->x2lamda(x[i],x[i],boxlo,h_inv);

      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          domain->lamda2x(x[i],x[i]);

      save_box_state();
    }
  }

  // clean up

  delete [] ops;

  // apply shrink-wrap boundary conditions

  if (domain->nonperiodic == 2) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->reset_box();
    if (domain->triclinic) domain->lamda2x(atom->nlocal);
  }

  // done if don't need to move atoms

  if (!move_atoms) return;

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc()
  //   in case box moved a long distance relative to atoms
  // use irregular() in case box moved a long distance relative to atoms

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  auto irregular = new Irregular(lmp);
  irregular->migrate_atoms(1);
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0)
    error->warning(FLERR,"Lost atoms via change_box: original {} "
                   "current {}", atom->natoms,natoms);
}

/* ----------------------------------------------------------------------
   parse optional parameters
------------------------------------------------------------------------- */

void ChangeBox::options(int narg, char **arg)
{
  if (narg < 0) utils::missing_cmd_args(FLERR, "change_box", error);

  scaleflag = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "change_box units", error);
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR, "Invalid change_box units argument: {}", arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Unknown change_box keyword: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   save current box state for converting atoms to lamda coords
------------------------------------------------------------------------- */

void ChangeBox::save_box_state()
{
  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];

  for (int i = 0; i < 6; i++)
    h_inv[i] = domain->h_inv[i];
}

/* ----------------------------------------------------------------------
   oldvol = box volume before dim3 changed
   newvol = box volume after dim3 changed
   reset box lengths of dim1/2 to preserve old volume
------------------------------------------------------------------------- */

void ChangeBox::volume_preserve(int dim1, int dim2, double oldvol)
{
  // invoke set_initial_box()
  // in case change by caller to dim3 was invalid or on shrink-wrapped dim

  domain->set_initial_box();

  // calculate newvol using boxlo/hi since xyz prd are not yet reset

  double newvol;
  if (domain->dimension == 2) {
    newvol = domain->boxhi[0] - domain->boxlo[0];
    newvol *= domain->boxhi[1] - domain->boxlo[1];
  } else {
    newvol = domain->boxhi[0] - domain->boxlo[0];
    newvol *= domain->boxhi[1] - domain->boxlo[1];
    newvol *= domain->boxhi[2] - domain->boxlo[2];
  }

  double scale = oldvol/newvol;
  double mid,delta;

  // change dim1 only

  if (dim2 < 0) {
    mid = 0.5 * (domain->boxlo[dim1] + domain->boxhi[dim1]);
    delta = domain->boxlo[dim1] - mid;
    domain->boxlo[dim1] = mid + scale*delta;
    delta = domain->boxhi[dim1] - mid;
    domain->boxhi[dim1] = mid + scale*delta;

  // change dim1 and dim2, keeping their relative aspect ratio constant
  // both are scaled by sqrt(scale)

  } else {
    mid = 0.5 * (domain->boxlo[dim1] + domain->boxhi[dim1]);
    delta = domain->boxlo[dim1] - mid;
    domain->boxlo[dim1] = mid + sqrt(scale)*delta;
    delta = domain->boxhi[dim1] - mid;
    domain->boxhi[dim1] = mid + sqrt(scale)*delta;
    mid = 0.5 * (domain->boxlo[dim2] + domain->boxhi[dim2]);
    delta = domain->boxlo[dim2] - mid;
    domain->boxlo[dim2] = mid + sqrt(scale)*delta;
    delta = domain->boxhi[dim2] - mid;
    domain->boxhi[dim2] = mid + sqrt(scale)*delta;
  }
}

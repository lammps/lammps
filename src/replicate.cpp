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

#include "stdlib.h"
#include "string.h"
#include "replicate.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define LB_FACTOR 1.1
#define EPSILON   1.0e-6

/* ---------------------------------------------------------------------- */

Replicate::Replicate(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Replicate::command(int narg, char **arg)
{
  int i,j,m,n;

  if (domain->box_exist == 0)
    error->all(FLERR,"Replicate command before simulation box is defined");
  if (narg != 3) error->all(FLERR,"Illegal replicate command");

  int me = comm->me;
  int nprocs = comm->nprocs;

  if (me == 0 && screen) fprintf(screen,"Replicating atoms ...\n");

  // nrep = total # of replications

  int nx = force->inumeric(FLERR,arg[0]);
  int ny = force->inumeric(FLERR,arg[1]);
  int nz = force->inumeric(FLERR,arg[2]);
  int nrep = nx*ny*nz;

  // error and warning checks

  if (nx <= 0 || ny <= 0 || nz <= 0)
    error->all(FLERR,"Illegal replicate command");
  if (domain->dimension == 2 && nz != 1)
    error->all(FLERR,"Cannot replicate 2d simulation in z dimension");
  if ((nx > 1 && domain->xperiodic == 0) ||
      (ny > 1 && domain->yperiodic == 0) ||
      (nz > 1 && domain->zperiodic == 0)) {
    if (comm->me == 0)
      error->warning(FLERR,"Replicating in a non-periodic dimension");
  }

  if (atom->nextra_grow || atom->nextra_restart || atom->nextra_store)
    error->all(FLERR,"Cannot replicate with fixes that store atom quantities");

  // maxtag = largest atom tag across all existing atoms

  tagint maxtag = 0;
  if (atom->tag_enable) {
    for (i = 0; i < atom->nlocal; i++) maxtag = MAX(atom->tag[i],maxtag);
    tagint maxtag_all;
    MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
    maxtag = maxtag_all;
  }

  // maxmol = largest molecule tag across all existing atoms

  tagint maxmol = 0;
  if (atom->molecule_flag) {
    for (i = 0; i < atom->nlocal; i++) maxmol = MAX(atom->molecule[i],maxmol);
    tagint maxmol_all;
    MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
    maxmol = maxmol_all;
  }

  // unmap existing atoms via image flags

  for (i = 0; i < atom->nlocal; i++)
    domain->unmap(atom->x[i],atom->image[i]);

  // communication buffer for all my atom's info
  // max_size = largest buffer needed by any proc
  // must do before new Atom class created,
  //   since size_restart() uses atom->nlocal

  int max_size;
  int send_size = atom->avec->size_restart();
  MPI_Allreduce(&send_size,&max_size,1,MPI_INT,MPI_MAX,world);

  double *buf;
  memory->create(buf,max_size,"replicate:buf");

  // old = original atom class
  // atom = new replicated atom class
  // if old atom style was hybrid, pass sub-style names to create_avec

  Atom *old = atom;
  atom = new Atom(lmp);
  atom->settings(old);

  int nstyles = 0;
  char **keywords = NULL;
  if (strcmp(old->atom_style,"hybrid") == 0) {
    AtomVecHybrid *avec_hybrid = (AtomVecHybrid *) old->avec;
    nstyles = avec_hybrid->nstyles;
    keywords = avec_hybrid->keywords;
  }
  atom->create_avec(old->atom_style,nstyles,keywords);

  // check that new system will not be too large
  // new tags cannot exceed MAXTAGINT
  // new system sizes cannot exceed MAXBIGINT

  if (atom->tag_enable) {
    bigint maxnewtag = maxtag + (nrep-1)*old->natoms;
    if (maxnewtag < 0 || maxnewtag >= MAXTAGINT)
      error->all(FLERR,"Replicated system atom IDs are too big");
  }

  if (nrep*old->natoms < 0 || nrep*old->natoms >= MAXBIGINT ||
      nrep*old->nbonds < 0 || nrep*old->nbonds >= MAXBIGINT ||
      nrep*old->nangles < 0 || nrep*old->nangles >= MAXBIGINT ||
      nrep*old->ndihedrals < 0 || nrep*old->ndihedrals >= MAXBIGINT ||
      nrep*old->nimpropers < 0 || nrep*old->nimpropers >= MAXBIGINT)
    error->all(FLERR,"Replicated system is too big");

  // assign atom and topology counts in new class from old one

  atom->natoms = old->natoms * nrep;
  atom->nbonds = old->nbonds * nrep;
  atom->nangles = old->nangles * nrep;
  atom->ndihedrals = old->ndihedrals * nrep;
  atom->nimpropers = old->nimpropers * nrep;

  atom->ntypes = old->ntypes;
  atom->nbondtypes = old->nbondtypes;
  atom->nangletypes = old->nangletypes;
  atom->ndihedraltypes = old->ndihedraltypes;
  atom->nimpropertypes = old->nimpropertypes;

  atom->bond_per_atom = old->bond_per_atom;
  atom->angle_per_atom = old->angle_per_atom;
  atom->dihedral_per_atom = old->dihedral_per_atom;
  atom->improper_per_atom = old->improper_per_atom;

  // store old simulation box

  int triclinic = domain->triclinic;
  double old_xprd = domain->xprd;
  double old_yprd = domain->yprd;
  double old_zprd = domain->zprd;
  double old_xy = domain->xy;
  double old_xz = domain->xz;
  double old_yz = domain->yz;

  // setup new simulation box

  domain->boxhi[0] = domain->boxlo[0] + nx*old_xprd;
  domain->boxhi[1] = domain->boxlo[1] + ny*old_yprd;
  domain->boxhi[2] = domain->boxlo[2] + nz*old_zprd;
  if (triclinic) {
    domain->xy *= ny;
    domain->xz *= nz;
    domain->yz *= nz;
  }

  // new problem setup using new box boundaries

  if (nprocs == 1) n = static_cast<int> (atom->natoms);
  else n = static_cast<int> (LB_FACTOR * atom->natoms / nprocs);

  atom->allocate_type_arrays();
  atom->avec->grow(n);
  n = atom->nmax;

  domain->print_box("  ");
  domain->set_initial_box();
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();

  // copy type arrays to new atom class

  if (atom->mass) {
    for (int itype = 1; itype <= atom->ntypes; itype++) {
      atom->mass_setflag[itype] = old->mass_setflag[itype];
      if (atom->mass_setflag[itype]) atom->mass[itype] = old->mass[itype];
    }
  }

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all replicated atoms will be owned even with round-off

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3],subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (domain->xperiodic) {
    if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
    if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += epsilon[0];
  }
  if (domain->yperiodic) {
    if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
    if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += epsilon[1];
  }
  if (domain->zperiodic) {
    if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
    if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += epsilon[2];
  }

  // loop over all procs
  // if this iteration of loop is me:
  //   pack my unmapped atom data into buf
  //   bcast it to all other procs
  // performs 3d replicate loop with while loop over atoms in buf
  //   x = new replicated position, remapped into simulation box
  //   unpack atom into new atom class from buf if I own it
  //   adjust tag, mol #, coord, topology info as needed

  AtomVec *old_avec = old->avec;
  AtomVec *avec = atom->avec;

  int ix,iy,iz;
  tagint atom_offset,mol_offset;
  imageint image;
  double x[3],lamda[3];
  double *coord;
  int tag_enable = atom->tag_enable;

  for (int iproc = 0; iproc < nprocs; iproc++) {
    if (me == iproc) {
      n = 0;
      for (i = 0; i < old->nlocal; i++) n += old_avec->pack_restart(i,&buf[n]);
    }
    MPI_Bcast(&n,1,MPI_INT,iproc,world);
    MPI_Bcast(buf,n,MPI_DOUBLE,iproc,world);

    for (ix = 0; ix < nx; ix++) {
      for (iy = 0; iy < ny; iy++) {
        for (iz = 0; iz < nz; iz++) {

          // while loop over one proc's atom list

          m = 0;
          while (m < n) {
            image = ((imageint) IMGMAX << IMG2BITS) |
              ((imageint) IMGMAX << IMGBITS) | IMGMAX;
            if (triclinic == 0) {
              x[0] = buf[m+1] + ix*old_xprd;
              x[1] = buf[m+2] + iy*old_yprd;
              x[2] = buf[m+3] + iz*old_zprd;
            } else {
              x[0] = buf[m+1] + ix*old_xprd + iy*old_xy + iz*old_xz;
              x[1] = buf[m+2] + iy*old_yprd + iz*old_yz;
              x[2] = buf[m+3] + iz*old_zprd;
            }
            domain->remap(x,image);
            if (triclinic) {
              domain->x2lamda(x,lamda);
              coord = lamda;
            } else coord = x;

            if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
                coord[1] >= sublo[1] && coord[1] < subhi[1] &&
                coord[2] >= sublo[2] && coord[2] < subhi[2]) {

              m += avec->unpack_restart(&buf[m]);

              i = atom->nlocal - 1;
              if (tag_enable)
                atom_offset = iz*ny*nx*maxtag + iy*nx*maxtag + ix*maxtag;
              else atom_offset = 0;
              mol_offset = iz*ny*nx*maxmol + iy*nx*maxmol + ix*maxmol;

              atom->x[i][0] = x[0];
              atom->x[i][1] = x[1];
              atom->x[i][2] = x[2];

              atom->tag[i] += atom_offset;
              atom->image[i] = image;

              if (atom->molecule_flag) {
                if (atom->molecule[i] > 0)
                  atom->molecule[i] += mol_offset;
                if (atom->avec->bonds_allow)
                  for (j = 0; j < atom->num_bond[i]; j++)
                    atom->bond_atom[i][j] += atom_offset;
                if (atom->avec->angles_allow)
                  for (j = 0; j < atom->num_angle[i]; j++) {
                    atom->angle_atom1[i][j] += atom_offset;
                    atom->angle_atom2[i][j] += atom_offset;
                    atom->angle_atom3[i][j] += atom_offset;
                  }
                if (atom->avec->dihedrals_allow)
                  for (j = 0; j < atom->num_dihedral[i]; j++) {
                    atom->dihedral_atom1[i][j] += atom_offset;
                    atom->dihedral_atom2[i][j] += atom_offset;
                    atom->dihedral_atom3[i][j] += atom_offset;
                    atom->dihedral_atom4[i][j] += atom_offset;
                  }
                if (atom->avec->impropers_allow)
                  for (j = 0; j < atom->num_improper[i]; j++) {
                    atom->improper_atom1[i][j] += atom_offset;
                    atom->improper_atom2[i][j] += atom_offset;
                    atom->improper_atom3[i][j] += atom_offset;
                    atom->improper_atom4[i][j] += atom_offset;
                  }
              }
            } else m += static_cast<int> (buf[m]);
          }
        }
      }
    }
  }

  // free communication buffer and old atom class

  memory->destroy(buf);
  delete old;

  // check that all atoms were assigned to procs

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",natoms);
  }

  if (natoms != atom->natoms)
    error->all(FLERR,"Replicate did not assign all atoms correctly");

  if (me == 0) {
    if (atom->nbonds) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
    }
    if (atom->nangles) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " angles\n",
                          atom->nangles);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " angles\n",
                           atom->nangles);
    }
    if (atom->ndihedrals) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " dihedrals\n",
                          atom->ndihedrals);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " dihedrals\n",
                           atom->ndihedrals);
    }
    if (atom->nimpropers) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " impropers\n",
                          atom->nimpropers);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " impropers\n",
                           atom->nimpropers);
    }
  }

  // check that atom IDs are valid

  atom->tag_check();

  // if molecular system or user-requested, create global mapping of atoms

  if (atom->molecular || atom->map_user) {
    atom->map_init();
    atom->map_set();
  }

  // create special bond lists for molecular systems

  if (atom->molecular) {
    Special special(lmp);
    special.build();
  }
}

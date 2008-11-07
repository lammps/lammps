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
#define MAXATOMS  0x7FFFFFFF
#define EPSILON   1.0e-6

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Replicate::Replicate(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void Replicate::command(int narg, char **arg)
{
  int i,j,m,n;

  if (domain->box_exist == 0)
    error->all("Replicate command before simulation box is defined");
  if (narg != 3) error->all("Illegal replicate command");

  int me = comm->me;
  int nprocs = comm->nprocs;

  if (me == 0 && screen) fprintf(screen,"Replicating atoms ...\n");

  // nrep = total # of replications

  int nx = atoi(arg[0]);
  int ny = atoi(arg[1]);
  int nz = atoi(arg[2]);
  int nrep = nx*ny*nz;

  // error and warning checks

  if (nx <= 0 || ny <= 0 || nz <= 0) error->all("Illegal replicate command");
  if (domain->dimension == 2 && nz != 1)
    error->all("Cannot replicate 2d simulation in z dimension");
  if ((nx > 1 && domain->xperiodic == 0) || 
      (ny > 1 && domain->yperiodic == 0) ||
      (nz > 1 && domain->zperiodic == 0)) 
    error->warning("Replicating in a non-periodic dimension");

  if (atom->nextra_grow || atom->nextra_restart || atom->nextra_store)
    error->all("Cannot replicate with fixes that store atom quantities");

  // maxtag = largest atom tag across all existing atoms

  int maxtag = 0;
  for (i = 0; i < atom->nlocal; i++) maxtag = MAX(atom->tag[i],maxtag);
  int maxtag_all;
  MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_INT,MPI_MAX,world);
  maxtag = maxtag_all;

  // maxmol = largest molecule tag across all existing atoms

  int maxmol = 0;
  if (atom->molecular) {
    for (i = 0; i < atom->nlocal; i++) maxmol = MAX(atom->molecule[i],maxmol);
    int maxmol_all;
    MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_INT,MPI_MAX,world);
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

  double *buf = 
    (double *) memory->smalloc(max_size*sizeof(double),"replicate:buf");

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

  // check that new problem size will not be too large
  // if N > 2^31, turn off tags
  // if molecular, N/Nbonds/etc cannot be > 2^31 else tags/counts invalid

  double rep = nrep;
  if (rep*old->natoms > MAXATOMS) atom->tag_enable = 0;

  if (atom->molecular) {
    if (rep*old->natoms > MAXATOMS || rep*old->nbonds > MAXATOMS ||
	rep*old->nangles > MAXATOMS || rep*old->ndihedrals > MAXATOMS ||
	rep*old->nimpropers > MAXATOMS)
      error->all("Too big a problem to replicate with molecular atom style");
  }

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
  comm->set_procs();
  domain->set_local_box();

  // copy type arrays to new atom class

  if (atom->mass) {
    for (int itype = 1; itype <= atom->ntypes; itype++) {
      atom->mass_setflag[itype] = old->mass_setflag[itype];
      if (atom->mass_setflag[itype]) atom->mass[itype] = old->mass[itype];
    }
  }

  if (atom->dipole) {
    for (int itype = 1; itype <= atom->ntypes; itype++) {
      atom->dipole_setflag[itype] = old->dipole_setflag[itype];
      if (atom->dipole_setflag[itype]) 
	atom->dipole[itype] = old->dipole[itype];
    }
  }

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all replicated atoms will be owned even with round-off

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
    if (comm->myloc[0] == 0) sublo[0] -= EPSILON;
    if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] += EPSILON;
  }
  if (domain->yperiodic) {
    if (comm->myloc[1] == 0) sublo[1] -= EPSILON;
    if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] += EPSILON;
  }
  if (domain->zperiodic) {
    if (comm->myloc[2] == 0) sublo[2] -= EPSILON;
    if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] += EPSILON;
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

  int ix,iy,iz,image,atom_offset,mol_offset;
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
	    image = (512 << 20) | (512 << 10) | 512;
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

	      if (atom->molecular) {
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
  } // end of proc loop

  // free communication buffer and old atom class

  memory->sfree(buf);
  delete old;

  // check that all atoms were assigned to procs

  double natoms;
  double rlocal = atom->nlocal;
  MPI_Allreduce(&rlocal,&natoms,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  %.15g atoms\n",natoms);
    if (logfile) fprintf(logfile,"  %.15g atoms\n",natoms);
  }

  if (natoms != atom->natoms)
    error->all("Replicate did not assign all atoms correctly");

  if (me == 0) {
    if (atom->nbonds) {
      if (screen) fprintf(screen,"  %d bonds\n",atom->nbonds);
      if (logfile) fprintf(logfile,"  %d bonds\n",atom->nbonds);
    }
    if (atom->nangles) {
      if (screen) fprintf(screen,"  %d angles\n",atom->nangles);
      if (logfile) fprintf(logfile,"  %d angles\n",atom->nangles);
    }
    if (atom->ndihedrals) {
      if (screen) fprintf(screen,"  %d dihedrals\n",atom->ndihedrals);
      if (logfile) fprintf(logfile,"  %d dihedrals\n",atom->ndihedrals);
    }
    if (atom->nimpropers) {
      if (screen) fprintf(screen,"  %d impropers\n",atom->nimpropers);
      if (logfile) fprintf(logfile,"  %d impropers\n",atom->nimpropers);
    }
  }

  // create global mapping and bond topology now that system is defined

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
  if (atom->molecular) {
    Special special(lmp);
    special.build();
  }
}

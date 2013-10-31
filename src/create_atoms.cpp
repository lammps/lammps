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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_atoms.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "random_park.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

#define BIG 1.0e30
#define EPSILON 1.0e-6

enum{BOX,REGION,SINGLE,RANDOM};

/* ---------------------------------------------------------------------- */

CreateAtoms::CreateAtoms(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void CreateAtoms::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Create_atoms command before simulation box is defined");
  if (modify->nfix_restart_peratom)
    error->all(FLERR,"Cannot create_atoms after "
               "reading restart file with per-atom info");

  // parse arguments

  if (narg < 2) error->all(FLERR,"Illegal create_atoms command");
  itype = force->inumeric(FLERR,arg[0]);
  if (itype <= 0 || itype > atom->ntypes)
    error->all(FLERR,"Invalid atom type in create_atoms command");

  int iarg;
  if (strcmp(arg[1],"box") == 0) {
    style = BOX;
    iarg = 2;
  } else if (strcmp(arg[1],"region") == 0) {
    style = REGION;
    if (narg < 3) error->all(FLERR,"Illegal create_atoms command");
    nregion = domain->find_region(arg[2]);
    if (nregion == -1) error->all(FLERR,
                                  "Create_atoms region ID does not exist");
    domain->regions[nregion]->init();
    iarg = 3;;
  } else if (strcmp(arg[1],"single") == 0) {
    style = SINGLE;
    if (narg < 5) error->all(FLERR,"Illegal create_atoms command");
    xone[0] = force->numeric(FLERR,arg[2]);
    xone[1] = force->numeric(FLERR,arg[3]);
    xone[2] = force->numeric(FLERR,arg[4]);
    iarg = 5;
  } else if (strcmp(arg[1],"random") == 0) {
    style = RANDOM;
    if (narg < 5) error->all(FLERR,"Illegal create_atoms command");
    nrandom = force->inumeric(FLERR,arg[2]);
    seed = force->inumeric(FLERR,arg[3]);
    if (strcmp(arg[4],"NULL") == 0) nregion = -1;
    else {
      nregion = domain->find_region(arg[4]);
      if (nregion == -1) error->all(FLERR,
                                    "Create_atoms region ID does not exist");
      domain->regions[nregion]->init();
    }
    iarg = 5;
  } else error->all(FLERR,"Illegal create_atoms command");

  // process optional keywords

  int scaleflag = 1;
  remapflag = 0;

  nbasis = domain->lattice->nbasis;
  basistype = new int[nbasis];
  for (int i = 0; i < nbasis; i++) basistype[i] = itype;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal create_atoms command");
      int ibasis = force->inumeric(FLERR,arg[iarg+1]);
      itype = force->inumeric(FLERR,arg[iarg+2]);
      if (ibasis <= 0 || ibasis > nbasis ||
          itype <= 0 || itype > atom->ntypes)
        error->all(FLERR,"Invalid basis setting in create_atoms command");
      basistype[ibasis-1] = itype;
      iarg += 3;
    } else if (strcmp(arg[iarg],"remap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_atoms command");
      if (strcmp(arg[iarg+1],"yes") == 0) remapflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) remapflag = 0;
      else error->all(FLERR,"Illegal create_atoms command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal create_atoms command");
      iarg += 2;
    } else error->all(FLERR,"Illegal create_atoms command");
  }

  // error checks

  if (style == RANDOM) {
    if (nrandom < 0) error->all(FLERR,"Illegal create_atoms command");
    if (seed <= 0) error->all(FLERR,"Illegal create_atoms command");
  }

  // demand non-none lattice be defined for BOX and REGION
  // else setup scaling for SINGLE and RANDOM
  // could use domain->lattice->lattice2box() to do conversion of
  //   lattice to box, but not consistent with other uses of units=lattice
  // triclinic remapping occurs in add_single()

  if (style == BOX || style == REGION) {
    if (nbasis == 0)
      error->all(FLERR,"Cannot create atoms with undefined lattice");
  } else if (scaleflag == 1) {
    xone[0] *= domain->lattice->xlattice;
    xone[1] *= domain->lattice->ylattice;
    xone[2] *= domain->lattice->zlattice;
  }

  // set bounds for my proc in sublo[3] & subhi[3]
  // if periodic and style = BOX or REGION, i.e. using lattice:
  //   should create exactly 1 atom when 2 images are both "on" the boundary
  //   either image may be slightly inside/outside true box due to round-off
  //   if I am lo proc, decrement lower bound by EPSILON
  //     this will insure lo image is created
  //   if I am hi proc, decrement upper bound by 2.0*EPSILON
  //     this will insure hi image is not created
  //   thus insertion box is EPSILON smaller than true box
  //     and is shifted away from true boundary
  //     which is where atoms are likely to be generated

  triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  if (triclinic == 0) {
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
  }

  if (style == BOX || style == REGION) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0]-1) subhi[0] -= 2.0*epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1]-1) subhi[1] -= 2.0*epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2]-1) subhi[2] -= 2.0*epsilon[2];
    }
  }

  // add atoms in one of 3 ways

  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;

  if (style == SINGLE) add_single();
  else if (style == RANDOM) add_random();
  else add_lattice();

  // invoke set_arrays() for fixes that need initialization of new atoms

  int nlocal = atom->nlocal;
  for (int m = 0; m < modify->nfix; m++) {
    Fix *fix = modify->fix[m];
    if (fix->create_attribute)
      for (int i = nlocal_previous; i < nlocal; i++)
        fix->set_arrays(i);
  }

  // clean up

  if (domain->lattice) delete [] basistype;

  // new total # of atoms

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
    error->all(FLERR,"Too many total atoms");

  // print status

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Created " BIGINT_FORMAT " atoms\n",
              atom->natoms-natoms_previous);
    if (logfile)
      fprintf(logfile,"Created " BIGINT_FORMAT " atoms\n",
              atom->natoms-natoms_previous);
  }

  // reset simulation now that more atoms are defined
  // add tags for newly created atoms if possible
  // if global map exists, reset it
  // change these to MAXTAGINT when allow tagint = bigint

  if (atom->natoms > MAXSMALLINT) {
    if (comm->me == 0) 
      error->warning(FLERR,"Total atom count exceeds ID limit, "
                     "atoms will not have individual IDs");
    atom->tag_enable = 0;
  }
  if (atom->natoms <= MAXSMALLINT) atom->tag_extend();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // if a molecular system, set nspecial to 0 for new atoms
  // NOTE: 31May12, don't think this is needed, avec->create_atom() does it

  //if (atom->molecular) {
  //  int **nspecial = atom->nspecial;
  //  for (int i = nlocal_previous; i < atom->nlocal; i++) {
  //    nspecial[i][0] = 0;
  //    nspecial[i][1] = 0;
  //    nspecial[i][2] = 0;
  //  }
  //}
}

/* ----------------------------------------------------------------------
   add single atom with coords at xone if it's in my sub-box
   if triclinic, xone is in lamda coords
------------------------------------------------------------------------- */

void CreateAtoms::add_single()
{
  // remap atom if requested

  if (remapflag) {
    tagint imagetmp = ((tagint) IMGMAX << IMG2BITS) |
      ((tagint) IMGMAX << IMGBITS) | IMGMAX;
    domain->remap(xone,imagetmp);
  }

  // if triclinic, convert to lamda coords (0-1)

  double lamda[3],*coord;
  if (triclinic) {
    domain->x2lamda(xone,lamda);
    coord = lamda;
  } else coord = xone;

  // if atom is in my subbox, create it

  if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
      coord[1] >= sublo[1] && coord[1] < subhi[1] &&
      coord[2] >= sublo[2] && coord[2] < subhi[2])
    atom->avec->create_atom(itype,xone);
}


/* ----------------------------------------------------------------------
   add Nrandom atoms at random locations
------------------------------------------------------------------------- */

void CreateAtoms::add_random()
{
  double xlo,ylo,zlo,xhi,yhi,zhi,zmid;
  double lamda[3],*coord;
  double *boxlo,*boxhi;

  // random number generator, same for all procs

  RanPark *random = new RanPark(lmp,seed);

  // bounding box for atom creation
  // in real units, even if triclinic
  // only limit bbox by region if its bboxflag is set (interior region)

  if (triclinic == 0) {
    xlo = domain->boxlo[0]; xhi = domain->boxhi[0];
    ylo = domain->boxlo[1]; yhi = domain->boxhi[1];
    zlo = domain->boxlo[2]; zhi = domain->boxhi[2];
    zmid = zlo + 0.5*(zhi-zlo);
  } else {
    xlo = domain->boxlo_bound[0]; xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1]; yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2]; zhi = domain->boxhi_bound[2];
    zmid = zlo + 0.5*(zhi-zlo);
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
  }

  if (nregion >= 0 && domain->regions[nregion]->bboxflag) {
    xlo = MAX(xlo,domain->regions[nregion]->extent_xlo);
    xhi = MIN(xhi,domain->regions[nregion]->extent_xhi);
    ylo = MAX(ylo,domain->regions[nregion]->extent_ylo);
    yhi = MIN(yhi,domain->regions[nregion]->extent_yhi);
    zlo = MAX(zlo,domain->regions[nregion]->extent_zlo);
    zhi = MIN(zhi,domain->regions[nregion]->extent_zhi);
  }

  // generate random positions for each new atom within bounding box
  // iterate until atom is within region and triclinic simulation box
  // if final atom position is in my subbox, create it

  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all(FLERR,"No overlap of box and region for create_atoms");

  int valid;
  for (int i = 0; i < nrandom; i++) {
    while (1) {
      xone[0] = xlo + random->uniform() * (xhi-xlo);
      xone[1] = ylo + random->uniform() * (yhi-ylo);
      xone[2] = zlo + random->uniform() * (zhi-zlo);
      if (domain->dimension == 2) xone[2] = zmid;

      valid = 1;
      if (nregion >= 0 &&
          domain->regions[nregion]->match(xone[0],xone[1],xone[2]) == 0)
        valid = 0;
      if (triclinic) {
        domain->x2lamda(xone,lamda);
        coord = lamda;
        if (coord[0] < boxlo[0] || coord[0] >= boxhi[0] ||
            coord[1] < boxlo[1] || coord[1] >= boxhi[1] ||
            coord[2] < boxlo[2] || coord[2] >= boxhi[2]) valid = 0;
      } else coord = xone;

      if (valid) break;
    }

    // if triclinic, coord is now in lamda units

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2])
      atom->avec->create_atom(itype,xone);
  }

  // clean-up

  delete random;
}

/* ----------------------------------------------------------------------
   add many atoms by looping over lattice
------------------------------------------------------------------------- */

void CreateAtoms::add_lattice()
{
  // convert 8 corners of my subdomain from box coords to lattice coords
  // for orthogonal, use corner pts of my subbox
  // for triclinic, use bounding box of my subbox
  // xyz min to max = bounding box around the domain corners in lattice space

  double bboxlo[3],bboxhi[3];

  if (triclinic == 0) {
    bboxlo[0] = domain->sublo[0]; bboxhi[0] = domain->subhi[0];
    bboxlo[1] = domain->sublo[1]; bboxhi[1] = domain->subhi[1];
    bboxlo[2] = domain->sublo[2]; bboxhi[2] = domain->subhi[2];
  } else domain->bbox(domain->sublo_lamda,domain->subhi_lamda,bboxlo,bboxhi);

  double xmin,ymin,zmin,xmax,ymax,zmax;
  xmin = ymin = zmin = BIG;
  xmax = ymax = zmax = -BIG;

  domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxlo[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxlo[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxlo[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxlo[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxhi[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxhi[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxhi[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);
  domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxhi[2],
                        xmin,ymin,zmin,xmax,ymax,zmax);

  // ilo:ihi,jlo:jhi,klo:khi = loop bounds for lattice overlap of my subbox
  // overlap = any part of a unit cell (face,edge,pt) in common with my subbox
  // in lattice space, subbox is a tilted box
  // but bbox of subbox is aligned with lattice axes
  // so ilo:khi unit cells should completely tile bounding box
  // decrement lo, increment hi to avoid round-off issues in lattice->bbox(),
  //   which can lead to missing atoms in rare cases
  // extra decrement of lo if min < 0, since static_cast(-1.5) = -1

  int ilo,ihi,jlo,jhi,klo,khi;
  ilo = static_cast<int> (xmin) - 1;
  jlo = static_cast<int> (ymin) - 1;
  klo = static_cast<int> (zmin) - 1;
  ihi = static_cast<int> (xmax) + 1;
  jhi = static_cast<int> (ymax) + 1;
  khi = static_cast<int> (zmax) + 1;

  if (xmin < 0.0) ilo--;
  if (ymin < 0.0) jlo--;
  if (zmin < 0.0) klo--;

  // iterate on 3d periodic lattice of unit cells using loop bounds
  // iterate on nbasis atoms in each unit cell
  // convert lattice coords to box coords
  // add atom if it meets all criteria

  double **basis = domain->lattice->basis;
  double x[3],lamda[3];
  double *coord;

  int i,j,k,m;
  for (k = klo; k <= khi; k++)
    for (j = jlo; j <= jhi; j++)
      for (i = ilo; i <= ihi; i++)
        for (m = 0; m < nbasis; m++) {

          x[0] = i + basis[m][0];
          x[1] = j + basis[m][1];
          x[2] = k + basis[m][2];

          // convert from lattice coords to box coords

          domain->lattice->lattice2box(x[0],x[1],x[2]);

          // if a region was specified, test if atom is in it

          if (style == REGION)
            if (!domain->regions[nregion]->match(x[0],x[1],x[2])) continue;

          // test if atom is in my subbox

          if (triclinic) {
            domain->x2lamda(x,lamda);
            coord = lamda;
          } else coord = x;

          if (coord[0] < sublo[0] || coord[0] >= subhi[0] ||
              coord[1] < sublo[1] || coord[1] >= subhi[1] ||
              coord[2] < sublo[2] || coord[2] >= subhi[2]) continue;

          // add the atom to my list of atoms

          atom->avec->create_atom(basistype[m],x);
        }
}

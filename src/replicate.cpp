// clang-format off
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

#include "replicate.h"

#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "special.h"
#include "label_map.h"

#include <cstring>

using namespace LAMMPS_NS;

#define LB_FACTOR 1.1
#define EPSILON   1.0e-6

/* ---------------------------------------------------------------------- */

Replicate::Replicate(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

void Replicate::command(int narg, char **arg)
{
  int i,j,m,n;

  if (domain->box_exist == 0)
    error->all(FLERR,"Replicate command before simulation box is defined");
  if (narg < 3 || narg > 4) error->all(FLERR,"Illegal replicate command");

  int me = comm->me;
  int nprocs = comm->nprocs;

  // nrep = total # of replications

  int nx = utils::inumeric(FLERR,arg[0],false,lmp);
  int ny = utils::inumeric(FLERR,arg[1],false,lmp);
  int nz = utils::inumeric(FLERR,arg[2],false,lmp);
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
    error->all(FLERR, "Illegal replication grid {}x{}x{}. All replications must be > 0",
               nx, ny, nz);

  int nrep = nx*ny*nz;
  if (me == 0)
    utils::logmesg(lmp, "Replication is creating a {}x{}x{} = {} times larger system...\n",
                   nx, ny, nz, nrep);

  int bbox_flag = 0;
  if (narg == 4)
    if (strcmp(arg[3],"bbox") == 0) bbox_flag = 1;

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

  // record wall time for atom replication

  MPI_Barrier(world);
  double time1 = platform::walltime();

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

  // check image flags maximum extent
  // only efficient small image flags compared to new system

  int _imagelo[3], _imagehi[3];
  _imagelo[0] = 0;
  _imagelo[1] = 0;
  _imagelo[2] = 0;
  _imagehi[0] = 0;
  _imagehi[1] = 0;
  _imagehi[2] = 0;

  if (bbox_flag) {

    for (i=0; i<atom->nlocal; ++i) {
      imageint image = atom->image[i];
      int xbox = (image & IMGMASK) - IMGMAX;
      int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (image >> IMG2BITS) - IMGMAX;

      if (xbox < _imagelo[0]) _imagelo[0] = xbox;
      if (ybox < _imagelo[1]) _imagelo[1] = ybox;
      if (zbox < _imagelo[2]) _imagelo[2] = zbox;

      if (xbox > _imagehi[0]) _imagehi[0] = xbox;
      if (ybox > _imagehi[1]) _imagehi[1] = ybox;
      if (zbox > _imagehi[2]) _imagehi[2] = zbox;
    }

    MPI_Allreduce(MPI_IN_PLACE, &(_imagelo[0]), 3, MPI_INT, MPI_MIN, world);
    MPI_Allreduce(MPI_IN_PLACE, &(_imagehi[0]), 3, MPI_INT, MPI_MAX, world);
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
  // also set atomKK for Kokkos version of Atom class

  Atom *old = atom;
  atomKK = nullptr;
  if (lmp->kokkos) atom = atomKK = new AtomKokkos(lmp);
  else atom = new Atom(lmp);

  atom->settings(old);

  // transfer molecule templates. needs to be done early for atom style template

  if (old->nmolecule) {
    atom->molecules = (Molecule **) memory->smalloc((old->nmolecule)*sizeof(Molecule *),
                                                    "atom::molecules");
    atom->nmolecule = old->nmolecule;
    for (i = 0; i < old->nmolecule; ++i)
      atom->molecules[i] = old->molecules[i];
    memory->sfree(old->molecules);
    old->molecules = nullptr;
    old->nmolecule = 0;
  }

  atom->create_avec(old->atom_style,old->avec->nargcopy,old->avec->argcopy,0);

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

  atom->nellipsoids = old->nellipsoids * nrep;
  atom->nlines = old->nlines * nrep;
  atom->ntris = old->ntris * nrep;
  atom->nbodies = old->nbodies * nrep;

  atom->ntypes = old->ntypes;
  atom->nbondtypes = old->nbondtypes;
  atom->nangletypes = old->nangletypes;
  atom->ndihedraltypes = old->ndihedraltypes;
  atom->nimpropertypes = old->nimpropertypes;

  atom->bond_per_atom = old->bond_per_atom;
  atom->angle_per_atom = old->angle_per_atom;
  atom->dihedral_per_atom = old->dihedral_per_atom;
  atom->improper_per_atom = old->improper_per_atom;

  atom->extra_bond_per_atom = old->extra_bond_per_atom;
  atom->extra_angle_per_atom = old->extra_angle_per_atom;
  atom->extra_dihedral_per_atom = old->extra_dihedral_per_atom;
  atom->extra_improper_per_atom = old->extra_improper_per_atom;
  atom->maxspecial = old->maxspecial;

  if (old->labelmapflag) {
    atom->add_label_map();
    atom->lmap->merge_lmap(old->lmap, Atom::ATOM);
    atom->lmap->merge_lmap(old->lmap, Atom::BOND);
    atom->lmap->merge_lmap(old->lmap, Atom::ANGLE);
    atom->lmap->merge_lmap(old->lmap, Atom::DIHEDRAL);
    atom->lmap->merge_lmap(old->lmap, Atom::IMPROPER);
  }

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

  // allocate atom arrays to size N, rounded up by AtomVec->DELTA

  bigint nbig = n;
  nbig = atom->avec->roundup(nbig);
  n = static_cast<int> (nbig);
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
  // ensures all replicated atoms will be owned even with round-off

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

  if (comm->layout != Comm::LAYOUT_TILED) {
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

  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] += epsilon[2];
    }
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

  if (bbox_flag) {

    // allgather size of buf on each proc

    n = 0;
    for (i = 0; i < old->nlocal; i++) n += old_avec->pack_restart(i,&buf[n]);

    int * size_buf_rnk;
    memory->create(size_buf_rnk, nprocs, "replicate:size_buf_rnk");

    MPI_Allgather(&n, 1, MPI_INT, size_buf_rnk, 1, MPI_INT, world);

    // size of buf_all

    int size_buf_all = 0;
    MPI_Allreduce(&n, &size_buf_all, 1, MPI_INT, MPI_SUM, world);

    if (me == 0) {
      auto mesg = fmt::format("  bounding box image = ({} {} {}) "
                              "to ({} {} {})\n",
                              _imagelo[0],_imagelo[1],_imagelo[2],
                              _imagehi[0],_imagehi[1],_imagehi[2]);
      mesg += fmt::format("  bounding box extra memory = {:.2f} MB\n",
                          (double)size_buf_all*sizeof(double)/1024/1024);
      utils::logmesg(lmp,mesg);
    }

    // rnk offsets

    int *disp_buf_rnk;
    memory->create(disp_buf_rnk, nprocs, "replicate:disp_buf_rnk");
    disp_buf_rnk[0] = 0;
    for (i = 1; i < nprocs; i++)
      disp_buf_rnk[i] = disp_buf_rnk[i-1] + size_buf_rnk[i-1];

    // allgather buf_all

    double * buf_all;
    memory->create(buf_all, size_buf_all, "replicate:buf_all");

    MPI_Allgatherv(buf,n,MPI_DOUBLE,buf_all,size_buf_rnk,disp_buf_rnk,
                   MPI_DOUBLE,world);

    // bounding box of original unwrapped system

    double _orig_lo[3], _orig_hi[3];
    if (triclinic) {
      _orig_lo[0] = domain->boxlo[0] +
        _imagelo[0] * old_xprd + _imagelo[1] * old_xy + _imagelo[2] * old_xz;
      _orig_lo[1] = domain->boxlo[1] +
        _imagelo[1] * old_yprd + _imagelo[2] * old_yz;
      _orig_lo[2] = domain->boxlo[2] + _imagelo[2] * old_zprd;

      _orig_hi[0] = domain->boxlo[0] +
        (_imagehi[0]+1) * old_xprd +
        (_imagehi[1]+1) * old_xy + (_imagehi[2]+1) * old_xz;
      _orig_hi[1] = domain->boxlo[1] +
        (_imagehi[1]+1) * old_yprd + (_imagehi[2]+1) * old_yz;
      _orig_hi[2] = domain->boxlo[2] + (_imagehi[2]+1) * old_zprd;
    } else {
      _orig_lo[0] = domain->boxlo[0] + _imagelo[0] * old_xprd;
      _orig_lo[1] = domain->boxlo[1] + _imagelo[1] * old_yprd;
      _orig_lo[2] = domain->boxlo[2] + _imagelo[2] * old_zprd;

      _orig_hi[0] = domain->boxlo[0] + (_imagehi[0]+1) * old_xprd;
      _orig_hi[1] = domain->boxlo[1] + (_imagehi[1]+1) * old_yprd;
      _orig_hi[2] = domain->boxlo[2] + (_imagehi[2]+1) * old_zprd;
    }

    double _lo[3], _hi[3];

    int num_replicas_added = 0;

    for (ix = 0; ix < nx; ix++) {
      for (iy = 0; iy < ny; iy++) {
        for (iz = 0; iz < nz; iz++) {

          // domain->remap() overwrites coordinates, so always recompute here

          if (triclinic) {
            _lo[0] = _orig_lo[0] + ix * old_xprd + iy * old_xy + iz * old_xz;
            _hi[0] = _orig_hi[0] + ix * old_xprd + iy * old_xy + iz * old_xz;

            _lo[1] = _orig_lo[1] + iy * old_yprd + iz * old_yz;
            _hi[1] = _orig_hi[1] + iy * old_yprd + iz * old_yz;

            _lo[2] = _orig_lo[2] + iz * old_zprd;
            _hi[2] = _orig_hi[2] + iz * old_zprd;
          } else {
            _lo[0] = _orig_lo[0] + ix * old_xprd;
            _hi[0] = _orig_hi[0] + ix * old_xprd;

            _lo[1] = _orig_lo[1] + iy * old_yprd;
            _hi[1] = _orig_hi[1] + iy * old_yprd;

            _lo[2] = _orig_lo[2] + iz * old_zprd;
            _hi[2] = _orig_hi[2] + iz * old_zprd;
          }

          // test if bounding box of shifted replica overlaps sub-domain of proc
          // if not, then skip testing atoms

          int xoverlap = 1;
          int yoverlap = 1;
          int zoverlap = 1;
          if (triclinic) {
            double _llo[3];
            domain->x2lamda(_lo,_llo);
            double _lhi[3];
            domain->x2lamda(_hi,_lhi);

            if (_llo[0] > (subhi[0] - EPSILON)
                || _lhi[0] < (sublo[0] + EPSILON) ) xoverlap = 0;
            if (_llo[1] > (subhi[1] - EPSILON)
                || _lhi[1] < (sublo[1] + EPSILON) ) yoverlap = 0;
            if (_llo[2] > (subhi[2] - EPSILON)
                || _lhi[2] < (sublo[2] + EPSILON) ) zoverlap = 0;
          } else {
            if (_lo[0] > (subhi[0] - EPSILON)
                || _hi[0] < (sublo[0] + EPSILON) ) xoverlap = 0;
            if (_lo[1] > (subhi[1] - EPSILON)
                || _hi[1] < (sublo[1] + EPSILON) ) yoverlap = 0;
            if (_lo[2] > (subhi[2] - EPSILON)
                || _hi[2] < (sublo[2] + EPSILON) ) zoverlap = 0;
          }

          int overlap = 0;
          if (xoverlap && yoverlap && zoverlap) overlap = 1;

          // if no overlap, test if bounding box wrapped back into new system

          if (!overlap) {

            // wrap back into cell

            imageint imagelo = ((imageint) IMGMAX << IMG2BITS) |
              ((imageint) IMGMAX << IMGBITS) | IMGMAX;
            domain->remap(&(_lo[0]), imagelo);
            int xboxlo = (imagelo & IMGMASK) - IMGMAX;
            int yboxlo = (imagelo >> IMGBITS & IMGMASK) - IMGMAX;
            int zboxlo = (imagelo >> IMG2BITS) - IMGMAX;

            imageint imagehi = ((imageint) IMGMAX << IMG2BITS) |
              ((imageint) IMGMAX << IMGBITS) | IMGMAX;
            domain->remap(&(_hi[0]), imagehi);
            int xboxhi = (imagehi & IMGMASK) - IMGMAX;
            int yboxhi = (imagehi >> IMGBITS & IMGMASK) - IMGMAX;
            int zboxhi = (imagehi >> IMG2BITS) - IMGMAX;

            if (triclinic) {
              double _llo[3];
              _llo[0] = _lo[0]; _llo[1] = _lo[1];  _llo[2] = _lo[2];
              domain->x2lamda(_llo,_lo);

              double _lhi[3];
              _lhi[0] = _hi[0]; _lhi[1] = _hi[1];  _lhi[2] = _hi[2];
              domain->x2lamda(_lhi,_hi);
            }

            // test all fragments for any overlap; ok to include false positives

            int _xoverlap1 = 0;
            int _xoverlap2 = 0;
            if (!xoverlap) {
              if (xboxlo < 0) {
                _xoverlap1 = 1;
                if (_lo[0] > (subhi[0] - EPSILON)) _xoverlap1 = 0;
              }

              if (xboxhi > 0) {
                _xoverlap2 = 1;
                if (_hi[0] < (sublo[0] + EPSILON)) _xoverlap2 = 0;
              }

              if (_xoverlap1 || _xoverlap2) xoverlap = 1;
            }

            int _yoverlap1 = 0;
            int _yoverlap2 = 0;
            if (!yoverlap) {
              if (yboxlo < 0) {
                _yoverlap1 = 1;
                if (_lo[1] > (subhi[1] - EPSILON)) _yoverlap1 = 0;
              }

              if (yboxhi > 0) {
                _yoverlap2 = 1;
                if (_hi[1] < (sublo[1] + EPSILON)) _yoverlap2 = 0;
              }

              if (_yoverlap1 || _yoverlap2) yoverlap = 1;
            }


            int _zoverlap1 = 0;
            int _zoverlap2 = 0;
            if (!zoverlap) {
              if (zboxlo < 0) {
                _zoverlap1 = 1;
                if (_lo[2] > (subhi[2] - EPSILON)) _zoverlap1 = 0;
              }

              if (zboxhi > 0) {
                _zoverlap2 = 1;
                if (_hi[2] < (sublo[2] + EPSILON)) _zoverlap2 = 0;
              }

              if (_zoverlap1 || _zoverlap2) zoverlap = 1;
            }

            // does either fragment overlap w/ sub-domain

            if (xoverlap && yoverlap && zoverlap) overlap = 1;
          }

          // while loop over one proc's atom list

          if (overlap) {
            num_replicas_added++;

            m = 0;
            while (m < size_buf_all) {
              image = ((imageint) IMGMAX << IMG2BITS) |
                ((imageint) IMGMAX << IMGBITS) | IMGMAX;
              if (triclinic == 0) {
                x[0] = buf_all[m+1] + ix*old_xprd;
                x[1] = buf_all[m+2] + iy*old_yprd;
                x[2] = buf_all[m+3] + iz*old_zprd;
              } else {
                x[0] = buf_all[m+1] + ix*old_xprd + iy*old_xy + iz*old_xz;
                x[1] = buf_all[m+2] + iy*old_yprd + iz*old_yz;
                x[2] = buf_all[m+3] + iz*old_zprd;
              }
              domain->remap(x,image);
              if (triclinic) {
                domain->x2lamda(x,lamda);
                coord = lamda;
              } else coord = x;

              if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
                  coord[1] >= sublo[1] && coord[1] < subhi[1] &&
                  coord[2] >= sublo[2] && coord[2] < subhi[2]) {

                m += avec->unpack_restart(&buf_all[m]);

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

                if (atom->molecular != Atom::ATOMIC) {
                  if (atom->molecule[i] > 0)
                    atom->molecule[i] += mol_offset;
                  if (atom->molecular == Atom::MOLECULAR) {
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
                }
              } else m += static_cast<int> (buf_all[m]);
            }
          } // if (overlap)

        }
      }
    }

    memory->destroy(size_buf_rnk);
    memory->destroy(disp_buf_rnk);
    memory->destroy(buf_all);

    int sum = 0;
    MPI_Reduce(&num_replicas_added, &sum, 1, MPI_INT, MPI_SUM, 0, world);
    double avg = (double) sum / nprocs;
    if (me == 0)
      utils::logmesg(lmp,"  average # of replicas added to proc = {:.2f} out "
                     "of {} ({:.2f}%)\n",avg,nx*ny*nz,avg/(nx*ny*nz)*100.0);
  } else {

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

                if (atom->molecular != Atom::ATOMIC) {
                  if (atom->molecule[i] > 0)
                    atom->molecule[i] += mol_offset;
                  if (atom->molecular == Atom::MOLECULAR) {
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
                }
              } else m += static_cast<int> (buf[m]);
            }
          }
        }
      }
    }
  } // if (bbox_flag)

  // free communication buffer and old atom class

  memory->destroy(buf);
  delete old;

  // check that all atoms were assigned to procs

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (me == 0) {
    utils::logmesg(lmp,"  {} atoms\n",natoms);
  }

  if (natoms != atom->natoms)
    error->all(FLERR,"Replicate did not assign all atoms correctly");

  if (me == 0) {
    const char *molstyle = "";
    if (atom->molecular == Atom::TEMPLATE) molstyle = "template ";
    if (atom->nbonds) {
      utils::logmesg(lmp,"  {} {}bonds\n",atom->nbonds,molstyle);
    }
    if (atom->nangles) {
      utils::logmesg(lmp,"  {} {}angles\n",atom->nangles,molstyle);
    }
    if (atom->ndihedrals) {
      utils::logmesg(lmp,"  {} {}dihedrals\n",atom->ndihedrals,molstyle);
    }
    if (atom->nimpropers) {
      utils::logmesg(lmp,"  {} {}impropers\n",atom->nimpropers,molstyle);
    }
  }

  // check that atom IDs are valid

  atom->tag_check();

  // create global mapping of atoms

  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }

  // create special bond lists for molecular systems

  if (atom->molecular == Atom::MOLECULAR) {
    Special special(lmp);
    special.build();
  }

  // total time

  MPI_Barrier(world);

  if (me == 0)
    utils::logmesg(lmp,"  replicate CPU = {:.3f} seconds\n",platform::walltime()-time1);
}

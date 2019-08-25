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

/* ----------------------------------------------------------------------
   Contributing authors: Lauren Abbott (Sandia)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "fix_mscg.h"
#include "mscg.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "region.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMSCG::FixMSCG(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix mscg command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix mscg command");

  me = comm->me;
  nprocs = comm->nprocs;
  if (nprocs > 1) error->all(FLERR,"Fix mscg does not yet support "
                             "parallel use via MPI");

  if (sizeof(tagint) != sizeof(smallint))
    error->all(FLERR,"Fix mscg must be used with 32-bit atom IDs");

  // initialize

  int natoms = atom->natoms;
  int ntypes = atom->ntypes;

  max_partners_bond = 4;
  max_partners_angle = 12;
  max_partners_dihedral = 36;
  nframes = n_frames = block_size = 0;
  range_flag = 0;
  name_flag = 0;
  f = NULL;

  type_names = new char*[natoms];
  for (int i = 0; i < natoms; i++) type_names[i] = new char[24];

  // parse remaining arguments

  int iarg = 4;
  while(iarg < narg) {
    if (strcmp(arg[iarg],"range") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix mscg command");
      if (strcmp(arg[iarg+1],"on") == 0)
        range_flag = 1;
      else if (strcmp(arg[iarg+1],"off") == 0)
        range_flag = 0;
      else error->all(FLERR,"Illegal fix mscg command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"name") == 0) {
      if (iarg+ntypes+1 > narg)
        error->all(FLERR,"Illegal fix mscg command");
      name_flag = 1;
      for (int i = 0; i < ntypes; i++) {
        iarg += 1;
        std::string str = arg[iarg];
        type_names[i] = strcat(strdup(str.c_str()),"\0");
      }
      iarg += 1;
    } else if (strcmp(arg[iarg],"max") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix mscg command");
      max_partners_bond = atoi(arg[iarg+1]);
      max_partners_angle = atoi(arg[iarg+2]);
      max_partners_dihedral = atoi(arg[iarg+3]);
      iarg += 4;
    } else error->all(FLERR,"Illegal fix mscg command");
  }

  if (name_flag == 0) {
    for (int i = 0; i < natoms; i++) {
      std::string str = std::to_string(i+1);
      type_names[i] = strcat(strdup(str.c_str()),"\0");
    }
  }
}

/* ---------------------------------------------------------------------- */

FixMSCG::~FixMSCG()
{
  memory->destroy(f);
}

/* ---------------------------------------------------------------------- */

int FixMSCG::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMSCG::post_constructor()
{
  if (domain->triclinic == 1)
    error->all(FLERR,"Fix mscg does not yet support triclinic geometries");

  // topology information
  // sort by atom id to send to mscg lib

  int natoms = atom->natoms;
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  int *num_angle = atom->num_angle;
  int **angle_atom1 = atom->angle_atom1;
  int **angle_atom3 = atom->angle_atom3;
  int *num_dihedral = atom->num_dihedral;
  int **dihedral_atom1 = atom->dihedral_atom1;
  int **dihedral_atom3 = atom->dihedral_atom3;
  int **dihedral_atom4 = atom->dihedral_atom4;
  double *prd_half = domain->prd_half;
  int i,ii,j,jj,jnum,k,l;

  n_cg_sites = natoms;
  n_cg_types = atom->ntypes;

  memory->grow(f,nlocal,3,"fix_mscg:f");
  f1d = new double[n_cg_sites*3]();
  x1d = new double[n_cg_sites*3]();
  cg_site_types = new int[n_cg_sites]();
  n_partners_bond = new unsigned[n_cg_sites]();
  n_partners_angle = new unsigned[n_cg_sites]();
  n_partners_dihedral = new unsigned[n_cg_sites]();
  partners_bond = new unsigned*[n_cg_sites];
  for (i = 0; i < n_cg_sites; i++)
    partners_bond[i] = new unsigned[1*max_partners_bond]();
  partners_angle = new unsigned*[n_cg_sites];
  for (i = 0; i < n_cg_sites; i++)
    partners_angle[i] = new unsigned[2*max_partners_angle]();
  partners_dihedral = new unsigned*[n_cg_sites];
  for (i = 0; i < n_cg_sites; i++)
    partners_dihedral[i] = new unsigned[3*max_partners_dihedral]();

  for (i = 0; i < 3; i++)
    box_half_lengths[i] = prd_half[i];

  for (i = 0; i < nlocal; i++) {
    cg_site_types[i] = 0;
    n_partners_bond[i] = 0;
    n_partners_angle[i] = 0;
    n_partners_dihedral[i] = 0;
  }

  for (ii = 0; ii < nlocal; ii++) {
    i = tag[ii];
    cg_site_types[i-1] = type[ii];

    jnum = num_bond[ii];
    for (jj = 0; jj < jnum; jj++) {
      j = bond_atom[ii][jj];
      if (n_partners_bond[i-1] >= max_partners_bond ||
          n_partners_bond[j-1] >= max_partners_bond)
        error->all(FLERR,"Bond list overflow, boost fix_mscg max");
      partners_bond[i-1][n_partners_bond[i-1]] = j-1;
      partners_bond[j-1][n_partners_bond[j-1]] = i-1;
      n_partners_bond[i-1]++;
      n_partners_bond[j-1]++;
    }

    jnum = num_angle[ii];
    for (jj = 0; jj < jnum; jj++) {
      j = angle_atom1[ii][jj];
      k = angle_atom3[ii][jj];
      if (n_partners_angle[j-1] >= max_partners_angle ||
          n_partners_angle[k-1] >= max_partners_angle)
        error->all(FLERR,"Angle list overflow, boost fix_mscg max");
      partners_angle[j-1][n_partners_angle[j-1]*2] = i-1;
      partners_angle[j-1][n_partners_angle[j-1]*2+1] = k-1;
      partners_angle[k-1][n_partners_angle[k-1]*2] = i-1;
      partners_angle[k-1][n_partners_angle[k-1]*2+1] = j-1;
      n_partners_angle[j-1]++;
      n_partners_angle[k-1]++;
    }

    jnum = num_dihedral[ii];
    for (jj = 0; jj < jnum; jj++) {
      j = dihedral_atom1[ii][jj];
      k = dihedral_atom3[ii][jj];
      l = dihedral_atom4[ii][jj];
      if (n_partners_dihedral[j-1] >= max_partners_dihedral ||
          n_partners_dihedral[l-1] >= max_partners_dihedral)
        error->all(FLERR,"Dihedral list overflow, boost fix_mscg max");
      partners_dihedral[j-1][n_partners_dihedral[j-1]*3] = i-1;
      partners_dihedral[j-1][n_partners_dihedral[j-1]*3+1] = k-1;
      partners_dihedral[j-1][n_partners_dihedral[j-1]*3+2] = l-1;
      partners_dihedral[l-1][n_partners_dihedral[l-1]*3] = k-1;
      partners_dihedral[l-1][n_partners_dihedral[l-1]*3+1] = i-1;
      partners_dihedral[l-1][n_partners_dihedral[l-1]*3+2] = j-1;
      n_partners_dihedral[j-1]++;
      n_partners_dihedral[l-1]++;
    }
  }

  // pass topology data to mscg code and run startup

  fprintf(screen,"Initializing MSCG with topology data ...\n");
  if (range_flag)
    mscg_struct = rangefinder_startup_part1(mscg_struct);
  else
    mscg_struct = mscg_startup_part1(mscg_struct);

  n_frames = get_n_frames(mscg_struct);
  block_size = get_block_size(mscg_struct);

  mscg_struct =
    setup_topology_and_frame(mscg_struct,n_cg_sites,n_cg_types,type_names,
                             cg_site_types,box_half_lengths);
  mscg_struct =
    set_bond_topology(mscg_struct,partners_bond,n_partners_bond);
  mscg_struct =
    set_angle_topology(mscg_struct,partners_angle,n_partners_angle);
  mscg_struct =
    set_dihedral_topology(mscg_struct,partners_dihedral,n_partners_dihedral);
  mscg_struct =
    generate_exclusion_topology(mscg_struct);

  if (range_flag)
    mscg_struct = rangefinder_startup_part2(mscg_struct);
  else
    mscg_struct = mscg_startup_part2(mscg_struct);
}

/* ---------------------------------------------------------------------- */

void FixMSCG::init()
{
  int nlocal = atom->nlocal;
  double **force = atom->f;
  int i;

  // forces are reset to 0 before pre_force, saved here
  // init called for each frame of dump in rerun command

  for (i = 0; i < nlocal; i++) {
    f[i][0] = force[i][0];
    f[i][1] = force[i][1];
    f[i][2] = force[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixMSCG::end_of_step()
{
  if (domain->triclinic == 1)
    error->all(FLERR,"Fix mscg does not yet support triclinic geometries");

  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  double **x = atom->x;
  int i,ii,j;

  // trajectory information

  for (ii = 0; ii < nlocal; ii++) {
    i = tag[ii];
    for (j = 0; j < 3; j++) {
      x1d[(i-1)*3+j] = x[ii][j];
      f1d[(i-1)*3+j] = f[ii][j];
    }
  }

  // pass x,f to mscg to update matrix

  nframes++;
  if (range_flag)
    mscg_struct = rangefinder_process_frame(mscg_struct,x1d,f1d);
  else
    mscg_struct = mscg_process_frame(mscg_struct,x1d,f1d);
}

/* ---------------------------------------------------------------------- */

void FixMSCG::post_run()
{
  // call mscg to solve matrix and generate output

  fprintf(screen,"Finalizing MSCG ...\n");

  if (nframes != n_frames)
    error->warning(FLERR,"Fix mscg n_frames is inconsistent with control.in");
  if (nframes % block_size != 0)
    error->warning(FLERR,"Fix mscg n_frames is not divisible by "
                   "block_size in control.in");

  if (range_flag)
    rangefinder_solve_and_output(mscg_struct);
  else
    mscg_solve_and_output(mscg_struct);
}

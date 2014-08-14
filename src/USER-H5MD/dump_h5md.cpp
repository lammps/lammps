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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "limits.h"
#include "ch5md.h"
#include "dump_h5md.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "output.h"
#include "error.h"
#include "force.h"
#include "memory.h"

using namespace LAMMPS_NS;

#define MYMIN(a,b) ((a) < (b) ? (a) : (b))
#define MYMAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DumpH5MD::DumpH5MD(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump xtc command");
  if (binary || compressed || multifile || multiproc)
    error->all(FLERR,"Invalid dump h5md filename");

  size_one = 3;
  sort_flag = 1;
  sortcol = 0;
  format_default = NULL;
  flush_flag = 0;
  unwrap_flag = 0;

  // allocate global array for atom coords

  bigint n = group->count(igroup);
  natoms = static_cast<int> (n);

  memory->create(coords,domain->dimension*natoms,"dump:coords");

  openfile();
  ntotal = 0;
}

/* ---------------------------------------------------------------------- */

DumpH5MD::~DumpH5MD()
{
  memory->destroy(coords);

}

/* ---------------------------------------------------------------------- */

void DumpH5MD::init_style()
{
  if (sort_flag == 0 || sortcol != 0)
    error->all(FLERR,"Dump h5md requires sorting by atom ID");
}

/* ---------------------------------------------------------------------- */

void DumpH5MD::openfile()
{
  char *group_name;
  int group_name_length;
  int dims[2];
  char *boundary[3];
  for (int i=0; i<3; i++) {
    boundary[i] = new char[9];
    if (domain->periodicity[i]==1) {
      strcpy(boundary[i], "periodic");
    } else {
      strcpy(boundary[i], "none");
    }
  }

  if (me == 0) {
    datafile = h5md_create_file(filename, "N/A", NULL, "lammps", "N/A");
    group_name_length = strlen(group->names[igroup]);
    group_name = new char[group_name_length];
    strcpy(group_name, group->names[igroup]);
    particles_data = h5md_create_particles_group(datafile, group_name);
    delete [] group_name;
    dims[0] = natoms;
    dims[1] = domain->dimension;
    particles_data.position = h5md_create_time_data(particles_data.group, "position", 2, dims, H5T_NATIVE_DOUBLE, NULL);
    h5md_create_box(&particles_data, dims[1], boundary, true, NULL);
  }
  for (int i=0; i<3; i++) {
    delete [] boundary[i];
  }

}

/* ---------------------------------------------------------------------- */

void DumpH5MD::write_header(bigint nbig)
{
  return;
}

/* ---------------------------------------------------------------------- */

void DumpH5MD::pack(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int dim=domain->dimension;

  m = n = 0;
  if (unwrap_flag == 1) {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        int ix = (image[i] & IMGMASK) - IMGMAX;
        int iy = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        int iz = (image[i] >> IMG2BITS) - IMGMAX;

	buf[m++] = (x[i][0] + ix * xprd);
	buf[m++] = (x[i][1] + iy * yprd);
	if (dim>2) buf[m++] = (x[i][2] + iz * zprd);
        ids[n++] = tag[i];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        buf[m++] = x[i][0];
        buf[m++] = x[i][1];
        if (dim>2) buf[m++] = x[i][2];
        ids[n++] = tag[i];
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpH5MD::write_data(int n, double *mybuf)
{
  // copy buf atom coords into global array

  int m = 0;
  int dim = domain->dimension;
  int k = dim*ntotal;
  for (int i = 0; i < n; i++) {
    for (int j=0; j<dim; j++) {
      coords[k++] = mybuf[m++];
    }
    ntotal++;
  }

  // if last chunk of atoms in this snapshot, write global arrays to file

  if (ntotal == natoms) {
    write_frame();
    ntotal = 0;
  }
}

/* ---------------------------------------------------------------------- */

int DumpH5MD::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"unwrap") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) unwrap_flag = 1;
    else if (strcmp(arg[1],"no") == 0) unwrap_flag = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpH5MD::write_frame()
{
  int local_step;
  double local_time;
  double edges[3];
  local_step = update->ntimestep;
  local_time = local_step * update->dt;
  h5md_append(particles_data.position, coords, local_step, local_time);
  edges[0] = boxxhi - boxxlo;
  edges[1] = boxyhi - boxylo;
  edges[2] = boxzhi - boxzlo;
  h5md_append(particles_data.box_edges, edges, local_step, local_time);
}


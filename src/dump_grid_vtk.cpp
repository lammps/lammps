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

#include "dump_grid_vtk.h"

#include "domain.h"
#include "error.h"
#include "memory.h"
#include "vtk_writer.h"

using namespace LAMMPS_NS;

enum { SCALAR, VECTOR };

/* ---------------------------------------------------------------------- */

DumpGridVTK::DumpGridVTK(LAMMPS *lmp, int narg, char **arg) :
    DumpGrid(lmp, narg, arg), xcoord(nullptr), ycoord(nullptr), zcoord(nullptr)
{
  if (binary || multiproc) error->all(FLERR, "Invalid dump grid/vtk filename: {}", filename);
  if (nfield != 1 && nfield != 3) error->all(FLERR, "Dump grid/vtk requires one or three fields\n");

  buffer_allow = 0;
  buffer_flag = 0;
  sort_flag = 1;
  sortcol = 0;

  mode = (nfield == 1) ? SCALAR : VECTOR;
}

/* ---------------------------------------------------------------------- */

DumpGridVTK::~DumpGridVTK()
{
  memory->destroy(xcoord);
  memory->destroy(ycoord);
  memory->destroy(zcoord);
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::init_style()
{
  DumpGrid::init_style();

  if (multifile == 0) error->all(FLERR, "Dump grid/vtk requires one snapshot per file");
  if (sort_flag == 0 || sortcol > 0) error->all(FLERR, "Dump grid/vtk requires sorting on IDs");
  if (domain->triclinic)
    error->all(FLERR,
               "Dump grid/vtk does not support triclinic simulation boxes, use dump grid instead");
  if (binary) error->all(FLERR, "Dump grid/vtk cannot write binary files\n");

  if (!xcoord) {
    memory->create(xcoord, nxgrid + 1, "dumpgridVTK:xcoord");
    memory->create(ycoord, nygrid + 1, "dumpgridVTK:ycoord");
    memory->create(zcoord, nzgrid + 1, "dumpgridVTK:zcoord");
  }
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_header(bigint /*ndump*/)
{
  if (me) return;

  // the grid data of one snapshot is collected by write_data() and only
  // written out by write_footer(), so all this has to do is prepare for it

  xyz_grid();

  values.clear();
  values.reserve((std::size_t) nxgrid * nygrid * nzgrid * nfield);
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_data(int n, double *mybuf)
{
  const std::size_t nvalues = (std::size_t) n * nfield;
  for (std::size_t m = 0; m < nvalues; m++) values.push_back(mybuf[m]);
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_footer()
{
  if (me) return;

  // grid cell coordinates bound the cells, so there is one more of them
  // than there are cells in each dimension

  const std::vector<double> xc(xcoord, xcoord + nxgrid + 1);
  const std::vector<double> yc(ycoord, ycoord + nygrid + 1);
  const std::vector<double> zc(zcoord, zcoord + nzgrid + 1);

  try {
    VTKWriter writer(VTKWriter::XML, false);
    writer.set_rectilinear_grid(xc, yc, zc);
    writer.add_cell_array((mode == SCALAR) ? "Scalar" : "Vector", nfield, values);
    writer.write(fp);
  } catch (VTKWriterException &e) {
    error->one(FLERR, "Cannot write dump grid/vtk file {}: {}", filename, e.what());
  }
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::xyz_grid()
{
  double *boxlo = domain->boxlo;

  double dx = domain->prd[0] / nxgrid;
  double dy = domain->prd[1] / nygrid;
  double dz = domain->prd[2] / nzgrid;

  for (int ix = 0; ix <= nxgrid; ix++) xcoord[ix] = boxlo[0] + ix * dx;
  for (int iy = 0; iy <= nygrid; iy++) ycoord[iy] = boxlo[1] + iy * dy;
  for (int iz = 0; iz <= nzgrid; iz++) zcoord[iz] = boxlo[2] + iz * dz;
}

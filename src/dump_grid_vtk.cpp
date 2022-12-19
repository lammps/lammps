/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

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

  xyz_grid();

  fprintf(fp, "<?xml version=\"1.0\"\?>\n");
  fprintf(fp, "<VTKFile type=\"RectilinearGrid\">\n");
  fprintf(fp,
          "<RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\" "
          "Origin=\"0 0 0\" Spacing=\"1 1 1\">\n",
          nxgrid, nygrid, nzgrid);
  fprintf(fp, "<Piece Extent=\"0 %d 0 %d 0 %d\">\n", nxgrid, nygrid, nzgrid);
  fprintf(fp, "<PointData>\n");
  fprintf(fp, "</PointData>\n");
  fprintf(fp, "<Coordinates>\n");

  // coords of center point of grid cells in each of xyz dimensions

  fprintf(fp, "<DataArray type=\"Float32\" format=\"ascii\">\n");
  for (int i = 0; i <= nxgrid; i++) fprintf(fp, "%g ", xcoord[i]);
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Float32\" format=\"ascii\">\n");
  for (int i = 0; i <= nygrid; i++) fprintf(fp, "%g ", ycoord[i]);
  fprintf(fp, "\n</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Float32\" format=\"ascii\">\n");
  for (int i = 0; i <= nzgrid; i++) fprintf(fp, "%g ", zcoord[i]);
  fprintf(fp, "\n</DataArray>\n");

  fprintf(fp, "</Coordinates>\n");

  fprintf(fp, "<CellData>\n");
  if (mode == SCALAR)
    fprintf(fp, "<DataArray type=\"Float32\" Name=\"Scalar\" format=\"ascii\">\n");
  else if (mode == VECTOR)
    fprintf(fp,
            "<DataArray type=\"Float32\" Name=\"Vector\" "
            "NumberOfComponents=\"3\" format=\"ascii\">\n");
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_footer()
{
  if (me) return;

  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</RectilinearGrid>\n");
  fprintf(fp, "</VTKFile>\n");
}

/* ---------------------------------------------------------------------- */

void DumpGridVTK::write_data(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < nfield; j++) {
      if (vtype[j] == Dump::INT)
        fprintf(fp, vformat[j], static_cast<int>(mybuf[m]));
      else if (vtype[j] == Dump::DOUBLE)
        fprintf(fp, vformat[j], mybuf[m]);
      else if (vtype[j] == Dump::BIGINT)
        fprintf(fp, vformat[j], static_cast<bigint>(mybuf[m]));
      m++;
    }
    fprintf(fp, "\n");
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

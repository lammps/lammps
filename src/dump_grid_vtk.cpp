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
#include "force.h"
#include "memory.h"

#include <cstring>
#include <string>
#include <utility>

using namespace LAMMPS_NS;

enum { SCALAR, VECTOR };
enum { VTKLEGACY, VTKXML };
enum { RECTILINEAR, IMAGE };

/* ---------------------------------------------------------------------- */

DumpGridVTK::DumpGridVTK(LAMMPS *lmp, int narg, char **arg) :
    DumpGrid(lmp, narg, arg), xcoord(nullptr), ycoord(nullptr), zcoord(nullptr)
{
  if (multiproc) error->all(FLERR, "Invalid dump grid/vtk filename: {}", filename);
  if (nfield != 1 && nfield != 3) error->all(FLERR, "Dump grid/vtk requires one or three fields\n");

  buffer_allow = 0;
  buffer_flag = 0;
  sort_flag = 1;
  sortcol = 0;

  mode = (nfield == 1) ? SCALAR : VECTOR;

  // the file name extension selects which of the file formats to write.
  // an unrecognized extension keeps the XML rectilinear grid default.

  vtkflavor = VTKXML;
  dataset = RECTILINEAR;
  precision_warned = 0;
  writeprec = VTKWriter::SINGLE;

  std::string fname(filename);
  auto dot = fname.find_last_of('.');
  if (dot != std::string::npos) {
    const std::string ext = fname.substr(dot);
    if (ext == ".vtk") {
      vtkflavor = VTKLEGACY;
      dataset = RECTILINEAR;
    } else if (ext == ".vti") {
      vtkflavor = VTKXML;
      dataset = IMAGE;
    }
  }
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

  // the image data format describes the grid by origin and spacing alone,
  // so the cell boundary coordinates are only needed for rectilinear grids

  if ((dataset == RECTILINEAR) && !xcoord) {
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

  if (dataset == RECTILINEAR) xyz_grid();

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

  try {
    VTKWriter writer((vtkflavor == VTKXML) ? VTKWriter::XML : VTKWriter::LEGACY, binary != 0,
                     writeprec);

    if (dataset == IMAGE) {

      // the grid cells all have the same size, so the grid can be described
      // by its origin and spacing alone

      const int dims[3] = {nxgrid + 1, nygrid + 1, nzgrid + 1};
      const double origin[3] = {domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};
      const double spacing[3] = {domain->prd[0] / nxgrid, domain->prd[1] / nygrid,
                                 domain->prd[2] / nzgrid};
      writer.set_image_data(dims, origin, spacing);

    } else {

      // grid cell coordinates bound the cells, so there is one more of them
      // than there are cells in each dimension

      writer.set_rectilinear_grid({xcoord, xcoord + nxgrid + 1}, {ycoord, ycoord + nygrid + 1},
                                  {zcoord, zcoord + nzgrid + 1});
    }

    // the collected data is donated to the writer, it is rebuilt for every
    // snapshot anyway

    writer.add_cell_array((mode == SCALAR) ? "Scalar" : "Vector", nfield, std::move(values));
    writer.write(fp);

    // grid coordinates are stored in single precision by default, warn once
    // if that no longer resolves them well enough.  with dump_modify double
    // yes nothing is tracked, so the warning is skipped automatically.

    const double maxcoord = writer.max_single_precision_value();
    const double resolution = VTKWriter::single_precision_resolution(maxcoord, force->angstrom);
    if (!precision_warned && (resolution > 0.0)) {
      precision_warned = 1;
      error->warning(FLERR,
                     "Dump grid/vtk writes grid coordinates in single precision, which "
                     "resolves the largest coordinate of this grid, {:.4g}, to only about {:.2g} "
                     "length units. Use dump_modify double yes if your analysis needs more "
                     "resolution than that.",
                     maxcoord, resolution);
    }
  } catch (VTKWriterException &e) {
    error->one(FLERR, "Cannot write dump grid/vtk file {}: {}", filename, e.what());
  }
}

/* ----------------------------------------------------------------------
   the VTK file name extensions collide with the LAMMPS convention of
   appending ".bin" to select binary output, so offer a keyword instead
------------------------------------------------------------------------- */

int DumpGridVTK::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "binary") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify binary", error);
    binary = utils::logical(FLERR, arg[1], false, lmp);
    return 2;
  }
  if (strcmp(arg[0], "double") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify double", error);
    writeprec = utils::logical(FLERR, arg[1], false, lmp) ? VTKWriter::DOUBLE : VTKWriter::SINGLE;
    return 2;
  }
  return 0;
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

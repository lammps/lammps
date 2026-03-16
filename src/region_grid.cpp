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

#include "region_grid.h"

#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "grid3d.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr int OFFSET = 16384;

/* ---------------------------------------------------------------------- */

RegGrid::RegGrid(LAMMPS *lmp, int narg, char **arg) :
    Region(lmp, narg, arg), gridref(nullptr), source_id(nullptr),
    grid_name(nullptr), data_name(nullptr), grid3d(nullptr), griddata(nullptr)
{
  // region ID grid gridref op value [region-options]

  if (narg < 5) utils::missing_cmd_args(FLERR, "region grid", error);

  varshape = 1;

  gridref = utils::strdup(arg[2]);

  if (strcmp(arg[3], "gt") == 0) compare_op = OP_GT;
  else if (strcmp(arg[3], "ge") == 0) compare_op = OP_GE;
  else if (strcmp(arg[3], "lt") == 0) compare_op = OP_LT;
  else if (strcmp(arg[3], "le") == 0) compare_op = OP_LE;
  else if (strcmp(arg[3], "eq") == 0) compare_op = OP_EQ;
  else if (strcmp(arg[3], "ne") == 0) compare_op = OP_NE;
  else error->all(FLERR, "Invalid comparison operator '{}' for region grid", arg[3]);

  threshold = utils::numeric(FLERR, arg[4], false, lmp);

  // parse gridref to extract source type, ID, grid name, data name, and column index
  // format: c_ID:gname:dname  or  f_ID:gname:dname  (optionally with [N])

  std::string ref(gridref);
  if (ref.size() < 3 || ref[1] != '_' || (ref[0] != 'c' && ref[0] != 'f'))
    error->all(FLERR, "Invalid grid reference '{}' for region grid", gridref);

  source_type = (ref[0] == 'c') ? COMPUTE_SOURCE : FIX_SOURCE;
  std::string rest = ref.substr(2);

  gridindex = 0;
  auto bracket_pos = rest.find('[');
  if (bracket_pos != std::string::npos) {
    auto end_pos = rest.find(']');
    if (end_pos == std::string::npos)
      error->all(FLERR, "Missing ']' in grid reference '{}'", gridref);
    gridindex = utils::inumeric(FLERR,
        rest.substr(bracket_pos + 1, end_pos - bracket_pos - 1), false, lmp);
    rest = rest.substr(0, bracket_pos);
  }

  auto pos1 = rest.find(':');
  if (pos1 == std::string::npos)
    error->all(FLERR, "Invalid grid reference '{}' - expected ID:gname:dname", gridref);
  auto pos2 = rest.find(':', pos1 + 1);
  if (pos2 == std::string::npos)
    error->all(FLERR, "Invalid grid reference '{}' - expected ID:gname:dname", gridref);

  source_id = utils::strdup(rest.substr(0, pos1));
  grid_name = utils::strdup(rest.substr(pos1 + 1, pos2 - pos1 - 1));
  data_name = utils::strdup(rest.substr(pos2 + 1));

  cmax = 6;
  contact = new Contact[cmax];
  tmax = 0;

  options(narg - 5, &arg[5]);
}

/* ---------------------------------------------------------------------- */

RegGrid::~RegGrid()
{
  delete[] gridref;
  delete[] source_id;
  delete[] grid_name;
  delete[] data_name;
}

/* ---------------------------------------------------------------------- */

void RegGrid::init()
{
  Region::init();
  resolve_grid_reference();
}

/* ---------------------------------------------------------------------- */

void RegGrid::resolve_grid_reference()
{
  int dim;

  if (source_type == COMPUTE_SOURCE) {
    auto *compute = modify->get_compute_by_id(source_id);
    if (!compute)
      error->all(FLERR, "Compute '{}' for region grid does not exist", source_id);

    igrid = compute->get_grid_by_name(grid_name, dim);
    if (igrid < 0)
      error->all(FLERR, "Compute '{}' does not provide grid '{}'", source_id, grid_name);
    if (dim != 3)
      error->all(FLERR, "Region grid requires a 3d grid from compute '{}'", source_id);

    idata = compute->get_griddata_by_name(igrid, data_name, ncol);
    if (idata < 0)
      error->all(FLERR, "Compute '{}' does not provide grid data '{}'", source_id, data_name);

    grid3d = (Grid3d *) compute->get_grid_by_index(igrid);
    griddata = compute->get_griddata_by_index(idata);

  } else {
    auto *fix = modify->get_fix_by_id(source_id);
    if (!fix)
      error->all(FLERR, "Fix '{}' for region grid does not exist", source_id);

    igrid = fix->get_grid_by_name(grid_name, dim);
    if (igrid < 0)
      error->all(FLERR, "Fix '{}' does not provide grid '{}'", source_id, grid_name);
    if (dim != 3)
      error->all(FLERR, "Region grid requires a 3d grid from fix '{}'", source_id);

    idata = fix->get_griddata_by_name(igrid, data_name, ncol);
    if (idata < 0)
      error->all(FLERR, "Fix '{}' does not provide grid data '{}'", source_id, data_name);

    grid3d = (Grid3d *) fix->get_grid_by_index(igrid);
    griddata = fix->get_griddata_by_index(idata);
  }

  if (!grid3d)
    error->all(FLERR, "Region grid could not resolve grid reference '{}'", gridref);

  grid3d->get_size(nx, ny, nz);
  grid3d->get_bounds_ghost(nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

  update_bbox();
}

/* ---------------------------------------------------------------------- */

void RegGrid::shape_update()
{
  if (source_type == COMPUTE_SOURCE) {
    auto *compute = modify->get_compute_by_id(source_id);
    if (!compute) return;

    if (!(compute->invoked_flag & Compute::INVOKED_PERGRID)) {
      compute->compute_pergrid();
      compute->invoked_flag |= Compute::INVOKED_PERGRID;
    }

    grid3d = (Grid3d *) compute->get_grid_by_index(igrid);
    griddata = compute->get_griddata_by_index(idata);

  } else {
    auto *fix = modify->get_fix_by_id(source_id);
    if (!fix) return;

    grid3d = (Grid3d *) fix->get_grid_by_index(igrid);
    griddata = fix->get_griddata_by_index(idata);
  }

  if (grid3d) {
    grid3d->get_size(nx, ny, nz);
    grid3d->get_bounds_ghost(nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);
    update_bbox();
  }
}

/* ---------------------------------------------------------------------- */

void RegGrid::update_bbox()
{
  if (!grid3d || !griddata) {
    bboxflag = 0;
    return;
  }

  bboxflag = 1;
  extent_xlo = domain->boxlo[0];
  extent_xhi = domain->boxhi[0];
  extent_ylo = domain->boxlo[1];
  extent_yhi = domain->boxhi[1];
  extent_zlo = domain->boxlo[2];
  extent_zhi = domain->boxhi[2];
}

/* ---------------------------------------------------------------------- */

int RegGrid::inside(double x, double y, double z)
{
  if (!grid3d || !griddata) return 0;

  double *boxlo = domain->boxlo;
  double dxinv = nx / domain->xprd;
  double dyinv = ny / domain->yprd;
  double dzinv = nz / domain->zprd;

  int ix = static_cast<int>((x - boxlo[0]) * dxinv + OFFSET) - OFFSET;
  int iy = static_cast<int>((y - boxlo[1]) * dyinv + OFFSET) - OFFSET;
  int iz = static_cast<int>((z - boxlo[2]) * dzinv + OFFSET) - OFFSET;

  if (ix < nxlo_out || ix > nxhi_out ||
      iy < nylo_out || iy > nyhi_out ||
      iz < nzlo_out || iz > nzhi_out) return 0;

  double value;
  if (ncol == 0) {
    auto vec3d = (double ***) griddata;
    value = vec3d[iz][iy][ix];
  } else {
    auto array3d = (double ****) griddata;
    int col = (gridindex > 0) ? gridindex - 1 : 0;
    value = array3d[iz][iy][ix][col];
  }

  return evaluate(value);
}

/* ---------------------------------------------------------------------- */

int RegGrid::surface_interior(double *x, double cutoff)
{
  if (!grid3d || !griddata) return 0;

  double *boxlo = domain->boxlo;
  double dx = domain->xprd / nx;
  double dy = domain->yprd / ny;
  double dz = domain->zprd / nz;

  int ix = static_cast<int>((x[0] - boxlo[0]) / dx + OFFSET) - OFFSET;
  int iy = static_cast<int>((x[1] - boxlo[1]) / dy + OFFSET) - OFFSET;
  int iz = static_cast<int>((x[2] - boxlo[2]) / dz + OFFSET) - OFFSET;

  if (!cell_inside(ix, iy, iz)) return 0;

  int n = 0;
  double delta;

  // -x face
  delta = x[0] - (boxlo[0] + ix * dx);
  if (delta < cutoff && !cell_inside(ix - 1, iy, iz)) {
    contact[n].r = delta;
    contact[n].delx = delta;
    contact[n].dely = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 0;
    n++;
  }

  // +x face
  delta = (boxlo[0] + (ix + 1) * dx) - x[0];
  if (delta < cutoff && !cell_inside(ix + 1, iy, iz)) {
    contact[n].r = delta;
    contact[n].delx = -delta;
    contact[n].dely = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 1;
    n++;
  }

  // -y face
  delta = x[1] - (boxlo[1] + iy * dy);
  if (delta < cutoff && !cell_inside(ix, iy - 1, iz)) {
    contact[n].r = delta;
    contact[n].dely = delta;
    contact[n].delx = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 2;
    n++;
  }

  // +y face
  delta = (boxlo[1] + (iy + 1) * dy) - x[1];
  if (delta < cutoff && !cell_inside(ix, iy + 1, iz)) {
    contact[n].r = delta;
    contact[n].dely = -delta;
    contact[n].delx = contact[n].delz = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 3;
    n++;
  }

  // -z face
  delta = x[2] - (boxlo[2] + iz * dz);
  if (delta < cutoff && !cell_inside(ix, iy, iz - 1)) {
    contact[n].r = delta;
    contact[n].delz = delta;
    contact[n].delx = contact[n].dely = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 4;
    n++;
  }

  // +z face
  delta = (boxlo[2] + (iz + 1) * dz) - x[2];
  if (delta < cutoff && !cell_inside(ix, iy, iz + 1)) {
    contact[n].r = delta;
    contact[n].delz = -delta;
    contact[n].delx = contact[n].dely = 0.0;
    contact[n].radius = 0;
    contact[n].iwall = 5;
    n++;
  }

  return n;
}

/* ---------------------------------------------------------------------- */

int RegGrid::surface_exterior(double *x, double cutoff)
{
  if (!grid3d || !griddata) return 0;

  double *boxlo = domain->boxlo;
  double dx = domain->xprd / nx;
  double dy = domain->yprd / ny;
  double dz = domain->zprd / nz;

  int ix = static_cast<int>((x[0] - boxlo[0]) / dx + OFFSET) - OFFSET;
  int iy = static_cast<int>((x[1] - boxlo[1]) / dy + OFFSET) - OFFSET;
  int iz = static_cast<int>((x[2] - boxlo[2]) / dz + OFFSET) - OFFSET;

  if (ix < nxlo_out || ix > nxhi_out ||
      iy < nylo_out || iy > nyhi_out ||
      iz < nzlo_out || iz > nzhi_out) return 0;

  if (cell_inside(ix, iy, iz)) return 0;

  double mindist = cutoff;
  double best_delx = 0.0, best_dely = 0.0, best_delz = 0.0;
  int found = 0;
  double delta;

  // -x face: neighbor at ix-1 is inside
  delta = x[0] - (boxlo[0] + ix * dx);
  if (delta < mindist && cell_inside(ix - 1, iy, iz)) {
    mindist = delta;
    best_delx = delta;
    best_dely = best_delz = 0.0;
    found = 1;
  }

  // +x face: neighbor at ix+1 is inside
  delta = (boxlo[0] + (ix + 1) * dx) - x[0];
  if (delta < mindist && cell_inside(ix + 1, iy, iz)) {
    mindist = delta;
    best_delx = -delta;
    best_dely = best_delz = 0.0;
    found = 1;
  }

  // -y face
  delta = x[1] - (boxlo[1] + iy * dy);
  if (delta < mindist && cell_inside(ix, iy - 1, iz)) {
    mindist = delta;
    best_dely = delta;
    best_delx = best_delz = 0.0;
    found = 1;
  }

  // +y face
  delta = (boxlo[1] + (iy + 1) * dy) - x[1];
  if (delta < mindist && cell_inside(ix, iy + 1, iz)) {
    mindist = delta;
    best_dely = -delta;
    best_delx = best_delz = 0.0;
    found = 1;
  }

  // -z face
  delta = x[2] - (boxlo[2] + iz * dz);
  if (delta < mindist && cell_inside(ix, iy, iz - 1)) {
    mindist = delta;
    best_delz = delta;
    best_delx = best_dely = 0.0;
    found = 1;
  }

  // +z face
  delta = (boxlo[2] + (iz + 1) * dz) - x[2];
  if (delta < mindist && cell_inside(ix, iy, iz + 1)) {
    mindist = delta;
    best_delz = -delta;
    best_delx = best_dely = 0.0;
    found = 1;
  }

  if (!found) return 0;

  contact[0].r = mindist;
  contact[0].delx = best_delx;
  contact[0].dely = best_dely;
  contact[0].delz = best_delz;
  contact[0].radius = 0;
  contact[0].iwall = 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

bool RegGrid::cell_inside(int ix, int iy, int iz)
{
  if (ix < nxlo_out || ix > nxhi_out ||
      iy < nylo_out || iy > nyhi_out ||
      iz < nzlo_out || iz > nzhi_out) return false;

  double value;
  if (ncol == 0) {
    auto vec3d = (double ***) griddata;
    value = vec3d[iz][iy][ix];
  } else {
    auto array3d = (double ****) griddata;
    int col = (gridindex > 0) ? gridindex - 1 : 0;
    value = array3d[iz][iy][ix][col];
  }

  return evaluate(value) == 1;
}

/* ---------------------------------------------------------------------- */

int RegGrid::evaluate(double value)
{
  switch (compare_op) {
    case OP_GT: return (value > threshold) ? 1 : 0;
    case OP_GE: return (value >= threshold) ? 1 : 0;
    case OP_LT: return (value < threshold) ? 1 : 0;
    case OP_LE: return (value <= threshold) ? 1 : 0;
    case OP_EQ: return (value == threshold) ? 1 : 0;
    case OP_NE: return (value != threshold) ? 1 : 0;
  }
  return 0;
}

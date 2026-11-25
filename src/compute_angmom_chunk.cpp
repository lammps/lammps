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

#include "compute_angmom_chunk.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

static constexpr double SINERTIA = 0.4;          // moment of inertia prefactor for sphere
static constexpr double LINERTIA = 1.0 / 12.0;   // moment of inertia prefactor for line

/* ---------------------------------------------------------------------- */

ComputeAngmomChunk::ComputeAngmomChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), com(nullptr),
    comall(nullptr), angmom(nullptr), angmomall(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute angmom/chunk command");

  array_flag = 1;
  size_array_cols = 3;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeAngmomChunk::init();
  ComputeAngmomChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeAngmomChunk::~ComputeAngmomChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(angmom);
  memory->destroy(angmomall);
}

/* ---------------------------------------------------------------------- */

void ComputeAngmomChunk::compute_array()
{
  ComputeChunk::compute_array();

  int i, index;
  double dx, dy, dz, massone;
  double unwrap[3];
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    angmom[i][0] = angmom[i][1] = angmom[i][2] = 0.0;
  }

  // compute COM for each chunk

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      domain->unmap(x[i], image[i], unwrap);
      massproc[index] += massone;
      com[index][0] += unwrap[0] * massone;
      com[index][1] += unwrap[1] * massone;
      com[index][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(massproc, masstotal, nchunk, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&com[0][0], &comall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);

  for (i = 0; i < nchunk; i++) {
    if (masstotal[i] > 0.0) {
      comall[i][0] /= masstotal[i];
      comall[i][1] /= masstotal[i];
      comall[i][2] /= masstotal[i];
    }
  }

  // compute angmom for each chunk
  // finite-size particles add their own (spin) angular momentum:
  //   ANGMOM-type (ellipsoid, superellipsoid, triangle, body) store it
  //   directly; OMEGA-type (sphere, line) carry an angular velocity, so
  //   L_spin = I_spin * omega. shape indices are checked before the sphere
  //   (radius) case, since finite-shape particles also carry a radius

  double **v = atom->v;
  double vunwrap[3];

  auto *avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto *avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto *avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto *avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  double *radius = atom->radius;
  double **omega = atom->omega;
  double **angmom_one = atom->angmom;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *body = atom->body;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      domain->unmap(x[i], v[i], image[i], mask[i], unwrap, vunwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      angmom[index][0] += massone * (dy * vunwrap[2] - dz * vunwrap[1]);
      angmom[index][1] += massone * (dz * vunwrap[0] - dx * vunwrap[2]);
      angmom[index][2] += massone * (dx * vunwrap[1] - dy * vunwrap[0]);

      if (avec_ellipsoid && ellipsoid[i] >= 0) {
        if (angmom_one) {
          angmom[index][0] += angmom_one[i][0];
          angmom[index][1] += angmom_one[i][1];
          angmom[index][2] += angmom_one[i][2];
        }
      } else if (avec_tri && tri[i] >= 0) {
        if (angmom_one) {
          angmom[index][0] += angmom_one[i][0];
          angmom[index][1] += angmom_one[i][1];
          angmom[index][2] += angmom_one[i][2];
        }
      } else if (avec_body && body[i] >= 0) {
        if (angmom_one) {
          angmom[index][0] += angmom_one[i][0];
          angmom[index][1] += angmom_one[i][1];
          angmom[index][2] += angmom_one[i][2];
        }
      } else if (avec_line && line[i] >= 0) {
        if (omega) {
          double length = avec_line->bonus[line[i]].length;
          angmom[index][2] += LINERTIA * massone * length * length * omega[i][2];
        }
      } else if (radius && radius[i] > 0.0) {
        if (omega) {
          double sphere = SINERTIA * massone * radius[i] * radius[i];
          angmom[index][0] += sphere * omega[i][0];
          angmom[index][1] += sphere * omega[i][1];
          angmom[index][2] += sphere * omega[i][2];
        }
      }
    }

  MPI_Allreduce(&angmom[0][0], &angmomall[0][0], 3 * nchunk, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeAngmomChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(angmom);
  memory->destroy(angmomall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "angmom/chunk:massproc");
  memory->create(masstotal, maxchunk, "angmom/chunk:masstotal");
  memory->create(com, maxchunk, 3, "angmom/chunk:com");
  memory->create(comall, maxchunk, 3, "angmom/chunk:comall");
  memory->create(angmom, maxchunk, 3, "angmom/chunk:angmom");
  memory->create(angmomall, maxchunk, 3, "angmom/chunk:angmomall");
  array = angmomall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeAngmomChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  return bytes;
}

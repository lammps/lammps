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

#include "compute_inertia_chunk.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "memory.h"

using namespace LAMMPS_NS;

static constexpr double SINERTIA = 0.4;    // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

ComputeInertiaChunk::ComputeInertiaChunk(LAMMPS *lmp, int narg, char **arg) :
    ComputeChunk(lmp, narg, arg), massproc(nullptr), masstotal(nullptr), com(nullptr),
    comall(nullptr), inertia(nullptr), inertiaall(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute inertia/chunk command");

  array_flag = 1;
  size_array_cols = 6;
  size_array_rows = 0;
  size_array_rows_variable = 1;
  extarray = 0;

  ComputeInertiaChunk::init();
  ComputeInertiaChunk::allocate();
}

/* ---------------------------------------------------------------------- */

ComputeInertiaChunk::~ComputeInertiaChunk()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(inertia);
  memory->destroy(inertiaall);
}

/* ---------------------------------------------------------------------- */

void ComputeInertiaChunk::compute_array()
{
  int i, j, index;
  double dx, dy, dz, massone;
  double unwrap[3];

  ComputeChunk::compute_array();
  int *ichunk = cchunk->ichunk;

  // zero local per-chunk values

  for (i = 0; i < nchunk; i++) {
    massproc[i] = 0.0;
    com[i][0] = com[i][1] = com[i][2] = 0.0;
    for (j = 0; j < 6; j++) inertia[i][j] = 0.0;
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

  // compute inertia tensor for each chunk
  // finite-size particles add their own moment of inertia about their center
  //   (parallel-axis theorem); MathExtra::inertia_* return the tensor in
  //   (Ixx,Iyy,Izz,Iyz,Ixz,Ixy) order, remapped here to the (Ixx,Iyy,Izz,
  //   Ixy,Iyz,Ixz) column order of this compute

  auto *avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  auto *avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  auto *avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  auto *avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  double *radius = atom->radius;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *body = atom->body;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      index = ichunk[i] - 1;
      if (index < 0) continue;
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - comall[index][0];
      dy = unwrap[1] - comall[index][1];
      dz = unwrap[2] - comall[index][2];
      inertia[index][0] += massone * (dy * dy + dz * dz);
      inertia[index][1] += massone * (dx * dx + dz * dz);
      inertia[index][2] += massone * (dx * dx + dy * dy);
      inertia[index][3] -= massone * dx * dy;
      inertia[index][4] -= massone * dy * dz;
      inertia[index][5] -= massone * dx * dz;

      // finite-size particle spin inertia as a Voigt 6-vector
      // (Ixx,Iyy,Izz,Iyz,Ixz,Ixy); zero for point particles
      // check the shape-defining bonus indices before the sphere (radius)
      // case: ellipsoid/superellipsoid/triangle/body particles also carry a
      // (bounding) radius, so the finite-sphere branch must only catch true
      // spheres

      double ivec[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      if (avec_ellipsoid && ellipsoid[i] >= 0) {
        if (atom->superellipsoid_flag)
          MathExtra::inertia_ellipsoid(avec_ellipsoid->bonus_super[ellipsoid[i]].inertia,
                                       avec_ellipsoid->bonus_super[ellipsoid[i]].quat, ivec);
        else
          MathExtra::inertia_ellipsoid(avec_ellipsoid->bonus[ellipsoid[i]].shape,
                                       avec_ellipsoid->bonus[ellipsoid[i]].quat, massone, ivec);
      } else if (avec_line && line[i] >= 0) {
        MathExtra::inertia_line(avec_line->bonus[line[i]].length, avec_line->bonus[line[i]].theta,
                                massone, ivec);
      } else if (avec_tri && tri[i] >= 0) {
        MathExtra::inertia_triangle(avec_tri->bonus[tri[i]].inertia, avec_tri->bonus[tri[i]].quat,
                                    massone, ivec);
      } else if (avec_body && body[i] >= 0) {
        MathExtra::inertia_ellipsoid(avec_body->bonus[body[i]].inertia,
                                     avec_body->bonus[body[i]].quat, ivec);
      } else if (radius && radius[i] > 0.0) {
        ivec[0] = ivec[1] = ivec[2] = SINERTIA * massone * radius[i] * radius[i];
      }
      inertia[index][0] += ivec[0];
      inertia[index][1] += ivec[1];
      inertia[index][2] += ivec[2];
      inertia[index][3] += ivec[5];
      inertia[index][4] += ivec[3];
      inertia[index][5] += ivec[4];
    }

  MPI_Allreduce(&inertia[0][0], &inertiaall[0][0], 6 * nchunk, MPI_DOUBLE, MPI_SUM, world);
}
/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeInertiaChunk::allocate()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(inertia);
  memory->destroy(inertiaall);
  maxchunk = nchunk;
  memory->create(massproc, maxchunk, "inertia/chunk:massproc");
  memory->create(masstotal, maxchunk, "inertia/chunk:masstotal");
  memory->create(com, maxchunk, 3, "inertia/chunk:com");
  memory->create(comall, maxchunk, 3, "inertia/chunk:comall");
  memory->create(inertia, maxchunk, 6, "inertia/chunk:inertia");
  memory->create(inertiaall, maxchunk, 6, "inertia/chunk:inertiaall");
  array = inertiaall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeInertiaChunk::memory_usage()
{
  double bytes = ComputeChunk::memory_usage();
  bytes += (bigint) maxchunk * 2 * sizeof(double);
  bytes += (double) maxchunk * 2 * 3 * sizeof(double);
  bytes += (double) maxchunk * 2 * 6 * sizeof(double);
  return bytes;
}

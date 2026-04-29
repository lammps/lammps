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

#include "fix_graphics_chunk.h"

#include "atom.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "graphics.h"
#include "image_objects.h"
#include "lattice.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "region.h"
#include "update.h"

#include <algorithm>
#include <cstring>
#include <unordered_map>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {

using namespace ImageObjects;

// Define vertices of an icosahedron to approximate a sphere of radius 1 around the origin.
// A and B are the normalized coordinates derived from the golden ratio:
// phi = (1 + sqrt(5)) / 2; A = 1 / sqrt(1 + phi^2); B = phi / sqrt(1 + phi^2)
constexpr double A = 0.5257311121191336;
constexpr double B = 0.8506508083520399;
const std::vector<vec3> icosahedron = {
    {-A, B, 0.0},  {A, B, 0.0},  {-A, -B, 0.0}, {A, -B, 0.0}, {0.0, -A, B},  {0.0, A, B},
    {0.0, -A, -B}, {0.0, A, -B}, {B, 0.0, -A},  {B, 0.0, A},  {-B, 0.0, -A}, {-B, 0.0, A},
};
constexpr int NUM_POINTS = 12;

// helpers to split clusters that cross PBC boundaries

struct DSU {
  std::vector<int> parent, rank;

  DSU(int n) : parent(n), rank(n, 0)
  {
    for (int i = 0; i < n; ++i) parent[i] = i;
  }

  int find(int x)
  {
    if (parent[x] != x) parent[x] = find(parent[x]);
    return parent[x];
  }

  void unite(int a, int b)
  {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (rank[a] < rank[b]) std::swap(a, b);
    parent[b] = a;
    if (rank[a] == rank[b]) ++rank[a];
  }
};

struct CellKey {
  bigint x, y, z;
  bool operator==(const CellKey &other) const
  {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct CellHash {
  size_t operator()(const CellKey &k) const noexcept
  {
    size_t h1 = std::hash<bigint>{}(k.x);
    size_t h2 = std::hash<bigint>{}(k.y);
    size_t h3 = std::hash<bigint>{}(k.z);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

double dist2(const double *const a, const double *const b)
{
  double dx = a[0] - b[0];
  double dy = a[1] - b[1];
  double dz = a[2] - b[2];
  return dx * dx + dy * dy + dz * dz;
}

CellKey cellOf(const double *const p, double cellSize)
{
  return {static_cast<bigint>(std::floor(p[0] / cellSize)),
          static_cast<bigint>(std::floor(p[1] / cellSize)),
          static_cast<bigint>(std::floor(p[2] / cellSize))};
}

std::vector<std::vector<int>> splitByDistance(const std::vector<int> &pts,
                                              const double *const *const x, double dist)
{
  const int n = static_cast<int>(pts.size());
  if (n == 0) return {};

  DSU dsu(n);
  const double d2 = dist * dist;

  std::unordered_map<CellKey, std::vector<int>, CellHash> grid;
  grid.reserve(n * 2);

  const double cellSize = dist;

  for (int i = 0; i < n; ++i) {
    CellKey c = cellOf(x[pts[i]], cellSize);

    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          CellKey nc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(nc);
          if (it == grid.end()) continue;

          for (int j : it->second) {
            if (dist2(x[pts[i]], x[pts[j]]) < d2) dsu.unite(i, j);
          }
        }
      }
    }

    grid[c].push_back(i);
  }

  std::unordered_map<int, std::vector<int>> groups;
  groups.reserve(n);

  for (int i = 0; i < n; ++i) groups[dsu.find(i)].push_back(pts[i]);

  std::vector<std::vector<int>> result;
  result.reserve(groups.size());
  for (auto &kv : groups) { result.push_back(std::move(kv.second)); }
  return result;
}
}    // namespace

/* ---------------------------------------------------------------------- */

FixGraphicsChunk::FixGraphicsChunk(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), chunkid(nullptr), atomrad(nullptr), id_chunk(nullptr), cchunk(nullptr),
    id_region(nullptr), region(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix graphics/chunk", error);

  maxchunk = -1;
  comm_forward = 2;

  // parse mandatory args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Fix graphics/chunk nevery value {} must be > 0", nevery);
  global_freq = nevery;
  dynamic_group_allow = 1;

  id_chunk = utils::strdup(arg[4]);
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(id_chunk));
  if (!cchunk)
    error->all(FLERR, 4, "Chunk/atom compute {} does not exist or is incorrect style for fix {}",
               id_chunk, style);

  // defaults
  numobjs = 0;
  alpha = 0.0;
  maxreplace = 100;
  mindist = 0.0;

  // fallback guess for default atom size
  radius = 0.5 * domain->lattice->xlattice;
  has_global_radius = false;
  smooth = true;
  clip = true;

  // parse optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "radius") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk radius", error);
      radius = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (radius <= 0.0) error->all(FLERR, iarg + 1, "Fix graphics/chunk radius value must be > 0");
      has_global_radius = true;
      iarg += 2;
    } else if (strcmp(arg[iarg], "alpha") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk alpha", error);
      alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (alpha < 0.0) error->all(FLERR, iarg + 1, "Fix graphics/chunk alpha value must be >= 0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk region", error);
      delete[] id_region;
      id_region = utils::strdup(arg[iarg + 1]);
      region = domain->get_region_by_id(id_region);
      if (!region)
        error->all(FLERR, iarg + 1, "Region {} for fix graphics/chunk does not exist", id_region);
      iarg += 2;
    } else if (strcmp(arg[iarg], "clip") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk clip", error);
      clip = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "maxreplace") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk maxreplace", error);
      maxreplace = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (maxreplace < 0)
        error->all(FLERR, iarg + 1, "Fix graphics/chunk maxreplace value must be >= 0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "mindist") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk mindist", error);
      mindist = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (mindist <= 0.0)
        error->all(FLERR, iarg + 1, "Fix graphics/chunk mindist value must be > 0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "shading") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk shading", error);
      if (strcmp(arg[iarg + 1], "smooth") == 0) {
        smooth = true;
      } else if (strcmp(arg[iarg + 1], "flat") == 0) {
        smooth = false;
      } else {
        error->all(FLERR, iarg + 1, "Unknown fix graphics/chunk shading setting {}", arg[iarg + 1]);
      }
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/chunk keyword {}", arg[iarg]);
    }
  }

  if (domain->dimension == 2)
    error->all(FLERR, "Fix graphics/chunk is currently not compatible with 2d systems");
}

/* ---------------------------------------------------------------------- */

FixGraphicsChunk::~FixGraphicsChunk()
{
  delete[] id_chunk;
  delete[] id_region;
  memory->destroy(chunkid);
  memory->destroy(atomrad);
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsChunk::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsChunk::init()
{
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(id_chunk));
  if (!cchunk)
    error->all(FLERR, Error::NOLASTLINE,
               "Chunk/atom compute {} does not exist or is incorrect style for fix {}", id_chunk,
               style);

  if (id_region) {
    region = domain->get_region_by_id(id_region);
    if (!region)
      error->all(FLERR, Error::NOLASTLINE, "Region {} for fix graphics/chunk does not exist",
                 id_region);
  }
}

/* ---------------------------------------------------------------------- */

void FixGraphicsChunk::setup(int)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsChunk::end_of_step()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
  imgobjs = nullptr;
  imgparms = nullptr;
  numobjs = 0;

  // invoke chunk/atom compute and get chunk assignments

  modify->clearstep_compute();

  // update region if necessary

  if (region) region->prematch();

  int nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int *const mask = atom->mask;
  const int *const type = atom->type;
  const double *const *const x = atom->x;
  const int *const ichunk = cchunk->ichunk;

  // transfer chunk-ID data to chunkid array and communicate it to ghost atoms

  if (maxchunk < atom->nmax) {
    maxchunk = atom->nmax;
    memory->grow(chunkid, maxchunk, "graphics/chunk:chunkid");
    memory->grow(atomrad, maxchunk, "graphics/chunk:atomrad");
  }

  // determine per-atom radius: use per-atom property if available,
  // otherwise get sigma from potential
  // else use fixed radius from input or its default

  int dim = 0;
  double **sigma = nullptr;
  if (force->pair) {
    sigma = (double **) force->pair->extract("sigma", dim);
    if (dim != 2) sigma = nullptr;
  }

  double maxradius, oneradius = 0.0;
  bool has_peratom_radius = (atom->radius != nullptr);
  bool has_pertype_radius = (sigma != nullptr);

  for (int i = 0; i < nlocal; ++i) {
    chunkid[i] = (mask[i] & groupbit) ? ichunk[i] : 0;
    atomrad[i] = radius;
    if (!has_global_radius) {
      if (has_peratom_radius)
        atomrad[i] = atom->radius[i];
      else if (has_pertype_radius)
        atomrad[i] = 0.5 * sigma[type[i]][type[i]];
    }
    // fall back to global default radius if current value <= 0.0
    if (atomrad[i] <= 0.0) atomrad[i] = radius;
    oneradius = std::max(oneradius, atomrad[i]);
  }
  comm->forward_comm(this);
  MPI_Allreduce(&oneradius, &maxradius, 1, MPI_DOUBLE, MPI_MAX, world);

  // gather points per chunk along with their atom types
  // chunks are 1-based, ichunk[i] == 0 means atom is not in any chunk

  std::vector<std::vector<int>> chunk_atoms(nchunk);
  for (int i = 0; i < nall; ++i) {
    // skip atoms not in fix group
    if (!(mask[i] & groupbit)) continue;

    // skip atoms not in region, if enabled
    if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

    // skip atoms without assigned chunk
    int ic = chunkid[i];
    if (ic < 1 || ic > nchunk) continue;

    chunk_atoms[ic - 1].push_back(i);
  }

  // use automatic minimum cluster splitting distance, if not set explicitly
  if (mindist <= 0.0) mindist = 4.0 * maxradius;

  std::vector<std::vector<int>> new_chunks;

  // split chunks that are too far apart.
  for (const auto &chunk : chunk_atoms) {
    auto split = splitByDistance(chunk, x, mindist);
    for (const auto &c : split) new_chunks.emplace_back(c);
  }
  chunk_atoms.clear();

  // build convex hulls for each chunk and collect all TRINORM objects

  struct ObjData {
    double type0, type1, type2;
    vec3 v0, v1, v2;
    vec3 n0, n1, n2;
  };

  const double *const sublo = domain->sublo;
  const double *const subhi = domain->subhi;

  std::vector<ObjData> all_objs;
  ImageObjects::ConvexHullObj hull;

  for (const auto &iatoms : new_chunks) {
    const auto natoms = iatoms.size();
    if (iatoms.empty()) continue;

    // create list of shifted points
    int num_points = 1;
    std::vector<vec3> pts;

    // if number of atoms in cluster is maxreplace or smaller we replace the atom
    // positions with vectices from an icosahedron scaled to radius.
    if ((int)natoms <= maxreplace) {
      num_points = NUM_POINTS;
      pts.reserve(num_points * natoms);
      for (const auto i : iatoms) {
        for (const auto &ico : icosahedron) {
          vec3 pos{atomrad[i] * ico + vec3{x[i][0], x[i][1], x[i][2]}};
          if (!clip || domain->inside(pos.data())) pts.push_back(pos);
        }
      }
    } else {

      // shift atom positions away from center by radius

      vec3 center{0.0, 0.0, 0.0};
      for (const auto i : iatoms) {
        center[0] += x[i][0];
        center[1] += x[i][1];
        center[2] += x[i][2];
      }
      center[0] /= double(natoms);
      center[1] /= double(natoms);
      center[2] /= double(natoms);

      pts.reserve(natoms);
      for (const auto i : iatoms) {
        vec3 extra{x[i][0] - center[0], x[i][1] - center[1], x[i][2] - center[2]};
        MathExtra::norm3(extra.data());
        vec3 pos{atomrad[i] * extra + vec3{x[i][0], x[i][1], x[i][2]}};
        if (!clip || domain->inside(pos.data())) pts.push_back(pos);
      }
    }

    // build convex hull
    hull.build(pts, smooth, alpha);

    const auto &tris = hull.get_triangles();
    const auto &norms = hull.get_normals();
    const auto &cidx = hull.get_color_indices();

    for (size_t t = 0; t < tris.size(); ++t) {

      // get index into original atoms array for access to atom type
      int ci0 = cidx[t][0] / NUM_POINTS;
      int ci1 = cidx[t][1] / NUM_POINTS;
      int ci2 = cidx[t][2] / NUM_POINTS;
      if ((ci0 < 0) || (ci0 >= (int) iatoms.size())) ci0 = 0;
      if ((ci1 < 0) || (ci1 >= (int) iatoms.size())) ci1 = 0;
      if ((ci2 < 0) || (ci2 >= (int) iatoms.size())) ci2 = 0;

      bool addme = true;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if ((tris[t][i][j] < sublo[j] - mindist) || (tris[t][i][j] > subhi[j] + mindist))
            addme = false;

      if (addme)
        all_objs.push_back({double(type[iatoms[ci0]]), double(type[iatoms[ci1]]),
                            double(type[iatoms[ci2]]), tris[t][0], tris[t][1], tris[t][2],
                            norms[t][0], norms[t][1], norms[t][2]});
    }
  }
  new_chunks.clear();

  // allocate and fill imgobjs and imgparms arrays

  numobjs = static_cast<int>(all_objs.size());
  if (numobjs > 0) {
    memory->create(imgobjs, numobjs, "fix_graphics_chunk:imgobjs");
    memory->create(imgparms, numobjs, 21, "fix_graphics_chunk:imgparms");

    for (int n = 0; n < numobjs; ++n) {
      const auto &od = all_objs[n];
      imgobjs[n] = Graphics::TRINORM;
      imgparms[n][0] = od.type0;
      imgparms[n][1] = od.type1;
      imgparms[n][2] = od.type2;
      imgparms[n][3] = od.v0[0];
      imgparms[n][4] = od.v0[1];
      imgparms[n][5] = od.v0[2];
      imgparms[n][6] = od.v1[0];
      imgparms[n][7] = od.v1[1];
      imgparms[n][8] = od.v1[2];
      imgparms[n][9] = od.v2[0];
      imgparms[n][10] = od.v2[1];
      imgparms[n][11] = od.v2[2];
      imgparms[n][12] = od.n0[0];
      imgparms[n][13] = od.n0[1];
      imgparms[n][14] = od.n0[2];
      imgparms[n][15] = od.n1[0];
      imgparms[n][16] = od.n1[1];
      imgparms[n][17] = od.n1[2];
      imgparms[n][18] = od.n2[0];
      imgparms[n][19] = od.n2[1];
      imgparms[n][20] = od.n2[2];
    }
  }

  modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsChunk::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; ++i) {
    int j = list[i];
    buf[m++] = chunkid[j];
    buf[m++] = atomrad[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsChunk::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;

  for (int i = first; i < last; ++i) {
    chunkid[i] = (int) buf[m++];
    atomrad[i] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsChunk::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}

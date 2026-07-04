/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Joel Clemmer (SNL), Dan Bolintineanu (SNL)
----------------------------------------------------------------------- */

#include "fix_surface_local.h"

#include "atom.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_move.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "stl_reader.h"
#include "update.h"

#include <cstring>
#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static constexpr double EPSILON = 0.001;
static constexpr int NBIN = 100;
static constexpr double BIG = 1.0e20;
static constexpr int MAXTRIPOINT = 24;
static constexpr int DELTA_CONNECT = 1024;
static constexpr int DELTA_RVOUS = 1024;    // must be >= 8

enum { INTERNAL = 0, EXTERNAL, UNCONNECTED };

static constexpr double FLATTHRESH =
    0.00015230484360876085;        // = 1.0-cos(MY_PI/180.0); = 1 degree
static constexpr int RVOUS = 1;    // 0 for irregular, 1 for all2all

enum { MOLTEMPLATE, STLFILE };
enum { LAYOUT_UNIFORM, LAYOUT_NONUNIFORM, LAYOUT_TILED };    // several files

// allocate space for static class variable

FixSurfaceLocal *FixSurfaceLocal::fptr;

/* ---------------------------------------------------------------------- */

FixSurfaceLocal::FixSurfaceLocal(LAMMPS *lmp, int narg, char **arg) :
    FixSurface(lmp, narg, arg), connect2d(nullptr), connect3d(nullptr), atom2connect(nullptr),
    tstr(nullptr), input_modes(nullptr), input_stypes(nullptr), input_smols(nullptr),
    input_sources(nullptr), avec_line(nullptr), avec_tri(nullptr), ipc(nullptr), tpc(nullptr),
    pool2d(nullptr), pool3d(nullptr), connect2atom(nullptr), neigh_p1(nullptr), neigh_p2(nullptr),
    neigh_e1(nullptr), neigh_e2(nullptr), neigh_e3(nullptr), neigh_c1(nullptr), neigh_c2(nullptr),
    neigh_c3(nullptr), points(nullptr), lines(nullptr), tris(nullptr), connect2dall(nullptr),
    connect3dall(nullptr), plines(nullptr), etris(nullptr), ctris(nullptr), pts(nullptr),
    bincount(nullptr), binfirst(nullptr)
{
  create_attribute = 1;

  dimension = domain->dimension;

  // process zero or more inputs
  // just store info for use in post_constructor()

  ninput = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "input") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix surface/local input", error);
      if (strcmp(arg[iarg + 1], "mol") == 0) {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix surface/local input mol", error);
        input_modes = (int *) memory->srealloc(input_modes, (ninput + 1) * sizeof(int),
                                               "surface/local:input_modes");
        input_sources = (char **) memory->srealloc(input_sources, (ninput + 1) * sizeof(char **),
                                                   "surface/local:input_sources");
        input_modes[ninput] = MOLTEMPLATE;
        int n = strlen(arg[iarg + 2]) + 1;
        char *sourceID = new char[n];
        strcpy(sourceID, arg[iarg + 2]);
        input_sources[ninput] = sourceID;
        iarg += 3;
      } else if (strcmp(arg[iarg + 1], "stl") == 0) {
        if (iarg + 5 > narg) utils::missing_cmd_args(FLERR, "fix surface/local input stl", error);
        input_modes = (int *) memory->srealloc(input_modes, (ninput + 1) * sizeof(int),
                                               "surface/local:input_modes");
        input_sources = (char **) memory->srealloc(input_sources, (ninput + 1) * sizeof(char **),
                                                   "surface/local:input_sources");
        input_stypes = (int *) memory->srealloc(input_stypes, (ninput + 1) * sizeof(int),
                                                "surface/local:input_stypes");
        input_smols = (tagint *) memory->srealloc(input_smols, (ninput + 1) * sizeof(tagint),
                                                  "surface/local:input_smols");
        input_modes[ninput] = STLFILE;
        int stype = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
        input_stypes[ninput] = stype;
        tagint smol = utils::tnumeric(FLERR, arg[iarg + 3], false, lmp);
        input_smols[ninput] = smol;
        int n = strlen(arg[iarg + 4]) + 1;
        char *sourceID = new char[n];
        strcpy(sourceID, arg[iarg + 4]);
        input_sources[ninput] = sourceID;
        iarg += 5;
      } else
        error->all(FLERR, iarg + 1, "Unknown fix surface/local input keyword: {}", arg[iarg + 1]);
    } else
      break;

    ninput++;
  }

  // optional command-line args
  // smaxtype overrides max surf type of input surfs
  // flat overrides FLATTHRESH of one degree

  flatthresh = FLATTHRESH;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "flat") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix surface/local flat", error);
      double flat = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if ((flat < 0.0) || (flat > 90.0))
        error->all(FLERR, iarg + 1, "Invalid value {} for fix surface/local flat keyword", flat);
      flatthresh = 1.0 - cos(MY_PI * flat / 180.0);
      iarg += 2;
    } else
      error->all(FLERR, iarg, "Unknown fix surface/local keyword: {}", arg[iarg]);
  }

  // error check

  if (dimension == 2) {
    avec_line = (AtomVecLine *) atom->style_match("line");
    if (!avec_line) error->all(FLERR, 2, "Fix surface/local requires atom style line");
  } else if (dimension == 3) {
    avec_tri = (AtomVecTri *) atom->style_match("tri");
    if (!avec_tri) error->all(FLERR, 2, "Fix surface/local requires atom style tri");
  }

  // initializations

  atom2connect = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(2);

  nlocal_connect = nghost_connect = nmax_connect = 0;
  connect2d = nullptr;
  connect3d = nullptr;
  pool2d = nullptr;
  pool3d = nullptr;
  connect2atom = nullptr;

  ipc = new MyPoolChunk<int>(1, MAXTRIPOINT, 6);
  tpc = new MyPoolChunk<int>(1, MAXTRIPOINT, 6);

  flag_complete = 0;
  epssq = -1.0;
}

/* ---------------------------------------------------------------------- */

FixSurfaceLocal::~FixSurfaceLocal()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id, 0);

  memory->destroy(atom2connect);
  memory->destroy(connect2atom);

  // return all connection vector memory to ipc and tpc allocators

  if (dimension == 2) {
    int nall = nlocal_connect + nghost_connect;
    for (int i = 0; i < nall; i++) {
      if (connect2d[i].neigh_p1) {
        tpc->put(pool2d[i].neigh_p1);
        ipc->put(pool2d[i].pwhich_p1);
        ipc->put(pool2d[i].nside_p1);
        ipc->put(pool2d[i].aflag_p1);
        ipc->put(pool2d[i].fflag_p1);
      }
      if (connect2d[i].neigh_p2) {
        tpc->put(pool2d[i].neigh_p2);
        ipc->put(pool2d[i].pwhich_p2);
        ipc->put(pool2d[i].nside_p2);
        ipc->put(pool2d[i].aflag_p2);
        ipc->put(pool2d[i].fflag_p2);
      }
    }
    memory->sfree(connect2d);
  } else {
    int nall = nlocal_connect + nghost_connect;
    for (int i = 0; i < nall; i++) {
      if (connect3d[i].neigh_e1) {
        tpc->put(pool3d[i].neigh_e1);
        ipc->put(pool3d[i].ewhich_e1);
        ipc->put(pool3d[i].nside_e1);
        ipc->put(pool3d[i].aflag_e1);
        ipc->put(pool3d[i].fflag_e1);
      }
      if (connect3d[i].neigh_e2) {
        tpc->put(pool3d[i].neigh_e2);
        ipc->put(pool3d[i].ewhich_e2);
        ipc->put(pool3d[i].nside_e2);
        ipc->put(pool3d[i].aflag_e2);
        ipc->put(pool3d[i].fflag_e2);
      }
      if (connect3d[i].neigh_e3) {
        tpc->put(pool3d[i].neigh_e3);
        ipc->put(pool3d[i].ewhich_e3);
        ipc->put(pool3d[i].nside_e3);
        ipc->put(pool3d[i].aflag_e3);
        ipc->put(pool3d[i].fflag_e3);
      }
      if (connect3d[i].neigh_c1) {
        tpc->put(pool3d[i].neigh_c1);
        ipc->put(pool3d[i].cwhich_c1);
        ipc->put(pool3d[i].nside_c1);
        ipc->put(pool3d[i].fflag_c1);
      }
      if (connect3d[i].neigh_c2) {
        tpc->put(pool3d[i].neigh_c2);
        ipc->put(pool3d[i].cwhich_c2);
        ipc->put(pool3d[i].nside_c2);
        ipc->put(pool3d[i].fflag_c2);
      }
      if (connect3d[i].neigh_c3) {
        tpc->put(pool3d[i].neigh_c3);
        ipc->put(pool3d[i].cwhich_c3);
        ipc->put(pool3d[i].nside_c3);
        ipc->put(pool3d[i].fflag_c3);
      }
    }
    memory->sfree(connect3d);
  }

  memory->sfree(pool2d);
  memory->sfree(pool3d);

  delete ipc;
  delete tpc;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceLocal::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;    // only needed for DEBUG tests
  return mask;
}

/* ----------------------------------------------------------------------
   one-time setup of distributed lines/tri and their connectivity
   must be done in post_constructor() for MOLTEMPLATE and STLFILE
     they add owned lines/tris to AtomVec class
     its grow() method makes a callback to grow_arrays() in this fix
     callback can't be invoked unless fix is fully instantiated
------------------------------------------------------------------------- */

void FixSurfaceLocal::post_constructor()
{
  // if line/tri particles already exist from data file, initialize their connectivity

  nlocal0 = 0;
  if (check_exist()) {
    if (comm->me == 0) utils::logmesg(lmp, "Connecting line/tri particles ...\n");

    if (dimension == 2)
      connectivity2d_local();
    else
      connectivity3d_local();
    nlocal0 = atom->nlocal;
  }

  // loop over instances of input keyword

  if (ninput) {
    npoints = maxpoints = 0;
    nlines = ntris = 0;
    points = nullptr;
    lines = nullptr;
    tris = nullptr;

    std::map<std::tuple<double, double, double, int>, int> *hash =
        new std::map<std::tuple<double, double, double, int>, int>();

    for (int i = 0; i < ninput; i++) {
      int mode = input_modes[i];
      char *sourceID = input_sources[i];

      if (comm->me == 0) {
        if (mode == MOLTEMPLATE)
          utils::logmesg(lmp, "Converting molecule file to line/tri particles ...\n");
        if (mode == STLFILE) utils::logmesg(lmp, "Reading STL file for triangle particles ...\n");
      }

      // read in lines/tris from appropriate source
      // each proc builds global data structs of all points/lines/tris

      if (mode == MOLTEMPLATE)
        extract_from_molecule(sourceID, hash, npoints, maxpoints, points, nlines, lines, ntris,
                              tris);
      if (mode == STLFILE) {
        int stype = input_stypes[i];
        tagint smol = input_smols[i];
        extract_from_stlfile(sourceID, stype, smol, hash, npoints, maxpoints, points, ntris, tris);
      }
    }

    delete hash;

    memory->sfree(input_modes);
    for (int i = 0; i < ninput; i++) delete[] input_sources[i];
    memory->sfree(input_sources);
    memory->sfree(input_stypes);
    memory->sfree(input_smols);

    // each proc infers global connectivity from global data structs
    // distribute lines/surfs across procs, based on center pt coords

    if (dimension == 2) {
      connectivity2d_global(npoints, nlines, lines, connect2dall, neigh_p1, neigh_p2);
      assign2d();
    } else {
      connectivity3d_global(npoints, ntris, tris, connect3dall, neigh_e1, neigh_e2, neigh_e3,
                            neigh_c1, neigh_c2, neigh_c3);
      assign3d();
    }

    // delete global data structs

    memory->sfree(points);
    memory->sfree(lines);
    memory->sfree(tris);
  }

  // check that line/triangle particles now exist

  if (!check_exist())
    error->all(FLERR, Error::NOLASTLINE, "Fix surface/local defines no line/triangle particles");

  // confirm all tri/line atoms are defined by this fix

  for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
    if (domain->dimension == 2) {
      if ((atom->line[i] >= 0) && (atom2connect[i] == -1))
        error->one(FLERR, Error::NOLASTLINE,
                   "Unassociated line atom {} identified by fix surface/local", atom->tag[i]);
    } else {
      if ((atom->tri[i] >= 0) && (atom2connect[i] == -1))
        error->one(FLERR, Error::NOLASTLINE,
                   "Unassociated tri atom {} identified by fix surface/local", atom->tag[i]);
    }
  }

  // set max size for comm of connection info
  // 2d = 2 end points, 4 vectors, each of length npmaxall
  // 3d = 3 edges, 4 vectors, each of length nemaxall
  //      plus 3 corner points, 2 vectors, each of length ncmaxall

  if (dimension == 2) {
    int nlocal = atom->nlocal;
    int *line = atom->line;
    int iconnect;

    int npmax = 0;
    for (int i = 0; i < nlocal; i++) {
      if (line[i] < 0) continue;
      iconnect = atom2connect[i];
      npmax = MAX(npmax, connect2d[iconnect].np1);
      npmax = MAX(npmax, connect2d[iconnect].np2);
    }

    int npmaxall;
    MPI_Allreduce(&npmax, &npmaxall, 1, MPI_INT, MPI_MAX, world);
    comm_border = comm_forward = 3 + 2 * 5 * npmaxall;

  } else if (dimension == 3) {
    int nlocal = atom->nlocal;
    int *tri = atom->tri;
    int iconnect;

    int nemax = 0;
    int ncmax = 0;
    for (int i = 0; i < nlocal; i++) {
      if (tri[i] < 0) continue;
      iconnect = atom2connect[i];
      nemax = MAX(nemax, connect3d[iconnect].ne1);
      nemax = MAX(nemax, connect3d[iconnect].ne2);
      nemax = MAX(nemax, connect3d[iconnect].ne3);
      ncmax = MAX(ncmax, connect3d[iconnect].nc1);
      ncmax = MAX(ncmax, connect3d[iconnect].nc2);
      ncmax = MAX(ncmax, connect3d[iconnect].nc3);
    }

    int nemaxall, ncmaxall;
    MPI_Allreduce(&nemax, &nemaxall, 1, MPI_INT, MPI_MAX, world);
    MPI_Allreduce(&ncmax, &ncmaxall, 1, MPI_INT, MPI_MAX, world);
    comm_border = comm_forward = 7 + 3 * 5 * nemaxall + 3 * 4 * ncmaxall;
  }

  // error checks on duplicate surfs or zero-size surfs

  if (dimension == 2)
    check2d();
  else
    check3d();

  // print stats on surfs and their connectivity

  if (dimension == 2)
    stats2d();
  else
    stats3d();
}

/* ----------------------------------------------------------------------
   error check that each surf is assigned to no more than one motion
   error check that connected surfs have same motion
   done at beginning of each run in case fix_modify move has changed motions
   do at setup() when connectivity for owned/ghost surfs is fully initialized
------------------------------------------------------------------------- */

void FixSurfaceLocal::setup(int /*vflag*/)
{
  // count instances of FixMove
  // just return if no FixMove commands
  // movelist = pointers to each instance

  auto movelist = modify->get_fix_by_style("^move");
  if (movelist.size() == 0) return;

  // check for inconsistent surf motion
  // error if a single surf is assigned to multiple motions
  // error if two surfs are connected but are assigned to different motions
  // check before every run, b/c fix_modify move can be used between runs

  if (dimension == 2) {

    int *line = atom->line;
    int *mask = atom->mask;
    tagint *tag = atom->tag;
    int nlocal = atom->nlocal;

    // check that each line is assigned to no more than one FixMove

    int mcount;
    int ecount = 0;

    for (int i = 0; i < nlocal; i++) {
      if (line[i] < 0) continue;

      mcount = 0;
      for (auto *ifix : movelist)
        if (mask[i] & ifix->groupbit) mcount++;
      if (mcount > 1) ecount++;
    }

    int selfcount;
    MPI_Allreduce(&ecount, &selfcount, 1, MPI_INT, MPI_SUM, world);

    if (selfcount)
      error->all(FLERR, Error::NOLASTLINE,
                 "Multiple fix surface/local motions assigned to {} lines", selfcount);

    // check connections to either endpoint of each line
    // avoid double counting by jtag < itag check

    tagint itag;
    int m, iline, jline, np1, np2, imotion, jmotion;
    ecount = 0;

    for (int i = 0; i < nlocal; i++) {
      if (line[i] < 0) continue;
      itag = tag[i];
      iline = atom2connect[i];

      imotion = -1;
      m = 0;
      for (auto *ifix : movelist) {
        if (mask[i] & ifix->groupbit) imotion = m;
        ++m;
      }

      np1 = connect2d[iline].np1;
      np2 = connect2d[iline].np2;
      for (int j = 0; j < np1; j++) {
        jline = atom->map(connect2d[iline].neigh_p1[j]);
        if (tag[jline] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jline] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }
      for (int j = 0; j < np2; j++) {
        jline = atom->map(connect2d[iline].neigh_p2[j]);
        if (tag[jline] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jline] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }
    }

    int paircount;
    MPI_Allreduce(&ecount, &paircount, 1, MPI_INT, MPI_SUM, world);

    if (paircount)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix surface/local motions assigned to "
                 "{} connected line pairs are inconsistent",
                 paircount);

  } else if (dimension == 3) {

    int *tri = atom->tri;
    int *mask = atom->mask;
    tagint *tag = atom->tag;
    int nlocal = atom->nlocal;

    // check that each tri is assigned to no more than one FixMove

    int mcount;
    int ecount = 0;

    for (int i = 0; i < nlocal; i++) {
      if (tri[i] < 0) continue;

      mcount = 0;
      for (auto *ifix : movelist)
        if (mask[i] & ifix->groupbit) mcount++;
      if (mcount > 1) ecount++;
    }

    int selfcount;
    MPI_Allreduce(&ecount, &selfcount, 1, MPI_INT, MPI_SUM, world);

    if (selfcount)
      error->all(FLERR, Error::NOLASTLINE, "Multiple fix surface/local motions assigned to {} tris",
                 selfcount);

    // check connections to each tri's 3 edges and 3 corner pts
    // avoid double counting by jtag < itag check
    // but edge and corner connections both contribute to error count

    tagint itag;
    int m, itri, jtri, ne1, ne2, ne3, nc1, nc2, nc3, imotion, jmotion;
    ecount = 0;

    for (int i = 0; i < nlocal; i++) {
      if (tri[i] < 0) continue;
      itag = tag[i];
      itri = atom2connect[i];

      imotion = -1;
      m = 0;
      for (auto *ifix : movelist) {
        if (mask[i] & ifix->groupbit) imotion = m;
        ++m;
      }

      ne1 = connect3d[itri].ne1;
      ne2 = connect3d[itri].ne2;
      ne3 = connect3d[itri].ne3;

      for (int j = 0; j < ne1; j++) {
        jtri = atom->map(connect3d[itri].neigh_e1[j]);
        if (tag[jtri] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jtri] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }
      for (int j = 0; j < ne2; j++) {
        jtri = atom->map(connect3d[itri].neigh_e2[j]);
        if (tag[jtri] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jtri] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }
      for (int j = 0; j < ne3; j++) {
        jtri = atom->map(connect3d[itri].neigh_e3[j]);
        if (tag[jtri] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jtri] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }

      nc1 = connect3d[itri].nc1;
      nc2 = connect3d[itri].nc2;
      nc3 = connect3d[itri].nc3;

      for (int j = 0; j < nc1; j++) {
        jtri = atom->map(connect3d[itri].neigh_c1[j]);
        if (tag[jtri] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jtri] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }
      for (int j = 0; j < nc2; j++) {
        jtri = atom->map(connect3d[itri].neigh_c2[j]);
        if (tag[jtri] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jtri] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }
      for (int j = 0; j < nc3; j++) {
        jtri = atom->map(connect3d[itri].neigh_c3[j]);
        if (tag[jtri] < itag) continue;
        jmotion = -1;
        m = 0;
        for (auto *ifix : movelist) {
          if (mask[jtri] & ifix->groupbit) jmotion = m;
          ++m;
        }
        if (imotion != jmotion) ecount++;
      }
    }

    int paircount;
    MPI_Allreduce(&ecount, &paircount, 1, MPI_INT, MPI_SUM, world);

    if (paircount)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix surface/local motions assigned to "
                 "{} connected tri pairs are inconsistent",
                 paircount);
  }
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixSurfaceLocal::grow_arrays(int nmax)
{
  memory->grow(atom2connect, nmax, "surface/local:atom2connect");
}

/* ---------------------------------------------------------------------- */

void FixSurfaceLocal::setup_pre_neighbor()
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, Error::NOLASTLINE, "Fix surface/local requires an atom map");

  // check ghost cutoff large enough to walk all possible connections, > 2x max radius
  //  not perfect (doesn't check multi), but at least a warning

  const double cutghost = MAX(neighbor->cutneighmax, comm->cutghostuser);
  if (cutghost < 2 * max_radius)
    error->warning(FLERR, "Maximum triangle diameter {} may be less than ghost cutoff {}",
                   2 * max_radius, cutghost);

  // one-time calculation of remaining fields in Connect2d/3d
  // cannot do until now, b/c need ghost connection info via border comm
  // then re-communicate new owned-particle connectivity to ghost particles

  if (!flag_complete) {
    if (dimension == 2)
      connectivity2d_complete();
    else
      connectivity3d_complete();
    flag_complete = 1;

    clear_bonus();
    comm->forward_comm(this);
  }

  pre_neighbor();
}

/* ----------------------------------------------------------------------
   only for DEBUGGING
   check that atom2connect is identical to line/tri
   check that atom2connect is consistent with connect2atom
------------------------------------------------------------------------- */

void FixSurfaceLocal::pre_neighbor()
{
  int n = atom->nlocal + atom->nghost;

  int *surf = atom->line;
  if (dimension == 3) surf = atom->tri;

  int count1 = 0;
  int count2 = 0;

  for (int i = 0; i < n; i++) {
    if (surf[i] < 0) continue;

    if (atom2connect[i] != surf[i]) {
      count1++;
      continue;
    }
    if (connect2atom[atom2connect[i]] != i) count2++;
  }

  int all1, all2;
  MPI_Allreduce(&count1, &all1, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&count2, &all2, 1, MPI_INT, MPI_SUM, world);

  if ((all1 || all2) && (comm->me == 0))
    error->warning(FLERR, "Fix surface/local atom2connect vector mis-match: {} {}: {}\n", all1,
                   all2, update->ntimestep);
}

/* ----------------------------------------------------------------------
   grow connect array by DELTA_CONNECT
------------------------------------------------------------------------- */

void FixSurfaceLocal::grow_connect()
{
  nmax_connect = nmax_connect / DELTA_CONNECT * DELTA_CONNECT;
  bigint newmax = (bigint) nmax_connect + DELTA_CONNECT;
  if (newmax > MAXSMALLINT)
    error->one(FLERR, Error::NOLASTLINE, "Too many fix surface/local connections");
  nmax_connect = newmax;

  memory->grow(connect2atom, nmax_connect, "surface/local:connect2atom");

  if (dimension == 2) {
    connect2d = (Connect2d *) memory->srealloc(connect2d, nmax_connect * sizeof(Connect2d),
                                               "surface/local:connect2d");
    pool2d =
        (Pool2d *) memory->srealloc(pool2d, nmax_connect * sizeof(Pool2d), "surface/local:pool2d");
  } else {
    connect3d = (Connect3d *) memory->srealloc(connect3d, nmax_connect * sizeof(Connect3d),
                                               "surface/local:connect3d");
    pool3d =
        (Pool3d *) memory->srealloc(pool3d, nmax_connect * sizeof(Pool3d), "surface/local:pool3d");
  }
}

/* ----------------------------------------------------------------------
   copy atom I to atom J with optional delflag
   I or J or both can be a line/triangle with connection info
------------------------------------------------------------------------- */

void FixSurfaceLocal::copy_arrays(int i, int j, int delflag)
{
  // if overwriting atom J via delflag and J is a line/tri with connection data:
  //   return J's vector memory to TCP before memcpy()
  //   shrink list of owned connections by copying last element to Kth location

  if (dimension == 2) {
    if (delflag && atom2connect[j] >= 0) {
      int k = atom2connect[j];
      if (connect2d[k].neigh_p1) {
        tpc->put(pool2d[k].neigh_p1);
        ipc->put(pool2d[k].pwhich_p1);
        ipc->put(pool2d[k].nside_p1);
        ipc->put(pool2d[k].aflag_p1);
        ipc->put(pool2d[k].fflag_p1);
      }
      if (connect2d[k].neigh_p2) {
        tpc->put(pool2d[k].neigh_p2);
        ipc->put(pool2d[k].pwhich_p2);
        ipc->put(pool2d[k].nside_p2);
        ipc->put(pool2d[k].aflag_p2);
        ipc->put(pool2d[k].fflag_p2);
      }
      memcpy(&connect2d[k], &connect2d[nlocal_connect - 1], sizeof(Connect2d));
      memcpy(&pool2d[k], &pool2d[nlocal_connect - 1], sizeof(Pool2d));
      int iatom = connect2atom[nlocal_connect - 1];
      atom2connect[iatom] = k;
      connect2atom[k] = iatom;
      nlocal_connect--;
    }
  } else {
    if (delflag && atom2connect[j] >= 0) {
      int k = atom2connect[j];
      if (connect3d[k].neigh_e1) {
        tpc->put(pool3d[k].neigh_e1);
        ipc->put(pool3d[k].ewhich_e1);
        ipc->put(pool3d[k].nside_e1);
        ipc->put(pool3d[k].aflag_e1);
        ipc->put(pool3d[k].fflag_e1);
      }
      if (connect3d[k].neigh_e2) {
        tpc->put(pool3d[k].neigh_e2);
        ipc->put(pool3d[k].ewhich_e2);
        ipc->put(pool3d[k].nside_e2);
        ipc->put(pool3d[k].aflag_e2);
        ipc->put(pool3d[k].fflag_e2);
      }
      if (connect3d[k].neigh_e3) {
        tpc->put(pool3d[k].neigh_e3);
        ipc->put(pool3d[k].ewhich_e3);
        ipc->put(pool3d[k].nside_e3);
        ipc->put(pool3d[k].aflag_e3);
        ipc->put(pool3d[k].fflag_e3);
      }
      if (connect3d[k].neigh_c1) {
        tpc->put(pool3d[k].neigh_c1);
        ipc->put(pool3d[k].cwhich_c1);
        ipc->put(pool3d[k].nside_c1);
        ipc->put(pool3d[k].fflag_c1);
      }
      if (connect3d[k].neigh_c2) {
        tpc->put(pool3d[k].neigh_c2);
        ipc->put(pool3d[k].cwhich_c2);
        ipc->put(pool3d[k].nside_c2);
        ipc->put(pool3d[k].fflag_c2);
      }
      if (connect3d[k].neigh_c3) {
        tpc->put(pool3d[k].neigh_c3);
        ipc->put(pool3d[k].cwhich_c3);
        ipc->put(pool3d[k].nside_c3);
        ipc->put(pool3d[k].fflag_c3);
      }
      memcpy(&connect3d[k], &connect3d[nlocal_connect - 1], sizeof(Connect3d));
      memcpy(&pool3d[k], &pool3d[nlocal_connect - 1], sizeof(Pool3d));
      int iatom = connect2atom[nlocal_connect - 1];
      atom2connect[iatom] = k;
      connect2atom[k] = iatom;
      nlocal_connect--;
    }
  }

  // if atom I has connection data, reset connect2atom[I] to loc J
  // do NOT do this if self-copy (I=J) since I's connection data is already deleted above

  if (atom2connect[i] >= 0 && i != j) connect2atom[atom2connect[i]] = j;
  atom2connect[j] = atom2connect[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixSurfaceLocal::set_arrays(int i)
{
  atom2connect[i] = -1;
}

/* ----------------------------------------------------------------------
   clear ghost info in connect array
   called before ghosts are recommunicated in comm and irregular
   return all vector memory to TCP
------------------------------------------------------------------------- */

void FixSurfaceLocal::clear_bonus()
{
  if (dimension == 2) {
    int nall = nlocal_connect + nghost_connect;
    for (int i = nlocal_connect; i < nall; i++) {
      if (connect2d[i].neigh_p1) {
        tpc->put(pool2d[i].neigh_p1);
        ipc->put(pool2d[i].pwhich_p1);
        ipc->put(pool2d[i].nside_p1);
        ipc->put(pool2d[i].aflag_p1);
        ipc->put(pool2d[i].fflag_p1);
      }
      if (connect2d[i].neigh_p2) {
        tpc->put(pool2d[i].neigh_p2);
        ipc->put(pool2d[i].pwhich_p2);
        ipc->put(pool2d[i].nside_p2);
        ipc->put(pool2d[i].aflag_p2);
        ipc->put(pool2d[i].fflag_p2);
      }
    }
  } else if (dimension == 3) {
    int nall = nlocal_connect + nghost_connect;
    for (int i = nlocal_connect; i < nall; i++) {
      if (connect3d[i].neigh_e1) {
        tpc->put(pool3d[i].neigh_e1);
        ipc->put(pool3d[i].ewhich_e1);
        ipc->put(pool3d[i].nside_e1);
        ipc->put(pool3d[i].aflag_e1);
        ipc->put(pool3d[i].fflag_e1);
      }
      if (connect3d[i].neigh_e2) {
        tpc->put(pool3d[i].neigh_e2);
        ipc->put(pool3d[i].ewhich_e2);
        ipc->put(pool3d[i].nside_e2);
        ipc->put(pool3d[i].aflag_e2);
        ipc->put(pool3d[i].fflag_e2);
      }
      if (connect3d[i].neigh_e3) {
        tpc->put(pool3d[i].neigh_e3);
        ipc->put(pool3d[i].ewhich_e3);
        ipc->put(pool3d[i].nside_e3);
        ipc->put(pool3d[i].aflag_e3);
        ipc->put(pool3d[i].fflag_e3);
      }
      if (connect3d[i].neigh_c1) {
        tpc->put(pool3d[i].neigh_c1);
        ipc->put(pool3d[i].cwhich_c1);
        ipc->put(pool3d[i].nside_c1);
        ipc->put(pool3d[i].fflag_c1);
      }
      if (connect3d[i].neigh_c2) {
        tpc->put(pool3d[i].neigh_c2);
        ipc->put(pool3d[i].cwhich_c2);
        ipc->put(pool3d[i].nside_c2);
        ipc->put(pool3d[i].fflag_c2);
      }
      if (connect3d[i].neigh_c3) {
        tpc->put(pool3d[i].neigh_c3);
        ipc->put(pool3d[i].cwhich_c3);
        ipc->put(pool3d[i].nside_c3);
        ipc->put(pool3d[i].fflag_c3);
      }
    }
  }

  nghost_connect = 0;
}

/* ----------------------------------------------------------------------
   pack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixSurfaceLocal::pack_border(int n, int *list, double *buf)
{
  int m, ic;

  m = 0;

  if (dimension == 2) {
    int np1, np2;

    for (int i = 0; i < n; i++) {
      int j = list[i];
      if (atom2connect[j] < 0)
        buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        ic = atom2connect[j];

        np1 = connect2d[ic].np1;
        np2 = connect2d[ic].np2;
        buf[m++] = ubuf(np1).d;
        buf[m++] = ubuf(np2).d;
        buf[m++] = ubuf(connect2d[ic].external_pt[0]).d;
        buf[m++] = ubuf(connect2d[ic].external_pt[1]).d;
        if (np1)
          for (int k = 0; k < np1; k++) {
            buf[m++] = ubuf(connect2d[ic].neigh_p1[k]).d;
            buf[m++] = ubuf(connect2d[ic].pwhich_p1[k]).d;
            buf[m++] = ubuf(connect2d[ic].nside_p1[k]).d;
            buf[m++] = ubuf(connect2d[ic].aflag_p1[k]).d;
            buf[m++] = ubuf(connect2d[ic].fflag_p1[k]).d;
          }
        if (np2)
          for (int k = 0; k < np2; k++) {
            buf[m++] = ubuf(connect2d[ic].neigh_p2[k]).d;
            buf[m++] = ubuf(connect2d[ic].pwhich_p2[k]).d;
            buf[m++] = ubuf(connect2d[ic].nside_p2[k]).d;
            buf[m++] = ubuf(connect2d[ic].aflag_p2[k]).d;
            buf[m++] = ubuf(connect2d[ic].fflag_p2[k]).d;
          }
      }
    }

  } else {
    int ne1, ne2, ne3, nc1, nc2, nc3;

    for (int i = 0; i < n; i++) {
      int j = list[i];
      if (atom2connect[j] < 0)
        buf[m++] = ubuf(0).d;
      else {
        buf[m++] = ubuf(1).d;
        ic = atom2connect[j];

        ne1 = connect3d[ic].ne1;
        ne2 = connect3d[ic].ne2;
        ne3 = connect3d[ic].ne3;
        buf[m++] = ubuf(ne1).d;
        buf[m++] = ubuf(ne2).d;
        buf[m++] = ubuf(ne3).d;
        buf[m++] = ubuf(connect3d[ic].external_edge[0]).d;
        buf[m++] = ubuf(connect3d[ic].external_edge[1]).d;
        buf[m++] = ubuf(connect3d[ic].external_edge[2]).d;
        if (ne1)
          for (int k = 0; k < ne1; k++) {
            buf[m++] = ubuf(connect3d[ic].neigh_e1[k]).d;
            buf[m++] = ubuf(connect3d[ic].ewhich_e1[k]).d;
            buf[m++] = ubuf(connect3d[ic].nside_e1[k]).d;
            buf[m++] = ubuf(connect3d[ic].aflag_e1[k]).d;
            buf[m++] = ubuf(connect3d[ic].fflag_e1[k]).d;
          }
        if (ne2)
          for (int k = 0; k < ne2; k++) {
            buf[m++] = ubuf(connect3d[ic].neigh_e2[k]).d;
            buf[m++] = ubuf(connect3d[ic].ewhich_e2[k]).d;
            buf[m++] = ubuf(connect3d[ic].nside_e2[k]).d;
            buf[m++] = ubuf(connect3d[ic].aflag_e2[k]).d;
            buf[m++] = ubuf(connect3d[ic].fflag_e2[k]).d;
          }
        if (ne3)
          for (int k = 0; k < ne3; k++) {
            buf[m++] = ubuf(connect3d[ic].neigh_e3[k]).d;
            buf[m++] = ubuf(connect3d[ic].ewhich_e3[k]).d;
            buf[m++] = ubuf(connect3d[ic].nside_e3[k]).d;
            buf[m++] = ubuf(connect3d[ic].aflag_e3[k]).d;
            buf[m++] = ubuf(connect3d[ic].fflag_e3[k]).d;
          }

        nc1 = connect3d[ic].nc1;
        nc2 = connect3d[ic].nc2;
        nc3 = connect3d[ic].nc3;
        buf[m++] = ubuf(nc1).d;
        buf[m++] = ubuf(nc2).d;
        buf[m++] = ubuf(nc3).d;
        buf[m++] = ubuf(connect3d[ic].external_cor[0]).d;
        buf[m++] = ubuf(connect3d[ic].external_cor[1]).d;
        buf[m++] = ubuf(connect3d[ic].external_cor[2]).d;
        if (nc1)
          for (int k = 0; k < nc1; k++) {
            buf[m++] = ubuf(connect3d[ic].neigh_c1[k]).d;
            buf[m++] = ubuf(connect3d[ic].cwhich_c1[k]).d;
            buf[m++] = ubuf(connect3d[ic].nside_c1[k]).d;
            buf[m++] = ubuf(connect3d[ic].fflag_c1[k]).d;
          }
        if (nc2)
          for (int k = 0; k < nc2; k++) {
            buf[m++] = ubuf(connect3d[ic].neigh_c2[k]).d;
            buf[m++] = ubuf(connect3d[ic].cwhich_c2[k]).d;
            buf[m++] = ubuf(connect3d[ic].nside_c2[k]).d;
            buf[m++] = ubuf(connect3d[ic].fflag_c2[k]).d;
          }
        if (nc3)
          for (int k = 0; k < nc3; k++) {
            buf[m++] = ubuf(connect3d[ic].neigh_c3[k]).d;
            buf[m++] = ubuf(connect3d[ic].cwhich_c3[k]).d;
            buf[m++] = ubuf(connect3d[ic].nside_c3[k]).d;
            buf[m++] = ubuf(connect3d[ic].fflag_c3[k]).d;
          }
      }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values for border communication at re-neighboring
------------------------------------------------------------------------- */

int FixSurfaceLocal::unpack_border(int n, int first, double *buf)
{
  int m, last, flag;

  m = 0;
  last = first + n;

  if (dimension == 2) {
    int np1, np2;

    for (int i = first; i < last; i++) {
      flag = (int) ubuf(buf[m++]).i;
      if (flag == 0)
        atom2connect[i] = -1;
      else {
        int j = nlocal_connect + nghost_connect;
        if (j == nmax_connect) grow_connect();

        np1 = (int) ubuf(buf[m++]).i;
        np2 = (int) ubuf(buf[m++]).i;
        connect2d[j].np1 = np1;
        connect2d[j].np2 = np2;

        connect2d[j].external_pt[0] = (int) ubuf(buf[m++]).i;
        connect2d[j].external_pt[1] = (int) ubuf(buf[m++]).i;

        if (np1) {
          connect2d[j].neigh_p1 = tpc->get(np1, pool2d[j].neigh_p1);
          connect2d[j].pwhich_p1 = ipc->get(np1, pool2d[j].pwhich_p1);
          connect2d[j].nside_p1 = ipc->get(np1, pool2d[j].nside_p1);
          connect2d[j].aflag_p1 = ipc->get(np1, pool2d[j].aflag_p1);
          connect2d[j].fflag_p1 = ipc->get(np1, pool2d[j].fflag_p1);
          for (int k = 0; k < np1; k++) {
            connect2d[j].neigh_p1[k] = (tagint) ubuf(buf[m++]).i;
            connect2d[j].pwhich_p1[k] = (int) ubuf(buf[m++]).i;
            connect2d[j].nside_p1[k] = (int) ubuf(buf[m++]).i;
            connect2d[j].aflag_p1[k] = (int) ubuf(buf[m++]).i;
            connect2d[j].fflag_p1[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect2d[j].neigh_p1 = nullptr;
          connect2d[j].pwhich_p1 = nullptr;
          connect2d[j].nside_p1 = nullptr;
          connect2d[j].aflag_p1 = nullptr;
          connect2d[j].fflag_p1 = nullptr;
        }

        if (np2) {
          connect2d[j].neigh_p2 = tpc->get(np1, pool2d[j].neigh_p2);
          connect2d[j].pwhich_p2 = ipc->get(np1, pool2d[j].pwhich_p2);
          connect2d[j].nside_p2 = ipc->get(np1, pool2d[j].nside_p2);
          connect2d[j].aflag_p2 = ipc->get(np1, pool2d[j].aflag_p2);
          connect2d[j].fflag_p2 = ipc->get(np1, pool2d[j].fflag_p2);
          for (int k = 0; k < np1; k++) {
            connect2d[j].neigh_p2[k] = (tagint) ubuf(buf[m++]).i;
            connect2d[j].pwhich_p2[k] = (int) ubuf(buf[m++]).i;
            connect2d[j].nside_p2[k] = (int) ubuf(buf[m++]).i;
            connect2d[j].aflag_p2[k] = (int) ubuf(buf[m++]).i;
            connect2d[j].fflag_p2[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect2d[j].neigh_p2 = nullptr;
          connect2d[j].pwhich_p2 = nullptr;
          connect2d[j].nside_p2 = nullptr;
          connect2d[j].aflag_p2 = nullptr;
          connect2d[j].fflag_p2 = nullptr;
        }

        connect2atom[j] = i;
        atom2connect[i] = j;
        nghost_connect++;
      }
    }

  } else {
    int ne1, ne2, ne3, nc1, nc2, nc3;

    for (int i = first; i < last; i++) {
      flag = (int) ubuf(buf[m++]).i;
      if (flag == 0)
        atom2connect[i] = -1;
      else {
        int j = nlocal_connect + nghost_connect;
        if (j == nmax_connect) grow_connect();

        ne1 = (int) ubuf(buf[m++]).i;
        ne2 = (int) ubuf(buf[m++]).i;
        ne3 = (int) ubuf(buf[m++]).i;
        connect3d[j].ne1 = ne1;
        connect3d[j].ne2 = ne2;
        connect3d[j].ne3 = ne3;

        connect3d[j].external_edge[0] = (int) ubuf(buf[m++]).i;
        connect3d[j].external_edge[1] = (int) ubuf(buf[m++]).i;
        connect3d[j].external_edge[2] = (int) ubuf(buf[m++]).i;

        if (ne1) {
          connect3d[j].neigh_e1 = tpc->get(ne1, pool3d[j].neigh_e1);
          connect3d[j].ewhich_e1 = ipc->get(ne1, pool3d[j].ewhich_e1);
          connect3d[j].nside_e1 = ipc->get(ne1, pool3d[j].nside_e1);
          connect3d[j].aflag_e1 = ipc->get(ne1, pool3d[j].aflag_e1);
          connect3d[j].fflag_e1 = ipc->get(ne1, pool3d[j].fflag_e1);
          for (int k = 0; k < ne1; k++) {
            connect3d[j].neigh_e1[k] = (tagint) ubuf(buf[m++]).i;
            connect3d[j].ewhich_e1[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].nside_e1[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].aflag_e1[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].fflag_e1[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect3d[j].neigh_e1 = nullptr;
          connect3d[j].ewhich_e1 = nullptr;
          connect3d[j].nside_e1 = nullptr;
          connect3d[j].aflag_e1 = nullptr;
          connect3d[j].fflag_e1 = nullptr;
        }

        if (ne2) {
          connect3d[j].neigh_e2 = tpc->get(ne2, pool3d[j].neigh_e2);
          connect3d[j].ewhich_e2 = ipc->get(ne2, pool3d[j].ewhich_e2);
          connect3d[j].nside_e2 = ipc->get(ne2, pool3d[j].nside_e2);
          connect3d[j].aflag_e2 = ipc->get(ne2, pool3d[j].aflag_e2);
          connect3d[j].fflag_e2 = ipc->get(ne2, pool3d[j].fflag_e2);
          for (int k = 0; k < ne2; k++) {
            connect3d[j].neigh_e2[k] = (tagint) ubuf(buf[m++]).i;
            connect3d[j].ewhich_e2[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].nside_e2[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].aflag_e2[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].fflag_e2[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect3d[j].neigh_e2 = nullptr;
          connect3d[j].ewhich_e2 = nullptr;
          connect3d[j].nside_e2 = nullptr;
          connect3d[j].aflag_e2 = nullptr;
          connect3d[j].fflag_e2 = nullptr;
        }

        if (ne3) {
          connect3d[j].neigh_e3 = tpc->get(ne3, pool3d[j].neigh_e3);
          connect3d[j].ewhich_e3 = ipc->get(ne3, pool3d[j].ewhich_e3);
          connect3d[j].nside_e3 = ipc->get(ne3, pool3d[j].nside_e3);
          connect3d[j].aflag_e3 = ipc->get(ne3, pool3d[j].aflag_e3);
          connect3d[j].fflag_e3 = ipc->get(ne3, pool3d[j].fflag_e3);
          for (int k = 0; k < ne3; k++) {
            connect3d[j].neigh_e3[k] = (tagint) ubuf(buf[m++]).i;
            connect3d[j].ewhich_e3[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].nside_e3[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].aflag_e3[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].fflag_e3[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect3d[j].neigh_e3 = nullptr;
          connect3d[j].ewhich_e3 = nullptr;
          connect3d[j].nside_e3 = nullptr;
          connect3d[j].aflag_e3 = nullptr;
          connect3d[j].fflag_e3 = nullptr;
        }

        nc1 = (int) ubuf(buf[m++]).i;
        nc2 = (int) ubuf(buf[m++]).i;
        nc3 = (int) ubuf(buf[m++]).i;
        connect3d[j].nc1 = nc1;
        connect3d[j].nc2 = nc2;
        connect3d[j].nc3 = nc3;

        connect3d[j].external_cor[0] = (int) ubuf(buf[m++]).i;
        connect3d[j].external_cor[1] = (int) ubuf(buf[m++]).i;
        connect3d[j].external_cor[2] = (int) ubuf(buf[m++]).i;

        if (nc1) {
          connect3d[j].neigh_c1 = tpc->get(nc1, pool3d[j].neigh_c1);
          connect3d[j].cwhich_c1 = ipc->get(nc1, pool3d[j].cwhich_c1);
          connect3d[j].nside_c1 = ipc->get(nc1, pool3d[j].nside_c1);
          connect3d[j].fflag_c1 = ipc->get(nc1, pool3d[j].fflag_c1);
          for (int k = 0; k < nc1; k++) {
            connect3d[j].neigh_c1[k] = (tagint) ubuf(buf[m++]).i;
            connect3d[j].cwhich_c1[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].nside_c1[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].fflag_c1[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect3d[j].neigh_c1 = nullptr;
          connect3d[j].cwhich_c1 = nullptr;
          connect3d[j].nside_c1 = nullptr;
          connect3d[j].fflag_c1 = nullptr;
        }

        if (nc2) {
          connect3d[j].neigh_c2 = tpc->get(nc2, pool3d[j].neigh_c2);
          connect3d[j].cwhich_c2 = ipc->get(nc2, pool3d[j].cwhich_c2);
          connect3d[j].nside_c2 = ipc->get(nc2, pool3d[j].nside_c2);
          connect3d[j].fflag_c2 = ipc->get(nc2, pool3d[j].fflag_c2);
          for (int k = 0; k < nc2; k++) {
            connect3d[j].neigh_c2[k] = (tagint) ubuf(buf[m++]).i;
            connect3d[j].cwhich_c2[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].nside_c2[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].fflag_c2[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect3d[j].neigh_c2 = nullptr;
          connect3d[j].cwhich_c2 = nullptr;
          connect3d[j].nside_c2 = nullptr;
          connect3d[j].fflag_c2 = nullptr;
        }

        if (nc3) {
          connect3d[j].neigh_c3 = tpc->get(nc3, pool3d[j].neigh_c3);
          connect3d[j].cwhich_c3 = ipc->get(nc3, pool3d[j].cwhich_c3);
          connect3d[j].nside_c3 = ipc->get(nc3, pool3d[j].nside_c3);
          connect3d[j].fflag_c3 = ipc->get(nc3, pool3d[j].fflag_c3);
          for (int k = 0; k < nc3; k++) {
            connect3d[j].neigh_c3[k] = (tagint) ubuf(buf[m++]).i;
            connect3d[j].cwhich_c3[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].nside_c3[k] = (int) ubuf(buf[m++]).i;
            connect3d[j].fflag_c3[k] = (int) ubuf(buf[m++]).i;
          }
        } else {
          connect3d[j].neigh_c3 = nullptr;
          connect3d[j].cwhich_c3 = nullptr;
          connect3d[j].nside_c3 = nullptr;
          connect3d[j].fflag_c3 = nullptr;
        }

        connect2atom[j] = i;
        atom2connect[i] = j;
        nghost_connect++;
      }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local connect array for exchange with another proc
------------------------------------------------------------------------- */

int FixSurfaceLocal::pack_exchange(int i, double *buf)
{
  int m = 0;

  if (dimension == 2) {
    int np1, np2;

    if (atom2connect[i] < 0)
      buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      int j = atom2connect[i];

      np1 = connect2d[j].np1;
      np2 = connect2d[j].np2;
      buf[m++] = ubuf(np1).d;
      buf[m++] = ubuf(np2).d;

      buf[m++] = ubuf(connect2d[j].external_pt[0]).d;
      buf[m++] = ubuf(connect2d[j].external_pt[1]).d;

      if (np1)
        for (int k = 0; k < np1; k++) {
          buf[m++] = ubuf(connect2d[j].neigh_p1[k]).d;
          buf[m++] = ubuf(connect2d[j].pwhich_p1[k]).d;
          buf[m++] = ubuf(connect2d[j].nside_p1[k]).d;
          buf[m++] = ubuf(connect2d[j].aflag_p1[k]).d;
          buf[m++] = ubuf(connect2d[j].fflag_p1[k]).d;
        }
      if (np2)
        for (int k = 0; k < np2; k++) {
          buf[m++] = ubuf(connect2d[j].neigh_p2[k]).d;
          buf[m++] = ubuf(connect2d[j].pwhich_p2[k]).d;
          buf[m++] = ubuf(connect2d[j].nside_p2[k]).d;
          buf[m++] = ubuf(connect2d[j].aflag_p2[k]).d;
          buf[m++] = ubuf(connect2d[j].fflag_p2[k]).d;
        }
    }

  } else {
    int ne1, ne2, ne3, nc1, nc2, nc3;

    if (atom2connect[i] < 0)
      buf[m++] = ubuf(0).d;
    else {
      buf[m++] = ubuf(1).d;
      int j = atom2connect[i];

      ne1 = connect3d[j].ne1;
      ne2 = connect3d[j].ne2;
      ne3 = connect3d[j].ne3;
      buf[m++] = ubuf(ne1).d;
      buf[m++] = ubuf(ne2).d;
      buf[m++] = ubuf(ne3).d;

      buf[m++] = ubuf(connect3d[j].external_edge[0]).d;
      buf[m++] = ubuf(connect3d[j].external_edge[1]).d;
      buf[m++] = ubuf(connect3d[j].external_edge[2]).d;

      if (ne1)
        for (int k = 0; k < ne1; k++) {
          buf[m++] = ubuf(connect3d[j].neigh_e1[k]).d;
          buf[m++] = ubuf(connect3d[j].ewhich_e1[k]).d;
          buf[m++] = ubuf(connect3d[j].nside_e1[k]).d;
          buf[m++] = ubuf(connect3d[j].aflag_e1[k]).d;
          buf[m++] = ubuf(connect3d[j].fflag_e1[k]).d;
        }
      if (ne2)
        for (int k = 0; k < ne2; k++) {
          buf[m++] = ubuf(connect3d[j].neigh_e2[k]).d;
          buf[m++] = ubuf(connect3d[j].ewhich_e2[k]).d;
          buf[m++] = ubuf(connect3d[j].nside_e2[k]).d;
          buf[m++] = ubuf(connect3d[j].aflag_e2[k]).d;
          buf[m++] = ubuf(connect3d[j].fflag_e2[k]).d;
        }
      if (ne3)
        for (int k = 0; k < ne3; k++) {
          buf[m++] = ubuf(connect3d[j].neigh_e3[k]).d;
          buf[m++] = ubuf(connect3d[j].ewhich_e3[k]).d;
          buf[m++] = ubuf(connect3d[j].nside_e3[k]).d;
          buf[m++] = ubuf(connect3d[j].aflag_e3[k]).d;
          buf[m++] = ubuf(connect3d[j].fflag_e3[k]).d;
        }

      nc1 = connect3d[j].nc1;
      nc2 = connect3d[j].nc2;
      nc3 = connect3d[j].nc3;
      buf[m++] = ubuf(nc1).d;
      buf[m++] = ubuf(nc2).d;
      buf[m++] = ubuf(nc3).d;

      buf[m++] = ubuf(connect3d[j].external_cor[0]).d;
      buf[m++] = ubuf(connect3d[j].external_cor[1]).d;
      buf[m++] = ubuf(connect3d[j].external_cor[2]).d;

      if (nc1)
        for (int k = 0; k < nc1; k++) {
          buf[m++] = ubuf(connect3d[j].neigh_c1[k]).d;
          buf[m++] = ubuf(connect3d[j].cwhich_c1[k]).d;
          buf[m++] = ubuf(connect3d[j].nside_c1[k]).d;
          buf[m++] = ubuf(connect3d[j].fflag_c1[k]).d;
        }
      if (nc2)
        for (int k = 0; k < nc2; k++) {
          buf[m++] = ubuf(connect3d[j].neigh_c2[k]).d;
          buf[m++] = ubuf(connect3d[j].cwhich_c2[k]).d;
          buf[m++] = ubuf(connect3d[j].nside_c2[k]).d;
          buf[m++] = ubuf(connect3d[j].fflag_c2[k]).d;
        }
      if (nc3)
        for (int k = 0; k < nc3; k++) {
          buf[m++] = ubuf(connect3d[j].neigh_c3[k]).d;
          buf[m++] = ubuf(connect3d[j].cwhich_c3[k]).d;
          buf[m++] = ubuf(connect3d[j].nside_c3[k]).d;
          buf[m++] = ubuf(connect3d[j].fflag_c3[k]).d;
        }
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local connect array from exchange with another proc
------------------------------------------------------------------------- */

int FixSurfaceLocal::unpack_exchange(int nlocal, double *buf)
{
  int flag;
  int m = 0;

  if (dimension == 2) {
    int np1, np2;

    flag = (int) ubuf(buf[m++]).i;
    if (flag == 0)
      atom2connect[nlocal] = -1;
    else {
      if (nlocal_connect == nmax_connect) grow_connect();

      np1 = (int) ubuf(buf[m++]).i;
      np2 = (int) ubuf(buf[m++]).i;
      connect2d[nlocal_connect].np1 = np1;
      connect2d[nlocal_connect].np2 = np2;

      connect2d[nlocal_connect].external_pt[0] = (int) ubuf(buf[m++]).i;
      connect2d[nlocal_connect].external_pt[1] = (int) ubuf(buf[m++]).i;

      if (np1) {
        connect2d[nlocal_connect].neigh_p1 = tpc->get(np1, pool2d[nlocal_connect].neigh_p1);
        connect2d[nlocal_connect].pwhich_p1 = ipc->get(np1, pool2d[nlocal_connect].pwhich_p1);
        connect2d[nlocal_connect].nside_p1 = ipc->get(np1, pool2d[nlocal_connect].nside_p1);
        connect2d[nlocal_connect].aflag_p1 = ipc->get(np1, pool2d[nlocal_connect].aflag_p1);
        connect2d[nlocal_connect].fflag_p1 = ipc->get(np1, pool2d[nlocal_connect].fflag_p1);
        for (int k = 0; k < np1; k++) {
          connect2d[nlocal_connect].neigh_p1[k] = (tagint) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].pwhich_p1[k] = (int) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].nside_p1[k] = (int) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].aflag_p1[k] = (int) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].fflag_p1[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect2d[nlocal_connect].neigh_p1 = nullptr;
        connect2d[nlocal_connect].pwhich_p1 = nullptr;
        connect2d[nlocal_connect].nside_p1 = nullptr;
        connect2d[nlocal_connect].aflag_p1 = nullptr;
        connect2d[nlocal_connect].fflag_p1 = nullptr;
      }

      if (np2) {
        connect2d[nlocal_connect].neigh_p2 = tpc->get(np2, pool2d[nlocal_connect].neigh_p2);
        connect2d[nlocal_connect].pwhich_p2 = ipc->get(np2, pool2d[nlocal_connect].pwhich_p2);
        connect2d[nlocal_connect].nside_p2 = ipc->get(np2, pool2d[nlocal_connect].nside_p2);
        connect2d[nlocal_connect].aflag_p2 = ipc->get(np2, pool2d[nlocal_connect].aflag_p2);
        connect2d[nlocal_connect].fflag_p2 = ipc->get(np2, pool2d[nlocal_connect].fflag_p2);
        for (int k = 0; k < np2; k++) {
          connect2d[nlocal_connect].neigh_p2[k] = (tagint) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].pwhich_p2[k] = (int) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].nside_p2[k] = (int) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].aflag_p2[k] = (int) ubuf(buf[m++]).i;
          connect2d[nlocal_connect].fflag_p2[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect2d[nlocal_connect].neigh_p2 = nullptr;
        connect2d[nlocal_connect].pwhich_p2 = nullptr;
        connect2d[nlocal_connect].nside_p2 = nullptr;
        connect2d[nlocal_connect].aflag_p2 = nullptr;
        connect2d[nlocal_connect].fflag_p2 = nullptr;
      }

      connect2atom[nlocal_connect] = nlocal;
      atom2connect[nlocal] = nlocal_connect++;
    }

  } else {
    int ne1, ne2, ne3, nc1, nc2, nc3;

    flag = (int) ubuf(buf[m++]).i;
    if (flag == 0)
      atom2connect[nlocal] = -1;
    else {
      if (nlocal_connect == nmax_connect) grow_connect();

      ne1 = (int) ubuf(buf[m++]).i;
      ne2 = (int) ubuf(buf[m++]).i;
      ne3 = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].ne1 = ne1;
      connect3d[nlocal_connect].ne2 = ne2;
      connect3d[nlocal_connect].ne3 = ne3;

      connect3d[nlocal_connect].external_edge[0] = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].external_edge[1] = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].external_edge[2] = (int) ubuf(buf[m++]).i;

      if (ne1) {
        connect3d[nlocal_connect].neigh_e1 = tpc->get(ne1, pool3d[nlocal_connect].neigh_e1);
        connect3d[nlocal_connect].ewhich_e1 = ipc->get(ne1, pool3d[nlocal_connect].ewhich_e1);
        connect3d[nlocal_connect].nside_e1 = ipc->get(ne1, pool3d[nlocal_connect].nside_e1);
        connect3d[nlocal_connect].aflag_e1 = ipc->get(ne1, pool3d[nlocal_connect].aflag_e1);
        connect3d[nlocal_connect].fflag_e1 = ipc->get(ne1, pool3d[nlocal_connect].fflag_e1);
        for (int k = 0; k < ne1; k++) {
          connect3d[nlocal_connect].neigh_e1[k] = (tagint) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].ewhich_e1[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].nside_e1[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].aflag_e1[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].fflag_e1[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect3d[nlocal_connect].neigh_e1 = nullptr;
        connect3d[nlocal_connect].ewhich_e1 = nullptr;
        connect3d[nlocal_connect].nside_e1 = nullptr;
        connect3d[nlocal_connect].aflag_e1 = nullptr;
        connect3d[nlocal_connect].fflag_e1 = nullptr;
      }

      if (ne2) {
        connect3d[nlocal_connect].neigh_e2 = tpc->get(ne2, pool3d[nlocal_connect].neigh_e2);
        connect3d[nlocal_connect].ewhich_e2 = ipc->get(ne2, pool3d[nlocal_connect].ewhich_e2);
        connect3d[nlocal_connect].nside_e2 = ipc->get(ne2, pool3d[nlocal_connect].nside_e2);
        connect3d[nlocal_connect].aflag_e2 = ipc->get(ne2, pool3d[nlocal_connect].aflag_e2);
        connect3d[nlocal_connect].fflag_e2 = ipc->get(ne2, pool3d[nlocal_connect].fflag_e2);
        for (int k = 0; k < ne2; k++) {
          connect3d[nlocal_connect].neigh_e2[k] = (tagint) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].ewhich_e2[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].nside_e2[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].aflag_e2[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].fflag_e2[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect3d[nlocal_connect].neigh_e2 = nullptr;
        connect3d[nlocal_connect].ewhich_e2 = nullptr;
        connect3d[nlocal_connect].nside_e2 = nullptr;
        connect3d[nlocal_connect].aflag_e2 = nullptr;
        connect3d[nlocal_connect].fflag_e2 = nullptr;
      }

      if (ne3) {
        connect3d[nlocal_connect].neigh_e3 = tpc->get(ne3, pool3d[nlocal_connect].neigh_e3);
        connect3d[nlocal_connect].ewhich_e3 = ipc->get(ne3, pool3d[nlocal_connect].ewhich_e3);
        connect3d[nlocal_connect].nside_e3 = ipc->get(ne3, pool3d[nlocal_connect].nside_e3);
        connect3d[nlocal_connect].aflag_e3 = ipc->get(ne3, pool3d[nlocal_connect].aflag_e3);
        connect3d[nlocal_connect].fflag_e3 = ipc->get(ne3, pool3d[nlocal_connect].fflag_e3);
        for (int k = 0; k < ne3; k++) {
          connect3d[nlocal_connect].neigh_e3[k] = (tagint) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].ewhich_e3[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].nside_e3[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].aflag_e3[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].fflag_e3[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect3d[nlocal_connect].neigh_e3 = nullptr;
        connect3d[nlocal_connect].ewhich_e3 = nullptr;
        connect3d[nlocal_connect].nside_e3 = nullptr;
        connect3d[nlocal_connect].aflag_e3 = nullptr;
        connect3d[nlocal_connect].fflag_e3 = nullptr;
      }

      nc1 = (int) ubuf(buf[m++]).i;
      nc2 = (int) ubuf(buf[m++]).i;
      nc3 = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].nc1 = nc1;
      connect3d[nlocal_connect].nc2 = nc2;
      connect3d[nlocal_connect].nc3 = nc3;

      connect3d[nlocal_connect].external_cor[0] = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].external_cor[1] = (int) ubuf(buf[m++]).i;
      connect3d[nlocal_connect].external_cor[2] = (int) ubuf(buf[m++]).i;

      if (nc1) {
        connect3d[nlocal_connect].neigh_c1 = tpc->get(nc1, pool3d[nlocal_connect].neigh_c1);
        connect3d[nlocal_connect].cwhich_c1 = ipc->get(nc1, pool3d[nlocal_connect].cwhich_c1);
        connect3d[nlocal_connect].nside_c1 = ipc->get(nc1, pool3d[nlocal_connect].nside_c1);
        connect3d[nlocal_connect].fflag_c1 = ipc->get(nc1, pool3d[nlocal_connect].fflag_c1);
        for (int k = 0; k < nc1; k++) {
          connect3d[nlocal_connect].neigh_c1[k] = (tagint) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].cwhich_c1[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].nside_c1[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].fflag_c1[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect3d[nlocal_connect].neigh_c1 = nullptr;
        connect3d[nlocal_connect].cwhich_c1 = nullptr;
        connect3d[nlocal_connect].nside_c1 = nullptr;
        connect3d[nlocal_connect].fflag_c1 = nullptr;
      }

      if (nc2) {
        connect3d[nlocal_connect].neigh_c2 = tpc->get(nc2, pool3d[nlocal_connect].neigh_c2);
        connect3d[nlocal_connect].cwhich_c2 = ipc->get(nc2, pool3d[nlocal_connect].cwhich_c2);
        connect3d[nlocal_connect].nside_c2 = ipc->get(nc2, pool3d[nlocal_connect].nside_c2);
        connect3d[nlocal_connect].fflag_c2 = ipc->get(nc2, pool3d[nlocal_connect].fflag_c2);
        for (int k = 0; k < nc2; k++) {
          connect3d[nlocal_connect].neigh_c2[k] = (tagint) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].cwhich_c2[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].nside_c2[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].fflag_c2[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect3d[nlocal_connect].neigh_c2 = nullptr;
        connect3d[nlocal_connect].cwhich_c2 = nullptr;
        connect3d[nlocal_connect].nside_c2 = nullptr;
        connect3d[nlocal_connect].fflag_c2 = nullptr;
      }

      if (nc3) {
        connect3d[nlocal_connect].neigh_c3 = tpc->get(nc3, pool3d[nlocal_connect].neigh_c3);
        connect3d[nlocal_connect].cwhich_c3 = ipc->get(nc3, pool3d[nlocal_connect].cwhich_c3);
        connect3d[nlocal_connect].nside_c3 = ipc->get(nc3, pool3d[nlocal_connect].nside_c3);
        connect3d[nlocal_connect].fflag_c3 = ipc->get(nc3, pool3d[nlocal_connect].fflag_c3);
        for (int k = 0; k < nc3; k++) {
          connect3d[nlocal_connect].neigh_c3[k] = (tagint) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].cwhich_c3[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].nside_c3[k] = (int) ubuf(buf[m++]).i;
          connect3d[nlocal_connect].fflag_c3[k] = (int) ubuf(buf[m++]).i;
        }
      } else {
        connect3d[nlocal_connect].neigh_c3 = nullptr;
        connect3d[nlocal_connect].cwhich_c3 = nullptr;
        connect3d[nlocal_connect].nside_c3 = nullptr;
        connect3d[nlocal_connect].fflag_c3 = nullptr;
      }

      connect2atom[nlocal_connect] = nlocal;
      atom2connect[nlocal] = nlocal_connect++;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   one-time forward comm of connectivity info from owned to ghost particles
   same pack/unpack operations as comm->borders() uses
------------------------------------------------------------------------- */

int FixSurfaceLocal::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                       int * /*pbc*/)
{
  return pack_border(n, list, buf);
}

/* ---------------------------------------------------------------------- */

void FixSurfaceLocal::unpack_forward_comm(int n, int first, double *buf)
{
  unpack_border(n, first, buf);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSurfaceLocal::memory_usage()
{
  double bytes = 0.0;

  if (dimension == 2) {
    bytes = nmax_connect * sizeof(Connect2d);
    bytes = nmax_connect * sizeof(Pool2d);
  } else {
    bytes = nmax_connect * sizeof(Connect3d);
    bytes = nmax_connect * sizeof(Pool3d);
  }

  bytes += atom->nmax * sizeof(int);      // atom2connect vector
  bytes += nmax_connect * sizeof(int);    // connect2atom vector

  return bytes;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for distributed connectivity builds in 2d or 3d
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create and initialize Connect2d info for owned lines
     only np1,np2 and neigh_p1,neigh_p2
     also atom2connect,connect2atom
   this must be done with INEXACT point matching
     done via Rvous algorithm on set of slightly overlapping bins
   inexact b/c once datafile is read, lines have only a center point, length, and theta
     thus same point calculated by 2 lines may be epsilon different
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity2d_local()
{
  int m, n;

  avec_line = (AtomVecLine *) atom->style_match("line");

  // calculate epssq = square of EPSILON fraction of minimum line length
  // for use in point_match()

  if (epssq < 0.0) epsilon_calculate();

  // nline = count of owned lines

  int *line = atom->line;
  int nlocal = atom->nlocal;

  int nline = 0;
  for (int i = 0; i < nlocal; i++)
    if (line[i] >= 0) nline++;

  // allocate connection info for owned lines
  // initialize atom2connect for both particles and lines

  nlocal_connect = nmax_connect = nline;
  nghost_connect = 0;
  grow_connect();

  for (int i = 0; i < nlocal; i++) {
    atom2connect[i] = line[i];
    if (line[i] < 0) continue;
    int j = line[i];
    connect2atom[j] = i;
  }

  // calculate current endpts of owned lines

  double **endpts;
  memory->create(endpts, nline, 4, "surface/local:endpts");
  calculate_endpts(endpts);

  // compute min/max extent of my line endpts
  // MPI_Allreduce for bbox of all lines

  double mylo[2], myhi[2];

  mylo[0] = mylo[1] = BIG;
  myhi[0] = myhi[1] = -BIG;

  m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    mylo[0] = MIN(mylo[0], endpts[m][0]);
    myhi[0] = MAX(myhi[0], endpts[m][0]);
    mylo[1] = MIN(mylo[1], endpts[m][1]);
    myhi[1] = MAX(myhi[1], endpts[m][1]);
    mylo[0] = MIN(mylo[0], endpts[m][2]);
    myhi[0] = MAX(myhi[0], endpts[m][2]);
    mylo[1] = MIN(mylo[1], endpts[m][3]);
    myhi[1] = MAX(myhi[1], endpts[m][3]);
    m++;
  }

  MPI_Allreduce(mylo, bboxlo, 2, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(myhi, bboxhi, 2, MPI_DOUBLE, MPI_MAX, world);

  // add 2*EPS to all 4 edges of bbox

  bboxlo[0] -= 2.0 * eps;
  bboxlo[1] -= 2.0 * eps;
  bboxhi[0] += 2.0 * eps;
  bboxhi[1] += 2.0 * eps;

  // conceptual binning of bbox by up to NBIN x NBIN bins
  // ensure bin size is not <= 2*EPS so that pt +/- EPS cannot overlap > 4 bins
  // nbin xy = # of bins in each dim

  nbinx = static_cast<int>((bboxhi[0] - bboxlo[0]) / (4.0 * eps));
  nbinx = MIN(nbinx, NBIN);
  nbinx = MAX(nbinx, 1);
  nbiny = static_cast<int>((bboxhi[1] - bboxlo[1]) / (4.0 * eps));
  nbiny = MIN(nbiny, NBIN);
  nbiny = MAX(nbiny, 1);

  nbins = nbinx * nbiny;

  invbinx = nbinx / (bboxhi[0] - bboxlo[0]);
  invbiny = nbiny / (bboxhi[1] - bboxlo[1]);

  // inbuf = list of datums to send to procs in Rvous decomposition of bins
  // every Pth bin is assigned to each proc
  // use overlap_bins_2d() to find all bins a line endpt overlaps within EPS
  // allows for matching pts in Rvous decomp which are up to EPS apart

  int me = comm->me;
  int nprocs = comm->nprocs;

  tagint *molecule = atom->molecule;
  tagint *tag = atom->tag;

  int *proclist = nullptr;
  InRvous *inbuf = nullptr;
  int ncount = 0;
  int maxcount = 0;

  int indices[4];

  m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;

    for (int ipoint = 0; ipoint < 2; ipoint++) {
      n = overlap_bins_2d(&endpts[m][2 * ipoint], eps, indices);

      if (ncount + n > maxcount) {
        maxcount += DELTA_RVOUS;
        memory->grow(proclist, maxcount, "fix/surface/local:proclist");
        inbuf =
            (InRvous *) memory->srealloc(inbuf, maxcount * sizeof(InRvous), "surface/local:inbuf");
      }

      for (int k = 0; k < n; k++) {
        proclist[ncount] = indices[k] % nprocs;
        inbuf[ncount].proc = me;
        inbuf[ncount].ibin = indices[k];
        inbuf[ncount].ilocal = i;
        inbuf[ncount].ipoint = ipoint;
        inbuf[ncount].mol = molecule[i];
        inbuf[ncount].atomID = tag[i];
        inbuf[ncount].x[0] = endpts[m][2 * ipoint];
        inbuf[ncount].x[1] = endpts[m][2 * ipoint + 1];
        inbuf[ncount].x[2] = 0.0;
        ncount++;
      }
    }

    m++;
  }

  memory->destroy(endpts);

  // perform rendezvous operation
  // each proc owns every Pth bin
  // receives all points in those bins

  char *buf;
  int nreturn = comm->rendezvous(RVOUS, ncount, (char *) inbuf, sizeof(InRvous), 0, proclist,
                                 point_match, 0, buf, sizeof(OutRvous), (void *) this);
  auto outbuf = (OutRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // loop over received Rvous datums
  // p12_counts = # of connections for each point on my lines
  // this will overcount (potentially by 4x) due to bins overlapping by EPS
  // datums do NOT include self connection

  int *p1_counts, *p2_counts;
  memory->create(p1_counts, nlocal_connect, "surface/local:p1_counts");
  memory->create(p2_counts, nlocal_connect, "surface/local:p2_counts");

  int ilocal, iline;

  for (int i = 0; i < nlocal_connect; i++) p1_counts[i] = p2_counts[i] = 0;

  for (int i = 0; i < nreturn; i++) {
    ilocal = outbuf[i].ilocal;
    iline = line[ilocal];
    if (outbuf[i].ipoint == 0)
      p1_counts[iline]++;
    else
      p2_counts[iline]++;
  }

  // allocate ragged tneigh12 vectors using p12_counts
  // this will overallocate because of bins overlapping by EPS
  //   will reallocate below when know exact size

  tagint **tneigh1, **tneigh2;
  memory->create_ragged(tneigh1, nlocal_connect, p1_counts, "surface/local:tneigh1");
  memory->create_ragged(tneigh2, nlocal_connect, p2_counts, "surface/local:tneigh2");

  // loop over received Rvous datums
  // add each atomID to tneigh12 but only if not already in the list
  // recalculate exact p12_counts

  for (int i = 0; i < nlocal_connect; i++) p1_counts[i] = p2_counts[i] = 0;

  int np, j;
  tagint atomID;
  tagint *neigh;

  for (int i = 0; i < nreturn; i++) {
    ilocal = outbuf[i].ilocal;
    iline = line[ilocal];
    if (outbuf[i].ipoint == 0) {
      atomID = outbuf[i].atomID;
      np = p1_counts[iline];
      neigh = tneigh1[iline];
      for (j = 0; j < np; j++)
        if (neigh[j] == atomID) break;
      if (j == np) {
        neigh[np] = atomID;
        p1_counts[iline]++;
      }
    } else {
      atomID = outbuf[i].atomID;
      np = p2_counts[iline];
      neigh = tneigh2[iline];
      for (j = 0; j < np; j++)
        if (neigh[j] == atomID) break;
      if (j == np) {
        neigh[np] = atomID;
        p2_counts[iline]++;
      }
    }
  }

  // set np1,np2 via exact p12_counts
  // likewise allocate all vectors within Connect2d
  // use ragged tneigh12 to set neigh_p12 within Connect2d
  // initialize other vectors to 0 for first-time borders comm
  //   they will be set correctly in connectivity2d_complete()

  for (int i = 0; i < nlocal_connect; i++) {
    connect2d[i].np1 = p1_counts[i];
    if (connect2d[i].np1) {
      connect2d[i].neigh_p1 = tpc->get(connect2d[i].np1, pool2d[i].neigh_p1);
      connect2d[i].pwhich_p1 = ipc->get(connect2d[i].np1, pool2d[i].pwhich_p1);
      connect2d[i].nside_p1 = ipc->get(connect2d[i].np1, pool2d[i].nside_p1);
      connect2d[i].aflag_p1 = ipc->get(connect2d[i].np1, pool2d[i].aflag_p1);
      connect2d[i].fflag_p1 = ipc->get(connect2d[i].np1, pool2d[i].fflag_p1);
      for (int j = 0; j < connect2d[i].np1; j++) {
        connect2d[i].neigh_p1[j] = tneigh1[i][j];
        connect2d[i].pwhich_p1[j] = 0;
        connect2d[i].nside_p1[j] = 0;
        connect2d[i].aflag_p1[j] = 0;
        connect2d[i].fflag_p1[j] = 0;
      }
    } else {
      connect2d[i].neigh_p1 = nullptr;
      connect2d[i].pwhich_p1 = nullptr;
      connect2d[i].nside_p1 = nullptr;
      connect2d[i].aflag_p1 = nullptr;
      connect2d[i].fflag_p1 = nullptr;
    }

    connect2d[i].np2 = p2_counts[i];
    if (connect2d[i].np2) {
      connect2d[i].neigh_p2 = tpc->get(connect2d[i].np2, pool2d[i].neigh_p2);
      connect2d[i].pwhich_p2 = ipc->get(connect2d[i].np2, pool2d[i].pwhich_p2);
      connect2d[i].nside_p2 = ipc->get(connect2d[i].np2, pool2d[i].nside_p2);
      connect2d[i].aflag_p2 = ipc->get(connect2d[i].np2, pool2d[i].aflag_p2);
      connect2d[i].fflag_p2 = ipc->get(connect2d[i].np2, pool2d[i].fflag_p2);
      for (int j = 0; j < connect2d[i].np2; j++) {
        connect2d[i].neigh_p2[j] = tneigh2[i][j];
        connect2d[i].pwhich_p2[j] = 0;
        connect2d[i].nside_p2[j] = 0;
        connect2d[i].aflag_p2[j] = 0;
        connect2d[i].fflag_p2[j] = 0;
      }
    } else {
      connect2d[i].neigh_p2 = nullptr;
      connect2d[i].pwhich_p2 = nullptr;
      connect2d[i].nside_p2 = nullptr;
      connect2d[i].aflag_p2 = nullptr;
      connect2d[i].fflag_p2 = nullptr;
    }
  }

  // clean up

  memory->sfree(outbuf);
  memory->destroy(p1_counts);
  memory->destroy(p2_counts);
  memory->destroy(tneigh1);
  memory->destroy(tneigh2);
}

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for owned tris
     only ne1,ne2,ne3 and neigh_e1,neigh_e2,neigh_e3
     only nc1,nc2,nc3 and neigh_c1,neigh_c2,neigh_c3
     also atom2connect,connect2atom
   this must be done with INEXACT point matching
     done via Rvous algorithm on set of slightly overlapping bins
   inexact b/c once datafile is read, tris have a center point, quaternion,
     and body-frame corner point displacements
     thus same point calculated by 2 tris may be epsilon different
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity3d_local()
{
  int k, m, n;

  avec_tri = (AtomVecTri *) atom->style_match("tri");

  // epssq = square of EPSILON fraction of minimum tri diameter
  // for use in point_match()

  if (epssq < 0.0) epsilon_calculate();

  // ntri = count of owned triangles

  int *tri = atom->tri;
  int nlocal = atom->nlocal;

  int ntri = 0;
  for (int i = 0; i < nlocal; i++)
    if (tri[i] >= 0) ntri++;

  // allocate connection info for owned triangles
  // initialize atom2connect for both particles and triangles

  nlocal_connect = nmax_connect = ntri;
  nghost_connect = 0;
  grow_connect();

  for (int i = 0; i < nlocal; i++) {
    atom2connect[i] = tri[i];
    if (tri[i] < 0) continue;
    int j = tri[i];
    connect2atom[j] = i;
  }

  // calculate current corners of owned tris

  double **corners;
  memory->create(corners, ntri, 9, "surface/local:corners");
  calculate_corners(corners);

  // compute min/max extent of my tri corners
  // MPI_Allreduce for bbox of all lines

  double mylo[3], myhi[3];

  mylo[0] = mylo[1] = mylo[2] = BIG;
  myhi[0] = myhi[1] = myhi[2] = -BIG;

  m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    int k = 0;
    for (int j = 0; j < 3; j++) {
      mylo[0] = MIN(mylo[0], corners[m][k]);
      myhi[0] = MAX(myhi[0], corners[m][k]);
      mylo[1] = MIN(mylo[1], corners[m][k + 1]);
      myhi[1] = MAX(myhi[1], corners[m][k + 1]);
      mylo[2] = MIN(mylo[2], corners[m][k + 2]);
      myhi[2] = MAX(myhi[2], corners[m][k + 2]);
      k += 3;
    }
    m++;
  }

  MPI_Allreduce(mylo, bboxlo, 3, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(myhi, bboxhi, 3, MPI_DOUBLE, MPI_MAX, world);

  // add 2*EPS to all 4 faces of bbox

  bboxlo[0] -= 2.0 * eps;
  bboxlo[1] -= 2.0 * eps;
  bboxlo[2] -= 2.0 * eps;
  bboxhi[0] += 2.0 * eps;
  bboxhi[1] += 2.0 * eps;
  bboxhi[2] += 2.0 * eps;

  // conceptual binning of bbox by up to NBIN x NBIN x NBIN bins
  // ensure bin size is not <= 2*EPS so that pt +/- EPS cannot overlap > 8 bins
  // nbin xyz = # of bins in each dim

  nbinx = static_cast<int>((bboxhi[0] - bboxlo[0]) / (4.0 * eps));
  nbinx = MIN(nbinx, NBIN);
  nbinx = MAX(nbinx, 1);
  nbiny = static_cast<int>((bboxhi[1] - bboxlo[1]) / (4.0 * eps));
  nbiny = MIN(nbiny, NBIN);
  nbiny = MAX(nbiny, 1);
  nbinz = static_cast<int>((bboxhi[2] - bboxlo[2]) / (4.0 * eps));
  nbinz = MIN(nbinz, NBIN);
  nbinz = MAX(nbinz, 1);

  nbins = nbinx * nbiny * nbinz;

  invbinx = nbinx / (bboxhi[0] - bboxlo[0]);
  invbiny = nbiny / (bboxhi[1] - bboxlo[1]);
  invbinz = nbinz / (bboxhi[2] - bboxlo[2]);

  // inbuf = list of datums to send to procs in Rvous decomposition of bins
  // every Pth bin is assigned to each proc
  // use overlap_bins_3d() to find all bins a tri corner pt overlaps within EPS
  // allows for matching pts in Rvous decomp which are up to EPS apart

  int me = comm->me;
  int nprocs = comm->nprocs;

  tagint *molecule = atom->molecule;
  tagint *tag = atom->tag;

  int *proclist = nullptr;
  InRvous *inbuf = nullptr;
  int ncount = 0;
  int maxcount = 0;

  int indices[8];

  m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;

    for (int ipoint = 0; ipoint < 3; ipoint++) {
      n = overlap_bins_3d(&corners[m][3 * ipoint], eps, indices);

      if (ncount + n > maxcount) {
        maxcount += DELTA_RVOUS;
        memory->grow(proclist, maxcount, "fix/surface/local:proclist");
        inbuf =
            (InRvous *) memory->srealloc(inbuf, maxcount * sizeof(InRvous), "surface/local:inbuf");
      }

      for (int k = 0; k < n; k++) {
        proclist[ncount] = indices[k] % nprocs;
        inbuf[ncount].proc = me;
        inbuf[ncount].ibin = indices[k];
        inbuf[ncount].ilocal = i;
        inbuf[ncount].ipoint = ipoint;
        inbuf[ncount].mol = molecule[i];
        inbuf[ncount].atomID = tag[i];
        inbuf[ncount].x[0] = corners[m][3 * ipoint];
        inbuf[ncount].x[1] = corners[m][3 * ipoint + 1];
        inbuf[ncount].x[2] = corners[m][3 * ipoint + 2];
        ncount++;
      }
    }

    m++;
  }

  memory->destroy(corners);

  // perform rendezvous operation
  // each proc owns every Pth bin
  // receives all points in those bins

  char *buf;
  int nreturn = comm->rendezvous(RVOUS, ncount, (char *) inbuf, sizeof(InRvous), 0, proclist,
                                 point_match, 0, buf, sizeof(OutRvous), (void *) this);
  auto outbuf = (OutRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // loop over received Rvous datums
  // n1/n2/n3_counts = # of connections for each corner point on my tris
  // this will overcount (potentially by 8x) due to bins overlapping by EPS
  // datums do NOT include self connection

  int *n1_counts, *n2_counts, *n3_counts;
  memory->create(n1_counts, nlocal_connect, "surface/local:n1_counts");
  memory->create(n2_counts, nlocal_connect, "surface/local:n2_counts");
  memory->create(n3_counts, nlocal_connect, "surface/local:n3_counts");

  int ilocal, iconnect;

  for (int i = 0; i < nlocal_connect; i++) n1_counts[i] = n2_counts[i] = n3_counts[i] = 0;

  for (int i = 0; i < nreturn; i++) {
    ilocal = outbuf[i].ilocal;
    iconnect = atom2connect[ilocal];
    if (outbuf[i].ipoint == 0)
      n1_counts[iconnect]++;
    else if (outbuf[i].ipoint == 1)
      n2_counts[iconnect]++;
    else
      n3_counts[iconnect]++;
  }

  // allocate ragged tneigh123 vectors using n1/n2/n3_counts
  // this will overallocate because of bins overlapping by EPS
  //   will reallocate below when know exact size of edge and corner neighs

  tagint **tneigh1, **tneigh2, **tneigh3;
  memory->create_ragged(tneigh1, nlocal_connect, n1_counts, "surface/local:tneigh1");
  memory->create_ragged(tneigh2, nlocal_connect, n2_counts, "surface/local:tneigh2");
  memory->create_ragged(tneigh3, nlocal_connect, n3_counts, "surface/local:tneigh3");

  // loop over received Rvous datums
  // add each atomID to tneigh123 but only if not already in the list
  // recalculate exact n1/n2/n3_counts

  for (int i = 0; i < nlocal_connect; i++) n1_counts[i] = n2_counts[i] = n3_counts[i] = 0;

  int j, np;
  tagint atomID;
  tagint *neigh;

  for (int i = 0; i < nreturn; i++) {
    ilocal = outbuf[i].ilocal;
    iconnect = atom2connect[ilocal];
    if (outbuf[i].ipoint == 0) {
      atomID = outbuf[i].atomID;
      np = n1_counts[iconnect];
      neigh = tneigh1[iconnect];
      for (j = 0; j < np; j++)
        if (neigh[j] == atomID) break;
      if (j == np) {
        neigh[np] = atomID;
        n1_counts[iconnect]++;
      }
    } else if (outbuf[i].ipoint == 1) {
      atomID = outbuf[i].atomID;
      np = n2_counts[iconnect];
      neigh = tneigh2[iconnect];
      for (int j = 0; j < np; j++)
        if (neigh[j] == atomID) break;
      if (j == np) {
        neigh[np] = atomID;
        n2_counts[iconnect]++;
      }
    } else if (outbuf[i].ipoint == 2) {
      atomID = outbuf[i].atomID;
      np = n3_counts[iconnect];
      neigh = tneigh3[iconnect];
      for (int j = 0; j < np; j++)
        if (neigh[j] == atomID) break;
      if (j == np) {
        neigh[np] = atomID;
        n3_counts[iconnect]++;
      }
    }
  }

  memory->sfree(outbuf);

  // now have exact list of all neighbor tris of each corner point
  // each neighbor is either an edge neighbor or corner neighbor
  // edge neighbor = appears in neigh lists of 2 adjacent corners
  // corner neighbor = appears only in neigh list of a single corner

  int n1, n2;
  tagint *neigh1, *neigh2;

  for (int i = 0; i < nlocal_connect; i++) {
    connect3d[i].ne1 = connect3d[i].ne2 = connect3d[i].ne3 = 0;
    connect3d[i].nc1 = connect3d[i].nc2 = connect3d[i].nc3 = 0;
  }

  for (int i = 0; i < nlocal_connect; i++) {

    // count edge neighbors first

    n1 = n1_counts[i];
    n2 = n2_counts[i];
    neigh1 = tneigh1[i];
    neigh2 = tneigh2[i];
    for (int j = 0; j < n1; j++) {
      for (int k = 0; k < n2; k++) {
        if (neigh2[k] == neigh1[j]) {
          connect3d[i].ne1++;
          break;
        }
      }
    }

    n1 = n2_counts[i];
    n2 = n3_counts[i];
    neigh1 = tneigh2[i];
    neigh2 = tneigh3[i];
    for (int j = 0; j < n1; j++) {
      for (int k = 0; k < n2; k++) {
        if (neigh2[k] == neigh1[j]) {
          connect3d[i].ne2++;
          break;
        }
      }
    }

    n1 = n3_counts[i];
    n2 = n1_counts[i];
    neigh1 = tneigh3[i];
    neigh2 = tneigh1[i];
    for (int j = 0; j < n1; j++) {
      for (int k = 0; k < n2; k++) {
        if (neigh2[k] == neigh1[j]) {
          connect3d[i].ne3++;
          break;
        }
      }
    }

    // corner neighbor count = all neighbors minus 2 sets of edge neighbors

    connect3d[i].nc1 = n1_counts[i] - connect3d[i].ne1 - connect3d[i].ne3;
    connect3d[i].nc2 = n2_counts[i] - connect3d[i].ne2 - connect3d[i].ne1;
    connect3d[i].nc3 = n3_counts[i] - connect3d[i].ne3 - connect3d[i].ne2;
  }

  // use exact n1/n2/n3_counts to allocate all vectors within Connect3d
  // neigh_e123 and neigh_c123 will be initialized below
  // initialize other vectors to 0 for first-time borders comm
  //   they will be set correctly in connectivity3d_complete()

  for (int i = 0; i < nlocal_connect; i++) {
    if (connect3d[i].ne1) {
      connect3d[i].neigh_e1 = tpc->get(connect3d[i].ne1, pool3d[i].neigh_e1);
      connect3d[i].ewhich_e1 = ipc->get(connect3d[i].ne1, pool3d[i].ewhich_e1);
      connect3d[i].nside_e1 = ipc->get(connect3d[i].ne1, pool3d[i].nside_e1);
      connect3d[i].aflag_e1 = ipc->get(connect3d[i].ne1, pool3d[i].aflag_e1);
      connect3d[i].fflag_e1 = ipc->get(connect3d[i].ne1, pool3d[i].fflag_e1);
      for (int j = 0; j < connect3d[i].ne1; j++) {
        connect3d[i].ewhich_e1[j] = 0;
        connect3d[i].nside_e1[j] = 0;
        connect3d[i].aflag_e1[j] = 0;
        connect3d[i].fflag_e1[j] = 0;
      }
    } else {
      connect3d[i].neigh_e1 = nullptr;
      connect3d[i].ewhich_e1 = nullptr;
      connect3d[i].nside_e1 = nullptr;
      connect3d[i].aflag_e1 = nullptr;
      connect3d[i].fflag_e1 = nullptr;
    }
    if (connect3d[i].ne2) {
      connect3d[i].neigh_e2 = tpc->get(connect3d[i].ne2, pool3d[i].neigh_e2);
      connect3d[i].ewhich_e2 = ipc->get(connect3d[i].ne2, pool3d[i].ewhich_e2);
      connect3d[i].nside_e2 = ipc->get(connect3d[i].ne2, pool3d[i].nside_e2);
      connect3d[i].aflag_e2 = ipc->get(connect3d[i].ne2, pool3d[i].aflag_e2);
      connect3d[i].fflag_e2 = ipc->get(connect3d[i].ne2, pool3d[i].fflag_e2);
      for (int j = 0; j < connect3d[i].ne2; j++) {
        connect3d[i].ewhich_e2[j] = 0;
        connect3d[i].nside_e2[j] = 0;
        connect3d[i].aflag_e2[j] = 0;
        connect3d[i].fflag_e2[j] = 0;
      }
    } else {
      connect3d[i].neigh_e2 = nullptr;
      connect3d[i].ewhich_e2 = nullptr;
      connect3d[i].nside_e2 = nullptr;
      connect3d[i].aflag_e2 = nullptr;
      connect3d[i].fflag_e2 = nullptr;
    }
    if (connect3d[i].ne3) {
      connect3d[i].neigh_e3 = tpc->get(connect3d[i].ne3, pool3d[i].neigh_e3);
      connect3d[i].ewhich_e3 = ipc->get(connect3d[i].ne3, pool3d[i].ewhich_e3);
      connect3d[i].nside_e3 = ipc->get(connect3d[i].ne3, pool3d[i].nside_e3);
      connect3d[i].aflag_e3 = ipc->get(connect3d[i].ne3, pool3d[i].aflag_e3);
      connect3d[i].fflag_e3 = ipc->get(connect3d[i].ne3, pool3d[i].fflag_e3);
      for (int j = 0; j < connect3d[i].ne3; j++) {
        connect3d[i].ewhich_e3[j] = 0;
        connect3d[i].nside_e3[j] = 0;
        connect3d[i].aflag_e3[j] = 0;
        connect3d[i].fflag_e3[j] = 0;
      }
    } else {
      connect3d[i].neigh_e3 = nullptr;
      connect3d[i].ewhich_e3 = nullptr;
      connect3d[i].nside_e3 = nullptr;
      connect3d[i].aflag_e3 = nullptr;
      connect3d[i].fflag_e3 = nullptr;
    }

    if (connect3d[i].nc1) {
      connect3d[i].neigh_c1 = tpc->get(connect3d[i].nc1, pool3d[i].neigh_c1);
      connect3d[i].cwhich_c1 = ipc->get(connect3d[i].nc1, pool3d[i].cwhich_c1);
      connect3d[i].nside_c1 = ipc->get(connect3d[i].nc1, pool3d[i].nside_c1);
      connect3d[i].fflag_c1 = ipc->get(connect3d[i].nc1, pool3d[i].fflag_c1);
      for (int j = 0; j < connect3d[i].nc1; j++) {
        connect3d[i].cwhich_c1[j] = 0;
        connect3d[i].nside_c1[j] = 0;
        connect3d[i].fflag_c1[j] = 0;
      }
    } else {
      connect3d[i].neigh_c1 = nullptr;
      connect3d[i].cwhich_c1 = nullptr;
      connect3d[i].nside_c1 = nullptr;
      connect3d[i].fflag_c1 = nullptr;
    }
    if (connect3d[i].nc2) {
      connect3d[i].neigh_c2 = tpc->get(connect3d[i].nc2, pool3d[i].neigh_c2);
      connect3d[i].cwhich_c2 = ipc->get(connect3d[i].nc2, pool3d[i].cwhich_c2);
      connect3d[i].nside_c2 = ipc->get(connect3d[i].nc2, pool3d[i].nside_c2);
      connect3d[i].fflag_c2 = ipc->get(connect3d[i].nc2, pool3d[i].fflag_c2);
      for (int j = 0; j < connect3d[i].nc2; j++) {
        connect3d[i].cwhich_c2[j] = 0;
        connect3d[i].nside_c2[j] = 0;
        connect3d[i].fflag_c2[j] = 0;
      }
    } else {
      connect3d[i].neigh_c2 = nullptr;
      connect3d[i].cwhich_c2 = nullptr;
      connect3d[i].nside_c2 = nullptr;
      connect3d[i].fflag_c2 = nullptr;
    }
    if (connect3d[i].nc3) {
      connect3d[i].neigh_c3 = tpc->get(connect3d[i].nc3, pool3d[i].neigh_c3);
      connect3d[i].cwhich_c3 = ipc->get(connect3d[i].nc3, pool3d[i].cwhich_c3);
      connect3d[i].nside_c3 = ipc->get(connect3d[i].nc3, pool3d[i].nside_c3);
      connect3d[i].fflag_c3 = ipc->get(connect3d[i].nc3, pool3d[i].fflag_c3);
      for (int j = 0; j < connect3d[i].nc3; j++) {
        connect3d[i].cwhich_c3[j] = 0;
        connect3d[i].nside_c3[j] = 0;
        connect3d[i].fflag_c3[j] = 0;
      }
    } else {
      connect3d[i].neigh_c3 = nullptr;
      connect3d[i].cwhich_c3 = nullptr;
      connect3d[i].nside_c3 = nullptr;
      connect3d[i].fflag_c3 = nullptr;
    }
  }

  // populate neigh_e123 vectors within Connect3d

  for (int i = 0; i < nlocal_connect; i++)
    connect3d[i].ne1 = connect3d[i].ne2 = connect3d[i].ne3 = 0;

  for (int i = 0; i < nlocal_connect; i++) {
    n1 = n1_counts[i];
    n2 = n2_counts[i];
    neigh1 = tneigh1[i];
    neigh2 = tneigh2[i];
    for (int j = 0; j < n1; j++) {
      for (int k = 0; k < n2; k++) {
        if (neigh2[k] == neigh1[j]) {
          connect3d[i].neigh_e1[connect3d[i].ne1++] = neigh1[j];
          break;
        }
      }
    }

    n1 = n2_counts[i];
    n2 = n3_counts[i];
    neigh1 = tneigh2[i];
    neigh2 = tneigh3[i];
    for (int j = 0; j < n1; j++) {
      for (int k = 0; k < n2; k++) {
        if (neigh2[k] == neigh1[j]) {
          connect3d[i].neigh_e2[connect3d[i].ne2++] = neigh1[j];
          break;
        }
      }
    }

    n1 = n3_counts[i];
    n2 = n1_counts[i];
    neigh1 = tneigh3[i];
    neigh2 = tneigh1[i];
    for (int j = 0; j < n1; j++) {
      for (int k = 0; k < n2; k++) {
        if (neigh2[k] == neigh1[j]) {
          connect3d[i].neigh_e3[connect3d[i].ne3++] = neigh1[j];
          break;
        }
      }
    }
  }

  // populate neigh_c123 vectors within Connect3d

  for (int i = 0; i < nlocal_connect; i++)
    connect3d[i].nc1 = connect3d[i].nc2 = connect3d[i].nc3 = 0;

  int flag;

  for (int i = 0; i < nlocal_connect; i++) {
    n = n1_counts[i];
    neigh = tneigh1[i];
    for (int j = 0; j < n; j++) {
      flag = 0;
      for (int k = 0; k < connect3d[i].ne1; k++) {
        if (neigh[j] == connect3d[i].neigh_e1[k]) {
          flag = 1;
          break;
        }
      }
      for (int k = 0; k < connect3d[i].ne3; k++) {
        if (neigh[j] == connect3d[i].neigh_e3[k]) {
          flag = 1;
          break;
        }
      }
      if (!flag) connect3d[i].neigh_c1[connect3d[i].nc1++] = neigh[j];
    }

    n = n2_counts[i];
    neigh = tneigh2[i];
    for (int j = 0; j < n; j++) {
      flag = 0;
      for (int k = 0; k < connect3d[i].ne2; k++) {
        if (neigh[j] == connect3d[i].neigh_e2[k]) {
          flag = 1;
          break;
        }
      }
      for (int k = 0; k < connect3d[i].ne1; k++) {
        if (neigh[j] == connect3d[i].neigh_e1[k]) {
          flag = 1;
          break;
        }
      }
      if (!flag) connect3d[i].neigh_c2[connect3d[i].nc2++] = neigh[j];
    }

    n = n3_counts[i];
    neigh = tneigh3[i];
    for (int j = 0; j < n; j++) {
      flag = 0;
      for (int k = 0; k < connect3d[i].ne3; k++) {
        if (neigh[j] == connect3d[i].neigh_e3[k]) {
          flag = 1;
          break;
        }
      }
      for (int k = 0; k < connect3d[i].ne2; k++) {
        if (neigh[j] == connect3d[i].neigh_e2[k]) {
          flag = 1;
          break;
        }
      }
      if (!flag) connect3d[i].neigh_c3[connect3d[i].nc3++] = neigh[j];
    }
  }

  // clean up

  memory->destroy(n1_counts);
  memory->destroy(n2_counts);
  memory->destroy(n3_counts);
  memory->destroy(tneigh1);
  memory->destroy(tneigh2);
  memory->destroy(tneigh3);
}

/* ----------------------------------------------------------------------
   callback from comm->rvous() for 2d or 3d
   operates in Rvous decomposition of bins
------------------------------------------------------------------------- */

int FixSurfaceLocal::point_match(int n, char *inbuf, int &rflag, int *&proclist, char *&outbuf,
                                 void *ptr)
{
  // access class data for epssq and bin count

  auto fslptr = (FixSurfaceLocal *) ptr;
  Memory *memory = fslptr->memory;
  Comm *comm = fslptr->comm;

  double epssq = fslptr->epssq;
  int nbins = fslptr->nbins;
  int nprocs = comm->nprocs;
  int me = comm->me;

  // nmine = # of bins I own
  // num[nmine] = count of datums in each bin
  // first[nmine] = index of 1st datum in each bin
  // next[n] = index to next datum in each bin

  int nmine = nbins / nprocs;
  if (me < nbins % nprocs) nmine++;

  int *num, *first, *next;
  memory->create(num, nmine, "surface/local:num");
  memory->create(first, nmine, "surface/local:first");
  memory->create(next, n, "surface/local:next");

  for (int i = 0; i < nmine; i++) num[i] = 0;
  for (int i = 0; i < nmine; i++) first[i] = -1;

  auto in = (InRvous *) inbuf;

  int ibin, whichbin;

  for (int i = n - 1; i >= 0; i--) {
    ibin = in[i].ibin;
    whichbin = ibin / nprocs;
    if (first[whichbin] < 0)
      next[i] = -1;
    else
      next[i] = first[whichbin];
    first[whichbin] = i;
    num[whichbin]++;
  }

  // double loop over datums in each bin to to identify point matches
  // match = 2 points within EPS distance of each other
  // add match to outbuf twice:
  //   once for each of the 2 points
  //   send each the atomID of other line/tri it matches

  proclist = nullptr;
  OutRvous *out = nullptr;
  int ncount = 0;
  int maxcount = 0;

  double dx, dy, dz, rsq;

  for (int ibin = 0; ibin < nmine; ibin++) {
    int i = first[ibin];

    while (i >= 0) {
      int j = first[ibin];
      while (j >= 0) {
        if (j == i) break;
        dx = in[i].x[0] - in[j].x[0];
        dy = in[i].x[1] - in[j].x[1];
        dz = in[i].x[2] - in[j].x[2];
        rsq = dx * dx + dy * dy + dz * dz;
        if (rsq < epssq && in[i].mol == in[j].mol) {
          if (ncount + 2 > maxcount) {
            maxcount += DELTA_RVOUS;
            memory->grow(proclist, maxcount, "surface/local:proclist");
            out = (OutRvous *) memory->srealloc(out, maxcount * sizeof(OutRvous),
                                                "surface/local:outbuf");
          }

          proclist[ncount] = in[i].proc;
          out[ncount].ilocal = in[i].ilocal;
          out[ncount].ipoint = in[i].ipoint;
          out[ncount].atomID = in[j].atomID;
          ncount++;

          proclist[ncount] = in[j].proc;
          out[ncount].ilocal = in[j].ilocal;
          out[ncount].ipoint = in[j].ipoint;
          out[ncount].atomID = in[i].atomID;
          ncount++;
        }
        j = next[j];
      }
      i = next[i];
    }
  }

  // clean up

  memory->destroy(num);
  memory->destroy(first);
  memory->destroy(next);

  // return values

  rflag = 2;
  outbuf = (char *) out;
  return ncount;
}

/* ----------------------------------------------------------------------
   compute current end points of my owned lines
------------------------------------------------------------------------- */

void FixSurfaceLocal::calculate_endpts(double **endpts)
{
  double length, theta, dx, dy;

  AtomVecLine::Bonus *bonus = avec_line->bonus;
  double **x = atom->x;
  int *line = atom->line;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    length = bonus[line[i]].length;
    theta = bonus[line[i]].theta;
    dx = 0.5 * length * cos(theta);
    dy = 0.5 * length * sin(theta);
    endpts[m][0] = x[i][0] - dx;
    endpts[m][1] = x[i][1] - dy;
    endpts[m][2] = x[i][0] + dx;
    endpts[m][3] = x[i][1] + dy;
    m++;
  }
}

/* ----------------------------------------------------------------------
   calculate all bins which pt +/- EPS overlaps with
   return N = # of overlapped bins
   return indices = list of bin IDs (0 to Nbins-1)
------------------------------------------------------------------------- */

int FixSurfaceLocal::overlap_bins_2d(double *pt, double eps, int *indices)
{
  int ilo = static_cast<int>((pt[0] - eps - bboxlo[0]) * invbinx);
  int ihi = static_cast<int>((pt[0] + eps - bboxlo[0]) * invbinx);
  int jlo = static_cast<int>((pt[1] - eps - bboxlo[1]) * invbiny);
  int jhi = static_cast<int>((pt[1] + eps - bboxlo[1]) * invbiny);

  ilo = MAX(ilo, 0);
  ihi = MIN(ihi, nbinx - 1);
  jlo = MAX(jlo, 0);
  jhi = MIN(jhi, nbiny - 1);

  if (ilo == ihi && jlo == jhi) {
    indices[0] = jlo * nbinx + ilo;
    return 1;
  }

  int n = 0;
  for (int j = jlo; j <= jhi; j++)
    for (int i = ilo; i <= ihi; i++) indices[n++] = j * nbinx + i;
  return n;
}

/* ----------------------------------------------------------------------
   compute current corner points of my owned triangles
------------------------------------------------------------------------- */

void FixSurfaceLocal::calculate_corners(double **corners)
{
  int ibonus;
  double p[3][3];
  double *corner;

  AtomVecTri::Bonus *bonus = avec_tri->bonus;
  double **x = atom->x;
  int *tri = atom->tri;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    ibonus = tri[i];
    corner = corners[m];
    MathExtra::quat_to_mat(bonus[ibonus].quat, p);
    MathExtra::matvec(p, bonus[ibonus].c1, &corner[0]);
    MathExtra::add3(x[i], &corner[0], &corner[0]);
    MathExtra::matvec(p, bonus[ibonus].c2, &corner[3]);
    MathExtra::add3(x[i], &corner[3], &corner[3]);
    MathExtra::matvec(p, bonus[ibonus].c3, &corner[6]);
    MathExtra::add3(x[i], &corner[6], &corner[6]);
    m++;
  }
}

/* ----------------------------------------------------------------------
   map point +/ EPS in both dims to all bins it overlaps with
   return N = # of overlapped bins
   return indices = list of bin IDs, each from 0 to Nbins-1
------------------------------------------------------------------------- */

int FixSurfaceLocal::overlap_bins_3d(double *pt, double eps, int *indices)
{
  int ilo = static_cast<int>((pt[0] - eps - bboxlo[0]) * invbinx);
  int ihi = static_cast<int>((pt[0] + eps - bboxlo[0]) * invbinx);
  int jlo = static_cast<int>((pt[1] - eps - bboxlo[1]) * invbiny);
  int jhi = static_cast<int>((pt[1] + eps - bboxlo[1]) * invbiny);
  int klo = static_cast<int>((pt[2] - eps - bboxlo[2]) * invbinz);
  int khi = static_cast<int>((pt[2] + eps - bboxlo[2]) * invbinz);

  ilo = MAX(ilo, 0);
  ihi = MIN(ihi, nbinx - 1);
  jlo = MAX(jlo, 0);
  jhi = MIN(jhi, nbiny - 1);
  klo = MAX(klo, 0);
  khi = MIN(khi, nbinz - 1);

  if (ilo == ihi && jlo == jhi && klo == khi) {
    indices[0] = klo * nbiny * nbinx + jlo * nbinx + ilo;
    return 1;
  }

  int n = 0;
  for (int k = klo; k <= khi; k++)
    for (int j = jlo; j <= jhi; j++)
      for (int i = ilo; i <= ihi; i++) indices[n++] = k * nbiny * nbinx + j * nbinx + i;
  return n;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for global surf connectivity build, assignment to procs
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   check if distributed line/tri particles already exist
------------------------------------------------------------------------- */

int FixSurfaceLocal::check_exist()
{
  int *surf = atom->line;
  if (dimension == 3) surf = atom->tri;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (surf[i] >= 0) flag = 1;

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);

  if (flagall) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   assign lines and their connectivity from global data structs to each proc
   logic for proc assign follows Atom::data_atoms(), as if read from data file
------------------------------------------------------------------------- */

void FixSurfaceLocal::assign2d()
{
  AtomVec *avec = atom->avec;
  avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
  if (!avec_line)
    error->all(FLERR, Error::NOLASTLINE, "Fix surface/local requires atom style line");

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic)
    epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3], subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0];
    subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1];
    subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2];
    subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0];
    subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1];
    subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2];
    subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0] - 1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1] - 1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2] - 1) subhi[2] += epsilon[2];
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

  // idmaxall = largest existing atom ID
  // set atom2connect = -1 for all current particles

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint idmax = 0;
  for (int i = nlocal0; i < nlocal; i++) {
    idmax = MAX(tag[i], idmax);
    atom2connect[i] = -1;
  }
  tagint idmaxall;
  MPI_Allreduce(&idmax, &idmaxall, 1, MPI_LMP_TAGINT, MPI_MAX, world);

  // loop over all global lines
  // compute their midpoints
  // keep lines belonging to this proc

  int num;
  imageint imagedata;
  double xmid[3], lamda[3];
  double *coord, *x1, *x2;
  int *global, *local;
  std::vector<std::string> svalues(5);

  for (int i = 0; i < nlines; i++) {

    // midpt of line

    x1 = points[lines[i].p1].x;
    x2 = points[lines[i].p2].x;
    xmid[0] = 0.5 * (x1[0] + x2[0]);
    xmid[1] = 0.5 * (x1[1] + x2[1]);
    xmid[2] = 0.5 * (x1[2] + x2[2]);

    imagedata = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    domain->remap(xmid, imagedata);
    if (triclinic) {
      domain->x2lamda(xmid, lamda);
      coord = lamda;
    } else
      coord = xmid;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] && coord[1] >= sublo[1] &&
        coord[1] < subhi[1] && coord[2] >= sublo[2] && coord[2] < subhi[2]) {

      if (lines[i].type < 1 || lines[i].type > atom->ntypes)
        error->all(FLERR, Error::NOLASTLINE, "Invalid line type {} in fix surface/local",
                   lines[i].type);

      // create a default particle of correct style

      avec->create_atom(lines[i].type, xmid);

      // change it to be a line

      int n = atom->nlocal - 1;
      atom->tag[n] = idmaxall + i + 1;
      atom->image[n] = imagedata;
      atom->molecule[n] = lines[i].mol;
      atom->line[n] = 0;
      atom->rmass[n] = 1.0;    // set line density = 1.0
      atom->radius[n] = 0.0;

      // set the bonus data for a line = end points
      // svalues[0] is not used here
      // used when caller of data_atom_bonus() is via read_data command

      svalues[1] = fmt::format("{:.16e}", x1[0]);
      svalues[2] = fmt::format("{:.16e}", x1[1]);
      svalues[3] = fmt::format("{:.16e}", x2[0]);
      svalues[4] = fmt::format("{:.16e}", x2[1]);

      avec_line->data_atom_bonus(n, svalues);

      // allocate all vectors in Connect2d
      // copy neigh_p12 from global connect2dall to local connect2d
      //   reset global line indices to new line IDs beyond idmaxall
      // initialize other vectors to 0 for first-time borders comm
      //   they will be set correctly in connectivity2d_complete()

      if (nlocal_connect == nmax_connect) grow_connect();
      memcpy(&connect2d[nlocal_connect], &connect2dall[i], sizeof(Connect2d));

      num = connect2d[nlocal_connect].np1;

      if (num) {
        connect2d[nlocal_connect].neigh_p1 = tpc->get(num, pool2d[nlocal_connect].neigh_p1);
        connect2d[nlocal_connect].pwhich_p1 = ipc->get(num, pool2d[nlocal_connect].pwhich_p1);
        connect2d[nlocal_connect].nside_p1 = ipc->get(num, pool2d[nlocal_connect].nside_p1);
        connect2d[nlocal_connect].aflag_p1 = ipc->get(num, pool2d[nlocal_connect].aflag_p1);
        connect2d[nlocal_connect].fflag_p1 = ipc->get(num, pool2d[nlocal_connect].fflag_p1);
        global = connect2dall[i].neigh_p1;
        local = connect2d[nlocal_connect].neigh_p1;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect2d[nlocal_connect].pwhich_p1[j] = 0;
          connect2d[nlocal_connect].nside_p1[j] = 0;
          connect2d[nlocal_connect].aflag_p1[j] = 0;
          connect2d[nlocal_connect].fflag_p1[j] = 0;
        }
      } else {
        connect2d[nlocal_connect].neigh_p1 = nullptr;
        connect2d[nlocal_connect].pwhich_p1 = nullptr;
        connect2d[nlocal_connect].nside_p1 = nullptr;
        connect2d[nlocal_connect].aflag_p1 = nullptr;
        connect2d[nlocal_connect].fflag_p1 = nullptr;
      }

      num = connect2d[nlocal_connect].np2;

      if (num) {
        connect2d[nlocal_connect].neigh_p2 = tpc->get(num, pool2d[nlocal_connect].neigh_p2);
        connect2d[nlocal_connect].pwhich_p2 = ipc->get(num, pool2d[nlocal_connect].pwhich_p2);
        connect2d[nlocal_connect].nside_p2 = ipc->get(num, pool2d[nlocal_connect].nside_p2);
        connect2d[nlocal_connect].aflag_p2 = ipc->get(num, pool2d[nlocal_connect].aflag_p2);
        connect2d[nlocal_connect].fflag_p2 = ipc->get(num, pool2d[nlocal_connect].fflag_p2);
        global = connect2dall[i].neigh_p2;
        local = connect2d[nlocal_connect].neigh_p2;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect2d[nlocal_connect].pwhich_p2[j] = 0;
          connect2d[nlocal_connect].nside_p2[j] = 0;
          connect2d[nlocal_connect].aflag_p2[j] = 0;
          connect2d[nlocal_connect].fflag_p2[j] = 0;
        }
      } else {
        connect2d[nlocal_connect].neigh_p2 = nullptr;
        connect2d[nlocal_connect].pwhich_p2 = nullptr;
        connect2d[nlocal_connect].nside_p2 = nullptr;
        connect2d[nlocal_connect].aflag_p2 = nullptr;
        connect2d[nlocal_connect].fflag_p2 = nullptr;
      }

      connect2atom[nlocal_connect] = n;
      atom2connect[n] = nlocal_connect;
      nlocal_connect++;
    }
  }

  memory->sfree(connect2dall);
  memory->destroy(neigh_p1);
  memory->destroy(neigh_p2);

  // set new total # of atoms and error check

  atom->natoms += nlines;
  atom->nlines += nlines;
  if ((atom->natoms < 0) || (atom->natoms >= MAXBIGINT))
    error->all(FLERR, Error::NOLASTLINE, "Too many total atoms");

  // recreate atom map that includes added lines

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
}

/* ----------------------------------------------------------------------
   assign tris and their connectivity from global data structs to each proc
------------------------------------------------------------------------- */

void FixSurfaceLocal::assign3d()
{
  AtomVec *avec = atom->avec;
  avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
  if (!avec_tri) error->all(FLERR, Error::NOLASTLINE, "Fix surface/local requires atom style tri");

  // set bounds for my proc
  // if periodic and I am lo/hi proc, adjust bounds by EPSILON
  // insures all data atoms will be owned even with round-off

  int triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic)
    epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  double sublo[3], subhi[3];
  if (triclinic == 0) {
    sublo[0] = domain->sublo[0];
    subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1];
    subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2];
    subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0];
    subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1];
    subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2];
    subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0] - 1) subhi[0] += epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1] - 1) subhi[1] += epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2] - 1) subhi[2] += epsilon[2];
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

  // idmaxall = largest existing atom ID
  // set atom2connect = -1 for all current particles

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint idmax = 0;
  for (int i = nlocal0; i < nlocal; i++) {
    idmax = MAX(tag[i], idmax);
    atom2connect[i] = -1;
  }
  tagint idmaxall;
  MPI_Allreduce(&idmax, &idmaxall, 1, MPI_LMP_TAGINT, MPI_MAX, world);

  // loop over all global triangles
  // compute their center points
  // keep lines belonging to this proc

  int num;
  imageint imagedata;
  double xmid[3], lamda[3];
  double *coord, *x1, *x2, *x3;
  int *global, *local;
  std::vector<std::string> svalues(10);

  for (int i = 0; i < ntris; i++) {

    // center pt of triangle

    x1 = points[tris[i].p1].x;
    x2 = points[tris[i].p2].x;
    x3 = points[tris[i].p3].x;
    xmid[0] = (x1[0] + x2[0] + x3[0]) / 3.0;
    xmid[1] = (x1[1] + x2[1] + x3[1]) / 3.0;
    xmid[2] = (x1[2] + x2[2] + x3[2]) / 3.0;

    imagedata = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

    domain->remap(xmid, imagedata);
    if (triclinic) {
      domain->x2lamda(xmid, lamda);
      coord = lamda;
    } else
      coord = xmid;

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] && coord[1] >= sublo[1] &&
        coord[1] < subhi[1] && coord[2] >= sublo[2] && coord[2] < subhi[2]) {

      if ((tris[i].type < 1) || (tris[i].type > atom->ntypes))
        error->all(FLERR, Error::NOLASTLINE, "Invalid triangle type {} in fix surface/local",
                   tris[i].type);

      // create a default particle of correct style

      avec->create_atom(tris[i].type, xmid);

      // change it to be a triangle

      int n = atom->nlocal - 1;
      atom->tag[n] = idmaxall + i + 1;
      atom->image[n] = imagedata;
      atom->molecule[n] = tris[i].mol;
      atom->tri[n] = 0;
      atom->rmass[n] = 1.0;    // set tri density = 1.0
      atom->radius[n] = 0.0;

      // set the bonus data for a tri = corner pts
      // svalues[0] is not used here
      // used when caller of data_atom_bonus() is via read_data command

      svalues[1] = fmt::format("{:.16e}", x1[0]);
      svalues[2] = fmt::format("{:.16e}", x1[1]);
      svalues[3] = fmt::format("{:.16e}", x1[2]);
      svalues[4] = fmt::format("{:.16e}", x2[0]);
      svalues[5] = fmt::format("{:.16e}", x2[1]);
      svalues[6] = fmt::format("{:.16e}", x2[2]);
      svalues[7] = fmt::format("{:.16e}", x3[0]);
      svalues[8] = fmt::format("{:.16e}", x3[1]);
      svalues[9] = fmt::format("{:.16e}", x3[2]);

      avec_tri->data_atom_bonus(n, svalues);

      // allocate all vectors in Connect3d
      // copy neigh_e123 and neigh_c123 from global connect3dall to local connect3d
      //   reset global tri indices to new tri IDs beyond idmaxall
      // initialize other vectors to 0 for first-time borders comm
      //   they will be set correctly in connectivity3d_complete()

      if (nlocal_connect == nmax_connect) grow_connect();
      memcpy(&connect3d[nlocal_connect], &connect3dall[i], sizeof(Connect3d));

      num = connect3d[nlocal_connect].ne1;
      if (num) {
        connect3d[nlocal_connect].neigh_e1 = tpc->get(num, pool3d[nlocal_connect].neigh_e1);
        connect3d[nlocal_connect].ewhich_e1 = ipc->get(num, pool3d[nlocal_connect].ewhich_e1);
        connect3d[nlocal_connect].nside_e1 = ipc->get(num, pool3d[nlocal_connect].nside_e1);
        connect3d[nlocal_connect].aflag_e1 = ipc->get(num, pool3d[nlocal_connect].aflag_e1);
        connect3d[nlocal_connect].fflag_e1 = ipc->get(num, pool3d[nlocal_connect].fflag_e1);
        global = connect3dall[i].neigh_e1;
        local = connect3d[nlocal_connect].neigh_e1;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect3d[nlocal_connect].ewhich_e1[j] = 0;
          connect3d[nlocal_connect].nside_e1[j] = 0;
          connect3d[nlocal_connect].aflag_e1[j] = 0;
          connect3d[nlocal_connect].fflag_e1[j] = 0;
        }
      } else {
        connect3d[nlocal_connect].neigh_e1 = nullptr;
        connect3d[nlocal_connect].ewhich_e1 = nullptr;
        connect3d[nlocal_connect].nside_e1 = nullptr;
        connect3d[nlocal_connect].aflag_e1 = nullptr;
        connect3d[nlocal_connect].fflag_e1 = nullptr;
      }

      num = connect3d[nlocal_connect].ne2;
      if (num) {
        connect3d[nlocal_connect].neigh_e2 = tpc->get(num, pool3d[nlocal_connect].neigh_e2);
        connect3d[nlocal_connect].ewhich_e2 = ipc->get(num, pool3d[nlocal_connect].ewhich_e2);
        connect3d[nlocal_connect].nside_e2 = ipc->get(num, pool3d[nlocal_connect].nside_e2);
        connect3d[nlocal_connect].aflag_e2 = ipc->get(num, pool3d[nlocal_connect].aflag_e2);
        connect3d[nlocal_connect].fflag_e2 = ipc->get(num, pool3d[nlocal_connect].fflag_e2);
        global = connect3dall[i].neigh_e2;
        local = connect3d[nlocal_connect].neigh_e2;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect3d[nlocal_connect].ewhich_e2[j] = 0;
          connect3d[nlocal_connect].nside_e2[j] = 0;
          connect3d[nlocal_connect].nside_e2[j] = 0;
          connect3d[nlocal_connect].aflag_e2[j] = 0;
          connect3d[nlocal_connect].fflag_e2[j] = 0;
        }
      } else {
        connect3d[nlocal_connect].neigh_e2 = nullptr;
        connect3d[nlocal_connect].ewhich_e2 = nullptr;
        connect3d[nlocal_connect].nside_e2 = nullptr;
        connect3d[nlocal_connect].aflag_e2 = nullptr;
        connect3d[nlocal_connect].fflag_e2 = nullptr;
      }

      num = connect3d[nlocal_connect].ne3;
      if (num) {
        connect3d[nlocal_connect].neigh_e3 = tpc->get(num, pool3d[nlocal_connect].neigh_e3);
        connect3d[nlocal_connect].ewhich_e3 = ipc->get(num, pool3d[nlocal_connect].ewhich_e3);
        connect3d[nlocal_connect].nside_e3 = ipc->get(num, pool3d[nlocal_connect].nside_e3);
        connect3d[nlocal_connect].aflag_e3 = ipc->get(num, pool3d[nlocal_connect].aflag_e3);
        connect3d[nlocal_connect].fflag_e3 = ipc->get(num, pool3d[nlocal_connect].fflag_e3);
        global = connect3dall[i].neigh_e3;
        local = connect3d[nlocal_connect].neigh_e3;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect3d[nlocal_connect].ewhich_e3[j] = 0;
          connect3d[nlocal_connect].nside_e3[j] = 0;
          connect3d[nlocal_connect].aflag_e3[j] = 0;
          connect3d[nlocal_connect].fflag_e3[j] = 0;
        }
      } else {
        connect3d[nlocal_connect].neigh_e3 = nullptr;
        connect3d[nlocal_connect].ewhich_e3 = nullptr;
        connect3d[nlocal_connect].nside_e3 = nullptr;
        connect3d[nlocal_connect].aflag_e3 = nullptr;
        connect3d[nlocal_connect].fflag_e3 = nullptr;
      }

      num = connect3d[nlocal_connect].nc1;
      if (num) {
        connect3d[nlocal_connect].neigh_c1 = tpc->get(num, pool3d[nlocal_connect].neigh_c1);
        connect3d[nlocal_connect].cwhich_c1 = ipc->get(num, pool3d[nlocal_connect].cwhich_c1);
        connect3d[nlocal_connect].nside_c1 = ipc->get(num, pool3d[nlocal_connect].nside_c1);
        connect3d[nlocal_connect].fflag_c1 = ipc->get(num, pool3d[nlocal_connect].fflag_c1);
        global = connect3dall[i].neigh_c1;
        local = connect3d[nlocal_connect].neigh_c1;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect3d[nlocal_connect].cwhich_c1[j] = 0;
          connect3d[nlocal_connect].nside_c1[j] = 0;
          connect3d[nlocal_connect].fflag_c1[j] = 0;
        }
      } else {
        connect3d[nlocal_connect].neigh_c1 = nullptr;
        connect3d[nlocal_connect].cwhich_c1 = nullptr;
        connect3d[nlocal_connect].nside_c1 = nullptr;
        connect3d[nlocal_connect].fflag_c1 = nullptr;
      }

      num = connect3d[nlocal_connect].nc2;
      if (num) {
        connect3d[nlocal_connect].neigh_c2 = tpc->get(num, pool3d[nlocal_connect].neigh_c2);
        connect3d[nlocal_connect].cwhich_c2 = ipc->get(num, pool3d[nlocal_connect].cwhich_c2);
        connect3d[nlocal_connect].nside_c2 = ipc->get(num, pool3d[nlocal_connect].nside_c2);
        connect3d[nlocal_connect].fflag_c2 = ipc->get(num, pool3d[nlocal_connect].fflag_c2);
        global = connect3dall[i].neigh_c2;
        local = connect3d[nlocal_connect].neigh_c2;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect3d[nlocal_connect].cwhich_c2[j] = 0;
          connect3d[nlocal_connect].nside_c2[j] = 0;
          connect3d[nlocal_connect].fflag_c2[j] = 0;
        }
      } else {
        connect3d[nlocal_connect].neigh_c2 = nullptr;
        connect3d[nlocal_connect].cwhich_c2 = nullptr;
        connect3d[nlocal_connect].nside_c2 = nullptr;
        connect3d[nlocal_connect].fflag_c2 = nullptr;
      }

      num = connect3d[nlocal_connect].nc3;
      if (num) {
        connect3d[nlocal_connect].neigh_c3 = tpc->get(num, pool3d[nlocal_connect].neigh_c3);
        connect3d[nlocal_connect].cwhich_c3 = ipc->get(num, pool3d[nlocal_connect].cwhich_c3);
        connect3d[nlocal_connect].nside_c3 = ipc->get(num, pool3d[nlocal_connect].nside_c3);
        connect3d[nlocal_connect].fflag_c3 = ipc->get(num, pool3d[nlocal_connect].fflag_c3);
        global = connect3dall[i].neigh_c3;
        local = connect3d[nlocal_connect].neigh_c3;
        for (int j = 0; j < num; j++) {
          local[j] = global[j] + idmaxall + 1;
          connect3d[nlocal_connect].cwhich_c3[j] = 0;
          connect3d[nlocal_connect].nside_c3[j] = 0;
          connect3d[nlocal_connect].fflag_c3[j] = 0;
        }
      } else {
        connect3d[nlocal_connect].neigh_c3 = nullptr;
        connect3d[nlocal_connect].cwhich_c3 = nullptr;
        connect3d[nlocal_connect].nside_c3 = nullptr;
        connect3d[nlocal_connect].fflag_c3 = nullptr;
      }

      connect2atom[nlocal_connect] = n;
      atom2connect[n] = nlocal_connect;
      nlocal_connect++;
    }
  }

  memory->sfree(connect3dall);
  memory->destroy(neigh_e1);
  memory->destroy(neigh_e2);
  memory->destroy(neigh_e3);
  memory->destroy(neigh_c1);
  memory->destroy(neigh_c2);
  memory->destroy(neigh_c3);

  // set new total # of atoms and error check

  atom->natoms += ntris;
  atom->ntris += ntris;
  if ((atom->natoms < 0) || (atom->natoms >= MAXBIGINT))
    error->all(FLERR, Error::NOLASTLINE, "Too many total atoms");

  // recreate atom map that includes added tris

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for final stages of surf connectivity build, for global or local
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   assign values to all other fields in Connect2d for line connections
   only np1,np2 and neigh_p12 were previously assigned
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity2d_complete()
{
  // epssq = square of EPSILON fraction of minimum line length
  // for use in same_point()

  if (epssq < 0.0) epsilon_calculate();

  // calculate endpts and normals for all owned+ghost line particles

  double ***endpts, **normals;
  memory->create(endpts, nlocal_connect + nghost_connect, 2, 3, "surface/local:endpts");
  memory->create(normals, nlocal_connect + nghost_connect, 3, "surface/local:normals");

  AtomVecLine *avec = (AtomVecLine *) atom->style_match("line");
  AtomVecLine::Bonus *bonus = avec->bonus;
  double **x = atom->x;
  int *line = atom->line;
  int nall = atom->nlocal + atom->nghost;

  int m, iconnect;
  double length, theta, dx, dy;
  double p12[3];
  double zunit[3] = {0.0, 0.0, 1.0};

  for (int i = 0; i < nall; i++) {
    if (line[i] < 0) continue;
    m = line[i];
    length = bonus[m].length;
    theta = bonus[m].theta;
    dx = 0.5 * length * cos(theta);
    dy = 0.5 * length * sin(theta);

    iconnect = atom2connect[i];

    endpts[iconnect][0][0] = x[i][0] - dx;
    endpts[iconnect][0][1] = x[i][1] - dy;
    endpts[iconnect][0][2] = 0.0;
    endpts[iconnect][1][0] = x[i][0] + dx;
    endpts[iconnect][1][1] = x[i][1] + dy;
    endpts[iconnect][1][2] = 0.0;

    MathExtra::sub3(endpts[iconnect][1], endpts[iconnect][0], p12);
    MathExtra::cross3(zunit, p12, normals[iconnect]);
    MathExtra::norm3(normals[iconnect]);
  }

  // set connect2d pwhich/nside/aflag for each end point of each line
  // see fsl.h file for an explanation of each vector in Connect2d
  // aflag is based on dot and cross product of 2 connected line normals
  //   cross product is either along +z or -z direction

  int nlocal = atom->nlocal;

  int i, j, jconnect;
  tagint jtag;

  for (int i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    iconnect = atom2connect[i];

    for (int m = 0; m < connect2d[iconnect].np1; m++) {
      jtag = connect2d[iconnect].neigh_p1[m];
      j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];

      int jwhich = -1;
      if (same_point(endpts[iconnect][0], endpts[jconnect][0]))
        jwhich = 0;
      else if (same_point(endpts[iconnect][0], endpts[jconnect][1]))
        jwhich = 1;
      point_connection2d(normals[iconnect], normals[jconnect], 0, jwhich, flatthresh,
                         connect2d[iconnect].fflag_p1[m], connect2d[iconnect].nside_p1[m],
                         connect2d[iconnect].aflag_p1[m]);
      connect2d[iconnect].pwhich_p1[m] = jwhich;
    }

    for (int m = 0; m < connect2d[iconnect].np2; m++) {
      jtag = connect2d[iconnect].neigh_p2[m];
      j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];

      int jwhich = -1;
      if (same_point(endpts[iconnect][1], endpts[jconnect][0]))
        jwhich = 0;
      else if (same_point(endpts[iconnect][1], endpts[jconnect][1]))
        jwhich = 1;
      point_connection2d(normals[iconnect], normals[jconnect], 1, jwhich, flatthresh,
                         connect2d[iconnect].fflag_p2[m], connect2d[iconnect].nside_p2[m],
                         connect2d[iconnect].aflag_p2[m]);
      connect2d[iconnect].pwhich_p2[m] = jwhich;
    }
  }

  memory->destroy(endpts);
  memory->destroy(normals);

  // determine whether pts are external based on connectivity

  for (i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    iconnect = atom2connect[i];
    connect2d[iconnect].external_pt[0] = INTERNAL;
    connect2d[iconnect].external_pt[1] = INTERNAL;
  }

  int n;
  for (int i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    iconnect = atom2connect[i];
    // external if there's a nonflat connection
    for (n = 0; n < connect2d[iconnect].np1; n++)
      if (connect2d[iconnect].fflag_p1[n] != FLAT) connect2d[iconnect].external_pt[0] = EXTERNAL;
    for (n = 0; n < connect2d[iconnect].np2; n++)
      if (connect2d[iconnect].fflag_p2[n] != FLAT) connect2d[iconnect].external_pt[1] = EXTERNAL;
    // or if unconnected on border
    if (connect2d[iconnect].np1 == 0) connect2d[iconnect].external_pt[0] = UNCONNECTED;
    if (connect2d[iconnect].np2 == 0) connect2d[iconnect].external_pt[1] = UNCONNECTED;
  }
}

/* ----------------------------------------------------------------------
   assign values to all other fields in Connect3d for tri connections
   only ne1,ne2,ne3 and neigh_e123 were previously assigned
   only nc1,nc2,nc3 and neigh_c123 were previously assigned
------------------------------------------------------------------------- */

void FixSurfaceLocal::connectivity3d_complete()
{
  // epssq = square of EPSILON fraction of minimum tri diameter
  // for use in same_point()

  if (epssq < 0.0) epsilon_calculate();

  // calculate corner pts and normals for all owned+ghost tri particles

  double ***cpts, **normals;
  memory->create(cpts, nlocal_connect + nghost_connect, 3, 3, "surface/local:cpts");
  memory->create(normals, nlocal_connect + nghost_connect, 3, "surface/local:normals");

  AtomVecTri *avec = (AtomVecTri *) atom->style_match("tri");
  AtomVecTri::Bonus *bonus = avec->bonus;
  double **x = atom->x;
  int *tri = atom->tri;
  int nall = atom->nlocal + atom->nghost;

  int m, iconnect;
  double p[3][3], p12[3], p13[3];

  for (int i = 0; i < nall; i++) {
    if (tri[i] < 0) continue;
    m = tri[i];
    iconnect = atom2connect[i];

    MathExtra::quat_to_mat(bonus[m].quat, p);
    MathExtra::matvec(p, bonus[m].c1, cpts[iconnect][0]);
    MathExtra::add3(x[i], cpts[iconnect][0], cpts[iconnect][0]);
    MathExtra::matvec(p, bonus[m].c2, cpts[iconnect][1]);
    MathExtra::add3(x[i], cpts[iconnect][1], cpts[iconnect][1]);
    MathExtra::matvec(p, bonus[m].c3, cpts[iconnect][2]);
    MathExtra::add3(x[i], cpts[iconnect][2], cpts[iconnect][2]);

    MathExtra::sub3(cpts[iconnect][1], cpts[iconnect][0], p12);
    MathExtra::sub3(cpts[iconnect][2], cpts[iconnect][0], p13);
    MathExtra::cross3(p12, p13, normals[iconnect]);
    MathExtra::norm3(normals[iconnect]);
  }

  // set connect3d edge ewhich/nside/aflag for each edge of each tri
  // see fsl.h file for an explanation of each edge vector in Connect3d
  // aflag is based on dot and cross product of 2 connected tri normals
  //   cross product is either along itri edge or in opposite dir

  int nlocal = atom->nlocal;

  int jconnect, jpfirst, jpsecond;
  tagint jtag;
  double iedge[3];

  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    iconnect = atom2connect[i];

    for (int m = 0; m < connect3d[iconnect].ne1; m++) {
      jtag = connect3d[iconnect].neigh_e1[m];
      int j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];

      jpfirst = jpsecond = -1;
      if (same_point(cpts[iconnect][0], cpts[jconnect][0]))
        jpfirst = 1;
      else if (same_point(cpts[iconnect][0], cpts[jconnect][1]))
        jpfirst = 2;
      else if (same_point(cpts[iconnect][0], cpts[jconnect][2]))
        jpfirst = 3;

      if (same_point(cpts[iconnect][1], cpts[jconnect][0]))
        jpsecond = 1;
      else if (same_point(cpts[iconnect][1], cpts[jconnect][1]))
        jpsecond = 2;
      else if (same_point(cpts[iconnect][1], cpts[jconnect][2]))
        jpsecond = 3;

      MathExtra::sub3(cpts[iconnect][1], cpts[iconnect][0], iedge);
      edge_connection3d(normals[iconnect], normals[jconnect], iedge, jpfirst, jpsecond,
                        flatthresh, connect3d[iconnect].fflag_e1[m],
                        connect3d[iconnect].ewhich_e1[m], connect3d[iconnect].nside_e1[m],
                        connect3d[iconnect].aflag_e1[m]);
    }

    for (int m = 0; m < connect3d[iconnect].ne2; m++) {
      jtag = connect3d[iconnect].neigh_e2[m];
      int j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];

      jpfirst = jpsecond = -1;
      if (same_point(cpts[iconnect][1], cpts[jconnect][0]))
        jpfirst = 1;
      else if (same_point(cpts[iconnect][1], cpts[jconnect][1]))
        jpfirst = 2;
      else if (same_point(cpts[iconnect][1], cpts[jconnect][2]))
        jpfirst = 3;

      if (same_point(cpts[iconnect][2], cpts[jconnect][0]))
        jpsecond = 1;
      else if (same_point(cpts[iconnect][2], cpts[jconnect][1]))
        jpsecond = 2;
      else if (same_point(cpts[iconnect][2], cpts[jconnect][2]))
        jpsecond = 3;

      MathExtra::sub3(cpts[iconnect][2], cpts[iconnect][1], iedge);
      edge_connection3d(normals[iconnect], normals[jconnect], iedge, jpfirst, jpsecond,
                        flatthresh, connect3d[iconnect].fflag_e2[m],
                        connect3d[iconnect].ewhich_e2[m], connect3d[iconnect].nside_e2[m],
                        connect3d[iconnect].aflag_e2[m]);
    }

    for (int m = 0; m < connect3d[iconnect].ne3; m++) {
      jtag = connect3d[iconnect].neigh_e3[m];
      int j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];

      jpfirst = jpsecond = -1;
      if (same_point(cpts[iconnect][2], cpts[jconnect][0]))
        jpfirst = 1;
      else if (same_point(cpts[iconnect][2], cpts[jconnect][1]))
        jpfirst = 2;
      else if (same_point(cpts[iconnect][2], cpts[jconnect][2]))
        jpfirst = 3;

      if (same_point(cpts[iconnect][0], cpts[jconnect][0]))
        jpsecond = 1;
      else if (same_point(cpts[iconnect][0], cpts[jconnect][1]))
        jpsecond = 2;
      else if (same_point(cpts[iconnect][0], cpts[jconnect][2]))
        jpsecond = 3;

      MathExtra::sub3(cpts[iconnect][0], cpts[iconnect][2], iedge);
      edge_connection3d(normals[iconnect], normals[jconnect], iedge, jpfirst, jpsecond,
                        flatthresh, connect3d[iconnect].fflag_e3[m],
                        connect3d[iconnect].ewhich_e3[m], connect3d[iconnect].nside_e3[m],
                        connect3d[iconnect].aflag_e3[m]);
    }
  }

  // set connect3d cwhich for each corner point of each tri
  // see fsl.h file for an explanation of each corner vector in Connect3d
  // aflag is based on dot of 2 connected tri normals
  //   only designates flat vs. non-flat (treated like concave)

  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    iconnect = atom2connect[i];

    for (int m = 0; m < connect3d[iconnect].nc1; m++) {
      jtag = connect3d[iconnect].neigh_c1[m];
      int j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];
      int cwhich = -1;
      if (same_point(cpts[iconnect][0], cpts[jconnect][0]))
        cwhich = 0;
      else if (same_point(cpts[iconnect][0], cpts[jconnect][1]))
        cwhich = 1;
      else if (same_point(cpts[iconnect][0], cpts[jconnect][2]))
        cwhich = 2;
      corner_connection3d(normals[iconnect], normals[jconnect], cwhich, flatthresh,
                          connect3d[iconnect].fflag_c1[m], connect3d[iconnect].nside_c1[m]);
      connect3d[iconnect].cwhich_c1[m] = cwhich;
    }

    for (int m = 0; m < connect3d[iconnect].nc2; m++) {
      jtag = connect3d[iconnect].neigh_c2[m];
      int j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];
      int cwhich = -1;
      if (same_point(cpts[iconnect][1], cpts[jconnect][0]))
        cwhich = 0;
      else if (same_point(cpts[iconnect][1], cpts[jconnect][1]))
        cwhich = 1;
      else if (same_point(cpts[iconnect][1], cpts[jconnect][2]))
        cwhich = 2;
      corner_connection3d(normals[iconnect], normals[jconnect], cwhich, flatthresh,
                          connect3d[iconnect].fflag_c2[m], connect3d[iconnect].nside_c2[m]);
      connect3d[iconnect].cwhich_c2[m] = cwhich;
    }

    for (int m = 0; m < connect3d[iconnect].nc3; m++) {
      jtag = connect3d[iconnect].neigh_c3[m];
      int j = atom->map(jtag);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Missing tri atom {} from surface", jtag);
      jconnect = atom2connect[j];
      int cwhich = -1;
      if (same_point(cpts[iconnect][2], cpts[jconnect][0]))
        cwhich = 0;
      else if (same_point(cpts[iconnect][2], cpts[jconnect][1]))
        cwhich = 1;
      else if (same_point(cpts[iconnect][2], cpts[jconnect][2]))
        cwhich = 2;
      corner_connection3d(normals[iconnect], normals[jconnect], cwhich, flatthresh,
                          connect3d[iconnect].fflag_c3[m], connect3d[iconnect].nside_c3[m]);
      connect3d[iconnect].cwhich_c3[m] = cwhich;
    }
  }

  // clean up

  memory->destroy(cpts);
  memory->destroy(normals);

  // determine whether edges/pts are external based on connectivity

  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    iconnect = atom2connect[i];
    for (int a = 0; a < 3; a++) connect3d[iconnect].external_edge[a] = INTERNAL;
  }

  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    iconnect = atom2connect[i];

    // external edge if there's a nonflat connection
    for (int n = 0; n < connect3d[iconnect].ne1; n++)
      if (connect3d[iconnect].fflag_e1[n] != FLAT) connect3d[iconnect].external_edge[0] = EXTERNAL;
    for (int n = 0; n < connect3d[iconnect].ne2; n++)
      if (connect3d[iconnect].fflag_e2[n] != FLAT) connect3d[iconnect].external_edge[1] = EXTERNAL;
    for (int n = 0; n < connect3d[iconnect].ne3; n++)
      if (connect3d[iconnect].fflag_e3[n] != FLAT) connect3d[iconnect].external_edge[2] = EXTERNAL;

    // or unconnected on border
    if (connect3d[iconnect].ne1 == 0) connect3d[iconnect].external_edge[0] = UNCONNECTED;
    if (connect3d[iconnect].ne2 == 0) connect3d[iconnect].external_edge[1] = UNCONNECTED;
    if (connect3d[iconnect].ne3 == 0) connect3d[iconnect].external_edge[2] = UNCONNECTED;

    // corners basically inherit from connected edges
    connect3d[iconnect].external_cor[0] =
        MAX(connect3d[iconnect].external_edge[0], connect3d[iconnect].external_edge[2]);
    connect3d[iconnect].external_cor[1] =
        MAX(connect3d[iconnect].external_edge[0], connect3d[iconnect].external_edge[1]);
    connect3d[iconnect].external_cor[2] =
        MAX(connect3d[iconnect].external_edge[1], connect3d[iconnect].external_edge[2]);
  }
}

/* ----------------------------------------------------------------------
   return 1 if two points are same within eps
   else return 0
------------------------------------------------------------------------- */

int FixSurfaceLocal::same_point(double *pt1, double *pt2)
{
  double dx = pt1[0] - pt2[0];
  double dy = pt1[1] - pt2[1];
  double dz = pt1[2] - pt2[2];
  double rsq = dx * dx + dy * dy + dz * dz;
  if (rsq < epssq) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   calculate eps = threshold for matching 2 points from different lines/tris
   stored as eps and epssq
------------------------------------------------------------------------- */

void FixSurfaceLocal::epsilon_calculate()
{
  double minlen = BIG;
  if (dimension == 2) {
    AtomVecLine::Bonus *bonus = avec_line->bonus;
    int *line = atom->line;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (line[i] >= 0) minlen = MIN(minlen, bonus[line[i]].length);

    // eps = EPSILON fraction of minimum line length

  } else if (dimension == 3) {
    double *radius = atom->radius;
    int *tri = atom->tri;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (tri[i] >= 0) minlen = MIN(minlen, radius[i]);
    minlen *= 2.0;

    // eps = EPSILON fraction of minimum tri diameter
  }

  MPI_Allreduce(&minlen, &eps, 1, MPI_DOUBLE, MPI_MIN, world);
  eps *= EPSILON;
  epssq = eps * eps;
}

/* ----------------------------------------------------------------------
   error checks on lines
   no zero-length lines
   don't check for duplcate lines since endpt coords are inexact
------------------------------------------------------------------------- */

void FixSurfaceLocal::check2d()
{
  // check for zero length lines

  AtomVecLine::Bonus *bonus = avec_line->bonus;
  int *line = atom->line;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    if (bonus[line[i]].length == 0.0) flag++;
  }

  int allflag;
  MPI_Allreduce(&flag, &allflag, 1, MPI_INT, MPI_SUM, world);
  if (allflag)
    error->all(FLERR, Error::NOLASTLINE, "Fix surface/local defines {} zero-length lines", allflag);
}

/* ----------------------------------------------------------------------
   error checks on tris
   no tris with a zero-length edge
   don't check for duplcate tris since corner point coords are inexact
------------------------------------------------------------------------- */

void FixSurfaceLocal::check3d()
{
  // check for zero length tri edges via corner points in bonus

  AtomVecTri::Bonus *bonus = avec_tri->bonus;
  int *tri = atom->tri;
  int nlocal = atom->nlocal;

  double *c1, *c2, *c3;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    c1 = bonus[tri[i]].c1;
    c2 = bonus[tri[i]].c2;
    c3 = bonus[tri[i]].c3;

    if (c1[0] == c2[0] && c1[1] == c2[1] && c1[2] == c2[2]) flag++;
    if (c2[0] == c3[0] && c2[1] == c3[1] && c2[2] == c3[2]) flag++;
    if (c3[0] == c1[0] && c3[1] == c1[1] && c3[2] == c1[2]) flag++;
  }

  int allflag;
  MPI_Allreduce(&flag, &allflag, 1, MPI_INT, MPI_SUM, world);
  if (allflag)
    error->all(FLERR, Error::NOLASTLINE, "Fix surface/local defines {} zero-length triangle edges",
               allflag);
}

/* ----------------------------------------------------------------------
   print stats on all lines and their connections
------------------------------------------------------------------------- */

void FixSurfaceLocal::stats2d()
{
  AtomVecLine::Bonus *bonus = avec_line->bonus;
  int *line = atom->line;
  int nlocal = atom->nlocal;

  int nlines = 0;
  int nconnect = 0;
  int nfree = 0;
  double partial_point = 0.0;
  double minsize = BIG;
  double maxsize = 0.0;

  int j;
  double size;

  for (int i = 0; i < nlocal; i++) {
    if (line[i] < 0) continue;
    j = atom2connect[i];
    nlines++;
    nconnect += connect2d[j].np1 + connect2d[j].np2;

    // free point requires no connections
    // same partial point tally is done by each line which shares a point

    if (connect2d[j].np1 == 0)
      nfree++;
    else
      partial_point += 1.0 / (connect2d[j].np1 + 1.0);
    if (connect2d[j].np2 == 0)
      nfree++;
    else
      partial_point += 1.0 / (connect2d[j].np2 + 1.0);

    size = bonus[line[i]].length;
    minsize = MIN(minsize, size);
    maxsize = MAX(maxsize, size);
  }

  int alllines, allpoints, allconnect, allfree;
  double allpartial_point, allminsize, allmaxsize;

  MPI_Allreduce(&nlines, &alllines, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nconnect, &allconnect, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nfree, &allfree, 1, MPI_INT, MPI_SUM, world);

  MPI_Allreduce(&partial_point, &allpartial_point, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&minsize, &allminsize, 1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(&maxsize, &allmaxsize, 1, MPI_DOUBLE, MPI_MAX, world);

  allconnect /= 2;
  allpoints = (int) (allpartial_point + EPSILON) + allfree;

  if (comm->me == 0) {
    utils::logmesg(lmp, "Fix surface/local line segment creation:\n");
    utils::logmesg(lmp, "  {} lines\n", alllines);
    utils::logmesg(lmp, "  {} line end points\n", allpoints);
    utils::logmesg(lmp, "  {} end point connections\n", allconnect);
    utils::logmesg(lmp, "  {} free end points\n", allfree);
    utils::logmesg(lmp, "  {} min line length\n", allminsize);
    utils::logmesg(lmp, "  {} max line length\n", allmaxsize);
  }

  max_radius = maxsize * 0.5;
  MPI_Bcast(&max_radius, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   print stats on all tris and their connections
------------------------------------------------------------------------- */

void FixSurfaceLocal::stats3d()
{
  double size, area;
  double p[3][3];
  double c1[3], c2[3], c3[3];
  double delta[3], edge12[3], edge13[3], cross[3];

  double **x = atom->x;
  AtomVecTri::Bonus *bonus = avec_tri->bonus;
  int *tri = atom->tri;
  int nlocal = atom->nlocal;

  int ntris = 0;
  int nconnect_edge = 0;
  int nconnect_corner = 0;
  int nfree_edge = 0;
  int nfree_corner = 0;
  double partial_edge = 0.0;
  double partial_corner = 0.0;
  double minedge = BIG;
  double maxedge = 0.0;
  double minarea = BIG;
  double maxarea = 0.0;

  int j, ibonus;

  for (int i = 0; i < nlocal; i++) {
    if (tri[i] < 0) continue;
    j = atom2connect[i];
    ntris++;

    nconnect_edge += connect3d[j].ne1 + connect3d[j].ne2 + connect3d[j].ne3;
    nconnect_corner += connect3d[j].nc1 + connect3d[j].nc2 + connect3d[j].nc3;

    // free edge requires no connections
    // same partial edge tally is done by each tri which shares an edge

    if (connect3d[j].ne1 == 0)
      nfree_edge++;
    else
      partial_edge += 1.0 / (connect3d[j].ne1 + 1.0);
    if (connect3d[j].ne2 == 0)
      nfree_edge++;
    else
      partial_edge += 1.0 / (connect3d[j].ne2 + 1.0);
    if (connect3d[j].ne3 == 0)
      nfree_edge++;
    else
      partial_edge += 1.0 / (connect3d[j].ne3 + 1.0);

    // free corner requires 2 adjacent edges also have no connections
    // same partial corner tally is done by each tri which shares a corner point

    if (connect3d[j].nc1 == 0 && (connect3d[j].ne3 == 0 && connect3d[j].ne1 == 0))
      nfree_corner++;
    else
      partial_corner += 1.0 / (connect3d[j].ne3 + connect3d[j].ne1 + connect3d[j].nc1 + 1.0);

    if (connect3d[j].nc2 == 0 && (connect3d[j].ne1 == 0 && connect3d[j].ne2 == 0))
      nfree_corner++;
    else
      partial_corner += 1.0 / (connect3d[j].ne1 + connect3d[j].ne2 + connect3d[j].nc2 + 1.0);

    if (connect3d[j].nc3 == 0 && (connect3d[j].ne2 == 0 && connect3d[j].ne3 == 0))
      nfree_corner++;
    else
      partial_corner += 1.0 / (connect3d[j].ne2 + connect3d[j].ne3 + connect3d[j].nc3 + 1.0);

    ibonus = tri[i];

    MathExtra::quat_to_mat(bonus[ibonus].quat, p);
    MathExtra::matvec(p, bonus[ibonus].c1, c1);
    MathExtra::add3(x[i], c1, c1);
    MathExtra::matvec(p, bonus[ibonus].c2, c2);
    MathExtra::add3(x[i], c2, c2);
    MathExtra::matvec(p, bonus[ibonus].c3, c3);
    MathExtra::add3(x[i], c3, c3);

    MathExtra::sub3(c1, c2, delta);
    size = MathExtra::len3(delta);
    minedge = MIN(minedge, size);
    maxedge = MAX(maxedge, size);
    MathExtra::sub3(c2, c3, delta);
    size = MathExtra::len3(delta);
    minedge = MIN(minedge, size);
    maxedge = MAX(maxedge, size);
    MathExtra::sub3(c3, c1, delta);
    size = MathExtra::len3(delta);
    minedge = MIN(minedge, size);
    maxedge = MAX(maxedge, size);

    MathExtra::sub3(c2, c1, edge12);
    MathExtra::sub3(c3, c1, edge13);
    MathExtra::cross3(edge12, edge13, cross);
    area = 0.5 * MathExtra::len3(cross);
    minarea = MIN(minarea, area);
    maxarea = MAX(maxarea, area);
  }

  int alltris, alledges, allpoints;
  int allconnect_edge, allconnect_corner, allfree_edge, allfree_corner;
  double allpartial_edge, allpartial_corner;
  double allminedge, allmaxedge, allminarea, allmaxarea;

  MPI_Allreduce(&ntris, &alltris, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nconnect_edge, &allconnect_edge, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nconnect_corner, &allconnect_corner, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nfree_edge, &allfree_edge, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&nfree_corner, &allfree_corner, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&partial_edge, &allpartial_edge, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&partial_corner, &allpartial_corner, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&minedge, &allminedge, 1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(&maxedge, &allmaxedge, 1, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&minarea, &allminarea, 1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(&maxarea, &allmaxarea, 1, MPI_DOUBLE, MPI_MAX, world);

  allconnect_edge /= 2;
  allconnect_corner /= 2;
  alledges = (int) (allpartial_edge + EPSILON) + allfree_edge;
  allpoints = (int) (allpartial_corner + EPSILON) + allfree_corner;

  if (comm->me == 0) {
    utils::logmesg(lmp, "Fix surface/local triangle creation:\n");
    utils::logmesg(lmp, "  {} tris\n", alltris);
    utils::logmesg(lmp, "  {} tri edges\n", alledges);
    utils::logmesg(lmp, "  {} tri corner points\n", allpoints);
    utils::logmesg(lmp, "  {} edge connections\n", allconnect_edge);
    utils::logmesg(lmp, "  {} corner point connections\n", allconnect_corner);
    utils::logmesg(lmp, "  {} free edges\n", allfree_edge);
    utils::logmesg(lmp, "  {} free corner points\n", allfree_corner);
    utils::logmesg(lmp, "  {} min edge length\n", allminedge);
    utils::logmesg(lmp, "  {} max edge length\n", allmaxedge);
    utils::logmesg(lmp, "  {} min tri area\n", allminarea);
    utils::logmesg(lmp, "  {} max tri area\n", allmaxarea);
  }

  max_radius = allmaxedge * 0.5;
  MPI_Bcast(&max_radius, 1, MPI_DOUBLE, 0, world);
}

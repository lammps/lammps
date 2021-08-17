// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// NOTE: allow for bin center to be variables for sphere/cylinder

#include "compute_chunk_atom.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_store.h"
#include "group.h"
#include "input.h"
#include "lattice.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <map>
#include <utility>

using namespace LAMMPS_NS;
using namespace MathConst;

enum{LOWER,CENTER,UPPER,COORD};
enum{BOX,LATTICE,REDUCED};
enum{NODISCARD,MIXED,YESDISCARD};
enum{ONCE,NFREQ,EVERY};              // used in several files
enum{LIMITMAX,LIMITEXACT};

#define IDMAX 1024*1024

/* ---------------------------------------------------------------------- */

ComputeChunkAtom::ComputeChunkAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  chunk_volume_vec(nullptr), coord(nullptr), ichunk(nullptr), chunkID(nullptr),
  cfvid(nullptr), idregion(nullptr), region(nullptr), cchunk(nullptr), fchunk(nullptr),
  varatom(nullptr), id_fix(nullptr), fixstore(nullptr), lockfix(nullptr), chunk(nullptr),
  exclude(nullptr), hash(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute chunk/atom command");

  peratom_flag = 1;
  scalar_flag = 1;
  extscalar = 0;
  size_peratom_cols = 0;
  create_attribute = 1;

  // chunk style and its args

  int iarg = 0;

  binflag = 0;
  ncoord = 0;
  cfvid = nullptr;

  if (strcmp(arg[3],"bin/1d") == 0) {
    binflag = 1;
    which = ArgInfo::BIN1D;
    ncoord = 1;
    iarg = 4;
    readdim(narg,arg,iarg,0);
    iarg += 3;
  } else if (strcmp(arg[3],"bin/2d") == 0) {
    binflag = 1;
    which = ArgInfo::BIN2D;
    ncoord = 2;
    iarg = 4;
    readdim(narg,arg,iarg,0);
    readdim(narg,arg,iarg+3,1);
    iarg += 6;
  } else if (strcmp(arg[3],"bin/3d") == 0) {
    binflag = 1;
    which = ArgInfo::BIN3D;
    ncoord = 3;
    iarg = 4;
    readdim(narg,arg,iarg,0);
    readdim(narg,arg,iarg+3,1);
    readdim(narg,arg,iarg+6,2);
    iarg += 9;

  } else if (strcmp(arg[3],"bin/sphere") == 0) {
    binflag = 1;
    which = ArgInfo::BINSPHERE;
    ncoord = 1;
    iarg = 4;
    if (iarg+6 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
    sorigin_user[0] = utils::numeric(FLERR,arg[iarg],false,lmp);
    sorigin_user[1] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    sorigin_user[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    sradmin_user = utils::numeric(FLERR,arg[iarg+3],false,lmp);
    sradmax_user = utils::numeric(FLERR,arg[iarg+4],false,lmp);
    nsbin = utils::inumeric(FLERR,arg[iarg+5],false,lmp);
    iarg += 6;
  } else if (strcmp(arg[3],"bin/cylinder") == 0) {
    binflag = 1;
    which = ArgInfo::BINCYLINDER;
    ncoord = 2;
    iarg = 4;
    readdim(narg,arg,iarg,0);
    iarg += 3;
    if (dim[0] == 0) {
      cdim1 = 1;
      cdim2 = 2;
    } else if (dim[0] == 1) {
      cdim1 = 0;
      cdim2 = 2;
    } else {
      cdim1 = 0;
      cdim2 = 1;
    }
    if (iarg+5 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
    corigin_user[dim[0]] = 0.0;
    corigin_user[cdim1] = utils::numeric(FLERR,arg[iarg],false,lmp);
    corigin_user[cdim2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    cradmin_user = utils::numeric(FLERR,arg[iarg+2],false,lmp);
    cradmax_user = utils::numeric(FLERR,arg[iarg+3],false,lmp);

    ncbin = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
    iarg += 5;

  } else if (strcmp(arg[3],"type") == 0) {
    which = ArgInfo::TYPE;
    iarg = 4;
  } else if (strcmp(arg[3],"molecule") == 0) {
    which = ArgInfo::MOLECULE;
    iarg = 4;

  } else {
    ArgInfo argi(arg[3]);

    which = argi.get_type();
    argindex = argi.get_index1();
    cfvid = argi.copy_name();

    if ((which == ArgInfo::UNKNOWN) || (which == ArgInfo::NONE)
        || (argi.get_dim() > 1))
      error->all(FLERR,"Illegal compute chunk/atom command");
    iarg = 4;
  }

  // optional args

  regionflag = 0;
  idregion = nullptr;
  nchunksetflag = 0;
  nchunkflag = EVERY;
  limit = 0;
  limitstyle = LIMITMAX;
  limitfirst = 0;
  idsflag = EVERY;
  compress = 0;
  int discardsetflag = 0;
  discard = MIXED;
  minflag[0] = LOWER;
  minflag[1] = LOWER;
  minflag[2] = LOWER;
  maxflag[0] = UPPER;
  maxflag[1] = UPPER;
  maxflag[2] = UPPER;
  scaleflag = LATTICE;
  pbcflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for compute chunk/atom does not exist");
      idregion = utils::strdup(arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"nchunk") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"once") == 0) nchunkflag = ONCE;
      else if (strcmp(arg[iarg+1],"every") == 0) nchunkflag = EVERY;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      nchunksetflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"limit") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      limit = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (limit < 0) error->all(FLERR,"Illegal compute chunk/atom command");
      if (limit && !compress) limitfirst = 1;
      iarg += 2;
      if (limit) {
        if (iarg > narg)
          error->all(FLERR,"Illegal compute chunk/atom command");
        if (strcmp(arg[iarg],"max") == 0) limitstyle = LIMITMAX;
        else if (strcmp(arg[iarg],"exact") == 0) limitstyle = LIMITEXACT;
        else error->all(FLERR,"Illegal compute chunk/atom command");
        iarg++;
      }
    } else if (strcmp(arg[iarg],"ids") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"once") == 0) idsflag = ONCE;
      else if (strcmp(arg[iarg+1],"nfreq") == 0) idsflag = NFREQ;
      else if (strcmp(arg[iarg+1],"every") == 0) idsflag = EVERY;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      else if (strcmp(arg[iarg+1],"no") == 0) compress = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) compress = 1;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"discard") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"mixed") == 0) discard = MIXED;
      else if (strcmp(arg[iarg+1],"no") == 0) discard = NODISCARD;
      else if (strcmp(arg[iarg+1],"yes") == 0) discard = YESDISCARD;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      discardsetflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"bound") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      int idim = 0;
      if (strcmp(arg[iarg+1],"x") == 0) idim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) idim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) idim = 2;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      minflag[idim] = COORD;
      if (strcmp(arg[iarg+2],"lower") == 0) minflag[idim] = LOWER;
      else minvalue[idim] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      maxflag[idim] = COORD;
      if (strcmp(arg[iarg+3],"upper") == 0) maxflag[idim] = UPPER;
      else maxvalue[idim] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg+1],"reduced") == 0) scaleflag = REDUCED;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"no") == 0) pbcflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) pbcflag = 1;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute chunk/atom command");
  }

  // set nchunkflag and discard to default values if not explicitly set
  // for binning style, also check in init() if simulation box is static,
  //   which sets nchunkflag = ONCE

  if (!nchunksetflag) {
    if (binflag) {
      if (scaleflag == REDUCED) nchunkflag = ONCE;
      else nchunkflag = EVERY;
    }
    if (which == ArgInfo::TYPE) nchunkflag = ONCE;
    if (which == ArgInfo::MOLECULE) {
      if (regionflag) nchunkflag = EVERY;
      else nchunkflag = ONCE;
    }
    if (compress) nchunkflag = EVERY;
  }

  if (!discardsetflag) {
    if (binflag) discard = MIXED;
    else discard = YESDISCARD;
  }

  // error checks

  if (which == ArgInfo::MOLECULE && !atom->molecule_flag)
    error->all(FLERR,"Compute chunk/atom molecule for non-molecular system");

  if (!binflag && discard == MIXED)
    error->all(FLERR,"Compute chunk/atom without bins "
               "cannot use discard mixed");
  if (which == ArgInfo::BIN1D && delta[0] <= 0.0)
    error->all(FLERR,"Illegal compute chunk/atom command");
  if (which == ArgInfo::BIN2D && (delta[0] <= 0.0 || delta[1] <= 0.0))
    error->all(FLERR,"Illegal compute chunk/atom command");
  if (which == ArgInfo::BIN2D && (dim[0] == dim[1]))
    error->all(FLERR,"Illegal compute chunk/atom command");
  if (which == ArgInfo::BIN3D &&
      (delta[0] <= 0.0 || delta[1] <= 0.0 || delta[2] <= 0.0))
      error->all(FLERR,"Illegal compute chunk/atom command");
  if (which == ArgInfo::BIN3D &&
      (dim[0] == dim[1] || dim[1] == dim[2] || dim[0] == dim[2]))
      error->all(FLERR,"Illegal compute chunk/atom command");
  if (which == ArgInfo::BINSPHERE) {
    if (domain->dimension == 2 && sorigin_user[2] != 0.0)
      error->all(FLERR,"Compute chunk/atom sphere z origin must be 0.0 for 2d");
    if (sradmin_user < 0.0 || sradmin_user >= sradmax_user || nsbin < 1)
      error->all(FLERR,"Illegal compute chunk/atom command");
  }
  if (which == ArgInfo::BINCYLINDER) {
    if (delta[0] <= 0.0)
      error->all(FLERR,"Illegal compute chunk/atom command");
    if (domain->dimension == 2 && dim[0] != 2)
      error->all(FLERR,"Compute chunk/atom cylinder axis must be z for 2d");
    if (cradmin_user < 0.0 || cradmin_user >= cradmax_user || ncbin < 1)
      error->all(FLERR,"Illegal compute chunk/atom command");
  }

  if (which == ArgInfo::COMPUTE) {
    int icompute = modify->find_compute(cfvid);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for compute chunk /atom does not exist");
    if (modify->compute[icompute]->peratom_flag == 0)
      error->all(FLERR,
                 "Compute chunk/atom compute does not calculate "
                 "per-atom values");
    if (argindex == 0 &&
        modify->compute[icompute]->size_peratom_cols != 0)
      error->all(FLERR,"Compute chunk/atom compute does not "
                 "calculate a per-atom vector");
    if (argindex && modify->compute[icompute]->size_peratom_cols == 0)
      error->all(FLERR,"Compute chunk/atom compute does not "
                 "calculate a per-atom array");
    if (argindex &&
        argindex > modify->compute[icompute]->size_peratom_cols)
      error->all(FLERR,"Compute chunk/atom compute array is "
                 "accessed out-of-range");
  }

  if (which == ArgInfo::FIX) {
    int ifix = modify->find_fix(cfvid);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for compute chunk/atom does not exist");
    if (modify->fix[ifix]->peratom_flag == 0)
      error->all(FLERR,"Compute chunk/atom fix does not calculate "
                 "per-atom values");
    if (argindex == 0 && modify->fix[ifix]->size_peratom_cols != 0)
      error->all(FLERR,
                 "Compute chunk/atom fix does not calculate a per-atom vector");
    if (argindex && modify->fix[ifix]->size_peratom_cols == 0)
      error->all(FLERR,
                 "Compute chunk/atom fix does not calculate a per-atom array");
    if (argindex && argindex > modify->fix[ifix]->size_peratom_cols)
      error->all(FLERR,"Compute chunk/atom fix array is accessed out-of-range");
  }

  if (which == ArgInfo::VARIABLE) {
    int ivariable = input->variable->find(cfvid);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for compute chunk/atom does not exist");
    if (input->variable->atomstyle(ivariable) == 0)
      error->all(FLERR,"Compute chunk/atom variable is not "
                 "atom-style variable");
  }

  // setup scaling

  if (binflag) {
    if (domain->triclinic == 1 && scaleflag != REDUCED)
      error->all(FLERR,"Compute chunk/atom for triclinic boxes "
                 "requires units reduced");
  }

  if (scaleflag == LATTICE) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  } else xscale = yscale = zscale = 1.0;

  // apply scaling factors and cylinder dims orthogonal to axis

  if (binflag) {
    double scale = 1.0;
    if (which == ArgInfo::BIN1D || which == ArgInfo::BIN2D
        || which == ArgInfo::BIN3D || which == ArgInfo::BINCYLINDER) {
      if (which == ArgInfo::BIN1D || which == ArgInfo::BINCYLINDER) ndim = 1;
      if (which == ArgInfo::BIN2D) ndim = 2;
      if (which == ArgInfo::BIN3D) ndim = 3;
      for (int idim = 0; idim < ndim; idim++) {
        if (dim[idim] == 0) scale = xscale;
        else if (dim[idim] == 1) scale = yscale;
        else if (dim[idim] == 2) scale = zscale;
        delta[idim] *= scale;
        invdelta[idim] = 1.0/delta[idim];
        if (originflag[idim] == COORD) origin[idim] *= scale;
        if (minflag[idim] == COORD) minvalue[idim] *= scale;
        if (maxflag[idim] == COORD) maxvalue[idim] *= scale;
      }
    } else if (which == ArgInfo::BINSPHERE) {
      sorigin_user[0] *= xscale;
      sorigin_user[1] *= yscale;
      sorigin_user[2] *= zscale;
      sradmin_user *= xscale;     // radii are scaled by xscale
      sradmax_user *= xscale;
    } else if (which == ArgInfo::BINCYLINDER) {
      if (dim[0] == 0) {
        corigin_user[cdim1] *= yscale;
        corigin_user[cdim2] *= zscale;
        cradmin_user *= yscale;     // radii are scaled by first non-axis dim
        cradmax_user *= yscale;
      } else if (dim[0] == 1) {
        corigin_user[cdim1] *= xscale;
        corigin_user[cdim2] *= zscale;
        cradmin_user *= xscale;
        cradmax_user *= xscale;
      } else {
        corigin_user[cdim1] *= xscale;
        corigin_user[cdim2] *= yscale;
        cradmin_user *= xscale;
        cradmax_user *= xscale;
      }
    }
  }

  // initialize chunk vector and per-chunk info

  nmax = 0;
  chunk = nullptr;
  nmaxint = -1;
  ichunk = nullptr;
  exclude = nullptr;

  nchunk = 0;
  chunk_volume_scalar = 1.0;
  chunk_volume_vec = nullptr;
  coord = nullptr;
  chunkID = nullptr;

  // computeflag = 1 if this compute might invoke another compute
  // during assign_chunk_ids()

  if (which == ArgInfo::COMPUTE || which == ArgInfo::FIX || which == ArgInfo::VARIABLE) computeflag = 1;
  else computeflag = 0;

  // other initializations

  invoked_setup = -1;
  invoked_ichunk = -1;

  id_fix = nullptr;
  fixstore = nullptr;

  if (compress) hash = new std::map<tagint,int>();
  else hash = nullptr;

  maxvar = 0;
  varatom = nullptr;

  lockcount = 0;
  lockfix = nullptr;

  if (which == ArgInfo::MOLECULE) molcheck = 1;
  else molcheck = 0;
}

/* ---------------------------------------------------------------------- */

ComputeChunkAtom::~ComputeChunkAtom()
{
  // check nfix in case all fixes have already been deleted

  if (id_fix && modify->nfix) modify->delete_fix(id_fix);
  delete [] id_fix;

  memory->destroy(chunk);
  memory->destroy(ichunk);
  memory->destroy(exclude);
  memory->destroy(chunk_volume_vec);
  memory->destroy(coord);
  memory->destroy(chunkID);

  delete [] idregion;
  delete [] cfvid;
  delete hash;

  memory->destroy(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeChunkAtom::init()
{
  // set and check validity of region

  if (regionflag) {
    int iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for compute chunk/atom does not exist");
    region = domain->regions[iregion];
  }

  // set compute,fix,variable

  if (which == ArgInfo::COMPUTE) {
    int icompute = modify->find_compute(cfvid);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for compute chunk/atom does not exist");
    cchunk = modify->compute[icompute];
  } else if (which == ArgInfo::FIX) {
    int ifix = modify->find_fix(cfvid);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for compute chunk/atom does not exist");
    fchunk = modify->fix[ifix];
  } else if (which == ArgInfo::VARIABLE) {
    int ivariable = input->variable->find(cfvid);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for compute chunk/atom does not exist");
    vchunk = ivariable;
  }

  // for style MOLECULE, check that no mol IDs exceed MAXSMALLINT
  // don't worry about group or optional region

  if (which == ArgInfo::MOLECULE) {
    tagint *molecule = atom->molecule;
    int nlocal = atom->nlocal;
    tagint maxone = -1;
    for (int i = 0; i < nlocal; i++)
      if (molecule[i] > maxone) maxone = molecule[i];
    tagint maxall;
    MPI_Allreduce(&maxone,&maxall,1,MPI_LMP_TAGINT,MPI_MAX,world);
    if (maxall > MAXSMALLINT)
      error->all(FLERR,"Molecule IDs too large for compute chunk/atom");
  }

  // for binning, if nchunkflag not already set, set it to ONCE or EVERY
  // depends on whether simulation box size is static or dynamic
  // reset invoked_setup if this is not first run and box just became static

  if (binflag && !nchunksetflag && !compress && scaleflag != REDUCED) {
    if (domain->box_change_size == 0) {
      if (nchunkflag == EVERY && invoked_setup >= 0) invoked_setup = -1;
      nchunkflag = ONCE;
    } else nchunkflag = EVERY;
  }

  // require nchunkflag = ONCE if idsflag = ONCE
  // b/c nchunk cannot change if chunk IDs are frozen
  // can't check until now since nchunkflag may have been adjusted in init()

  if (idsflag == ONCE && nchunkflag != ONCE)
    error->all(FLERR,"Compute chunk/atom ids once but nchunk is not once");

  // create/destroy fix STORE for persistent chunk IDs as needed
  // need to do this if idsflag = ONCE or locks will be used by other commands
  // need to wait until init() so that fix command(s) are in place
  //   they increment lockcount if they lock this compute
  // fixstore ID = compute-ID + COMPUTE_STORE, fix group = compute group
  // fixstore initializes all values to 0.0

  if ((idsflag == ONCE || lockcount) && !fixstore) {
    id_fix = utils::strdup(id + std::string("_COMPUTE_STORE"));
    fixstore = (FixStore *) modify->add_fix(fmt::format("{} {} STORE peratom 1 1",
                                                        id_fix, group->names[igroup]));
  }

  if ((idsflag != ONCE && !lockcount) && fixstore) {
    modify->delete_fix(id_fix);
    fixstore = nullptr;
  }
}

/* ----------------------------------------------------------------------
   invoke setup_chunks and/or compute_ichunk if only done ONCE
   so that nchunks or chunk IDs are assigned when this compute was specified
     as opposed to first times compute_peratom() or compute_ichunk() is called
------------------------------------------------------------------------- */

void ComputeChunkAtom::setup()
{
  if (nchunkflag == ONCE) setup_chunks();
  if (idsflag == ONCE) compute_ichunk();
  else invoked_ichunk = -1;
}

/* ----------------------------------------------------------------------
   only called by classes that use per-atom computes in standard way
     dump, variable, thermo output, other computes, etc
   not called by fix chunk or compute chunk commands
     they invoke setup_chunks() and compute_ichunk() directly
------------------------------------------------------------------------- */

void ComputeChunkAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow floating point chunk vector if necessary

  if (atom->nmax > nmax) {
    memory->destroy(chunk);
    nmax = atom->nmax;
    memory->create(chunk,nmax,"chunk/atom:chunk");
    vector_atom = chunk;
  }

  setup_chunks();
  compute_ichunk();

  // copy integer indices into floating-point chunk vector

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) chunk[i] = ichunk[i];
}


/* ----------------------------------------------------------------------
   to return the number of chunks, we first need to make certain
   that compute_peratom() has been called.
------------------------------------------------------------------------- */
double ComputeChunkAtom::compute_scalar()
{
  if (invoked_peratom != update->ntimestep)
    compute_peratom();
  invoked_scalar = update->ntimestep;

  return (scalar = nchunk);
}

/* ----------------------------------------------------------------------
   set lock, so that nchunk will not change from startstep to stopstep
   called by fix for duration of time it requires lock
   OK if called by multiple fix commands
     error if all callers do not have same duration
     last caller holds the lock, so it can also unlock
   stopstep can be positive for final step of finite-size time window
   or can be -1 for infinite-size time window
------------------------------------------------------------------------- */

void ComputeChunkAtom::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  if (lockfix == nullptr) {
    lockfix = fixptr;
    lockstart = startstep;
    lockstop = stopstep;
    return;
  }

  if (startstep != lockstart || stopstep != lockstop)
    error->all(FLERR,"Two fix commands using "
               "same compute chunk/atom command in incompatible ways");

  // set lock to last calling Fix, since it will be last to unlock()

  lockfix = fixptr;
}

/* ----------------------------------------------------------------------
   unset lock
   can only be done by fix command that holds the lock
------------------------------------------------------------------------- */

void ComputeChunkAtom::unlock(Fix *fixptr)
{
  if (fixptr != lockfix) return;
  lockfix = nullptr;
}

/* ----------------------------------------------------------------------
   assign chunk IDs from 1 to Nchunk to every atom, or 0 if not in chunk
------------------------------------------------------------------------- */

void ComputeChunkAtom::compute_ichunk()
{
  int i;

  // skip if already done on this step

  if (invoked_ichunk == update->ntimestep) return;

  // if old IDs persist via storage in fixstore, then just retrieve them
  // yes if idsflag = ONCE, and already done once
  //   or if idsflag = NFREQ and lock is in place and are on later timestep
  // else proceed to recalculate per-atom chunk assignments

  const int nlocal = atom->nlocal;
  int restore = 0;
  if (idsflag == ONCE && invoked_ichunk >= 0) restore = 1;
  if (idsflag == NFREQ && lockfix && update->ntimestep > lockstart) restore = 1;

  if (restore) {
    invoked_ichunk = update->ntimestep;
    double *vstore = fixstore->vstore;
    for (i = 0; i < nlocal; i++) ichunk[i] = static_cast<int> (vstore[i]);
    return;
  }

  // assign chunk IDs to atoms
  // will exclude atoms not in group or in optional region
  // already invoked if this is same timestep as last setup_chunks()
  // however, when between runs or using rerun, we need it again.

  if ((update->ntimestep > invoked_setup) || (invoked_ichunk < 0)) assign_chunk_ids();

  invoked_ichunk = update->ntimestep;

  // compress chunk IDs via hash of the original uncompressed IDs
  // also apply discard rule except for binning styles which already did

  if (compress) {
    if (binflag) {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (hash->find(ichunk[i]) == hash->end()) exclude[i] = 1;
        else ichunk[i] = hash->find(ichunk[i])->second;
      }
    } else if (discard == NODISCARD) {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (hash->find(ichunk[i]) == hash->end()) ichunk[i] = nchunk;
        else ichunk[i] = hash->find(ichunk[i])->second;
      }
    } else {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (hash->find(ichunk[i]) == hash->end()) exclude[i] = 1;
        else ichunk[i] = hash->find(ichunk[i])->second;
      }
    }

  // else if no compression apply discard rule by itself

  } else {
    if (discard == NODISCARD) {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (ichunk[i] < 1 || ichunk[i] > nchunk) ichunk[i] = nchunk;;
      }
    } else {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (ichunk[i] < 1 || ichunk[i] > nchunk) exclude[i] = 1;
      }
    }
  }

  // set ichunk = 0 for excluded atoms
  // this should set any ichunk values which have not yet been set

  for (i = 0; i < nlocal; i++)
    if (exclude[i]) ichunk[i] = 0;

  // if newly calculated IDs need to persist, store them in fixstore
  // yes if idsflag = ONCE or idsflag = NFREQ and lock is in place

  if (idsflag == ONCE || (idsflag == NFREQ && lockfix)) {
    double *vstore = fixstore->vstore;
    for (i = 0; i < nlocal; i++) vstore[i] = ichunk[i];
  }

  // one-time check if which = MOLECULE and
  // any chunks do not contain all atoms in the molecule

  if (molcheck) {
    check_molecules();
    molcheck = 0;
  }
}

/* ----------------------------------------------------------------------
   setup chunks
   return nchunk = # of chunks
     all atoms will be assigned a chunk ID from 1 to Nchunk, or 0
   also setup any internal state needed to quickly assign atoms to chunks
   called from compute_peratom() and also directly from
     fix chunk and compute chunk commands
------------------------------------------------------------------------- */

int ComputeChunkAtom::setup_chunks()
{
  if (invoked_setup == update->ntimestep) return nchunk;

  // check if setup needs to be done
  // no if lock is in place
  // no if nchunkflag = ONCE, and already done once
  // otherwise yes
  // even if no, check if need to re-compute bin volumes
  //   so that fix ave/chunk can do proper density normalization

  int flag = 0;
  if (lockfix) flag = 1;
  if (nchunkflag == ONCE && invoked_setup >= 0) flag = 1;

  if (flag) {
    if (binflag && scaleflag == REDUCED && domain->box_change_size)
      bin_volumes();
    return nchunk;
  }

  invoked_setup = update->ntimestep;

  // assign chunk IDs to atoms
  // will exclude atoms not in group or in optional region
  // for binning styles, need to setup bins and their volume first
  //   else chunk_volume_scalar = entire box volume
  // IDs are needed to scan for max ID and for compress()

  if (binflag) {
    if (which == ArgInfo::BIN1D || which == ArgInfo::BIN2D
        || which == ArgInfo::BIN3D)
      nchunk = setup_xyz_bins();
    else if (which == ArgInfo::BINSPHERE) nchunk = setup_sphere_bins();
    else if (which == ArgInfo::BINCYLINDER) nchunk = setup_cylinder_bins();
    bin_volumes();
  } else {
    chunk_volume_scalar = domain->xprd * domain->yprd;
    if (domain->dimension == 3) chunk_volume_scalar *= domain->zprd;
  }

  assign_chunk_ids();

  // set nchunk for chunk styles other than binning
  // for styles other than TYPE, scan for max ID

  if (which == ArgInfo::TYPE) nchunk = atom->ntypes;
  else if (!binflag) {

    int nlocal = atom->nlocal;
    int hi = -1;
    for (int i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      if (ichunk[i] > hi) hi = ichunk[i];
    }

    MPI_Allreduce(&hi,&nchunk,1,MPI_INT,MPI_MAX,world);
    if (nchunk <= 0) nchunk = 1;
  }

  // apply limit setting as well as compression of chunks with no atoms
  // if limit is set, there are 3 cases:
  //   no compression, limit specified before compression, or vice versa

  if (limit && !binflag) {
    if (!compress) {
      if (limitstyle == LIMITMAX) nchunk = MIN(nchunk,limit);
      else if (limitstyle == LIMITEXACT) nchunk = limit;
    } else if (limitfirst) {
      nchunk = MIN(nchunk,limit);
    }
  }

  if (compress) compress_chunk_ids();

  if (limit && !binflag && compress) {
    if (limitstyle == LIMITMAX) nchunk = MIN(nchunk,limit);
    else if (limitstyle == LIMITEXACT) nchunk = limit;
  }

  return nchunk;
}

/* ----------------------------------------------------------------------
   assign chunk IDs for all atoms, via ichunk vector
   except excluded atoms, their chunk IDs are set to 0 later
   also set exclude vector to 0/1 for all atoms
     excluded atoms are those not in group or in optional region
   called from compute_ichunk() and setup_chunks()
------------------------------------------------------------------------- */

void ComputeChunkAtom::assign_chunk_ids()
{
  int i;

  // grow integer chunk index vector if necessary

  if (atom->nmax > nmaxint) {
    memory->destroy(ichunk);
    memory->destroy(exclude);
    nmaxint = atom->nmax;
    memory->create(ichunk,nmaxint,"chunk/atom:ichunk");
    memory->create(exclude,nmaxint,"chunk/atom:exclude");
  }

  // update region if necessary

  if (regionflag) region->prematch();

  // exclude = 1 if atom is not assigned to a chunk
  // exclude atoms not in group or not in optional region

  double **x = atom->x;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  if (regionflag) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit &&
          region->match(x[i][0],x[i][1],x[i][2])) exclude[i] = 0;
      else exclude[i] = 1;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) exclude[i] = 0;
      else exclude[i] = 1;
    }
  }

  // set ichunk to style value for included atoms
  // binning styles apply discard rule, others do not yet

  if (binflag) {
    if (which == ArgInfo::BIN1D) atom2bin1d();
    else if (which == ArgInfo::BIN2D) atom2bin2d();
    else if (which == ArgInfo::BIN3D) atom2bin3d();
    else if (which == ArgInfo::BINSPHERE) atom2binsphere();
    else if (which == ArgInfo::BINCYLINDER) atom2bincylinder();

  } else if (which == ArgInfo::TYPE) {
    int *type = atom->type;
    for (i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      ichunk[i] = type[i];
    }

  } else if (which == ArgInfo::MOLECULE) {
    tagint *molecule = atom->molecule;
    for (i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      ichunk[i] = static_cast<int> (molecule[i]);
    }

  } else if (which == ArgInfo::COMPUTE) {
    if (!(cchunk->invoked_flag & Compute::INVOKED_PERATOM)) {
      cchunk->compute_peratom();
      cchunk->invoked_flag |= Compute::INVOKED_PERATOM;
    }

    if (argindex == 0) {
      double *vec = cchunk->vector_atom;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (vec[i]);
      }
    } else {
      double **array = cchunk->array_atom;
      int argm1 = argindex - 1;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (array[i][argm1]);
      }
    }

  } else if (which == ArgInfo::FIX) {
    if (update->ntimestep % fchunk->peratom_freq)
      error->all(FLERR,"Fix used in compute chunk/atom not "
                 "computed at compatible time");

    if (argindex == 0) {
      double *vec = fchunk->vector_atom;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (vec[i]);
      }
    } else {
      double **array = fchunk->array_atom;
      int argm1 = argindex - 1;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (array[i][argm1]);
      }
    }

  } else if (which == ArgInfo::VARIABLE) {
    if (atom->nmax > maxvar) {
      maxvar = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom,maxvar,"chunk/atom:varatom");
    }

    input->variable->compute_atom(vchunk,igroup,varatom,1,0);
    for (i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      ichunk[i] = static_cast<int> (varatom[i]);
    }
  }
}

/* ----------------------------------------------------------------------
   compress chunk IDs currently assigned to atoms across all processors
     by removing those with no atoms assigned
   current assignment excludes atoms not in group or in optional region
   current Nchunk = max ID
   operation:
     use hash to store list of populated IDs that I own
     add new IDs to populated lists communicated from all other procs
     final hash has global list of populated ideas
   reset Nchunk = length of global list
   called by setup_chunks() when setting Nchunk
   remapping of chunk IDs to smaller Nchunk occurs later in compute_ichunk()
------------------------------------------------------------------------- */

void ComputeChunkAtom::compress_chunk_ids()
{
  hash->clear();

  // put my IDs into hash

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;
    if (hash->find(ichunk[i]) == hash->end()) (*hash)[ichunk[i]] = 0;
  }

  // n = # of my populated IDs
  // nall = n summed across all procs

  int n = hash->size();
  bigint nbone = n;
  bigint nball;
  MPI_Allreduce(&nbone,&nball,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // create my list of populated IDs

  int *list = nullptr;
  memory->create(list,n,"chunk/atom:list");

  n = 0;
  std::map<tagint,int>::iterator pos;
  for (pos = hash->begin(); pos != hash->end(); ++pos)
    list[n++] = pos->first;

  // if nall < 1M, just allgather all ID lists on every proc
  // else perform ring comm
  // add IDs from all procs to my hash

  if (nball <= IDMAX) {

    // setup for allgatherv

    int nprocs = comm->nprocs;
    int nall = nball;
    int *recvcounts,*displs,*listall;
    memory->create(recvcounts,nprocs,"chunk/atom:recvcounts");
    memory->create(displs,nprocs,"chunk/atom:displs");
    memory->create(listall,nall,"chunk/atom:listall");

    MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

    displs[0] = 0;
    for (int iproc = 1; iproc < nprocs; iproc++)
      displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

    // allgatherv acquires list of populated IDs from all procs

    MPI_Allgatherv(list,n,MPI_INT,listall,recvcounts,displs,MPI_INT,world);

    // add all unique IDs in listall to my hash

    for (int i = 0; i < nall; i++)
      if (hash->find(listall[i]) == hash->end()) (*hash)[listall[i]] = 0;

    // clean up

    memory->destroy(recvcounts);
    memory->destroy(displs);
    memory->destroy(listall);

  } else {
    comm->ring(n,sizeof(int),list,1,idring,nullptr,(void *)this,0);
  }

  memory->destroy(list);

  // nchunk = length of hash containing populated IDs from all procs

  nchunk = hash->size();

  // reset hash value of each original chunk ID to ordered index
  //   ordered index = new compressed chunk ID (1 to Nchunk)
  //   leverages fact that map stores keys in ascending order
  // also allocate and set chunkID = list of original chunk IDs
  //   used by fix ave/chunk and compute property/chunk

  memory->destroy(chunkID);
  memory->create(chunkID,nchunk,"chunk/atom:chunkID");

  n = 0;
  for (pos = hash->begin(); pos != hash->end(); ++pos) {
    chunkID[n] = pos->first;
    (*hash)[pos->first] = ++n;
  }
}

/* ----------------------------------------------------------------------
   callback from comm->ring()
   cbuf = list of N chunk IDs from another proc
   loop over the list, add each to my hash
   hash ends up storing all unique IDs across all procs
------------------------------------------------------------------------- */

void ComputeChunkAtom::idring(int n, char *cbuf, void *ptr)
{
  ComputeChunkAtom *cptr = (ComputeChunkAtom *)ptr;
  tagint *list = (tagint *) cbuf;
  std::map<tagint,int> *hash = cptr->hash;
  for (int i = 0; i < n; i++) (*hash)[list[i]] = 0;
}

/* ----------------------------------------------------------------------
   one-time check for which = MOLECULE to check
     if each chunk contains all atoms in the molecule
   issue warning if not
   note that this check is without regard to discard rule
   if discard == NODISCARD, there is no easy way to check that all
     atoms in an out-of-bounds molecule were added to a chunk,
     some could have been excluded by group or region, others not
------------------------------------------------------------------------- */

void ComputeChunkAtom::check_molecules()
{
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  int flag = 0;

  if (!compress) {
    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] > 0 && molecule[i] <= nchunk &&
          ichunk[i] == 0) flag = 1;
    }
  } else {
    int molid;
    for (int i = 0; i < nlocal; i++) {
      molid = static_cast<int> (molecule[i]);
      if (hash->find(molid) != hash->end() && ichunk[i] == 0) flag = 1;
    }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,
                   "One or more chunks do not contain all atoms in molecule");
}

/* ----------------------------------------------------------------------
   setup xyz spatial bins and their extent and coordinates
   return nbins = # of bins, will become # of chunks
   called from setup_chunks()
------------------------------------------------------------------------- */

int ComputeChunkAtom::setup_xyz_bins()
{
  int i,j,k,m,n,idim;
  double lo,hi,coord1,coord2;

  // lo = bin boundary immediately below boxlo or minvalue
  // hi = bin boundary immediately above boxhi or maxvalue
  // allocate and initialize arrays based on new bin count

  double binlo[3],binhi[3];
  if (scaleflag == REDUCED) {
    binlo[0] = domain->boxlo_lamda[0];
    binlo[1] = domain->boxlo_lamda[1];
    binlo[2] = domain->boxlo_lamda[2];
    binhi[0] = domain->boxhi_lamda[0];
    binhi[1] = domain->boxhi_lamda[1];
    binhi[2] = domain->boxhi_lamda[2];
  } else {
    binlo[0] = domain->boxlo[0];
    binlo[1] = domain->boxlo[1];
    binlo[2] = domain->boxlo[2];
    binhi[0] = domain->boxhi[0];
    binhi[1] = domain->boxhi[1];
    binhi[2] = domain->boxhi[2];
  }

  if (minflag[0] == COORD) binlo[0] = minvalue[0];
  if (minflag[1] == COORD) binlo[1] = minvalue[1];
  if (minflag[2] == COORD) binlo[2] = minvalue[2];
  if (maxflag[0] == COORD) binhi[0] = maxvalue[0];
  if (maxflag[1] == COORD) binhi[1] = maxvalue[1];
  if (maxflag[2] == COORD) binhi[2] = maxvalue[2];

  int nbins = 1;

  for (m = 0; m < ndim; m++) {
    idim = dim[m];
    if (originflag[m] == LOWER) origin[m] = binlo[idim];
    else if (originflag[m] == UPPER) origin[m] = binhi[idim];
    else if (originflag[m] == CENTER)
      origin[m] = 0.5 * (binlo[idim] + binhi[idim]);

    if (origin[m] < binlo[idim]) {
      n = static_cast<int> ((binlo[idim] - origin[m]) * invdelta[m]);
      lo = origin[m] + n*delta[m];
    } else {
      n = static_cast<int> ((origin[m] - binlo[idim]) * invdelta[m]);
      lo = origin[m] - n*delta[m];
      if (lo > binlo[idim]) lo -= delta[m];
    }
    if (origin[m] < binhi[idim]) {
      n = static_cast<int> ((binhi[idim] - origin[m]) * invdelta[m]);
      hi = origin[m] + n*delta[m];
      if (hi < binhi[idim]) hi += delta[m];
    } else {
      n = static_cast<int> ((origin[m] - binhi[idim]) * invdelta[m]);
      hi = origin[m] - n*delta[m];
    }

    if (lo > hi) error->all(FLERR,"Invalid bin bounds in compute chunk/atom");

    offset[m] = lo;
    nlayers[m] = static_cast<int> ((hi-lo) * invdelta[m] + 0.5);
    nbins *= nlayers[m];
  }

  // allocate and set bin coordinates

  memory->destroy(coord);
  memory->create(coord,nbins,ndim,"chunk/atom:coord");

  if (ndim == 1) {
    for (i = 0; i < nlayers[0]; i++)
      coord[i][0] = offset[0] + (i+0.5)*delta[0];
  } else if (ndim == 2) {
    m = 0;
    for (i = 0; i < nlayers[0]; i++) {
      coord1 = offset[0] + (i+0.5)*delta[0];
      for (j = 0; j < nlayers[1]; j++) {
        coord[m][0] = coord1;
        coord[m][1] = offset[1] + (j+0.5)*delta[1];
        m++;
      }
    }
  } else if (ndim == 3) {
    m = 0;
    for (i = 0; i < nlayers[0]; i++) {
      coord1 = offset[0] + (i+0.5)*delta[0];
      for (j = 0; j < nlayers[1]; j++) {
        coord2 = offset[1] + (j+0.5)*delta[1];
        for (k = 0; k < nlayers[2]; k++) {
          coord[m][0] = coord1;
          coord[m][1] = coord2;
          coord[m][2] = offset[2] + (k+0.5)*delta[2];
          m++;
        }
      }
    }
  }

  return nbins;
}

/* ----------------------------------------------------------------------
   setup spherical spatial bins and their single coordinate
   return nsphere = # of bins, will become # of chunks
   called from setup_chunks()
------------------------------------------------------------------------- */

int ComputeChunkAtom::setup_sphere_bins()
{
  // convert sorigin_user to sorigin
  // sorigin,srad are always in box units, for orthogonal or triclinic domains
  // lamda2x works for either orthogonal or triclinic

  if (scaleflag == REDUCED) {
    domain->lamda2x(sorigin_user,sorigin);
    sradmin = sradmin_user * (domain->boxhi[0]-domain->boxlo[0]);
    sradmax = sradmax_user * (domain->boxhi[0]-domain->boxlo[0]);
  } else {
    sorigin[0] = sorigin_user[0];
    sorigin[1] = sorigin_user[1];
    sorigin[2] = sorigin_user[2];
    sradmin = sradmin_user;
    sradmax = sradmax_user;
  }

  // if pbcflag set, sradmax must be < 1/2 box in any periodic dim
  // treat orthogonal and triclinic the same
  // check every time bins are created

  if (pbcflag) {
    double *prd_half = domain->prd_half;
    int *periodicity = domain->periodicity;
    int flag = 0;
    if (periodicity[0] && sradmax > prd_half[0]) flag = 1;
    if (periodicity[1] && sradmax > prd_half[1]) flag = 1;
    if (domain->dimension == 3 &&
        periodicity[2] && sradmax > prd_half[2]) flag = 1;
    if (flag)
      error->all(FLERR,"Compute chunk/atom bin/sphere radius "
                 "is too large for periodic box");
  }

  sinvrad = nsbin / (sradmax-sradmin);

  // allocate and set bin coordinates
  // coord = midpt of radii for a spherical shell

  memory->destroy(coord);
  memory->create(coord,nsbin,1,"chunk/atom:coord");

  double rlo,rhi;

  for (int i = 0; i < nsbin; i++) {
    rlo = sradmin + i * (sradmax-sradmin) / nsbin;
    rhi = sradmin + (i+1) * (sradmax-sradmin) / nsbin;
    if (i == nsbin-1) rhi = sradmax;
    coord[i][0] = 0.5 * (rlo+rhi);
  }

  return nsbin;
}

/* ----------------------------------------------------------------------
   setup cylindrical spatial bins and their single coordinate
   return nsphere = # of bins, will become # of chunks
   called from setup_chunks()
------------------------------------------------------------------------- */

int ComputeChunkAtom::setup_cylinder_bins()
{
  // setup bins along cylinder axis
  // ncplane = # of axis bins

  ncplane = setup_xyz_bins();

  // convert corigin_user to corigin
  // corigin is always in box units, for orthogonal or triclinic domains
  // lamda2x works for either orthogonal or triclinic

  if (scaleflag == REDUCED) {
    domain->lamda2x(corigin_user,corigin);
    cradmin = cradmin_user * (domain->boxhi[cdim1]-domain->boxlo[cdim1]);
    cradmax = cradmax_user * (domain->boxhi[cdim1]-domain->boxlo[cdim1]);
  } else {
    corigin[cdim1] = corigin_user[cdim1];
    corigin[cdim2] = corigin_user[cdim2];
    cradmin = cradmin_user;
    cradmax = cradmax_user;
  }

  // if pbcflag set, sradmax must be < 1/2 box in any periodic non-axis dim
  // treat orthogonal and triclinic the same
  // check every time bins are created

  if (pbcflag) {
    double *prd_half = domain->prd_half;
    int *periodicity = domain->periodicity;
    int flag = 0;
    if (periodicity[cdim1] && sradmax > prd_half[cdim1]) flag = 1;
    if (periodicity[cdim2] && sradmax > prd_half[cdim2]) flag = 1;
    if (flag)
      error->all(FLERR,"Compute chunk/atom bin/cylinder radius "
                 "is too large for periodic box");
  }

  cinvrad = ncbin / (cradmax-cradmin);

  // allocate and set radial bin coordinates
  // radial coord = midpt of radii for a cylindrical shell
  // axiscoord = saved bin coords along cylndrical axis
  // radcoord = saved bin coords in radial direction

  double **axiscoord = coord;
  memory->create(coord,ncbin,1,"chunk/atom:coord");
  double **radcoord = coord;

  double rlo,rhi;

  for (int i = 0; i < ncbin; i++) {
    rlo = cradmin + i * (cradmax-cradmin) / ncbin;
    rhi = cradmin + (i+1) * (cradmax-cradmin) / ncbin;
    if (i == ncbin-1) rhi = cradmax;
    coord[i][0] = 0.5 * (rlo+rhi);
  }

  // create array of combined coords for all bins

  memory->create(coord,ncbin*ncplane,2,"chunk/atom:coord");
  int m = 0;
  for (int i = 0; i < ncbin; i++)
    for (int j = 0; j < ncplane; j++) {
      coord[m][0] = radcoord[i][0];
      coord[m][1] = axiscoord[j][0];
      m++;
    }
  memory->destroy(axiscoord);
  memory->destroy(radcoord);

  return ncbin*ncplane;
}

/* ----------------------------------------------------------------------
   calculate chunk volumes = bin volumes
   scalar if all bins have same volume
   vector if per-bin volumes are different
------------------------------------------------------------------------- */

void ComputeChunkAtom::bin_volumes()
{
  if (which == ArgInfo::BIN1D || which == ArgInfo::BIN2D
      || which == ArgInfo::BIN3D) {
    if (domain->dimension == 3)
      chunk_volume_scalar = domain->xprd * domain->yprd * domain->zprd;
    else chunk_volume_scalar = domain->xprd * domain->yprd;
    double *prd;
    if (scaleflag == REDUCED) prd = domain->prd_lamda;
    else prd = domain->prd;
    for (int m = 0; m < ndim; m++)
      chunk_volume_scalar *= delta[m]/prd[dim[m]];

  } else if (which == ArgInfo::BINSPHERE) {
    memory->destroy(chunk_volume_vec);
    memory->create(chunk_volume_vec,nchunk,"chunk/atom:chunk_volume_vec");
    double rlo,rhi,vollo,volhi;
    for (int i = 0; i < nchunk; i++) {
      rlo = sradmin + i * (sradmax-sradmin) / nsbin;
      rhi = sradmin + (i+1) * (sradmax-sradmin) / nsbin;
      if (i == nchunk-1) rhi = sradmax;
      vollo = 4.0/3.0 * MY_PI * rlo*rlo*rlo;
      volhi = 4.0/3.0 * MY_PI * rhi*rhi*rhi;
      chunk_volume_vec[i] = volhi - vollo;
    }

  } else if (which == ArgInfo::BINCYLINDER) {
    memory->destroy(chunk_volume_vec);
    memory->create(chunk_volume_vec,nchunk,"chunk/atom:chunk_volume_vec");

    // slabthick = delta of bins along cylinder axis

    double *prd;
    if (scaleflag == REDUCED) prd = domain->prd_lamda;
    else prd = domain->prd;
    double slabthick = domain->prd[dim[0]] * delta[0]/prd[dim[0]];

    // area lo/hi of concentric circles in radial direction

    int iradbin;
    double rlo,rhi,arealo,areahi;
    for (int i = 0; i < nchunk; i++) {
      iradbin = i / ncplane;
      rlo = cradmin + iradbin * (cradmax-cradmin) / ncbin;
      rhi = cradmin + (iradbin+1) * (cradmax-cradmin) / ncbin;
      if (iradbin == ncbin-1) rhi = cradmax;
      arealo = MY_PI * rlo*rlo;
      areahi = MY_PI * rhi*rhi;
      chunk_volume_vec[i] = (areahi-arealo) * slabthick;
    }
  }
}

/* ----------------------------------------------------------------------
   assign each atom to a 1d spatial bin (layer)
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2bin1d()
{
  int i,ibin;
  double *boxlo,*boxhi,*prd;
  double xremap;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int nlayer1m1 = nlayers[0] - 1;
  int periodicity = domain->periodicity[idim];

  if (periodicity) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords
  // apply discard rule

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;

    xremap = x[i][idim];
    if (periodicity) {
      if (xremap < boxlo[idim]) xremap += prd[idim];
      if (xremap >= boxhi[idim]) xremap -= prd[idim];
    }

    ibin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
    if (xremap < offset[0]) ibin--;

    if (discard == MIXED) {
      if (!minflag[idim]) ibin = MAX(ibin,0);
      else if (ibin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[idim]) ibin = MIN(ibin,nlayer1m1);
      else if (ibin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      ibin = MAX(ibin,0);
      ibin = MIN(ibin,nlayer1m1);
    } else if (ibin < 0 || ibin > nlayer1m1) {
      exclude[i] = 1;
      continue;
    }

    ichunk[i] = ibin+1;
  }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);
}

/* ----------------------------------------------------------------------
   assign each atom to a 2d spatial bin (pencil)
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2bin2d()
{
  int i,ibin,i1bin,i2bin;
  double *boxlo,*boxhi,*prd;
  double xremap,yremap;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int jdim = dim[1];
  int nlayer1m1 = nlayers[0] - 1;
  int nlayer2m1 = nlayers[1] - 1;
  int *periodicity = domain->periodicity;

  if (periodicity[idim] || periodicity[jdim]) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords
  // apply discard rule

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;

    xremap = x[i][idim];
    if (periodicity[idim]) {
      if (xremap < boxlo[idim]) xremap += prd[idim];
      if (xremap >= boxhi[idim]) xremap -= prd[idim];
    }

    i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
    if (xremap < offset[0]) i1bin--;

    if (discard == MIXED) {
      if (!minflag[idim]) i1bin = MAX(i1bin,0);
      else if (i1bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[idim]) i1bin = MIN(i1bin,nlayer1m1);
      else if (i1bin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i1bin = MAX(i1bin,0);
      i1bin = MIN(i1bin,nlayer1m1);
    } else if (i1bin < 0 || i1bin > nlayer1m1) {
      exclude[i] = 1;
      continue;
    }

    yremap = x[i][jdim];
    if (periodicity[jdim]) {
      if (yremap < boxlo[jdim]) yremap += prd[jdim];
      if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
    }

    i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
    if (yremap < offset[1]) i2bin--;

    if (discard == MIXED) {
      if (!minflag[jdim]) i2bin = MAX(i2bin,0);
      else if (i2bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[jdim]) i2bin = MIN(i2bin,nlayer2m1);
      else if (i2bin > nlayer2m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i2bin = MAX(i2bin,0);
      i2bin = MIN(i2bin,nlayer2m1);
    } else if (i2bin < 0 || i2bin > nlayer2m1) {
      exclude[i] = 1;
      continue;
    }

    ibin = i1bin*nlayers[1] + i2bin;
    ichunk[i] = ibin+1;
  }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);
}

/* ----------------------------------------------------------------------
   assign each atom to a 3d spatial bin (brick)
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2bin3d()
{
  int i,ibin,i1bin,i2bin,i3bin;
  double *boxlo,*boxhi,*prd;
  double xremap,yremap,zremap;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int jdim = dim[1];
  int kdim = dim[2];
  int nlayer1m1 = nlayers[0] - 1;
  int nlayer2m1 = nlayers[1] - 1;
  int nlayer3m1 = nlayers[2] - 1;
  int *periodicity = domain->periodicity;

  if (periodicity[idim] || periodicity[jdim] || periodicity[kdim]) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords
  // apply discard rule

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;

    xremap = x[i][idim];
    if (periodicity[idim]) {
      if (xremap < boxlo[idim]) xremap += prd[idim];
      if (xremap >= boxhi[idim]) xremap -= prd[idim];
    }

    i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
    if (xremap < offset[0]) i1bin--;

    if (discard == MIXED) {
      if (!minflag[idim]) i1bin = MAX(i1bin,0);
      else if (i1bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[idim]) i1bin = MIN(i1bin,nlayer1m1);
      else if (i1bin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i1bin = MAX(i1bin,0);
      i1bin = MIN(i1bin,nlayer1m1);
    } else if (i1bin < 0 || i1bin > nlayer1m1) {
      exclude[i] = 1;
      continue;
    }

    yremap = x[i][jdim];
    if (periodicity[jdim]) {
      if (yremap < boxlo[jdim]) yremap += prd[jdim];
      if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
    }

    i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
    if (yremap < offset[1]) i2bin--;

    if (discard == MIXED) {
      if (!minflag[jdim]) i2bin = MAX(i2bin,0);
      else if (i2bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[jdim]) i2bin = MIN(i2bin,nlayer2m1);
      else if (i2bin > nlayer2m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i2bin = MAX(i2bin,0);
      i2bin = MIN(i2bin,nlayer2m1);
    } else if (i2bin < 0 || i2bin > nlayer2m1) {
      exclude[i] = 1;
      continue;
    }

    zremap = x[i][kdim];
    if (periodicity[kdim]) {
      if (zremap < boxlo[kdim]) zremap += prd[kdim];
      if (zremap >= boxhi[kdim]) zremap -= prd[kdim];
    }

    i3bin = static_cast<int> ((zremap - offset[2]) * invdelta[2]);
    if (zremap < offset[2]) i3bin--;

    if (discard == MIXED) {
      if (!minflag[kdim]) i3bin = MAX(i3bin,0);
      else if (i3bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[kdim]) i3bin = MIN(i3bin,nlayer3m1);
      else if (i3bin > nlayer3m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i3bin = MAX(i3bin,0);
      i3bin = MIN(i3bin,nlayer3m1);
    } else if (i3bin < 0 || i3bin > nlayer3m1) {
      exclude[i] = 1;
      continue;
    }

    ibin = i1bin*nlayers[1]*nlayers[2] + i2bin*nlayers[2] + i3bin;
    ichunk[i] = ibin+1;
  }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);
}

/* ----------------------------------------------------------------------
   assign each atom to a spherical bin
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2binsphere()
{
  int i,ibin;
  double dx,dy,dz,r;
  double xremap,yremap,zremap;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;
  double *prd_half = domain->prd_half;
  int *periodicity = domain->periodicity;

  // remap each atom's relevant coords back into box via PBC if necessary
  // apply discard rule based on rmin and rmax

  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;

    xremap = x[i][0];
    if (periodicity[0]) {
      while (xremap < boxlo[0]) {xremap += prd[0];}
      while (xremap >= boxhi[0]) {xremap -= prd[0];}
    }
    yremap = x[i][1];
    if (periodicity[1]) {
      while (yremap < boxlo[1]) {yremap += prd[1];}
      while (yremap >= boxhi[1]) {yremap -= prd[1];}
    }
    zremap = x[i][2];
    if (periodicity[2]) {
      while (zremap < boxlo[2]) {zremap += prd[2];}
      while (zremap >= boxhi[2]) {zremap -= prd[2];}
    }

    dx = xremap - sorigin[0];
    dy = yremap - sorigin[1];
    dz = zremap - sorigin[2];

    // if requested, apply PBC to distance from sphere center
    // treat orthogonal and triclinic the same
    //   with dx,dy,dz = lengths independent of each other
    // so do not use domain->minimum_image() which couples for triclinic

    if (pbcflag) {
      if (periodicity[0]) {
        while (fabs(dx) > prd_half[0]) {
          if (dx < 0.0) dx += prd[0];
          else dx -= prd[0];
        }
      }
      if (periodicity[1]) {
        while (fabs(dy) > prd_half[1]) {
          if (dy < 0.0) dy += prd[1];
          else dy -= prd[1];
        }
      }
      if (periodicity[2]) {
        while (fabs(dz) > prd_half[2]) {
          if (dz < 0.0) dz += prd[2];
          else dz -= prd[2];
        }
      }
    }

    r = sqrt(dx*dx + dy*dy + dz*dz);

    ibin = static_cast<int> ((r - sradmin) * sinvrad);
    if (r < sradmin) ibin--;

    if (discard == MIXED || discard == NODISCARD) {
      ibin = MAX(ibin,0);
      ibin = MIN(ibin,nchunk-1);
    } else if (ibin < 0 || ibin >= nchunk) {
      exclude[i] = 1;
      continue;
    }

    ichunk[i] = ibin+1;
  }
}

/* ----------------------------------------------------------------------
   assign each atom to a cylindrical bin
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2bincylinder()
{
  int i,rbin,kbin;
  double d1,d2,r;
  double remap1,remap2;

  // first use atom2bin1d() to bin all atoms along cylinder axis

  atom2bin1d();

  // now bin in radial direction
  // kbin = bin along cylinder axis
  // rbin = bin in radial direction

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;
  double *prd_half = domain->prd_half;
  int *periodicity = domain->periodicity;

  // remap each atom's relevant coords back into box via PBC if necessary
  // apply discard rule based on rmin and rmax

  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;
    kbin = ichunk[i] - 1;

    remap1 = x[i][cdim1];
    if (periodicity[cdim1]) {
      if (remap1 < boxlo[cdim1]) remap1 += prd[cdim1];
      if (remap1 >= boxhi[cdim1]) remap1 -= prd[cdim1];
    }
    remap2 = x[i][cdim2];
    if (periodicity[cdim2]) {
      if (remap2 < boxlo[cdim2]) remap2 += prd[cdim2];
      if (remap2 >= boxhi[cdim2]) remap2 -= prd[cdim2];
    }

    d1 = remap1 - corigin[cdim1];
    d2 = remap2 - corigin[cdim2];

    // if requested, apply PBC to distance from cylinder axis
    // treat orthogonal and triclinic the same
    //   with d1,d2 = lengths independent of each other

    if (pbcflag) {
      if (periodicity[cdim1]) {
        if (fabs(d1) > prd_half[cdim1]) {
          if (d1 < 0.0) d1 += prd[cdim1];
          else d1 -= prd[cdim1];
        }
      }
      if (periodicity[cdim2]) {
        if (fabs(d2) > prd_half[cdim2]) {
          if (d2 < 0.0) d2 += prd[cdim2];
          else d2 -= prd[cdim2];
        }
      }
    }

    r = sqrt(d1*d1 + d2*d2);

    rbin = static_cast<int> ((r - cradmin) * cinvrad);
    if (r < cradmin) rbin--;

    if (discard == MIXED || discard == NODISCARD) {
      rbin = MAX(rbin,0);
      rbin = MIN(rbin,ncbin-1);
    } else if (rbin < 0 || rbin >= ncbin) {
      exclude[i] = 1;
      continue;
    }

    // combine axis and radial bin indices to set ichunk

    ichunk[i] = rbin*ncplane + kbin + 1;
  }
}

/* ----------------------------------------------------------------------
   process args for one dimension of binning info
   set dim, originflag, origin, delta
------------------------------------------------------------------------- */

void ComputeChunkAtom::readdim(int narg, char **arg, int iarg, int idim)
{
  if (narg < iarg+3) error->all(FLERR,"Illegal compute chunk/atom command");
  if (strcmp(arg[iarg],"x") == 0) dim[idim] = 0;
  else if (strcmp(arg[iarg],"y") == 0) dim[idim] = 1;
  else if (strcmp(arg[iarg],"z") == 0) dim[idim] = 2;
  else error->all(FLERR,"Illegal compute chunk/atom command");

  if (dim[idim] == 2 && domain->dimension == 2)
    error->all(FLERR,"Cannot use compute chunk/atom bin z for 2d model");

  if (strcmp(arg[iarg+1],"lower") == 0) originflag[idim] = LOWER;
  else if (strcmp(arg[iarg+1],"center") == 0) originflag[idim] = CENTER;
  else if (strcmp(arg[iarg+1],"upper") == 0) originflag[idim] = UPPER;
  else originflag[idim] = COORD;
  if (originflag[idim] == COORD)
    origin[idim] = utils::numeric(FLERR,arg[iarg+1],false,lmp);

  delta[idim] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
   just set chunkID to 0 for new atom
------------------------------------------------------------------------- */

void ComputeChunkAtom::set_arrays(int i)
{
  if (!fixstore) return;
  double *vstore = fixstore->vstore;
  vstore[i] = 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays and per-chunk arrays
   note: nchunk is actually 0 until first call
------------------------------------------------------------------------- */

double ComputeChunkAtom::memory_usage()
{
  double bytes = 2*MAX(nmaxint,0) * sizeof(int);   // ichunk,exclude
  bytes += (double)nmax * sizeof(double);                  // chunk
  bytes += (double)ncoord*nchunk * sizeof(double);         // coord
  if (compress) bytes += (double)nchunk * sizeof(int);     // chunkID
  return bytes;
}

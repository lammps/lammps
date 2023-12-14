// clang-format off
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

#include "fix_momentum_chunk.h"

#include "atom.h"
#include "compute.h"
#include "compute_chunk_atom.h"
#include "compute_com_chunk.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "modify.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Contributing author: Jiang Xiao (Hong Kong Polytechnic University)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

FixMomentumChunk::FixMomentumChunk(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  cchunk(nullptr), ccom(nullptr), cvcm(nullptr), comega(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix momentum/chunk command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix momentum/chunk command");

  id_chunk = arg[4];
  int icompute = modify->find_compute(id_chunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for fix momentum/chunk");

  id_com.clear();
  id_vcm.clear();
  id_omega.clear();

  linear = angular = rescale = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"linear") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix momentum command");
      linear = 1;
      xflag = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      yflag = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      zflag = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"angular") == 0) {
      angular = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"rescale") == 0) {
      rescale = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix momentum/chunk command");
  }

  if (linear == 0 && angular == 0)
    error->all(FLERR,"Illegal fix momentum/chunk command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1)
      error->all(FLERR,"Illegal fix momentum/chunk command");

  dynamic_group_allow = 0;
}

FixMomentumChunk::~FixMomentumChunk()
{
  if (!id_com.empty()) modify->delete_compute(id_com);
  if (!id_vcm.empty()) modify->delete_compute(id_vcm);
  if (!id_omega.empty()) modify->delete_compute(id_omega);
}

/* ---------------------------------------------------------------------- */

int FixMomentumChunk::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMomentumChunk::init()
{
  // current indices for idchunk and idcom

  int icompute = modify->find_compute(id_chunk);
  if (icompute < 0)
    error->all(FLERR,"Chunk/atom compute does not exist for fix momentum/chunk");
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->compute[icompute]);
  if (strcmp(cchunk->style,"chunk/atom") != 0)
    error->all(FLERR,"Fix momentum/chunk does not use chunk/atom compute");

  // create computes dependent on chunks

  id_com = id + id_chunk + "_com";
  icompute = modify->find_compute(id_com);
  if (icompute >= 0) modify->delete_compute(id_com);
  auto cmd = fmt::format("{} {} com/chunk {}",id_com,group->names[igroup],id_chunk);
  modify->add_compute(cmd);
  icompute = modify->find_compute(id_com);
  ccom = dynamic_cast<ComputeCOMChunk *>(modify->compute[icompute]);

  id_vcm = id + id_chunk + "_vcm";
  icompute = modify->find_compute(id_vcm);
  if (icompute >= 0) modify->delete_compute(id_vcm);
  cmd = fmt::format("{} {} vcm/chunk {}",id_vcm,group->names[igroup],id_chunk);
  modify->add_compute(cmd);
  icompute = modify->find_compute(id_vcm);
  cvcm = modify->compute[icompute];

  id_omega = id + id_chunk + "_omega";
  icompute = modify->find_compute(id_omega);
  if (icompute >= 0) modify->delete_compute(id_omega);
  cmd = fmt::format("{} {} omega/chunk {}",id_omega,group->names[igroup],id_chunk);
  modify->add_compute(cmd);
  icompute = modify->find_compute(id_omega);
  comega = modify->compute[icompute];
}

/* ---------------------------------------------------------------------- */

void FixMomentumChunk::end_of_step()
{
  // calculate per-chunk properties.
  // this will also trigger a compute/update of the chunks if needed.

  ccom->compute_array();
  cvcm->compute_array();
  comega->compute_array();

  nchunk = cchunk->nchunk;
  int *ichunk = cchunk->ichunk;
  double **com = ccom->array;
  double **vcm = cvcm->array;
  double **omega = comega->array;


  // apply removing translational and rotational velocity from atoms in each chunk

  double **v = atom->v;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  // compute per-chunk kinetic energy before momentum removal, if needed

  double *ke_chunk_old,*ke_chunk_new,*ke_chunk_local,*factor;
  if (rescale) {
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    ke_chunk_old = new double[nchunk];
    ke_chunk_local = new double[nchunk];
    memset(ke_chunk_local,0,nchunk*sizeof(double));

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int m = ichunk[i]-1;
        if (m < 0) continue;

        if (rmass)
          ke_chunk_local[m] += rmass[i] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
        else
          ke_chunk_local[m] +=  mass[type[i]] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      }
    }
    MPI_Allreduce(ke_chunk_local,ke_chunk_old,nchunk,MPI_DOUBLE,MPI_SUM,world);
  }

  if (linear) {

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int m = ichunk[i]-1;
        if (m < 0) continue;
        if (xflag) v[i][0] -= vcm[m][0];
        if (yflag) v[i][1] -= vcm[m][1];
        if (zflag) v[i][2] -= vcm[m][2];
      }
    }
  }

  if (angular) {

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    double **x = atom->x;
    imageint *image = atom->image;
    double dx,dy,dz;
    double unwrap[3];

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int m = ichunk[i]-1;
        if (m < 0) continue;
        domain->unmap(x[i],image[i],unwrap);
        dx = unwrap[0] - com[m][0];
        dy = unwrap[1] - com[m][1];
        dz = unwrap[2] - com[m][2];
        v[i][0] -= omega[m][1]*dz - omega[m][2]*dy;
        v[i][1] -= omega[m][2]*dx - omega[m][0]*dz;
        v[i][2] -= omega[m][0]*dy - omega[m][1]*dx;
      }
    }
  }

  // compute kinetic energy after momentum removal, if needed

  if (rescale) {
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    ke_chunk_new = new double[nchunk];
    factor = new double[nchunk];
    memset(ke_chunk_local,0,nchunk*sizeof(double));

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int m = ichunk[i]-1;
        if (m < 0) continue;

        if (rmass)
          ke_chunk_local[m] += rmass[i] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
        else
          ke_chunk_local[m] +=  mass[type[i]] *
            (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      }
    }
    MPI_Allreduce(ke_chunk_local,ke_chunk_new,nchunk,MPI_DOUBLE,MPI_SUM,world);

    // get scaling factors

    for (int m = 0; m < nchunk; ++m)
      factor[m] = (ke_chunk_new[0] > 0.0) ? sqrt(ke_chunk_old[m]/ke_chunk_new[m]) : 1.0;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        int m = ichunk[i]-1;
        if (m < 0) continue;
        v[i][0] *= factor[m];
        v[i][1] *= factor[m];
        v[i][2] *= factor[m];
      }
    }
    delete[] factor;
    delete[] ke_chunk_local;
    delete[] ke_chunk_old;
    delete[] ke_chunk_new;
  }
}

/* ---------------------------------------------------------------------- */

void FixMomentumChunk::post_run()
{
  modify->delete_compute(id_com);
  modify->delete_compute(id_vcm);
  modify->delete_compute(id_omega);
  id_com.clear();
  id_vcm.clear();
  id_omega.clear();
}

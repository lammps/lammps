/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "chunk_wrap.h"
#include <mpi.h>
#include <cstring>
#include <string>
#include "atom.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "modify.h"

using namespace LAMMPS_NS;

ChunkWrap::ChunkWrap(LAMMPS *lmp, const char *chunkid, const char *wrapstyle)
  : Pointers(lmp), imgdiff(NULL), x(NULL), image(NULL)
{
  maxchunk = -1;
  wrapflag = UNWRAP;
  int n = strlen(chunkid) + 1;
  idchunk = new char[n];
  strcpy(idchunk,chunkid);

  if (strcmp(wrapstyle,"unwrap") == 0) {
    wrapflag = UNWRAP;
  } else if (strcmp(wrapstyle,"atom") == 0) {
    wrapflag = ATOM;
  } else if (strcmp(wrapstyle,"chunk") == 0) {
    wrapflag = CHUNK;
  } else if (strcmp(wrapstyle,"com") == 0) {
    wrapflag = COM;
  } else {
    std::string errmsg = "Unknown chunk wrap style: ";
    errmsg += wrapstyle;
    error->all(FLERR,errmsg.c_str());
  }
}

ChunkWrap::~ChunkWrap()
{
  memory->destroy(imgdiff);
}

void ChunkWrap::init()
{
  x = atom->x;
  image = atom->image;

  // imageint value for (0,0,0)

  const imageint zeroimg = ((imageint) IMGMAX)
    | (((imageint) IMGMAX) << IMGBITS)
    | (((imageint) IMGMAX) << IMG2BITS);

  if ((wrapflag == CHUNK) || (wrapflag == COM)) {
    int icompute = modify->find_compute(idchunk);
    if (icompute < 0)
      error->all(FLERR,"Chunk/atom compute does not exist for wrap");
    ComputeChunkAtom *cca = (ComputeChunkAtom *) modify->compute[icompute];
    nchunk = cca->setup_chunks();
    cca->compute_ichunk();

    if (nchunk+1 > maxchunk) {
      memory->destroy(imgdiff);
      maxchunk = nchunk+1;
      memory->create(imgdiff, maxchunk, "chunkwrap:imgdiff");
    }

    if (wrapflag == CHUNK) {
      // find atom with the smallest atom id for each chunk,
      // and store the imageflag difference.

      tagint *localid, *refid;
      memory->create(localid, nchunk+1, "chunkwrap:localid");
      memory->create(refid, nchunk+1, "chunkwrap:refid");
      for (int i=0; i < nchunk+1; ++i) localid[i] = MAXTAGINT;

      const int nlocal = atom->nlocal;
      tagint *tag = atom->tag;
      ichunk = cca->ichunk;
      for (int i=0; i < nlocal; ++i) {
        const int mychunk = ichunk[i];
        if (mychunk > 0) localid[mychunk] = MIN(localid[mychunk],tag[i]);
      }
      MPI_Allreduce(localid,refid,nchunk+1,MPI_LMP_TAGINT,MPI_MIN,world);
      memory->destroy(localid);

      imageint *localimg;
      memory->create(localimg,nchunk+1,"chunkwrap:localimg");
      localimg[0] = 0;
      for (int i=1; i < nchunk+1; ++i) {
        const int idx = atom->map(refid[i]);
        if ((idx >= 0) && (idx < nlocal)) {
          localimg[i] = image[idx] - zeroimg;
        } else localimg[i] = 0;
      }
      MPI_Allreduce(localimg,imgdiff,nchunk+1,MPI_LMP_IMAGEINT,MPI_SUM,world);
      memory->destroy(refid);
      memory->destroy(localimg);

    } else if (wrapflag == COM) {

      // get per-chunk center of mass from unwrapped coordinates, compute
      // the imageflags for that position and store the difference to (0,0,0)

      double *masslocal, *masstotal, **comlocal, **comtotal;
      memory->create(masslocal, nchunk+1, "chunkwrap:masslocal");
      memory->create(masstotal, nchunk+1, "chunkwrap:masstotal");
      memory->create(comlocal, nchunk+1, 3, "chunkwrap:comlocal");
      memory->create(comtotal, nchunk+1, 3, "chunkwrap:comtotal");
      memset(masslocal,0,static_cast<size_t>(nchunk+1)*sizeof(double));
      memset(&comlocal[0][0],0,3*static_cast<size_t>(nchunk+1)*sizeof(double));

      int *type = atom->type;
      double *rmass = atom->rmass;
      double *mass = atom->mass;
      const int nlocal = atom->nlocal;
      double massone, unwrap[3];
      imageint myimage;
      ichunk = cca->ichunk;
      for (int i=0; i < nlocal; ++i) {
        const int mychunk = ichunk[i];
        if (mychunk > 0) {
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          domain->unmap(x[i],image[i],unwrap);
          comlocal[mychunk][0] += unwrap[0] * massone;
          comlocal[mychunk][1] += unwrap[1] * massone;
          comlocal[mychunk][2] += unwrap[2] * massone;
          masslocal[mychunk] += massone;
        }
      }
      MPI_Allreduce(masslocal,masstotal,nchunk+1,MPI_DOUBLE,MPI_SUM,world);
      MPI_Allreduce(&comlocal[0][0],&comtotal[0][0],3*(nchunk+1),MPI_DOUBLE,MPI_SUM,world);
      memory->destroy(masslocal);
      memory->destroy(comlocal);
      for (int i=1; i < nchunk+1; ++i) {
        myimage = zeroimg;
        if (masstotal[i] > 0.0) {
          comtotal[i][0] /= masstotal[i];
          comtotal[i][1] /= masstotal[i];
          comtotal[i][2] /= masstotal[i];
          domain->remap(&comtotal[i][0],myimage);
        }
        imgdiff[i] = myimage - zeroimg;
      }
      memory->destroy(masstotal);
      memory->destroy(comtotal);
    }
  }
}

void ChunkWrap::wrap(int idx, double *out)
{
  if (wrapflag == ATOM) {
    out[0] = x[idx][0];
    out[1] = x[idx][1];
    out[2] = x[idx][2];
  } else if (wrapflag == UNWRAP) {
    domain->unmap(x[idx],image[idx],out);
  } else if ((wrapflag == CHUNK) || (wrapflag == COM)) {
    int id = ichunk[idx];
    if (id > 0) {
      imageint chunkimg = image[idx] - imgdiff[id];
      domain->unmap(x[idx],chunkimg,out);
    }
  } else { // we should never get here.
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
  }
}

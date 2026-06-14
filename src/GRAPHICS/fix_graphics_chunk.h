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

#ifdef FIX_CLASS
// clang-format off
FixStyle(graphics/chunk,FixGraphicsChunk);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_CHUNK_H
#define LMP_FIX_GRAPHICS_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {
class ComputeChunkAtom;
class Region;

class FixGraphicsChunk : public Fix {
 public:
  FixGraphicsChunk(class LAMMPS *, int, char **);
  ~FixGraphicsChunk() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  int image(int *&, double **&) override;

 private:
  double radius;               // global atom radius for hull inflation
  double alpha;                // triangulation algorithm parameter determining curvature
  bool has_global_radius;      // override auto-determined radius with global value
  bool smooth;                 // smooth vs flat shading
  bool clip;                   // clip point cloud to box
  int maxreplace;              // largest cluster size where atoms are replaced by icosahedra
  int maxchunk;                // size of array to store chunk-IDs for local and ghost atoms
  int *chunkid;                // copy of chunk-IDs
  double *atomrad;             // copy of atomradius
  double mindist;              // minimum distance between positions in chunks
  char *id_chunk;              // compute chunk/atom ID
  ComputeChunkAtom *cchunk;    // pointer to chunk compute
  char *id_region;             // region ID
  Region *region;              // pointer to region

  int numobjs;
  int *imgobjs;
  double **imgparms;

};
}    // namespace LAMMPS_NS
#endif
#endif

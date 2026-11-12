// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
------------------------------------------------------------------------- */

#include "pair_body_rounded_polyhedron_granular_omp.h"
#include "pair_body_rounded_polyhedron_omp_impl.h"

#include "granular_model.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace Granular_NS;

namespace LAMMPS_NS {
template class PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedronGranular>;
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedronGranularOMP::PairBodyRoundedPolyhedronGranularOMP(LAMMPS *lmp) :
    PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedronGranular>(lmp), models_stale(1)
{
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedronGranularOMP::~PairBodyRoundedPolyhedronGranularOMP()
{
  destroy_models_thr();
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranularOMP::destroy_models_thr()
{
  for (auto &models_t : models_thr)
    for (auto *model : models_t) delete model;
  models_thr.clear();
}

/* ----------------------------------------------------------------------
   the models of mixed pairs of types are made in init_one(), after
   init_style(), so the copies per thread are made at the next compute()
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranularOMP::init_style()
{
  PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedronGranular>::init_style();
  models_stale = 1;
}

/* ----------------------------------------------------------------------
   copies of all models for each thread, made in the same way as the models,
   see PairBodyRoundedPolyhedronGranular::clone_model()
------------------------------------------------------------------------- */

void PairBodyRoundedPolyhedronGranularOMP::init_thr(int nthreads)
{
  if (!models_stale && ((int) models_thr.size() == nthreads)) return;

  destroy_models_thr();
  models_thr.resize(nthreads);
  for (int t = 0; t < nthreads; t++)
    for (int m = 0; m < (int) models.size(); m++) models_thr[t].push_back(clone_model(m));
  models_stale = 0;
}

/* ----------------------------------------------------------------------
   the copy of the model of a pair of types for the calling thread
------------------------------------------------------------------------- */

GranularModel *PairBodyRoundedPolyhedronGranularOMP::contact_model(int itype, int jtype)
{
  int tid = 0;
#if defined(_OPENMP)
  tid = omp_get_thread_num();
#endif
  if (tid >= (int) models_thr.size()) return type_model(itype, jtype);
  return models_thr[tid][model_index[itype][jtype]];
}

/* ---------------------------------------------------------------------- */

double PairBodyRoundedPolyhedronGranularOMP::memory_usage()
{
  double bytes = PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedronGranular>::memory_usage();
  bytes += (double) models_thr.size() * models.size() * sizeof(GranularModel);
  return bytes;
}

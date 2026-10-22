#include "pair_pace_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"

#include <cstring>
#include <exception>

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

namespace LAMMPS_NS {
struct ACEImpl {
  ACEImpl() : basis_set(nullptr), ace(nullptr) {}
  ~ACEImpl()
  {
    delete basis_set;
    delete ace;
  }
  ACECTildeBasisSet *basis_set;
  ACERecursiveEvaluator *ace;
};
}    // namespace LAMMPS_NS

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;

// Do not need the list of elements since it is used by PairPACE::coeff only

/* ---------------------------------------------------------------------- */

PairPACEOMP::PairPACEOMP(LAMMPS *lmp) : PairPACE(lmp), ThrOMP(lmp, THR_PAIR), ace_thr()
{
  suffix_flag |= Suffix::OMP;

  // Create and resize vector of per-thread ACE evaluators
  const int nthreads = comm->nthreads;
  ace_thr.resize(nthreads, nullptr);
  basis_thr.resize(nthreads, nullptr);
}

/* ---------------------------------------------------------------------- */

PairPACEOMP::~PairPACEOMP()
{

  // Clean up resources PairPACEOMP owns
  for (ACERecursiveEvaluator *ptr : ace_thr) { delete ptr; }
  ace_thr.clear();
  for (ACECTildeBasisSet *ptr : basis_thr) { delete ptr; }
  basis_thr.clear();
}

/* ---------------------------------------------------------------------- */

void PairPACEOMP::coeff(int narg, char **arg)
{
  // Use base class to create one ACERecursiveEvaluator object
  PairPACE::coeff(narg, arg);
  // Build per-thread ACERecursiveEvaluator
  build_thread_evaluators();
}

/* ---------------------------------------------------------------------- */

void PairPACEOMP::build_thread_evaluators()
{
  const int nthreads = comm->nthreads;

  // Clean up any existing evaluators
  for (ACERecursiveEvaluator *&ptr : ace_thr) {
    delete ptr;
    ptr = nullptr;
  }
  // Clean up existing basis
  for (ACECTildeBasisSet *&ptr : basis_thr) {
    delete ptr;
    ptr = nullptr;
  }

  // Copy configuration from the base ACERecursiveEvaluator
  // Need a copy of both the evaluator and the basis to ensure no memory is shared across threads
  for (int t = 0; t < nthreads; ++t) {
    auto *ace = new ACERecursiveEvaluator();
    ACECTildeBasisSet *basis = new ACECTildeBasisSet(*aceimpl->basis_set);
    ace->set_recursive(recursive);
    ace->set_basis(*basis, 1);
    ace->element_type_mapping = aceimpl->ace->element_type_mapping;
    ace_thr[t] = ace;
    basis_thr[t] = basis;
  }
}

/* ---------------------------------------------------------------------- */

void PairPACEOMP::compute(int eflag, int vflag)
{
  if (copymode)
    ev_init(eflag, vflag, 0);
  else
    ev_init(eflag, vflag, 1);

  int i, ii, inum, jnum;
  int *ilist, *numneigh, **firstneigh;
  double **x = atom->x;
  int *type = atom->type;

  // number of atoms in cell
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;
  const int nthreads = comm->nthreads;

  // inum: length of the neighborlists list
  inum = list->inum;

  // ilist: list of "i" atoms for which neighbor lists exist
  ilist = list->ilist;

  //numneigh: the length of each these neighbor list
  numneigh = list->numneigh;

  // the pointer to the list of neighbors of "i"
  firstneigh = list->firstneigh;

  // corerep_factor can be shared as each i atom is handled by exactly one thread
  if (flag_corerep_factor && atom->nlocal > nmax_corerep) {
    memory->destroy(corerep_factor);
    nmax_corerep = atom->nlocal;
    memory->create(corerep_factor, nmax_corerep, "pace/atom:corerep_factor");
    //zeroify array
    memset(corerep_factor, 0, nmax_corerep * sizeof(*corerep_factor));
  }

  //determine the maximum number of neighbours
  int max_jnum = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jnum = numneigh[i];
    if (jnum > max_jnum) max_jnum = jnum;
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(                                               \
    eflag, eflag_either, vflag, vflag_either, eatom, vatom, x, type, nlocal, nall, newton_pair, \
    inum, ilist, numneigh, firstneigh, max_jnum, fix, ace_thr, scale, flag_corerep_factor,      \
    corerep_factor, nthreads)
#endif
  {

    int i, j, ii, jj, jnum;
    int *jlist;
    double delx, dely, delz, evdwl;
    double fij[3];

    // thread setup
    // ifrom, ito = i-atom thread limits
    // tid = thread id
    int ifrom, ito, tid;
    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    // per-thread force view
    double *const *const f_thr = thr->get_f();

    // per-thread ACE evaluator with its own basis
    ACERecursiveEvaluator *ace = ace_thr[tid];

    // per-thread memory resize
    ace->resize_neighbours_cache(max_jnum);

    for (ii = ifrom; ii < ito; ii++) {
      i = ilist[ii];
      const int itype = type[i];

      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];

      jlist = firstneigh[i];
      jnum = numneigh[i];

      try {
        ace->compute_atom(i, x, type, jnum, jlist);
      } catch (std::exception &e) {
        error->one(FLERR, e.what());
      }

      if (flag_corerep_factor) corerep_factor[i] = 1 - ace->ace_fcut;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        delx = x[j][0] - xtmp;
        dely = x[j][1] - ytmp;
        delz = x[j][2] - ztmp;

        fij[0] = scale[itype][itype] * ace->neighbours_forces(jj, 0);
        fij[1] = scale[itype][itype] * ace->neighbours_forces(jj, 1);
        fij[2] = scale[itype][itype] * ace->neighbours_forces(jj, 2);

        f_thr[i][0] += fij[0];
        f_thr[i][1] += fij[1];
        f_thr[i][2] += fij[2];
        f_thr[j][0] -= fij[0];
        f_thr[j][1] -= fij[1];
        f_thr[j][2] -= fij[2];

        // tally per-atom virial contribution
        if (vflag_either)
          ev_tally_xyz_thr(this, i, j, nlocal, newton_pair, 0.0, 0.0, fij[0], fij[1], fij[2], -delx,
                           -dely, -delz, thr);
      }    // loop over jj

      // tally energy contribution
      if (eflag_either) {
        // evdwl = energy of atom i
        evdwl = scale[itype][itype] * ace->e_atom;
        ev_tally_full_thr(this, i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0, thr);
      }
    }    // loop over ii

    // reduce thread's contributions into global arrays
    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

double PairPACEOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairPACE::memory_usage();

  return bytes;
}

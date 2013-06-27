/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_lj_long_tip4p_long_omp.h"
#include "atom.h"
#include "comm.h"
#include "math_vector.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "domain.h"

#include "suffix.h"
using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJLongTIP4PLongOMP::PairLJLongTIP4PLongOMP(LAMMPS *lmp) :
  PairLJLongTIP4PLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 1;
  cut_respa = NULL;
  newsite_thr = NULL;
  hneigh_thr = NULL;
  tip4pflag = dispersionflag = 1;
  no_virial_fdotr_compute = 1;
  single_enable = 0;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

PairLJLongTIP4PLongOMP::~PairLJLongTIP4PLongOMP()
{
  memory->destroy(hneigh_thr);
  memory->destroy(newsite_thr);
}

/* ---------------------------------------------------------------------- */

void PairLJLongTIP4PLongOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  // reallocate hneigh_thr & newsite_thr if necessary
  // initialize hneigh_thr[0] to -1 on steps when reneighboring occured
  // initialize hneigh_thr[2] to 0 every step
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh_thr);
    memory->create(hneigh_thr,nmax,"pair:hneigh_thr");
    memory->destroy(newsite_thr);
    memory->create(newsite_thr,nmax,"pair:newsite_thr");
  }

  int i;
  // tag entire list as completely invalid after a neighbor
  // list update, since that can change the order of atoms.
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh_thr[i].a = -1;

  // indicate that the coordinates for the M point need to
  // be updated. this needs to be done in every step.
  for (i = 0; i < nall; i++) hneigh_thr[i].t = 0;

  const int order1 = ewald_order&(1<<1);
  const int order6 = ewald_order&(1<<6);

  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (order6) {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,1,1>(ifrom, ito, thr);
                else eval<1,1,0,0,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,1,1>(ifrom, ito, thr);
                else eval<1,0,0,0,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,1,1>(ifrom, ito, thr);
              else eval<0,0,0,0,0,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,1,1>(ifrom, ito, thr);
                else eval<1,1,0,1,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,1,1>(ifrom, ito, thr);
                else eval<1,0,0,1,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,1,1>(ifrom, ito, thr);
              else eval<0,0,0,1,0,1,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,1,1>(ifrom, ito, thr);
                else eval<1,1,0,0,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,1,1>(ifrom, ito, thr);
                else eval<1,0,0,0,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,1,1,1>(ifrom, ito, thr);
              else eval<0,0,0,0,1,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,1,1>(ifrom, ito, thr);
                else eval<1,1,0,1,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,1,1>(ifrom, ito, thr);
                else eval<1,0,0,1,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,1,1>(ifrom, ito, thr);
              else eval<0,0,0,1,1,1,1>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,0,1>(ifrom, ito, thr);
                else eval<1,1,0,0,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,0,1>(ifrom, ito, thr);
                else eval<1,0,0,0,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,0,1>(ifrom, ito, thr);
              else eval<0,0,0,0,0,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,0,1>(ifrom, ito, thr);
                else eval<1,1,0,1,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,0,1>(ifrom, ito, thr);
                else eval<1,0,0,1,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,0,1>(ifrom, ito, thr);
              else eval<0,0,0,1,0,0,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,0,1>(ifrom, ito, thr);
                else eval<1,1,0,0,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,0,1>(ifrom, ito, thr);
                else eval<1,0,0,0,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,1,0,1>(ifrom, ito, thr);
              else eval<0,0,0,0,1,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,0,1>(ifrom, ito, thr);
                else eval<1,1,0,1,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,0,1>(ifrom, ito, thr);
                else eval<1,0,0,1,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,0,1>(ifrom, ito, thr);
              else eval<0,0,0,1,1,0,1>(ifrom, ito, thr);
            }
          }
        }
      }
    } else {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,1,0>(ifrom, ito, thr);
                else eval<1,1,0,0,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,1,0>(ifrom, ito, thr);
                else eval<1,0,0,0,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,1,0>(ifrom, ito, thr);
              else eval<0,0,0,0,0,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,1,0>(ifrom, ito, thr);
                else eval<1,1,0,1,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,1,0>(ifrom, ito, thr);
                else eval<1,0,0,1,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,1,0>(ifrom, ito, thr);
              else eval<0,0,0,1,0,1,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,1,0>(ifrom, ito, thr);
                else eval<1,1,0,0,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,1,0>(ifrom, ito, thr);
                else eval<1,0,0,0,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,1,1,0>(ifrom, ito, thr);
              else eval<0,0,0,0,1,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,1,0>(ifrom, ito, thr);
                else eval<1,1,0,1,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,1,0>(ifrom, ito, thr);
                else eval<1,0,0,1,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,1,0>(ifrom, ito, thr);
              else eval<0,0,0,1,1,1,0>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,0,0>(ifrom, ito, thr);
                else eval<1,1,0,0,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,0,0>(ifrom, ito, thr);
                else eval<1,0,0,0,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,0,0>(ifrom, ito, thr);
              else eval<0,0,0,0,0,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,0,0>(ifrom, ito, thr);
                else eval<1,1,0,1,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,0,0>(ifrom, ito, thr);
                else eval<1,0,0,1,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,0,0>(ifrom, ito, thr);
              else eval<0,0,0,1,0,0,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,0,0>(ifrom, ito, thr);
                else eval<1,1,0,0,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,0,0>(ifrom, ito, thr);
                else eval<1,0,0,0,1,0,0>(ifrom, ito, thr);
              }
            } else {
            if (force->newton_pair) eval<0,0,1,0,1,0,0>(ifrom, ito, thr);
              else eval<0,0,0,0,1,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,0,0>(ifrom, ito, thr);
                else eval<1,1,0,1,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,0,0>(ifrom, ito, thr);
                else eval<1,0,0,1,1,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,0,0>(ifrom, ito, thr);
              else eval<0,0,0,1,1,0,0>(ifrom, ito, thr);
            }
          }
        }
      }
    } 

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

void PairLJLongTIP4PLongOMP::compute_inner()
{
   // reallocate hneigh_thr & newsite_thr if necessary
  // initialize hneigh_thr[0] to -1 on steps when reneighboring occured
  // initialize hneigh_thr[2] to 0 every step
  const int nall = atom->nlocal + atom->nghost;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh_thr);
    memory->create(hneigh_thr,nmax,"pair:hneigh_thr");
    memory->destroy(newsite_thr);
    memory->create(newsite_thr,nmax,"pair:newsite_thr");
  }

  int i;
  // tag entire list as completely invalid after a neighbor
  // list update, since that can change the order of atoms.
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh_thr[i].a = -1;

  // indicate that the coordinates for the M point need to
  // be updated. this needs to be done in every step.
  for (i = 0; i < nall; i++) hneigh_thr[i].t = 0;

  const int nthreads = comm->nthreads;
  const int inum = listinner->inum;
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(0, 0, nall, 0, 0, thr);
    eval_inner(ifrom, ito, thr);

  }  // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

void PairLJLongTIP4PLongOMP::compute_middle()
{

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = listmiddle->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(0, 0, nall, 0, 0, thr);
    eval_middle(ifrom, ito, thr);

  }  // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

void PairLJLongTIP4PLongOMP::compute_outer(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;
  const int order1 = ewald_order&(1<<1);
  const int order6 = ewald_order&(1<<6);

  const int nall = atom->nlocal + atom->nghost;

  // reallocate hneigh_thr & newsite_thr if necessary
  // initialize hneigh_thr[0] to -1 on steps when reneighboring occured
  // initialize hneigh_thr[2] to 0 every step

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh_thr);
    memory->create(hneigh_thr,nmax,"pair:hneigh_thr");
    memory->destroy(newsite_thr);
    memory->create(newsite_thr,nmax,"pair:newsite_thr");
  }

  int i;
  // tag entire list as completely invalid after a neighbor
  // list update, since that can change the order of atoms.
  if (neighbor->ago == 0) {
    for (i = 0; i < nall; i++) hneigh_thr[i].a = -1;
    // indicate that the coordinates for the M point need to
    // be updated. this needs to be done only if neighbor list
    // has been updated in compute_outer
    for (i = 0; i < nall; i++) hneigh_thr[i].t = 0;
  }

  const int nthreads = comm->nthreads;
  const int inum = listouter->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (order6) {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,1,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,1,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,1,1>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,0,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,1,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,0,1>(ifrom, ito, thr);
            }
          }
        }
      }
    } else {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,1,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,1,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,1,0>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,0,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,0,0>(ifrom, ito, thr);
              }
            } else {
            if (force->newton_pair) eval_outer<0,0,1,0,1,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,0,0>(ifrom, ito, thr);
            }
          }
        }
      }
    } 

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template < const int EVFLAG, const int EFLAG,
           const int NEWTON_PAIR, const int CTABLE, const int LJTABLE, const int ORDER1, const int ORDER6 >
void PairLJLongTIP4PLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const double cut_coulsqplus = (cut_coul+2.0*qdist)*(cut_coul+2.0*qdist);

  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int n,vlist[6];
  int key;
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double r,r2inv,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double xiM[3],xjM[3],fO[3],fH[3],fd[3],v[6],xH1[3],xH2[3];// f1[3];
  dbl3_t x1,x2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;

  int ni;
  double  *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    if (itype == typeO) {
      if (hneigh_thr[i].a < 0) {
        hneigh_thr[i].a = iH1 = atom->map(atom->tag[i] + 1);
        hneigh_thr[i].b = iH2 = atom->map(atom->tag[i] + 2);
        hneigh_thr[i].t = 1;
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
      } else {
        iH1 = hneigh_thr[i].a;
        iH2 = hneigh_thr[i].b;
        if (hneigh_thr[i].t == 0) {
          hneigh_thr[i].t = 1;
          compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
        }
      }
      x1 = newsite_thr[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    offseti = offset[itype];
    lj1i = lj1[itype]; lj2i = lj2[itype]; lj3i = lj3[itype]; lj4i = lj4[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      ni = sbmask(j);
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
        
      if (rsq < cut_ljsq[itype][jtype]) {			// lj
        r2inv = 1.0/rsq;
       	if (ORDER6) {					// long-range lj
          if (!LJTABLE || rsq <= tabinnerdispsq) {
	    register double rn = r2inv*r2inv*r2inv;
	    register double x2 = g2*rsq, a2 = 1.0/x2;
	    x2 = a2*exp(-x2)*lj4i[jtype];
	    if (ni == 0) {
	      forcelj =
	        (rn*=rn)*lj1i[jtype]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	      if (EFLAG)
	        evdwl = rn*lj3i[jtype]-g6*((a2+1.0)*a2+0.5)*x2;
	    }
	    else {					// special case
	      register double f = special_lj[ni], t = rn*(1.0-f);
	      forcelj = f*(rn *= rn)*lj1i[jtype]-
	        g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[jtype];
	      if (EFLAG) 
	        evdwl = f*rn*lj3i[jtype]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[jtype];
	    }
          }
          else {                                        // table real space
            register union_int_float_t disp_t;
            disp_t.f = rsq;
            register const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            register double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            register double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              forcelj = (rn*=rn)*lj1i[jtype]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[jtype];
              if (EFLAG) evdwl = rn*lj3i[jtype]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[jtype];
            }
            else {					// special case
              register double f = special_lj[ni], t = rn*(1.0-f);
              forcelj = f*(rn *= rn)*lj1i[jtype]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[jtype]+t*lj2i[jtype];
              if (EFLAG) evdwl = f*rn*lj3i[jtype]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[jtype]+t*lj4i[jtype];
            }
          }
	}
	else {						// cut lj
	  register double rn = r2inv*r2inv*r2inv;
	  if (ni == 0) {
	    forcelj = rn*(rn*lj1i[jtype]-lj2i[jtype]);
	    if (EFLAG) evdwl = rn*(rn*lj3i[jtype]-lj4i[jtype])-offseti[jtype];
	  }
	  else {					// special case
	    register double f = special_lj[ni];
	    forcelj = f*rn*(rn*lj1i[jtype]-lj2i[jtype]);
	    if (EFLAG)
	      evdwl = f * (rn*(rn*lj3i[jtype]-lj4i[jtype])-offseti[jtype]);
	  }
        }

        forcelj *= r2inv;
	f[i][0] += delx*forcelj;
	f[i][1] += dely*forcelj;
	f[i][2] += delz*forcelj;
	f[j][0] -= delx*forcelj;
	f[j][1] -= dely*forcelj;
	f[j][2] -= delz*forcelj;

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal, /* newton_pair = */ 1,
				 evdwl,0.0,forcelj,delx,dely,delz,thr);
      }

      
      // adjust rsq and delxyz for off-site O charge(s)
      // ADDITIONAL REQEUST REQUIRED HERE!!!!!

      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) { 
	  if (jtype == typeO) {
            if (hneigh_thr[j].a < 0) {
              hneigh_thr[j].a = jH1 = atom->map(atom->tag[j] + 1);
              hneigh_thr[j].b = jH2 = atom->map(atom->tag[j] + 2);
              hneigh_thr[j].t = 1;
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
            } else {
              jH1 = hneigh_thr[j].a;
              jH2 = hneigh_thr[j].b;
              if (hneigh_thr[j].t == 0) {
                hneigh_thr[j].t = 1;
                compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
              }
            }
            x2 = newsite_thr[j];
	  } else x2 = x[j];
	  delx = x1.x - x2.x;
	  dely = x1.y - x2.y;
	  delz = x1.z - x2.z;
	  rsq = delx*delx + dely*dely + delz*delz;
        }

	// test current rsq against cutoff and compute Coulombic force
      
        if (rsq < cut_coulsq && ORDER1) {
	  r2inv = 1.0 / rsq;
	  if (!CTABLE || rsq <= tabinnersq) {
	    r = sqrt(rsq);
	    grij = g_ewald * r;
	    expm2 = exp(-grij*grij);
	    t = 1.0 / (1.0 + EWALD_P*grij);
	    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	    prefactor = qqrd2e * qtmp*q[j]/r;
	    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	    if (factor_coul < 1.0) {
	      forcecoul -= (1.0-factor_coul)*prefactor; 
	    }
	  } else {
	    union_int_float_t rsq_lookup;
	    rsq_lookup.f = rsq;
	    itable = rsq_lookup.i & ncoulmask;
	    itable >>= ncoulshiftbits;
	    fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
	    table = ftable[itable] + fraction*dftable[itable];
	    forcecoul = qtmp*q[j] * table;
	    if (factor_coul < 1.0) {
	      table = ctable[itable] + fraction*dctable[itable];
	      prefactor = qtmp*q[j] * table;
	      forcecoul -= (1.0-factor_coul)*prefactor;
	    }
	  }

	  cforce = forcecoul * r2inv;

	  //if (evflag) ev_tally(i,j,nlocal,newton_pair,
          //		       evdwl,0.0,cforce,delx,dely,delz);

	  // if i,j are not O atoms, force is applied directly
	  // if i or j are O atoms, force is on fictitious atom & partitioned
	  // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
	  // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
	  // preserves total force and torque on water molecule
	  // virial = sum(r x F) where each water's atoms are near xi and xj
	  // vlist stores 2,4,6 atoms whose forces contribute to virial

	  n = 0;
          key = 0;

	  if (itype != typeO) {
	    f[i][0] += delx * cforce;
	    f[i][1] += dely * cforce;
	    f[i][2] += delz * cforce;

            if (EVFLAG) {
              v[0] = x[i].x * delx * cforce;
              v[1] = x[i].y * dely * cforce;
              v[2] = x[i].z * delz * cforce;
              v[3] = x[i].x * dely * cforce;
              v[4] = x[i].x * delz * cforce;
              v[5] = x[i].y * delz * cforce;
            }
          vlist[n++] = i;

	  } else {
            key += 1;
            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];

	    if (EVFLAG) {
	      domain->closest_image(&x[i].x,&x[iH1].x,xH1);
	      domain->closest_image(&x[i].x,&x[iH2].x,xH2);

	      v[0] = x[i].x*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] = x[i].y*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] = x[i].z*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] = x[i].x*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] = x[i].x*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] = x[i].y*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
	    }
	    vlist[n++] = i;
	    vlist[n++] = iH1;
	    vlist[n++] = iH2;
  	  }

	  if (jtype != typeO) {
	    f[j][0] -= delx * cforce;
	    f[j][1] -= dely * cforce;
	    f[j][2] -= delz * cforce;

	    if (EVFLAG) {
	      v[0] -= x[j].x * delx * cforce;
	      v[1] -= x[j].y * dely * cforce;
	      v[2] -= x[j].z * delz * cforce;
	      v[3] -= x[j].x * dely * cforce;
	      v[4] -= x[j].x * delz * cforce;
	      v[5] -= x[j].y * delz * cforce;
            }
	    vlist[n++] = j;

	  } else {
            key += 2;

	    fd[0] = -delx*cforce;
	    fd[1] = -dely*cforce;
	    fd[2] = -delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2]; 

	    f[j][0] += fO[0];
	    f[j][1] += fO[1];
	    f[j][2] += fO[2];

	    f[jH1][0] += fH[0];
	    f[jH1][1] += fH[1];
	    f[jH1][2] += fH[2];

	    f[jH2][0] += fH[0];
	    f[jH2][1] += fH[1];
	    f[jH2][2] += fH[2];

	    if (EVFLAG) {
	      domain->closest_image(&x[j].x,&x[jH1].x,xH1);
	      domain->closest_image(&x[j].x,&x[jH2].x,xH2);

	      v[0] += x[j].x*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] += x[j].y*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] += x[j].z*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] += x[j].x*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] += x[j].x*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] += x[j].y*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
            }
      	    vlist[n++] = j;
	    vlist[n++] = jH1;
	    vlist[n++] = jH2;
	  }

	  if (EFLAG) {
	    if (!CTABLE || rsq <= tabinnersq)
	      ecoul = prefactor*erfc;
	    else {
	      table = etable[itable] + fraction*detable[itable];
	      ecoul = qtmp*q[j] * table;
	    }
	    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
 
          if (EVFLAG) ev_tally_list_thr(this,key,vlist,v,ecoul,alpha,thr);
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJLongTIP4PLongOMP::eval_inner(int iifrom, int iito, ThrData * const thr)
{
  double r, rsq, r2inv, forcecoul = 0.0, forcelj, cforce, fpair;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const int newton_pair = force->newton_pair;
  const double cut_coulsqplus = (cut_coul+2.0*qdist)*(cut_coul+2.0*qdist);

  const double cut_out_on = cut_respa[0];
  const double cut_out_off = cut_respa[1];

  const double cut_out_diff = cut_out_off - cut_out_on;
  const double cut_out_on_sq = cut_out_on*cut_out_on;
  const double cut_out_off_sq = cut_out_off*cut_out_off;

  int *jneigh, *jneighn, typei, typej, ni;
  const int order1 = (ewald_order|(ewald_off^-1))&(1<<1);
  double qri;
  vector xi, d;

  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double xiM[3],xjM[3],fO[3],fH[3],fd[3],v[6],xH1[3],xH2[3];// f1[3];
  dbl3_t x1,x2;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double  *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;

  ilist = listinner->ilist;
  numneigh = listinner->numneigh;
  firstneigh = listinner->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    if (itype == typeO) {
      if (hneigh_thr[i].a < 0) {
        hneigh_thr[i].a = iH1 = atom->map(atom->tag[i] + 1);
        hneigh_thr[i].b = iH2 = atom->map(atom->tag[i] + 2);
        hneigh_thr[i].t = 1;
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
      } else {
        iH1 = hneigh_thr[i].a;
        iH2 = hneigh_thr[i].b;
        if (hneigh_thr[i].t == 0) {
          hneigh_thr[i].t = 1;
          compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
        }
      }
      x1 = newsite_thr[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    offseti = offset[itype];
    lj1i = lj1[itype]; lj2i = lj2[itype]; lj3i = lj3[itype]; lj4i = lj4[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      ni = sbmask(j);
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
        
      if (rsq < cut_ljsq[itype][jtype] && rsq < cut_out_off_sq ) {  // lj
        r2inv = 1.0/rsq;				
	register double rn = r2inv*r2inv*r2inv;
	if (ni == 0) forcelj = rn*(rn*lj1i[jtype]-lj2i[jtype]);
	else {					// special case
	  register double f = special_lj[ni];
	  forcelj = f*rn*(rn*lj1i[jtype]-lj2i[jtype]);
	 }

        if (rsq > cut_out_on_sq) {                        // switching
          register double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          forcelj  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
        }

        forcelj *= r2inv;
	f[i][0] += delx*forcelj;
	f[i][1] += dely*forcelj;
	f[i][2] += delz*forcelj;
	f[j][0] -= delx*forcelj;
	f[j][1] -= dely*forcelj;
	f[j][2] -= delz*forcelj;
      }

      
      // adjust rsq and delxyz for off-site O charge(s)
      // ADDITIONAL REQEUST REQUIRED HERE!!!!!

      if (rsq < cut_coulsqplus && order1) {
        if (itype == typeO || jtype == typeO) { 
	  if (jtype == typeO) {
            if (hneigh_thr[j].a < 0) {
              hneigh_thr[j].a = jH1 = atom->map(atom->tag[j] + 1);
              hneigh_thr[j].b = jH2 = atom->map(atom->tag[j] + 2);
              hneigh_thr[j].t = 1;
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
            } else {
              jH1 = hneigh_thr[j].a;
              jH2 = hneigh_thr[j].b;
              if (hneigh_thr[j].t == 0) {
                hneigh_thr[j].t = 1;
                compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
              }
            }
            x2 = newsite_thr[j];
	  } else x2 = x[j];
	  delx = x1.x - x2.x;
	  dely = x1.y - x2.y;
	  delz = x1.z - x2.z;
	  rsq = delx*delx + dely*dely + delz*delz;
        }

	// test current rsq against cutoff and compute Coulombic force
      
        if (rsq < cut_coulsq && rsq < cut_out_off_sq) {
	  r2inv = 1.0 / rsq;
          qri = qqrd2e*qtmp;
          if (ni == 0) forcecoul = qri*q[j]*sqrt(r2inv);
          else {
            forcecoul = qri*q[j]*sqrt(r2inv)*special_coul[ni];
          }

          if (rsq > cut_out_on_sq) {                        // switching
            register double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
            forcecoul  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
          }

	  cforce = forcecoul * r2inv;

	  //if (evflag) ev_tally(i,j,nlocal,newton_pair,
          //		       evdwl,0.0,cforce,delx,dely,delz);

	  // if i,j are not O atoms, force is applied directly
	  // if i or j are O atoms, force is on fictitious atom & partitioned
	  // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
	  // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
	  // preserves total force and torque on water molecule
	  // virial = sum(r x F) where each water's atoms are near xi and xj
	  // vlist stores 2,4,6 atoms whose forces contribute to virial

	  if (itype != typeO) {
	    f[i][0] += delx * cforce;
	    f[i][1] += dely * cforce;
	    f[i][2] += delz * cforce;

	  } else {
            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];
          }

	  if (jtype != typeO) {
	    f[j][0] -= delx * cforce;
	    f[j][1] -= dely * cforce;
	    f[j][2] -= delz * cforce;

	  } else {
	    fd[0] = -delx*cforce;
	    fd[1] = -dely*cforce;
	    fd[2] = -delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2]; 

	    f[j][0] += fO[0];
	    f[j][1] += fO[1];
	    f[j][2] += fO[2];

	    f[jH1][0] += fH[0];
	    f[jH1][1] += fH[1];
	    f[jH1][2] += fH[2];

	    f[jH2][0] += fH[0];
	    f[jH2][1] += fH[1];
	    f[jH2][2] += fH[2];
          }
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJLongTIP4PLongOMP::eval_middle(int iifrom, int iito, ThrData * const thr)
{
  double r, rsq, r2inv, forcecoul,forcelj, cforce, fpair;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  const double cut_coulsqplus = (cut_coul+2.0*qdist)*(cut_coul+2.0*qdist);
  const int newton_pair = force->newton_pair;
  const int order1 = (ewald_order|(ewald_off^-1))&(1<<1);

  const double cut_in_off = cut_respa[0];
  const double cut_in_on = cut_respa[1];
  const double cut_out_on = cut_respa[2];
  const double cut_out_off = cut_respa[3];

  const double cut_in_diff = cut_in_on - cut_in_off;
  const double cut_out_diff = cut_out_off - cut_out_on;
  const double cut_in_off_sq = cut_in_off*cut_in_off;
  const double cut_in_on_sq = cut_in_on*cut_in_on;
  const double cut_out_on_sq = cut_out_on*cut_out_on;
  const double cut_out_off_sq = cut_out_off*cut_out_off;

  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double xiM[3],xjM[3],fO[3],fH[3],fd[3],v[6],xH1[3],xH2[3];// f1[3];
  dbl3_t x1,x2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double qri;

  int ni;
  double  *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;

  ilist = listmiddle->ilist;
  numneigh = listmiddle->numneigh;
  firstneigh = listmiddle->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    if (itype == typeO) {
      if (hneigh_thr[i].a < 0) {
        hneigh_thr[i].a = iH1 = atom->map(atom->tag[i] + 1);
        hneigh_thr[i].b = iH2 = atom->map(atom->tag[i] + 2);
        hneigh_thr[i].t = 1;
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
      } else {
        iH1 = hneigh_thr[i].a;
        iH2 = hneigh_thr[i].b;
        if (hneigh_thr[i].t == 0) {
          hneigh_thr[i].t = 1;
          compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
        }
      }
      x1 = newsite_thr[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    offseti = offset[itype];
    lj1i = lj1[itype]; lj2i = lj2[itype]; lj3i = lj3[itype]; lj4i = lj4[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      ni = sbmask(j);
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
        
      if (rsq < cut_ljsq[itype][jtype] && rsq >= cut_in_off_sq && rsq <= cut_out_off_sq ) {  // lj
        r2inv = 1.0/rsq;				
	register double rn = r2inv*r2inv*r2inv;
	if (ni == 0) forcelj = rn*(rn*lj1i[jtype]-lj2i[jtype]);
	else {					// special case
	  register double f = special_lj[ni];
	  forcelj = f*rn*(rn*lj1i[jtype]-lj2i[jtype]);
	 }

        if (rsq < cut_in_on_sq) {                                // switching
          register double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
          forcelj  *= rsw*rsw*(3.0 - 2.0*rsw);
        }
        if (rsq > cut_out_on_sq) {
          register double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          forcelj  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
        }

        forcelj *= r2inv;
	f[i][0] += delx*forcelj;
	f[i][1] += dely*forcelj;
	f[i][2] += delz*forcelj;
	f[j][0] -= delx*forcelj;
	f[j][1] -= dely*forcelj;
	f[j][2] -= delz*forcelj;
      }

      
      // adjust rsq and delxyz for off-site O charge(s)
      // ADDITIONAL REQEUST REQUIRED HERE!!!!!

      if (rsq < cut_coulsqplus && order1) {
        if (itype == typeO || jtype == typeO) { 
	  if (jtype == typeO) {
            if (hneigh_thr[j].a < 0) {
              hneigh_thr[j].a = jH1 = atom->map(atom->tag[j] + 1);
              hneigh_thr[j].b = jH2 = atom->map(atom->tag[j] + 2);
              hneigh_thr[j].t = 1;
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
            } else {
              jH1 = hneigh_thr[j].a;
              jH2 = hneigh_thr[j].b;
              if (hneigh_thr[j].t == 0) {
                hneigh_thr[j].t = 1;
                compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
              }
            }
            x2 = newsite_thr[j];
	  } else x2 = x[j];
	  delx = x1.x - x2.x;
	  dely = x1.y - x2.y;
	  delz = x1.z - x2.z;
	  rsq = delx*delx + dely*dely + delz*delz;
        }

	// test current rsq against cutoff and compute Coulombic force
      
        if (rsq < cut_coulsq &&  rsq >= cut_in_off_sq && rsq <= cut_out_off_sq) {
	  r2inv = 1.0 / rsq;
          qri = qqrd2e*qtmp;
          if (ni == 0) forcecoul = qri*q[j]*sqrt(r2inv);
          else {
            forcecoul = qri*q[j]*sqrt(r2inv)*special_coul[ni];
          }

          if (rsq < cut_in_on_sq) {                                // switching
            register double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
            forcecoul  *= rsw*rsw*(3.0 - 2.0*rsw);
          }
          if (rsq > cut_out_on_sq) {
            register double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
            forcecoul  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
          }

	  cforce = forcecoul * r2inv;

	  //if (evflag) ev_tally(i,j,nlocal,newton_pair,
          //		       evdwl,0.0,cforce,delx,dely,delz);

	  // if i,j are not O atoms, force is applied directly
	  // if i or j are O atoms, force is on fictitious atom & partitioned
	  // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
	  // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
	  // preserves total force and torque on water molecule
	  // virial = sum(r x F) where each water's atoms are near xi and xj
	  // vlist stores 2,4,6 atoms whose forces contribute to virial

	  if (itype != typeO) {
	    f[i][0] += delx * cforce;
	    f[i][1] += dely * cforce;
	    f[i][2] += delz * cforce;

	  } else {
            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];
          }

	  if (jtype != typeO) {
	    f[j][0] -= delx * cforce;
	    f[j][1] -= dely * cforce;
	    f[j][2] -= delz * cforce;

	  } else {
	    fd[0] = -delx*cforce;
	    fd[1] = -dely*cforce;
	    fd[2] = -delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2]; 

	    f[j][0] += fO[0];
	    f[j][1] += fO[1];
	    f[j][2] += fO[2];

	    f[jH1][0] += fH[0];
	    f[jH1][1] += fH[1];
	    f[jH1][2] += fH[2];

	    f[jH2][0] += fH[0];
	    f[jH2][1] += fH[1];
	    f[jH2][2] += fH[2];
          }
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template < const int EVFLAG, const int EFLAG,
           const int NEWTON_PAIR, const int CTABLE, const int LJTABLE, const int ORDER1, const int ORDER6 >
void PairLJLongTIP4PLongOMP::eval_outer(int iifrom, int iito, ThrData * const thr)
{
  double evdwl,ecoul,fvirial,fpair;
  evdwl = ecoul = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const double cut_coulsqplus = (cut_coul+2.0*qdist)*(cut_coul+2.0*qdist);

  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int n,vlist[6];
  int key;
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double fraction,table;
  double r,r2inv,forcecoul,forcelj,cforce, respa_coul, respa_lj, frespa;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double xiM[3],xjM[3],fO[3],fH[3],fd[3],v[6],xH1[3],xH2[3];// f1[3];
  dbl3_t x1,x2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq,qri;
  int respa_flag;
  
  int ni;
  double  *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;

  const double cut_in_off = cut_respa[2];
  const double cut_in_on = cut_respa[3];
  
  const double cut_in_diff = cut_in_on - cut_in_off;
  const double cut_in_off_sq = cut_in_off*cut_in_off;
  const double cut_in_on_sq = cut_in_on*cut_in_on;

  ilist = listouter->ilist;
  numneigh = listouter->numneigh;
  firstneigh = listouter->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    qri = qtmp*qqrd2e;
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    if (itype == typeO) {
      if (hneigh_thr[i].a < 0) {
        hneigh_thr[i].a = iH1 = atom->map(atom->tag[i] + 1);
        hneigh_thr[i].b = iH2 = atom->map(atom->tag[i] + 2);
        hneigh_thr[i].t = 1;
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
      } else {
        iH1 = hneigh_thr[i].a;
        iH2 = hneigh_thr[i].b;
        if (hneigh_thr[i].t == 0) {
          hneigh_thr[i].t = 1;
          compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
        }
      }
      x1 = newsite_thr[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    offseti = offset[itype];
    lj1i = lj1[itype]; lj2i = lj2[itype]; lj3i = lj3[itype]; lj4i = lj4[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      ni = sbmask(j);
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
        
      respa_coul = 0;
      respa_lj = 0;
      if (rsq < cut_ljsq[itype][jtype]) {			// lj
        frespa = 1.0;                                       // check whether and how to compute respa corrections
        respa_flag = rsq < cut_in_on_sq ? 1 : 0;
        if (respa_flag && (rsq > cut_in_off_sq)) {
          register double rsw = (sqrt(rsq)-cut_in_off)/cut_in_diff;
          frespa = 1-rsw*rsw*(3.0-2.0*rsw);
        }

        r2inv = 1.0/rsq;
        register double rn = r2inv*r2inv*r2inv;
        if (respa_flag) respa_lj = ni == 0 ?                 // correct for respa
            frespa*rn*(rn*lj1i[jtype]-lj2i[jtype]) :
            frespa*rn*(rn*lj1i[jtype]-lj2i[jtype])*special_lj[ni];
        if (ORDER6) {                                        // long-range form
          if (!ndisptablebits || rsq <= tabinnerdispsq) {
            register double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*lj4i[jtype];
            if (ni == 0) {
              forcelj =
                (rn*=rn)*lj1i[jtype]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq-respa_lj;
              if (EFLAG) evdwl = rn*lj3i[jtype]-g6*((a2+1.0)*a2+0.5)*x2;
            }
            else {                                        // correct for special
              register double f = special_lj[ni], t = rn*(1.0-f);
              forcelj = f*(rn *= rn)*lj1i[jtype]-
                g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[jtype]-respa_lj;
              if (EFLAG)
                evdwl = f*rn*lj3i[jtype]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[jtype];
            }
          }
          else {						// table real space
            register union_int_float_t disp_t;
            disp_t.f = rsq;
            register const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            register double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            if (ni == 0) {
              forcelj = (rn*=rn)*lj1i[jtype]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[jtype]-respa_lj;
              if (EFLAG) evdwl = rn*lj3i[jtype]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[jtype];
            }
            else {					// special case
              register double f = special_lj[ni], t = rn*(1.0-f);
              forcelj = f*(rn *= rn)*lj1i[jtype]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[jtype]+t*lj2i[jtype]-respa_lj;
              if (EFLAG) evdwl = f*rn*lj3i[jtype]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[jtype]+t*lj4i[jtype];
            }
          }
        }
        else {                                                // cut form
          if (ni == 0) {
            forcelj = rn*(rn*lj1i[jtype]-lj2i[jtype])-respa_lj;
            if (EFLAG) evdwl = rn*(rn*lj3i[jtype]-lj4i[jtype])-offseti[jtype];
          }
          else {                                        // correct for special
            register double f = special_lj[ni];
            forcelj = f*rn*(rn*lj1i[jtype]-lj2i[jtype])-respa_lj;
            if (EFLAG)
              evdwl = f*(rn*(rn*lj3i[jtype]-lj4i[jtype])-offseti[jtype]);
          }
        }

        forcelj *= r2inv;
	f[i][0] += delx*forcelj;
	f[i][1] += dely*forcelj;
	f[i][2] += delz*forcelj;
	f[j][0] -= delx*forcelj;
	f[j][1] -= dely*forcelj;
	f[j][2] -= delz*forcelj;
      
        if (EVFLAG) {
          fvirial = forcelj + respa_lj*r2inv;
          ev_tally_thr(this,i,j,nlocal,/*newton_pair = */ 1,
		       evdwl,0.0,fvirial,delx,dely,delz, thr);
        }
      }

      
      // adjust rsq and delxyz for off-site O charge(s)
      // ADDITIONAL REQEUST REQUIRED HERE!!!!!

      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) { 
	  if (jtype == typeO) {
            if (hneigh_thr[j].a < 0) {
              hneigh_thr[j].a = jH1 = atom->map(atom->tag[j] + 1);
              hneigh_thr[j].b = jH2 = atom->map(atom->tag[j] + 2);
              hneigh_thr[j].t = 1;
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
            } else {
              jH1 = hneigh_thr[j].a;
              jH2 = hneigh_thr[j].b;
              if (hneigh_thr[j].t == 0) {
                hneigh_thr[j].t = 1;
                compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
              }
            }
            x2 = newsite_thr[j];
	  } else x2 = x[j];
	  delx = x1.x - x2.x;
	  dely = x1.y - x2.y;
	  delz = x1.z - x2.z;
	  rsq = delx*delx + dely*dely + delz*delz;
        }

	// test current rsq against cutoff and compute Coulombic force
        if ((rsq < cut_coulsq) && ORDER1) {

          frespa = 1.0;                                       // check whether and how to compute respa corrections
          respa_flag = rsq < cut_in_on_sq ? 1 : 0;
          if (respa_flag && (rsq > cut_in_off_sq)) {
            register double rsw = (sqrt(rsq)-cut_in_off)/cut_in_diff;
            frespa = 1-rsw*rsw*(3.0-2.0*rsw);
          }

          r2inv = 1.0 / rsq;
          if (!CTABLE || rsq <= tabinnersq) {        // series real space
            register double r = sqrt(rsq), s = qri*q[j];
            if (respa_flag)                                // correct for respa
              respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
            register double x = g_ewald*r, t = 1.0/(1.0+EWALD_P*x);
            if (ni == 0) {
              s *= g_ewald*exp(-x*x);
              forcecoul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-respa_coul;
              if (EFLAG) ecoul = t;
            }
            else {                                        // correct for special
              r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
              forcecoul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r-respa_coul;
              if (EFLAG) ecoul = t-r;
            }
          }                                                // table real space
          else {
            if (respa_flag) {
              register double r = sqrt(rsq), s = qri*q[j];
              respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
            }
            register union_int_float_t t;
            t.f = rsq;
            register const int k = (t.i & ncoulmask) >> ncoulshiftbits;
            register double f = (t.f-rtable[k])*drtable[k], qiqj = qtmp*q[j];
            if (ni == 0) {
              forcecoul = qiqj*(ftable[k]+f*dftable[k]);
              if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]);
            }
            else {                                        // correct for special
              t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
              forcecoul = qiqj*(ftable[k]+f*dftable[k]-t.f);
              if (EFLAG) {
                t.f = (1.0-special_coul[ni])*(ptable[k]+f*dptable[k]);
                ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
              }
            }
          } 

          cforce = forcecoul * r2inv;
          fvirial = (forcecoul + respa_coul) * r2inv;

          // if i,j are not O atoms, force is applied directly
	  // if i or j are O atoms, force is on fictitious atom & partitioned
	  // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
	  // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
	  // preserves total force and torque on water molecule
	  // virial = sum(r x F) where each water's atoms are near xi and xj
	  // vlist stores 2,4,6 atoms whose forces contribute to virial

	  n = 0;
          key = 0;

	  if (itype != typeO) {
	    f[i][0] += delx * cforce;
            f[i][1] += dely * cforce;
	    f[i][2] += delz * cforce;

            if (EVFLAG) {
              v[0] = x[i].x * delx * fvirial;
              v[1] = x[i].y * dely * fvirial;
              v[2] = x[i].z * delz * fvirial;
              v[3] = x[i].x * dely * fvirial;
              v[4] = x[i].x * delz * fvirial;
              v[5] = x[i].y * delz * fvirial;
            }
          vlist[n++] = i;

	  } else {
            key += 1;
            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];

	    if (EVFLAG) {

              fd[0] = delx*fvirial;
              fd[1] = dely*fvirial;
              fd[2] = delz*fvirial;

              fO[0] = fd[0]*(1 - alpha);
              fO[1] = fd[1]*(1 - alpha);
              fO[2] = fd[2]*(1 - alpha);

              fH[0] = 0.5 * alpha * fd[0];
              fH[1] = 0.5 * alpha * fd[1];
              fH[2] = 0.5 * alpha * fd[2];

	      domain->closest_image(&x[i].x,&x[iH1].x,xH1);
	      domain->closest_image(&x[i].x,&x[iH2].x,xH2);

	      v[0] = x[i].x*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] = x[i].y*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] = x[i].z*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] = x[i].x*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] = x[i].x*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] = x[i].y*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
	    }
	    vlist[n++] = i;
	    vlist[n++] = iH1;
	    vlist[n++] = iH2;
  	  }

	  if (jtype != typeO) {
	    f[j][0] -= delx * cforce;
	    f[j][1] -= dely * cforce;
	    f[j][2] -= delz * cforce;

	    if (EVFLAG) {
	      v[0] -= x[j].x * delx * fvirial;
	      v[1] -= x[j].y * dely * fvirial;
	      v[2] -= x[j].z * delz * fvirial;
	      v[3] -= x[j].x * dely * fvirial;
	      v[4] -= x[j].x * delz * fvirial;
	      v[5] -= x[j].y * delz * fvirial;
            }
	    vlist[n++] = j;

	  } else {
            key += 2;

	    fd[0] = -delx*cforce;
	    fd[1] = -dely*cforce;
	    fd[2] = -delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2]; 

	    f[j][0] += fO[0];
	    f[j][1] += fO[1];
	    f[j][2] += fO[2];

	    f[jH1][0] += fH[0];
	    f[jH1][1] += fH[1];
	    f[jH1][2] += fH[2];

	    f[jH2][0] += fH[0];
	    f[jH2][1] += fH[1];
	    f[jH2][2] += fH[2];

	    if (EVFLAG) {

	      fd[0] = -delx*fvirial;
	      fd[1] = -dely*fvirial;
	      fd[2] = -delz*fvirial;

              fO[0] = fd[0]*(1 - alpha);
              fO[1] = fd[1]*(1 - alpha);
              fO[2] = fd[2]*(1 - alpha);

              fH[0] = 0.5 * alpha * fd[0];
              fH[1] = 0.5 * alpha * fd[1];
              fH[2] = 0.5 * alpha * fd[2]; 

	      domain->closest_image(&x[j].x,&x[jH1].x,xH1);
	      domain->closest_image(&x[j].x,&x[jH2].x,xH2);

	      v[0] += x[j].x*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] += x[j].y*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] += x[j].z*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] += x[j].x*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] += x[j].x*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] += x[j].y*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
            }
      	    vlist[n++] = j;
	    vlist[n++] = jH1;
	    vlist[n++] = jH2;
	  }
 
          if (EVFLAG) ev_tally_list_thr(this,key,vlist,v,ecoul,alpha,thr);
        }
      }
    }
  }
}


/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairLJLongTIP4PLongOMP::compute_newsite_thr(const dbl3_t &xO,
                                                const dbl3_t &xH1,
                                                const dbl3_t &xH2,
                                                dbl3_t &xM) const
{
  double delx1 = xH1.x - xO.x;
  double dely1 = xH1.y - xO.y;
  double delz1 = xH1.z - xO.z;
  domain->minimum_image(delx1,dely1,delz1);

  double delx2 = xH2.x - xO.x;
  double dely2 = xH2.y - xO.y;
  double delz2 = xH2.z - xO.z;
  domain->minimum_image(delx2,dely2,delz2);

  const double prefac = alpha * 0.5;
  xM.x = xO.x + prefac * (delx1 + delx2);
  xM.y = xO.y + prefac * (dely1 + dely2);
  xM.z = xO.z + prefac * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

double PairLJLongTIP4PLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJLongTIP4PLong::memory_usage();
  return bytes;
}

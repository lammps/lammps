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

/* ----------------------------------------------------------------------
   OPT version: Wayne Mitchell (Loyola University New Orleans)
------------------------------------------------------------------------- */

#include "pair_lj_long_coul_long_opt.h"

#include "atom.h"
#include "force.h"
#include "math_extra.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathExtra;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJLongCoulLongOpt::PairLJLongCoulLongOpt(LAMMPS *lmp) : PairLJLongCoulLong(lmp)
{
  respa_enable = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLongOpt::compute(int eflag, int vflag)
{

  ev_init(eflag,vflag);
  int order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);

  if (order6) {
    if (order1) {
      if (!ndisptablebits) {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,0,0,1,1>();
              else return eval<1,1,0,0,0,1,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,0,1,1>();
              else return eval<1,0,0,0,0,1,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,0,1,1>();
            else return eval<0,0,0,0,0,1,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,0,1,1>();
              else return eval<1,1,0,1,0,1,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,0,1,1>();
              else return eval<1,0,0,1,0,1,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,0,1,1>();
            else return eval<0,0,0,1,0,1,1>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,0,1,1,1>();
              else return eval<1,1,0,0,1,1,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,1,1,1>();
              else return eval<1,0,0,0,1,1,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,1,1,1>();
            else return eval<0,0,0,0,1,1,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,1,1,1>();
              else return eval<1,1,0,1,1,1,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,1,1,1>();
              else return eval<1,0,0,1,1,1,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,1,1,1>();
            else return eval<0,0,0,1,1,1,1>();
          }
        }
      }
    } else {
      if (!ndisptablebits) {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,0,0,0,1>();
              else return eval<1,1,0,0,0,0,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,0,0,1>();
              else return eval<1,0,0,0,0,0,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,0,0,1>();
            else return eval<0,0,0,0,0,0,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,0,0,1>();
              else return eval<1,1,0,1,0,0,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,0,0,1>();
              else return eval<1,0,0,1,0,0,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,0,0,1>();
            else return eval<0,0,0,1,0,0,1>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,0,1,0,1>();
              else return eval<1,1,0,0,1,0,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,1,0,1>();
              else return eval<1,0,0,0,1,0,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,1,0,1>();
            else return eval<0,0,0,0,1,0,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,1,0,1>();
              else return eval<1,1,0,1,1,0,1>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,1,0,1>();
              else return eval<1,0,0,1,1,0,1>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,1,0,1>();
            else return eval<0,0,0,1,1,0,1>();
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
              if (force->newton_pair) return eval<1,1,1,0,0,1,0>();
              else return eval<1,1,0,0,0,1,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,0,1,0>();
              else return eval<1,0,0,0,0,1,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,0,1,0>();
            else return eval<0,0,0,0,0,1,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,0,1,0>();
              else return eval<1,1,0,1,0,1,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,0,1,0>();
              else return eval<1,0,0,1,0,1,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,0,1,0>();
            else return eval<0,0,0,1,0,1,0>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,0,1,1,0>();
              else return eval<1,1,0,0,1,1,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,1,1,0>();
              else return eval<1,0,0,0,1,1,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,1,1,0>();
            else return eval<0,0,0,0,1,1,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,1,1,0>();
              else return eval<1,1,0,1,1,1,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,1,1,0>();
              else return eval<1,0,0,1,1,1,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,1,1,0>();
            else return eval<0,0,0,1,1,1,0>();
          }
        }
      }
    } else {
      if (!ndisptablebits) {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,0,0,0,0>();
              else return eval<1,1,0,0,0,0,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,0,0,0>();
              else return eval<1,0,0,0,0,0,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,0,0,0>();
            else return eval<0,0,0,0,0,0,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,0,0,0>();
              else return eval<1,1,0,1,0,0,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,0,0,0>();
              else return eval<1,0,0,1,0,0,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,0,0,0>();
            else return eval<0,0,0,1,0,0,0>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,0,1,0,0>();
              else return eval<1,1,0,0,1,0,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,0,1,0,0>();
              else return eval<1,0,0,0,1,0,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,0,1,0,0>();
            else return eval<0,0,0,0,1,0,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval<1,1,1,1,1,0,0>();
              else return eval<1,1,0,1,1,0,0>();
            } else {
              if (force->newton_pair) return eval<1,0,1,1,1,0,0>();
              else return eval<1,0,0,1,1,0,0>();
            }
          } else {
            if (force->newton_pair) return eval<0,0,1,1,1,0,0>();
            else return eval<0,0,0,1,1,0,0>();
          }
        }
      }
    }
  }
}

void PairLJLongCoulLongOpt::compute_outer(int eflag, int vflag)
{

  ev_init(eflag,vflag);
  int order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);

  if (order6) {
    if (order1) {
      if (!ndisptablebits) {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,0,0,1,1>();
              else return eval_outer<1,1,0,0,0,1,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,0,1,1>();
              else return eval_outer<1,0,0,0,0,1,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,0,1,1>();
            else return eval_outer<0,0,0,0,0,1,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,0,1,1>();
              else return eval_outer<1,1,0,1,0,1,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,0,1,1>();
              else return eval_outer<1,0,0,1,0,1,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,0,1,1>();
            else return eval_outer<0,0,0,1,0,1,1>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,0,1,1,1>();
              else return eval_outer<1,1,0,0,1,1,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,1,1,1>();
              else return eval_outer<1,0,0,0,1,1,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,1,1,1>();
            else return eval_outer<0,0,0,0,1,1,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,1,1,1>();
              else return eval_outer<1,1,0,1,1,1,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,1,1,1>();
              else return eval_outer<1,0,0,1,1,1,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,1,1,1>();
            else return eval_outer<0,0,0,1,1,1,1>();
          }
        }
      }
    } else {
      if (!ndisptablebits) {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,0,0,0,1>();
              else return eval_outer<1,1,0,0,0,0,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,0,0,1>();
              else return eval_outer<1,0,0,0,0,0,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,0,0,1>();
            else return eval_outer<0,0,0,0,0,0,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,0,0,1>();
              else return eval_outer<1,1,0,1,0,0,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,0,0,1>();
              else return eval_outer<1,0,0,1,0,0,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,0,0,1>();
            else return eval_outer<0,0,0,1,0,0,1>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,0,1,0,1>();
              else return eval_outer<1,1,0,0,1,0,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,1,0,1>();
              else return eval_outer<1,0,0,0,1,0,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,1,0,1>();
            else return eval_outer<0,0,0,0,1,0,1>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,1,0,1>();
              else return eval_outer<1,1,0,1,1,0,1>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,1,0,1>();
              else return eval_outer<1,0,0,1,1,0,1>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,1,0,1>();
            else return eval_outer<0,0,0,1,1,0,1>();
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
              if (force->newton_pair) return eval_outer<1,1,1,0,0,1,0>();
              else return eval_outer<1,1,0,0,0,1,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,0,1,0>();
              else return eval_outer<1,0,0,0,0,1,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,0,1,0>();
            else return eval_outer<0,0,0,0,0,1,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,0,1,0>();
              else return eval_outer<1,1,0,1,0,1,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,0,1,0>();
              else return eval_outer<1,0,0,1,0,1,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,0,1,0>();
            else return eval_outer<0,0,0,1,0,1,0>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,0,1,1,0>();
              else return eval_outer<1,1,0,0,1,1,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,1,1,0>();
              else return eval_outer<1,0,0,0,1,1,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,1,1,0>();
            else return eval_outer<0,0,0,0,1,1,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,1,1,0>();
              else return eval_outer<1,1,0,1,1,1,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,1,1,0>();
              else return eval_outer<1,0,0,1,1,1,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,1,1,0>();
            else return eval_outer<0,0,0,1,1,1,0>();
          }
        }
      }
    } else {
      if (!ndisptablebits) {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,0,0,0,0>();
              else return eval_outer<1,1,0,0,0,0,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,0,0,0>();
              else return eval_outer<1,0,0,0,0,0,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,0,0,0>();
            else return eval_outer<0,0,0,0,0,0,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,0,0,0>();
              else return eval_outer<1,1,0,1,0,0,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,0,0,0>();
              else return eval_outer<1,0,0,1,0,0,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,0,0,0>();
            else return eval_outer<0,0,0,1,0,0,0>();
          }
        }
      } else {
        if (!ncoultablebits) {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,0,1,0,0>();
              else return eval_outer<1,1,0,0,1,0,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,0,1,0,0>();
              else return eval_outer<1,0,0,0,1,0,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,0,1,0,0>();
            else return eval_outer<0,0,0,0,1,0,0>();
          }
        } else {
          if (evflag) {
            if (eflag) {
              if (force->newton_pair) return eval_outer<1,1,1,1,1,0,0>();
              else return eval_outer<1,1,0,1,1,0,0>();
            } else {
              if (force->newton_pair) return eval_outer<1,0,1,1,1,0,0>();
              else return eval_outer<1,0,0,1,1,0,0>();
            }
          } else {
            if (force->newton_pair) return eval_outer<0,0,1,1,1,0,0>();
            else return eval_outer<0,0,0,1,1,0,0>();
          }
        }
      }
    }
  }
}


template < const int EVFLAG, const int EFLAG,
           const int NEWTON_PAIR, const int CTABLE, const int LJTABLE, const int ORDER1, const int ORDER6 >
void PairLJLongCoulLongOpt::eval()
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;

  double **x = atom->x, *x0 = x[0];
  double **f = atom->f, *f0 = f[0], *fi = f0;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  int i, j;
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  double qi = 0.0, qri = 0.0;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double xi[3], d[3];

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {                        // loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (ORDER1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    offseti = offset[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    memcpy(xi, x0+(i+(i<<1)), 3*sizeof(double));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {                        // loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;

      { double *xj = x0+(j+(j<<1));
        d[0] = xi[0] - xj[0];                                // pair vector
        d[1] = xi[1] - xj[1];
        d[2] = xi[2] - xj[2]; }

      if ((rsq = dot3(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;

      if (ORDER1 && (rsq < cut_coulsq)) {                // coulombic
        if (!CTABLE || rsq <= tabinnersq) {        // series real space
          double r = sqrt(rsq), x = g_ewald*r;
          double s = qri*q[j], t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
            if (EFLAG) ecoul = t;
          }
          else {                                        // special case
            r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r;
            if (EFLAG) ecoul = t-r;
          }
        }                                                // table real space
        else {
          union_int_float_t t;
          t.f = rsq;
          const int k = (t.i & ncoulmask)>>ncoulshiftbits;
          double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+f*dftable[k]);
            if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]);
          }
          else {                                        // special case
            t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
            force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
            if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
          }
        }
      }
      else force_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {                        // lj
        if (ORDER6) {                                        // long-range lj
          if (!LJTABLE || rsq <= tabinnerdispsq) {               // series real space
            double rn = r2inv*r2inv*r2inv;
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*lj4i[typej];
            if (ni == 0) {
              force_lj =
              (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
              if (EFLAG)
                evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            }
            else {                                        // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-
              g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej];
              if (EFLAG)
                evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
            }
          }
          else {                        // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              force_lj = (rn*=rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej];
              if (EFLAG) evdwl = rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej];
            }
            else {                  // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]+t*lj2i[typej];
              if (EFLAG) evdwl = f*rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej]+t*lj4i[typej];
            }
          }
        }
        else {                                                // cut lj
          double rn = r2inv*r2inv*r2inv;
          if (ni == 0) {
            force_lj = rn*(rn*lj1i[typej]-lj2i[typej]);
            if (EFLAG) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
          }
          else {                                        // special case
            double f = special_lj[ni];
            force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej]);
            if (EFLAG)
              evdwl = f * (rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      }
      else force_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

      if (NEWTON_PAIR || j < nlocal) {
        double *fj = f0+(j+(j<<1)), f;
        fi[0] += f = d[0]*fpair; fj[0] -= f;
        fi[1] += f = d[1]*fpair; fj[1] -= f;
        fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
        fi[0] += d[0]*fpair;
        fi[1] += d[1]*fpair;
        fi[2] += d[2]*fpair;
      }

      if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                           evdwl,ecoul,fpair,d[0],d[1],d[2]);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ---------------------------------------------------------------------- */

template < const int EVFLAG, const int EFLAG,
           const int NEWTON_PAIR, const int CTABLE, const int LJTABLE, const int ORDER1, const int ORDER6 >
void PairLJLongCoulLongOpt::eval_outer()
{
  double evdwl,ecoul,fvirial,fpair;
  evdwl = ecoul = 0.0;

  double **x = atom->x, *x0 = x[0];
  double **f = atom->f, *f0 = f[0], *fi = f0;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j;
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni, respa_flag;
  double qi = 0.0, qri = 0.0;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double respa_lj = 0.0, respa_coul = 0.0, frespa = 0.0;
  double xi[3], d[3];

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {                        // loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (ORDER1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    offseti = offset[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    memcpy(xi, x0+(i+(i<<1)), 3*sizeof(double));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {                        // loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;

      { double *xj = x0+(j+(j<<1));
        d[0] = xi[0] - xj[0];                                // pair vector
        d[1] = xi[1] - xj[1];
        d[2] = xi[2] - xj[2]; }

      if ((rsq = dot3(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;

      frespa = 1.0;                                       // check whether and how to compute respa corrections
      respa_coul = 0;
      respa_lj = 0;
      respa_flag = rsq < cut_in_on_sq ? 1 : 0;
      if (respa_flag && (rsq > cut_in_off_sq)) {
        double rsw = (sqrt(rsq)-cut_in_off)/cut_in_diff;
        frespa = 1-rsw*rsw*(3.0-2.0*rsw);
      }

      if (ORDER1 && (rsq < cut_coulsq)) {                // coulombic
        if (!CTABLE || rsq <= tabinnersq) {        // series real space
          double r = sqrt(rsq), s = qri*q[j];
          if (respa_flag)                                // correct for respa
            respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
          double x = g_ewald*r, t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-respa_coul;
            if (EFLAG) ecoul = t;
          }
          else {                                        // correct for special
            r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r-respa_coul;
            if (EFLAG) ecoul = t-r;
          }
        }                                                // table real space
        else {
          if (respa_flag) {
            double r = sqrt(rsq), s = qri*q[j];
            respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
          }
          union_int_float_t t;
          t.f = rsq;
          const int k = (t.i & ncoulmask) >> ncoulshiftbits;
          double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+f*dftable[k]);
            if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]);
          }
          else {                                        // correct for special
            t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
            force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
            if (EFLAG) {
              t.f = (1.0-special_coul[ni])*(ptable[k]+f*dptable[k]);
              ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
            }
          }
        }
      }

      else force_coul = respa_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {                        // lennard-jones
        double rn = r2inv*r2inv*r2inv;
        if (respa_flag) respa_lj = ni == 0 ?                 // correct for respa
            frespa*rn*(rn*lj1i[typej]-lj2i[typej]) :
            frespa*rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
        if (ORDER6) {                                        // long-range form
          if (!LJTABLE || rsq <= tabinnerdispsq) {
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*lj4i[typej];
            if (ni == 0) {
              force_lj =
                (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq-respa_lj;
              if (EFLAG) evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            }
            else {                                        // correct for special
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-
                g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej]-respa_lj;
              if (EFLAG)
                evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
            }
          }
          else {                        // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              force_lj = (rn*=rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]-respa_lj;
              if (EFLAG) evdwl = rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej];
            }
            else {                  // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]+t*lj2i[typej]-respa_lj;
              if (EFLAG) evdwl = f*rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej]+t*lj4i[typej];
            }
          }
        }
        else {                                                // cut form
          if (ni == 0) {
            force_lj = rn*(rn*lj1i[typej]-lj2i[typej])-respa_lj;
            if (EFLAG) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
          }
          else {                                        // correct for special
            double f = special_lj[ni];
            force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej])-respa_lj;
            if (EFLAG)
              evdwl = f*(rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      }
      else force_lj = respa_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

      if (NEWTON_PAIR || j < nlocal) {
        double *fj = f0+(j+(j<<1)), f;
        fi[0] += f = d[0]*fpair; fj[0] -= f;
        fi[1] += f = d[1]*fpair; fj[1] -= f;
        fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
        fi[0] += d[0]*fpair;
        fi[1] += d[1]*fpair;
        fi[2] += d[2]*fpair;
      }

      if (EVFLAG) {
        fvirial = (force_coul + force_lj + respa_coul + respa_lj)*r2inv;
        ev_tally(i,j,nlocal,newton_pair,
                 evdwl,ecoul,fvirial,d[0],d[1],d[2]);
      }
    }
  }
}

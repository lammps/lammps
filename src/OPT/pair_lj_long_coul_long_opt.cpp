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

/* ----------------------------------------------------------------------
   OPT version: Wayne Mitchell (Loyola University New Orleans)
------------------------------------------------------------------------- */

#include "pair_lj_long_coul_long_opt.h"

#include "atom.h"
#include "ewald_const.h"
#include "force.h"
#include "math_extra.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairLJLongCoulLongOpt::PairLJLongCoulLongOpt(LAMMPS *lmp) : PairLJLongCoulLong(lmp)
{
  respa_enable = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLongOpt::compute(int eflag, int vflag)
{

  ev_init(eflag,vflag);
  int order1 = ewald_order & EWALD_COUL, order6 = ewald_order & EWALD_DISP;

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
  int order1 = ewald_order & EWALD_COUL, order6 = ewald_order & EWALD_DISP;

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

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  int i, j, ii, jj, typei, typej, ni;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double qi = 0.0, qri = 0.0;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double xi[3], d[3];

  for (ii = 0; ii < inum; ii++) {                          // loop over my atoms
    i = ilist[ii];
    if (ORDER1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    typei = type[i];
    offseti = offset[typei];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                          // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
      r2inv = 1.0/rsq;

      if (ORDER1 && (rsq < cut_coulsq)) {                // coulombic
        if (!CTABLE || rsq <= tabinnersq) {        // series real space
          double r = sqrt(rsq);
          double grij = g_ewald*r;
          double t = 1.0/(1.0+EWALD_P*grij);
          double erfc_poly = ((((t*A5+A4)*t+A3)*t+A2)*t+A1);
          double pre = qri*q[j];                            // qqrd2e * qi * qj
          if (ni == 0) {
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre;
            if (EFLAG) ecoul = t;
          } else {                                        // special case
            double adjust = pre*(1.0-special_coul[ni])/r;
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre-adjust;
            if (EFLAG) ecoul = t-adjust;
          }
        } else { // table real space
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          const int k = (rsq_lookup.i & ncoulmask)>>ncoulshiftbits;
          double fraction = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]);
            if (EFLAG) ecoul = qiqj*(etable[k]+fraction*detable[k]);
          } else {                                        // special case
            rsq_lookup.f = (1.0-special_coul[ni])*(ctable[k]+fraction*dctable[k]);
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]-(double)rsq_lookup.f);
            if (EFLAG) ecoul = qiqj*(etable[k]+fraction*detable[k]-(double)rsq_lookup.f);
          }
        }
      }
      else force_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {                        // lj
        if (ORDER6) {                                        // long-range lj
          if (!LJTABLE || rsq <= tabinnerdispsq) {               // series real space
            double r6inv = r2inv*r2inv*r2inv;
            double r12inv = r6inv*r6inv;
            double gr2 = g2*rsq, a2 = 1.0/gr2;
            double expterm = a2*exp(-gr2)*lj4i[typej];
            double g6term = g6*((a2+1.0)*a2+0.5)*expterm;
            double g8term = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*expterm*rsq;
            if (ni == 0) {
              force_lj = r12inv*lj1i[typej]-g8term;
              if (EFLAG) evdwl = r12inv*lj3i[typej]-g6term;
            } else {                                        // special case
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_lj = factor*r12inv*lj1i[typej]-g8term+t*lj2i[typej];
              if (EFLAG) evdwl = factor*r12inv*lj3i[typej]-g6term+t*lj4i[typej];
            }
          } else {                        // table real space
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            const int disp_k = (rsq_lookup.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double ftable_disp = fdisptable[disp_k]+f_disp*dfdisptable[disp_k];
            double etable_disp = edisptable[disp_k]+f_disp*dedisptable[disp_k];
            double r6inv = r2inv*r2inv*r2inv;
            double r12inv = r6inv*r6inv;
            if (ni == 0) {
              force_lj = r12inv*lj1i[typej]-ftable_disp*lj4i[typej];
              if (EFLAG) evdwl = r12inv*lj3i[typej]-etable_disp*lj4i[typej];
            } else {                  // special case
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_lj = factor*r12inv*lj1i[typej]-ftable_disp*lj4i[typej]+t*lj2i[typej];
              if (EFLAG) evdwl = factor*r12inv*lj3i[typej]-etable_disp*lj4i[typej]+t*lj4i[typej];
            }
          }
        } else {                                                // cut lj
          double r6inv = r2inv*r2inv*r2inv;
          if (ni == 0) {
            force_lj = r6inv*(r6inv*lj1i[typej]-lj2i[typej]);
            if (EFLAG) evdwl = r6inv*(r6inv*lj3i[typej]-lj4i[typej])-offseti[typej];
          } else {                                        // special case
            double factor = special_lj[ni];
            force_lj = factor*r6inv*(r6inv*lj1i[typej]-lj2i[typej]);
            if (EFLAG)
              evdwl = factor * (r6inv*(r6inv*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      }
      else force_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      f[i][0] += fdx;
      f[i][1] += fdy;
      f[i][2] += fdz;
      if (NEWTON_PAIR || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
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

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, ii, jj, typei, typej, ni, respa_flag;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
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

  for (ii = 0; ii < inum; ii++) {                          // loop over my atoms
    i = ilist[ii];
    if (ORDER1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    typei = type[i];
    offseti = offset[typei];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                          // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
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
          double r = sqrt(rsq), pre = qri*q[j];
          if (respa_flag)                                // correct for respa
            respa_coul = ni == 0 ? frespa*pre/r : frespa*pre/r*special_coul[ni];
          double grij = g_ewald*r, t = 1.0/(1.0+EWALD_P*grij);
          double erfc_poly = ((((t*A5+A4)*t+A3)*t+A2)*t+A1);
          if (ni == 0) {
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre-respa_coul;
            if (EFLAG) ecoul = t;
          } else {                                        // correct for special
            double adjust = pre*(1.0-special_coul[ni])/r;
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre-adjust-respa_coul;
            if (EFLAG) ecoul = t-adjust;
          }
        } else {                                            // table real space
          if (respa_flag) {
            double r = sqrt(rsq), pre = qri*q[j];
            respa_coul = ni == 0 ? frespa*pre/r : frespa*pre/r*special_coul[ni];
          }
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          const int k = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
          double fraction = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]);
            if (EFLAG) ecoul = qiqj*(etable[k]+fraction*detable[k]);
          }
          else {                                        // correct for special
            rsq_lookup.f = (1.0-special_coul[ni])*(ctable[k]+fraction*dctable[k]);
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]-(double)rsq_lookup.f);
            if (EFLAG) {
              rsq_lookup.f = (1.0-special_coul[ni])*(ptable[k]+fraction*dptable[k]);
              ecoul = qiqj*(etable[k]+fraction*detable[k]-(double)rsq_lookup.f);
            }
          }
        }
      }

      else force_coul = respa_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {                        // lennard-jones
        double r6inv = r2inv*r2inv*r2inv;
        if (respa_flag) respa_lj = ni == 0 ?                 // correct for respa
            frespa*r6inv*(r6inv*lj1i[typej]-lj2i[typej]) :
            frespa*r6inv*(r6inv*lj1i[typej]-lj2i[typej])*special_lj[ni];
        if (ORDER6) {                                        // long-range form
          if (!LJTABLE || rsq <= tabinnerdispsq) {
            double r12inv = r6inv*r6inv;
            double gr2 = g2*rsq, a2 = 1.0/gr2;
            double expterm = a2*exp(-gr2)*lj4i[typej];
            double g6term = g6*((a2+1.0)*a2+0.5)*expterm;
            double g8term = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*expterm*rsq;
            if (ni == 0) {
              force_lj = r12inv*lj1i[typej]-g8term-respa_lj;
              if (EFLAG) evdwl = r12inv*lj3i[typej]-g6term;
            } else {                                        // correct for special
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_lj = factor*r12inv*lj1i[typej]-g8term+t*lj2i[typej]-respa_lj;
              if (EFLAG) evdwl = factor*r12inv*lj3i[typej]-g6term+t*lj4i[typej];
            }
          } else {                        // table real space
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            const int disp_k = (rsq_lookup.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double ftable_disp = fdisptable[disp_k]+f_disp*dfdisptable[disp_k];
            double etable_disp = edisptable[disp_k]+f_disp*dedisptable[disp_k];
            double r12inv = r6inv*r6inv;
            if (ni == 0) {
              force_lj = r12inv*lj1i[typej]-ftable_disp*lj4i[typej]-respa_lj;
              if (EFLAG) evdwl = r12inv*lj3i[typej]-etable_disp*lj4i[typej];
            } else {                  // special case
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_lj = factor*r12inv*lj1i[typej]-ftable_disp*lj4i[typej]+t*lj2i[typej]-respa_lj;
              if (EFLAG) evdwl = factor*r12inv*lj3i[typej]-etable_disp*lj4i[typej]+t*lj4i[typej];
            }
          }
        } else {                                                // cut form
          if (ni == 0) {
            force_lj = r6inv*(r6inv*lj1i[typej]-lj2i[typej])-respa_lj;
            if (EFLAG) evdwl = r6inv*(r6inv*lj3i[typej]-lj4i[typej])-offseti[typej];
          } else {                                        // correct for special
            double factor = special_lj[ni];
            force_lj = factor*r6inv*(r6inv*lj1i[typej]-lj2i[typej])-respa_lj;
            if (EFLAG)
              evdwl = factor*(r6inv*(r6inv*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      }
      else force_lj = respa_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      f[i][0] += fdx;
      f[i][1] += fdy;
      f[i][2] += fdz;
      if (NEWTON_PAIR || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }

      if (EVFLAG) {
        fvirial = (force_coul + force_lj + respa_coul + respa_lj)*r2inv;
        ev_tally(i,j,nlocal,newton_pair,
                 evdwl,ecoul,fvirial,d[0],d[1],d[2]);
      }
    }
  }
}

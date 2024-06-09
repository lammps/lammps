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

#include "pair_hybrid_molecular.h"

#include "atom.h"
#include "error.h"
#include "neigh_request.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybridMolecular::PairHybridMolecular(LAMMPS *lmp) : PairHybridOverlay(lmp) {}

/* ----------------------------------------------------------------------
  modify neighbor list requests
------------------------------------------------------------------------- */

void PairHybridMolecular::init_style()
{
  if (!atom->molecule_flag)
    error->all(FLERR, "Pair style hybrid/molecular requires atom attribute molecule");
  if (manybody_flag)
    error->all(FLERR, "Pair style hybrid/molecular is not compatible with manybody potentials");

  PairHybridOverlay::init_style();

  // modify neighbor list requests

  bool first = true;
  for (auto &request : neighbor->get_pair_requests()) {
    if (first) {
      request->set_molskip(NeighRequest::INTRA);
      first = false;
    } else {
      request->set_molskip(NeighRequest::INTER);
    }
  }
  born_matrix_enable = 0;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHybridMolecular::init_one(int i, int j)
{
  // if I,J is not set explicitly:
  // perform mixing only if I,I sub-style = J,J sub-style
  // plus I,I and J,J need the same number of substyles

  if (setflag[i][j] == 0) {
    if (nmap[i][i] != nmap[j][j])
      error->one(FLERR,"All pair coeffs are not set");
    int num = 0;
    for (int k = 0; k < nmap[i][i]; ++k) {
      for (int l = 0; l < nmap[j][j]; ++l) {
        if (map[i][i][k] == map[j][j][l]) {
          map[i][j][num] = map[i][i][k];
          ++num;
          nmap[i][j] = num;
        }
      }
    }
    if (nmap[i][i] != nmap[i][j])
      error->one(FLERR,"All pair coeffs are not set");
  }
  nmap[j][i] = nmap[i][j];

  // call init/mixing for all sub-styles of I,J
  // set cutsq in sub-style just as Pair::init() does via call to init_one()
  // set cutghost for I,J and J,I just as sub-style does
  // sum tail corrections for I,J
  // return max cutoff of all sub-styles assigned to I,J
  // if no sub-styles assigned to I,J (pair_coeff none), cutmax = 0.0 returned

  double cutmax = 0.0;
  cutghost[i][j] = cutghost[j][i] = 0.0;
  if (tail_flag) etail_ij = ptail_ij = 0.0;

  for (int k = 0; k < nmap[i][j]; k++) {
    map[j][i][k] = map[i][j][k];
    double cut = styles[map[i][j][k]]->init_one(i,j);
    if (styles[map[i][j][k]]->did_mix) did_mix = true;
    styles[map[i][j][k]]->cutsq[i][j] = styles[map[i][j][k]]->cutsq[j][i] = cut*cut;
    if (styles[map[i][j][k]]->ghostneigh)
      cutghost[i][j] = cutghost[j][i] = MAX(cutghost[i][j],styles[map[i][j][k]]->cutghost[i][j]);
    if (tail_flag) {
      etail_ij += styles[map[i][j][k]]->etail_ij;
      ptail_ij += styles[map[i][j][k]]->ptail_ij;
    }
    cutmax = MAX(cutmax,cut);

    int istyle;
    for (istyle = 0; istyle < nstyles; istyle++)
      if (styles[istyle] == styles[map[i][j][k]]) break;

    if (styles[istyle]->trim_flag) {

      if (cut > cutmax_style[istyle]) {
        cutmax_style[istyle] = cut;

        for (auto &request : neighbor->get_pair_requests()) {
          if (styles[istyle] == request->get_requestor()) {
            request->set_cutoff(cutmax_style[istyle]);
            break;
          }
        }
      }
    }
  }
  return cutmax;
}

/* ----------------------------------------------------------------------
   call sub-style to compute single interaction
   error if sub-style does not support single() call
   since overlay could have multiple sub-styles, sum results explicitly
------------------------------------------------------------------------- */

double PairHybridMolecular::single(int i, int j, int itype, int jtype, double rsq,
                                   double factor_coul, double factor_lj, double &fforce)
{
  if (nmap[itype][jtype] == 0)
    error->one(FLERR,"Invoked pair single() on sub-style none");

  double fone;
  fforce = 0.0;
  double esum = 0.0;
  if (nmap[itype][jtype] == 2) {
    int m = 0;
    if (atom->molecule[i] != atom->molecule[j]) m = 1;
    const int mystyle = map[itype][jtype][m];
    if (rsq < styles[mystyle]->cutsq[itype][jtype]) {
      if (styles[mystyle]->single_enable == 0)
        error->one(FLERR,"Pair hybrid/molecular sub-style {} does not support single() call",
                   keywords[mystyle]);

      if ((special_lj[mystyle] != nullptr) || (special_coul[mystyle] != nullptr))
        error->one(FLERR,"Pair hybrid/molecular single() calls do not support per sub-style "
                   "special bond values");

      esum += styles[mystyle]->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fone);
      fforce += fone;
    }
  }

  if (single_extra) copy_svector(itype,jtype);
  return esum;
}

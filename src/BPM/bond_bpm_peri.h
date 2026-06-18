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

#ifdef BOND_CLASS
// clang-format off
BondStyle(bpm/peri,BondBPMPeri);
// clang-format on
#else

#ifndef LMP_BOND_BPM_PERI_H
#define LMP_BOND_BPM_PERI_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMPeri : public BondBPM {
 public:
  BondBPMPeri(class LAMMPS *);
  ~BondBPMPeri() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void settings(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, double, int, int, double &) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  // peridynamic constitutive models (internal ids; lower-case keywords in input)
  enum { PMB, LPS, VES, EPS };
  // tags selecting which per-atom array the next forward_comm carries
  enum { COMM_SMIN, COMM_THETA, COMM_WVOLUME, COMM_TDNORM };

  int *model;             // per bond type: constitutive model id
  double *c;              // PMB micromodulus c = 18 K / (pi delta^4)
  double *kbulk, *gshear; // LPS/VES/EPS bulk and shear moduli (per bond type)
  double *lambda, *tau;   // VES viscoelastic modulus ratio and relaxation time
  double *yieldstress;    // EPS yield stress (per bond type)
  double *cut;            // horizon delta (per bond type)
  double *s00, *alpha;    // critical-stretch bond-break parameters (#984 rule)

  // per-atom property/atom storage: vfrac (user-supplied nodal volume), the
  // critical-stretch bookkeeping s0 (diagnostic) / smin (break state), and the
  // EPS accumulated plastic multiplier lambdaValue (public diagnostic)
  char *id_fix_property_peri;
  int index_vfrac, index_s0, index_smin, index_lambda;
  int index_vinter;       // reference interaction volume (for damage diagnostics)
  int index_dtheta;       // per-step dilatation theta exposed for visualization

  // per-step scratch for the break bookkeeping (committed to smin/s0 each step)
  double *smin_new, *s0_new;
  int nmax;

  // state-based (LPS/VES/EPS) machinery: weighted volume (static) and dilatation
  // theta (per step), both ghost-communicated; allocated only when in use
  int state_based;        // 1 if any bond type uses a state-based model
  int wvolume_setup;      // 1 once the static weighted volume has been computed
  double kbulk_rep;       // bulk modulus for the per-atom volumetric energy
  double *wvolume, *theta;
  double *winv;           // 1/wvolume, precomputed each step (avoids per-bond divides)
  int commflag;           // selects the array packed by the next forward_comm

  // EPS plasticity: per-atom deviatoric force-state norm (ghost-communicated),
  // plus a representative yield/shear for the per-atom return-mapping
  int plastic;            // 1 if any bond type uses EPS
  double *tdnorm;
  double pointwise_yield, gshear_rep;

  void allocate();
  void store_data() override;
  void compute_wvolume();
  void compute_dilatation();
  void compute_plastic_state();
};

}    // namespace LAMMPS_NS

#endif
#endif

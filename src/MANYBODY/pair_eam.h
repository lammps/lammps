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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eam,PairEAM);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_H
#define LMP_PAIR_EAM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEAM : public Pair {
 public:
  friend class FixSemiGrandCanonicalMC;    // Alex Stukowski option

  // public variables so ATC package can access them

  double cutmax;

  // potentials as array data

  int nrho, nr;
  int nfrho, nrhor, nz2r;
  double **frho, **rhor, **z2r;
  int *type2frho, **type2rhor, **type2z2r;

  // potentials in spline form used for force computation

  double dr, rdr, drho, rdrho, rhomax, rhomin;
  double ***rhor_spline, ***frho_spline, ***z2r_spline;

  PairEAM(class LAMMPS *);
  ~PairEAM() override;
  void compute(int, int) override;
  double compute_atomic_energy(int, NeighList *) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void swap_eam(double *, double **) override;

 protected:
  int nmax;    // allocated size of per-atom arrays
  double cutforcesq;
  double **scale;
  bigint embedstep;    // timestep, the embedding term was computed
  int he_flag;         // 1 if eam/he format: embedding table starts at rhomin
                       // (usually negative) with two-sided linear extrapolation

  // which potential file format coeff() reads:
  // one funcfl file per element (eam), or one setfl file for all
  // elements (eam/alloy), or one Finnis-Sinclair setfl variant file
  // with per-element-pair densities (eam/fs and eam/he)

  enum EAMFileFormat { FUNCFL, SETFL, FS };
  EAMFileFormat fileformat;

  int exceeded_rhomax;    // global flag for whether rho[i] has exceeded rhomax
                          // on a step energy is computed - 0 = no, 1 = yes

  // per-atom arrays

  double *rho, *fp;
  int *numforce;

  // potentials as file data

  struct Funcfl {
    char *file;
    int nrho, nr;
    double drho, dr, cut, mass;
    double *frho, *rhor, *zr;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, **rhor, ***z2r;
  };
  Setfl *setfl;

  struct Fs {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, ***rhor, ***z2r;
  };
  Fs *fs;

  virtual void allocate();
  virtual void array2spline();
  void interpolate(int, double, double *, double **);

  // coeff(), read_file(), and file2array() dispatch on fileformat
  // to the per-format implementations below

  virtual void read_file(char *);
  virtual void file2array();

  void coeff_funcfl(int, char **);
  void coeff_mapped(int, char **);
  void read_funcfl(char *);
  void read_setfl(char *);
  void read_fs(char *);
  void file2array_funcfl();
  void file2array_setfl();
  void file2array_fs();
  void delete_setfl();
  void delete_fs();

  // embedding spline table index m and fractional offset p for density rho_i
  // classic tables start at 0 and clamp at the table ends;
  // eam/he tables start at rhomin and allow extrapolation past both ends.
  // the HE template parameter lets hot loops hoist the format selection
  // out of the loop; the non-template overload dispatches on he_flag.

  template <int HE> inline void embedding_index(double rho_i, int &m, double &p) const
  {
    if (HE) {
      p = (rho_i - rhomin) * rdrho + 1.0;
      m = static_cast<int>(p);
      if (m < 2) {
        m = 2;
      } else if (m > nrho - 1) {
        m = nrho - 1;
      }
      p -= m;
      if (p < -1.0) {
        p = -1.0;
      } else if (p > 1.0) {
        p = 1.0;
      }
    } else {
      p = rho_i * rdrho + 1.0;
      m = static_cast<int>(p);
      m = MAX(1, MIN(m, nrho - 1));
      p -= m;
      p = MIN(p, 1.0);
    }
  }

  inline void embedding_index(double rho_i, int &m, double &p) const
  {
    if (he_flag)
      embedding_index<1>(rho_i, m, p);
    else
      embedding_index<0>(rho_i, m, p);
  }

  // embedding energy evaluation loop of compute(), templated on the
  // table format so the format test is hoisted out of the loop

  template <int HE> void compute_embedding(int, int &);
};

}    // namespace LAMMPS_NS

#endif
#endif

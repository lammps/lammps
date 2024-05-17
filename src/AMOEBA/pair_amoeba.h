/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(amoeba,PairAmoeba);
// clang-format on
#else

#ifndef LMP_PAIR_AMOEBA_H
#define LMP_PAIR_AMOEBA_H

#include "lmpfftsettings.h"    // IWYU pragma: export
#include "pair.h"

namespace LAMMPS_NS {

#define SBBITS15 29
#define NEIGHMASK15 0x1FFFFFFF

class PairAmoeba : public Pair {
 public:
  PairAmoeba(class LAMMPS *);
  ~PairAmoeba() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void finish() override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  void reset_grid() override;

  void pack_forward_grid(int, void *, int, int *) override;
  void unpack_forward_grid(int, void *, int, int *) override;
  void pack_reverse_grid(int, void *, int, int *) override;
  void unpack_reverse_grid(int, void *, int, int *) override;

  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;
  double memory_usage() override;

 protected:
  int nmax;                // allocation for owned+ghost
  int cfstyle, crstyle;    // style of forward/reverse comm operations
  int nualt;
  double electric;
  double rotate[3][3];    // rotation matrix

  bool amoeba;               // which force field: amoeba == true, hippo == false
  std::string mystyle;       // text label for style
  int first_flag;            // 1 before first init_style()
  int first_flag_compute;    // 1 before first call to compute()
  int optlevel;

  // turn on/off components of force field

  int hal_flag, repulse_flag, qxfer_flag;
  int disp_rspace_flag, disp_kspace_flag;
  int polar_rspace_flag, polar_kspace_flag;
  int mpole_rspace_flag, mpole_kspace_flag;
  int bond_flag, angle_flag, dihedral_flag, improper_flag;
  int urey_flag, pitorsion_flag, bitorsion_flag;

  // DEBUG timers

  double time_init, time_hal, time_repulse, time_disp;
  double time_mpole, time_induce, time_polar, time_qxfer;

  double time_mpole_rspace, time_mpole_kspace;
  double time_direct_rspace, time_direct_kspace;
  double time_mutual_rspace, time_mutual_kspace;
  double time_polar_rspace, time_polar_kspace;
  double time_grid_uind, time_fphi_uind;

  // energy/virial components

  double ehal, erepulse, edisp, epolar, empole, eqxfer;
  double virhal[6], virrepulse[6], virdisp[6], virpolar[6], virmpole[6], virqxfer[6];

  // scalar values defined in force-field file

  char *forcefield;    // FF name
  double am_dielectric;

  int opbendtype, vdwtype;
  int radius_rule, radius_type, radius_size, epsilon_rule;

  double bond_cubic, bond_quartic;
  double angle_cubic, angle_quartic, angle_pentic, angle_sextic;
  double opbend_cubic, opbend_quartic, opbend_pentic, opbend_sextic;
  double torsion_unit;

  int poltyp;

  double special_hal[8];
  double special_repel[8];
  double special_disp[8];
  double special_mpole[8];
  double special_polar_pscale[8];
  double special_polar_piscale[8];
  double special_polar_wscale[8];

  double polar_dscale, polar_uscale;

  // scalar values defined in keyfile

  double dhal, ghal;

  double vdwcut, vdwtaper;
  double repcut, reptaper;
  double dispcut, disptaper;
  double mpolecut, mpoletaper;
  double ctrncut, ctrntaper;

  double ewaldcut;
  double dewaldcut;
  double usolvcut;

  int use_ewald, use_dewald;

  int use_pred;
  int politer, polpred;
  int pcgprec, pcgguess;
  double pcgpeek;
  int tcgnab, optorder;
  int maxualt;
  double poleps;
  double udiag;

  int aeewald_key, apewald_key, adewald_key;
  int pmegrid_key, dpmegrid_key;

  // types and classes

  int n_amtype;       // # of defined AMOEBA types, 1-N
  int n_amclass;      // # of defined AMOEBA classes, 1-N
  int max_amtype;     // allocation length of per-type data
  int max_amclass;    // allocation length of per-class data

  int *amtype_defined;     // 1 if type was defined in FF file
  int *amclass_defined;    // 1 if class was defined in FF file
  int *amtype2class;       // amt2c[i] = class which type I belongs to

  // static per-atom properties, must persist as atoms migrate

  int index_amtype, index_amgroup, index_redID;
  int index_xyzaxis, index_polaxe, index_pval;

  int *amtype;     // AMOEBA type, 1 to N_amtype
  int *amgroup;    // AMOEBA polarization group, 1 to Ngroup

  char *id_pole, *id_udalt, *id_upalt;
  class FixStoreAtom *fixpole;     // stores pole = multipole components
  class FixStoreAtom *fixudalt;    // stores udalt = induced dipole history
  class FixStoreAtom *fixupalt;    // stores upalt = induced dipole history

  // static per-type properties defined in force-field file

  int *atomic_num;    // atomic number
  int *valence;       // valence (# of possible bonds)
  double *am_mass;    // atomic weight
  double *am_q;       // charge
  double **am_mu;     // dipole moment

  double *polarity;    // for polar
  double *pdamp;       // for polar
  double *thole;       // for polar
  double *dirdamp;     // for polar
  int *npolgroup;      // # of other types in polarization group, per-type
  int **polgroup;      // list of other types in polarization group, per-type

  double *sizpr, *dmppr, *elepr;

  // multipole frame info for each amtype, read from PRM file

  int *nmultiframe;                 // # of frames for each type
  int **mpaxis;                     // polaxe values
  int **xpole, **ypole, **zpole;    // other types in xyz dirs for multipole frame
  double ***fpole;                  // 13 values from file
                                    // 0 = monopole, same as q
                                    // 1,2,3 = 3 dipole components
                                    // 4-12 = 9 quadrupole components

  // static per-class properties defined in force-field file

  double *vdwl_eps;          // Vdwl epsilon for each class of atom
  double *vdwl_sigma;        // Vdwl sigma for each class of atom
  double *kred;              // fraction that H atoms move towards bonded atom
                             // used in Vdwl, 0.0 if not H atom
  double *csix, *adisp;      // used in dispersion
  double *chgct, *dmpct;     // used in charge transfer
  double *pcore, *palpha;    // for multipole

  int **vdwl_class_pair;      // Vdwl iclass/jclass for pair of classes
  double *vdwl_eps_pair;      // Vdwl epsilon for pair of classes
  double *vdwl_sigma_pair;    // Vdwl sigma for pair of classes
  int nvdwl_pair;             // # of pairwise Vdwl entries in file
  int max_vdwl_pair;          // size of allocated data for pairwise Vdwl

  // vectors and arrays of small size

  double *copt, *copm;    // 0:optorder in length
  double *gear, *aspc;

  double *a_ualt, *ap_ualt;                     // maxualt*(maxualt+1)/2 in length
  double *b_ualt, *bp_ualt;                     // maxualt in length
  double **c_ualt, **cp_ualt;                   // maxualt x maxualt in size
                                                // indices NOT flipped vs Fortran
  double *bpred, *bpredp, *bpreds, *bpredps;    // maxualt in length

  double vmsave[6];    // multipole virial saved to use in polar

  double csixpr;    // square of csix for all atoms

  // params common to pairwise terms

  double off2, cut2;
  double c0, c1, c2, c3, c4, c5;

  // Vdwl hal params - only for AMOEBA

  double **radmin, **epsilon;
  double **radmin4, **epsilon4;

  // peratom values computed each step
  // none of them persist with atoms
  // some of them need communication to ghosts

  double **rpole;    // multipole, comm to ghosts

  int *xaxis2local, *yaxis2local, *zaxis2local;    // xyz axis IDs -> local indices
                                                   // just for owned atoms
                                                   // set to self if not defined

  int *red2local;    // local indices of ired IDs, computed for owned and ghost
  double **xred;     // altered coords for H atoms for Vdwl, comm to ghosts

  double **tq;    // torque from pairwise multipole, reverse comm from ghosts

  double **uind, **uinp;    // computed by induce, comm to ghosts
  double **udirp;
  double **rsd, **rsdp;    // used by induce, comm to ghosts

  double **field, **fieldp;    // used by induce, reverse comm from ghosts
  double ***uopt, ***uoptp;    // Nlocal x Optorder+1 x 3 arrays

  double **ufld, **dufld;    // used by polar, reverse comm from ghosts
  double **zrsd, **zrsdp;    // used by induce, reverse comm from ghosts

  double ***uad, ***uap, ***ubd, ***ubp;    // used by TCG (not for now)

  double ***fopt, ***foptp;    // computed in induce, used by polar, if OPT
                               // Nlocal x optorder x 10

  double *poli;
  double **conj, **conjp;
  double **vec, **vecp;
  double **udir, **usum, **usump;

  double **fuind, **fuinp;
  double **fdip_phi1, **fdip_phi2, **fdip_sum_phi;
  double **dipfield1, **dipfield2;

  double **fphid, **fphip;
  double **fphidp, **cphidp;

  // derived local neighbor lists

  int *numneigh_dipole;         // number of dipole neighs for each atom
  int **firstneigh_dipole;      // ptr to each atom's dipole neigh indices
  MyPage<int> *ipage_dipole;    // pages of neighbor indices for dipole neighs

  double **firstneigh_dipdip;      // ptr to each atom's dip/dip values
  MyPage<double> *dpage_dipdip;    // pages of dip/dip values for dipole neighs

  int *numneigh_precond;         // number of precond neighs for each atom
  int **firstneigh_precond;      // ptr to each atom's precond neigh indices
  MyPage<int> *ipage_precond;    // pages of neighbor indices for precond neighs

  double **firstneigh_pcpc;      // ptr to each atom's pc/pc values
  MyPage<double> *dpage_pcpc;    // pages of pc/pc values for precond neighs

  // KSpace data
  // in indices = owned portion of grid in spatial decomp
  // out indices = in + ghost grid cells
  // fft indices = owned portion of grid in FFT decomp

  int nefft1, nefft2, nefft3;    // for electrostatic PME operations
  int ndfft1, ndfft2, ndfft3;    // for dispersion PME operations

  int bseorder;      // for electrostatics
  int bsporder;      // for polarization
  int bsdorder;      // for dispersion
  int bsordermax;    // max of 3 bsorder values

  double aewald;     // current Ewald alpha
  double aeewald;    // for electrostatics
  double apewald;    // for polarization
  double adewald;    // for dispersion

  double *bsmod1, *bsmod2, *bsmod3;    // B-spline module along abc axes
                                       // set to max of any nfft1,nfft2,nfft3

  double ***thetai1, ***thetai2, ***thetai3;    // B-spline coeffs along abc axes
                                                // Nlocal x max bsorder x 4

  int **igrid;    // grid indices for each owned particle, Nlocal x 3

  double **bsbuild;    // used internally in bsplgen, max-bsorder x max-bsorder
                       // indices ARE flipped vs Fortran

  // Kspace data for induce and polar

  double *qfac;        // convoulution pre-factors
  double *gridfft1;    // copy of p_kspace FFT grid

  double **cmp, **fmp;    // Cartesian and fractional multipoles
  double **cphi, **fphi;

  double *_moduli_array;    // buffers for moduli
  double *_moduli_bsarray;
  int _nfft_max;

  // params for current KSpace solve and FFT being worked on

  int nfft1, nfft2, nfft3;    // size of FFT
  int bsorder;                // stencil size
  double recip[3][3];         // indices NOT flipped vs Fortran
  double ctf[10][10];         // indices NOT flipped vs Fortran
  double ftc[10][10];         // indices NOT flipped vs Fortran

  class AmoebaConvolution *m_kspace;    // multipole KSpace
  class AmoebaConvolution *p_kspace;    // polar KSpace
  class AmoebaConvolution *pc_kspace;
  class AmoebaConvolution *d_kspace;    // dispersion KSpace
  class AmoebaConvolution *i_kspace;    // induce KSpace
  class AmoebaConvolution *ic_kspace;

  // FFT grid size factors

  int nfactors;    // # of factors
  int *factors;    // list of possible factors (2,3,5)

  // components of force field

  void hal();

  virtual void repulsion();
  void damprep(double, double, double, double, double, double, double, double, int, double, double,
               double *);

  void dispersion();
  virtual void dispersion_real();
  void dispersion_kspace();

  void multipole();
  virtual void multipole_real();
  void multipole_kspace();

  void polar();
  void polar_energy();
  virtual void polar_real();
  virtual void polar_kspace();
  void damppole(double, int, double, double, double *, double *, double *);

  virtual void induce();
  void ulspred();
  virtual void ufield0c(double **, double **);
  void uscale0b(int, double **, double **, double **, double **);
  void dfield0c(double **, double **);
  virtual void umutual1(double **, double **);
  virtual void umutual2b(double **, double **);
  void udirect1(double **);
  virtual void udirect2b(double **, double **);
  void dampmut(double, double, double, double *);
  void dampdir(double, double, double, double *, double *);
  void cholesky(int, double *, double *);

  void charge_transfer();

  // KSpace methods

  void lattice();
  void moduli();
  void bspline(double, int, double *);
  void dftmod(double *, double *, int, int);
  void bspline_fill();
  void bsplgen(double, double **);
  void cmp_to_fmp(double **, double **);
  void cart_to_frac();
  void fphi_to_cphi(double **, double **);
  void frac_to_cart();

  void grid_mpole(double **, FFT_SCALAR ***);
  void fphi_mpole(FFT_SCALAR ***, double **);
  void grid_uind(double **, double **, FFT_SCALAR ****);
  virtual void fphi_uind(FFT_SCALAR ****, double **, double **, double **);
  void grid_disp(FFT_SCALAR ***);

  void kewald();
  void kewald_parallel(int, int, int, int, int &, int &, int &, int &, int &, int &, int &, int &,
                       int &, int &, int &, int &, int &, int &, int &, int &, int &, int &);
  double ewaldcof(double);
  int factorable(int);

  double final_accuracy_mpole();
  double rms(int km, double prd, bigint natoms, double g_ewald, double q2);
  double two_charge_force;

  // debug methods

  FILE *fp_uind;
  void dump6(FILE *, const char *, double, double **, double **);

  // functions in pair_amoeba.cpp

  void allocate();
  void print_settings();

  void initialize_vdwl();
  void allocate_vdwl();
  void deallocate_vdwl();

  void initialize_smallsize();
  void allocate_smallsize();
  void deallocate_smallsize();

  void assign_groups();
  void pbc_xred();
  void precond_neigh();
  void choose(int);
  void mix();
  void zero_energy_force_virial();
  void grow_local();

  // functions in amoeba_utils.cpp

  void kmpole();

  void chkpole(int);
  void rotmat(int);
  void rotsite(int);

  void add_onefive_neighbors();
  void find_hydrogen_neighbors();
  void find_multipole_neighbors();

  void torque2force(int, double *, double *, double *, double *, double **);

  // functions in file_amoeba.cpp

  void set_defaults();
  void read_prmfile(char *);
  void read_keyfile(char *);

  void initialize_type_class();
  void allocate_type_class(int, int);
  void deallocate_type_class();

  void file_ffield(const std::vector<std::string> &, int);
  void file_literature(const std::vector<std::string> &, int);
  void file_atomtype(const std::vector<std::string> &, int);
  void file_vdwl(const std::vector<std::string> &, int);
  void file_vdwl_pair(const std::vector<std::string> &, int);
  void file_bstretch(const std::vector<std::string> &, int);
  void file_sbend(const std::vector<std::string> &, int);
  void file_abend(const std::vector<std::string> &, int);
  void file_pauli(const std::vector<std::string> &, int);
  void file_dispersion(const std::vector<std::string> &, int);
  void file_ub(const std::vector<std::string> &, int);
  void file_outplane(const std::vector<std::string> &, int);
  void file_torsion(const std::vector<std::string> &, int);
  void file_pitorsion(const std::vector<std::string> &, int);
  void file_multipole(const std::vector<std::string> &, int);
  void file_charge_penetration(const std::vector<std::string> &, int);
  void file_dippolar(const std::vector<std::string> &, int);
  void file_charge_transfer(const std::vector<std::string> &, int);

  // inline function for neighbor list unmasking

  inline int sbmask15(int j) const { return j >> SBBITS15 & 7; }
};
}    // namespace LAMMPS_NS
#endif
#endif

#ifndef LMP_MLIAP_SO3_H
#define LMP_MLIAP_SO3_H

#include "pointers.h"

namespace LAMMPS_NS {

class MLIAP_SO3 : protected Pointers {

 public:
  MLIAP_SO3(LAMMPS *, double vrcut, int vlmax, int vnmax, double valpha);
  MLIAP_SO3(LAMMPS *lmp) : Pointers(lmp){};

  ~MLIAP_SO3();

  void init();

  double memory_usage();

  int ncoeff;
  double *m_plist_r;
  double *m_dplist_r;

 private:
  double alloc_init, alloc_arrays;
  int *m_ellpl1, *m_ellm1;
  double *m_pfac;
  int m_pfac_l1, m_pfac_l2;
  double *m_dfac0, *m_dfac1, *m_dfac2, *m_dfac3, *m_dfac4, *m_dfac5;
  int m_dfac_l1, m_dfac_l2;
  double m_rcut, m_alpha;
  int m_lmax, m_nmax, m_Nmax;
  double *m_g_array, *m_w, *m_rootpq;
  int *m_idxu_block, *m_idxylm;
  int m_idxu_count, m_idxy_count;
  int m_numYlms;

  double *m_clisttot_r, *m_clisttot_i;

  double *m_rip_array, *m_rip_darray;
  double *m_sbes_array, *m_sbes_darray;
  int m_init_arrays;

  double *m_plist_i;
  double *m_clist_r, *m_clist_i;
  double *m_ulist_r, *m_ulist_i;

  double *m_Ylms_r, *m_Ylms_i, *m_dYlm_r, *m_dYlm_i;
  double *m_dplist_i, *m_dclist_r, *m_dclist_i, *m_tempdp_r;

 public:
  void spectrum(int nlocal, int *numneighs, int *jelems, double *wjelem, double **rij, int nmax,
                int lmax, double rcut, double alpha, int ncoefs);
  void spectrum_dxdr(int nlocal, int *numneighs, int *jelems, double *wjelem, double **rij,
                     int nmax, int lmax, double rcut, double alpha, bigint npairs, int ncoefs);

 private:
  double Cosine(double Rij, double Rc);
  double CosinePrime(double Rij, double Rc);
  double compute_sfac(double r, double rcut);
  double compute_dsfac(double r, double rcut);
  void get_sbes_array(int nlocal, int *numneighs, double **rij, int lmax, double rcut,
                      double alpha);
  void get_rip_array(int nlocal, int *numneighs, double **rij, int nmax, int lmax, double alpha);
  void init_arrays(int nlocal, int ncoefs);
  void init_garray(int nmax, int lmax, double rcut, double alpha, double *w, int lw1,
                   double *g_array, int lg2);

  void compute_uarray_recursive(double x, double y, double z, double r, int twol, double *ulist_r,
                                double *ulist_i, int *idxu_block, double *rootpqarray);
  void compute_ncoeff();

  int get_sum(int istart, int iend, int id, int imult);

  void compute_dpidrj(int nmax, int lmax, double *clisttot_r, double *clisttot_i, int lctot2,
                      double *dclist_r, double *dclist_i, int ldcli2, int ldcli3, double *dplist_r,
                      int dpli2);

  double compute_g(double r, int n, int nmax, double rcut, double *w, int lw1);
  double phi(double r, int alpha, double rcut);
  void compute_pi(int nmax, int lmax, double *clisttot_r, double *clisttot_i, int lcl2,
                  double *plist_r, double *plist_i, int lpl2, int indpl);

  void compute_W(int nmax, double *arr);
};

}    // namespace LAMMPS_NS

#endif

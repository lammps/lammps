/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This file is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_QMMM_XTB_ADAPTER_H
#define LMP_QMMM_XTB_ADAPTER_H

extern "C" {

// The adapter owns one GFN1-xTB or GFN2-xTB calculator on the calling MPI rank.  LAMMPS
// invokes it on rank zero and broadcasts the results, because libxtb uses
// shared-memory parallelism internally.
int lammps_qmmm_xtb_create(int nqm, const int *atomic_numbers, const double *qm_xyz_bohr,
                           int method, int charge, int uhf, double accuracy, int maxiter,
                           double electronic_temperature);

int lammps_qmmm_xtb_calculate(int nqm, const double *qm_xyz_bohr, int npoint,
                              const double *point_xyz_bohr, const double *point_charge,
                              const int *point_atomic_numbers, double mm_hardness,
                              const double *mm_shift_hartree, const double *image_response_hartree,
                              double *energy_hartree, double *qm_gradient_hartree_bohr,
                              double *mulliken_charge, double *point_gradient_hartree_bohr);

void lammps_qmmm_xtb_destroy();
}

#endif

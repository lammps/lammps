//
// Created by charlie sievers on 6/21/18.
//

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(dynamical_matrix/kk,DynamicalMatrixKokkos);
CommandStyle(dynamical_matrix/kk/device,DynamicalMatrixKokkos);
CommandStyle(dynamical_matrix/kk/host,DynamicalMatrixKokkos);
// clang-format on
#else

#ifndef LMP_DYNAMICAL_MATRIX_KOKKOS_H
#define LMP_DYNAMICAL_MATRIX_KOKKOS_H

#include "dynamical_matrix.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class DynamicalMatrixKokkos : public DynamicalMatrix {
 public:
  DynamicalMatrixKokkos(class LAMMPS *);
  virtual ~DynamicalMatrixKokkos();
  void command(int, char **);
  void setup();

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    f(i,0) += f_merge_copy(i,0);
    f(i,1) += f_merge_copy(i,1);
    f(i,2) += f_merge_copy(i,2);
  }

 protected:
  void update_force();
  void force_clear();
  DAT::t_f_array f_merge_copy,f;


};
}    // namespace LAMMPS_NS

#endif    //LMP_DYNAMICAL_MATRIX_KOKKOS_H
#endif

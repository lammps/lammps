// variable_kokkos.h
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

#ifndef LMP_VARIABLE_KOKKOS_H
#define LMP_VARIABLE_KOKKOS_H

#include "variable.h"
#include "kokkos_type.h"
#include <vector>

namespace LAMMPS_NS {

struct VarOpcode {
  int opcode;
  double value;
  int nstride;
  /// Index into packed per-atom aux buffer (fix/compute/custom columns); -1 if unused
  int aux_slot;

  KOKKOS_INLINE_FUNCTION
  VarOpcode() : opcode(0), value(0.0), nstride(1), aux_slot(-1) {}
};

class VariableKokkos : public Variable {
 public:
  VariableKokkos(class LAMMPS *);
  ~VariableKokkos() override;

  // Kokkos implementation of atom-style variables (delegates ATOMFILE to base class)
  void compute_atom(int, int, double *, int, int) override;

  // Bytecode buffers (explicit HostSpace: OpenMP default device has no HostMirror typedef)
  using OpcodeDeviceView =
      Kokkos::View<VarOpcode *, Kokkos::LayoutRight, typename KKDevice<LMPDeviceType>::value>;
  using OpcodeHostView = Kokkos::View<VarOpcode *, Kokkos::LayoutRight, Kokkos::HostSpace>;

 private:
  class AtomKokkos *atomKK;

  struct CodeInfo {
    OpcodeDeviceView d_code;
    OpcodeHostView h_code;
    int max_stack;
    int length;
    // Host column bases for ATOMARRAY nodes not mapped to core AtomKokkos views (e.g. f_fix[i])
    std::vector<const double *> aux_col_bases;
    std::vector<int> aux_nstride;
    Kokkos::View<double *, Kokkos::LayoutRight, typename KKDevice<LMPDeviceType>::value> d_aux_packed;
  };

  std::vector<CodeInfo> compiled_vars;

  // Compiles the Tree AST into a flat RPN bytecode array
  void compile_tree(Tree *tree, CodeInfo &info);
  void flatten_ast(Tree *tree, CodeInfo &info, int &current_stack, int &max_stack);
  
  // Identifies which atom array a raw pointer from the base Variable class corresponds to
  int map_atom_array(double *ptr);
  int map_int_array(int *ptr);
  int map_bigint_array(bigint *ptr);
};

} // namespace LAMMPS_NS

#endif // !LMP_VARIABLE_KOKKOS_H

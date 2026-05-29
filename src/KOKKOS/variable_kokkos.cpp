// variable_kokkos.cpp
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

#include "variable_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "group.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "universe.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

// Mirror the base class enums to interpret the Tree->type field
namespace {
enum {
  DONE, ADD, SUBTRACT, MULTIPLY, DIVIDE, CARAT, MODULO, UNARY,
  NOT, EQ, NE, LT, LE, GT, GE, AND, OR, XOR,
  SQRT, EXP, LN, LOG, ABS, SIN, COS, TAN, ASIN, ACOS, ATAN, ATAN2,
  RANDOM, NORMAL, CEIL, FLOOR, ROUND, TERNARY,
  RAMP, STAGGER, LOGFREQ, LOGFREQ2, LOGFREQ3, STRIDE, STRIDE2,
  VDISPLACE, SWIGGLE, CWIGGLE, SIGN, GMASK, RMASK,
  GRMASK, IS_ACTIVE, IS_DEFINED, IS_AVAILABLE, IS_FILE, EXTRACT_SETTING,
  PYWRAPPER,
  VALUE, ATOMARRAY, TYPEARRAY, INTARRAY, BIGINTARRAY, VECTORARRAY
};

// Kokkos-specific opcodes for identifying mapped atom properties on the GPU
enum {
  OP_KOKKOS_X = 1000, OP_KOKKOS_Y, OP_KOKKOS_Z,
  OP_KOKKOS_VX, OP_KOKKOS_VY, OP_KOKKOS_VZ,
  OP_KOKKOS_FX, OP_KOKKOS_FY, OP_KOKKOS_FZ,
  OP_KOKKOS_Q, OP_KOKKOS_MASS, OP_KOKKOS_RMASS, OP_KOKKOS_TYPE, OP_KOKKOS_ID,
  // Packed generic per-atom double column (fix/compute/custom), see CodeInfo::d_aux_packed
  OP_KOKKOS_ATOMARRAY_AUX = 1100
};
}

/* ---------------------------------------------------------------------- */

VariableKokkos::VariableKokkos(LAMMPS *lmp) : Variable(lmp)
{
  atomKK = dynamic_cast<AtomKokkos *>(atom);
}

/* ---------------------------------------------------------------------- */

VariableKokkos::~VariableKokkos()
{
  // Cleanup is handled automatically by Kokkos Views and std::vector
}

/* ----------------------------------------------------------------------
   Flatten the dynamically allocated AST into a contiguous Postfix array
------------------------------------------------------------------------- */

void VariableKokkos::flatten_ast(Tree *tree, CodeInfo &info, int &current_stack, int &max_stack)
{
  if (!tree) return;

  // Post-order traversal for RPN execution (evaluate children first)
  if (tree->first) {
    flatten_ast(tree->first, info, current_stack, max_stack);
  }
  if (tree->second) {
    flatten_ast(tree->second, info, current_stack, max_stack);
  }
  if (tree->nextra) {
    for (int i = 0; i < tree->nextra; i++) {
      flatten_ast(tree->extra[i], info, current_stack, max_stack);
    }
  }

  // Ensure capacity in the Host mirror view
  if (info.length >= info.h_code.extent(0)) {
    int new_size = info.h_code.extent(0) == 0 ? 64 : info.h_code.extent(0) * 2;
    Kokkos::resize(info.h_code, new_size);
  }

  VarOpcode &op = info.h_code(info.length++);
  op.opcode = tree->type;
  op.nstride = tree->nstride;
  op.value = 0.0;
  op.aux_slot = -1;

  // Map data values and array accesses
  if (tree->type == VALUE) {
    op.value = tree->value;
    current_stack++;
  }
  else if (tree->type == ATOMARRAY) {
    int mapped_op = map_atom_array(tree->array);
    if (mapped_op >= 0) {
      op.opcode = mapped_op;
    } else {
      if (!tree->array)
        error->all(FLERR, "VariableKokkos encountered null double ATOMARRAY");
      op.opcode = OP_KOKKOS_ATOMARRAY_AUX;
      op.aux_slot = static_cast<int>(info.aux_col_bases.size());
      info.aux_col_bases.push_back(tree->array);
      const int str = tree->nstride > 0 ? tree->nstride : 1;
      info.aux_nstride.push_back(str);
    }
    current_stack++;
  }
  else if (tree->type == TYPEARRAY) {
    int mapped_op = map_atom_array(tree->array);
    if (mapped_op >= 0) op.opcode = mapped_op;
    else error->all(FLERR, "VariableKokkos encountered unmapped TYPEARRAY");
    current_stack++;
  }
  else if (tree->type == INTARRAY) {
    int mapped_op = map_int_array(tree->iarray);
    if (mapped_op >= 0) op.opcode = mapped_op;
    else error->all(FLERR, "VariableKokkos encountered unmapped INTARRAY");
    current_stack++;
  }
  else if (tree->type == BIGINTARRAY) {
    int mapped_op = map_bigint_array(tree->barray);
    if (mapped_op >= 0) op.opcode = mapped_op;
    else error->all(FLERR, "VariableKokkos encountered unmapped BIGINTARRAY");
    current_stack++;
  } else {
    const int operand_count = (tree->first ? 1 : 0) + (tree->second ? 1 : 0) + tree->nextra;

    // All non-leaf operator/function nodes consume their operands and push
    // a single result onto the RPN stack.
    if (operand_count > 0) {
      if (current_stack < operand_count)
        error->all(FLERR, "VariableKokkos encountered invalid AST stack underflow");
      current_stack += 1 - operand_count;
    }
  }

  if (current_stack > max_stack) max_stack = current_stack;
}

/* ---------------------------------------------------------------------- */

void VariableKokkos::compile_tree(Tree *tree, CodeInfo &info)
{
  info.aux_col_bases.clear();
  info.aux_nstride.clear();
  info.length = 0;
  info.max_stack = 0;
  int current_stack = 0;

  flatten_ast(tree, info, current_stack, info.max_stack);

  // Resize device view and deep copy the compiled bytecode
  Kokkos::resize(info.d_code, info.length);
  Kokkos::deep_copy(
    Kokkos::subview(info.d_code, std::make_pair(0, info.length)),
    Kokkos::subview(info.h_code, std::make_pair(0, info.length))
  );
}

/* ----------------------------------------------------------------------
   Identify raw pointers to core atom arrays to convert them to Kokkos ops
------------------------------------------------------------------------- */

int VariableKokkos::map_atom_array(double *ptr)
{
  if (ptr == &atom->x[0][0]) return OP_KOKKOS_X;
  if (ptr == &atom->x[0][1]) return OP_KOKKOS_Y;
  if (ptr == &atom->x[0][2]) return OP_KOKKOS_Z;
  if (ptr == &atom->v[0][0]) return OP_KOKKOS_VX;
  if (ptr == &atom->v[0][1]) return OP_KOKKOS_VY;
  if (ptr == &atom->v[0][2]) return OP_KOKKOS_VZ;
  if (ptr == &atom->f[0][0]) return OP_KOKKOS_FX;
  if (ptr == &atom->f[0][1]) return OP_KOKKOS_FY;
  if (ptr == &atom->f[0][2]) return OP_KOKKOS_FZ;
  if (ptr == atom->q) return OP_KOKKOS_Q;
  if (ptr == atom->mass) return OP_KOKKOS_MASS;
  if (ptr == atom->rmass) return OP_KOKKOS_RMASS;
  return -1;
}

int VariableKokkos::map_int_array(int *ptr)
{
  if (ptr == atom->type) return OP_KOKKOS_TYPE;
  if (sizeof(tagint) == sizeof(smallint) && (tagint*)ptr == atom->tag) return OP_KOKKOS_ID;
  return -1;
}

int VariableKokkos::map_bigint_array(bigint *ptr)
{
  if (sizeof(tagint) == sizeof(bigint) && (tagint*)ptr == atom->tag) return OP_KOKKOS_ID;
  return -1;
}

/* ----------------------------------------------------------------------
   Functor to evaluate the flattened RPN bytecode on the Device
------------------------------------------------------------------------- */

template<int SUMFLAG>
struct EvalAtomVarFunctor {
  using DAT = ArrayTypes<LMPDeviceType>;

  VariableKokkos::OpcodeDeviceView code;
  int code_len;

  Kokkos::View<double *, Kokkos::LayoutRight, typename KKDevice<LMPDeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      result;
  int stride;
  int groupbit;

  typename DAT::t_kkfloat_1d_3_lr_randomread x;
  typename DAT::t_kkfloat_1d_3_randomread v;
  typename DAT::t_kkacc_1d_3_randomread f;
  typename DAT::t_kkfloat_1d_randomread q;
  typename DAT::t_kkfloat_1d_randomread mass;
  typename DAT::t_kkfloat_1d_randomread rmass;
  typename DAT::t_int_1d_randomread type;
  typename DAT::t_tagint_1d_randomread tag;
  typename DAT::t_int_1d_randomread mask;

  Kokkos::View<const double *, Kokkos::LayoutRight, typename KKDevice<LMPDeviceType>::value> aux_packed;
  int aux_nlocal;

  EvalAtomVarFunctor(VariableKokkos::OpcodeDeviceView _code, int _code_len,
                       Kokkos::View<double *, Kokkos::LayoutRight, typename KKDevice<LMPDeviceType>::value,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>> _result,
                       int _stride, int _groupbit, AtomKokkos *atomKK,
                       Kokkos::View<const double *, Kokkos::LayoutRight,
                                    typename KKDevice<LMPDeviceType>::value> _aux_packed,
                       int _aux_nlocal)
      : code(_code),
        code_len(_code_len),
        result(_result),
        stride(_stride),
        groupbit(_groupbit),
        aux_packed(_aux_packed),
        aux_nlocal(_aux_nlocal)
  {
    x = atomKK->k_x.view<LMPDeviceType>();
    v = atomKK->k_v.view<LMPDeviceType>();
    f = atomKK->k_f.view<LMPDeviceType>();
    q = atomKK->k_q.view<LMPDeviceType>();
    mass = atomKK->k_mass.view<LMPDeviceType>();
    rmass = atomKK->k_rmass.view<LMPDeviceType>();
    type = atomKK->k_type.view<LMPDeviceType>();
    tag = atomKK->k_tag.view<LMPDeviceType>();
    mask = atomKK->k_mask.view<LMPDeviceType>();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (!(mask(i) & groupbit)) {
      if (!SUMFLAG) result(i * stride) = 0.0;
      return;
    }

    double stack[64]; // Fixed bounds for hardware registers. Assumes max_stack <= 64.
    int sp = 0;

    for (int pc = 0; pc < code_len; ++pc) {
      const VarOpcode& op = code(pc);

      switch (op.opcode) {
        case VALUE: stack[sp++] = op.value; break;

        // Kokkos Mapped Arrays
        case OP_KOKKOS_X: stack[sp++] = x(i, 0); break;
        case OP_KOKKOS_Y: stack[sp++] = x(i, 1); break;
        case OP_KOKKOS_Z: stack[sp++] = x(i, 2); break;
        case OP_KOKKOS_VX: stack[sp++] = v(i, 0); break;
        case OP_KOKKOS_VY: stack[sp++] = v(i, 1); break;
        case OP_KOKKOS_VZ: stack[sp++] = v(i, 2); break;
        case OP_KOKKOS_FX: stack[sp++] = f(i, 0); break;
        case OP_KOKKOS_FY: stack[sp++] = f(i, 1); break;
        case OP_KOKKOS_FZ: stack[sp++] = f(i, 2); break;
        case OP_KOKKOS_Q: stack[sp++] = q(i); break;
        case OP_KOKKOS_TYPE: stack[sp++] = (double)type(i); break;
        case OP_KOKKOS_ID: stack[sp++] = (double)tag(i); break;
        case OP_KOKKOS_RMASS: stack[sp++] = rmass(i); break;
        case OP_KOKKOS_MASS: stack[sp++] = mass(type(i)); break;

        case OP_KOKKOS_ATOMARRAY_AUX:
          stack[sp++] = aux_packed(op.aux_slot * aux_nlocal + i);
          break;

        // Binary Math
        case ADD: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = a + b;
          break;
        }
        case SUBTRACT: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = a - b;
          break;
        }
        case MULTIPLY: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = a * b;
          break;
        }
        case DIVIDE: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = a / b;
          break;
        }
        case MODULO: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = fmod(a, b);
          break;
        }
        case CARAT: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = pow(a, b);
          break;
        }

        // Logic
        case EQ: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a == b) ? 1.0 : 0.0;
          break;
        }
        case NE: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a != b) ? 1.0 : 0.0;
          break;
        }
        case LT: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a < b) ? 1.0 : 0.0;
          break;
        }
        case LE: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a <= b) ? 1.0 : 0.0;
          break;
        }
        case GT: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a > b) ? 1.0 : 0.0;
          break;
        }
        case GE: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a >= b) ? 1.0 : 0.0;
          break;
        }
        case AND: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a != 0.0 && b != 0.0) ? 1.0 : 0.0;
          break;
        }
        case OR: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a != 0.0 || b != 0.0) ? 1.0 : 0.0;
          break;
        }
        case XOR: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] =
              ((a == 0.0 && b != 0.0) || (a != 0.0 && b == 0.0)) ? 1.0 : 0.0;
          break;
        }

        // Unary Math
        case UNARY: stack[sp-1] = -stack[sp-1]; break;
        case NOT: stack[sp-1] = (stack[sp-1] == 0.0) ? 1.0 : 0.0; break;
        case SQRT: stack[sp-1] = sqrt(stack[sp-1]); break;
        case EXP: stack[sp-1] = exp(stack[sp-1]); break;
        case LN: stack[sp-1] = log(stack[sp-1]); break;
        case LOG: stack[sp-1] = log10(stack[sp-1]); break;
        case ABS: stack[sp-1] = fabs(stack[sp-1]); break;
        case SIN: stack[sp-1] = sin(stack[sp-1]); break;
        case COS: stack[sp-1] = cos(stack[sp-1]); break;
        case TAN: stack[sp-1] = tan(stack[sp-1]); break;
        case ASIN: stack[sp-1] = asin(stack[sp-1]); break;
        case ACOS: stack[sp-1] = acos(stack[sp-1]); break;
        case ATAN: stack[sp-1] = atan(stack[sp-1]); break;

        // Complex / Multi-pop
        case ATAN2: {
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = atan2(a, b);
          break;
        }
        case TERNARY: {
          double extra = stack[--sp];
          double b = stack[--sp];
          double a = stack[--sp];
          stack[sp++] = (a != 0.0) ? b : extra;
          break;
        }
      }
    }

    if (SUMFLAG) result(i * stride) += stack[0];
    else result(i * stride) = stack[0];
  }
};

/* ----------------------------------------------------------------------
   compute result of atom-style variable evaluation directly on device
------------------------------------------------------------------------- */

void VariableKokkos::compute_atom(int ivar, int igroup, double *result_ptr, int stride, int sumflag)
{
  if (ivar < 0 || ivar >= maxvar) return;
  if (style[ivar] != ATOM) {
    // Fall back to CPU evaluation for ATOMFILE or other unhandled styles
    Variable::compute_atom(ivar, igroup, result_ptr, stride, sumflag);
    return;
  }

  if (eval_in_progress[ivar]) print_var_error(FLERR, "has a circular dependency", ivar);
  eval_in_progress[ivar] = 1;

  try {
    // Input is instantiated before create(), so `atom` is not AtomKokkos yet at ctor time.
    if (!atomKK) atomKK = dynamic_cast<AtomKokkos *>(atom);
    if (!atomKK) {
      eval_in_progress[ivar] = 0;
      Variable::compute_atom(ivar, igroup, result_ptr, stride, sumflag);
      return;
    }

    // The base variable class parses the string into a Tree AST every time.
    // We mimic this, then flatten the tree to Kokkos bytecode.
    Tree *tree = nullptr;
    treetype = ATOM;
    evaluate(data[ivar][0], &tree, ivar);
    collapse_tree(tree);

    if (ivar >= compiled_vars.size()) compiled_vars.resize(ivar + 1);
    CodeInfo &info = compiled_vars[ivar];
    compile_tree(tree, info);

    free_tree(tree);

    if (info.max_stack > 64) {
      error->all(FLERR, "VariableKokkos encountered AST depth exceeding device stack allocation");
    }

    if (result_ptr) {
      // Sync necessary atom data to the device before evaluating
      atomKK->sync(ExecutionSpaceFromDevice<LMPDeviceType>::space,
                   X_MASK | V_MASK | F_MASK | Q_MASK | TYPE_MASK | MASK_MASK);

      int nlocal = atomKK->nlocal;
      int groupbit = group->bitmask[igroup];

      const int na = static_cast<int>(info.aux_col_bases.size());
      if (na > 0) {
        Kokkos::resize(info.d_aux_packed, na * nlocal);
        auto h_aux = Kokkos::create_mirror_view(info.d_aux_packed);
        for (int a = 0; a < na; ++a) {
          const double *col = info.aux_col_bases[a];
          const int str = info.aux_nstride[a] > 0 ? info.aux_nstride[a] : 1;
          for (int i = 0; i < nlocal; ++i) h_aux(a * nlocal + i) = col[i * str];
        }
        Kokkos::deep_copy(info.d_aux_packed, h_aux);
      } else {
        Kokkos::resize(info.d_aux_packed, 0);
      }

      Kokkos::View<const double *, Kokkos::LayoutRight, typename KKDevice<LMPDeviceType>::value> aux_read(
          info.d_aux_packed);

      // Wrap host buffer as device-accessible unmanaged view (MemoryTraits last).
      Kokkos::View<double *, Kokkos::LayoutRight, typename KKDevice<LMPDeviceType>::value,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          d_result(result_ptr, nlocal * stride);

      if (sumflag == 0) {
        EvalAtomVarFunctor<0> functor(info.d_code, info.length, d_result, stride, groupbit, atomKK, aux_read,
                                      nlocal);
        Kokkos::parallel_for("VariableKokkos::compute_atom", nlocal, functor);
      } else {
        EvalAtomVarFunctor<1> functor(info.d_code, info.length, d_result, stride, groupbit, atomKK, aux_read,
                                      nlocal);
        Kokkos::parallel_for("VariableKokkos::compute_atom", nlocal, functor);
      }

      // Ensure GPU completes execution
      Kokkos::fence();
    }
  } catch (...) {
    eval_in_progress[ivar] = 0;
    throw;
  }

  eval_in_progress[ivar] = 0;
}

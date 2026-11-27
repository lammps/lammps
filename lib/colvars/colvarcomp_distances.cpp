// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>

#include "colvardeps.h"
#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvar_rotation_derivative.h"

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
#include "cuda/colvarcomp_distance_kernel.h"
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)


colvar::distance::distance()
{
  set_function_type("distance");
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
}


int colvar::distance::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  if (!group1 || !group2) {
    return error_code | COLVARS_INPUT_ERROR;
  }

  error_code |= init_total_force_params(conf);

  return error_code;
}


void colvar::distance::calc_value()
{
  dist_v = boundary_conditions.position_distance(group1->center_of_mass(), group2->center_of_mass());
  x.real_value = dist_v.norm();
}


void colvar::distance::calc_gradients()
{
  cvm::rvector const u = dist_v.unit();
  group1->set_weighted_gradient(-1.0 * u);
  group2->set_weighted_gradient(       u);
}


void colvar::distance::calc_force_invgrads()
{
  group1->read_total_forces();
  if (is_enabled(f_cvc_one_site_total_force)) {
    ft.real_value = -1.0 * (group1->total_force() * dist_v.unit());
  } else {
    group2->read_total_forces();
    ft.real_value = 0.5 * ((group2->total_force() - group1->total_force()) * dist_v.unit());
  }
}


void colvar::distance::calc_Jacobian_derivative()
{
  jd.real_value = x.real_value ? (2.0 / x.real_value) : 0.0;
}



colvar::distance_vec::distance_vec()
{
  set_function_type("distanceVec");
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_3vector);
}


void colvar::distance_vec::calc_value()
{
  x.rvector_value =
      boundary_conditions.position_distance(group1->center_of_mass(), group2->center_of_mass());
}


void colvar::distance_vec::calc_gradients()
{
  // gradients are not stored: a 3x3 matrix for each atom would be
  // needed to store just the identity matrix
}


void colvar::distance_vec::apply_force(colvarvalue const &force)
{
  if (!group1->noforce)
    group1->apply_force(-1.0 * force.rvector_value);

  if (!group2->noforce)
    group2->apply_force(       force.rvector_value);
}


cvm::real colvar::distance_vec::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  return (boundary_conditions.position_distance(x1.rvector_value, x2.rvector_value)).norm2();
}


colvarvalue colvar::distance_vec::dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return 2.0 * boundary_conditions.position_distance(x2.rvector_value, x1.rvector_value);
}


colvarvalue colvar::distance_vec::dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return distance_vec::dist2_lgrad(x2, x1);
}


void colvar::distance_vec::wrap(colvarvalue & /* x_unwrapped */) const {}


colvar::distance_z::distance_z()
{
  set_function_type("distanceZ");
  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
  provide(f_cvc_periodic);
  x.type(colvarvalue::type_scalar);
}


int colvar::distance_z::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  main = parse_group(conf, "main");
  ref1 = parse_group(conf, "ref");
  // this group is optional
  ref2 = parse_group(conf, "ref2", true);

  if ( ref2 ) {
    cvmodule->log("Using axis joining the centers of mass of groups \"ref\" and \"ref2\"\n");
    fixed_axis = false;
    if (key_lookup(conf, "axis"))
      cvmodule->log("Warning: explicit axis definition will be ignored!\n");
  } else {
    if (get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0))) {
      if (axis.norm2() == 0.0) {
        error_code |= cvmodule->error("Axis vector is zero!", COLVARS_INPUT_ERROR);
      }
      if (axis.norm2() != 1.0) {
        axis = axis.unit();
        cvmodule->log("The normalized axis is: "+cvm::to_str(axis)+".\n");
      }
    }
    fixed_axis = true;
  }

  error_code |= init_total_force_params(conf);

  return error_code;
}


void colvar::distance_z::calc_value()
{
  cvm::rvector const M = main->center_of_mass();
  cvm::rvector const R1 = ref1->center_of_mass();
  if (fixed_axis) {
    dist_v = boundary_conditions.position_distance(R1, M);
  } else {
    cvm::rvector const R2 = ref2->center_of_mass();
    cvm::rvector const C = 0.5 * (R1 + R2);
    dist_v = boundary_conditions.position_distance(C, M);
    axis = boundary_conditions.position_distance(R1, R2);
    axis_norm = axis.norm();
    axis = axis.unit();
  }
  x.real_value = axis * dist_v;
  wrap(x);
}


void colvar::distance_z::calc_gradients()
{
  main->set_weighted_gradient( axis );

  if (fixed_axis) {
    ref1->set_weighted_gradient(-1.0 * axis);
  } else {
    // Perpendicular term: (M - C) - Â(Â·(M - C)) = (M - C) - Â x
    cvm::rvector const perp_term = dist_v - axis * x.real_value;
    // ∂x/∂R₁ = -common_term/‖A‖ - (1/2)Â
    cvm::rvector const grad_R1 = -perp_term / axis_norm - 0.5 * axis;
    // ∂x/∂R₂ = common_term/‖A‖ - (1/2)Â
    cvm::rvector const grad_R2 = perp_term / axis_norm - 0.5 * axis;
    ref1->set_weighted_gradient(grad_R1);
    ref2->set_weighted_gradient(grad_R2);
  }
}


void colvar::distance_z::calc_force_invgrads()
{
  main->read_total_forces();

  if (fixed_axis && !is_enabled(f_cvc_one_site_total_force)) {
    ref1->read_total_forces();
    ft.real_value = 0.5 * ((main->total_force() - ref1->total_force()) * axis);
  } else {
    ft.real_value = main->total_force() * axis;
  }
}


void colvar::distance_z::calc_Jacobian_derivative()
{
  jd.real_value = 0.0;
}



colvar::distance_xy::distance_xy()
{
  set_function_type("distanceXY");
  provide(f_cvc_periodic, false); // Disable inherited distance_z flag
  init_as_distance();
}


void colvar::distance_xy::calc_value()
{
  dist_v = boundary_conditions.position_distance(ref1->center_of_mass(), main->center_of_mass());
  if (!fixed_axis) {
    v12 = boundary_conditions.position_distance(ref1->center_of_mass(), ref2->center_of_mass());
    axis_norm = v12.norm();
    axis = v12.unit();
  }

  dist_v_ortho = dist_v - (dist_v * axis) * axis;
  x.real_value = dist_v_ortho.norm();
}


void colvar::distance_xy::calc_gradients()
{
  // Intermediate quantity (r_P3 / r_12 where P is the projection
  // of 3(main) on the plane orthogonal to 12, containing 1 (ref1))
  cvm::real A;
  cvm::real x_inv;

  if (x.real_value == 0.0) return;
  x_inv = 1.0 / x.real_value;

  if (fixed_axis) {
    ref1->set_weighted_gradient(-1.0 * x_inv * dist_v_ortho);
    main->set_weighted_gradient(       x_inv * dist_v_ortho);
  } else {
    v13 = boundary_conditions.position_distance(ref1->center_of_mass(), main->center_of_mass());
    A = (dist_v * axis) / axis_norm;

    ref1->set_weighted_gradient( (A - 1.0) * x_inv * dist_v_ortho);
    ref2->set_weighted_gradient( -A        * x_inv * dist_v_ortho);
    main->set_weighted_gradient(      1.0  * x_inv * dist_v_ortho);
  }
}


void colvar::distance_xy::calc_force_invgrads()
{
  main->read_total_forces();

  if (fixed_axis && !is_enabled(f_cvc_one_site_total_force)) {
    ref1->read_total_forces();
    ft.real_value = 0.5 / x.real_value * ((main->total_force() - ref1->total_force()) * dist_v_ortho);
  } else {
    ft.real_value = 1.0 / x.real_value * main->total_force() * dist_v_ortho;
  }
}


void colvar::distance_xy::calc_Jacobian_derivative()
{
  jd.real_value = x.real_value ? (1.0 / x.real_value) : 0.0;
}



colvar::distance_dir::distance_dir()
{
  set_function_type("distanceDir");
  enable(f_cvc_com_based);
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_unit3vector);
}


void colvar::distance_dir::calc_value()
{
  dist_v =
      boundary_conditions.position_distance(group1->center_of_mass(), group2->center_of_mass());
  x.rvector_value = dist_v.unit();
}


void colvar::distance_dir::calc_gradients()
{
  // gradients are computed on the fly within apply_force()
  // Note: could be a problem if a future bias relies on gradient
  // calculations...
  // TODO in new deps system: remove dependency of biasing force to gradient?
  // That way we could tell apart an explicit gradient dependency
}


void colvar::distance_dir::apply_force(colvarvalue const &force)
{
  // remove the radial force component
  cvm::real const iprod = force.rvector_value * x.rvector_value;
  cvm::rvector const force_tang = force.rvector_value - iprod * x.rvector_value;

  if (!group1->noforce) {
    group1->apply_force(-1.0 / dist_v.norm() * force_tang);
  }
  if (!group2->noforce) {
    group2->apply_force( 1.0 / dist_v.norm() * force_tang);
  }
}


cvm::real colvar::distance_dir::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  return x1.dist2(x2);
}


colvarvalue colvar::distance_dir::dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return x1.dist2_grad(x2);
}


colvarvalue colvar::distance_dir::dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return x2.dist2_grad(x1);
}


void colvar::distance_dir::wrap(colvarvalue & /* x_unwrapped */) const {}



colvar::distance_inv::distance_inv()
{
  set_function_type("distanceInv");
  init_as_distance();
}


int colvar::distance_inv::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  get_keyval(conf, "exponent", exponent, exponent);
  if (exponent % 2) {
    error_code |=
        cvmodule->error("Error: odd exponent provided, can only use even ones.\n", COLVARS_INPUT_ERROR);
  }
  if (exponent <= 0) {
    error_code |= cvmodule->error("Error: negative or zero exponent provided.\n", COLVARS_INPUT_ERROR);
  }

  if (cvm::atom_group::overlap(*group1, *group2) != 0) {
    error_code |= cvmodule->error("Error: group1 and group2 have some atoms in common: this is not "
                            "allowed for distanceInv.\n",
                            COLVARS_INPUT_ERROR);
  }

  if (is_enabled(f_cvc_debug_gradient)) {
    cvmodule->log("Warning: debugGradients will not give correct results "
             "for distanceInv, because its value and gradients are computed "
             "simultaneously.\n");
  }

  return error_code;
}


void colvar::distance_inv::calc_value()
{
#define CALL_KERNEL() do {                             \
  const int factor = -1*(exponent/2);                  \
  for (size_t i = 0; i < group1->size(); ++i) {        \
    const cvm::atom_pos pos1(group1->pos_x(i),         \
                             group1->pos_y(i),         \
                             group1->pos_z(i));        \
    cvm::rvector g1(0, 0, 0);                          \
    for (size_t j = 0; j < group2->size(); ++j) {      \
      const cvm::atom_pos pos2(group2->pos_x(j),       \
                               group2->pos_y(j),       \
                               group2->pos_z(j));      \
      cvm::rvector dv;                                 \
      dv = boundary_conditions.position_distance(pos1, pos2); \
      cvm::real const d2 = dv.norm2();                                      \
      cvm::real const dinv = cvm::integer_power(d2, factor);                \
      x.real_value += dinv;                                                 \
      cvm::rvector const dsumddv = factor * dinv/d2 * 2.0 * dv;             \
      g1 += -1.0 * dsumddv;             \
      group2->grad_x(j) += dsumddv.x;   \
      group2->grad_y(j) += dsumddv.y;   \
      group2->grad_z(j) += dsumddv.z;   \
    }                                   \
    group1->grad_x(i) += g1.x; \
    group1->grad_y(i) += g1.y; \
    group1->grad_z(i) += g1.z; \
  }                            \
} while (0);
  x.real_value = 0.0;
  CALL_KERNEL();

  x.real_value *= 1.0 / cvm::real(group1->size() * group2->size());
  x.real_value = cvm::pow(x.real_value, -1.0/cvm::real(exponent));

  cvm::real const dxdsum = (-1.0/(cvm::real(exponent))) *
    cvm::integer_power(x.real_value, exponent+1) /
    cvm::real(group1->size() * group2->size());
  const size_t num_grads_xyz_group1 = 3 * group1->size();
  for (size_t i = 0; i < num_grads_xyz_group1; ++i) {
    group1->grad_x(i) = group1->grad_x(i) * dxdsum;
  }
  const size_t num_grads_xyz_group2 = 3 * group2->size();
  for (size_t i = 0; i < num_grads_xyz_group2; ++i) {
    group2->grad_x(i) = group2->grad_x(i) * dxdsum;
  }
#undef CALL_KERNEL
}


void colvar::distance_inv::calc_gradients()
{
}



colvar::distance_pairs::distance_pairs()
{
  set_function_type("distancePairs");
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_vector);
}


int colvar::distance_pairs::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");
  x.vector1d_value.resize(group1->size() * group2->size());

  return error_code;
}


void colvar::distance_pairs::calc_value()
{
  x.vector1d_value.resize(group1->size() * group2->size());
#define CALL_KERNEL() do {                                         \
  for (size_t i1 = 0; i1 < group1->size(); ++i1) {                 \
    const cvm::atom_pos pos1(group1->pos_x(i1),                    \
                             group1->pos_y(i1),                    \
                             group1->pos_z(i1));                   \
    cvm::rvector g1(0, 0, 0);                                      \
    for (size_t i2 = 0; i2 < group2->size(); i2++) {               \
      const cvm::atom_pos pos2(group2->pos_x(i2),                  \
                               group2->pos_y(i2),                  \
                               group2->pos_z(i2));                 \
      const cvm::rvector dv = boundary_conditions.position_distance(pos1, pos2); \
      cvm::real const d = dv.norm();                               \
      x.vector1d_value[i1*group2->size() + i2] = d;                \
      const cvm::rvector g2 = dv.unit();                           \
      g1 += -g2;                                                   \
      /*group2->grad_x(i2) += g2.x;                                  \
      group2->grad_y(i2) += g2.y;                                  \
      group2->grad_z(i2) += g2.z;*/                                  \
    }                                                              \
    /*group1->grad_x(i1) += g1.x;                                    \
    group1->grad_y(i1) += g1.y;                                    \
    group1->grad_z(i1) += g1.z;*/                                    \
  }                                                                \
} while (0);
  CALL_KERNEL();
#undef CALL_KERNEL
}


void colvar::distance_pairs::calc_gradients()
{
  // will be calculated on the fly in apply_force()
}


void colvar::distance_pairs::apply_force(colvarvalue const &force)
{
#define CALL_KERNEL() do {                                              \
  auto group1_force_obj = group1->get_group_force_object();             \
  auto group2_force_obj = group2->get_group_force_object();             \
  for (size_t i1 = 0; i1 < group1->size(); i1++) {                      \
    const cvm::atom_pos pos1(group1->pos_x(i1),                         \
                             group1->pos_y(i1),                         \
                             group1->pos_z(i1));                        \
    cvm::rvector f1(0, 0, 0);                                           \
    for (size_t i2 = 0; i2 < group2->size(); i2++) {                    \
      const cvm::atom_pos pos2(group2->pos_x(i2),                       \
                               group2->pos_y(i2),                       \
                               group2->pos_z(i2));                      \
      const cvm::rvector dv = boundary_conditions.position_distance(pos1, pos2); \
      cvm::real const d = dv.norm();                                    \
      x.vector1d_value[i1*group2->size() + i2] = d;                     \
      const cvm::rvector f2 = force[i1*group2->size() + i2] * dv.unit();\
      f1 += -f2;                                                        \
      group2_force_obj.add_atom_force(i2, f2);                          \
    }                                                                   \
    group1_force_obj.add_atom_force(i1, f1);                            \
  }                                                                     \
} while (0);
  CALL_KERNEL();
#undef CALL_KERNEL
}


cvm::real colvar::distance_pairs::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  return (x1.vector1d_value - x2.vector1d_value).norm2();
}


colvarvalue colvar::distance_pairs::dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return 2.0 * (x1.vector1d_value - x2.vector1d_value);
}


colvarvalue colvar::distance_pairs::dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return distance_pairs::dist2_lgrad(x1, x2);
}


void colvar::distance_pairs::wrap(colvarvalue & /* x_unwrapped */) const {}



colvar::dipole_magnitude::dipole_magnitude()
{
  set_function_type("dipoleMagnitude");
  x.type(colvarvalue::type_scalar);
}


int colvar::dipole_magnitude::init(std::string const &conf)
{
  int error_code = cvc::init(conf);
  atoms = parse_group(conf, "atoms");
  if (!atoms) error_code |= COLVARS_INPUT_ERROR;
  return error_code;
}


void colvar::dipole_magnitude::calc_value()
{
  cvm::atom_pos const atomsCom = atoms->center_of_mass();
  atoms->calc_dipole(atomsCom);
  dipoleV = atoms->dipole();
  x.real_value = dipoleV.norm();
}


void colvar::dipole_magnitude::calc_gradients()
{
  cvm::real const aux1 = atoms->total_charge/atoms->total_mass;
  cvm::atom_pos const dipVunit = dipoleV.unit();
  for (size_t i = 0; i < atoms->size(); ++i) {
    const cvm::rvector grad = (atoms->charge(i) - aux1 * atoms->mass(i)) * dipVunit;
    atoms->grad_x(i) = grad.x;
    atoms->grad_y(i) = grad.y;
    atoms->grad_z(i) = grad.z;
  }
}



colvar::gyration::gyration()
{
  set_function_type("gyration");
  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  init_as_distance();
}


int colvar::gyration::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  atoms = parse_group(conf, "atoms");

  if (atoms->b_user_defined_fit) {
    cvmodule->log("WARNING: explicit fitting parameters were provided for atom group \"atoms\".\n");
  } else {
    atoms->enable(f_ag_center);
    std::vector<cvm::atom_pos> ref_pos_aos{cvm::atom_pos(0, 0, 0)};
    atoms->set_ref_pos_from_aos(ref_pos_aos);
    atoms->fit_gradients.assign(3 * atoms->size(), 0);
  }

  return error_code;
}


void colvar::gyration::calc_value()
{
  x.real_value = 0.0;
  const size_t num_pos_xyz = 3 * atoms->size();
  for (size_t i = 0; i < num_pos_xyz; ++i) {
    x.real_value += atoms->pos_x(i) * atoms->pos_x(i);
  }
  x.real_value = cvm::sqrt(x.real_value / cvm::real(atoms->size()));
}


void colvar::gyration::calc_gradients()
{
  cvm::real const drdx = 1.0/(cvm::real(atoms->size()) * x.real_value);
  const size_t num_pos_xyz = 3 * atoms->size();
  for (size_t i = 0; i < num_pos_xyz; ++i) {
    atoms->grad_x(i) = drdx * atoms->pos_x(i);
  }
}


void colvar::gyration::calc_force_invgrads()
{
  atoms->read_total_forces();

  cvm::real const dxdr = 1.0/x.real_value;
  ft.real_value = 0.0;
  const size_t num_pos_xyz = 3 * atoms->size();
  for (size_t i = 0; i < num_pos_xyz; ++i) {
    ft.real_value += atoms->pos_x(i) * atoms->total_force_x(i);
  }
  ft.real_value *= dxdr;
}


void colvar::gyration::calc_Jacobian_derivative()
{
  jd = x.real_value ? (3.0 * cvm::real(atoms->size()) - 4.0) / x.real_value : 0.0;
}



colvar::inertia::inertia()
{
  set_function_type("inertia");
}


void colvar::inertia::calc_value()
{
  x.real_value = 0.0;
  const size_t num_pos_xyz = 3 * atoms->size();
  for (size_t i = 0; i < num_pos_xyz; ++i) {
    x.real_value += atoms->pos_x(i) * atoms->pos_x(i);
  }
}


void colvar::inertia::calc_gradients()
{
  const size_t num_pos_xyz = 3 * atoms->size();
  for (size_t i = 0; i < num_pos_xyz; ++i) {
    atoms->grad_x(i) = 2.0 * atoms->pos_x(i);
  }
}



colvar::inertia_z::inertia_z()
{
  set_function_type("inertiaZ");
}


int colvar::inertia_z::init(std::string const &conf)
{
  int error_code = inertia::init(conf);
  if (get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0))) {
    if (axis.norm2() == 0.0) {
      error_code |= cvmodule->error("Axis vector is zero!", COLVARS_INPUT_ERROR);
    }
    if (axis.norm2() != 1.0) {
      axis = axis.unit();
      cvmodule->log("The normalized axis is: "+cvm::to_str(axis)+".\n");
    }
  }
  return error_code;
}


void colvar::inertia_z::calc_value()
{
  x.real_value = 0.0;
  for (size_t i = 0; i < atoms->size(); ++i) {
    const cvm::atom_pos pos(atoms->pos_x(i),
                            atoms->pos_y(i),
                            atoms->pos_z(i));
    cvm::real const iprod = pos * axis;
    x.real_value += iprod * iprod;
  }
}


void colvar::inertia_z::calc_gradients()
{
  for (size_t i = 0; i < atoms->size(); ++i) {
    const cvm::atom_pos pos(atoms->pos_x(i),
                            atoms->pos_y(i),
                            atoms->pos_z(i));
    cvm::real const iprod = pos * axis;
    const cvm::rvector grad = 2.0 * iprod * axis;
    atoms->grad_x(i) = grad.x;
    atoms->grad_y(i) = grad.y;
    atoms->grad_z(i) = grad.z;
  }
}



colvar::rmsd::rmsd()
{
  set_function_type("rmsd");
  init_as_distance();
  provide(f_cvc_inv_gradient);
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  d_ref_pos_soa = nullptr;
  d_permutation_msds = nullptr;
  d_tbcounts = nullptr;
  h_rmsd = nullptr;
  h_best_perm_index = nullptr;
  d_ft = nullptr;
  h_ft = nullptr;
  d_jd = nullptr;
  h_jd = nullptr;
  d_tbcount_ft = nullptr;
  d_tbcount_jd = nullptr;
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}

colvar::rmsd::~rmsd() {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  colvarproxy* p = cvmodule->proxy;
  if (colvar::rmsd::has_gpu_implementation()) {
    p->deallocate_device(&d_ref_pos_soa);
    p->deallocate_device(&d_permutation_msds);
    p->deallocate_device(&d_tbcounts);
    p->deallocate_host(&h_rmsd);
    p->deallocate_host(&h_best_perm_index);
    p->deallocate_device(&d_ft);
    p->deallocate_device(&d_jd);
    p->deallocate_host(&h_ft);
    p->deallocate_host(&h_jd);
    p->deallocate_device(&d_tbcount_ft);
    p->deallocate_device(&d_tbcount_jd);
  }
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}

bool colvar::rmsd::has_gpu_implementation() const {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  const colvarproxy* p = cvmodule->proxy;
  if (p->get_smp_mode() == colvarproxy_smp::smp_mode_t::gpu){
    return true;
  } else {
    return false;
  }
#else
  return false;
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
}


int colvar::rmsd::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  atoms = parse_group(conf, "atoms");

  if (!atoms || atoms->size() == 0) {
    return error_code | COLVARS_INPUT_ERROR;
  }

  bool b_Jacobian_derivative = true;
  if (atoms->fitting_group != NULL && b_Jacobian_derivative) {
    cvmodule->log("The option \"fittingGroup\" (alternative group for fitting) was enabled: "
              "Jacobian derivatives of the RMSD will not be calculated.\n");
    b_Jacobian_derivative = false;
  }
  if (b_Jacobian_derivative) provide(f_cvc_Jacobian);

  std::vector<cvm::atom_pos> ref_pos;
  // the following is a simplified version of the corresponding atom group options;
  // we need this because the reference coordinates defined inside the atom group
  // may be used only for fitting, and even more so if fitting_group is used
  if (get_keyval(conf, "refPositions", ref_pos, ref_pos)) {
    cvmodule->log("Using reference positions from configuration file to calculate the variable.\n");
    if (ref_pos.size() != atoms->size()) {
      error_code |= cvmodule->error("Error: the number of reference positions provided (" +
                                   cvm::to_str(ref_pos.size()) +
                                   ") does not match the number of atoms of group \"atoms\" (" +
                                   cvm::to_str(atoms->size()) + ").\n",
                               COLVARS_INPUT_ERROR);
    }
  } else { // Only look for ref pos file if ref positions not already provided
    std::string ref_pos_file;
    if (get_keyval(conf, "refPositionsFile", ref_pos_file, std::string(""))) {

      if (ref_pos.size()) {
        error_code |= cvmodule->error("Error: cannot specify \"refPositionsFile\" and "
                                 "\"refPositions\" at the same time.\n",
                                 COLVARS_INPUT_ERROR);
      }

      std::string ref_pos_col;
      double ref_pos_col_value=0.0;

      if (get_keyval(conf, "refPositionsCol", ref_pos_col, std::string(""))) {
        // if provided, use PDB column to select coordinates
        bool found = get_keyval(conf, "refPositionsColValue", ref_pos_col_value, 0.0);
        if (found && ref_pos_col_value==0.0) {
          error_code |= cvmodule->error("Error: refPositionsColValue, "
                                   "if provided, must be non-zero.\n",
                                   COLVARS_INPUT_ERROR);
        }
      }

      ref_pos.resize(atoms->size());

      cvmodule->load_coords(ref_pos_file.c_str(), &ref_pos, atoms,
                       ref_pos_col, ref_pos_col_value);
    } else {
      error_code |= cvmodule->error(
          "Error: no reference positions for RMSD; use either refPositions of refPositionsFile.",
          COLVARS_INPUT_ERROR);
    }
  }

  if (ref_pos.size() != atoms->size()) {
    error_code |=
        cvmodule->error("Error: found " + cvm::to_str(ref_pos.size()) +
                       " reference positions for RMSD; expected " + cvm::to_str(atoms->size()),
                   COLVARS_INPUT_ERROR);
  }

  if (atoms->b_user_defined_fit) {
    cvmodule->log("WARNING: explicit fitting parameters were provided for atom group \"atoms\".\n");
  } else {
    // Default: fit everything
    cvmodule->log("Enabling \"centerToReference\" and \"rotateToReference\", to minimize RMSD before calculating it as a variable: "
              "if this is not the desired behavior, disable them explicitly within the \"atoms\" block.\n");
    atoms->enable(f_ag_center);
    atoms->enable(f_ag_rotate);
    // default case: reference positions for calculating the rmsd are also those used
    // for fitting
    atoms->set_ref_pos_from_aos(ref_pos);
    atoms->center_ref_pos();

    cvmodule->log("This is a standard minimum RMSD, derivatives of the optimal rotation "
              "will not be computed as they cancel out in the gradients.");
    atoms->disable(f_ag_fit_gradients);
  }
  atoms->setup_rotation_derivative();

  error_code |= init_permutation(ref_pos, conf);

  num_ref_pos = ref_pos.size();
  ref_pos_soa = cvm::atom_group::pos_aos_to_soa(ref_pos);
  if (has_gpu_implementation()) {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
    colvarproxy* p = cvmodule->proxy;
    error_code |= p->reallocate_device(&d_ref_pos_soa, 3 * num_ref_pos);
    error_code |= p->copy_HtoD(ref_pos_soa.data(), d_ref_pos_soa, 3 * num_ref_pos);
    error_code |= p->reallocate_device(&d_permutation_msds, n_permutations);
    error_code |= p->reallocate_device(&d_tbcounts, n_permutations);
    error_code |= p->reallocate_host(&h_rmsd, 1);
    error_code |= p->reallocate_host(&h_best_perm_index, 1);
    error_code |= p->clear_device_array(d_permutation_msds, n_permutations);
    error_code |= p->clear_device_array(d_tbcounts, n_permutations);
    error_code |= p->reallocate_device(&d_ft, 1);
    error_code |= p->clear_device_array(d_ft, 1);
    error_code |= p->reallocate_device(&d_jd, 1);
    error_code |= p->clear_device_array(d_jd, 1);
    error_code |= p->reallocate_host(&h_ft, 1);
    error_code |= p->reallocate_host(&h_jd, 1);
    error_code |= p->allocate_device(&d_tbcount_ft, 1);
    error_code |= p->allocate_device(&d_tbcount_jd, 1);
    error_code |= p->clear_device_array(d_tbcount_ft, 1);
    error_code |= p->clear_device_array(d_tbcount_jd, 1);
    disable(f_cvc_require_cpu_buffers);
#endif // (COLVARS_CUDA) || defined (COLVARS_GPU)
  }
  return error_code;
}


int colvar::rmsd::init_permutation(std::vector<cvm::atom_pos>& ref_pos, std::string const &conf)
{
  int error_code = COLVARS_OK;
  std::string perm_conf;
  size_t pos = 0; // current position in config string
  n_permutations = 1;

  while (key_lookup(conf, "atomPermutation", &perm_conf, &pos)) {
    cvmodule->cite_feature("Symmetry-adapted RMSD");
    std::vector<size_t> perm;
    if (perm_conf.size()) {
      std::istringstream is(perm_conf);
      size_t index;
      while (is >> index) {
        std::vector<int> const &ids = atoms->ids();
        size_t const ia = std::find(ids.begin(), ids.end(), index-1) - ids.begin();
        if (ia == atoms->size()) {
          error_code |= cvmodule->error("Error: atom id " + cvm::to_str(index) +
                                       " is not a member of group \"atoms\".",
                                   COLVARS_INPUT_ERROR);
        }
        if (std::find(perm.begin(), perm.end(), ia) != perm.end()) {
          error_code |= cvmodule->error("Error: atom id " + cvm::to_str(index) +
                                       " is mentioned more than once in atomPermutation list.",
                                   COLVARS_INPUT_ERROR);
        }
        perm.push_back(ia);
      }
      if (perm.size() != atoms->size()) {
        error_code |= cvmodule->error(
            "Error: symmetry permutation in input contains " + cvm::to_str(perm.size()) +
                " indices, but group \"atoms\" contains " + cvm::to_str(atoms->size()) + " atoms.",
            COLVARS_INPUT_ERROR);
      }
      cvmodule->log("atomPermutation = " + cvm::to_str(perm));
      n_permutations++;
      // Record a copy of reference positions in new order
      for (size_t ia = 0; ia < atoms->size(); ia++) {
        ref_pos.push_back(ref_pos[perm[ia]]);
      }
    }
  }

  permutation_msds.assign(n_permutations, 0);

  return error_code;
}


void colvar::rmsd::calc_value()
{
  // rotational-translational fit is handled by the atom group
  permutation_msds.assign(n_permutations, 0);
  // x.real_value = 0.0;
  for (size_t ip = 0; ip < n_permutations; ip++) {
    const auto ref_pos_x = ref_pos_soa.begin() + ip * atoms->size();
    const auto ref_pos_y = ref_pos_x + num_ref_pos;
    const auto ref_pos_z = ref_pos_y + num_ref_pos;
    for (size_t ia = 0; ia < atoms->size(); ia++) {
      const cvm::real diff_x = atoms->pos_x(ia) - ref_pos_x[ia];
      const cvm::real diff_y = atoms->pos_y(ia) - ref_pos_y[ia];
      const cvm::real diff_z = atoms->pos_z(ia) - ref_pos_z[ia];
      permutation_msds[ip] += diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
    }
  }

  best_perm_index = std::min_element(
    permutation_msds.begin(), permutation_msds.end()) - permutation_msds.begin();

  x.real_value = permutation_msds[best_perm_index];
  x.real_value /= cvm::real(atoms->size()); // MSD
  x.real_value = cvm::sqrt(x.real_value);
}


void colvar::rmsd::calc_gradients()
{
  cvm::real const drmsddx2 = (x.real_value > 0.0) ?
    0.5 / (x.real_value * cvm::real(atoms->size())) :
    0.0;

  // Use the appropriate symmetry permutation of reference positions to calculate gradients
  // size_t const start = atoms->size() * best_perm_index;
  const auto ref_pos_x = ref_pos_soa.begin() + best_perm_index * atoms->size();
  const auto ref_pos_y = ref_pos_x + num_ref_pos;
  const auto ref_pos_z = ref_pos_y + num_ref_pos;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    const cvm::atom_pos pos_ia(
      atoms->pos_x(ia), atoms->pos_y(ia), atoms->pos_z(ia));
    const cvm::atom_pos ref_pos_ia(
      ref_pos_x[ia], ref_pos_y[ia], ref_pos_z[ia]);
    const cvm::rvector grad = (drmsddx2 * 2.0 * (pos_ia - ref_pos_ia));
    atoms->grad_x(ia) = grad.x;
    atoms->grad_y(ia) = grad.y;
    atoms->grad_z(ia) = grad.z;
  }
}


void colvar::rmsd::calc_force_invgrads()
{
  atoms->read_total_forces();
  ft.real_value = 0.0;

  // Note: gradient square norm is 1/N_atoms
  const size_t num_grad_xyz = 3 * atoms->size();
  for (size_t i = 0; i < num_grad_xyz; i++) {
    ft.real_value += atoms->grad_x(i) * atoms->total_force_x(i);
  }
  ft.real_value *= atoms->size();
}


void colvar::rmsd::calc_Jacobian_derivative()
{
  // divergence of the rotated coordinates (including only derivatives of the rotation matrix)
  cvm::real rotation_term = 0.0;

  // The rotation term only applies is coordinates are rotated
  if (atoms->is_enabled(f_ag_rotate)) {

    // gradient of the rotation matrix
    cvm::matrix2d<cvm::rvector> grad_rot_mat(3, 3);
    // gradients of products of 2 quaternion components
    cvm::rvector g11, g22, g33, g01, g02, g03, g12, g13, g23;
    std::array<cvm::rvector, 4> dq;
    atoms->rot_deriv->prepare_derivative(rotation_derivative_dldq::use_dq);
    const auto ref_pos_x = ref_pos_soa.begin() + best_perm_index * atoms->size();
    const auto ref_pos_y = ref_pos_x + num_ref_pos;
    const auto ref_pos_z = ref_pos_y + num_ref_pos;
    for (size_t ia = 0; ia < atoms->size(); ia++) {

      // Gradient of optimal quaternion wrt current Cartesian position
      atoms->rot_deriv->calc_derivative_wrt_group1<false, true, false>(ia, nullptr, &dq);

      g11 = 2.0 * (atoms->rot.q)[1]*dq[1];
      g22 = 2.0 * (atoms->rot.q)[2]*dq[2];
      g33 = 2.0 * (atoms->rot.q)[3]*dq[3];
      g01 = (atoms->rot.q)[0]*dq[1] + (atoms->rot.q)[1]*dq[0];
      g02 = (atoms->rot.q)[0]*dq[2] + (atoms->rot.q)[2]*dq[0];
      g03 = (atoms->rot.q)[0]*dq[3] + (atoms->rot.q)[3]*dq[0];
      g12 = (atoms->rot.q)[1]*dq[2] + (atoms->rot.q)[2]*dq[1];
      g13 = (atoms->rot.q)[1]*dq[3] + (atoms->rot.q)[3]*dq[1];
      g23 = (atoms->rot.q)[2]*dq[3] + (atoms->rot.q)[3]*dq[2];

      // Gradient of the rotation matrix wrt current Cartesian position
      grad_rot_mat[0][0] = -2.0 * (g22 + g33);
      grad_rot_mat[1][0] =  2.0 * (g12 + g03);
      grad_rot_mat[2][0] =  2.0 * (g13 - g02);
      grad_rot_mat[0][1] =  2.0 * (g12 - g03);
      grad_rot_mat[1][1] = -2.0 * (g11 + g33);
      grad_rot_mat[2][1] =  2.0 * (g01 + g23);
      grad_rot_mat[0][2] =  2.0 * (g02 + g13);
      grad_rot_mat[1][2] =  2.0 * (g23 - g01);
      grad_rot_mat[2][2] = -2.0 * (g11 + g22);

      const cvm::atom_pos y(ref_pos_x[ia], ref_pos_y[ia], ref_pos_z[ia]);

      for (size_t alpha = 0; alpha < 3; alpha++) {
        for (size_t beta = 0; beta < 3; beta++) {
          rotation_term += grad_rot_mat[beta][alpha][alpha] * y[beta];
        // Note: equation was derived for inverse rotation (see colvars paper)
        // so here the matrix is transposed
        // (eq would give   divergence += grad_rot_mat[alpha][beta][alpha] * y[beta];)
        }
      }
    }
  }

  // The translation term only applies is coordinates are centered
  cvm::real translation_term = atoms->is_enabled(f_ag_center) ? 3.0 : 0.0;

  jd.real_value = x.real_value > 0.0 ?
    (3.0 * atoms->size() - 1.0 - translation_term - rotation_term) / x.real_value :
    0.0;
}

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
int colvar::rmsd::add_calc_value_node(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  std::vector<cudaGraphNode_t> calc_value_rmsd(n_permutations, nullptr);
  auto& gpu_buffers = atoms->get_gpu_atom_group()->get_gpu_buffers();
  error_code |= colvars_gpu::calc_value_rmsd(
    d_ref_pos_soa, gpu_buffers.d_atoms_pos,
    d_permutation_msds, permutation_msds.data(),
    atoms->size(), n_permutations, num_ref_pos,
    d_tbcounts, calc_value_rmsd, graph, {});
  for (int i = 0; i < n_permutations; ++i) {
    const std::string node_name = name + "_calc_value_permute_" + cvm::to_str(i);
    nodes_map[node_name] = calc_value_rmsd[i];
  }
  return error_code;
}

int colvar::rmsd::calc_value_after_gpu() {
  best_perm_index = std::min_element(
    permutation_msds.begin(), permutation_msds.end()) - permutation_msds.begin();
  x.real_value = permutation_msds[best_perm_index];
  x.real_value /= cvm::real(atoms->size()); // MSD
  x.real_value = cvm::sqrt(x.real_value);
  // Save the results to host-pinned memory for further gradients calculation
  *h_best_perm_index = best_perm_index;
  *h_rmsd = x.real_value;
  return COLVARS_OK;
}

int colvar::rmsd::add_calc_gradients_node(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  auto& gpu_buffers = atoms->get_gpu_atom_group()->get_gpu_buffers();
  cudaGraphNode_t calc_gradients_rmsd;
  error_code |= colvars_gpu::calc_gradients_rmsd(
    h_rmsd, h_best_perm_index, d_ref_pos_soa,
    gpu_buffers.d_atoms_pos, gpu_buffers.d_atoms_grad,
    atoms->size(), num_ref_pos, calc_gradients_rmsd, graph, {});
  const std::string node_name = name + "_calc_gradients";
  nodes_map[node_name] = calc_gradients_rmsd;
  return error_code;
}

int colvar::rmsd::add_calc_force_invgrads_node(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  auto& gpu_buffers = atoms->get_gpu_atom_group()->get_gpu_buffers();
  // If I explicitly read the total force to gpu_buffer.d_atoms_total_force
  // there would be a race condition if the underlyding atom group is used
  // by multiple other CVCs, so I have to fuse the reading of total force
  // with the Jacobian derivative calculation, and left
  // gpu_buffer.d_atoms_total_force unchanged.
  cudaGraphNode_t calc_force_invgrads;
  colvarproxy *p = cvmodule->proxy;
  const auto& rot = atoms->get_gpu_atom_group()->get_rot_gpu();
  error_code |= colvars_gpu::calc_force_invgrads_rmsd(
    atoms->is_enabled(f_ag_rotate),
    gpu_buffers.d_atoms_index,
    p->proxy_atoms_total_forces_gpu(),
    rot.get_q(), gpu_buffers.d_atoms_grad, d_ft, h_ft,
    atoms->size(), p->get_atom_ids()->size(), d_tbcount_ft,
    calc_force_invgrads, graph, {});
  const std::string node_name = name + "_calc_force_invgrads";
  nodes_map[node_name] = calc_force_invgrads;
  return error_code;
}

int colvar::rmsd::calc_force_invgrads_after_gpu() {
  ft.real_value = *h_ft;
  return COLVARS_OK;
}

int colvar::rmsd::add_calc_Jacobian_derivative_node(
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  const auto& rot = atoms->get_gpu_atom_group()->get_rot_gpu();
  auto* rot_deriv = atoms->get_gpu_atom_group()->get_rot_deriv_gpu();
  if (atoms->is_enabled(f_ag_rotate)) {
    if (rot_deriv == nullptr) {
      return cvmodule->error("BUG: rot_deriv is null.\n");
    }
    error_code |= rot_deriv->add_prepare_derivative_nodes(
      rotation_derivative_dldq::use_dq, graph, nodes_map);
  }
  cudaGraphNode_t calc_Jacobian_derivative;
  std::vector<cudaGraphNode_t> dependencies;
  if (atoms->is_enabled(f_ag_rotate)) {
    error_code |= colvars_gpu::prepare_dependencies(
      {{"prepare_rotation_derivative", true}}, dependencies, nodes_map);
  }
  error_code |= colvars_gpu::calc_Jacobian_derivative_rmsd(
    atoms->is_enabled(f_ag_center), atoms->is_enabled(f_ag_rotate),
    h_best_perm_index, d_ref_pos_soa, rot.get_q(),
    rot_deriv, h_rmsd, d_jd, h_jd, atoms->size(),
    num_ref_pos, d_tbcount_jd, calc_Jacobian_derivative,
    graph, dependencies);
  return error_code;
}

int colvar::rmsd::calc_Jacobian_derivative_after_gpu() {
  jd.real_value = *h_jd;
  return COLVARS_OK;
}
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)


colvar::eigenvector::eigenvector()
{
  set_function_type("eigenvector");
  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  x.type(colvarvalue::type_scalar);
}


int colvar::eigenvector::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  atoms = parse_group(conf, "atoms");
  if (!atoms || atoms->size() == 0) {
    return error_code | COLVARS_INPUT_ERROR;
  }

  {
    bool const b_inline = get_keyval(conf, "refPositions", ref_pos, ref_pos);

    if (b_inline) {
      cvmodule->log("Using reference positions from input file.\n");
      if (ref_pos.size() != atoms->size()) {
        error_code |= cvmodule->error("Error: reference positions do not "
                                 "match the number of requested atoms.\n",
                                 COLVARS_INPUT_ERROR);
      }
    }

    std::string file_name;
    if (get_keyval(conf, "refPositionsFile", file_name)) {

      if (b_inline) {
        error_code |= cvmodule->error(
            "Error: refPositions and refPositionsFile cannot be specified at the same time.\n",
            COLVARS_INPUT_ERROR);
      }
      std::string file_col;
      double file_col_value=0.0;
      if (get_keyval(conf, "refPositionsCol", file_col, std::string(""))) {
        // use PDB flags if column is provided
        bool found = get_keyval(conf, "refPositionsColValue", file_col_value, 0.0);
        if (found && file_col_value == 0.0) {
          error_code |= cvmodule->error("Error: refPositionsColValue, if provided, must be non-zero.\n",
                                   COLVARS_INPUT_ERROR);
        }
      }

      ref_pos.resize(atoms->size());
      cvmodule->load_coords(file_name.c_str(), &ref_pos, atoms,
                       file_col, file_col_value);
    }
  }

  if (ref_pos.size() == 0) {
    error_code |=
        cvmodule->error("Error: reference positions were not provided.\n", COLVARS_INPUT_ERROR);
  }

  if (ref_pos.size() != atoms->size()) {
    error_code |= cvmodule->error("Error: reference positions do not "
                             "match the number of requested atoms.\n",
                             COLVARS_INPUT_ERROR);
  }

  // save for later the geometric center of the provided positions (may not be the origin)
  cvm::rvector ref_pos_center(0.0, 0.0, 0.0);
  for (size_t i = 0; i < atoms->size(); i++) {
    ref_pos_center += ref_pos[i];
  }
  ref_pos_center *= 1.0 / atoms->size();

  if (atoms->b_user_defined_fit) {
    cvmodule->log("WARNING: explicit fitting parameters were provided for atom group \"atoms\".\n");
  } else {
    // default: fit everything
    cvmodule->log("Enabling \"centerToReference\" and \"rotateToReference\", to minimize RMSD before calculating the vector projection: "
              "if this is not the desired behavior, disable them explicitly within the \"atoms\" block.\n");
    atoms->enable(f_ag_center);
    atoms->enable(f_ag_rotate);
    atoms->enable(f_ag_fit_gradients);
    atoms->set_ref_pos_from_aos(ref_pos);
    atoms->center_ref_pos();
  }
  atoms->setup_rotation_derivative();

  {
    bool const b_inline = get_keyval(conf, "vector", eigenvec, eigenvec);
    // now load the eigenvector
    if (b_inline) {
      cvmodule->log("Using vector components from input file.\n");
      if (eigenvec.size() != atoms->size()) {
        error_code |= cvmodule->error("Error: vector components do not "
                                 "match the number of requested atoms->\n", COLVARS_INPUT_ERROR);
      }
    }

    std::string file_name;
    if (get_keyval(conf, "vectorFile", file_name)) {

      if (b_inline) {
        error_code |=
            cvmodule->error("Error: vector and vectorFile cannot be specified at the same time.\n",
                       COLVARS_INPUT_ERROR);
      }

      std::string file_col;
      double file_col_value=0.0;
      if (get_keyval(conf, "vectorCol", file_col, std::string(""))) {
        // use PDB flags if column is provided
        bool found = get_keyval(conf, "vectorColValue", file_col_value, 0.0);
        if (found && file_col_value==0.0) {
          error_code |= cvmodule->error("Error: vectorColValue, if provided, must be non-zero.\n",
                                   COLVARS_INPUT_ERROR);
        }
      }

      eigenvec.resize(atoms->size());
      cvmodule->load_coords(file_name.c_str(), &eigenvec, atoms,
                       file_col, file_col_value);
    }
  }

  if (!ref_pos.size() || !eigenvec.size()) {
    error_code |= cvmodule->error("Error: both reference coordinates and eigenvector must be defined.\n",
                             COLVARS_INPUT_ERROR);
  }

  cvm::atom_pos eig_center(0.0, 0.0, 0.0);
  for (size_t eil = 0; eil < atoms->size(); eil++) {
    eig_center += eigenvec[eil];
  }
  eig_center *= 1.0 / atoms->size();
  cvmodule->log("Geometric center of the provided vector: "+cvm::to_str(eig_center)+"\n");

  bool b_difference_vector = false;
  get_keyval(conf, "differenceVector", b_difference_vector, false);

  if (b_difference_vector) {

    if (atoms->is_enabled(f_ag_center)) {
      // both sets should be centered on the origin for fitting
      for (size_t i = 0; i < atoms->size(); i++) {
        eigenvec[i] -= eig_center;
        ref_pos[i]  -= ref_pos_center;
      }
    }
    if (atoms->is_enabled(f_ag_rotate)) {
      atoms->rot.calc_optimal_rotation(eigenvec, ref_pos);
      const auto rot_mat = atoms->rot.matrix();
      for (size_t i = 0; i < atoms->size(); i++) {
        eigenvec[i] = rot_mat * eigenvec[i];
      }
    }
    cvmodule->log("\"differenceVector\" is on: subtracting the reference positions from the provided vector: v = x_vec - x_ref.\n");
    for (size_t i = 0; i < atoms->size(); i++) {
      eigenvec[i] -= ref_pos[i];
    }
    if (atoms->is_enabled(f_ag_center)) {
      // bring back the ref positions to where they were
      for (size_t i = 0; i < atoms->size(); i++) {
        ref_pos[i] += ref_pos_center;
      }
    }

  } else {
    cvmodule->log("Centering the provided vector to zero.\n");
    for (size_t i = 0; i < atoms->size(); i++) {
      eigenvec[i] -= eig_center;
    }
  }

  // eigenvec_invnorm2 is used when computing inverse gradients
  eigenvec_invnorm2 = 0.0;
  for (size_t ein = 0; ein < atoms->size(); ein++) {
    eigenvec_invnorm2 += eigenvec[ein].norm2();
  }
  eigenvec_invnorm2 = 1.0 / eigenvec_invnorm2;

  // Vector normalization overrides the default normalization for differenceVector
  bool normalize = false;
  get_keyval(conf, "normalizeVector", normalize, normalize);

  if (normalize) {
    cvmodule->log("Normalizing the vector so that |v| = 1.\n");
    for (size_t i = 0; i < atoms->size(); i++) {
      eigenvec[i] *= cvm::sqrt(eigenvec_invnorm2);
    }
    eigenvec_invnorm2 = 1.0;
  } else if (b_difference_vector) {
    cvmodule->log("Normalizing the vector so that the norm of the projection |v ⋅ (x_vec - x_ref)| = 1.\n");
    for (size_t i = 0; i < atoms->size(); i++) {
      eigenvec[i] *= eigenvec_invnorm2;
    }
    eigenvec_invnorm2 = 1.0/eigenvec_invnorm2;
  } else {
    cvmodule->log("The norm of the vector is |v| = "+
             cvm::to_str(1.0/cvm::sqrt(eigenvec_invnorm2))+".\n");
  }

  return error_code;
}


void colvar::eigenvector::calc_value()
{
  x.real_value = 0.0;
  for (size_t i = 0; i < atoms->size(); i++) {
    const cvm::atom_pos pos_i(
      atoms->pos_x(i), atoms->pos_y(i), atoms->pos_z(i));
    x.real_value += (pos_i - ref_pos[i]) * eigenvec[i];
  }
}


void colvar::eigenvector::calc_gradients()
{
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    atoms->grad_x(ia) = eigenvec[ia].x;
    atoms->grad_y(ia) = eigenvec[ia].y;
    atoms->grad_z(ia) = eigenvec[ia].z;
  }
}


void colvar::eigenvector::calc_force_invgrads()
{
  atoms->read_total_forces();
  ft.real_value = 0.0;
  const size_t num_grad_xyz = 3 * atoms->size();
  for (size_t i = 0; i < num_grad_xyz; i++) {
    ft.real_value += atoms->grad_x(i) * atoms->total_force_x(i);
  }
  ft.real_value *= eigenvec_invnorm2;
}


void colvar::eigenvector::calc_Jacobian_derivative()
{
  // gradient of the rotation matrix
  cvm::matrix2d<cvm::rvector> grad_rot_mat(3, 3);
  cvm::quaternion &quat0 = atoms->rot.q;

  // gradients of products of 2 quaternion components
  cvm::rvector g11, g22, g33, g01, g02, g03, g12, g13, g23;

  cvm::real sum = 0.0;

  std::array<cvm::rvector, 4> dq_1;
  atoms->rot_deriv->prepare_derivative(rotation_derivative_dldq::use_dq);
  for (size_t ia = 0; ia < atoms->size(); ia++) {

    // Gradient of optimal quaternion wrt current Cartesian position
    // trick: d(R^-1)/dx = d(R^t)/dx = (dR/dx)^t
    // we can just transpose the derivatives of the direct matrix
    atoms->rot_deriv->calc_derivative_wrt_group1<false, true, false>(ia, nullptr, &dq_1);

    g11 = 2.0 * quat0[1]*dq_1[1];
    g22 = 2.0 * quat0[2]*dq_1[2];
    g33 = 2.0 * quat0[3]*dq_1[3];
    g01 = quat0[0]*dq_1[1] + quat0[1]*dq_1[0];
    g02 = quat0[0]*dq_1[2] + quat0[2]*dq_1[0];
    g03 = quat0[0]*dq_1[3] + quat0[3]*dq_1[0];
    g12 = quat0[1]*dq_1[2] + quat0[2]*dq_1[1];
    g13 = quat0[1]*dq_1[3] + quat0[3]*dq_1[1];
    g23 = quat0[2]*dq_1[3] + quat0[3]*dq_1[2];

    // Gradient of the inverse rotation matrix wrt current Cartesian position
    // (transpose of the gradient of the direct rotation)
    grad_rot_mat[0][0] = -2.0 * (g22 + g33);
    grad_rot_mat[0][1] =  2.0 * (g12 + g03);
    grad_rot_mat[0][2] =  2.0 * (g13 - g02);
    grad_rot_mat[1][0] =  2.0 * (g12 - g03);
    grad_rot_mat[1][1] = -2.0 * (g11 + g33);
    grad_rot_mat[1][2] =  2.0 * (g01 + g23);
    grad_rot_mat[2][0] =  2.0 * (g02 + g13);
    grad_rot_mat[2][1] =  2.0 * (g23 - g01);
    grad_rot_mat[2][2] = -2.0 * (g11 + g22);

    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        sum += grad_rot_mat[i][j][i] * eigenvec[ia][j];
      }
    }
  }

  jd.real_value = sum * cvm::sqrt(eigenvec_invnorm2);
}



colvar::cartesian::cartesian()
{
  set_function_type("cartesian");
  x.type(colvarvalue::type_vector);
  disable(f_cvc_explicit_gradient);
}


int colvar::cartesian::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  atoms = parse_group(conf, "atoms");

  bool use_x, use_y, use_z;
  get_keyval(conf, "useX", use_x, true);
  get_keyval(conf, "useY", use_y, true);
  get_keyval(conf, "useZ", use_z, true);

  axes.clear();
  if (use_x) axes.push_back(0);
  if (use_y) axes.push_back(1);
  if (use_z) axes.push_back(2);

  if (axes.size() == 0) {
    error_code |=
        cvmodule->error("Error: a \"cartesian\" component was defined with all three axes disabled.\n",
                   COLVARS_INPUT_ERROR);
  }

  // Don't try to access atoms if creation of the atom group failed
  if (atoms) x.vector1d_value.resize(atoms->size() * axes.size());

  return error_code;
}


void colvar::cartesian::calc_value()
{
  size_t const dim = axes.size();
  size_t ia, j;
  for (ia = 0; ia < atoms->size(); ia++) {
    for (j = 0; j < dim; j++) {
      const cvm::atom_pos pos_ia(
        atoms->pos_x(ia), atoms->pos_y(ia), atoms->pos_z(ia));
      x.vector1d_value[dim*ia + j] = pos_ia[axes[j]];
    }
  }
}


void colvar::cartesian::calc_gradients()
{
  // we're not using the "grad" member of each
  // atom object, because it only can represent the gradient of a
  // scalar colvar
}


void colvar::cartesian::apply_force(colvarvalue const &force)
{
  size_t const dim = axes.size();
  size_t ia, j;
  if (!atoms->noforce) {
    cvm::rvector f;
    auto ag_force = atoms->get_group_force_object();
    for (ia = 0; ia < atoms->size(); ia++) {
      for (j = 0; j < dim; j++) {
        f[axes[j]] = force.vector1d_value[dim*ia + j];
      }
      ag_force.add_atom_force(ia, f);
    }
  }
}


cvm::real colvar::cartesian::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  return (x1.vector1d_value - x2.vector1d_value).norm2();
}


colvarvalue colvar::cartesian::dist2_lgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return 2.0 * (x1.vector1d_value - x2.vector1d_value);
}


colvarvalue colvar::cartesian::dist2_rgrad(colvarvalue const &x1, colvarvalue const &x2) const
{
  return cartesian::dist2_lgrad(x1, x2);
}


void colvar::cartesian::wrap(colvarvalue & /* x_unwrapped */) const {}

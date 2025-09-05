// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvar_rotation_derivative.h"

struct colvar::orientation::rotation_derivative_impl_: public rotation_derivative {
public:
  rotation_derivative_impl_(colvar::orientation* orientation_cvc):
   rotation_derivative(
    orientation_cvc->rot,
    orientation_cvc->ref_pos_soa,
    orientation_cvc->shifted_pos_soa,
    orientation_cvc->num_ref_pos,
    orientation_cvc->num_shifted_pos) {}
};


colvar::orientation::orientation()
{
  set_function_type("orientation");
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_quaternion);
}


colvar::orientation::~orientation() {
  num_ref_pos = 0;
}


int colvar::orientation::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  atoms = parse_group(conf, "atoms");
  if (!atoms || atoms->size() == 0) {
    return error_code | COLVARS_INPUT_ERROR;
  }
  std::vector<cvm::atom_pos> ref_pos;
  ref_pos.reserve(atoms->size());

  if (get_keyval(conf, "refPositions", ref_pos, ref_pos)) {
    cvm::log("Using reference positions from input file.\n");
    if (ref_pos.size() != atoms->size()) {
      return cvm::error("Error: reference positions do not "
                        "match the number of requested atoms.\n", COLVARS_INPUT_ERROR);
    }
  }

  {
    std::string file_name;
    if (get_keyval(conf, "refPositionsFile", file_name)) {

      std::string file_col;
      double file_col_value=0.0;
      if (get_keyval(conf, "refPositionsCol", file_col, std::string(""))) {
        // use PDB flags if column is provided
        bool found = get_keyval(conf, "refPositionsColValue", file_col_value, 0.0);
        if (found && file_col_value==0.0) {
          return cvm::error("Error: refPositionsColValue, "
                            "if provided, must be non-zero.\n", COLVARS_INPUT_ERROR);
        }
      }

      ref_pos.resize(atoms->size());
      error_code |= cvm::load_coords(file_name.c_str(), &ref_pos, atoms,
                       file_col, file_col_value);
    }
  }

  if (error_code != COLVARS_OK) return error_code;

  if (!ref_pos.size()) {
    return cvm::error("Error: must define a set of "
                      "reference coordinates.\n", COLVARS_INPUT_ERROR);
  }

  cvm::rvector ref_cog(0.0, 0.0, 0.0);
  size_t i;
  for (i = 0; i < ref_pos.size(); i++) {
    ref_cog += ref_pos[i];
  }
  ref_cog /= cvm::real(ref_pos.size());
  cvm::log("Centering the reference coordinates on the origin by subtracting "
           "the center of geometry at "+
           cvm::to_str(-1.0 * ref_cog)+"; it is "
           "assumed that each atom is the closest "
           "periodic image to the center of geometry.\n");
  for (i = 0; i < ref_pos.size(); i++) {
    ref_pos[i] -= ref_cog;
  }

  get_keyval(conf, "closestToQuaternion", ref_quat, cvm::quaternion(1.0, 0.0, 0.0, 0.0));

  // If the debug gradients feature is active, debug the rotation gradients
  // (note that this won't be active for the orientation CVC itself, because
  // colvardeps prevents the flag's activation)
  rot.b_debug_gradients = is_enabled(f_cvc_debug_gradient);
  ref_pos_soa = cvm::atom_group::pos_aos_to_soa(ref_pos);
  num_ref_pos = ref_pos.size();
  num_shifted_pos = atoms->size();
  rot_deriv_impl = std::unique_ptr<rotation_derivative_impl_>(new rotation_derivative_impl_(this));
  return error_code;
}


void colvar::orientation::calc_value()
{
  atoms_cog = atoms->center_of_geometry();
  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());

  if ((rot.q).inner(ref_quat) >= 0.0) {
    x.quaternion_value = rot.q;
  } else {
    x.quaternion_value = -1.0 * rot.q;
  }
}


void colvar::orientation::calc_gradients()
{
  // gradients have already been calculated and stored within the
  // member object "rot"; we're not using the "grad" member of each
  // atom object, because it only can represent the gradient of a
  // scalar colvar
}


void colvar::orientation::apply_force(colvarvalue const &force)
{
  cvm::quaternion const &FQ = force.quaternion_value;

  if (!atoms->noforce) {
    const cvm::real sign = (rot.q).inner(ref_quat) >= 0.0 ? 1.0 : -1.0;
    rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
    cvm::vector1d<cvm::rvector> dq0_2;
    auto ag_force = atoms->get_group_force_object();
    for (size_t ia = 0; ia < atoms->size(); ia++) {
      rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
      const auto f_ia = sign * (FQ[0] * dq0_2[0] +
                        FQ[1] * dq0_2[1] +
                        FQ[2] * dq0_2[2] +
                        FQ[3] * dq0_2[3]);
      ag_force.add_atom_force(ia, f_ia);
    }
  }
}


cvm::real colvar::orientation::dist2(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return x1.quaternion_value.dist2(x2);
}


colvarvalue colvar::orientation::dist2_lgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return x1.quaternion_value.dist2_grad(x2);
}


colvarvalue colvar::orientation::dist2_rgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return x2.quaternion_value.dist2_grad(x1);
}


void colvar::orientation::wrap(colvarvalue & /* x_unwrapped */) const {}



colvar::orientation_angle::orientation_angle()
{
  set_function_type("orientationAngle");
  init_as_angle();
  enable(f_cvc_explicit_gradient);
}


void colvar::orientation_angle::calc_value()
{
  atoms_cog = atoms->center_of_geometry();
  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());

  if ((rot.q).q0 >= 0.0) {
    x.real_value = (180.0/PI) * 2.0 * cvm::acos((rot.q).q0);
  } else {
    x.real_value = (180.0/PI) * 2.0 * cvm::acos(-1.0 * (rot.q).q0);
  }
}


void colvar::orientation_angle::calc_gradients()
{
  const cvm::real sign = (rot.q).q0 >= 0.0 ? 1.0 : -1.0;
  cvm::real const dxdq0 =
    sign * ( ((rot.q).q0 * (rot.q).q0 < 1.0) ?
      ((180.0 / PI) * (-2.0) / cvm::sqrt(1.0 - ((rot.q).q0 * (rot.q).q0))) :
      0.0 );

  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
    const cvm::rvector g = dxdq0 * dq0_2[0];
    atoms->grad_x(ia) = g.x;
    atoms->grad_y(ia) = g.y;
    atoms->grad_z(ia) = g.z;
  }
}


void colvar::orientation_angle::apply_force(colvarvalue const &force)
{
  cvc::apply_force(force);
}


cvm::real colvar::orientation_angle::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  return cvc::dist2(x1, x2);
}


colvarvalue colvar::orientation_angle::dist2_lgrad(colvarvalue const &x1,
                                                   colvarvalue const &x2) const
{
  return cvc::dist2_lgrad(x1, x2);
}


colvarvalue colvar::orientation_angle::dist2_rgrad(colvarvalue const &x1,
                                                   colvarvalue const &x2) const
{
  return cvc::dist2_rgrad(x1, x2);
}


void colvar::orientation_angle::wrap(colvarvalue & /* x_unwrapped */) const {}



colvar::orientation_proj::orientation_proj()
{
  set_function_type("orientationProj");
  enable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);
}


void colvar::orientation_proj::calc_value()
{
  atoms_cog = atoms->center_of_geometry();
  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());
  x.real_value = 2.0 * (rot.q).q0 * (rot.q).q0 - 1.0;
}


void colvar::orientation_proj::calc_gradients()
{
  cvm::real const dxdq0 = 2.0 * 2.0 * (rot.q).q0;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
    const cvm::rvector g = dxdq0 * dq0_2[0];
    atoms->grad_x(ia) = g.x;
    atoms->grad_y(ia) = g.y;
    atoms->grad_z(ia) = g.z;
  }
}



colvar::tilt::tilt()
{
  set_function_type("tilt");
  x.type(colvarvalue::type_scalar);
  enable(f_cvc_explicit_gradient);
  init_scalar_boundaries(-1.0, 1.0);
}


int colvar::tilt::init(std::string const &conf)
{
  int error_code = orientation_proj::init(conf);

  get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0));
  if (axis.norm2() != 1.0) {
    axis /= axis.norm();
    cvm::log("Normalizing rotation axis to "+cvm::to_str(axis)+".\n");
  }

  return error_code;
}


void colvar::tilt::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());

  x.real_value = rot.cos_theta(axis);
}


void colvar::tilt::calc_gradients()
{
  cvm::quaternion const dxdq = rot.dcos_theta_dq(axis);

  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
    cvm::rvector grad(0, 0, 0);
    for (size_t iq = 0; iq < 4; iq++) {
      grad += (dxdq[iq] * dq0_2[iq]);
    }
    atoms->grad_x(ia) = grad.x;
    atoms->grad_y(ia) = grad.y;
    atoms->grad_z(ia) = grad.z;
  }
}



colvar::spin_angle::spin_angle()
{
  set_function_type("spinAngle");
  init_as_periodic_angle();
  enable(f_cvc_explicit_gradient);
}


void colvar::spin_angle::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());

  x.real_value = rot.spin_angle(axis);
  wrap(x);
}


void colvar::spin_angle::calc_gradients()
{
  cvm::quaternion const dxdq = rot.dspin_angle_dq(axis);

  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
    cvm::rvector grad(0, 0, 0);
    for (size_t iq = 0; iq < 4; iq++) {
      grad += (dxdq[iq] * dq0_2[iq]);
    }
    atoms->grad_x(ia) = grad.x;
    atoms->grad_y(ia) = grad.y;
    atoms->grad_z(ia) = grad.z;
  }
}



colvar::euler_phi::euler_phi()
{
  set_function_type("eulerPhi");
  init_as_periodic_angle();
  enable(f_cvc_explicit_gradient);
}


void colvar::euler_phi::calc_value()
{
  atoms_cog = atoms->center_of_geometry();
  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());

  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real tmp_y = 2 * (q0 * q1 + q2 * q3);
  const cvm::real tmp_x = 1 - 2 * (q1 * q1 + q2 * q2);
  x.real_value = cvm::atan2(tmp_y, tmp_x) * (180.0/PI);
}


void colvar::euler_phi::calc_gradients()
{
  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real denominator = (2 * q0 * q1 + 2 * q2 * q3) * (2 * q0 * q1 + 2 * q2 * q3) + (-2 * q1 * q1 - 2 * q2 * q2 + 1) * (-2 * q1 * q1 - 2 * q2 * q2 + 1);
  const cvm::real dxdq0 = (180.0/PI) * 2 * q1 * (-2 * q1 * q1 - 2 * q2 * q2 + 1) / denominator;
  const cvm::real dxdq1 = (180.0/PI) * (2 * q0 * (-2 * q1 * q1 - 2 * q2 * q2 + 1) - 4 * q1 * (-2 * q0 * q1 - 2 * q2 * q3)) / denominator;
  const cvm::real dxdq2 = (180.0/PI) * (-4 * q2 * (-2 * q0 * q1 - 2 * q2 * q3) + 2 * q3 * (-2 * q1 * q1 - 2 * q2 * q2 + 1)) / denominator;
  const cvm::real dxdq3 = (180.0/PI) * 2 * q2 * (-2 * q1 * q1 - 2 * q2 * q2 + 1) / denominator;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
    const cvm::rvector grad = (dxdq0 * dq0_2[0]) +
                              (dxdq1 * dq0_2[1]) +
                              (dxdq2 * dq0_2[2]) +
                              (dxdq3 * dq0_2[3]);
    atoms->grad_x(ia) = grad.x;
    atoms->grad_y(ia) = grad.y;
    atoms->grad_z(ia) = grad.z;
  }
}



colvar::euler_psi::euler_psi()
{
  set_function_type("eulerPsi");
  init_as_periodic_angle();
  enable(f_cvc_explicit_gradient);
}


void colvar::euler_psi::calc_value()
{
  atoms_cog = atoms->center_of_geometry();
  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());

  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real tmp_y = 2 * (q0 * q3 + q1 * q2);
  const cvm::real tmp_x = 1 - 2 * (q2 * q2 + q3 * q3);
  x.real_value = cvm::atan2(tmp_y, tmp_x) * (180.0/PI);
}


void colvar::euler_psi::calc_gradients()
{
  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real denominator = (2 * q0 * q3 + 2 * q1 * q2) * (2 * q0 * q3 + 2 * q1 * q2) + (-2 * q2 * q2 - 2 * q3 * q3 + 1) * (-2 * q2 * q2 - 2 * q3 * q3 + 1);
  const cvm::real dxdq0 = (180.0/PI) * 2 * q3 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) / denominator;
  const cvm::real dxdq1 = (180.0/PI) * 2 * q2 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) / denominator;
  const cvm::real dxdq2 = (180.0/PI) * (2 * q1 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) - 4 * q2 * (-2 * q0 * q3 - 2 * q1 * q2)) / denominator;
  const cvm::real dxdq3 = (180.0/PI) * (2 * q0 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) - 4 * q3 * (-2 * q0 * q3 - 2 * q1 * q2)) / denominator;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
    const cvm::rvector grad = (dxdq0 * dq0_2[0]) +
                              (dxdq1 * dq0_2[1]) +
                              (dxdq2 * dq0_2[2]) +
                              (dxdq3 * dq0_2[3]);
    atoms->grad_x(ia) = grad.x;
    atoms->grad_y(ia) = grad.y;
    atoms->grad_z(ia) = grad.z;
  }
}



colvar::euler_theta::euler_theta()
{
  set_function_type("eulerTheta");
  init_as_angle();
  enable(f_cvc_explicit_gradient);
}


void colvar::euler_theta::calc_value()
{
  atoms_cog = atoms->center_of_geometry();
  shifted_pos_soa = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation_soa(ref_pos_soa, shifted_pos_soa, num_ref_pos, atoms->size());

  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  x.real_value = cvm::asin(2 * (q0 * q2 - q3 * q1)) * (180.0/PI);
}


void colvar::euler_theta::calc_gradients()
{
  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real denominator = cvm::sqrt(1 - (2 * q0 * q2 - 2 * q1 * q3) * (2 * q0 * q2 - 2 * q1 * q3));
  const cvm::real dxdq0 = (180.0/PI) * 2 * q2 / denominator;
  const cvm::real dxdq1 = (180.0/PI) * -2 * q3 / denominator;
  const cvm::real dxdq2 = (180.0/PI) * 2 * q0 / denominator;
  const cvm::real dxdq3 = (180.0/PI) * -2 * q1 / denominator;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2<false, true, false>(ia, nullptr, &dq0_2);
    const cvm::rvector grad = (dxdq0 * dq0_2[0]) +
                              (dxdq1 * dq0_2[1]) +
                              (dxdq2 * dq0_2[2]) +
                              (dxdq3 * dq0_2[3]);
    atoms->grad_x(ia) = grad.x;
    atoms->grad_y(ia) = grad.y;
    atoms->grad_z(ia) = grad.z;
  }
}

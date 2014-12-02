/// -*- c++ -*-

#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"



colvar::orientation::orientation(std::string const &conf)
  : cvc(conf)
{
  function_type = "orientation";
  parse_group(conf, "atoms", atoms);
  atom_groups.push_back(&atoms);
  x.type(colvarvalue::type_quaternion);

  ref_pos.reserve(atoms.size());

  if (get_keyval(conf, "refPositions", ref_pos, ref_pos)) {
    cvm::log("Using reference positions from input file.\n");
    if (ref_pos.size() != atoms.size()) {
      cvm::fatal_error("Error: reference positions do not "
                        "match the number of requested atoms.\n");
    }
  }

  {
    std::string file_name;
    if (get_keyval(conf, "refPositionsFile", file_name)) {

      std::string file_col;
        double file_col_value;
      if (get_keyval(conf, "refPositionsCol", file_col, std::string(""))) {
        // use PDB flags if column is provided
        bool found = get_keyval(conf, "refPositionsColValue", file_col_value, 0.0);
        if (found && !file_col_value)
          cvm::fatal_error("Error: refPositionsColValue, "
                            "if provided, must be non-zero.\n");
      } else {
        // if not, use atom indices
        atoms.create_sorted_ids();
      }
      ref_pos.resize(atoms.size());
      cvm::load_coords(file_name.c_str(), ref_pos, atoms.sorted_ids, file_col, file_col_value);
    }
  }

  if (!ref_pos.size()) {
    cvm::fatal_error("Error: must define a set of "
                      "reference coordinates.\n");
  }


  cvm::log("Centering the reference coordinates: it is "
            "assumed that each atom is the closest "
            "periodic image to the center of geometry.\n");
  cvm::rvector cog(0.0, 0.0, 0.0);
  size_t i;
  for (i = 0; i < ref_pos.size(); i++) {
    cog += ref_pos[i];
  }
  cog /= cvm::real(ref_pos.size());
  for (i = 0; i < ref_pos.size(); i++) {
    ref_pos[i] -= cog;
  }

  get_keyval(conf, "closestToQuaternion", ref_quat, cvm::quaternion(1.0, 0.0, 0.0, 0.0));

  // initialize rot member data
  if (!atoms.noforce) {
    rot.request_group2_gradients(atoms.size());
  }

}


colvar::orientation::orientation()
  : cvc()
{
  function_type = "orientation";
  x.type(colvarvalue::type_quaternion);
}


void colvar::orientation::calc_value()
{
  // atoms.reset_atoms_data();
  // atoms.read_positions();

  atoms_cog = atoms.center_of_geometry();

  rot.calc_optimal_rotation(ref_pos, atoms.positions_shifted(-1.0 * atoms_cog));

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

  if (!atoms.noforce) {
    for (size_t ia = 0; ia < atoms.size(); ia++) {
      for (size_t i = 0; i < 4; i++) {
        atoms[ia].apply_force(FQ[i] * rot.dQ0_2[ia][i]);
      }
    }
  }
}



colvar::orientation_angle::orientation_angle(std::string const &conf)
  : orientation(conf)
{
  function_type = "orientation_angle";
  x.type(colvarvalue::type_scalar);
}


colvar::orientation_angle::orientation_angle()
  : orientation()
{
  function_type = "orientation_angle";
  x.type(colvarvalue::type_scalar);
}


void colvar::orientation_angle::calc_value()
{
  // atoms.reset_atoms_data();
  // atoms.read_positions();

  atoms_cog = atoms.center_of_geometry();

  rot.calc_optimal_rotation(ref_pos, atoms.positions_shifted(-1.0 * atoms_cog));

  if ((rot.q).q0 >= 0.0) {
    x.real_value = (180.0/PI) * 2.0 * std::acos((rot.q).q0);
  } else {
    x.real_value = (180.0/PI) * 2.0 * std::acos(-1.0 * (rot.q).q0);
  }
}


void colvar::orientation_angle::calc_gradients()
{
  cvm::real const dxdq0 =
    ( ((rot.q).q0 * (rot.q).q0 < 1.0) ?
      ((180.0 / PI) * (-2.0) / std::sqrt(1.0 - ((rot.q).q0 * (rot.q).q0))) :
      0.0 );

  for (size_t ia = 0; ia < atoms.size(); ia++) {
    atoms[ia].grad = (dxdq0 * (rot.dQ0_2[ia])[0]);
  }
  if (b_debug_gradients) {
    cvm::log("Debugging orientationAngle component gradients:\n");
    debug_gradients(atoms);
  }
}


void colvar::orientation_angle::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;
  if (!atoms.noforce) {
    atoms.apply_colvar_force(fw);
  }
}



colvar::orientation_proj::orientation_proj(std::string const &conf)
  : orientation(conf)
{
  function_type = "orientation_proj";
  x.type(colvarvalue::type_scalar);
}


colvar::orientation_proj::orientation_proj()
  : orientation()
{
  function_type = "orientation_proj";
  x.type(colvarvalue::type_scalar);
}


void colvar::orientation_proj::calc_value()
{
  atoms_cog = atoms.center_of_geometry();
  rot.calc_optimal_rotation(ref_pos, atoms.positions_shifted(-1.0 * atoms_cog));
  x.real_value = 2.0 * (rot.q).q0 * (rot.q).q0 - 1.0;
}


void colvar::orientation_proj::calc_gradients()
{
  cvm::real const dxdq0 = 2.0 * 2.0 * (rot.q).q0;
  for (size_t ia = 0; ia < atoms.size(); ia++) {
    atoms[ia].grad = (dxdq0 * (rot.dQ0_2[ia])[0]);
  }
  if (b_debug_gradients) {
    cvm::log("Debugging orientationProj component gradients:\n");
    debug_gradients(atoms);
  }
}


void colvar::orientation_proj::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;

  if (!atoms.noforce) {
    atoms.apply_colvar_force(fw);
  }
}



colvar::tilt::tilt(std::string const &conf)
  : orientation(conf)
{
  function_type = "tilt";

  get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0));

  if (axis.norm2() != 1.0) {
    axis /= axis.norm();
    cvm::log("Normalizing rotation axis to "+cvm::to_str(axis)+".\n");
  }

  x.type(colvarvalue::type_scalar);
}


colvar::tilt::tilt()
  : orientation()
{
  function_type = "tilt";
  x.type(colvarvalue::type_scalar);
}


void colvar::tilt::calc_value()
{
  // atoms.reset_atoms_data();
  // atoms.read_positions();

  atoms_cog = atoms.center_of_geometry();

  rot.calc_optimal_rotation(ref_pos, atoms.positions_shifted(-1.0 * atoms_cog));

  x.real_value = rot.cos_theta(axis);
}


void colvar::tilt::calc_gradients()
{
  cvm::quaternion const dxdq = rot.dcos_theta_dq(axis);

  for (size_t ia = 0; ia < atoms.size(); ia++) {
    atoms[ia].grad = cvm::rvector(0.0, 0.0, 0.0);
    for (size_t iq = 0; iq < 4; iq++) {
      atoms[ia].grad += (dxdq[iq] * (rot.dQ0_2[ia])[iq]);
    }
  }

  if (b_debug_gradients) {
    cvm::log("Debugging tilt component gradients:\n");
    debug_gradients(atoms);
  }
}


void colvar::tilt::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;

  if (!atoms.noforce) {
    atoms.apply_colvar_force(fw);
  }
}



colvar::spin_angle::spin_angle(std::string const &conf)
  : orientation(conf)
{
  function_type = "spin_angle";

  get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0));

  if (axis.norm2() != 1.0) {
    axis /= axis.norm();
    cvm::log("Normalizing rotation axis to "+cvm::to_str(axis)+".\n");
  }

  period = 360.0;
  b_periodic = true;
  x.type(colvarvalue::type_scalar);
}


colvar::spin_angle::spin_angle()
  : orientation()
{
  function_type = "spin_angle";
  period = 360.0;
  b_periodic = true;
  x.type(colvarvalue::type_scalar);
}


void colvar::spin_angle::calc_value()
{
  // atoms.reset_atoms_data();
  // atoms.read_positions();

  atoms_cog = atoms.center_of_geometry();

  rot.calc_optimal_rotation(ref_pos, atoms.positions_shifted(-1.0 * atoms_cog));

  x.real_value = rot.spin_angle(axis);
  this->wrap(x);
}


void colvar::spin_angle::calc_gradients()
{
  cvm::quaternion const dxdq = rot.dspin_angle_dq(axis);

  for (size_t ia = 0; ia < atoms.size(); ia++) {
    atoms[ia].grad = cvm::rvector(0.0, 0.0, 0.0);
    for (size_t iq = 0; iq < 4; iq++) {
      atoms[ia].grad += (dxdq[iq] * (rot.dQ0_2[ia])[iq]);
    }
  }
}


void colvar::spin_angle::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;

  if (!atoms.noforce) {
    atoms.apply_colvar_force(fw);
  }
}

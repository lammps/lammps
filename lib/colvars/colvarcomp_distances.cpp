// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"



colvar::distance::distance(std::string const &conf)
  : cvc(conf)
{
  function_type = "distance";
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  init_total_force_params(conf);
}


colvar::distance::distance()
  : cvc()
{
  function_type = "distance";
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
}


void colvar::distance::calc_value()
{
  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    dist_v = group2->center_of_mass() - group1->center_of_mass();
  } else {
    dist_v = cvm::position_distance(group1->center_of_mass(),
                                    group2->center_of_mass());
  }
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


void colvar::distance::apply_force(colvarvalue const &force)
{
  if (!group1->noforce)
    group1->apply_colvar_force(force.real_value);

  if (!group2->noforce)
    group2->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(distance)



colvar::distance_vec::distance_vec(std::string const &conf)
  : distance(conf)
{
  function_type = "distance_vec";
  enable(f_cvc_com_based);
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_3vector);
}


colvar::distance_vec::distance_vec()
  : distance()
{
  function_type = "distance_vec";
  enable(f_cvc_com_based);
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_3vector);
}


void colvar::distance_vec::calc_value()
{
  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    x.rvector_value = group2->center_of_mass() - group1->center_of_mass();
  } else {
    x.rvector_value = cvm::position_distance(group1->center_of_mass(),
                                             group2->center_of_mass());
  }
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


cvm::real colvar::distance_vec::dist2(colvarvalue const &x1,
                                      colvarvalue const &x2) const
{
  return (cvm::position_distance(x1.rvector_value, x2.rvector_value)).norm2();
}


colvarvalue colvar::distance_vec::dist2_lgrad(colvarvalue const &x1,
                                              colvarvalue const &x2) const
{
  return 2.0 * cvm::position_distance(x2.rvector_value, x1.rvector_value);
}


colvarvalue colvar::distance_vec::dist2_rgrad(colvarvalue const &x1,
                                              colvarvalue const &x2) const
{
  return 2.0 * cvm::position_distance(x2.rvector_value, x1.rvector_value);
}



colvar::distance_z::distance_z(std::string const &conf)
  : cvc(conf)
{
  function_type = "distance_z";
  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
  x.type(colvarvalue::type_scalar);

  // TODO detect PBC from MD engine (in simple cases)
  // and then update period in real time
  if (period != 0.0)
    enable(f_cvc_periodic);

  if ((wrap_center != 0.0) && (period == 0.0)) {
    cvm::error("Error: wrapAround was defined in a distanceZ component,"
                " but its period has not been set.\n");
    return;
  }

  main = parse_group(conf, "main");
  ref1 = parse_group(conf, "ref");
  // this group is optional
  ref2 = parse_group(conf, "ref2", true);

  if ( ref2 ) {
    cvm::log("Using axis joining the centers of mass of groups \"ref\" and \"ref2\"");
    fixed_axis = false;
    if (key_lookup(conf, "axis"))
      cvm::log("Warning: explicit axis definition will be ignored!");
  } else {
    if (get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0))) {
      if (axis.norm2() == 0.0) {
        cvm::error("Axis vector is zero!");
        return;
      }
      if (axis.norm2() != 1.0) {
        axis = axis.unit();
        cvm::log("The normalized axis is: "+cvm::to_str(axis)+".\n");
      }
    }
    fixed_axis = true;
  }

  init_total_force_params(conf);

}


colvar::distance_z::distance_z()
{
  function_type = "distance_z";
  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
  x.type(colvarvalue::type_scalar);
}


void colvar::distance_z::calc_value()
{
  if (fixed_axis) {
    if (!is_enabled(f_cvc_pbc_minimum_image)) {
      dist_v = main->center_of_mass() - ref1->center_of_mass();
    } else {
      dist_v = cvm::position_distance(ref1->center_of_mass(),
                                      main->center_of_mass());
    }
  } else {

    if (!is_enabled(f_cvc_pbc_minimum_image)) {
      dist_v = main->center_of_mass() -
               (0.5 * (ref1->center_of_mass() + ref2->center_of_mass()));
      axis = ref2->center_of_mass() - ref1->center_of_mass();
    } else {
      dist_v = cvm::position_distance(0.5 * (ref1->center_of_mass() +
                                             ref2->center_of_mass()),
                                      main->center_of_mass());
      axis = cvm::position_distance(ref1->center_of_mass(),
                                    ref2->center_of_mass());
    }
    axis_norm = axis.norm();
    axis = axis.unit();
  }
  x.real_value = axis * dist_v;
  this->wrap(x);
}


void colvar::distance_z::calc_gradients()
{
  main->set_weighted_gradient( axis );

  if (fixed_axis) {
    ref1->set_weighted_gradient(-1.0 * axis);
  } else {
    if (!is_enabled(f_cvc_pbc_minimum_image)) {
      ref1->set_weighted_gradient( 1.0 / axis_norm *
                                   (main->center_of_mass() - ref2->center_of_mass() -
                                   x.real_value * axis ));
      ref2->set_weighted_gradient( 1.0 / axis_norm *
                                   (ref1->center_of_mass() - main->center_of_mass() +
                                   x.real_value * axis ));
    } else {
      ref1->set_weighted_gradient( 1.0 / axis_norm * (
        cvm::position_distance(ref2->center_of_mass(),
                               main->center_of_mass()) - x.real_value * axis ));
      ref2->set_weighted_gradient( 1.0 / axis_norm * (
        cvm::position_distance(main->center_of_mass(),
                               ref1->center_of_mass()) + x.real_value * axis ));
    }
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


void colvar::distance_z::apply_force(colvarvalue const &force)
{
  if (!ref1->noforce)
    ref1->apply_colvar_force(force.real_value);

  if (ref2 && !ref2->noforce)
    ref2->apply_colvar_force(force.real_value);

  if (!main->noforce)
    main->apply_colvar_force(force.real_value);
}


// Differences should always be wrapped around 0 (ignoring wrap_center)
cvm::real colvar::distance_z::dist2(colvarvalue const &x1,
                                    colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) {
    cvm::real shift = cvm::floor(diff/period + 0.5);
    diff -= shift * period;
  }
  return diff * diff;
}


colvarvalue colvar::distance_z::dist2_lgrad(colvarvalue const &x1,
                                            colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) {
    cvm::real shift = cvm::floor(diff/period + 0.5);
    diff -= shift * period;
  }
  return 2.0 * diff;
}


colvarvalue colvar::distance_z::dist2_rgrad(colvarvalue const &x1,
                                            colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) {
    cvm::real shift = cvm::floor(diff/period + 0.5);
    diff -= shift * period;
  }
  return (-2.0) * diff;
}


void colvar::distance_z::wrap(colvarvalue &x_unwrapped) const
{
  if (!is_enabled(f_cvc_periodic)) {
    // don't wrap if the period has not been set
    return;
  }
  cvm::real shift =
    cvm::floor((x_unwrapped.real_value - wrap_center) / period + 0.5);
  x_unwrapped.real_value -= shift * period;
}



colvar::distance_xy::distance_xy(std::string const &conf)
  : distance_z(conf)
{
  function_type = "distance_xy";
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
}


colvar::distance_xy::distance_xy()
  : distance_z()
{
  function_type = "distance_xy";
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
}


void colvar::distance_xy::calc_value()
{
  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    dist_v = main->center_of_mass() - ref1->center_of_mass();
  } else {
    dist_v = cvm::position_distance(ref1->center_of_mass(),
                                    main->center_of_mass());
  }
  if (!fixed_axis) {
    if (!is_enabled(f_cvc_pbc_minimum_image)) {
      v12 = ref2->center_of_mass() - ref1->center_of_mass();
    } else {
      v12 = cvm::position_distance(ref1->center_of_mass(),
                                   ref2->center_of_mass());
    }
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
    if (!is_enabled(f_cvc_pbc_minimum_image)) {
      v13 = main->center_of_mass() - ref1->center_of_mass();
    } else {
      v13 = cvm::position_distance(ref1->center_of_mass(),
                                   main->center_of_mass());
    }
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


void colvar::distance_xy::apply_force(colvarvalue const &force)
{
  if (!ref1->noforce)
    ref1->apply_colvar_force(force.real_value);

  if (ref2 && !ref2->noforce)
    ref2->apply_colvar_force(force.real_value);

  if (!main->noforce)
    main->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(distance_xy)



colvar::distance_dir::distance_dir(std::string const &conf)
  : distance(conf)
{
  function_type = "distance_dir";
  enable(f_cvc_com_based);
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_unit3vector);
}


colvar::distance_dir::distance_dir()
  : distance()
{
  function_type = "distance_dir";
  enable(f_cvc_com_based);
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_unit3vector);
}


void colvar::distance_dir::calc_value()
{
  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    dist_v = group2->center_of_mass() - group1->center_of_mass();
  } else {
    dist_v = cvm::position_distance(group1->center_of_mass(),
                                    group2->center_of_mass());
  }
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

  if (!group1->noforce)
    group1->apply_force(-1.0 * force_tang);

  if (!group2->noforce)
    group2->apply_force(       force_tang);
}


cvm::real colvar::distance_dir::dist2(colvarvalue const &x1,
                                      colvarvalue const &x2) const
{
  return (x1.rvector_value - x2.rvector_value).norm2();
}


colvarvalue colvar::distance_dir::dist2_lgrad(colvarvalue const &x1,
                                              colvarvalue const &x2) const
{
  return colvarvalue((x1.rvector_value - x2.rvector_value), colvarvalue::type_unit3vectorderiv);
}


colvarvalue colvar::distance_dir::dist2_rgrad(colvarvalue const &x1,
                                              colvarvalue const &x2) const
{
  return colvarvalue((x2.rvector_value - x1.rvector_value), colvarvalue::type_unit3vectorderiv);
}



colvar::distance_inv::distance_inv(std::string const &conf)
  : cvc(conf)
{
  function_type = "distance_inv";
  init_as_distance();

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  get_keyval(conf, "exponent", exponent, 6);
  if (exponent%2) {
    cvm::error("Error: odd exponent provided, can only use even ones.\n");
    return;
  }
  if (exponent <= 0) {
    cvm::error("Error: negative or zero exponent provided.\n");
    return;
  }

  for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
    for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
      if (ai1->id == ai2->id) {
        cvm::error("Error: group1 and group2 have some atoms in common: this is not allowed for distanceInv.\n");
        return;
      }
    }
  }

  if (is_enabled(f_cvc_debug_gradient)) {
    cvm::log("Warning: debugGradients will not give correct results "
             "for distanceInv, because its value and gradients are computed "
             "simultaneously.\n");
  }
}


void colvar::distance_inv::calc_value()
{
  x.real_value = 0.0;
  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
      for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
        cvm::rvector const dv = ai2->pos - ai1->pos;
        cvm::real const d2 = dv.norm2();
        cvm::real const dinv = cvm::integer_power(d2, -1*(exponent/2));
        x.real_value += dinv;
        cvm::rvector const dsumddv = -1.0*(exponent/2) * dinv/d2 * 2.0 * dv;
        ai1->grad += -1.0 * dsumddv;
        ai2->grad +=        dsumddv;
      }
    }
  } else {
    for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
      for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
        cvm::rvector const dv = cvm::position_distance(ai1->pos, ai2->pos);
        cvm::real const d2 = dv.norm2();
        cvm::real const dinv = cvm::integer_power(d2, -1*(exponent/2));
        x.real_value += dinv;
        cvm::rvector const dsumddv = -1.0*(exponent/2) * dinv/d2 * 2.0 * dv;
        ai1->grad += -1.0 * dsumddv;
        ai2->grad +=        dsumddv;
      }
    }
  }

  x.real_value *= 1.0 / cvm::real(group1->size() * group2->size());
  x.real_value = cvm::pow(x.real_value, -1.0/cvm::real(exponent));

  cvm::real const dxdsum = (-1.0/(cvm::real(exponent))) *
    cvm::integer_power(x.real_value, exponent+1) /
    cvm::real(group1->size() * group2->size());
  for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
    ai1->grad *= dxdsum;
  }
  for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
    ai2->grad *= dxdsum;
  }
}


void colvar::distance_inv::calc_gradients()
{
}


void colvar::distance_inv::apply_force(colvarvalue const &force)
{
  if (!group1->noforce)
    group1->apply_colvar_force(force.real_value);

  if (!group2->noforce)
    group2->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(distance_inv)



colvar::distance_pairs::distance_pairs(std::string const &conf)
  : cvc(conf)
{
  function_type = "distance_pairs";

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  x.type(colvarvalue::type_vector);
  disable(f_cvc_explicit_gradient);
  x.vector1d_value.resize(group1->size() * group2->size());
}


colvar::distance_pairs::distance_pairs()
{
  function_type = "distance_pairs";
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_vector);
}


void colvar::distance_pairs::calc_value()
{
  x.vector1d_value.resize(group1->size() * group2->size());

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    size_t i1, i2;
    for (i1 = 0; i1 < group1->size(); i1++) {
      for (i2 = 0; i2 < group2->size(); i2++) {
        cvm::rvector const dv = (*group2)[i2].pos - (*group1)[i1].pos;
        cvm::real const d = dv.norm();
        x.vector1d_value[i1*group2->size() + i2] = d;
        (*group1)[i1].grad = -1.0 * dv.unit();
        (*group2)[i2].grad =  dv.unit();
      }
    }
  } else {
    size_t i1, i2;
    for (i1 = 0; i1 < group1->size(); i1++) {
      for (i2 = 0; i2 < group2->size(); i2++) {
        cvm::rvector const dv = cvm::position_distance((*group1)[i1].pos,
                                                       (*group2)[i2].pos);
        cvm::real const d = dv.norm();
        x.vector1d_value[i1*group2->size() + i2] = d;
        (*group1)[i1].grad = -1.0 * dv.unit();
        (*group2)[i2].grad =  dv.unit();
      }
    }
  }
}


void colvar::distance_pairs::calc_gradients()
{
  // will be calculated on the fly in apply_force()
}


void colvar::distance_pairs::apply_force(colvarvalue const &force)
{
  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    size_t i1, i2;
    for (i1 = 0; i1 < group1->size(); i1++) {
      for (i2 = 0; i2 < group2->size(); i2++) {
        cvm::rvector const dv = (*group2)[i2].pos - (*group1)[i1].pos;
        (*group1)[i1].apply_force(force[i1*group2->size() + i2] * (-1.0) * dv.unit());
        (*group2)[i2].apply_force(force[i1*group2->size() + i2] * dv.unit());
      }
    }
  } else {
    size_t i1, i2;
    for (i1 = 0; i1 < group1->size(); i1++) {
      for (i2 = 0; i2 < group2->size(); i2++) {
        cvm::rvector const dv = cvm::position_distance((*group1)[i1].pos,
                                                       (*group2)[i2].pos);
        (*group1)[i1].apply_force(force[i1*group2->size() + i2] * (-1.0) * dv.unit());
        (*group2)[i2].apply_force(force[i1*group2->size() + i2] * dv.unit());
      }
    }
  }
}



colvar::dipole_magnitude::dipole_magnitude(std::string const &conf)
  : cvc(conf)
{
  function_type = "dipole_magnitude";
  atoms = parse_group(conf, "atoms");
  init_total_force_params(conf);
  x.type(colvarvalue::type_scalar);
}


colvar::dipole_magnitude::dipole_magnitude(cvm::atom const &a1)
{
  atoms = new cvm::atom_group(std::vector<cvm::atom>(1, a1));
  register_atom_group(atoms);
  x.type(colvarvalue::type_scalar);
}


colvar::dipole_magnitude::dipole_magnitude()
{
  function_type = "dipole_magnitude";
  x.type(colvarvalue::type_scalar);
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

  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    ai->grad = (ai->charge - aux1*ai->mass) * dipVunit;
  }
}


void colvar::dipole_magnitude::apply_force(colvarvalue const &force)
{
  if (!atoms->noforce) {
    atoms->apply_colvar_force(force.real_value);
  }
}


simple_scalar_dist_functions(dipole_magnitude)



colvar::gyration::gyration(std::string const &conf)
  : cvc(conf)
{
  function_type = "gyration";
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  atoms = parse_group(conf, "atoms");

  if (atoms->b_user_defined_fit) {
    cvm::log("WARNING: explicit fitting parameters were provided for atom group \"atoms\".");
  } else {
    atoms->b_center = true;
    atoms->ref_pos.assign(1, cvm::atom_pos(0.0, 0.0, 0.0));
    atoms->fit_gradients.assign(atoms->size(), cvm::rvector(0.0, 0.0, 0.0));
  }
}


void colvar::gyration::calc_value()
{
  x.real_value = 0.0;
  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    x.real_value += (ai->pos).norm2();
  }
  x.real_value = cvm::sqrt(x.real_value / cvm::real(atoms->size()));
}


void colvar::gyration::calc_gradients()
{
  cvm::real const drdx = 1.0/(cvm::real(atoms->size()) * x.real_value);
  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    ai->grad = drdx * ai->pos;
  }
}


void colvar::gyration::calc_force_invgrads()
{
  atoms->read_total_forces();

  cvm::real const dxdr = 1.0/x.real_value;
  ft.real_value = 0.0;

  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    ft.real_value += dxdr * ai->pos * ai->total_force;
  }
}


void colvar::gyration::calc_Jacobian_derivative()
{
  jd = x.real_value ? (3.0 * cvm::real(atoms->size()) - 4.0) / x.real_value : 0.0;
}


void colvar::gyration::apply_force(colvarvalue const &force)
{
  if (!atoms->noforce)
    atoms->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(gyration)



colvar::inertia::inertia(std::string const &conf)
  : gyration(conf)
{
  function_type = "inertia";
  init_as_distance();
}


void colvar::inertia::calc_value()
{
  x.real_value = 0.0;
  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    x.real_value += (ai->pos).norm2();
  }
}


void colvar::inertia::calc_gradients()
{
  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    ai->grad = 2.0 * ai->pos;
  }
}


void colvar::inertia::apply_force(colvarvalue const &force)
{
  if (!atoms->noforce)
    atoms->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(inertia_z)



colvar::inertia_z::inertia_z(std::string const &conf)
  : inertia(conf)
{
  function_type = "inertia_z";
  init_as_distance();
  if (get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0))) {
    if (axis.norm2() == 0.0) {
      cvm::error("Axis vector is zero!", INPUT_ERROR);
      return;
    }
    if (axis.norm2() != 1.0) {
      axis = axis.unit();
      cvm::log("The normalized axis is: "+cvm::to_str(axis)+".\n");
    }
  }
}


void colvar::inertia_z::calc_value()
{
  x.real_value = 0.0;
  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    cvm::real const iprod = ai->pos * axis;
    x.real_value += iprod * iprod;
  }
}


void colvar::inertia_z::calc_gradients()
{
  for (cvm::atom_iter ai = atoms->begin(); ai != atoms->end(); ai++) {
    ai->grad = 2.0 * (ai->pos * axis) * axis;
  }
}


void colvar::inertia_z::apply_force(colvarvalue const &force)
{
  if (!atoms->noforce)
    atoms->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(inertia)




colvar::rmsd::rmsd(std::string const &conf)
  : cvc(conf)
{
  function_type = "rmsd";
  init_as_distance();

  provide(f_cvc_inv_gradient);

  atoms = parse_group(conf, "atoms");

  if (!atoms || atoms->size() == 0) {
    cvm::error("Error: \"atoms\" must contain at least 1 atom to compute RMSD.");
    return;
  }

  bool b_Jacobian_derivative = true;
  if (atoms->fitting_group != NULL && b_Jacobian_derivative) {
    cvm::log("The option \"fittingGroup\" (alternative group for fitting) was enabled: "
              "Jacobian derivatives of the RMSD will not be calculated.\n");
    b_Jacobian_derivative = false;
  }
  if (b_Jacobian_derivative) provide(f_cvc_Jacobian);

  // the following is a simplified version of the corresponding atom group options;
  // we need this because the reference coordinates defined inside the atom group
  // may be used only for fitting, and even more so if fitting_group is used
  if (get_keyval(conf, "refPositions", ref_pos, ref_pos)) {
    cvm::log("Using reference positions from configuration file to calculate the variable.\n");
    if (ref_pos.size() != atoms->size()) {
      cvm::error("Error: the number of reference positions provided ("+
                  cvm::to_str(ref_pos.size())+
                  ") does not match the number of atoms of group \"atoms\" ("+
                  cvm::to_str(atoms->size())+").\n");
      return;
    }
  } else { // Only look for ref pos file if ref positions not already provided
    std::string ref_pos_file;
    if (get_keyval(conf, "refPositionsFile", ref_pos_file, std::string(""))) {

      if (ref_pos.size()) {
        cvm::error("Error: cannot specify \"refPositionsFile\" and "
                          "\"refPositions\" at the same time.\n");
        return;
      }

      std::string ref_pos_col;
      double ref_pos_col_value=0.0;

      if (get_keyval(conf, "refPositionsCol", ref_pos_col, std::string(""))) {
        // if provided, use PDB column to select coordinates
        bool found = get_keyval(conf, "refPositionsColValue", ref_pos_col_value, 0.0);
        if (found && ref_pos_col_value==0.0) {
          cvm::error("Error: refPositionsColValue, "
                     "if provided, must be non-zero.\n");
          return;
        }
      }

      ref_pos.resize(atoms->size());

      cvm::load_coords(ref_pos_file.c_str(), &ref_pos, atoms,
                       ref_pos_col, ref_pos_col_value);
    } else {
      cvm::error("Error: no reference positions for RMSD; use either refPositions of refPositionsFile.");
      return;
    }
  }

  if (ref_pos.size() != atoms->size()) {
    cvm::error("Error: found " + cvm::to_str(ref_pos.size()) +
                    " reference positions for RMSD; expected " + cvm::to_str(atoms->size()));
    return;
  }

  if (atoms->b_user_defined_fit) {
    cvm::log("WARNING: explicit fitting parameters were provided for atom group \"atoms\".");
  } else {
    // Default: fit everything
    cvm::log("Enabling \"centerReference\" and \"rotateReference\", to minimize RMSD before calculating it as a variable: "
              "if this is not the desired behavior, disable them explicitly within the \"atoms\" block.\n");
    atoms->b_center = true;
    atoms->b_rotate = true;
    // default case: reference positions for calculating the rmsd are also those used
    // for fitting
    atoms->ref_pos = ref_pos;
    atoms->center_ref_pos();

    cvm::log("This is a standard minimum RMSD, derivatives of the optimal rotation "
              "will not be computed as they cancel out in the gradients.");
    atoms->disable(f_ag_fit_gradients);

    // request the calculation of the derivatives of the rotation defined by the atom group
    atoms->rot.request_group1_gradients(atoms->size());
    // request derivatives of optimal rotation wrt reference coordinates for Jacobian:
    // this is only required for ABF, but we do both groups here for better caching
    atoms->rot.request_group2_gradients(atoms->size());
  }
}


void colvar::rmsd::calc_value()
{
  // rotational-translational fit is handled by the atom group

  x.real_value = 0.0;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    x.real_value += ((*atoms)[ia].pos - ref_pos[ia]).norm2();
  }
  x.real_value /= cvm::real(atoms->size()); // MSD
  x.real_value = cvm::sqrt(x.real_value);
}


void colvar::rmsd::calc_gradients()
{
  cvm::real const drmsddx2 = (x.real_value > 0.0) ?
    0.5 / (x.real_value * cvm::real(atoms->size())) :
    0.0;

  for (size_t ia = 0; ia < atoms->size(); ia++) {
    (*atoms)[ia].grad = (drmsddx2 * 2.0 * ((*atoms)[ia].pos - ref_pos[ia]));
  }
}


void colvar::rmsd::apply_force(colvarvalue const &force)
{
  if (!atoms->noforce)
    atoms->apply_colvar_force(force.real_value);
}


void colvar::rmsd::calc_force_invgrads()
{
  atoms->read_total_forces();
  ft.real_value = 0.0;

  // Note: gradient square norm is 1/N_atoms

  for (size_t ia = 0; ia < atoms->size(); ia++) {
    ft.real_value += (*atoms)[ia].grad * (*atoms)[ia].total_force;
  }
  ft.real_value *= atoms->size();
}


void colvar::rmsd::calc_Jacobian_derivative()
{
  // divergence of the rotated coordinates (including only derivatives of the rotation matrix)
  cvm::real rotation_term = 0.0;

  // The rotation term only applies is coordinates are rotated
  if (atoms->b_rotate) {

    // gradient of the rotation matrix
    cvm::matrix2d<cvm::rvector> grad_rot_mat(3, 3);
    // gradients of products of 2 quaternion components
    cvm::rvector g11, g22, g33, g01, g02, g03, g12, g13, g23;
    for (size_t ia = 0; ia < atoms->size(); ia++) {

      // Gradient of optimal quaternion wrt current Cartesian position
      cvm::vector1d<cvm::rvector> &dq = atoms->rot.dQ0_1[ia];

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

      cvm::atom_pos &y = ref_pos[ia];

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
  cvm::real translation_term = atoms->b_center ? 3.0 : 0.0;

  jd.real_value = x.real_value > 0.0 ?
    (3.0 * atoms->size() - 1.0 - translation_term - rotation_term) / x.real_value :
    0.0;
}


simple_scalar_dist_functions(rmsd)



colvar::eigenvector::eigenvector(std::string const &conf)
  : cvc(conf)
{
  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  function_type = "eigenvector";
  x.type(colvarvalue::type_scalar);

  atoms = parse_group(conf, "atoms");

  {
    bool const b_inline = get_keyval(conf, "refPositions", ref_pos, ref_pos);

    if (b_inline) {
      cvm::log("Using reference positions from input file.\n");
      if (ref_pos.size() != atoms->size()) {
        cvm::error("Error: reference positions do not "
                   "match the number of requested atoms.\n");
        return;
      }
    }

    std::string file_name;
    if (get_keyval(conf, "refPositionsFile", file_name)) {

      if (b_inline) {
        cvm::error("Error: refPositions and refPositionsFile cannot be specified at the same time.\n");
        return;
      }

      std::string file_col;
      double file_col_value=0.0;
      if (get_keyval(conf, "refPositionsCol", file_col, std::string(""))) {
        // use PDB flags if column is provided
        bool found = get_keyval(conf, "refPositionsColValue", file_col_value, 0.0);
        if (found && file_col_value==0.0) {
          cvm::error("Error: refPositionsColValue, "
                            "if provided, must be non-zero.\n");
          return;
        }
      }

      ref_pos.resize(atoms->size());
      cvm::load_coords(file_name.c_str(), &ref_pos, atoms,
                       file_col, file_col_value);
    }
  }

  if (ref_pos.size() == 0) {
    cvm::error("Error: reference positions were not provided.\n", INPUT_ERROR);
    return;
  }

  if (ref_pos.size() != atoms->size()) {
    cvm::error("Error: reference positions do not "
               "match the number of requested atoms.\n", INPUT_ERROR);
    return;
  }

  // save for later the geometric center of the provided positions (may not be the origin)
  cvm::rvector ref_pos_center(0.0, 0.0, 0.0);
  for (size_t i = 0; i < atoms->size(); i++) {
    ref_pos_center += ref_pos[i];
  }
  ref_pos_center *= 1.0 / atoms->size();

  if (atoms->b_user_defined_fit) {
    cvm::log("WARNING: explicit fitting parameters were provided for atom group \"atoms\".\n");
  } else {
    // default: fit everything
    cvm::log("Enabling \"centerReference\" and \"rotateReference\", to minimize RMSD before calculating the vector projection: "
              "if this is not the desired behavior, disable them explicitly within the \"atoms\" block.\n");
    atoms->b_center = true;
    atoms->b_rotate = true;
    atoms->ref_pos = ref_pos;
    atoms->center_ref_pos();
    atoms->disable(f_ag_fit_gradients); // cancel out if group is fitted on itself
                                        // and cvc is translationally invariant

    // request the calculation of the derivatives of the rotation defined by the atom group
    atoms->rot.request_group1_gradients(atoms->size());
    // request derivatives of optimal rotation wrt reference coordinates for Jacobian:
    // this is only required for ABF, but we do both groups here for better caching
    atoms->rot.request_group2_gradients(atoms->size());
  }

  {
    bool const b_inline = get_keyval(conf, "vector", eigenvec, eigenvec);
    // now load the eigenvector
    if (b_inline) {
      cvm::log("Using vector components from input file.\n");
      if (eigenvec.size() != atoms->size()) {
        cvm::error("Error: vector components do not "
                          "match the number of requested atoms->\n");
        return;
      }
    }

    std::string file_name;
    if (get_keyval(conf, "vectorFile", file_name)) {

      if (b_inline) {
        cvm::error("Error: vector and vectorFile cannot be specified at the same time.\n");
        return;
      }

      std::string file_col;
      double file_col_value=0.0;
      if (get_keyval(conf, "vectorCol", file_col, std::string(""))) {
        // use PDB flags if column is provided
        bool found = get_keyval(conf, "vectorColValue", file_col_value, 0.0);
        if (found && file_col_value==0.0) {
          cvm::error("Error: vectorColValue, if provided, must be non-zero.\n");
          return;
        }
      }

      eigenvec.resize(atoms->size());
      cvm::load_coords(file_name.c_str(), &eigenvec, atoms,
                       file_col, file_col_value);
    }
  }

  if (!ref_pos.size() || !eigenvec.size()) {
    cvm::error("Error: both reference coordinates"
                      "and eigenvector must be defined.\n");
    return;
  }

  cvm::atom_pos eig_center(0.0, 0.0, 0.0);
  for (size_t eil = 0; eil < atoms->size(); eil++) {
    eig_center += eigenvec[eil];
  }
  eig_center *= 1.0 / atoms->size();
  cvm::log("Geometric center of the provided vector: "+cvm::to_str(eig_center)+"\n");

  bool b_difference_vector = false;
  get_keyval(conf, "differenceVector", b_difference_vector, false);

  if (b_difference_vector) {

    if (atoms->b_center) {
      // both sets should be centered on the origin for fitting
      for (size_t i = 0; i < atoms->size(); i++) {
        eigenvec[i] -= eig_center;
        ref_pos[i]  -= ref_pos_center;
      }
    }
    if (atoms->b_rotate) {
      atoms->rot.calc_optimal_rotation(eigenvec, ref_pos);
      for (size_t i = 0; i < atoms->size(); i++) {
        eigenvec[i] = atoms->rot.rotate(eigenvec[i]);
      }
    }
    cvm::log("\"differenceVector\" is on: subtracting the reference positions from the provided vector: v = v - x0.\n");
    for (size_t i = 0; i < atoms->size(); i++) {
      eigenvec[i] -= ref_pos[i];
    }
    if (atoms->b_center) {
      // bring back the ref positions to where they were
      for (size_t i = 0; i < atoms->size(); i++) {
        ref_pos[i] += ref_pos_center;
      }
    }

  } else {
    cvm::log("Centering the provided vector to zero.\n");
    for (size_t i = 0; i < atoms->size(); i++) {
      eigenvec[i] -= eig_center;
    }
  }

  // cvm::log("The first three components(v1x, v1y, v1z) of the resulting vector are: "+cvm::to_str (eigenvec[0])+".\n");

  // for inverse gradients
  eigenvec_invnorm2 = 0.0;
  for (size_t ein = 0; ein < atoms->size(); ein++) {
    eigenvec_invnorm2 += eigenvec[ein].norm2();
  }
  eigenvec_invnorm2 = 1.0 / eigenvec_invnorm2;

  if (b_difference_vector) {
    cvm::log("\"differenceVector\" is on: normalizing the vector.\n");
    for (size_t i = 0; i < atoms->size(); i++) {
      eigenvec[i] *= eigenvec_invnorm2;
    }
  } else {
    cvm::log("The norm of the vector is |v| = "+cvm::to_str(eigenvec_invnorm2)+".\n");
  }
}


void colvar::eigenvector::calc_value()
{
  x.real_value = 0.0;
  for (size_t i = 0; i < atoms->size(); i++) {
    x.real_value += ((*atoms)[i].pos - ref_pos[i]) * eigenvec[i];
  }
}


void colvar::eigenvector::calc_gradients()
{
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    (*atoms)[ia].grad = eigenvec[ia];
  }
}


void colvar::eigenvector::apply_force(colvarvalue const &force)
{
  if (!atoms->noforce)
    atoms->apply_colvar_force(force.real_value);
}


void colvar::eigenvector::calc_force_invgrads()
{
  atoms->read_total_forces();
  ft.real_value = 0.0;

  for (size_t ia = 0; ia < atoms->size(); ia++) {
    ft.real_value += eigenvec_invnorm2 * (*atoms)[ia].grad *
      (*atoms)[ia].total_force;
  }
}


void colvar::eigenvector::calc_Jacobian_derivative()
{
  // gradient of the rotation matrix
  cvm::matrix2d<cvm::rvector> grad_rot_mat(3, 3);
  cvm::quaternion &quat0 = atoms->rot.q;

  // gradients of products of 2 quaternion components
  cvm::rvector g11, g22, g33, g01, g02, g03, g12, g13, g23;

  cvm::real sum = 0.0;

  for (size_t ia = 0; ia < atoms->size(); ia++) {

    // Gradient of optimal quaternion wrt current Cartesian position
    // trick: d(R^-1)/dx = d(R^t)/dx = (dR/dx)^t
    // we can just transpose the derivatives of the direct matrix
    cvm::vector1d<cvm::rvector> &dq_1 = atoms->rot.dQ0_1[ia];

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


simple_scalar_dist_functions(eigenvector)



colvar::cartesian::cartesian(std::string const &conf)
  : cvc(conf)
{
  function_type = "cartesian";

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
    cvm::error("Error: a \"cartesian\" component was defined with all three axes disabled.\n");
    return;
  }

  x.type(colvarvalue::type_vector);
  disable(f_cvc_explicit_gradient);
  x.vector1d_value.resize(atoms->size() * axes.size());
}


void colvar::cartesian::calc_value()
{
  size_t const dim = axes.size();
  size_t ia, j;
  for (ia = 0; ia < atoms->size(); ia++) {
    for (j = 0; j < dim; j++) {
      x.vector1d_value[dim*ia + j] = (*atoms)[ia].pos[axes[j]];
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
    for (ia = 0; ia < atoms->size(); ia++) {
      for (j = 0; j < dim; j++) {
        f[axes[j]] = force.vector1d_value[dim*ia + j];
      }
      (*atoms)[ia].apply_force(f);
    }
  }
}


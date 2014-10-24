/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarcomp.h"

#include <cmath>


colvar::angle::angle(std::string const &conf)
  : cvc(conf)
{
  function_type = "angle";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  parse_group(conf, "group1", group1);
  parse_group(conf, "group2", group2);
  parse_group(conf, "group3", group3);
  atom_groups.push_back(&group1);
  atom_groups.push_back(&group2);
  atom_groups.push_back(&group3);
  if (get_keyval(conf, "oneSiteSystemForce", b_1site_force, false)) {
    cvm::log("Computing system force on group 1 only");
  }
  x.type(colvarvalue::type_scalar);
}


colvar::angle::angle(cvm::atom const &a1,
                      cvm::atom const &a2,
                      cvm::atom const &a3)
  : group1(std::vector<cvm::atom> (1, a1)),
    group2(std::vector<cvm::atom> (1, a2)),
    group3(std::vector<cvm::atom> (1, a3))
{
  function_type = "angle";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  b_1site_force = false;
  atom_groups.push_back(&group1);
  atom_groups.push_back(&group2);
  atom_groups.push_back(&group3);

  x.type(colvarvalue::type_scalar);
}


colvar::angle::angle()
{
  function_type = "angle";
  x.type(colvarvalue::type_scalar);
}


void colvar::angle::calc_value()
{
  group1.read_positions();
  group2.read_positions();
  group3.read_positions();

  cvm::atom_pos const g1_pos = group1.center_of_mass();
  cvm::atom_pos const g2_pos = group2.center_of_mass();
  cvm::atom_pos const g3_pos = group3.center_of_mass();

  r21  = cvm::position_distance(g2_pos, g1_pos);
  r21l = r21.norm();
  r23  = cvm::position_distance(g2_pos, g3_pos);
  r23l = r23.norm();

  cvm::real     const cos_theta = (r21*r23)/(r21l*r23l);

  x.real_value = (180.0/PI) * std::acos(cos_theta);
}


void colvar::angle::calc_gradients()
{
  size_t i;
  cvm::real const cos_theta = (r21*r23)/(r21l*r23l);
  cvm::real const dxdcos = -1.0 / std::sqrt(1.0 - cos_theta*cos_theta);

  dxdr1 = (180.0/PI) * dxdcos *
    (1.0/r21l) * ( r23/r23l + (-1.0) * cos_theta * r21/r21l );

  dxdr3 = (180.0/PI) * dxdcos *
    (1.0/r23l) * ( r21/r21l + (-1.0) * cos_theta * r23/r23l );

  for (i = 0; i < group1.size(); i++) {
    group1[i].grad = (group1[i].mass/group1.total_mass) *
      (dxdr1);
  }

  for (i = 0; i < group2.size(); i++) {
    group2[i].grad = (group2[i].mass/group2.total_mass) *
      (dxdr1 + dxdr3) * (-1.0);
  }

  for (i = 0; i < group3.size(); i++) {
    group3[i].grad = (group3[i].mass/group3.total_mass) *
      (dxdr3);
  }
}

void colvar::angle::calc_force_invgrads()
{
  // This uses a force measurement on groups 1 and 3 only
  // to keep in line with the implicit variable change used to
  // evaluate the Jacobian term (essentially polar coordinates
  // centered on group2, which means group2 is kept fixed
  // when propagating changes in the angle)

  if (b_1site_force) {
    group1.read_system_forces();
    cvm::real norm_fact = 1.0 / dxdr1.norm2();
    ft.real_value = norm_fact * dxdr1 * group1.system_force();
  } else {
    group1.read_system_forces();
    group3.read_system_forces();
    cvm::real norm_fact = 1.0 / (dxdr1.norm2() + dxdr3.norm2());
    ft.real_value = norm_fact * ( dxdr1 * group1.system_force()
                                + dxdr3 * group3.system_force());
  }
  return;
}

void colvar::angle::calc_Jacobian_derivative()
{
  // det(J) = (2 pi) r^2 * sin(theta)
  // hence Jd = cot(theta)
  const cvm::real theta = x.real_value * PI / 180.0;
  jd = PI / 180.0 * (theta != 0.0 ? std::cos(theta) / std::sin(theta) : 0.0);
}


void colvar::angle::apply_force(colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force(force.real_value);

  if (!group2.noforce)
    group2.apply_colvar_force(force.real_value);

  if (!group3.noforce)
    group3.apply_colvar_force(force.real_value);
}




colvar::dihedral::dihedral(std::string const &conf)
  : cvc(conf)
{
  function_type = "dihedral";
  period = 360.0;
  b_periodic = true;
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  if (get_keyval(conf, "oneSiteSystemForce", b_1site_force, false)) {
    cvm::log("Computing system force on group 1 only");
  }
  parse_group(conf, "group1", group1);
  parse_group(conf, "group2", group2);
  parse_group(conf, "group3", group3);
  parse_group(conf, "group4", group4);
  atom_groups.push_back(&group1);
  atom_groups.push_back(&group2);
  atom_groups.push_back(&group3);
  atom_groups.push_back(&group4);

  x.type(colvarvalue::type_scalar);
}


colvar::dihedral::dihedral(cvm::atom const &a1,
                            cvm::atom const &a2,
                            cvm::atom const &a3,
                            cvm::atom const &a4)
  : group1(std::vector<cvm::atom> (1, a1)),
    group2(std::vector<cvm::atom> (1, a2)),
    group3(std::vector<cvm::atom> (1, a3)),
    group4(std::vector<cvm::atom> (1, a4))
{
  if (cvm::debug())
    cvm::log("Initializing dihedral object from atom groups.\n");

  function_type = "dihedral";
  period = 360.0;
  b_periodic = true;
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  b_1site_force = false;

  atom_groups.push_back(&group1);
  atom_groups.push_back(&group2);
  atom_groups.push_back(&group3);
  atom_groups.push_back(&group4);

  x.type(colvarvalue::type_scalar);

  if (cvm::debug())
    cvm::log("Done initializing dihedral object from atom groups.\n");
}


colvar::dihedral::dihedral()
{
  function_type = "dihedral";
  period = 360.0;
  b_periodic = true;
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  x.type(colvarvalue::type_scalar);
}


void colvar::dihedral::calc_value()
{
  group1.read_positions();
  group2.read_positions();
  group3.read_positions();
  group4.read_positions();

  cvm::atom_pos const g1_pos = group1.center_of_mass();
  cvm::atom_pos const g2_pos = group2.center_of_mass();
  cvm::atom_pos const g3_pos = group3.center_of_mass();
  cvm::atom_pos const g4_pos = group4.center_of_mass();

  // Usual sign convention: r12 = r2 - r1
  r12 = cvm::position_distance(g1_pos, g2_pos);
  r23 = cvm::position_distance(g2_pos, g3_pos);
  r34 = cvm::position_distance(g3_pos, g4_pos);

  cvm::rvector const n1 = cvm::rvector::outer(r12, r23);
  cvm::rvector const n2 = cvm::rvector::outer(r23, r34);

  cvm::real const cos_phi = n1 * n2;
  cvm::real const sin_phi = n1 * r34 * r23.norm();

  x.real_value = (180.0/PI) * std::atan2(sin_phi, cos_phi);
  this->wrap(x);
}


void colvar::dihedral::calc_gradients()
{
  cvm::rvector A = cvm::rvector::outer(r12, r23);
  cvm::real   rA = A.norm();
  cvm::rvector B = cvm::rvector::outer(r23, r34);
  cvm::real   rB = B.norm();
  cvm::rvector C = cvm::rvector::outer(r23, A);
  cvm::real   rC = C.norm();

  cvm::real const cos_phi = (A*B)/(rA*rB);
  cvm::real const sin_phi = (C*B)/(rC*rB);

  cvm::rvector f1, f2, f3;

  rB = 1.0/rB;
  B *= rB;

  if (std::fabs(sin_phi) > 0.1) {
    rA = 1.0/rA;
    A *= rA;
    cvm::rvector const dcosdA = rA*(cos_phi*A-B);
    cvm::rvector const dcosdB = rB*(cos_phi*B-A);
    rA = 1.0;

    cvm::real const K = (1.0/sin_phi) * (180.0/PI);

	f1 = K * cvm::rvector::outer(r23, dcosdA);
	f3 = K * cvm::rvector::outer(dcosdB, r23);
	f2 = K * (cvm::rvector::outer(dcosdA, r12)
		   +  cvm::rvector::outer(r34, dcosdB));
  }
  else {
    rC = 1.0/rC;
    C *= rC;
    cvm::rvector const dsindC = rC*(sin_phi*C-B);
    cvm::rvector const dsindB = rB*(sin_phi*B-C);
    rC = 1.0;

    cvm::real    const K = (-1.0/cos_phi) * (180.0/PI);

    f1.x = K*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
              - r23.x*r23.y*dsindC.y
              - r23.x*r23.z*dsindC.z);
    f1.y = K*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
              - r23.y*r23.z*dsindC.z
              - r23.y*r23.x*dsindC.x);
    f1.z = K*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
              - r23.z*r23.x*dsindC.x
              - r23.z*r23.y*dsindC.y);

    f3 = cvm::rvector::outer(dsindB, r23);
    f3 *= K;

    f2.x = K*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
              +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
              +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
              +dsindB.z*r34.y - dsindB.y*r34.z);
    f2.y = K*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
              +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
              +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
              +dsindB.x*r34.z - dsindB.z*r34.x);
    f2.z = K*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
              +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
              +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
              +dsindB.y*r34.x - dsindB.x*r34.y);
  }

  size_t i;
  for (i = 0; i < group1.size(); i++)
    group1[i].grad = (group1[i].mass/group1.total_mass) * (-f1);

  for (i = 0; i < group2.size(); i++)
    group2[i].grad = (group2[i].mass/group2.total_mass) * (-f2 + f1);

  for (i = 0; i < group3.size(); i++)
	group3[i].grad = (group3[i].mass/group3.total_mass) * (-f3 + f2);

  for (i = 0; i < group4.size(); i++)
    group4[i].grad = (group4[i].mass/group4.total_mass) * (f3);
}


void colvar::dihedral::calc_force_invgrads()
{
  cvm::rvector const u12 = r12.unit();
  cvm::rvector const u23 = r23.unit();
  cvm::rvector const u34 = r34.unit();

  cvm::real const d12 = r12.norm();
  cvm::real const d34 = r34.norm();

  cvm::rvector const cross1 = (cvm::rvector::outer(u23, u12)).unit();
  cvm::rvector const cross4 = (cvm::rvector::outer(u23, u34)).unit();

  cvm::real const dot1 = u23 * u12;
  cvm::real const dot4 = u23 * u34;

  cvm::real const fact1 = d12 * std::sqrt(1.0 - dot1 * dot1);
  cvm::real const fact4 = d34 * std::sqrt(1.0 - dot4 * dot4);

  group1.read_system_forces();
  if ( b_1site_force ) {
    // This is only measuring the force on group 1
    ft.real_value = PI/180.0 * fact1 * (cross1 * group1.system_force());
  } else {
    // Default case: use groups 1 and 4
    group4.read_system_forces();
    ft.real_value = PI/180.0 * 0.5 * (fact1 * (cross1 * group1.system_force())
				      + fact4 * (cross4 * group4.system_force()));
  }
}


void colvar::dihedral::calc_Jacobian_derivative()
{
  // With this choice of inverse gradient ("internal coordinates"), Jacobian correction is 0
  jd = 0.0;
}


void colvar::dihedral::apply_force(colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force(force.real_value);

  if (!group2.noforce)
    group2.apply_colvar_force(force.real_value);

  if (!group3.noforce)
    group3.apply_colvar_force(force.real_value);

  if (!group4.noforce)
    group4.apply_colvar_force(force.real_value);
}



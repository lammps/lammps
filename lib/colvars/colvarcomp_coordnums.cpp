// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <cmath>

#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvaratoms.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



template<bool calculate_gradients>
cvm::real colvar::coordnum::switching_function(cvm::real const &r0,
                                               int const &en,
                                               int const &ed,
                                               cvm::atom &A1,
                                               cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);
  cvm::real const l2 = diff.norm2()/(r0*r0);

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = std::pow(l2, en2);
  cvm::real const xd = std::pow(l2, ed2);
  cvm::real const func = (1.0-xn)/(1.0-xd);

  if (calculate_gradients) {
    cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2) - func*ed2*(xd/l2))*(-1.0);
    cvm::rvector const dl2dx = (2.0/(r0*r0))*diff;
    A1.grad += (-1.0)*dFdl2*dl2dx;
    A2.grad +=        dFdl2*dl2dx;
  }

  return func;
}


template<bool calculate_gradients>
cvm::real colvar::coordnum::switching_function(cvm::rvector const &r0_vec,
                                               int const &en,
                                               int const &ed,
                                               cvm::atom &A1,
                                               cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);
  cvm::rvector const scal_diff(diff.x/r0_vec.x, diff.y/r0_vec.y, diff.z/r0_vec.z);
  cvm::real const l2 = scal_diff.norm2();

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = std::pow(l2, en2);
  cvm::real const xd = std::pow(l2, ed2);
  cvm::real const func = (1.0-xn)/(1.0-xd);

  if (calculate_gradients) {
    cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2) - func*ed2*(xd/l2))*(-1.0);
    cvm::rvector const dl2dx((2.0/(r0_vec.x*r0_vec.x))*diff.x,
                             (2.0/(r0_vec.y*r0_vec.y))*diff.y,
                             (2.0/(r0_vec.z*r0_vec.z))*diff.z);
    A1.grad += (-1.0)*dFdl2*dl2dx;
    A2.grad +=        dFdl2*dl2dx;
  }
  return func;
}


colvar::coordnum::coordnum(std::string const &conf)
  : cvc(conf), b_anisotropic(false), b_group2_center_only(false)
{
  function_type = "coordnum";
  x.type(colvarvalue::type_scalar);

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  if (group1->b_dummy) {
    cvm::error("Error: only group2 is allowed to be a dummy atom\n");
    return;
  }

  bool const b_isotropic = get_keyval(conf, "cutoff", r0,
                                      cvm::real(4.0 * cvm::unit_angstrom()));

  if (get_keyval(conf, "cutoff3", r0_vec, cvm::rvector(4.0 * cvm::unit_angstrom(),
                                                       4.0 * cvm::unit_angstrom(),
                                                       4.0 * cvm::unit_angstrom()))) {
    if (b_isotropic) {
      cvm::error("Error: cannot specify \"cutoff\" and \"cutoff3\" at the same time.\n",
                 INPUT_ERROR);
      return;
    }

    b_anisotropic = true;
    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;
  }

  get_keyval(conf, "expNumer", en, int(6) );
  get_keyval(conf, "expDenom", ed, int(12));

  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponents provided, can only use even ones.\n", INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

  get_keyval(conf, "group2CenterOnly", b_group2_center_only, group2->b_dummy);
}


colvar::coordnum::coordnum()
  : b_anisotropic(false), b_group2_center_only(false)
{
  function_type = "coordnum";
  x.type(colvarvalue::type_scalar);
}


void colvar::coordnum::calc_value()
{
  x.real_value = 0.0;

  if (b_group2_center_only) {

    // create a fake atom to hold the group2 com coordinates
    cvm::atom group2_com_atom;
    group2_com_atom.pos = group2->center_of_mass();

    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        x.real_value += switching_function<false>(r0_vec, en, ed, *ai1, group2_com_atom);
    } else {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        x.real_value += switching_function<false>(r0, en, ed, *ai1, group2_com_atom);
    }

  } else {

    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
          x.real_value += switching_function<false>(r0_vec, en, ed, *ai1, *ai2);
        }
    } else {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
          x.real_value += switching_function<false>(r0, en, ed, *ai1, *ai2);
        }
    }
  }
}


void colvar::coordnum::calc_gradients()
{
  if (b_group2_center_only) {

    // create a fake atom to hold the group2 com coordinates
    cvm::atom group2_com_atom;
    group2_com_atom.pos = group2->center_of_mass();


    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        switching_function<true>(r0_vec, en, ed, *ai1, group2_com_atom);
    } else {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        switching_function<true>(r0, en, ed, *ai1, group2_com_atom);
    }

    group2->set_weighted_gradient(group2_com_atom.grad);

  } else {

    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
          switching_function<true>(r0_vec, en, ed, *ai1, *ai2);
        }
    } else {
      for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++)
        for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
          switching_function<true>(r0, en, ed, *ai1, *ai2);
        }
    }
  }
}


void colvar::coordnum::apply_force(colvarvalue const &force)
{
  if (!group1->noforce)
    group1->apply_colvar_force(force.real_value);

  if (!group2->noforce)
    group2->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(coordnum)



// h_bond member functions

colvar::h_bond::h_bond(std::string const &conf)
  : cvc(conf)
{
  if (cvm::debug())
    cvm::log("Initializing h_bond object.\n");

  function_type = "h_bond";
  x.type(colvarvalue::type_scalar);

  int a_num, d_num;
  get_keyval(conf, "acceptor", a_num, -1);
  get_keyval(conf, "donor",    d_num, -1);

  if ( (a_num == -1) || (d_num == -1) ) {
    cvm::error("Error: either acceptor or donor undefined.\n");
    return;
  }

  cvm::atom acceptor = cvm::atom(a_num);
  cvm::atom donor    = cvm::atom(d_num);
  register_atom_group(new cvm::atom_group);
  atom_groups[0]->add_atom(acceptor);
  atom_groups[0]->add_atom(donor);

  get_keyval(conf, "cutoff",   r0, (3.3 * cvm::unit_angstrom()));
  get_keyval(conf, "expNumer", en, 6);
  get_keyval(conf, "expDenom", ed, 8);

  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponents provided, can only use even ones.\n");
    return;
  }

  if (cvm::debug())
    cvm::log("Done initializing h_bond object.\n");
}


colvar::h_bond::h_bond(cvm::atom const &acceptor,
                       cvm::atom const &donor,
                       cvm::real r0_i, int en_i, int ed_i)
  : r0(r0_i), en(en_i), ed(ed_i)
{
  function_type = "h_bond";
  x.type(colvarvalue::type_scalar);

  register_atom_group(new cvm::atom_group);
  atom_groups[0]->add_atom(acceptor);
  atom_groups[0]->add_atom(donor);
}


colvar::h_bond::h_bond()
  : cvc()
{
  function_type = "h_bond";
  x.type(colvarvalue::type_scalar);
}


colvar::h_bond::~h_bond()
{
  delete atom_groups[0];
}


void colvar::h_bond::calc_value()
{
  x.real_value = colvar::coordnum::switching_function<false>(r0, en, ed, (*atom_groups[0])[0], (*atom_groups[0])[1]);
}


void colvar::h_bond::calc_gradients()
{
  colvar::coordnum::switching_function<true>(r0, en, ed, (*atom_groups[0])[0], (*atom_groups[0])[1]);
}


void colvar::h_bond::apply_force(colvarvalue const &force)
{
  (atom_groups[0])->apply_colvar_force(force);
}


simple_scalar_dist_functions(h_bond)



colvar::selfcoordnum::selfcoordnum(std::string const &conf)
  : cvc(conf)
{
  function_type = "selfcoordnum";
  x.type(colvarvalue::type_scalar);

  group1 = parse_group(conf, "group1");

  get_keyval(conf, "cutoff", r0, cvm::real(4.0 * cvm::unit_angstrom()));
  get_keyval(conf, "expNumer", en, int(6) );
  get_keyval(conf, "expDenom", ed, int(12));

  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponents provided, can only use even ones.\n");
    return;
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }
}


colvar::selfcoordnum::selfcoordnum()
{
  function_type = "selfcoordnum";
  x.type(colvarvalue::type_scalar);
}


void colvar::selfcoordnum::calc_value()
{
  x.real_value = 0.0;
  for (size_t i = 0; i < group1->size() - 1; i++) {
    for (size_t j = i + 1; j < group1->size(); j++) {
      x.real_value += colvar::coordnum::switching_function<false>(r0, en, ed, (*group1)[i], (*group1)[j]);
    }
  }
}


void colvar::selfcoordnum::calc_gradients()
{
  for (size_t i = 0; i < group1->size() - 1; i++) {
    for (size_t j = i + 1; j < group1->size(); j++) {
      colvar::coordnum::switching_function<true>(r0, en, ed, (*group1)[i], (*group1)[j]);
    }
  }
}

void colvar::selfcoordnum::apply_force(colvarvalue const &force)
{
  if (!group1->noforce) {
    group1->apply_colvar_force(force.real_value);
  }
}


simple_scalar_dist_functions(selfcoordnum)


// groupcoordnum member functions
colvar::groupcoordnum::groupcoordnum(std::string const &conf)
  : distance(conf), b_anisotropic(false)
{
  function_type = "groupcoordnum";
  x.type(colvarvalue::type_scalar);

  // group1 and group2 are already initialized by distance()
  if (group1->b_dummy || group2->b_dummy) {
    cvm::error("Error: neither group can be a dummy atom\n");
    return;
  }

  bool const b_scale = get_keyval(conf, "cutoff", r0,
                                  cvm::real(4.0 * cvm::unit_angstrom()));

  if (get_keyval(conf, "cutoff3", r0_vec,
                 cvm::rvector(4.0, 4.0, 4.0), parse_silent)) {

    if (b_scale) {
      cvm::error("Error: cannot specify \"scale\" and "
                       "\"scale3\" at the same time.\n");
      return;
    }
    b_anisotropic = true;
    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;
  }

  get_keyval(conf, "expNumer", en, int(6) );
  get_keyval(conf, "expDenom", ed, int(12));

  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponents provided, can only use even ones.\n");
    return;
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

}


colvar::groupcoordnum::groupcoordnum()
  : b_anisotropic(false)
{
  function_type = "groupcoordnum";
  x.type(colvarvalue::type_scalar);
}


template<bool calculate_gradients>
cvm::real colvar::groupcoordnum::switching_function(cvm::real const &r0,
                                                    int const &en,
                                                    int const &ed,
                                                    cvm::atom &A1,
                                                    cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);
  cvm::real const l2 = diff.norm2()/(r0*r0);

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = std::pow(l2, en2);
  cvm::real const xd = std::pow(l2, ed2);
  cvm::real const func = (1.0-xn)/(1.0-xd);

  if (calculate_gradients) {
    cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2) - func*ed2*(xd/l2))*(-1.0);
    cvm::rvector const dl2dx = (2.0/(r0*r0))*diff;
    A1.grad += (-1.0)*dFdl2*dl2dx;
    A2.grad +=        dFdl2*dl2dx;
  }

  return func;
}


#if 0  // AMG: I don't think there's any reason to support anisotropic,
       //      and I don't have those flags below in calc_value, but
       //      if I need them, I'll also need to uncomment this method
template<bool calculate_gradients>
cvm::real colvar::groupcoordnum::switching_function(cvm::rvector const &r0_vec,
                                                    int const &en,
                                                    int const &ed,
                                                    cvm::atom &A1,
                                                    cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);
  cvm::rvector const scal_diff(diff.x/r0_vec.x, diff.y/r0_vec.y, diff.z/r0_vec.z);
  cvm::real const l2 = scal_diff.norm2();

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = std::pow(l2, en2);
  cvm::real const xd = std::pow(l2, ed2);
  cvm::real const func = (1.0-xn)/(1.0-xd);

  if (calculate_gradients) {
    cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2) - func*ed2*(xd/l2))*(-1.0);
    cvm::rvector const dl2dx((2.0/(r0_vec.x*r0_vec.x))*diff.x,
                             (2.0/(r0_vec.y*r0_vec.y))*diff.y,
                             (2.0/(r0_vec.z*r0_vec.z))*diff.z);
    A1.grad += (-1.0)*dFdl2*dl2dx;
    A2.grad +=        dFdl2*dl2dx;
  }
  return func;
}
#endif


void colvar::groupcoordnum::calc_value()
{

  // create fake atoms to hold the com coordinates
  cvm::atom group1_com_atom;
  cvm::atom group2_com_atom;
  group1_com_atom.pos = group1->center_of_mass();
  group2_com_atom.pos = group2->center_of_mass();

  x.real_value = coordnum::switching_function<false>(r0, en, ed,
                                                     group1_com_atom, group2_com_atom);
}


void colvar::groupcoordnum::calc_gradients()
{
  cvm::atom group1_com_atom;
  cvm::atom group2_com_atom;
  group1_com_atom.pos = group1->center_of_mass();
  group2_com_atom.pos = group2->center_of_mass();

  coordnum::switching_function<true>(r0, en, ed, group1_com_atom, group2_com_atom);
  group1->set_weighted_gradient(group1_com_atom.grad);
  group2->set_weighted_gradient(group2_com_atom.grad);

}


void colvar::groupcoordnum::apply_force(colvarvalue const &force)
{
  if (!group1->noforce)
    group1->apply_colvar_force(force.real_value);

  if (!group2->noforce)
    group2->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(groupcoordnum)

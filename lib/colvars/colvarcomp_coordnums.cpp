/// -*- c++ -*-

#include <cmath>

#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvaratoms.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



template<bool calculate_gradients>
cvm::real colvar::coordnum::switching_function (cvm::real const &r0,
                                                int const &en,
                                                int const &ed,
                                                cvm::atom &A1,
                                                cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance (A1.pos, A2.pos);
  cvm::real const l2 = diff.norm2()/(r0*r0);

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = std::pow (l2, en2);
  cvm::real const xd = std::pow (l2, ed2);
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
cvm::real colvar::coordnum::switching_function (cvm::rvector const &r0_vec,
                                                int const &en,
                                                int const &ed,
                                                cvm::atom &A1,
                                                cvm::atom &A2)
{
  cvm::rvector const diff = cvm::position_distance (A1.pos, A2.pos);
  cvm::rvector const scal_diff (diff.x/r0_vec.x, diff.y/r0_vec.y, diff.z/r0_vec.z);
  cvm::real const l2 = scal_diff.norm2();

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = std::pow (l2, en2);
  cvm::real const xd = std::pow (l2, ed2);
  cvm::real const func = (1.0-xn)/(1.0-xd);

  if (calculate_gradients) {
    cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2) - func*ed2*(xd/l2))*(-1.0);
    cvm::rvector const dl2dx ((2.0/(r0_vec.x*r0_vec.x))*diff.x,
                              (2.0/(r0_vec.y*r0_vec.y))*diff.y,
                              (2.0/(r0_vec.z*r0_vec.z))*diff.z);
    A1.grad += (-1.0)*dFdl2*dl2dx;
    A2.grad +=        dFdl2*dl2dx;
  }
  return func;
}



colvar::coordnum::coordnum (std::string const &conf)
  : distance (conf), b_anisotropic (false), b_group2_center_only (false)
{
  function_type = "coordnum";
  x.type (colvarvalue::type_scalar);

  // group1 and group2 are already initialized by distance()
  if (group1.b_dummy)
    cvm::fatal_error ("Error: only group2 is allowed to be a dummy atom\n");


  // need to specify this explicitly because the distance() constructor
  // has set it to true
  b_inverse_gradients = false;

  bool const b_scale = get_keyval (conf, "cutoff", r0,
                                   cvm::real (4.0 * cvm::unit_angstrom()));

  if (get_keyval (conf, "cutoff3", r0_vec,
                  cvm::rvector (4.0, 4.0, 4.0), parse_silent)) {

    if (b_scale)
      cvm::fatal_error ("Error: cannot specify \"scale\" and "
                        "\"scale3\" at the same time.\n");
    b_anisotropic = true;
    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;
  }

  get_keyval (conf, "expNumer", en, int (6) );
  get_keyval (conf, "expDenom", ed, int (12));

  if ( (en%2) || (ed%2) ) {
    cvm::fatal_error ("Error: odd exponents provided, can only use even ones.\n");
  }

  get_keyval (conf, "group2CenterOnly", b_group2_center_only, group2.b_dummy);
}


colvar::coordnum::coordnum()
  : b_anisotropic (false), b_group2_center_only (false)
{
  function_type = "coordnum";
  x.type (colvarvalue::type_scalar);
}


void colvar::coordnum::calc_value()
{
  x.real_value = 0.0;

  if (b_group2_center_only) {

    // create a fake atom to hold the group2 com coordinates
    cvm::atom group2_com_atom;
    group2_com_atom.pos = group2.center_of_mass();

    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        x.real_value += switching_function<false> (r0_vec, en, ed, *ai1, group2_com_atom);
    } else {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        x.real_value += switching_function<false> (r0, en, ed, *ai1, group2_com_atom);
    }

  } else {

    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        for (cvm::atom_iter ai2 = group2.begin(); ai2 != group2.end(); ai2++) {
          x.real_value += switching_function<false> (r0_vec, en, ed, *ai1, *ai2);
        }
    } else {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        for (cvm::atom_iter ai2 = group2.begin(); ai2 != group2.end(); ai2++) {
          x.real_value += switching_function<false> (r0, en, ed, *ai1, *ai2);
        }
    }
  }
}


void colvar::coordnum::calc_gradients()
{
  if (b_group2_center_only) {

    // create a fake atom to hold the group2 com coordinates
    cvm::atom group2_com_atom;
    group2_com_atom.pos = group2.center_of_mass();


    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        switching_function<true> (r0_vec, en, ed, *ai1, group2_com_atom);
    } else {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        switching_function<true> (r0, en, ed, *ai1, group2_com_atom);
    }

    group2.set_weighted_gradient (group2_com_atom.grad);

  } else {

    if (b_anisotropic) {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        for (cvm::atom_iter ai2 = group2.begin(); ai2 != group2.end(); ai2++) {
          switching_function<true> (r0_vec, en, ed, *ai1, *ai2);
        }
    } else {
      for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++)
        for (cvm::atom_iter ai2 = group2.begin(); ai2 != group2.end(); ai2++) {
          switching_function<true> (r0, en, ed, *ai1, *ai2);
        }
    }
  }

  //   if (cvm::debug()) {
  //     for (size_t i = 0; i < group1.size(); i++) {
  //       cvm::log ("atom["+cvm::to_str (group1[i].id+1)+"] gradient: "+
  //                 cvm::to_str (group1[i].grad)+"\n");
  //     }

  //     for (size_t i = 0; i < group2.size(); i++) {
  //       cvm::log ("atom["+cvm::to_str (group2[i].id+1)+"] gradient: "+
  //                 cvm::to_str (group2[i].grad)+"\n");
  //     }
  //   }
}

void colvar::coordnum::apply_force (colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force (force.real_value);

  if (!group2.noforce)
    group2.apply_colvar_force (force.real_value);
}



// h_bond member functions

colvar::h_bond::h_bond (std::string const &conf)
  : cvc (conf)
{
  if (cvm::debug())
    cvm::log ("Initializing h_bond object.\n");

  function_type = "h_bond";
  x.type (colvarvalue::type_scalar);

  int a_num, d_num;
  get_keyval (conf, "acceptor", a_num, -1);
  get_keyval (conf, "donor",    d_num, -1);

  if ( (a_num == -1) || (d_num == -1) ) {
    cvm::fatal_error ("Error: either acceptor or donor undefined.\n");
  }

  acceptor = cvm::atom (a_num);
  donor    = cvm::atom (d_num);
  atom_groups.push_back (new cvm::atom_group);
  atom_groups[0]->add_atom (acceptor);
  atom_groups[0]->add_atom (donor);

  get_keyval (conf, "cutoff",   r0, (3.3 * cvm::unit_angstrom()));
  get_keyval (conf, "expNumer", en, 6);
  get_keyval (conf, "expDenom", ed, 8);

  if ( (en%2) || (ed%2) ) {
    cvm::fatal_error ("Error: odd exponents provided, can only use even ones.\n");
  }

  if (cvm::debug())
    cvm::log ("Done initializing h_bond object.\n");
}


colvar::h_bond::h_bond (cvm::atom const &acceptor_i,
                        cvm::atom const &donor_i,
                        cvm::real r0_i, int en_i, int ed_i)
  : acceptor (acceptor_i),
    donor (donor_i),
    r0 (r0_i), en (en_i), ed (ed_i)
{
  function_type = "h_bond";
  x.type (colvarvalue::type_scalar);

  atom_groups.push_back (new cvm::atom_group);
  atom_groups[0]->add_atom (acceptor);
  atom_groups[0]->add_atom (donor);
}

colvar::h_bond::h_bond()
  : cvc ()
{
  function_type = "h_bond";
  x.type (colvarvalue::type_scalar);
}


colvar::h_bond::~h_bond()
{
  for (unsigned int i=0; i<atom_groups.size(); i++) {
    delete atom_groups[i];
  }
}


void colvar::h_bond::calc_value()
{
  x.real_value = colvar::coordnum::switching_function<false> (r0, en, ed, acceptor, donor);
}


void colvar::h_bond::calc_gradients()
{
  colvar::coordnum::switching_function<true> (r0, en, ed, acceptor, donor);
  (*atom_groups[0])[0].grad = acceptor.grad;
  (*atom_groups[0])[1].grad = donor.grad;
}


void colvar::h_bond::apply_force (colvarvalue const &force)
{
  acceptor.apply_force (force.real_value * acceptor.grad);
  donor.apply_force    (force.real_value * donor.grad);
}




colvar::selfcoordnum::selfcoordnum (std::string const &conf)
 : distance (conf, false)
{
  function_type = "selfcoordnum";
  x.type (colvarvalue::type_scalar);

  // group1 is already initialized by distance()

  // need to specify this explicitly because the distance() constructor
  // has set it to true
  b_inverse_gradients = false;

  get_keyval (conf, "cutoff", r0, cvm::real (4.0 * cvm::unit_angstrom()));
  get_keyval (conf, "expNumer", en, int (6) );
  get_keyval (conf, "expDenom", ed, int (12));

  if ( (en%2) || (ed%2) ) {
    cvm::fatal_error ("Error: odd exponents provided, can only use even ones.\n");
  }
}


colvar::selfcoordnum::selfcoordnum()
{
  function_type = "selfcoordnum";
  x.type (colvarvalue::type_scalar);
}


void colvar::selfcoordnum::calc_value()
{
  x.real_value = 0.0;

  for (size_t i = 0; i < group1.size() - 1; i++)
    for (size_t j = i + 1; j < group1.size(); j++)
      x.real_value += colvar::coordnum::switching_function<false> (r0, en, ed, group1[i], group1[j]);
}


void colvar::selfcoordnum::calc_gradients()
{
  for (size_t i = 0; i < group1.size() - 1; i++)
    for (size_t j = i + 1; j < group1.size(); j++)
      colvar::coordnum::switching_function<true> (r0, en, ed, group1[i], group1[j]);
}

void colvar::selfcoordnum::apply_force (colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force (force.real_value);
}


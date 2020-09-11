// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvaratoms.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



template<int flags>
cvm::real colvar::coordnum::switching_function(cvm::real const &r0,
                                               cvm::rvector const &r0_vec,
                                               int en,
                                               int ed,
                                               cvm::atom &A1,
                                               cvm::atom &A2,
                                               bool **pairlist_elem,
                                               cvm::real pairlist_tol)
{
  if ((flags & ef_use_pairlist) && !(flags & ef_rebuild_pairlist)) {
    bool const within = **pairlist_elem;
    (*pairlist_elem)++;
    if (!within) {
      return 0.0;
    }
  }

  cvm::rvector const r0sq_vec(r0_vec.x*r0_vec.x,
                              r0_vec.y*r0_vec.y,
                              r0_vec.z*r0_vec.z);

  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);

  cvm::rvector const scal_diff(diff.x/((flags & ef_anisotropic) ?
                                       r0_vec.x : r0),
                               diff.y/((flags & ef_anisotropic) ?
                                       r0_vec.y : r0),
                               diff.z/((flags & ef_anisotropic) ?
                                       r0_vec.z : r0));
  cvm::real const l2 = scal_diff.norm2();

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = cvm::integer_power(l2, en2);
  cvm::real const xd = cvm::integer_power(l2, ed2);
  //The subtraction and division stretches the function back to the range of [0,1] from [pairlist_tol,1]
  cvm::real const func = (((1.0-xn)/(1.0-xd)) - pairlist_tol) / (1.0-pairlist_tol);

  if (flags & ef_rebuild_pairlist) {
    //Particles just outside of the cutoff also are considered if they come near.
    **pairlist_elem = (func > (-pairlist_tol * 0.5)) ? true : false;
    (*pairlist_elem)++;
  }
  //If the value is too small, we need to exclude it, rather than let it contribute to the sum or the gradients.
  if (func < 0)
    return 0;

  if (flags & ef_gradients) {
    //This is the old, completely correct expression for dFdl2:
    //cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2) -
    //                                        func*ed2*(xd/l2))*(-1.0);
    //This can become:
    //cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2)*(1.0-xn)/(1.0-xn) -
    //                                        func*ed2*(xd/l2))*(-1.0);
    //Recognizing that func = (1.0-xn)/(1.0-xd), we can group together the "func" and get a version of dFdl2 that is 0
    //when func=0, which lets us skip this gradient calculation when func=0.
    cvm::real const dFdl2 = func * ((ed2*xd/((1.0-xd)*l2)) - (en2*xn/((1.0-xn)*l2)));
    cvm::rvector const dl2dx((2.0/((flags & ef_anisotropic) ? r0sq_vec.x :
                                   r0*r0)) * diff.x,
                             (2.0/((flags & ef_anisotropic) ? r0sq_vec.y :
                                   r0*r0)) * diff.y,
                             (2.0/((flags & ef_anisotropic) ? r0sq_vec.z :
                                   r0*r0)) * diff.z);
    A1.grad += (-1.0)*dFdl2*dl2dx;
    A2.grad +=        dFdl2*dl2dx;
  }

  return func;
}


colvar::coordnum::coordnum(std::string const &conf)
  : cvc(conf), b_anisotropic(false), pairlist(NULL)

{
  function_type = "coordnum";
  x.type(colvarvalue::type_scalar);

  colvarproxy *proxy = cvm::main()->proxy;

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  if (group1 == NULL || group2 == NULL) {
    cvm::error("Error: failed to initialize atom groups.\n",
                INPUT_ERROR);
    return;
  }

  if (int atom_number = cvm::atom_group::overlap(*group1, *group2)) {
    cvm::error("Error: group1 and group2 share a common atom (number: " +
               cvm::to_str(atom_number) + ")\n", INPUT_ERROR);
    return;
  }

  if (group1->b_dummy) {
    cvm::error("Error: only group2 is allowed to be a dummy atom\n",
               INPUT_ERROR);
    return;
  }

  bool const b_isotropic = get_keyval(conf, "cutoff", r0,
                                      cvm::real(4.0 * proxy->angstrom_value));

  if (get_keyval(conf, "cutoff3", r0_vec,
                 cvm::rvector(4.0 * proxy->angstrom_value,
                              4.0 * proxy->angstrom_value,
                              4.0 * proxy->angstrom_value))) {
    if (b_isotropic) {
      cvm::error("Error: cannot specify \"cutoff\" and \"cutoff3\" "
                 "at the same time.\n",
                 INPUT_ERROR);
      return;
    }

    b_anisotropic = true;
    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;
  }

  get_keyval(conf, "expNumer", en, 6);
  get_keyval(conf, "expDenom", ed, 12);

  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
               INPUT_ERROR);
  }

  if ( (en <= 0) || (ed <= 0) ) {
    cvm::error("Error: negative exponent(s) provided.\n",
               INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

  get_keyval(conf, "group2CenterOnly", b_group2_center_only, group2->b_dummy);

  get_keyval(conf, "tolerance", tolerance, 0.0);
  if (tolerance > 0) {
    get_keyval(conf, "pairListFrequency", pairlist_freq, 100);
    if ( ! (pairlist_freq > 0) ) {
      cvm::error("Error: non-positive pairlistfrequency provided.\n",
                 INPUT_ERROR);
      return; // and do not allocate the pairlists below
    }
    if (b_group2_center_only) {
      pairlist = new bool[group1->size()];
    }
    else {
      pairlist = new bool[group1->size() * group2->size()];
    }
  }

  init_scalar_boundaries(0.0, b_group2_center_only ? group1->size() :
                         group1->size() * group2->size());
}


colvar::coordnum::~coordnum()
{
  if (pairlist != NULL) {
    delete [] pairlist;
  }
}


template<int flags> void colvar::coordnum::main_loop(bool **pairlist_elem)
{
  if (b_group2_center_only) {
    cvm::atom group2_com_atom;
    group2_com_atom.pos = group2->center_of_mass();
    for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
      x.real_value += switching_function<flags>(r0, r0_vec, en, ed,
                                                *ai1, group2_com_atom,
                                                pairlist_elem,
                                                tolerance);
    }
    if (b_group2_center_only) {
      group2->set_weighted_gradient(group2_com_atom.grad);
    }
  } else {
    for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
      for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
        x.real_value += switching_function<flags>(r0, r0_vec, en, ed,
                                                  *ai1, *ai2,
                                                  pairlist_elem,
                                                  tolerance);
      }
    }
  }
}


template<int compute_flags> int colvar::coordnum::compute_coordnum()
{
  bool const use_pairlist = (pairlist != NULL);
  bool const rebuild_pairlist = (pairlist != NULL) &&
    (cvm::step_relative() % pairlist_freq == 0);

  bool *pairlist_elem = use_pairlist ? pairlist : NULL;

  if (b_anisotropic) {

    if (use_pairlist) {
      if (rebuild_pairlist) {
        int const flags = compute_flags | ef_anisotropic | ef_use_pairlist |
          ef_rebuild_pairlist;
        main_loop<flags>(&pairlist_elem);
      } else {
        int const flags = compute_flags | ef_anisotropic | ef_use_pairlist;
        main_loop<flags>(&pairlist_elem);
      }

    } else {

      int const flags = compute_flags | ef_anisotropic;
      main_loop<flags>(NULL);
    }

  } else {

    if (use_pairlist) {

      if (rebuild_pairlist) {
        int const flags = compute_flags | ef_use_pairlist | ef_rebuild_pairlist;
        main_loop<flags>(&pairlist_elem);
      } else {
        int const flags = compute_flags | ef_use_pairlist;
        main_loop<flags>(&pairlist_elem);
      }

    } else {

      int const flags = compute_flags;
      main_loop<flags>(NULL);
    }
  }

  return COLVARS_OK;
}


void colvar::coordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {
    compute_coordnum<ef_gradients>();
  } else {
    compute_coordnum<ef_null>();
  }
}


void colvar::coordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
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
  init_scalar_boundaries(0.0, 1.0);

  colvarproxy *proxy = cvm::main()->proxy;

  int a_num = -1, d_num = -1;
  get_keyval(conf, "acceptor", a_num, a_num);
  get_keyval(conf, "donor",    d_num, a_num);

  if ( (a_num == -1) || (d_num == -1) ) {
    cvm::error("Error: either acceptor or donor undefined.\n");
    return;
  }

  cvm::atom acceptor = cvm::atom(a_num);
  cvm::atom donor    = cvm::atom(d_num);
  register_atom_group(new cvm::atom_group);
  atom_groups[0]->add_atom(acceptor);
  atom_groups[0]->add_atom(donor);

  get_keyval(conf, "cutoff",   r0, (3.3 * proxy->angstrom_value));
  get_keyval(conf, "expNumer", en, 6);
  get_keyval(conf, "expDenom", ed, 8);

  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
               INPUT_ERROR);
  }

  if ( (en <= 0) || (ed <= 0) ) {
    cvm::error("Error: negative exponent(s) provided.\n",
               INPUT_ERROR);
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
  init_scalar_boundaries(0.0, 1.0);

  register_atom_group(new cvm::atom_group);
  atom_groups[0]->add_atom(acceptor);
  atom_groups[0]->add_atom(donor);
}


void colvar::h_bond::calc_value()
{
  int const flags = coordnum::ef_null;
  cvm::rvector const r0_vec(0.0); // TODO enable the flag?
  x.real_value =
    coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                        (*atom_groups[0])[0],
                                        (*atom_groups[0])[1],
                                        NULL, 0.0);
}


void colvar::h_bond::calc_gradients()
{
  int const flags = coordnum::ef_gradients;
  cvm::rvector const r0_vec(0.0); // TODO enable the flag?
  coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                      (*atom_groups[0])[0],
                                      (*atom_groups[0])[1],
                                      NULL, 0.0);
}


void colvar::h_bond::apply_force(colvarvalue const &force)
{
  (atom_groups[0])->apply_colvar_force(force);
}


simple_scalar_dist_functions(h_bond)



colvar::selfcoordnum::selfcoordnum(std::string const &conf)
  : cvc(conf), pairlist(NULL)
{
  function_type = "selfcoordnum";
  x.type(colvarvalue::type_scalar);

  colvarproxy *proxy = cvm::main()->proxy;

  group1 = parse_group(conf, "group1");

  get_keyval(conf, "cutoff", r0, cvm::real(4.0 * proxy->angstrom_value));
  get_keyval(conf, "expNumer", en, 6);
  get_keyval(conf, "expDenom", ed, 12);


  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
               INPUT_ERROR);
  }

  if ( (en <= 0) || (ed <= 0) ) {
    cvm::error("Error: negative exponent(s) provided.\n",
               INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

  get_keyval(conf, "tolerance", tolerance, 0.0);
  if (tolerance > 0) {
    get_keyval(conf, "pairListFrequency", pairlist_freq, 100);
    if ( ! (pairlist_freq > 0) ) {
      cvm::error("Error: non-positive pairlistfrequency provided.\n",
                 INPUT_ERROR);
      return;
    }
    pairlist = new bool[(group1->size()-1) * (group1->size()-1)];
  }

  init_scalar_boundaries(0.0, (group1->size()-1) * (group1->size()-1));
}


colvar::selfcoordnum::~selfcoordnum()
{
  if (pairlist != NULL) {
    delete [] pairlist;
  }
}


template<int compute_flags> int colvar::selfcoordnum::compute_selfcoordnum()
{
  cvm::rvector const r0_vec(0.0); // TODO enable the flag?

  bool const use_pairlist = (pairlist != NULL);
  bool const rebuild_pairlist = (pairlist != NULL) &&
    (cvm::step_relative() % pairlist_freq == 0);

  bool *pairlist_elem = use_pairlist ? pairlist : NULL;
  size_t i = 0, j = 0;
  size_t const n = group1->size();

  // Always isotropic (TODO: enable the ellipsoid?)

  if (use_pairlist) {

    if (rebuild_pairlist) {
      int const flags = compute_flags | coordnum::ef_use_pairlist |
        coordnum::ef_rebuild_pairlist;
      for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
          x.real_value +=
            coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                (*group1)[i],
                                                (*group1)[j],
                                                &pairlist_elem,
                                                tolerance);
        }
      }
    } else {
      int const flags = compute_flags | coordnum::ef_use_pairlist;
      for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
          x.real_value +=
            coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                (*group1)[i],
                                                (*group1)[j],
                                                &pairlist_elem,
                                                tolerance);
        }
      }
    }

  } else { // if (use_pairlist) {

    int const flags = compute_flags | coordnum::ef_null;
    for (i = 0; i < n - 1; i++) {
      for (j = i + 1; j < n; j++) {
        x.real_value +=
          coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                              (*group1)[i],
                                              (*group1)[j],
                                              &pairlist_elem,
                                              tolerance);
      }
    }
  }

  return COLVARS_OK;
}


void colvar::selfcoordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {
    compute_selfcoordnum<coordnum::ef_gradients>();
  } else {
    compute_selfcoordnum<coordnum::ef_null>();
  }
}


void colvar::selfcoordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
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
  init_scalar_boundaries(0.0, 1.0);

  colvarproxy *proxy = cvm::main()->proxy;

  // group1 and group2 are already initialized by distance()
  if (group1->b_dummy || group2->b_dummy) {
    cvm::error("Error: neither group can be a dummy atom\n");
    return;
  }

  bool const b_scale = get_keyval(conf, "cutoff", r0,
                                  cvm::real(4.0 * proxy->angstrom_value));

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

  get_keyval(conf, "expNumer", en, 6);
  get_keyval(conf, "expDenom", ed, 12);

  if ( (en%2) || (ed%2) ) {
    cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
               INPUT_ERROR);
  }

  if ( (en <= 0) || (ed <= 0) ) {
    cvm::error("Error: negative exponent(s) provided.\n",
               INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

}


void colvar::groupcoordnum::calc_value()
{
  // create fake atoms to hold the com coordinates
  cvm::atom group1_com_atom;
  cvm::atom group2_com_atom;
  group1_com_atom.pos = group1->center_of_mass();
  group2_com_atom.pos = group2->center_of_mass();
  if (b_anisotropic) {
    int const flags = coordnum::ef_anisotropic;
    x.real_value = coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                       group1_com_atom,
                                                       group2_com_atom,
                                                       NULL, 0.0);
  } else {
    int const flags = coordnum::ef_null;
    x.real_value = coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                       group1_com_atom,
                                                       group2_com_atom,
                                                       NULL, 0.0);
  }
}


void colvar::groupcoordnum::calc_gradients()
{
  cvm::atom group1_com_atom;
  cvm::atom group2_com_atom;
  group1_com_atom.pos = group1->center_of_mass();
  group2_com_atom.pos = group2->center_of_mass();

  if (b_anisotropic) {
    int const flags = coordnum::ef_gradients | coordnum::ef_anisotropic;
    coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                        group1_com_atom,
                                        group2_com_atom,
                                        NULL, 0.0);
  } else {
    int const flags = coordnum::ef_gradients;
    coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                        group1_com_atom,
                                        group2_com_atom,
                                        NULL, 0.0);
  }

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

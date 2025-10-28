// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"


colvar::alpha_angles::alpha_angles()
{
  set_function_type("alpha");
  enable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_scalar);
  colvarproxy *proxy = cvm::main()->proxy;
  r0 = proxy->angstrom_to_internal(3.3);
}


int colvar::alpha_angles::init(std::string const &conf)
{
  int error_code = cvc::init(conf);
  if (error_code != COLVARS_OK) return error_code;

  std::string segment_id;
  std::vector<int> residues;

  bool b_use_index_groups = false;
  cvm::atom_group group_CA, group_N, group_O;

  std::string residues_conf = "";
  std::string prefix;

  // residueRange is mandatory for the topology-based case
  if (key_lookup(conf, "residueRange", &residues_conf)) {
    if (residues_conf.size()) {
      std::istringstream is(residues_conf);
      int initial, final;
      char dash;
      if ( (is >> initial) && (initial > 0) &&
          (is >> dash) && (dash == '-') &&
          (is >> final) && (final > 0) ) {
        for (int rnum = initial; rnum <= final; rnum++) {
          residues.push_back(rnum);
        }
      }
    } else {
      return cvm::error("Error: no residues defined in \"residueRange\".\n", COLVARS_INPUT_ERROR);
    }

    if (residues.size() < 5) {
      return cvm::error("Error: not enough residues defined in \"residueRange\".\n", COLVARS_INPUT_ERROR);
    }
    get_keyval(conf, "psfSegID", segment_id, std::string("MAIN"));

  } else {
    b_use_index_groups = true;
    get_keyval(conf, "prefix", prefix, "alpha_");

    // Not all groups are mandatory, parse silently
    {
      // These lines must be in its own scope to ensure RAII
      auto modify_group_CA = group_CA.get_atom_modifier();
      auto modify_group_N = group_N.get_atom_modifier();
      auto modify_group_O = group_O.get_atom_modifier();
      modify_group_CA.add_index_group(prefix + "CA", true);
      modify_group_N.add_index_group(prefix + "N", true);
      modify_group_O.add_index_group(prefix + "O", true);
    }

    int na = group_CA.size();
    int nn = group_N.size();
    int no = group_O.size();
    if ((nn != 0 || no != 0) && (nn != no)) {
      return cvm::error("Error: If either is provided, atom groups " + prefix + "N and " + prefix + "O must have the same number of atoms.",
                        COLVARS_INPUT_ERROR);
    }
    if (nn != 0 && na != 0 && nn != na) {
      return cvm::error("Error: If both are provided, atom groups " + prefix + "N and " + prefix + "CA must have the same number of atoms.",
                        COLVARS_INPUT_ERROR);
    }
  }

  std::string const &sid    = segment_id;
  std::vector<int> const &r = residues;


  get_keyval(conf, "hBondCoeff", hb_coeff, hb_coeff);
  if ((hb_coeff < 0.0) || (hb_coeff > 1.0)) {
    return cvm::error("Error: hBondCoeff must be defined between 0 and 1.\n", COLVARS_INPUT_ERROR);
  }


  get_keyval(conf, "angleRef", theta_ref, theta_ref);
  get_keyval(conf, "angleTol", theta_tol, theta_tol);

  if (hb_coeff < 1.0) {
    if (b_use_index_groups) {
      if (group_CA.size() < 5) {
        return cvm::error("Not enough atoms (" + cvm::to_str(group_CA.size()) + ") in index group \"" + prefix + "CA\"",
                          COLVARS_INPUT_ERROR);
      }
      for (size_t i = 0; i < group_CA.size()-2; i++) {
        // Note: the angle constructor constructs copies of the atom objects
        theta.push_back(new colvar::angle(group_CA[i],
                                          group_CA[i+1],
                                          group_CA[i+2]));
        register_atom_group(theta.back()->atom_groups[0]);
        register_atom_group(theta.back()->atom_groups[1]);
        register_atom_group(theta.back()->atom_groups[2]);
      }
    } else {
      colvarproxy* const p = cvm::main()->proxy;
      for (size_t i = 0; i < residues.size()-2; i++) {
        theta.push_back(
          new colvar::angle(
            cvm::atom_group::init_atom_from_proxy(p, r[i  ], "CA", sid),
            cvm::atom_group::init_atom_from_proxy(p, r[i+1], "CA", sid),
            cvm::atom_group::init_atom_from_proxy(p, r[i+2], "CA", sid)));
        register_atom_group(theta.back()->atom_groups[0]);
        register_atom_group(theta.back()->atom_groups[1]);
        register_atom_group(theta.back()->atom_groups[2]);
      }
    }

  } else {
    cvm::log("The hBondCoeff specified will disable the Calpha-Calpha-Calpha angle terms.\n");
  }

  {
    get_keyval(conf, "hBondCutoff",   r0, r0);
    get_keyval(conf, "hBondExpNumer", en, en);
    get_keyval(conf, "hBondExpDenom", ed, ed);

    if (hb_coeff > 0.0) {
      colvarproxy* const p = cvm::main()->proxy;
      if (b_use_index_groups) {
        if (group_N.size() < 5) {
          return cvm::error("Not enough atoms (" + cvm::to_str(group_N.size()) + ") in index group \"" + prefix + "N\"",
                            COLVARS_INPUT_ERROR);
        }
        for (size_t i = 0; i < group_N.size()-4; i++) {
          // Note: we need to call the atom copy constructor here because
          // the h_bond constructor does not make copies of the provided atoms
          hb.push_back(
            new colvar::h_bond(cvm::atom_group::init_atom_from_proxy(p,group_O[i]),
                               cvm::atom_group::init_atom_from_proxy(p,group_N[i+4]),
                               r0, en, ed));
          register_atom_group(hb.back()->atom_groups[0]);
        }
      } else {
        for (size_t i = 0; i < residues.size()-4; i++) {
          hb.push_back(
            new colvar::h_bond(cvm::atom_group::init_atom_from_proxy(p,r[i  ], "O",  sid),
                               cvm::atom_group::init_atom_from_proxy(p,r[i+4], "N",  sid),
                               r0, en, ed));
          register_atom_group(hb.back()->atom_groups[0]);
        }
      }
    } else {
      cvm::log("The hBondCoeff specified will disable the hydrogen bond terms.\n");
    }
  }

  return error_code;
}


colvar::alpha_angles::~alpha_angles()
{
  while (theta.size() != 0) {
    delete theta.back();
    theta.pop_back();
  }
  while (hb.size() != 0) {
    delete hb.back();
    hb.pop_back();
  }
  // Our references to atom groups have become invalid now that children cvcs are deleted
  atom_groups.clear();
}


void colvar::alpha_angles::calc_value()
{
  x.real_value = 0.0;

  if (theta.size()) {

    cvm::real const theta_norm =
      (1.0-hb_coeff) / cvm::real(theta.size());

    for (size_t i = 0; i < theta.size(); i++) {

      (theta[i])->calc_value();

      cvm::real const t = ((theta[i])->value().real_value-theta_ref)/theta_tol;
      cvm::real const f = ( (1.0 - (t*t)) /
                            (1.0 - (t*t*t*t)) );

      x.real_value += theta_norm * f;

      if (cvm::debug())
        cvm::log("Calpha-Calpha angle no. "+cvm::to_str(i+1)+" in \""+
                  this->name+"\" has a value of "+
                  (cvm::to_str((theta[i])->value().real_value))+
                  " degrees, f = "+cvm::to_str(f)+".\n");
    }
  }

  if (hb.size()) {

    cvm::real const hb_norm =
      hb_coeff / cvm::real(hb.size());

    for (size_t i = 0; i < hb.size(); i++) {
      (hb[i])->calc_value();
      x.real_value += hb_norm * (hb[i])->value().real_value;
      if (cvm::debug())
        cvm::log("Hydrogen bond no. "+cvm::to_str(i+1)+" in \""+
                  this->name+"\" has a value of "+
                  (cvm::to_str((hb[i])->value().real_value))+".\n");
    }
  }
}


void colvar::alpha_angles::calc_gradients()
{
  size_t i;
  for (i = 0; i < theta.size(); i++)
    (theta[i])->calc_gradients();

  for (i = 0; i < hb.size(); i++)
    (hb[i])->calc_gradients();
}


void colvar::alpha_angles::collect_gradients(std::vector<int> const &atom_ids, std::vector<cvm::rvector> &atomic_gradients)
{
  cvm::real cvc_coeff = sup_coeff * cvm::real(sup_np) * cvm::integer_power(value().real_value, sup_np-1);

  if (theta.size()) {
    cvm::real const theta_norm = (1.0-hb_coeff) / cvm::real(theta.size());

    for (size_t i = 0; i < theta.size(); i++) {
      cvm::real const t = ((theta[i])->value().real_value-theta_ref)/theta_tol;
      cvm::real const f = ( (1.0 - (t*t)) /
                            (1.0 - (t*t*t*t)) );
      cvm::real const dfdt =
        1.0/(1.0 - (t*t*t*t)) *
        ( (-2.0 * t) + (-1.0*f)*(-4.0 * (t*t*t)) );

      // Coefficient of this CVC's gradient in the colvar gradient, times coefficient of this
      // angle's gradient in the CVC's gradient
      cvm::real const coeff = cvc_coeff * theta_norm * dfdt * (1.0/theta_tol);

      for (size_t j = 0; j < theta[i]->atom_groups.size(); j++) {
        auto &ag = *(theta[i]->atom_groups[j]);
        for (size_t k = 0; k < ag.size(); k++) {
          size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                      ag.id(k)) - atom_ids.begin();
          atomic_gradients[a] += coeff * cvm::rvector(ag.grad_x(k), ag.grad_y(k), ag.grad_z(k));
        }
      }
    }
  }

  if (hb.size()) {

    cvm::real const hb_norm = hb_coeff / cvm::real(hb.size());

    for (size_t i = 0; i < hb.size(); i++) {
      // Coefficient of this CVC's gradient in the colvar gradient, times coefficient of this
      // hbond's gradient in the CVC's gradient
      cvm::real const coeff = cvc_coeff * 0.5 * hb_norm;

      for (size_t j = 0; j < hb[i]->atom_groups.size(); j++) {
        auto &ag = *(hb[i]->atom_groups[j]);
        for (size_t k = 0; k < ag.size(); k++) {
          size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                      ag.id(k)) - atom_ids.begin();
          atomic_gradients[a] += coeff * cvm::rvector(ag.grad_x(k), ag.grad_y(k), ag.grad_z(k));
        }
      }
    }
  }
}


void colvar::alpha_angles::apply_force(colvarvalue const &force)
{

  if (theta.size()) {

    cvm::real const theta_norm =
      (1.0-hb_coeff) / cvm::real(theta.size());

    for (size_t i = 0; i < theta.size(); i++) {

      cvm::real const t = ((theta[i])->value().real_value-theta_ref)/theta_tol;
      cvm::real const f = ( (1.0 - (t*t)) /
                            (1.0 - (t*t*t*t)) );

      cvm::real const dfdt =
        1.0/(1.0 - (t*t*t*t)) *
        ( (-2.0 * t) + (-1.0*f)*(-4.0 * (t*t*t)) );

      (theta[i])->apply_force(theta_norm *
                               dfdt * (1.0/theta_tol) *
                               force.real_value );
    }
  }

  if (hb.size()) {

    cvm::real const hb_norm =
      hb_coeff / cvm::real(hb.size());

    for (size_t i = 0; i < hb.size(); i++) {
      (hb[i])->apply_force(0.5 * hb_norm * force.real_value);
    }
  }
}




//////////////////////////////////////////////////////////////////////
// dihedral principal component
//////////////////////////////////////////////////////////////////////

colvar::dihedPC::dihedPC()
{
  set_function_type("dihedPC");
  // Supported through references to atom groups of children cvcs
  enable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_scalar);
}


int colvar::dihedPC::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  if (cvm::debug())
    cvm::log("Initializing dihedral PC object.\n");

  bool b_use_index_groups = false;
  std::string segment_id;
  std::vector<int> residues;
  size_t n_residues;
  std::string residues_conf = "";
  std::string prefix;
  cvm::atom_group group_CA, group_N, group_C;

  // residueRange is mandatory for the topology-based case
  if (key_lookup(conf, "residueRange", &residues_conf)) {
    if (residues_conf.size()) {
      std::istringstream is(residues_conf);
      int initial, final;
      char dash;
      if ( (is >> initial) && (initial > 0) &&
          (is >> dash) && (dash == '-') &&
          (is >> final) && (final > 0) ) {
        for (int rnum = initial; rnum <= final; rnum++) {
          residues.push_back(rnum);
        }
      }
    } else {
      return cvm::error("Error: no residues defined in \"residueRange\".\n", COLVARS_INPUT_ERROR);
    }
    n_residues = residues.size();
    get_keyval(conf, "psfSegID", segment_id, std::string("MAIN"));

  } else {

    b_use_index_groups = true;
    get_keyval(conf, "prefix", prefix, "dihed_");

    // All three groups are required
    {
      auto modify_group_CA = group_CA.get_atom_modifier();
      auto modify_group_N = group_N.get_atom_modifier();
      auto modify_group_C = group_C.get_atom_modifier();
      modify_group_CA.add_index_group(prefix + "CA");
      modify_group_N.add_index_group(prefix + "N");
      modify_group_C.add_index_group(prefix + "C");
    }
    int na = group_CA.size();
    int nn = group_N.size();
    int nc = group_C.size();
    if ((nn != na || na != nc)) {
      return cvm::error("Error: atom groups " + prefix + "N, " + prefix + "CA, and " + prefix +
                        "C must have the same number of atoms.", COLVARS_INPUT_ERROR);
    }
    n_residues = nn;
  }
  if (n_residues < 2) {
    error_code |=
      cvm::error("Error: dihedralPC requires at least two residues.\n", COLVARS_INPUT_ERROR);
  }

  std::string const &sid    = segment_id;
  std::vector<int> const &r = residues;

  std::string vecFileName;
  if (get_keyval(conf, "vectorFile", vecFileName, vecFileName)) {
    int vecNumber;
    get_keyval(conf, "vectorNumber", vecNumber, 0);
    if (vecNumber < 1) {
      error_code |=
          cvm::error("A positive value of vectorNumber is required.", COLVARS_INPUT_ERROR);
    }

    std::istream &vecFile =
      cvm::main()->proxy->input_stream(vecFileName,
                                       "dihedral PCA vector file");
    if (!vecFile) {
      return COLVARS_INPUT_ERROR;
    }

    // TODO: adapt to different formats by setting this flag
    // bool eigenvectors_as_columns = true;
    // if (eigenvectors_as_columns) {
      // Carma-style dPCA file
      std::string line;
      cvm::real c;
      while (vecFile.good()) {
        getline(vecFile, line);
        if (line.length() < 2) break;
        std::istringstream ls(line);
        for (int i=0; i<vecNumber; i++) ls >> c;
        coeffs.push_back(c);
      }
    /* } else { // Uncomment this when different formats are recognized
      // Eigenvectors as lines
      // Skip to the right line
      for (int i = 1; i<vecNumber; i++)
        vecFile.ignore(999999, '\n');

      if (!vecFile.good()) {
        cvm::error("Error reading dihedral PCA vector file " + vecFileName);
      }

      std::string line;
      getline(vecFile, line);
      std::istringstream ls(line);
      cvm::real c;
      while (ls.good()) {
        ls >> c;
        coeffs.push_back(c);
      }
    }
 */
    cvm::main()->proxy->close_input_stream(vecFileName);

  } else {
    get_keyval(conf, "vector", coeffs, coeffs);
  }

  if ( coeffs.size() != 4 * (n_residues - 1)) {
    error_code |= cvm::error("Error: wrong number of coefficients: " + cvm::to_str(coeffs.size()) +
                             ". Expected " + cvm::to_str(4 * (n_residues - 1)) +
                             " (4 coeffs per residue, minus one residue).\n",
                             COLVARS_INPUT_ERROR);
  }
  colvarproxy* const p = cvm::main()->proxy;
  for (size_t i = 0; i < n_residues-1; i++) {
    // Psi
    if (b_use_index_groups) {
      theta.push_back(new colvar::dihedral( group_N[i],
                                            group_CA[i],
                                            group_C[i],
                                            group_N[i+1]));
    } else {
      theta.push_back(
        new colvar::dihedral(
          cvm::atom_group::init_atom_from_proxy(p,r[i  ], "N", sid),
          cvm::atom_group::init_atom_from_proxy(p,r[i  ], "CA", sid),
          cvm::atom_group::init_atom_from_proxy(p,r[i  ], "C", sid),
          cvm::atom_group::init_atom_from_proxy(p,r[i+1], "N", sid)));
    }
    if (cvm::get_error()) {
      return cvm::get_error();
    }
    register_atom_group(theta.back()->atom_groups[0]);
    register_atom_group(theta.back()->atom_groups[1]);
    register_atom_group(theta.back()->atom_groups[2]);
    register_atom_group(theta.back()->atom_groups[3]);
    // Phi (next res)
    if (b_use_index_groups) {
      theta.push_back(new colvar::dihedral(group_C[i],
                                           group_N[i+1],
                                           group_CA[i+1],
                                           group_C[i+1]));
    } else {
      theta.push_back(
        new colvar::dihedral(cvm::atom_group::init_atom_from_proxy(p,r[i  ], "C", sid),
                             cvm::atom_group::init_atom_from_proxy(p,r[i+1], "N", sid),
                             cvm::atom_group::init_atom_from_proxy(p,r[i+1], "CA", sid),
                             cvm::atom_group::init_atom_from_proxy(p,r[i+1], "C", sid)));
    }
    if (cvm::get_error()) {
      return cvm::get_error();
    }
    register_atom_group(theta.back()->atom_groups[0]);
    register_atom_group(theta.back()->atom_groups[1]);
    register_atom_group(theta.back()->atom_groups[2]);
    register_atom_group(theta.back()->atom_groups[3]);
  }

  if (cvm::debug())
    cvm::log("Done initializing dihedPC object.\n");

  return error_code;
}


colvar::dihedPC::~dihedPC()
{
  while (theta.size() != 0) {
    delete theta.back();
    theta.pop_back();
  }
  // Our references to atom groups have become invalid now that children cvcs are deleted
  atom_groups.clear();
}


void colvar::dihedPC::calc_value()
{
  x.real_value = 0.0;
  for (size_t i = 0; i < theta.size(); i++) {
    theta[i]->calc_value();
    cvm::real const t = (PI / 180.) * theta[i]->value().real_value;
    x.real_value += coeffs[2*i  ] * cvm::cos(t)
                  + coeffs[2*i+1] * cvm::sin(t);
  }
}


void colvar::dihedPC::calc_gradients()
{
  for (size_t i = 0; i < theta.size(); i++) {
    theta[i]->calc_gradients();
  }
}


void colvar::dihedPC::collect_gradients(std::vector<int> const &atom_ids, std::vector<cvm::rvector> &atomic_gradients)
{
  cvm::real cvc_coeff = sup_coeff * cvm::real(sup_np) * cvm::integer_power(value().real_value, sup_np-1);
  for (size_t i = 0; i < theta.size(); i++) {
    cvm::real const t = (PI / 180.) * theta[i]->value().real_value;
    cvm::real const dcosdt = - (PI / 180.) * cvm::sin(t);
    cvm::real const dsindt =   (PI / 180.) * cvm::cos(t);
    // Coefficient of this dihedPC's gradient in the colvar gradient, times coefficient of this
    // dihedral's gradient in the dihedPC's gradient
    cvm::real const coeff = cvc_coeff * (coeffs[2*i] * dcosdt + coeffs[2*i+1] * dsindt);

    for (size_t j = 0; j < theta[i]->atom_groups.size(); j++) {
      auto &ag = *(theta[i]->atom_groups[j]);
      for (size_t k = 0; k < ag.size(); k++) {
        size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                    ag.id(k)) - atom_ids.begin();
        atomic_gradients[a] += coeff * cvm::rvector(ag.grad_x(k), ag.grad_y(k), ag.grad_z(k));
      }
    }
  }
}


void colvar::dihedPC::apply_force(colvarvalue const &force)
{
  for (size_t i = 0; i < theta.size(); i++) {
    cvm::real const t = (PI / 180.) * theta[i]->value().real_value;
    cvm::real const dcosdt = - (PI / 180.) * cvm::sin(t);
    cvm::real const dsindt =   (PI / 180.) * cvm::cos(t);

    theta[i]->apply_force((coeffs[2*i  ] * dcosdt +
                           coeffs[2*i+1] * dsindt) * force);
  }
}

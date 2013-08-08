#include "WeakEquationSchrodinger.h"
#include "Material.h"
#include <string>
#include <iostream>
#include <fstream>

namespace ATC {

//==============================================================
//  Class WeakEquationSchrodinger
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationSchrodinger::WeakEquationSchrodinger()
  : WeakEquation(EIGENVALUE_PDE,ELECTRON_WAVEFUNCTION,1) // w=0 inhomo soln
{}
//--------------------------------------------------------------
//  Destructor
//---------------------------------------------------------------------
WeakEquationSchrodinger::~WeakEquationSchrodinger(void)
{}

//---------------------------------------------------------------------
void WeakEquationSchrodinger::BB_tangent_coefficients(
  const FieldName field,
  const FIELD_MATS & fields,
  const Material* material,
  DENS_MAT &coefs) const
{
  material->inv_effective_mass(fields, coefs);// scaled by 1/2 hbar^2
}
//---------------------------------------------------------------------

void WeakEquationSchrodinger::NN_tangent_coefficients(
  const FieldName field,
  const FIELD_MATS & fields,
  const Material* material,
  DENS_MAT & V) const
{
  material->band_edge_potential(fields,V);
  FIELD_MATS::const_iterator phiField = fields.find(ELECTRIC_POTENTIAL);
  const DENS_MAT &  phi = phiField->second;
  V -= phi; // phi in volts equals |e|*phi in [eV]'s
}

};

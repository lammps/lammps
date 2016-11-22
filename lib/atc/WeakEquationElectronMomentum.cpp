#include "WeakEquationElectronMomentum.h"
#include "Material.h"
#include "LammpsInterface.h"

namespace ATC {

//==============================================================
//  Class WeakEquationElectronMomentum
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationElectronMomentum::WeakEquationElectronMomentum()
  : WeakEquation(STATIC_PDE,ELECTRON_VELOCITY,3)
{}

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationElectronMomentum::~WeakEquationElectronMomentum()
{}

void WeakEquationElectronMomentum::convection(const FIELD_MATS &fields,
                                              const Material * material,
                                              DENS_MAT_VEC & flux) const
{
  // set up mass density
  FIELD_MATS::const_iterator nField = fields.find(ELECTRON_DENSITY);
  const DENS_MAT & n = nField->second;
  DENS_MAT nMe(n.nRows(),n.nCols()); 
  material->inv_effective_mass(fields,nMe);
  nMe = n.div_by_element(nMe);

  // set up velocity and flux
  FIELD_MATS::const_iterator vField = fields.find(ELECTRON_VELOCITY);
  const DENS_MAT & velocity = vField->second;
  const CLON_VEC u(velocity,CLONE_COL,0);
  const CLON_VEC v(velocity,CLONE_COL,1);
  const CLON_VEC w(velocity,CLONE_COL,2);
  flux[0] = velocity; 
  flux[1] = velocity;
  flux[2] = velocity;
  CLON_VEC nuu(flux[0],CLONE_COL,0);
  CLON_VEC nuv(flux[1],CLONE_COL,0);
  CLON_VEC nuw(flux[2],CLONE_COL,0);
  CLON_VEC nvu(flux[0],CLONE_COL,1);
  CLON_VEC nvv(flux[1],CLONE_COL,1);
  CLON_VEC nvw(flux[2],CLONE_COL,1);
  CLON_VEC nwu(flux[0],CLONE_COL,2);
  CLON_VEC nwv(flux[1],CLONE_COL,2);
  CLON_VEC nww(flux[2],CLONE_COL,2);
  
  for (int i = 0; i < n.nRows(); i++) {
    // tensor product of velocities
    nuu(i) *= nMe(i,0)*u(i);
    nuv(i) *= nMe(i,0)*v(i);
    nuw(i) *= nMe(i,0)*w(i);
    nvu(i) *= nMe(i,0)*u(i);
    nvv(i) *= nMe(i,0)*v(i);
    nvw(i) *= nMe(i,0)*w(i);
    nwu(i) *= nMe(i,0)*u(i);
    nwv(i) *= nMe(i,0)*v(i);
    nww(i) *= nMe(i,0)*w(i);
  }
};

//---------------------------------------------------------------------
//   compute mass density
//---------------------------------------------------------------------
void WeakEquationElectronMomentum::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & density ) const
{
  material->electron_mass_density(fields, density);
}

//--------------------------------------------------------------
void WeakEquationElectronMomentum::B_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT_VEC &flux) const
{
  convection(fields,material,flux);
}

//==============================================================
//  Class WeakEquationElectronMomentumDDM
//==============================================================

//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
WeakEquationElectronMomentumDDM::WeakEquationElectronMomentumDDM()
  : WeakEquationElectronMomentum()
{
  DENS_MAT dummy;
  _dnCp_.reserve(nsd_);
  for (int i = 0; i < nsd_; i++)
    _dnCp_.push_back(dummy);

  _electricForce_.reserve(nsd_);
  for (int i = 0; i < nsd_; i++)
    _electricForce_.push_back(dummy);
}

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
WeakEquationElectronMomentumDDM::~WeakEquationElectronMomentumDDM()
{}

void WeakEquationElectronMomentumDDM::thermal_stress(const FIELD_MATS &fields,
                                                     const GRAD_FIELD_MATS &gradFields,
                                                     const Material * material,
                                                     DENS_MAT &flux) const
{
  GRAD_FIELD_MATS::const_iterator dtField = gradFields.find(ELECTRON_TEMPERATURE);
  const DENS_MAT_VEC & DTe = dtField->second;

  CLON_VEC tsx(flux,CLONE_COL,0);
  CLON_VEC tsy(flux,CLONE_COL,1);
  CLON_VEC tsz(flux,CLONE_COL,2);

  // ith velocity component has thermal stress of
  // d_i n * Cp * Te
  DENS_MAT nCp(DTe[0].nRows(),DTe[0].nCols());  
  material->electron_heat_capacity(fields,nCp);
  nCp *= 2./3.;  // correction to capacity account for convection
  
  tsx += nCp.mult_by_element(DTe[0]);
  tsy += nCp.mult_by_element(DTe[1]);
  tsz += nCp.mult_by_element(DTe[2]);
  
  FIELD_MATS::const_iterator tField = fields.find(ELECTRON_TEMPERATURE);
  const DENS_MAT & Te = tField->second;
  
  material->D_electron_heat_capacity(fields,gradFields,_dnCp_);
  for (int i = 0; i < nsd_; i++)
    _dnCp_[i] *= 2./3.; // correction to capacity account for convection
  tsx += Te.mult_by_element(_dnCp_[0]);
  tsy += Te.mult_by_element(_dnCp_[1]);
  tsz += Te.mult_by_element(_dnCp_[2]);
}

//---------------------------------------------------------------------
//   compute mass density
//---------------------------------------------------------------------
void WeakEquationElectronMomentumDDM::M_integrand(
  const FIELD_MATS &fields,
  const Material * material,
  DENS_MAT & density ) const
{
  material->electron_drag_velocity_coefficient(fields, density);
}

//--------------------------------------------------------------
bool WeakEquationElectronMomentumDDM::N_integrand(
  const FIELD_MATS &fields,
  const GRAD_FIELD_MATS &grad_fields,
  const Material * material,
  DENS_MAT &flux) const
{
  FIELD_MATS::const_iterator vField = fields.find(ELECTRON_VELOCITY);
  const DENS_MAT & velocity = vField->second;
  flux.reset(velocity.nRows(),velocity.nCols());
  thermal_stress(fields, grad_fields, material, flux);
  material->electric_displacement(fields, grad_fields, _electricForce_);


  FIELD_MATS::const_iterator nField = fields.find(ELECTRON_DENSITY);
  const DENS_MAT & n = nField->second;

  
  CLON_VEC tsx(flux,CLONE_COL,0);
  CLON_VEC tsy(flux,CLONE_COL,1);
  CLON_VEC tsz(flux,CLONE_COL,2);
  tsx += n.mult_by_element(_electricForce_[0]);
  tsy += n.mult_by_element(_electricForce_[1]);
  tsz += n.mult_by_element(_electricForce_[2]);
  return true;
}

}; // END namespace

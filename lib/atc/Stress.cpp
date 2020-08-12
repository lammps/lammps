#include "Stress.h"
#include "CauchyBorn.h"
#include "CBLattice.h"
#include "CbLjCut.h" 
#include "CbLjSmoothLinear.h" 
#include "CbEam.h" 
#include "ATC_Error.h"
#include "LammpsInterface.h"
#include "VoigtOperations.h"
#include <iostream>

using ATC_Utility::command_line;
using ATC_Utility::str2dbl;
using voigt3::voigt_idx1;
using voigt3::voigt_idx2;
using voigt3::to_voigt_unsymmetric;
using voigt3::from_voigt_unsymmetric;
using voigt3::to_voigt;
using voigt3::from_voigt;
using std::stringstream;
using std::vector;
using std::string;
using std::fstream;

namespace ATC {
//=============================================================================
// extracts a stress at an integration point 
// Note: Utility function: not in header
//=============================================================================
DENS_MAT extract_stress(const DENS_MAT_VEC &sigma, INDEX ip=0)
{
  DENS_MAT s(3,3,false);
  for (int j=0; j<3; j++) for (int i=0; i<3; i++) s(i,j) = sigma[i](ip,j);
  return s;
}
//=============================================================================
// computes the pressure from the stress at the first quadrature point (in atm)
// Note: Utility function: not in header
//=============================================================================
double compute_pressure(const DENS_MAT_VEC &sigma, const DENS_MAT &F)
{
  // pressure in units (mass-velocity^2)/Volume (LAMMPS real)
  double p = (sigma[0](0,0) + sigma[1](0,1) + sigma[2](0,2)) * (1.0/3.0);
  p *= 1.0e14/6.0221415;            // convert from units real to Pa 
  p *= 1.0/101235.0;                // convert from Pa to ATM
  return p * pow(det(F), -1.0/3.0); // convert from PK2 to Cauchy stress
}
//=============================================================================
// extracts the deformation gradient at a quadrature point, q
// Note: Utility function: not in header
//=============================================================================
void deformation_gradient(const DENS_MAT_VEC &du, INDEX q, MATRIX &F)
{
  F.reset(du.size(), du.size(), false);
  for (INDEX j=0; j<F.nCols(); j++) {
    for (INDEX i=0; i<F.nRows(); i++) F(i,j) = du[j](q,i);
    F(j,j) += 1.0;
  } 
}

//=============================================================================

// E = 1/2 stress*strain for linear elastic models
//=============================================================================
  void Stress::elastic_energy(const FIELD_MATS & /* fields */,
                            const GRAD_FIELD_MATS &gradFields,
                            DENS_MAT &energy) const
{
 int nRows = ( ((gradFields.find(DISPLACEMENT))->second)[0]).nRows();
 energy.reset(nRows,1);
 ATC::LammpsInterface::instance()->print_msg("WARNING: returning dummy elastic energy");
 
}

//=============================================================================
// isotropic linear elastic
//=============================================================================
StressLinearElastic::StressLinearElastic(fstream &fileId) 
  : StressCubicElastic(), E_(0), nu_(0), mu_(0), lambda_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if (line[0] == "end") {
      mu_ = E_/(2.0+2.0*nu_);
      lambda_ = mu_*nu_ / (0.5 - nu_);
      StressCubicElastic::c11_ = E_*(1-nu_)/(1+nu_)/(1-2*nu_);
      StressCubicElastic::c12_ = E_*nu_    /(1+nu_)/(1-2*nu_);
      StressCubicElastic::c44_ = E_/(1+nu_)/2;
      if (nu_ < 0.0 || nu_ > 1.0) 
        throw ATC_Error("bad linear elastic constants");   
      if (lambda_ < 0.0 || mu_ < 0.0)
        throw ATC_Error("bad continuum material parameter");
      return;
    }
    else if (line[0]=="modulus")        E_ = str2dbl(line[1]);
    else if (line[0]=="possions_ratio") nu_ = str2dbl(line[1]);
    else throw ATC_Error( "unrecognized material function");
  }
}
//=============================================================================
// compute the stress at N integration points from the displacement gradient
// T_{ij} = 1/2*C_{ijkl}* (u_{k,l} + u_{l,k})
//=============================================================================
  void StressLinearElastic::stress(const FIELD_MATS      & /* fields */,
                                 const GRAD_FIELD_MATS &gradFields,
                                 DENS_MAT_VEC &sigma)
{
  GRAD_FIELD_MATS::const_iterator du_itr = gradFields.find(DISPLACEMENT);
  const DENS_MAT_VEC &du = du_itr->second;

  CLON_VEC uxx(du[0],CLONE_COL,0);
  CLON_VEC uxy(du[1],CLONE_COL,0);
  CLON_VEC uxz(du[2],CLONE_COL,0);
  CLON_VEC uyx(du[0],CLONE_COL,1);
  CLON_VEC uyy(du[1],CLONE_COL,1);
  CLON_VEC uyz(du[2],CLONE_COL,1);
  CLON_VEC uzx(du[0],CLONE_COL,2);
  CLON_VEC uzy(du[1],CLONE_COL,2);
  CLON_VEC uzz(du[2],CLONE_COL,2);

  const INDEX N = uxx.size();          // # of integration pts
  sigma.assign(3, DENS_MAT(N,3));

  // precompute the pressure and copy to the diagonal
  column(sigma[0],0) = (uxx + uyy + uzz)*(-lambda_);
  column(sigma[1],1) = column(sigma[0],0);
  column(sigma[2],2) = column(sigma[0],0);

  column(sigma[0],0) -= 2.0*mu_*uxx;
  column(sigma[0],1) = (uxy + uyx)*(-mu_);
  column(sigma[0],2) = (uxz + uzx)*(-mu_);
  column(sigma[1],0) = column(sigma[0],1);
  column(sigma[1],1) -= 2.0*mu_*uyy;
  column(sigma[1],2) = (uyz + uzy)*(-mu_);
  column(sigma[2],0) = column(sigma[0],2);
  column(sigma[2],1) = column(sigma[1],2);
  column(sigma[2],2) -= 2.0*mu_*uzz;
}
//=============================================================================
// cubic elastic
//=============================================================================
StressCubicElastic::StressCubicElastic(fstream &fileId) 
  : c11_(0), c12_(0), c44_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if      (line[0]=="end") return;
    else if (line[0]=="c11") c11_ = str2dbl(line[1]);
    else if (line[0]=="c12") c12_ = str2dbl(line[1]);
    else if (line[0]=="c44") c44_ = str2dbl(line[1]);
    else throw ATC_Error( "unrecognized material function"); 
  }
}
//---------------------------------------------------------------------------
// compute the stress at N integration points from the displacement gradient
// T_{ij} = 1/2*C_{ijkl}*(u_{k,l} + u_{l,k}) 
//---------------------------------------------------------------------------
  void StressCubicElastic::stress(const FIELD_MATS      & /* fields */,
                                const GRAD_FIELD_MATS &gradFields,
                                DENS_MAT_VEC  &sigma)  
{
  GRAD_FIELD_MATS::const_iterator du_itr = gradFields.find(DISPLACEMENT);
  const DENS_MAT_VEC &du = du_itr->second;
  CLON_VEC uxx(du[0],CLONE_COL,0);
  CLON_VEC uxy(du[1],CLONE_COL,0);
  CLON_VEC uxz(du[2],CLONE_COL,0);
  CLON_VEC uyx(du[0],CLONE_COL,1);
  CLON_VEC uyy(du[1],CLONE_COL,1);
  CLON_VEC uyz(du[2],CLONE_COL,1);
  CLON_VEC uzx(du[0],CLONE_COL,2);
  CLON_VEC uzy(du[1],CLONE_COL,2);
  CLON_VEC uzz(du[2],CLONE_COL,2);

  const INDEX N = uxx.size();          // # of integration pts
  sigma.assign(3, DENS_MAT(N,3));

  const double c12 = c12_;
  const double c11 = c11_;
  const double c44 = c44_;

  // scaling: stress must return (-) stress
  column(sigma[0],0) = -c11*uxx - c12*(uyy+uzz);
  column(sigma[1],1) = -c11*uyy - c12*(uxx+uzz);
  column(sigma[2],2) = -c11*uzz - c12*(uxx+uyy);
  column(sigma[0],1) = -c44*(uxy+uyx);
  column(sigma[1],0) = column(sigma[0],1);
  column(sigma[0],2) = -c44*(uxz+uzx);
  column(sigma[2],0) = column(sigma[0],2);
  column(sigma[1],2) = -c44*(uyz+uzy);
  column(sigma[2],1) = column(sigma[1],2);

}
//---------------------------------------------------------------------------
// compute the elastic energy at N integration points from displacement gradient
// E = 1/8*C_{ijkl}* (u_{k,l} + u_{l,k})* (u_{i,j} + u_{j,i})*rho ?
//   = 1/2 (4 c44 (u12^2 + u13^2 + u23^2) + 2 c12 (u11 u22 + u11 u33 + u22 u33) 
//        + c11 (u11^2 + u22^2 + u33^2))
//---------------------------------------------------------------------------
  void StressCubicElastic::elastic_energy(const FIELD_MATS      & /* fields */,
                                        const GRAD_FIELD_MATS &gradFields,
                                        DENS_MAT  &energy) const
{
  GRAD_FIELD_MATS::const_iterator du_itr = gradFields.find(DISPLACEMENT);
  const DENS_MAT_VEC &du = du_itr->second;
  CLON_VEC uxx(du[0],CLONE_COL,0);
  CLON_VEC uxy(du[1],CLONE_COL,0);
  CLON_VEC uxz(du[2],CLONE_COL,0);
  CLON_VEC uyx(du[0],CLONE_COL,1);
  CLON_VEC uyy(du[1],CLONE_COL,1);
  CLON_VEC uyz(du[2],CLONE_COL,1);
  CLON_VEC uzx(du[0],CLONE_COL,2);
  CLON_VEC uzy(du[1],CLONE_COL,2);
  CLON_VEC uzz(du[2],CLONE_COL,2);

  CLON_VEC E(energy,CLONE_COL,0);

  const double c12 = c12_;
  const double c11 = c11_;
  const double c44 = c44_;
  //double scale = (ATC::LammpsInterface::instance()->mvv2e());
  for (INDEX gp=0; gp<du.front().nRows(); gp++) {
    double u11 = uxx(gp); 
    double u22 = uyy(gp); 
    double u33 = uzz(gp); 
    double u12 = 0.5*(uxy(gp)+uyx(gp)); 
    double u13 = 0.5*(uxz(gp)+uzx(gp));
    double u23 = 0.5*(uyz(gp)+uzy(gp));
    double EE  = 0.5* (4.0*c44*(u12*u12 + u13*u13 + u23*u23) 
                     + 2.0*c12*(u11*u22 + u11*u33 + u22*u33) 
                         + c11*(u11*u11 + u22*u22 + u33*u33));
    
    E(gp) = EE;
  }
}

void StressCubicElastic::set_tangent(void) 
{
  C_.reset(6,6);
  C_(0,0)=C_(1,1)=C_(2,2)                        =c11_;
  C_(0,1)=C_(1,0)=C_(1,2)=C_(2,1)=C_(0,2)=C_(2,0)=c12_;
  C_(3,3)=C_(4,4)=C_(5,5)                        =c44_;
}

//=============================================================================
// damped cubic elastic
//=============================================================================
StressCubicElasticDamped::StressCubicElasticDamped(fstream &fileId) 
  : StressCubicElastic(), gamma_(0)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file");
  vector<string> line;
  while(fileId.good()) {
    command_line(fileId, line);
    if      (line[0]=="end") return;
    else if (line[0]=="c11") StressCubicElastic::c11_ = str2dbl(line[1]);
    else if (line[0]=="c12") StressCubicElastic::c12_ = str2dbl(line[1]);
    else if (line[0]=="c44") StressCubicElastic::c44_ = str2dbl(line[1]);
    else if (line[0]=="gamma") gamma_ = str2dbl(line[1]);
    else throw ATC_Error( "unrecognized material function");
  }
}
//---------------------------------------------------------------------------
// compute the stress at N integration points 
//---------------------------------------------------------------------------

void StressCubicElasticDamped::stress(const FIELD_MATS      &fields,
                                const GRAD_FIELD_MATS &gradFields,
                                DENS_MAT_VEC  &sigma)  
{
  StressCubicElastic::stress(fields,gradFields,sigma);
  GRAD_FIELD_MATS::const_iterator dv_itr = gradFields.find(VELOCITY);
  const DENS_MAT_VEC &dv = dv_itr->second;
  CLON_VEC vxx(dv[0],CLONE_COL,0);
  CLON_VEC vxy(dv[1],CLONE_COL,0);
  CLON_VEC vxz(dv[2],CLONE_COL,0);
  CLON_VEC vyx(dv[0],CLONE_COL,1);
  CLON_VEC vyy(dv[1],CLONE_COL,1);
  CLON_VEC vyz(dv[2],CLONE_COL,1);
  CLON_VEC vzx(dv[0],CLONE_COL,2);
  CLON_VEC vzy(dv[1],CLONE_COL,2);
  CLON_VEC vzz(dv[2],CLONE_COL,2);

  // scaling: stress must return (-) stress
  column(sigma[0],0) += -gamma_*vxx;
  column(sigma[1],1) += -gamma_*vyy;
  column(sigma[2],2) += -gamma_*vzz;
  column(sigma[0],1) += -0.5*gamma_*(vxy+vyx);
  column(sigma[1],0) += column(sigma[0],1);
  column(sigma[0],2) += -0.5*gamma_*(vxz+vzx);
  column(sigma[2],0) += column(sigma[0],2);
  column(sigma[1],2) += -0.5*gamma_*(vyz+vzy);
  column(sigma[2],1) += column(sigma[1],2);
}
//==============================================================================
// cauchy born model
//==============================================================================
StressCauchyBorn::StressCauchyBorn(fstream &fileId, CbData &cb)
  :  cblattice_(NULL),
     potential_(NULL),
     makeLinear_(false),
     cubicMat_(NULL),
     initialized_(false),
     fixed_temperature_(0.),
     cbdata_(cb)
{
  if (!fileId.is_open()) throw ATC_Error("cannot open material file"); 
  while(fileId.good()) {
    // reads a line from the material file
    vector<string> line;
    command_line(fileId, line);
    if (line.empty()) continue;   // skip blank lines
    else if (line[0]=="end") {
      delete cblattice_;
      if (!potential_) throw ATC_Error( "no potential defined"); 
      cblattice_ = new CBLattice(cbdata_.cell_vectors, cbdata_.basis_vectors); 
      return;
    }
    else if (line[0] == "pair_style") {
      if (line[1] == "lj/cut") {      // Lennard-Jones w/ cutoff radius
        if (line.size()<3) throw(ATC_Error("no lj/cut cutoff radius"));
        const double rc = str2dbl(line[2]);
        while (!fileId.eof()) { // find next pair_coeff command
          command_line(fileId, line); 
          if (line.size() && line[0]=="pair_coeff") break;
        }
        if (line[0] != "pair_coeff" || line.size() != 3) {
          throw(ATC_Error("lj/cut needs 2 coefficients"));
        }
        delete potential_;
        potential_ = new CbLjCut(str2dbl(line[1]), str2dbl(line[2]), rc);
      }
      else if (line[1] == "lj/smooth/linear") { // Lennard-Jones w/ cutoff radius and smoothed
        if (line.size()<3) throw(ATC_Error("no lj/smooth/linear cutoff radius"));
        const double rc = str2dbl(line[2]);
        while (!fileId.eof()) { // find next pair_coeff command
          command_line(fileId, line);
          if (line.size() && line[0]=="pair_coeff") break;
        }
        if (line[0] != "pair_coeff" || line.size() != 3) {
          throw(ATC_Error("lj/smooth/linear needs 2 coefficients"));
        }
        delete potential_;
        potential_ = new CbLjSmoothLinear(str2dbl(line[1]), str2dbl(line[2]), rc);
      }
      else if (line[1] == "eam") {      // Embedded atom method potential 
        delete potential_;
        potential_ = new CbEam();
      }
      else throw (ATC_Error("Invalid pair style"));
    }
    else if (line[0] == "linear") makeLinear_ = true;
    else if (line[0] == "temperature" && line.size() == 2) {
      fixed_temperature_ = str2dbl(line[1]);
    }
    else if (line[0]=="material" || line[0]=="stress") /* ignore this */; 
    else throw ATC_Error( "Unrecognized Cauchy-Born parameter: "+line[0]+".");
  }

}
//==============================================================================
//* default destructor - delete potential and lattice
//==============================================================================
StressCauchyBorn::~StressCauchyBorn()
{
  if (potential_) delete potential_;
  if (cblattice_) delete cblattice_;
  if (cubicMat_) delete cubicMat_;
}

//==============================================================================
// initialize
//==============================================================================
void StressCauchyBorn::initialize(void) 
{
  if (!initialized_) {
    if (makeLinear_) linearize();

    stringstream ss;
    double k = stiffness()*cbdata_.e2mvv; 
    double m = cbdata_.atom_mass;
    double w0 = sqrt(k*m);
    ss << "CB stiffness: " << stiffness() << " Einstein freq: " << w0;
    ATC::LammpsInterface::instance()->print_msg_once(ss.str());

    initialized_ = true;
  }
}

//==============================================================================
// compute the bond stiffness consistent with the einstein freq
//==============================================================================
double StressCauchyBorn::stiffness(void) const
{
  AtomCluster vac;
  cblattice_->atom_cluster(eye<double>(3,3), potential_->cutoff_radius(), vac);
  DENS_MAT k = vac.force_constants(0,potential_);
  return k(0,0);
}
//==============================================================================
// compute the stress at N integration points from the displacement gradient
//==============================================================================
void StressCauchyBorn::stress(const FIELD_MATS &fields, 
                              const GRAD_FIELD_MATS &gradFields, 
                              DENS_MAT_VEC &sigma)
{
  if (cubicMat_) {
    cubicMat_->stress(fields, gradFields, sigma);
    return;
  }
  FIELD_MATS::const_iterator               temp = fields.find(TEMPERATURE);
  GRAD_FIELD_MATS::const_iterator disp_gradient = gradFields.find(DISPLACEMENT);

  // Scaling factor - scale by atomic volume and energy conversion.
  // negative because stress must return (-) stress
  const double fact = -cbdata_.inv_atom_volume * cbdata_.e2mvv;
  const DENS_MAT_VEC &du(disp_gradient->second);
  const INDEX num_integration_pts = du.front().nRows(); 
  const INDEX nsd = du.size();
  DENS_MAT F(nsd,nsd);                // displacement gradient 
  bool temp_varies = (temp != fields.end());
  sigma.assign(nsd, DENS_MAT(num_integration_pts, nsd));
  StressAtIP S(sigma); // wrapper for quadrature points.
  AtomCluster vac;

  for (INDEX gp=0; gp<num_integration_pts; gp++) {
    // Sets the quadrature point to be computed.
    S.set_quadrature_point(gp);
    // Get displacement gradient and construct a virtual atom cluster.
    deformation_gradient(du, gp, F);
    // Generates the atom cluster, given the deformation gradient.
    cblattice_->atom_cluster(F, potential_->cutoff_radius(), vac);
    // Get temperature (assume 0K if no temperature field is present).
    const double T = (temp_varies ? temp->second[gp] : fixed_temperature_);
    // Computes the cauchy-born stresses.
    const StressArgs args(vac, potential_, cbdata_.boltzmann, cbdata_.hbar, T);
    cb_stress(args, S);

    // copy symmetric part of stress and scale by V0
    for (INDEX i=0; i<nsd; i++) { 
      S(i,i) *= fact;  
      for (INDEX j=i+1; j<nsd; j++)  S(j,i)=(S(i,j)*=fact);
    }
  }
}

//==============================================================================
// Computes free (T>0)/potential(T=0) energy density.  [mvv/L^3]
//==============================================================================
void StressCauchyBorn::elastic_energy(const FIELD_MATS &fields, 
                                      const GRAD_FIELD_MATS &gradFields,
                                      DENS_MAT &energy) const
{
  if (cubicMat_) {
    cubicMat_->elastic_energy(fields, gradFields, energy);
    return;
  }
  FIELD_MATS::const_iterator               temp = fields.find(TEMPERATURE);
  GRAD_FIELD_MATS::const_iterator disp_gradient = gradFields.find(DISPLACEMENT);
  const DENS_MAT_VEC &du(disp_gradient->second);
  DENS_MAT F(du.size(),du.size());
  AtomCluster vac;

  for (INDEX gp=0; gp<du.front().nRows(); gp++) {
    deformation_gradient(du, gp, F);
    cblattice_->atom_cluster(F, potential_->cutoff_radius(), vac);
    double T = (temp!=fields.end() ? temp->second[gp] : fixed_temperature_);
    energy[gp] = cb_energy(StressArgs(vac, potential_, cbdata_.boltzmann, cbdata_.hbar, T));
  }
  // Scaling factor - scale by atomic volume and energy conversion.
  
  energy *= cbdata_.inv_atom_volume * cbdata_.e2mvv;
}

//==============================================================================
// Computes entropic energy density.  [mvv/L^3]
//==============================================================================
void StressCauchyBorn::entropic_energy(const FIELD_MATS &fields, 
                                       const GRAD_FIELD_MATS &gradFields,
                                       DENS_MAT &energy) const
{
  FIELD_MATS::const_iterator               temp = fields.find(TEMPERATURE);
  GRAD_FIELD_MATS::const_iterator disp_gradient = gradFields.find(DISPLACEMENT);
  const DENS_MAT_VEC &du(disp_gradient->second);
  DENS_MAT F(du.size(),du.size());
  AtomCluster vac;

  for (INDEX gp=0; gp<du.front().nRows(); gp++) {
    deformation_gradient(du, gp, F);
    cblattice_->atom_cluster(F, potential_->cutoff_radius(), vac);
    double T = (temp!=fields.end() ? temp->second[gp] : fixed_temperature_);
    energy[gp] = cb_entropic_energy(StressArgs(vac, potential_, cbdata_.boltzmann, cbdata_.hbar, T));
  }
  // Scaling factor - scale by atomic volume and energy conversion.
  energy *= cbdata_.inv_atom_volume * cbdata_.e2mvv;
}

//==============================================================================
// creates a linearization for a deformation gradient
//==============================================================================
void StressCauchyBorn::linearize(MATRIX *F)
{
  if (cubicMat_) delete cubicMat_;
  DENS_MAT C;   
  if (F) tangent(*F, C);
  else   tangent(eye<double>(3,3), C);
  cubicMat_ = new StressCubicElastic(C(0,0), C(0,1), C(3,3));

  stringstream ss;
  double c11 = C(0,0)/cbdata_.e2mvv;
  double c12 = C(0,1)/cbdata_.e2mvv;
  double c44 = C(3,3)/cbdata_.e2mvv;
  ss << "created cubic stress function:"
     << "\n   lammps         ATC units" 
     << "\n   c11=" << c11 << " " << C(0,0)
     << "\n   c12=" << c12 << " " << C(0,1)
     << "\n   c44=" << c44 << " " << C(3,3);
  ATC::LammpsInterface::instance()->print_msg_once(ss.str());
}
//==============================================================================
// sets C as the material tangent modulus, given deformation gradient F
//==============================================================================
// Note: C is dS/dC which is 1/2 dS/dF_sym

void StressCauchyBorn::tangent(const MATRIX &F, MATRIX &C)  const
{
  if (cubicMat_) {
    cubicMat_->tangent(F,C);
    return;
  }
  elasticity_tensor(F,C);
}
//==============================================================================
// 1st elasticity tensor : B = dP/dF = C F F + S I ( 9 x 9 in Voigt notation)
// 2nd elasticity tensor : C = dS/dE ( 6 x 6 in Voigt notation)
//==============================================================================
DENS_VEC StressCauchyBorn::elasticity_tensor(const VECTOR &Fv, MATRIX &C, const ElasticityTensorType type)  const
{
   DENS_MAT F;
   if (Fv.nRows()==9) { F = from_voigt_unsymmetric(Fv); }
   else               { F = from_voigt(Fv); }
   return elasticity_tensor(F, C,type);
}
DENS_VEC StressCauchyBorn::elasticity_tensor(const MATRIX &F, MATRIX &C, const ElasticityTensorType type)  const
{
  double T = 0; 
  AtomCluster vac;
  cblattice_->atom_cluster(F, potential_->cutoff_radius(), vac);
  if (vac.size() < 4) throw ATC_Error("StressCauchyBorn::second_elasticity_tensor cluster does not have sufficient atoms");
  const StressArgs args(vac, potential_, cbdata_.boltzmann, cbdata_.hbar, T);
  // if using EAM potential, calculate embedding function and derivatives
  bool hasEAM = potential_->terms.embedding;
  double embed_p  = 0;
  double embed_pp = 0;
  if (hasEAM) {
    double e_density = cb_electron_density(args);
    embed_p  = potential_->F_p(e_density); // "F" in usual EAM symbology
    embed_pp = potential_->F_pp(e_density);
  }


  int size = 6;
  if (type == FIRST_ELASTICITY_TENSOR) { size = 9; }
  DENS_VEC Z(size), S(size), Zfp(size);
  Zfp = 0;
  C.reset(size,size);
  for (INDEX a=0; a<vac.size(); a++) {
    const DENS_VEC &Ra = vac.R(a);
    if (type == FIRST_ELASTICITY_TENSOR) { 
      DENS_VEC ra = F*Ra;
      for (INDEX i=0; i<size; i++) { Z(i)=ra(voigt_idx1[i])*Ra(voigt_idx2[i]); }
    }
    else {
      for (INDEX i=0; i<size; i++) { Z(i)=Ra(voigt_idx1[i])*Ra(voigt_idx2[i]); }
    }
    double d = vac.bond_length(a);
    double rinv   = 1.0/d;
    double phi_r  = potential_->phi_r(d);  // computes phi'
    double phi_rr = potential_->phi_rr(d); // computes phi''
    double fact1 = 0.5*phi_r*rinv;  // 1/2 see Philips
    double fact2 = 0.5*(phi_rr - phi_r*rinv) * rinv*rinv; 
    if (hasEAM) {
      double rho_r  = potential_->rho_r(d);  // computes rho'
      double rho_rr = potential_->rho_rr(d); // computes rho''
      fact1 += embed_p*rho_r*rinv;
      fact2 += embed_p*(rho_rr - rho_r*rinv) * rinv*rinv; 
      Zfp += Z*(rho_r*rinv);
    }
    for (INDEX i=0; i<size; i++) {
      S(i) += fact1*Z(i);
      for (INDEX j=0; j<size; j++) { 
        C(i,j) += fact2*Z(i)*Z(j); 
      }
    }
    if (type == FIRST_ELASTICITY_TENSOR) { 
      for (INDEX i=0; i<9; i++) {
        for (INDEX j=0; j<9; j++) { 
          if ( voigt_idx1[i] == voigt_idx1[j] )  { // \delta_ik S_JL
            C(i,j) += fact1*Ra(voigt_idx2[i])*Ra(voigt_idx2[j]);
          }
        }
      }
    }
  }
  if (hasEAM) {
    for (INDEX i=0; i<6; i++) {
      for (INDEX j=0; j<6; j++) { 
        C(i,j) += embed_pp*Zfp(i)*Zfp(j);  
      }
    }
  }
  double s = cbdata_.inv_atom_volume * cbdata_.e2mvv;  
  S *= s;
  C *= s;
  return S;
}
}// end atc namespace

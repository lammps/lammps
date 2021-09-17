#ifndef CBEAM_H
#define CBEAM_H

#include <iostream>
#include <cstdlib>
#include "CbPotential.h"
#include "LammpsInterface.h"
#include "MANYBODY/pair_eam.h"

namespace ATC
{

  /**
   *  @class  CbEam
   *  @brief  Class for computing Cauchy-Born quantities for an Embeded-Atom Method material
   *          (A factor of one-half is already included to split the
   *           bond energy between atoms)
   */

  class CbEam : public CbPotential
  {
  public:
    //! Constructor
    CbEam(void) : CbPotential(Interactions(PAIRWISE,EAM)) {
      // get pointer to lammps' pair_eam object
      lammps_eam = ATC::LammpsInterface::instance()->pair_eam();
      nrho  = &lammps_eam->nrho;
      nr    = &lammps_eam->nr;
      nfrho = &lammps_eam->nfrho;
      nrhor = &lammps_eam->nrhor;
      nz2r  = &lammps_eam->nz2r;
      type2frho = lammps_eam->type2frho;
      type2rhor = lammps_eam->type2rhor;
      type2z2r  = lammps_eam->type2z2r;
      dr    = &lammps_eam->dr;
      rdr   = &lammps_eam->rdr;
      drho  = &lammps_eam->drho;
      rdrho = &lammps_eam->rdrho;
      rhor_spline = &lammps_eam->rhor_spline;
      frho_spline = &lammps_eam->frho_spline;
      z2r_spline  = &lammps_eam->z2r_spline;
      cutmax = &lammps_eam->cutmax;
    }

    //! Returns the cutoff readius of the EAM potential functions rho and z2r.
    double cutoff_radius() const { return *cutmax; }
    //! Returns the EAM pair energy
    double phi(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*z2r_spline)[type2z2r[1][1]][m];
      double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      return z2/r;
    }
    //! Returns the first derivative of the pair energy.
    double phi_r(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*z2r_spline)[type2z2r[1][1]][m];
      double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      double z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
      return (1.0/r)*(z2p-z2/r);
    }
    //! Returns the second derivative of the pair energy.
    double phi_rr(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*z2r_spline)[type2z2r[1][1]][m];
      double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      double z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
      double z2pp = (*rdr)*(2.0*coeff[0]*p + coeff[1]);
      return (1.0/r)*(z2pp-2.0*z2p/r+2.0*z2/(r*r));
    }
    //! Returns the third derivative of the pair energy.
    double phi_rrr(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*z2r_spline)[type2z2r[1][1]][m];
      double z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      double z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
      double z2pp = (*rdr)*(2.0*coeff[0]*p + coeff[1]);
      double z2ppp = (*rdr)*(*rdr)*2.0*coeff[0];
      return (1.0/r)*(z2ppp-3.0*z2pp/r+6.0*z2p/(r*r)-6.0*z2/(r*r*r));
    }
    //! Returns the EAM atomic charge density.
    double rho(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*rhor_spline)[type2rhor[1][1]][m];
      return (((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]);
    }
   //! Returns the first derivative of the atomic charge density.
    double rho_r(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*rhor_spline)[type2rhor[1][1]][m];
      return ((coeff[0]*p + coeff[1])*p + coeff[2]);
    }
    //! Returns the second derivative of the atomic charge density.
    double rho_rr(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*rhor_spline)[type2rhor[1][1]][m];
      return ((*rdr)*(2.0*coeff[0]*p + coeff[1]));
    }
    //! Returns the third derivative of the atomic charge density.
    double rho_rrr(const double &r) const
    {
      double p = r*(*rdr) + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,(*nr)-1);
      p -= m;
      p = MIN(p,1.0);
      // for now, assume itype = jtype = 1
      double *coeff = (*rhor_spline)[type2rhor[1][1]][m];
      return ((*rdr)*(*rdr)*2.0*coeff[0]);
    }
    //! Returns the EAM embedding energy.
    double F(const double &p) const
    {
      double q = p*(*rdrho) + 1.0;
      int m = static_cast<int> (q);
      m = MIN(m,(*nrho)-1);
      q -= m;
      q = MIN(q,1.0);
      // for now, assume itype = 1
      double *coeff = (*frho_spline)[type2frho[1]][m];
      return (((coeff[3]*q + coeff[4])*q + coeff[5])*q + coeff[6]);
    }
   //! Returns the first derivative of the embedding energy.
    double F_p(const double &p) const
    {
      double q = p*(*rdrho) + 1.0;
      int m = static_cast<int> (q);
      m = MIN(m,(*nrho)-1);
      q -= m;
      q = MIN(q,1.0);
      // for now, assume itype = 1
      double *coeff = (*frho_spline)[type2frho[1]][m];
      return ((coeff[0]*q + coeff[1])*q + coeff[2]);
    }
    //! Returns the second derivative of the atomic charge density.
    double F_pp(const double &p) const
    {
      double q = p*(*rdrho) + 1.0;
      int m = static_cast<int> (q);
      m = MIN(m,(*nrho)-1);
      q -= m;
      q = MIN(q,1.0);
      // for now, assume itype = 1
      double *coeff = (*frho_spline)[type2frho[1]][m];
      return ((*rdrho)*(2.0*coeff[0]*q + coeff[1]));
    }
    //! Returns the third derivative of the atomic charge density.
    double F_ppp(const double &p) const
    {
      double q = p*(*rdrho) + 1.0;
      int m = static_cast<int> (q);
      m = MIN(m,(*nrho)-1);
      q -= m;
      q = MIN(q,1.0);
      // for now, assume itype = 1
      double *coeff = (*frho_spline)[type2frho[1]][m];
      return ((*rdrho)*(*rdrho)*2.0*coeff[0]);
    }
    int *nrho,*nr,*nfrho,*nrhor,*nz2r;
    int *type2frho,**type2rhor,**type2z2r;
    double *cutmax;
    double *dr,*rdr,*drho,*rdrho;
    double ****rhor_spline,****frho_spline,****z2r_spline;
    LAMMPS_NS::PairEAM* lammps_eam;
  };
}
#endif

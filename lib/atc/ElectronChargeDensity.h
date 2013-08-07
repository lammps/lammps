#ifndef ELECTRON_DENSITY_H
#define ELECTRON_DENSITY_H

#include <map>
#include <string>

using std::map;
using std::string;

#include "ATC_TypeDefs.h"
#include "Function.h"

const double tol = 1.e-8;

namespace ATC {

  /**
   *  @class  ElectronChargeDensity 
   *  @brief  Base class for models of extrinsic electric charges 
   */

  class ElectronChargeDensity
  {
    public:
      ElectronChargeDensity()   {};
      virtual ~ElectronChargeDensity() {};
      virtual bool electron_charge_density(const FIELD_MATS &fields,
                                           DENS_MAT &flux) const { return false; };

      
      virtual void D_electron_charge_density(const FieldName fieldName, 
                                             const FIELD_MATS &fields,
                                             DENS_MAT &flux) const 
        { throw ATC_Error("Charge density D_electron_charge_density unimplemented function");}

      virtual void band_edge_potential(const FIELD_MATS &fields,
                                       DENS_MAT &density) const
        { throw ATC_Error("Charge density band_edge_potential unimplemented function");}
  };
  //-----------------------------------------------------------------------
  /**
   *  @class  ElectronChargeDensityInterpolation
   *  @brief  Class for models of electron charge density as a tabular function of electric potential
   */

  class ElectronChargeDensityInterpolation : public ElectronChargeDensity
  {
    public:
      ElectronChargeDensityInterpolation(fstream &matfile,map<string,double> & parameters);
      virtual ~ElectronChargeDensityInterpolation() {};
      virtual bool electron_charge_density(const FIELD_MATS &fields,
                                           DENS_MAT &flux) const
      {
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        const DENS_MAT & phi = phi_field->second;
        int nNodes  = phi.nRows();
        flux.reset(nNodes,1,false); 
        for (int i = 0; i < nNodes; i++) {  // a mapping of a vector
          flux(i,0) = n_.f(phi(i,0));  
        }
        flux *= -1.;
        return true;
      }; 
      virtual void D_electron_charge_density(const FieldName field,
                                             const FIELD_MATS &fields,
                                             DENS_MAT &coef) const
      {
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        const DENS_MAT & phi = phi_field->second;
        int nNodes  = phi.nRows();
        coef.reset(nNodes,1,false); 
        for (int i = 0; i < nNodes; i++) { 
          coef(i,0) = n_.dfdt(phi(i,0)); 
          coef(i,0) = n_.dfdt(phi(i,0)); 
        }
        coef *= -1.;
      }
    private:
      InterpolationFunction n_;
  };
  //-----------------------------------------------------------------------
  /**
   *  @class  ElectronChargeDensityLinear
   *  @brief  Class for models of electron charge density proportional to electric potential
   */

  class ElectronChargeDensityLinear : public ElectronChargeDensity
  {
    public:
      ElectronChargeDensityLinear(fstream &matfile,map<string,double> & parameters);
      virtual ~ElectronChargeDensityLinear() {};
      virtual bool electron_charge_density(const FIELD_MATS &fields,
                                           DENS_MAT &flux) const
      {
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        flux = phi_field->second;
        flux *= -C_;
        return true;
      }; 
      virtual void D_electron_charge_density(const FieldName field,
                                             const FIELD_MATS &fields,
                                             DENS_MAT &coef) const
      {
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        const DENS_MAT & phi = phi_field->second;
        int nNodes  = phi.nRows();
        coef.reset(nNodes,1,false);
        coef = -C_;
      }
    private:
      double C_;
  };
  //-----------------------------------------------------------------------
  /**
   *  @class  ElectronChargeDensityExponential
   *  @brief  Class for models of electron charge density dependent on difference between electric potential and the Fermi level n = n_i exp ( (phi-E_i) / kB T)
   */ 

  class ElectronChargeDensityExponential : public ElectronChargeDensity
  {
    public:
      ElectronChargeDensityExponential(fstream &matfile,map<string,double> & parameters);
      virtual ~ElectronChargeDensityExponential() {};

      double n(const double phi, double T)  const
      {
        return -intrinsicConcentration_*exp((phi-intrinsicEnergy_)/(kBeV_*T)); 
      }
      double dndphi(const double phi, double T)  const
      {
        return n(phi,T)/(kBeV_*T);
      }
      virtual bool electron_charge_density(const FIELD_MATS &fields, 
                     DENS_MAT &density) const
      {
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        FIELD_MATS::const_iterator T_field   = fields.find(TEMPERATURE);
        double T = 300;
        bool hasTref = (referenceTemperature_ > 0 );
        const DENS_MAT & phi = phi_field->second;
        int nNodes = phi.nRows();
        density.resize(nNodes,1);
        if (hasTref) {
          T = referenceTemperature_;
          for (int i = 0; i < nNodes; i++) { 
            density(i,0) = n(phi(i,0),T); }
        }
        else {
          const DENS_MAT & temp = T_field->second;
          for (int i = 0; i < nNodes; i++) { 
            density(i,0) = n(phi(i,0),temp(i,0)); }
        }
        density *= -1.;
        return true;
      }; 
      virtual void D_electron_charge_density(const FieldName field,
                                             const FIELD_MATS &fields,
                                             DENS_MAT &coef) const
      {
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        FIELD_MATS::const_iterator T_field   = fields.find(TEMPERATURE);
        double T = 300;
        bool hasTref = (referenceTemperature_ > 0 );
        const DENS_MAT & phi = phi_field->second;
        int nNodes = phi.nRows();
        coef.resize(nNodes,1);
        if (hasTref) {
          T = referenceTemperature_;
          for (int i = 0; i < nNodes; i++) { 
            coef(i,0) = dndphi(phi(i,0),T); }
        }
        else {
          const DENS_MAT & temp = T_field->second;
          for (int i = 0; i < nNodes; i++) { 
            coef(i,0) = dndphi(phi(i,0),temp(i,0)); }
        }
        coef *= -1.;
      }; 
    protected:
      double intrinsicConcentration_,intrinsicEnergy_;
      double referenceTemperature_;
  };

  //-----------------------------------------------------------------------
  /**
   *  @class  ElectronChargeDensityFermiDirac
   *  @brief  Class for models of electron charge density based on Fermi-Dirac statistics 
   */ 

  class ElectronChargeDensityFermiDirac : public ElectronChargeDensity
  {
    public:
      ElectronChargeDensityFermiDirac(fstream &matfile,map<string,double> & parameters);
      virtual ~ElectronChargeDensityFermiDirac() {};
      double fermi_dirac(const double E, const double T) const
      {
        double f = 1.0;
        if      (T > 0)   f = 1.0 / ( exp((E-Ef_)/kBeV_/T)+1.0 );
        else if (E > Ef_) f = 0;
        return f;
      }; 
      virtual bool electron_charge_density(const FIELD_MATS &fields,
                                           DENS_MAT &density) const
      {
        // psi : the inhomogeneous solution
        FIELD_MATS::const_iterator psi_field = fields.find(ELECTRON_WAVEFUNCTION);

        const DENS_MAT & psi = psi_field->second;
        FIELD_MATS::const_iterator psis_field = fields.find(ELECTRON_WAVEFUNCTIONS);
        // if (psis_field==fields.end())
        //throw ATC_Error("Wavefunctions not defined");
        const DENS_MAT & psis = psis_field->second;
        FIELD_MATS::const_iterator E_field = fields.find(ELECTRON_WAVEFUNCTION_ENERGIES);
        const DENS_MAT & Es = E_field->second;
        FIELD_MATS::const_iterator T_field   = fields.find(ELECTRON_TEMPERATURE);
        const DENS_MAT & Ts = T_field->second;
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        const DENS_MAT & phi = phi_field->second;

        int nNodes = psi.nRows();
        density.reset(nNodes,1);  
        double T  = referenceTemperature_;
        int count = 0;
        for (int i = 0; i < nNodes; i++) {
          if (!hasReferenceTemperature_) { T = Ts(i,0); }
          int j = 0;
          for (j = 0; j < psis.nCols(); j++) {
            double E = Es(j,0); // Eigenvalue  
            double f = fermi_dirac(E,T);
            if (f < tol) break;
            else  count++;
            density(i,0) -= psis(i,j)*psis(i,j)*f;  // < 0
          }
          if (donorIonization_) { 
            double E = -1.0* phi(i,0);// units [eV] E = - |e| phi
            if ( E + Eb_ > Ef_+Ed_) density(i,0) += Nd_;  // > 0
          }
        }
        return true;
      }; 
      virtual void D_electron_charge_density(const FieldName fieldName, 
                                             const FIELD_MATS &fields,
                                             DENS_MAT &coef) const
      {
        FIELD_MATS::const_iterator phi_field = fields.find(ELECTRIC_POTENTIAL);
        const DENS_MAT & phi = phi_field->second;
        int nNodes  = phi.nRows();
        coef.reset(nNodes,1,false); 
      }

      virtual void band_edge_potential(const FIELD_MATS &fields,
                                       DENS_MAT &density) const
      {
        FIELD_MATS::const_iterator p_field   = fields.find(ELECTRIC_POTENTIAL);
        const DENS_MAT & phi = p_field->second;
        int nNodes  = phi.nRows();
        density.reset(nNodes,1,false);  
        density = Eb_;
      };

    protected:
      double Ef_;
      double referenceTemperature_;
      double Ed_, Nd_;
      double Eb_;
      bool hasReferenceTemperature_, donorIonization_;
  };
}

#endif 



#ifndef PHYSICS_MODEL_TWO_TEMPERATURE_H
#define PHYSICS_MODEL_TWO_TEMPERATURE_H

// included headers
#include "PhysicsModel.h"

namespace ATC{

class PhysicsModelTwoTemperature : public PhysicsModel {

  /** system:
      \dot T_e = diffusion_e + exchange_e
      \dot T_p = diffusion_p + exchange_p
  */

  
 public:
  
  // constructor (take material parameter/s)
  PhysicsModelTwoTemperature(string matFileName,
                             ATC_Transfer * atcTransfer);

  // destructor
  virtual ~PhysicsModelTwoTemperature();

  /** checks materials for necessary interfaces */
  virtual void initialize(void);

  virtual void get_num_fields(map<FieldName,int> & fieldSizes, 
                              Array2D<bool> & fieldMask) const
  { 
    fieldSizes[TEMPERATURE] = 1; fieldSizes[ELECTRON_TEMPERATURE] = 1;
    fieldMask(ELECTRON_TEMPERATURE,FLUX)   = true;
    fieldMask(ELECTRON_TEMPERATURE,SOURCE) = true;
  }

  /** energy */
  virtual void E_integrand(const Array<FieldName> & mask, 
                           const FIELDS &fields,
                           const GRAD_FIELDS &grad_fields,
                           FIELDS &energy,
                           const int matIndex = 0) const;

  /** capacity that used to form the mass matrix */
  virtual void M_integrand(const Array<FieldName> & mask, 
                           const FIELDS &fields,
                           FIELDS &flux,
                           const int matIndex = 0) const;

  /** this model has a B weighted integrand */
  virtual bool has_B_integrand() const {return true;}

  /** flux that is integrated with Grad N as its weight */
  virtual void B_integrand(const Array2D<bool> & mask, 
                           const FIELDS &fields,
                           const GRAD_FIELDS &grad_fields,
                           GRAD_FIELDS &flux,
                           const int matIndex = 0) const;

  /** this model has a N weighted integrand */
  virtual bool has_N_integrand() const {return true;}

  /** flux that is integrated with N as its weight */
  virtual void N_integrand(const Array2D<bool> &mask, 
                           const FIELDS &fields, 
                           const GRAD_FIELDS &grad_fields,
                           FIELDS &flux,
                           const int matIndex = 0) const;
};
}; // end namespace
#endif

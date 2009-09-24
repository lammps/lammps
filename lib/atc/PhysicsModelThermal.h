#ifndef PHYSICS_MODEL_THERMAL_H
#define PHYSICS_MODEL_THERMAL_H

#include "PhysicsModel.h"

namespace ATC{

  class PhysicsModelThermal : public PhysicsModel {

  public:

    // constructor (take material parameter/s)
    PhysicsModelThermal(std::string matFileName,
                        ATC_Transfer * atcTransfer);

    // destructor
    virtual ~PhysicsModelThermal();

    /** checks materials for necessary interfaces */
    virtual void initialize(void);

    virtual void get_num_fields(map<FieldName,int> & fieldSizes, Array2D<bool> & fieldMask) const
    { 
      fieldSizes[TEMPERATURE] = 1;
      fieldMask(TEMPERATURE,FLUX) = true;
    }

    /** heat capacity */
    virtual void M_integrand(const Array<FieldName> & mask, 
                             const FIELDS &fields,
                             FIELDS &capacity,
                             const int matIndex = 0) const;
    /** energy */
    virtual void E_integrand(const Array<FieldName> & mask, 
                             const FIELDS &fields,
                             const GRAD_FIELDS &grad_fields,
                             FIELDS &energy,
                             const int matIndex = 0) const;

    /** this model has a B weighted integrand */
    virtual bool has_B_integrand() const {return true;};

    /** flux that is integrated with Grad N as its weight */
    virtual void B_integrand(const Array2D<bool> & mask,
                             const FIELDS & fields,
                             const GRAD_FIELDS & grad_fields,
                             GRAD_FIELDS & flux,
                             const int matIndex = 0) const;

  };

};
#endif

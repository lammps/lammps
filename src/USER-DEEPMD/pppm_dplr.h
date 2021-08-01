#ifdef KSPACE_CLASS

KSpaceStyle(pppm/dplr,PPPMDPLR)

#else

#ifndef LMP_PPPM_DPLR_H
#define LMP_PPPM_DPLR_H

#ifdef HIGH_PREC
#define FLOAT_PREC double
#else
#define FLOAT_PREC float
#endif

#include "pppm.h"
#include <iostream>
#include <vector>

namespace LAMMPS_NS {

  class PPPMDPLR : public PPPM {
public:
#ifdef OLD_LMP_PPPM
    PPPMDPLR(class LAMMPS *, int, char **);
#else
    PPPMDPLR(class LAMMPS *);
#endif
    virtual ~PPPMDPLR () {};
    void init();
    const std::vector<double > & get_fele() const {return fele;};
protected:
    virtual void compute(int, int);
    virtual void fieldforce_ik();
    virtual void fieldforce_ad();    
private:
    std::vector<double > fele;
  };

}

#endif
#endif


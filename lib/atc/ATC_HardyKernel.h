/** ATC_HardyKernel: Hardy smoothing  */
#ifndef ATC_HARDY_KERNEL_H
#define ATC_HARDY_KERNEL_H

#include "LammpsInterface.h"
#include "MatrixLibrary.h"


namespace ATC {


  class ATC_HardyKernel {

  public:
  
    // constructor
    ATC_HardyKernel(int nparameters, double* parameters);
    // destructor
    virtual ~ATC_HardyKernel() {};
    // function value
    virtual double value(DENS_VEC& x_atom)=0; 
    // bond function value via quadrature
    virtual double bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2);
    // localization-volume intercepts for bond calculation
    virtual void bond_intercepts(DENS_VEC& xa, 
                 DENS_VEC& xb, double &lam1, double &lam2);
    double inv_vol(void) { return invVol_; }

  protected:
    double invVol_; // normalization factor
    double Rc_, invRc_; // cutoff radius
    int nsd_ ; // number of dimensions 

    /** pointer to lammps interface class */
    LammpsInterface * lammpsInterface_;

    /** periodicity flags and lengths */
    int periodicity[3];
    double box_bounds[2][3];
    double box_length[3];

  };

  /** a step with spherical support */
  class ATC_HardyKernelStep : public ATC_HardyKernel {

  public:
    // constructor
    ATC_HardyKernelStep(int nparameters, double* parameters);
    // destructor
    virtual ~ATC_HardyKernelStep() {};
    // function value
    double value(DENS_VEC& x_atom);
    // bond function value
    virtual double bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2)
      {return lam2 -lam1;}
  };

  /** a step with rectangular support suitable for a rectangular grid */
  class ATC_HardyKernelCell : public ATC_HardyKernel {

  public:
    // constructor
    ATC_HardyKernelCell(int nparameters, double* parameters); 
    // destructor
    virtual ~ATC_HardyKernelCell() {};
    // function value
    virtual double value(DENS_VEC& x_atom);
    // bond function value
    virtual double bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2)
      {return lam2 -lam1;}
    // bond intercept values : origin is the node position
    void bond_intercepts(DENS_VEC& xa,
                 DENS_VEC& xb, double &lam1, double &lam2);

  protected:
    double hx, hy, hz;
    DENS_VEC cellBounds_;
  };

  class ATC_HardyKernelCubicSphere : public ATC_HardyKernel {

  public:
    // constructor
    ATC_HardyKernelCubicSphere(int nparameters, double* parameters);
    // destructor
    virtual ~ATC_HardyKernelCubicSphere() {};
    // function value
    virtual double value(DENS_VEC& x_atom);
  };

  class ATC_HardyKernelQuarticSphere : public ATC_HardyKernel {

  public:
    // constructor
    ATC_HardyKernelQuarticSphere(int nparameters, double* parameters);
    // destructor
    virtual ~ATC_HardyKernelQuarticSphere() {};
    // function value
    virtual double value(DENS_VEC& x_atom);
  };


  class ATC_HardyKernelCubicCyl : public ATC_HardyKernel {

  public:
    // constructor
    ATC_HardyKernelCubicCyl(int nparameters, double* parameters);
    // destructor
    virtual ~ATC_HardyKernelCubicCyl() {};
    // function value
    virtual double value(DENS_VEC& x_atom) ;
  };

  class ATC_HardyKernelQuarticCyl : public ATC_HardyKernel {

  public:
  
    // constructor
    ATC_HardyKernelQuarticCyl(int nparameters, double* parameters);
    // destructor
    virtual ~ATC_HardyKernelQuarticCyl() {};
    // function value
    virtual double value(DENS_VEC& x_atom);
  };
};
#endif

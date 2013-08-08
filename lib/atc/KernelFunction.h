/** KernelFunction: Hardy smoothing  */
#ifndef KERNEL_FUNCTION_H
#define KERNEL_FUNCTION_H

#include <set>
#include "LammpsInterface.h"
#include "MatrixLibrary.h"

namespace ATC {

  /**
   *  @class  KernelFunctionMgr
   *  @brief  Base class for managing kernels 
   */  
  class KernelFunctionMgr {
  public:
   /** Static instance of this class */
    static KernelFunctionMgr * instance();
    class KernelFunction* function(char** arg, int nargs);
  protected:
   KernelFunctionMgr() {};
   ~KernelFunctionMgr();
  private:
    static KernelFunctionMgr * myInstance_;
    set<KernelFunction*> pointerSet_;
  };

  /**
   *  @class  KernelFunction 
   *  @brief  Base class for kernels for atom-continuum transfer    
   */  

  class KernelFunction {

  public:
  
    // constructor
    KernelFunction(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunction() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const =0 ;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const =0 ;
    // bond function value via quadrature
    virtual double bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2) const;
    // localization-volume intercepts for bond calculation
    virtual void bond_intercepts(DENS_VEC& xa, 
                 DENS_VEC& xb, double &lam1, double &lam2) const;
    virtual bool node_contributes(DENS_VEC node) const;
    virtual bool in_support(DENS_VEC node) const;
    double inv_vol(void) const { return invVol_; }
    virtual double dimensionality_factor(void) const { return 1.; }

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

  /**
   *  @class  KernelFunctionStep 
   *  @brief  Class for defining kernel function of a step with spherical support
   */  
  
  class KernelFunctionStep : public KernelFunction {

  public:
    // constructor
    KernelFunctionStep(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionStep() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;
    // bond function value
    virtual double bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2) const
      { return lam2-lam1; }
  };

  /**
   *  @class  KernelFunctionCell 
   *  @brief  Class for defining kernel function of a step with rectangular support
   *          suitable for a rectangular grid
   */  
  
  class KernelFunctionCell : public KernelFunction {

  public:
    // constructor
    KernelFunctionCell(int nparameters, double* parameters); 
    // destructor
    virtual ~KernelFunctionCell() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
    // bond function value
    virtual double bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2) const
      {return lam2 -lam1;}
    // bond intercept values : origin is the node position
    void bond_intercepts(DENS_VEC& xa, DENS_VEC& xb, 
                         double &lam1, double &lam2) const;
    bool node_contributes(DENS_VEC node) const;
    bool in_support(DENS_VEC dx) const;

  protected:
    double hx, hy, hz;
    DENS_VEC cellBounds_;
  };

  /**
   *  @class  KernelFunctionCubicSphere
   *  @brief  Class for defining kernel function of a cubic with spherical support
   */

  class KernelFunctionCubicSphere : public KernelFunction {

  public:
    // constructor
    KernelFunctionCubicSphere(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionCubicSphere() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
  };

  /**
   *  @class  KernelFunctionQuarticSphere
   *  @brief  Class for defining kernel function of a quartic with spherical support
   */

  class KernelFunctionQuarticSphere : public KernelFunction {

  public:
    // constructor
    KernelFunctionQuarticSphere(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionQuarticSphere() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
  };

  /**
   *  @class  KernelFunctionCubicCyl
   *  @brief  Class for defining kernel function of a cubic with cylindrical support
   */

  class KernelFunctionCubicCyl : public KernelFunction {

  public:
    // constructor
    KernelFunctionCubicCyl(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionCubicCyl() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
    virtual double dimensionality_factor(void) const { return 0.5; }
  };

  /**
   *  @class  KernelFunctionQuarticCyl
   *  @brief  Class for defining kernel function of a quartic with cylindrical support
   */

  class KernelFunctionQuarticCyl : public KernelFunction {

  public:
  
    // constructor
    KernelFunctionQuarticCyl(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionQuarticCyl() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
    virtual double dimensionality_factor(void) const { return 0.5; }
  };

  /**
   *  @class  KernelFunctionCubicBar
   *  @brief  Class for defining kernel function of a cubic with 1-dimensional (bar) support 
   */

  class KernelFunctionCubicBar : public KernelFunction {

  public:
    // constructor
    KernelFunctionCubicBar(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionCubicBar() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
    virtual double dimensionality_factor(void) const { return 0.25; }
  };

  /**
   *  @class  KernelFunctionCubicBar
   *  @brief  Class for defining kernel function of a cubic with 1-dimensional (bar) support 
   */

  class KernelFunctionLinearBar : public KernelFunction {

  public:
    // constructor
    KernelFunctionLinearBar(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionLinearBar() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
    virtual double dimensionality_factor(void) const { return 0.25; }
  };

  /**
   *  @class  KernelFunctionQuarticBar
   *  @brief  Class for defining kernel function of a quartic with 1-dimensional (bar) support
   */

  class KernelFunctionQuarticBar : public KernelFunction {

  public:

    // constructor
    KernelFunctionQuarticBar(int nparameters, double* parameters);
    // destructor
    virtual ~KernelFunctionQuarticBar() {};
    // function value
    virtual double value(DENS_VEC& x_atom) const;
    // function derivative 
    virtual void derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const;     
    virtual double dimensionality_factor(void) const { return 0.25; }
  };

};
#endif

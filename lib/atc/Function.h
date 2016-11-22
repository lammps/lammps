#ifndef XT_FUNCTION_H
#define XT_FUNCTION_H

#include <math.h>
#include <string>
#include <set>
#include <cstdlib>
#include <iostream>

#include "Array.h"
#include "ATC_TypeDefs.h"

namespace ATC {

  /**
   *  @class Function
   *  @brief Base class for functions of fields, space and time
   */

  class Function {
  public:
    Function(int nargs, char** args);
    virtual ~Function(void) {};

    /** name */
    const std::string & tag() { return tag_;}

    /** depdendencies */
    virtual inline ARG_NAMES args(void) {ARG_NAMES names; return names;};

    /** (1st) derivative of function wrt to a field */
    virtual inline double dfd(FieldName field, ARGS& args ) {return 0.0;};
    virtual inline void   dfd(FieldName field, ARGS& args, DENS_MAT vals ) {};

    // addl: d2fd2(field1, field2, args), linearization(), grad_args

  protected:
    /** tag : name of function */
    std::string tag_;
  };

  /**
   *  @class Function_Mgr
   *  @brief Base class that constructs and returns UXT_Function objects based on requests
   */
  class Function_Mgr {
  public:
    /** Static instance of this class */
    static Function_Mgr * instance();

    Function* function(char ** arg, int nargs);
    Function* copy_function(Function* other);
   protected:
    Function_Mgr() {};
    ~Function_Mgr();
   private:
    static Function_Mgr * myInstance_;
    /** set to store all generated objects for later deletion */
    std::set<Function * > pointerSet_;
  };


  /** 
   *  @class LinearFieldFunction
   *  @brief Class for functions returning values linear a given field
   */

  class LinearFieldFunction : public Function {
  public:
    LinearFieldFunction(int nargs, char** args);
    virtual ~LinearFieldFunction(void) {};
  
    inline double f(double* u, double* x, double t) {return c1_*u[0]-c0_;}
    inline double dfd(FieldName field, ARGS& args) {return c1_;}

    private :
      double c0_,c1_;
  };

  /**
   *  @class UXT_Function
   *  @brief Base class for functions of fields, space and time
   */

  class UXT_Function {
  public:
    UXT_Function(int nargs, double* args);
    virtual ~UXT_Function(void) {};

    const std::string & tag() { return tag_;}

    /** function value */
    virtual inline double f(double * u, double* x, double t) {return 0.0;};
    /** derivative of function wrt to field */
    virtual inline double dfdu(double * u, double* x, double t) {return 0.0;};

  protected:
    /** tag : name of function */
    std::string tag_;
  };

  /**
   *  @class UXT_Function_Mgr
   *  @brief Base class that constructs and returns UXT_Function objects based on requests
   */


  class UXT_Function_Mgr {
  public:
    /** Static instance of this class */
    static UXT_Function_Mgr * instance();

    UXT_Function* function(std::string & type, int nargs, double * arg);
    UXT_Function* function(char ** arg, int nargs);
    UXT_Function* linear_function(double c0, double c1);
    UXT_Function* copy_UXT_function(UXT_Function* other);
   protected:
    UXT_Function_Mgr() {};
    ~UXT_Function_Mgr();
   private:
    static UXT_Function_Mgr * myInstance_;
    /** set to store all generated objects for later deletion */
    std::set<UXT_Function * > pointerSet_;
  };


  /** 
   *  @class ScalarLinearFunction
   *  @brief Class for functions returning values linear in space 
   */

  class ScalarLinearFunction : public UXT_Function {
  public:
    ScalarLinearFunction(int nargs, double* args);
    virtual ~ScalarLinearFunction(void) {};
  
    //inline double f(double* u, double* x, double t) {return c1_*(u[0]-c0_);}

    inline double f(double* u, double* x, double t) {return c1_*u[0]+c0_;}
    inline double dfdu(double* u, double* x, double t) {return c1_;}

    private :
      double c0_,c1_;
  };

  /**
   *  @class XT_Function
   *  @brief Base class for functions based on space and time variables 
   */

  class XT_Function {
  public:
    XT_Function(int nargs, double* args);
    virtual ~XT_Function(void) {};

    const std::string & tag() { return tag_;}

    /** function value */
    virtual inline double f(double* x, double t) {return 0.0;};
    /** time derivative of function */
    virtual inline double dfdt(double* x, double t) {return 0.0;};
    /** 2nd time derivative of function */
    virtual inline double ddfdt(double* x, double t) {return 0.0;};
    /** 3rd time derivative of function */
    virtual inline double dddfdt(double* x, double t) {return 0.0;};

  protected:
    /** mask : masks x,y,z dependence, x0 : origin */
    double mask[3], x0[3];
    /** tag : name of function */
    std::string tag_;
  };

  /**
   *  @class XT_Function_Mgr
   *  @brief Base class that constructs and returns XT_Function objects based on requests
   */

  class XT_Function_Mgr {
  public:
    /** Static instance of this class */
    static XT_Function_Mgr * instance();

    XT_Function* function(std::string & type, int nargs, double * arg);
    XT_Function* function(char ** arg, int nargs);
    XT_Function* constant_function(double c);
    XT_Function* copy_XT_function(XT_Function* other);
   protected:
    XT_Function_Mgr() {};
    ~XT_Function_Mgr();
   private:
    static XT_Function_Mgr * myInstance_;
    /** set to store all generated objects for later deletion */
    std::set<XT_Function * > pointerSet_;

  };

  //------------------------------------------------------------------------
  // derived classes
  //------------------------------------------------------------------------

  /**
   *  @class ConstantFunction
   *  @brief Class for functions returning constant values 
   */

  class ConstantFunction : public XT_Function {
  public:
    ConstantFunction(int nargs, double* args);
    ConstantFunction(double arg);
    virtual ~ConstantFunction(void) {};
  
    inline double f(double* x, double t) 
    {return C0;};

    private :
      double C0;
  };

  /** 
   *  @class LinearFunction
   *  @brief Class for functions returning values linear in space 
   */

  class LinearFunction : public XT_Function {
  public:
    LinearFunction(int nargs, double* args);
    virtual ~LinearFunction(void) {};
  
    double f(double* x, double t) 
      {return mask[0]*(x[0]-x0[0])+mask[1]*(x[1]-x0[1])+mask[2]*(x[2]-x0[2]) + C0;};

    private :
      double C0;
  };
  /** 
   *  @class PiecewiseLinearFunction
   *  @brief Class for functions returning values piecewise linear in space 
   *  along given direction
   */

  class PiecewiseLinearFunction : public XT_Function {
  public:
    PiecewiseLinearFunction(int nargs, double* args);
    virtual ~PiecewiseLinearFunction(void) {};
  
    double f(double* x, double t) ;

    private :
      Array<double> xi; 
      Array<double> fi; 
  };

  /** 
   *  @class LinearTemporalRamp
   *  @brief Class for functions returning values linear in space and time 
   */

  class LinearTemporalRamp : public XT_Function {
  public:
    LinearTemporalRamp(int nargs, double* args);
    ~LinearTemporalRamp(void) {};
  
    double f(double* x, double t);
    double dfdt(double* x, double t);
    
    protected :
      double mask_slope[3];
      double C0_initial, C0_slope;
    
  };

  /**
   *  @class QuadraticFunction
   *  @brief Class for functions returning values quadratic in space
   */

  class QuadraticFunction : public XT_Function {
  public:
    QuadraticFunction(int nargs, double* args);
    virtual ~QuadraticFunction(void) {};
  
    inline double f(double* x, double t) 
      {return 
       C2[0]*(x[0]-x0[0])*(x[0]-x0[0])+
       C2[1]*(x[1]-x0[1])*(x[1]-x0[1])+
       C2[2]*(x[2]-x0[2])*(x[2]-x0[2])+
       2.0*C2[3]*(x[0]-x0[0])*(x[1]-x0[1]) +
       2.0*C2[4]*(x[0]-x0[0])*(x[2]-x0[2]) +
       2.0*C2[5]*(x[1]-x0[1])*(x[2]-x0[2]) +
       mask[0]*(x[0]-x0[0])+mask[1]*(x[1]-x0[1])+mask[2]*(x[2]-x0[2]) + C0;};

    private :
      double C0, C2[6]; // C2 1:xx 2:yy 3:zz 4:xy|yx 5:xz|zx 6:yz|zy
  };

  /**
   *  @class SineFunction 
   *  @brief Class for functions returning values sinusoidally varying in space and time
   */

  class SineFunction : public XT_Function {
  public:
    SineFunction(int nargs, double* args);
    virtual ~SineFunction(void){};

    inline double f(double* x, double t) 
      {return  C*sin( mask[0]*(x[0]-x0[0])
                     +mask[1]*(x[1]-x0[1])
                     +mask[2]*(x[2]-x0[2]) - w*t) + C0;};

    private :
      double C, C0, w;
  };

  /**
   *  @class GaussianFunction
   *  @brief Class for functions returning values according to a Gaussian distribution in space
   */

  class GaussianFunction : public XT_Function {
  public:
    GaussianFunction(int nargs, double* args);
    virtual ~GaussianFunction(void){};

    // 1/(2 pi \sigma)^(n/2) exp(-1/2 x.x/\sigma^2 ) for n = dimension
    inline double f(double* x, double t) 
      {return  C*exp(-(mask[0]*(x[0]-x0[0])*(x[0]-x0[0])
                      +mask[1]*(x[1]-x0[1])*(x[1]-x0[1])
                      +mask[2]*(x[2]-x0[2])*(x[2]-x0[2]))/tau/tau) + C0;};

    protected:
      double tau, C, C0;
  };

  /**
   *  @class GaussianTemporalRamp 
   *  @brief Class for functions returning values according to a Gaussian distribution in space and linearly in time
   */

  class GaussianTemporalRamp : public GaussianFunction {
  public:
    GaussianTemporalRamp(int nargs, double* args);
    virtual ~GaussianTemporalRamp(void){};

    double f(double* x, double t);
    double dfdt(double* x, double t);

  protected:
    double tau_initial, tau_slope;
    double C_initial, C_slope;
    double C0_initial, C0_slope;
  };
  
  /**
   *  @class TemporalRamp 
   *  @brief Class for functions returning values constant in space and varying linearly in time
   */

  class TemporalRamp : public XT_Function {
  public:
    TemporalRamp(int nargs, double* args);
    virtual ~TemporalRamp(void) {};
  
    inline double f(double* x, double t) 
    {return f_initial + slope*t;};

    inline double dfdt(double* x, double t) 
    {return slope;};

    private :
      double f_initial, slope;
  };

  /**
   *  @class RadialPower 
   *  @brief Class for functions returning values based on distance from a fix point raised to a specified power 
   */

  class RadialPower : public XT_Function {
  public:
    RadialPower(int nargs, double* args);
    virtual ~RadialPower(void) {};
    
    inline double f(double* x, double t) 
    {
      double dx = x[0]-x0[0]; double dy = x[1]-x0[1]; double dz = x[2]-x0[2];
      double r = mask[0]*dx*dx+mask[1]*dy*dy+mask[2]*dz*dz; r = sqrt(r);
      return C0*pow(r,n);
    };
    
  private :
    double C0, n;
  };

  /**
   *  @class InterpolationFunction
   *  @brief Base class for interpolation functions 
   */

  class InterpolationFunction {
  public:
    InterpolationFunction(void) : npts_(0) {};
    virtual ~InterpolationFunction(void) {};

    /** populate data */
    void initialize(int npts,std::fstream &fileId, double coef = 1.);

    /** function value */
    double f(const double t) const;
    /** derivative of function */
    double dfdt(const double t) const;

  protected:
    double coordinate(double x, 
      double & f0, double & fp0, double & f1, double & fp1, double & inv_dx) const;
    int npts_;
    Array<double> xs_;
    Array<double> fs_;
    Array<double> fps_;
  };

}
#endif

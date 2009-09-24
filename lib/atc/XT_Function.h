#ifndef XT_FUNCTION_H
#define XT_FUNCTION_H

#include <math.h>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

namespace ATC {
  //------------------------------------------------------------------------
  // base class
  //------------------------------------------------------------------------
  class XT_Function {
  public:
    XT_Function(int nargs, double* args);
    ~XT_Function(void) {};

    const string & get_tag() { return tag;}

    /** function value */
    virtual inline double f(double* x, double t) {return 0.0;};
    /** time derivative of function */
    virtual inline double dfdt(double* x, double t) {return 0.0;};
    /** 2nd time derivative of funtion */
    virtual inline double ddfdt(double* x, double t) {return 0.0;};
    /** 3rd time derivative of funtion */
    virtual inline double dddfdt(double* x, double t) {return 0.0;};

  protected:
    /** mask : masks x,y,z dependence, x0 : origin */
    double mask[3], x0[3];
    /** tag : name of function */
    string tag;
  };

  //------------------------------------------------------------------------
  // manager
  //------------------------------------------------------------------------
  class XT_Function_Mgr {
  public:
    /** Static instance of this class */
    static XT_Function_Mgr * instance();

    XT_Function* get_function(string & type, int nargs, double * arg);
    XT_Function* get_function(char ** arg, int nargs);
    XT_Function* get_constant_function(double c);
    XT_Function* copy_XT_function(XT_Function* other);
   protected:
    XT_Function_Mgr() {};
    ~XT_Function_Mgr();
    /** Vector of ptrs to functions created thus far */
    std::vector<XT_Function *> xtFunctionVec_;
   private:
    static XT_Function_Mgr * myInstance_;
  };

  //------------------------------------------------------------------------
  // derived classes
  //------------------------------------------------------------------------

  /** a constant  */
  class ConstantFunction : public XT_Function {
  public:
    ConstantFunction(int nargs, double* args);
    ~ConstantFunction(void) {};
  
    inline double f(double* x, double t) 
      {return C0;};

    private :
      double C0;
  };

  /** a linear in x-y-z */
  class LinearFunction : public XT_Function {
  public:
    LinearFunction(int nargs, double* args);
    ~LinearFunction(void) {};
  
    inline double f(double* x, double t) 
      {return mask[0]*(x[0]-x0[0])+mask[1]*(x[1]-x0[1])+mask[2]*(x[2]-x0[2]) + C0;};

    private :
      double C0;
  };

  /** a quadratic in x-y-z */
  class QuadraticFunction : public XT_Function {
  public:
    QuadraticFunction(int nargs, double* args);
    ~QuadraticFunction(void) {};
  
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

  class SineFunction : public XT_Function {
  public:
    SineFunction(int nargs, double* args);
    ~SineFunction(void);

    inline double f(double* x, double t) 
      {return  C*sin( mask[0]*(x[0]-x0[0])
                     +mask[1]*(x[1]-x0[1])
                     +mask[2]*(x[2]-x0[2]) - w*t);};

    private :
      double C, w;
  };

  /** a spatial gaussian function */
  class GaussianFunction : public XT_Function {
  public:
    GaussianFunction(int nargs, double* args);
    ~GaussianFunction(void){};

    // 1/(2 pi \sigma)^(n/2) exp(-1/2 x.x/\sigma^2 ) for n = dimension
    inline double f(double* x, double t) 
      {return  C*exp(-(mask[0]*(x[0]-x0[0])*(x[0]-x0[0])
                      +mask[1]*(x[1]-x0[1])*(x[1]-x0[1])
                      +mask[2]*(x[2]-x0[2])*(x[2]-x0[2]))/tau/tau) + C0;};

    protected:
      double tau, C, C0;
  };

  /** a spatial gaussian function that has variables ramp up in time */
  class GaussianTemporalRamp : public GaussianFunction {
  public:
    GaussianTemporalRamp(int nargs, double* args);
    ~GaussianTemporalRamp(void){};

    double f(double* x, double t);
    double dfdt(double* x, double t);

  protected:
    double tau_initial, tau_slope;
    double C_initial, C_slope;
    double C0_initial, C0_slope;
  };
  
  /** a ramp in time */
  class TemporalRamp : public XT_Function {
  public:
    TemporalRamp(int nargs, double* args);
    ~TemporalRamp(void) {};
  
    inline double f(double* x, double t) 
    {return f_initial + slope*t;};

    inline double dfdt(double* x, double t) 
      {return slope;};

    private :
      double f_initial, slope;
  };


}
#endif

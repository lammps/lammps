#include "Function.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"
#include <sstream>

namespace ATC {

  //====================================================================
  //  UXT_Function
  //===================================================================
  UXT_Function::UXT_Function(int narg, double* args) { }
  //====================================================================
  //  UXT_Function_Mgr
  //====================================================================
  UXT_Function_Mgr * UXT_Function_Mgr::myInstance_ = NULL;
  // -----------------------------------------------------------------
  //  instance()
  // -----------------------------------------------------------------
  UXT_Function_Mgr * UXT_Function_Mgr::instance()
  {
    if (myInstance_ == NULL) {
      myInstance_ = new UXT_Function_Mgr();
    }
    return myInstance_;
  }

  // Destructor
  UXT_Function_Mgr::~UXT_Function_Mgr()
  {
    // Delete all functions created using "new"
    set<UXT_Function * >::iterator it;
    for (it = pointerSet_.begin(); it != pointerSet_.end(); it++)
      if (*it) delete *it;
  }

  // add user function into the if statement and assign returnFunction to it
  UXT_Function* UXT_Function_Mgr::function(string & type, int nargs, double * args)
  {
    UXT_Function * returnFunction;
    if      (type=="linear") {
      returnFunction = new ScalarLinearFunction(nargs,args);
    }
    else
      throw ATC_Error("Bad user function name");

    pointerSet_.insert(returnFunction);

    return returnFunction;
  }

  // add user function into the if statement and assign returnFunction to it
  UXT_Function* UXT_Function_Mgr::function(char ** args, int nargs)
  {
    string type = args[0];
    int narg = nargs -1;
    double dargs[narg];
    for (int i = 0; i < narg; ++i) dargs[i] = atof(args[i+1]);
  
    return function(type, narg, dargs);
  }

  // add constant function
  UXT_Function* UXT_Function_Mgr::linear_function(double c0, double c1)
  {
    double args[2] = {c0,c1};
    UXT_Function * returnFunction = new ScalarLinearFunction(2,args);
    pointerSet_.insert(returnFunction);
    return (returnFunction);
  }

  UXT_Function* UXT_Function_Mgr::copy_UXT_function(UXT_Function* other)
  {
    string tag = other->tag();
  
    UXT_Function * returnFunction = NULL;
    if (tag=="linear") {
      ScalarLinearFunction * other_cast = (ScalarLinearFunction*) other;
      returnFunction = new ScalarLinearFunction(*other_cast);
    }
    pointerSet_.insert(returnFunction);
    return returnFunction;
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  ScalarLinearFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  ScalarLinearFunction::ScalarLinearFunction(int narg, double* args) 
    : UXT_Function(narg,args)
  {
    tag_ = "linear";
    c0_ = args[0];
    c1_ = args[1];
    stringstream ss;
    ss << "created function : " << c0_ << " + " << c1_ << "*u";
    ATC::LammpsInterface::instance()->print_msg_once(ss.str());
  }

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  XT_Function
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  XT_Function::XT_Function(int narg, double* args) 
  {

    if (narg > 5 ) { 
      x0[0] = args[0]; 
      x0[1] = args[1];
      x0[2] = args[2];
      mask[0] = args[3];
      mask[1] = args[4];
      mask[2] = args[5];
    } 
    else {
      x0[0] = 0.0; 
      x0[1] = 0.0;
      x0[2] = 0.0;
      mask[0] = 0.0;
      mask[1] = 0.0;
      mask[2] = 0.0;
    }
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  XT_Function_Mgr
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
XT_Function_Mgr * XT_Function_Mgr::myInstance_ = NULL;

// -----------------------------------------------------------------
//  instance()
// -----------------------------------------------------------------
  XT_Function_Mgr * XT_Function_Mgr::instance()
  {
    if (myInstance_ == NULL) {
      myInstance_ = new XT_Function_Mgr();
    }
    return myInstance_;
  }

  // Destructor
  XT_Function_Mgr::~XT_Function_Mgr()
  {
    // Delete all functions created using "new"
    set<XT_Function * >::iterator it;
    for (it = pointerSet_.begin(); it != pointerSet_.end(); it++)
      if (*it) delete *it;
  }

  // add user function into the if statement and assign returnFunction to it
  XT_Function* XT_Function_Mgr::function(string & type, int nargs, double * args)
  {
    XT_Function * returnFunction;
    if      (type=="constant") {
      returnFunction = new ConstantFunction(nargs,args);
    }
    else if (type=="temporal_ramp") {
      returnFunction = new TemporalRamp(nargs,args);
    }
    else if (type=="linear")
      returnFunction = new LinearFunction(nargs,args);
    else if (type=="piecewise_linear")
      returnFunction = new PiecewiseLinearFunction(nargs,args);
    else if (type=="linear_temporal_ramp")
      returnFunction = new LinearTemporalRamp(nargs,args);
    else if (type=="quadratic")
      returnFunction = new QuadraticFunction(nargs,args);
    else if (type=="sine")
      returnFunction = new SineFunction(nargs,args);
    else if (type=="gaussian")
      returnFunction = new GaussianFunction(nargs,args);
    else if (type=="gaussian_temporal_ramp")
      returnFunction = new GaussianTemporalRamp(nargs,args);
    else if (type=="radial_power")
      returnFunction = new RadialPower(nargs,args);
    else
      throw ATC_Error("Bad user function name");

    pointerSet_.insert(returnFunction);

    return returnFunction;
  }

  // add user function into the if statement and assign returnFunction to it
  XT_Function* XT_Function_Mgr::function(char ** args, int nargs)
  {
    string type = args[0];
    int narg = nargs -1;
    double dargs[narg];
    for (int i = 0; i < narg; ++i) dargs[i] = atof(args[i+1]);
  
    return function(type, narg, dargs);
  }

  // add constant function
  XT_Function* XT_Function_Mgr::constant_function(double c)
  {
    XT_Function * returnFunction = new ConstantFunction(c);
    pointerSet_.insert(returnFunction);
    return (returnFunction);
  }

  XT_Function* XT_Function_Mgr::copy_XT_function(XT_Function* other)
  {
    string tag = other->tag();
  
    XT_Function * returnFunction = NULL;
    if (tag=="linear") {
      LinearFunction * other_cast = (LinearFunction*) other;
      returnFunction = new LinearFunction(*other_cast);
    }
    else if (tag=="piecewise_linear") {
      PiecewiseLinearFunction * other_cast = (PiecewiseLinearFunction*) other;
      returnFunction = new PiecewiseLinearFunction(*other_cast);
    }
    else if (tag=="quadratic") {
      QuadraticFunction * other_cast = (QuadraticFunction*) other;
      returnFunction = new QuadraticFunction(*other_cast);
    }
    else if (tag=="sine") {
      SineFunction * other_cast = (SineFunction*) other;
      returnFunction = new SineFunction(*other_cast);
    }
    else if (tag=="gaussian") {
      GaussianFunction * other_cast = (GaussianFunction*) other;
      returnFunction = new GaussianFunction(*other_cast);
    }
    else if (tag=="gaussian_temporal_ramp") {
      GaussianTemporalRamp * other_cast = (GaussianTemporalRamp*) other;
      returnFunction = new GaussianTemporalRamp(*other_cast);
    }
    else if (tag=="temporal_ramp") {
      TemporalRamp * other_cast = (TemporalRamp*) other;
      returnFunction = new TemporalRamp(*other_cast);
    }
    else if (tag=="radial_power") {
      RadialPower * other_cast = (RadialPower*) other;
      returnFunction = new RadialPower(*other_cast);
    }
    pointerSet_.insert(returnFunction);
    return returnFunction;
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  ConstantFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  ConstantFunction::ConstantFunction(int narg, double* args) 
    : XT_Function(narg,args),
      C0(args[0])
  {
    tag_ = "constant";
  }
  //--------------------------------------------------------------------
  ConstantFunction::ConstantFunction(double arg) 
    : XT_Function(1,&arg),
      C0(arg)
  {
    tag_ = "constant";
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  LinearFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  LinearFunction::LinearFunction(int narg, double* args) 
    : XT_Function(narg,args)
  {
    C0 = args[6];
    tag_ = "linear";
    stringstream ss;
    ss << "created function : " << C0 << " + " << mask[0] << "(x-"<< x0[0] << ")+"<< mask[1] << "(y-"<<x0[1]<<")+"<<mask[2]<<"(z-"<<x0[2] << ")";
    ATC::LammpsInterface::instance()->print_msg_once(ss.str());
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  PiecewiseLinearFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  PiecewiseLinearFunction::PiecewiseLinearFunction(int narg, double* args) 
    : XT_Function(narg,args)
  {
    int i=0, idx = 6, n = (narg-idx)/2;
    xi.reset(n);
    fi.reset(n);
    while (idx < narg) {
      xi(i)   = args[idx++];
      fi(i++) = args[idx++];
    }
    tag_ = "piecewise_linear";
  }
  double PiecewiseLinearFunction::f(double * x, double t)
  {

    double s = mask[0]*(x[0]-x0[0])+mask[1]*(x[1]-x0[1])+mask[2]*(x[2]-x0[2]);
    int index = xi.index(s);

    if      (index < 0)              return fi(0);
    else if (index >= xi.size()-1 )  return fi(xi.size()-1);
    else {
      double f = fi(index) 
        + (fi(index+1)-fi(index))*(s-xi(index))/(xi(index+1)-xi(index));
      return f;
    }
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  LinearTemporalRamp
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  LinearTemporalRamp::LinearTemporalRamp(int narg, double* args)
    : XT_Function(narg,args)
  {
    double mask_final[3];
    mask_final[0] = args[6];
    mask_final[1] = args[7];
    mask_final[2] = args[8];
    C0_initial = args[9];
    double C0_final = args[10];
    double delta_t = args[11];

    for (int i = 0; i < 3; i++)
      mask_slope[i] = (mask_final[i] - mask[i])/delta_t;
    C0_slope = (C0_initial - C0_final)/delta_t;
  }

  double LinearTemporalRamp::f(double* x, double t) {
    double slope[3];
    for (int i = 0; i < 3; i++)
      slope[i] = mask[i] + mask_slope[i]*t;
    double C0 = C0_initial + C0_slope*t;
    return slope[0]*(x[0]-x0[0])+slope[1]*(x[1]-x0[1])+slope[2]*(x[2]-x0[2]) + C0;
  }

  double LinearTemporalRamp::dfdt(double* x, double t) {
    return mask_slope[0]*(x[0]-x0[0])+mask_slope[1]*(x[1]-x0[1])+mask_slope[2]*(x[2]-x0[2]) + C0_slope;
  }

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  QuadraticFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  QuadraticFunction::QuadraticFunction(int narg, double* args) 
    : XT_Function(narg,args)
  {
    C0 = args[6];
    C2[0] = args[7];
    C2[1] = args[8];
    C2[2] = args[9];
    C2[3] = args[10];
    C2[4] = args[11];
    C2[5] = args[12];
    tag_ = "quadratic";
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  SineFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  SineFunction::SineFunction(int narg, double* args) 
    : XT_Function(narg,args)
  {
    C = args[6];
    w = args[7];
    C0 = args[8];
    tag_ = "sine";
    stringstream ss;
    ss << "created function : " << C << " sin( " << mask[0] << "(x-"<< x0[0] << ")+"<< mask[1] << "(y-"<<x0[1]<<")+"<<mask[2]<<"(z-"<<x0[2] << ") - " << w << "t ) + " << C0;
    ATC::LammpsInterface::instance()->print_msg_once(ss.str());
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  GaussianFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  GaussianFunction::GaussianFunction(int narg, double* args) 
    : XT_Function(narg,args)
  {
    tau = args[6];
    C   = args[7];
    C0  = args[8];
    tag_ = "gaussian";
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  GaussianTemporalRamp
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  GaussianTemporalRamp::GaussianTemporalRamp(int narg, double* args) 
    : GaussianFunction(narg,args)
  {
    tau_initial = args[9];
    C_initial   = args[10];
    C0_initial  = args[11];
    double delta_t = args[12];

    tau_slope = (tau - tau_initial)/delta_t;
    C_slope   = (C - C_initial)/delta_t;
    C0_slope  = (C0 - C0_initial)/delta_t;

    tag_ = "gaussian_temporal_ramp";
  }
  double GaussianTemporalRamp::f(double* x, double t) {
    tau = tau_initial + tau_slope*t;
    C   = C_initial + C_slope*t;
    C0  = C0_initial + C0_slope*t;
    return GaussianFunction::f(x,t);
  }
  double GaussianTemporalRamp::dfdt(double* x, double t) {
    tau = tau_initial + tau_slope*t;
    C   = C_initial + C_slope*t;
    C0  = C0_initial + C0_slope*t;
    double dfdt = 0.;
    dfdt += C_slope*exp(-(mask[0]*(x[0]-x0[0])*(x[0]-x0[0])
                          +mask[1]*(x[1]-x0[1])*(x[1]-x0[1])
                          +mask[2]*(x[2]-x0[2])*(x[2]-x0[2]))
                        /tau/tau);
    dfdt += C*exp(2.*tau_slope*(mask[0]*(x[0]-x0[0])*(x[0]-x0[0])
                                +mask[1]*(x[1]-x0[1])*(x[1]-x0[1])
                                +mask[2]*(x[2]-x0[2])*(x[2]-x0[2]))
                  /tau/tau/tau);
    dfdt += C0_slope;
    return dfdt;
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  TemporalRamp
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  TemporalRamp::TemporalRamp(int narg, double* args) 
    : XT_Function(narg,args)
  {
    f_initial = args[0];
    double f_final = args[1];
    double delta_t = args[2];
    slope = (f_final - f_initial)/delta_t;
    tag_ = "temporal_ramp";
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  RadialPower
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  RadialPower::RadialPower(int narg, double* args) 
    : XT_Function(narg,args)
  {
    C0 = args[6];
    n = args[7];
    tag_ = "radial_power";
  }

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  InterpolationFunction
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  void InterpolationFunction::initialize(int npts, fstream &fileId, double coef)
  {  // read data
     npts_ = npts;
     xs_.reset(npts);
     fs_.reset(npts);
     fps_.reset(npts);
     double x,f,fp;
     int i = 0;
     while(fileId.good() && i < npts) {
       fileId >> x >> f >> fp;
       xs_(i)=x;
       fs_(i)=coef*f;
       fps_(i)=coef*fp;
       i++;
     }
     // scale tangents 
     double dx, dx0 = xs_(1)-xs_(0); 
     for (int i = 0; i < npts_ ; i++) {
       if      (i   == 0)     { dx = xs_(1)-xs_(0); }
       else if (i+1 == npts_) { dx = xs_(npts_-1)-xs_(npts_-2); }
       else { dx= 0.5*(xs_(i+1)-xs_(i-1)); }
       if (abs(dx-dx0) > 1.e-8) throw ATC_Error("InterpolationFunction::initialize non-uniform data spacing not handled currently");
       fps_(i) *= dx;
     }
     // options: calculate / adjust tangents for monotonicity

  }

  double InterpolationFunction::coordinate(double x, 
    double & f0, double & fp0, double & f1, double & fp1, double & inv_dx ) const
  {
    int i0 = xs_.index(x);
    int i1 = i0+1;
    if      (i0 < 0) {
      double x0 = xs_(0),   x1 = xs_(1);
      inv_dx = 1./(x1-x0); 
      fp0 = fp1 = fps_(0);
      f1 = fs_(0);
      f0 = fp0*(x-xs_(0))+f1;
      return 0;
    }
    else if (i1 >= npts_) {
      double x0 = xs_(npts_-2),   x1 = xs_(npts_-1);
      inv_dx = 1./(x1-x0); 
      fp0 = fp1 = fps_(i0);
      f0 = fs_(i0);
      f1 = fp0*(x-xs_(i0))+f0;
      return 1;
    }
    else {
      double x0 = xs_(i0),   x1 = xs_(i1);
      inv_dx = 1./(x1-x0); 
      f0  = fs_ (i0); f1  = fs_ (i1);
      fp0 = fps_(i0); fp1 = fps_(i1);
      double t = (x-x0)*inv_dx;
      return t;
    }
  }
  double InterpolationFunction::f(const double x) const
  {
    double f0,fp0,f1,fp1,inv_dx;
    double t = coordinate(x,f0,fp0,f1,fp1,inv_dx);
    double t2 = t*t; 
    double t3 = t*t2;
    double h00 = 2*t3 - 3*t2     + 1;
    double h10 =   t3 - 2*t2 + t;
    double h01 =-2*t3 + 3*t2;
    double h11 =   t3 -   t2;
    double f = h00 * f0 + h10 * fp0 + h01 * f1 + h11 * fp1;
    return f;
  }
  double InterpolationFunction::dfdt(const double x) const
  {
    double f0,fp0,f1,fp1,inv_dx;
    double t = coordinate(x,f0,fp0,f1,fp1,inv_dx);
    double t2 = t*t; 
    double d00 = 6*t2 - 6*t;
    double d10 = 3*t2 - 4*t + 1;
    double d01 =-6*t2 + 6*t;
    double d11 = 3*t2 - 2*t;
    double fp = d00 * f0 + d10 * fp0 + d01 * f1 + d11 * fp1;
    fp *= inv_dx;
    return fp;
  }
}

#include "XT_Function.h"
#include "ATC_Error.h"

namespace ATC {

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //  XT_Function
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  XT_Function::XT_Function(int narg, double* args) 
  {
// NOTE need to adjust scaling to match input nodal coordinates
    if (narg > 5 ) { // NOTE kludge for temporal only functions 
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
    for (int i = 0; i < xtFunctionVec_.size(); i++) {
      delete xtFunctionVec_[i];
    }
  }

  // add user function into the if statement and assign returnFunction to it
  XT_Function* XT_Function_Mgr::get_function(string & type, int nargs, double * args)
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
    else if (type=="quadratic")
      returnFunction = new QuadraticFunction(nargs,args);
    else if (type=="sine")
      returnFunction = new SineFunction(nargs,args);
    else if (type=="gaussian")
      returnFunction = new GaussianFunction(nargs,args);
    else if (type=="gaussian_temporal_ramp")
      returnFunction = new GaussianTemporalRamp(nargs,args);
    else
      throw ATC_Error(0,"Bad user function name");

    xtFunctionVec_.push_back(returnFunction);

    return returnFunction;
  }

  // add user function into the if statement and assign returnFunction to it
  XT_Function* XT_Function_Mgr::get_function(char ** args, int nargs)
  {
    string type = args[0];
    int narg = nargs -1;
    double dargs[narg];
    for (int i = 0; i < narg; ++i) dargs[i] = atof(args[i+1]);
  
    return get_function(type, narg, dargs);
  }

  // add constant function
  XT_Function* XT_Function_Mgr::get_constant_function(double c)
  {
    double args[1] = {c}; // NOTE kludge
    XT_Function * returnFunction = new ConstantFunction(1,args);
    xtFunctionVec_.push_back(returnFunction);
    return (returnFunction);
  }

  XT_Function* XT_Function_Mgr::copy_XT_function(XT_Function* other)
  {
    string tag = other->get_tag();
  
    XT_Function * returnFunction = NULL;
    if (tag=="linear") {
      LinearFunction * other_cast = (LinearFunction*) other;
      returnFunction = new LinearFunction(*other_cast);
    }
    if (tag=="quadratic") {
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
    tag = "constant";
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
    tag = "linear";
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
    tag = "quadratic";
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
    tag = "sine";
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
    tag = "gaussian";
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

    tag = "gaussian_temporal_ramp";
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
    tag = "temporal_ramp";
  }
}

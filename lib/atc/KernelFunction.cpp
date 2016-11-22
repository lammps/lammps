#include "KernelFunction.h"
#include "math.h"
#include <vector>
#include "ATC_Error.h"
#include "Quadrature.h"
#include "Utility.h"


using namespace std;
using namespace ATC_Utility;

static const double tol = 1.0e-8;

static const int line_ngauss = 10; 
static double line_xg[line_ngauss], line_wg[line_ngauss];

namespace ATC {

  
  //========================================================================
  //  KernelFunctionMgr
  //========================================================================
  KernelFunctionMgr * KernelFunctionMgr::myInstance_ = NULL;
  //------------------------------------------------------------------------
  //  instance
  //------------------------------------------------------------------------
  KernelFunctionMgr * KernelFunctionMgr::instance()
  {
    if (myInstance_ == NULL) {
      myInstance_ = new KernelFunctionMgr();
    }
    return myInstance_;
  }
  //------------------------------------------------------------------------
  // get function from args
  //------------------------------------------------------------------------
  KernelFunction* KernelFunctionMgr::function(char ** arg, int narg)
  {
    /*! \page man_kernel_function fix_modify AtC kernel 
      \section syntax
      fix_modify AtC kernel <type> <parameters>
      - type (keyword) = step, cell, cubic_bar, cubic_cylinder, cubic_sphere, 
                         quartic_bar, quartic_cylinder, quartic_sphere \n
      - parameters :\n
      step = radius (double) \n
      cell = hx, hy, hz (double) or h (double) \n
      cubic_bar = half-width (double) \n
      cubic_cylinder = radius (double) \n
      cubic_sphere = radius (double) \n
      quartic_bar = half-width (double) \n
      quartic_cylinder = radius (double) \n
      quartic_sphere = radius (double) \n
      \section examples
      fix_modify AtC kernel cell 1.0 1.0 1.0
      fix_modify AtC kernel quartic_sphere 10.0
      \section description
      
      \section restrictions
      Must be used with the hardy AtC fix \n
      For bar kernel types, half-width oriented along x-direction \n
      For cylinder kernel types, cylindrical axis is assumed to be in z-direction \n
      ( see \ref man_fix_atc )
      \section related
      \section default
      No default
    */
    int argIdx = 0;
    KernelFunction * ptr = NULL;
    char* type = arg[argIdx++];
    if (strcmp(type,"step")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff radius
      ptr = new KernelFunctionStep(1,parameters);
    }
    else if (strcmp(type,"cell")==0) {
      double parameters[3];
      parameters[0] = parameters[1] = parameters[2] = atof(arg[argIdx++]);
      if (narg > argIdx) { // L_x, L_y, L_z
        for (int i = 1; i < 3; i++) { parameters[i] = atof(arg[argIdx++]); }
      }
      ptr = new KernelFunctionCell(2,parameters);
    }
    else if (strcmp(type,"cubic_bar")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff half-length 
      ptr = new KernelFunctionCubicBar(1,parameters);
    }
    else if (strcmp(type,"linear_bar")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff half-length 
      ptr = new KernelFunctionLinearBar(1,parameters);
    }
    else if (strcmp(type,"cubic_cylinder")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff radius
      ptr = new KernelFunctionCubicCyl(1,parameters);
    }
    else if (strcmp(type,"cubic_sphere")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff radius
      ptr = new KernelFunctionCubicSphere(1,parameters);
    }
    else if (strcmp(type,"quartic_bar")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff half-length 
      ptr = new KernelFunctionQuarticBar(1,parameters);
    }
    else if (strcmp(type,"quartic_cylinder")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff radius
      ptr = new KernelFunctionQuarticCyl(1,parameters);
    }
    else if (strcmp(type,"quartic_sphere")==0) {
      double parameters[1] = {atof(arg[argIdx])}; // cutoff radius
      ptr = new KernelFunctionQuarticSphere(1,parameters);
    }
    pointerSet_.insert(ptr);
    return ptr;
  }

  // Destructor
  KernelFunctionMgr::~KernelFunctionMgr()
  {
     set<KernelFunction * >::iterator it;
    for (it = pointerSet_.begin(); it != pointerSet_.end(); it++)
      if (*it) delete *it;
  }

  //------------------------------------------------------------------------
  //  KernelFunction
  //------------------------------------------------------------------------
  // constructor
  KernelFunction::KernelFunction(int nparameters, double* parameters):
    Rc_(0),invRc_(0),nsd_(3),
    lammpsInterface_(LammpsInterface::instance())
  { 
    Rc_ = parameters[0];
    invRc_ = 1.0/Rc_;
    Rc_ = parameters[0];
    invRc_ = 1.0/Rc_;
    invVol_ = 1.0/(4.0/3.0*Pi_*pow(Rc_,3));

    ATC::Quadrature::instance()->set_line_quadrature(line_ngauss,line_xg,line_wg);

    // get periodicity and box bounds/lengths
    lammpsInterface_->box_periodicity(periodicity[0],
                                          periodicity[1],periodicity[2]);
    lammpsInterface_->box_bounds(box_bounds[0][0],box_bounds[1][0],
                                     box_bounds[0][1],box_bounds[1][1],
                                     box_bounds[0][2],box_bounds[1][2]);
    for (int k = 0; k < 3; k++) {
      box_length[k] = box_bounds[1][k] - box_bounds[0][k]; 
    } 
  }

  // does an input node's kernel intersect bonds on this processor
  bool KernelFunction::node_contributes(DENS_VEC node) const
  {
    DENS_VEC ghostNode = node;
    bool contributes = true;
    bool ghostContributes = lammpsInterface_->nperiodic();
    double kernel_bounds[2][3];
    lammpsInterface_->sub_bounds(kernel_bounds[0][0],kernel_bounds[1][0],
                                     kernel_bounds[0][1],kernel_bounds[1][1],
                                     kernel_bounds[0][2],kernel_bounds[1][2]);
    for (int i=0; i<3; ++i) {
      if (i < nsd_) {
        kernel_bounds[0][i] -= (Rc_+lammpsInterface_->pair_cutoff());
        kernel_bounds[1][i] += (Rc_+lammpsInterface_->pair_cutoff());
        contributes = contributes && (node(i) >= kernel_bounds[0][i]) 
                                  && (node(i) <  kernel_bounds[1][i]);
        if (periodicity[i]) {
          if (node[i] <= box_bounds[0][i] + box_length[i]/2) {
            ghostNode[i] += box_length[i];
          } else {
            ghostNode[i] -= box_length[i];
          }
          ghostContributes = ghostContributes 
                          && ((ghostNode(i) >= kernel_bounds[0][i]) || 
                                   (node(i) >= kernel_bounds[0][i]))
                          && ((ghostNode(i) <  kernel_bounds[1][i]) ||
                                   (node(i) <  kernel_bounds[1][i]));
        }
      }
      if (!(contributes || ghostContributes)) break;
    }
    return true;
  }

  bool KernelFunction::in_support(DENS_VEC dx) const
  {
    if (dx.norm() > Rc_) return false;
    return true;
  }

  // bond function value via quadrature
  double KernelFunction::bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2) const
  {
    DENS_VEC xab(nsd_), q(nsd_);
    double lamg;
    double bhsum=0.0;
    xab = xa - xb;
    for (int i = 0; i < line_ngauss; i++) {
      lamg=0.5*((lam2-lam1)*line_xg[i]+(lam2+lam1));
      q = lamg*xab + xb;
      double locg_value=this->value(q);
      bhsum+=locg_value*line_wg[i];
    }
    return 0.5*(lam2-lam1)*bhsum;
  }

  // localization-volume intercepts for bond calculation 
  // bond intercept values assuming spherical support
  void KernelFunction::bond_intercepts(DENS_VEC& xa,
                 DENS_VEC& xb, double &lam1, double &lam2) const
  {
    if (nsd_ == 2) {// for cylinders, axis is always z! 
      const int iaxis = 2;
      xa[iaxis] = 0.0;
      xb[iaxis] = 0.0;
    } else if (nsd_ == 1) {// for bars, analysis is 1D in x
      xa[1] = 0.0;
      xa[2] = 0.0;
      xb[1] = 0.0;
      xb[2] = 0.0;
    }
    lam1 = lam2 = -1;
    double ra_n = invRc_*xa.norm(); // lambda = 1
    double rb_n = invRc_*xb.norm(); // lambda = 0
    bool a_in = (ra_n <= 1.0);
    bool b_in = (rb_n <= 1.0);
    if (a_in && b_in) {
      lam1 = 0.0; 
      lam2 = 1.0; 
      return;
    }
    DENS_VEC xab = xa - xb;
    double rab_n = invRc_*xab.norm();
    double a = rab_n*rab_n; // always at least an interatomic distance
    double b = 2.0*invRc_*invRc_*xab.dot(xb);
    double c = rb_n*rb_n - 1.0;
    double discrim = b*b - 4.0*a*c; // discriminant
    if (discrim < 0) return; // no intersection
    // num recipes:
    double s1, s2;
    if (b < 0.0) {
      double aux = -0.5*(b-sqrt(discrim));
      s1 = c/aux;
      s2 = aux/a;
    } 
    else {
      double aux = -0.5*(b+sqrt(discrim));
      s1 = aux/a;
      s2 = c/aux;
    }
    if (a_in && !b_in) {
      lam1 = s1;
      lam2 = 1.0;
    } 
    else if (!a_in && b_in) { 
      lam1 = 0.0;
      lam2 = s2; 
    }
    else {
      if (s1 >= 0.0 && s2 <= 1.0) { 
        lam1 = s1; 
        lam2 = s2; 
      }
    }
  }

  //------------------------------------------------------------------------
  // constructor
  KernelFunctionStep::KernelFunctionStep
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters) 
  {
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
        }
      }
    } 
  }

  // function value
  double KernelFunctionStep::value(DENS_VEC& x_atom) const
  {
    double rn=invRc_*x_atom.norm();
    if (rn <= 1.0) { return 1.0; }
    else           { return 0.0; }
  }
  
  // function derivative value
  void KernelFunctionStep::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }

  //------------------------------------------------------------------------
  /** a step with rectangular support suitable for a rectangular grid */
  // constructor
  KernelFunctionCell::KernelFunctionCell
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters) 
  {
    hx = parameters[0];
    hy = parameters[1];
    hz = parameters[2];
    invVol_ = 1.0/8.0/(hx*hy*hz);
    cellBounds_.reset(6);
    cellBounds_(0) = -hx;
    cellBounds_(1) =  hx;
    cellBounds_(2) = -hy;
    cellBounds_(3) =  hy;
    cellBounds_(4) = -hz;
    cellBounds_(5) =  hz;
      
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (parameters[k] > 0.5*box_length[k]) {
          throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
        }
      }
    }
  }

  // does an input node's kernel intersect bonds on this processor
  bool KernelFunctionCell::node_contributes(DENS_VEC node) const
  {
    DENS_VEC ghostNode = node;
    bool contributes = true;
    bool ghostContributes = lammpsInterface_->nperiodic();
    double kernel_bounds[2][3];
    lammpsInterface_->sub_bounds(kernel_bounds[0][0],kernel_bounds[1][0],
                                     kernel_bounds[0][1],kernel_bounds[1][1],
                                     kernel_bounds[0][2],kernel_bounds[1][2]);
    for (int i=0; i<3; ++i) {
      kernel_bounds[0][i] -= (cellBounds_(i*2+1) +
                              lammpsInterface_->pair_cutoff());
      kernel_bounds[1][i] += (cellBounds_(i*2+1) + 
                              lammpsInterface_->pair_cutoff());
      contributes = contributes && (node(i) >= kernel_bounds[0][i]) 
                                && (node(i) <  kernel_bounds[1][i]);
      if (periodicity[i]) {
        if (node[i] <= box_bounds[0][i] + box_length[i]/2) {
          ghostNode[i] += box_length[i];
        } else {
          ghostNode[i] -= box_length[i];
        }
        ghostContributes = ghostContributes 
                        && ((ghostNode(i) >= kernel_bounds[0][i]) ||
                                 (node(i) >= kernel_bounds[0][i]))
                        && ((ghostNode(i) <  kernel_bounds[1][i]) ||
                                 (node(i) <  kernel_bounds[1][i]));
      }
      if (!(contributes || ghostContributes)) break;
    }
    return true;
  }

  bool KernelFunctionCell::in_support(DENS_VEC dx) const
  {
    if (dx(0) < cellBounds_(0)
     || dx(0) > cellBounds_(1)
     || dx(1) < cellBounds_(2)
     || dx(1) > cellBounds_(3)
     || dx(2) < cellBounds_(4)
     || dx(2) > cellBounds_(5) ) return false;
    return true;
  }

  // function value
  double KernelFunctionCell::value(DENS_VEC& x_atom) const
  {
    if ((cellBounds_(0) <= x_atom(0)) && (x_atom(0) < cellBounds_(1)) 
     && (cellBounds_(2) <= x_atom(1)) && (x_atom(1) < cellBounds_(3)) 
     && (cellBounds_(4) <= x_atom(2)) && (x_atom(2) < cellBounds_(5))) { 
      return 1.0; 
    } 
    else { 
      return 0.0; 
    }
  }

  // function derivative value
  void KernelFunctionCell::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }
 
  // bond intercept values for rectangular region : origin is the node position
  void KernelFunctionCell::bond_intercepts(DENS_VEC& xa,
                 DENS_VEC& xb, double &lam1, double &lam2) const
  {
    lam1 = 0.0; // start
    lam2 = 1.0; // end

    bool a_in = (value(xa) > 0.0);
    bool b_in = (value(xb) > 0.0);

    // (1) both in, no intersection needed
    if (a_in && b_in) {
      return;
    }
    // (2) for one in & one out -> one plane intersection
    // determine the points of intersection between the line joining
    // atoms a and b and the bounding planes of the localization volume
    else if (a_in || b_in) {
      DENS_VEC xab = xa - xb;
      for (int i = 0; i < nsd_; i++) {
        // check if segment is parallel to face
        if (fabs(xab(i)) > tol) {
          for (int j = 0; j < 2; j++) {
            double s = (cellBounds_(2*i+j) - xb(i))/xab(i);
            // check if between a & b
            if (s >= 0 && s <= 1) {
              bool in_bounds = false; 
              DENS_VEC x = xb + s*xab;
              if (i == 0) { 
                if ((cellBounds_(2) <= x(1)) && (x(1) <= cellBounds_(3))
                 && (cellBounds_(4) <= x(2)) && (x(2) <= cellBounds_(5))) {
                  in_bounds = true;
                }
              }
              else if (i == 1) { 
                if ((cellBounds_(0) <= x(0)) && (x(0) <= cellBounds_(1))
                 && (cellBounds_(4) <= x(2)) && (x(2) <= cellBounds_(5))) {
                  in_bounds = true;
                }
              }
              else if (i == 2) { 
                if ((cellBounds_(0) <= x(0)) && (x(0) <= cellBounds_(1))
                 && (cellBounds_(2) <= x(1)) && (x(1) <= cellBounds_(3))) {
                  in_bounds = true;
                }
              }
              if (in_bounds) {
                if (a_in) { lam1 = s;}
                else      { lam2 = s;}
                return;
              }
            }
          }
        } 
      }
      throw ATC_Error("logic failure in HardyKernel Cell for single intersection\n");
    }
    // (3) both out -> corner intersection
    else {
      lam2 = lam1; // default to no intersection
      DENS_VEC xab = xa - xb;
      double ss[6] = {-1,-1,-1,-1,-1,-1};
      int is = 0;
      for (int i = 0; i < nsd_; i++) {
        // check if segment is parallel to face
        if (fabs(xab(i)) > tol) {
          for (int j = 0; j < 2; j++) {
            double s = (cellBounds_(2*i+j) - xb(i))/xab(i);
            // check if between a & b
            if (s >= 0 && s <= 1) {
              // check if in face
              DENS_VEC x = xb + s*xab;
              if (i == 0) { 
                if ((cellBounds_(2) <= x(1)) && (x(1) <= cellBounds_(3))
                 && (cellBounds_(4) <= x(2)) && (x(2) <= cellBounds_(5))) {
                  ss[is++] = s;
                }
              }
              else if (i == 1) { 
                if ((cellBounds_(0) <= x(0)) && (x(0) <= cellBounds_(1))
                 && (cellBounds_(4) <= x(2)) && (x(2) <= cellBounds_(5))) {
                  ss[is++] = s;
                }
              }
              else if (i == 2) { 
                if ((cellBounds_(0) <= x(0)) && (x(0) <= cellBounds_(1))
                 && (cellBounds_(2) <= x(1)) && (x(1) <= cellBounds_(3))) {
                  ss[is++] = s;
                }
              }
            }
          }
        } 
      }
      if (is == 1) {
      // intersection occurs at a box edge - leave lam1 = lam2
      }
      else if (is == 2) {
        lam1 = min(ss[0],ss[1]);
        lam2 = max(ss[0],ss[1]);
      }
      else if (is == 3) {
      // intersection occurs at a box vertex - leave lam1 = lam2
      }
      else {
        if (is != 0) throw ATC_Error("logic failure in HardyKernel Cell for corner intersection\n");
      }
    }
  }

  //------------------------------------------------------------------------
  // constructor
  KernelFunctionCubicSphere::KernelFunctionCubicSphere
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters) 
  {
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double KernelFunctionCubicSphere::value(DENS_VEC& x_atom) const
  {
     double r=x_atom.norm();
     double rn=r/Rc_;
     if (rn < 1.0) { return 5.0*(1.0-3.0*rn*rn+2.0*rn*rn*rn); }
     else          { return 0.0; }
  } 

  // function derivative value
  void KernelFunctionCubicSphere::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }

  //------------------------------------------------------------------------
  // constructor
  KernelFunctionQuarticSphere::KernelFunctionQuarticSphere
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters) 
  {
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double KernelFunctionQuarticSphere::value(DENS_VEC& x_atom) const 
  {
     double r=x_atom.norm();
     double rn=r/Rc_;
     if (rn < 1.0) { return 35.0/8.0*pow((1.0-rn*rn),2); }
     else          { return 0.0; }
  } 

  // function derivative value
  void KernelFunctionQuarticSphere::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }

  //------------------------------------------------------------------------
  // constructor
  KernelFunctionCubicCyl::KernelFunctionCubicCyl
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters)
  {
    nsd_ = 2;
    double Lz = box_length[2];
    invVol_ = 1.0/(Pi_*pow(Rc_,2)*Lz);
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double KernelFunctionCubicCyl::value(DENS_VEC& x_atom) const
  {
     double r=sqrt(pow(x_atom(0),2)+pow(x_atom(1),2));
     double rn=r/Rc_;
     if (rn < 1.0) { return 10.0/3.0*(1.0-3.0*rn*rn+2.0*rn*rn*rn); }
     else          { return 0.0; }
  } 

  // function derivative value
  void KernelFunctionCubicCyl::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }

  //------------------------------------------------------------------------
  // constructor
  KernelFunctionQuarticCyl::KernelFunctionQuarticCyl
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters)
  {
    nsd_ = 2;
    double Lz = box_length[2];
    invVol_ = 1.0/(Pi_*pow(Rc_,2)*Lz);
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double KernelFunctionQuarticCyl::value(DENS_VEC& x_atom) const
  {
    double r=sqrt(pow(x_atom(0),2)+pow(x_atom(1),2));
    double rn=r/Rc_;
    if (rn < 1.0) { return 3.0*pow((1.0-rn*rn),2); }
    else          { return 0.0; }
  } 

  // function derivative value
  void KernelFunctionQuarticCyl::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }
  //------------------------------------------------------------------------
  // constructor
  KernelFunctionCubicBar::KernelFunctionCubicBar
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters)
  {
    // Note: Bar is assumed to be oriented in the x(0) direction
    nsd_ = 1;
    double Ly = box_length[1];
    double Lz = box_length[2];
    invVol_ = 1.0/(2*Rc_*Ly*Lz);
    if ((bool) periodicity[0]) {
      if (Rc_ > 0.5*box_length[0]) {
        throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
      };
    };
  }

  // function value
  double KernelFunctionCubicBar::value(DENS_VEC& x_atom) const
  {
     double r=abs(x_atom(0));
     double rn=r/Rc_;
     if (rn < 1.0) { return 2.0*(1.0-3.0*rn*rn+2.0*rn*rn*rn); }
     else          { return 0.0; }
  } 

  // function derivative value
  void KernelFunctionCubicBar::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }

  //------------------------------------------------------------------------
  // constructor
  KernelFunctionLinearBar::KernelFunctionLinearBar
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters)
  {
    // Note: Bar is assumed to be oriented in the z(0) direction
    double Lx = box_length[0];
    double Ly = box_length[1];
    invVol_ = 1.0/(Lx*Ly*Rc_);
    if ((bool) periodicity[2]) {
      if (Rc_ > 0.5*box_length[2]) {
        throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
      };
    };
  }

  // function value
  double KernelFunctionLinearBar::value(DENS_VEC& x_atom) const
  {
     double r=abs(x_atom(2));
     double rn=r/Rc_;
     if (rn < 1.0) { return 1.0-rn; }
     else          { return 0.0; }
  } 

  // function derivative value
  void KernelFunctionLinearBar::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
     deriv(0) = 0.0;
     deriv(1) = 0.0;
     double r=abs(x_atom(2));
     double rn=r/Rc_;
     if (rn < 1.0 && x_atom(2) <= 0.0) { deriv(2) = 1.0/Rc_; }
     else if (rn < 1.0 && x_atom(2) > 0.0) { deriv(2) = -1.0/Rc_; }
     else { deriv(2) = 0.0; }
  }

  //------------------------------------------------------------------------
  // constructor
  KernelFunctionQuarticBar::KernelFunctionQuarticBar
    (int nparameters, double* parameters): 
    KernelFunction(nparameters, parameters)
  {
    // Note: Bar is assumed to be oriented in the x(0) direction
    nsd_ = 1;
    double Ly = box_length[1];
    double Lz = box_length[2];
    invVol_ = 1.0/(2*Rc_*Ly*Lz);
    if ((bool) periodicity[0]) {
      if (Rc_ > 0.5*box_length[0]) {
        throw ATC_Error("Size of localization volume is too large for periodic boundary condition");
      };
    };
  }

  // function value
  double KernelFunctionQuarticBar::value(DENS_VEC& x_atom) const
  {
     double r=abs(x_atom(0));
     double rn=r/Rc_;
     // if (rn < 1.0) { return 5.0/2.0*(1.0-6*rn*rn+8*rn*rn*rn-3*rn*rn*rn*rn); } - alternative quartic
     if (rn < 1.0) { return 15.0/8.0*pow((1.0-rn*rn),2); }
     else          { return 0.0; }
  } 

  // function derivative value
  void KernelFunctionQuarticBar::derivative(const DENS_VEC& x_atom, DENS_VEC& deriv) const
  {
     deriv.reset(nsd_);
  }
};

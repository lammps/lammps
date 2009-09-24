#include "ATC_HardyKernel.h"
#include "math.h"
#include <iostream> //for debugging purposes; take this out after I'm done
#include <vector>
#include "ATC_Error.h"
#include "Quadrature.h"

using namespace std;

static const double Pi = 4.0*atan(1.0);
static const double tol = 1.0e-8;

namespace ATC {


  //------------------------------------------------------------------------
  // constructor
  ATC_HardyKernel::ATC_HardyKernel(int nparameters, double* parameters):
    lammpsInterface_(LammpsInterface::instance()),
    Rc_(0),invRc_(0),nsd_(3)
  { 
    Rc_ = parameters[0];
    invRc_ = 1.0/Rc_;
    Rc_ = parameters[0];
    invRc_ = 1.0/Rc_;
    invVol_ = 1.0/(4.0/3.0*Pi*pow(Rc_,3));

    set_line_quadrature(line_ngauss,line_xg,line_wg);

    // get periodicity and box bounds/lengths
    lammpsInterface_->get_box_periodicity(periodicity[0],
                                          periodicity[1],periodicity[2]);
    lammpsInterface_->get_box_bounds(box_bounds[0][0],box_bounds[1][0],
                                     box_bounds[0][1],box_bounds[1][1],
                                     box_bounds[0][2],box_bounds[1][2]);
    for (int k = 0; k < 3; k++) {
      box_length[k] = box_bounds[1][k] - box_bounds[0][k]; 
    } 
  };

  // bond function value via quadrature
  double ATC_HardyKernel::bond(DENS_VEC& xa, DENS_VEC&xb, double lam1, double lam2)
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
  void ATC_HardyKernel::bond_intercepts(DENS_VEC& xa,
                 DENS_VEC& xb, double &lam1, double &lam2) 
  {
    if (nsd_ == 2) {// for cylinders, axis is always z! 
      const int iaxis = 2;
      xa[iaxis] = 0.0;
      xb[iaxis] = 0.0;
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
  };

  //------------------------------------------------------------------------
  // constructor
  ATC_HardyKernelStep::ATC_HardyKernelStep
    (int nparameters, double* parameters): 
    ATC_HardyKernel(nparameters, parameters) 
  {
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error(0,"Size of localization volume is too large for periodic boundary condition");
        };
      };
    }; 
  }

  // function value
  double ATC_HardyKernelStep::value(DENS_VEC& x_atom) 
  {
    double rn=invRc_*x_atom.norm();
    if (rn <= 1.0) { return 1.0; }
    else           { return 0.0; }
  };

  //------------------------------------------------------------------------
  /** a step with rectangular support suitable for a rectangular grid */
  // constructor
  ATC_HardyKernelCell::ATC_HardyKernelCell
    (int nparameters, double* parameters): 
    ATC_HardyKernel(nparameters, parameters) 
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
          throw ATC_Error(0,"Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double ATC_HardyKernelCell::value(DENS_VEC& x_atom) 
  {
    if ((cellBounds_(0) <= x_atom(0)) && (x_atom(0) < cellBounds_(1)) 
     && (cellBounds_(2) <= x_atom(1)) && (x_atom(1) < cellBounds_(3)) 
     && (cellBounds_(4) <= x_atom(2)) && (x_atom(2) < cellBounds_(5))) { 
      return 1.0; 
    } 
    else { 
      return 0.0; 
    }
  };

  // bond intercept values for rectangular region : origin is the node position
  void ATC_HardyKernelCell::bond_intercepts(DENS_VEC& xa,
                 DENS_VEC& xb, double &lam1, double &lam2) 
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
      throw ATC_Error(0,"logic failure in HardyKernel Cell for single intersection\n");
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
        if (is != 0) throw ATC_Error(0,"logic failure in HardyKernel Cell for corner intersection\n");
      }
    }
  }

  //------------------------------------------------------------------------
  // constructor
  ATC_HardyKernelCubicSphere::ATC_HardyKernelCubicSphere
    (int nparameters, double* parameters): 
    ATC_HardyKernel(nparameters, parameters) 
  {
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error(0,"Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double ATC_HardyKernelCubicSphere::value(DENS_VEC& x_atom) 
  {
     double r=x_atom.norm();
     double rn=r/Rc_;
     if (rn < 1.0) { return 5.0*(1.0-3.0*rn*rn+2.0*rn*rn*rn); }
     else          { return 0.0; }
  } 
  //------------------------------------------------------------------------
  // constructor
  ATC_HardyKernelQuarticSphere::ATC_HardyKernelQuarticSphere
    (int nparameters, double* parameters): 
    ATC_HardyKernel(nparameters, parameters) 
  {
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error(0,"Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double ATC_HardyKernelQuarticSphere::value(DENS_VEC& x_atom) 
  {
     double r=x_atom.norm();
     double rn=r/Rc_;
     if (rn < 1.0) { return 35.0/8.0*pow((1.0-rn*rn),2); }
     else          { return 0.0; }
  } 
  //------------------------------------------------------------------------
  // constructor
  ATC_HardyKernelCubicCyl::ATC_HardyKernelCubicCyl
    (int nparameters, double* parameters): 
    ATC_HardyKernel(nparameters, parameters)
  {
    nsd_ = 2;
    double Lz = box_length[2];
    invVol_ = 1.0/(Pi*pow(Rc_,2)*Lz);
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error(0,"Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double ATC_HardyKernelCubicCyl::value(DENS_VEC& x_atom) 
  {
     double r=sqrt(pow(x_atom(0),2)+pow(x_atom(1),2));
     double rn=r/Rc_;
     if (rn < 1.0) { return 10.0/3.0*(1.0-3.0*rn*rn+2.0*rn*rn*rn); }
     else          { return 0.0; }
  } 


  //------------------------------------------------------------------------
  // constructor
  ATC_HardyKernelQuarticCyl::ATC_HardyKernelQuarticCyl
    (int nparameters, double* parameters): 
    ATC_HardyKernel(nparameters, parameters)
  {
    nsd_ = 2;
    double Lz = box_length[2];
    invVol_ = 1.0/(Pi*pow(Rc_,2)*Lz);
    for (int k = 0; k < nsd_; k++ ) {
      if ((bool) periodicity[k]) {
        if (Rc_ > 0.5*box_length[k]) {
          throw ATC_Error(0,"Size of localization volume is too large for periodic boundary condition");
        };
      };
    };
  }

  // function value
  double ATC_HardyKernelQuarticCyl::value(DENS_VEC& x_atom) 
  {
    double r=sqrt(pow(x_atom(0),2)+pow(x_atom(1),2));
    double rn=r/Rc_;
    if (rn < 1.0) { return 3.0*pow((1.0-rn*rn),2); }
    else          { return 0.0; }
  } 
};

// ATC header files
#include "ATC_Error.h"
#include "FE_Mesh.h"
#include "FE_Element.h"
#include "FE_Interpolate.h"
#include "LinearSolver.h"
#include "PolynomialSolver.h"
#include "Utility.h"

// Other headers
#include "math.h"

using ATC_Utility::dbl_geq;
using ATC_Utility::det3;
using std::vector;

namespace ATC {

static const int localCoordinatesMaxIterations_ = 40;
static const double localCoordinatesTolerance = 1.e-09;


  // =============================================================
  //   class FE_Element
  // =============================================================
  FE_Element::FE_Element(const int nSD,
                         int numFaces,
                         int numNodes,
                         int numFaceNodes,
                         int numNodes1d)
    : nSD_(nSD),
      numFaces_(numFaces),
      numNodes_(numNodes),
      numFaceNodes_(numFaceNodes),
      numNodes1d_(numNodes1d),
      tolerance_(localCoordinatesTolerance),
      projectionGuess_(COORDINATE_ALIGNED)
  {
    feInterpolate_ = NULL;
  }

  FE_Element::~FE_Element()
  {
    if (feInterpolate_) delete feInterpolate_;
  }

  int FE_Element::num_ips() const
  { 
    return feInterpolate_->num_ips(); 
  }

  int FE_Element::num_face_ips() const
  { 
    return feInterpolate_->num_face_ips(); 
  }

  void FE_Element::face_coordinates(const DENS_MAT &eltCoords, 
                                    const int faceID,
                                    DENS_MAT & faceCoords) const
  {
    faceCoords.reset(nSD_, numFaceNodes_);

    for (int inode=0; inode < numFaceNodes_; inode++) 
    {
      int id = localFaceConn_(faceID,inode);
      for (int isd=0; isd<nSD_; isd++) {
        faceCoords(isd,inode) = eltCoords(isd,id);
      }
    }
  }

  void FE_Element::mapping(const int inode, vector<int> &mapping) const
  {
    for (int iSD=0; iSD<nSD_; ++iSD) {
      mapping[iSD] = static_cast<int>((localCoords_(iSD,inode)+1)/2*
                                      (numNodes1d_-1));
    }
  }

  DENS_VEC FE_Element::local_coords_1d() const
  {
    DENS_VEC localCoords1d(numNodes1d_);
    for (int inode1d=0; inode1d<numNodes1d_; ++inode1d) {
       localCoords1d(inode1d) = (double(inode1d)/double(numNodes1d_-1))*2 - 1;
    }
    return localCoords1d;
  }

  void FE_Element::centroid(const DENS_MAT &eltCoords,
                                  DENS_VEC &centroid) const
  {
    centroid.reset(nSD_);
    for (int i = 0; i < nSD_; i++) {
      centroid(i) = eltCoords.row_mean(i);
    }
  }
 
  // -------------------------------------------------------------
  //  generic conversion from global to local coordinates using
  //  Newton's method to solve the nonliear equation that arises
  // -------------------------------------------------------------
  bool FE_Element::local_coordinates(const DENS_MAT &eltCoords,
                                     const DENS_VEC &x,
                                     DENS_VEC &xi) const
  {
    // initial guess
    DENS_VEC xiGuess(nSD_);
    this->initial_local_coordinates(eltCoords,x,xiGuess);
    // clip out-of-range guesses
    if (fabs(xiGuess(0)) > 1.0) xiGuess(0) = 0.;
    if (fabs(xiGuess(1)) > 1.0) xiGuess(1) = 0.;
    if (fabs(xiGuess(2)) > 1.0) xiGuess(2) = 0.;
    
    // iteratively solve the equation by calculating the global 
    // position of the guess and bringing the difference between it
    // and the actual global position of x to zero
    //
    // uses Newton's method
    DENS_VEC N(numNodes_);
    DENS_MAT dNdr(nSD_,numNodes_);

    DENS_VEC xGuess(nSD_);
    DENS_VEC xDiff(nSD_);
    DENS_MAT eltCoordsT = transpose(eltCoords);
    int count = 0;
    bool converged = false;
    while (count < localCoordinatesMaxIterations_) {
      feInterpolate_->compute_N_dNdr(xiGuess,N,dNdr);
      xGuess = N*eltCoordsT;
      xDiff  = xGuess-x; 
      // determine if the guess is close enough.
      // if it is, take it and run
      // if not, use Newton's method to update the guess
      if (!dbl_geq(abs(xDiff(0)),tolerance_) && 
          !dbl_geq(abs(xDiff(1)),tolerance_) && 
          !dbl_geq(abs(xDiff(2)),tolerance_)) {
        converged = true;
        xi = xiGuess;
        break;
      } else {
        xiGuess = xiGuess - transpose(inv(dNdr*eltCoordsT))*xDiff;
      }
      ++count;
    }
    return converged;
  }

  // -------------------------------------------------------------
  //  guess for initial local coordinates
  // -------------------------------------------------------------
  void FE_Element::initial_local_coordinates(const DENS_MAT &eltCoords,
                                             const DENS_VEC &x,
                                             DENS_VEC &xi) const
  {
    xi.reset(nSD_);
    if      (projectionGuess_ == COORDINATE_ALIGNED) {
      double min=0; double max=0;
      for (int i=0; i<nSD_; ++i) {
        bounds_in_dim(eltCoords,i,min,max);
        xi(i) = 2.0*(x(i)-min)/(max-min) - 1.0;
      }
    }
    else if (projectionGuess_ == CENTROID_LINEARIZED) { 
      DENS_VEC xi0(nSD_); xi0 = 0;
      DENS_VEC x0(nSD_), dx(nSD_);
      centroid(eltCoords,x0); 
      dx = x - x0; 
      vector<DENS_VEC> ts; ts.reserve(nSD_);
      tangents(eltCoords,xi0,ts);
      DENS_VEC & t1 = ts[0];
      DENS_VEC & t2 = ts[1];
      DENS_VEC & t3 = ts[2];
      double J = det3(t1.ptr(),t2.ptr(),t3.ptr());
      double J1 = det3(dx.ptr(),t2.ptr(),t3.ptr());
      double J2 = det3(t1.ptr(),dx.ptr(),t3.ptr());
      double J3 = det3(t1.ptr(),t2.ptr(),dx.ptr());
      xi(0) = J1/J;
      xi(1) = J2/J;
      xi(2) = J3/J;
    }
    else if (projectionGuess_ == TWOD_ANALYTIC) { 
      // assume x-y planar and HEX8
      double x0 = x(0), y0 = x(1);
      double X[4] = {eltCoords(0,0),eltCoords(0,1),eltCoords(0,2),eltCoords(0,3)};
      double Y[4] = {eltCoords(1,0),eltCoords(1,1),eltCoords(1,2),eltCoords(1,3)};
      double c[3]={0,0,0};
      c[0] = y0*X[1] - y0*X[2] - y0*X[3] + y0*X[4] - x0*Y[1] + (X[2]*Y[1])*0.5 + (X[3]*Y[1])*0.5 + x0*Y[2] - (X[1]*Y[2])*0.5 - (X[4]*Y[2])*0.5 + x0*Y[3] - (X[1]*Y[3])*0.5 - (X[4]*Y[3])*0.5 - x0*Y[4] + (X[2]*Y[4])*0.5 + (X[3]*Y[4])*0.5;
      c[1] = -(y0*X[1]) + y0*X[2] - y0*X[3] + y0*X[4] + x0*Y[1] - X[2]*Y[1] - x0*Y[2] + X[1]*Y[2] + x0*Y[3] - X[4]*Y[3] - x0*Y[4] + X[3]*Y[4];
      c[2] = (X[2]*Y[1])*0.5 - (X[3]*Y[1])*0.5 - (X[1]*Y[2])*0.5 + (X[4]*Y[2])*0.5 + (X[1]*Y[3])*0.5 - (X[4]*Y[3])*0.5 - (X[2]*Y[4])*0.5 + (X[3]*Y[4])*0.5;
      double xi2[2]={0,0};
      int nroots = solve_quadratic(c,xi2);
      if (nroots == 0) throw ATC_Error("no real roots in 2D analytic projection guess");
      double xi1[2]={0,0};
      xi1[0] = (4*x0 - X[1] + xi2[0]*X[1] - X[2] + xi2[0]*X[2] - X[3] - xi2[0]*X[3] - X[4] - xi2[0]*X[4])/(-X[1] + xi2[0]*X[1] + X[2] - xi2[0]*X[2] + X[3] + xi2[0]*X[3] - X[4] - xi2[0]*X[4]);
      xi1[1] = (4*x0 - X[1] + xi2[1]*X[1] - X[2] + xi2[1]*X[2] - X[3] - xi2[1]*X[3] - X[4] - xi2[1]*X[4])/(-X[1] + xi2[1]*X[1] + X[2] - xi2[1]*X[2] + X[3] + xi2[1]*X[3] - X[4] - xi2[1]*X[4]);
      // choose which one gives back x
      xi(0) = xi1[0]; 
      xi(1) = xi2[0]; 
      xi(2) = 0; 
    }
  }

 
  bool FE_Element::range_check(const DENS_MAT &eltCoords,
                               const DENS_VEC &x) const
  {
    double min=0; double max=0;
    for (int i=0; i<nSD_; ++i) {
      bounds_in_dim(eltCoords,i,min,max);
      if (!dbl_geq(x(i),min) || !dbl_geq(max,x(i))) return false;
    }
    return true;
  }


  // -------------------------------------------------------------
  //   Note: Only works for convex elements with planar faces with 
  //   outward normals
  // -------------------------------------------------------------
  bool FE_Element::contains_point(const DENS_MAT &eltCoords,
                                  const DENS_VEC &x) const
  {
    if (! range_check(eltCoords,x) ) return false;
    DENS_MAT faceCoords;
    DENS_VEC normal;
    normal.reset(nSD_);
    DENS_VEC faceToPoint;
    double dot;
    bool inside = true;
    for (int faceID=0; faceID<numFaces_; ++faceID) {
      face_coordinates(eltCoords, faceID, faceCoords);
      feInterpolate_->face_normal(faceCoords, 
                                  0,
                                  normal);
      faceToPoint = x - column(faceCoords, 0);
      dot = normal.dot(faceToPoint);
      if (dbl_geq(dot,0)) {
        inside = false;
        break;
      }
    }
    return inside;
  }

  // -------------------------------------------------------------
  //   returns the minimum and maximum values of an element in the
  //   specified dimension
  // -------------------------------------------------------------
  void FE_Element::bounds_in_dim(const DENS_MAT &eltCoords, const int dim,
                                 double &min, double &max) const
  {
    int it;
    // iterate over all nodes
    min = eltCoords(dim,0);
    it = 1;
    while (it < numNodes_) {
      if (dbl_geq(min,eltCoords(dim,it))) {
        if (dbl_geq(eltCoords(dim,it),min)) {
          ++it;
        } else {
          // set min to this node's coord in the specified dim, if it's 
          // smaller than the value previously stored 
          min = eltCoords(dim,it);
        }
      } else {
        ++it;
      }
    }
    max = eltCoords(dim,0);
    it = 1;
    while (it < numNodes_) {
      if (dbl_geq(max,eltCoords(dim,it))) {
        ++it;
      } else {
        // same, except max/larger
        max = eltCoords(dim,it);
      }
    }
  }
  
  // -------------------------------------------------------------
  //   shape_function calls should stay generic at all costs
  // -------------------------------------------------------------
  void FE_Element::shape_function(const VECTOR &xi,
                                  DENS_VEC &N) const
  {
    feInterpolate_->shape_function(xi, N);
  }

  void FE_Element::shape_function(const DENS_MAT eltCoords,
                                  const VECTOR &x,
                                  DENS_VEC &N)
  {
    DENS_VEC xi;
    local_coordinates(eltCoords, x, xi);

    feInterpolate_->shape_function(xi, N);
  }

  void FE_Element::shape_function(const DENS_MAT eltCoords,
                                  const VECTOR &x,
                                  DENS_VEC &N,
                                  DENS_MAT &dNdx)
  {
    DENS_VEC xi;
    local_coordinates(eltCoords, x, xi);

    feInterpolate_->shape_function(eltCoords, xi, N, dNdx);
  }  

  void FE_Element::shape_function_derivatives(const DENS_MAT eltCoords,
                                  const VECTOR &x,
                                  DENS_MAT &dNdx)
  {
    DENS_VEC xi;
    local_coordinates(eltCoords, x, xi);

    feInterpolate_->shape_function_derivatives(eltCoords, xi, dNdx);
  }  


  void FE_Element::shape_function(const DENS_MAT eltCoords,
                                  DENS_MAT &N,
                                  vector<DENS_MAT> &dN,
                                  DIAG_MAT &weights)
  {
    feInterpolate_->shape_function(eltCoords, N, dN, weights);
  }

  void FE_Element::face_shape_function(const DENS_MAT &eltCoords,
                                       const int faceID,
                                       DENS_MAT &N,
                                       DENS_MAT &n,
                                       DIAG_MAT &weights)
  {
    DENS_MAT faceCoords;
    face_coordinates(eltCoords, faceID, faceCoords); 
    
    feInterpolate_->face_shape_function(eltCoords, faceCoords, faceID,
                                        N, n, weights);
  }

  void FE_Element::face_shape_function(const DENS_MAT &eltCoords,
                                       const int faceID,
                                       DENS_MAT &N,
                                       vector<DENS_MAT> &dN,
                                       vector<DENS_MAT> &Nn,
                                       DIAG_MAT &weights)
  {
    DENS_MAT faceCoords;
    face_coordinates(eltCoords, faceID, faceCoords); 
    
    feInterpolate_->face_shape_function(eltCoords, faceCoords, faceID,
                                        N, dN, Nn, weights);
  }

  double FE_Element::face_normal(const DENS_MAT &eltCoords,
                                 const int faceID,
                                 int ip,
                                 DENS_VEC &normal) 
  {
    DENS_MAT faceCoords;
    face_coordinates(eltCoords, faceID, faceCoords);

    double J = feInterpolate_->face_normal(faceCoords, ip, normal);
    return J;
  }

  void FE_Element::tangents(const DENS_MAT &eltCoords,
                            const DENS_VEC & localCoords,
                            vector<DENS_VEC> &tangents,
                            const bool normalize)  const
  {
    feInterpolate_->tangents(eltCoords,localCoords,tangents,normalize);
  }
  
  
  // =============================================================
  //   class FE_ElementHex
  // =============================================================
  FE_ElementHex::FE_ElementHex(int numNodes, 
                               int numFaceNodes,
                               int numNodes1d)
    : FE_Element(3,                 // number of spatial dimensions
                 6,                 // number of faces
                 numNodes,
                 numFaceNodes,
                 numNodes1d)
  {
    //       3 --- 2
    //      /|    /|
    //     / 0 --/ 1     y
    //    7 --- 6 /      |
    //    |     |/       |
    //    4 --- 5         -----x
    //                  /
    //                 /
    //                z
  
    // Basic properties of element:
    vol_ = 8.0;
    faceArea_ = 4.0;
    
    // Order-specific information:
    if (numNodes != 8  &&
        numNodes != 20 &&
        numNodes != 27) {
      throw ATC_Error("Unrecognized interpolation order specified " 
                      "for element class: \n"                       
                      "  element only knows how to construct lin "  
                        "and quad elments.");
    }

    localCoords_.resize(nSD_,numNodes_);
    localFaceConn_ = Array2D<int>(numFaces_,numFaceNodes_);
    
    // Matrix of local nodal coordinates
    localCoords_(0,0) = -1; localCoords_(0,4) = -1; 
    localCoords_(1,0) = -1; localCoords_(1,4) = -1; 
    localCoords_(2,0) = -1; localCoords_(2,4) =  1;
    // 
    localCoords_(0,1) =  1; localCoords_(0,5) =  1; 
    localCoords_(1,1) = -1; localCoords_(1,5) = -1; 
    localCoords_(2,1) = -1; localCoords_(2,5) =  1;
    //
    localCoords_(0,2) =  1; localCoords_(0,6) =  1; 
    localCoords_(1,2) =  1; localCoords_(1,6) =  1; 
    localCoords_(2,2) = -1; localCoords_(2,6) =  1;
    //
    localCoords_(0,3) = -1; localCoords_(0,7) = -1; 
    localCoords_(1,3) =  1; localCoords_(1,7) =  1; 
    localCoords_(2,3) = -1; localCoords_(2,7) =  1;
    if (numNodes >= 20) {
      // only for quads
      localCoords_(0,8) =  0;  localCoords_(0,14) =  1; 
      localCoords_(1,8) = -1;  localCoords_(1,14) =  1; 
      localCoords_(2,8) = -1;  localCoords_(2,14) =  0;
      // 
      localCoords_(0,9) =  1;  localCoords_(0,15) = -1; 
      localCoords_(1,9) =  0;  localCoords_(1,15) =  1; 
      localCoords_(2,9) = -1;  localCoords_(2,15) =  0;
      //
      localCoords_(0,10) =  0; localCoords_(0,16) =  0; 
      localCoords_(1,10) =  1; localCoords_(1,16) = -1; 
      localCoords_(2,10) = -1; localCoords_(2,16) =  1;
      //
      localCoords_(0,11) = -1; localCoords_(0,17) =  1; 
      localCoords_(1,11) =  0; localCoords_(1,17) =  0; 
      localCoords_(2,11) = -1; localCoords_(2,17) =  1;
      //
      localCoords_(0,12) = -1; localCoords_(0,18) =  0; 
      localCoords_(1,12) = -1; localCoords_(1,18) =  1; 
      localCoords_(2,12) =  0; localCoords_(2,18) =  1;
      //
      localCoords_(0,13) =  1; localCoords_(0,19) = -1; 
      localCoords_(1,13) = -1; localCoords_(1,19) =  0; 
      localCoords_(2,13) =  0; localCoords_(2,19) =  1;
      if (numNodes >= 27) {
        // only for quads
        localCoords_(0,20) =  0; localCoords_(0,24) =  1; 
        localCoords_(1,20) =  0; localCoords_(1,24) =  0; 
        localCoords_(2,20) =  0; localCoords_(2,24) =  0;
        //
        localCoords_(0,21) =  0; localCoords_(0,25) =  0; 
        localCoords_(1,21) =  0; localCoords_(1,25) = -1; 
        localCoords_(2,21) = -1; localCoords_(2,25) =  0;
        //
        localCoords_(0,22) =  0; localCoords_(0,26) =  0; 
        localCoords_(1,22) =  0; localCoords_(1,26) =  1; 
        localCoords_(2,22) =  1; localCoords_(2,26) =  0;
        //
        localCoords_(0,23) = -1; 
        localCoords_(1,23) =  0; 
        localCoords_(2,23) =  0;
      }
    }

    // Matrix of local face connectivity
    // -x                    // +x
    localFaceConn_(0,0) = 0; localFaceConn_(1,0) = 1; 
    localFaceConn_(0,1) = 4; localFaceConn_(1,1) = 2; 
    localFaceConn_(0,2) = 7; localFaceConn_(1,2) = 6; 
    localFaceConn_(0,3) = 3; localFaceConn_(1,3) = 5; 
    if (numNodes >= 20) {
      localFaceConn_(0,4) = 12; localFaceConn_(1,4) = 9; 
      localFaceConn_(0,5) = 19; localFaceConn_(1,5) = 14; 
      localFaceConn_(0,6) = 15; localFaceConn_(1,6) = 17; 
      localFaceConn_(0,7) = 11; localFaceConn_(1,7) = 13; 
      if (numNodes >= 27) {
        localFaceConn_(0,8) = 23; localFaceConn_(1,8) = 24;
      }
    }

    // -y                    // +y
    localFaceConn_(2,0) = 0; localFaceConn_(3,0) = 3; 
    localFaceConn_(2,1) = 1; localFaceConn_(3,1) = 7; 
    localFaceConn_(2,2) = 5; localFaceConn_(3,2) = 6; 
    localFaceConn_(2,3) = 4; localFaceConn_(3,3) = 2; 
    if (numNodes >= 20) {
      localFaceConn_(2,4) = 8;  localFaceConn_(3,4) = 15; 
      localFaceConn_(2,5) = 13; localFaceConn_(3,5) = 18; 
      localFaceConn_(2,6) = 16; localFaceConn_(3,6) = 14; 
      localFaceConn_(2,7) = 12; localFaceConn_(3,7) = 10; 
      if (numNodes >= 27) {
        localFaceConn_(2,8) = 25; localFaceConn_(3,8) = 26;
      }
    }
 
    // -z                    // +z
    localFaceConn_(4,0) = 0; localFaceConn_(5,0) = 4; 
    localFaceConn_(4,1) = 3; localFaceConn_(5,1) = 5; 
    localFaceConn_(4,2) = 2; localFaceConn_(5,2) = 6; 
    localFaceConn_(4,3) = 1; localFaceConn_(5,3) = 7; 
    if (numNodes >= 20) {
      localFaceConn_(4,4) = 8;  localFaceConn_(5,4) = 16; 
      localFaceConn_(4,5) = 11; localFaceConn_(5,5) = 17; 
      localFaceConn_(4,6) = 10; localFaceConn_(5,6) = 18; 
      localFaceConn_(4,7) = 9;  localFaceConn_(5,7) = 19; 
      if (numNodes >= 27) {
        localFaceConn_(4,8) = 21; localFaceConn_(5,8) = 22;
      }
    }

    if (numNodes == 8) {
      feInterpolate_ = new FE_InterpolateCartLin(this);
    } else if (numNodes == 20) {
      feInterpolate_ = new FE_InterpolateCartSerendipity(this);
    } else if (numNodes == 27) {
      feInterpolate_ = new FE_InterpolateCartLagrange(this);
    }

    // determine alignment and skewness to see which guess we should use
    
  }

  FE_ElementHex::~FE_ElementHex()
  {
    // Handled by base class
  }

  void FE_ElementHex::set_quadrature(FeIntQuadrature type) 
  { 
    feInterpolate_->set_quadrature(HEXA,type);
  }
  
  bool FE_ElementHex::contains_point(const DENS_MAT &eltCoords,
                                     const DENS_VEC &x) const
  {
    if (! range_check(eltCoords,x) ) return false;

    DENS_VEC xi;
    bool converged = local_coordinates(eltCoords,x,xi); 
    if (!converged) return false;
    for (int i=0; i<nSD_; ++i) {
      if (!dbl_geq(1.0,abs(xi(i)))) return false;
    }
    return true;
  }

  // =============================================================
  //   class FE_ElementRect
  // =============================================================
  FE_ElementRect::FE_ElementRect()
    : FE_ElementHex(8,4,2)
  {
    // Handled by hex class
  }

  FE_ElementRect::~FE_ElementRect()
  {
    // Handled by base class
  }

  bool FE_ElementRect::contains_point(const DENS_MAT &eltCoords,
                                      const DENS_VEC &x) const
  {
    return range_check(eltCoords,x);
  }

  // much faster than the unstructured method
  bool FE_ElementRect::local_coordinates(const DENS_MAT &eltCoords,
                                         const DENS_VEC &x,
                                         DENS_VEC &xi) const
  {
    xi.reset(nSD_);
    double min = 0.0;
    double max = 0.0;
    for (int iSD=0; iSD<nSD_; ++iSD) {
      min = eltCoords(iSD,0);
      max = eltCoords(iSD,6);
      xi(iSD) = 2.0*(x(iSD)-min)/(max-min) - 1.0;
    }
    return true;
  }
  

  // =============================================================
  //   class FE_ElementTet
  // =============================================================
  FE_ElementTet::FE_ElementTet(int numNodes,
                               int numFaceNodes,
                               int numNodes1d)
    : FE_Element(3,                 // number of spatial dimensions
                 4,                 // number of faces
                 numNodes,
                 numFaceNodes,
                 numNodes1d)
  {
    // t
    // ^
    // |
    // |
    //                          s
    // 3                     .7
    // |\```---           .
    // |  \     ```--2 .
    // |    \       .| 
    // |      \  .   |
    // |      . \    |
    // |   .      \  |
    // |.___________\|  -------> r
    // 0             1  
    //
    // (This is as dictated by the EXODUSII standard.)
    //
    // The face opposite point 1 has r = 0,
    //                         2 has s = 0,
    //                         3 has t = 0,
    //                         0 has u = 0.

    // Basic properties of element:
    vol_ = 1.0/6.0; // local volume
    faceArea_ = 1.0/2.0;
 
    // Order-specific information:
    if (numNodes != 4 && numNodes != 10) {
      throw ATC_Error("Unrecognized interpolation order specified " 
                      "for element class: \n"                       
                      "  element only knows how to construct lin "  
                        "and quad elments.");
    }
    
    localCoords_.resize(nSD_+1, numNodes_);
    localFaceConn_ = Array2D<int>(numFaces_,numFaceNodes_);
    
    // Matrix of local nodal coordinates
    //
    // Remember, there's actually another coordinate too (u), coming
    // out from the third non-normal face. But we can deal with it on
    // the fly; these coordinates are barycentric, such that
    // r + s + t + u = 1 always. As such we only need to deal with r,
    // s, and t and the computations become easy.
    //
    // The first three axes correspond to x, y, and z (essentially),
    // for the canonical element.
    
    // Everyone gets these nodes...
    localCoords_(0,0) = 0; localCoords_(0,2) = 0;
    localCoords_(1,0) = 0; localCoords_(1,2) = 1; 
    localCoords_(2,0) = 0; localCoords_(2,2) = 0;
    localCoords_(3,0) = 1; localCoords_(3,2) = 0;
    //
    localCoords_(0,1) = 1; localCoords_(0,3) = 0; 
    localCoords_(1,1) = 0; localCoords_(1,3) = 0; 
    localCoords_(2,1) = 0; localCoords_(2,3) = 1; 
    localCoords_(3,1) = 0; localCoords_(3,3) = 0;
    if (numNodes >= 10) {
      // ...quads get even more!
      localCoords_(0,4) = 0.5; localCoords_(0,5) = 0.5;
      localCoords_(1,4) = 0.0; localCoords_(1,5) = 0.5; 
      localCoords_(2,4) = 0.0; localCoords_(2,5) = 0.0;
      localCoords_(3,4) = 0.5; localCoords_(3,5) = 0.0;
      //
      localCoords_(0,6) = 0.0; localCoords_(0,7) = 0.0; 
      localCoords_(1,6) = 0.5; localCoords_(1,7) = 0.0; 
      localCoords_(2,6) = 0.0; localCoords_(2,7) = 0.5; 
      localCoords_(3,6) = 0.5; localCoords_(3,7) = 0.5;
      //
      localCoords_(0,8) = 0.5; localCoords_(0,9) = 0.0; 
      localCoords_(1,8) = 0.0; localCoords_(1,9) = 0.5; 
      localCoords_(2,8) = 0.5; localCoords_(2,9) = 0.5; 
      localCoords_(3,8) = 0.0; localCoords_(3,9) = 0.0;
    }
  
    // Matrix of local face connectivity:
    // ...opposite point 0,     ...opposite point 2, 
    localFaceConn_(0,0) = 1; localFaceConn_(2,0) = 0;
    localFaceConn_(0,1) = 2; localFaceConn_(2,1) = 1;
    localFaceConn_(0,2) = 3; localFaceConn_(2,2) = 3;
   
    // ...opposite point 1,     ...opposite point 3.
    localFaceConn_(1,0) = 2; localFaceConn_(3,0) = 0; 
    localFaceConn_(1,1) = 0; localFaceConn_(3,1) = 2; 
    localFaceConn_(1,2) = 3; localFaceConn_(3,2) = 1; 

    feInterpolate_ = new FE_InterpolateSimpLin(this);
  }

  FE_ElementTet::~FE_ElementTet()
  {
    // Handled by base class
  }

  void FE_ElementTet::set_quadrature(FeIntQuadrature type) 
  { 
    feInterpolate_->set_quadrature(TETRA,type);
  }

  bool FE_ElementTet::local_coordinates(const DENS_MAT &eltCoords,
                                        const DENS_VEC &x,
                                        DENS_VEC &xi) const
  {
    DENS_MAT T(nSD_, numNodes_-1);
    DENS_VEC r(nSD_);
    for (int iSD=0; iSD<nSD_; ++iSD) {
      for (int inode=1; inode<numNodes_; ++inode) {
        T(iSD, inode-1) = eltCoords(iSD, inode) -
                          eltCoords(iSD, 0);
      }
      r(iSD) = x(iSD) - eltCoords(iSD, 0);
    }
    MultMv(inv(T), r, xi, false, 1.0, 0.0);
    return true;
  }

  bool FE_ElementTet::contains_point(const DENS_MAT &eltCoords,
                                     const DENS_VEC &x) const
  {
    if (! range_check(eltCoords,x) ) return false;
    DENS_VEC xi(nSD_);
    bool converged = local_coordinates(eltCoords, x, xi);
    if (! converged) return false;
    double sum = 0.0;
    bool inside = true;
    for (int iSD = 0; iSD < nSD_; ++iSD) {
      if (dbl_geq(xi(iSD),1.0) || dbl_geq(0.0,xi(iSD))) {
        inside = false;
        break;
      }
      sum += xi(iSD);
    }
    if (dbl_geq(sum,1.0)) inside = false;
    return inside;
  }


}; // namespace ATC

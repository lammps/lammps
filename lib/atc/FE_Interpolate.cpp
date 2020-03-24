// ATC header files
#include "ATC_Error.h"
#include "FE_Element.h"
#include "FE_Interpolate.h"
#include "FE_Quadrature.h"

// Other headers
#include <cmath>

using std::map;
using std::vector;

namespace ATC {

  FE_Interpolate::FE_Interpolate(FE_Element *feElement)
    : feElement_(feElement),
      nSD_(feElement->num_dims())
  {
    // Nothing to do here
  }

  FE_Interpolate::~FE_Interpolate()
  {
    if (!feQuadList_.empty()) {
      map<FeIntQuadrature,FE_Quadrature *>::iterator qit;
      for (qit  = feQuadList_.begin();
           qit != feQuadList_.end(); ++qit) {
        delete (qit->second);
      }
    }
  }

  void FE_Interpolate::set_quadrature(FeEltGeometry geo,
                                      FeIntQuadrature quad) 
  {
    if (feQuadList_.count(quad) == 0) {
      feQuad_ = new FE_Quadrature(geo,quad);
      feQuadList_[quad] = feQuad_;
    } else {
      feQuad_ = feQuadList_[quad];
    }
    precalculate_shape_functions();
  }

  void FE_Interpolate::precalculate_shape_functions()
  {
    int numEltNodes = feElement_->num_elt_nodes();
    int numFaces = feElement_->num_faces();
    int numFaceNodes = feElement_->num_face_nodes();
    
    int numIPs = feQuad_->numIPs;
    DENS_MAT &ipCoords = feQuad_->ipCoords;
    int numFaceIPs = feQuad_->numFaceIPs;
    vector<DENS_MAT> &ipFaceCoords = feQuad_->ipFaceCoords;
    DENS_MAT &ipFace2DCoords = feQuad_->ipFace2DCoords;

    // Compute elemental shape functions at ips
    N_.reset(numIPs,numEltNodes);
    dNdr_.assign(numIPs,DENS_MAT(nSD_,numEltNodes));
    for (int ip = 0; ip < numIPs; ip++) {
      CLON_VEC thisIP = column(ipCoords,ip);
      CLON_VEC thisN = row(N_,ip);
      DENS_MAT &thisdNdr = dNdr_[ip];
      compute_N(thisIP,thisN);
      compute_N_dNdr(thisIP,thisN,thisdNdr);
    }

    // Compute face shape functions at ip's
    NFace_.assign(numFaces,DENS_MAT(numFaceIPs,numEltNodes));
    dNdrFace_.resize(numFaces);
    for (int f = 0; f < numFaces; f++) {
      dNdrFace_[f].assign(numIPs,DENS_MAT(nSD_,numEltNodes));
    }
    for (int f = 0; f < numFaces; f++) {
      for (int ip = 0; ip < numFaceIPs; ip++) {
        CLON_VEC thisIP = column(ipFaceCoords[f],ip);
        CLON_VEC thisN = row(NFace_[f],ip);
        DENS_MAT &thisdNdr = dNdrFace_[f][ip];
        compute_N_dNdr(thisIP,thisN,thisdNdr);
      }
    }
  
    // Compute 2D face shape function derivatives
    dNdrFace2D_.assign(numFaceIPs,DENS_MAT(nSD_-1,numFaceNodes));
    for (int ip = 0; ip < numFaceIPs; ip++) {
      CLON_VEC thisIP = column(ipFace2DCoords,ip);
      DENS_MAT &thisdNdr = dNdrFace2D_[ip];
      compute_dNdr(thisIP,
                   numFaceNodes,nSD_-1,feElement_->face_area(),
                   thisdNdr);
    }
  }

  //-----------------------------------------------------------------
  // shape function value at a particular point given local coordinates
  //-----------------------------------------------------------------
  void FE_Interpolate::shape_function(const VECTOR &xi,
                                      DENS_VEC &N)
  {
    int numEltNodes = feElement_->num_elt_nodes();
    N.resize(numEltNodes);
    compute_N(xi,N);
  }

  void FE_Interpolate::shape_function(const DENS_MAT &eltCoords,
                                      const VECTOR &xi,
                                      DENS_VEC &N,
                                      DENS_MAT &dNdx)
  {
    int numEltNodes = feElement_->num_elt_nodes();
    N.resize(numEltNodes);
    DENS_MAT dNdr(nSD_,numEltNodes,false);
    compute_N_dNdr(xi,N,dNdr);
    
    DENS_MAT eltCoordsT = transpose(eltCoords);
    DENS_MAT dxdr, drdx;
    
    dxdr = dNdr*eltCoordsT;
    drdx = inv(dxdr);
    dNdx = drdx*dNdr;
  } 

  void FE_Interpolate::shape_function_derivatives(const DENS_MAT &eltCoords,
                                      const VECTOR &xi,
                                      DENS_MAT &dNdx)
  {
    int numEltNodes = feElement_->num_elt_nodes();
    DENS_MAT dNdr(nSD_,numEltNodes,false);
    
    DENS_VEC N(numEltNodes); 
    compute_N_dNdr(xi,N,dNdr);
    DENS_MAT eltCoordsT = transpose(eltCoords);
    DENS_MAT dxdr, drdx;
    dxdr = dNdr*eltCoordsT; // tangents or Jacobian matrix
    
    drdx = inv(dxdr);
    dNdx = drdx*dNdr; // dN/dx = dN/dxi (dx/dxi)^-1
  } 

  void FE_Interpolate::tangents(const DENS_MAT &eltCoords,
                                const VECTOR &xi,
                                DENS_MAT &dxdr) const
  {
    int numEltNodes = feElement_->num_elt_nodes();
    DENS_MAT dNdr(nSD_,numEltNodes,false);
    
    DENS_VEC N(numEltNodes); 
    compute_N_dNdr(xi,N,dNdr);
//dNdr.print("dNdr");
    DENS_MAT eltCoordsT = transpose(eltCoords);
//eltCoordsT.print("elt coords");
    dxdr = dNdr*eltCoordsT;
//dxdr.print("dxdr");
  } 

  void FE_Interpolate::tangents(const DENS_MAT &eltCoords,
                                const VECTOR &xi,
                                vector<DENS_VEC> & dxdxi,
                                const bool normalize) const
  {
     DENS_MAT dxdr;
     tangents(eltCoords,xi,dxdr);
//dxdr.print("dxdr-post");
     dxdxi.resize(nSD_);
     //for (int i = 0; i < nSD_; ++i) dxdxi[i] = CLON_VEC(dxdr,CLONE_COL,i); 
     for (int i = 0; i < nSD_; ++i) {
       dxdxi[i].resize(nSD_);
       for (int j = 0; j < nSD_; ++j) {
         dxdxi[i](j) = dxdr(i,j);
       }
     }
//dxdxi[0].print("t1");
//dxdxi[1].print("t2");
//dxdxi[2].print("t3");
     if (normalize) {
       for (int j = 0; j < nSD_; ++j) {
         double norm = 0;
         VECTOR & t = dxdxi[j];
         for (int i = 0; i < nSD_; ++i) norm += t(i)*t(i);
         norm = 1./sqrt(norm);
         for (int i = 0; i < nSD_; ++i) t(i) *= norm;
       }
     }
  }

  // -------------------------------------------------------------
  //   shape_function values at nodes
  // -------------------------------------------------------------
  void FE_Interpolate::shape_function(const DENS_MAT &eltCoords,
                                      DENS_MAT &N,
                                      vector<DENS_MAT> &dN,
                                      DIAG_MAT &weights)
  {
    int numEltNodes = feElement_->num_elt_nodes();
    
    // Transpose eltCoords
    DENS_MAT eltCoordsT(transpose(eltCoords));

    // Shape functions are simply the canonical element values
    N = N_;

    // Set sizes of matrices and vectors
    if ((int)dN.size() != nSD_) dN.resize(nSD_);
    for (int isd = 0; isd < nSD_; isd++) 
      dN[isd].resize(feQuad_->numIPs,numEltNodes);
    weights.resize(feQuad_->numIPs,feQuad_->numIPs);

    // Create some temporary matrices:
    // Jacobian matrix: [dx/dr dy/ds dz/dt | dx/ds ... ]
    DENS_MAT dxdr, drdx, dNdx;

    // Loop over integration points
    for (int ip = 0; ip < feQuad_->numIPs; ip++) {
      // Compute dx/dxi matrix
      dxdr = dNdr_[ip]*eltCoordsT;
      drdx = inv(dxdr);

      // Compute dNdx and fill dN matrix
      dNdx = drdx * dNdr_[ip];
      for (int isd = 0; isd < nSD_; isd++)
        for (int inode = 0; inode < numEltNodes; inode++) 
          dN[isd](ip,inode) = dNdx(isd,inode);

      // Compute jacobian determinant of dxdr at this ip
      double J = dxdr(0,0) * (dxdr(1,1)*dxdr(2,2) - dxdr(2,1)*dxdr(1,2))
               - dxdr(0,1) * (dxdr(1,0)*dxdr(2,2) - dxdr(1,2)*dxdr(2,0))
               + dxdr(0,2) * (dxdr(1,0)*dxdr(2,1) - dxdr(1,1)*dxdr(2,0));

      // Compute ip weight
      weights(ip,ip) = feQuad_->ipWeights(ip)*J; 
    }
  
  }

  //-----------------------------------------------------------------
  // shape functions on a given face
  //-----------------------------------------------------------------
  void FE_Interpolate::face_shape_function(const DENS_MAT &eltCoords,
                                           const DENS_MAT &faceCoords,
                                           const int faceID,
                                           DENS_MAT &N,
                                           DENS_MAT &n,
                                           DIAG_MAT &weights)
  {
    int numFaceIPs = feQuad_->numFaceIPs;

    // Transpose eltCoords
    DENS_MAT eltCoordsT = transpose(eltCoords);

    // Shape functions are simply the canonical element values
    N = NFace_[faceID];

    // Create some temporary matrices:
    // Jacobian matrix: [dx/dr dy/ds dz/dt | dx/ds ... ]
    DENS_MAT dxdr, drdx, dNdx;

    // Loop over integration points
    DENS_VEC normal(nSD_);
    n.resize(nSD_,numFaceIPs);
    weights.resize(numFaceIPs,numFaceIPs);
    for (int ip = 0; ip < numFaceIPs; ip++) {
      // Compute 2d jacobian determinant of dxdr at this ip
      double J = face_normal(faceCoords,ip,normal);

      // Copy normal at integration point
      for (int isd = 0; isd < nSD_; isd++) {
        n(isd,ip) = normal(isd);
      }

      // Compute ip weight
      weights(ip,ip) = feQuad_->ipFaceWeights(ip)*J; 
    }
  
  }

  void FE_Interpolate::face_shape_function(const DENS_MAT &eltCoords,
                                           const DENS_MAT &faceCoords,
                                           const int faceID,
                                           DENS_MAT &N,
                                           vector<DENS_MAT> &dN,
                                           vector<DENS_MAT> &Nn,
                                           DIAG_MAT &weights)
  {
    int numEltNodes = feElement_->num_elt_nodes();
    int numFaceIPs = feQuad_->numFaceIPs;

    // Transpose eltCoords
    DENS_MAT eltCoordsT = transpose(eltCoords);

    // Shape functions are simply the canonical element values
    N = NFace_[faceID];

    // Set sizes of matrices and vectors
    if ((int)dN.size() != nSD_) dN.resize(nSD_);
    if ((int)Nn.size() != nSD_) Nn.resize(nSD_);
    for (int isd = 0; isd < nSD_; isd++) {
      dN[isd].resize(numFaceIPs,numEltNodes);
      Nn[isd].resize(numFaceIPs,numEltNodes);
    }
    weights.resize(numFaceIPs,numFaceIPs);

    // Create some temporary matrices:
    // Jacobian matrix: [dx/dr dy/ds dz/dt | dx/ds ... ]
    DENS_MAT dxdr, drdx, dNdx;
    DENS_VEC normal(nSD_);

    // Loop over integration points
    for (int ip = 0; ip < numFaceIPs; ip++) {
      // Compute dx/dxi matrix
      dxdr = dNdrFace_[faceID][ip] * eltCoordsT;
      drdx = inv(dxdr);

      // Compute 2d jacobian determinant of dxdr at this ip
      double J = face_normal(faceCoords,ip,normal);

      // Compute dNdx and fill dN matrix
      dNdx = drdx * dNdrFace_[faceID][ip];
      for (int isd = 0; isd < nSD_; isd++) {
        for (int inode = 0; inode < numEltNodes; inode++) {
          dN[isd](ip,inode) = dNdx(isd,inode);
          Nn[isd](ip,inode) = N(ip,inode)*normal(isd);
        }
      }

      // Compute ip weight
      weights(ip,ip) = feQuad_->ipFaceWeights(ip)*J; 
    }
  
  }

  // -------------------------------------------------------------
  //   face normal
  // -------------------------------------------------------------
  double FE_Interpolate::face_normal(const DENS_MAT &faceCoords,
                                     int ip,
                                     DENS_VEC &normal) 
  {
    // Compute dx/dr for deformed element
    DENS_MAT faceCoordsT = transpose(faceCoords);
    DENS_MAT dxdr = dNdrFace2D_[ip]*faceCoordsT;

    // Normal vector from cross product, hardcoded for 3D, sad
    normal(0) = dxdr(0,1)*dxdr(1,2) - dxdr(0,2)*dxdr(1,1);
    normal(1) = dxdr(0,2)*dxdr(1,0) - dxdr(0,0)*dxdr(1,2);
    normal(2) = dxdr(0,0)*dxdr(1,1) - dxdr(0,1)*dxdr(1,0);

    // Jacobian (3D)
    double J = sqrt(normal(0)*normal(0) + 
                    normal(1)*normal(1) + 
                    normal(2)*normal(2));
    double invJ = 1.0/J;
    normal(0) *= invJ;
    normal(1) *= invJ;
    normal(2) *= invJ;
    return J;
  }
  
  int FE_Interpolate::num_ips() const
  { 
    return feQuad_->numIPs; 
  }

  int FE_Interpolate::num_face_ips() const
  { 
    return feQuad_->numFaceIPs; 
  }

  
  /*********************************************************
   * Class FE_InterpolateCartLagrange
   *
   * For computing Lagrange shape functions using Cartesian
   * coordinate systems (all quads/hexes fall under this
   * category, and any elements derived by degenerating 
   * them). Not to be used for serendipity elements, which
   * should be implemented for SPEED.
   *
   *********************************************************/
  FE_InterpolateCartLagrange::FE_InterpolateCartLagrange(
                                          FE_Element *feElement)
    : FE_Interpolate(feElement)
  {
    set_quadrature(HEXA,GAUSS2);
  }

  FE_InterpolateCartLagrange::~FE_InterpolateCartLagrange()
  {
    // Handled by base class
  }

  void FE_InterpolateCartLagrange::compute_N(const VECTOR &point,
                                             VECTOR &N)
  {
    // *** see comments for compute_N_dNdr ***
    const DENS_VEC &localCoords1d = feElement_->local_coords_1d();
    int numEltNodes = feElement_->num_elt_nodes();
    int numEltNodes1d = feElement_->num_elt_nodes_1d();
    
    DENS_MAT lagrangeTerms(nSD_,numEltNodes1d);
    DENS_MAT lagrangeDenom(nSD_,numEltNodes1d);
    lagrangeTerms = 1.0;
    lagrangeDenom = 1.0;
    for (int iSD = 0; iSD < nSD_; ++iSD) {
      for (int inode = 0; inode < numEltNodes1d; ++inode) {
        for (int icont = 0; icont < numEltNodes1d; ++icont) {
          if (inode != icont) {
            lagrangeDenom(iSD,inode) *= (localCoords1d(inode) - 
                                         localCoords1d(icont));
            lagrangeTerms(iSD,inode) *= (point(iSD)-localCoords1d(icont));
          }
        }
      }
    }
    for (int iSD=0; iSD<nSD_; ++iSD) {
      for (int inode=0; inode<numEltNodes1d; ++inode) {
        lagrangeTerms(iSD,inode) /= lagrangeDenom(iSD,inode);
      }
    }

    N = 1.0;
    vector<int> mapping(nSD_);
    for (int inode=0; inode<numEltNodes; ++inode) {
      feElement_->mapping(inode,mapping);
      for (int iSD=0; iSD<nSD_; ++iSD) {
        N(inode) *= lagrangeTerms(iSD,mapping[iSD]);
      }
    }
  }

  // Sort of a test-ride for a generic version that can be used for
  // faces too. The only thing that's not "generic" is localCoords,
  // which very magically works in both cases.
  void FE_InterpolateCartLagrange::compute_dNdr(const VECTOR &point,
                                                const int numNodes,
                                                const int nD,
                                                const double,
                                                DENS_MAT &dNdr)
  {
    // *** see comments for compute_N_dNdr ***
    const DENS_VEC &localCoords1d = feElement_->local_coords_1d();
    int numEltNodes1d = feElement_->num_elt_nodes_1d();
    
    DENS_MAT lagrangeTerms(nD,numEltNodes1d);
    DENS_MAT lagrangeDenom(nD,numEltNodes1d);
    DENS_MAT lagrangeDeriv(nD,numEltNodes1d);
    lagrangeDenom = 1.0;
    lagrangeTerms = 1.0;
    lagrangeDeriv = 0.0;
    DENS_VEC productRuleVec(numEltNodes1d);
    productRuleVec = 1.0;
    for (int iSD = 0; iSD < nD; ++iSD) {
      for (int inode = 0; inode < numEltNodes1d; ++inode) {
        for (int icont = 0; icont < numEltNodes1d; ++icont) {
          if (inode != icont) {
            lagrangeTerms(iSD,inode) *= (point(iSD)-localCoords1d(icont));
            lagrangeDenom(iSD,inode) *= (localCoords1d(inode) - 
                                         localCoords1d(icont));
            for (int dcont=0; dcont<numEltNodes1d; ++dcont) {
              if (inode == dcont) {
                productRuleVec(dcont) = 0.0;
              } else if (icont == dcont) {
              } else {
                productRuleVec(dcont) *= (point(iSD)-localCoords1d(icont));
              }
            }
          }
        }
        for (int dcont=0; dcont<numEltNodes1d; ++dcont) {
          lagrangeDeriv(iSD,inode) += productRuleVec(dcont);
        }
        productRuleVec = 1.0;
      }
    }
    for (int iSD=0; iSD<nD; ++iSD) {
      for (int inode=0; inode<numEltNodes1d; ++inode) {
        lagrangeTerms(iSD,inode) /= lagrangeDenom(iSD,inode);
        lagrangeDeriv(iSD,inode) /= lagrangeDenom(iSD,inode);
      }
    }
    
    dNdr = 1.0;
    vector<int> mapping(nD);
    for (int inode=0; inode<numNodes; ++inode) {
      feElement_->mapping(inode,mapping);
      for (int iSD=0; iSD<nD; ++iSD) {
        for (int dSD=0; dSD<nD; ++dSD) {
          if (iSD == dSD) {
            dNdr(dSD,inode) *= lagrangeDeriv(iSD,mapping[iSD]);
          } else {
            dNdr(dSD,inode) *= lagrangeTerms(iSD,mapping[iSD]);
          }
        }
      }
    }
  }

  void FE_InterpolateCartLagrange::compute_N_dNdr(const VECTOR &point,
                                                  VECTOR &N,
                                                  DENS_MAT &dNdr) const
  {
    // Required data from element class
    const DENS_VEC &localCoords1d = feElement_->local_coords_1d();
    int numEltNodes = feElement_->num_elt_nodes();
    int numEltNodes1d = feElement_->num_elt_nodes_1d();
    
    // lagrangeTerms stores the numerator for the various Lagrange polynomials
    // in one dimension, that will be used to produce the three dimensional
    // shape functions
    DENS_MAT lagrangeTerms(nSD_,numEltNodes1d);
    // lagrangeDenom stores the denominator. Stored separately to reduce
    // redundancy, because it will be used for the shape functions and derivs
    DENS_MAT lagrangeDenom(nSD_,numEltNodes1d);
    // lagrangeDeriv stores the numerator for the derivative of the Lagrange
    // polynomials
    DENS_MAT lagrangeDeriv(nSD_,numEltNodes1d);
    // Terms/Denom are products, Deriv will be a sum, so initialize as such:
    lagrangeTerms = 1.0;
    lagrangeDenom = 1.0;
    lagrangeDeriv = 0.0;
    // the derivative requires use of the product rule; to store the prodcuts
    // which make up the terms produced by the product rule, we'll use this
    // vector
    DENS_VEC productRuleVec(numEltNodes1d);
    productRuleVec = 1.0;
    for (int iSD = 0; iSD < nSD_; ++iSD) {
      for (int inode = 0; inode < numEltNodes1d; ++inode) {
        for (int icont = 0; icont < numEltNodes1d; ++icont) {
          if (inode != icont) {
            // each dimension and each 1d node per dimension has a 
            // contribution from all nodes besides the current node
            lagrangeTerms(iSD,inode) *= (point(iSD)-localCoords1d(icont));
            lagrangeDenom(iSD,inode) *= (localCoords1d(inode) - 
                                         localCoords1d(icont));
            // complciated; each sum produced by the product rule has one
            // "derivative", and the rest are just identical to the terms 
            // above
            for (int dcont=0; dcont<numEltNodes1d; ++dcont) {
              if (inode == dcont) {
                // skip this term, derivative is 0
                productRuleVec(dcont) = 0.0;
              } else if (icont == dcont) {
                // no numerator contribution, derivative is 1
              } else {
                // part of the "constant"
                productRuleVec(dcont) *= (point(iSD)-localCoords1d(icont));
              }
            }
          }
        }
        // sum the terms produced by the product rule and store in Deriv
        for (int dcont=0; dcont<numEltNodes1d; ++dcont) {
          lagrangeDeriv(iSD,inode) += productRuleVec(dcont);
        }
        productRuleVec = 1.0;
      }
    }
    // divide by denom
    for (int iSD=0; iSD<nSD_; ++iSD) {
      for (int inode=0; inode<numEltNodes1d; ++inode) {
        lagrangeTerms(iSD,inode) /= lagrangeDenom(iSD,inode);
        lagrangeDeriv(iSD,inode) /= lagrangeDenom(iSD,inode);
      }
    }
    
    N = 1.0;
    dNdr = 1.0;
    // mapping returns the 1d nodes in each dimension that should be multiplied
    // to achieve the shape functions in 3d
    vector<int> mapping(nSD_);
    for (int inode=0; inode<numEltNodes; ++inode) {
      feElement_->mapping(inode,mapping);
      for (int iSD=0; iSD<nSD_; ++iSD) {
        N(inode) *= lagrangeTerms(iSD,mapping[iSD]);
        for (int dSD=0; dSD<nSD_; ++dSD) {
          // only use Deriv for the dimension in which we're taking the
          // derivative, because the rest is essentially a "constant"
          if (iSD == dSD) {
            dNdr(dSD,inode) *= lagrangeDeriv(iSD,mapping[iSD]);
          } else {
            dNdr(dSD,inode) *= lagrangeTerms(iSD,mapping[iSD]);
          }
        }
      }
    }
  }

  /*********************************************************
   * Class FE_InterpolateCartLin
   *
   * For computing linear shape functions using Cartesian
   * coordinate systems (all quads/hexes fall under this
   * category, and any elements derived by degenerating 
   * them).
   *
   *********************************************************/
  FE_InterpolateCartLin::FE_InterpolateCartLin(
                                          FE_Element *feElement)
    : FE_Interpolate(feElement)
  {
    set_quadrature(HEXA,GAUSS2);
  }

  FE_InterpolateCartLin::~FE_InterpolateCartLin()
  {
    // Handled by base class
  }
  
  void FE_InterpolateCartLin::compute_N(const VECTOR &point,
                                        VECTOR &N)
  {
    // *** see comments for compute_N_dNdr ***
    const DENS_MAT &localCoords = feElement_->local_coords();
    double invVol = 1.0/(feElement_->vol());
    int numEltNodes = feElement_->num_elt_nodes();
    
    for (int inode = 0; inode < numEltNodes; ++inode) {
      N(inode) = invVol;
      for (int isd = 0; isd < nSD_; ++isd) {
        N(inode) *= (1.0 + point(isd)*localCoords(isd,inode));
      }
    }
  }

  // Sort of a test-ride for a generic version that can be used for
  // faces too. The only thing that's not "generic" is localCoords,
  // which very magically works in both cases.
  void FE_InterpolateCartLin::compute_dNdr(const VECTOR &point,
                                           const int numNodes,
                                           const int nD,
                                           const double vol,
                                           DENS_MAT &dNdr)
  {
    // *** see comments for compute_N_dNdr ***
    const DENS_MAT &localCoords = feElement_->local_coords();
    double invVol = 1.0/vol;
    
    for (int inode = 0; inode < numNodes; ++inode) {
      for (int idr = 0; idr < nD; ++idr) {
        dNdr(idr,inode) = invVol;
      }
      for (int id = 0; id < nD; ++id) {
        for (int idr = 0; idr < nD; ++idr) {
          if (id == idr) dNdr(idr,inode) *= localCoords(id,inode);
          else dNdr(idr,inode) *= 1.0 + 
                                  point(id)*localCoords(id,inode);
        }
      }
    }
  }

  void FE_InterpolateCartLin::compute_N_dNdr(const VECTOR &point,
                                             VECTOR &N,
                                             DENS_MAT &dNdr) const
  {
    // Required data from element class
    const DENS_MAT &localCoords = feElement_->local_coords();
    double invVol = 1.0/(feElement_->vol());
    int numEltNodes = feElement_->num_elt_nodes();
    
    // Fill in for each node
    for (int inode = 0; inode < numEltNodes; ++inode) {
      // Intiialize shape function and derivatives
      N(inode) = invVol;
      for (int idr = 0; idr < nSD_; ++idr) {
        dNdr(idr,inode) = invVol;
      }
      for (int isd = 0; isd < nSD_; ++isd) {
        // One term for each dimension
        N(inode) *= (1.0 + point(isd)*localCoords(isd,inode));
        // One term for each dimension, only deriv in deriv's dimension
        for (int idr = 0; idr < nSD_; ++idr) {
          if (isd == idr) dNdr(idr,inode) *= localCoords(isd,inode);
          else dNdr(idr,inode) *= 1.0 + 
                                  point(isd)*localCoords(isd,inode);
        }
      }
    }
  }

  
  /*********************************************************
   * Class FE_InterpolateCartSerendipity
   *
   * For computing shape functions for quadratic serendipity
   * elements, implemented for SPEED.
   * 
   *********************************************************/
  FE_InterpolateCartSerendipity::FE_InterpolateCartSerendipity(
                                   FE_Element *feElement)
    : FE_Interpolate(feElement)
  {
    set_quadrature(HEXA,GAUSS2);
  }

  FE_InterpolateCartSerendipity::~FE_InterpolateCartSerendipity()
  {
    // Handled by base class
  }
  
  void FE_InterpolateCartSerendipity::compute_N(const VECTOR &point,
                                                VECTOR &N)
  {
    // *** see comments for compute_N_dNdr ***
    const DENS_MAT &localCoords = feElement_->local_coords();
    double invVol = 1.0/(feElement_->vol());
    int numEltNodes = feElement_->num_elt_nodes();
    
    for (int inode = 0; inode < numEltNodes; ++inode) {
      N(inode) = invVol;
      for (int isd = 0; isd < nSD_; ++isd) {
        if (((inode == 8  || inode == 10 || inode == 16 || inode == 18) && 
             (isd == 0))                                                   || 
            ((inode == 9  || inode == 11 || inode == 17 || inode == 19) &&
             (isd == 1))                                                   || 
            ((inode == 12 || inode == 13 || inode == 14 || inode == 15) &&
             (isd == 2))) {
          N(inode) *= (1.0 - pow(point(isd),2))*2;
        } else {
          N(inode) *= (1.0 + point(isd)*localCoords(isd,inode));
        }
      }
      if (inode < 8) {
        N(inode) *= (point(0)*localCoords(0,inode) +
                     point(1)*localCoords(1,inode) +
                     point(2)*localCoords(2,inode) - 2);
      }
    }
  }

  // Sort of a test-ride for a generic version that can be used for
  // faces too. The only thing that's not "generic" is localCoords,
  // which very magically works in both cases.
  void FE_InterpolateCartSerendipity::compute_dNdr(const VECTOR &point,
                                                   const int numNodes,
                                                   const int nD,
                                                   const double vol,
                                                   DENS_MAT &dNdr)
  {
    // *** see comments for compute_N_dNdr ***
    const DENS_MAT &localCoords = feElement_->local_coords();
    double invVol = 1.0/vol;
    bool serendipityNode = false;
    double productRule1 = 0.0;
    double productRule2 = 0.0;
    
    if (nD != 3 && nD != 2) {
      ATC_Error("Serendipity dNdr calculations are too hard-wired to do "
                "what you want them to. Only 2D and 3D currently work.");
    }

    for (int inode = 0; inode < numNodes; ++inode) {
      for (int idr = 0; idr < nD; ++idr) {
        dNdr(idr,inode) = invVol;
      }
      for (int id = 0; id < nD; ++id) {
        for (int idr = 0; idr < nD; ++idr) {
          // identify nodes/dims differently if 3d or 2d case
          if (nD == 3) {
            serendipityNode =
             (((inode == 8  || inode == 10 || inode == 16 || inode == 18) && 
               (id == 0))                                                    || 
              ((inode == 9  || inode == 11 || inode == 17 || inode == 19) &&
               (id == 1))                                                    || 
              ((inode == 12 || inode == 13 || inode == 14 || inode == 15) &&
               (id == 2)));
          } else if (nD == 2) {
            serendipityNode =
             (((inode == 4 || inode == 6) && (id == 0))  || 
              ((inode == 5 || inode == 7) && (id == 1)));
          }
          if (serendipityNode) {
            if (id == idr) {
              dNdr(idr,inode) *= point(id)*(-4);
            } else {
              dNdr(idr,inode) *= (1.0 - pow(point(id),2))*2;
            }
          } else {
            if (id == idr) {
              dNdr(idr,inode) *= localCoords(id,inode);
            } else {
              dNdr(idr,inode) *= (1.0 + point(id)*localCoords(id,inode));
            }
          }
        }
      }
      for (int idr = 0; idr < nD; ++idr) {
        if (inode < 8) {
          // final corner contribution slightly different for 3d and 2d cases
          if (nD == 3) {
            productRule2 = (point(0)*localCoords(0,inode) +
                            point(1)*localCoords(1,inode) +
                            point(2)*localCoords(2,inode) - 2);
          } else if (nD == 2) {
            productRule2 = (point(0)*localCoords(0,inode) +
                            point(1)*localCoords(1,inode) - 1);
          }
          productRule1 = dNdr(idr,inode) * 
                         (1 + point(idr)*localCoords(idr,inode));
          productRule2 *= dNdr(idr,inode);
          dNdr(idr,inode) = productRule1 + productRule2;
        }
      }
    }
  }

  void FE_InterpolateCartSerendipity::compute_N_dNdr(const VECTOR &point,
                                                     VECTOR &N,
                                                     DENS_MAT &dNdr) const
  {
    // Required data from element class
    const DENS_MAT &localCoords = feElement_->local_coords();
    double invVol = 1.0/(feElement_->vol());
    int numEltNodes = feElement_->num_elt_nodes();

    // Will store terms for product rule derivative for dNdr
    double productRule1;
    double productRule2;

    // Fill in for each node
    for (int inode = 0; inode < numEltNodes; ++inode) {
      // Initialize shape functions and derivatives
      N(inode) = invVol;
      for (int idr = 0; idr < nSD_; ++idr) {
        dNdr(idr,inode) = invVol;
      }
      // Add components from each dimension
      for (int isd = 0; isd < nSD_; ++isd) {
        for (int idr = 0; idr < nSD_; ++idr) {
          // Check to see if the node is NOT a corner node, and if its
          // "0-coordinate" is in the same dimension as the one we're currently
          // iterating over. If that's the case, we want to contribute to its
          // shape functions and derivatives in a modified way:
          if (((inode == 8  || inode == 10 || inode == 16 || inode == 18) && 
               (isd == 0))                                                   || 
              ((inode == 9  || inode == 11 || inode == 17 || inode == 19) &&
               (isd == 1))                                                   || 
              ((inode == 12 || inode == 13 || inode == 14 || inode == 15) &&
               (isd == 2))) {
            // If the 1d shape function dimension matches the derivative
            // dimension...
            if (isd == idr) {
              // contribute to N; sloppy, but this is the easiest way to get
              // N to work right without adding extra, arguably unnecessary 
              // loops, while also computing the shape functions
              N(inode) *= (1.0 - pow(point(isd),2))*2;
              // contribute to dNdr with the derivative of this shape function
              // contribution
              dNdr(idr,inode) *= point(isd)*(-4);
            } else {
              // otherwise, just use the "constant" contribution to the deriv
              dNdr(idr,inode) *= (1.0 - pow(point(isd),2))*2;
            }
          } else {
            // non-serendipity style contributions
            if (isd == idr) {
              N(inode) *= (1.0 + point(isd)*localCoords(isd,inode));
              dNdr(idr,inode) *= localCoords(isd,inode);
            } else {
              dNdr(idr,inode) *= (1.0 + point(isd)*localCoords(isd,inode));
            }
          }
        }
      }
      // serendipity corner nodes require more extra handling
      if (inode < 8) {
        N(inode) *= (point(0)*localCoords(0,inode) +
                     point(1)*localCoords(1,inode) +
                     point(2)*localCoords(2,inode) - 2);
      }
      for (int idr = 0; idr < nSD_; ++idr) {
        if (inode < 8) {
          productRule1 = dNdr(idr,inode) * 
                         (1 + point(idr)*localCoords(idr,inode));
          productRule2 = dNdr(idr,inode) * (point(0)*localCoords(0,inode) +
                                            point(1)*localCoords(1,inode) +
                                            point(2)*localCoords(2,inode) - 2);
          dNdr(idr,inode) = productRule1 + productRule2;
        }
      }
    }
  }


  /*********************************************************
   * Class FE_InterpolateSimpLin
   * 
   * For computing linear shape functions of simplices,
   * which are rather different from those computed
   * in Cartesian coordinates.
   *
   * Note: degenerating quads/hexes can yield simplices
   *       as well, but this class is for computing these
   *       shape functions _natively_, in their own
   *       triangular/tetrahedral coordinate systems.
   *
   *********************************************************/
  FE_InterpolateSimpLin::FE_InterpolateSimpLin(
                                        FE_Element *feElement) 
    : FE_Interpolate(feElement)
  {
    set_quadrature(TETRA,GAUSS2);
  }

  FE_InterpolateSimpLin::~FE_InterpolateSimpLin()
  {
    // Handled by base class
  }

  void FE_InterpolateSimpLin::compute_N(const VECTOR &point,
                                        VECTOR &N)
  {
    int numEltNodes = feElement_->num_elt_nodes();
    
    // Fill in for each node
    for (int inode = 0; inode < numEltNodes; ++inode) {
      if (inode == 0) {
        // Fill N...the ips are serving as proxies for "dimensions"
        // since we're in tetrahedral coordinates, except that
        //   0th node = 3rd "dimension" (u or O_o)
        //   1st node = 0th "dimension" (x or r)
        //   2nd node = 1st "dimension" (y or s)
        //   3rd node = 3nd "dimension" (z or t)
        // and remember that u = 1 - r - s - t for tet coords
        N(inode) = 1;
        for (int icont = 0; icont < nSD_; ++icont) {
          N(inode) -= point(icont);
        }
      } else {
        N(inode) = point(inode-1);
      }
    }
  }

  void FE_InterpolateSimpLin::compute_dNdr(const VECTOR &,
                                           const int numNodes,
                                           const int nD,
                                           const double,
                                           DENS_MAT &dNdr)
  {
    // Fill in for each node
    for (int inode = 0; inode < numNodes; ++inode) {
      // Fill dNdr_; we want 1 if the dimension of derivative
      // and variable within N correspond. That is, if N == r,
      // we want the 0th dimension to contain (d/dr)r = 1. Of
      // course, (d/di)r = 0 forall i != r, so we need that as
      // well. This is a bit elusively complicated. Also, the 0th
      // integration point contains the term u = 1 - r - s - t.
      // (which map to x, y, and z). Therefore, the derivative in
      // each dimension are -1.
      //
      // The idea is similar for 2 dimensions, which this can 
      // handle as well.
      for (int idr = 0; idr < nD; ++idr) {
        if (inode == 0) {
          dNdr(idr,inode) = -1;
        } else {
          dNdr(idr,inode) = (inode == (idr + 1)) ? 1 : 0;
        }
      }
    }
  }

  void FE_InterpolateSimpLin::compute_N_dNdr(const VECTOR &point,
                                             VECTOR &N,
                                             DENS_MAT &dNdr) const
  {
    int numEltNodes = feElement_->num_elt_nodes();
    
    // Fill in for each node
    for (int inode = 0; inode < numEltNodes; ++inode) {
      // Fill N...
      if (inode == 0) {
        N(inode) = 1;
        for (int icont = 0; icont < nSD_; ++icont) {
          N(inode) -= point(icont);
        }
      } else {
        N(inode) = point(inode-1);
      }
      // Fill dNdr...
      for (int idr = 0; idr < nSD_; ++idr) {
        if (inode == 0) {
          dNdr(idr,inode) = -1;
        } else {
          dNdr(idr,inode) = (inode == (idr + 1)) ? 1 : 0;
        }
      }
    }
  }


} // namespace ATC


// ATC header files
#include "ATC_Error.h"
#include "FE_Element.h"
#include "FE_Mesh.h"

// Other headers
#include "math.h"

namespace ATC {

  // -------------------------------------------------------------
  // -------------------------------------------------------------
  //   class FE_Element
  // -------------------------------------------------------------
  // -------------------------------------------------------------

  FE_Element::FE_Element(FE_Mesh * feMesh, 
                         int nSD,
                         int numEltNodes, 
                         int numIPs,
                         int numFaces,
                         int numFaceNodes,
                         int numFaceIPs)
    : feMesh_(feMesh),
      nSD_(nSD),
      numEltNodes_(numEltNodes),
      numIPs_(numIPs),
      numFaces_(numFaces),
      numFaceNodes_(numFaceNodes),
      numFaceIPs_(numFaceIPs),
      localCoords_(nSD,numEltNodes),
      ipCoords_(nSD,numIPs),
      N_(numIPs,numEltNodes),
      ipWeights_(numIPs),
      localFaceConn_(numFaces,numFaceNodes),
      ipFace2DCoords_(nSD-1,numFaceIPs),
      ipFaceCoords_(numFaces),
      //   localFaceCoords_(nSD-1,numFaceNodes),
      NFace_(numFaces),
      ipFaceWeights_(numFaceIPs)
  {
    for (int f = 0; f < numFaces; f++) {
      ipFaceCoords_[f].reset(nSD,numFaceIPs);
      NFace_[f].reset(numFaceIPs,numEltNodes);
    }
  }

  FE_Element::~FE_Element()
  {
    // Nothing else to do
  }


  // -------------------------------------------------------------
  // -------------------------------------------------------------
  //   class FE_ElementHex
  // -------------------------------------------------------------
  // -------------------------------------------------------------

  FE_ElementHex::FE_ElementHex(FE_Mesh * feMesh)
    : FE_Element(feMesh, 3, 8, 8, 6, 4, 4)
  {

    // Matrix of local nodal coordinates
    localCoords_(0,0) = -1; localCoords_(1,0) = -1; localCoords_(2,0) = -1;
    localCoords_(0,1) = +1; localCoords_(1,1) = -1; localCoords_(2,1) = -1;
    localCoords_(0,2) = +1; localCoords_(1,2) = +1; localCoords_(2,2) = -1;
    localCoords_(0,3) = -1; localCoords_(1,3) = +1; localCoords_(2,3) = -1;
    localCoords_(0,4) = -1; localCoords_(1,4) = -1; localCoords_(2,4) = +1;
    localCoords_(0,5) = +1; localCoords_(1,5) = -1; localCoords_(2,5) = +1;
    localCoords_(0,6) = +1; localCoords_(1,6) = +1; localCoords_(2,6) = +1;
    localCoords_(0,7) = -1; localCoords_(1,7) = +1; localCoords_(2,7) = +1;


    //       3 --- 2
    //      /|    /|
    //     / 0 --/ 1     y
    //    7 --- 6 /      |
    //    |     |/       |
    //    4 --- 5         ---> x
    //                  /
    //                 /
    //                z
  
    // Matrix of local face connectivity
    // -x
    localFaceConn_(0,0) = 0; 
    localFaceConn_(0,1) = 4; 
    localFaceConn_(0,2) = 7; 
    localFaceConn_(0,3) = 3; 
 
    // +x
    localFaceConn_(1,0) = 1; 
    localFaceConn_(1,1) = 2; 
    localFaceConn_(1,2) = 6; 
    localFaceConn_(1,3) = 5; 

    // -y
    localFaceConn_(2,0) = 0; 
    localFaceConn_(2,1) = 1; 
    localFaceConn_(2,2) = 5; 
    localFaceConn_(2,3) = 4; 
 
    // +y
    localFaceConn_(3,0) = 2; 
    localFaceConn_(3,1) = 3; 
    localFaceConn_(3,2) = 7; 
    localFaceConn_(3,3) = 6; 
 
    // -z
    localFaceConn_(4,0) = 0; 
    localFaceConn_(4,1) = 3; 
    localFaceConn_(4,2) = 2; 
    localFaceConn_(4,3) = 1; 

    // +z
    localFaceConn_(5,0) = 4; 
    localFaceConn_(5,1) = 5; 
    localFaceConn_(5,2) = 6; 
    localFaceConn_(5,3) = 7; 

    /*
      localFaceCoords_(0,0) = -1; localFaceCoords_(1,0) = -1;
      localFaceCoords_(0,1) = +1; localFaceCoords_(1,1) = -1; 
      localFaceCoords_(0,2) = +1; localFaceCoords_(1,2) = +1; 
      localFaceCoords_(0,3) = -1; localFaceCoords_(1,3) = +1;
    */

    set_quadrature(GAUSSIAN_QUADRATURE);

  }

  //-----------------------------------------------------------------
  void FE_ElementHex::set_quadrature(int quadrature) 
  {
    double a = 1./sqrt(3.);
    if (quadrature == NODAL_QUADRATURE) {
           a = 1.0;
    }
    // Matrix of integration point locations  (Gaussian) & follows local conn
    ipCoords_(0,0) = -a; ipCoords_(1,0) = -a; ipCoords_(2,0) = -a;
    ipCoords_(0,1) = +a; ipCoords_(1,1) = -a; ipCoords_(2,1) = -a;
    ipCoords_(0,2) = +a; ipCoords_(1,2) = +a; ipCoords_(2,2) = -a;
    ipCoords_(0,3) = -a; ipCoords_(1,3) = +a; ipCoords_(2,3) = -a;
    ipCoords_(0,4) = -a; ipCoords_(1,4) = -a; ipCoords_(2,4) = +a;
    ipCoords_(0,5) = +a; ipCoords_(1,5) = -a; ipCoords_(2,5) = +a;
    ipCoords_(0,6) = +a; ipCoords_(1,6) = +a; ipCoords_(2,6) = +a;
    ipCoords_(0,7) = -a; ipCoords_(1,7) = +a; ipCoords_(2,7) = +a;

    // Compute shape functions at ip's
    for (int ip = 0; ip < numIPs_; ip++) {
      double r = ipCoords_(0, ip);
      double s = ipCoords_(1, ip);
      double t = ipCoords_(2, ip);
      for (int iNode = 0; iNode < numEltNodes_; iNode++) {
        double rI = localCoords_(0, iNode);
        double sI = localCoords_(1, iNode);
        double tI = localCoords_(2, iNode);
        N_(ip, iNode) = 0.125 * (1 + r*rI) * (1 + s*sI) * (1 + t*tI);
      }
    }

    // Integration point weights
    ipWeights_ = 1.0;

    // integration points by face
    ipFace2DCoords_(0,0)=-a; ipFace2DCoords_(1,0)=-a;
    ipFace2DCoords_(0,1)= a; ipFace2DCoords_(1,1)=-a; 
    ipFace2DCoords_(0,2)= a; ipFace2DCoords_(1,2)= a; 
    ipFace2DCoords_(0,3)=-a; ipFace2DCoords_(1,3)= a; 

    ipFaceCoords_[0](0,0)=-1;ipFaceCoords_[0](1,0)=-a;ipFaceCoords_[0](2,0)=-a;
    ipFaceCoords_[0](0,1)=-1;ipFaceCoords_[0](1,1)= a;ipFaceCoords_[0](2,1)=-a; 
    ipFaceCoords_[0](0,2)=-1;ipFaceCoords_[0](1,2)= a;ipFaceCoords_[0](2,2)= a; 
    ipFaceCoords_[0](0,3)=-1;ipFaceCoords_[0](1,3)=-a;ipFaceCoords_[0](2,3)= a; 

    ipFaceCoords_[1](0,0)= 1;ipFaceCoords_[1](1,0)=-a;ipFaceCoords_[1](2,0)=-a;
    ipFaceCoords_[1](0,1)= 1;ipFaceCoords_[1](1,1)= a;ipFaceCoords_[1](2,1)=-a; 
    ipFaceCoords_[1](0,2)= 1;ipFaceCoords_[1](1,2)= a;ipFaceCoords_[1](2,2)= a; 
    ipFaceCoords_[1](0,3)= 1;ipFaceCoords_[1](1,3)=-a;ipFaceCoords_[1](2,3)= a; 

    ipFaceCoords_[2](0,0)=-a;ipFaceCoords_[2](1,0)=-1;ipFaceCoords_[2](2,0)=-a;
    ipFaceCoords_[2](0,1)=-a;ipFaceCoords_[2](1,1)=-1;ipFaceCoords_[2](2,1)= a; 
    ipFaceCoords_[2](0,2)= a;ipFaceCoords_[2](1,2)=-1;ipFaceCoords_[2](2,2)= a; 
    ipFaceCoords_[2](0,3)= a;ipFaceCoords_[2](1,3)=-1;ipFaceCoords_[2](2,3)=-a; 

    ipFaceCoords_[3](0,0)=-a;ipFaceCoords_[3](1,0)= 1;ipFaceCoords_[3](2,0)=-a;
    ipFaceCoords_[3](0,1)=-a;ipFaceCoords_[3](1,1)= 1;ipFaceCoords_[3](2,1)= a; 
    ipFaceCoords_[3](0,2)= a;ipFaceCoords_[3](1,2)= 1;ipFaceCoords_[3](2,2)= a; 
    ipFaceCoords_[3](0,3)= a;ipFaceCoords_[3](1,3)= 1;ipFaceCoords_[2](2,3)=-a; 

    ipFaceCoords_[4](0,0)=-a;ipFaceCoords_[4](1,0)=-a;ipFaceCoords_[4](2,0)=-1;
    ipFaceCoords_[4](0,1)= a;ipFaceCoords_[4](1,1)=-a;ipFaceCoords_[4](2,1)=-1; 
    ipFaceCoords_[4](0,2)= a;ipFaceCoords_[4](1,2)= a;ipFaceCoords_[4](2,2)=-1; 
    ipFaceCoords_[4](0,3)=-a;ipFaceCoords_[4](1,3)= a;ipFaceCoords_[4](2,3)=-1; 

    ipFaceCoords_[5](0,0)=-a;ipFaceCoords_[5](1,0)=-a;ipFaceCoords_[5](2,0)= 1;
    ipFaceCoords_[5](0,1)= a;ipFaceCoords_[5](1,1)=-a;ipFaceCoords_[5](2,1)= 1; 
    ipFaceCoords_[5](0,2)= a;ipFaceCoords_[5](1,2)= a;ipFaceCoords_[5](2,2)= 1; 
    ipFaceCoords_[5](0,3)=-a;ipFaceCoords_[5](1,3)= a;ipFaceCoords_[5](2,3)= 1; 

    // Compute shape functions at ip's
    for (int f = 0; f < numFaces_; f++) {
      for (int ip = 0; ip < numFaceIPs_; ip++) {
        double r = ipFaceCoords_[f](0, ip);
        double s = ipFaceCoords_[f](1, ip);
        double t = ipFaceCoords_[f](2, ip);
        for (int iNode = 0; iNode < numEltNodes_; iNode++) {
          double rI = localCoords_(0, iNode);
          double sI = localCoords_(1, iNode);
          double tI = localCoords_(2, iNode);
          NFace_[f](ip, iNode) = 0.125 * (1 + r*rI) * (1 + s*sI) * (1 + t*tI);
        }
      }
    }

    // integration points
    ipFaceWeights_ = 1.0;
  }

  FE_ElementHex::~FE_ElementHex()
  {
    // Nothing to do here
  }

  // -------------------------------------------------------------
  //   shape_function
  // -------------------------------------------------------------
  void FE_ElementHex::shape_function(const int eltID,
                                     DENS_MAT &N,
                                     vector<DENS_MAT> &dN,
                                     DIAG_MAT &weights)
  {
    // Get element node coordinates from mesh
    DENS_MAT xCoords;
    feMesh_->element_coordinates(eltID, xCoords);

    // Transpose xCoords
    DENS_MAT xCoordsT(transpose(xCoords));

    // Shape functions are simply the canonical element values
    N = N_;

    // Set sizes of matrices and vectors
    if ((int)dN.size() != nSD_) dN.resize(nSD_);
    for (int isd = 0; isd < nSD_; isd++) dN[isd].resize(numIPs_, numEltNodes_);
    weights.resize(numIPs_,numIPs_);

    // Create some temporary matrices:

    // Shape function deriv.'s w.r.t. element coords
    DENS_MAT dNdr(nSD_, numEltNodes_, false);

    // Jacobian matrix: [dx/dr dy/ds dz/dt | dx/ds ... ]
    DENS_MAT dxdr, drdx, dNdx;

    // Loop over integration points
    for (int ip = 0; ip < numIPs_; ip++) {
      const double &r = ipCoords_(0,ip);
      const double &s = ipCoords_(1,ip);
      const double &t = ipCoords_(2,ip);

      // Fill dNdr
      for (int inode = 0; inode < numEltNodes_; inode++) {
        const double &rI = localCoords_(0,inode); 
        const double &sI = localCoords_(1,inode); 
        const double &tI = localCoords_(2,inode);
        dNdr(0,inode) = 0.125 * rI * (1 + s*sI) * (1 + t*tI);
        dNdr(1,inode) = 0.125 * sI * (1 + r*rI) * (1 + t*tI);
        dNdr(2,inode) = 0.125 * tI * (1 + r*rI) * (1 + s*sI);
      }

      // Compute dx/dxi matrix
      dxdr = dNdr * xCoordsT;
      drdx = inv(dxdr);

      // Compute dNdx and fill dN matrix
      dNdx = drdx*dNdr;
      for (int isd = 0; isd < nSD_; isd++)
        for (int inode = 0; inode < numEltNodes_; inode++) 
          dN[isd](ip,inode) = dNdx(isd,inode);

      // Compute jacobian determinant of dxdr at this ip
      double J = dxdr(0,0) * ( dxdr(1,1)*dxdr(2,2) - dxdr(2,1)*dxdr(1,2) )
        - dxdr(0,1) * ( dxdr(1,0)*dxdr(2,2) - dxdr(1,2)*dxdr(2,0) )
        + dxdr(0,2) * ( dxdr(1,0)*dxdr(2,1) - dxdr(1,1)*dxdr(2,0) );

      // Compute ip weight
      weights(ip,ip) = ipWeights_(ip)*J; 
    }
  
  }


  //-----------------------------------------------------------------
  void FE_ElementHex::shape_function(const int eltID,
                                     const VECTOR &xi,
                                     DENS_VEC &N)
  {
    N.resize(numEltNodes_);
    for (int inode = 0; inode < numEltNodes_; ++inode) 
    {
      double shp = 0.125;
      for (int i = 0; i < nSD_; ++i) shp *= (1 + xi(i)*localCoords_(i, inode));
      N(inode) = shp;
    }
  }

  //-----------------------------------------------------------------
  void FE_ElementHex::shape_function(const int eltID,
                                     const VECTOR &xi,
                                     DENS_VEC &N,
                                     DENS_MAT &dNdx)
  {
    N.resize(numEltNodes_);
    
    // Get element node coordinates from mesh
    DENS_MAT xCoords;
    feMesh_->element_coordinates(eltID, xCoords);
    DENS_MAT xCoordsT = transpose(xCoords);
    DENS_MAT dNdr(3,8,false), dxdr(3,3,false), drdx;
    double r = xi(0);
    double s = xi(1);
    double t = xi(2);
    for (int inode = 0; inode < numEltNodes_; ++inode) 
    {
      double rI = localCoords_(0,inode); 
      double sI = localCoords_(1,inode); 
      double tI = localCoords_(2,inode);

      // Shape function
      N(inode) = 0.125 * (1.0 + r*rI) * (1.0 + s*sI) * (1.0 + t*tI);

      // Shape function derivative wrt (r,s,t)
      dNdr(0,inode) = 0.125 * rI * (1.0 + s*sI) * (1.0 + t*tI);
      dNdr(1,inode) = 0.125 * sI * (1.0 + r*rI) * (1.0 + t*tI);
      dNdr(2,inode) = 0.125 * tI * (1.0 + r*rI) * (1.0 + s*sI);
    }
    
    // Derivatives wrt (x,y,z)
    dxdr = dNdr * xCoordsT;
    drdx = inv(dxdr);
    dNdx = drdx*dNdr;
  }  

  //-----------------------------------------------------------------
  void FE_ElementHex::face_shape_function(const PAIR & face,
                                          DENS_MAT &N,
                                          vector<DENS_MAT> &dN,
                                          vector<DENS_MAT> &Nn,
                                          DIAG_MAT &weights)
  {
    int eltID = face.first;
    int local_faceID = face.second;
    //  const double * face_normal = face.normal();
    //  double * normal = face_normal; 
    DENS_VEC normal(nSD_);

    // Get element node coordinates from mesh
    DENS_MAT xCoords;
    feMesh_->element_coordinates(eltID, xCoords);

    // Transpose xCoords
    DENS_MAT xCoordsT = transpose(xCoords);

    // Shape functions are simply the canonical element values
    N.reset(numFaceIPs_, numEltNodes_);
    N = NFace_[local_faceID];

    // Set sizes of matrices and vectors
    if ((int)dN.size() != nSD_) dN.resize(nSD_);
    for (int isd = 0; isd < nSD_; isd++) dN[isd].reset(numFaceIPs_, numEltNodes_);
    weights.reset(numFaceIPs_,numFaceIPs_);

    // Create some temporary matrices:

    // Shape function deriv.'s w.r.t. element coords
    DENS_MAT dNdr(nSD_, numEltNodes_);

    // Jacobian matrix: [dx/dr dy/ds dz/dt | dx/ds ... ]
    DENS_MAT dxdr, drdx, dNdx;

    // Loop over integration points
    for (int ip = 0; ip < numFaceIPs_; ip++) {
      double r = ipFaceCoords_[local_faceID](0,ip);
      double s = ipFaceCoords_[local_faceID](1,ip);
      double t = ipFaceCoords_[local_faceID](2,ip);

      // Fill dNdr
      for (int inode = 0; inode < numEltNodes_; inode++) {
        double rI = localCoords_(0,inode); 
        double sI = localCoords_(1,inode); 
        double tI = localCoords_(2,inode);
        dNdr(0,inode) = 0.125 * rI * (1 + s*sI) * (1 + t*tI);
        dNdr(1,inode) = 0.125 * sI * (1 + r*rI) * (1 + t*tI);
        dNdr(2,inode) = 0.125 * tI * (1 + r*rI) * (1 + s*sI);
      }

      // Compute dx/dxi matrix
      dxdr = dNdr * xCoordsT;
      drdx = inv(dxdr);

      // Compute 2d jacobian determinant of dxdr at this ip
      double J = face_normal(face, ip, normal);

      // Compute dNdx and fill dN matrix
      dNdx = drdx*dNdr;
      for (int isd = 0; isd < nSD_; isd++) {
        for (int inode = 0; inode < numEltNodes_; inode++) {
          dN[isd](ip,inode) = dNdx(isd,inode);
          Nn[isd](ip,inode) = N(ip,inode)*normal(isd);
        }
      }

      // Compute ip weight
      weights(ip,ip) = ipFaceWeights_(ip)*J; 
    }
  
  }

  //-----------------------------------------------------------------
  void FE_ElementHex::face_shape_function(const PAIR & face,
                                          DENS_MAT &N,
                                          DENS_MAT &n,
                                          DIAG_MAT &weights)
  {
    int eltID = face.first;
    int local_faceID = face.second;
    //  const double * face_normal = face.normal();
    //  double * normal = face_normal; 
    DENS_VEC normal(nSD_);

    // Get element node coordinates from mesh
    DENS_MAT xCoords;
    feMesh_->element_coordinates(eltID, xCoords);

    // Transpose xCoords
    DENS_MAT xCoordsT = transpose(xCoords);

    // Shape functions are simply the canonical element values
    N.reset(numFaceIPs_, numEltNodes_);
    N = NFace_[local_faceID];

    weights.reset(numFaceIPs_,numFaceIPs_);

    // Create some temporary matrices:

    // Shape function deriv.'s w.r.t. element coords
    DENS_MAT dNdr(nSD_, numEltNodes_);

    // Jacobian matrix: [dx/dr dy/ds dz/dt | dx/ds ... ]
    DENS_MAT dxdr, drdx, dNdx;

    // Loop over integration points
    for (int ip = 0; ip < numFaceIPs_; ip++) {
      double r = ipFaceCoords_[local_faceID](0,ip);
      double s = ipFaceCoords_[local_faceID](1,ip);
      double t = ipFaceCoords_[local_faceID](2,ip);

      // Compute 2d jacobian determinant of dxdr at this ip
      double J = face_normal(face, ip, normal);

      // Copy normal at integration point
      for (int isd = 0; isd < nSD_; isd++) {
        n(isd,ip) = normal(isd);
      }

      // Compute ip weight
      weights(ip,ip) = ipFaceWeights_(ip)*J; 
    }
  
  }

  // -------------------------------------------------------------
  //   face_normal
  // -------------------------------------------------------------
  double FE_ElementHex::face_normal(const PAIR & face,
                                    int ip,
                                    DENS_VEC &normal) 
  {
    // Get element node coordinates from mesh
    DENS_MAT xCoords;
    feMesh_->face_coordinates(face, xCoords); 

    // Transpose xCoords
    DENS_MAT xCoordsT = transpose(xCoords);

    double r = ipFace2DCoords_(0,ip);
    double s = ipFace2DCoords_(1,ip);
    DENS_MAT dNdr(nSD_-1, numFaceNodes_);
    for (int inode = 0; inode < numFaceNodes_; inode++) {
      double rI = localCoords_(0,inode); // NOTE a little dangerous
      double sI = localCoords_(1,inode); 
      dNdr(0,inode) = 0.25 * rI * (1 + s*sI);
      dNdr(1,inode) = 0.25 * sI * (1 + r*rI);
    }
    DENS_MAT dxdr(2,3);

    // Compute dx/dxi matrix
    dxdr = dNdr * xCoordsT;
    normal(0) =  dxdr(0,1)*dxdr(1,2) - dxdr(0,2)*dxdr(1,1) ;
    normal(1) =  dxdr(0,2)*dxdr(1,0) - dxdr(0,0)*dxdr(1,2) ;
    normal(2) =  dxdr(0,0)*dxdr(1,1) - dxdr(0,1)*dxdr(1,0) ;

    double J = sqrt(normal(0)*normal(0)+normal(1)*normal(1)+normal(2)*normal(2));
    //  double inv_J = normal_sense(face)/J;
    double inv_J = 1.0/J;
    normal(0) *= inv_J;
    normal(1) *= inv_J;
    normal(2) *= inv_J;
    return J;
  }


}; // namespace ATC_Transfer

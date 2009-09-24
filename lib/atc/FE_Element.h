#ifndef FE_ELEMENT_H
#define FE_ELEMENT_H

// ATC_Transfer headers
#include "MatrixLibrary.h"
#include "Array2D.h"
#include "ATC_TypeDefs.h"

using namespace std;

namespace ATC {

  // Forward declarations
  class FE_Mesh;

  /**
   *  @class  FE_Element
   *  @brief  Base class for a finite element holding info for canonical element
   */
 
  class FE_Element {
   
  public:
    FE_Element(FE_Mesh *feMesh, 
               int nSD,
               int numEltNodes, 
               int numIPs,
               int numFaces,
               int numFaceNodes,
               int numFaceIPs);
    virtual ~FE_Element();

    /** get number of element nodes */
    int num_elt_nodes() { return numEltNodes_; }
  
    /** get number of integration points */
    int num_ips() { return numIPs_; }

    /** get number of faces */
    int num_faces() { return numFaces_; }
  
    /** get number of face nodes */
    int num_face_nodes() { return numFaceNodes_; }
  
    /** get number of integration points */
    int num_face_ips() { return numFaceIPs_; }

    /** cannonical face numbering */
    const Array2D<int> & local_face_conn(void) const {return localFaceConn_ ;}

    /** 
     *  compute shape functions at all ip's:
     *  indexed: N(ip,node)
     *           dN[nsd](ip,node)
     *           weights(ip)
     */
    virtual void shape_function(const int eltID,
                                DENS_MAT &N,
                                vector<DENS_MAT> &dN, 
                                DIAG_MAT &weights)=0;

    /** 
     *  compute shape functions at a single point, given the local coordinates
     *  indexed: N(node)
     *           dN[nsd](node)
     */
    virtual void shape_function(const int eltID,
                                const VECTOR &xi,
                                DENS_VEC &N) = 0;

    virtual void shape_function(const int eltID,
                                const VECTOR &xi,
                                DENS_VEC &N,
                                DENS_MAT &dN) = 0;

    /** 
     *  compute shape functions at all face ip's:
     *  indexed: N(ip,node)
     *           dN[nsd](ip,node)
     *           Nn[nsd](ip,node)
     *           weights(ip)
     */
    virtual void face_shape_function(const PAIR       &face,
                                     DENS_MAT         &N,
                                     vector<DENS_MAT> &dN,
                                     vector<DENS_MAT> &Nn,
                                     DIAG_MAT         &weights) = 0;

    virtual void face_shape_function(const PAIR & face,
                                          DENS_MAT &N,
                                          DENS_MAT &n,
                                          DIAG_MAT &weights) = 0;

    virtual double face_normal(const PAIR & face,
                               const int ip,
                               DENS_VEC &normal) = 0;

    enum FE_ElementQuadrature { GAUSSIAN_QUADRATURE, NODAL_QUADRATURE};
    virtual void set_quadrature(int quadrature_type) = 0;

  protected:

    // Mesh
    FE_Mesh * feMesh_;

    // Number of spatial dimensions
    int nSD_;

    // Number of element nodes
    int numEltNodes_;

    // Number of integration points
    int numIPs_;

    // Number of faces
    int numFaces_;

    // Number of face nodes
    int numFaceNodes_;

    // Number of face integration points
    int numFaceIPs_;

    // local coords of nodes: localCoords_(isd, ip)
    DENS_MAT localCoords_;

    // quadrature scheme: localCoords_(isd, ip)
    DENS_MAT ipCoords_; // local coordinates

    // matrix of shape functions at ip's: N_(ip, node)
    DENS_MAT N_;

    // integration point weights: ipWeights_(ip)
    DENS_VEC ipWeights_; // local coordinates

    // local face numbering
    Array2D<int>  localFaceConn_;

    // quadrature scheme: localCoords_(isd, ip)
    vector<DENS_MAT> ipFaceCoords_; 
    DENS_MAT ipFace2DCoords_; 

    // matrix of shape functions at ip's: N_(ip, node)
    vector<DENS_MAT> NFace_;

    // integration point weights: ipWeights_(ip)
    DENS_VEC ipFaceWeights_; 

  };



  /**
   *  @class  FE_ElementHex
   *  @author Greg Wagner
   *  @brief  3D, linear 8-node hex element with 2x2x2 quadrature
   */
  class FE_ElementHex : public FE_Element {

  public:
    FE_ElementHex(FE_Mesh * feMesh);
    // Dump state info to disk for later restart
    void write_restart(FILE *);
    ~FE_ElementHex();

  protected:
    virtual void shape_function(const int eltID,
                                DENS_MAT &N,
                                vector<DENS_MAT> &dN,
                                DiagonalMatrix<double> &weights);

    /** 
     *  compute shape functions at a single point, given the local coordinates
     *  indexed: N(node)
     *           dN[nsd](node)
     */
    virtual void shape_function(const int eltID,
                                const VECTOR &xi,
                                DENS_VEC &N);

    virtual void shape_function(const int eltID,
                                const VECTOR &xi,
                                DENS_VEC &N,
				DENS_MAT &dN);

    virtual void face_shape_function(const PAIR &face,
                                     DENS_MAT &N,
                                     vector<DENS_MAT> &dN,
                                     vector<DENS_MAT> &Nn,
                                     DiagonalMatrix<double> &weights);

    virtual void face_shape_function(const PAIR & face,
                                          DENS_MAT &N,
                                          DENS_MAT &n,
                                          DIAG_MAT &weights); 

    virtual double face_normal(const PAIR &face, const int ip, DENS_VEC &normal);

    virtual void set_quadrature(int quadrature_type);
  };

}; // namespace ATC_Transfer

#endif // FE_ELEMENT_H

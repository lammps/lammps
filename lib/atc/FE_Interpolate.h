#ifndef FE_INTERPOLATE_H
#define FE_INTERPOLATE_H

// ATC_Transfer headers
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

using namespace std;

namespace ATC {

  // Forward declarations
  struct FE_Quadrature;
  class FE_Element;

  /**
   *  @class  FE_Interpolate
   *  @brief  Base class for computing shape functions, nested inside FE_Element
   */
  class FE_Interpolate {

  public:
    FE_Interpolate(FE_Element *feElement);

    virtual ~FE_Interpolate();
    
    /** compute the quadrature for a given element type */
    void set_quadrature(FeEltGeometry geo, FeIntQuadrature quad);

    /** store N_ and dNdr_ (and eventually dNdrFace_) the given
     *  quadrature */
    virtual void precalculate_shape_functions();

    /** compute the shape functions at the given point; 
     *  we use CLON_VECs */
    virtual void compute_N(const VECTOR &point,
                           VECTOR &N) = 0;
 
    // middle args get no names because they won't be used by some
    // child classes.
    /** compute the shape function derivatives at the given point; 
     *  generic, so can work for 2D or 3D case */
    virtual void compute_dNdr(const VECTOR &point,
                              const int,
                              const int,
                              const double,
                              DENS_MAT &dNdr) = 0;

    /** compute both shape functions and derivatives, using the
     *  fastest, shortcuttiest methods available for batch use */
    virtual void compute_N_dNdr(const VECTOR &point,
                                VECTOR &N,
                                DENS_MAT &dNdr) const = 0;

    /** compute shape functions at a single point, given the local 
     *    coordinates; eltCoords needed for derivatives:
     *  indexed: N(node)
     *           dN[nsd](node) */
    virtual void shape_function(const VECTOR &xi,
                                DENS_VEC &N);

    virtual void shape_function(const DENS_MAT &eltCoords,
                                const VECTOR &xi,
                                DENS_VEC &N,
                                DENS_MAT &dNdx);
    virtual void shape_function_derivatives(const DENS_MAT &eltCoords,
                                const VECTOR &xi,
                                DENS_MAT &dNdx);

    /** compute shape functions at all ip's:
     *  indexed: N(ip,node)
     *           dN[nsd](ip,node)
     *           weights(ip) */
    virtual void shape_function(const DENS_MAT &eltCoords,
                                DENS_MAT &N,
                                vector<DENS_MAT> &dN, 
                                DIAG_MAT &weights);

    /** compute shape functions at all face ip's:
     *  indexed: N(ip,node)
     *           dN[nsd](ip,node)
     *           n(ip,node)/Nn[nsd](ip,node)
     *           weights(ip) */
    virtual void face_shape_function(const DENS_MAT &eltCoords,
                                     const DENS_MAT &faceCoords,
                                     const int faceID,
                                     DENS_MAT &N,
                                     DENS_MAT &n,
                                     DIAG_MAT &weights);

    virtual void face_shape_function(const DENS_MAT &eltCoords,
                                     const DENS_MAT &faceCoords,
                                     const int faceID,
                                     DENS_MAT &N,
                                     vector<DENS_MAT> &dN,
                                     vector<DENS_MAT> &Nn,
                                     DIAG_MAT &weights);

    /** compute unit normal vector for a face:
     *  indexed: N(ip,node)
     *           dN[nsd](ip,node)
     *           n(ip,node)/Nn[nsd](ip,node)
     *           weights(ip) */
    virtual double face_normal(const DENS_MAT &faceCoords,
                               int ip,
                               DENS_VEC &normal);

    /** compute tangents to coordinate lines
     *  indexed:             */
    virtual void tangents(const DENS_MAT &eltCoords,
                          const VECTOR &localCoords,
                          DENS_MAT &dxdr) const;
    virtual void tangents(const DENS_MAT &eltCoords,
                          const VECTOR &localCoords,
                          vector<DENS_VEC> &t,
                          const bool normalize = false) const;
    /** get number of ips in scheme */
    int num_ips() const;

    /** get number of ips per face */
    int num_face_ips() const;

  protected:
    // owner element
    FE_Element *feElement_;

    // Number of dimensions
    int nSD_;

    map<FeIntQuadrature,FE_Quadrature *> feQuadList_;
    FE_Quadrature *feQuad_;

    // matrix of shape functions at ip's: N_(ip, node)
    DENS_MAT N_;

    
    vector<DENS_MAT> dNdr_;

    // matrix of shape functions at ip's: N_(ip, node)
    vector<DENS_MAT> NFace_;

    // matrix of generic face shape function derivatives
    vector<vector<DENS_MAT> > dNdrFace_;

    // matrix of generic face shape function derivatives
    vector<DENS_MAT> dNdrFace2D_;

  };


  /**********************************************
   * Class for linear interpolation functions with
   * Cartesian coordinates
   **********************************************/
  class FE_InterpolateCartLagrange : public FE_Interpolate {

  public:
    FE_InterpolateCartLagrange(FE_Element *feElement);

    virtual ~FE_InterpolateCartLagrange();

    virtual void compute_N(const VECTOR &point,
                           VECTOR &N);
  
    virtual void compute_dNdr(const VECTOR &point,
                              const int numNodes,
                              const int nD,
                              const double,
                              DENS_MAT &dNdr);

    virtual void compute_N_dNdr(const VECTOR &point,
                                VECTOR &N,
                                DENS_MAT &dNdr) const;

  };


  /**********************************************
   * Class for linear elements with
   * Cartesian coordinates
   **********************************************/
  class FE_InterpolateCartLin : public FE_Interpolate {

  public:
    FE_InterpolateCartLin(FE_Element *feElement);

    virtual ~FE_InterpolateCartLin();

    virtual void compute_N(const VECTOR &point,
                           VECTOR &N);
  
    virtual void compute_dNdr(const VECTOR &point,
                              const int numNodes,
                              const int nD,
                              const double vol,
                              DENS_MAT &dNdr);

    virtual void compute_N_dNdr(const VECTOR &point,
                                VECTOR &N,
                                DENS_MAT &dNdr) const;

  };


  /**********************************************
   * Class for quadratic serendipity elements with
   * Cartesian coordinates
   **********************************************/
  class FE_InterpolateCartSerendipity : public FE_Interpolate {

  public:
    FE_InterpolateCartSerendipity(FE_Element *feElement);

    virtual ~FE_InterpolateCartSerendipity();

    virtual void compute_N(const VECTOR &point,
                           VECTOR &N);
  
    virtual void compute_dNdr(const VECTOR &point,
                              const int numNodes,
                              const int nD,
                              const double vol,
                              DENS_MAT &dNdr);

    virtual void compute_N_dNdr(const VECTOR &point,
                                VECTOR &N,
                                DENS_MAT &dNdr) const;

  };


  /**********************************************
   * Class for linear interpolation functions with
   * volumetric coordinates
   **********************************************/
  class FE_InterpolateSimpLin : public FE_Interpolate {

  public:
    // "Simp"ly overrides all standard shape function methods
    FE_InterpolateSimpLin(FE_Element *feElement);

    virtual ~FE_InterpolateSimpLin();

    virtual void compute_N(const VECTOR &point,
                           VECTOR &N);
  
    // middle args get no names because they won't be used by some
    // child classes.
    virtual void compute_dNdr(const VECTOR &,
                              const int,
                              const int,
                              const double,
                              DENS_MAT &dNdr);

    virtual void compute_N_dNdr(const VECTOR &point,
                                VECTOR &N,
                                DENS_MAT &dNdr) const;

  };


}; // namespace ATC

#endif // FE_INTERPOLATE_H

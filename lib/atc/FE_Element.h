#ifndef FE_ELEMENT_H
#define FE_ELEMENT_H

#include <vector>
#include <string>

// ATC_Transfer headers
#include "MatrixLibrary.h"
#include "Array2D.h"
#include "ATC_TypeDefs.h"

namespace ATC {

  enum ProjectionGuessType {
    COORDINATE_ALIGNED=0, 
    CENTROID_LINEARIZED,
    TWOD_ANALYTIC};

  // Forward declarations
  class FE_Interpolate;

  /**
   *  @class  FE_Element
   *  @brief  Base class for a finite element holding info for canonical element
   */
 
  class FE_Element {
   
  public:
    
    ///////////////////////////////////////////////////////////////////////////
    //
    // CONSTRUCTOR AND DESTRUCTOR 

    FE_Element(const int nSD,
               int numFaces,  
               int numNodes,
               int numFaceNodes,
               int numNodes1d);

    virtual ~FE_Element();
    
    ///////////////////////////////////////////////////////////////////////////
    //
    // GETTERS

    /** get number of spatial dimensions (almost always 3) */
    int num_dims() { return nSD_; }

    /** get number of element nodes */
    int num_elt_nodes() { return numNodes_; }
  
    /** get number of element nodes */
    int num_elt_nodes_1d() { return numNodes1d_; }
  
    /** get number of faces */
    int num_faces() { return numFaces_; }
  
    /** get number of face nodes */
    int num_face_nodes() { return numFaceNodes_; }
 
    // Getters for FE_Interpoate to have access to coordinates and connectivity
    /** get canonical coordinates */
    const DENS_MAT &local_coords() const { return localCoords_; }

    /** get canonical coordinates in 1d */
    DENS_VEC local_coords_1d() const; 
    
    /** get canonical connectivity of nodes and faces */
    const Array2D<int> &local_face_conn() const { return localFaceConn_; }

    /** return volume of the element */
    double vol() const { return vol_; }

    /** return area of a face */
    double face_area() const { return faceArea_; }

    // the following two are pass-throughs to the interpolate class, and
    // can thus only be declared in the class body (or else the
    // interpolate class is "incomplete" and cannot be referenced)
    
    /** get number of integration points */
    int num_ips() const;

    /** get number of integration points */
    int num_face_ips() const;
  
    /** order of interpolation */
    int order() const {return numNodes1d_;} 

    /** compute the quadrature for a given element type */
    virtual void set_quadrature(FeIntQuadrature type) = 0; 

    /** return the set of 1d nodes that correspond to this node in 3d space */
    void mapping(const int inode, std::vector<int> &mapping) const;

    /** extract face coordinates from element coordinates */
    void face_coordinates(const DENS_MAT &eltCoords, 
                          const int faceID,
                          DENS_MAT &faceCoords) const;

    /** set initial guess type for point in element search */ 
    void set_projection_guess(ProjectionGuessType type) 
       { projectionGuess_ = type;}
    ///////////////////////////////////////////////////////////////////////////
    //
    // GENERIC ELEMENT COMPUTATIONS

    /** compute local coordinates from global */
    virtual bool local_coordinates(const DENS_MAT &eltCoords,
                                   const DENS_VEC &x,
                                   DENS_VEC &xi) const;

    /** location of local coordinates (0,0,0) */
    virtual void centroid(const DENS_MAT &eltCoords,
                          DENS_VEC & centroid) const;

    /** test if a specified element actually contains the given point */
    virtual bool contains_point(const DENS_MAT &eltCoords,
                                const DENS_VEC &x) const;

    /** check if element bounding box contains the given point */
    bool range_check(const DENS_MAT &eltCoords, const DENS_VEC & x) const;

    /** get the min and max coordinate of any point in an element in a 
     *  dimension */
    void bounds_in_dim(const DENS_MAT &eltCoords, const int dim,
                       double &min, double &max) const;

    ///////////////////////////////////////////////////////////////////////////
    //
    //PASS-THROUGHS TO INTERPOLATE CLASS
    virtual void shape_function(const VECTOR & xi,
                                DENS_VEC &N) const;

    /** 
     *  compute shape functions at all ip's:
     *  indexed: N(ip,node)
     *           dN[nsd](ip,node)
     *           weights(ip)
     */
    virtual void shape_function(const DENS_MAT eltCoords,
                                DENS_MAT &N,
                                std::vector<DENS_MAT> &dN, 
                                DIAG_MAT &weights);

    /** 
     *  compute shape functions and derivatives at a single point, 
     *    given the point and the element that contains it
     *  indexed: N(node)
     */
    virtual void shape_function(const DENS_MAT eltCoords,
                                const VECTOR &x,
                                DENS_VEC &N);

    /** 
     *  compute shape functions and derivatives at a single point, 
     *    given the point and the element that contains it
     *  indexed: N(node)
     *           dNdx(ip,nSD)
     */
    virtual void shape_function(const DENS_MAT eltCoords,
                                const VECTOR &x,
                                DENS_VEC &N,
                                DENS_MAT &dNdx);

    /** 
     *  compute shape functions and derivatives at a single point, 
     *    given the point and the element that contains it
     *  indexed: 
     *           dNdx(ip,nSD)
     */
    virtual void shape_function_derivatives(const DENS_MAT eltCoords,
                                const VECTOR &x,
                                DENS_MAT &dNdx);

    /** 
     *  compute shape functions at all face ip's:
     *  indexed: N(ip,node)
     *           n[nsd](ip,node)
     *           weights(ip)
     */
    virtual void face_shape_function(const DENS_MAT &eltCoords,
                                     const int faceID,
                                     DENS_MAT &N,
                                     DENS_MAT &n,
                                     DIAG_MAT &weights);

    /** 
     *  compute shape functions at all face ip's:
     *  indexed: N(ip,node)
     *           dN[nsd](ip,node)
     *           Nn[nsd](ip,node)
     *           weights(ip)
     */
    virtual void face_shape_function(const DENS_MAT &eltCoords,
                                     const int faceID,
                                     DENS_MAT &N,
                                     std::vector<DENS_MAT> &dN,
                                     std::vector<DENS_MAT> &Nn,
                                     DIAG_MAT &weights);

    /** 
     *  compute normal vector from the specified face
     *  indexed: normal(nSD)
     */
    virtual double face_normal(const DENS_MAT &eltCoords,
                               const int faceID,
                               int ip,
                               DENS_VEC &normal); 
  
    /** 
     *  compute tangents to local coordinates
     *  indexed: 
     */
    virtual void tangents(const DENS_MAT &eltCoords,
                          const DENS_VEC &x,
                          std::vector<DENS_VEC> & tangents,
                          const bool normalize=false) const; 

  protected:

    ///////////////////////////////////////////////////////////////////////////
    //
    // HELPERS
    
    /** 
     *  generate the appropriate interpolation class 
     */
    FE_Interpolate *interpolate_factory(std::string interpolateType);

    /** initial guess for local coordinates */
    virtual void initial_local_coordinates(const DENS_MAT &eltCoords,
                                           const DENS_VEC &x,
                                           DENS_VEC &xiInitial) const;
    
    ///////////////////////////////////////////////////////////////////////////
    //
    // PROTECTED MEMBERS
    
    // Currently used interpolation class
    FE_Interpolate *feInterpolate_;

    // Number of spatial dimensions
    int nSD_;
    // Number of faces, used for generic contains_point
    int numFaces_;

    // Number of element nodes
    int numNodes_;
    // Number of face nodes
    int numFaceNodes_;
    // Number of nodes in one dimension
    int numNodes1d_; 


    // local coords of nodes: localCoords_(isd, ip)
    DENS_MAT localCoords_;

    // local face numbering
    Array2D<int> localFaceConn_;

    // volume of canonical element
    double vol_;

    // area of faces of canonical element
    double faceArea_;

    /** tolerance used in solving Newton's method for local coordinates */
    double tolerance_;
    ProjectionGuessType projectionGuess_;
    
  };


  /**
   *  @class  FE_ElementHex
   *  @author Sean Laguna
   *  @brief  3D, linear 8-node hex element
   */
  class FE_ElementHex : public FE_Element {

  public:
    FE_ElementHex(int numNodes,
                  int numFaceNodes,
                  int numNodes1d);

    // Dump state info to disk for later restart (unimplemented)
    void write_restart(FILE *);

    ~FE_ElementHex();

    void set_quadrature(FeIntQuadrature type); 

    bool contains_point(const DENS_MAT &eltCoords,
                        const DENS_VEC &x) const;

  };


  /**
   *  @class  FE_ElementRect
   *  @author Greg Wagner, amended by Sean Laguna
   *  @brief  3D, linear 8-node rectilinear hex element
   */
  class FE_ElementRect : public FE_ElementHex {

  public:
    FE_ElementRect();

    // Dump state info to disk for later restart (unimplemented)
    void write_restart(FILE *);

    ~FE_ElementRect();
    
    bool local_coordinates(const DENS_MAT &eltCoords,
                           const DENS_VEC &x,
                           DENS_VEC &xi) const;

   protected:
    virtual bool contains_point(const DENS_MAT &eltCoords,
                                const DENS_VEC &x) const;
  
  };


  /**
   *  @class  FE_ElementTet
   *  @author Aaron Gable & Sean Laguna
   *  @brief  3D, linear 4-node tetrahedral element
   */
  class FE_ElementTet : public FE_Element {

  public:
    FE_ElementTet(int numNodes,
                  int numFaceNodes,
                  int numNodes1d);

    // Dump state info to disk for later restart (unimplemented)
    void write_restart(FILE *);
    
    ~FE_ElementTet();

    void set_quadrature(FeIntQuadrature type); 

    bool local_coordinates(const DENS_MAT &eltCoords,
                           const DENS_VEC &x,
                           DENS_VEC &xi) const;

    bool contains_point(const DENS_MAT &eltCoords,
                        const DENS_VEC &x) const;

  };

}; // namespace ATC

#endif // FE_ELEMENT_H

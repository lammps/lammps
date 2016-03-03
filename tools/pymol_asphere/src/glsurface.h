/***************************************************************************
                                 glsurface.h
                               W. Michael Brown
                             -------------------

  Store and manipulate surfaces as OpenGL primitives
  
    begin                : Sun Jun 8 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#ifndef GLSURFACE_H
#define GLSURFACE_H

#include "miscm.h"
#include "cartesian.h"
#include "colors.h"
//#include "gridnum.h"

#include <vector>
#include <algorithm>
#include <fstream>
using namespace std;

/// Stores the ID numbers of the vertices making a triangle
struct Triangle {
  /// Triangle position
  unsigned point[3];
};

/// Vertex information used for all primitives
struct Vertex {
  /// True if normal is being set for vertex
  bool valid_normal;    
  /// Cartesian coords
  cPt cpt;
  // Used for surface area calculations
  vectorPt normal;
  /// Color
  colorPt color;
  /// For coloring surface
  char marker;          
  /// For transparency (0 transparent, 1 opaque)
  double transparency;
};

/// ID numbers for vertices making a line
struct GLline {
  unsigned points[2];
};

// This stuff is used for removing duplicate lines by sorting;
bool operator < (const GLline &one, const GLline &two);
bool operator == (const GLline &one, const GLline &two);

/// Line Strip (loop) primitive (stores vertex indices)
struct LineStrip {
  /// True for loop, false for line_strip
  bool loop;
  vector<unsigned> line;
};

/// Sphere primitive
struct Sphere {
  unsigned i;
  double radius;
};

/// Class for storage and I/O of a set of OpenGL primitives
class GLSurface {
public: 
  GLSurface();
	~GLSurface();
  /// Reserve vector space for specified number of vertices and triangles
  void reserve(unsigned num_v, unsigned num_t);
  /// Delete all primitives
  void clear();
    
  /// Add a triangle
  void addtriangle(Triangle &t);
  /// Add a vertex
  void addvertex(Vertex &v);

  /// Add many triangles (via swap)
  /** \note triangles are removed from input vector */
  void addtriangles(vector<Triangle> &t);
  /// Add many lines (via swap)
  /** \note lines are removed from input vector */
  void addlinestrips(vector<LineStrip> &l);
  /// Add a mesh (via swap)
  /** \note lines are removed from input vector */
  void addxyzmesh(vector<GLline> &x,vector<GLline> &y,vector<GLline> &z);
  /// Add a line
  void addline(unsigned v1, unsigned v2);
  /// Add a sphere based on a current vertex and a radius
  void add_sphere(unsigned index, double radius) {
    Sphere s;
    s.i=index;
    s.radius=radius;
    spheres.push_back(s);
  }
  
  /// Add an ellipsoid with specified resolution
  /** The number of triangles is equal to twice the resolution squared **/
  void add_ellipsoid(const cPt &center, const cPt &rad,
                       const Quaternion &rot, const colorPt &color,
                       const double alpha,const unsigned resolution);
  /// Add an ellipsoid with resolution of 10
  /** \sa add_ellipsoid **/
  void add_ellipsoid(const cPt &center, const cPt &rad,
                       const Quaternion &rot, const colorPt &color,
                       const double alpha);

  /// Calculate the unit normals for any vertices marked with invalid normals
  void calculatenormals();
  /// Normalize normals.
  void normalize();
  /// Flip normals
  void flipnormals();
  /// Calculate the Area of the Surface
  double surfacearea();
  /// Calculate the Area of a Triangle
  double triarea(Triangle &t);
  
  /// Return the number of spheres
  unsigned size_spheres() { return spheres.size(); }
  /// Return the number of vertices
  unsigned size_vertices() { return vertices.size(); }
  
  /// Solid color a surface
  void color(const colorPt &color);
  /// Color a surface by using a gradient
  void colorbygradient(const colorPt &start, const colorPt &end);
  /// Color a surface based on a marker
  void colorbymarker(const colorPt &zero, const colorPt &nonzero);
  /// Color a surface based on interpolation of values on a grid
  /** \param grid Grid with appropriate interpolation scheme set
    * \param neg Color for the value at and below minv
    * \param mid Color for the value at midv
    * \param pos Color for the value at and above posv **/
  //void colorbygrid(GridNUM<double> &grid,const colorPt &neg,
  //                 const colorPt &mid, const colorPt &pos, double minv,
  //                 double midv, double maxv);
  /// Set the transparency of the entire surface (0-1 [1=no transparency])
  void set_transparency(double alpha);
  
  /// Write out triangles as BEGIN,TRIANGLES,END primitives
  void writetris(ofstream &out, const string &objname);
  /// Write out triangles as TRIANGLE primitives
  void writetris_surf(ofstream &out, const string &objname);
  /// Write out triangles as POINTS
  void writetris_points(ofstream &out, const string &objname);
  /// Write out vertices as spheres
  void write_vspheres(ofstream &out, const string &objname, double radius);
  /// Write out triangles as Mesh
  void writetris_mesh(ofstream &out, const string &objname);

  /// Write out triangles as an xyzmesh
  void writexyzmesh(ofstream &out,const string &objname);

  void writelines(ofstream &out, const string &objname);
  void writelines(vector<GLline> &l, ofstream &out, const string &objname);
  /// Write out Line Strips (Loops)
  void writelinestrips(ofstream &out, const string &objname);
  void writespheres(ofstream &out, const string &objname);
  
  /// Write header for Python scripts for rendering in PyMol
  void writepymolheader(ofstream &out);
  /// Write tail for Python script for rendering in Pymol
  void writepymoltail(ofstream &out, const string &objname);
private:
  vector<Triangle> triangles;
  vector<Vertex> vertices;
  vector<GLline> gllines;
  vector<Sphere> spheres;
  
  // For xyzmesh
  vector<GLline> xmesh;
  vector<GLline> ymesh;
  vector<GLline> zmesh;
  
  vector<LineStrip> linestrips;

  void add_super_ellipsoid(const cPt &cen, const cPt &rad,const Quaternion &q,
                           const double n, const double e, const double u1,
                           const double u2, const double v1, const double v2,
                           const unsigned u_segs, const unsigned v_segs,
                           const colorPt &color, const double alpha);
  void SQE_helper(Vertex &ver, const cPt &rad, const double u,
                  const double v, const double n, const double e);
};

#endif


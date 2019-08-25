/***************************************************************************
                                glsurface.cpp
                               W. Michael Brown
                             -------------------

  Store and manipulate surfaces as OpenGL primitives

    begin                : Sun Jun 8 2003
    copyright            : (C) 2003 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

#include "glsurface.h"

GLSurface::GLSurface(){
}
GLSurface::~GLSurface(){
}

// Operations on GLline structure
bool operator < (const GLline &one, const GLline &two) {
  unsigned minone, mintwo, maxone, maxtwo;

  if (one.points[0]<one.points[1]) {
    minone=one.points[0];
    maxone=one.points[1];
  } else {
    minone=one.points[1];
    maxone=one.points[0];
  }
  if (two.points[0]<two.points[1]) {
    mintwo=two.points[0];
    maxtwo=two.points[1];
  } else {
    mintwo=two.points[1];
    maxtwo=two.points[0];
  }

  if (minone<mintwo)
    return true;
  else if (minone==mintwo)
    if (maxone<maxtwo)
      return true;

  return false;
}

bool operator == (const GLline &one, const GLline &two) {
  if (one.points[0]==two.points[0] && one.points[1]==two.points[1])
    return true;
  if (one.points[0]==two.points[1] && one.points[1]==two.points[0])
    return true;
  return false;
}

void GLSurface::clear() {
  gllines.clear();
  linestrips.clear();
  triangles.clear();
  vertices.clear();
  xmesh.clear();
  ymesh.clear();
  zmesh.clear();
}

// Reserve vector space for specified number of vertices and triangles
void GLSurface::reserve(unsigned num_v, unsigned num_t) {
  vertices.reserve(num_v); triangles.reserve(num_t);
}

// Add a triangle
void GLSurface::addtriangle(Triangle &t) {
  triangles.push_back(t);
}

// Add a vertex
void GLSurface::addvertex(Vertex &v) {
  vertices.push_back(v);
}

// Add many triangles (via swap)
void GLSurface::addtriangles(vector<Triangle> &t) {
  triangles.swap(t);
}

// Add many lines (via swap)
void GLSurface::addlinestrips(vector<LineStrip> &l) {
  linestrips.swap(l);
}

// Add a mesh (via swap)
void GLSurface::addxyzmesh(vector<GLline> &x,vector<GLline> &y,
                           vector<GLline> &z) {
  xmesh.swap(x); ymesh.swap(y); zmesh.swap(z);
}

// Preconditions:
//   Triangles all loaded with correct vertex ordering
//   Vertices all loaded
//   Vertices that need to be calculated are marked as invalid
// Postconditions:
//   For any normals that are marked as invalid, calculate normals by
//     averaging normals for each triangle
void GLSurface::calculatenormals() {
  unsigned i;

  // Set all invalid normals to zero
  for (i=0; i<vertices.size(); i++) 
    if (vertices[i].valid_normal==false)
      vertices[i].normal=vectorPt(0,0,0);

  for (i=0; i<triangles.size(); i++) {
    vectorPt vec1, vec2, normal;
		unsigned id0, id1, id2;
    id0=triangles[i].point[0];
		id1=triangles[i].point[1];
		id2=triangles[i].point[2];
    vec1=vertices[id1].cpt-vertices[id0].cpt;
    vec2=vertices[id2].cpt-vertices[id0].cpt;
		normal[X]=vec1[Z]*vec2[Y]-vec1[Y]*vec2[Z];
		normal[Y]=vec1[X]*vec2[Z]-vec1[Z]*vec2[X];
		normal[Z]=vec1[Y]*vec2[X]-vec1[X]*vec2[Y];
    // Only do vertex averaging if the normal has not been
    //   explicitly set
    if (vertices[id0].valid_normal==false)
      vertices[id0].normal+=normal;
    if (vertices[id1].valid_normal==false)
      vertices[id1].normal+=normal;
    if (vertices[id2].valid_normal==false)
      vertices[id2].normal+=normal;
  }

  // Mark all normals as valid
  for (i=0; i<vertices.size(); i++)
    if (vertices[i].valid_normal==false) {
     vertices[i].valid_normal=true;
     double no=vertices[i].normal.norm();
     if (no==0)
       vertices[i].normal=ORIGIN;
     else
       vertices[i].normal/=no*-1;
  }
  //error->note[31] << "Calculated normals: " << flush;

  return;
}
    
// Normalize normals.
void GLSurface::normalize() {
  unsigned i;

  for (i=0; i<vertices.size(); i++) {
    vertices[i].normal.normalize();
    if (am::not_nan(vertices[i].normal.x)==false) {
        cerr << "Bad normal!"; // cerr fix me
      vertices[i].normal=vectorPt(0,0,0);
    }
	}
  //error->note[31] << "Normalized normals: " << flush;
}

// Flip normals
void GLSurface::flipnormals() {
  unsigned i;

  for (i=0; i<vertices.size(); i++)
    vertices[i].normal*=-1.0;
  //error->note[31] << "Flipped normals: " << flush;
}

// Calculate the Area of the Surface
double GLSurface::surfacearea() {
  double sa=0;

  for (unsigned i=0; i<triangles.size(); i++) 
    sa+=triarea(triangles[i]);

  return sa;
}

double GLSurface::triarea(Triangle &t) {
  double a=vertices[t.point[0]].cpt.dist(vertices[t.point[1]].cpt);
  double b=vertices[t.point[1]].cpt.dist(vertices[t.point[2]].cpt);
  double c=vertices[t.point[0]].cpt.dist(vertices[t.point[2]].cpt);
  double d=a+b+c;
  
  return (sqrt(d*(d-a)*(d-b)*(d-c))/4);
}

// Color surfaces
void GLSurface::color(const colorPt &color) {
  unsigned i;

  for (i=0; i<vertices.size(); i++)
		vertices[i].color=color;
}
 
void GLSurface::colorbygradient(const colorPt &start, const colorPt &end) {
  colorPt step, current;
  unsigned i;
  
  step=(end-start)/double(vertices.size());
  current=start;
  for (i=0; i<vertices.size(); i++) {
		vertices[i].color=current;
		current+=step;
	}
}

void GLSurface::colorbymarker(const colorPt &zero, const colorPt &nonzero) {
  unsigned i;

  for (i=0; i<vertices.size(); i++) {
    if (vertices[i].marker==0)
      vertices[i].color=zero;
    else
      vertices[i].color=nonzero;
  }
  return;
}

// Color a surface based on interpolation of values on a grid
// void GLSurface::colorbygrid(GridNUM<double> &grid,const colorPt &neg,
//                             const colorPt &mid, const colorPt &pos, double minv,
//                             double midv, double maxv) {
//   Colors c;
//   for (unsigned i=0; i<vertices.size(); i++)
//     vertices[i].color=c.gradient(minv,midv,maxv,grid.interp(vertices[i].cpt),
//                                  neg,mid,pos);
// }    

void GLSurface::set_transparency(double alpha) {
	for (unsigned i=0; i<vertices.size(); i++)
		vertices[i].transparency=alpha;
}

// Write out triangles as BEGIN,TRIANGLES,END primitives
void GLSurface::writetris(ofstream &out, const string &objname) {
  unsigned i,j,k;
  int index;

  writepymolheader(out);
  out << "BEGIN, TRIANGLES," << endl;
  for (i=0; i<triangles.size(); i++) {
    // Output coords, colors, and normals of points
    for (j=0; j<3; j++) {
      index=triangles[i].point[j];
//      if (vertices[index].transparency!=1)
        out << "ALPHA," << vertices[index].transparency << ",";
      out << "COLOR,";
      for (k=X; k<=Z; k++)
        out << vertices[index].color[k] << ",";
      if (vertices[index].valid_normal) {
        out << "NORMAL,";
        for (k=RED; k<=BLUE; k++)
          out << vertices[index].normal[k] << ",";
      }
      out << "VERTEX,";
      for (k=X; k<=Z; k++)
        out << vertices[index].cpt[k] << ",";
    }
    out << endl;
  }
  out << "END";
  writepymoltail(out,objname);
}

// Write out triangles as TRIANGLE primitives
void GLSurface::writetris_surf(ofstream &out, const string &objname) {
  unsigned i,j,k;
  int index;

  writepymolheader(out);
  for (i=0; i<triangles.size(); i++) {
//    if (vertices[triangles[i].point[0]].transparency!=1)
      out << "ALPHA," << vertices[triangles[i].point[0]].transparency << ",";
    out << "TRIANGLE,";
    // Output coords of points
    for (j=0; j<3; j++) {
      index=triangles[i].point[j];
      for (k=X; k<=Z; k++)
        out << vertices[index].cpt[k] << ",";
    }
    // Output normals of points
    for (j=0; j<3; j++) {
      index=triangles[i].point[j];
      for (k=X; k<=Z; k++)
        out << vertices[index].normal[k] << ",";
    }
    // Output colors of points
    for (j=0; j<3; j++) {
      index=triangles[i].point[j];
      for (k=RED; k<=BLUE; k++)
        out << vertices[index].color[k] << ",";
    }
    out << endl;
  }
  writepymoltail(out,objname);
}

// Write out triangles as POINTS
void GLSurface::writetris_points(ofstream &out, const string &objname) {
  unsigned i,k;

  writepymolheader(out);
  out << "BEGIN, POINTS," << endl;
  for (i=0; i<vertices.size(); i++) {
    // Output coords, colors, and normals of points
//    if (vertices[i].transparency!=1)
      out << "ALPHA," << vertices[i].transparency << ",";
    out << "COLOR,";
    for (k=X; k<=Z; k++)
      out << vertices[i].color[k] << ",";
    if (vertices[i].valid_normal) {
      out << "NORMAL,";
      for (k=RED; k<=BLUE; k++)
        out << vertices[i].normal[k] << ",";
    }
    out << "VERTEX,";
    for (k=X; k<=Z; k++)
      out << vertices[i].cpt[k] << ",";
    out << endl;
  }
  out << "END\n";
  writepymoltail(out,objname);
}

// Write out triangles as POINTS
void GLSurface::write_vspheres(ofstream &out, const string &objname,
                               double radius) {
  unsigned i,k;

  writepymolheader(out);
  for (i=0; i<vertices.size(); i++) {
    // Output coords, colors, and normals of points
    if (vertices[i].transparency==0)
      continue;
      out << "ALPHA," << vertices[i].transparency << ",";
    out << "COLOR,";
    for (k=X; k<=Z; k++)
      out << vertices[i].color[k] << ",";
    out << "SPHERE,";
    for (k=X; k<=Z; k++)
      out << vertices[i].cpt[k] << ",";
    out << radius << ",\n";
  }
  writepymoltail(out,objname);
}

// Write out triangles as POINTS
void GLSurface::writespheres(ofstream &out, const string &objname) {
  unsigned i,k;

  writepymolheader(out);
  for (i=0; i<spheres.size(); i++) {
    // Output coords, colors, and normals of points
    if (vertices[i].transparency==0)
      continue;
      out << "ALPHA," << vertices[spheres[i].i].transparency << ",";
    out << "COLOR,";
    for (k=X; k<=Z; k++)
      out << vertices[spheres[i].i].color[k] << ",";
    out << "SPHERE,";
    for (k=X; k<=Z; k++)
      out << vertices[spheres[i].i].cpt[k] << ",";
    out << spheres[i].radius << ",\n";
  }
  writepymoltail(out,objname);
}

// Write out triangles as mesh
// OPTIMIZE TO REMOVE IDENTICAL LINES
void GLSurface::writetris_mesh(ofstream &out, const string &objname) {
  unsigned i,j,k;
  int index;

  writepymolheader(out);
  for (i=0; i<triangles.size(); i++) {
    out << "BEGIN, LINE_LOOP," << endl;
    // Output coords, colors, and normals of points
    for (j=0; j<3; j++) {
      index=triangles[i].point[j];
//      if (vertices[index].transparency!=1)
        out << "ALPHA," << vertices[index].transparency << ",";
      out << "COLOR,";
      for (k=X; k<=Z; k++)
        out << vertices[index].color[k] << ",";
      if (vertices[index].valid_normal) {
        out << "NORMAL,";
        for (k=RED; k<=BLUE; k++)
          out << vertices[index].normal[k] << ",";
      }
      out << "VERTEX,";
      for (k=X; k<=Z; k++)
        out << vertices[index].cpt[k] << ",";
    }
    out << "END,";
    out << endl;
  }
  writepymoltail(out,objname);
}
  
// Write out Line Strips (Loops)
void GLSurface::writelinestrips(ofstream &out, const string &objname) {
  unsigned i,j,k, size;
  int index;

  writepymolheader(out);
  //out << "LINEWIDTH, 2," << endl;
  for (i=0; i<linestrips.size(); i++) {
    size=linestrips[i].line.size();
    if (size==1)
      out << "BEGIN, POINTS," << endl;
    else if (size==2)
      out << "BEGIN, LINES," << endl;
    else if (linestrips[i].loop)
      out << "BEGIN, LINE_LOOP," << endl;
    else
      out << "BEGIN, LINE_STRIP," << endl;
    // Output coords, colors, and normals of points
    for (j=0; j<linestrips[i].line.size(); j++) {
      index=linestrips[i].line[j];
//      if (vertices[index].transparency!=1)
        out << "ALPHA," << vertices[index].transparency << ",";
      out << "COLOR,";
      for (k=X; k<=Z; k++)
        out << vertices[index].color[k] << ",";
      if (vertices[index].valid_normal) {
        out << "NORMAL,";
        for (k=RED; k<=BLUE; k++)
          out << vertices[index].normal[k] << ",";
      }
      out << "VERTEX,";
      for (k=X; k<=Z; k++)
        out << vertices[index].cpt[k] << ",";
    }
    out << "END,";
    out << endl;
  }
  writepymoltail(out,objname);
}

// Add a line
void GLSurface::addline(unsigned v1, unsigned v2) {
  GLline line;
  line.points[0]=v1;
  line.points[1]=v2;
  gllines.push_back(line);
}

// Write out GLlines
void GLSurface::writelines(ofstream &out, const string &objname) {
  writelines(gllines, out, objname);
}

void GLSurface::writelines(vector<GLline> &l, ofstream &out,
													 const string &objname) {
  unsigned i,j,k;
  int index;

  writepymolheader(out);
  //out << "LINEWIDTH, 2," << endl;
  out << "BEGIN, LINES," << endl;
  for (i=0; i<l.size(); i++) 
    for (j=0; j<2; j++) {
      index=l[i].points[j];
//      if (vertices[index].transparency!=1)
        out << "ALPHA," << vertices[index].transparency << ",";
      out << "COLOR,";
      for (k=X; k<=Z; k++)
        out << vertices[index].color[k] << ",";
      if (vertices[i].valid_normal) {
        out << "NORMAL,";
        for (k=RED; k<=BLUE; k++)
          out << vertices[index].normal[k] << ",";
      }
      out << "VERTEX,";
      for (k=X; k<=Z; k++)
        out << vertices[index].cpt[k] << ",";
      out << endl;
    }
  out << "END,\n";
  writepymoltail(out,objname);
}

// Write a set of triangles as an xyzmesh - Requires surface was created on
//    a grid of 'resolution' with minimum latice point at 'min'
void GLSurface::writexyzmesh(ofstream &out, const string &objname) {
  writelines(xmesh,out,objname+"x");
  writelines(ymesh,out,objname+"y");
  writelines(zmesh,out,objname+"z");
  return;
}

void GLSurface::writepymolheader(ofstream &out) {
  out << "from pymol.cgo import *\n"
      << "from pymol import cmd\n"
      << "obj = [\n";
}

void GLSurface::writepymoltail(ofstream &out, const string &objname) {
  out << "]\n"
      << "cmd.load_cgo(obj,'" << objname << "')\n\n";
}

// Helper function to set the xyz and the normal for a ellipsoid vertex
void GLSurface::SQE_helper(Vertex &ver, const cPt &rad, const double u,
                           const double v, const double n,
                           const double e) {
  double sqCun=am::sign(cos(u))*powf(fabs(cos(u)),n);
  double sqSun=am::sign(sin(u))*powf(fabs(sin(u)),n);
  double sqCve=am::sign(cos(v))*powf(fabs(cos(v)),e);
  double sqSve=am::sign(sin(v))*powf(fabs(sin(v)),e);
  double sqCu2n=am::sign(cos(u))*powf(fabs(cos(u)),2.0-n);
  double sqSu2n=am::sign(sin(u))*powf(fabs(sin(u)),2.0-n);
  double sqCv2e=am::sign(cos(v))*powf(fabs(cos(v)),2.0-e);
  double sqSv2e=am::sign(sin(v))*powf(fabs(sin(v)),2.0-e);
      
  ver.cpt=rad;
  ver.cpt.x*=sqCun*sqCve;
  ver.cpt.y*=sqCun*sqSve;
  ver.cpt.z*=sqSun;
  ver.normal.x=sqCu2n*sqCv2e/rad.x;
  ver.normal.y=sqCu2n*sqSv2e/rad.y;
  ver.normal.z=sqSu2n/rad.z;
}

void GLSurface::add_super_ellipsoid(const cPt &cen, const cPt &rad,
                                    const Quaternion &q,
                                    const double n, const double e,
                                    const double u1, const double u2,
                                    const double v1, const double v2,
                                    const unsigned u_segs,
                                    const unsigned v_segs,
                                    const colorPt &color, const double alpha){
  RotMat rot(q);
  Vertex ver;
  ver.color=color;
  ver.transparency=alpha;
  ver.valid_normal=true;
  double dU=(u2-u1)/u_segs;
  double dV=(v2-v1)/v_segs;
  double U=u1;
  for (unsigned i=0; i<u_segs; i++) {
    double V=v1;
    for (unsigned j=0; j<v_segs; j++) {
      unsigned offset=vertices.size();

      SQE_helper(ver,rad,U,V,n,e);
      ver.cpt=rot*ver.cpt+cen;
      ver.normal=rot*ver.normal;
      addvertex(ver);
      SQE_helper(ver,rad,U+dU,V,n,e);
      ver.cpt=rot*ver.cpt+cen;
      ver.normal=rot*ver.normal;
      addvertex(ver);
      SQE_helper(ver,rad,U+dU,V+dV,n,e);
      ver.cpt=rot*ver.cpt+cen;
      ver.normal=rot*ver.normal;
      addvertex(ver);
      SQE_helper(ver,rad,U,V+dV,n,e);
      ver.cpt=rot*ver.cpt+cen;
      ver.normal=rot*ver.normal;
      addvertex(ver);

      Triangle t;
      t.point[0]=offset;
      t.point[1]=offset+1;
      t.point[2]=offset+3;
      addtriangle(t);
      t.point[0]=offset+1;
      t.point[1]=offset+2;
      t.point[2]=offset+3;
      addtriangle(t);

      V+=dV;
    }
    U+=dU;
  }      
}   

void GLSurface::add_ellipsoid(const cPt &center, const cPt &rad,
                              const Quaternion &rot, const colorPt &color,
                              const double alpha, const unsigned res) {
  assert(res!=0);
  add_super_ellipsoid(center,rad,rot,1.0,1.0,HALFPI*-1.0,HALFPI,PI*-1.0,PI,
                      res,res,color,alpha);
}

void GLSurface::add_ellipsoid(const cPt &center, const cPt &rad,
                              const Quaternion &rot, const colorPt &color,
                              const double alpha) {
  add_super_ellipsoid(center,rad,rot,1.0,1.0,HALFPI*-1.0,HALFPI,PI*-1.0,PI,
                      5,5,color,alpha);
}

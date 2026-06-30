// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "surf_extra.h"

#include "math_extra.h"

static constexpr double EPSILON = 1e-14;

using namespace MathExtra;

namespace SurfExtra {

/* ----------------------------------------------------------------------
   compute nearest point between sphere I and line segment J
   return 0 if no contact, 1 if pt is interior to line segment,
     -1/-2 if pt = line end point 1/2
   if contact, return:
     pt = point on line segment
     r = vector from pt to sphere center
     rsq = squared length of r
   based on geometry.cpp::distsq_point_line() in SPARTA
------------------------------------------------------------------------- */

int overlap_sphere_line(double *xsphere, double radius, double *p1, double *p2,
                        double *pt, double *r, double &rsq)
{
  double a[3],b[3];

  // A = vector from P1 to Xsphere
  // B = vector from P1 to P2

  MathExtra::sub3(xsphere,p1,a);
  MathExtra::sub3(p2,p1,b);

  // alpha = fraction of distance from P1 to P2 that P is located at
  // P = projected point on infinite line that is nearest to Xsphere center
  // alpha can be any value

  double alpha = MathExtra::dot3(a,b) / MathExtra::lensq3(b);

  // pt = point on line segment that is nearest to Xsphere center
  // if alpha <= 0.0, pt = P1, ptflag = -1
  // if alpha >= 1.0, pt = P2, ptflag = -2
  // else pt = P1 + alpha*(P2-P1), ptflag = 1

  int ptflag;
  if (alpha <= 0.0) {
    ptflag = -1;
    pt[0] = p1[0];
    pt[1] = p1[1];
    pt[2] = p1[2];
  } else if (alpha >= 1.0) {
    ptflag = -2;
    pt[0] = p2[0];
    pt[1] = p2[1];
    pt[2] = p2[2];
  } else {
    ptflag = 1;
    pt[0] = p1[0] + alpha*b[0];
    pt[1] = p1[1] + alpha*b[1];
    pt[2] = p1[2] + alpha*b[2];
  }

  // R = vector from nearest pt on line to Xsphere center
  // return ptflag if len(R) < sphere radius
  // else no contact, return 0

  double radsq = radius * radius;
  MathExtra::sub3(xsphere,pt,r);

  rsq = MathExtra::lensq3(r);
  if (rsq < radsq) return ptflag;
  return 0;
}

/* ----------------------------------------------------------------------
   compute nearest point between sphere I and triangle J
   return 0 if no contact, 1 if pt is interior to triangle,
     -1/-2/-3 if pt on tri edges, -4/-5/-6 if pt = tri corners 1/2/3
   if contact, return:
     pt = point on triangle
     r = vector from pt to sphere center
     rsq = squared length of r
   based on geometry.cpp::distsq_point_tri() in SPARTA
------------------------------------------------------------------------- */

int overlap_sphere_tri(double *xsphere, double radius,
                       double *p1, double *p2, double *p3, double *norm,
                       double *pt, double *r, double &rsq)
{
  int e12flag,e23flag,e31flag,o12flag,o23flag,o31flag;
  int esum,osum,lineflag;
  double dot, xnorm;
  double a[3],point[3],edge[3],pvec[3],xproduct[3];

  // A = vector from P1 to Xsphere
  // alpha = distance from Xsphere to triangle plane
  // pt = projection of Xsphere onto triangle plane

  MathExtra::sub3(xsphere,p1,a);
  double alpha = MathExtra::dot3(a,norm);
  pt[0] = xsphere[0] - alpha*norm[0];
  pt[1] = xsphere[1] - alpha*norm[1];
  pt[2] = xsphere[2] - alpha*norm[2];

  // test if projected point is inside triangle
  // inside = interior + boundary of tri
  // edge = edge vector of triangle
  // pvec = vector from triangle vertex to projected point
  // xproduct = cross product of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   (and pvec isn't parallel with edge),
  //   projected point is outside tri

  int inside = 1;
  e12flag = e23flag = e31flag = 0;
  o12flag = o23flag = o31flag = 0;

  MathExtra::sub3(p2,p1,edge);
  MathExtra::sub3(pt,p1,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  dot = MathExtra::dot3(xproduct,norm);
  xnorm = MathExtra::lensq3(xproduct);
  if (dot <= 0.0 && xnorm > EPSILON) {
    o12flag = 1;
    if (dot == 0.0) e12flag = 1;
    else inside = 0;
  }

  MathExtra::sub3(p3,p2,edge);
  MathExtra::sub3(pt,p2,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  dot = MathExtra::dot3(xproduct,norm);
  xnorm = MathExtra::lensq3(xproduct);
  if (dot <= 0.0 && xnorm > EPSILON) {
    o23flag = 1;
    if (dot == 0.0) e23flag = 1;
    else inside = 0;
  }

  MathExtra::sub3(p1,p3,edge);
  MathExtra::sub3(pt,p3,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  dot = MathExtra::dot3(xproduct,norm);
  xnorm = MathExtra::lensq3(xproduct);
  if (dot <= 0.0 && xnorm > EPSILON) {
    o31flag = 1;
    if (dot == 0.0) e31flag = 1;
    else inside = 0;
  }

  // projected point is inside tri = interior or boundary of tri
  // the projected point is the nearest triangle point to sphere center
  // set ptflag = 1 for interior
  // set ptflag = -1,-2,-3 for 3 edges E12,E23,E31
  // set ptflag = -4,-5,-6 for 3 corner pts P1,P2,P3

  int flag = 0;
  if (inside) {
    flag = 1;
    esum = e12flag + e23flag + e31flag;
    if (esum) {
      if (esum == 1) {
        if (e12flag) flag = -1;
        else if (e23flag) flag = -2;
        else flag = -3;
      } else {
        if (!e12flag) flag = -6;
        else if (!e23flag) flag = -4;
        else flag = -5;
      }
    }

  // projected point is outside tri
  // reset pt to be nearest triangle point to sphere center
  // if osum = 1, projected point is outside one edge
  //   nearest tri pt must be on that edge
  // if osum = 2, projected pt is outside two edges
  //   nearest tri pt could be interior to either edge or on both (corner)
  // set ptflag = -1,-2,-3 if nearest tri pt is interior to an edge
  // set ptflag = -4,-5,-6 if nearest tri pt is a corner pt

  } else {
    osum = o12flag + o23flag + o31flag;

    // projected point is outside a single edge
    // 3 possible cases:
    // nearest tri pt is interior to edge
    // or nearest tri pt = either corner pt of edge

    if (osum == 1) {
      if (o12flag) {
        lineflag = SurfExtra::nearest_point_line(xsphere,p1,p2,pt);
        if (lineflag == 1) flag = -1;
        else if (lineflag == -1) flag = -4;
        else flag = -5;
      } else if (o23flag) {
        lineflag = SurfExtra::nearest_point_line(xsphere,p2,p3,pt);
        if (lineflag == 1) flag = -2;
        else if (lineflag == -1) flag = -5;
        else flag = -6;
      } else {
        lineflag = SurfExtra::nearest_point_line(xsphere,p3,p1,pt);
        if (lineflag == 1) flag = -3;
        else if (lineflag == -1) flag = -6;
        else flag = -4;
      }

    // projected point is outside two edges
    // 3 possible cases:
    // nearest tri pt is interior to either of the two edges
    //   (only possible if angle between 2 edges is obtuse)
    // or nearest tri pt = corner point between the two edges

    } else {
      if (!o12flag) {
        lineflag = SurfExtra::nearest_point_line(xsphere,p2,p3,pt);
        if (lineflag == 1) flag = -2;
        else {
          lineflag = SurfExtra::nearest_point_line(xsphere,p3,p1,pt);
          if (lineflag == 1) flag = -3;
          else {
            flag = -6;
            pt[0] = p3[0];
            pt[1] = p3[1];
            pt[2] = p3[2];
          }
        }
      } else if (!o23flag) {
        lineflag = SurfExtra::nearest_point_line(xsphere,p3,p1,pt);
        if (lineflag == 1) flag = -3;
        else {
          lineflag = SurfExtra::nearest_point_line(xsphere,p1,p2,pt);
          if (lineflag == 1) flag = -1;
          else {
            flag = -4;
            pt[0] = p1[0];
            pt[1] = p1[1];
            pt[2] = p1[2];
          }
        }
      } else {
        lineflag = SurfExtra::nearest_point_line(xsphere,p1,p2,pt);
        if (lineflag == 1) flag = -1;
        else {
          lineflag = SurfExtra::nearest_point_line(xsphere,p2,p3,pt);
          if (lineflag == 1) flag = -2;
          else {
            flag = -5;
            pt[0] = p2[0];
            pt[1] = p2[1];
            pt[2] = p2[2];
          }
        }
      }
    }
  }

  // R = vector from nearest pt on triangle to Xsphere center
  // return flag if len(R) < sphere radius
  // else no contact, return 0

  double radsq = radius * radius;
  MathExtra::sub3(xsphere,pt,r);
  rsq = MathExtra::lensq3(r);

  if (rsq < radsq) return flag;
  return 0;
}

/* ----------------------------------------------------------------------
   compute nearest point between point X and line segment P1 to P2
   return pt = nearest point within line segment
   return 1 if pt is interior to line segment
   return -1/-2 if pt = line segment end point 1/2
   based on geometry.cpp::distsq_point_line() in SPARTA
------------------------------------------------------------------------- */

int nearest_point_line(double *x, double *p1, double *p2, double *pt)
{
  double a[3],b[3];

  // A = vector from P1 to X
  // B = vector from P1 to P2

  MathExtra::sub3(x,p1,a);
  MathExtra::sub3(p2,p1,b);

  // alpha = fraction of distance from P1 to P2 that P is located at
  // P = projected point on infinite line that is nearest to X
  // alpha can be any value

  double alpha = MathExtra::dot3(a,b) / MathExtra::lensq3(b);

  // pt = point on line segment that is nearest to X
  // if alpha <= 0.0, pt = P1, ptflag = -1
  // if alpha >= 1.0, pt = P2, ptflag = -2
  // else pt = P1 + alpha*(P2-P1), ptflag = 1

  int ptflag;
  if (alpha <= 0.0) {
    ptflag = -1;
    pt[0] = p1[0];
    pt[1] = p1[1];
    pt[2] = p1[2];
  } else if (alpha >= 1.0) {
    ptflag = -2;
    pt[0] = p2[0];
    pt[1] = p2[1];
    pt[2] = p2[2];
  } else {
    ptflag = 1;
    pt[0] = p1[0] + alpha*b[0];
    pt[1] = p1[1] + alpha*b[1];
    pt[2] = p1[2] + alpha*b[2];
  }

  return ptflag;
}

/* ---------------------------------------------------------------------- */

}

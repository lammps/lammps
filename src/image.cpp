/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Nathan Fabian (Sandia)
------------------------------------------------------------------------- */

#include "image.h"
#include <mpi.h>
#include <cmath>
#include <cctype>
#include <cstring>
#include "math_extra.h"
#include "random_mars.h"
#include "math_const.h"
#include "error.h"
#include "force.h"
#include "memory.h"

#ifdef LAMMPS_JPEG
#include <jpeglib.h>
#endif

#ifdef LAMMPS_PNG
#include <png.h>
#include <zlib.h>
#include <csetjmp>
#include "version.h"
#endif

using namespace LAMMPS_NS;
using namespace MathConst;

#define NCOLORS 140
#define NELEMENTS 109
#define EPSILON 1.0e-6

enum{NUMERIC,MINVALUE,MAXVALUE};
enum{CONTINUOUS,DISCRETE,SEQUENTIAL};
enum{ABSOLUTE,FRACTIONAL};
enum{NO,YES};

/* ---------------------------------------------------------------------- */

Image::Image(LAMMPS *lmp, int nmap_caller) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // defaults for 3d viz

  width = height = 512;
  theta = 60.0 * MY_PI/180.0;
  phi = 30.0 * MY_PI/180.0;
  zoom = 1.0;
  persp = 0.0;
  shiny = 1.0;
  ssao = NO;

  up[0] = 0.0;
  up[1] = 0.0;
  up[2] = 1.0;

  // colors

  ncolors = 0;
  username = NULL;
  userrgb = NULL;

  boxcolor = color2rgb("yellow");
  background[0] = background[1] = background[2] = 0;

  // define nmap colormaps, all with default settings

  nmap = nmap_caller;
  maps = new ColorMap*[nmap];
  for (int i = 0; i < nmap; i++)
    maps[i] = new ColorMap(lmp,this);

  // static parameters

  FOV = MY_PI/6.0;              // 30 degrees
  ambientColor[0] = 0.0;
  ambientColor[1] = 0.0;
  ambientColor[2] = 0.0;

  keyLightPhi = -MY_PI4;        // -45 degrees
  keyLightTheta = MY_PI/6.0;    // 30 degrees
  keyLightColor[0] = 0.9;
  keyLightColor[1] = 0.9;
  keyLightColor[2] = 0.9;

  fillLightPhi = MY_PI/6.0;     // 30 degrees
  fillLightTheta = 0;
  fillLightColor[0] = 0.45;
  fillLightColor[1] = 0.45;
  fillLightColor[2] = 0.45;

  backLightPhi = MY_PI;         // 180 degrees
  backLightTheta = MY_PI/12.0;  // 15 degrees
  backLightColor[0] = 0.9;
  backLightColor[1] = 0.9;
  backLightColor[2] = 0.9;

  random = NULL;
}

/* ---------------------------------------------------------------------- */

Image::~Image()
{
  for (int i = 0; i < nmap; i++) delete maps[i];
  delete [] maps;

  for (int i = 0; i < ncolors; i++) delete [] username[i];
  memory->sfree(username);
  memory->destroy(userrgb);

  memory->destroy(depthBuffer);
  memory->destroy(surfaceBuffer);
  memory->destroy(imageBuffer);
  memory->destroy(depthcopy);
  memory->destroy(surfacecopy);
  memory->destroy(rgbcopy);

  if (random) delete random;
}

/* ----------------------------------------------------------------------
   allocate image and depth buffers
   called after image size is set
------------------------------------------------------------------------- */

void Image::buffers()
{
  npixels = width * height;
  memory->create(depthBuffer,npixels,"image:depthBuffer");
  memory->create(surfaceBuffer,2*npixels,"image:surfaceBuffer");
  memory->create(imageBuffer,3*npixels,"image:imageBuffer");
  memory->create(depthcopy,npixels,"image:depthcopy");
  memory->create(surfacecopy,2*npixels,"image:surfacecopy");
  memory->create(rgbcopy,3*npixels,"image:rgbcopy");
}

/* ----------------------------------------------------------------------
   reset view parameters
   called once if view is STATIC
   called before every render if view is DYNAMIC
------------------------------------------------------------------------- */

void Image::view_params(double boxxlo, double boxxhi, double boxylo,
                        double boxyhi, double boxzlo, double boxzhi)
{
  // camDir points at the camera, view direction = -camDir

  camDir[0] = sin(theta)*cos(phi);
  camDir[1] = sin(theta)*sin(phi);
  camDir[2] = cos(theta);

  // normalize up vector

  if (up[0] == 0.0 && up[1] == 0.0 && up[2] == 0.0)
    error->all(FLERR,"Invalid image up vector");
  MathExtra::norm3(up);

  // adjust camDir by epsilon if camDir and up are parallel
  // do this by tweaking view direction, not up direction
  // try to insure continuous images as changing view passes thru up
  // sufficient to handle common cases where theta = 0 or 180 is degenerate?

  double dot = MathExtra::dot3(up,camDir);
  if (fabs(dot) > 1.0-EPSILON) {
    if (theta == 0.0) {
      camDir[0] = sin(EPSILON)*cos(phi);
      camDir[1] = sin(EPSILON)*sin(phi);
      camDir[2] = cos(EPSILON);
    } else if (theta == MY_PI) {
      camDir[0] = sin(theta-EPSILON)*cos(phi);
      camDir[1] = sin(theta-EPSILON)*sin(phi);
      camDir[2] = cos(theta-EPSILON);
    } else {
      camDir[0] = sin(theta+EPSILON)*cos(phi);
      camDir[1] = sin(theta+EPSILON)*sin(phi);
      camDir[2] = cos(theta+EPSILON);
    }
  }

  // camUp = camDir x (Up x camDir)

  MathExtra::cross3(up,camDir,camRight);
  MathExtra::norm3(camRight);
  MathExtra::cross3(camDir,camRight,camUp);
  if (camUp[0] == 0.0 && camUp[1] == 0.0 && camUp[2] == 0.0)
    error->all(FLERR,"Invalid image up vector");
  MathExtra::norm3(camUp);

  // zdist = camera distance = function of zoom & bounding box
  // camPos = camera position = function of camDir and zdist

  double delx = 2.0*(boxxhi-boxxlo);
  double dely = 2.0*(boxyhi-boxylo);
  double delz = 2.0*(boxzhi-boxzlo);
  double maxdel = MAX(delx,dely);
  maxdel = MAX(maxdel,delz);

  zdist = maxdel;
  zdist /= tan(FOV);
  zdist += 0.5 * (delx*camDir[0] + dely*camDir[1] + delz*camDir[2]);
  zdist /= zoom;

  camPos[0] = camDir[0] * zdist;
  camPos[1] = camDir[1] * zdist;
  camPos[2] = camDir[2] * zdist;

  // light directions in terms of -camDir = z

  keyLightDir[0] = cos(keyLightTheta) * sin(keyLightPhi);
  keyLightDir[1] = sin(keyLightTheta);
  keyLightDir[2] = cos(keyLightTheta) * cos(keyLightPhi);

  fillLightDir[0] = cos(fillLightTheta) * sin(fillLightPhi);
  fillLightDir[1] = sin(fillLightTheta);
  fillLightDir[2] = cos(fillLightTheta) * cos(fillLightPhi);

  backLightDir[0] = cos(backLightTheta) * sin(backLightPhi);
  backLightDir[1] = sin(backLightTheta);
  backLightDir[2] = cos(backLightTheta) * cos(backLightPhi);

  keyHalfDir[0] = 0 + keyLightDir[0];
  keyHalfDir[1] = 0 + keyLightDir[1];
  keyHalfDir[2] = 1 + keyLightDir[2];
  MathExtra::norm3(keyHalfDir);

  // adjust shinyness of the reflection

  specularHardness = 16.0 * shiny;
  specularIntensity = shiny;

  // adjust strength of the SSAO

  if (ssao) {
    if (!random) random = new RanMars(lmp,seed+me);
    SSAORadius = maxdel * 0.05 * ssaoint;
    SSAOSamples = static_cast<int> (8.0 + 32.0*ssaoint);
    SSAOJitter = MY_PI / 12;
    ambientColor[0] = 0.5;
    ambientColor[1] = 0.5;
    ambientColor[2] = 0.5;
  }

  // param for rasterizing spheres

  tanPerPixel = -(maxdel / (double) height);
}

/* ----------------------------------------------------------------------
   initialize image to background color and depth buffer
   no need to init surfaceBuffer, since will be based on depth
------------------------------------------------------------------------- */

void Image::clear()
{
  int red = background[0];
  int green = background[1];
  int blue = background[2];

  int ix,iy;
  for (iy = 0; iy < height; iy ++)
    for (ix = 0; ix < width; ix ++) {
      imageBuffer[iy * width * 3 + ix * 3 + 0] = red;
      imageBuffer[iy * width * 3 + ix * 3 + 1] = green;
      imageBuffer[iy * width * 3 + ix * 3 + 2] = blue;
      depthBuffer[iy * width + ix] = -1;
    }
}

/* ----------------------------------------------------------------------
   merge image from each processor into one composite image
   done pixel by pixel, respecting depth buffer
   hi procs send to lo procs, cascading down logarithmically
------------------------------------------------------------------------- */

void Image::merge()
{
  MPI_Request requests[3];

  int nhalf = 1;
  while (nhalf < nprocs) nhalf *= 2;
  nhalf /= 2;

  while (nhalf) {
    if (me < nhalf && me+nhalf < nprocs) {
      MPI_Irecv(rgbcopy,npixels*3,MPI_BYTE,me+nhalf,0,world,&requests[0]);
      MPI_Irecv(depthcopy,npixels,MPI_DOUBLE,me+nhalf,0,world,&requests[1]);
      if (ssao)
        MPI_Irecv(surfacecopy,npixels*2,MPI_DOUBLE,
                  me+nhalf,0,world,&requests[2]);
      if (ssao) MPI_Waitall(3,requests,MPI_STATUS_IGNORE);
      else MPI_Waitall(2,requests,MPI_STATUS_IGNORE);

      for (int i = 0; i < npixels; i++) {
        if (depthBuffer[i] < 0 || (depthcopy[i] >= 0 &&
                                   depthcopy[i] < depthBuffer[i])) {
          depthBuffer[i] = depthcopy[i];
          imageBuffer[i*3+0] = rgbcopy[i*3+0];
          imageBuffer[i*3+1] = rgbcopy[i*3+1];
          imageBuffer[i*3+2] = rgbcopy[i*3+2];
          if (ssao) {
            surfaceBuffer[i*2+0] = surfacecopy[i*2+0];
            surfaceBuffer[i*2+1] = surfacecopy[i*2+1];
          }
        }
      }

    } else if (me >= nhalf && me < 2*nhalf) {
      MPI_Send(imageBuffer,npixels*3,MPI_BYTE,me-nhalf,0,world);
      MPI_Send(depthBuffer,npixels,MPI_DOUBLE,me-nhalf,0,world);
      if (ssao) MPI_Send(surfaceBuffer,npixels*2,MPI_DOUBLE,me-nhalf,0,world);
    }

    nhalf /= 2;
  }

  // extra SSAO enhancement
  // bcast full image to all procs
  // each works on subset of pixels
  // gather result back to proc 0

  if (ssao) {
    MPI_Bcast(imageBuffer,npixels*3,MPI_BYTE,0,world);
    MPI_Bcast(surfaceBuffer,npixels*2,MPI_DOUBLE,0,world);
    MPI_Bcast(depthBuffer,npixels,MPI_DOUBLE,0,world);
    compute_SSAO();
    int pixelPart = height/nprocs * width*3;
    MPI_Gather(imageBuffer+me*pixelPart,pixelPart,MPI_BYTE,
               rgbcopy,pixelPart,MPI_BYTE,0,world);
    writeBuffer = rgbcopy;
  } else {
    writeBuffer = imageBuffer;
  }
}

/* ----------------------------------------------------------------------
   draw simulation bounding box as 12 cylinders
------------------------------------------------------------------------- */

void Image::draw_box(double (*corners)[3], double diameter)
{
  draw_cylinder(corners[0],corners[1],boxcolor,diameter,3);
  draw_cylinder(corners[2],corners[3],boxcolor,diameter,3);
  draw_cylinder(corners[0],corners[2],boxcolor,diameter,3);
  draw_cylinder(corners[1],corners[3],boxcolor,diameter,3);
  draw_cylinder(corners[0],corners[4],boxcolor,diameter,3);
  draw_cylinder(corners[1],corners[5],boxcolor,diameter,3);
  draw_cylinder(corners[2],corners[6],boxcolor,diameter,3);
  draw_cylinder(corners[3],corners[7],boxcolor,diameter,3);
  draw_cylinder(corners[4],corners[5],boxcolor,diameter,3);
  draw_cylinder(corners[6],corners[7],boxcolor,diameter,3);
  draw_cylinder(corners[4],corners[6],boxcolor,diameter,3);
  draw_cylinder(corners[5],corners[7],boxcolor,diameter,3);
}

/* ----------------------------------------------------------------------
   draw XYZ axes in red/green/blue
   axes = 4 end points
------------------------------------------------------------------------- */

void Image::draw_axes(double (*axes)[3], double diameter)
{
  draw_cylinder(axes[0],axes[1],color2rgb("red"),diameter,3);
  draw_cylinder(axes[0],axes[2],color2rgb("green"),diameter,3);
  draw_cylinder(axes[0],axes[3],color2rgb("blue"),diameter,3);
}

/* ----------------------------------------------------------------------
   draw sphere at x with surfaceColor and diameter
   render pixel by pixel onto image plane with depth buffering
------------------------------------------------------------------------- */

void Image::draw_sphere(double *x, double *surfaceColor, double diameter)
{
  int ix,iy;
  double projRad;
  double xlocal[3],surface[3];
  double depth;

  xlocal[0] = x[0] - xctr;
  xlocal[1] = x[1] - yctr;
  xlocal[2] = x[2] - zctr;

  double xmap = MathExtra::dot3(camRight,xlocal);
  double ymap = MathExtra::dot3(camUp,xlocal);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(xlocal,camDir);

  double radius = 0.5*diameter;
  double radsq = radius*radius;
  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist :
    -tanPerPixel / zoom;
  double pixelRadiusFull = radius / pixelWidth;
  int pixelRadius = static_cast<int> (pixelRadiusFull + 0.5) + 1;

  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int> (xf);
  int yc = static_cast<int> (yf);
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  for (iy = yc - pixelRadius; iy <= yc + pixelRadius; iy++) {
    for (ix = xc - pixelRadius; ix <= xc + pixelRadius; ix++) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;

      surface[1] = ((iy - yc) - height_error) * pixelWidth;
      surface[0] = ((ix - xc) - width_error) * pixelWidth;
      projRad = surface[0]*surface[0] + surface[1]*surface[1];

      // outside the sphere in the projected image

      if (projRad > radsq) continue;
      surface[2] = sqrt(radsq - projRad);
      depth = dist - surface[2];

      surface[0] /= radius;
      surface[1] /= radius;
      surface[2] /= radius;

      draw_pixel (ix, iy, depth, surface, surfaceColor);
    }
  }
}

/* ----------------------------------------------------------------------
   draw axis oriented cube at x with surfaceColor and diameter in size
   render pixel by pixel onto image plane with depth buffering
------------------------------------------------------------------------- */

void Image::draw_cube(double *x, double *surfaceColor, double diameter)
{
  double xlocal[3],surface[3],normal[3];
  double t,tdir[3];
  double depth;

  xlocal[0] = x[0] - xctr;
  xlocal[1] = x[1] - yctr;
  xlocal[2] = x[2] - zctr;

  double xmap = MathExtra::dot3(camRight,xlocal);
  double ymap = MathExtra::dot3(camUp,xlocal);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(xlocal,camDir);

  double radius = 0.5*diameter;
  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist :
    -tanPerPixel / zoom;

  double halfWidth = diameter;
  double pixelHalfWidthFull = halfWidth / pixelWidth;
  int pixelHalfWidth = static_cast<int> (pixelHalfWidthFull + 0.5);

  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int> (xf);
  int yc = static_cast<int> (yf);
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  for (int iy = yc - pixelHalfWidth; iy <= yc + pixelHalfWidth; iy ++) {
    for (int ix = xc - pixelHalfWidth; ix <= xc + pixelHalfWidth; ix ++) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;

      double sy = ((iy - yc) - height_error) * pixelWidth;
      double sx = ((ix - xc) - width_error) * pixelWidth;
      surface[0] = camRight[0] * sx + camUp[0] * sy;
      surface[1] = camRight[1] * sx + camUp[1] * sy;
      surface[2] = camRight[2] * sx + camUp[2] * sy;

      // iterate through each of the 6 axis-oriented planes of the box
      // only render up to 3 which are facing the camera
      // these checks short circuit a dot product, testing for > 0

      for (int dim = 0; dim < 3; dim ++) {
        if (camDir[dim] > 0) {          // positive faces camera
          t = (radius - surface[dim]) / camDir[dim];
          normal[0] = camRight[dim];
          normal[1] = camUp[dim];
          normal[2] = camDir[dim];
        } else if (camDir[dim] < 0) {   // negative faces camera
          t = -(radius + surface[dim]) / camDir[dim];
          normal[0] = -camRight[dim];
          normal[1] = -camUp[dim];
          normal[2] = -camDir[dim];
        }
        if (camDir[dim] != 0) {
          tdir[0] = camDir[0] * t;
          tdir[1] = camDir[1] * t;
          tdir[2] = camDir[2] * t;

          bool xin = ((surface[0]+tdir[0]) >= -radius) &&
            ((surface[0]+tdir[0]) <= radius);
          bool yin = ((surface[1]+tdir[1]) >= -radius) &&
            ((surface[1]+tdir[1]) <= radius);
          bool zin = ((surface[2]+tdir[2]) >= -radius) &&
            ((surface[2]+tdir[2]) <= radius);

          switch (dim) {
          case 0:
            if (yin & zin) {
              depth = dist - t;
              draw_pixel (ix, iy, depth, normal, surfaceColor);
            }
            break;
          case 1:
            if (xin & zin) {
              depth = dist - t;
              draw_pixel (ix, iy, depth, normal, surfaceColor);
            }
            break;
          case 2:
            if (xin & yin) {
              depth = dist - t;
              draw_pixel (ix, iy, depth, normal, surfaceColor);
            }
            break;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   draw cylinder from x to y with surfaceColor and diameter
   render pixel by pixel onto image plane with depth buffering
   if sflag = 0, draw no end spheres
   if sflag = 1, draw 1st end sphere
   if sflag = 2, draw 2nd end sphere
   if sflag = 3, draw both end spheres
------------------------------------------------------------------------- */

void Image::draw_cylinder(double *x, double *y,
                          double *surfaceColor, double diameter, int sflag)
{
  double surface[3], normal[3];
  double mid[3],xaxis[3],yaxis[3],zaxis[3];
  double camLDir[3], camLRight[3], camLUp[3];
  double zmin, zmax;

  if (sflag % 2) draw_sphere(x,surfaceColor,diameter);
  if (sflag/2) draw_sphere(y,surfaceColor,diameter);

  double radius = 0.5*diameter;
  double radsq = radius*radius;

  zaxis[0] = y[0] - x[0];
  zaxis[1] = y[1] - x[1];
  zaxis[2] = y[2] - x[2];

  double rasterWidth = fabs(MathExtra::dot3(zaxis, camRight)) + diameter;
  double rasterHeight = fabs(MathExtra::dot3(zaxis, camUp)) + diameter;

  mid[0] = (y[0] + x[0]) * 0.5 - xctr;
  mid[1] = (y[1] + x[1]) * 0.5 - yctr;
  mid[2] = (y[2] + x[2]) * 0.5 - zctr;

  double len = MathExtra::len3(zaxis);
  MathExtra::scale3(1.0/len,zaxis);
  len *= 0.5;
  zmax = len;
  zmin = -len;

  double xmap = MathExtra::dot3(camRight,mid);
  double ymap = MathExtra::dot3(camUp,mid);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(mid,camDir);

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist :
    -tanPerPixel / zoom;

  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int> (xf);
  int yc = static_cast<int> (yf);
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  double pixelHalfWidthFull = (rasterWidth * 0.5) / pixelWidth;
  double pixelHalfHeightFull = (rasterHeight * 0.5) / pixelWidth;
  int pixelHalfWidth = static_cast<int> (pixelHalfWidthFull + 0.5);
  int pixelHalfHeight = static_cast<int> (pixelHalfHeightFull + 0.5);

  if (zaxis[0] == camDir[0] && zaxis[1] == camDir[1] && zaxis[2] == camDir[2])
    return;
  if (zaxis[0] == -camDir[0] && zaxis[1] == -camDir[1] &&
      zaxis[2] == -camDir[2]) return;

  MathExtra::cross3(zaxis,camDir,yaxis);
  MathExtra::norm3(yaxis);
  MathExtra::cross3(yaxis,zaxis,xaxis);
  MathExtra::norm3(xaxis);

  camLDir[0] = MathExtra::dot3(camDir,xaxis);
  camLDir[1] = 0.0;
  camLDir[2] = MathExtra::dot3(camDir,zaxis);

  camLRight[0] = MathExtra::dot3(camRight,xaxis);
  camLRight[1] = MathExtra::dot3(camRight,yaxis);
  camLRight[2] = MathExtra::dot3(camRight,zaxis);
  MathExtra::norm3(camLRight);

  camLUp[0] = MathExtra::dot3(camUp,xaxis);
  camLUp[1] = MathExtra::dot3(camUp,yaxis);
  camLUp[2] = MathExtra::dot3(camUp,zaxis);
  MathExtra::norm3(camLUp);

  double a = camLDir[0] * camLDir[0];

  for (int iy = yc - pixelHalfHeight; iy <= yc + pixelHalfHeight; iy ++) {
    for (int ix = xc - pixelHalfWidth; ix <= xc + pixelHalfWidth; ix ++) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;

      double sy = ((iy - yc) - height_error) * pixelWidth;
      double sx = ((ix - xc) - width_error) * pixelWidth;
      surface[0] = camLRight[0] * sx + camLUp[0] * sy;
      surface[1] = camLRight[1] * sx + camLUp[1] * sy;
      surface[2] = camLRight[2] * sx + camLUp[2] * sy;

      double b = 2 * camLDir[0] * surface[0];
      double c = surface[0] * surface[0] + surface[1] * surface[1] - radsq;

      double partial = b*b - 4*a*c;
      if (partial < 0) continue;
      partial = sqrt (partial);

      double t = (-b + partial) / (2*a);
      double t2 = (-b - partial) / (2*a);
      if (t2 > t) { t = t2; }

      surface[0] += t * camLDir[0];
      surface[1] += t * camLDir[1];
      surface[2] += t * camLDir[2];

      if (surface[2] > zmax || surface[2] < zmin) continue;

      // convert surface into the surface normal

      normal[0] = surface[0] / radius;
      normal[1] = surface[1] / radius;
      normal[2] = 0.0;

      // in camera space

      surface[0] = MathExtra::dot3 (normal, camLRight);
      surface[1] = MathExtra::dot3 (normal, camLUp);
      surface[2] = MathExtra::dot3 (normal, camLDir);

      double depth = dist - t;
      draw_pixel (ix, iy, depth, surface, surfaceColor);
    }
  }
}

/* ----------------------------------------------------------------------
   draw triangle with 3 corner points x,y,z and surfaceColor
------------------------------------------------------------------------- */

void Image::draw_triangle(double *x, double *y, double *z, double *surfaceColor)
{
  double d1[3], d1len, d2[3], d2len, normal[3], invndotd;
  double xlocal[3], ylocal[3], zlocal[3];
  double surface[3];
  double depth;

  xlocal[0] = x[0] - xctr;
  xlocal[1] = x[1] - yctr;
  xlocal[2] = x[2] - zctr;
  ylocal[0] = y[0] - xctr;
  ylocal[1] = y[1] - yctr;
  ylocal[2] = y[2] - zctr;
  zlocal[0] = z[0] - xctr;
  zlocal[1] = z[1] - yctr;
  zlocal[2] = z[2] - zctr;

  MathExtra::sub3 (xlocal, ylocal, d1);
  d1len = MathExtra::len3 (d1);
  MathExtra::scale3 (1.0 / d1len, d1);
  MathExtra::sub3 (zlocal, ylocal, d2);
  d2len = MathExtra::len3 (d2);
  MathExtra::scale3 (1.0 / d2len, d2);

  MathExtra::cross3 (d1, d2, normal);
  MathExtra::norm3 (normal);
  invndotd = 1.0 / MathExtra::dot3(normal, camDir);

  // invalid triangle (parallel)

  if (invndotd == 0) return;

  double r[3],u[3];

  r[0] = MathExtra::dot3(camRight,xlocal);
  r[1] = MathExtra::dot3(camRight,ylocal);
  r[2] = MathExtra::dot3(camRight,zlocal);

  u[0] = MathExtra::dot3(camUp,xlocal);
  u[1] = MathExtra::dot3(camUp,ylocal);
  u[2] = MathExtra::dot3(camUp,zlocal);

  double rasterLeft = r[0] - MIN(r[0],MIN(r[1],r[2]));
  double rasterRight = MAX(r[0],MAX(r[1],r[2])) - r[0];
  double rasterDown = u[0] - MIN(u[0],MIN(u[1],u[2]));
  double rasterUp = MAX(u[0],MAX(u[1],u[2])) - u[0];

  double xmap = MathExtra::dot3(camRight,xlocal);
  double ymap = MathExtra::dot3(camUp,xlocal);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(xlocal,camDir);

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist :
    -tanPerPixel / zoom;

  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int> (xf);
  int yc = static_cast<int> (yf);
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  double pixelLeftFull = rasterLeft / pixelWidth;
  double pixelRightFull = rasterRight / pixelWidth;
  double pixelDownFull = rasterDown / pixelWidth;
  double pixelUpFull = rasterUp / pixelWidth;
  int pixelLeft = static_cast<int> (pixelLeftFull + 0.5);
  int pixelRight = static_cast<int> (pixelRightFull + 0.5);
  int pixelDown = static_cast<int> (pixelDownFull + 0.5);
  int pixelUp = static_cast<int> (pixelUpFull + 0.5);

  for (int iy = yc - pixelDown; iy <= yc + pixelUp; iy ++) {
    for (int ix = xc - pixelLeft; ix <= xc + pixelRight; ix ++) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;

      double sy = ((iy - yc) - height_error) * pixelWidth;
      double sx = ((ix - xc) - width_error) * pixelWidth;
      surface[0] = camRight[0] * sx + camUp[0] * sy;
      surface[1] = camRight[1] * sx + camUp[1] * sy;
      surface[2] = camRight[2] * sx + camUp[2] * sy;

      double t = -MathExtra::dot3(normal,surface) * invndotd;

      // test inside the triangle

      double p[3];
      p[0] = xlocal[0] + surface[0] + camDir[0] * t;
      p[1] = xlocal[1] + surface[1] + camDir[1] * t;
      p[2] = xlocal[2] + surface[2] + camDir[2] * t;

      double s1[3], s2[3], s3[3];
      double c1[3], c2[3];

      MathExtra::sub3 (zlocal, xlocal, s1);
      MathExtra::sub3 (ylocal, xlocal, s2);
      MathExtra::sub3 (p, xlocal, s3);
      MathExtra::cross3 (s1, s2, c1);
      MathExtra::cross3 (s1, s3, c2);
      if (MathExtra::dot3 (c1, c2) <= 0) continue;

      MathExtra::sub3 (xlocal, ylocal, s1);
      MathExtra::sub3 (zlocal, ylocal, s2);
      MathExtra::sub3 (p, ylocal, s3);
      MathExtra::cross3 (s1, s2, c1);
      MathExtra::cross3 (s1, s3, c2);
      if (MathExtra::dot3 (c1, c2) <= 0) continue;

      MathExtra::sub3 (ylocal, zlocal, s1);
      MathExtra::sub3 (xlocal, zlocal, s2);
      MathExtra::sub3 (p, zlocal, s3);
      MathExtra::cross3 (s1, s2, c1);
      MathExtra::cross3 (s1, s3, c2);
      if (MathExtra::dot3 (c1, c2) <= 0) continue;

      double cNormal[3];
      cNormal[0] = MathExtra::dot3(camRight, normal);
      cNormal[1] = MathExtra::dot3(camUp, normal);
      cNormal[2] = MathExtra::dot3(camDir, normal);

      depth = dist - t;
      draw_pixel(ix,iy,depth,cNormal,surfaceColor);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Image::draw_pixel(int ix, int iy, double depth,
                           double *surface, double *surfaceColor)
{
  double diffuseKey,diffuseFill,diffuseBack,specularKey;
  if (depth < 0 || (depthBuffer[ix + iy*width] >= 0 &&
                    depth >= depthBuffer[ix + iy*width])) return;
  depthBuffer[ix + iy*width] = depth;

  // store only the tangent relative to the camera normal (0,0,-1)

  surfaceBuffer[0 + ix * 2 + iy*width * 2] = surface[1];
  surfaceBuffer[1 + ix * 2 + iy*width * 2] = -surface[0];

  diffuseKey = saturate(MathExtra::dot3(surface, keyLightDir));
  diffuseFill = saturate(MathExtra::dot3(surface, fillLightDir));
  diffuseBack = saturate(MathExtra::dot3(surface, backLightDir));
  specularKey = pow(saturate(MathExtra::dot3(surface, keyHalfDir)),
                    specularHardness) * specularIntensity;

  double c[3];
  c[0] = surfaceColor[0] * ambientColor[0];
  c[1] = surfaceColor[1] * ambientColor[1];
  c[2] = surfaceColor[2] * ambientColor[2];

  c[0] += surfaceColor[0] * keyLightColor[0] * diffuseKey;
  c[1] += surfaceColor[1] * keyLightColor[1] * diffuseKey;
  c[2] += surfaceColor[2] * keyLightColor[2] * diffuseKey;

  c[0] += keyLightColor[0] * specularKey;
  c[1] += keyLightColor[1] * specularKey;
  c[2] += keyLightColor[2] * specularKey;

  c[0] += surfaceColor[0] * fillLightColor[0] * diffuseFill;
  c[1] += surfaceColor[1] * fillLightColor[1] * diffuseFill;
  c[2] += surfaceColor[2] * fillLightColor[2] * diffuseFill;

  c[0] += surfaceColor[0] * backLightColor[0] * diffuseBack;
  c[1] += surfaceColor[1] * backLightColor[1] * diffuseBack;
  c[2] += surfaceColor[2] * backLightColor[2] * diffuseBack;

  c[0] = saturate(c[0]);
  c[1] = saturate(c[1]);
  c[2] = saturate(c[2]);

  imageBuffer[0 + ix*3 + iy*width*3] = static_cast<int>(c[0] * 255.0);
  imageBuffer[1 + ix*3 + iy*width*3] = static_cast<int>(c[1] * 255.0);
  imageBuffer[2 + ix*3 + iy*width*3] = static_cast<int>(c[2] * 255.0);
}

/* ---------------------------------------------------------------------- */

void Image::compute_SSAO()
{
  // used for rasterizing the spheres

  double delTheta = 2.0*MY_PI / SSAOSamples;

  // typical neighborhood value for shading

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel :
        -tanPerPixel / zoom;
  int pixelRadius = (int) trunc (SSAORadius / pixelWidth + 0.5);

  int x,y,s;
  int hPart = height / nprocs;
  int index = me * hPart * width;
  for (y = me * hPart; y < (me + 1) * hPart; y ++) {
    for (x = 0; x < width; x ++, index ++) {
      double cdepth = depthBuffer[index];
      if (cdepth < 0) { continue; }

      double sx = surfaceBuffer[index * 2 + 0];
      double sy = surfaceBuffer[index * 2 + 1];
      double sin_t = -sqrt(sx*sx + sy*sy);

      double mytheta = random->uniform() * SSAOJitter;
      double ao = 0.0;

      for (s = 0; s < SSAOSamples; s ++) {
        double hx = cos(mytheta);
        double hy = sin(mytheta);
        mytheta += delTheta;

        // multiply by z cross surface tangent
        // so that dot (aka cos) works here

        double scaled_sin_t = sin_t * (hx*sy + hy*sx);

        // Bresenham's line algorithm to march over depthBuffer

        int dx = static_cast<int> (hx * pixelRadius);
        int dy = static_cast<int> (hy * pixelRadius);
        int ex = x + dx;
        if (ex < 0) { ex = 0; } if (ex >= width) { ex = width - 1; }
        int ey = y + dy;
        if (ey < 0) { ey = 0; } if (ey >= height) { ey = height - 1; }
        double delta;
        int small, large;
        double lenIncr;
        if (fabs(hx) > fabs(hy)) {
          small = (hx > 0) ? 1 : -1;
          large = (hy > 0) ? width : -width;
          delta = fabs(hy / hx);
        } else {
          small = (hy > 0) ? width : -width;
          large = (hx > 0) ? 1 : -1;
          delta = fabs(hx / hy);
        }
        lenIncr = sqrt (1 + delta * delta) * pixelWidth;

        // initialize with one step
        // because the center point doesn't need testing

        int end = ex + ey * width;
        int ind = index + small;
        double len = lenIncr;
        double err = delta;
        if (err >= 1.0) {
          ind += large;
          err -= 1.0;
        }

        double minPeak = -1;
        double peakLen = 0.0;
        int stepsTaken = 1;
        while ((small > 0 && ind <= end) || (small < 0 && ind >= end)) {
          if (ind < 0 || ind >= (width*height)) {
            break;
          }

          // cdepth - depthBuffer B/C we want it in the negative z direction

          if (minPeak < 0 || (depthBuffer[ind] >= 0 &&
                              depthBuffer[ind] < minPeak)) {
            minPeak = depthBuffer[ind];
            peakLen = len;
          }
          ind += small;
          len += lenIncr;
          err += delta;
          if (err >= 1.0) {
            ind += large;
            err -= 1.0;
          }
          stepsTaken ++;
        }

        if (peakLen > 0) {
          double h = atan ((cdepth - minPeak) / peakLen);
          ao += saturate(sin (h) - scaled_sin_t);
        } else {
          ao += saturate(-scaled_sin_t);
        }
      }
      ao /= (double)SSAOSamples;

      double c[3];
      c[0] = (double) (*(unsigned char *) &imageBuffer[index * 3 + 0]);
      c[1] = (double) (*(unsigned char *) &imageBuffer[index * 3 + 1]);
      c[2] = (double) (*(unsigned char *) &imageBuffer[index * 3 + 2]);
      c[0] *= (1.0 - ao);
      c[1] *= (1.0 - ao);
      c[2] *= (1.0 - ao);
      imageBuffer[index * 3 + 0] = (int) c[0];
      imageBuffer[index * 3 + 1] = (int) c[1];
      imageBuffer[index * 3 + 2] = (int) c[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void Image::write_JPG(FILE *fp)
{
#ifdef LAMMPS_JPEG
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo,fp);
  cinfo.image_width = width;
  cinfo.image_height = height;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo,85,TRUE);
  jpeg_start_compress(&cinfo,TRUE);

  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer = (JSAMPROW)
      &writeBuffer[(cinfo.image_height - 1 - cinfo.next_scanline) * 3 * width];
    jpeg_write_scanlines(&cinfo,&row_pointer,1);
  }

  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
#else
  LMP_UNUSED_PARAM(fp);
#endif
}

/* ---------------------------------------------------------------------- */

void Image::write_PNG(FILE *fp)
{
#ifdef LAMMPS_PNG
  png_structp png_ptr;
  png_infop info_ptr;

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) return;

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr, NULL);
    return;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    return;
  }

  png_init_io(png_ptr, fp);
  png_set_compression_level(png_ptr,Z_BEST_COMPRESSION);
  png_set_IHDR(png_ptr,info_ptr,width,height,8,PNG_COLOR_TYPE_RGB,
    PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);

  png_text text_ptr[2];
  memset(text_ptr,0,2*sizeof(png_text));

  char key0[]  = "Software";
  char text0[] = "LAMMPS " LAMMPS_VERSION;
  char key1[]  = "Description";
  char text1[] = "Dump image snapshot";
  text_ptr[0].key = key0;
  text_ptr[0].text = text0;
  text_ptr[1].key = key1;
  text_ptr[1].text = text1;
  text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
  text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;

  png_set_text(png_ptr,info_ptr,text_ptr,1);
  png_write_info(png_ptr,info_ptr);

  png_bytep *row_pointers = new png_bytep[height];
  for (int i=0; i < height; ++i)
    row_pointers[i] = (png_bytep) &writeBuffer[(height-i-1)*3*width];

  png_write_image(png_ptr, row_pointers);
  png_write_end(png_ptr, info_ptr);

  png_destroy_write_struct(&png_ptr, &info_ptr);
  delete[] row_pointers;
#else
  LMP_UNUSED_PARAM(fp);
#endif
}

/* ---------------------------------------------------------------------- */

void Image::write_PPM(FILE *fp)
{
  fprintf(fp,"P6\n%d %d\n255\n",width,height);

  int y;
  for (y = height-1; y >= 0; y--)
    fwrite(&writeBuffer[y*width*3],3,width,fp);
}

/* ----------------------------------------------------------------------
   return static/dynamic status of color map index
------------------------------------------------------------------------- */

int Image::map_dynamic(int index)
{
  return maps[index]->dynamic;
}

/* ----------------------------------------------------------------------
   redefine properties of the color map index
   return 1 if any error in args, else return 0
------------------------------------------------------------------------- */

int Image::map_reset(int index, int narg, char **arg)
{
  return maps[index]->reset(narg,arg);
}

/* ----------------------------------------------------------------------
   set min/max bounds of dynamic color map index
------------------------------------------------------------------------- */

int Image::map_minmax(int index, double mindynamic, double maxdynamic)
{
  return maps[index]->minmax(mindynamic,maxdynamic);
}

/* ----------------------------------------------------------------------
   return 3-vector color corresponding to value from color map index
------------------------------------------------------------------------- */

double *Image::map_value2color(int index, double value)
{
  return maps[index]->value2color(value);
}

/* ----------------------------------------------------------------------
   add a new color to username and userrgb
   redefine RGB values in userrgb if name already exists
   return 1 if RGB values are invalid, else return 0
------------------------------------------------------------------------- */

int Image::addcolor(char *name, double r, double g, double b)
{
  int icolor;
  for (icolor = 0; icolor < ncolors; icolor++)
    if (strcmp(name,username[icolor]) == 0) break;

  if (icolor == ncolors) {
    username = (char **)
      memory->srealloc(username,(ncolors+1)*sizeof(char *),"image:username");
    memory->grow(userrgb,ncolors+1,3,"image:userrgb");
    ncolors++;
  }

  int n = strlen(name) + 1;
  username[icolor] = new char[n];
  strcpy(username[icolor],name);

  if (r < 0.0 || r > 1.0 || g < 0.0 || g > 1.0 || b < 0.0 || b > 1.0)
    return 1;

  userrgb[icolor][0] = r;
  userrgb[icolor][1] = g;
  userrgb[icolor][2] = b;

  return 0;
}

/* ----------------------------------------------------------------------
   if index > 0, return ptr to index-1 color from rgb
   if index < 0, return ptr to -index-1 color from userrgb
   if index = 0, search the 2 lists of color names for the string color
   search user-defined color names first, then the list of NCOLORS names
   return a pointer to the 3 floating point RGB values or NULL if didn't find
------------------------------------------------------------------------- */

double *Image::color2rgb(const char *color, int index)
{
  static const char *name[NCOLORS] = {
    "aliceblue",
    "antiquewhite",
    "aqua",
    "aquamarine",
    "azure",
    "beige",
    "bisque",
    "black",
    "blanchedalmond",
    "blue",
    "blueviolet",
    "brown",
    "burlywood",
    "cadetblue",
    "chartreuse",
    "chocolate",
    "coral",
    "cornflowerblue",
    "cornsilk",
    "crimson",
    "cyan",
    "darkblue",
    "darkcyan",
    "darkgoldenrod",
    "darkgray",
    "darkgreen",
    "darkkhaki",
    "darkmagenta",
    "darkolivegreen",
    "darkorange",
    "darkorchid",
    "darkred",
    "darksalmon",
    "darkseagreen",
    "darkslateblue",
    "darkslategray",
    "darkturquoise",
    "darkviolet",
    "deeppink",
    "deepskyblue",
    "dimgray",
    "dodgerblue",
    "firebrick",
    "floralwhite",
    "forestgreen",
    "fuchsia",
    "gainsboro",
    "ghostwhite",
    "gold",
    "goldenrod",
    "gray",
    "green",
    "greenyellow",
    "honeydew",
    "hotpink",
    "indianred",
    "indigo",
    "ivory",
    "khaki",
    "lavender",
    "lavenderblush",
    "lawngreen",
    "lemonchiffon",
    "lightblue",
    "lightcoral",
    "lightcyan",
    "lightgoldenrodyellow",
    "lightgreen",
    "lightgrey",
    "lightpink",
    "lightsalmon",
    "lightseagreen",
    "lightskyblue",
    "lightslategray",
    "lightsteelblue",
    "lightyellow",
    "lime",
    "limegreen",
    "linen",
    "magenta",
    "maroon",
    "mediumaquamarine",
    "mediumblue",
    "mediumorchid",
    "mediumpurple",
    "mediumseagreen",
    "mediumslateblue",
    "mediumspringgreen",
    "mediumturquoise",
    "mediumvioletred",
    "midnightblue",
    "mintcream",
    "mistyrose",
    "moccasin",
    "navajowhite",
    "navy",
    "oldlace",
    "olive",
    "olivedrab",
    "orange",
    "orangered",
    "orchid",
    "palegoldenrod",
    "palegreen",
    "paleturquoise",
    "palevioletred",
    "papayawhip",
    "peachpuff",
    "peru",
    "pink",
    "plum",
    "powderblue",
    "purple",
    "red",
    "rosybrown",
    "royalblue",
    "saddlebrown",
    "salmon",
    "sandybrown",
    "seagreen",
    "seashell",
    "sienna",
    "silver",
    "skyblue",
    "slateblue",
    "slategray",
    "snow",
    "springgreen",
    "steelblue",
    "tan",
    "teal",
    "thistle",
    "tomato",
    "turquoise",
    "violet",
    "wheat",
    "white",
    "whitesmoke",
    "yellow",
    "yellowgreen"
  };

  static double rgb[NCOLORS][3] = {
    {240/255.0, 248/255.0, 255/255.0},
    {250/255.0, 235/255.0, 215/255.0},
    {0/255.0, 255/255.0, 255/255.0},
    {127/255.0, 255/255.0, 212/255.0},
    {240/255.0, 255/255.0, 255/255.0},
    {245/255.0, 245/255.0, 220/255.0},
    {255/255.0, 228/255.0, 196/255.0},
    {0/255.0, 0/255.0, 0/255.0},
    {255/255.0, 255/255.0, 205/255.0},
    {0/255.0, 0/255.0, 255/255.0},
    {138/255.0, 43/255.0, 226/255.0},
    {165/255.0, 42/255.0, 42/255.0},
    {222/255.0, 184/255.0, 135/255.0},
    {95/255.0, 158/255.0, 160/255.0},
    {127/255.0, 255/255.0, 0/255.0},
    {210/255.0, 105/255.0, 30/255.0},
    {255/255.0, 127/255.0, 80/255.0},
    {100/255.0, 149/255.0, 237/255.0},
    {255/255.0, 248/255.0, 220/255.0},
    {220/255.0, 20/255.0, 60/255.0},
    {0/255.0, 255/255.0, 255/255.0},
    {0/255.0, 0/255.0, 139/255.0},
    {0/255.0, 139/255.0, 139/255.0},
    {184/255.0, 134/255.0, 11/255.0},
    {169/255.0, 169/255.0, 169/255.0},
    {0/255.0, 100/255.0, 0/255.0},
    {189/255.0, 183/255.0, 107/255.0},
    {139/255.0, 0/255.0, 139/255.0},
    {85/255.0, 107/255.0, 47/255.0},
    {255/255.0, 140/255.0, 0/255.0},
    {153/255.0, 50/255.0, 204/255.0},
    {139/255.0, 0/255.0, 0/255.0},
    {233/255.0, 150/255.0, 122/255.0},
    {143/255.0, 188/255.0, 143/255.0},
    {72/255.0, 61/255.0, 139/255.0},
    {47/255.0, 79/255.0, 79/255.0},
    {0/255.0, 206/255.0, 209/255.0},
    {148/255.0, 0/255.0, 211/255.0},
    {255/255.0, 20/255.0, 147/255.0},
    {0/255.0, 191/255.0, 255/255.0},
    {105/255.0, 105/255.0, 105/255.0},
    {30/255.0, 144/255.0, 255/255.0},
    {178/255.0, 34/255.0, 34/255.0},
    {255/255.0, 250/255.0, 240/255.0},
    {34/255.0, 139/255.0, 34/255.0},
    {255/255.0, 0/255.0, 255/255.0},
    {220/255.0, 220/255.0, 220/255.0},
    {248/255.0, 248/255.0, 255/255.0},
    {255/255.0, 215/255.0, 0/255.0},
    {218/255.0, 165/255.0, 32/255.0},
    {128/255.0, 128/255.0, 128/255.0},
    {0/255.0, 128/255.0, 0/255.0},
    {173/255.0, 255/255.0, 47/255.0},
    {240/255.0, 255/255.0, 240/255.0},
    {255/255.0, 105/255.0, 180/255.0},
    {205/255.0, 92/255.0, 92/255.0},
    {75/255.0, 0/255.0, 130/255.0},
    {255/255.0, 240/255.0, 240/255.0},
    {240/255.0, 230/255.0, 140/255.0},
    {230/255.0, 230/255.0, 250/255.0},
    {255/255.0, 240/255.0, 245/255.0},
    {124/255.0, 252/255.0, 0/255.0},
    {255/255.0, 250/255.0, 205/255.0},
    {173/255.0, 216/255.0, 230/255.0},
    {240/255.0, 128/255.0, 128/255.0},
    {224/255.0, 255/255.0, 255/255.0},
    {250/255.0, 250/255.0, 210/255.0},
    {144/255.0, 238/255.0, 144/255.0},
    {211/255.0, 211/255.0, 211/255.0},
    {255/255.0, 182/255.0, 193/255.0},
    {255/255.0, 160/255.0, 122/255.0},
    {32/255.0, 178/255.0, 170/255.0},
    {135/255.0, 206/255.0, 250/255.0},
    {119/255.0, 136/255.0, 153/255.0},
    {176/255.0, 196/255.0, 222/255.0},
    {255/255.0, 255/255.0, 224/255.0},
    {0/255.0, 255/255.0, 0/255.0},
    {50/255.0, 205/255.0, 50/255.0},
    {250/255.0, 240/255.0, 230/255.0},
    {255/255.0, 0/255.0, 255/255.0},
    {128/255.0, 0/255.0, 0/255.0},
    {102/255.0, 205/255.0, 170/255.0},
    {0/255.0, 0/255.0, 205/255.0},
    {186/255.0, 85/255.0, 211/255.0},
    {147/255.0, 112/255.0, 219/255.0},
    {60/255.0, 179/255.0, 113/255.0},
    {123/255.0, 104/255.0, 238/255.0},
    {0/255.0, 250/255.0, 154/255.0},
    {72/255.0, 209/255.0, 204/255.0},
    {199/255.0, 21/255.0, 133/255.0},
    {25/255.0, 25/255.0, 112/255.0},
    {245/255.0, 255/255.0, 250/255.0},
    {255/255.0, 228/255.0, 225/255.0},
    {255/255.0, 228/255.0, 181/255.0},
    {255/255.0, 222/255.0, 173/255.0},
    {0/255.0, 0/255.0, 128/255.0},
    {253/255.0, 245/255.0, 230/255.0},
    {128/255.0, 128/255.0, 0/255.0},
    {107/255.0, 142/255.0, 35/255.0},
    {255/255.0, 165/255.0, 0/255.0},
    {255/255.0, 69/255.0, 0/255.0},
    {218/255.0, 112/255.0, 214/255.0},
    {238/255.0, 232/255.0, 170/255.0},
    {152/255.0, 251/255.0, 152/255.0},
    {175/255.0, 238/255.0, 238/255.0},
    {219/255.0, 112/255.0, 147/255.0},
    {255/255.0, 239/255.0, 213/255.0},
    {255/255.0, 239/255.0, 213/255.0},
    {205/255.0, 133/255.0, 63/255.0},
    {255/255.0, 192/255.0, 203/255.0},
    {221/255.0, 160/255.0, 221/255.0},
    {176/255.0, 224/255.0, 230/255.0},
    {128/255.0, 0/255.0, 128/255.0},
    {255/255.0, 0/255.0, 0/255.0},
    {188/255.0, 143/255.0, 143/255.0},
    {65/255.0, 105/255.0, 225/255.0},
    {139/255.0, 69/255.0, 19/255.0},
    {250/255.0, 128/255.0, 114/255.0},
    {244/255.0, 164/255.0, 96/255.0},
    {46/255.0, 139/255.0, 87/255.0},
    {255/255.0, 245/255.0, 238/255.0},
    {160/255.0, 82/255.0, 45/255.0},
    {192/255.0, 192/255.0, 192/255.0},
    {135/255.0, 206/255.0, 235/255.0},
    {106/255.0, 90/255.0, 205/255.0},
    {112/255.0, 128/255.0, 144/255.0},
    {255/255.0, 250/255.0, 250/255.0},
    {0/255.0, 255/255.0, 127/255.0},
    {70/255.0, 130/255.0, 180/255.0},
    {210/255.0, 180/255.0, 140/255.0},
    {0/255.0, 128/255.0, 128/255.0},
    {216/255.0, 191/255.0, 216/255.0},
    {253/255.0, 99/255.0, 71/255.0},
    {64/255.0, 224/255.0, 208/255.0},
    {238/255.0, 130/255.0, 238/255.0},
    {245/255.0, 222/255.0, 179/255.0},
    {255/255.0, 255/255.0, 255/255.0},
    {245/255.0, 245/255.0, 245/255.0},
    {255/255.0, 255/255.0, 0/255.0},
    {154/255.0, 205/255.0, 50/255.0}
  };

  if (index > 0) {
    if (index > NCOLORS) return NULL;
    return rgb[index-1];
  }
  if (index < 0) {
    if (-index > ncolors) return NULL;
    return userrgb[-index-1];
  }

  if (color) {
    for (int i = 0; i < ncolors; i++)
      if (strcmp(color,username[i]) == 0) return userrgb[i];
    for (int i = 0; i < NCOLORS; i++)
      if (strcmp(color,name[i]) == 0) return rgb[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   return number of default colors
------------------------------------------------------------------------- */

int Image::default_colors()
{
  return NCOLORS;
}

/* ----------------------------------------------------------------------
   search the list of element names for the string element
   return a pointer to the 3 floating point RGB values
   this list is used by AtomEye and is taken from its Mendeleyev.c file
------------------------------------------------------------------------- */

double *Image::element2color(char *element)
{
  static const char *name[NELEMENTS] = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"
  };

  static double rgb[NELEMENTS][3] = {
    {0.8, 0.8, 0.8},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.7, 0.7, 0.7},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.9, 0.4, 0},
    {0.35, 0.35, 0.35},
    {0.2, 0.2, 0.8},
    {0.8, 0.2, 0.2},
    {0.7, 0.85, 0.45},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6, 0.6, 0.6},
    {0.6, 0.6, 0.7},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6901960784, 0.768627451, 0.8705882353},
    {0.1, 0.7, 0.3},
    {0.95, 0.9, 0.2},
    {0.15, 0.5, 0.1},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.5, 0.5, 0.5},
    {0.8, 0.8, 0.7},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0, 0.8, 0},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.5176470588, 0.5764705882, 0.6529411765},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.257254902, 0.2666666667, 0.271372549},
    {0.95, 0.7900735294, 0.01385869565},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.9, 0, 1},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {1, 1, 0.3},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.5, 0.08, 0.12},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.5, 0.1, 0.5},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.8, 0.8, 0},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {1, 0.8431372549, 0},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.9, 0.8, 0},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.8, 0.2, 0.2},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.1, 0.7, 0.3},
    {0.1, 0.3, 0.7},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.9, 0.8, 0},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6431372549, 0.6666666667, 0.6784313725}
  };

  for (int i = 0; i < NELEMENTS; i++)
    if (strcmp(element,name[i]) == 0) return rgb[i];
  return NULL;
}

/* ----------------------------------------------------------------------
   search the list of element names for the string element
   return a pointer to the 3 floating point RGB values
   this list is used by AtomEye and is taken from its Mendeleyev.c file
------------------------------------------------------------------------- */

double Image::element2diam(char *element)
{
  static const char *name[NELEMENTS] = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"
  };

  static double diameter[NELEMENTS] = {
    0.35, 1.785, 1.45, 1.05, 0.85, 0.72, 0.65, 0.6, 0.5, 1.5662,
    1.8, 1.5, 1.4255, 1.07, 1, 1, 1, 1.8597, 2.2, 1.8,
    1.6, 1.4, 1.51995, 1.44225, 1.4, 1.43325, 1.35, 1.35, 1.278, 1.35,
    1.3, 1.25, 1.15, 1.15, 1.15, 2.0223, 2.35, 2, 1.8, 1.55,
    1.6504, 1.3872, 1.35, 1.3, 1.35, 1.4, 1.6, 1.55, 1.55, 1.45,
    1.45, 1.4, 1.4, 2.192, 2.6, 2.15, 1.95, 1.85, 1.85, 1.85,
    1.85, 1.85, 1.85, 1.8, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
    1.75, 1.55, 1.6529, 1.5826, 1.35, 1.3, 1.35, 1.35, 1.35, 1.5,
    1.9, 1.8, 1.6, 1.9, 1.6, 1.0, 1.0, 2.15, 1.95, 1.8,
    1.8, 1.75, 1.75, 1.75, 1.75, 1.0, 1.0, 1.6, 1.6, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.6, 1.0, 1.0, 1.0, 1.0
  };

  for (int i = 0; i < NELEMENTS; i++)
    if (strcmp(element,name[i]) == 0) return diameter[i];
  return 0.0;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ColorMap class
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

ColorMap::ColorMap(LAMMPS *lmp, Image *caller) : Pointers(lmp)
{
  image = caller;

  // default color map

  dynamic = 1;

  mlo = MINVALUE;
  mhi = MAXVALUE;
  mstyle = CONTINUOUS;
  mrange = FRACTIONAL;

  nentry = 2;
  mentry = new MapEntry[nentry];
  mentry[0].single = MINVALUE;
  mentry[0].color = image->color2rgb("blue");
  mentry[1].single = MAXVALUE;
  mentry[1].color = image->color2rgb("red");
}

/* ---------------------------------------------------------------------- */

ColorMap::~ColorMap()
{
  delete [] mentry;
}

/* ----------------------------------------------------------------------
   redefine color map
   args = lo hi style delta N entry1 entry2 ... entryN as defined by caller
   return 1 if any error in args, else return 0
------------------------------------------------------------------------- */

int ColorMap::reset(int narg, char **arg)
{
  if (!islower(arg[0][0])) {
    mlo = NUMERIC;
    mlovalue = force->numeric(FLERR,arg[0]);
  } else if (strcmp(arg[0],"min") == 0) mlo = MINVALUE;
  else return 1;

  if (!islower(arg[1][0])) {
    mhi = NUMERIC;
    mhivalue = force->numeric(FLERR,arg[1]);
  } else if (strcmp(arg[1],"max") == 0) mhi = MAXVALUE;
  else return 1;

  if (mlo == NUMERIC && mhi == NUMERIC && mlovalue >= mhivalue) return 1;

  if (mlo == MINVALUE || mhi == MAXVALUE) dynamic = 1;
  else dynamic = 0;

  if (strlen(arg[2]) != 2) return 1;
  if (arg[2][0] == 'c') mstyle = CONTINUOUS;
  else if (arg[2][0] == 'd') mstyle = DISCRETE;
  else if (arg[2][0] == 's') mstyle = SEQUENTIAL;
  else return 1;
  if (arg[2][1] == 'a') mrange = ABSOLUTE;
  else if (arg[2][1] == 'f') mrange = FRACTIONAL;
  else return 1;

  if (mstyle == SEQUENTIAL) {
    mbinsize = force->numeric(FLERR,arg[3]);
    if (mbinsize <= 0.0) return 1;
    mbinsizeinv = 1.0/mbinsize;
  }

  nentry = force->inumeric(FLERR,arg[4]);
  if (nentry < 1) return 1;
  delete [] mentry;
  mentry = new MapEntry[nentry];

  int expandflag = 0;

  int n = 5;
  for (int i = 0; i < nentry; i++) {
    if (mstyle == CONTINUOUS) {
      if (n+2 > narg) return 1;
      if (!islower(arg[n][0])) {
        mentry[i].single = NUMERIC;
        mentry[i].svalue = force->numeric(FLERR,arg[n]);
      } else if (strcmp(arg[n],"min") == 0) mentry[i].single = MINVALUE;
      else if (strcmp(arg[n],"max") == 0) mentry[i].single = MAXVALUE;
      else return 1;
      mentry[i].color = image->color2rgb(arg[n+1]);
      n += 2;
    } else if (mstyle == DISCRETE) {
      if (n+3 > narg) return 1;
      if (!islower(arg[n][0])) {
        mentry[i].lo = NUMERIC;
        mentry[i].lvalue = force->numeric(FLERR,arg[n]);
      } else if (strcmp(arg[n],"min") == 0) mentry[i].lo = MINVALUE;
      else if (strcmp(arg[n],"max") == 0) mentry[i].lo = MAXVALUE;
      else return 1;
      if (!islower(arg[n+1][0])) {
        mentry[i].hi = NUMERIC;
        mentry[i].hvalue = force->numeric(FLERR,arg[n+1]);
      } else if (strcmp(arg[n+1],"min") == 0) mentry[i].hi = MINVALUE;
      else if (strcmp(arg[n+1],"max") == 0) mentry[i].hi = MAXVALUE;
      else return 1;
      mentry[i].color = image->color2rgb(arg[n+2]);
      n += 3;
    } else if (mstyle == SEQUENTIAL) {
      // NOTE: this is unfinished code, not sure how useful it is
      // idea is to allow a list of colors to be specified with a single arg
      // problem is that sequential colors in ALL are not very different
      // e.g. ALL or USER or ALL5:10 or USER1:10:2
      // current code is just 1st nentry values of ALL or USER
      // need to comment out error check in DumpImage::modify_param()
      //   for amap check on (narg < n) to get it to work
      // need to add extra logic here to check not accessing undefined colors
      if (i == 0) {
        if (n+1 > narg) return 1;
        if (strcmp(arg[n],"ALL") == 0) expandflag = 1;
        if (strcmp(arg[n],"USER") == 0) expandflag = 2;
      }
      if (expandflag == 0) {
        if (n+1 > narg) return 1;
        mentry[i].color = image->color2rgb(arg[n]);
      } else if (expandflag == 1) {
        mentry[i].color = image->color2rgb(NULL,i+1);
      } else if (expandflag == 2) {
        mentry[i].color = image->color2rgb(NULL,-(i+1));
      }
      n += 1;
    }
    if (mentry[i].color == NULL) return 1;
  }

  if (mstyle == CONTINUOUS) {
    if (nentry < 2) return 1;
    if (mentry[0].single != MINVALUE || mentry[nentry-1].single != MAXVALUE)
      return 1;
    for (int i = 2; i < nentry-1; i++)
      if (mentry[i].svalue <= mentry[i-1].svalue) return 1;
  } else if (mstyle == DISCRETE) {
    if (nentry < 1) return 1;
    if (mentry[nentry-1].lo != MINVALUE || mentry[nentry-1].hi != MAXVALUE)
      return 1;
  } else if (mstyle == SEQUENTIAL) {
    if (nentry < 1) return 1;
  }

  // one-time call to minmax if color map is static

  if (!dynamic) return minmax(mlovalue,mhivalue);

  return 0;
}

/* ----------------------------------------------------------------------
   set explicit values for all min/max settings in color map
     from min/max dynamic values
   set lo/hi current and lvalue/hvalue entries that are MIN/MAX VALUE
   called only once if mlo/mhi != MIN/MAX VALUE, else called repeatedly
   return 1 = error if any values now overlap incorrectly with dynamic bounds
   else return 0
------------------------------------------------------------------------- */

int ColorMap::minmax(double mindynamic, double maxdynamic)
{
  if (mlo == MINVALUE) locurrent = mindynamic;
  else locurrent = mlovalue;
  if (mhi == MAXVALUE) hicurrent = maxdynamic;
  else hicurrent = mhivalue;
  if (locurrent > hicurrent) return 1;

  if (mstyle == CONTINUOUS) {
    if (mrange == ABSOLUTE) mentry[0].svalue = locurrent;
    else mentry[0].svalue = 0.0;
    if (mrange == ABSOLUTE) mentry[nentry-1].svalue = hicurrent;
    else mentry[nentry-1].svalue = 1.0;

    // error in ABSOLUTE mode if new lo/hi current cause
    // first/last entry to become lo > hi with adjacent entry

    if (mrange == ABSOLUTE) {
      if (mentry[0].svalue > mentry[1].svalue) return 1;
      if (mentry[nentry-2].svalue > mentry[nentry-1].svalue) return 1;
    }

  // OK if new lo/hi current cause an entry to have lo > hi,
  // since last entry will always be a match

  } else if (mstyle == DISCRETE) {
    for (int i = 0; i < nentry; i++) {
      if (mentry[i].lo == MINVALUE) {
        if (mrange == ABSOLUTE) mentry[i].lvalue = locurrent;
        else mentry[i].lvalue = 0.0;
      }
      if (mentry[i].hi == MAXVALUE) {
        if (mrange == ABSOLUTE) mentry[i].hvalue = hicurrent;
        else mentry[i].hvalue = 1.0;
      }
    }
  }

  return 0;
}

/* ----------------------------------------------------------------------
   convert value into an RGB color via color map
   return pointer to 3-vector
------------------------------------------------------------------------- */

double *ColorMap::value2color(double value)
{
  double lo;//,hi;

  value = MAX(value,locurrent);
  value = MIN(value,hicurrent);

  if (mrange == FRACTIONAL) {
    if (locurrent == hicurrent) value = 0.0;
    else value = (value-locurrent) / (hicurrent-locurrent);
    lo = 0.0;
    //hi = 1.0;
  } else {
    lo = locurrent;
    //hi = hicurrent;
  }

  if (mstyle == CONTINUOUS) {
    for (int i = 0; i < nentry-1; i++)
      if (value >= mentry[i].svalue && value <= mentry[i+1].svalue) {
        double fraction = (value-mentry[i].svalue) /
          (mentry[i+1].svalue-mentry[i].svalue);
        interpolate[0] = mentry[i].color[0] +
          fraction*(mentry[i+1].color[0]-mentry[i].color[0]);
        interpolate[1] = mentry[i].color[1] +
          fraction*(mentry[i+1].color[1]-mentry[i].color[1]);
        interpolate[2] = mentry[i].color[2] +
          fraction*(mentry[i+1].color[2]-mentry[i].color[2]);
        return interpolate;
      }
  } else if (mstyle == DISCRETE) {
    for (int i = 0; i < nentry; i++)
      if (value >= mentry[i].lvalue && value <= mentry[i].hvalue)
        return mentry[i].color;
  } else {
    int ibin = static_cast<int> ((value-lo) * mbinsizeinv);
    return mentry[ibin%nentry].color;
  }

  return NULL;
}

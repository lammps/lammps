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

/* ----------------------------------------------------------------------
   Contributing authors: Nathan Fabian (Sandia)
                     and Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include "image.h"

#include "domain.h"
#include "error.h"
#include "image_objects.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "random_mars.h"
#include "version.h"

#include <cctype>
#include <cmath>
#include <cstring>

#ifdef LAMMPS_JPEG
#include <jpeglib.h>
#endif

#ifdef LAMMPS_PNG
#include <png.h>
#include <zlib.h>
#include <csetjmp>
#endif

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::MY_PI;
using MathConst::MY_PI4;

// clang-format on
namespace {
constexpr int NCOLORS = 140;
constexpr int NELEMENTS = 109;
constexpr double EPSILON = 1.0e-6;
constexpr double TRANS_DELTA = 0.01;

enum { NUMERIC, MINVALUE, MAXVALUE };
enum { CONTINUOUS, DISCRETE, SEQUENTIAL };
enum { ABSOLUTE, FRACTIONAL };
enum { NO, YES };

struct TGAHeader {
  unsigned char idlength;
  unsigned char colormaptype;
  unsigned char datatypecode;
  unsigned char colormaporigin[2];
  unsigned char colormaplength[2];
  unsigned char colormapdepth;
  unsigned char x_origin[2];
  unsigned char y_origin[2];
  unsigned char width[2];
  unsigned char height[2];
  unsigned char bitsperpixel;
  unsigned char imagedescriptor;
};

struct TGAFooter {
  unsigned int extoffset;
  unsigned int devoffset;
  char signature[18];
};

////////////////////////////////////////////////////////////////////////
// the following regular Bayer threshold matrix can be created for any
// power of 2 ranks with the following python code.
// See https://en.wikipedia.org/wiki/Ordered_dithering
//
// import numpy as np
//
// def bayer_matrix(n: int) -> np.ndarray:
//     """Generate an n x n Bayer (ordered dither) matrix for n a power of 2."""
//     if n == 1:
//         return np.array([[0]], dtype=int)
//
//     m = n // 2
//     B = bayer_matrix(m)
//
//     # Standard Bayer recursion
//     return np.block([
//         [4 * B + 0, 4 * B + 2],
//         [4 * B + 3, 4 * B + 1]
//     ])
//
// # set rank
// n=16
//
// # normalized matrix from recursion
// matrix = (bayer_matrix(n) + 0.5) / (n * n)
//
// print("constexpr int TRANK = %d;" % n)
// print("constexpr double transthresh[TRANK][TRANK] = {")
// nx, ny = matrix.shape
// for iy in range(0,ny):
//     print("{")
//     for ix in range(0,nx):
//         if ix < nx -1:
//             print(f"{matrix[ix][iy]:12.9f},")
//         else:
//             print(f"{matrix[ix][iy]:12.9f}")
//     if iy < ny - 1:
//         print("},")
//     else:
//         print("}")
//
// print("};")
////////////////////////////////////////////////////////////////////////
constexpr int TRANK = 16;
constexpr double transthresh[TRANK][TRANK] = {
    {0.001953125, 0.751953125, 0.189453125, 0.939453125, 0.048828125, 0.798828125, 0.236328125,
     0.986328125, 0.013671875, 0.763671875, 0.201171875, 0.951171875, 0.060546875, 0.810546875,
     0.248046875, 0.998046875},
    {0.501953125, 0.251953125, 0.689453125, 0.439453125, 0.548828125, 0.298828125, 0.736328125,
     0.486328125, 0.513671875, 0.263671875, 0.701171875, 0.451171875, 0.560546875, 0.310546875,
     0.748046875, 0.498046875},
    {0.126953125, 0.876953125, 0.064453125, 0.814453125, 0.173828125, 0.923828125, 0.111328125,
     0.861328125, 0.138671875, 0.888671875, 0.076171875, 0.826171875, 0.185546875, 0.935546875,
     0.123046875, 0.873046875},
    {0.626953125, 0.376953125, 0.564453125, 0.314453125, 0.673828125, 0.423828125, 0.611328125,
     0.361328125, 0.638671875, 0.388671875, 0.576171875, 0.326171875, 0.685546875, 0.435546875,
     0.623046875, 0.373046875},
    {0.033203125, 0.783203125, 0.220703125, 0.970703125, 0.017578125, 0.767578125, 0.205078125,
     0.955078125, 0.044921875, 0.794921875, 0.232421875, 0.982421875, 0.029296875, 0.779296875,
     0.216796875, 0.966796875},
    {0.533203125, 0.283203125, 0.720703125, 0.470703125, 0.517578125, 0.267578125, 0.705078125,
     0.455078125, 0.544921875, 0.294921875, 0.732421875, 0.482421875, 0.529296875, 0.279296875,
     0.716796875, 0.466796875},
    {0.158203125, 0.908203125, 0.095703125, 0.845703125, 0.142578125, 0.892578125, 0.080078125,
     0.830078125, 0.169921875, 0.919921875, 0.107421875, 0.857421875, 0.154296875, 0.904296875,
     0.091796875, 0.841796875},
    {0.658203125, 0.408203125, 0.595703125, 0.345703125, 0.642578125, 0.392578125, 0.580078125,
     0.330078125, 0.669921875, 0.419921875, 0.607421875, 0.357421875, 0.654296875, 0.404296875,
     0.591796875, 0.341796875},
    {0.009765625, 0.759765625, 0.197265625, 0.947265625, 0.056640625, 0.806640625, 0.244140625,
     0.994140625, 0.005859375, 0.755859375, 0.193359375, 0.943359375, 0.052734375, 0.802734375,
     0.240234375, 0.990234375},
    {0.509765625, 0.259765625, 0.697265625, 0.447265625, 0.556640625, 0.306640625, 0.744140625,
     0.494140625, 0.505859375, 0.255859375, 0.693359375, 0.443359375, 0.552734375, 0.302734375,
     0.740234375, 0.490234375},
    {0.134765625, 0.884765625, 0.072265625, 0.822265625, 0.181640625, 0.931640625, 0.119140625,
     0.869140625, 0.130859375, 0.880859375, 0.068359375, 0.818359375, 0.177734375, 0.927734375,
     0.115234375, 0.865234375},
    {0.634765625, 0.384765625, 0.572265625, 0.322265625, 0.681640625, 0.431640625, 0.619140625,
     0.369140625, 0.630859375, 0.380859375, 0.568359375, 0.318359375, 0.677734375, 0.427734375,
     0.615234375, 0.365234375},
    {0.041015625, 0.791015625, 0.228515625, 0.978515625, 0.025390625, 0.775390625, 0.212890625,
     0.962890625, 0.037109375, 0.787109375, 0.224609375, 0.974609375, 0.021484375, 0.771484375,
     0.208984375, 0.958984375},
    {0.541015625, 0.291015625, 0.728515625, 0.478515625, 0.525390625, 0.275390625, 0.712890625,
     0.462890625, 0.537109375, 0.287109375, 0.724609375, 0.474609375, 0.521484375, 0.271484375,
     0.708984375, 0.458984375},
    {0.166015625, 0.916015625, 0.103515625, 0.853515625, 0.150390625, 0.900390625, 0.087890625,
     0.837890625, 0.162109375, 0.912109375, 0.099609375, 0.849609375, 0.146484375, 0.896484375,
     0.083984375, 0.833984375},
    {0.666015625, 0.416015625, 0.603515625, 0.353515625, 0.650390625, 0.400390625, 0.587890625,
     0.337890625, 0.662109375, 0.412109375, 0.599609375, 0.349609375, 0.646484375, 0.396484375,
     0.583984375, 0.333984375}};

// function to apply bilinear scaling to pixmap
void scale_pixmap(int ow, int oh, const unsigned char *opix, int nw, int nh, unsigned char *npix)
{
  double x_ratio = (double) (ow - 1) / nw;
  double y_ratio = (double) (oh - 1) / nh;

  for (int i = 0; i < nh; i++) {
    for (int j = 0; j < nw; j++) {
      int x = (int) (x_ratio * j);
      int y = (int) (y_ratio * i);
      double x_diff = (x_ratio * j) - x;
      double y_diff = (y_ratio * i) - y;

      // Get the four neighboring pixels in the original pixmap
      int offs = y * 3 * ow + 3 * x;
      unsigned char a[3] = {opix[offs], opix[offs + 1], opix[offs + 2]};
      offs = y * 3 * ow + 3 * (x + 1);
      unsigned char b[3] = {opix[offs], opix[offs + 1], opix[offs + 2]};
      offs = (y + 1) * 3 * ow + 3 * x;
      unsigned char c[3] = {opix[offs], opix[offs + 1], opix[offs + 2]};
      offs = (y + 1) * 3 * ow + 3 * (x + 1);
      unsigned char d[3] = {opix[offs], opix[offs + 1], opix[offs + 2]};

      // interpolate R, G, and B channels separately with bilinear scaling
      npix[i * 3 * nw + 3 * j] =
          (unsigned char) (a[0] * (1 - x_diff) * (1 - y_diff) + b[0] * (x_diff) * (1 - y_diff) +
                           c[0] * (y_diff) * (1 - x_diff) + d[0] * (x_diff * y_diff));

      npix[i * 3 * nw + 3 * j + 1] =
          (unsigned char) (a[1] * (1 - x_diff) * (1 - y_diff) + b[1] * (x_diff) * (1 - y_diff) +
                           c[1] * (y_diff) * (1 - x_diff) + d[1] * (x_diff * y_diff));

      npix[i * 3 * nw + 3 * j + 2] =
          (unsigned char) (a[2] * (1 - x_diff) * (1 - y_diff) + b[2] * (x_diff) * (1 - y_diff) +
                           c[2] * (y_diff) * (1 - x_diff) + d[2] * (x_diff * y_diff));
    }
  }
}

// convert an RGB color to YUV colorspace
void rgb2yuv(int *rgb, int *yuv)
{
  yuv[0] = static_cast<int>((0.299 * rgb[0]) + (0.587 * rgb[1]) + (0.114 * rgb[2]));
  yuv[1] = static_cast<int>(-(0.14713 * rgb[0]) - (0.28886 * rgb[1]) + (0.436 * rgb[2]));
  yuv[2] = static_cast<int>((0.615 * rgb[0]) - (0.51499 * rgb[1]) - (0.10001 * rgb[2]));
}

// convert XPM-like bitmap to 8-bit pixmap
void xpm2pix(int width, int height, const char *xpm, unsigned char *pix, double *fg, double *bg)
{
  for (int j = height; j > 0; --j) {
    for (int i = 0; i < width; ++i) {
      double *color = (xpm[(j - 1) * width + i] == '#') ? fg : bg;
      int idx = (height - j) * width * 3 + 3 * i;
      pix[idx] = static_cast<unsigned char>(color[0] * 255.0);
      pix[idx + 1] = static_cast<unsigned char>(color[1] * 255.0);
      pix[idx + 2] = static_cast<unsigned char>(color[2] * 255.0);
    }
  }
}
// clang-format off

// the following are 32x32 pixel XPM-like bitmaps of the letters X, Y, and Z for the axis labels.
constexpr char letter_x[] = {
  "                                "
  "                                "
  "    #####             #####     "
  "     #####           #####      "
  "      #####         #####       "
  "      #####         #####       "
  "       #####       #####        "
  "        #####     #####         "
  "        #####     #####         "
  "         #####   #####          "
  "          ##### #####           "
  "          ##### #####           "
  "           #########            "
  "            #######             "
  "            #######             "
  "             #####              "
  "            ######              "
  "            #######             "
  "           ########             "
  "          ##########            "
  "          ##### #####           "
  "         #####   ####           "
  "         ####    #####          "
  "        #####     #####         "
  "       #####       #####        "
  "       ####        #####        "
  "      #####         #####       "
  "     #####           #####      "
  "     ####            #####      "
  "    #####             #####     "
  "   #####               #####    "
  "                                "};

constexpr char letter_y[] = {
  "                                "
  "                                "
  "    #####              #####    "
  "     #####            #####     "
  "      ####            ####      "
  "      #####          #####      "
  "       #####        #####       "
  "        ####        ####        "
  "        #####      #####        "
  "         #####    #####         "
  "          ####    ####          "
  "          #####  #####          "
  "           ##########           "
  "            ########            "
  "            ########            "
  "             ######             "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "              ####              "
  "                                "};

constexpr char letter_z[] = {
  "                                "
  "                                "
  "    #######################     "
  "    #######################     "
  "    #######################     "
  "                     ######     "
  "                     #####      "
  "                    #####       "
  "                   #####        "
  "                  ######        "
  "                  #####         "
  "                 #####          "
  "                #####           "
  "               #####            "
  "              ######            "
  "              #####             "
  "             #####              "
  "            #####               "
  "           #####                "
  "           #####                "
  "          #####                 "
  "         #####                  "
  "        #####                   "
  "       ######                   "
  "       #####                    "
  "      #####                     "
  "     #####                      "
  "    #####                       "
  "    ########################    "
  "    ########################    "
  "    ########################    "
  "                                "};

}    // namespace

/* ---------------------------------------------------------------------- */

Image::Image(LAMMPS *lmp, int nmap_caller) :
  Pointers(lmp), maps(nullptr), depthBuffer(nullptr), surfaceBuffer(nullptr), depthcopy(nullptr),
  surfacecopy(nullptr), imageBuffer(nullptr), rgbcopy(nullptr), writeBuffer(nullptr),
  recvcounts(nullptr), displs(nullptr), username(nullptr), userrgb(nullptr), random(nullptr)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // defaults for 3d viz

  width = height = 512;
  theta = 60.0 * DEG2RAD;
  phi = 30.0 * DEG2RAD;
  zoom = 1.0;
  shiny = 1.0;
  ssao = NO;
  fsaa = NO;

  up[0] = 0.0;
  up[1] = 0.0;
  up[2] = 1.0;

  // colors

  ncolors = 0;
  boxcolor = color2rgb("yellow");
  background[0] = background[1] = background[2] = 0.0;
  background2[0] = background2[1] = background2[2] = -1.0;

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
}

/* ---------------------------------------------------------------------- */

Image::~Image()
{
  for (int i = 0; i < nmap; i++) delete maps[i];
  delete[] maps;

  for (int i = 0; i < ncolors; i++) delete [] username[i];
  memory->sfree(username);
  memory->destroy(userrgb);

  memory->destroy(depthBuffer);
  memory->destroy(surfaceBuffer);
  memory->destroy(imageBuffer);
  memory->destroy(depthcopy);
  memory->destroy(surfacecopy);
  memory->destroy(rgbcopy);

  delete random;

  memory->destroy(recvcounts);
  memory->destroy(displs);
}

/* ----------------------------------------------------------------------
   allocate image and depth buffers
   called after image size is set
------------------------------------------------------------------------- */

void Image::buffers()
{
  memory->destroy(depthBuffer);
  memory->destroy(surfaceBuffer);
  memory->destroy(imageBuffer);
  memory->destroy(depthcopy);
  memory->destroy(surfacecopy);
  memory->destroy(rgbcopy);

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
  // try to ensure continuous images as changing view passes thru up
  // sufficient to handle common cases where theta = 0 or 180 is degenerate?

  double dot = MathExtra::dot3(up,camDir);
  if (fabs(dot) > (1.0 - EPSILON)) {
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
    SSAOSamples = static_cast<int>(8.0 + 32.0*ssaoint);
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
   create background gradient, if background2[0] is >= 0
     otherwise use single background color
------------------------------------------------------------------------- */

void Image::clear()
{
  int red   = background[0];
  int green = background[1];
  int blue  = background[2];

  if (background2[0] < 0.0) {
    for (int iy = 0; iy < height; ++iy) {
      for (int ix = 0; ix < width; ++ix) {
        imageBuffer[iy * width * 3 + ix * 3 + 0] = red;
        imageBuffer[iy * width * 3 + ix * 3 + 1] = green;
        imageBuffer[iy * width * 3 + ix * 3 + 2] = blue;
        depthBuffer[iy * width + ix] = -1;
      }
    }
  } else {
    for (int iy = 0; iy < height; ++iy) {
      double fraction = (double) iy / (double) height;
      red   = static_cast<int>(fraction * background2[0] + (1.0 - fraction) * background[0]);
      green = static_cast<int>(fraction * background2[1] + (1.0 - fraction) * background[1]);
      blue  = static_cast<int>(fraction * background2[2] + (1.0 - fraction) * background[2]);
      for (int ix = 0; ix < width; ++ix) {
        imageBuffer[iy * width * 3 + ix * 3 + 0] = red;
        imageBuffer[iy * width * 3 + ix * 3 + 1] = green;
        imageBuffer[iy * width * 3 + ix * 3 + 2] = blue;
        depthBuffer[iy * width + ix] = -1;
      }
    }
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
        if (depthBuffer[i] < 0 || (depthcopy[i] >= 0 && depthcopy[i] < depthBuffer[i])) {
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
  // MPI_Gather() result back to proc 0
  // use Gatherv() if subset of pixels is not the same size on every proc

  if (ssao) {
    MPI_Bcast(imageBuffer,npixels*3,MPI_BYTE,0,world);
    MPI_Bcast(surfaceBuffer,npixels*2,MPI_DOUBLE,0,world);
    MPI_Bcast(depthBuffer,npixels,MPI_DOUBLE,0,world);
    compute_SSAO();

    int pixelstart = 3 * static_cast<int>(1.0*me/nprocs * npixels);
    int pixelstop = 3 * static_cast<int>(1.0*(me+1)/nprocs * npixels);
    int mypixels = pixelstop - pixelstart;

    if (npixels % nprocs == 0) {
      MPI_Gather(imageBuffer+pixelstart,mypixels,MPI_BYTE,
                 rgbcopy,mypixels,MPI_BYTE,0,world);

    } else {
      if (recvcounts == nullptr) {
        memory->create(recvcounts,nprocs,"image:recvcounts");
        memory->create(displs,nprocs,"image:displs");
        MPI_Allgather(&mypixels,1,MPI_INT,recvcounts,1,MPI_INT,world);
        displs[0] = 0;
        for (int i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
      }

      MPI_Gatherv(imageBuffer+pixelstart,mypixels,MPI_BYTE,
                  rgbcopy,recvcounts,displs,MPI_BYTE,0,world);
    }

    writeBuffer = rgbcopy;
  } else {
    writeBuffer = imageBuffer;
  }

  // scale down image for antialiasing. can be done in place with simple averaging
  if (fsaa) {
    for (int h=0; h < height; h += 2) {
      for (int w=0; w < width; w +=2) {
        int idx1 = 3*width*h + 3*w;
        int idx2 = 3*width*h + 3*(w+1);
        int idx3 = 3*width*(h+1) + 3*w;
        int idx4 = 3*width*(h+1) + 3*(w+1);

        int out = 3*(width/2)*(h/2) + 3*(w/2);
        for (int i=0; i < 3; ++i) {
          writeBuffer[out+i] = (unsigned char) (0.25*((int)writeBuffer[idx1+i]
                                                      +(int)writeBuffer[idx2+i]
                                                      +(int)writeBuffer[idx3+i]
                                                      +(int)writeBuffer[idx4+i]));
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   draw simulation bounding box as 12 cylinders
------------------------------------------------------------------------- */

void Image::draw_box(double (*corners)[3], double diameter, double opacity)
{
  if ((diameter <= 0.0) || (opacity <= 0.0)) return;  // nothing to do

  draw_cylinder(corners[0],corners[1],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[2],corners[3],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[0],corners[2],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[1],corners[3],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[0],corners[4],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[1],corners[5],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[2],corners[6],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[3],corners[7],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[4],corners[5],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[6],corners[7],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[4],corners[6],boxcolor,diameter,3,opacity);
  draw_cylinder(corners[5],corners[7],boxcolor,diameter,3,opacity);
}

/* ----------------------------------------------------------------------
   draw XYZ axes with X, Y, and Z labels in red/green/blue
   axes = 4 end points
   user arrow image object to draw axes as arrows

   labels are created from 32x32 size bitmaps stored as string in an XPM-like format
   convert the bitmap into an RGB pixmap with foreground and transparent background
      choose text and background color to be somewhat similar to minimize artifacts from scaling
   switch colors from white/silver to black/darkgray depending on the luminance of the background
   offset labels in every direction to avoid them being obscured by the arrows
   scale letter pixmaps based on the smaller of image width or height
------------------------------------------------------------------------- */

void Image::draw_axes(double (*axes)[3], double diameter, double opacity)
{
  if ((diameter <= 0.0) || (opacity <= 0.0)) return;  // nothing to do

  // draw arrows

  const double radius = 0.5 * diameter;
  draw_sphere(axes[0], color2rgb("gray"), radius, opacity);
  ImageObjects::ArrowObj arrow;
  arrow.draw(this, color2rgb("red"), axes[0], axes[1], radius, opacity);
  arrow.draw(this, color2rgb("green"), axes[0], axes[2], radius, opacity);
  if (domain->dimension == 3)
    arrow.draw(this, color2rgb("blue"), axes[0], axes[3], radius, opacity);

  // adjust size of labels based on image size,
  // with FSAA active, width and height are doubled; adjust the scale factor accordingly

  double scale = static_cast<double>(MIN(width, height)) / 1440.0;
  if (fsaa) scale *= 0.5;

  // determine color of labels

  double *fontcolor = color2rgb("white");
  double *backcolor = color2rgb("silver");
  int bgyuv[3];
  rgb2yuv(background, bgyuv);
  if (bgyuv[0] > 192) {    // switch to black text only for very bright backgrounds
    fontcolor = color2rgb("black");
    backcolor = color2rgb("darkgray");
  }

  // convert bitmap of letters to pixmap and scale/draw.

  unsigned char rgbbuffer[32 * 32 * 3];
  double shiftedpos[3];
  constexpr double DIROFFS = 0.05;
  xpm2pix(32, 32, letter_x, rgbbuffer, fontcolor, backcolor);
  shiftedpos[0] = axes[1][0] + DIROFFS * (axes[1][0] - axes[0][0]);
  shiftedpos[1] = axes[1][1] + radius;
  shiftedpos[2] = axes[1][2] - radius;    // moving in lower z-direction reduces overlap for X
  draw_pixmap(shiftedpos, 32, 32, rgbbuffer, backcolor, scale, opacity);

  xpm2pix(32, 32, letter_y, rgbbuffer, fontcolor, backcolor);
  shiftedpos[0] = axes[2][0] + radius;
  shiftedpos[1] = axes[2][1] + DIROFFS * (axes[2][1] - axes[0][1]);
  shiftedpos[2] = axes[2][2] + radius;
  draw_pixmap(shiftedpos, 32, 32, rgbbuffer, backcolor, scale, opacity);

  if (domain->dimension == 3) {
    xpm2pix(32, 32, letter_z, rgbbuffer, fontcolor, backcolor);
    shiftedpos[0] = axes[3][0] + radius;
    shiftedpos[1] = axes[3][1] + radius;
    shiftedpos[2] = axes[3][2] + DIROFFS * (axes[3][2] - axes[0][2]);
    draw_pixmap(shiftedpos, 32, 32, rgbbuffer, backcolor, scale, opacity);
  }
}

/* ----------------------------------------------------------------------
   scale and add pixmap centered at location x into image with depth buffering
   background color indicates transparency and pixels in that color are skipped
------------------------------------------------------------------------- */

void Image::draw_pixmap(const double *x, int pixwidth, int pixheight, const unsigned char *pixmap,
                        double *transcolor, double scale, double opacity)
{
  // nothing to do
  if (!pixmap || (pixwidth == 0) || (pixheight == 0) || (scale <= 0.0) || (opacity <= 0.0)) return;

  double xlocal[3] = {x[0] - xctr, x[1] - yctr, x[2] - zctr};
  double xmap = MathExtra::dot3(camRight,xlocal);
  double ymap = MathExtra::dot3(camUp,xlocal);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(xlocal,camDir);

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist : -tanPerPixel / zoom;
  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int>(xf);
  int yc = static_cast<int>(yf);

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  // convert back to non-FSAA image coordinates, so we can re-use the pixmap drawing code
  if (fsaa) {
    xc /= 2;
    yc /= 2;
  }

  draw_pixmap(xc, yc, pixwidth, pixheight, pixmap, transcolor, scale, opacity, dist);
}

/* ----------------------------------------------------------------------
   scale and add pixmap centered at location xc, yc in image coordinates to image
   background color indicates transparency and pixels in that color are skipped
------------------------------------------------------------------------- */

void Image::draw_pixmap(int xc, int yc, int pixwidth, int pixheight, const unsigned char *pixmap,
                        double *transcolor, double scale, double opacity, double dist)
{
  if (opacity <= 0.0) return;  // nothing to do

  const unsigned char *mypixmap = pixmap;
  unsigned char *npixmap = nullptr;

  // adjust scale factor and image location for FSAA
  if (fsaa) {
    scale *= 2.0;
    xc *= 2;
    yc *= 2;
  }

  // only scale as much as needed.
  if (scale != 1.0) {
    int nwidth = std::lround(scale * pixwidth + 0.5);
    int nheight = std::lround(scale * pixheight + 0.5);
    npixmap = new unsigned char[3*nwidth*nheight];
    scale_pixmap(pixwidth, pixheight, pixmap, nwidth, nheight, npixmap);
    mypixmap = npixmap;
    pixwidth = nwidth;
    pixheight = nheight;
  }

  int ylo = yc;
  int xlo = xc;
  double normal[3] = {0.0, 0.0, 1.0};

  for (int j = 0; j < pixheight; ++j) {
    for (int i = 0; i < pixwidth; ++i) {
      int iy = ylo + j - pixheight/2;
      int ix = xlo + i - pixwidth/2;
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;
      if ((opacity < 1.0) && (transthresh[ix % TRANK][iy % TRANK] > opacity))
        continue;

      // get color of pixel at x/y position of pixmap

      double pixelcolor[3];
      int offs = 3*j*pixwidth + 3*i;
      pixelcolor[0] = (double)mypixmap[offs] / 255.0;
      pixelcolor[1] = (double)mypixmap[offs + 1] / 255.0;
      pixelcolor[2] = (double)mypixmap[offs + 2] / 255.0;

      // check for transparency color and skip if it matches
      // we allow a few steps difference for each channel to account
      // for rounding errors and reduce "bleeding" from interpolation

      if ((fabs(pixelcolor[0] - transcolor[0]) < TRANS_DELTA) &&
          (fabs(pixelcolor[1] - transcolor[1]) < TRANS_DELTA) &&
          (fabs(pixelcolor[2] - transcolor[2]) < TRANS_DELTA)) continue;

      draw_pixel(ix, iy, dist, normal, pixelcolor);
    }
  }
  delete[] npixmap;
}

/* ----------------------------------------------------------------------
   draw sphere at x with surfaceColor and diameter
   render pixel by pixel onto image plane with depth buffering
------------------------------------------------------------------------- */

void Image::draw_sphere(const double *x, const double *surfaceColor, double diameter,
                        double opacity)
{
  if ((diameter <= 0.0) || (opacity <= 0.0)) return;  // nothing to do

  double xlocal[3];

  xlocal[0] = x[0] - xctr;
  xlocal[1] = x[1] - yctr;
  xlocal[2] = x[2] - zctr;

  double xmap = MathExtra::dot3(camRight,xlocal);
  double ymap = MathExtra::dot3(camUp,xlocal);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(xlocal,camDir);

  double radius = 0.5*diameter;
  double radsq = radius*radius;
  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist : -tanPerPixel / zoom;
  double pixelRadiusFull = radius / pixelWidth;
  int pixelRadius = std::lround(pixelRadiusFull) + 1;

  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int>(xf);
  int yc = static_cast<int>(yf);
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  for (int iy = yc - pixelRadius; iy <= yc + pixelRadius; iy++) {
    for (int ix = xc - pixelRadius; ix <= xc + pixelRadius; ix++) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;
      if ((opacity < 1.0) && (transthresh[ix % TRANK][iy % TRANK] > opacity))
        continue;

      double surface[3];
      surface[1] = ((iy - yc) - height_error) * pixelWidth;
      surface[0] = ((ix - xc) - width_error) * pixelWidth;
      double projRad = surface[0]*surface[0] + surface[1]*surface[1];

      // outside the sphere in the projected image

      if (projRad > radsq) continue;
      surface[2] = sqrt(radsq - projRad);
      double depth = dist - surface[2];

      surface[0] /= radius;
      surface[1] /= radius;
      surface[2] /= radius;

      draw_pixel(ix, iy, depth, surface, surfaceColor);
    }
  }
}

/* ----------------------------------------------------------------------
   draw axis oriented cube at x with surfaceColor and diameter in size
   render pixel by pixel onto image plane with depth buffering
------------------------------------------------------------------------- */

void Image::draw_cube(const double *x, const double *surfaceColor, double diameter, double opacity)
{
  if ((diameter <= 0.0) || (opacity <= 0.0)) return;  // nothing to do

  double xlocal[3],surface[3];
  double normal[3] = {0.0, 0.0, 1.0};
  double t = 1.0;
  double tdir[3] = {0.5, 0.5, 0.0};
  double depth;

  xlocal[0] = x[0] - xctr;
  xlocal[1] = x[1] - yctr;
  xlocal[2] = x[2] - zctr;

  double xmap = MathExtra::dot3(camRight,xlocal);
  double ymap = MathExtra::dot3(camUp,xlocal);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(xlocal,camDir);

  double radius = 0.5*diameter;
  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist : -tanPerPixel / zoom;

  double halfWidth = diameter;
  double pixelHalfWidthFull = halfWidth / pixelWidth;
  int pixelHalfWidth = std::lround(pixelHalfWidthFull);

  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int>(xf);
  int yc = static_cast<int>(yf);
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  for (int iy = yc - pixelHalfWidth; iy <= yc + pixelHalfWidth; ++iy) {
    for (int ix = xc - pixelHalfWidth; ix <= xc + pixelHalfWidth; ++ix) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;
      if ((opacity < 1.0) && (transthresh[ix % TRANK][iy % TRANK] > opacity))
        continue;

      double sy = ((iy - yc) - height_error) * pixelWidth;
      double sx = ((ix - xc) - width_error) * pixelWidth;
      surface[0] = camRight[0] * sx + camUp[0] * sy;
      surface[1] = camRight[1] * sx + camUp[1] * sy;
      surface[2] = camRight[2] * sx + camUp[2] * sy;

      // iterate through each of the 6 axis-oriented planes of the box
      // only render up to 3 which are facing the camera
      // these checks short circuit a dot product, testing for > 0

      for (int dim = 0; dim < 3; ++dim) {
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
              draw_pixel(ix, iy, depth, normal, surfaceColor);
            }
            break;
          case 1:
            if (xin & zin) {
              depth = dist - t;
              draw_pixel(ix, iy, depth, normal, surfaceColor);
            }
            break;
          case 2:
            if (xin & yin) {
              depth = dist - t;
              draw_pixel(ix, iy, depth, normal, surfaceColor);
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

void Image::draw_cylinder(const double *x, const double *y,
                          const double *surfaceColor, double diameter, int sflag, double opacity)
{
  if ((diameter <= 0.0) || (opacity <= 0.0)) return;  // nothing to do

  double mid[3],xaxis[3],yaxis[3],zaxis[3];
  double camLDir[3], camLRight[3], camLUp[3];
  double zmin, zmax;

  if (sflag % 2) draw_sphere(x,surfaceColor,diameter,opacity);
  if (sflag / 2) draw_sphere(y,surfaceColor,diameter,opacity);
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
  if (len == 0.0) return;       // nothing left to do
  MathExtra::scale3(1.0/len,zaxis);
  len *= 0.5;
  zmax = len;
  zmin = -len;

  double xmap = MathExtra::dot3(camRight,mid);
  double ymap = MathExtra::dot3(camUp,mid);
  double dist = MathExtra::dot3(camPos,camDir) - MathExtra::dot3(mid,camDir);

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist : -tanPerPixel / zoom;

  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int>(xf);
  int yc = static_cast<int>(yf);
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  double pixelHalfWidthFull = (rasterWidth * 0.5) / pixelWidth;
  double pixelHalfHeightFull = (rasterHeight * 0.5) / pixelWidth;
  int pixelHalfWidth = std::lround(pixelHalfWidthFull);
  int pixelHalfHeight = std::lround(pixelHalfHeightFull);

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

  for (int iy = yc - pixelHalfHeight; iy <= yc + pixelHalfHeight; ++iy) {
    for (int ix = xc - pixelHalfWidth; ix <= xc + pixelHalfWidth; ++ix) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;
      if ((opacity < 1.0) && (transthresh[ix % TRANK][iy % TRANK] > opacity))
        continue;

      double surface[3], normal[3];
      double sy = ((iy - yc) - height_error) * pixelWidth;
      double sx = ((ix - xc) - width_error) * pixelWidth;
      surface[0] = camLRight[0] * sx + camLUp[0] * sy;
      surface[1] = camLRight[1] * sx + camLUp[1] * sy;
      surface[2] = camLRight[2] * sx + camLUp[2] * sy;

      double b = 2 * camLDir[0] * surface[0];
      double c = surface[0] * surface[0] + surface[1] * surface[1] - radsq;

      double partial = b*b - 4*a*c;
      if ((partial < 0.0) || (a == 0.0)) continue;
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

      surface[0] = MathExtra::dot3(normal, camLRight);
      surface[1] = MathExtra::dot3(normal, camLUp);
      surface[2] = MathExtra::dot3(normal, camLDir);

      double depth = dist - t;
      draw_pixel(ix, iy, depth, surface, surfaceColor);
    }
  }
}

/* ----------------------------------------------------------------------
   draw triangle with 3 corner points x,y,z, surfaceColor
------------------------------------------------------------------------- */

void Image::draw_triangle(const double *x, const double *y, const double *z,
                          const double *surfaceColor, const double opacity)
{
  if (opacity <= 0.0) return;  // nothing to do

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

  MathExtra::sub3(xlocal, ylocal, d1);
  d1len = MathExtra::len3(d1);
  if (d1len == 0.0) return;     // zero length of triangle side
  MathExtra::scale3(1.0 / d1len, d1);

  MathExtra::sub3(zlocal, ylocal, d2);
  d2len = MathExtra::len3(d2);
  if (d2len == 0.0) return;     // zero length of triangle side
  MathExtra::scale3(1.0 / d2len, d2);

  MathExtra::cross3(d1, d2, normal);
  MathExtra::norm3(normal);
  invndotd = MathExtra::dot3(normal, camDir);

  // triangle parallel to camera and thus invisible
  if (invndotd == 0.0) return;
  invndotd = 1.0 / invndotd;

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

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist : -tanPerPixel / zoom;
  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int>(floor(xf));
  int yc = static_cast<int>(floor(yf));
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  double pixelLeftFull = rasterLeft / pixelWidth;
  double pixelRightFull = rasterRight / pixelWidth;
  double pixelDownFull = rasterDown / pixelWidth;
  double pixelUpFull = rasterUp / pixelWidth;
  int pixelLeft = static_cast<int>(ceil(pixelLeftFull));
  int pixelRight = static_cast<int>(ceil(pixelRightFull));
  int pixelDown = static_cast<int>(ceil(pixelDownFull));
  int pixelUp = static_cast<int>(ceil(pixelUpFull));

  for (int iy = yc - pixelDown; iy <= yc + pixelUp; ++iy) {
    for (int ix = xc - pixelLeft; ix <= xc + pixelRight; ++ix) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;
      if ((opacity < 1.0) && (transthresh[ix % TRANK][iy % TRANK] > opacity))
        continue;

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

      MathExtra::sub3(zlocal, xlocal, s1);
      MathExtra::sub3(ylocal, xlocal, s2);
      MathExtra::sub3(p, xlocal, s3);
      MathExtra::cross3(s1, s2, c1);
      MathExtra::cross3(s1, s3, c2);
      if (MathExtra::dot3(c1, c2) < 0.0) continue;

      MathExtra::sub3(xlocal, ylocal, s1);
      MathExtra::sub3(zlocal, ylocal, s2);
      MathExtra::sub3(p, ylocal, s3);
      MathExtra::cross3(s1, s2, c1);
      MathExtra::cross3(s1, s3, c2);
      if (MathExtra::dot3(c1, c2) < 0.0) continue;

      MathExtra::sub3(ylocal, zlocal, s1);
      MathExtra::sub3(xlocal, zlocal, s2);
      MathExtra::sub3(p, zlocal, s3);
      MathExtra::cross3(s1, s2, c1);
      MathExtra::cross3(s1, s3, c2);
      if (MathExtra::dot3(c1, c2) < 0.0) continue;

      double cNormal[3];
      cNormal[0] = MathExtra::dot3(camRight, normal);
      cNormal[1] = MathExtra::dot3(camUp, normal);
      cNormal[2] = MathExtra::dot3(camDir, normal);

      depth = dist - t;
      draw_pixel(ix,iy,depth,cNormal,surfaceColor);
    }
  }
}

/* ----------------------------------------------------------------------
   draw triangle with 3 corner points x,y,z
   3 normals nx,ny,nz and 3 colors cx,cy,cz, one per corner
   normal and color for each pixel are interpolated using barycentric coords
------------------------------------------------------------------------- */

void Image::draw_trinorm(const double *x, const double *y, const double *z,
                         const double *nx, const double *ny, const double *nz,
                         const double *cx, const double *cy, const double *cz,
                         const double opacity)
{
  if (opacity <= 0.0) return;  // nothing to do

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

  MathExtra::sub3(xlocal, ylocal, d1);
  d1len = MathExtra::len3(d1);
  if (d1len == 0.0) return;     // zero length of triangle side
  MathExtra::scale3(1.0 / d1len, d1);

  MathExtra::sub3(zlocal, ylocal, d2);
  d2len = MathExtra::len3(d2);
  if (d2len == 0.0) return;     // zero length of triangle side
  MathExtra::scale3(1.0 / d2len, d2);

  MathExtra::cross3(d1, d2, normal);
  MathExtra::norm3(normal);
  invndotd = MathExtra::dot3(normal, camDir);

  // triangle parallel to camera and thus invisible
  if (invndotd == 0.0) return;
  invndotd = 1.0 / invndotd;

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

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel * dist : -tanPerPixel / zoom;
  double xf = xmap / pixelWidth;
  double yf = ymap / pixelWidth;
  int xc = static_cast<int>(floor(xf));
  int yc = static_cast<int>(floor(yf));
  double width_error = xf - xc;
  double height_error = yf - yc;

  // shift 0,0 to screen center (vs lower left)

  xc += width / 2;
  yc += height / 2;

  double pixelLeftFull = rasterLeft / pixelWidth;
  double pixelRightFull = rasterRight / pixelWidth;
  double pixelDownFull = rasterDown / pixelWidth;
  double pixelUpFull = rasterUp / pixelWidth;
  int pixelLeft = static_cast<int>(ceil(pixelLeftFull));
  int pixelRight = static_cast<int>(ceil(pixelRightFull));
  int pixelDown = static_cast<int>(ceil(pixelDownFull));
  int pixelUp = static_cast<int>(ceil(pixelUpFull));

  // precompute for barycentric coordinates

  double v0[3], v1[3];
  MathExtra::sub3(ylocal, xlocal, v0);
  MathExtra::sub3(zlocal, xlocal, v1);
  double d00 = MathExtra::dot3(v0, v0);
  double d01 = MathExtra::dot3(v0, v1);
  double d11 = MathExtra::dot3(v1, v1);
  double denom = d00 * d11 - d01 * d01;
  if (denom == 0.0) return;    // degenerate triangle
  double inv_denom = 1.0 / denom;

  for (int iy = yc - pixelDown; iy <= yc + pixelUp; ++iy) {
    for (int ix = xc - pixelLeft; ix <= xc + pixelRight; ++ix) {
      if (iy < 0 || iy >= height || ix < 0 || ix >= width) continue;
      if ((opacity < 1.0) && (transthresh[ix % TRANK][iy % TRANK] > opacity))
        continue;

      double sy = ((iy - yc) - height_error) * pixelWidth;
      double sx = ((ix - xc) - width_error) * pixelWidth;
      surface[0] = camRight[0] * sx + camUp[0] * sy;
      surface[1] = camRight[1] * sx + camUp[1] * sy;
      surface[2] = camRight[2] * sx + camUp[2] * sy;

      double t = -MathExtra::dot3(normal,surface) * invndotd;

      // compute point on triangle plane

      double p[3];
      p[0] = xlocal[0] + surface[0] + camDir[0] * t;
      p[1] = xlocal[1] + surface[1] + camDir[1] * t;
      p[2] = xlocal[2] + surface[2] + camDir[2] * t;

      // compute barycentric coordinates

      double v2[3];
      MathExtra::sub3(p, xlocal, v2);
      double d20 = MathExtra::dot3(v2, v0);
      double d21 = MathExtra::dot3(v2, v1);
      double lambda_y = (d11 * d20 - d01 * d21) * inv_denom;
      double lambda_z = (d00 * d21 - d01 * d20) * inv_denom;
      double lambda_x = 1.0 - lambda_y - lambda_z;

      // point outside triangle if any barycentric coordinate is negative

      if ((lambda_x < 0.0) || (lambda_y < 0.0) || (lambda_z < 0.0)) continue;

      // interpolate normal from per-vertex normals

      double inormal[3];
      inormal[0] = lambda_x * nx[0] + lambda_y * ny[0] + lambda_z * nz[0];
      inormal[1] = lambda_x * nx[1] + lambda_y * ny[1] + lambda_z * nz[1];
      inormal[2] = lambda_x * nx[2] + lambda_y * ny[2] + lambda_z * nz[2];
      MathExtra::norm3(inormal);

      // interpolate color from per-vertex colors

      double icolor[3];
      icolor[0] = lambda_x * cx[0] + lambda_y * cy[0] + lambda_z * cz[0];
      icolor[1] = lambda_x * cx[1] + lambda_y * cy[1] + lambda_z * cz[1];
      icolor[2] = lambda_x * cx[2] + lambda_y * cy[2] + lambda_z * cz[2];

      // transform interpolated normal to camera space

      double cNormal[3];
      cNormal[0] = MathExtra::dot3(camRight, inormal);
      cNormal[1] = MathExtra::dot3(camUp, inormal);
      cNormal[2] = MathExtra::dot3(camDir, inormal);

      depth = dist - t;
      draw_pixel(ix, iy, depth, cNormal, icolor);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Image::draw_pixel(int ix, int iy, double depth, const double *surface,
                       const double *surfaceColor)
{
  if (!std::isfinite(depth)) return; // reject pixels with invalid depth buffer values
  if (!surfaceColor) return;         // reject pixels with an invalid color

  double diffuseKey,diffuseFill,diffuseBack,specularKey;
  if (depth < 0 || (depthBuffer[ix + iy*width] >= 0 && depth >= depthBuffer[ix + iy*width])) return;
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

  double pixelWidth = (tanPerPixel > 0) ? tanPerPixel : -tanPerPixel / zoom;
  int pixelRadius = (int) trunc (SSAORadius / pixelWidth + 0.5);

  // each proc is assigned a subset of contiguous pixels from the full image
  // pixels are contiguous in x (columns within a row), then by row
  // index = pixels from 0 to npixel-1
  // x = column # from 0 to width-1
  // y = row # from 0 to height-1

  int pixelstart = static_cast<int>(1.0*me/nprocs * npixels);
  int pixelstop = static_cast<int>(1.0*(me+1)/nprocs * npixels);

  // fill buffer with random numbers to avoid race conditions
  auto *uniform = new double[pixelstop - pixelstart];
  for (int i = 0; i < pixelstop - pixelstart; ++i) uniform[i] = random->uniform();

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int index = pixelstart; index < pixelstop; index++) {
    int x = index % width;
    int y = index / width;

    double cdepth = depthBuffer[index];
    if (cdepth < 0) continue;

    double sx = surfaceBuffer[index * 2 + 0];
    double sy = surfaceBuffer[index * 2 + 1];
    double sin_t = -sqrt(sx*sx + sy*sy);

    double mytheta = uniform[index - pixelstart] * SSAOJitter;
    double ao = 0.0;

    for (int s = 0; s < SSAOSamples; ++s) {
      double hx = cos(mytheta);
      double hy = sin(mytheta);
      mytheta += delTheta;

      // multiply by z cross surface tangent
      // so that dot (aka cos) works here

      double scaled_sin_t = sin_t * (hx*sy + hy*sx);

      // Bresenham's line algorithm to march over depthBuffer

      int dx = static_cast<int>(hx * pixelRadius);
      int dy = static_cast<int>(hy * pixelRadius);
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
  delete[] uniform;
}

/* ---------------------------------------------------------------------- */

void Image::write_JPG(FILE *fp)
{
#ifdef LAMMPS_JPEG
  if (!fp) return;

  const int aafactor = fsaa ? 2 : 1;
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo,fp);
  cinfo.image_width = width/aafactor;
  cinfo.image_height = height/aafactor;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo,85,TRUE);
  jpeg_start_compress(&cinfo,TRUE);

  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer = (JSAMPROW)
      &writeBuffer[(cinfo.image_height - 1 - cinfo.next_scanline) * 3 * cinfo.image_width];
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
  if (!fp) return;

  const int aafactor = fsaa ? 2 : 1;
  const int pngwidth = width/aafactor;
  const int pngheight = height/aafactor;
  png_structp png_ptr;
  png_infop info_ptr;

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

  if (!png_ptr) return;

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr, nullptr);
    return;
  }

  if (setjmp(png_jmpbuf(png_ptr))) { // NOLINT
    png_destroy_write_struct(&png_ptr, &info_ptr);
    return;
  }

  png_init_io(png_ptr, fp);
  png_set_compression_level(png_ptr,Z_BEST_SPEED);
  png_set_IHDR(png_ptr,info_ptr,pngwidth,pngheight,8,PNG_COLOR_TYPE_RGB,
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

  auto *row_pointers = new png_bytep[pngheight];
  for (int i=0; i < pngheight; ++i)
    row_pointers[i] = (png_bytep) &writeBuffer[(pngheight-i-1)*3*pngwidth];

  png_write_image(png_ptr, row_pointers);
  png_write_end(png_ptr, info_ptr);

  png_destroy_write_struct(&png_ptr, &info_ptr);
  delete[] row_pointers;
#else
  LMP_UNUSED_PARAM(fp);
#endif
}

/* ---------------------------------------------------------------------- */

void Image::write_TGA(FILE *fp, bool compressed)
{
  if (!fp) return;

  const int aafactor = fsaa ? 2 : 1;
  const int tgaheight = height/aafactor;
  const int tgawidth = width/aafactor;

  TGAHeader header;
  memset(&header, 0, sizeof(header));
  header.width[0] = static_cast<unsigned char>(tgawidth & 0xFF);
  header.width[1] = static_cast<unsigned char>((tgawidth & 0xFF00) >> 8);
  header.height[0] = static_cast<unsigned char>(tgaheight & 0xFF);
  header.height[1] = static_cast<unsigned char>((tgaheight & 0xFF00) >> 8);
  header.bitsperpixel = 3 * 8 * sizeof(unsigned char);

  if (compressed) {
    header.datatypecode = 10;    // RLE compressed RGB
    fwrite(&header, sizeof(header), 1, fp);

    unsigned char *pix;
    unsigned char old[3];
    unsigned char len;
    for (int i=0; i < tgaheight; ++i) {
      len = 0;
      for (int j = 0; j < tgawidth; ++j) {
        pix = &writeBuffer[i*3*tgawidth + 3*j];
        if (len == 0) {
          old[0] = pix[0];
          old[1] = pix[1];
          old[2] = pix[2];
        }

        if (memcmp(old, pix, 3) || (len == 127) || (j == tgawidth - 1)) {
          if (j != tgawidth - 1) --len;
          len |= 0x80U;
          fputc(len, fp);
          // TGA stores RGB as BGR
          fputc(old[2], fp);
          fputc(old[1], fp);
          fputc(old[0], fp);
          old[0] = pix[0];
          old[1] = pix[1];
          old[2] = pix[2];
          len = 1;
        } else ++len;
      }
    }
  } else {
    header.datatypecode = 2;    // uncompressed RGB
    fwrite(&header, sizeof(header), 1, fp);

    unsigned char *pix;
    for (int i=0; i < tgaheight; ++i) {
      for (int j = 0; j < tgawidth; ++j) {
        pix = &writeBuffer[i*3*tgawidth + 3*j];
        // TGA stores RGB as BGR
        fputc(pix[2], fp);
        fputc(pix[1], fp);
        fputc(pix[0], fp);
      }
    }
  }
  TGAFooter footer{0, 0, "TRUEVISION-XFILE."};
  fwrite(&footer, sizeof(footer), 1, fp);
}

/* ---------------------------------------------------------------------- */

void Image::write_PPM(FILE *fp)
{
  if (!fp) return;

  const int aafactor = fsaa ? 2 : 1;
  const int ppmheight = height/aafactor;
  const int ppmwidth = width/aafactor;

  fprintf(fp,"P6\n%d %d\n",ppmwidth,ppmheight);
  fprintf(fp,"# CREATOR: dump image\n# SOFTWARE: LAMMPS version %s\n255\n", LAMMPS_VERSION);

  for (int y = ppmheight-1; y >= 0; y--)
    fwrite(&writeBuffer[y*ppmwidth*3],3,ppmwidth,fp);
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
   get min/max bounds of dynamic color map index and return 1 if dynamic
------------------------------------------------------------------------- */

int Image::map_info(int index, double &min, double &max)
{
  return maps[index]->info(min, max);
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

  if (r < 0.0 || r > 1.0) return 1;
  if (g < 0.0 || g > 1.0) return 2;
  if (b < 0.0 || b > 1.0) return 3;

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
   return a pointer to the 3 floating point RGB values or nullptr if didn't find
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
    {69/255.0, 69/255.0, 69/255.0},
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
    if (index > NCOLORS) return nullptr;
    return rgb[index-1];
  }
  if (index < 0) {
    if (-index > ncolors) return nullptr;
    return userrgb[-index-1];
  }

  if (color) {
    if (strcmp(color,"none") == 0) return nullptr;
    for (int i = 0; i < ncolors; i++)
      if (strcmp(color,username[i]) == 0) return userrgb[i];
    for (int i = 0; i < NCOLORS; i++)
      if (strcmp(color,name[i]) == 0) return rgb[i];
  }
  return nullptr;
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
    {0.6, 0.6, 0.8},
    {0.6, 0.6, 0.7},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.6901960784, 0.768627451, 0.8705882353},
    {0.1, 0.7, 0.3},
    {0.95, 0.9, 0.2},
    {0.15, 0.5, 0.1},
    {0.6431372549, 0.6666666667, 0.6784313725},
    {0.8, 0.5, 0.5},
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
  return nullptr;
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

  locurrent = mlo = MINVALUE;
  hicurrent = mhi = MAXVALUE;
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
   return > 0 if any error in args, else return 0
   return value is position of failed arg, i.e. array index+1
------------------------------------------------------------------------- */

int ColorMap::reset(int narg, char **arg)
{
  if (utils::is_double(arg[0])) {
    mlo = NUMERIC;
    mlovalue = utils::numeric(FLERR,arg[0],false,lmp);
  } else if (strcmp(arg[0],"min") == 0) mlo = MINVALUE;
  else return 1;

  if (utils::is_double(arg[1])) {
    mhi = NUMERIC;
    mhivalue = utils::numeric(FLERR,arg[1],false,lmp);
  } else if (strcmp(arg[1],"max") == 0) mhi = MAXVALUE;
  else return 2;

  if ((mlo == NUMERIC) && (mhi == NUMERIC) && (mlovalue >= mhivalue)) return 1;

  if ((mlo == MINVALUE) || (mhi == MAXVALUE)) dynamic = 1;
  else dynamic = 0;

  if (strlen(arg[2]) != 2) return 3;
  if (arg[2][0] == 'c') mstyle = CONTINUOUS;
  else if (arg[2][0] == 'd') mstyle = DISCRETE;
  else if (arg[2][0] == 's') mstyle = SEQUENTIAL;
  else return 3;
  if (arg[2][1] == 'a') mrange = ABSOLUTE;
  else if (arg[2][1] == 'f') mrange = FRACTIONAL;
  else return 3;

  if (mstyle == SEQUENTIAL) {
    mbinsize = utils::numeric(FLERR,arg[3],false,lmp);
    if (mbinsize <= 0.0) return 4;
    mbinsizeinv = 1.0/mbinsize;
  }

  nentry = utils::inumeric(FLERR,arg[4],false,lmp);
  if (nentry < 1) return 5;
  delete [] mentry;
  mentry = new MapEntry[nentry];
  mentry[0].svalue = 0.0;

  int n = 5;
  for (int i = 0; i < nentry; i++) {
    if (mstyle == CONTINUOUS) {
      if (n+2 > narg) return n;
      if (utils::is_double(arg[n])) {
        mentry[i].single = NUMERIC;
        mentry[i].svalue = utils::numeric(FLERR,arg[n],false,lmp);
      } else if (strcmp(arg[n],"min") == 0) mentry[i].single = MINVALUE;
      else if (strcmp(arg[n],"max") == 0) mentry[i].single = MAXVALUE;
      else return n+1;
      mentry[i].color = image->color2rgb(arg[n+1]);
      n += 2;
    } else if (mstyle == DISCRETE) {
      if (n+3 > narg) return n+1;
      if (utils::is_double(arg[n])) {
        mentry[i].lo = NUMERIC;
        mentry[i].lvalue = utils::numeric(FLERR,arg[n],false,lmp);
      } else if (strcmp(arg[n],"min") == 0) mentry[i].lo = MINVALUE;
      else if (strcmp(arg[n],"max") == 0) mentry[i].lo = MAXVALUE;
      else return n+1;
      if (utils::is_double(arg[n+1])) {
        mentry[i].hi = NUMERIC;
        mentry[i].hvalue = utils::numeric(FLERR,arg[n+1],false,lmp);
      } else if (strcmp(arg[n+1],"min") == 0) mentry[i].hi = MINVALUE;
      else if (strcmp(arg[n+1],"max") == 0) mentry[i].hi = MAXVALUE;
      else return n+2;
      mentry[i].color = image->color2rgb(arg[n+2]);
      n += 3;
    } else if (mstyle == SEQUENTIAL) {
      if (n+1 > narg) return n+1;
      mentry[i].color = image->color2rgb(arg[n]);
      n += 1;
    }
    if (mentry[i].color == nullptr) return n;
  }

  if (mstyle == CONTINUOUS) {
    if (nentry < 2) return 5;
    if (mentry[0].single != MINVALUE)
      return 6;
    if (mentry[nentry-1].single != MAXVALUE)
      return 4 + nentry*2;
    for (int i = 2; i < nentry-1; i++)
      if (mentry[i].svalue <= mentry[i-1].svalue) return 4 + 2*(i+1);
  } else if (mstyle == DISCRETE) {
    if (nentry < 1) return 5;
    if (mentry[nentry-1].lo != MINVALUE)
      return 3 + nentry*3;
    if (mentry[nentry-1].hi != MAXVALUE)
      return 4 + nentry*3;
  } else if (mstyle == SEQUENTIAL) {
    if (nentry < 1) return 5;
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

int ColorMap::info(double &min, double &max)
{
  min = locurrent;
  max = hicurrent;
  return dynamic;
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
    int ibin = static_cast<int>((value-lo) * mbinsizeinv);
    return mentry[ibin%nentry].color;
  }

  // always return a non-NULL pointer
  return mentry[0].color;
}

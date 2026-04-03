/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_GRAPHICS_H
#define LMP_GRAPHICS_H

// common definitions and declarations for graphics support in LAMMPS

namespace LAMMPS_NS::Graphics {
enum {
  NONE,
  SPHERE,      // a single sphere with radius provided
  LINE,        // a cylinder with diameter given through fflag2
  TRI,         // a surface mesh as triangles or cylinder mesh based on fflag1, fflag2 sets diameter
  CYLINDER,    // a cylinder with diameter given by fix, fflag1 choose caps, fflag2 adjusts diameter
  TRIANGLE,    // a regular triangle, no settings apply
  TRINORM,     // a curved triangle with surface normals and a color for each corner
  BOND,        // two connected cylinders with bond diameter, colored by atom types, fflag1 sets cap
  ARROW,       // a cylinder with a conical tip and a flat cap at the bottom
  CONE,        // a truncated cone with flat caps, fflag1 sets caps
  PIXMAP       // a pointer to a pixmap buffer at x,y,z location
};    // used by some Body and Fix child classes

// definitions for rendering caps and sides of a truncated cone
enum {
  CONE_TOP = 1 << 0,                            // draw top cap of cone/cylinder
  CONE_BOT = 1 << 1,                            // draw bottom cap of cone/cylinder
  CONE_SIDE = 1 << 2,                           // draw side of cone/cylinder
  CONE_ALL = CONE_TOP | CONE_BOT | CONE_SIDE    // all of the above
};
}    // namespace LAMMPS_NS::Graphics
#endif

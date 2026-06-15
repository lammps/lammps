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

#ifdef FIX_CLASS
// clang-format off
FixStyle(graphics/objects,FixGraphicsObjects);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_OBJECTS_H
#define LMP_FIX_GRAPHICS_OBJECTS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGraphicsObjects : public Fix {
 public:
  FixGraphicsObjects(class LAMMPS *, int, char **);
  ~FixGraphicsObjects() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;

  int image(int *&, double **&) override;

 protected:
  int varflag;
  int numobjs;
  int *imgobjs;
  double **imgparms;

  struct SphereItem {
    int style;
    int type;
    double pos[3];
    double diameter;
    char *xstr, *ystr, *zstr, *dstr;
    int xvar, yvar, zvar, dvar;
  };

  struct CylinderItem {
    int style;
    int type;
    double pos1[3];
    double pos2[3];
    double diameter;
    char *x1str, *y1str, *z1str, *x2str, *y2str, *z2str, *dstr;
    int x1var, y1var, z1var, x2var, y2var, z2var, dvar;
  };

  struct ArrowItem {
    int style;
    int type;
    double tip[3];
    double bot[3];
    double diameter;
    double ratio;
    char *x1str, *y1str, *z1str, *x2str, *y2str, *z2str, *dstr;
    int x1var, y1var, z1var, x2var, y2var, z2var, dvar;
  };

  struct ConeItem {
    int style;
    int type;
    double top[3];
    double bot[3];
    double topdiam;
    double botdiam;
    char *x1str, *y1str, *z1str, *x2str, *y2str, *z2str, *d1str, *d2str;
    int x1var, y1var, z1var, x2var, y2var, z2var, d1var, d2var;
    int sides;
  };

  struct ProgbarItem {
    int style;
    int type1;
    int type2;
    int dim;
    int tics;
    double pos[3];
    double diameter;
    double length;
    double progress;
    char *pstr;
    int pvar;
  };

  union GraphicsItem {
    GraphicsItem() = delete;
    GraphicsItem(const SphereItem &s) : sphere(s) {}
    GraphicsItem(const CylinderItem &c) : cylinder(c) {}
    GraphicsItem(const ArrowItem &a) : arrow(a) {}
    GraphicsItem(const ConeItem &c) : cone(c) {}
    GraphicsItem(const ProgbarItem &p) : progbar(p) {}

    int style;
    SphereItem sphere;
    CylinderItem cylinder;
    ArrowItem arrow;
    ConeItem cone;
    ProgbarItem progbar;
  };

  std::vector<GraphicsItem> items;
};
}    // namespace LAMMPS_NS
#endif
#endif

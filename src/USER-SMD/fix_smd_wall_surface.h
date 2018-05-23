/* -*- c++ -*- ----------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(smd/wall_surface,FixSMDWallSurface)

#else

#ifndef LMP_FIX_SMD_WALL_SURFACE_H
#define LMP_FIX_SMD_WALL_SURFACE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSMDWallSurface: public Fix {

public:
        FixSMDWallSurface(class LAMMPS *, int, char **);
        virtual ~FixSMDWallSurface();
        int setmask();
        void init();
        void setup(int);
        void min_setup(int);

        int count_words(const char *line);
        void read_triangles(int pass);

private:
        int first;                 // flag for first time initialization
        double sublo[3], subhi[3];   // epsilon-extended proc sub-box for adding atoms;
        char *filename;
        int wall_particle_type;
        int wall_molecule_id;
};
}

#endif
#endif


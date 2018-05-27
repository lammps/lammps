/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the USER-SMD package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

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

#ifdef FIX_CLASS

FixStyle(smd/move_tri_surf,FixSMDMoveTriSurf)

#else

#ifndef LMP_FIX_SMD_INTEGRATE_TRIANGULAR_SURFACE_H
#define LMP_FIX_SMD_INTEGRATE_TRIANGULAR_SURFACE_H

#include <Eigen/Eigen>
#include "fix.h"

namespace LAMMPS_NS {

class FixSMDMoveTriSurf: public Fix {
public:
        FixSMDMoveTriSurf(class LAMMPS *, int, char **);
        ~FixSMDMoveTriSurf();
        int setmask();
        virtual void init();
        virtual void initial_integrate(int);
        void reset_dt();
        int pack_forward_comm(int, int *, double *, int, int *);
        void unpack_forward_comm(int, int, double *);

protected:
        double dtv;
        bool linearFlag, wiggleFlag, rotateFlag;
        double vx, vy, vz;
        Eigen::Vector3d rotation_axis, origin;
        double rotation_period;
        Eigen::Matrix3d u_cross, uxu;
        double wiggle_travel, wiggle_max_travel, wiggle_direction;
};

}

#endif
#endif

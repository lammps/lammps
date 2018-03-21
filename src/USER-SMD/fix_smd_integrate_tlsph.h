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

FixStyle(smd/integrate_tlsph,FixSMDIntegrateTlsph)

#else

#ifndef LMP_FIX_SMD_INTEGRATE_TLSPH_H
#define LMP_FIX_SMD_INTEGRATE_TLSPH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSMDIntegrateTlsph: public Fix {
        friend class Neighbor;
        friend class PairTlsph;
public:
    FixSMDIntegrateTlsph(class LAMMPS *, int, char **);
    virtual ~FixSMDIntegrateTlsph() {
    }
    int setmask();
    virtual void init();
    virtual void initial_integrate(int);
    virtual void final_integrate();
    virtual void reset_dt();

protected:
    double dtv, dtf, vlimit, vlimitsq;
    int mass_require;
    bool xsphFlag;

    class Pair *pair;
};

}

#endif
#endif

/* ERROR/WARNING messages:

 E: Illegal ... command

 Self-explanatory.  Check the input script syntax and compare to the
 documentation for the command.  You can use -echo screen as a
 command-line option when running LAMMPS to see the offending line.

 */

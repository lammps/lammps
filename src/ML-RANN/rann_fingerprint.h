/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*  ----------------------------------------------------------------------
   Contributing authors: Christopher Barrett (MSU) barrett@me.msstate.edu
                              Doyl Dickel (MSU) doyl@me.msstate.edu
    ----------------------------------------------------------------------*/
/*
“The research described and the resulting data presented herein, unless
otherwise noted, was funded under PE 0602784A, Project T53 "Military
Engineering Applied Research", Task 002 under Contract No. W56HZV-17-C-0095,
managed by the U.S. Army Combat Capabilities Development Command (CCDC) and
the Engineer Research and Development Center (ERDC).  The work described in
this document was conducted at CAVS, MSU.  Permission was granted by ERDC
to publish this information. Any opinions, findings and conclusions or
recommendations expressed in this material are those of the author(s) and
do not necessarily reflect the views of the United States Army.​”

DISTRIBUTION A. Approved for public release; distribution unlimited. OPSEC#4918

----------------*/

#ifndef LMP_RANN_FINGERPRINT_H
#define LMP_RANN_FINGERPRINT_H

#include <string>
#include <vector>

namespace LAMMPS_NS {
class PairRANN;
namespace RANN {
  class Fingerprint {
   public:
    Fingerprint(PairRANN *);
    virtual ~Fingerprint() {}

    virtual bool parse_values(std::string, std::vector<std::string>) { return false; }
    virtual void write_values(FILE *) {}

    virtual void init(int *, int) {}
    virtual void allocate() {}

    void init_screen(int);

    //no screen,no spin
    virtual void compute_fingerprint(double *, double *, double *, double *, int, int, double *,
                                     double *, double *, int *, int, int *)
    {
    }
    //screen
    virtual void compute_fingerprint(double *, double *, double *, double *, double *, double *,
                                     double *, double *, double *, double *, double *, bool *, int,
                                     int, double *, double *, double *, int *, int, int *)
    {
    }
    //spin
    virtual void compute_fingerprint(double *, double *, double *, double *, double *, double *,
                                     double *, int, int, double *, double *, double *, int *, int,
                                     int *)
    {
    }
    //spin,screen
    virtual void compute_fingerprint(double *, double *, double *, double *, double *, double *,
                                     double *, double *, double *, double *, double *, double *,
                                     double *, double *, bool *, int, int, double *, double *,
                                     double *, int *, int, int *)
    {
    }

    virtual int get_length() { return 0; };
    virtual double cutofffunction(double, double, double);
    virtual void generate_rinvssqrttable();
    bool spin;
    bool screen;
    int n_body_type;    //i-j vs. i-j-k vs. i-j-k-l, etc.
    bool empty;
    bool fullydefined;
    int startingneuron;
    int id;    //based on ordering of fingerprints listed for i-j in potential file
    const char *style;
    int *atomtypes;
    double *rinvsqrttable;
    double rc;
    PairRANN *pair;
  };
}    // namespace RANN
}    // namespace LAMMPS_NS

#endif /* RANN_FINGERPRINT_H_ */

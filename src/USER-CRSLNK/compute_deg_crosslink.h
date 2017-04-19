/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Amin Aramoon, aaramoo1@jhu.edu

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(crosslink,ComputeCrosslink)

#else

#ifndef LMP_COMPUTE_CROSS_H
#define LMP_COMPUTE_CROSS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCrosslink: public Compute {
public:
	ComputeCrosslink(class LAMMPS *, int, char **);
	virtual ~ComputeCrosslink();
	void init() {
	}
	void setup() {
	}
	double compute_scalar();
	int pack_reverse_comm(int, int, double *);
	void unpack_reverse_comm(int, int *, double *);

protected:
	void count_bond();
	int *bondcount;
	int act_type, max_act, crs_type, max_crs, is_invoked;
	int is_crs_control, is_act_control;
	bigint nbond_pre, natoms;
	unsigned int capacity_total, capacity_cur;

	void get_max();

}
;

}

#endif
#endif

/* ERROR/WARNING messages:

 E: Illegal ... command

 Self-explanatory.  Check the input script syntax and compare to the
 documentation for the command.  You can use -echo screen as a
 command-line option when running LAMMPS to see the offending line.

 */

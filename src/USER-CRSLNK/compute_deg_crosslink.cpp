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

#include "mpi.h"
#include "string.h"
#include "atom.h"
#include "compute_deg_crosslink.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "error.h"
#include "neighbor.h"
#include "memory.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCrosslink::ComputeCrosslink(LAMMPS *lmp, int narg, char **arg) :
		Compute(lmp, narg, arg) {
	int iarg = 3;
	if (narg > 3) {
		for (iarg = 3; iarg < narg; iarg += 3) {
			if (strcmp(arg[iarg], "active") == 0) {
				act_type = force->inumeric(FLERR, arg[iarg + 1]);
				max_act = force->inumeric(FLERR, arg[iarg + 2]);
			} else if (strcmp(arg[iarg], "crosslink") == 0) {
				crs_type = force->inumeric(FLERR, arg[iarg + 1]);
				max_crs = force->inumeric(FLERR, arg[iarg + 2]);
			} else
				error->all(FLERR, "Illegal input command in compute crosslink");
		}
	} else
		error->all(FLERR, "Illegal input command in compute crosslink");

	scalar_flag = 1;
	comm_reverse = 1;
	vector_flag = 0;
	extscalar = 0;
	nbond_pre = atom->nbonds;
	capacity_total = 0;
	natoms = atom->natoms;
	is_invoked = 0;
	capacity_cur = 0;
	bondcount = NULL;
	is_crs_control = is_act_control = 0;

}

/* ---------------------------------------------------------------------- */

ComputeCrosslink::~ComputeCrosslink() {
	memory->destroy(bondcount);
}

/* ---------------------------------------------------------------------- */

double ComputeCrosslink::compute_scalar() {
	/*invoked_scalar = update->ntimestep;
	 int i1, i2, i;
	 if (is_invoked == 0 || natoms != atom->natoms) { // first call
	 natoms = atom->natoms;
	 nbond_pre = atom->nbonds;
	 get_max();
	 is_invoked = 1;
	 } else {
	 ntaken += (atom->nbonds - nbond_pre);
	 nbond_pre = atom->nbonds;
	 }
	 double t_1 = (double) ((int) ntaken);
	 double t_2 = (double) navailable;
	 scalar = t_1 / t_2;
	 return scalar;*/

	int nlocal = atom->nlocal;
	int *type = atom->type;

	count_bond();

	unsigned int present_act = 0;
	unsigned int present_crs = 0;

	if (is_invoked) {

		int present = 0;
		if (is_act_control)
			for (int i = 0; i < nlocal; i++) {
				if (type[i] == act_type) {
					present += max_act - bondcount[i];
				}
			}
		else
			for (int i = 0; i < nlocal; i++) {
				if (type[i] == crs_type) {
					present += max_crs - bondcount[i];
				}
			}

		capacity_cur = 0;

		MPI_Allreduce(&present, &capacity_cur, 1, MPI_INT,
		MPI_SUM, world);

		//printf("%d\t%d\t%d\n", present, capacity_cur, capacity_total);

		scalar = ((double) (capacity_total - capacity_cur)) / capacity_total;
		return scalar;

	}

	for (int i = 0; i < nlocal; i++) {
		if (type[i] == crs_type) {
			present_crs += max_crs - bondcount[i];
		}
		if (type[i] == act_type) {
			present_act += max_act - bondcount[i];
		}
	}

	capacity_cur = capacity_total = 0;

	int capacity_act, capacity_crs;
	capacity_act = capacity_crs = 0;

	MPI_Allreduce(&present_act, &capacity_act, 1, MPI_INT,
	MPI_SUM, world);
	MPI_Allreduce(&present_crs, &capacity_crs, 1, MPI_INT,
	MPI_SUM, world);

	if (capacity_act < capacity_crs) {
		is_act_control = 1;
		is_crs_control = 0;
		capacity_cur = capacity_total = capacity_act;
	} else {
		is_act_control = 0;
		is_crs_control = 1;
		capacity_cur = capacity_total = capacity_crs;
	}

	is_invoked = 1;

	//printf("%d\t%d\t%d\n", 0, capacity_cur, capacity_total);

	scalar = ((double) (capacity_total - capacity_cur)) / capacity_total;
	return scalar;

}

/* ---------------------------------------------------------------------- */

void ComputeCrosslink::get_max() {
	/*
	 int *type = atom->type;
	 int nlocal = atom->nlocal;
	 int nghost = atom->nghost;
	 int nall = nlocal + nghost;
	 int *num_bond = new int[nlocal + nghost];
	 for (i = 0; i < nall; ++i)
	 num_bond[i] = 0;
	 int **bondlist = neighbor->bondlist;
	 int nbondlist = neighbor->nbondlist;
	 int newton_bond = force->newton_bond;
	 ntaken = 0;

	 for (i = 0; i < nbondlist; i++) {
	 i1 = bondlist[i][0];
	 i2 = bondlist[i][1];
	 if ((type[i1] == act_type && type[i2] == crs_type)
	 || (type[i2] == act_type && type[i1] == crs_type))
	 ntaken++;
	 else {
	 num_bond[i1]++;
	 num_bond[i2]++;
	 }
	 }

	 int temp = ntaken;
	 MPI_Allreduce(&temp, &ntaken, 1, MPI_LONG,
	 MPI_SUM, world);

	 int taken_active = 0;
	 int taken_crslinking = 0;

	 for (i = 0; i < nlocal; i++) {
	 int t = type[i];
	 int nbnd = num_bond[i];
	 if (t == act_type) {
	 taken_active += (max_act - nbnd);
	 } else if (t == crs_type) {
	 taken_crslinking += (max_crs - nbnd);
	 }
	 }

	 // passed here

	 long potential_act_all = 0;
	 long potential_crs_all = 0;

	 MPI_Allreduce(&taken_active, &potential_act_all, 1, MPI_LONG,
	 MPI_SUM, world);
	 MPI_Allreduce(&taken_crslinking, &potential_crs_all, 1, MPI_LONG,
	 MPI_SUM, world);

	 // passed here

	 navailable =
	 (potential_act_all > potential_crs_all) ?
	 potential_crs_all : potential_act_all;

	 delete[] num_bond;
	 */

}

/* ---------------------------------------------------------------------- */

int ComputeCrosslink::pack_reverse_comm(int n, int first, double *buf) {
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++)
		buf[m++] = ubuf(bondcount[i]).d;
	return m;
}

/* ---------------------------------------------------------------------- */

void ComputeCrosslink::unpack_reverse_comm(int n, int *list, double *buf) {
	int i, j, m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		bondcount[j] += (int) ubuf(buf[m++]).i;
	}
}

/* ---------------------------------------------------------------------- */

void ComputeCrosslink::count_bond() {
	int i, j, m;

	int *num_bond = atom->num_bond;
	int **bond_type = atom->bond_type;
	tagint **bond_atom = atom->bond_atom;
	int nlocal = atom->nlocal;
	int nghost = atom->nghost;
	int nall = nlocal + nghost;
	int newton_bond = force->newton_bond;

	int nmax = atom->nmax;
	if (!bondcount) {
		bondcount = new int[nmax];
	}

	for (i = 0; i < nall; i++)
		bondcount[i] = 0;

	for (i = 0; i < nlocal; i++)
		for (j = 0; j < num_bond[i]; j++) {
			bondcount[i]++;
			if (newton_bond) {
				m = atom->map(bond_atom[i][j]);
				if (m < 0)
					continue;
				bondcount[m]++;
			}
		}

	// if newton_bond is set, need to sum bondcount

	if (newton_bond)
		comm->reverse_comm_compute(this);

}

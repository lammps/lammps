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

#include "stdlib.h"
#include "atom_vec_meso.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecMeso::AtomVecMeso(LAMMPS *lmp, int narg, char **arg) :
	AtomVec(lmp, narg, arg) {
	molecular = 0;
	mass_type = 1;

	comm_x_only = 0; // we communicate not only x forward but also vest ...
	comm_f_only = 0; // we also communicate de and drho in reverse direction
	size_forward = 8; // 3 + rho + e + vest[3], that means we may only communicate 5 in hybrid
	size_reverse = 5; // 3 + drho + de
	size_border = 12; // 6 + rho + e + vest[3] + cv
	size_velocity = 3;
	size_data_atom = 8;
	size_data_vel = 4;
	xcol_data = 6;

	atom->e_flag = 1;
	atom->rho_flag = 1;
	atom->cv_flag = 1;
	atom->vest_flag = 1;
}

/* ----------------------------------------------------------------------
 grow atom arrays
 n = 0 grows arrays by DELTA
 n > 0 allocates arrays to size n
 ------------------------------------------------------------------------- */

void AtomVecMeso::grow(int n) {
	if (n == 0)
		nmax += DELTA;
	else
		nmax = n;
	atom->nmax = nmax;
	if (nmax < 0 || nmax > MAXSMALLINT)
		error->one(FLERR,"Per-processor system is too big");

	tag = memory->grow(atom->tag, nmax, "atom:tag");
	type = memory->grow(atom->type, nmax, "atom:type");
	mask = memory->grow(atom->mask, nmax, "atom:mask");
	image = memory->grow(atom->image, nmax, "atom:image");
	x = memory->grow(atom->x, nmax, 3, "atom:x");
	v = memory->grow(atom->v, nmax, 3, "atom:v");
	f = memory->grow(atom->f, nmax*comm->nthreads, 3, "atom:f");

	rho = memory->grow(atom->rho, nmax, "atom:rho");
	drho = memory->grow(atom->drho, nmax*comm->nthreads, "atom:drho");
	e = memory->grow(atom->e, nmax, "atom:e");
	de = memory->grow(atom->de, nmax*comm->nthreads, "atom:de");
	vest = memory->grow(atom->vest, nmax, 3, "atom:vest");
	cv = memory->grow(atom->cv, nmax, "atom:cv");

	if (atom->nextra_grow)
		for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
			modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
 reset local array ptrs
 ------------------------------------------------------------------------- */

void AtomVecMeso::grow_reset() {
	tag = atom->tag;
	type = atom->type;
	mask = atom->mask;
	image = atom->image;
	x = atom->x;
	v = atom->v;
	f = atom->f;
	rho = atom->rho;
	drho = atom->drho;
	e = atom->e;
	de = atom->de;
	vest = atom->vest;
	cv = atom->cv;
}

/* ---------------------------------------------------------------------- */

void AtomVecMeso::copy(int i, int j, int delflag) {
	//printf("in AtomVecMeso::copy\n");
	tag[j] = tag[i];
	type[j] = type[i];
	mask[j] = mask[i];
	image[j] = image[i];
	x[j][0] = x[i][0];
	x[j][1] = x[i][1];
	x[j][2] = x[i][2];
	v[j][0] = v[i][0];
	v[j][1] = v[i][1];
	v[j][2] = v[i][2];

	rho[j] = rho[i];
	drho[j] = drho[i];
	e[j] = e[i];
	de[j] = de[i];
	cv[j] = cv[i];
	vest[j][0] = vest[i][0];
	vest[j][1] = vest[i][1];
	vest[j][2] = vest[i][2];

	if (atom->nextra_grow)
		for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
			modify->fix[atom->extra_grow[iextra]]->copy_arrays(i, j);
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_comm_hybrid(int n, int *list, double *buf) {
	//printf("in AtomVecMeso::pack_comm_hybrid\n");
	int i, j, m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		buf[m++] = rho[j];
		buf[m++] = e[j];
		buf[m++] = vest[j][0];
		buf[m++] = vest[j][1];
		buf[m++] = vest[j][2];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::unpack_comm_hybrid(int n, int first, double *buf) {
	//printf("in AtomVecMeso::unpack_comm_hybrid\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		rho[i] = buf[m++];
		e[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_border_hybrid(int n, int *list, double *buf) {
	//printf("in AtomVecMeso::pack_border_hybrid\n");
	int i, j, m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		buf[m++] = rho[j];
		buf[m++] = e[j];
		buf[m++] = vest[j][0];
		buf[m++] = vest[j][1];
		buf[m++] = vest[j][2];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::unpack_border_hybrid(int n, int first, double *buf) {
	//printf("in AtomVecMeso::unpack_border_hybrid\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		rho[i] = buf[m++];
		e[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_reverse_hybrid(int n, int first, double *buf) {
	//printf("in AtomVecMeso::pack_reverse_hybrid\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		buf[m++] = drho[i];
		buf[m++] = de[i];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::unpack_reverse_hybrid(int n, int *list, double *buf) {
	//printf("in AtomVecMeso::unpack_reverse_hybrid\n");
	int i, j, m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		drho[j] += buf[m++];
		de[j] += buf[m++];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_comm(int n, int *list, double *buf, int pbc_flag,
		int *pbc) {
	//printf("in AtomVecMeso::pack_comm\n");
	int i, j, m;
	double dx, dy, dz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0];
			buf[m++] = x[j][1];
			buf[m++] = x[j][2];
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	} else {
		if (domain->triclinic == 0) {
			dx = pbc[0] * domain->xprd;
			dy = pbc[1] * domain->yprd;
			dz = pbc[2] * domain->zprd;
		} else {
			dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
			dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
			dz = pbc[2] * domain->zprd;
		}
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0] + dx;
			buf[m++] = x[j][1] + dy;
			buf[m++] = x[j][2] + dz;
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_comm_vel(int n, int *list, double *buf, int pbc_flag,
		int *pbc) {
	//printf("in AtomVecMeso::pack_comm_vel\n");
	int i, j, m;
	double dx, dy, dz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0];
			buf[m++] = x[j][1];
			buf[m++] = x[j][2];
			buf[m++] = v[j][0];
			buf[m++] = v[j][1];
			buf[m++] = v[j][2];
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	} else {
		if (domain->triclinic == 0) {
			dx = pbc[0] * domain->xprd;
			dy = pbc[1] * domain->yprd;
			dz = pbc[2] * domain->zprd;
		} else {
			dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
			dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
			dz = pbc[2] * domain->zprd;
		}
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0] + dx;
			buf[m++] = x[j][1] + dy;
			buf[m++] = x[j][2] + dz;
			buf[m++] = v[j][0];
			buf[m++] = v[j][1];
			buf[m++] = v[j][2];
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecMeso::unpack_comm(int n, int first, double *buf) {
	//printf("in AtomVecMeso::unpack_comm\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		rho[i] = buf[m++];
		e[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

void AtomVecMeso::unpack_comm_vel(int n, int first, double *buf) {
	//printf("in AtomVecMeso::unpack_comm_vel\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		v[i][0] = buf[m++];
		v[i][1] = buf[m++];
		v[i][2] = buf[m++];
		rho[i] = buf[m++];
		e[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_reverse(int n, int first, double *buf) {
	//printf("in AtomVecMeso::pack_reverse\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		buf[m++] = f[i][0];
		buf[m++] = f[i][1];
		buf[m++] = f[i][2];
		buf[m++] = drho[i];
		buf[m++] = de[i];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecMeso::unpack_reverse(int n, int *list, double *buf) {
	//printf("in AtomVecMeso::unpack_reverse\n");
	int i, j, m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		f[j][0] += buf[m++];
		f[j][1] += buf[m++];
		f[j][2] += buf[m++];
		drho[j] += buf[m++];
		de[j] += buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_border(int n, int *list, double *buf, int pbc_flag,
		int *pbc) {
	//printf("in AtomVecMeso::pack_border\n");
	int i, j, m;
	double dx, dy, dz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0];
			buf[m++] = x[j][1];
			buf[m++] = x[j][2];
			buf[m++] = tag[j];
			buf[m++] = type[j];
			buf[m++] = mask[j];
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = cv[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	} else {
		if (domain->triclinic == 0) {
			dx = pbc[0] * domain->xprd;
			dy = pbc[1] * domain->yprd;
			dz = pbc[2] * domain->zprd;
		} else {
			dx = pbc[0];
			dy = pbc[1];
			dz = pbc[2];
		}
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0] + dx;
			buf[m++] = x[j][1] + dy;
			buf[m++] = x[j][2] + dz;
			buf[m++] = tag[j];
			buf[m++] = type[j];
			buf[m++] = mask[j];
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = cv[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::pack_border_vel(int n, int *list, double *buf, int pbc_flag,
		int *pbc) {
	//printf("in AtomVecMeso::pack_border_vel\n");
	int i, j, m;
	double dx, dy, dz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0];
			buf[m++] = x[j][1];
			buf[m++] = x[j][2];
			buf[m++] = tag[j];
			buf[m++] = type[j];
			buf[m++] = mask[j];
			buf[m++] = v[j][0];
			buf[m++] = v[j][1];
			buf[m++] = v[j][2];
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = cv[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	} else {
		if (domain->triclinic == 0) {
			dx = pbc[0] * domain->xprd;
			dy = pbc[1] * domain->yprd;
			dz = pbc[2] * domain->zprd;
		} else {
			dx = pbc[0];
			dy = pbc[1];
			dz = pbc[2];
		}
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0] + dx;
			buf[m++] = x[j][1] + dy;
			buf[m++] = x[j][2] + dz;
			buf[m++] = tag[j];
			buf[m++] = type[j];
			buf[m++] = mask[j];
			buf[m++] = v[j][0];
			buf[m++] = v[j][1];
			buf[m++] = v[j][2];
			buf[m++] = rho[j];
			buf[m++] = e[j];
			buf[m++] = cv[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecMeso::unpack_border(int n, int first, double *buf) {
	//printf("in AtomVecMeso::unpack_border\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		if (i == nmax)
			grow(0);
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		tag[i] = static_cast<int> (buf[m++]);
		type[i] = static_cast<int> (buf[m++]);
		mask[i] = static_cast<int> (buf[m++]);
		rho[i] = buf[m++];
		e[i] = buf[m++];
		cv[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

void AtomVecMeso::unpack_border_vel(int n, int first, double *buf) {
	//printf("in AtomVecMeso::unpack_border_vel\n");
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		if (i == nmax)
			grow(0);
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		tag[i] = static_cast<int> (buf[m++]);
		type[i] = static_cast<int> (buf[m++]);
		mask[i] = static_cast<int> (buf[m++]);
		v[i][0] = buf[m++];
		v[i][1] = buf[m++];
		v[i][2] = buf[m++];
		rho[i] = buf[m++];
		e[i] = buf[m++];
		cv[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
	}
}

/* ----------------------------------------------------------------------
 pack data for atom I for sending to another proc
 xyz must be 1st 3 values, so comm::exchange() can test on them
 ------------------------------------------------------------------------- */

int AtomVecMeso::pack_exchange(int i, double *buf) {
	//printf("in AtomVecMeso::pack_exchange\n");
	int m = 1;
	buf[m++] = x[i][0];
	buf[m++] = x[i][1];
	buf[m++] = x[i][2];
	buf[m++] = v[i][0];
	buf[m++] = v[i][1];
	buf[m++] = v[i][2];
	buf[m++] = tag[i];
	buf[m++] = type[i];
	buf[m++] = mask[i];
	buf[m++] = image[i];
	buf[m++] = rho[i];
	buf[m++] = e[i];
	buf[m++] = cv[i];
	buf[m++] = vest[i][0];
	buf[m++] = vest[i][1];
	buf[m++] = vest[i][2];

	if (atom->nextra_grow)
		for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
			m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i, &buf[m]);

	buf[0] = m;
	return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecMeso::unpack_exchange(double *buf) {
	//printf("in AtomVecMeso::unpack_exchange\n");
	int nlocal = atom->nlocal;
	if (nlocal == nmax)
		grow(0);

	int m = 1;
	x[nlocal][0] = buf[m++];
	x[nlocal][1] = buf[m++];
	x[nlocal][2] = buf[m++];
	v[nlocal][0] = buf[m++];
	v[nlocal][1] = buf[m++];
	v[nlocal][2] = buf[m++];
	tag[nlocal] = static_cast<int> (buf[m++]);
	type[nlocal] = static_cast<int> (buf[m++]);
	mask[nlocal] = static_cast<int> (buf[m++]);
	image[nlocal] = static_cast<int> (buf[m++]);
	rho[nlocal] = buf[m++];
	e[nlocal] = buf[m++];
	cv[nlocal] = buf[m++];
	vest[nlocal][0] = buf[m++];
	vest[nlocal][1] = buf[m++];
	vest[nlocal][2] = buf[m++];

	if (atom->nextra_grow)
		for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
			m += modify->fix[atom->extra_grow[iextra]]-> unpack_exchange(nlocal,
					&buf[m]);

	atom->nlocal++;
	return m;
}

/* ----------------------------------------------------------------------
 size of restart data for all atoms owned by this proc
 include extra data stored by fixes
 ------------------------------------------------------------------------- */

int AtomVecMeso::size_restart() {
	int i;

	int nlocal = atom->nlocal;
	int n = 17 * nlocal; // 11 + rho + e + cv + vest[3]

	if (atom->nextra_restart)
		for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
			for (i = 0; i < nlocal; i++)
				n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

	return n;
}

/* ----------------------------------------------------------------------
 pack atom I's data for restart file including extra quantities
 xyz must be 1st 3 values, so that read_restart can test on them
 molecular types may be negative, but write as positive
 ------------------------------------------------------------------------- */

int AtomVecMeso::pack_restart(int i, double *buf) {
	int m = 1;
	buf[m++] = x[i][0];
	buf[m++] = x[i][1];
	buf[m++] = x[i][2];
	buf[m++] = tag[i];
	buf[m++] = type[i];
	buf[m++] = mask[i];
	buf[m++] = image[i];
	buf[m++] = v[i][0];
	buf[m++] = v[i][1];
	buf[m++] = v[i][2];
	buf[m++] = rho[i];
	buf[m++] = e[i];
	buf[m++] = cv[i];
	buf[m++] = vest[i][0];
	buf[m++] = vest[i][1];
	buf[m++] = vest[i][2];

	if (atom->nextra_restart)
		for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
			m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i, &buf[m]);

	buf[0] = m;
	return m;
}

/* ----------------------------------------------------------------------
 unpack data for one atom from restart file including extra quantities
 ------------------------------------------------------------------------- */

int AtomVecMeso::unpack_restart(double *buf) {
	int nlocal = atom->nlocal;
	if (nlocal == nmax) {
		grow(0);
		if (atom->nextra_store)
			memory->grow(atom->extra, nmax, atom->nextra_store, "atom:extra");
	}

	int m = 1;
	x[nlocal][0] = buf[m++];
	x[nlocal][1] = buf[m++];
	x[nlocal][2] = buf[m++];
	tag[nlocal] = static_cast<int> (buf[m++]);
	type[nlocal] = static_cast<int> (buf[m++]);
	mask[nlocal] = static_cast<int> (buf[m++]);
	image[nlocal] = static_cast<int> (buf[m++]);
	v[nlocal][0] = buf[m++];
	v[nlocal][1] = buf[m++];
	v[nlocal][2] = buf[m++];
	rho[nlocal] = buf[m++];
	e[nlocal] = buf[m++];
	cv[nlocal] = buf[m++];
	vest[nlocal][0] = buf[m++];
	vest[nlocal][1] = buf[m++];
	vest[nlocal][2] = buf[m++];

	double **extra = atom->extra;
	if (atom->nextra_store) {
		int size = static_cast<int> (buf[0]) - m;
		for (int i = 0; i < size; i++)
			extra[nlocal][i] = buf[m++];
	}

	atom->nlocal++;
	return m;
}

/* ----------------------------------------------------------------------
 create one atom of itype at coord
 set other values to defaults
 ------------------------------------------------------------------------- */

void AtomVecMeso::create_atom(int itype, double *coord) {
	int nlocal = atom->nlocal;
	if (nlocal == nmax)
		grow(0);

	tag[nlocal] = 0;
	type[nlocal] = itype;
	x[nlocal][0] = coord[0];
	x[nlocal][1] = coord[1];
	x[nlocal][2] = coord[2];
	mask[nlocal] = 1;
	image[nlocal] = (512 << 20) | (512 << 10) | 512;
	v[nlocal][0] = 0.0;
	v[nlocal][1] = 0.0;
	v[nlocal][2] = 0.0;
	rho[nlocal] = 0.0;
	e[nlocal] = 0.0;
	cv[nlocal] = 1.0;
	vest[nlocal][0] = 0.0;
	vest[nlocal][1] = 0.0;
	vest[nlocal][2] = 0.0;
	de[nlocal] = 0.0;
	drho[nlocal] = 0.0;

	atom->nlocal++;
}

/* ----------------------------------------------------------------------
 unpack one line from Atoms section of data file
 initialize other atom quantities
 ------------------------------------------------------------------------- */

void AtomVecMeso::data_atom(double *coord, int imagetmp, char **values) {
	int nlocal = atom->nlocal;
	if (nlocal == nmax)
		grow(0);

	tag[nlocal] = atoi(values[0]);
	if (tag[nlocal] <= 0)
		error->one(FLERR,"Invalid atom ID in Atoms section of data file");

	type[nlocal] = atoi(values[1]);
	if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
		error->one(FLERR,"Invalid atom type in Atoms section of data file");

	rho[nlocal] = atof(values[2]);
	e[nlocal] = atof(values[3]);
	cv[nlocal] = atof(values[4]);

	x[nlocal][0] = coord[0];
	x[nlocal][1] = coord[1];
	x[nlocal][2] = coord[2];

	//printf("rho=%f, e=%f, cv=%f, x=%f\n", rho[nlocal], e[nlocal], cv[nlocal], x[nlocal][0]);

	image[nlocal] = imagetmp;

	mask[nlocal] = 1;
	v[nlocal][0] = 0.0;
	v[nlocal][1] = 0.0;
	v[nlocal][2] = 0.0;

	vest[nlocal][0] = 0.0;
	vest[nlocal][1] = 0.0;
	vest[nlocal][2] = 0.0;

	de[nlocal] = 0.0;
	drho[nlocal] = 0.0;

	atom->nlocal++;
}

/* ----------------------------------------------------------------------
 unpack hybrid quantities from one line in Atoms section of data file
 initialize other atom quantities for this sub-style
 ------------------------------------------------------------------------- */

int AtomVecMeso::data_atom_hybrid(int nlocal, char **values) {

	rho[nlocal] = atof(values[0]);
	e[nlocal] = atof(values[1]);
	cv[nlocal] = atof(values[2]);

	return 3;
}

/* ----------------------------------------------------------------------
 return # of bytes of allocated memory
 ------------------------------------------------------------------------- */

bigint AtomVecMeso::memory_usage() {
	bigint bytes = 0;

	if (atom->memcheck("tag"))
		bytes += memory->usage(tag, nmax);
	if (atom->memcheck("type"))
		bytes += memory->usage(type, nmax);
	if (atom->memcheck("mask"))
		bytes += memory->usage(mask, nmax);
	if (atom->memcheck("image"))
		bytes += memory->usage(image, nmax);
	if (atom->memcheck("x"))
		bytes += memory->usage(x, nmax, 3);
	if (atom->memcheck("v"))
		bytes += memory->usage(v, nmax, 3);
	if (atom->memcheck("f"))
		bytes += memory->usage(f, nmax*comm->nthreads, 3);
	if (atom->memcheck("rho"))
		bytes += memory->usage(rho, nmax);
	if (atom->memcheck("drho"))
		bytes += memory->usage(drho, nmax*comm->nthreads);
	if (atom->memcheck("e"))
		bytes += memory->usage(e, nmax);
	if (atom->memcheck("de"))
		bytes += memory->usage(de, nmax*comm->nthreads);
	if (atom->memcheck("cv"))
		bytes += memory->usage(cv, nmax);
	if (atom->memcheck("vest"))
		bytes += memory->usage(vest, nmax);

	return bytes;
}

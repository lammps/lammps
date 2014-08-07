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

/* ----------------------------------------------------------------------
   Contributing author: Alexander Stukowski (LLNL), alex@stukowski.com
   see LLNL copyright notice at bottom of file
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 * File history of changes:
 * 25-Oct-10 - AS: First code version.
 * 17-Feb-11 - AS: Several optimizations (introduced MEAM2Body struct).
 * 25-Mar-11 - AS: Fixed calculation of per-atom virial stress.
 * 11-Apr-11 - AS: Adapted code to new memory management of LAMMPS.
 * 24-Sep-11 - AS: Adapted code to new interface of Error::one() function.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_meam_spline.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMEAMSpline::PairMEAMSpline(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;

  Uprime_values = NULL;
  nmax = 0;
  maxNeighbors = 0;
  twoBodyInfo = NULL;

  comm_forward = 1;
  comm_reverse = 0;
}

/* ---------------------------------------------------------------------- */

PairMEAMSpline::~PairMEAMSpline()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;

  delete[] twoBodyInfo;
  memory->destroy(Uprime_values);

  if(allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairMEAMSpline::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag, vflag);
  else evflag = vflag_fdotr =
         eflag_global = vflag_global = eflag_atom = vflag_atom = 0;

  double cutforcesq = cutoff*cutoff;

  // Grow per-atom array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(Uprime_values);
    nmax = atom->nmax;
    memory->create(Uprime_values,nmax,"pair:Uprime");
  }

  double** const x = atom->x;
  double** forces = atom->f;
  int nlocal = atom->nlocal;
  bool newton_pair = force->newton_pair;

  int inum_full = listfull->inum;
  int* ilist_full = listfull->ilist;
  int* numneigh_full = listfull->numneigh;
  int** firstneigh_full = listfull->firstneigh;

  // Determine the maximum number of neighbors a single atom has

  int newMaxNeighbors = 0;
  for(int ii = 0; ii < inum_full; ii++) {
    int jnum = numneigh_full[ilist_full[ii]];
    if(jnum > newMaxNeighbors) newMaxNeighbors = jnum;
  }

  // Allocate array for temporary bond info

  if(newMaxNeighbors > maxNeighbors) {
    maxNeighbors = newMaxNeighbors;
    delete[] twoBodyInfo;
    twoBodyInfo = new MEAM2Body[maxNeighbors];
  }

  // Sum three-body contributions to charge density and
  // compute embedding energies

  for(int ii = 0; ii < inum_full; ii++) {
    int i = ilist_full[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    int* jlist = firstneigh_full[i];
    int jnum = numneigh_full[i];
    double rho_value = 0;
    int numBonds = 0;
    MEAM2Body* nextTwoBodyInfo = twoBodyInfo;

    for(int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      double jdelx = x[j][0] - xtmp;
      double jdely = x[j][1] - ytmp;
      double jdelz = x[j][2] - ztmp;
      double rij_sq = jdelx*jdelx + jdely*jdely + jdelz*jdelz;

      if(rij_sq < cutforcesq) {
        double rij = sqrt(rij_sq);
        double partial_sum = 0;

        nextTwoBodyInfo->tag = j;
        nextTwoBodyInfo->r = rij;
        nextTwoBodyInfo->f = f.eval(rij, nextTwoBodyInfo->fprime);
        nextTwoBodyInfo->del[0] = jdelx / rij;
        nextTwoBodyInfo->del[1] = jdely / rij;
        nextTwoBodyInfo->del[2] = jdelz / rij;

        for(int kk = 0; kk < numBonds; kk++) {
          const MEAM2Body& bondk = twoBodyInfo[kk];
          double cos_theta = (nextTwoBodyInfo->del[0]*bondk.del[0] +
                              nextTwoBodyInfo->del[1]*bondk.del[1] +
                              nextTwoBodyInfo->del[2]*bondk.del[2]);
          partial_sum += bondk.f * g.eval(cos_theta);
        }

        rho_value += nextTwoBodyInfo->f * partial_sum;
        rho_value += rho.eval(rij);

        numBonds++;
        nextTwoBodyInfo++;
      }
    }

    // Compute embedding energy and its derivative

    double Uprime_i;
    double embeddingEnergy = U.eval(rho_value, Uprime_i) - zero_atom_energy;
    Uprime_values[i] = Uprime_i;
    if(eflag) {
      if(eflag_global) eng_vdwl += embeddingEnergy;
      if(eflag_atom) eatom[i] += embeddingEnergy;
    }

    double forces_i[3] = {0, 0, 0};

    // Compute three-body contributions to force

    for(int jj = 0; jj < numBonds; jj++) {
      const MEAM2Body bondj = twoBodyInfo[jj];
      double rij = bondj.r;
      int j = bondj.tag;

      double f_rij_prime = bondj.fprime;
      double f_rij = bondj.f;

      double forces_j[3] = {0, 0, 0};

      MEAM2Body const* bondk = twoBodyInfo;
      for(int kk = 0; kk < jj; kk++, ++bondk) {
        double rik = bondk->r;

        double cos_theta = (bondj.del[0]*bondk->del[0] +
                            bondj.del[1]*bondk->del[1] +
                            bondj.del[2]*bondk->del[2]);
        double g_prime;
        double g_value = g.eval(cos_theta, g_prime);
        double f_rik_prime = bondk->fprime;
        double f_rik = bondk->f;

        double fij = -Uprime_i * g_value * f_rik * f_rij_prime;
        double fik = -Uprime_i * g_value * f_rij * f_rik_prime;

        double prefactor = Uprime_i * f_rij * f_rik * g_prime;
        double prefactor_ij = prefactor / rij;
        double prefactor_ik = prefactor / rik;
        fij += prefactor_ij * cos_theta;
        fik += prefactor_ik * cos_theta;

        double fj[3], fk[3];

        fj[0] = bondj.del[0] * fij - bondk->del[0] * prefactor_ij;
        fj[1] = bondj.del[1] * fij - bondk->del[1] * prefactor_ij;
        fj[2] = bondj.del[2] * fij - bondk->del[2] * prefactor_ij;
        forces_j[0] += fj[0];
        forces_j[1] += fj[1];
        forces_j[2] += fj[2];

        fk[0] = bondk->del[0] * fik - bondj.del[0] * prefactor_ik;
        fk[1] = bondk->del[1] * fik - bondj.del[1] * prefactor_ik;
        fk[2] = bondk->del[2] * fik - bondj.del[2] * prefactor_ik;
        forces_i[0] -= fk[0];
        forces_i[1] -= fk[1];
        forces_i[2] -= fk[2];

        int k = bondk->tag;
        forces[k][0] += fk[0];
        forces[k][1] += fk[1];
        forces[k][2] += fk[2];

        if(evflag) {
          double delta_ij[3];
          double delta_ik[3];
          delta_ij[0] = bondj.del[0] * rij;
          delta_ij[1] = bondj.del[1] * rij;
          delta_ij[2] = bondj.del[2] * rij;
          delta_ik[0] = bondk->del[0] * rik;
          delta_ik[1] = bondk->del[1] * rik;
          delta_ik[2] = bondk->del[2] * rik;
          ev_tally3(i, j, k, 0.0, 0.0, fj, fk, delta_ij, delta_ik);
        }
      }

      forces[i][0] -= forces_j[0];
      forces[i][1] -= forces_j[1];
      forces[i][2] -= forces_j[2];
      forces[j][0] += forces_j[0];
      forces[j][1] += forces_j[1];
      forces[j][2] += forces_j[2];
    }

    forces[i][0] += forces_i[0];
    forces[i][1] += forces_i[1];
    forces[i][2] += forces_i[2];
  }

  // Communicate U'(rho) values

  comm->forward_comm_pair(this);

  int inum_half = listhalf->inum;
  int* ilist_half = listhalf->ilist;
  int* numneigh_half = listhalf->numneigh;
  int** firstneigh_half = listhalf->firstneigh;

  // Compute two-body pair interactions

  for(int ii = 0; ii < inum_half; ii++) {
    int i = ilist_half[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    int* jlist = firstneigh_half[i];
    int jnum = numneigh_half[i];

    for(int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      double jdel[3];
      jdel[0] = x[j][0] - xtmp;
      jdel[1] = x[j][1] - ytmp;
      jdel[2] = x[j][2] - ztmp;
      double rij_sq = jdel[0]*jdel[0] + jdel[1]*jdel[1] + jdel[2]*jdel[2];

      if(rij_sq < cutforcesq) {
        double rij = sqrt(rij_sq);

        double rho_prime;
        rho.eval(rij, rho_prime);
        double fpair = rho_prime * (Uprime_values[i] + Uprime_values[j]);

        double pair_pot_deriv;
        double pair_pot = phi.eval(rij, pair_pot_deriv);
        fpair += pair_pot_deriv;

        // Divide by r_ij to get forces from gradient

        fpair /= rij;

        forces[i][0] += jdel[0]*fpair;
        forces[i][1] += jdel[1]*fpair;
        forces[i][2] += jdel[2]*fpair;
        forces[j][0] -= jdel[0]*fpair;
        forces[j][1] -= jdel[1]*fpair;
        forces[j][2] -= jdel[2]*fpair;
        if (evflag) ev_tally(i, j, nlocal, newton_pair,
                             pair_pot, 0.0, -fpair, jdel[0], jdel[1], jdel[2]);
      }
    }
  }

  if(vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairMEAMSpline::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMEAMSpline::settings(int narg, char **arg)
{
  if(narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMEAMSpline::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // for now, only allow single element

  if (nelements > 1)
    error->all(FLERR,
               "Pair meam/spline only supports single element potentials");

  // read potential file

  read_file(arg[2]);

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

#define MAXLINE 1024

void PairMEAMSpline::read_file(const char* filename)
{
        if(comm->me == 0) {
                FILE *fp = force->open_potential(filename);
                if(fp == NULL) {
                        char str[1024];
                        sprintf(str,"Cannot open spline MEAM potential file %s", filename);
                        error->one(FLERR,str);
                }

                // Skip first line of file.
                char line[MAXLINE];
                fgets(line, MAXLINE, fp);

                // Parse spline functions.
                phi.parse(fp, error);
                rho.parse(fp, error);
                U.parse(fp, error);
                f.parse(fp, error);
                g.parse(fp, error);

                fclose(fp);
        }

        // Transfer spline functions from master processor to all other processors.
        phi.communicate(world, comm->me);
        rho.communicate(world, comm->me);
        f.communicate(world, comm->me);
        U.communicate(world, comm->me);
        g.communicate(world, comm->me);

        // Calculate 'zero-point energy' of single atom in vacuum.
        zero_atom_energy = U.eval(0.0);

        // Determine maximum cutoff radius of all relevant spline functions.
        cutoff = 0.0;
        if(phi.cutoff() > cutoff) cutoff = phi.cutoff();
        if(rho.cutoff() > cutoff) cutoff = rho.cutoff();
        if(f.cutoff() > cutoff) cutoff = f.cutoff();

        // Set LAMMPS pair interaction flags.
        for(int i = 1; i <= atom->ntypes; i++) {
                for(int j = 1; j <= atom->ntypes; j++) {
                        setflag[i][j] = 1;
                        cutsq[i][j] = cutoff;
                }
        }

        //phi.writeGnuplot("phi.gp", "Phi(r)");
        //rho.writeGnuplot("rho.gp", "Rho(r)");
        //f.writeGnuplot("f.gp", "f(r)");
        //U.writeGnuplot("U.gp", "U(rho)");
        //g.writeGnuplot("g.gp", "g(x)");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairMEAMSpline::init_style()
{
        if(force->newton_pair == 0)
                error->all(FLERR,"Pair style meam/spline requires newton pair on");

        // Need both full and half neighbor list.
        int irequest_full = neighbor->request(this);
        neighbor->requests[irequest_full]->id = 1;
        neighbor->requests[irequest_full]->half = 0;
        neighbor->requests[irequest_full]->full = 1;
        int irequest_half = neighbor->request(this);
        neighbor->requests[irequest_half]->id = 2;
        neighbor->requests[irequest_half]->half = 0;
        neighbor->requests[irequest_half]->half_from_full = 1;
        neighbor->requests[irequest_half]->otherlist = irequest_full;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
------------------------------------------------------------------------- */
void PairMEAMSpline::init_list(int id, NeighList *ptr)
{
        if(id == 1) listfull = ptr;
        else if(id == 2) listhalf = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
double PairMEAMSpline::init_one(int i, int j)
{
        return cutoff;
}

/* ---------------------------------------------------------------------- */

int PairMEAMSpline::pack_forward_comm(int n, int *list, double *buf, 
                                      int pbc_flag, int *pbc)
{
        int* list_iter = list;
        int* list_iter_end = list + n;
        while(list_iter != list_iter_end)
                *buf++ = Uprime_values[*list_iter++];
        return n;
}

/* ---------------------------------------------------------------------- */

void PairMEAMSpline::unpack_forward_comm(int n, int first, double *buf)
{
        memcpy(&Uprime_values[first], buf, n * sizeof(buf[0]));
}

/* ---------------------------------------------------------------------- */

int PairMEAMSpline::pack_reverse_comm(int n, int first, double *buf)
{
        return 0;
}

/* ---------------------------------------------------------------------- */

void PairMEAMSpline::unpack_reverse_comm(int n, int *list, double *buf)
{
}

/* ----------------------------------------------------------------------
   Returns memory usage of local atom-based arrays
------------------------------------------------------------------------- */
double PairMEAMSpline::memory_usage()
{
        return nmax * sizeof(double);        // The Uprime_values array.
}


/// Parses the spline knots from a text file.
void PairMEAMSpline::SplineFunction::parse(FILE* fp, Error* error)
{
        char line[MAXLINE];

        // Parse number of spline knots.
        fgets(line, MAXLINE, fp);
        int n = atoi(line);
        if(n < 2)
                error->one(FLERR,"Invalid number of spline knots in MEAM potential file");

        // Parse first derivatives at beginning and end of spline.
        fgets(line, MAXLINE, fp);
        double d0 = atof(strtok(line, " \t\n\r\f"));
        double dN = atof(strtok(NULL, " \t\n\r\f"));
        init(n, d0, dN);

        // Skip line.
        fgets(line, MAXLINE, fp);

        // Parse knot coordinates.
        for(int i=0; i<n; i++) {
                fgets(line, MAXLINE, fp);
                double x, y, y2;
                if(sscanf(line, "%lg %lg %lg", &x, &y, &y2) != 3) {
                        error->one(FLERR,"Invalid knot line in MEAM potential file");
                }
                setKnot(i, x, y);
        }

        prepareSpline(error);
}

/// Calculates the second derivatives at the knots of the cubic spline.
void PairMEAMSpline::SplineFunction::prepareSpline(Error* error)
{
        xmin = X[0];
        xmax = X[N-1];

        isGridSpline = true;
        h = (xmax-xmin)/(N-1);
        hsq = h*h;

        double* u = new double[N];
        Y2[0] = -0.5;
        u[0] = (3.0/(X[1]-X[0])) * ((Y[1]-Y[0])/(X[1]-X[0]) - deriv0);
        for(int i = 1; i <= N-2; i++) {
                double sig = (X[i]-X[i-1]) / (X[i+1]-X[i-1]);
                double p = sig * Y2[i-1] + 2.0;
                Y2[i] = (sig - 1.0) / p;
                u[i] = (Y[i+1]-Y[i]) / (X[i+1]-X[i]) - (Y[i]-Y[i-1])/(X[i]-X[i-1]);
                u[i] = (6.0 * u[i]/(X[i+1]-X[i-1]) - sig*u[i-1])/p;

                if(fabs(h*i+xmin - X[i]) > 1e-8)
                        isGridSpline = false;
        }

        double qn = 0.5;
        double un = (3.0/(X[N-1]-X[N-2])) * (derivN - (Y[N-1]-Y[N-2])/(X[N-1]-X[N-2]));
        Y2[N-1] = (un - qn*u[N-2]) / (qn * Y2[N-2] + 1.0);
        for(int k = N-2; k >= 0; k--) {
                Y2[k] = Y2[k] * Y2[k+1] + u[k];
        }

        delete[] u;

#if !SPLINE_MEAM_SUPPORT_NON_GRID_SPLINES
        if(!isGridSpline)
                error->one(FLERR,"Support for MEAM potentials with non-uniform cubic splines has not been enabled in the MEAM potential code. Set SPLINE_MEAM_SUPPORT_NON_GRID_SPLINES in pair_spline_meam.h to 1 to enable it");
#endif

        // Shift the spline to X=0 to speed up interpolation.
        for(int i = 0; i < N; i++) {
                Xs[i] = X[i] - xmin;
#if !SPLINE_MEAM_SUPPORT_NON_GRID_SPLINES
                if(i < N-1) Ydelta[i] = (Y[i+1]-Y[i])/h;
                Y2[i] /= h*6.0;
#endif
        }
        xmax_shifted = xmax - xmin;
}

/// Broadcasts the spline function parameters to all processors.
void PairMEAMSpline::SplineFunction::communicate(MPI_Comm& world, int me)
{
        MPI_Bcast(&N, 1, MPI_INT, 0, world);
        MPI_Bcast(&deriv0, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&derivN, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&xmin, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&xmax, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&xmax_shifted, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&isGridSpline, 1, MPI_INT, 0, world);
        MPI_Bcast(&h, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&hsq, 1, MPI_DOUBLE, 0, world);
        if(me != 0) {
                X = new double[N];
                Xs = new double[N];
                Y = new double[N];
                Y2 = new double[N];
                Ydelta = new double[N];
        }
        MPI_Bcast(X, N, MPI_DOUBLE, 0, world);
        MPI_Bcast(Xs, N, MPI_DOUBLE, 0, world);
        MPI_Bcast(Y, N, MPI_DOUBLE, 0, world);
        MPI_Bcast(Y2, N, MPI_DOUBLE, 0, world);
        MPI_Bcast(Ydelta, N, MPI_DOUBLE, 0, world);
}

/// Writes a Gnuplot script that plots the spline function.
///
/// This function is for debugging only!
void PairMEAMSpline::SplineFunction::writeGnuplot(const char* filename, const char* title) const
{
        FILE* fp = fopen(filename, "w");
        fprintf(fp, "#!/usr/bin/env gnuplot\n");
        if(title) fprintf(fp, "set title \"%s\"\n", title);
        double tmin = X[0] - (X[N-1] - X[0]) * 0.05;
        double tmax = X[N-1] + (X[N-1] - X[0]) * 0.05;
        double delta = (tmax - tmin) / (N*200);
        fprintf(fp, "set xrange [%f:%f]\n", tmin, tmax);
        fprintf(fp, "plot '-' with lines notitle, '-' with points notitle pt 3 lc 3\n");
        for(double x = tmin; x <= tmax+1e-8; x += delta) {
                double y = eval(x);
                fprintf(fp, "%f %f\n", x, y);
        }
        fprintf(fp, "e\n");
        for(int i = 0; i < N; i++) {
                fprintf(fp, "%f %f\n", X[i], Y[i]);
        }
        fprintf(fp, "e\n");
        fclose(fp);
}

/* ----------------------------------------------------------------------
 * Spline-based Modified Embedded Atom method (MEAM) potential routine.
 *
 * Copyright (2011) Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Alexander Stukowski (<alex@stukowski.com>).
 * LLNL-CODE-525797 All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free
 * Software Foundation) version 2, dated June 1991.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the
 * GNU General Public License for more details.
 *
 * Our Preamble Notice
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the
 * Lawrence Livermore National Laboratory under Contract No.
 * DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process,
 * or services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
------------------------------------------------------------------------- */

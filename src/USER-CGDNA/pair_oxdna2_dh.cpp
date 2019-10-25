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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxdna2_dh.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOxdna2Dh::PairOxdna2Dh(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairOxdna2Dh::~PairOxdna2Dh()
{
  if (allocated) {

    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(qeff_dh_pf);

    memory->destroy(kappa_dh);
    memory->destroy(b_dh);
    memory->destroy(cut_dh_ast);
    memory->destroy(cutsq_dh_ast);
    memory->destroy(cut_dh_c);
    memory->destroy(cutsq_dh_c);

  }
}

/* ----------------------------------------------------------------------
    compute vector COM-sugar-phosphate backbone interaction site in oxDNA2
------------------------------------------------------------------------- */
void PairOxdna2Dh::compute_interaction_sites(double e1[3],
  double e2[3], double r[3])
{
  double d_cs_x=-0.34, d_cs_y=+0.3408;

  r[0] = d_cs_x*e1[0] + d_cs_y*e2[0];
  r[1] = d_cs_x*e1[1] + d_cs_y*e2[1];
  r[2] = d_cs_x*e1[2] + d_cs_y*e2[2];
}

/* ----------------------------------------------------------------------
   compute function for oxDNA pair interactions
   s=sugar-phosphate backbone site, b=base site, st=stacking site
------------------------------------------------------------------------- */

void PairOxdna2Dh::compute(int eflag, int vflag)
{
  double delf[3],delta[3],deltb[3]; // force, torque increment;;
  double rtmp_s[3],delr[3];
  double evdwl,fpair,factor_lj;
  double r,rsq,rinv;
  // vectors COM-backbone sites in lab frame
  double ra_cs[3],rb_cs[3];

  // quaternions and Cartesian unit vectors in lab frame
  double *qa,ax[3],ay[3],az[3];
  double *qb,bx[3],by[3],bz[3];
  double *special_lj = force->special_lj;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *alist,*blist,*numneigh,**firstneigh;

  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  int a,b,ia,ib,anum,bnum,atype,btype;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  anum = list->inum;
  alist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over pair interaction neighbors of my atoms

  for (ia = 0; ia < anum; ia++) {

    a = alist[ia];
    atype = type[a];

    qa=bonus[a].quat;
    MathExtra::q_to_exyz(qa,ax,ay,az);

    // vector COM-backbone site a
    compute_interaction_sites(ax,ay,ra_cs);

    rtmp_s[0] = x[a][0] + ra_cs[0];
    rtmp_s[1] = x[a][1] + ra_cs[1];
    rtmp_s[2] = x[a][2] + ra_cs[2];

    blist = firstneigh[a];
    bnum = numneigh[a];

    for (ib = 0; ib < bnum; ib++) {

      b = blist[ib];
      factor_lj = special_lj[sbmask(b)]; // = 0 for nearest neighbors
      b &= NEIGHMASK;
      btype = type[b];

      qb=bonus[b].quat;
      MathExtra::q_to_exyz(qb,bx,by,bz);

      // vector COM-backbone site b
      compute_interaction_sites(bx,by,rb_cs);

      // vector backbone site b to a
      delr[0] = rtmp_s[0] - x[b][0] - rb_cs[0];
      delr[1] = rtmp_s[1] - x[b][1] - rb_cs[1];
      delr[2] = rtmp_s[2] - x[b][2] - rb_cs[2];
      rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];

      if (rsq <= cutsq_dh_c[atype][btype]) {

        r = sqrt(rsq);
        rinv = 1.0/r;

        if (r <= cut_dh_ast[atype][btype]) {

          fpair = qeff_dh_pf[atype][btype] * exp(-kappa_dh[atype][btype] * r) *
                  (kappa_dh[atype][btype] + rinv) * rinv * rinv;

          if (eflag) {
            evdwl = qeff_dh_pf[atype][btype] * exp(-kappa_dh[atype][btype]*r) * rinv;
          }

        }
        else {

          fpair = 2.0 * b_dh[atype][btype] * (cut_dh_c[atype][btype] - r) * rinv;

          if (eflag) {
            evdwl = b_dh[atype][btype] * (r - cut_dh_c[atype][btype]) * (r - cut_dh_c[atype][btype]);
          }

        }

        // knock out nearest-neighbor interaction between adjacent backbone sites
        fpair *= factor_lj;
        evdwl *= factor_lj;

        delf[0] = delr[0] * fpair;
        delf[1] = delr[1] * fpair;
        delf[2] = delr[2] * fpair;

        // apply force and torque to each of 2 atoms

        if (newton_pair || a < nlocal) {

          f[a][0] += delf[0];
          f[a][1] += delf[1];
          f[a][2] += delf[2];

          MathExtra::cross3(ra_cs,delf,delta);

          torque[a][0] += delta[0];
          torque[a][1] += delta[1];
          torque[a][2] += delta[2];

        }

        if (newton_pair || b < nlocal) {

          f[b][0] -= delf[0];
          f[b][1] -= delf[1];
          f[b][2] -= delf[2];

          MathExtra::cross3(rb_cs,delf,deltb);

          torque[b][0] -= deltb[0];
          torque[b][1] -= deltb[1];
          torque[b][2] -= deltb[2];

        }


        // increment energy and virial
        // NOTE: The virial is calculated on the 'molecular' basis.
        // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

        if (evflag) ev_tally_xyz(a,b,nlocal,newton_pair,evdwl,0.0,
            delf[0],delf[1],delf[2],x[a][0]-x[b][0],x[a][1]-x[b][1],x[a][2]-x[b][2]);

      }

    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOxdna2Dh::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(kappa_dh,n+1,n+1,"pair:kappa_dh");
  memory->create(qeff_dh_pf,n+1,n+1,"pair:qeff_dh_pf");

  memory->create(b_dh,n+1,n+1,"pair:b_dh");
  memory->create(cut_dh_ast,n+1,n+1,"pair:cut_dh_ast");
  memory->create(cutsq_dh_ast,n+1,n+1,"pair:cutsq_dh_ast");
  memory->create(cut_dh_c,n+1,n+1,"pair:cut_dh_c");
  memory->create(cutsq_dh_c,n+1,n+1,"pair:cutsq_dh_c");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOxdna2Dh::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdna2Dh::coeff(int narg, char **arg)
{
  int count;

  if (narg != 5) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/dh");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  count = 0;

  double T, rhos_dh_one, qeff_dh_one;

  T = force->numeric(FLERR,arg[2]);
  rhos_dh_one = force->numeric(FLERR,arg[3]);
  qeff_dh_one  = force->numeric(FLERR,arg[4]);

  double lambda_dh_one, kappa_dh_one, qeff_dh_pf_one;
  double b_dh_one, cut_dh_ast_one, cut_dh_c_one;

  // Debye length and inverse Debye length

  /*
  NOTE:
    The numerical factor is the Debye length in s.u.
    lambda(T = 300 K = 0.1) =
    sqrt(eps_0 * eps_r * k_B * T/(2 * N_A * e^2 * 1000 mol/m^3))
          * 1/oxDNA_energy_unit
    (see B. Snodin et al., J. Chem. Phys. 142, 234901 (2015).)

  We use
    eps_0 = vacuum permittivity = 8.854187817e-12 F/m
    eps_r = relative permittivity of water = 80
    k_B = Boltzmann constant = 1.3806485279e-23 J/K
    T = absolute temperature = 300 K
    N_A = Avogadro constant = 6.02214085774e23 / mol
    e = elementary charge = 1.6021766208e-19 C
    oxDNA_length_unit = 8.518e-10 m
  */

  lambda_dh_one = 0.3616455075438555*sqrt(T/0.1/rhos_dh_one);
  kappa_dh_one = 1.0/lambda_dh_one;

  // prefactor in DH interaction containing qeff^2

  /*
    NOTE:
      The numerical factor is
      qeff_dh_pf = e^2/(4 * pi * eps_0 * eps_r)
                    * 1/(oxDNA_energy_unit * oxDNA_length_unit)
      (see B. Snodin et al., J. Chem. Phys. 142, 234901 (2015).)

    In addition to the above units we use
      oxDNA_energy_unit = 4.142e-20 J
  */

  qeff_dh_pf_one = 0.08173808693529228*qeff_dh_one*qeff_dh_one;

  // smoothing parameters - determined through continuity and differentiability

  cut_dh_ast_one = 3.0*lambda_dh_one;

  b_dh_one = -(exp(-cut_dh_ast_one/lambda_dh_one) * qeff_dh_pf_one * qeff_dh_pf_one *
      (cut_dh_ast_one + lambda_dh_one) * (cut_dh_ast_one + lambda_dh_one))/
      (-4.0 * cut_dh_ast_one * cut_dh_ast_one * cut_dh_ast_one *
      lambda_dh_one * lambda_dh_one * qeff_dh_pf_one);

  cut_dh_c_one =  cut_dh_ast_one * (qeff_dh_pf_one*cut_dh_ast_one +
      3.0*qeff_dh_pf_one * lambda_dh_one)/
      (qeff_dh_pf_one * (cut_dh_ast_one+lambda_dh_one));

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {

      kappa_dh[i][j] = kappa_dh_one;
      qeff_dh_pf[i][j] = qeff_dh_pf_one;
      b_dh[i][j] = b_dh_one;
      cut_dh_ast[i][j] = cut_dh_ast_one;
      cut_dh_c[i][j] = cut_dh_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/dh");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairOxdna2Dh::init_style()
{
  int irequest;

  // request regular neighbor lists

  irequest = neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use regular
------------------------------------------------------------------------- */

void PairOxdna2Dh::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  if (id  > 0) error->all(FLERR,"Respa not supported");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOxdna2Dh::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    error->all(FLERR,"Coefficient mixing not defined in oxDNA");
  }
  if (offset_flag) {
    error->all(FLERR,"Offset not supported in oxDNA");
  }

  kappa_dh[j][i] = kappa_dh[i][j];
  qeff_dh_pf[j][i] = qeff_dh_pf[i][j];

  b_dh[j][i] = b_dh[i][j];

  cut_dh_ast[j][i] = cut_dh_ast[i][j];
  cut_dh_c[j][i] = cut_dh_c[i][j];

  cutsq_dh_ast[i][j] = cut_dh_ast[i][j]*cut_dh_ast[i][j];
  cutsq_dh_ast[j][i] = cutsq_dh_ast[i][j];

  cutsq_dh_c[i][j] = cut_dh_c[i][j]*cut_dh_c[i][j];
  cutsq_dh_c[j][i] = cutsq_dh_c[i][j];

  // set the master list distance cutoff
  return cut_dh_c[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdna2Dh::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {

        fwrite(&kappa_dh[i][j],sizeof(double),1,fp);
        fwrite(&qeff_dh_pf[i][j],sizeof(double),1,fp);
        fwrite(&b_dh[i][j],sizeof(double),1,fp);
        fwrite(&cut_dh_ast[i][j],sizeof(double),1,fp);
        fwrite(&cut_dh_c[i][j],sizeof(double),1,fp);

    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdna2Dh::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {

          utils::sfread(FLERR,&kappa_dh[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&qeff_dh_pf[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&b_dh[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut_dh_ast[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut_dh_c[i][j],sizeof(double),1,fp,NULL,error);

        }

        MPI_Bcast(&kappa_dh[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&qeff_dh_pf[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b_dh[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_dh_ast[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_dh_c[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOxdna2Dh::write_restart_settings(FILE *fp)
{
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOxdna2Dh::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOxdna2Dh::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\
         %g %g\
         %g %g %g\
         \n",i,
        kappa_dh[i][i],qeff_dh_pf[i][i],
        b_dh[i][i],cut_dh_ast[i][i],cut_dh_c[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOxdna2Dh::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d\
         %g %g\
         %g %g %g\
         \n",i,j,
        kappa_dh[i][j],qeff_dh_pf[i][j],
        b_dh[i][j],cut_dh_ast[i][j],cut_dh_c[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairOxdna2Dh::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str,"kappa_dh") == 0) return (void *) kappa_dh;
  if (strcmp(str,"qeff_dh_pf") == 0) return (void *) qeff_dh_pf;
  if (strcmp(str,"b_dh") == 0) return (void *) b_dh;
  if (strcmp(str,"cut_dh_ast") == 0) return (void *) cut_dh_ast;
  if (strcmp(str,"cut_dh_c") == 0) return (void *) cut_dh_c;

  return NULL;
}

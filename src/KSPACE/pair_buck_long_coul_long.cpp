// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Pieter J. in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "pair_buck_long_coul_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "info.h"
#include "kspace.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairBuckLongCoulLong::PairBuckLongCoulLong(LAMMPS *lmp) :
    Pair(lmp), cut_buck(nullptr), cut_buck_read(nullptr), cut_bucksq(nullptr),
    buck_a_read(nullptr), buck_a(nullptr), buck_c_read(nullptr), buck_c(nullptr), buck1(nullptr),
    buck2(nullptr), buck_rho_read(nullptr), buck_rho(nullptr), rhoinv(nullptr), offset(nullptr)
{
  dispersionflag = ewaldflag = pppmflag = 1;
  respa_enable = 1;
  writedata = 1;
  ftable = nullptr;
  fdisptable = nullptr;
  cut_respa = nullptr;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::options(char **arg, int mask)
{
  if (!*arg) error->all(FLERR,"Illegal pair_style buck/long/coul/long command");

  // "long" : treat this interaction with the long-range (k-space) solver
  // "cut"  : plain real-space cutoff (default, leaves both flags clear)
  // "off"  : interaction disabled

  if (strcmp(*arg,"long") == 0) ewald_order |= mask;
  else if (strcmp(*arg,"off") == 0) ewald_off |= mask;
  else if (strcmp(*arg,"cut") != 0)
    error->all(FLERR,"Illegal pair_style buck/long/coul/long command");
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLong::settings(int narg, char **arg)
{
  if (narg != 3 && narg != 4) error->all(FLERR,"Illegal pair_style command");

  ewald_order = 0;
  ewald_off = 0;

  options(arg, EWALD_DISP);
  options(++arg, EWALD_COUL);

  if (!comm->me && ewald_order == (EWALD_COUL | EWALD_DISP))
    error->warning(FLERR,"Using largest cutoff for buck/long/coul/long");
  if (!*(++arg))
    error->all(FLERR,"Cutoffs missing in pair_style buck/long/coul/long");
  if (!((ewald_order^ewald_off) & EWALD_DISP))
    dispersionflag = 0;
  if (ewald_off & EWALD_DISP)
    error->all(FLERR,"LJ6 off not supported in pair_style buck/long/coul/long");
  if (!((ewald_order^ewald_off) & EWALD_COUL))
    error->all(FLERR,
               "Coulomb cut not supported in pair_style buck/long/coul/coul");
  cut_buck_global = utils::numeric(FLERR,*(arg++),false,lmp);
  if (narg == 4 && ((ewald_order & (EWALD_COUL|EWALD_DISP)) == (EWALD_COUL|EWALD_DISP)))
    error->all(FLERR,"Only one cutoff allowed when requesting all long");
  if (narg == 4) cut_coul = utils::numeric(FLERR,*arg,false,lmp);
  else cut_coul = cut_buck_global;

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_buck[i][j] = cut_buck_global;
  }
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairBuckLongCoulLong::~PairBuckLongCoulLong()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_buck_read);
    memory->destroy(cut_buck);
    memory->destroy(cut_bucksq);
    memory->destroy(buck_a_read);
    memory->destroy(buck_a);
    memory->destroy(buck_c_read);
    memory->destroy(buck_c);
    memory->destroy(buck_rho_read);
    memory->destroy(buck_rho);
    memory->destroy(buck1);
    memory->destroy(buck2);
    memory->destroy(rhoinv);
    memory->destroy(offset);
  }
  if (ftable) free_tables();
  if (fdisptable) free_disp_tables();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_buck_read,n+1,n+1,"pair:cut_buck_read");
  memory->create(cut_buck,n+1,n+1,"pair:cut_buck");
  memory->create(cut_bucksq,n+1,n+1,"pair:cut_bucksq");
  memory->create(buck_a_read,n+1,n+1,"pair:buck_a_read");
  memory->create(buck_a,n+1,n+1,"pair:buck_a");
  memory->create(buck_c_read,n+1,n+1,"pair:buck_c_read");
  memory->create(buck_c,n+1,n+1,"pair:buck_c");
  memory->create(buck_rho_read,n+1,n+1,"pair:buck_rho_read");
  memory->create(buck_rho,n+1,n+1,"pair:buck_rho");
  memory->create(buck1,n+1,n+1,"pair:buck1");
  memory->create(buck2,n+1,n+1,"pair:buck2");
  memory->create(rhoinv,n+1,n+1,"pair:rhoinv");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   extract protected data from object
------------------------------------------------------------------------- */

void *PairBuckLongCoulLong::extract(const char *id, int &dim)
{
  // per-type-pair arrays are 2d, scalars are 0d

  dim = 2;
  if (strcmp(id,"B") == 0) return (void *) buck_c;

  dim = 0;
  if (strcmp(id,"ewald_order") == 0) return (void *) &ewald_order;
  if (strcmp(id,"ewald_cut") == 0) return (void *) &cut_coul;
  if (strcmp(id,"ewald_mix") == 0) return (void *) &mix_flag;
  if (strcmp(id,"cut_coul") == 0) return (void *) &cut_coul;
  if (strcmp(id,"cut_vdwl") == 0) return (void *) &cut_buck_global;
  // "cut_LJ" is a deprecated alias for "cut_vdwl", kept for backward
  // compatibility; remove after a suitable deprecation period
  if (strcmp(id,"cut_LJ") == 0) return (void *) &cut_buck_global;
  return nullptr;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,*(arg++),1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,*(arg++),1,atom->ntypes,jlo,jhi,error);

  double buck_a_one = utils::numeric(FLERR,*(arg++),false,lmp);
  double buck_rho_one = utils::numeric(FLERR,*(arg++),false,lmp);
  double buck_c_one = utils::numeric(FLERR,*(arg++),false,lmp);

  double cut_buck_one = cut_buck_global;
  if (narg == 6) cut_buck_one = utils::numeric(FLERR,*(arg++),false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      buck_a_read[i][j] = buck_a_one;
      buck_c_read[i][j] = buck_c_one;
      buck_rho_read[i][j] = buck_rho_one;
      cut_buck_read[i][j] = cut_buck_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::init_style()
{
  // require an atom style with charge defined

  if (!atom->q_flag && (ewald_order & EWALD_COUL))
    error->all(FLERR,
               "Invoking coulombic in pair style buck/long/coul/long "
               "requires atom attribute q");

  // ensure use of KSpace long-range solver, set two g_ewalds

  if (force->kspace == nullptr)
    error->all(FLERR,"Pair style requires a KSpace style");
  if (ewald_order & EWALD_COUL) g_ewald = force->kspace->g_ewald;
  if (ewald_order & EWALD_DISP) g_ewald_6 = force->kspace->g_ewald_6;

  // set rRESPA cutoffs

  if (utils::strmatch(update->integrate_style,"^respa") &&
      (dynamic_cast<Respa *>(update->integrate))->level_inner >= 0)
    cut_respa = (dynamic_cast<Respa *>(update->integrate))->cutoff;
  else cut_respa = nullptr;

  // setup force tables

  if (ncoultablebits && (ewald_order & EWALD_COUL)) init_tables(cut_coul,cut_respa);
  if (ndisptablebits && (ewald_order & EWALD_DISP)) init_tables_disp(cut_buck_global);

  // request regular or rRESPA neighbor lists if neighrequest_flag != 0

  if (force->kspace->neighrequest_flag) {
    int list_style = NeighConst::REQ_DEFAULT;

    if (update->whichflag == 1 && utils::strmatch(update->integrate_style, "^respa")) {
      auto *respa = dynamic_cast<Respa *>(update->integrate);
      if (respa->level_inner >= 0) list_style = NeighConst::REQ_RESPA_INOUT;
      if (respa->level_middle >= 0) list_style = NeighConst::REQ_RESPA_ALL;
    }
    neighbor->add_request(this, list_style);
  }

  cut_coulsq = cut_coul * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBuckLongCoulLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "All pair coeffs are not set. Status\n" + Info::get_pair_coeff_status(lmp));

  if (ewald_order & EWALD_DISP) cut_buck[i][j] = cut_buck_global;
  else cut_buck[i][j] = cut_buck_read[i][j];
  buck_a[i][j] = buck_a_read[i][j];
  buck_c[i][j] = buck_c_read[i][j];
  buck_rho[i][j] = buck_rho_read[i][j];

  double cut = MAX(cut_buck[i][j],cut_coul);
  cutsq[i][j] = cut*cut;
  cut_bucksq[i][j] = cut_buck[i][j] * cut_buck[i][j];

  buck1[i][j] = buck_a[i][j]/buck_rho[i][j];
  buck2[i][j] = 6.0*buck_c[i][j];
  rhoinv[i][j] = 1.0/buck_rho[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_buck[i][j],cut_coul) < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  if (offset_flag && (cut_buck[i][j] > 0.0)) {
    double rexp = exp(-cut_buck[i][j]/buck_rho[i][j]);
    offset[i][j] = buck_a[i][j]*rexp - buck_c[i][j]/pow(cut_buck[i][j],6.0);
  } else offset[i][j] = 0.0;

  cutsq[j][i] = cutsq[i][j];
  cut_bucksq[j][i] = cut_bucksq[i][j];
  buck_a[j][i] = buck_a[i][j];
  buck_c[j][i] = buck_c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&buck_a_read[i][j],sizeof(double),1,fp);
        fwrite(&buck_rho_read[i][j],sizeof(double),1,fp);
        fwrite(&buck_c_read[i][j],sizeof(double),1,fp);
        fwrite(&cut_buck_read[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&buck_a_read[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&buck_rho_read[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&buck_c_read[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_buck_read[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&buck_a_read[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&buck_rho_read[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&buck_c_read[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_buck_read[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_buck_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
  fwrite(&ewald_order,sizeof(int),1,fp);
  fwrite(&dispersionflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_buck_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&ncoultablebits,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tabinner,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&ewald_order,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&dispersionflag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_buck_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&ewald_order,1,MPI_INT,0,world);
  MPI_Bcast(&dispersionflag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,
            buck_a_read[i][i],buck_rho_read[i][i],buck_c_read[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (ewald_order & EWALD_DISP) {
      fprintf(fp,"%d %d %g %g\n",i,j,
              buck_a_read[i][j],buck_rho_read[i][j]);
      } else {
        fprintf(fp,"%d %d %g %g %g\n",i,j,
                buck_a_read[i][j],buck_rho_read[i][j],buck_c_read[i][j]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute pair interactions
------------------------------------------------------------------------- */

void PairBuckLongCoulLong::compute(int eflag, int vflag)
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, ii, jj, typei, typej, ni;
  int order1 = ewald_order & EWALD_COUL, order6 = ewald_order & EWALD_DISP;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double qi = 0.0, qri = 0.0, *cutsqi, *cut_bucksqi,
         *buck1i, *buck2i, *buckai, *buckci, *rhoinvi, *offseti;
  double r, rsq, r2inv, force_coul, force_buck;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double xi[3], d[3], fi[3];

  for (ii = 0; ii < inum; ii++) {                          // loop over my atoms
    i = ilist[ii];
    if (order1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    typei = type[i];
    offseti = offset[typei];
    buck1i = buck1[typei];
    buck2i = buck2[typei];
    buckai = buck_a[typei];
    buckci = buck_c[typei];
    rhoinvi = rhoinv[typei];
    cutsqi = cutsq[typei];
    cut_bucksqi = cut_bucksq[typei];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    fi[0] = fi[1] = fi[2] = 0.0;
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                          // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq)) {                // coulombic
        if (!ncoultablebits || rsq <= tabinnersq) {        // series real space
          double grij = g_ewald*r;
          double t = 1.0/(1.0+EWALD_P*grij);
          double erfc_poly = ((((t*A5+A4)*t+A3)*t+A2)*t+A1);
          double pre = qri*q[j];                            // qqrd2e * qi * qj
          if (ni == 0) {
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre;
            if (eflag) ecoul = t;
          } else {                                        // special case
            double adjust = pre*(1.0-special_coul[ni])/r;
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre-adjust;
            if (eflag) ecoul = t-adjust;
          }
        } else {                                             // table real space
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          const int k = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
          double fraction = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]);
            if (eflag) ecoul = qiqj*(etable[k]+fraction*detable[k]);
          } else {                                        // special case
            rsq_lookup.f = (1.0-special_coul[ni])*(ctable[k]+fraction*dctable[k]);
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]-(double)rsq_lookup.f);
            if (eflag) ecoul = qiqj*(etable[k]+fraction*detable[k]-(double)rsq_lookup.f);
          }
        }
      } else force_coul = ecoul = 0.0;

      if (rsq < cut_bucksqi[typej]) {                        // buckingham
        double r6inv = r2inv*r2inv*r2inv;
        double expr = exp(-r*rhoinvi[typej]);                // A exp(-r/rho) factor
        if (order6) {                                        // long-range
          if (!ndisptablebits || rsq <= tabinnerdispsq) {
            double gr2 = g2*rsq, a2 = 1.0/gr2;
            double expterm = a2*exp(-gr2)*buckci[typej];      // damped 1/r^6 reciprocal term
            double g6term = g6*((a2+1.0)*a2+0.5)*expterm;
            double g8term = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*expterm*rsq;
            if (ni == 0) {
              force_buck = r*expr*buck1i[typej]-g8term;
              if (eflag) evdwl = expr*buckai[typej]-g6term;
            } else {                                        // special case
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_buck = factor*r*expr*buck1i[typej]-g8term+t*buck2i[typej];
              if (eflag) evdwl = factor*expr*buckai[typej]-g6term+t*buckci[typej];
            }
          } else {                                              //table real space
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            const int disp_k = (rsq_lookup.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double ftable_disp = fdisptable[disp_k]+f_disp*dfdisptable[disp_k];
            double etable_disp = edisptable[disp_k]+f_disp*dedisptable[disp_k];
            if (ni == 0) {
              force_buck = r*expr*buck1i[typej]-ftable_disp*buckci[typej];
              if (eflag) evdwl = expr*buckai[typej]-etable_disp*buckci[typej];
            } else {                                             //special case
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_buck = factor*r*expr*buck1i[typej]-ftable_disp*buckci[typej]+t*buck2i[typej];
              if (eflag) evdwl = factor*expr*buckai[typej]-etable_disp*buckci[typej]+t*buckci[typej];
            }
          }
        } else {                                                // cut
          if (ni == 0) {
            force_buck = r*expr*buck1i[typej]-r6inv*buck2i[typej];
            if (eflag) evdwl = expr*buckai[typej]-r6inv*buckci[typej]-offseti[typej];
          } else {                                        // special case
            double factor = special_lj[ni];
            force_buck = factor*(r*expr*buck1i[typej]-r6inv*buck2i[typej]);
            if (eflag)
              evdwl = factor*(expr*buckai[typej]-r6inv*buckci[typej]-offseti[typej]);
          }
        }
      } else force_buck = evdwl = 0.0;

      fpair = (force_coul+force_buck)*r2inv;

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      fi[0] += fdx;                        // force update
      fi[1] += fdy;
      fi[2] += fdz;
      if (newton_pair || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,ecoul,fpair,d[0],d[1],d[2]);
    }
    f[i][0] += fi[0];
    f[i][1] += fi[1];
    f[i][2] += fi[2];
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLong::compute_inner()
{
  double r, rsq, r2inv, force_coul = 0.0, force_buck, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];

  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  int i, j, ii, jj, typei, typej, ni;
  int order1 = (ewald_order | ~ewald_off) & EWALD_COUL;
  int inum = list->inum_inner;
  int *ilist = list->ilist_inner;
  int *numneigh = list->numneigh_inner;
  int **firstneigh = list->firstneigh_inner;
  double qri, *cut_bucksqi, *buck1i, *buck2i, *rhoinvi;
  double xi[3], d[3], fi[3];

  for (ii = 0; ii < inum; ii++) {                          // loop over my atoms
    i = ilist[ii];
    if (order1) qri = qqrd2e*q[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    fi[0] = fi[1] = fi[2] = 0.0;
    typei = type[i];
    cut_bucksqi = cut_bucksq[typei];
    buck1i = buck1[typei];
    buck2i = buck2[typei];
    rhoinvi = rhoinv[typei];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                          // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      if (rsq >= cut_out_off_sq) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq))                        // coulombic
        force_coul = ni == 0 ?
          qri*q[j]/r : qri*q[j]/r*special_coul[ni];

      typej = type[j];
      if (rsq < cut_bucksqi[typej]) {                          // buckingham
        double r6inv = r2inv*r2inv*r2inv,
                        expr = exp(-r*rhoinvi[typej]);
        force_buck = ni == 0 ?
          (r*expr*buck1i[typej]-r6inv*buck2i[typej]) :
          (r*expr*buck1i[typej]-r6inv*buck2i[typej])*special_lj[ni];
      } else force_buck = 0.0;

      fpair = (force_coul + force_buck) * r2inv;

      if (rsq > cut_out_on_sq) {                        // switching
        double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
        fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      fi[0] += fdx;                       // force update
      fi[1] += fdy;
      fi[2] += fdz;
      if (newton_pair || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }
    }
    f[i][0] += fi[0];
    f[i][1] += fi[1];
    f[i][2] += fi[2];
  }
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLong::compute_middle()
{
  double r, rsq, r2inv, force_coul = 0.0, force_buck, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  int i, j, ii, jj, typei, typej, ni;
  int order1 = (ewald_order | ~ewald_off) & EWALD_COUL;
  int inum = list->inum_middle;
  int *ilist = list->ilist_middle;
  int *numneigh = list->numneigh_middle;
  int **firstneigh = list->firstneigh_middle;
  double qri, *cut_bucksqi, *buck1i, *buck2i, *rhoinvi;
  double xi[3], d[3], fi[3];

  for (ii = 0; ii < inum; ii++) {                          // loop over my atoms
    i = ilist[ii];
    if (order1) qri = qqrd2e*q[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    fi[0] = fi[1] = fi[2] = 0.0;
    typei = type[i];
    cut_bucksqi = cut_bucksq[typei];
    buck1i = buck1[typei];
    buck2i = buck2[typei];
    rhoinvi = rhoinv[typei];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                          // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      if (rsq >= cut_out_off_sq) continue;
      if (rsq <= cut_in_off_sq) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq))                        // coulombic
        force_coul = ni == 0 ?
          qri*q[j]/r : qri*q[j]/r*special_coul[ni];

      typej = type[j];
      if (rsq < cut_bucksqi[typej]) {                          // buckingham
        double r6inv = r2inv*r2inv*r2inv,
                        expr = exp(-r*rhoinvi[typej]);
        force_buck = ni == 0 ?
          (r*expr*buck1i[typej]-r6inv*buck2i[typej]) :
          (r*expr*buck1i[typej]-r6inv*buck2i[typej])*special_lj[ni];
      } else force_buck = 0.0;

      fpair = (force_coul + force_buck) * r2inv;

      if (rsq < cut_in_on_sq) {                                // switching
        double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
        fpair  *= rsw*rsw*(3.0 - 2.0*rsw);
      }
      if (rsq > cut_out_on_sq) {
        double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
        fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      fi[0] += fdx;                       // force update
      fi[1] += fdy;
      fi[2] += fdz;
      if (newton_pair || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }
    }
    f[i][0] += fi[0];
    f[i][1] += fi[1];
    f[i][2] += fi[2];
  }
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLong::compute_outer(int eflag, int vflag)
{
  double evdwl,ecoul,fpair,fvirial;
  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, ii, jj, typei, typej, ni, respa_flag;
  int order1 = ewald_order & EWALD_COUL, order6 = ewald_order & EWALD_DISP;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  double qi = 0.0, qri = 0.0, *cutsqi, *cut_bucksqi,
         *buck1i, *buck2i, *buckai, *buckci, *rhoinvi, *offseti;
  double r, rsq, r2inv, force_coul, force_buck;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double respa_buck = 0.0, respa_coul = 0.0, frespa = 0.0;
  double xi[3], d[3], fi[3];

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  for (ii = 0; ii < inum; ii++) {                          // loop over my atoms
    i = ilist[ii];
    if (order1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    typei = type[i];
    offseti = offset[typei];
    buck1i = buck1[typei];
    buck2i = buck2[typei];
    buckai = buck_a[typei];
    buckci = buck_c[typei];
    rhoinvi = rhoinv[typei];
    cutsqi = cutsq[typei];
    cut_bucksqi = cut_bucksq[typei];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    fi[0] = fi[1] = fi[2] = 0.0;
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                          // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      frespa = 1.0;      //check whether and how to compute respa corrections
      respa_coul = 0.0;
      respa_buck = 0.0;
      respa_flag = rsq < cut_in_on_sq ? 1 : 0;
      if (respa_flag && (rsq > cut_in_off_sq)) {
        double rsw = (r-cut_in_off)/cut_in_diff;
        frespa = 1-rsw*rsw*(3.0-2.0*rsw);
      }

      if (order1 && (rsq < cut_coulsq)) {                // coulombic
        if (!ncoultablebits || rsq <= tabinnersq) {        // series real space
          double pre = qri*q[j];                            // qqrd2e * qi * qj
          if (respa_flag)                                // correct for respa
            respa_coul = ni == 0 ? frespa*pre/r : frespa*pre/r*special_coul[ni];
          double grij = g_ewald*r, t = 1.0/(1.0+EWALD_P*grij);
          double erfc_poly = ((((t*A5+A4)*t+A3)*t+A2)*t+A1);
          if (ni == 0) {
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre-respa_coul;
            if (eflag) ecoul = t;
          } else {                                        // correct for special
            double adjust = pre*(1.0-special_coul[ni])/r;
            pre *= g_ewald*exp(-grij*grij);
            t *= erfc_poly*pre/grij;
            force_coul = t+EWALD_F*pre-adjust-respa_coul;
            if (eflag) ecoul = t-adjust;
          }
        } else {                                         // table real space
          if (respa_flag) {
            double pre = qri*q[j];
            respa_coul = ni == 0 ? frespa*pre/r : frespa*pre/r*special_coul[ni];
          }
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          const int k = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
          double fraction = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]);
            if (eflag) ecoul = qiqj*(etable[k]+fraction*detable[k]);
          } else {                                        // correct for special
            rsq_lookup.f = (1.0-special_coul[ni])*(ctable[k]+fraction*dctable[k]);
            force_coul = qiqj*(ftable[k]+fraction*dftable[k]-(double)rsq_lookup.f);
            if (eflag) {
              rsq_lookup.f = (1.0-special_coul[ni])*(ptable[k]+fraction*dptable[k]);
              ecoul = qiqj*(etable[k]+fraction*detable[k]-(double)rsq_lookup.f);
            }
          }
        }
      } else force_coul = respa_coul = ecoul = 0.0;

      if (rsq < cut_bucksqi[typej]) {                        // buckingham
        double r6inv = r2inv*r2inv*r2inv,
                        expr = exp(-r*rhoinvi[typej]);
        if (respa_flag) respa_buck = ni == 0 ?                 // correct for respa
            frespa*(r*expr*buck1i[typej]-r6inv*buck2i[typej]) :
            frespa*(r*expr*buck1i[typej]-r6inv*buck2i[typej])*special_lj[ni];
        if (order6) {                                        // long-range form
          if (!ndisptablebits || rsq <= tabinnerdispsq) {
            double gr2 = g2*rsq, a2 = 1.0/gr2;
            double expterm = a2*exp(-gr2)*buckci[typej];
            double g6term = g6*((a2+1.0)*a2+0.5)*expterm;
            double g8term = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*expterm*rsq;
            if (ni == 0) {
              force_buck = r*expr*buck1i[typej]-g8term-respa_buck;
              if (eflag) evdwl = expr*buckai[typej]-g6term;
            } else {                                        // correct for special
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_buck = factor*r*expr*buck1i[typej]-g8term+t*buck2i[typej]-respa_buck;
              if (eflag) evdwl = factor*expr*buckai[typej]-g6term+t*buckci[typej];
            }
          } else {          // table real space
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            const int disp_k = (rsq_lookup.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double ftable_disp = fdisptable[disp_k]+f_disp*dfdisptable[disp_k];
            double etable_disp = edisptable[disp_k]+f_disp*dedisptable[disp_k];
            if (ni == 0) {
              force_buck = r*expr*buck1i[typej]-ftable_disp*buckci[typej]-respa_buck;
              if (eflag) evdwl = expr*buckai[typej]-etable_disp*buckci[typej];
            } else {                             //special case
              double factor = special_lj[ni], t = r6inv*(1.0-factor);
              force_buck = factor*r*expr*buck1i[typej]-ftable_disp*buckci[typej]+t*buck2i[typej]-respa_buck;
              if (eflag) evdwl = factor*expr*buckai[typej]-etable_disp*buckci[typej]+t*buckci[typej];
            }
          }
        } else {                                                // cut form
          if (ni == 0) {
            force_buck = r*expr*buck1i[typej]-r6inv*buck2i[typej]-respa_buck;
            if (eflag)
              evdwl = expr*buckai[typej]-r6inv*buckci[typej]-offseti[typej];
          } else {                                        // correct for special
            double factor = special_lj[ni];
            force_buck = factor*(r*expr*buck1i[typej]-r6inv*buck2i[typej])-respa_buck;
            if (eflag)
              evdwl = factor*(expr*buckai[typej]-r6inv*buckci[typej]-offseti[typej]);
          }
        }
      } else force_buck = respa_buck = evdwl = 0.0;

      fpair = (force_coul+force_buck)*r2inv;

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      fi[0] += fdx;                        // force update
      fi[1] += fdy;
      fi[2] += fdz;
      if (newton_pair || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }

      if (evflag) {
        fvirial = (force_coul + force_buck + respa_coul + respa_buck)*r2inv;
        ev_tally(i,j,nlocal,newton_pair,
                 evdwl,ecoul,fvirial,d[0],d[1],d[2]);
      }
    }
    f[i][0] += fi[0];
    f[i][1] += fi[1];
    f[i][2] += fi[2];
  }
}

/* ---------------------------------------------------------------------- */

double PairBuckLongCoulLong::single(int i, int j, int itype, int jtype,
                            double rsq, double factor_coul, double factor_buck,
                            double &fforce)
{
  double r, r2inv, r6inv, force_coul, force_buck;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2, *q = atom->q;

  r = sqrt(rsq);
  r2inv = 1.0/rsq;
  double eng = 0.0;

  if ((ewald_order & EWALD_COUL) && (rsq < cut_coulsq)) {     // coulombic
    if (!ncoultablebits || rsq <= tabinnersq) {                // series real space
      double grij = g_ewald*r;
      double pre = force->qqrd2e*q[i]*q[j];
      double t = 1.0/(1.0+EWALD_P*grij);
      double erfc_poly = ((((t*A5+A4)*t+A3)*t+A2)*t+A1);
      double adjust = pre*(1.0-factor_coul)/r;
      pre *= g_ewald*exp(-grij*grij);
      t *= erfc_poly*pre/grij;
      force_coul = t+EWALD_F*pre-adjust;
      eng += t-adjust;
    } else {                                                // table real space
      union_int_float_t rsq_lookup;
      rsq_lookup.f = rsq;
      const int k = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
      double fraction = (rsq-rtable[k])*drtable[k], qiqj = q[i]*q[j];
      rsq_lookup.f = (1.0-factor_coul)*(ctable[k]+fraction*dctable[k]);
      force_coul = qiqj*(ftable[k]+fraction*dftable[k]-(double)rsq_lookup.f);
      eng += qiqj*(etable[k]+fraction*detable[k]-(double)rsq_lookup.f);
    }
  } else force_coul = 0.0;

  if (rsq < cut_bucksq[itype][jtype]) {                        // buckingham
    double expr = factor_buck*exp(-sqrt(rsq)*rhoinv[itype][jtype]);
    r6inv = r2inv*r2inv*r2inv;
    if (ewald_order & EWALD_DISP) {                      // long-range
      double gr2 = g2*rsq, a2 = 1.0/gr2;
      double expterm = a2*exp(-gr2)*buck_c[itype][jtype];
      double t = r6inv*(1.0-factor_buck);
      double g6term = g6*((a2+1.0)*a2+0.5)*expterm;
      double g8term = g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*expterm*rsq;
      force_buck = buck1[itype][jtype]*r*expr-g8term+t*buck2[itype][jtype];
      eng += buck_a[itype][jtype]*expr-g6term+t*buck_c[itype][jtype];
    } else {                                                // cut
      force_buck =
        factor_buck*(buck1[itype][jtype]*r*expr-buck2[itype][jtype]*r6inv);
      eng += buck_a[itype][jtype]*expr-
        factor_buck*(buck_c[itype][jtype]*r6inv-offset[itype][jtype]);
    }
  } else force_buck = 0.0;

  fforce = (force_coul+force_buck)*r2inv;
  return eng;
}

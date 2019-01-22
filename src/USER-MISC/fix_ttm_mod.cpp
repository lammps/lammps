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
   Contributing authors: (in addition to authors of original fix ttm)
   Sergey Starikov (Joint Institute for High Temperatures of RAS)
   Vasily Pisarev (Joint Institute for High Temperatures of RAS)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "fix_ttm_mod.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define MAXLINE 1024

static const char cite_fix_ttm_mod[] =
  "fix ttm/mod command:\n\n"
  "@article{Pisarev2014,\n"
  "author = {Pisarev, V. V. and Starikov, S. V.},\n"
  "title = {{Atomistic simulation of ion track formation in UO2.}},\n"
  "journal = {J.~Phys.:~Condens.~Matter},\n"
  "volume = {26},\n"
  "number = {47},\n"
  "pages = {475401},\n"
  "year = {2014}\n"
  "}\n\n"
  "@article{Norman2013,\n"
  "author = {Norman, G. E. and Starikov, S. V. and Stegailov, V. V. and Saitov, I. M. and Zhilyaev, P. A.},\n"
  "title = {{Atomistic Modeling of Warm Dense Matter in the Two-Temperature State}},\n"
  "journal = {Contrib.~Plasm.~Phys.},\n"
  "number = {2},\n"
  "volume = {53},\n"
  "pages = {129--139},\n"
  "year = {2013}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixTTMMod::FixTTMMod(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_ttm_mod);

  if (narg < 9) error->all(FLERR,"Illegal fix ttm/mod command");
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  seed = force->inumeric(FLERR,arg[3]);
  if (seed <= 0)
    error->all(FLERR,"Invalid random number seed in fix ttm/mod command");

  FILE *fpr_2 = force->open_potential(arg[4]);
  if (fpr_2 == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",arg[4]);
    error->all(FLERR,str);
  }

  nxnodes = force->inumeric(FLERR,arg[5]);
  nynodes = force->inumeric(FLERR,arg[6]);
  nznodes = force->inumeric(FLERR,arg[7]);
  if (nxnodes <= 0 || nynodes <= 0 || nznodes <= 0)
    error->all(FLERR,"Fix ttm/mod number of nodes must be > 0");

  FILE *fpr = force->open_potential(arg[8]);
  if (fpr == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",arg[8]);
    error->all(FLERR,str);
  }

  nfileevery = force->inumeric(FLERR,arg[9]);
  if (nfileevery > 0) {
    if (narg != 11) error->all(FLERR,"Illegal fix ttm/mod command");
    MPI_Comm_rank(world,&me);
    if (me == 0) {
      fp = fopen(arg[10],"w");
      if (fp == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open fix ttm/mod file %s",arg[10]);
        error->one(FLERR,str);
      }
    }
  }
  char linee[MAXLINE];
  double tresh_d;
  int tresh_i;
  // C0 (metal)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  esheat_0 = tresh_d;
  // C1 (metal*10^3)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  esheat_1 = tresh_d;
  // C2 (metal*10^6)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  esheat_2 = tresh_d;
  // C3 (metal*10^9)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  esheat_3 = tresh_d;
  // C4 (metal*10^12)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  esheat_4 = tresh_d;
  // C_limit
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  C_limit = tresh_d;
  //Temperature damping factor
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  T_damp = tresh_d;
  // rho_e
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  electronic_density = tresh_d;
  //thermal_diffusion
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  el_th_diff = tresh_d;
  // gamma_p
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  gamma_p = tresh_d;
  // gamma_s
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  gamma_s = tresh_d;
  // v0
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  v_0 = tresh_d;
  // average intensity of pulse (source of energy) (metal units)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  intensity = tresh_d;
  // coordinate of 1st surface in x-direction (in box units) - constant
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%d",&tresh_i);
  surface_l = tresh_i;
  // coordinate of 2nd surface in x-direction (in box units) - constant
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%d",&tresh_i);
  surface_r = tresh_i;
  // skin_layer = intensity is reduced (I=I0*exp[-x/skin_layer])
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%d",&tresh_i);
  skin_layer = tresh_i;
  // width of pulse (picoseconds)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  width = tresh_d;
  // factor of electronic pressure (PF) Pe = PF*Ce*Te
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  pres_factor = tresh_d;
  // effective free path of electrons (angstrom)
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  free_path = tresh_d;
  // ionic density (ions*angstrom^{-3})
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  ionic_density = tresh_d;
  // if movsur = 0: surface is freezed
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%d",&tresh_i);
  movsur = tresh_i;
  // electron_temperature_min
  fgets(linee,MAXLINE,fpr_2);
  fgets(linee,MAXLINE,fpr_2);
  sscanf(linee,"%lg",&tresh_d);
  electron_temperature_min = tresh_d;
  fclose(fpr_2);
  //t_surface is determined by electronic temperature (not constant)
  t_surface_l = surface_l;
  mult_factor = intensity;
  duration = 0.0;
  v_0_sq = v_0*v_0;
  surface_double = double(t_surface_l)*(domain->xprd/nxnodes);
  if ((C_limit+esheat_0) < 0.0)
    error->all(FLERR,"Fix ttm/mod electronic_specific_heat must be >= 0.0");
  if (electronic_density <= 0.0)
    error->all(FLERR,"Fix ttm/mod electronic_density must be > 0.0");
  if (gamma_p < 0.0) error->all(FLERR,"Fix ttm/mod gamma_p must be >= 0.0");
  if (gamma_s < 0.0) error->all(FLERR,"Fix ttm/mod gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all(FLERR,"Fix ttm/mod v_0 must be >= 0.0");
  if (ionic_density <= 0.0) error->all(FLERR,"Fix ttm/mod ionic_density must be > 0.0");
  if (surface_l < 0) error->all(FLERR,"Surface coordinates must be >= 0");
  if (surface_l >= surface_r) error->all(FLERR, "Left surface coordinate must be less than right surface coordinate");
  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
  // allocate per-type arrays for force prefactors
  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];
  // allocate 3d grid variables
  total_nnodes = nxnodes*nynodes*nznodes;
  memory->create(nsum,nxnodes,nynodes,nznodes,"ttm/mod:nsum");
  memory->create(nsum_all,nxnodes,nynodes,nznodes,"ttm/mod:nsum_all");
  memory->create(T_initial_set,nxnodes,nynodes,nznodes,"ttm/mod:T_initial_set");
  memory->create(sum_vsq,nxnodes,nynodes,nznodes,"ttm/mod:sum_vsq");
  memory->create(sum_mass_vsq,nxnodes,nynodes,nznodes,"ttm/mod:sum_mass_vsq");
  memory->create(sum_vsq_all,nxnodes,nynodes,nznodes,"ttm/mod:sum_vsq_all");
  memory->create(sum_mass_vsq_all,nxnodes,nynodes,nznodes,
                 "ttm/mod:sum_mass_vsq_all");
  memory->create(T_electron_old,nxnodes,nynodes,nznodes,"ttm/mod:T_electron_old");
  memory->create(T_electron_first,nxnodes,nynodes,nznodes,"ttm/mod:T_electron_first");
  memory->create(T_electron,nxnodes,nynodes,nznodes,"ttm/mod:T_electron");
  memory->create(net_energy_transfer,nxnodes,nynodes,nznodes,
                 "ttm/mod:net_energy_transfer");
  memory->create(net_energy_transfer_all,nxnodes,nynodes,nznodes,
                 "ttm/mod:net_energy_transfer_all");
  flangevin = NULL;
  grow_arrays(atom->nmax);
  // zero out the flangevin array
  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0;
    flangevin[i][1] = 0;
    flangevin[i][2] = 0;
  }
  atom->add_callback(0);
  atom->add_callback(1);
  // set initial electron temperatures from user input file
  if (me == 0) read_initial_electron_temperatures(fpr);
  MPI_Bcast(&T_electron[0][0][0],total_nnodes,MPI_DOUBLE,0,world);
  fclose(fpr);
}

/* ---------------------------------------------------------------------- */

FixTTMMod::~FixTTMMod()
{
  if (nfileevery && me == 0) fclose(fp);
  delete random;
  delete [] gfactor1;
  delete [] gfactor2;
  memory->destroy(nsum);
  memory->destroy(nsum_all);
  memory->destroy(T_initial_set);
  memory->destroy(sum_vsq);
  memory->destroy(sum_mass_vsq);
  memory->destroy(sum_vsq_all);
  memory->destroy(sum_mass_vsq_all);
  memory->destroy(T_electron_first);
  memory->destroy(T_electron_old);
  memory->destroy(T_electron);
  memory->destroy(flangevin);
  memory->destroy(net_energy_transfer);
  memory->destroy(net_energy_transfer_all);
}

/* ---------------------------------------------------------------------- */

int FixTTMMod::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::init()
{
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix ttm/mod with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use non-periodic boundares with fix ttm/mod");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix ttm/mod with triclinic box");
  // set force prefactors
  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - gamma_p / force->ftm2v;
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        net_energy_transfer_all[ixnode][iynode][iznode] = 0;
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    post_force_setup(vflag);
  } else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa_setup(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nynodes;
  double gamma1,gamma2;
  // apply damping and thermostat to all atoms in fix group
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;
      if (T_electron[ixnode][iynode][iznode] < 0)
        error->all(FLERR,"Electronic temperature dropped below zero");
      double tsqrt = sqrt(T_electron[ixnode][iynode][iznode]);
      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;
      if (ixnode >= surface_l){
        if (ixnode < surface_r){
          flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
          flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
          flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
          double x_surf = dx*double(surface_l)+dx;
          double x_at = x[i][0] - domain->boxlo[0];
          int right_xnode = ixnode + 1;
          int right_ynode = iynode + 1;
          int right_znode = iznode + 1;
          if (right_xnode == nxnodes) right_xnode = 0;
          if (right_ynode == nynodes) right_ynode = 0;
          if (right_znode == nznodes) right_znode = 0;
          int left_xnode = ixnode - 1;
          int left_ynode = iynode - 1;
          int left_znode = iznode - 1;
          if (left_xnode == -1) left_xnode = nxnodes - 1;
          if (left_ynode == -1) left_ynode = nynodes - 1;
          if (left_znode == -1) left_znode = nznodes - 1;
          double T_i = T_electron[ixnode][iynode][iznode];
          double T_ir = T_electron[right_xnode][iynode][iznode];
          double T_iu = T_electron[ixnode][right_ynode][iznode];
          double T_if = T_electron[ixnode][iynode][right_znode];
          double C_i = el_properties(T_electron[ixnode][iynode][iznode]).el_heat_capacity;
          double C_ir = el_properties(T_electron[right_xnode][iynode][iznode]).el_heat_capacity;
          double C_iu = el_properties(T_electron[ixnode][right_ynode][iznode]).el_heat_capacity;
          double C_if = el_properties(T_electron[ixnode][iynode][right_znode]).el_heat_capacity;
          double diff_x = (x_at - x_surf)*(x_at - x_surf);
          diff_x = pow(diff_x,0.5);
          double len_factor = diff_x/(diff_x+free_path);
          if (movsur == 1){
            if (x_at >= x_surf){
              flangevin[i][0] -= pres_factor/ionic_density*((C_ir*T_ir*free_path/(diff_x+free_path)/(diff_x+free_path)) +
                                                            (len_factor/dx)*(C_ir*T_ir-C_i*T_i));
              flangevin[i][1] -= pres_factor/ionic_density/dy*(C_iu*T_iu-C_i*T_i);
              flangevin[i][2] -= pres_factor/ionic_density/dz*(C_if*T_if-C_i*T_i);
            }
          } else {
            flangevin[i][0] -= pres_factor/ionic_density/dx*(C_ir*T_ir-C_i*T_i);
            flangevin[i][1] -= pres_factor/ionic_density/dy*(C_iu*T_iu-C_i*T_i);
            flangevin[i][2] -= pres_factor/ionic_density/dz*(C_if*T_if-C_i*T_i);
          }
          f[i][0] += flangevin[i][0];
          f[i][1] += flangevin[i][1];
          f[i][2] += flangevin[i][2];
        }
      }
      if (movsur == 1){
        if (ixnode < surface_l){
          t_surface_l = ixnode;
        }
      }
    }
  }
  MPI_Allreduce(&t_surface_l,&surface_l,1,MPI_INT,MPI_MIN,world);
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::post_force_setup(int /*vflag*/)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  // apply langevin forces that have been stored from previous run
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::post_force_respa_setup(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force_setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::reset_dt()
{
  for (int i = 1; i <= atom->ntypes; i++)
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixTTMMod::read_initial_electron_temperatures(FILE *fpr)
{
  char line[MAXLINE];
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_initial_set[ixnode][iynode][iznode] = 0;
  // read initial electron temperature values from file
  int ixnode,iynode,iznode;
  double T_tmp;
  while (1) {
    if (fgets(line,MAXLINE,fpr) == NULL) break;
    sscanf(line,"%d %d %d %lg",&ixnode,&iynode,&iznode,&T_tmp);
    if (T_tmp < 0.0) error->one(FLERR,"Fix ttm/mod electron temperatures must be >= 0.0");
    T_electron[ixnode][iynode][iznode] = T_tmp;
    T_initial_set[ixnode][iynode][iznode] = 1;
  }
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        if (T_initial_set[ixnode][iynode][iznode] == 0)
          error->one(FLERR,"Initial temperatures not all set in fix ttm/mod");
}

/* ---------------------------------------------------------------------- */

el_heat_capacity_thermal_conductivity FixTTMMod::el_properties(double T_e)
{
  el_heat_capacity_thermal_conductivity properties;
  double T_temp = T_e/1000.0, T_reduced = T_damp*T_temp;
  double T2 = T_temp*T_temp;
  double T3 = T2*T_temp;
  double T4 = T3*T_temp;
  double poly = esheat_0 + esheat_1*T_temp + esheat_2*T2 + esheat_3*T3 + esheat_4*T4;
  properties.el_heat_capacity = electronic_density*(poly*exp(-T_reduced*T_reduced) + C_limit); // heat capacity
  properties.el_thermal_conductivity = el_th_diff*properties.el_heat_capacity; // thermal conductivity
  return properties;
}
double FixTTMMod::el_sp_heat_integral(double T_e)
{
  double T_temp = T_e/1000.0, T_reduced = T_damp*T_temp;
  if (T_damp != 0)
    return electronic_density*(MY_PIS*(3*esheat_4/pow(T_damp,5)+2*esheat_2/pow(T_damp,3)+4*esheat_0/T_damp)*erf(T_reduced)+
                               4*esheat_3/pow(T_damp,4)+4*esheat_1/T_damp/T_damp-
                               ((6*esheat_4*T_temp+4*esheat_3)/pow(T_damp,4)+
                                (4*esheat_1+4*esheat_4*pow(T_temp,3)+4*esheat_3*T_temp*T_temp+4*esheat_2*T_temp)/T_damp/T_damp)*exp(-T_reduced*T_reduced))*125.0+electronic_density*C_limit*T_e;
  else
    return electronic_density*((esheat_0 + C_limit)*T_e + esheat_1*T_temp*T_e/2.0 + esheat_2*T_temp*T_temp*T_e/3.0 + esheat_3*pow(T_temp,3)*T_e/4.0 + esheat_4*pow(T_temp,4)*T_e/5.0);
}
void FixTTMMod::end_of_step()
{
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (movsur == 1){
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++){
          double TTT = T_electron[ixnode][iynode][iznode];
          if (TTT > 0){
            if (ixnode < t_surface_l)
              t_surface_l = ixnode;
          }
        }
  }
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        net_energy_transfer[ixnode][iynode][iznode] = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;
      if (ixnode >= t_surface_l){
        if (ixnode < surface_r)
          net_energy_transfer[ixnode][iynode][iznode] +=
            (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
             flangevin[i][2]*v[i][2]);
      }
    }
  MPI_Allreduce(&net_energy_transfer[0][0][0],
                &net_energy_transfer_all[0][0][0],
                total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;
  double el_specific_heat = 0.0;
  double el_thermal_conductivity = el_properties(electron_temperature_min).el_thermal_conductivity;
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        {
          if (el_properties(T_electron[ixnode][iynode][iznode]).el_thermal_conductivity > el_thermal_conductivity)
            el_thermal_conductivity = el_properties(T_electron[ixnode][iynode][iznode]).el_thermal_conductivity;
          if (el_specific_heat > 0.0)
            {
              if ((T_electron[ixnode][iynode][iznode] > 0.0) && (el_properties(T_electron[ixnode][iynode][iznode]).el_heat_capacity < el_specific_heat))
                el_specific_heat = el_properties(T_electron[ixnode][iynode][iznode]).el_heat_capacity;
            }
          else if (T_electron[ixnode][iynode][iznode] > 0.0) el_specific_heat = el_properties(T_electron[ixnode][iynode][iznode]).el_heat_capacity;
        }
  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve
  int num_inner_timesteps = 1;
  double inner_dt = update->dt;
  double stability_criterion = 0.0;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_electron_first[ixnode][iynode][iznode] =
          T_electron[ixnode][iynode][iznode];
  do {
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          T_electron[ixnode][iynode][iznode] =
            T_electron_first[ixnode][iynode][iznode];

    stability_criterion = 1.0 -
      2.0*inner_dt/el_specific_heat *
      (el_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    if (stability_criterion < 0.0) {
      inner_dt = 0.25*el_specific_heat /
        (el_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    }
    num_inner_timesteps = static_cast<unsigned int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm/mod",0);
    for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps;
         ith_inner_timestep++) {
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            T_electron_old[ixnode][iynode][iznode] =
              T_electron[ixnode][iynode][iznode];
      // compute new electron T profile
      duration = duration + inner_dt;
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            int right_xnode = ixnode + 1;
            int right_ynode = iynode + 1;
            int right_znode = iznode + 1;
            if (right_xnode == nxnodes) right_xnode = 0;
            if (right_ynode == nynodes) right_ynode = 0;
            if (right_znode == nznodes) right_znode = 0;
            int left_xnode = ixnode - 1;
            int left_ynode = iynode - 1;
            int left_znode = iznode - 1;
            if (left_xnode == -1) left_xnode = nxnodes - 1;
            if (left_ynode == -1) left_ynode = nynodes - 1;
            if (left_znode == -1) left_znode = nznodes - 1;
            double skin_layer_d = double(skin_layer);
            double ixnode_d = double(ixnode);
            double surface_d = double(t_surface_l);
            mult_factor = 0.0;
            if (duration < width){
              if (ixnode >= t_surface_l) mult_factor = (intensity/(dx*skin_layer_d))*exp((-1.0)*(ixnode_d - surface_d)/skin_layer_d);
            }
            if (ixnode < t_surface_l) net_energy_transfer_all[ixnode][iynode][iznode] = 0.0;
            double cr_vac = 1;
            if (T_electron_old[ixnode][iynode][iznode] == 0) cr_vac = 0;
            double cr_v_l_x = 1;
            if (T_electron_old[left_xnode][iynode][iznode] == 0) cr_v_l_x = 0;
            double cr_v_r_x = 1;
            if (T_electron_old[right_xnode][iynode][iznode] == 0) cr_v_r_x = 0;
            double cr_v_l_y = 1;
            if (T_electron_old[ixnode][left_ynode][iznode] == 0) cr_v_l_y = 0;
            double cr_v_r_y = 1;
            if (T_electron_old[ixnode][right_ynode][iznode] == 0) cr_v_r_y = 0;
            double cr_v_l_z = 1;
            if (T_electron_old[ixnode][iynode][left_znode] == 0) cr_v_l_z = 0;
            double cr_v_r_z = 1;
            if (T_electron_old[ixnode][iynode][right_znode] == 0) cr_v_r_z = 0;
            if (cr_vac != 0) {
              T_electron[ixnode][iynode][iznode] =
                T_electron_old[ixnode][iynode][iznode] +
                inner_dt/el_properties(T_electron_old[ixnode][iynode][iznode]).el_heat_capacity *
                ((cr_v_r_x*el_properties(T_electron_old[ixnode][iynode][iznode]/2.0+T_electron_old[right_xnode][iynode][iznode]/2.0).el_thermal_conductivity*
                  (T_electron_old[right_xnode][iynode][iznode]-T_electron_old[ixnode][iynode][iznode])/dx -
                  cr_v_l_x*el_properties(T_electron_old[ixnode][iynode][iznode]/2.0+T_electron_old[left_xnode][iynode][iznode]/2.0).el_thermal_conductivity*
                  (T_electron_old[ixnode][iynode][iznode]-T_electron_old[left_xnode][iynode][iznode])/dx)/dx +
                 (cr_v_r_y*el_properties(T_electron_old[ixnode][iynode][iznode]/2.0+T_electron_old[ixnode][right_ynode][iznode]/2.0).el_thermal_conductivity*
                  (T_electron_old[ixnode][right_ynode][iznode]-T_electron_old[ixnode][iynode][iznode])/dy -
                  cr_v_l_y*el_properties(T_electron_old[ixnode][iynode][iznode]/2.0+T_electron_old[ixnode][left_ynode][iznode]/2.0).el_thermal_conductivity*
                  (T_electron_old[ixnode][iynode][iznode]-T_electron_old[ixnode][left_ynode][iznode])/dy)/dy +
                 (cr_v_r_z*el_properties(T_electron_old[ixnode][iynode][iznode]/2.0+T_electron_old[ixnode][iynode][right_znode]/2.0).el_thermal_conductivity*
                  (T_electron_old[ixnode][iynode][right_znode]-T_electron_old[ixnode][iynode][iznode])/dz -
                  cr_v_l_z*el_properties(T_electron_old[ixnode][iynode][iznode]/2.0+T_electron_old[ixnode][iynode][left_znode]/2.0).el_thermal_conductivity*
                  (T_electron_old[ixnode][iynode][iznode]-T_electron_old[ixnode][iynode][left_znode])/dz)/dz);
              T_electron[ixnode][iynode][iznode]+=inner_dt/el_properties(T_electron[ixnode][iynode][iznode]).el_heat_capacity*
                (mult_factor -
                 net_energy_transfer_all[ixnode][iynode][iznode]/del_vol);
            }
            else T_electron[ixnode][iynode][iznode] =
                   T_electron_old[ixnode][iynode][iznode];
            if ((T_electron[ixnode][iynode][iznode] > 0.0) && (T_electron[ixnode][iynode][iznode] < electron_temperature_min))
              T_electron[ixnode][iynode][iznode] = T_electron[ixnode][iynode][iznode] + 0.5*(electron_temperature_min - T_electron[ixnode][iynode][iznode]);

            if (el_properties(T_electron[ixnode][iynode][iznode]).el_thermal_conductivity > el_thermal_conductivity)
              el_thermal_conductivity = el_properties(T_electron[ixnode][iynode][iznode]).el_thermal_conductivity;
            if ((T_electron[ixnode][iynode][iznode] > 0.0) && (el_properties(T_electron[ixnode][iynode][iznode]).el_heat_capacity < el_specific_heat))
              el_specific_heat = el_properties(T_electron[ixnode][iynode][iznode]).el_heat_capacity;
          }
    }
    stability_criterion = 1.0 -
      2.0*inner_dt/el_specific_heat *
      (el_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));

  } while (stability_criterion < 0.0);
  // output nodal temperatures for current timestep
  if ((nfileevery) && !(update->ntimestep % nfileevery)) {
    // compute atomic Ta for each grid point
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          nsum[ixnode][iynode][iznode] = 0;
          nsum_all[ixnode][iynode][iznode] = 0;
          sum_vsq[ixnode][iynode][iznode] = 0.0;
          sum_mass_vsq[ixnode][iynode][iznode] = 0.0;
          sum_vsq_all[ixnode][iynode][iznode] = 0.0;
          sum_mass_vsq_all[ixnode][iynode][iznode] = 0.0;
        }
    double massone;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
        double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
        double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
        int ixnode = static_cast<int>(xscale*nxnodes);
        int iynode = static_cast<int>(yscale*nynodes);
        int iznode = static_cast<int>(zscale*nznodes);
        while (ixnode > nxnodes-1) ixnode -= nxnodes;
        while (iynode > nynodes-1) iynode -= nynodes;
        while (iznode > nznodes-1) iznode -= nznodes;
        while (ixnode < 0) ixnode += nxnodes;
        while (iynode < 0) iynode += nynodes;
        while (iznode < 0) iznode += nznodes;
        double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        nsum[ixnode][iynode][iznode] += 1;
        sum_vsq[ixnode][iynode][iznode] += vsq;
        sum_mass_vsq[ixnode][iynode][iznode] += massone*vsq;
      }
    MPI_Allreduce(&nsum[0][0][0],&nsum_all[0][0][0],total_nnodes,
                  MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&sum_vsq[0][0][0],&sum_vsq_all[0][0][0],total_nnodes,
                  MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&sum_mass_vsq[0][0][0],&sum_mass_vsq_all[0][0][0],
                  total_nnodes,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&t_surface_l,&surface_l,
                  1,MPI_INT,MPI_MIN,world);
    if (me == 0) {
      fprintf(fp,BIGINT_FORMAT,update->ntimestep);
      double T_a;
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            T_a = 0;
            if (nsum_all[ixnode][iynode][iznode] > 0){
              T_a = sum_mass_vsq_all[ixnode][iynode][iznode]/
                (3.0*force->boltz*nsum_all[ixnode][iynode][iznode]/force->mvv2e);
              if (movsur == 1){
                if (T_electron[ixnode][iynode][iznode]==0.0) T_electron[ixnode][iynode][iznode] = electron_temperature_min;
              }
            }
            fprintf(fp," %f",T_a);
          }
      fprintf(fp,"\t");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            fprintf(fp,"%f ",T_electron[ixnode][iynode][iznode]);
      fprintf(fp,"\n");
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of 3d grid
------------------------------------------------------------------------- */

double FixTTMMod::memory_usage()
{
  double bytes = 0.0;
  bytes += 5*total_nnodes * sizeof(int);
  bytes += 14*total_nnodes * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixTTMMod::grow_arrays(int ngrow)
{
  memory->grow(flangevin,ngrow,3,"ttm/mod:flangevin");
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTMMod::compute_vector(int n)
{
  double e_energy = 0.0;
  double transfer_energy = 0.0;
  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        e_energy += el_sp_heat_integral(T_electron[ixnode][iynode][iznode])*del_vol;
        transfer_energy +=
          net_energy_transfer_all[ixnode][iynode][iznode]*update->dt;
      }
  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTMMod::write_restart(FILE *fp)
{
  double *rlist;
  memory->create(rlist,nxnodes*nynodes*nznodes+1,"ttm/mod:rlist");
  int n = 0;
  rlist[n++] = seed;
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        rlist[n++] = T_electron[ixnode][iynode][iznode];
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),n,fp);
  }
  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTMMod::restart(char *buf)
{
  int n = 0;
  double *rlist = (double *) buf;
  // the seed must be changed from the initial seed
  seed = static_cast<int> (0.5*rlist[n++]);
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_electron[ixnode][iynode][iznode] = rlist[n++];
  delete random;
  random = new RanMars(lmp,seed+comm->me);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTTMMod::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = flangevin[i][0];
  buf[2] = flangevin[i][1];
  buf[3] = flangevin[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTTMMod::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;
  // skip to Nth set of extra values
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;
  flangevin[nlocal][0] = extra[nlocal][m++];
  flangevin[nlocal][1] = extra[nlocal][m++];
  flangevin[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTTMMod::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTTMMod::size_restart(int /*nlocal*/)
{
  return 4;
}

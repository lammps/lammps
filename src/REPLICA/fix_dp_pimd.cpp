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
   Package      FixDPPimd
   Purpose      A Path Integral Molecular Dynamics Package developed by DeepModeling community
   Authors      Yifan Li (mail_liyifan@163.com, yifanl@princeton.edu)

   Updated      Jul-02-2021
------------------------------------------------------------------------- */

#include "fix_dp_pimd.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "universe.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "utils.h"
#include "timer.h"
#include "atom.h"
#include "compute.h"
#include "modify.h"
#include "domain.h"
#include "update.h"
#include "math_const.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{PIMD,NMPIMD,CMD};
enum{physical, normal};
enum{baoab};
enum{SVR, PILE_L, PILE_G};
enum{nve, nvt, nph, npt};
enum{MSTI, SCTI};

#define INVOKED_SCALAR 1

/* ---------------------------------------------------------------------- */

FixDPPimd::FixDPPimd(LAMMPS *lmp, int narg, char **arg) : 
  Fix(lmp, narg, arg),
  random(nullptr), c_pe(nullptr), c_press(nullptr)
{
  method       = NMPIMD;
  fmmode       = physical;
  integrator   = baoab;
  thermostat   = PILE_L;
  ensemble     = nvt;
  fmass        = 1.0;
  temp         = 298.15;
  baoab_temp   = 298.15;
  sp           = 1.0;
  harmonicflag = 0;
  omega        = 0.0;
  tiflag       = 0;
  timethod     = MSTI;
  lambda       = 0.0;
  pextflag     = 0;

  for(int i=3; i<narg-1; i+=2)
  {
    if(strcmp(arg[i],"method")==0)
    {
      if(strcmp(arg[i+1],"pimd")==0) method=PIMD;
      else if(strcmp(arg[i+1],"nmpimd")==0) method=NMPIMD;
      else if(strcmp(arg[i+1],"cmd")==0) method=CMD;
      else error->universe_all(FLERR,"Unknown method parameter for fix pimd");
    }

    else if(strcmp(arg[i], "integrator")==0)
    {
      if(strcmp(arg[i+1], "baoab")==0) integrator=baoab;
      else error->universe_all(FLERR, "Unknown integrator parameter for fix pimd. Only baoab integrator is supported!");
    }

    else if(strcmp(arg[i], "ensemble")==0)
    {
      if(strcmp(arg[i+1], "nve")==0) ensemble=nve;
      else if(strcmp(arg[i+1], "nvt")==0) ensemble=nvt;
      else if(strcmp(arg[i+1], "nph")==0) {ensemble=nph; pextflag=1;}
      else if(strcmp(arg[i+1], "npt")==0) {ensemble=npt; pextflag=1;}
      else error->universe_all(FLERR, "Unknown ensemble parameter for fix pimd. Only nve and nvt ensembles are supported!");
    }

    else if(strcmp(arg[i],"fmass")==0)
    {
      fmass = atof(arg[i+1]);
      if(fmass<0.0 || fmass>1.0) error->universe_all(FLERR,"Invalid fmass value for fix pimd");
    }

    else if(strcmp(arg[i], "fmmode")==0)
    {
      if(strcmp(arg[i+1], "physical")==0) fmmode=physical;
      else if(strcmp(arg[i+1], "normal")==0) fmmode=normal;
      else error->universe_all(FLERR, "Unknown fictitious mass mode for fix pimd. Only physical mass and normal mode mass are supported!");
    }

    else if(strcmp(arg[i],"scale")==0)
    {
      pilescale = atof(arg[i+1]);
      if(pilescale<0.0) error->universe_all(FLERR,"Invalid pile scale value for fix pimd");
    }

    else if(strcmp(arg[i],"sp")==0)
    {
      sp = atof(arg[i+1]);
      if(fmass<0.0) error->universe_all(FLERR,"Invalid sp value for fix pimd");
    }

    else if(strcmp(arg[i],"temp")==0)
    {
      temp = atof(arg[i+1]);
      if(temp<0.0) error->universe_all(FLERR,"Invalid temp value for fix pimd");
    } 

    else if(strcmp(arg[i], "press")==0)
    {
      Pext = atof(arg[i+1]);
      if(Pext<0.0) error->universe_all(FLERR,"Invalid press value for fix pimd");
    }

    else if(strcmp(arg[i], "taup")==0)
    {
      tau_p = atof(arg[i+1]);
      if(tau_p<=0.0) error->universe_all(FLERR, "Invalid tau_p value for fix pimd");
    }

    else if(strcmp(arg[i], "thermostat")==0)
    {
      if(strcmp(arg[i+1],"PILE_G")==0) 
      {
        thermostat = PILE_G;
        seed = atoi(arg[i+2]);
        i++;
      }
      else if(strcmp(arg[i+1], "SVR")==0)
      {
        thermostat = SVR;
        seed = atoi(arg[i+2]);
        i++;
      }
      else if(strcmp(arg[i+1],"PILE_L")==0) 
      {
        thermostat = PILE_L;
        seed = atoi(arg[i+2]);
        i++;
      }
      else error->universe_all(FLERR,"Unknown thermostat parameter for fix pimd");
    }

    else if(strcmp(arg[i], "tau")==0)
    {
      tau = atof(arg[i+1]);
    }
  
    else if(strcmp(arg[i], "ti")==0)
    {
      tiflag = 1;
      if(strcmp(arg[i+1], "MSTI")==0)  timethod = MSTI;
      else if(strcmp(arg[i+1], "SCTI")==0)  timethod = SCTI;
      else error->universe_all(FLERR, "Unknown method parameter for thermodynamic integration");
      lambda = atof(arg[i+2]);
      i++;
    }
 
    else if(strcmp(arg[i], "model")==0)
    {
      harmonicflag = 1;
      omega = atof(arg[i+1]);
      if(omega<0) error->universe_all(FLERR,"Invalid model frequency value for fix pimd");
    }
    else error->universe_all(arg[i],i+1,"Unknown keyword for fix pimd");
  }

  // initialize Marsaglia RNG with processor-unique seed

  if(integrator==baoab)
  {
    baoab_temp = temp;
    random = new RanMars(lmp, seed + universe->me);
  }
  
  /* Initiation */

  max_nsend = 0;
  tag_send = nullptr;
  buf_send = nullptr;

  max_nlocal = 0;
  buf_recv = nullptr;
  buf_beads = nullptr;

  coords_send = coords_recv = nullptr;
  forces_send = forces_recv = nullptr;
  nsend = nrecv = 0;
  tags_send = nullptr;
  coords = nullptr;
  forces = nullptr;
  size_plan = 0;
  plan_send = plan_recv = nullptr;

  xc = fc = nullptr;
  xf = 0.0;
  t_vir = t_cv = 0.0;
  total_spring_energy = 0.0;
  t_prim = 0.0;

  for(int i=0; i<9; i++) virial[i] = 0.0;

  tote = totke = totenthalpy = 0.0;
  ke_bead = 0.0;

  dfdl = 0.0;
  x_scaled = nullptr;

  M_x2xp = M_xp2x = nullptr; // M_f2fp = M_fp2f = nullptr;
  lam = nullptr;
  mode_index = nullptr;
  x_unwrap = nullptr;

  mass = nullptr;

  array_atom = nullptr;

  gamma = 0.0;
  c1 = 0.0;
  c2 = 0.0;

  restart_peratom = 1;
  peratom_flag    = 1;
  peratom_freq    = 1;

  global_freq = 1;
  thermo_energy = 1;
  vector_flag = 1;
  size_vector = 12;
  scalar_flag = 1;
  extvector   = 1;
  comm_forward = 3;

  atom->add_callback(0); // Call LAMMPS to allocate memory for per-atom array
  atom->add_callback(1); // Call LAMMPS to re-assign restart-data for per-atom array


  // some initilizations

  baoab_ready = false;

  r1 = 0.0;
  r2 = 0.0;

  id_pe = new char[8];
  strcpy(id_pe, "pimd_pe");
  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;

  id_press = new char[12];
  strcpy(id_press, "pimd_press");
  newarg = new char*[5];
  newarg[0] = id_press;
  newarg[1] = (char*) "all";
  newarg[2] = (char*) "pressure";
  newarg[3] = (char*) "thermo_temp";
  newarg[4] = (char*) "virial";
  modify->add_compute(5, newarg);
  delete [] newarg;
  
  domain->set_global_box();
 
  //FILE *frand;
  //std::string fname = "rand_";
  //fname += std::to_string(universe->iworld);
  //fname += ".txt";
  //frand = fopen(fname.c_str(), "w");
}

/* ---------------------------------------------------------------------- */

FixDPPimd::~FixDPPimd()
{
  delete _omega_k;
  delete baoab_c, baoab_s;
  if(integrator==baoab)
  {
    delete random;
  }

  if(thermostat==PILE_L)
  {
    delete tau_k ,c1_k, c2_k;
  }
  //fclose(frand);
}

/* ---------------------------------------------------------------------- */

int FixDPPimd::setmask()
{
  int mask = 0;
  //mask |= PRE_EXCHANGE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::end_of_step()
{
  compute_totke();
  compute_vir();
  compute_vir_();
  // compute_p_prim();
  compute_p_cv();
  compute_tote();
  if(pextflag) compute_totenthalpy();

  if(update->ntimestep % 10000 == 0)
  {
  if(universe->me==0) printf("This is the end of step %ld.\n", update->ntimestep);
  }

}

/* ---------------------------------------------------------------------- */

void FixDPPimd::init()
{
  if (atom->map_style == 0)
    error->all(FLERR,"Fix pimd requires an atom map, see atom_modify");

  if(universe->me==0 && screen) fprintf(screen,"Fix pimd initializing Path-Integral ...\n");
  // fprintf(stdout, "Fix pimd initilizing Path-Integral ...\n");

  // prepare the constants

  np = universe->nworlds;
  inverse_np = 1.0 / np;

  /* The first solution for the force constant, using SI units

  const double Boltzmann = 1.3806488E-23;    // SI unit: J/K
  const double Plank     = 6.6260755E-34;    // SI unit: m^2 kg / s

  double hbar = Plank / ( 2.0 * MY_PI ) * sp;
  double beta = 1.0 / ( Boltzmann * input.nh_temp);

  // - P / ( beta^2 * hbar^2)   SI unit: s^-2
  double _fbond = -1.0 / (beta*beta*hbar*hbar) * input.nbeads;

  // convert the units: s^-2 -> (kcal/mol) / (g/mol) / (A^2)
  fbond = _fbond * 4.184E+26;

  */

  /* The current solution, using LAMMPS internal real units */

  const double Boltzmann = force->boltz;
  const double Plank     = force->hplanck;

  double hbar   = Plank / ( 2.0 * MY_PI );
  double beta   = 1.0 / (Boltzmann * temp);
  double _fbond = 1.0 * np*np / (beta*beta*hbar*hbar) ;

  omega_np = np / (hbar * beta) * sqrt(force->mvv2e);
  fbond = _fbond * force->mvv2e;

  beta_np = 1.0 / force->boltz / baoab_temp / np;

  if(universe->me==0)
    printf("Fix pimd -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", fbond);

  if(integrator==baoab)   
  {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = 1./3 * dtv2 * dtv;
  }
  else
  {
    error->universe_all(FLERR,"Unknown integrator parameter for fix pimd");
  }

  comm_init();

  mass = new double [atom->ntypes+1];

  if(method==CMD || method==NMPIMD) nmpimd_init();
  else for(int i=1; i<=atom->ntypes; i++) mass[i] = atom->mass[i] / np * fmass;

  if(integrator==baoab)
  {
    if(!baoab_ready)
    {
      baoab_init();
    }
    // fprintf(stdout, "baoab thermostat initialized!\n");
  }
  else error->universe_all(FLERR,"Unknown integrator parameter for fix pimd");

  if(pextflag)
  {
    W = (3*atom->natoms) * tau_p * tau_p / beta_np; // consistent with the definition in i-Pi
    //W = 4 * tau_p * tau_p / beta_np;
    //printf("N=%d, tau_p=%f, beta=%f, W=%f\n", atom->natoms, tau_p, beta_np, W);
    Vcoeff = 1.0;
    vw = 0.0;
  }

  // initialize compute pe 
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];
  
  // initialize compute press
  int ipress = modify->find_compute(id_press);
  c_press = modify->compute[ipress];

  t_prim = t_vir = t_cv = p_prim = p_vir = p_cv = 0.0;

  if(universe->me==0) fprintf(screen, "Fix pimd successfully initialized!\n");
}

void FixDPPimd::setup(int vflag)
{
    if(method==NMPIMD)
    {
      nmpimd_fill(atom->v);
      comm_exec(atom->v);
      nmpimd_transform(buf_beads, atom->v, M_x2xp[universe->iworld]);
    }
  if(universe->me==0 && screen) fprintf(screen,"Setting up Path-Integral ...\n");
  if(universe->me==0) printf("Setting up Path-Integral ...\n");
  post_force(vflag);
  compute_p_cv();
  //compute_p_vir();
  end_of_step();
  //fprintf(stdout, "virial=%.8e.\n", virial[0]+virial[4]+virial[8]);
  //fprintf(stdout, "vir=%.8e.\n", vir);
  c_pe->addstep(update->ntimestep+1); 
  c_press->addstep(update->ntimestep+1);
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  //fprintf(stdout, "%.8e, %.8e, %.8e, %.8e, %.8e, %.8e\n", boxlo[0], boxlo[1], boxlo[2], boxhi[0], boxhi[1], boxhi[2]);

  //fprintf(stdout, "x=%.4e.\n", atom->x[0][0]);  
  vol_ = domain->xprd * domain->yprd * domain->zprd;
  
/*
  // Make the initial force zero so that it matches i-PI. 
  // Only for debug purpose.
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  tagint *tag = atom->tag;

  // printf("start qc_step, x:\n");
  for(int i=0; i<nlocal; i++)
  {
    // printf("%ld  ", tag[i]);
    for(int j=0; j<3; j++)
    {
      // printf("%.8e  ", x[i][j]);
      f[i][j] = 0.0;
    }
    // printf("\n");
  }
  // printf("\n");
*/
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::initial_integrate(int /*vflag*/)
{
  // unmap the atom coordinates and image flags so that the ring polymer is not wrapped
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  for(int i=0; i<nlocal; i++)
  {
    domain->unmap(x[i], image[i]);
  }
  
  if(integrator==baoab)
  {
    if(pextflag) press_v_step();
    b_step();
    if(method==NMPIMD)
    {
      nmpimd_fill(atom->x);
      comm_exec(atom->x);
      nmpimd_transform(buf_beads, atom->x, M_x2xp[universe->iworld]);
    }
  }

  else
  {
    error->universe_all(FLERR,"Unknown integrator parameter for fix pimd");
  }
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::post_integrate()
{
  if(integrator==baoab)
  {
    qc_step();
    a_step();
    if(ensemble==nvt || ensemble==npt)
    {
      o_step();
      if(pextflag) press_o_step();
    }
    else if(ensemble==nve || ensemble==nph)
    {
    
    }
    else
    {
      error->universe_all(FLERR, "Unknown ensemble parameter for fix pimd. Only nve and nvt are supported!\n");
    }
    qc_step();
    a_step();

    compute_spring_energy();

    if(method==NMPIMD)
    {
      nmpimd_fill(atom->x);
      comm_exec(atom->x);
      nmpimd_transform(buf_beads, atom->x, M_xp2x[universe->iworld]);
    }

    int nlocal = atom->nlocal;
    double **x = atom->x;
    imageint *image = atom->image;

    // remap the atom coordinates and image flags so that all the atoms are in the box and domain->pbc() does not change their coordinates
    for(int i=0; i<nlocal; i++)
    {
      domain->remap(x[i], image[i]);
    }
  
  }

  else
  {
    error->universe_all(FLERR, "Unknown integrator parameter for fix pimd");
  }
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::final_integrate()
{
  if(integrator==baoab)
  {
    // compute_totke();
    compute_p_cv();
    if(pextflag) press_v_step();
    b_step();
  }
  else
  {
    error->universe_all(FLERR,"Unknown integrator parameter for fix pimd");
  }
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::post_force(int /*flag*/)
{
  // unmap the atom coordinates and image flags so that the ring polymer is not wrapped
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  for(int i=0; i<nlocal; i++)
  {
    domain->unmap(x[i], image[i]);
  }

  inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  comm_exec(atom->x);
  comm_coords();
  comm_forces();
  compute_xc();
  compute_fc();
  compute_vir();
  compute_vir_();
  compute_t_prim();
  compute_t_vir();
  compute_pote();
  

  // transform the force into normal mode representation
  if(method==NMPIMD)
  {
    nmpimd_fill(atom->f);
    comm_exec(atom->f);
    nmpimd_transform(buf_beads, atom->f, M_x2xp[universe->iworld]);
  }
  c_pe->addstep(update->ntimestep+1); 
  c_press->addstep(update->ntimestep+1); 
}

/* ----------------------------------------------------------------------
   Langevin thermostat, BAOAB integrator
------------------------------------------------------------------------- */

void FixDPPimd::baoab_init()
{
  //fprintf(stdout, "baoab_temp=%.2f.\n", baoab_temp);
  double KT = force->boltz * baoab_temp;
  double beta = 1.0 / KT;
  double hbar = force->hplanck / (2.0 * MY_PI);
  _omega_np = np / beta / hbar;
  double _omega_np_dt_half = _omega_np * update->dt * 0.5;

  _omega_k = new double[np];
  baoab_c = new double[np];
  baoab_s = new double[np];
  if(fmmode==physical){
    for (int i=0; i<np; i++)
    {
      _omega_k[i] = _omega_np * sqrt(lam[i]); 
      baoab_c[i] = cos(sqrt(lam[i])*_omega_np_dt_half);
      baoab_s[i] = sin(sqrt(lam[i])*_omega_np_dt_half);
    }
  }
  else if(fmmode==normal){
    for (int i=0; i<np; i++)
    {
      _omega_k[i] = _omega_np; 
      baoab_c[i] = cos(_omega_np_dt_half);
      baoab_s[i] = sin(_omega_np_dt_half);
    }
  }
  if(tau > 0) gamma = 1.0 / tau;
  else gamma = np / beta / hbar;
  c1 = exp(-gamma * update->dt); // tau is the damping time of the centroid mode.
  c2 = sqrt(1.0 - c1 * c1); // note that c1 and c2 here only works for the centroid mode.

  if(thermostat == PILE_L || thermostat == PILE_G)
  {
    std::string out = "\nInitializing PI Langevin equation thermostat...\n";
    out += "Bead ID    |    omega    |    tau    |    c1    |    c2\n"; 
    //if(universe->iworld==0) fprintf(stdout, "Initializing PILE_L thermostat.\n");
    tau_k = new double[np];
    c1_k = new double[np];
    c2_k = new double[np];
    tau_k[0] = tau; c1_k[0] = c1; c2_k[0] = c2;
    for(int i=1; i<np; i++)
    {
      tau_k[i] = 0.5 / pilescale / _omega_k[i];
      c1_k[i] = exp(-1.0 * update->dt / tau_k[i]);
      c2_k[i] = sqrt(1.0 - c1_k[i] * c1_k[i]);
    }
    for(int i=0; i<np; i++)
    {
      out += fmt::format("    {:d}     {:.8e} {:.8e} {:.8e} {:.8e}\n", i, _omega_k[i], tau_k[i], c1_k[i], c2_k[i]);
    }
    if(thermostat == PILE_L) out += "PILE_L thermostat successfully initialized!\n";
    else if(thermostat == PILE_G) out += "PILE_G thermostat successfully initialized!\n";
    out += "\n";
    utils::logmesg(lmp, out);
  }


  baoab_ready = true;
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::b_step()
{
  int n = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for(int i=0; i<n; i++)
  {
    double dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }

}

/* ---------------------------------------------------------------------- */

void FixDPPimd::qc_step(){
  printf("\nstart qc_step, vol = %.8e h = (%.8e %.8e %.8e)\n", domain->xprd*domain->yprd*domain->zprd, domain->xprd, domain->yprd, domain->zprd);
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  tagint *tag = atom->tag;

  printf("start qc_step, x:\n");
  for(int i=0; i<nlocal; i++)
  {
    printf("%ld  ", tag[i]);
    for(int j=0; j<3; j++)
    {
      printf("%.8e  ", x[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  if(!pextflag) {
    if(universe->iworld == 0)
    {

      //fprintf(stdout, "executing qc_step, iworld=%ld.\n", universe->iworld);
      for(int i=0; i<nlocal; i++)
      {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
      //fprintf(stdout, "iworld=%lld, x=%.6f.\n", universe->iworld, x[0][0]);
    }
  }
  else if(pextflag) {
    double expq = exp(dtv * vw);
    printf("\nin qc_step, expq = %.15e\n\n", expq);
    //double expv = exp(-(1. + 1./atom->natoms) * dtv * vw);
    double expv = exp(-dtv * vw);
    if(universe->iworld == 0)
    {

      printf("in qc_step, v:\n");
      for(int i=0; i<nlocal; i++)
      {
        printf("%ld  ", tag[i]);
        for(int j=0; j<3; j++)
        {
          printf("%.8e  ", v[i][j]);
        }
        printf("\n");
      }
      printf("\n");

      for(int i=0; i<nlocal; i++)
      {
        for(int j=0; j<3; j++)
        {
          x[i][j] = expq * x[i][j] + (expq - expv) / 2. / vw * v[i][j];
          v[i][j] = expv * v[i][j];
        } 
      }
    }
    domain->xprd *= expq;
    domain->yprd *= expq;
    domain->zprd *= expq;

    domain->boxlo[0] = -0.5*domain->xprd;
    domain->boxlo[1] = -0.5*domain->yprd;
    domain->boxlo[2] = -0.5*domain->zprd;
    domain->boxhi[0] = 0.5*domain->xprd;
    domain->boxhi[1] = 0.5*domain->yprd;
    domain->boxhi[2] = 0.5*domain->zprd;

    domain->set_global_box();
    domain->set_local_box();

  }

  printf("end qc_step, vol = %.8e h = (%.8e %.8e %.8e)\n\n", domain->xprd*domain->yprd*domain->zprd, domain->xprd, domain->yprd, domain->zprd);
  printf("end qc_step, x:\n");
  for(int i=0; i<nlocal; i++)
  {
    printf("%ld  ", tag[i]);
    for(int j=0; j<3; j++)
    {
      printf("%.8e  ", x[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void FixDPPimd::a_step(){
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double x0, x1, x2, v0, v1, v2; // three components of x[i] and v[i]

  if(universe->iworld != 0)
  {
    for(int i=0; i<n; i++)
    {
      x0 = x[i][0]; x1 = x[i][1]; x2 = x[i][2];
      v0 = v[i][0]; v1 = v[i][1]; v2 = v[i][2];
      x[i][0] = baoab_c[universe->iworld] * x0 + 1./_omega_k[universe->iworld] * baoab_s[universe->iworld] * v0;
      x[i][1] = baoab_c[universe->iworld] * x1 + 1./_omega_k[universe->iworld] * baoab_s[universe->iworld] * v1;
      x[i][2] = baoab_c[universe->iworld] * x2 + 1./_omega_k[universe->iworld] * baoab_s[universe->iworld] * v2;
      v[i][0] = -1.*_omega_k[universe->iworld] * baoab_s[universe->iworld] * x0 + baoab_c[universe->iworld] * v0;
      v[i][1] = -1.*_omega_k[universe->iworld] * baoab_s[universe->iworld] * x1 + baoab_c[universe->iworld] * v1;
      v[i][2] = -1.*_omega_k[universe->iworld] * baoab_s[universe->iworld] * x2 + baoab_c[universe->iworld] * v2;
    }
  }

}

/* ---------------------------------------------------------------------- */
void FixDPPimd::svr_step(MPI_Comm which)
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double beta_np = 1.0 / force->boltz / baoab_temp / np * force->mvv2e;

  // compute bead kinetic energy
  double ke_0 = 0.0, ke_total = 0.0;
  for(int i=0; i<nlocal; i++) for(int j=0; j<3; j++) ke_0 += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j];
  MPI_Allreduce(&ke_0, &ke_total, 1, MPI_DOUBLE, MPI_SUM, which);

  // compute alpha
  double noise_ = 0.0, noise_total = 0.0, ksi0_ = 0.0, ksi_ = 0.0;
  for(int i=0; i<atom->natoms; i++) 
  {
    for(int j=0; j<3; j++) 
    {
      ksi_ = random->gaussian();
      if(i==0 && j==0 && universe->iworld==0) ksi0_ = ksi_;
      noise_ += ksi_ * ksi_;
    }
  }
  MPI_Allreduce(&noise_, &noise_total, 1, MPI_DOUBLE, MPI_SUM, which);
  //MPI_Bcast(&ksi0_, 1, MPI_DOUBLE, 0, which);
  
  if(universe->me == 0)
  {
    alpha2 = c1 + (1.0 - c1) * (noise_total) / 2 / beta_np / ke_total + 2 * ksi0_ * sqrt(c1 * (1.0 - c1) / 2 / beta_np / ke_total);
    sgn_ = ksi0_ + sqrt(2 * beta_np * ke_total * c1 / (1.0 - c1));
    // sgn = sgn_ / abs(sgn_);
    if(sgn_<0) sgn = -1.0;
    else sgn = 1.0;
    alpha = sgn * sqrt(alpha2);
  }
  //fprintf(stdout, "iworld = %d, ke_total = %.6e, ksi0_ = %.6e, noise_total = %.6e, c1 = %.6e, alpha2 = %.6e.\n", universe->iworld, ke_total, ksi0_, noise_total, c1, alpha2);

  // broadcast alpha to the other processes in this world world
  MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, which);

  //fprintf(stdout, "iworld = %d, me = %d, alpha = %.6e.\n", universe->iworld, universe->me, alpha);

  // scale the velocities
  for(int i=0; i<nlocal; i++)
  {
    for(int j=0; j<3; j++)
    {
      atom->v[i][j] *= alpha;
    }
  }

}

void FixDPPimd::press_v_step()
{
  printf("iworld = %d, start press_v_step, pw = %.8e.\n", universe->iworld, vw*W);
  int nlocal = atom->nlocal;
  double **f = atom->f;
  double **v = atom->v;
  int *type = atom->type; 
  double volume = domain->xprd * domain->yprd * domain->zprd;
  vw += dtv * 3 * (volume * np * (p_cv - Pext) + Vcoeff / beta_np) / W;
  printf("iworld = %d, p_cv = %.6e, Pext = %.6e, beta_np = %.6e, W = %.6e.\n", universe->iworld, p_cv, Pext, beta_np, W);
  printf("iworld = %d, after adding kinetic part, pw = %.8e.\n", universe->iworld, vw*W);
  if(universe->iworld==0){
    double dvw = 0.0;
    for(int i = 0; i < nlocal; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        dvw += dtv2 * f[i][j] * v[i][j] / W + dtv3 * f[i][j] * f[i][j] / mass[type[i]] / W;
      }
    }
    MPI_Allreduce(&dvw, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
    vw += dvw;
  }
  MPI_Bcast(&vw, 1, MPI_DOUBLE, 0, universe->uworld);
  printf("iworld = %d, ending press_v_step, pw = %.8e.\n", universe->iworld, vw*W);
}


void FixDPPimd::press_o_step()
{
  r1 = random->gaussian();
  vw = c1 * vw + c2 * sqrt(1. / W / beta_np) * r1;
}

void FixDPPimd::o_step()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double beta_np = 1.0 / force->boltz / baoab_temp / np * force->mvv2e;
  if(thermostat == PILE_L)
  {
    for(int i=0; i<nlocal; i++)
    {
      r1 = random->gaussian();
      r2 = random->gaussian();
      r3 = random->gaussian();
      //r1 = r2 = r3 = 0.0;
      //char* rns;
      //sprintf(rns, "%.6e %.6e %.6e\n", r1, r2, r3); 
      //fwrite(rns, sizeof(char), sizeof(rns), frand);
      //fprintf(frand, "%ld %d %.6e %.6e %.6e\n", update->ntimestep, i, r1, r2, r3);
      atom->v[i][0] = c1_k[universe->iworld] * atom->v[i][0] + c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r1; 
      atom->v[i][1] = c1_k[universe->iworld] * atom->v[i][1] + c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r2;
      atom->v[i][2] = c1_k[universe->iworld] * atom->v[i][2] + c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r3;
    }
      //fprintf(stdout, "iworld = %d, after o, v = %.8e.\n", universe->iworld, atom->v[0][0]);
  }
  else if(thermostat == SVR)
  {
    svr_step(universe->uworld);
  }
  else if(thermostat == PILE_G)
  {
    if(universe->iworld == 0)
    {
      svr_step(world);
    }
    else
    {
      for(int i=0; i<nlocal; i++)
      {
        r1 = random->gaussian();
        r2 = random->gaussian();
        r3 = random->gaussian();
        atom->v[i][0] = c1_k[universe->iworld] * atom->v[i][0] + c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r1; 
        atom->v[i][1] = c1_k[universe->iworld] * atom->v[i][1] + c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r2;
        atom->v[i][2] = c1_k[universe->iworld] * atom->v[i][2] + c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r3;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Normal Mode PIMD
------------------------------------------------------------------------- */

void FixDPPimd::nmpimd_init()
{
  memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");

  lam = (double*) memory->smalloc(sizeof(double)*np, "FixDPPimd::lam");

  // Set up  eigenvalues

  lam[0] = 0.0;
  if(np%2==0) lam[np-1] = 4.0;

  for(int i=2; i<=np/2; i++)
  {
    lam[2*i-3] = lam[2*i-2] = 2.0 * (1.0 - 1.0 *cos(2.0*MY_PI*(i-1)/np));
  }

  // Set up eigenvectors for non-degenerated modes

  for(int i=0; i<np; i++)
  {
    M_x2xp[0][i] = 1.0 / sqrt(np);
    if(np%2==0) M_x2xp[np-1][i] = 1.0 / sqrt(np) * pow(-1.0, i);
  }

  // Set up eigenvectors for degenerated modes

  for(int i=0; i<(np-1)/2; i++) for(int j=0; j<np; j++)
  {
    M_x2xp[2*i+1][j] =   sqrt(2.0) * cos ( 2.0 * MY_PI * (i+1) * j / np) / sqrt(np);
    M_x2xp[2*i+2][j] = - sqrt(2.0) * sin ( 2.0 * MY_PI * (i+1) * j / np) / sqrt(np);
  }

  // Set up Ut

  for(int i=0; i<np; i++)
    for(int j=0; j<np; j++)
    {
      M_xp2x[i][j] = M_x2xp[j][i];
    }

  // Set up masses

  int iworld = universe->iworld;

  for(int i=1; i<=atom->ntypes; i++)
  {
    mass[i] = atom->mass[i];

    if(iworld)
    {
      if(fmmode==physical) { mass[i] *= 1.0; }
      else if(fmmode==normal) { mass[i] *= lam[iworld]; }
//      mass[i] *= lam[iworld];
      mass[i] *= fmass;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::nmpimd_fill(double **ptr)
{
  comm_ptr = ptr;
  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::nmpimd_transform(double** src, double** des, double *vector)
{
  int n = atom->nlocal;
  int m = 0;

  //fprintf(stdout, "starting, src=%.6f %.6f, des=%.6f %.6f, vec=%.6f %.6f\n", src[0][0], src[1][0], des[0][0], des[1][0], vector[0], vector[1]);

  for(int i=0; i<n; i++) for(int d=0; d<3; d++)
  {
    des[i][d] = 0.0;
    for(int j=0; j<np; j++) { des[i][d] += (src[j][m] * vector[j]); }
    m++;
  }
  //fprintf(stdout, "ending, src=%.6f %.6f, des=%.6f %.6f, vec=%.6f %.6f\n", src[0][0], src[1][0], des[0][0], des[1][0], vector[0], vector[1]);
}

/* ----------------------------------------------------------------------
   Comm operations
------------------------------------------------------------------------- */

void FixDPPimd::comm_init()
{
  if(size_plan)
  {
    delete [] plan_send;
    delete [] plan_recv;
  }

  if(method == PIMD)
  {
    size_plan = 2;
    plan_send = new int [2];
    plan_recv = new int [2];
    mode_index = new int [2];

    int rank_last = universe->me - comm->nprocs;
    int rank_next = universe->me + comm->nprocs;
    if(rank_last<0) rank_last += universe->nprocs;
    if(rank_next>=universe->nprocs) rank_next -= universe->nprocs;

    plan_send[0] = rank_next; plan_send[1] = rank_last;
    plan_recv[0] = rank_last; plan_recv[1] = rank_next;

    mode_index[0] = 0; mode_index[1] = 1;
    x_last = 1; x_next = 0;
  }
  else
  {
    size_plan = np - 1;
    plan_send = new int [size_plan];
    plan_recv = new int [size_plan];
    mode_index = new int [size_plan];

    for(int i=0; i<size_plan; i++)
    {
      plan_send[i] = universe->me + comm->nprocs * (i+1);
      if(plan_send[i]>=universe->nprocs) plan_send[i] -= universe->nprocs;

      plan_recv[i] = universe->me - comm->nprocs * (i+1);
      if(plan_recv[i]<0) plan_recv[i] += universe->nprocs;

      mode_index[i]=(universe->iworld+i+1)%(universe->nworlds);
    }

    x_next = (universe->iworld+1+universe->nworlds)%(universe->nworlds);
    x_last = (universe->iworld-1+universe->nworlds)%(universe->nworlds);
  }

  if(buf_beads)
  {
    for(int i=0; i<np; i++) if(buf_beads[i]) delete [] buf_beads[i];
    delete [] buf_beads;
  }

  buf_beads = new double* [np];
  for(int i=0; i<np; i++) buf_beads[i] = nullptr;
  
  if(coords)
  {
    for(int i=0; i<np; i++) if(coords[i]) delete [] coords[i];
    delete [] coords;
  }

  if(forces)
  {
    for(int i=0; i<np; i++) if(forces[i]) delete [] forces[i];
    delete [] forces;
  }
  
  coords = new double* [np];
  for(int i=0; i<np; i++) coords[i] = nullptr;

  forces = new double* [np];
  for(int i=0; i<np; i++) forces[i] = nullptr;
  
  if(x_scaled)
  {
    for(int i=0; i<np; i++) if(x_scaled[i]) delete [] x_scaled[i];
    delete [] x_scaled;
  }

  x_scaled = new double* [np];
  for(int i=0; i<np; i++) x_scaled[i] = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::comm_exec(double **ptr)
{
  int nlocal = atom->nlocal;

  if(nlocal > max_nlocal)
  {
    max_nlocal = nlocal+200;
    int size = sizeof(double) * max_nlocal * 3;
    buf_recv = (double*) memory->srealloc(buf_recv, size, "FixDPPimd:x_recv");

    for(int i=0; i<np; i++)
      buf_beads[i] = (double*) memory->srealloc(buf_beads[i], size, "FixDPPimd:x_beads[i]");
  }

  // copy local positions

  memcpy(buf_beads[universe->iworld], &(ptr[0][0]), sizeof(double)*nlocal*3);

  // go over comm plans

  for(int iplan = 0; iplan<size_plan; iplan++)
  {
    // sendrecv nlocal

    int nsend;

    MPI_Sendrecv( &(nlocal), 1, MPI_INT, plan_send[iplan], 0,
                  &(nsend),  1, MPI_INT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // allocate arrays

    if(nsend > max_nsend)
    {
      max_nsend = nsend+200;
      tag_send = (tagint*) memory->srealloc(tag_send, sizeof(tagint)*max_nsend, "FixDPPimd:tag_send");
      buf_send = (double*) memory->srealloc(buf_send, sizeof(double)*max_nsend*3, "FixDPPimd:x_send");
    }

    // send tags

    MPI_Sendrecv( atom->tag, nlocal, MPI_LMP_TAGINT, plan_send[iplan], 0,
                  tag_send,  nsend,  MPI_LMP_TAGINT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // wrap positions

    double *wrap_ptr = buf_send;
    int ncpy = sizeof(double)*3;

    for(int i=0; i<nsend; i++)
    {
      int index = atom->map(tag_send[i]);

      if(index<0)
      {
        char error_line[256];

        sprintf(error_line, "Atom " TAGINT_FORMAT " is missing at world [%d] "
                "rank [%d] required by  rank [%d] (" TAGINT_FORMAT ", "
                TAGINT_FORMAT ", " TAGINT_FORMAT ").\n", tag_send[i],
                universe->iworld, comm->me, plan_recv[iplan],
                atom->tag[0], atom->tag[1], atom->tag[2]);

        error->universe_one(FLERR,error_line);
      }

      memcpy(wrap_ptr, ptr[index], ncpy);
      wrap_ptr += 3;
    }

    // sendrecv x

    MPI_Sendrecv( buf_send, nsend*3,  MPI_DOUBLE, plan_recv[iplan], 0,
                  buf_recv, nlocal*3, MPI_DOUBLE, plan_send[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // copy x

    memcpy(buf_beads[mode_index[iplan]], buf_recv, sizeof(double)*nlocal*3);
  }
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::comm_coords()
{
  int nlocal = atom->nlocal;

  // assign memory for arrays
  int size_coords = sizeof(double) * nlocal * 3;
  int size_tags;// = sizeof(tagint) * nlocal;
  coords_recv = (double*) memory->srealloc(coords_recv, size_coords, "FixDPPimd:coords_recv");
  for(int i=0; i<np; i++)
  {
    coords[i] = (double*) memory->srealloc(coords[i], size_coords, "FixDPPimd:coords[i]");
  }
  
  // copy local positions and tags
  memcpy(coords[universe->iworld], &(atom->x[0][0]), size_coords);

  // traversing over all the other worlds
  for(int dworld=1; dworld<=np-1; dworld++)
  {
      // send the tags and coords to the process proc_send
      // receive the tags and coords from the process proc_recv
    int proc_send = (universe->me + dworld * comm->nprocs) % universe->nprocs; 
    int proc_recv = (universe->me - dworld * comm->nprocs + universe->nprocs) % universe->nprocs;
    int world_recv = (int)(proc_recv / comm->nprocs);
    
    // determine the number of atoms to be sent to and received from the other worlds
    MPI_Sendrecv(&(nlocal), 1, MPI_INT, proc_send, 0, 
                 &(nsend), 1, MPI_INT, proc_recv, 0, 
                 universe->uworld, MPI_STATUS_IGNORE);
    nrecv = nlocal;

    size_coords = sizeof(double) * nsend * 3;
    size_tags = sizeof(tagint) * nsend;

    coords_send = (double*) memory->srealloc(coords_send, size_coords, "FixDPPimd:coords_send");
    tags_send = (tagint*) memory->srealloc(tags_send, size_tags, "FixDPPimd:tags_send");

    MPI_Sendrecv(atom->tag, nlocal, MPI_LMP_TAGINT, proc_send, 0,
                 tags_send, nsend, MPI_LMP_TAGINT, proc_recv, 0,
                 universe->uworld, MPI_STATUS_IGNORE);

    // wrap positions
    double *wrap_ptr = coords_send;
    int ncpy = sizeof(double)*3;

    for(int i=0; i<nsend; i++)
    {
      int index = atom->map(tags_send[i]);
      if(index < 0)
      {
        char error_line[256];

        sprintf(error_line, "Atom " TAGINT_FORMAT " is missing at world [%d] "
                "rank [%d] required by  rank [%d] (" TAGINT_FORMAT ", "
                TAGINT_FORMAT ", " TAGINT_FORMAT ").\n", tags_send[i],
                universe->iworld, comm->me, proc_recv,
                atom->tag[0], atom->tag[1], atom->tag[2]);

        error->universe_one(FLERR,error_line);
      }    
      memcpy(wrap_ptr, atom->x[index], ncpy);
      wrap_ptr += 3;
    }
    MPI_Sendrecv(coords_send, nsend*3, MPI_DOUBLE, proc_recv, 0,
                coords_recv, nrecv*3, MPI_DOUBLE, proc_send, 0,
                universe->uworld, MPI_STATUS_IGNORE);

    memcpy(coords[world_recv], coords_recv, sizeof(double)*nlocal*3);          
  }
}

void FixDPPimd::comm_forces()
{
  int nlocal = atom->nlocal;

  // assign memory for arrays
  int size_forces = sizeof(double) * nlocal * 3;
  int size_tags;// = sizeof(tagint) * nlocal;
  forces_recv = (double*) memory->srealloc(forces_recv, size_forces, "FixDPPimd:forces_recv");
  for(int i=0; i<np; i++)
  {
    forces[i] = (double*) memory->srealloc(forces[i], size_forces, "FixDPPimd:forces[i]");
  }
  
  // copy local positions and tags
  memcpy(forces[universe->iworld], &(atom->f[0][0]), size_forces);

  // traversing over all the other worlds
  for(int dworld=1; dworld<=np-1; dworld++)
  {
      // send the tags and forces to the process proc_send
      // receive the tags and forces from the process proc_recv
    int proc_send = (universe->me + dworld * comm->nprocs) % universe->nprocs; 
    int proc_recv = (universe->me - dworld * comm->nprocs + universe->nprocs) % universe->nprocs;
    int world_recv = (int)(proc_recv / comm->nprocs);
    
    // determine the number of atoms to be sent to and received from the other worlds
    MPI_Sendrecv(&(nlocal), 1, MPI_INT, proc_send, 0, 
                 &(nsend), 1, MPI_INT, proc_recv, 0, 
                 universe->uworld, MPI_STATUS_IGNORE);
    nrecv = nlocal;

    size_forces = sizeof(double) * nsend * 3;
    size_tags = sizeof(tagint) * nsend;

    forces_send = (double*) memory->srealloc(forces_send, size_forces, "FixDPPimd:forces_send");
    tags_send = (tagint*) memory->srealloc(tags_send, size_tags, "FixDPPimd:tags_send");

    MPI_Sendrecv(atom->tag, nlocal, MPI_LMP_TAGINT, proc_send, 0,
                 tags_send, nsend, MPI_LMP_TAGINT, proc_recv, 0,
                 universe->uworld, MPI_STATUS_IGNORE);

    // wrap positions
    double *wrap_ptr = forces_send;
    int ncpy = sizeof(double)*3;

    for(int i=0; i<nsend; i++)
    {
      int index = atom->map(tags_send[i]);
      if(index < 0)
      {
        char error_line[256];

        sprintf(error_line, "Atom " TAGINT_FORMAT " is missing at world [%d] "
                "rank [%d] required by  rank [%d] (" TAGINT_FORMAT ", "
                TAGINT_FORMAT ", " TAGINT_FORMAT ").\n", tags_send[i],
                universe->iworld, comm->me, proc_recv,
                atom->tag[0], atom->tag[1], atom->tag[2]);

        error->universe_one(FLERR,error_line);
      }    
      memcpy(wrap_ptr, atom->f[index], ncpy);
      wrap_ptr += 3;
    }
    MPI_Sendrecv(forces_send, nsend*3, MPI_DOUBLE, proc_recv, 0,
                forces_recv, nrecv*3, MPI_DOUBLE, proc_send, 0,
                universe->uworld, MPI_STATUS_IGNORE);

    memcpy(forces[world_recv], forces_recv, sizeof(double)*nlocal*3);          
  }
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::compute_xc()
{
  int nlocal = atom->nlocal;
  xc = (double*) memory->srealloc(xc, sizeof(double) * nlocal * 3, "FixDPPimd:xc");
  for(int i=0; i<nlocal; i++)
  {
    xc[3*i] = xc[3*i+1] = xc[3*i+2] = 0.0;
    for(int j=0; j<np; j++)
    {
      xc[3*i] += coords[j][3*i];
      xc[3*i+1] += coords[j][3*i+1];
      xc[3*i+2] += coords[j][3*i+2];
    }
    xc[3*i] /= np;
    xc[3*i+1] /= np;
    xc[3*i+2] /= np;
  } 
}

void FixDPPimd::compute_fc()
{
  int nlocal = atom->nlocal;
  fc = (double*) memory->srealloc(fc, sizeof(double) * nlocal * 3, "FixDPPimd:fc");
  for(int i=0; i<nlocal; i++)
  {
    fc[3*i] = fc[3*i+1] = fc[3*i+2] = 0.0;
    for(int j=0; j<np; j++)
    {
      fc[3*i] += forces[j][3*i];
      fc[3*i+1] += forces[j][3*i+1];
      fc[3*i+2] += forces[j][3*i+2];
    }
    // fc[3*i] /= np;
    // fc[3*i+1] /= np;
    // fc[3*i+2] /= np;
  } 
}

void FixDPPimd::compute_vir_()
{
  int nlocal = atom->nlocal;
  xf = vir_ = xcf = centroid_vir = 0.0;
  for(int i=0; i<nlocal; i++)
  {
    for(int j=0; j<3; j++)
    {
      xf += atom->x[i][j] * atom->f[i][j];
      xcf += (atom->x[i][j] - xc[3*i+j]) * atom->f[i][j];
    }
  }
  MPI_Allreduce(&xf, &vir_, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&xcf, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

void FixDPPimd::compute_vir()
{
  double volume = domain->xprd * domain->yprd * domain->zprd;
  c_press->compute_vector();
  virial[0] = c_press->vector[0]*volume;
  virial[4] = c_press->vector[1]*volume;
  virial[8] = c_press->vector[2]*volume;
  virial[1] = c_press->vector[3]*volume;
  virial[2] = c_press->vector[4]*volume;
  virial[5] = c_press->vector[5]*volume;
  //int nlocal = atom->nlocal;  
  //xf = vir = xcfc = centroid_vir = 0.0;
  //for(int i=0; i<nlocal; i++)
  //{
    //for(int j=0; j<3; j++)
    //{
      //xf += atom->x[i][j] * atom->f[i][j];
      //xcfc += xc[3*i+j] * fc[3*i+j];
    //}
  //}

  //MPI_Allreduce(&xcfc, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, world);
  printf("computing vir, vir:\n");
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      printf("%.8e  ", virial[3*i+j]);
    }
    printf("\n");
  }
  printf("\n");
  vir=(virial[0]+virial[4]+virial[8]);
  MPI_Allreduce(&vir,&vir,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  //printf("iworld=%d, vir=%.4e.\n", universe->iworld, vir);
}
/* ---------------------------------------------------------------------- */

void FixDPPimd::compute_xscaled()
{
  int nlocal = atom->nlocal;
  for(int i=0; i<np; i++)
  {
    x_scaled[i] = (double*) memory->srealloc(x_scaled[i], sizeof(double) * nlocal * 3, "FixDPPimd:x_scaled[i]");
  }
  for(int i=0; i<np; i++)
  {
    for(int j=0; j<nlocal; j++)
    {
    x_scaled[i][3*j] = lambda * coords[i][3*j] + (1.0 - lambda) * xc[3*j];
    x_scaled[i][3*j+1] = lambda * coords[i][3*j+1] + (1.0 - lambda) * xc[3*j+1];
    x_scaled[i][3*j+2] = lambda * coords[i][3*j+2] + (1.0 - lambda) * xc[3*j+2];
    }
  }
}

/* ---------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Compute centroid-virial kinetic energy estimator
------------------------------------------------------------------------- */

void FixDPPimd::compute_t_vir()
{
  t_vir = -0.5 / np * vir_;
  t_cv = 1.5 * atom->natoms * force->boltz * temp - 0.5 / np * centroid_vir;
}

/* ----------------------------------------------------------------------
   Compute primitive kinetic energy estimator
------------------------------------------------------------------------- */

void FixDPPimd::compute_t_prim()
{
  // fprintf(stdout, "in compute_t_prim, me = %d, N = %d, np = %d, force->boltz = %2.8f, temp = %2.8f, total_spring_energy = %2.8e.\n", universe->me, atom->natoms, np, force->boltz, temp, total_spring_energy);
  t_prim = 1.5 * atom->natoms * np * force->boltz * temp - total_spring_energy;
}

void FixDPPimd::compute_p_prim()
{
  //p_prim = atom->natoms * force->boltz * temp * inv_volume - 1.0 / 1.5 * inv_volume * total_spring_energy;
  //p_prim = atom->natoms * force->boltz * temp * inv_volume - 1.0 / 1.5 * inv_volume * total_spring_energy + 1.0 / 3 / np * inv_volume * vir;
  p_prim = atom->natoms * np * force->boltz * temp * inv_volume - 1.0 / 1.5 * inv_volume * total_spring_energy;
}

void FixDPPimd::compute_p_cv()
{
  //p_cv = 2. / 3.  * inv_volume / np * totke - 1. / 3. / np * inv_volume * centroid_vir; 
  if(universe->iworld == 0)
  {
    p_cv = 1. / 3.  * inv_volume  * (2. * ke_bead - 1. * centroid_vir + 1. * vir) / force->nktv2p / np; 
  }
  MPI_Bcast(&p_cv, 1, MPI_DOUBLE, 0, universe->uworld);
  fprintf(stdout, "in compute_p_cv, iworld = %d, ke_bead = %.8e, centroid_vir = %.8e, vir = %.8e, p_cv = %.8e.\n", universe->iworld, ke_bead, centroid_vir, vir, p_cv);
  // fprintf(stdout, "iworld = %d.\n", universe->iworld);
}

void FixDPPimd::compute_p_vir()
{
  //inv_volume = 1. / (domain->xprd * domain->yprd * domain->zprd);
  //inv_volume = 1. / vol_;
  //p_vir = 2.  / 3 * inv_volume * totke + 1. / 3 * inv_volume * vir * force->nktv2p;
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::compute_totke()
{
  double kine = 0.0;
  totke = 0.0;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  //double *_mass = atom->mass;
  for(int i=0; i<nlocal; i++)
  {
    for(int j=0; j<3; j++)
    {
      kine += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j];
    }
  }
  MPI_Allreduce(&kine, &ke_bead, 1, MPI_DOUBLE, MPI_SUM, world);
  ke_bead *= force->mvv2e;
  MPI_Allreduce(&kine, &totke, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  totke *= force->mvv2e / np;

  c_press->compute_scalar();
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::compute_pote()
{
  double pot_energy_partition = 0.0;
  pote = 0.0;
  c_pe->compute_scalar();
  pot_energy_partition = c_pe->scalar;
  pot_energy_partition /= np;
  MPI_Allreduce(&pot_energy_partition, &pote, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::compute_spring_energy()
{
  spring_energy = 0.0;

  double **x = atom->x;
  double* _mass = atom->mass;
  int* type = atom->type;
  int nlocal = atom->nlocal;

/*
  double* xlast = buf_beads[x_last];
  double* xnext = buf_beads[x_next];

  for(int i=0; i<nlocal; i++)
  {
    double delx1 = xlast[0] - x[i][0];
    double dely1 = xlast[1] - x[i][1];
    double delz1 = xlast[2] - x[i][2];
    xlast += 3;
    domain->minimum_image(delx1, dely1, delz1);

    double delx2 = xnext[0] - x[i][0];
    double dely2 = xnext[1] - x[i][1];
    double delz2 = xnext[2] - x[i][2];
    xnext += 3;
    domain->minimum_image(delx2, dely2, delz2);

    double ff = fbond * _mass[type[i]];

    double dx = delx1+delx2;
    double dy = dely1+dely2;
    double dz = delz1+delz2;

    spring_energy += -ff * (delx1*delx1+dely1*dely1+delz1*delz1+delx2*delx2+dely2*dely2+delz2*delz2);
  }
  MPI_Allreduce(&spring_energy, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  total_spring_energy *= 0.25;
  total_spring_energy /= np;
*/
  for(int i=0; i<nlocal; i++)
  {
    spring_energy += 0.5 * _mass[type[i]] * fbond * lam[universe->iworld] * (x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2]); 
  }
  //fprintf(stdout, "iworld=%d, _mass=%.2e, fbond=%.2e, lam=%.2e, x=%.2e, se=%.2e.\n", universe->iworld, _mass[type[0]], fbond, lam[universe->iworld], x[0][0], spring_energy);
  MPI_Allreduce(&spring_energy, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  total_spring_energy /= np;
}

/* ---------------------------------------------------------------------- */

void FixDPPimd::compute_tote()
{
  // tote = totke + hope;
  tote = totke + pote + total_spring_energy;
  //printf("totke=%f.\n", totke);
}

void FixDPPimd::compute_totenthalpy()
{
  double volume = domain->xprd * domain->yprd * domain->zprd;
  totenthalpy = tote + 0.5*W*vw*vw + Pext * volume - Vcoeff/beta_np * log(volume);
  //totenthalpy = tote + 0.5*W*vw*vw + Pext * volume ;
  //totenthalpy = tote + 0.5*W*vw*vw + Pext * vol_ ;
  //printf("vol=%f, enth=%f.\n", volume, totenthalpy);
}

/* ---------------------------------------------------------------------- */

double FixDPPimd::compute_vector(int n)
{
  //if(n==0) { return totke; }
  if(n==0) { return ke_bead; }
  if(n==1) { return spring_energy; }
  //if(n==1) { return atom->v[0][0]; }
  if(n==2) { return pote; }
  //if(n==2) { return atom->x[0][0]; }
  //if(n==3) { if(!pextflag) {return tote;} else {return totenthalpy;} }
  //if(n==3) { return totenthalpy; }
  if(n==3) { return tote; }
  if(n==4) { return t_prim; }
  //if(n==4) { printf("returning vol_=%f\n", vol_);  return vol_; }
  //if(n==4) { return vol_; }
  //if(n==5) { return domain->xprd; }
  //if(n==6) { return domain->yprd; }
  if(n==5) { return t_vir; }
  if(n==6) { return t_cv; }
  if(n==7) { return p_prim; }
  if(n==8) { return p_vir; }
  if(n==9) { return p_cv; }
  //if(pextflag) size_vector = 11;
  if(n==10) {return vw;}
  if(n==11) {return totenthalpy;}
  return 0.0;
}

double FixDPPimd::compute_scalar()
{
  return vol_;
}

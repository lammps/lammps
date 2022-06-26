/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Package      FixPIMD
   Purpose      Quantum Path Integral Algorithm for Quantum Chemistry
   Copyright    Voth Group @ University of Chicago
   Authors      Chris Knight & Yuxing Peng (yuxing at uchicago.edu)

   Updated      Oct-01-2011
   Version      1.0

   Updated      Jun-07-2022
   Yifan Li @ Princeton University (yifanl0716@gmail.com)
   Added components:
   - Multi-processor parallelism for each bead
   - White-noise Langevin thermostat
   - Bussi-Zykova-Parrinello barostat (isotropic and anisotropic)
   - Several quantum estimators
   Futher plans:
   - Triclinic barostat
------------------------------------------------------------------------- */

#include "fix_pimd.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "universe.h"
#include "utils.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { PIMD, NMPIMD, CMD };
enum { physical, normal };
enum { baoab, obabo };
enum { ISO, ANISO, TRICLINIC };
enum { PILE_L, NHC };
enum { MTTK, BZP };
// char* Barostats[] = {"MTTK", "BZP"};
std::map<int, std::string> Barostats {{MTTK, "MTTK"}, {BZP, "BZP"}};
enum { nve, nvt, nph, npt };

/* ---------------------------------------------------------------------- */

FixPIMD::FixPIMD(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), random(nullptr), c_pe(nullptr), c_press(nullptr) {
  time_integrate = 1;
  max_nsend = 0;
  tag_send = nullptr;
  buf_send = nullptr;

  max_nlocal = 0;
  buf_recv = nullptr;
  buf_beads = nullptr;

  size_plan = 0;
  plan_send = plan_recv = nullptr;

  M_x2xp = M_xp2x = M_f2fp = M_fp2f = nullptr;
  lam = nullptr;
  mode_index = nullptr;

  mass = nullptr;

  array_atom = nullptr;
  nhc_eta = nullptr;
  nhc_eta_dot = nullptr;
  nhc_eta_dotdot = nullptr;
  nhc_eta_mass = nullptr;

  method = PIMD;
  ensemble      = nvt;
  fmass = 1.0;
  nhc_temp = 298.15;
  nhc_nchain = 2;
  sp = 1.0;
  temp          = 298.15;
  Lan_temp      = 298.15;
  tau           = 1.0;
  tau_p         = 1.0;
  Pext          = 1.0;
  tstat_flag    = 1;
  pstat_flag    = 0;
  mapflag       = 1;
  removecomflag = 1;
  pstyle        = ISO;

  for (int i = 3; i < narg - 1; i += 2) {
    if (strcmp(arg[i], "method") == 0) {
      if (strcmp(arg[i + 1], "pimd") == 0)
        method = PIMD;
      else if (strcmp(arg[i + 1], "nmpimd") == 0)
        method = NMPIMD;
      else if (strcmp(arg[i + 1], "cmd") == 0)
        method = CMD;
      else
        error->universe_all(FLERR, "Unknown method parameter for fix pimd");
    }
    else if(strcmp(arg[i], "integrator")==0)
    {
      if(strcmp(arg[i+1], "obabo")==0) integrator=obabo;
      else if(strcmp(arg[i+1], "baoab")==0) integrator=baoab;
      else error->universe_all(FLERR, "Unknown integrator parameter for fix pimd. Only obabo and baoab integrators are supported!");
    }
    else if(strcmp(arg[i], "ensemble")==0)
    {
      if(strcmp(arg[i+1], "nve")==0) { ensemble = nve; tstat_flag = 0; pstat_flag = 0; }
      else if(strcmp(arg[i+1], "nvt")==0) { ensemble = nvt; tstat_flag = 1; pstat_flag = 0; }
      else if(strcmp(arg[i+1], "nph")==0) { ensemble = nph; tstat_flag = 0; pstat_flag = 1; }
      else if(strcmp(arg[i+1], "npt")==0) { ensemble = npt; tstat_flag = 1; pstat_flag = 1; }
      else error->universe_all(FLERR, "Unknown ensemble parameter for fix pimd. Only nve and nvt ensembles are supported!");
    }
      else if (strcmp(arg[i], "fmass") == 0) {
      fmass = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (fmass < 0.0 || fmass > 1.0)
        error->universe_all(FLERR, "Invalid fmass value for fix pimd");
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
    else if (strcmp(arg[i], "sp") == 0) {
      sp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (fmass < 0.0) error->universe_all(FLERR, "Invalid sp value for fix pimd");
    } else if (strcmp(arg[i], "temp") == 0) {
      nhc_temp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (nhc_temp < 0.0) error->universe_all(FLERR, "Invalid temp value for fix pimd");
    } 
    else if(strcmp(arg[i], "thermostat")==0)
    {
      if(strcmp(arg[i+1],"PILE_L")==0) 
      {
        thermostat = PILE_L;
        seed = atoi(arg[i+2]);
        i++;
      }
      else if(strcmp(arg[i+1],"NHC")==0)
      {
        thermostat = NHC;
      }
    }

    else if(strcmp(arg[i], "tau")==0)
    {
      tau = atof(arg[i+1]);
    }
  
    else if(strcmp(arg[i], "press")==0)
    {
      Pext = atof(arg[i+1]);
      if(Pext<0.0) error->universe_all(FLERR,"Invalid press value for fix pimd");
    }

    else if(strcmp(arg[i], "barostat")==0)
    {
      if(strcmp(arg[i+1],"MTTK")==0) 
      {
        barostat = MTTK;
      }
      else if(strcmp(arg[i+1],"BZP")==0)
      {
        barostat = BZP;
      }
      else error->universe_all(FLERR,"Unknown barostat parameter for fix pimd");
    }

    else if(strcmp(arg[i], "iso")==0)
    {
      pstyle = ISO;
      i--;
    }

    else if(strcmp(arg[i], "aniso")==0)
    {
      pstyle = ANISO;
      i--;
    }

    else if(strcmp(arg[i], "taup")==0)
    {
      tau_p = atof(arg[i+1]);
      if(tau_p<=0.0) error->universe_all(FLERR, "Invalid tau_p value for fix pimd");
    }
    else if(strcmp(arg[i], "fixcom")==0)
    {
      if(strcmp(arg[i+1], "yes")==0) removecomflag = 1;
      else if(strcmp(arg[i+1], "no")==0) removecomflag = 0;
    }

    else if(strcmp(arg[i], "map")==0)
    {
      if(strcmp(arg[i+1], "yes")==0) mapflag = 1;
      else if(strcmp(arg[i+1], "no")==0) mapflag = 0;
    }

    else if (strcmp(arg[i], "nhc") == 0) {
      nhc_nchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
      if (nhc_nchain < 2) error->universe_all(FLERR, "Invalid nhc value for fix pimd");
    } else
      error->universe_all(FLERR, fmt::format("Unknown keyword {} for fix pimd", arg[i]));
  }

  /* Initiation */

  size_peratom_cols = 12 * nhc_nchain + 3;

  nhc_offset_one_1 = 3 * nhc_nchain;
  nhc_offset_one_2 = 3 * nhc_nchain + 3;
  nhc_size_one_1 = sizeof(double) * nhc_offset_one_1;
  nhc_size_one_2 = sizeof(double) * nhc_offset_one_2;

  restart_peratom = 1;
  peratom_flag = 1;
  peratom_freq = 1;

  global_freq = 1;
  vector_flag = 1;
  if(pstyle==ISO) {size_vector = 15;}
  else if(pstyle==ANISO) {size_vector = 23;}
  extvector = 1;
  comm_forward = 3;

  atom->add_callback(Atom::GROW);       // Call LAMMPS to allocate memory for per-atom array
  atom->add_callback(Atom::RESTART);    // Call LAMMPS to re-assign restart-data for per-atom array

  grow_arrays(atom->nmax);

  // some initilizations

  nhc_ready = false;

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

  vol0 = domain->xprd * domain->yprd * domain->zprd;

  fixedpoint[0] = 0.5*(domain->boxlo[0]+domain->boxhi[0]);
  fixedpoint[1] = 0.5*(domain->boxlo[1]+domain->boxhi[1]);
  fixedpoint[2] = 0.5*(domain->boxlo[2]+domain->boxhi[2]);

  // initialize Marsaglia RNG with processor-unique seed

  if(integrator==baoab || integrator==obabo)
  {
    Lan_temp = temp;
    random = new RanMars(lmp, seed + universe->me);
  }
  
}

/* ---------------------------------------------------------------------- */

FixPIMD::~FixPIMD()
{
  // delete[] mass;
  // atom->delete_callback(id, Atom::GROW);
  // atom->delete_callback(id, Atom::RESTART);

  // memory->destroy(M_x2xp);
  // memory->destroy(M_xp2x);
  // memory->destroy(M_f2fp);
  // memory->destroy(M_fp2f);
  // memory->sfree(lam);
  // memory->sfree(buf_beads);

  // delete[] buf_beads;
  // delete[] plan_send;
  // delete[] plan_recv;
  // delete[] mode_index;

  // memory->sfree(tag_send);
  // memory->sfree(buf_send);
  // memory->sfree(buf_recv);

  // memory->destroy(array_atom);
  // memory->destroy(nhc_eta);
  // memory->destroy(nhc_eta_dot);
  // memory->destroy(nhc_eta_dotdot);
  // memory->destroy(nhc_eta_mass);
}

/* ---------------------------------------------------------------------- */
int FixPIMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::init()
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix pimd requires an atom map, see atom_modify");

  if (universe->me == 0 && universe->uscreen)
    fprintf(universe->uscreen, "Fix pimd initializing Path-Integral ...\n");

  // prepare the constants

  masstotal = group->mass(igroup);
  np = universe->nworlds;
  inverse_np = 1.0 / np;

  const double Boltzmann = force->boltz;
  const double Planck = force->hplanck;

  hbar = Planck / (2.0 * MY_PI);
  kBT = force->boltz * temp;
  double beta = 1.0 / (Boltzmann * temp);
  double _fbond = 1.0 * np * np / (beta * beta * hbar * hbar);

  omega_np = np / (hbar * beta) * sqrt(force->mvv2e);
  beta_np = 1.0 / force->boltz / temp / np;
  fbond = -_fbond * force->mvv2e;

  if (universe->me == 0)
    printf("Fix pimd -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", fbond);

  if(integrator==obabo)
  {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = 1./3 * dtv2 * dtv * force->ftm2v;
  }
  else if(integrator==baoab)
  {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = 1./3 * dtv2 * dtv * force->ftm2v;
  }
  else
  {
    error->universe_all(FLERR,"Unknown integrator parameter for fix pimd");
  }

  comm_init();

  mass = new double[atom->ntypes + 1];

  if (method == CMD || method == NMPIMD)
    nmpimd_init();
  else
    for (int i = 1; i <= atom->ntypes; i++) mass[i] = atom->mass[i] / np * fmass;

  if (!nhc_ready && thermostat == NHC) nhc_init();
  Langevin_init();
  if (pstat_flag) baro_init();

  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  int ipress = modify->find_compute(id_press);
  c_press = modify->compute[ipress];

  t_prim = t_vir = t_cv = p_prim = p_vir = p_cv = p_md = 0.0;

  if(universe->me==0) fprintf(screen, "Fix pimd successfully initialized!\n");
}

/* ---------------------------------------------------------------------- */

void FixPIMD::setup(int vflag)
{
  if(universe->me==0) printf("Setting up Path-Integral ...\n");
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  double **x = atom->x;
  imageint *image = atom->image;
  if(mapflag){
    for(int i=0; i<nlocal; i++)
    {
      domain->unmap(x[i], image[i]);
    }
  }
  if(method==NMPIMD)
  {
    nmpimd_fill(atom->x);
    comm_exec(atom->x);
    nmpimd_transform(buf_beads, atom->x, M_x2xp[universe->iworld]);
  }
  compute_spring_energy();
  if(method==NMPIMD)
  {
    nmpimd_fill(atom->x);
    comm_exec(atom->x);
    nmpimd_transform(buf_beads, atom->x, M_xp2x[universe->iworld]);
  }
  if(method==NMPIMD)
  {
    nmpimd_fill(atom->v);
    comm_exec(atom->v);
    nmpimd_transform(buf_beads, atom->v, M_x2xp[universe->iworld]);
  }
  update_x_unwrap();
  compute_xc();
  if(mapflag)
  {
    for(int i=0; i<nlocal; i++)
    {
      domain->unmap_inv(x[i], image[i]);
    }
  }
  post_force(vflag);
  compute_totke();
  compute_pote();
  if(pstyle==ANISO) compute_stress_tensor();
  end_of_step();
  c_pe->addstep(update->ntimestep+1); 
  c_press->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::initial_integrate(int /*vflag*/)
{
  if(ensemble == nvt && thermostat == NHC)
  {
    nhc_update_v();
    nhc_update_x();
  }
  else{
    int nlocal = atom->nlocal;
    tagint *tag = atom->tag;
    double **x = atom->x;
    imageint *image = atom->image;
    if(mapflag){
      for(int i=0; i<nlocal; i++)
      {
        domain->unmap(x[i], image[i]);
      }
    }
    if(integrator==obabo)
    {
      if(tstat_flag)
      {
        o_step();
        if(removecomflag) remove_com_motion();
        if(pstat_flag) press_o_step();
      }
      compute_totke();
      compute_p_cv();
      if(pstat_flag) 
      {
        press_v_step();
      }
      b_step();
      if(removecomflag) remove_com_motion();
      if(method==NMPIMD)
      {
        nmpimd_fill(atom->x);
        comm_exec(atom->x);
        nmpimd_transform(buf_beads, atom->x, M_x2xp[universe->iworld]);
      }
      qc_step();
      a_step();
      qc_step();
      a_step();
    }
    else if(integrator==baoab)
    {
      compute_totke();
      compute_p_cv();
      if(pstat_flag) 
      {
        press_v_step();
      }
      b_step();
      if(removecomflag) remove_com_motion();
      if(method==NMPIMD)
      {
        nmpimd_fill(atom->x);
        comm_exec(atom->x);
        nmpimd_transform(buf_beads, atom->x, M_x2xp[universe->iworld]);
      }
      qc_step();
      a_step();
      if(tstat_flag)
      {
        o_step();
        if(removecomflag) remove_com_motion();
        if(pstat_flag) press_o_step();
      }
      qc_step();
      a_step();      
    }
    else
    {
      error->universe_all(FLERR,"Unknown integrator parameter for fix pimd");
    }
    compute_spring_energy();
    if(method==NMPIMD)
    {
      nmpimd_fill(atom->x);
      comm_exec(atom->x);
      nmpimd_transform(buf_beads, atom->x, M_xp2x[universe->iworld]);
    }
    if(mapflag){
      for(int i=0; i<nlocal; i++)
      {
        domain->unmap_inv(x[i], image[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::final_integrate()
{
  if(ensemble == nvt && thermostat == NHC)
  {
    nhc_update_v();
  }
  else{
    if(pstat_flag) 
    {
      compute_totke();
      compute_p_cv();
      press_v_step();
    }
    b_step();
    if(removecomflag) remove_com_motion();
    if(integrator==obabo)
    {
      if(tstat_flag)
      {
        o_step();
        if(removecomflag) remove_com_motion();
        if(pstat_flag) press_o_step();
      }
    }
    else if(integrator==baoab)
    {
    
    }
    else
    {
      error->universe_all(FLERR,"Unknown integrator parameter for fix pimd");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::post_force(int /*flag*/)
{
  if(ensemble == nvt && thermostat == NHC)
  {
    for (int i = 0; i < atom->nlocal; i++)
      for (int j = 0; j < 3; j++) atom->f[i][j] /= np;

    comm_exec(atom->x);
    spring_force();

    if (method == CMD || method == NMPIMD) {
      /* forward comm for the force on ghost atoms */

      nmpimd_fill(atom->f);

      /* inter-partition comm */

      comm_exec(atom->f);

      /* normal-mode transform */

      nmpimd_transform(buf_beads, atom->f, M_f2fp[universe->iworld]);
    }
  }
  else{
   int nlocal = atom->nlocal;
   tagint *tag = atom->tag;
   double **x = atom->x;
   imageint *image = atom->image;
   if(mapflag){
     for(int i=0; i<nlocal; i++)
     {
       domain->unmap(x[i], image[i]);
     }
   }
   comm_exec(atom->x);
   update_x_unwrap();
   compute_xc();
   if(mapflag)
   {
     for(int i=0; i<nlocal; i++)
     {
       domain->unmap_inv(x[i], image[i]);
     }
   }
   compute_vir();
   compute_vir_();
   compute_t_prim();
   compute_t_vir();
   compute_pote();
   if(method==NMPIMD)
   {
     nmpimd_fill(atom->f);
     comm_exec(atom->f);
     nmpimd_transform(buf_beads, atom->f, M_x2xp[universe->iworld]);
   }
   c_pe->addstep(update->ntimestep+1); 
   c_press->addstep(update->ntimestep+1); 
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::end_of_step()
{
  compute_totke();
  inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  compute_p_prim();
  compute_p_cv();
  compute_tote();
  if(pstat_flag) compute_totenthalpy();

  if(update->ntimestep % 10000 == 0)
  {
  if(universe->me==0) printf("This is the end of step %ld.\n", update->ntimestep);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::update_x_unwrap()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  // delete x_unwrap;
  // memory->sfree(x_unwrap);
  x_unwrap = new double[nlocal*3];
  for(int i=0; i<nlocal; i++)
  { 
    for(int j=0; j<3; j++)
    {
      x_unwrap[3*i+j] = x[i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_xc()
{
  int natoms = atom->natoms;
  comm_exec(atom->x);
  int nlocal = atom->nlocal;
  // if(xc) delete xc;
  xc = new double[nlocal*3];
  for(int i=0; i<nlocal; i++)
  {
    xc[3*i] = xc[3*i+1] = xc[3*i+2] = 0.0;
    for(int j=0; j<np; j++)
    {
      xc[3*i] += buf_beads[j][3*i+0];
      xc[3*i+1] += buf_beads[j][3*i+1];
      xc[3*i+2] += buf_beads[j][3*i+2];
    }
    xc[3*i] /= np;
    xc[3*i+1] /= np;
    xc[3*i+2] /= np;
  } 
}

/* ----------------------------------------------------------------------
   Nose-Hoover Chains
------------------------------------------------------------------------- */

void FixPIMD::nhc_init()
{
  double tau = 1.0 / omega_np;
  double KT = force->boltz * nhc_temp;

  double mass0 = KT * tau * tau;
  int max = 3 * atom->nlocal;

  for (int i = 0; i < max; i++) {
    for (int ichain = 0; ichain < nhc_nchain; ichain++) {
      nhc_eta[i][ichain] = 0.0;
      nhc_eta_dot[i][ichain] = 0.0;
      nhc_eta_dot[i][ichain] = 0.0;
      nhc_eta_dotdot[i][ichain] = 0.0;
      nhc_eta_mass[i][ichain] = mass0;
      if ((method == CMD || method == NMPIMD) && universe->iworld == 0)
        ;
      else
        nhc_eta_mass[i][ichain] *= fmass;
    }

    nhc_eta_dot[i][nhc_nchain] = 0.0;

    for (int ichain = 1; ichain < nhc_nchain; ichain++)
      nhc_eta_dotdot[i][ichain] = (nhc_eta_mass[i][ichain - 1] * nhc_eta_dot[i][ichain - 1] *
                                       nhc_eta_dot[i][ichain - 1] * force->mvv2e -
                                   KT) /
          nhc_eta_mass[i][ichain];
  }

  // Zero NH acceleration for CMD

  if (method == CMD && universe->iworld == 0)
    for (int i = 0; i < max; i++)
      for (int ichain = 0; ichain < nhc_nchain; ichain++) nhc_eta_dotdot[i][ichain] = 0.0;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nhc_update_x()
{
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;

  if (method == CMD || method == NMPIMD) {
    nmpimd_fill(atom->v);
    comm_exec(atom->v);

    /* borrow the space of atom->f to store v in cartisian */

    v = atom->f;
    nmpimd_transform(buf_beads, v, M_xp2x[universe->iworld]);
  }

  for (int i = 0; i < n; i++) {
    x[i][0] += dtv * v[i][0];
    x[i][1] += dtv * v[i][1];
    x[i][2] += dtv * v[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nhc_update_v()
{
  int n = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for (int i = 0; i < n; i++) {
    double dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }

  t_sys = 0.0;
  if (method == CMD && universe->iworld == 0) return;

  double expfac;
  int nmax = 3 * atom->nlocal;
  double KT = force->boltz * nhc_temp;
  double kecurrent, t_current;

  double dthalf = 0.5 * update->dt;
  double dt4 = 0.25 * update->dt;
  double dt8 = 0.125 * update->dt;

  for (int i = 0; i < nmax; i++) {
    int iatm = i / 3;
    int idim = i % 3;

    double *vv = v[iatm];

    kecurrent = mass[type[iatm]] * vv[idim] * vv[idim] * force->mvv2e;
    t_current = kecurrent / force->boltz;

    double *eta = nhc_eta[i];
    double *eta_dot = nhc_eta_dot[i];
    double *eta_dotdot = nhc_eta_dotdot[i];

    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for (int ichain = nhc_nchain - 1; ichain > 0; ichain--) {
      expfac = exp(-dt8 * eta_dot[ichain + 1]);
      eta_dot[ichain] *= expfac;
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    expfac = exp(-dt8 * eta_dot[1]);
    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    // Update particle velocities half-step

    double factor_eta = exp(-dthalf * eta_dot[0]);
    vv[idim] *= factor_eta;

    t_current *= (factor_eta * factor_eta);
    kecurrent = force->boltz * t_current;
    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for (int ichain = 0; ichain < nhc_nchain; ichain++) eta[ichain] += dthalf * eta_dot[ichain];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    for (int ichain = 1; ichain < nhc_nchain; ichain++) {
      expfac = exp(-dt8 * eta_dot[ichain + 1]);
      eta_dot[ichain] *= expfac;
      eta_dotdot[ichain] =
          (nhc_eta_mass[i][ichain - 1] * eta_dot[ichain - 1] * eta_dot[ichain - 1] - KT) /
          nhc_eta_mass[i][ichain];
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    t_sys += t_current;
  }

  t_sys /= nmax;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::b_step()
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

void FixPIMD::qc_step(){
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  tagint *tag = atom->tag;
  double oldlo, oldhi;
  if(!pstat_flag) {
    if(universe->iworld == 0)
    {
      for(int i=0; i<nlocal; i++)
      {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }
  else{
    if(universe->iworld == 0)
    {
      double expp[3], expq[3];
      if(pstyle == ISO) {vw[1] = vw[0]; vw[2] = vw[0];}
      for(int j=0; j<3; j++)
      {
        expq[j] = exp(dtv * vw[j]);
        expp[j] = exp(-dtv * vw[j]);
      }
      if(barostat == BZP)
      {
        for(int i=0; i<nlocal; i++)
        {
          for(int j=0; j<3; j++)
          {
            x[i][j] = expq[j] * x[i][j] + (expq[j] - expp[j]) / 2. / vw[j] * v[i][j];
            v[i][j] = expp[j] * v[i][j];
          } 
        }
        oldlo = domain->boxlo[0];
        oldhi = domain->boxhi[0];
      
        domain->boxlo[0] = (oldlo-fixedpoint[0])*expq[0] + fixedpoint[0];
        domain->boxhi[0] = (oldhi-fixedpoint[0])*expq[0] + fixedpoint[0];

        oldlo = domain->boxlo[1];
        oldhi = domain->boxhi[1];
        domain->boxlo[1] = (oldlo-fixedpoint[1])*expq[1] + fixedpoint[1];
        domain->boxhi[1] = (oldhi-fixedpoint[1])*expq[1] + fixedpoint[1];
  
        oldlo = domain->boxlo[2];
        oldhi = domain->boxhi[2];
        domain->boxlo[2] = (oldlo-fixedpoint[2])*expq[2] + fixedpoint[2];
        domain->boxhi[2] = (oldhi-fixedpoint[2])*expq[2] + fixedpoint[2];
      }
    }    
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&domain->boxlo[0], 3, MPI_DOUBLE, 0, universe->uworld);
    MPI_Bcast(&domain->boxhi[0], 3, MPI_DOUBLE, 0, universe->uworld);
    domain->set_global_box();
    domain->set_local_box();
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::a_step(){
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
      x[i][0] = Lan_c[universe->iworld] * x0 + 1./_omega_k[universe->iworld] * Lan_s[universe->iworld] * v0;
      x[i][1] = Lan_c[universe->iworld] * x1 + 1./_omega_k[universe->iworld] * Lan_s[universe->iworld] * v1;
      x[i][2] = Lan_c[universe->iworld] * x2 + 1./_omega_k[universe->iworld] * Lan_s[universe->iworld] * v2;
      v[i][0] = -1.*_omega_k[universe->iworld] * Lan_s[universe->iworld] * x0 + Lan_c[universe->iworld] * v0;
      v[i][1] = -1.*_omega_k[universe->iworld] * Lan_s[universe->iworld] * x1 + Lan_c[universe->iworld] * v1;
      v[i][2] = -1.*_omega_k[universe->iworld] * Lan_s[universe->iworld] * x2 + Lan_c[universe->iworld] * v2;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::remove_com_motion(){
  if(universe->iworld == 0)
  {
  // double **x = atom->x;
    double **v = atom->v;
  int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (dynamic)  masstotal = group->mass(igroup);
    double vcm[3];
    group->vcm(igroup,masstotal,vcm);    
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] -= vcm[0];
        v[i][1] -= vcm[1];
        v[i][2] -= vcm[2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::baro_init()
{
  if(pstyle == ISO) {W = 3 * (atom->natoms) * tau_p * tau_p * np * kBT;} // consistent with the definition in i-Pi
  // printf("tau_p = %.6e np = %.6e kBT = %.6e W = %.6e\n", tau_p, np, kBT, W);}
  else if(pstyle == ANISO) {W = atom->natoms * tau_p * tau_p * np * kBT;}
  Vcoeff = 1.0;
  std::string out = fmt::format("\nInitializing PIMD {:s} barostat...\n", Barostats[barostat]);
  out += fmt::format("The barostat mass is W = {:.16e}\n", W);
  utils::logmesg(lmp, out);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::press_v_step()
{
  int nlocal = atom->nlocal;
  double **f = atom->f;
  double **v = atom->v;
  int *type = atom->type; 
  volume = domain->xprd * domain->yprd * domain->zprd;

  if(pstyle == ISO)
  {
    if(barostat == BZP)
    {
      vw[0] += dtv * 3 * (volume * np * (p_cv - Pext) / force->nktv2p + Vcoeff / beta_np) / W;
      if(universe->iworld==0)
      {
        double dvw_proc = 0.0, dvw = 0.0;
        for(int i = 0; i < nlocal; i++)
        {
          for(int j = 0; j < 3; j++)
          {
            dvw_proc += dtv2 * f[i][j] * v[i][j] / W + dtv3 * f[i][j] * f[i][j] / mass[type[i]] / W;
          }
        }
        MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
        vw[0] += dvw;
      }
      MPI_Barrier(universe->uworld);
      MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
    }
    else if(barostat == MTTK)
    {
        mtk_term1 = 2. / atom->natoms * totke / 3;
        f_omega = (volume * np * (p_md - Pext) + mtk_term1) / W;
        vw[0] += 0.5 * dtv * f_omega;
    }
  }
  else if(pstyle == ANISO)
  {
    compute_stress_tensor();
    for(int ii=0; ii<3; ii++)
    {
      vw[ii] += dtv * (volume * np * (stress_tensor[ii] - Pext) / force->nktv2p + Vcoeff / beta_np) / W;
      if(universe->iworld==0)
      {
        double dvw_proc = 0.0, dvw = 0.0;
        for(int i = 0; i < nlocal; i++)
        {
          dvw_proc += dtv2 * f[i][ii] * v[i][ii] / W + dtv3 * f[i][ii] * f[i][ii] / mass[type[i]] / W;
        }
        MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
        vw[ii] += dvw;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::press_o_step()
{
  if(pstyle==ISO)
  {
    if(universe->me==0)
    {
      r1 = random->gaussian();
      vw[0] = c1 * vw[0] + c2 * sqrt(1. / W / beta_np) * r1;
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
  }
  else if(pstyle==ANISO)
  {
    if(universe->me==0)
    {
      r1 = random->gaussian();
      r2 = random->gaussian();
      r3 = random->gaussian();
      vw[0] = c1 * vw[0] + c2 * sqrt(1. / W / beta_np) * r1;
      vw[1] = c1 * vw[1] + c2 * sqrt(1. / W / beta_np) * r2;
      vw[2] = c1 * vw[2] + c2 * sqrt(1. / W / beta_np) * r3;
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw, 3, MPI_DOUBLE, 0, universe->uworld);    
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::Langevin_init()
{
  double KT = force->boltz * temp;
  double beta = 1.0 / KT;
  _omega_np = np / beta / hbar;
  double _omega_np_dt_half = _omega_np * update->dt * 0.5;

  _omega_k = new double[np];
  Lan_c = new double[np];
  Lan_s = new double[np];
  if(fmmode==physical){
    for (int i=0; i<np; i++)
    {
      _omega_k[i] = _omega_np * sqrt(lam[i]); 
      Lan_c[i] = cos(sqrt(lam[i])*_omega_np_dt_half);
      Lan_s[i] = sin(sqrt(lam[i])*_omega_np_dt_half);
    }
  }
  else if(fmmode==normal){
    for (int i=0; i<np; i++)
    {
      _omega_k[i] = _omega_np; 
      Lan_c[i] = cos(_omega_np_dt_half);
      Lan_s[i] = sin(_omega_np_dt_half);
    }
  }
  if(tau > 0) gamma = 1.0 / tau;
  else gamma = np / beta / hbar;
  
  if(integrator==obabo) c1 = exp(-gamma * 0.5 * update->dt); // tau is the damping time of the centroid mode.
  else if(integrator==baoab) c1 = exp(-gamma * update->dt); 
  else error->universe_all(FLERR, "Unknown integrator parameter for fix pimd. Only obabo and baoab integrators is supported!");

  c2 = sqrt(1.0 - c1 * c1); // note that c1 and c2 here only works for the centroid mode.

  if( thermostat == PILE_L )
  {
    std::string out = "\nInitializing PI Langevin equation thermostat...\n";
    out += "Bead ID    |    omega    |    tau    |    c1    |    c2\n"; 
    tau_k = new double[np];
    c1_k = new double[np];
    c2_k = new double[np];
    tau_k[0] = tau; c1_k[0] = c1; c2_k[0] = c2;
    for(int i=1; i<np; i++)
    {
      tau_k[i] = 0.5 / pilescale / _omega_k[i];
      if(integrator==obabo) c1_k[i] = exp(-0.5 * update->dt / tau_k[i]);
      else if(integrator==baoab) c1_k[i] = exp(-1.0 * update->dt / tau_k[i]);
      else error->universe_all(FLERR, "Unknown integrator parameter for fix pimd. Only obabo and baoab integrators is supported!");
      c2_k[i] = sqrt(1.0 - c1_k[i] * c1_k[i]);
    }
    for(int i=0; i<np; i++)
    {
      out += fmt::format("    {:d}     {:.8e} {:.8e} {:.8e} {:.8e}\n", i, _omega_k[i], tau_k[i], c1_k[i], c2_k[i]);
    }
    if(thermostat == PILE_L) out += "PILE_L thermostat successfully initialized!\n";
    out += "\n";
    utils::logmesg(lmp, out);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::o_step()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double beta_np = 1.0 / force->boltz / Lan_temp / np * force->mvv2e;
  if(thermostat == PILE_L)
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

/* ----------------------------------------------------------------------
   Normal Mode PIMD
------------------------------------------------------------------------- */

void FixPIMD::nmpimd_init()
{
  if(ensemble == nvt && thermostat == NHC)
  {
    memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
    memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");
    memory->create(M_f2fp, np, np, "fix_feynman:M_f2fp");
    memory->create(M_fp2f, np, np, "fix_feynman:M_fp2f");
  
    lam = (double *) memory->smalloc(sizeof(double) * np, "FixPIMD::lam");
  
    // Set up  eigenvalues
  
    lam[0] = 0.0;
    if (np % 2 == 0) lam[np - 1] = 4.0 * np;
  
    for (int i = 2; i <= np / 2; i++) {
      lam[2 * i - 3] = lam[2 * i - 2] = 2.0 * np * (1.0 - 1.0 * cos(2.0 * MY_PI * (i - 1) / np));
    }
  
    // Set up eigenvectors for non-degenerated modes
  
    for (int i = 0; i < np; i++) {
      M_x2xp[0][i] = 1.0 / np;
      if (np % 2 == 0) M_x2xp[np - 1][i] = 1.0 / np * pow(-1.0, i);
    }
  
    // Set up eigenvectors for degenerated modes
  
    for (int i = 0; i < (np - 1) / 2; i++)
      for (int j = 0; j < np; j++) {
        M_x2xp[2 * i + 1][j] = sqrt(2.0) * cos(2.0 * MY_PI * (i + 1) * j / np) / np;
        M_x2xp[2 * i + 2][j] = -sqrt(2.0) * sin(2.0 * MY_PI * (i + 1) * j / np) / np;
      }
  
    // Set up Ut
  
    for (int i = 0; i < np; i++)
      for (int j = 0; j < np; j++) {
        M_xp2x[i][j] = M_x2xp[j][i] * np;
        M_f2fp[i][j] = M_x2xp[i][j] * np;
        M_fp2f[i][j] = M_xp2x[i][j];
      }
  
    // Set up masses
  
    int iworld = universe->iworld;
  
    for (int i = 1; i <= atom->ntypes; i++) {
      mass[i] = atom->mass[i];
  
      if (iworld) {
        mass[i] *= lam[iworld];
        mass[i] *= fmass;
      }
    }
  }
  else
  {
    memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
    memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");

    lam = (double*) memory->smalloc(sizeof(double)*np, "FixPIMD::lam");

    // Set up  eigenvalues
    for(int i=0; i<np; i++){
      double sin_tmp = sin(i * MY_PI / np);
      lam[i] = 4 * sin_tmp * sin_tmp;
    }

    // Set up eigenvectors for degenerated modes
    for(int j=0; j<np; j++){
      for(int i=1; i<int(np/2) + 1; i++) 
      {
        M_x2xp[i][j] =   sqrt(2.0) * cos ( 2.0 * MY_PI * double(i) * double(j) / double(np)) / sqrt(np);
      }
      for(int i=int(np/2)+1; i<np; i++)
      {
        M_x2xp[i][j] =   sqrt(2.0) * sin ( 2.0 * MY_PI * double(i) * double(j) / double(np)) / sqrt(np);
      }
    }

    // Set up eigenvectors for non-degenerated modes
    for(int i=0; i<np; i++)
    {
      M_x2xp[0][i] = 1.0 / sqrt(np);
      if(np%2==0) M_x2xp[np/2][i] = 1.0 / sqrt(np) * pow(-1.0, i);
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
        mass[i] *= fmass;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nmpimd_fill(double **ptr)
{
  comm_ptr = ptr;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nmpimd_transform(double **src, double **des, double *vector)
{
  int n = atom->nlocal;
  int m = 0;

  for (int i = 0; i < n; i++)
    for (int d = 0; d < 3; d++) {
      des[i][d] = 0.0;
      for (int j = 0; j < np; j++) { des[i][d] += (src[j][m] * vector[j]); }
      m++;
    }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::spring_force()
{
  spring_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *_mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double *xlast = buf_beads[x_last];
  double *xnext = buf_beads[x_next];

  for (int i = 0; i < nlocal; i++) {
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

    double dx = delx1 + delx2;
    double dy = dely1 + dely2;
    double dz = delz1 + delz2;

    f[i][0] -= (dx) *ff;
    f[i][1] -= (dy) *ff;
    f[i][2] -= (dz) *ff;

    spring_energy += (dx * dx + dy * dy + dz * dz);
  }
}

/* ----------------------------------------------------------------------
   Comm operations
------------------------------------------------------------------------- */

void FixPIMD::comm_init()
{
  if (size_plan) {
    delete[] plan_send;
    delete[] plan_recv;
  }

  if (method == PIMD) {
    size_plan = 2;
    plan_send = new int[2];
    plan_recv = new int[2];
    mode_index = new int[2];

    int rank_last = universe->me - comm->nprocs;
    int rank_next = universe->me + comm->nprocs;
    if (rank_last < 0) rank_last += universe->nprocs;
    if (rank_next >= universe->nprocs) rank_next -= universe->nprocs;

    plan_send[0] = rank_next;
    plan_send[1] = rank_last;
    plan_recv[0] = rank_last;
    plan_recv[1] = rank_next;

    mode_index[0] = 0;
    mode_index[1] = 1;
    x_last = 1;
    x_next = 0;
  } else {
    size_plan = np - 1;
    plan_send = new int[size_plan];
    plan_recv = new int[size_plan];
    mode_index = new int[size_plan];

    for (int i = 0; i < size_plan; i++) {
      plan_send[i] = universe->me + comm->nprocs * (i + 1);
      if (plan_send[i] >= universe->nprocs) plan_send[i] -= universe->nprocs;

      plan_recv[i] = universe->me - comm->nprocs * (i + 1);
      if (plan_recv[i] < 0) plan_recv[i] += universe->nprocs;

      mode_index[i] = (universe->iworld + i + 1) % (universe->nworlds);
    }

    x_next = (universe->iworld + 1 + universe->nworlds) % (universe->nworlds);
    x_last = (universe->iworld - 1 + universe->nworlds) % (universe->nworlds);
  }

  if (buf_beads) {
    for (int i = 0; i < np; i++) delete[] buf_beads[i];
    delete[] buf_beads;
  }

  buf_beads = new double *[np];
  for (int i = 0; i < np; i++) buf_beads[i] = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::comm_exec(double **ptr)
{
  int nlocal = atom->nlocal;

  if (nlocal > max_nlocal) {
    max_nlocal = nlocal + 200;
    int size = sizeof(double) * max_nlocal * 3;
    buf_recv = (double *) memory->srealloc(buf_recv, size, "FixPIMD:x_recv");

    for (int i = 0; i < np; i++)
      buf_beads[i] = (double *) memory->srealloc(buf_beads[i], size, "FixPIMD:x_beads[i]");
  }

  // copy local positions

  memcpy(buf_beads[universe->iworld], &(ptr[0][0]), sizeof(double) * nlocal * 3);

  // go over comm plans

  for (int iplan = 0; iplan < size_plan; iplan++) {
    // sendrecv nlocal

    int nsend;

    MPI_Sendrecv(&(nlocal), 1, MPI_INT, plan_send[iplan], 0, &(nsend), 1, MPI_INT, plan_recv[iplan],
                 0, universe->uworld, MPI_STATUS_IGNORE);

    // allocate arrays

    if (nsend > max_nsend) {
      max_nsend = nsend + 200;
      tag_send =
          (tagint *) memory->srealloc(tag_send, sizeof(tagint) * max_nsend, "FixPIMD:tag_send");
      buf_send =
          (double *) memory->srealloc(buf_send, sizeof(double) * max_nsend * 3, "FixPIMD:x_send");
    }

    // send tags

    MPI_Sendrecv(atom->tag, nlocal, MPI_LMP_TAGINT, plan_send[iplan], 0, tag_send, nsend,
                 MPI_LMP_TAGINT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // wrap positions

    double *wrap_ptr = buf_send;
    int ncpy = sizeof(double) * 3;

    for (int i = 0; i < nsend; i++) {
      int index = atom->map(tag_send[i]);

      if (index < 0) {
        auto mesg = fmt::format("Atom {} is missing at world [{}] rank [{}] "
                                "required by rank [{}] ({}, {}, {}).\n",
                                tag_send[i], universe->iworld, comm->me, plan_recv[iplan],
                                atom->tag[0], atom->tag[1], atom->tag[2]);
        error->universe_one(FLERR, mesg);
      }

      memcpy(wrap_ptr, ptr[index], ncpy);
      wrap_ptr += 3;
    }

    // sendrecv x

    MPI_Sendrecv(buf_send, nsend * 3, MPI_DOUBLE, plan_recv[iplan], 0, buf_recv, nlocal * 3,
                 MPI_DOUBLE, plan_send[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // copy x

    memcpy(buf_beads[mode_index[iplan]], buf_recv, sizeof(double) * nlocal * 3);
  }
}

/* ---------------------------------------------------------------------- */

int FixPIMD::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = comm_ptr[j][0];
    buf[m++] = comm_ptr[j][1];
    buf[m++] = comm_ptr[j][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    comm_ptr[i][0] = buf[m++];
    comm_ptr[i][1] = buf[m++];
    comm_ptr[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   Memory operations
------------------------------------------------------------------------- */

double FixPIMD::memory_usage()
{
  return (double) atom->nmax * size_peratom_cols * sizeof(double);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::grow_arrays(int nmax)
{
  if (nmax == 0) return;
  int count = nmax * 3;

  memory->grow(array_atom, nmax, size_peratom_cols, "FixPIMD::array_atom");
  memory->grow(nhc_eta, count, nhc_nchain, "FixPIMD::nh_eta");
  memory->grow(nhc_eta_dot, count, nhc_nchain + 1, "FixPIMD::nh_eta_dot");
  memory->grow(nhc_eta_dotdot, count, nhc_nchain, "FixPIMD::nh_eta_dotdot");
  memory->grow(nhc_eta_mass, count, nhc_nchain, "FixPIMD::nh_eta_mass");
}

/* ---------------------------------------------------------------------- */

void FixPIMD::copy_arrays(int i, int j, int /*delflag*/)
{
  int i_pos = i * 3;
  int j_pos = j * 3;

  memcpy(nhc_eta[j_pos], nhc_eta[i_pos], nhc_size_one_1);
  memcpy(nhc_eta_dot[j_pos], nhc_eta_dot[i_pos], nhc_size_one_2);
  memcpy(nhc_eta_dotdot[j_pos], nhc_eta_dotdot[i_pos], nhc_size_one_1);
  memcpy(nhc_eta_mass[j_pos], nhc_eta_mass[i_pos], nhc_size_one_1);
}

/* ---------------------------------------------------------------------- */

int FixPIMD::pack_exchange(int i, double *buf)
{
  int offset = 0;
  int pos = i * 3;

  memcpy(buf + offset, nhc_eta[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_dot[pos], nhc_size_one_2);
  offset += nhc_offset_one_2;
  memcpy(buf + offset, nhc_eta_dotdot[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_mass[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::unpack_exchange(int nlocal, double *buf)
{
  int offset = 0;
  int pos = nlocal * 3;

  memcpy(nhc_eta[pos], buf + offset, nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos], buf + offset, nhc_size_one_2);
  offset += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], buf + offset, nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos], buf + offset, nhc_size_one_1);
  offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::pack_restart(int i, double *buf)
{
  int offset = 0;
  int pos = i * 3;
  // pack buf[0] this way because other fixes unpack it
  buf[offset++] = size_peratom_cols + 1;

  memcpy(buf + offset, nhc_eta[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_dot[pos], nhc_size_one_2);
  offset += nhc_offset_one_2;
  memcpy(buf + offset, nhc_eta_dotdot[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_mass[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;

  return size_peratom_cols + 1;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  int pos = nlocal * 3;

  memcpy(nhc_eta[pos], extra[nlocal] + m, nhc_size_one_1);
  m += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos], extra[nlocal] + m, nhc_size_one_2);
  m += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], extra[nlocal] + m, nhc_size_one_1);
  m += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos], extra[nlocal] + m, nhc_size_one_1);
  m += nhc_offset_one_1;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::maxsize_restart()
{
  return size_peratom_cols + 1;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::size_restart(int /*nlocal*/)
{
  return size_peratom_cols + 1;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_vir_()
{
  int nlocal = atom->nlocal;
  xf = vir_ = xcf = centroid_vir = 0.0;
  for(int i=0; i<nlocal; i++)
  {
    for(int j=0; j<3; j++)
    {
      xf += x_unwrap[3*i+j] * atom->f[i][j];
      xcf += (x_unwrap[3*i+j] - xc[3*i+j]) * atom->f[i][j];
    }
  }
  MPI_Allreduce(&xf, &vir_, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&xcf, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  if(pstyle == ANISO){
    for(int i=0; i<6; i++) c_vir_tensor[i] = 0.0;
    for(int i=0; i<nlocal; i++){
      c_vir_tensor[0] += (x_unwrap[3*i+0] - xc[3*i+0]) * atom->f[i][0];
      c_vir_tensor[1] += (x_unwrap[3*i+1] - xc[3*i+1]) * atom->f[i][1];
      c_vir_tensor[2] += (x_unwrap[3*i+2] - xc[3*i+2]) * atom->f[i][2];
      c_vir_tensor[3] += (x_unwrap[3*i+0] - xc[3*i+0]) * atom->f[i][1];
      c_vir_tensor[4] += (x_unwrap[3*i+0] - xc[3*i+0]) * atom->f[i][2];
      c_vir_tensor[5] += (x_unwrap[3*i+1] - xc[3*i+1]) * atom->f[i][2];
    }
    MPI_Allreduce(MPI_IN_PLACE, &c_vir_tensor, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_vir()
{
  volume = domain->xprd * domain->yprd * domain->zprd;
  c_press->compute_vector();
  virial[0] = c_press->vector[0]*volume;
  virial[1] = c_press->vector[1]*volume;
  virial[2] = c_press->vector[2]*volume;
  virial[3] = c_press->vector[3]*volume;
  virial[4] = c_press->vector[4]*volume;
  virial[5] = c_press->vector[5]*volume;
  for(int i=0; i<6; i++) virial[i] /= universe->procs_per_world[universe->iworld];
  vir=(virial[0]+virial[1]+virial[2]);
  MPI_Allreduce(&vir,&vir,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  MPI_Allreduce(MPI_IN_PLACE, &virial[0], 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_stress_tensor()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  if(universe->iworld == 0){
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    for(int i=0; i<6; i++) ke_tensor[i] = 0.0;
    for(int i=0; i<nlocal; i++){
      ke_tensor[0] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][0] * force->mvv2e;
      ke_tensor[1] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][1] * force->mvv2e;
      ke_tensor[2] += 0.5 * mass[type[i]] * atom->v[i][2] * atom->v[i][2] * force->mvv2e;
      ke_tensor[3] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][1] * force->mvv2e;
      ke_tensor[4] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][2] * force->mvv2e;
      ke_tensor[5] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][2] * force->mvv2e;
    }
    // MPI_Allreduce(&ke_tensor, &ke_tensor, 6, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, &ke_tensor, 6, MPI_DOUBLE, MPI_SUM, world);
    for(int i=0; i<6; i++) 
    {
      stress_tensor[i] = inv_volume * ((2*ke_tensor[i] - c_vir_tensor[i]) * force->nktv2p + virial[i]) / np;
      // printf("i = %d, c_vir = %.6f, virial = %.6f stress = %.6f\n", i, c_vir_tensor[i], virial[i], stress_tensor[i]);
    }
  }
  MPI_Bcast(&stress_tensor, 6, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_totke()
{
  kine = 0.0;
  totke = ke_bead = 0.0;
  int nlocal = atom->nlocal;
  int *type = atom->type;
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

void FixPIMD::compute_spring_energy()
{
  spring_energy = 0.0;
  total_spring_energy = se_bead = 0.0;

  double **x = atom->x;
  double* _mass = atom->mass;
  int* type = atom->type;
  int nlocal = atom->nlocal;

  for(int i=0; i<nlocal; i++)
  {
    spring_energy += 0.5 * _mass[type[i]] * fbond * lam[universe->iworld] * (x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2]); 
  }
  MPI_Allreduce(&spring_energy, &se_bead, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&spring_energy, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  total_spring_energy /= np;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_pote()
{
  pe_bead = 0.0;
  pot_energy_partition = 0.0;
  pote = 0.0;
  c_pe->compute_scalar();
  pe_bead = c_pe->scalar;
  pot_energy_partition = pe_bead / universe->procs_per_world[universe->iworld];
  MPI_Allreduce(&pot_energy_partition, &pote, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  pote /= np;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_tote()
{
  tote = totke + pote + total_spring_energy;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_t_prim()
{
  t_prim = 1.5 * atom->natoms * np * force->boltz * temp - total_spring_energy;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_t_vir()
{
  t_vir = -0.5 / np * vir_;
  t_cv = 1.5 * atom->natoms * force->boltz * temp - 0.5 / np * centroid_vir;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_p_prim()
{
  p_prim = atom->natoms * np * force->boltz * temp * inv_volume - 1.0 / 1.5 * inv_volume * total_spring_energy;
  p_prim *= force->nktv2p;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_p_cv()
{
  inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_md = 2. / 3 * inv_volume * ((totke - total_spring_energy) * force->nktv2p + 0.5 * vir / np) ;
  if(universe->iworld == 0)
  {
    p_cv = 1. / 3.  * inv_volume  * ((2. * ke_bead  - 1. * centroid_vir) * force->nktv2p + 1. * vir) / np; 
  }
  MPI_Bcast(&p_cv, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::compute_totenthalpy()
{
  volume = domain->xprd * domain->yprd * domain->zprd;
  if(barostat == BZP)  
  {
    if(pstyle == ISO)
    {
      totenthalpy = tote + 0.5*W*vw[0]*vw[0]/np + Pext * volume / force->nktv2p - Vcoeff * kBT * log(volume);
    }
    else if(pstyle == ANISO)
    {
      totenthalpy = tote + 0.5*W*vw[0]*vw[0]/np + 0.5*W*vw[1]*vw[1]/np + 0.5*W*vw[2]*vw[2]/np + Pext * volume / force->nktv2p - Vcoeff * kBT * log(volume);
    }
  }
  else if(barostat == MTTK)  totenthalpy = tote + 1.5*W*vw[0]*vw[0]/np + Pext * (volume - vol0);
}

/* ---------------------------------------------------------------------- */

double FixPIMD::compute_vector(int n)
{
  if(n==0) { return ke_bead; }
  if(n==1) { return se_bead; }
  if(n==2) { return pe_bead; }
  if(n==3) { return tote; }
  if(n==4) { return t_prim; }
  if(n==5) { return t_vir; }
  if(n==6) { return t_cv; }
  if(n==7) { return p_prim; }
  if(n==8) { return p_md; }
  if(n==9) { return p_cv; }
  if(n==10) {return totenthalpy;}
  if(pstyle == ISO){
    if(n==11) { return vw[0]; }
    if(n==12) { 
      if(barostat == BZP) {  return 0.5*W*vw[0]*vw[0]; }
      else if(barostat == MTTK) {  return 1.5*W*vw[0]*vw[0]; }
    }
    if(n==13) { volume = domain->xprd * domain->yprd * domain->zprd; return np * Pext * volume / force->nktv2p; }
    if(n==14) { volume = domain->xprd * domain->yprd * domain->zprd; return - Vcoeff * np * kBT * log(volume); }
  }
  else if(pstyle==ANISO){
    if(n>10 && n<=13) return vw[n-11];
    if(n==14) return 0.5*W*vw[0]*vw[0]+0.5*W*vw[1]*vw[1]+0.5*W*vw[2]*vw[2];
    if(n>14 && n<21) return stress_tensor[n-15];
    if(n==21) { volume = domain->xprd * domain->yprd * domain->zprd; return np * Pext * volume / force->nktv2p; }
    if(n==22) { volume = domain->xprd * domain->yprd * domain->zprd; return - Vcoeff * np * kBT * log(volume); }
  }
  return 0.0;
}

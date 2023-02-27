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
   Package      FixPIMDLangevin
   Purpose      Quantum Path Integral Algorithm for Quantum Chemistry

   Yifan Li @ Princeton University (yifanl0716@gmail.com)
   Added components:
   - Multi-processor parallelism for each bead
   - White-noise Langevin thermostat
   - Bussi-Zykova-Parrinello barostat (isotropic and anisotropic)
   - Several quantum estimators
   Futher plans:
   - Triclinic barostat
------------------------------------------------------------------------- */

#include "fix_pimd_langevin.h"

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
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { PIMD, NMPIMD, CMD };
enum { PHYSICAL, NORMAL };
enum { BAOAB, OBABO };
enum { ISO, ANISO, TRICLINIC };
enum { PILE_L };
enum { MTTK, BZP };
// char* Barostats[] = {"MTTK", "BZP"};

static std::map<int, std::string> Barostats{{MTTK, "MTTK"}, {BZP, "BZP"}};
enum { NVE, NVT, NPH, NPT };
enum { SINGLE_PROC, MULTI_PROC };

/* ---------------------------------------------------------------------- */

FixPIMDLangevin::FixPIMDLangevin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), random(nullptr), c_pe(nullptr), c_press(nullptr)
{
  time_integrate = 1;
  tagsend = tagrecv = nullptr;
  bufsend = bufrecv = nullptr;
  bufsendall = bufrecvall = nullptr;
  bufsorted = bufsortedall = nullptr;
  outsorted = buftransall = nullptr;

  ntotal = 0;
  maxlocal = maxunwrap = maxxc = 0;
  bufbeads = nullptr;
  x_unwrap = xc = nullptr;
  xcall = nullptr;
  counts = nullptr;

  sizeplan = 0;
  plansend = planrecv = nullptr;

  M_x2xp = M_xp2x = M_f2fp = M_fp2f = nullptr;
  lam = nullptr;
  modeindex = nullptr;

  mass = nullptr;

  method = PIMD;
  ensemble = NVT;
  integrator = OBABO;
  thermostat = PILE_L;
  barostat = BZP;
  fmass = 1.0;
  sp = 1.0;
  temp = 298.15;
  Lan_temp = 298.15;
  tau = 1.0;
  tau_p = 1.0;
  Pext = 1.0;
  pilescale = 1.0;
  tstat_flag = 1;
  pstat_flag = 0;
  mapflag = 1;
  removecomflag = 1;
  fmmode = PHYSICAL;
  pstyle = ISO;
  totenthalpy = 0.0;

  int seed = -1;

  for (int i = 0; i < 6; i++) p_flag[i] = 0;

  for (int i = 3; i < narg - 1; i += 2) {
    if (strcmp(arg[i], "method") == 0) {
      if (strcmp(arg[i + 1], "pimd") == 0)
        method = PIMD;
      else if (strcmp(arg[i + 1], "nmpimd") == 0)
        method = NMPIMD;
      else if (strcmp(arg[i + 1], "cmd") == 0)
        method = CMD;
      else
        error->universe_all(FLERR, "Unknown method parameter for fix pimd/langevin");
    } else if (strcmp(arg[i], "integrator") == 0) {
      if (strcmp(arg[i + 1], "obabo") == 0)
        integrator = OBABO;
      else if (strcmp(arg[i + 1], "baoab") == 0)
        integrator = BAOAB;
      else
        error->universe_all(
            FLERR,
            "Unknown integrator parameter for fix pimd/langevin. Only obabo and baoab "
            "integrators are supported!");
    } else if (strcmp(arg[i], "ensemble") == 0) {
      if (strcmp(arg[i + 1], "nve") == 0) {
        ensemble = NVE;
        tstat_flag = 0;
        pstat_flag = 0;
      } else if (strcmp(arg[i + 1], "nvt") == 0) {
        ensemble = NVT;
        tstat_flag = 1;
        pstat_flag = 0;
      } else if (strcmp(arg[i + 1], "nph") == 0) {
        ensemble = NPH;
        tstat_flag = 0;
        pstat_flag = 1;
      } else if (strcmp(arg[i + 1], "npt") == 0) {
        ensemble = NPT;
        tstat_flag = 1;
        pstat_flag = 1;
      } else
        error->universe_all(FLERR,
                            "Unknown ensemble parameter for fix pimd/langevin. Only nve and nvt "
                            "ensembles are supported!");
    } else if (strcmp(arg[i], "fmass") == 0) {
      fmass = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (fmass < 0.0 || fmass > 1.0)
        error->universe_all(FLERR, "Invalid fmass value for fix pimd/langevin");
    } else if (strcmp(arg[i], "fmmode") == 0) {
      if (strcmp(arg[i + 1], "physical") == 0)
        fmmode = PHYSICAL;
      else if (strcmp(arg[i + 1], "normal") == 0)
        fmmode = NORMAL;
      else
        error->universe_all(
            FLERR,
            "Unknown fictitious mass mode for fix pimd/langevin. Only physical mass and "
            "normal mode mass are supported!");
    } else if (strcmp(arg[i], "scale") == 0) {
      pilescale = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (pilescale < 0.0)
        error->universe_all(FLERR, "Invalid pile scale value for fix pimd/langevin");
    } else if (strcmp(arg[i], "temp") == 0) {
      temp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (temp < 0.0) error->universe_all(FLERR, "Invalid temp value for fix pimd/langevin");
    } else if (strcmp(arg[i], "lj") == 0) {
      lj_epsilon = utils::numeric(FLERR, arg[i + 1], false, lmp);
      lj_sigma = utils::numeric(FLERR, arg[i + 2], false, lmp);
      lj_mass = utils::numeric(FLERR, arg[i + 3], false, lmp);
      other_planck = utils::numeric(FLERR, arg[i + 4], false, lmp);
      i += 3;
    } else if (strcmp(arg[i], "thermostat") == 0) {
      if (strcmp(arg[i + 1], "PILE_L") == 0) {
        thermostat = PILE_L;
        seed = utils::inumeric(FLERR, arg[i + 2], false, lmp);
        i++;
      }
    } else if (strcmp(arg[i], "tau") == 0) {
      tau = utils::numeric(FLERR, arg[i + 1], false, lmp);
    } else if (strcmp(arg[i], "press") == 0) {
      Pext = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (Pext < 0.0) error->universe_all(FLERR, "Invalid press value for fix pimd/langevin");
    } else if (strcmp(arg[i], "barostat") == 0) {
      if (strcmp(arg[i + 1], "MTTK") == 0) {
        barostat = MTTK;
      } else if (strcmp(arg[i + 1], "BZP") == 0) {
        barostat = BZP;
      } else
        error->universe_all(FLERR, "Unknown barostat parameter for fix pimd/langevin");
    } else if (strcmp(arg[i], "iso") == 0) {
      pstyle = ISO;
      i--;
    } else if (strcmp(arg[i], "aniso") == 0) {
      pstyle = ANISO;
      i--;
    } else if (strcmp(arg[i], "taup") == 0) {
      tau_p = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (tau_p <= 0.0) error->universe_all(FLERR, "Invalid tau_p value for fix pimd/langevin");
    } else if (strcmp(arg[i], "fixcom") == 0) {
      if (strcmp(arg[i + 1], "yes") == 0)
        removecomflag = 1;
      else if (strcmp(arg[i + 1], "no") == 0)
        removecomflag = 0;
    } else if (strcmp(arg[i], "map") == 0) {
      if (strcmp(arg[i + 1], "yes") == 0)
        mapflag = 1;
      else if (strcmp(arg[i + 1], "no") == 0)
        mapflag = 0;
    } else {
      error->universe_all(FLERR, fmt::format("Unknown keyword {} for fix {}", arg[i], style));
    }
  }

  /* Initiation */

  global_freq = 1;
  vector_flag = 1;
  size_vector = 11;
  extvector = 1;

  // some initilizations

  id_pe = utils::strdup(std::string(id) + "_pimd_pe");
  modify->add_compute(std::string(id_pe) + " all pe");

  id_press = utils::strdup(std::string(id) + "_pimd_press");
  modify->add_compute(std::string(id_press) + " all pressure thermo_temp virial");

  vol0 = domain->xprd * domain->yprd * domain->zprd;

  fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
  fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
  fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);

  // initialize Marsaglia RNG with processor-unique seed

  if (integrator == BAOAB || integrator == OBABO) {
    Lan_temp = temp;
    random = new RanMars(lmp, seed + universe->me);
  }

  me = comm->me;
  nprocs = comm->nprocs;
  if (nprocs == 1)
    cmode = SINGLE_PROC;
  else
    cmode = MULTI_PROC;

  nprocs_universe = universe->nprocs;
  nreplica = universe->nworlds;
  ireplica = universe->iworld;

  int *iroots = new int[nreplica];
  MPI_Group uworldgroup, rootgroup;

  for (int i = 0; i < nreplica; i++) iroots[i] = universe->root_proc[i];
  MPI_Comm_group(universe->uworld, &uworldgroup);
  MPI_Group_incl(uworldgroup, nreplica, iroots, &rootgroup);
  MPI_Comm_create(universe->uworld, rootgroup, &rootworld);
  if (rootgroup != MPI_GROUP_NULL) MPI_Group_free(&rootgroup);
  if (uworldgroup != MPI_GROUP_NULL) MPI_Group_free(&uworldgroup);
  delete[] iroots;

  ntotal = atom->natoms;
  if (atom->nmax > maxlocal) reallocate();
  if (atom->nmax > maxunwrap) reallocate_x_unwrap();
  if (atom->nmax > maxxc) reallocate_xc();
  memory->create(xcall, ntotal * 3, "FixPIMDLangevin:xcall");

  if (cmode == SINGLE_PROC) {
    memory->create(bufsorted, ntotal, 3, "FixPIMDLangevin:bufsorted");
    memory->create(outsorted, ntotal, 3, "FixPIMDLangevin:outsorted");
    memory->create(bufsortedall, nreplica * ntotal, 3, "FixPIMDLangevin:bufsortedall");
    memory->create(buftransall, nreplica * ntotal, 3, "FixPIMDLangevin:buftransall");
    memory->create(counts, nreplica, "FixPIMDLangevin:counts");
    memory->create(displacements, nreplica, "FixPIMDLangevin:displacements");
  }

  if ((cmode == MULTI_PROC) && (counts == nullptr)) {
    memory->create(bufsendall, ntotal, 3, "FixPIMDLangevin:bufsendall");
    memory->create(bufrecvall, ntotal, 3, "FixPIMDLangevin:bufrecvall");
    memory->create(tagsendall, ntotal, "FixPIMDLangevin:tagsendall");
    memory->create(tagrecvall, ntotal, "FixPIMDLangevin:tagrecvall");
    memory->create(counts, nprocs, "FixPIMDLangevin:counts");
    memory->create(displacements, nprocs, "FixPIMDLangevin:displacements");
  }
}

/* ---------------------------------------------------------------------- */

FixPIMDLangevin::~FixPIMDLangevin()
{
  modify->delete_compute(id_pe);
  modify->delete_compute(id_press);
  delete[] id_pe;
  delete[] id_press;
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixPIMDLangevin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  // mask |= POST_NEIGHBOR;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::init()
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "fix pimd/langevin requires an atom map, see atom_modify");

  if (universe->me == 0 && universe->uscreen)
    fprintf(universe->uscreen, "fix pimd/langevin initializing Path-Integral ...\n");

  // prepare the constants

  masstotal = group->mass(igroup);
  np = universe->nworlds;
  inverse_np = 1.0 / np;

  double planck;
  if (strcmp(update->unit_style, "lj") == 0) {
    double planck_star = sqrt(lj_epsilon) * sqrt(atom->mass[0]) * lj_sigma;
    planck = other_planck / planck_star;
  } else {
    planck = force->hplanck;
  }

  hbar = planck / (2.0 * MY_PI);
  kBT = force->boltz * temp;
  double beta = 1.0 / (force->boltz * temp);
  double _fbond = 1.0 * np * np / (beta * beta * hbar * hbar);

  omega_np = np / (hbar * beta) * sqrt(force->mvv2e);
  beta_np = 1.0 / force->boltz / temp * inverse_np;
  fbond = _fbond * force->mvv2e;

  if ((universe->me == 0) && (universe->uscreen))
    fprintf(universe->uscreen, "fix pimd/langevin -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", fbond);

  if (integrator == OBABO) {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = THIRD * dtv2 * dtv * force->ftm2v;
  } else if (integrator == BAOAB) {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = THIRD * dtv2 * dtv * force->ftm2v;
  } else {
    error->universe_all(FLERR, "Unknown integrator parameter for fix pimd/langevin");
  }

  comm_init();

  mass = new double[atom->ntypes + 1];

  nmpimd_init();

  Langevin_init();
  if (pstat_flag) baro_init();

  c_pe = modify->get_compute_by_id(id_pe);
  c_press = modify->get_compute_by_id(id_press);

  t_prim = t_vir = t_cv = p_prim = p_vir = p_cv = p_md = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::setup(int vflag)
{
  if (universe->me == 0) printf("Setting up Path-Integral ...\n");
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap(x[i], image[i]);
  }

  if (method == NMPIMD) {
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
  }
  collect_xc();
  compute_spring_energy();
  if (method == NMPIMD) {
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_xp2x[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_xp2x[universe->iworld]);
  }
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap_inv(x[i], image[i]);
  }

  if (method == NMPIMD) {
    inter_replica_comm(v);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, v, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, v, M_x2xp[universe->iworld]);
  }
  if (universe->me == 0 && screen) fprintf(screen, "Setting up Path-Integral ...\n");
  if (universe->me == 0) printf("Setting up Path-Integral ...\n");
  post_force(vflag);
  compute_totke();
  end_of_step();
  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::initial_integrate(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap(x[i], image[i]);
  }
  if (integrator == OBABO) {
    if (tstat_flag) {
      o_step();
      if (pstat_flag) press_o_step();
    }
    if (pstat_flag) {
      compute_totke();
      compute_p_cv();
      press_v_step();
    }
    b_step();
    if (method == NMPIMD) {
      inter_replica_comm(x);
      if (cmode == SINGLE_PROC)
        nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
      else if (cmode == MULTI_PROC)
        nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
    }
    qc_step();
    a_step();
    qc_step();
    a_step();
  } else if (integrator == BAOAB) {
    if (pstat_flag) {
      compute_totke();
      compute_p_cv();
      press_v_step();
    }
    b_step();
    if (method == NMPIMD) {
      inter_replica_comm(x);
      if (cmode == SINGLE_PROC)
        nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
      else if (cmode == MULTI_PROC)
        nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
    }
    qc_step();
    a_step();
    if (tstat_flag) {
      o_step();
      if (pstat_flag) press_o_step();
    }
    qc_step();
    a_step();
  } else {
    error->universe_all(FLERR, "Unknown integrator parameter for fix pimd/langevin");
  }
  collect_xc();
  compute_spring_energy();

  if (method == NMPIMD) {
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_xp2x[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_xp2x[universe->iworld]);
  }

  if (mapflag) {
    for (int i = 0; i < nlocal; i++) { domain->unmap_inv(x[i], image[i]); }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::final_integrate()
{
  if (pstat_flag) {
    compute_totke();
    compute_p_cv();
    press_v_step();
  }
  b_step();
  if (integrator == OBABO) {
    if (tstat_flag) {
      o_step();
      if (pstat_flag) press_o_step();
    }
  } else if (integrator == BAOAB) {

  } else {
    error->universe_all(FLERR, "Unknown integrator parameter for fix pimd/langevin");
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::post_force(int /*flag*/)
{
  if (atom->nmax > maxunwrap) reallocate_x_unwrap();
  if (atom->nmax > maxxc) reallocate_xc();
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  imageint *image = atom->image;
  tagint *tag = atom->tag;
  for (int i = 0; i < nlocal; i++) {
    x_unwrap[i][0] = x[i][0];
    x_unwrap[i][1] = x[i][1];
    x_unwrap[i][2] = x[i][2];
  }
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) { domain->unmap(x_unwrap[i], image[i]); }
  }
  for (int i = 0; i < nlocal; i++) {
    xc[i][0] = xcall[3 * (tag[i] - 1) + 0];
    xc[i][1] = xcall[3 * (tag[i] - 1) + 1];
    xc[i][2] = xcall[3 * (tag[i] - 1) + 2];
  }

  compute_vir();
  compute_vir_();
  compute_pote();
  if (method == NMPIMD) {
    inter_replica_comm(f);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, f, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, f, M_x2xp[universe->iworld]);
  }
  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::end_of_step()
{
  compute_totke();
  compute_p_cv();
  compute_tote();
  if (pstat_flag) compute_totenthalpy();

  if (update->ntimestep % 10000 == 0) {
    if (universe->me == 0) printf("This is the end of step %ld.\n", update->ntimestep);
  }
}

void FixPIMDLangevin::collect_xc()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  tagint *tag = atom->tag;
  if (ireplica == 0) {
    if (cmode == SINGLE_PROC) {
      for (int i = 0; i < nlocal; i++) {
        xcall[3 * i + 0] = xcall[3 * i + 1] = xcall[3 * i + 2] = 0.0;
      }
    } else if (cmode == MULTI_PROC) {
      for (int i = 0; i < ntotal; i++) {
        xcall[3 * i + 0] = xcall[3 * i + 1] = xcall[3 * i + 2] = 0.0;
      }
    }

    for (int i = 0; i < nlocal; i++) {
      xcall[3 * (tag[i] - 1) + 0] = x[i][0] / sqrt(np);
      xcall[3 * (tag[i] - 1) + 1] = x[i][1] / sqrt(np);
      xcall[3 * (tag[i] - 1) + 2] = x[i][2] / sqrt(np);
    }

    if (cmode == MULTI_PROC) {
      MPI_Allreduce(MPI_IN_PLACE, xcall, ntotal * 3, MPI_DOUBLE, MPI_SUM, world);
    }
  }
  MPI_Bcast(xcall, ntotal * 3, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::b_step()
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
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::qc_step()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double oldlo, oldhi;

  if (!pstat_flag) {
    if (universe->iworld == 0) {
      for (int i = 0; i < nlocal; i++) {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  } else {
    if (universe->iworld == 0) {
      double expp[3], expq[3];
      if (pstyle == ISO) {
        vw[1] = vw[0];
        vw[2] = vw[0];
      }
      for (int j = 0; j < 3; j++) {
        expq[j] = exp(dtv * vw[j]);
        expp[j] = exp(-dtv * vw[j]);
      }
      if (barostat == BZP) {
        for (int i = 0; i < nlocal; i++) {
          for (int j = 0; j < 3; j++) {
            x[i][j] = expq[j] * x[i][j] + (expq[j] - expp[j]) / 2. / vw[j] * v[i][j];
            v[i][j] = expp[j] * v[i][j];
          }
        }
        oldlo = domain->boxlo[0];
        oldhi = domain->boxhi[0];

        domain->boxlo[0] = (oldlo - fixedpoint[0]) * expq[0] + fixedpoint[0];
        domain->boxhi[0] = (oldhi - fixedpoint[0]) * expq[0] + fixedpoint[0];

        oldlo = domain->boxlo[1];
        oldhi = domain->boxhi[1];
        domain->boxlo[1] = (oldlo - fixedpoint[1]) * expq[1] + fixedpoint[1];
        domain->boxhi[1] = (oldhi - fixedpoint[1]) * expq[1] + fixedpoint[1];

        oldlo = domain->boxlo[2];
        oldhi = domain->boxhi[2];
        domain->boxlo[2] = (oldlo - fixedpoint[2]) * expq[2] + fixedpoint[2];
        domain->boxhi[2] = (oldhi - fixedpoint[2]) * expq[2] + fixedpoint[2];
      }
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&domain->boxlo[0], 3, MPI_DOUBLE, 0, universe->uworld);
    MPI_Bcast(&domain->boxhi[0], 3, MPI_DOUBLE, 0, universe->uworld);
    domain->set_global_box();
    domain->set_local_box();
  }
  volume = domain->xprd * domain->yprd * domain->zprd;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::a_step()
{
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double x0, x1, x2, v0, v1, v2;    // three components of x[i] and v[i]

  if (universe->iworld != 0) {
    for (int i = 0; i < n; i++) {
      x0 = x[i][0];
      x1 = x[i][1];
      x2 = x[i][2];
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
      x[i][0] = Lan_c[universe->iworld] * x0 +
          1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v0;
      x[i][1] = Lan_c[universe->iworld] * x1 +
          1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v1;
      x[i][2] = Lan_c[universe->iworld] * x2 +
          1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v2;
      v[i][0] = -1.0 * _omega_k[universe->iworld] * Lan_s[universe->iworld] * x0 +
          Lan_c[universe->iworld] * v0;
      v[i][1] = -1.0 * _omega_k[universe->iworld] * Lan_s[universe->iworld] * x1 +
          Lan_c[universe->iworld] * v1;
      v[i][2] = -1.0 * _omega_k[universe->iworld] * Lan_s[universe->iworld] * x2 +
          Lan_c[universe->iworld] * v2;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::baro_init()
{
  vw[0] = vw[1] = vw[2] = vw[3] = vw[4] = vw[5] = 0.0;
  if (pstyle == ISO) {
    W = 3 * (atom->natoms) * tau_p * tau_p * np * kBT;
  }    // consistent with the definition in i-Pi
  else if (pstyle == ANISO) {
    W = atom->natoms * tau_p * tau_p * np * kBT;
  }
  Vcoeff = 1.0;
  std::string out = fmt::format("\nInitializing PIMD {:s} barostat...\n", Barostats[barostat]);
  out += fmt::format("The barostat mass is W = {:.16e}\n", W);
  utils::logmesg(lmp, out);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::press_v_step()
{
  int nlocal = atom->nlocal;
  double **f = atom->f;
  double **v = atom->v;
  int *type = atom->type;
  volume = domain->xprd * domain->yprd * domain->zprd;

  if (pstyle == ISO) {
    if (barostat == BZP) {
      vw[0] += dtv * 3 * (volume * np * (p_cv - Pext) / force->nktv2p + Vcoeff / beta_np) / W;
      if (universe->iworld == 0) {
        double dvw_proc = 0.0, dvw = 0.0;
        for (int i = 0; i < nlocal; i++) {
          for (int j = 0; j < 3; j++) {
            dvw_proc += dtv2 * f[i][j] * v[i][j] / W + dtv3 * f[i][j] * f[i][j] / mass[type[i]] / W;
          }
        }
        MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
        vw[0] += dvw;
      }
      MPI_Barrier(universe->uworld);
      MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
    } else if (barostat == MTTK) {
      mtk_term1 = 2. / atom->natoms * totke / 3;
      f_omega = (volume * np * (p_md - Pext) + mtk_term1) / W;
      vw[0] += 0.5 * dtv * f_omega;
    }
  } else if (pstyle == ANISO) {
    compute_stress_tensor();
    for (int ii = 0; ii < 3; ii++) {
      vw[ii] +=
          dtv * (volume * np * (stress_tensor[ii] - Pext) / force->nktv2p + Vcoeff / beta_np) / W;
      if (universe->iworld == 0) {
        double dvw_proc = 0.0, dvw = 0.0;
        for (int i = 0; i < nlocal; i++) {
          dvw_proc +=
              dtv2 * f[i][ii] * v[i][ii] / W + dtv3 * f[i][ii] * f[i][ii] / mass[type[i]] / W;
        }
        MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
        vw[ii] += dvw;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::press_o_step()
{
  if (pstyle == ISO) {
    if (universe->me == 0) {
      r1 = random->gaussian();
      vw[0] = c1 * vw[0] + c2 * sqrt(1.0 / W / beta_np) * r1;
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
  } else if (pstyle == ANISO) {
    if (universe->me == 0) {
      r1 = random->gaussian();
      r2 = random->gaussian();
      r3 = random->gaussian();
      vw[0] = c1 * vw[0] + c2 * sqrt(1.0 / W / beta_np) * r1;
      vw[1] = c1 * vw[1] + c2 * sqrt(1.0 / W / beta_np) * r2;
      vw[2] = c1 * vw[2] + c2 * sqrt(1.0 / W / beta_np) * r3;
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw, 3, MPI_DOUBLE, 0, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::Langevin_init()
{
  double beta = 1.0 / kBT;
  _omega_np = np / beta / hbar;
  double _omega_np_dt_half = _omega_np * update->dt * 0.5;

  _omega_k = new double[np];
  Lan_c = new double[np];
  Lan_s = new double[np];
  if (fmmode == PHYSICAL) {
    for (int i = 0; i < np; i++) {
      _omega_k[i] = _omega_np * sqrt(lam[i]);
      Lan_c[i] = cos(sqrt(lam[i]) * _omega_np_dt_half);
      Lan_s[i] = sin(sqrt(lam[i]) * _omega_np_dt_half);
    }
  } else if (fmmode == NORMAL) {
    for (int i = 0; i < np; i++) {
      _omega_k[i] = _omega_np;
      Lan_c[i] = cos(_omega_np_dt_half);
      Lan_s[i] = sin(_omega_np_dt_half);
    }
  }
  if (tau > 0)
    gamma = 1.0 / tau;
  else
    gamma = np / beta / hbar;

  if (integrator == OBABO)
    c1 = exp(-gamma * 0.5 * update->dt);    // tau is the damping time of the centroid mode.
  else if (integrator == BAOAB)
    c1 = exp(-gamma * update->dt);
  else
    error->universe_all(FLERR,
                        "Unknown integrator parameter for fix pimd/langevin. Only obabo and "
                        "baoab integrators are supported!");

  c2 = sqrt(1.0 - c1 * c1);    // note that c1 and c2 here only works for the centroid mode.

  if (thermostat == PILE_L) {
    std::string out = "\nInitializing PI Langevin equation thermostat...\n";
    out += "Bead ID    |    omega    |    tau    |    c1    |    c2\n";
    tau_k = new double[np];
    c1_k = new double[np];
    c2_k = new double[np];
    tau_k[0] = tau;
    c1_k[0] = c1;
    c2_k[0] = c2;
    for (int i = 1; i < np; i++) {
      tau_k[i] = 0.5 / pilescale / _omega_k[i];
      if (integrator == OBABO)
        c1_k[i] = exp(-0.5 * update->dt / tau_k[i]);
      else if (integrator == BAOAB)
        c1_k[i] = exp(-1.0 * update->dt / tau_k[i]);
      else
        error->universe_all(FLERR,
                            "Unknown integrator parameter for fix pimd/langevin. Only obabo and "
                            "baoab integrators are supported!");
      c2_k[i] = sqrt(1.0 - c1_k[i] * c1_k[i]);
    }
    for (int i = 0; i < np; i++) {
      out += fmt::format("    {:d}     {:.8e} {:.8e} {:.8e} {:.8e}\n", i, _omega_k[i], tau_k[i],
                         c1_k[i], c2_k[i]);
    }
    if (thermostat == PILE_L) out += "PILE_L thermostat successfully initialized!\n";
    out += "\n";
    utils::logmesg(lmp, out);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::o_step()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double beta_np = 1.0 / force->boltz / Lan_temp * inverse_np * force->mvv2e;
  if (thermostat == PILE_L) {
    for (int i = 0; i < nlocal; i++) {
      r1 = random->gaussian();
      r2 = random->gaussian();
      r3 = random->gaussian();
      atom->v[i][0] = c1_k[universe->iworld] * atom->v[i][0] +
          c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r1;
      atom->v[i][1] = c1_k[universe->iworld] * atom->v[i][1] +
          c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r2;
      atom->v[i][2] = c1_k[universe->iworld] * atom->v[i][2] +
          c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * r3;
    }
  }
}

/* ----------------------------------------------------------------------
   Normal Mode PIMD
------------------------------------------------------------------------- */

void FixPIMDLangevin::nmpimd_init()
{
  memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");

  lam = (double *) memory->smalloc(sizeof(double) * np, "FixPIMDLangevin::lam");

  // Set up  eigenvalues
  for (int i = 0; i < np; i++) {
    double sin_tmp = sin(i * MY_PI / np);
    lam[i] = 4 * sin_tmp * sin_tmp;
  }

  // Set up eigenvectors for degenerated modes
  for (int j = 0; j < np; j++) {
    for (int i = 1; i < int(np / 2) + 1; i++) {
      M_x2xp[i][j] = sqrt(2.0) * cos(2.0 * MY_PI * double(i) * double(j) / double(np)) / sqrt(np);
    }
    for (int i = int(np / 2) + 1; i < np; i++) {
      M_x2xp[i][j] = sqrt(2.0) * sin(2.0 * MY_PI * double(i) * double(j) / double(np)) / sqrt(np);
    }
  }

  // Set up eigenvectors for non-degenerated modes
  for (int i = 0; i < np; i++) {
    M_x2xp[0][i] = 1.0 / sqrt(np);
    if (np % 2 == 0) M_x2xp[np / 2][i] = 1.0 / sqrt(np) * pow(-1.0, i);
  }

  // Set up Ut
  for (int i = 0; i < np; i++)
    for (int j = 0; j < np; j++) { M_xp2x[i][j] = M_x2xp[j][i]; }

  // Set up masses
  int iworld = universe->iworld;
  for (int i = 1; i <= atom->ntypes; i++) {
    mass[i] = atom->mass[i];
    if (iworld) {
      if (fmmode == PHYSICAL) {
        mass[i] *= 1.0;
      } else if (fmmode == NORMAL) {
        mass[i] *= lam[iworld];
      }
      mass[i] *= fmass;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::nmpimd_transform(double **src, double **des, double *vector)
{
  if (cmode == SINGLE_PROC) {
    for (int i = 0; i < ntotal; i++) {
      for (int d = 0; d < 3; d++) {
        bufsorted[i][d] = 0.0;
        for (int j = 0; j < nreplica; j++) {
          bufsorted[i][d] += src[j * ntotal + i][d] * vector[j];
        }
      }
    }
    for (int i = 0; i < ntotal; i++) {
      tagint tagtmp = atom->tag[i];
      for (int d = 0; d < 3; d++) { des[i][d] = bufsorted[tagtmp - 1][d]; }
    }
  } else if (cmode == MULTI_PROC) {
    int n = atom->nlocal;
    int m = 0;

    for (int i = 0; i < n; i++)
      for (int d = 0; d < 3; d++) {
        des[i][d] = 0.0;
        for (int j = 0; j < np; j++) { des[i][d] += (src[j][m] * vector[j]); }
        m++;
      }
  }
}

/* ----------------------------------------------------------------------
   Comm operations
------------------------------------------------------------------------- */

void FixPIMDLangevin::comm_init()
{
  // if(me == 0){
  if (sizeplan) {
    delete[] plansend;
    delete[] planrecv;
  }

  sizeplan = np - 1;
  plansend = new int[sizeplan];
  planrecv = new int[sizeplan];
  modeindex = new int[sizeplan];
  for (int i = 0; i < sizeplan; i++) {
    int isend, irecv;
    isend = ireplica + i + 1;
    if (isend >= nreplica) isend -= nreplica;
    irecv = ireplica - (i + 1);
    if (irecv < 0) irecv += nreplica;
    plansend[i] = universe->root_proc[isend];
    planrecv[i] = universe->root_proc[irecv];
    modeindex[i] = irecv;
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::reallocate_xc()
{
  maxxc = atom->nmax;
  memory->destroy(xc);
  memory->create(xc, maxxc, 3, "FixPIMDLangevin:xc");
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::reallocate_x_unwrap()
{
  maxunwrap = atom->nmax;
  memory->destroy(x_unwrap);
  memory->create(x_unwrap, maxunwrap, 3, "FixPIMDLangevin:x_unwrap");
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::reallocate()
{
  maxlocal = atom->nmax;
  memory->destroy(bufsend);
  memory->destroy(bufrecv);
  memory->destroy(tagsend);
  memory->destroy(tagrecv);
  memory->destroy(bufbeads);
  memory->create(bufsend, maxlocal, 3, "FixPIMDLangevin:bufsend");
  memory->create(bufrecv, maxlocal, 3, "FixPIMDLangevin:bufrecv");
  memory->create(tagsend, maxlocal, "FixPIMDLangevin:tagsend");
  memory->create(tagrecv, maxlocal, "FixPIMDLangevin:tagrecv");
  memory->create(bufbeads, nreplica, maxlocal * 3, "FixPIMDLangevin:bufrecv");
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::inter_replica_comm(double **ptr)
{
  MPI_Request requests[2];
  MPI_Status statuses[2];
  if (atom->nmax > maxlocal) reallocate();
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int i, m;

  // copy local values
  for (i = 0; i < nlocal; i++) {
    bufbeads[ireplica][3 * i + 0] = ptr[i][0];
    bufbeads[ireplica][3 * i + 1] = ptr[i][1];
    bufbeads[ireplica][3 * i + 2] = ptr[i][2];
  }

  // communicate values from the other beads
  if (cmode == SINGLE_PROC) {
    m = 0;
    for (i = 0; i < nlocal; i++) {
      tagint tagtmp = atom->tag[i];
      bufsorted[tagtmp - 1][0] = ptr[i][0];
      bufsorted[tagtmp - 1][1] = ptr[i][1];
      bufsorted[tagtmp - 1][2] = ptr[i][2];
      m++;
    }
    MPI_Allgather(&m, 1, MPI_INT, counts, 1, MPI_INT, universe->uworld);
    for (i = 0; i < nreplica; i++) { counts[i] *= 3; }
    displacements[0] = 0;
    for (i = 0; i < nreplica - 1; i++) { displacements[i + 1] = displacements[i] + counts[i]; }
    MPI_Allgatherv(bufsorted[0], 3 * m, MPI_DOUBLE, bufsortedall[0], counts, displacements,
                   MPI_DOUBLE, universe->uworld);
  } else if (cmode == MULTI_PROC) {
    m = 0;
    for (i = 0; i < nlocal; i++) {
      tagsend[m] = tag[i];
      bufsend[m][0] = ptr[i][0];
      bufsend[m][1] = ptr[i][1];
      bufsend[m][2] = ptr[i][2];
      m++;
    }
    MPI_Gather(&m, 1, MPI_INT, counts, 1, MPI_INT, 0, world);
    displacements[0] = 0;
    for (i = 0; i < nprocs - 1; i++) { displacements[i + 1] = displacements[i] + counts[i]; }
    MPI_Gatherv(tagsend, m, MPI_LMP_TAGINT, tagsendall, counts, displacements, MPI_LMP_TAGINT, 0,
                world);
    for (i = 0; i < nprocs; i++) { counts[i] *= 3; }
    for (i = 0; i < nprocs - 1; i++) { displacements[i + 1] = displacements[i] + counts[i]; }
    MPI_Gatherv(bufsend[0], 3 * m, MPI_DOUBLE, bufsendall[0], counts, displacements, MPI_DOUBLE, 0,
                world);
    for (int iplan = 0; iplan < sizeplan; iplan++) {
      if (me == 0) {
        MPI_Irecv(bufrecvall[0], 3 * ntotal, MPI_DOUBLE, planrecv[iplan], 0, universe->uworld,
                  &requests[0]);
        MPI_Irecv(tagrecvall, ntotal, MPI_LMP_TAGINT, planrecv[iplan], 0, universe->uworld,
                  &requests[1]);
        MPI_Send(bufsendall[0], 3 * ntotal, MPI_DOUBLE, plansend[iplan], 0, universe->uworld);
        MPI_Send(tagsendall, ntotal, MPI_LMP_TAGINT, plansend[iplan], 0, universe->uworld);
        MPI_Waitall(2, requests, statuses);
      }
      MPI_Bcast(tagrecvall, ntotal, MPI_LMP_TAGINT, 0, world);
      MPI_Bcast(bufrecvall[0], 3 * ntotal, MPI_DOUBLE, 0, world);
      for (i = 0; i < ntotal; i++) {
        m = atom->map(tagrecvall[i]);
        if (m < 0 || m >= nlocal) continue;
        bufbeads[modeindex[iplan]][3 * m + 0] = bufrecvall[i][0];
        bufbeads[modeindex[iplan]][3 * m + 1] = bufrecvall[i][1];
        bufbeads[modeindex[iplan]][3 * m + 2] = bufrecvall[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_vir_()
{
  int nlocal = atom->nlocal;
  xf = vir_ = xcf = centroid_vir = 0.0;
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < 3; j++) {
      xf += x_unwrap[i][j] * atom->f[i][j];
      xcf += (x_unwrap[i][j] - xc[i][j]) * atom->f[i][j];
    }
  }
  MPI_Allreduce(&xf, &vir_, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&xcf, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  if (pstyle == ANISO) {
    for (int i = 0; i < 6; i++) c_vir_tensor[i] = 0.0;
    for (int i = 0; i < nlocal; i++) {
      c_vir_tensor[0] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][0];
      c_vir_tensor[1] += (x_unwrap[i][1] - xc[i][1]) * atom->f[i][1];
      c_vir_tensor[2] += (x_unwrap[i][2] - xc[i][2]) * atom->f[i][2];
      c_vir_tensor[3] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][1];
      c_vir_tensor[4] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][2];
      c_vir_tensor[5] += (x_unwrap[i][1] - xc[i][1]) * atom->f[i][2];
    }
    MPI_Allreduce(MPI_IN_PLACE, &c_vir_tensor, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_vir()
{
  volume = domain->xprd * domain->yprd * domain->zprd;
  c_press->compute_vector();
  virial[0] = c_press->vector[0] * volume;
  virial[1] = c_press->vector[1] * volume;
  virial[2] = c_press->vector[2] * volume;
  virial[3] = c_press->vector[3] * volume;
  virial[4] = c_press->vector[4] * volume;
  virial[5] = c_press->vector[5] * volume;
  for (int i = 0; i < 6; i++) virial[i] /= universe->procs_per_world[universe->iworld];
  double vir_bead = (virial[0] + virial[1] + virial[2]);
  MPI_Allreduce(&vir_bead, &vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(MPI_IN_PLACE, &virial[0], 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_stress_tensor()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  if (universe->iworld == 0) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    for (int i = 0; i < 6; i++) ke_tensor[i] = 0.0;
    for (int i = 0; i < nlocal; i++) {
      ke_tensor[0] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][0] * force->mvv2e;
      ke_tensor[1] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][1] * force->mvv2e;
      ke_tensor[2] += 0.5 * mass[type[i]] * atom->v[i][2] * atom->v[i][2] * force->mvv2e;
      ke_tensor[3] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][1] * force->mvv2e;
      ke_tensor[4] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][2] * force->mvv2e;
      ke_tensor[5] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][2] * force->mvv2e;
    }
    MPI_Allreduce(MPI_IN_PLACE, &ke_tensor, 6, MPI_DOUBLE, MPI_SUM, world);
    for (int i = 0; i < 6; i++) {
      stress_tensor[i] =
          inv_volume * ((2 * ke_tensor[i] - c_vir_tensor[i]) * force->nktv2p + virial[i]) / np;
    }
  }
  MPI_Bcast(&stress_tensor, 6, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_totke()
{
  kine = 0.0;
  totke = ke_bead = 0.0;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < 3; j++) { kine += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j]; }
  }
  kine *= force->mvv2e;
  MPI_Allreduce(&kine, &ke_bead, 1, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_spring_energy()
{
  spring_energy = 0.0;
  total_spring_energy = se_bead = 0.0;

  double **x = atom->x;
  double *_mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    spring_energy += 0.5 * _mass[type[i]] * fbond * lam[universe->iworld] *
        (x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]);
  }
  MPI_Allreduce(&spring_energy, &se_bead, 1, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_pote()
{
  pe_bead = 0.0;
  pot_energy_partition = 0.0;
  pote = 0.0;
  c_pe->compute_scalar();
  pe_bead = c_pe->scalar;
  pot_energy_partition = pe_bead / universe->procs_per_world[universe->iworld];
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_tote()
{
  tote = totke + pote + total_spring_energy;
}

/* ---------------------------------------------------------------------- */
/*
void FixPIMDLangevin::compute_t_prim()
{
  t_prim = 1.5 * atom->natoms * np * force->boltz * temp - total_spring_energy;
}
*/
/* ---------------------------------------------------------------------- */
/*
void FixPIMDLangevin::compute_t_vir()
{
  t_vir = -0.5 * inverse_np * vir_;
  t_cv = 1.5 * atom->natoms * force->boltz * temp - 0.5 * inverse_np * centroid_vir;
}
*/
/* ---------------------------------------------------------------------- */
/*
void FixPIMDLangevin::compute_p_prim()
{
  p_prim = atom->natoms * np * force->boltz * temp * inv_volume - 1.0 / 1.5 * inv_volume * total_spring_energy;
  p_prim *= force->nktv2p;
}
*/
/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_p_cv()
{
  inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  if (universe->iworld == 0) {
    p_cv = THIRD * inv_volume * ((2.0 * ke_bead - centroid_vir) * force->nktv2p + vir) / np;
  }
  MPI_Bcast(&p_cv, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_totenthalpy()
{
  volume = domain->xprd * domain->yprd * domain->zprd;
  if (barostat == BZP) {
    if (pstyle == ISO) {
      totenthalpy = tote + 0.5 * W * vw[0] * vw[0] * inverse_np + Pext * volume / force->nktv2p -
          Vcoeff * kBT * log(volume);
    } else if (pstyle == ANISO) {
      totenthalpy = tote + 0.5 * W * vw[0] * vw[0] * inverse_np +
          0.5 * W * vw[1] * vw[1] * inverse_np + 0.5 * W * vw[2] * vw[2] * inverse_np +
          Pext * volume / force->nktv2p - Vcoeff * kBT * log(volume);
    }
  } else if (barostat == MTTK)
    totenthalpy = tote + 1.5 * W * vw[0] * vw[0] * inverse_np + Pext * (volume - vol0);
}

/* ---------------------------------------------------------------------- */

double FixPIMDLangevin::compute_vector(int n)
{
  if (n == 0) return ke_bead;
  if (n == 1) return se_bead;
  if (n == 2) return pe_bead;
  if (n == 3) return tote;
  // if(n==3) { return W*vw[0]; }
  if (!pstat_flag) {
    if (n == 4) return t_prim;
    if (n == 5) return t_vir;
    if (n == 6) return t_cv;
  } else if (pstat_flag) {
    if (pstyle == ISO) {
      if (barostat == BZP) {
        if (n == 4) return 0.5 * W * vw[0] * vw[0];
      } else if (barostat == MTTK) {
        if (n == 4) return 1.5 * W * vw[0] * vw[0];
      }
      if (n == 5) {
        volume = domain->xprd * domain->yprd * domain->zprd;
        return np * Pext * volume / force->nktv2p;
      }
      if (n == 6) {
        volume = domain->xprd * domain->yprd * domain->zprd;
        return -Vcoeff * np * kBT * log(volume);
      }
      if (n == 7) return totenthalpy;
      if (n == 8) return p_cv;
      if (n == 9) return total_spring_energy;
    } else if (pstyle == ANISO) {
    }
  }

  /*

  if(n==7) { return p_prim; }
  if(n==8) { return p_md; }
  if(n==9) { return p_cv; }
  if(n==10) {return totenthalpy;
  if(pstyle == ISO){
    if(n==11) { return vw[0]; }
    if(n==12) {
      if(barostat == BZP) {  return 0.5*W*vw[0]*vw[0]; }
      else if(barostat == MTTK) {  return 1.5*W*vw[0]*vw[0]; }
    }
    if(n==13) { volume = domain->xprd * domain->yprd * domain->zprd; return np * Pext * volume / force->nktv2p; }
    if(n==14) { volume = domain->xprd * domain->yprd * domain->zprd;
    // printf("Vcoeff = %.6e np = %d kBT = %.6e logV = %.6e\n", Vcoeff, np, kBT, log(volume));
    return - Vcoeff * np * kBT * log(volume); }
  }
  else if(pstyle==ANISO){
    if(n>10 && n<=13) return vw[n-11];
    if(n==14) return 0.5*W*vw[0]*vw[0]+0.5*W*vw[1]*vw[1]+0.5*W*vw[2]*vw[2];
    if(n>14 && n<21) return stress_tensor[n-15];
    if(n==21) { volume = domain->xprd * domain->yprd * domain->zprd; return np * Pext * volume / force->nktv2p; }
    if(n==22) { volume = domain->xprd * domain->yprd * domain->zprd; return - Vcoeff * np * kBT * log(volume); }
  }  */
  return 0.0;
}

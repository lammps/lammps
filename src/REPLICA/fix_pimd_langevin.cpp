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
   Purpose      Path Integral Molecular Dynamics with Langevin Thermostat

   Yifan Li @ Princeton University (yifanl0716@gmail.com)
   Current Features:
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
#include "kspace.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;
using MathConst::MY_PI;
using MathConst::MY_SQRT2;
using MathConst::THIRD;
using MathSpecial::powint;

namespace {
std::map<int, std::string> Barostats{{FixPIMDLangevin::MTTK, "MTTK"},
                                     {FixPIMDLangevin::BZP, "BZP"}};
std::map<int, std::string> Ensembles{{FixPIMDLangevin::NVE, "NVE"},
                                     {FixPIMDLangevin::NVT, "NVT"},
                                     {FixPIMDLangevin::NPH, "NPH"},
                                     {FixPIMDLangevin::NPT, "NPT"}};
}    // namespace

/* ---------------------------------------------------------------------- */

FixPIMDLangevin::FixPIMDLangevin(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDNVE(lmp, narg, arg, true), tau_k(nullptr),
    c1_k(nullptr), c2_k(nullptr), random(nullptr)
{

  ensemble = NVT;
  thermostat = PILE_L;
  barostat = BZP;
  Lan_temp = 298.15;
  tau = 1.0;
  tau_p = 1.0;
  Pext = 1.0;
  pdim = 0;
  pilescale = 1.0;
  tstat_flag = 1;
  pstat_flag = 0;
  pstyle = ISO;
  totenthalpy = 0.0;

  seed = -1;

  for (int i = 0; i < 6; i++) {
    p_flag[i] = 0;
    p_target[i] = 0.0;
  }

  // process keywords

  for (int i = 3; i < narg;) {
    if (!parse_keyword(narg, arg, i))
      error->universe_all(FLERR, fmt::format("Unknown keyword {} for fix {}", arg[i], style));
  }

  if (pstat_flag && !pdim)
    error->universe_all(
        FLERR, fmt::format("Must use pressure coupling with {} ensemble", Ensembles[ensemble]));
  if (!pstat_flag && pdim)
    error->universe_all(
        FLERR, fmt::format("Must not use pressure coupling with {} ensemble", Ensembles[ensemble]));

  if (method == PIMD && pstat_flag)
    error->universe_all(FLERR,
                        "Pressure control has not been supported for method pimd yet. Please set "
                        "method to nmpimd.");

  if (method == PIMD && fmmode == NORMAL)
    error->universe_all(
        FLERR, "Normal mode mass is not supported for method pimd. Please set method to nmpimd.");

  /* Initiation */

  if (!pstat_flag) {
    size_vector = 10;
  } else if (pstat_flag) {
    if (pstyle == ISO) {
      size_vector = 15;
    } else if (pstyle == ANISO) {
      size_vector = 17;
    }
  }
  extvector = 1;
  kt = force->boltz * temp;
  if (pstat_flag) FixPIMDLangevin::baro_init();

  finish_constructor_setup();

  vol0 = domain->xprd * domain->yprd * domain->zprd;

  fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
  fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
  fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);
  if (pstat_flag) p_hydro = (p_target[0] + p_target[1] + p_target[2]) / pdim;

  // initialize Marsaglia RNG with processor-unique seed

  if (tstat_flag) {
    if (integrator == BAOAB || integrator == OBABO) {
      Lan_temp = temp;
      random = new RanMars(lmp, seed + universe->me);
    }
  }
}

/* ---------------------------------------------------------------------- */

bool FixPIMDLangevin::parse_keyword(int narg, char **arg, int &i)
{
  if (i + 2 > narg)
    utils::missing_cmd_args(FLERR, fmt::format("fix {} {}", style, arg[i]), error);
  // The bosonic argument filter replaces its own keyword with an empty string.
  if (arg[i][0] == '\0') {
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "method") == 0) {
    if (strcmp(arg[i + 1], "nmpimd") == 0)
      method = NMPIMD;
    else if (strcmp(arg[i + 1], "pimd") == 0)
      method = PIMD;
    else
      error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
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
                          fmt::format("Unknown ensemble parameter for fix {}. Only nve, nvt, "
                                      "nph, and npt ensembles are supported!",
                                      style));
  } else if (strcmp(arg[i], "scale") == 0) {
    if (method == PIMD)
      error->universe_all(
          FLERR,
          "The scale parameter of the PILE_L thermostat is not supported for method pimd. Delete "
          "scale parameter if you do want to use method pimd.");
    pilescale = utils::numeric(FLERR, arg[i + 1], false, lmp);
    if (pilescale < 0.0)
      error->universe_all(FLERR, fmt::format("Invalid PILE_L scale value for fix {}", style));
  } else if (strcmp(arg[i], "thermostat") == 0) {
    if (strcmp(arg[i + 1], "PILE_L") == 0) {
      if (i + 3 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} thermostat", style), error);
      thermostat = PILE_L;
      seed = utils::inumeric(FLERR, arg[i + 2], false, lmp);
      i++;
    } else
      error->universe_all(FLERR,
                          fmt::format("Unknown thermostat parameter {} for fix {}; only PILE_L is supported",
                                      arg[i + 1], style));
  } else if (strcmp(arg[i], "tau") == 0) {
    tau = utils::numeric(FLERR, arg[i + 1], false, lmp);
  } else if (strcmp(arg[i], "barostat") == 0) {
    if (strcmp(arg[i + 1], "MTTK") == 0) {
      barostat = MTTK;
    } else if (strcmp(arg[i + 1], "BZP") == 0) {
      barostat = BZP;
    } else
      error->universe_all(FLERR, fmt::format("Unknown barostat parameter for fix {}", style));
  } else if (strcmp(arg[i], "iso") == 0) {
    pstyle = ISO;
    p_flag[0] = p_flag[1] = p_flag[2] = 1;
    Pext = utils::numeric(FLERR, arg[i + 1], false, lmp);
    p_target[0] = p_target[1] = p_target[2] = Pext;
    pdim = 3;
  } else if (strcmp(arg[i], "aniso") == 0) {
    pstyle = ANISO;
    p_flag[0] = p_flag[1] = p_flag[2] = 1;
    Pext = utils::numeric(FLERR, arg[i + 1], false, lmp);
    p_target[0] = p_target[1] = p_target[2] = Pext;
    pdim = 3;
  } else if (strcmp(arg[i], "x") == 0) {
    pstyle = ANISO;
    p_flag[0] = 1;
    p_target[0] = utils::numeric(FLERR, arg[i + 1], false, lmp);
    pdim++;
  } else if (strcmp(arg[i], "y") == 0) {
    pstyle = ANISO;
    p_flag[1] = 1;
    p_target[1] = utils::numeric(FLERR, arg[i + 1], false, lmp);
    pdim++;
  } else if (strcmp(arg[i], "z") == 0) {
    pstyle = ANISO;
    p_flag[2] = 1;
    p_target[2] = utils::numeric(FLERR, arg[i + 1], false, lmp);
    pdim++;
  } else if (strcmp(arg[i], "taup") == 0) {
    tau_p = utils::numeric(FLERR, arg[i + 1], false, lmp);
    if (tau_p <= 0.0)
      error->universe_all(FLERR, fmt::format("Invalid tau_p value for fix {}", style));
  } else if (strcmp(arg[i], "fixcom") == 0) {
    if (strcmp(arg[i + 1], "yes") == 0)
      removecomflag = 1;
    else if (strcmp(arg[i + 1], "no") == 0)
      removecomflag = 0;
  } else if (strcmp(arg[i], "lj") == 0 || strcmp(arg[i], "removecom") == 0) {
    // These base options are not part of the Langevin command syntax.
    error->universe_all(FLERR, fmt::format("Unknown keyword {} for fix {}", arg[i], style));
  } else {
    return FixPIMDNVE::parse_keyword(narg, arg, i);
  }
  i += 2;
  return true;
}

/* ---------------------------------------------------------------------- */

FixPIMDLangevin::~FixPIMDLangevin()
{
  delete random;
  delete[] tau_k;
  delete[] c1_k;
  delete[] c2_k;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::setup_subclass_state()
{
  langevin_init();
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::b_step()
{
  if (pstat_flag) {
    compute_totke();
    compute_p_cv();
    press_v_step();
  }
  FixPIMDNVE::b_step();
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::end_of_step()
{
  FixPIMDNVE::end_of_step();
  if (pstat_flag) compute_totenthalpy();
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::qc_step()
{
  // used for NMPIMD
  // evolve the centroid mode
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double oldlo, oldhi;

  if (!pstat_flag) {
    FixPIMDNVE::qc_step();
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
          if (mask[i] & groupbit) {
            for (int j = 0; j < 3; j++) {
              if (p_flag[j]) {
                x[i][j] = expq[j] * x[i][j] + (expq[j] - expp[j]) / 2. / vw[j] * v[i][j];
                v[i][j] = expp[j] * v[i][j];
              } else {
                x[i][j] += dtv * v[i][j];
              }
            }
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
    if (force->kspace) force->kspace->setup();
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::baro_init()
{
  vw[0] = vw[1] = vw[2] = vw[3] = vw[4] = vw[5] = 0.0;
  if (pstyle == ISO) {
    W = 3 * (group->count(igroup)) * tau_p * tau_p * np * kt;
  }    // consistent with the definition in i-Pi
  else if (pstyle == ANISO) {
    W = group->count(igroup) * tau_p * tau_p * np * kt;
  }
  Vcoeff = 1.0;
  std::string out = fmt::format("\nInitializing PIMD {:s} barostat...\n", Barostats[barostat]);
  out += fmt::format("  The barostat mass is W = {:.16e}\n", W);
  if (universe->me == 0) utils::logmesg(lmp, out);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::press_v_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **f = atom->f;
  double **v = atom->v;
  int *type = atom->type;
  double volume = domain->xprd * domain->yprd * domain->zprd;

  if (pstyle == ISO) {
    if (barostat == BZP) {
      vw[0] += dtv * 3 * (volume * np * (p_cv - p_hydro) / force->nktv2p + Vcoeff / beta_np) / W;
      if (universe->iworld == 0) {
        double dvw_proc = 0.0, dvw = 0.0;
        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            for (int j = 0; j < 3; j++) {
              dvw_proc += dtv2 * f[i][j] * v[i][j] / W + dtv3 * f[i][j] * f[i][j] / mass[type[i]] / W;
            }
          }
        }
        MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
        vw[0] += dvw;
      }
      MPI_Barrier(universe->uworld);
      MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
    } else if (barostat == MTTK) {
      double mtk_term1 = 2.0 / group->count(igroup) * totke / 3.0;
      vw[0] += 0.5 * dtv * (volume * np * (p_md - p_hydro) + mtk_term1) / W;
    }
  } else if (pstyle == ANISO) {
    compute_stress_tensor();
    for (int ii = 0; ii < 3; ii++) {
      if (p_flag[ii]) {
        vw[ii] += dtv *
            (volume * np * (stress_tensor[ii] - p_hydro) / force->nktv2p + Vcoeff / beta_np) / W;
        if (universe->iworld == 0) {
          double dvw_proc = 0.0, dvw = 0.0;
          for (int i = 0; i < nlocal; i++) {
            if (mask[i] & groupbit) {
              dvw_proc +=
                  dtv2 * f[i][ii] * v[i][ii] / W + dtv3 * f[i][ii] * f[i][ii] / mass[type[i]] / W;
            }
          }
          MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
          vw[ii] += dvw;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::press_o_step()
{
  if (pstyle == ISO) {
    if (universe->me == 0) vw[0] = c1 * vw[0] + c2 * sqrt(1.0 / W / beta_np) * random->gaussian();
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
  } else if (pstyle == ANISO) {
    if (universe->me == 0) {
      for (int ii = 0; ii < 3; ii++) {
        if (p_flag[ii]) vw[ii] = c1 * vw[ii] + c2 * sqrt(1.0 / W / beta_np) * random->gaussian();
      }
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw, 3, MPI_DOUBLE, 0, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::langevin_init()
{
  double beta = 1.0 / kt;
  const double _omega_np = np / beta / hbar;
  const double _omega_np_dt_half = _omega_np * update->dt * 0.5;

  if (method == NMPIMD) {
    if (fmmode == PHYSICAL) {
      for (int i = 0; i < np; i++) {
        _omega_k[i] = _omega_np * sqrt(lam[i]) / sqrt(fmass);
        Lan_c[i] = cos(sqrt(lam[i]) * _omega_np_dt_half);
        Lan_s[i] = sin(sqrt(lam[i]) * _omega_np_dt_half);
      }
    } else if (fmmode == NORMAL) {
      for (int i = 0; i < np; i++) {
        _omega_k[i] = _omega_np / sqrt(fmass);
        Lan_c[i] = cos(_omega_np_dt_half);
        Lan_s[i] = sin(_omega_np_dt_half);
      }
    } else {
      error->universe_all(FLERR, "Unknown fmmode setting; only physical and normal are supported!");
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
                        fmt::format("Unknown integrator parameter for fix {}. Only obabo and baoab "
                                    "integrators are supported!",
                                    style));

  c2 = sqrt(1.0 - c1 * c1);    // note that c1 and c2 here only works for the centroid mode.

  if (thermostat == PILE_L) {
    std::string out = "Initializing PI Langevin equation thermostat...\n";
    out += "  Bead ID    |    omega    |    tau    |    c1    |    c2\n";
    if (method == NMPIMD) {
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
                              fmt::format("Unknown integrator parameter for fix {}. Only obabo and "
                                          "baoab integrators are supported!",
                                          style));
        c2_k[i] = sqrt(1.0 - c1_k[i] * c1_k[i]);
      }
      for (int i = 0; i < np; i++) {
        out += fmt::format("      {:d}     {:.8e} {:.8e} {:.8e} {:.8e}\n", i, _omega_k[i], tau_k[i],
                           c1_k[i], c2_k[i]);
      }
    } else if (method == PIMD) {
      for (int i = 0; i < np; i++) {
        out += fmt::format("      {:d}     {:.8e} {:.8e} {:.8e} {:.8e}\n", i,
                           _omega_np / sqrt(fmass), tau, c1, c2);
      }
    }
    if (thermostat == PILE_L) out += "  PILE_L thermostat successfully initialized!\n";
    out += "\n";
    if (universe->me == 0) utils::logmesg(lmp, out);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::o_step()
{
  if (!tstat_flag) return;

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double beta_np = 1.0 / force->boltz / Lan_temp * inverse_np * force->mvv2e;
  if (thermostat == PILE_L) {
    if (method == NMPIMD) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          atom->v[i][0] = c1_k[universe->iworld] * atom->v[i][0] +
              c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][1] = c1_k[universe->iworld] * atom->v[i][1] +
              c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][2] = c1_k[universe->iworld] * atom->v[i][2] +
              c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
        }
      }
    } else if (method == PIMD) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          atom->v[i][0] =
              c1 * atom->v[i][0] + c2 * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][1] =
              c1 * atom->v[i][1] + c2 * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][2] =
              c1 * atom->v[i][2] + c2 * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
        }
      }
    }
  }
  if (removecomflag) remove_com_motion();
  if (pstat_flag) press_o_step();
}

/* ----------------------------------------------------------------------
   Normal Mode PIMD
   ------------------------------------------------------------------------- */

void FixPIMDLangevin::remove_com_motion()
{
  if (method == NMPIMD) {
    FixPIMDNVE::remove_com_motion();
  } else if (method == PIMD) {
    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (dynamic) masstotal = group->mass(igroup);
    double vcm[3];
    group->vcm(igroup, masstotal, vcm);
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] -= vcm[0];
        v[i][1] -= vcm[1];
        v[i][2] -= vcm[2];
      }
    }
  } else {
    error->all(
        FLERR,
        fmt::format("Unknown method for fix {}. Only nmpimd and pimd are supported!", style));
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_cvir()
{
  FixPIMDNVE::compute_cvir();
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  if (pstyle == ANISO) {
    for (int i = 0; i < 6; i++) c_vir_tensor[i] = 0.0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        c_vir_tensor[0] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][0];
        c_vir_tensor[1] += (x_unwrap[i][1] - xc[i][1]) * atom->f[i][1];
        c_vir_tensor[2] += (x_unwrap[i][2] - xc[i][2]) * atom->f[i][2];
        c_vir_tensor[3] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][1];
        c_vir_tensor[4] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][2];
        c_vir_tensor[5] += (x_unwrap[i][1] - xc[i][1]) * atom->f[i][2];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &c_vir_tensor, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_stress_tensor()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  if (universe->iworld == 0) {
    double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    for (int i = 0; i < 6; i++) ke_tensor[i] = 0.0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        ke_tensor[0] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][0] * force->mvv2e;
        ke_tensor[1] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][1] * force->mvv2e;
        ke_tensor[2] += 0.5 * mass[type[i]] * atom->v[i][2] * atom->v[i][2] * force->mvv2e;
        ke_tensor[3] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][1] * force->mvv2e;
        ke_tensor[4] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][2] * force->mvv2e;
        ke_tensor[5] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][2] * force->mvv2e;
      }
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

void FixPIMDLangevin::compute_totenthalpy()
{
  double volume = domain->xprd * domain->yprd * domain->zprd;
  if (barostat == BZP) {
    if (pstyle == ISO) {
      totenthalpy = tote + 0.5 * W * vw[0] * vw[0] * inverse_np + p_hydro * volume / force->nktv2p -
          Vcoeff * kt * log(volume);
    } else if (pstyle == ANISO) {
      totenthalpy = tote + 0.5 * W * vw[0] * vw[0] * inverse_np +
          0.5 * W * vw[1] * vw[1] * inverse_np + 0.5 * W * vw[2] * vw[2] * inverse_np +
          p_hydro * volume / force->nktv2p - Vcoeff * kt * log(volume);
    }
  } else if (barostat == MTTK)
    totenthalpy = tote + 1.5 * W * vw[0] * vw[0] * inverse_np + p_hydro * (volume - vol0);
}

double FixPIMDLangevin::compute_subclass_vector(int n) const
{
  if (!pstat_flag) return 0.0;

  const double volume = domain->xprd * domain->yprd * domain->zprd;

  if (pstyle == ISO) {
    if (n == 0) return vw[0];
    if (n == 1) {
      if (barostat == BZP) return 0.5 * W * vw[0] * vw[0];
      if (barostat == MTTK) return 1.5 * W * vw[0] * vw[0];
      return 0.0;
    }
    if (n == 2) return np * Pext * volume / force->nktv2p;
    if (n == 3) return -Vcoeff * np * kt * log(volume);
    if (n == 4) return totenthalpy;
  } else if (pstyle == ANISO) {
    if (n == 0) return vw[0];
    if (n == 1) return vw[1];
    if (n == 2) return vw[2];
    if (n == 3) return 0.5 * W * (vw[0] * vw[0] + vw[1] * vw[1] + vw[2] * vw[2]);
    if (n == 4) return np * Pext * volume / force->nktv2p;
    if (n == 5) return -Vcoeff * np * kt * log(volume);
    if (n == 6) return totenthalpy;
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixPIMDLangevin::base_restart_size() const
{
  return 6;
}

/* ---------------------------------------------------------------------- */

int FixPIMDLangevin::pack_base_restart(double *list) const
{
  int n = 0;
  for (int i = 0; i < 6; i++) list[n++] = vw[i];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixPIMDLangevin::unpack_base_restart(const double *list)
{
  int n = 0;
  for (int i = 0; i < 6; i++) vw[i] = list[n++];
  return n;
}

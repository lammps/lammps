/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "gran_sub_mod_tangential.h"

#include "error.h"
#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_normal.h"
#include "granular_model.h"
#include "math_extra.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathExtra;

static constexpr double EPSILON = 1e-10;

/* ----------------------------------------------------------------------
   Default model
------------------------------------------------------------------------- */

GranSubModTangential::GranSubModTangential(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  allow_synchronization = 0;
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModTangentialNone::GranSubModTangentialNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  allow_synchronization = 1;
}

/* ----------------------------------------------------------------------
   Linear model with no history
------------------------------------------------------------------------- */

GranSubModTangentialLinearNoHistory::GranSubModTangentialLinearNoHistory(GranularModel *gm,
                                                                         LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 2;
  size_history = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::coeffs_to_local()
{
  k = 0.0;    // No tangential stiffness with no history
  xt = coeffs[0];
  mu = coeffs[1];

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal linear no history tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::calculate_forces()
{
  // classic pair gran/hooke (no history)
  damp = xt * gm->damping_model->get_damp_prefactor();

  double vrel = gm->vrel;
  double *vtr = gm->vtr;
  double *fs = gm->fs;
  double Fscrit = mu * gm->normal_model->get_fncrit();
  double fsmag = damp * vrel;

  double Ft;
  if (vrel != 0.0)
    Ft = MIN(Fscrit, fsmag) / vrel;
  else
    Ft = 0.0;

  scale3(-Ft, vtr, fs);
}

/* ----------------------------------------------------------------------
   Linear model with history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistory::GranSubModTangentialLinearHistory(GranularModel *gm,
                                                                     LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  allow_synchronization = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || xt < 0.0 || mu < 0.0) error->all(FLERR, "Illegal linear tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::calculate_forces()
{
  // Note: this is the same as the base Mindlin calculation except k isn't scaled by contact radius
  double magfs, rsht, shrmag;
  double hist_increment[3], fdamp[3], vtr2[3];
  int frame_update = 0;

  double *nx = gm->nx;
  double *nx_unrotated = gm->nx_unrotated;
  double *vtr = gm->vtr;
  double *fs = gm->fs;
  double dt = gm->dt;
  double *history = &gm->history[history_index];
  int history_update = gm->history_update;

  damp = xt * gm->damping_model->get_damp_prefactor();
  double Fscrit = gm->normal_model->get_fncrit() * mu;

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = dot3(history, nx);
    frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);

    if (frame_update) rotate_rescale_vec(history, nx);

    // update history, tangential force using velocities at half step
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    scale3(dt, vtr, hist_increment);
    add3(history, hist_increment, history);

    if(gm->synchronized_verlet == 1) {
      rsht = dot3(history, nx_unrotated);
      frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);
      //Second projection to nx (t+\Delta t)
      if (frame_update) rotate_rescale_vec(history, nx_unrotated);
    }
  }

  // tangential forces = history + tangential velocity damping
  scale3(-k, history, fs);
  //Rotating vtr for damping term in nx direction
  if (frame_update && gm->synchronized_verlet == 1) {
    copy3(vtr, vtr2);
    rotate_rescale_vec(vtr2, nx_unrotated);
  } else {
    copy3(vtr, vtr2);
  }
  scale3(-damp, vtr2, fdamp);
  add3(fs, fdamp, fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      // Rescale shear force
      scale3(Fscrit / magfs, fs);

      // Set shear to elastic component of rescaled force
      //  has extra factor of k that is then removed
      sub3(fs, fdamp, history);
      scale3(-1.0 / k, history);
    } else {
      zero3(fs);
    }
  }
}


//**********************************************
// Linear history w/ static
//**********************************************

GranSubModTangentialLinearHistoryStatic::GranSubModTangentialLinearHistoryStatic(GranularModel *gm, LAMMPS *lmp) :
		GranSubModTangentialLinearHistory(gm, lmp)
{
  num_coeffs = 4;
  size_history = 4;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistoryStatic::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu_static = coeffs[2];
  mu_dynamic = coeffs[3];

  if (k < 0.0 || xt < 0.0 || mu_static < 0.0 || mu_dynamic < 0.0)
    error->all(FLERR, "Illegal linear tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistoryStatic::calculate_forces()
{
  // Note: this is the same as the base Mindlin calculation except k isn't scaled by contact_radius
  double magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl, temp_array[3];
  int frame_update = 0;

  damp = xt * gm->damping_model->damp_prefactor;

  double Fncrit_static = gm->normal_model->Fncrit * mu_static;
  double Fncrit_dynamic = gm->normal_model->Fncrit * mu_dynamic;
  double *history = & gm->history[history_index];
  double Fscrit_static = gm->normal_model->Fncrit * mu_static;
  double Fscrit_dynamic = gm->normal_model->Fncrit * mu_dynamic;
  double Fscrit;
  int dynamic;

  dynamic = history[3];

  if (dynamic) {
	  Fscrit = gm->normal_model->Fncrit * mu_dynamic;
  } else {
	  Fscrit = gm->normal_model->Fncrit * mu_static;
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    frame_update = (fabs(rsht) * k) > (EPSILON * Fscrit);

    if (frame_update) {
      shrmag = len3(history);

      // projection
      scale3(rsht, gm->nx, temp_array);
      sub3(history, temp_array, history);

      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history, tangential force
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    scale3(gm->dt, gm->vtr, temp_array);
    add3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  scale3(-k, history, gm->fs);
  scale3(damp, gm->vtr, temp_array);
  sub3(gm->fs, temp_array, gm->fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
	shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);
      scale3(-1.0 / k, history);
      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
    if (!dynamic) history[3] = 1; // If force exceeds Fcrit_static,
  }  					          // switch to dynamic case  
  else {
	  if (dynamic) history[3] = 0; //If force drops below Fcrit_dynamic,
	} 	  	  	  	  	  	  	   //switch back to static case  
}


/* ----------------------------------------------------------------------
   Linear model with history from pair gran/hooke/history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistoryClassic::GranSubModTangentialLinearHistoryClassic(
    GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangentialLinearHistory(gm, lmp)
{
  // Whether contact radius scales normal force (as in Hertz)
  contact_radius_flag = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistoryClassic::calculate_forces()
{
  double magfs, rsht, shrmag;
  double hist_increment[3], fdamp[3];

  double *nx = gm->nx;
  double *vtr = gm->vtr;
  double *fs = gm->fs;
  double dt = gm->dt;
  double contact_radius = gm->contact_radius;
  double *history = &gm->history[history_index];
  int history_update = gm->history_update;

  damp = xt * gm->damping_model->get_damp_prefactor();
  double Fscrit = gm->normal_model->get_fncrit() * mu;

  // update history
  if (history_update) {
    scale3(dt, vtr, hist_increment);
    add3(history, hist_increment, history);
  }

  shrmag = len3(history);

  // rotate shear displacements
  if (history_update) {
    rsht = dot3(history, nx);
    scale3(rsht, nx, hist_increment);
    sub3(history, hist_increment, history);
  }

  // tangential forces = history + tangential velocity damping
  // classic model can only set contact_radius_flag through hertz
  if (contact_radius_flag)
    scale3(-k * contact_radius, history, fs);
  else
    scale3(-k, history, fs);

  // damping force, note that damp automatically has a factor
  //   of contact radius with hertz (sets viscoelastic damping)
  //   but not with hooke (sets mass_velocity damping)

  scale3(-damp, vtr, fdamp);
  add3(fs, fdamp, fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    if (shrmag != 0.0) {
      // Rescale shear force
      scale3(Fscrit / magfs, fs);

      // Set shear to elastic component of rescaled force
      //  has extra factor of kt (+ contact radius)
      sub3(fs, fdamp, history);

      // Remove extra prefactors from shear history
      if (contact_radius_flag)
        scale3(-1.0 / (k * contact_radius), history);
      else
        scale3(-1.0 / k, history);
    } else {
      zero3(fs);
    }
  }
}

/* ----------------------------------------------------------------------
   Mindlin from pair gran/hertz/history
------------------------------------------------------------------------- */

GranSubModTangentialMindlinClassic::GranSubModTangentialMindlinClassic(GranularModel *gm,
                                                                       LAMMPS *lmp) :
    GranSubModTangentialLinearHistoryClassic(gm, lmp)
{
  contact_radius_flag = 1;    // Sets gran/hertz/history behavior
}

/* ----------------------------------------------------------------------
   Mindlin model
------------------------------------------------------------------------- */

GranSubModTangentialMindlin::GranSubModTangentialMindlin(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  mindlin_force = 0;
  mindlin_rescale = 0;
  contact_radius_flag = 1;
  allow_synchronization = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k == -1) {
    if (!gm->normal_model->get_material_properties())
      error->all(FLERR,
                 "Must either specify tangential stiffness or material properties for normal model "
                 "for the Mindlin tangential style");

    double Emod = gm->normal_model->get_emod();
    double poiss = gm->normal_model->get_poiss();

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu < 0.0) error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1)
    coeffs[0] = -1;
  else
    coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::calculate_forces()
{
  double k_scaled, magfs, rsht, shrmag;
  double hist_increment[3], fdamp[3], vtr2[3];
  int frame_update = 0;

  double *nx = gm->nx;
  double *nx_unrotated = gm->nx_unrotated;
  double *vtr = gm->vtr;
  double *fs = gm->fs;
  double dt = gm->dt;
  double contact_radius = gm->contact_radius;
  double *history = &gm->history[history_index];
  int history_update = gm->history_update;

  damp = xt * gm->damping_model->get_damp_prefactor();
  double Fscrit = gm->normal_model->get_fncrit() * mu;

  k_scaled = k * gm->contact_radius;

  // on unloading, rescale the shear displacements/force
  if (mindlin_rescale)
    if (contact_radius < history[3]) scale3(contact_radius / history[3], history);

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = dot3(history, nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) rotate_rescale_vec(history, nx);

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * dt, vtr, hist_increment);
    } else {
      scale3(dt, vtr, hist_increment);
    }
    add3(history, hist_increment, history);

    if (mindlin_rescale) history[3] = contact_radius;

    if (gm->synchronized_verlet == 1) {
      // second projection to full step normal
      rsht = dot3(history, nx_unrotated);
      if (mindlin_force) {
        frame_update = fabs(rsht) > (EPSILON * Fscrit);
      } else {
        frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
      }
      if (frame_update) rotate_rescale_vec(history, nx_unrotated);
    }
  }

  // tangential forces = history + tangential velocity damping

  if (!mindlin_force) {
    scale3(-k_scaled, history, fs);
  } else {
    copy3(history, fs);
  }


  // Rotating vtr for damping term in nx direction
  if (frame_update && gm->synchronized_verlet) {
    copy3(vtr, vtr2);
    rotate_rescale_vec(vtr2, nx_unrotated);
  } else {
    copy3(vtr, vtr2);
  }
  scale3(-damp, vtr2, fdamp);
  add3(fs, fdamp, fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      // Rescale shear force
      scale3(Fscrit / magfs, fs);

      // Set shear to elastic component of rescaled force
      //  may have extra factor of k_scaled that is then removed
      sub3(fs, fdamp, history);
      if (!mindlin_force)
        scale3(-1.0 / k_scaled, history);
    } else {
      zero3(fs);
    }
  }
}

/*-----------------------------------------------------------------------
 * Mindlin with static friction coefficient
 */

GranSubModTangentialMindlinStatic::GranSubModTangentialMindlinStatic(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialMindlin(gm, lmp)
{
  num_coeffs = 4;
  size_history = 4;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu_static = coeffs[2];
  mu_dynamic = coeffs[3];

  if (k == -1) {
    if (!gm->normal_model->material_properties)
      error->all(FLERR, "Must either specify tangential stiffness or material properties for normal model for the Mindlin tangential style");

    double Emod = gm->normal_model->Emod;
    double poiss = gm->normal_model->poiss;

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu_static < 0.0 || mu_dynamic < 0.0)
    error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1) coeffs[0] = -1;
  else coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  damp = xt * gm->damping_model->damp_prefactor;

  double *history = & gm->history[history_index];
  double Fscrit_static = gm->normal_model->Fncrit * mu_static;
  double Fscrit_dynamic = gm->normal_model->Fncrit * mu_dynamic;
  double Fscrit;
  int dynamic;

  k_scaled = k * gm->contact_radius;
  dynamic = history[3];

  if (dynamic) {
	  Fscrit = gm->normal_model->Fncrit * mu_dynamic;
  } else {
	  Fscrit = gm->normal_model->Fncrit * mu_static;
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, gm->nx, temp_array);
      sub3(history, temp_array, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * gm->dt, gm->vtr, temp_array);
    } else {
      scale3(gm->dt, gm->vtr, temp_array);
    }
    add3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  scale3(-damp, gm->vtr, gm->fs);

  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(gm->fs, temp_array, gm->fs);
  } else {
    add3(gm->fs, history, gm->fs);
  }


  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);

      if (!mindlin_force) scale3(-1.0 / k_scaled, history);

      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
    if (!dynamic) history[3] = 1; // If force exceeds Fcrit_static,
    					          // switch to dynamic case
  }
  else{
	  if (dynamic) history[3] = 0; //If force drops below Fcrit_dynamic,
	  	  	  	  	  	  	  	   //switch back to static case
  }
}

/*-----------------------------------------------------------------------
 * Mindlin with static friction coefficient
 */

GranSubModTangentialMindlinStatic::GranSubModTangentialMindlinStatic(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialMindlin(gm, lmp)
{
  num_coeffs = 4;
  size_history = 4;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu_static = coeffs[2];
  mu_dynamic = coeffs[3];

  if (k == -1) {
    if (!gm->normal_model->material_properties)
      error->all(FLERR, "Must either specify tangential stiffness or material properties for normal model for the Mindlin tangential style");

    double Emod = gm->normal_model->Emod;
    double poiss = gm->normal_model->poiss;

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu_static < 0.0 || mu_dynamic < 0.0)
    error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1) coeffs[0] = -1;
  else coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  damp = xt * gm->damping_model->damp_prefactor;

  double *history = & gm->history[history_index];
  double Fscrit_static = gm->normal_model->Fncrit * mu_static;
  double Fscrit_dynamic = gm->normal_model->Fncrit * mu_dynamic;
  double Fscrit;
  int dynamic;

  k_scaled = k * gm->contact_radius;
  dynamic = history[3];

  if (dynamic) {
	Fscrit = gm->normal_model->Fncrit * mu_dynamic;
  }
  else {
	Fscrit = gm->normal_model->Fncrit * mu_static;
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, gm->nx, temp_array);
      sub3(history, temp_array, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * gm->dt, gm->vtr, temp_array);
    } else {
      scale3(gm->dt, gm->vtr, temp_array);
    }
    add3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  scale3(-damp, gm->vtr, gm->fs);

  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(gm->fs, temp_array, gm->fs);
  } else {
    add3(gm->fs, history, gm->fs);
  }


  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);

      if (!mindlin_force) scale3(-1.0 / k_scaled, history);

      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
    if (!dynamic) history[3] = 1; // If force exceeds Fcrit_static,
    					          // switch to dynamic case
  }
  else{
	  if (dynamic) history[3] = 0; //If force drops below Fcrit_dynamic,
	  	  	  	  	  	  	  	   //switch back to static case
  }
}

/*-----------------------------------------------------------------------
 * Mindlin with static friction coefficient
 */

GranSubModTangentialMindlinStatic::GranSubModTangentialMindlinStatic(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialMindlin(gm, lmp)
{
  num_coeffs = 4;
  size_history = 4;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu_static = coeffs[2];
  mu_dynamic = coeffs[3];

  if (k == -1) {
    if (!gm->normal_model->material_properties)
      error->all(FLERR, "Must either specify tangential stiffness or material properties for normal model for the Mindlin tangential style");

    double Emod = gm->normal_model->Emod;
    double poiss = gm->normal_model->poiss;

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu_static < 0.0 || mu_dynamic < 0.0)
    error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1) coeffs[0] = -1;
  else coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  damp = xt * gm->damping_model->damp_prefactor;

  double *history = & gm->history[history_index];
  double Fscrit_static = gm->normal_model->Fncrit * mu_static;
  double Fscrit_dynamic = gm->normal_model->Fncrit * mu_dynamic;
  double Fscrit;
  int dynamic;

  k_scaled = k * gm->area;
  dynamic = history[3];

  if (dynamic) {
	Fscrit = gm->normal_model->Fncrit * mu_dynamic;
  }
  else {
	Fscrit = gm->normal_model->Fncrit * mu_static;
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, gm->nx, temp_array);
      sub3(history, temp_array, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * gm->dt, gm->vtr, temp_array);
    } else {
      scale3(gm->dt, gm->vtr, temp_array);
    }
    add3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  scale3(-damp, gm->vtr, gm->fs);

  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(gm->fs, temp_array, gm->fs);
  } else {
    add3(gm->fs, history, gm->fs);
  }


  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);

      if (!mindlin_force)
        scale3(-1.0 / k_scaled, history);

      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
    if (!dynamic) history[3] = 1; // If force exceeds Fcrit_static,
    					          // switch to dynamic case
  } else {
	  if (dynamic) history[3] = 0; //If force drops below Fcrit_dynamic,
	}  	  	  	  	  	  	  	   //switch back to static case  
}

/*-----------------------------------------------------------------------
 * Mindlin with static friction coefficient
 */

GranSubModTangentialMindlinStatic::GranSubModTangentialMindlinStatic(GranularModel *gm, LAMMPS *lmp) : GranSubModTangentialMindlin(gm, lmp)
{
  num_coeffs = 4;
  size_history = 4;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu_static = coeffs[2];
  mu_dynamic = coeffs[3];

  if (k == -1) {
    if (!gm->normal_model->material_properties)
      error->all(FLERR, "Must either specify tangential stiffness or material properties for normal model for the Mindlin tangential style");

    double Emod = gm->normal_model->Emod;
    double poiss = gm->normal_model->poiss;

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu_static < 0.0 || mu_dynamic < 0.0)
    error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1) coeffs[0] = -1;
  else coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlinStatic::calculate_forces()
{
  double k_scaled, magfs, magfs_inv, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  damp = xt * gm->damping_model->damp_prefactor;

  double *history = & gm->history[history_index];
  double Fscrit_static = gm->normal_model->Fncrit * mu_static;
  double Fscrit_dynamic = gm->normal_model->Fncrit * mu_dynamic;
  double Fscrit;
  int dynamic;

  k_scaled = k * gm->contact_radius;
  dynamic = history[3];

  if (dynamic) {
	  Fscrit = gm->normal_model->Fncrit * mu_dynamic;
  } else {
	  Fscrit = gm->normal_model->Fncrit * mu_static;
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (gm->history_update) {
    rsht = dot3(history, gm->nx);
    if (mindlin_force) {
      frame_update = fabs(rsht) > (EPSILON * Fscrit);
    } else {
      frame_update = (fabs(rsht) * k_scaled) > (EPSILON * Fscrit);
    }

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, gm->nx, temp_array);
      sub3(history, temp_array, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      scale3(-k_scaled * gm->dt, gm->vtr, temp_array);
    } else {
      scale3(gm->dt, gm->vtr, temp_array);
    }
    add3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  scale3(-damp, gm->vtr, gm->fs);

  if (!mindlin_force) {
    scale3(k_scaled, history, temp_array);
    sub3(gm->fs, temp_array, gm->fs);
  } else {
    add3(gm->fs, history, gm->fs);
  }


  // rescale frictional displacements and forces if needed
  magfs = len3(gm->fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      magfs_inv = 1.0 / magfs;
      scale3(Fscrit * magfs_inv, gm->fs, history);
      scale3(damp, gm->vtr, temp_array);
      add3(history, temp_array, history);

      if (!mindlin_force)
        scale3(-1.0 / k_scaled, history);

      scale3(Fscrit * magfs_inv, gm->fs);
    } else {
      zero3(gm->fs);
    }
    if (!dynamic) history[3] = 1; // If force exceeds Fcrit_static,
    					          // switch to dynamic case
  } else {
	  if (dynamic) history[3] = 0; //If force drops below Fcrit_dynamic,
	}  	  	  	  	  	  	  	   //switch back to static case  
}

/* ----------------------------------------------------------------------
   Mindlin force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinForce::GranSubModTangentialMindlinForce(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  mindlin_force = 1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescale::GranSubModTangentialMindlinRescale(GranularModel *gm,
                                                                       LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescaleForce::GranSubModTangentialMindlinRescaleForce(GranularModel *gm,
                                                                                 LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_force = 1;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

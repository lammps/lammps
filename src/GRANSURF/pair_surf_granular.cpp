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

#include "pair_surf_granular.h"

#include "atom.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "fix_surface_local.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathExtra;

enum {NONE, LINE, TRI};

/* ---------------------------------------------------------------------- */

PairSurfGranular::PairSurfGranular(LAMMPS *lmp) : PairGranular(lmp)
{
  single_enable = 0;

  emax = 0;
  endpts = nullptr;
  cmax = 0;
  corners = nullptr;
}

/* ---------------------------------------------------------------------- */

PairSurfGranular::~PairSurfGranular()
{
  memory->destroy(endpts);
  memory->destroy(corners);
}

/* ---------------------------------------------------------------------- */

void PairSurfGranular::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,inum,jnum,itype,jtype;
  int isphere,itri,jflag,kflag,otherflag;
  double radsphere,rsq,dr[3],contact[3],ds[3],vs[3];
  double factor_couple,factor_lj,mi,mj,meff;
  double *forces, *torquesi, *torquesj, dq;
  double omega0[3] = {0.0, 0.0, 0.0};

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool touchflag = false;
  const bool history_update = update->setupflag == 0;

  class GranularModel* model;

  for (int i = 0; i < nmodels; i++)
    models_list[i]->history_update = history_update;

  ev_init(eflag,vflag);

  // if just reneighbored:
  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body
  // forward comm mass_rigid so have it for ghost lines
  // also grab current line connectivity info from FixSurfaceLocal

  if (neighbor->ago == 0) {
    if (fix_rigid) {
      int tmp;
      int *body = (int *) fix_rigid->extract("body",tmp);
      double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
      if (atom->nmax > nmax) {
        memory->destroy(mass_rigid);
        nmax = atom->nmax;
        memory->create(mass_rigid,nmax,"surf/granular:mass_rigid");
      }
      int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++)
        if (body[i] >= 0)
          mass_rigid[i] = mass_body[body[i]];
        else
          mass_rigid[i] = 0.0;
      comm->forward_comm(this);
    }

    connect2d = fsl->connect2d;
    connect3d = fsl->connect3d;
  }

  // pre-calculate current end pts of owned+ghost lines
  // only once per reneighbor if surfs not moving

  if (surfmoveflag || neighbor->ago == 0) {
    if (style == LINE) calculate_endpts();
    if (style == TRI) calculate_corners();
  }

  // loop over neighbors of my atoms
  // I is always sphere, J is always line/tri

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *line = atom->line;
  int *tri = atom->tri;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double *special_lj = force->special_lj;
  double *heatflow, *temperature;
  if (heat_flag) {
    heatflow = atom->heatflow;
    temperature = atom->temperature;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  if (use_history) {
    firsttouch = fix_history->firstflag;
    firsthistory = fix_history->firstvalue;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      jtype = type[j];
      model = models_list[types_indices[itype][jtype]];

      // Reset model and copy initial geometric data
      model->xi = x[i];
      model->xj = x[j];
      model->radi = radius[i];
      model->radj = radius[j];
      if (use_history) model->touch = touch[jj];

      touchflag = model->check_contact();

      if (!touchflag) {
        // unset non-touching neighbors
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // sanity check that neighbor list is built correctly

      if ((style == LINE) && (line[i] >= 0 || line[j] < 0))
        error->one(FLERR,"Pair surf/granular iteraction is invalid");

      if ((style == TRI) && (tri[i] >= 0 || tri[j] < 0))
        error->one(FLERR,"Pair surf/granular iteraction is invalid");

      // check for overlap of sphere and line segment/triangle
      // for line:
      //   jflag = 0 for no overlap, 1 for interior line pt, -1/-2 for end pts
      // for tri:
      //   jflag = 0 for no overlap, 1 for interior line pt,
      //     -1/-2/-3 for 3 edges, -4/-5/-6 for 3 corner pts
      // if no overlap, just continue
      // for overlap, also return:
      //   contact = nearest point on line/tri to sphere center
      //   dr = vector from contact pt to sphere center
      //   rsq = squared length of dr

      rsq = model->rsq;
      if (style == LINE) jflag = overlap_sphere_line(i,j,contact,dr,rsq);
      if (style == TRI) jflag = overlap_sphere_tri(i,j,contact,dr,rsq);

      if (!jflag) {
        // unset non-touching neighbors
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // if any history is needed
      if (use_history) touch[jj] = 1;

      // if contact = line end pt or tri edge/corner:
      // check overlap status of adjacent line to the end pt or tri(s) to edge/corner
      // otherflag = 0/1 for this/other line/tri performs calculation

      if (jflag < 0) {
        if (style == LINE) {
          otherflag = endpt_neigh_check(i,j,jflag);
          if (otherflag) continue;
        } else if (style == TRI) {
          if (jflag >= -3) {
            otherflag = edge_neigh_check(i,j,jflag);
            if (otherflag) continue;
          } else {
            otherflag = corner_neigh_check(i,j,jflag);
            if (otherflag) continue;
          }
        }
      }

      // calculate new data
      // ds = vector from line/tri center to contact pt
      // vs = velocity of contact pt on line/tri, translation + rotation
      // omega for tri was set from angmom by calculate_corners()

      model->rsq = rsq;

      ds[0] = contact[0] - x[j][0];
      ds[1] = contact[1] - x[j][1];
      ds[2] = contact[2] - x[j][2];

      vs[0] = v[j][0] + (omega[j][1] * ds[2] - omega[j][2] * ds[1]);
      vs[1] = v[j][1] + (omega[j][2] * ds[0] - omega[j][0] * ds[2]);
      vs[2] = v[j][2] + (omega[j][0] * ds[1] - omega[j][1] * ds[0]);

      model->dx[0] = dr[0];
      model->dx[1] = dr[1];
      model->dx[2] = dr[2];

      // NOTE: add logic to persist history history if contact has changed

      // NOTE: add logic to check for coupled contacts and weight them

      factor_couple = 1.0;

      // meff = effective mass of sphere and line/tri
      // if I or J is part of rigid body,use body mass
      // if line/tri is not part of rigidbody assume infinite mass

      meff = rmass[i];
      if (fix_rigid) {
        if (mass_rigid[i] > 0.0) meff = mass_rigid[i];
        if (mass_rigid[j] > 0.0) {
          mj = mass_rigid[j];
          meff = meff * mj / (meff + mj);
        }
      }

      // Copy additional information and prepare force calculations
      model->meff = meff;
      model->vi = v[i];
      model->vj = vs;
      model->omegai = omega[i];
      model->omegaj = omega0;

      if (use_history) {
        history = &allhistory[size_history * jj];
        model->history = history;
      }

      if (heat_flag) {
        model->Ti = temperature[i];
        model->Tj = temperature[j];
      }

      model->calculate_forces();

      // need to add support coupled contacts
      // is this just multiplying forces (+torques?) by factor_couple?

      forces = model->forces;
      torquesi = model->torquesi;
      torquesj = model->torquesj;

      // apply forces & torques
      scale3(factor_lj, forces);
      add3(f[i], forces, f[i]);

      scale3(factor_lj, torquesi);
      add3(torque[i], torquesi, torque[i]);

      if (force->newton_pair || j < nlocal) {
        sub3(f[j], forces, f[j]);
        scale3(factor_lj, torquesj);
        add3(torque[j], torquesj, torque[j]);
      }

      if (heat_flag) {
        dq = model->dq;
        heatflow[i] += dq;
        if (force->newton_pair || j < nlocal) heatflow[j] -= dq;
      }

      if (evflag) {
        ev_tally_xyz(i,j,nlocal,force->newton_pair,
          0.0,0.0,forces[0],forces[1],forces[2],model->dx[0],model->dx[1],model->dx[2]);
      }
    }
  }

  // NOTE: should there be virial contributions from boundary tris?

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSurfGranular::init_style()
{
  // error and warning checks

  avecline = (AtomVecLine *) atom->style_match("line");
  avectri = (AtomVecTri *) atom->style_match("tri");

  if (avecline) style = LINE;
  if (avectri) style = TRI;

  if (style == NONE)
    error->all(FLERR,"Pair surf/granular requires atom style line or tri");

  if (!atom->radius_flag || !atom->rmass_flag || !atom->omega_flag)
    error->all(FLERR, "Pair surf/granular requires atom attributes radius, rmass, omega");
  if (!force->newton_pair)
    error->all(FLERR,"Pair style surf/granular requires newton pair on");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair surf/granular requires ghost atoms store velocity");

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR,"Heat conduction in pair surf/granularular requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR,"Heat conduction in pair surf/granularular requires atom style with heatflow property");
  }

  // allocate history and initialize models
  class GranularModel* model;
  int size_max[NSUBMODELS] = {0};
  for (int n = 0; n < nmodels; n++) {
    model = models_list[n];
    model->contact_type = SURFACE;

    if (model->beyond_contact) {
      beyond_contact = 1;
      use_history = 1; // Need to track if in contact
    }
    if (model->size_history != 0) use_history = 1;

    for (int i = 0; i < NSUBMODELS; i++)
      if (model->sub_models[i]->size_history > size_max[i])
        size_max[i] = model->sub_models[i]->size_history;

    if (model->nondefault_history_transfer) nondefault_history_transfer = 1;
  }

  size_history = 0;
  if (use_history) {
    for (int i = 0; i < NSUBMODELS; i++) size_history += size_max[i];

    // Ensure size history is at least 1 to avoid errors in fix neigh/history
    // This could occur if normal model is beyond_contact but no other entries are required
    // E.g. JKR + linear_nohistory
    size_history = MAX(size_history, 1);
  }

  for (int n = 0; n < nmodels; n++) {
    model = models_list[n];
    int next_index = 0;
    for (int i = 0; i < NSUBMODELS; i++) {
      model->sub_models[i]->history_index = next_index;
      next_index += size_max[i];
    }

    if (model->beyond_contact)
      error->all(FLERR, "Beyond contact models not currenty supported");
  }

  // need a granular neighbor list

  if (use_history)
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_ONESIDED |
                          NeighConst::REQ_HISTORY);
  else
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_ONESIDED);

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (use_history && fix_history == nullptr) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->replace_fix(id_dummy, fmt::format("{} all NEIGH_HISTORY {} onesided", id_history, size_history),1));
    fix_history->pair = this;
  } else if (use_history) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id(id_history));
    if (!fix_history) error->all(FLERR,"Could not find pair fix neigh history ID");
  }

  // set ptr to FixSurfaceLocal for surf connectivity info

  fsl = nullptr;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style,"surface/local") == 0) {
      if (fsl)
        error->all(FLERR,"Pair surf/granular requires single fix surface/local");
      fsl = (FixSurfaceLocal *) modify->fix[m];
    }
  }
  if (!fsl) error->all(FLERR,"Pair surf/granular requires a fix surface/local");

  // surfmoveflag = 1 if surfs may move at every step
  // yes if fix move exists and its group includes lines
  // NOTE: are there other conditions, like fix deform or fix npt?

  surfmoveflag = 0;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style,"move") == 0) {
      int groupbit = modify->fix[m]->groupbit;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;
      int flag = 0;
      if (style == LINE) {
        int *line = atom->line;
        for (int i = 0; i < nlocal; i++) {
          if (line[i] < 0) continue;
          if (mask[i] & groupbit) flag = 1;
        }
      } else {
        int *tri = atom->tri;
        for (int i = 0; i < nlocal; i++) {
          if (tri[i] < 0) continue;
          if (mask[i] & groupbit) flag = 1;
        }
      }
      int any;
      MPI_Allreduce(&flag,&any,1,MPI_INT,MPI_SUM,world);
      if (any) surfmoveflag = 1;
    }
  }

  // check for FixFreeze and set freeze_group_bit

  auto fixlist = modify->get_fix_by_style("^freeze");
  if (fixlist.size() == 0)
    freeze_group_bit = 0;
  else if (fixlist.size() > 1)
    error->all(FLERR, "Only one fix freeze command at a time allowed");
  else
    freeze_group_bit = fixlist.front()->groupbit;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = nullptr;
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->rigid_flag) {
      if (fix_rigid)
        error->all(FLERR, "Only one fix rigid command at a time allowed");
      else
        fix_rigid = ifix;
    }
  }

  // check for FixPour and FixDeposit so can extract particle radii

  auto pours = modify->get_fix_by_style("^pour");
  auto deps = modify->get_fix_by_style("^deposit");

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic
  // lines/tris cannot be frozen

  int itype;
  for (int i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    for (auto &ipour : pours) {
      itype = i;
      double maxrad = *((double *) ipour->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
    for (auto &idep : deps) {
      itype = i;
      double maxrad = *((double *) idep->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
  }

  double *radius = atom->radius;
  int *line = atom->line;
  int *tri = atom->tri;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if ((style == LINE) && line[i] >= 0)
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);
    else if ((style == TRI) && tri[i] >= 0)
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);
    else {
      if (mask[i] & freeze_group_bit)
        onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
      else
        onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);
    }
  }

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairSurfGranular::memory_usage()
{
  double bytes = nmax * sizeof(double);
  if (style == LINE) bytes = emax * 4 * sizeof(double);        // endpts array for line particles
  if (style == TRI) bytes = emax * 12 * sizeof(double);        // corners array for tri particles
  return bytes;
}

/* ----------------------------------------------------------------------
   compute current end points of owned and ghost lines
   nothing computed for particles that are not lines
------------------------------------------------------------------------- */

void PairSurfGranular::calculate_endpts()
{
  int i,m;
  double length,theta,dx,dy;
  double *endpt;

  // realloc endpts array if necssary

  if (fsl->nmax_connect > emax) {
    memory->destroy(endpts);
    emax = fsl->nmax_connect;
    memory->create(endpts,emax,4,"surf/granular:endpts");
  }

  AtomVecLine::Bonus *bonus = avecline->bonus;
  double **x = atom->x;
  int *line = atom->line;
  int n = atom->nlocal + atom->nghost;

  for (i = 0; i < n; i++) {
    if (line[i] < 0) continue;
    m = line[i];
    length = bonus[m].length;
    theta = bonus[m].theta;
    dx = 0.5*length*cos(theta);
    dy = 0.5*length*sin(theta);
    endpt = endpts[m];
    endpt[0] = x[i][0] - dx;
    endpt[1] = x[i][1] - dy;
    endpt[2] = x[i][0] + dx;
    endpt[3] = x[i][1] + dy;
  }
}

/* ----------------------------------------------------------------------
   compute nearest point between sphere I and line segment J
   return 0 if no contact, 1 if pt is interior to line segment,
     -1/-2 if pt = line end point 1/2
   if contact, return:
     pt = point on line segment
     r = vector from pt to sphere center
     rsq = squared length of r
   based on geometry.cpp::distsq_point_line() in SPARTA
------------------------------------------------------------------------- */

int PairSurfGranular::overlap_sphere_line(int i, int j, double *pt,
                                           double *r, double &rsq)
{
  double p1[3],p2[3];
  double a[3],b[3];

  // P1,P2 = end points of line segment

  double *endpt = endpts[atom->line[j]];

  p1[0] = endpt[0];
  p1[1] = endpt[1];
  p1[2] = 0.0;
  p2[0] = endpt[2];
  p2[1] = endpt[3];
  p2[2] = 0.0;

  // A = vector from P1 to Xsphere
  // B = vector from P1 to P2

  double *xsphere = atom->x[i];
  MathExtra::sub3(xsphere,p1,a);
  MathExtra::sub3(p2,p1,b);

  // alpha = fraction of distance from P1 to P2 that P is located at
  // P = projected point on infinite line that is nearest to Xsphere center
  // alpha can be any value

  double alpha = MathExtra::dot3(a,b) / MathExtra::lensq3(b);

  // pt = point on line segment that is nearest to Xsphere center
  // if alpha <= 0.0, pt = P1, ptflag = -1
  // if alpha >= 1.0, pt = P2, ptflag = -2
  // else pt = P1 + alpha*(P2-P1), ptflag = 1

  int ptflag;
  if (alpha <= 0.0) {
    ptflag = -1;
    pt[0] = p1[0];
    pt[1] = p1[1];
    pt[2] = p1[2];
  } else if (alpha >= 1.0) {
    ptflag = -2;
    pt[0] = p2[0];
    pt[1] = p2[1];
    pt[2] = p2[2];
  } else {
    ptflag = 1;
    pt[0] = p1[0] + alpha*b[0];
    pt[1] = p1[1] + alpha*b[1];
    pt[2] = p1[2] + alpha*b[2];
  }

  // R = vector from nearest pt on line to Xsphere center
  // return ptflag if len(R) < sphere radius
  // else no contact, return 0

  double radsq = atom->radius[i] * atom->radius[i];
  MathExtra::sub3(xsphere,pt,r);
  rsq = MathExtra::lensq3(r);
  if (rsq < radsq) return ptflag;
  return 0;
}

/* ----------------------------------------------------------------------
   check overlap status of sphere I with line J versus any neighbor lines K
   I overlaps J at jflag = -1,-2 for two end points
   return 0 if this line J performs computation
   return 1 if some other line K performs computation
------------------------------------------------------------------------- */

int PairSurfGranular::endpt_neigh_check(int i, int j, int jflag)
{
  // ncheck = # of neighbor lines to check
  // neighs = indices of neighbor lines (including self)

  int ncheck;
  int *neighs;

  int jc = atom->line[j];
  if (jflag == -1) {
    if (connect2d[jc].np1 == 1) return 0;
    ncheck = connect2d[jc].np1;
    neighs = connect2d[jc].neigh_p1;
  } else if (jflag == -2) {
    if (connect2d[jc].np2 == 1) return 0;
    ncheck = connect2d[jc].np2;
    neighs = connect2d[jc].neigh_p2;
  }

  // check overlap with each neighbor line
  // if any line has interior overlap, another line computes
  // if all lines have endpt overlap, line with lowest index computes
  // kflag = overlap status with neighbor line
  // kflag = 1, interior overlap
  // kflag = 0, no overlap, should not be possible
  // kflag < 0, overlap at endpt

  tagint *tag = atom->tag;

  int k,kflag;
  double rsq;
  double dr[3],contact[3];

  int linemin = tag[j];

  for (int m = 0; m < ncheck; m++) {
    if (neighs[m] == tag[j]) continue;     // skip self line
    k = atom->map(neighs[m]);
    if (k < 0) error->one(FLERR,"Pair surf/granular neighbor line is missing");
    kflag = overlap_sphere_line(i,k,contact,dr,rsq);
    if (kflag > 0) return 1;
    if (kflag == 0) error->one(FLERR,"Fix surface/global neighbor line overlap is invalid");
    linemin = MIN(linemin,tag[k]);
  }

  if (tag[j] == linemin) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   compute current corner points and current norm of N triangles
   also compute omega from angmom
   N = nlocal or nlocal+nghost atoms
   nothing computed for particles that are not tris
------------------------------------------------------------------------- */

void PairSurfGranular::calculate_corners()
{
  int i,m;
  double ex[3],ey[3],ez[3],p[3][3];
  double *corner;

  // realloc corners array if necssary

  if (fsl->nmax_connect > cmax) {
    memory->destroy(corners);
    cmax = fsl->nmax_connect;
    memory->create(corners,cmax,12,"surf/granular:corners");
  }

  AtomVecTri::Bonus *bonus = avectri->bonus;
  double **x = atom->x;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  int *tri = atom->tri;
  int n = atom->nlocal + atom->nghost;

  for (int i = 0; i < n; i++) {
    if (tri[i] < 0) continue;
    m = tri[i];
    corner = corners[m];
    MathExtra::quat_to_mat(bonus[m].quat,p);
    MathExtra::matvec(p,bonus[m].c1,&corner[0]);
    MathExtra::add3(x[i],&corner[0],&corner[0]);
    MathExtra::matvec(p,bonus[m].c2,&corner[3]);
    MathExtra::add3(x[i],&corner[3],&corner[3]);
    MathExtra::matvec(p,bonus[m].c3,&corner[6]);
    MathExtra::add3(x[i],&corner[6],&corner[6]);
    corners2norm(corner,&corner[9]);

    // omega from angmom of tri particles

    if (angmom[i][0] == 0.0 && angmom[i][1] == 0.0 && angmom[i][2] == 0.0) {
      omega[i][0] = omega[i][1] = omega[i][2] = 0.0;
      continue;
    }
    MathExtra::q_to_exyz(bonus[m].quat,ex,ey,ez);
    MathExtra::angmom_to_omega(angmom[i],ex,ey,ez,
                               bonus[m].inertia,omega[i]);
  }
}

/* ----------------------------------------------------------------------
   compute norm of a triangle based on its 3 corner pts
------------------------------------------------------------------------- */

void PairSurfGranular::corners2norm(double *corners, double *norm)
{
  double p12[3],p13[3];

  MathExtra::sub3(&corners[0],&corners[3],p12);
  MathExtra::sub3(&corners[0],&corners[6],p13);
  MathExtra::cross3(p12,p13,norm);
  MathExtra::norm3(norm);
}

/* ----------------------------------------------------------------------
   compute nearest point between sphere I and triangle J
   return 0 if no contact, 1 if pt is interior to triangle,
     -1/-2/-3 if pt on tri edges, -4/-5/-6 if pt = tri corners 1/2/3
   if contact, return:
     pt = point on triangle
     r = vector from pt to sphere center
     rsq = squared length of r
   based on geometry.cpp::distsq_point_tri() in SPARTA
------------------------------------------------------------------------- */

int PairSurfGranular::overlap_sphere_tri(int i, int j, double *pt,
                                                double *r, double &rsq)
{
  int e12flag,e23flag,e31flag,o12flag,o23flag,o31flag;
  int esum,osum,lineflag;
  double dot;
  double *p1,*p2,*p3,*norm;
  double a[3],point[3],edge[3],pvec[3],xproduct[3];

  // P1,P2,P3 = corner points of triangle
  // norm = current norm of triangle

  double *corner = corners[atom->tri[j]];

  p1 = &corner[0];
  p2 = &corner[3];
  p3 = &corner[6];
  norm = &corner[9];

  // A = vector from P1 to Xsphere

  double *xsphere = atom->x[i];
  MathExtra::sub3(xsphere,p1,a);

  // pt = projected point on infinite triangle plane

  double alpha = MathExtra::dot3(a,norm);
  pt[0] = xsphere[0] - alpha*norm[0];
  pt[1] = xsphere[1] - alpha*norm[1];
  pt[2] = xsphere[2] - alpha*norm[2];

  // test if projected point is inside triangle
  // inside = interior + boundary of tri
  // edge = edge vector of triangle
  // pvec = vector from triangle vertex to projected point
  // xproduct = cross product of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   projected point is outside tri
  // NOTE: worry about round-off for pt being on edge or corner?

  int inside = 1;
  e12flag = e23flag = e31flag = 0;
  o12flag = o23flag = o31flag = 0;

  MathExtra::sub3(p2,p1,edge);
  MathExtra::sub3(pt,p1,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  dot = MathExtra::dot3(xproduct,norm);
  if (dot <= 0.0) {
    o12flag = 1;
    if (dot == 0.0) e12flag = 1;
    else inside = 0;
  }

  MathExtra::sub3(p3,p2,edge);
  MathExtra::sub3(pt,p2,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  dot = MathExtra::dot3(xproduct,norm);
  if (dot <= 0.0) {
    o23flag = 1;
    if (dot == 0.0) e23flag = 2;
    else inside = 0;
  }

  MathExtra::sub3(p1,p3,edge);
  MathExtra::sub3(pt,p3,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  dot = MathExtra::dot3(xproduct,norm);
  if (dot <= 0.0) {
    o31flag = 1;
    if (dot == 0.0) e31flag = 3;
    else inside = 0;
  }

  // projected point is inside tri = interior or boundary
  // set ptflag = 1 for interior
  // set ptflag = -1,-2,-3 for 3 edges E12,E23,E31
  // set ptflag = -4,-5,-6 for 3 corner pts P1,P2,P3

  int flag = 0;
  if (inside) {
    flag = 1;
    esum = e12flag + e23flag + e31flag;
    if (esum) {
      if (esum == 1) {
        if (e12flag) flag = -1;
        else if (e23flag) flag = -2;
        else flag = -3;
      } else {
        if (!e12flag) flag = -6;
        else if (!e23flag) flag = -4;
        else flag = -5;
      }
    }

  // projected point is outside tri
  // reset pt to nearest point to tri center
  // set ptflag = -1,-2,-3 if pt on edges
  // set ptflag = -4,-5,-6 if pt = corner pts

  } else {
    osum = o12flag + o23flag + o31flag;
    if (osum == 1) {
      if (o12flag) {
        lineflag = nearest_point_line(xsphere,p1,p2,pt);
        if (lineflag == 1) flag = -1;
        else if (lineflag == -1) flag = -4;
        else flag = -5;
      } else if (o23flag) {
        lineflag = nearest_point_line(xsphere,p2,p3,pt);
        if (lineflag == 1) flag = -2;
        else if (lineflag == -1) flag = -5;
        else flag = -6;
      } else {
        lineflag = nearest_point_line(xsphere,p3,p1,pt);
        if (lineflag == 1) flag = -3;
        else if (lineflag == -1) flag = -6;
        else flag = -4;
      }
    } else {
      if (!o12flag) {
        flag = -6;
        pt[0] = p3[0];
        pt[1] = p3[1];
        pt[2] = p3[2];
      } else if (!o23flag) {
        flag = -4;
        pt[0] = p1[0];
        pt[1] = p1[1];
        pt[2] = p1[2];
      } else {
        flag = -5;
        pt[0] = p2[0];
        pt[1] = p2[1];
        pt[2] = p2[2];
      }
    }
  }

  // test if point is exactly a corner pt
  // if so, reset ptwhich to corner pt

  /*
  if (pt[0] == p1[0] && pt[1] == p1[1] && pt[2] == p1[2]) flag == -4;
  else if (pt[0] == p2[0] && pt[1] == p2[1] && pt[2] == p2[2]) flag == -5;
  else if (pt[0] == p3[0] && pt[1] == p3[1] && pt[2] == p3[2]) flag == -6;
  */

  // R = vector from nearest pt on line to Xsphere center
  // return flag if len(R) < sphere radius
  // else no contact, return 0

  double radsq = atom->radius[i] * atom->radius[i];
  MathExtra::sub3(xsphere,pt,r);
  rsq = MathExtra::lensq3(r);
  if (rsq < radsq) return flag;
  return 0;
}

/* ----------------------------------------------------------------------
   compute nearest point between point X and line segment P1 to P2
   return pt = nearest point within line segment
   return 1 if pt is interior to line segment
   return -1/-2 if pt = line segment end point 1/2
   based on geometry.cpp::distsq_point_line() in SPARTA
------------------------------------------------------------------------- */

int PairSurfGranular::nearest_point_line(double *x,
                                                double *p1, double *p2,
                                                double *pt)
{
  double a[3],b[3];

  // A = vector from P1 to X
  // B = vector from P1 to P2

  MathExtra::sub3(x,p1,a);
  MathExtra::sub3(p2,p1,b);

  // alpha = fraction of distance from P1 to P2 that P is located at
  // P = projected point on infinite line that is nearest to X
  // alpha can be any value

  double alpha = MathExtra::dot3(a,b) / MathExtra::lensq3(b);

  // pt = point on line segment that is nearest to X
  // if alpha <= 0.0, pt = P1, ptflag = -1
  // if alpha >= 1.0, pt = P2, ptflag = -2
  // else pt = P1 + alpha*(P2-P1), ptflag = 1

  int ptflag;
  if (alpha <= 0.0) {
    ptflag = -1;
    pt[0] = p1[0];
    pt[1] = p1[1];
    pt[2] = p1[2];
  } else if (alpha >= 1.0) {
    ptflag = -2;
    pt[0] = p2[0];
    pt[1] = p2[1];
    pt[2] = p2[2];
  } else {
    ptflag = 1;
    pt[0] = p1[0] + alpha*b[0];
    pt[1] = p1[1] + alpha*b[1];
    pt[2] = p1[2] + alpha*b[2];
  }

  return ptflag;
}

/* ----------------------------------------------------------------------
   check overlap status of sphere I with tri J versus any neighbor tris K
   I overlaps J at jflag = -1,-2,-3 for three edges
   return 0 if this tri J performs computation
   return 1 if other tri K performs computation
------------------------------------------------------------------------- */

int PairSurfGranular::edge_neigh_check(int i, int j, int jflag)
{
  // ncheck = # of neighbor tris to check
  // neighs = indices of neighbor tris (including self)

  int ncheck;
  int *neighs;

  int jc = atom->tri[j];
  if (jflag == -1) {
    if (connect3d[jc].ne1 == 1) return 0;
    ncheck = connect3d[jc].ne1;
    neighs = connect3d[jc].neigh_e1;
  } else if (jflag == -2) {
    if (connect3d[jc].ne2 == 1) return 0;
    ncheck = connect3d[jc].ne2;
    neighs = connect3d[jc].neigh_e2;
  } else if (jflag == -3) {
    if (connect3d[jc].ne3 == 1) return 0;
    ncheck = connect3d[jc].ne3;
    neighs = connect3d[jc].neigh_e3;
  }

  // check overlap with each neighbor tri
  // if any tri has interior overlap, another tri computes
  // if all tris have edge overlap, tri with lowest ID computes
  // kflag = overlap status with neighbor tri
  // kflag = 1, interior overlap
  // kflag = 0, no overlap, should not be possible
  // kflag < 0, overlap at edge (overlap at corner pt should not be possible)

  tagint *tag = atom->tag;

  int k,kflag;
  double rsq;
  double dr[3],contact[3];

  tagint trimin = tag[j];

  for (int m = 0; m < ncheck; m++) {
    if (neighs[m] == tag[j]) continue;     // skip self tri
    k = atom->map(neighs[m]);
    if (k < 0) error->one(FLERR,"Pair surf/granular neighbor tri is missing");
    kflag = overlap_sphere_tri(i,k,contact,dr,rsq);
    if (kflag > 0) return 1;
    if (kflag == 0) error->one(FLERR,"Pair surf/granular neighbor tri overlap is invalid");
    trimin = MIN(trimin,k);
  }

  if (tag[j] == trimin) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check overlap status of sphere I with tri J versus any neighbor tris K
   I overlaps J at jflag = -4,-5,-6 for three corners
   return 0 if this tri J performs computation
   return 1 if some other tri K performs computation
------------------------------------------------------------------------- */

int PairSurfGranular::corner_neigh_check(int i, int j, int jflag)
{
  // ncheck = # of neighbor tris to check
  // neighs = indices of neighbor tris (including self)

  int ncheck;
  int *neighs;

  int jc = atom->tri[j];
  if (jflag == -4) {
    if (connect3d[jc].nc1 == 1) return 0;
    ncheck = connect3d[jc].nc1;
    neighs = connect3d[jc].neigh_c1;
  } else if (jflag == -5) {
    if (connect3d[jc].nc2 == 1) return 0;
    ncheck = connect3d[jc].nc2;
    neighs = connect3d[jc].neigh_c2;
  } else if (jflag == -6) {
    if (connect3d[jc].nc3 == 1) return 0;
    ncheck = connect3d[jc].nc3;
    neighs = connect3d[jc].neigh_c3;
  }

  // check overlap with each neighbor tri
  // if any tri has interior or edge overlap, another tri computes
  // if all tris have corner pt overlap, tri with lowest ID computes
  // kflag = overlap status with neighbor tri
  // kflag = 1, interior overlap
  // kflag = 0, no overlap, should not be possible
  // kflag = -1/-2/-3, overlap at edge
  // kflag = -4/-5/-6, overlap at corner pt

  tagint *tag = atom->tag;

  int k,kflag;
  double rsq;
  double dr[3],contact[3];

  tagint trimin = tag[j];

  for (int m = 0; m < ncheck; m++) {
    if (neighs[m] == tag[j]) continue;     // skip self tri
    k = atom->map(neighs[m]);
    if (k < 0) error->one(FLERR,"Pair surf/granular neighbor tri is missing");
    kflag = overlap_sphere_tri(i,k,contact,dr,rsq);
    if (kflag > 0) return 1;
    if (kflag == 0) error->one(FLERR,"Pair surf/granular neighbor tri overlap is invalid");
    if (kflag >= -3) return 1;
    trimin = MIN(trimin,tag[k]);
  }

  if (tag[j] == trimin) return 0;
  return 1;
}

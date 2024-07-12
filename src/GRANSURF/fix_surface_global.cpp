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

#include "fix_surface_global.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "input.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "stl_reader.h"
#include "update.h"
#include "variable.h"

#include <map>
#include <tuple>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace Granular_NS;
using namespace MathConst;
using namespace MathExtra;

enum{SPHERE,LINE,TRI};           // also in DumpImage
enum{NONE,LINEAR,WIGGLE,ROTATE,VARIABLE};

#define DELTA 128

/* ---------------------------------------------------------------------- */

FixSurfaceGlobal::FixSurfaceGlobal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), tstr(nullptr)
{
  if (narg < 11) error->all(FLERR,"Illegal fix surface/global command");

  if (!atom->omega_flag) error->all(FLERR,"Fix surface/global requires atom attribute omega");
  if (!atom->radius_flag) error->all(FLERR,"Fix surface/global requires atom attribute radius");

  // set interaction style
  // disable bonded/history option for now
  
  model = new GranularModel(lmp);
  model->contact_type = SURFGLOBAL;

  heat_flag = 0;
  int classic_flag = 1;
  if (strcmp(arg[4], "granular") == 0)  classic_flag = 0;

  // wall/particle coefficients

  int iarg;
  if (classic_flag) {
    iarg = model->define_classic_model(arg, 4, narg);

    if (iarg < narg) {
      if (strcmp(arg[iarg],"limit_damping") == 0) {
        model->limit_damping = 1;
        iarg += 1;
      }
    }

  } else {
    iarg = 5;
    iarg = model->add_sub_model(arg, iarg, narg, NORMAL);

    while (iarg < narg) {
      if (strcmp(arg[iarg], "damping") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, DAMPING);
      } else if (strcmp(arg[iarg], "tangential") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, TANGENTIAL);
      } else if (strcmp(arg[iarg], "rolling") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, ROLLING);
      } else if (strcmp(arg[iarg], "twisting") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, TWISTING);
      } else if (strcmp(arg[iarg], "heat") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, HEAT);
        heat_flag = 1;
      } else if (strcmp(arg[iarg],"limit_damping") == 0) {
        model->limit_damping = 1;
        iarg += 1;
      } else {
        break;
      }
    }
  }

  // define default damping sub model if unspecified, takes no args
  
  if (!model->damping_model) model->construct_sub_model("viscoelastic", DAMPING);
  model->init();

  size_history = model->size_history;
  if (model->beyond_contact) size_history += 1; //Need to track if particle is touching
  if (size_history == 0) use_history = 0;
  else use_history = 1;

  // optional args

  int scaleflag = 0;
  int Twall_defined = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix surface/global command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix surface/global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix surface/global command");
      if (utils::strmatch(arg[iarg+1], "^v_")) {
        tstr = utils::strdup(arg[iarg+1] + 2);
      } else {
        Twall = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      }
      Twall_defined = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix surface/global command");
  }

  if (heat_flag != Twall_defined)
    error->all(FLERR, "Must define wall temperature with heat model");

  // initializations

  dimension = domain->dimension;

  mstyle = NONE;
  points_lastneigh = nullptr;
  points_original = nullptr;
  xsurf_original = nullptr;

  connect2d = nullptr;
  connect3d = nullptr;
  plist = nullptr;
  elist = nullptr;
  clist = nullptr;

  nmax = 0;
  mass_rigid = nullptr;

  fix_rigid = nullptr;
  fix_history = nullptr;

  list = new NeighList(lmp);
  if (use_history) {
    listhistory = new NeighList(lmp);
    zeroes = new double[size_history];
    for (int i = 0; i < size_history; i++) zeroes[i] = 0.0;
  } else {
    listhistory = nullptr;
    zeroes = nullptr;
  }

  imax = 0;
  imflag = nullptr;
  imdata = nullptr;

  firsttime = 1;

  // setup scale factors for possible fix modify move settings

  if (scaleflag) {
    double xscale = domain->lattice->xlattice;
    double yscale = domain->lattice->ylattice;
    double zscale = domain->lattice->zlattice;
  } else xscale = yscale = zscale = 1.0;

  // define points/lines/tris and their connectivity
  // done via molecule template ID or STL file
  // check if arg = valid molecute template ID, else treat as STL file
  
  int imol = atom->find_molecule(arg[3]);
  if (imol >= 0) extract_from_molecules(arg[3]);
  else extract_from_stlfile(arg[3]);

  if (dimension == 2) connectivity2d_global();
  else connectivity3d_global();

  nsurf = nlines;
  if (dimension == 3) nsurf = ntris;

  surface_attributes();
}

/* ---------------------------------------------------------------------- */

FixSurfaceGlobal::~FixSurfaceGlobal()
{
  memory->sfree(points);
  memory->sfree(lines);
  memory->sfree(tris);

  memory->destroy(points_lastneigh);
  memory->destroy(points_original);
  memory->destroy(xsurf_original);

  memory->sfree(connect2d);
  memory->sfree(connect3d);
  memory->destroy(plist);
  memory->destroy(elist);
  memory->destroy(clist);

  memory->destroy(xsurf);
  memory->destroy(vsurf);
  memory->destroy(omegasurf);
  memory->destroy(radsurf);

  memory->destroy(mass_rigid);

  delete list;
  delete listhistory;
  delete [] zeroes;
  delete[] tstr;

  if (use_history)
    modify->delete_fix("NEIGH_HISTORY_SURFACE_GLOBAL_" + std::to_string(instance_me));

  memory->destroy(imflag);
  memory->destroy(imdata);
}

/* ----------------------------------------------------------------------
   create Fix needed for storing shear history if needed
   must be done in post_constructor()
------------------------------------------------------------------------- */

void FixSurfaceGlobal::post_constructor()
{
  if (use_history) {
    auto cmd = fmt::format("NEIGH_HISTORY_SURFACE_GLOBAL_" + std::to_string(instance_me) + " all NEIGH_HISTORY {}",size_history);
    fix_history = dynamic_cast<FixNeighHistory *>(modify->add_fix(cmd));
  } else
    fix_history = nullptr;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceGlobal::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::init()
{
  dt = update->dt;
  triggersq = 0.25 * neighbor->skin * neighbor->skin;

  // check for compatible heat conduction atom style

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR, "Heat conduction in fix surface/global requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR, "Heat conduction in fix surface/global requires atom style with heatflow property");
  }

  // define history indices

  int next_index = 0;
  if (model->beyond_contact) //next_index = 1;
    error->all(FLERR, "Beyond contact models not currenty supported");

  for (int i = 0; i < NSUBMODELS; i++) {
    model->sub_models[i]->history_index = next_index;
    next_index += model->sub_models[i]->size_history;
  }

  // one-time setup and allocation of neighbor list
  // wait until now, so neighbor settings have been made

  if (firsttime) {
    firsttime = 0;
    int pgsize = neighbor->pgsize;
    int oneatom = neighbor->oneatom;
    list->setup_pages(pgsize,oneatom);
    list->grow(atom->nmax,atom->nmax);

    if (use_history) {
      listhistory->setup_pages(pgsize,oneatom);
      listhistory->grow(atom->nmax,atom->nmax);
    }
  }

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) error->all(FLERR, "Variable {} for fix surface/global does not exist", tstr);
    if (! input->variable->equalstyle(tvar))
      error->all(FLERR, "Variable {} for fix surface/global must be an equal style variable", tstr);
  }
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   move surfaces via fix_modify setting
   similar to fix move operations
------------------------------------------------------------------------- */

void FixSurfaceGlobal::initial_integrate(int vflag)
{
  double ddotr;
  double a[3],b[3],c[3],d[3],disp[3],p12[3],p13[3];
  double *pt,*p1,*p2,*p3;

  double delta = (update->ntimestep - time_origin) * dt;

  // for rotate by right-hand rule around omega:
  // P = point = vector = point of rotation
  // R = vector = axis of rotation
  // w = omega of rotation (from period)
  // X0 = xoriginal = initial coord of atom
  // R0 = runit = unit vector for R
  // D = X0 - P = vector from P to X0
  // C = (D dot R0) R0 = projection of atom coord onto R line
  // A = D - C = vector from R line to X0
  // B = R0 cross A = vector perp to A in plane of rotation
  // A,B define plane of circular rotation around R line
  // X = P + C + A cos(w*dt) + B sin(w*dt)
  // V = w R0 cross (A cos(w*dt) + B sin(w*dt))

  if (mstyle == ROTATE) {
    double arg = omega_rotate * delta;
    double cosine = cos(arg);
    double sine = sin(arg);

    for (int i = 0; i < npoints; i++) {
      d[0] = points_original[i][0] - rpoint[0];
      d[1] = points_original[i][1] - rpoint[1];
      d[2] = points_original[i][2] - rpoint[2];

      ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
      c[0] = ddotr*runit[0];
      c[1] = ddotr*runit[1];
      c[2] = ddotr*runit[2];
      a[0] = d[0] - c[0];
      a[1] = d[1] - c[1];
      a[2] = d[2] - c[2];
      b[0] = runit[1]*a[2] - runit[2]*a[1];
      b[1] = runit[2]*a[0] - runit[0]*a[2];
      b[2] = runit[0]*a[1] - runit[1]*a[0];
      disp[0] = a[0]*cosine  + b[0]*sine;
      disp[1] = a[1]*cosine  + b[1]*sine;
      disp[2] = a[2]*cosine  + b[2]*sine;

      pt = points[i].x;
      pt[0] = rpoint[0] + c[0] + disp[0];
      pt[1] = rpoint[1] + c[1] + disp[1];
      pt[2] = rpoint[2] + c[2] + disp[2];
    }

    for (int i = 0; i < nsurf; i++) {
      d[0] = xsurf_original[i][0] - rpoint[0];
      d[1] = xsurf_original[i][1] - rpoint[1];
      d[2] = xsurf_original[i][2] - rpoint[2];
      ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
      c[0] = ddotr*runit[0];
      c[1] = ddotr*runit[1];
      c[2] = ddotr*runit[2];
      a[0] = d[0] - c[0];
      a[1] = d[1] - c[1];
      a[2] = d[2] - c[2];
      b[0] = runit[1]*a[2] - runit[2]*a[1];
      b[1] = runit[2]*a[0] - runit[0]*a[2];
      b[2] = runit[0]*a[1] - runit[1]*a[0];
      disp[0] = a[0]*cosine  + b[0]*sine;
      disp[1] = a[1]*cosine  + b[1]*sine;
      disp[2] = a[2]*cosine  + b[2]*sine;

      xsurf[i][0] = rpoint[0] + c[0] + disp[0];
      xsurf[i][1] = rpoint[1] + c[1] + disp[1];
      xsurf[i][2] = rpoint[2] + c[2] + disp[2];
      vsurf[i][0] = omega_rotate * (runit[1]*disp[2] - runit[2]*disp[1]);
      vsurf[i][1] = omega_rotate * (runit[2]*disp[0] - runit[0]*disp[2]);
      vsurf[i][2] = omega_rotate * (runit[0]*disp[1] - runit[1]*disp[0]);
    }

    if (dimension == 3) {
      for (int i = 0; i < nsurf; i++) {
        p1 = points[tris[i].p1].x;
        p2 = points[tris[i].p2].x;
        p3 = points[tris[i].p3].x;
        MathExtra::sub3(p1,p2,p12);
        MathExtra::sub3(p1,p3,p13);
        MathExtra::cross3(p12,p13,tris[i].norm);
        MathExtra::norm3(tris[i].norm);
      }
    }
  }

  // trigger reneighbor if any point has moved skin/2 distance

  double dx,dy,dz,rsq;

  int triggerflag = 0;

  if (mstyle != NONE) {
    for (int i = 0; i < npoints; i++) {
      pt = points[i].x;
      dx = pt[0] - points_lastneigh[i][0];
      dy = pt[1] - points_lastneigh[i][1];
      dz = pt[2] - points_lastneigh[i][2];
      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq > triggersq) {
        triggerflag = 1;
        break;
      }
    }
  }

  if (triggerflag) next_reneighbor = update->ntimestep;
}

/* ----------------------------------------------------------------------
   build neighbor list for sphere/surf interactions
   I = sphere, J = surf
   similar to methods in neigh_gran.cpp
------------------------------------------------------------------------- */

void FixSurfaceGlobal::pre_neighbor()
{
  int i,j,m,n,nn,dnum,dnumbytes;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,rsq,radsum,cutsq;
  int *neighptr,*touchptr;
  double *shearptr;

  int *npartner;
  tagint **partner;
  double **shearpartner;
  int **firsttouch;
  double **firstshear;
  MyPage<int> *ipage_touch;
  MyPage<double> *dpage_shear;

  double **x = atom->x;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double skin = neighbor->skin;

  list->grow(nlocal,nall);
  if (use_history) listhistory->grow(nlocal,nall);

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  if (use_history) {
    fix_history->nlocal_neigh = nlocal;
    //npartner = fix_history->npartner;
    //partner = fix_history->partner;
    //shearpartner = fix_history->shearpartner;
    firsttouch = fix_history->firstflag;
    firstshear = fix_history->firstvalue;
    //ipage_touch = listhistory->ipage;
    //dpage_shear = listhistory->dpage;
    //dnum = listhistory->dnum;
    dnumbytes = dnum * sizeof(double);
  }

  // store current point positions for future neighbor trigger check
  // check is performed in intitial_integrate()

  if (mstyle != NONE) {
    for (i = 0; i < npoints; i++) {
      points_lastneigh[i][0] = points[i].x[0];
      points_lastneigh[i][1] = points[i].x[1];
      points_lastneigh[i][2] = points[i].x[2];
    }
  }

  int inum = 0;
  ipage->reset();
  if (use_history) {
    ipage_touch->reset();
    dpage_shear->reset();
  }

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    if (use_history) {
      nn = 0;
      touchptr = ipage_touch->vget();
      shearptr = dpage_shear->vget();
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // for now, loop over all surfs
    // NOTE: use a more sophisticated neighbor check

    for (j = 0; j < nsurf; j++) {
      delx = xtmp - xsurf[j][0];
      dely = ytmp - xsurf[j][1];
      delz = ztmp - xsurf[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radsurf[j] + skin;
      cutsq = radsum*radsum;
      if (rsq <= cutsq) {
        neighptr[n] = j;

        if (use_history) {
          if (rsq < radsum*radsum) {
            for (m = 0; m < npartner[i]; m++)
              if (partner[i][m] == j) break;
            if (m < npartner[i]) {
              touchptr[n] = 1;
              memcpy(&shearptr[nn],&shearpartner[i][dnum*m],dnumbytes);
              nn += dnum;
            } else {
              touchptr[n] = 0;
              memcpy(&shearptr[nn],zeroes,dnumbytes);
              nn += dnum;
            }
          } else {
            touchptr[n] = 0;
            memcpy(&shearptr[nn],zeroes,dnumbytes);
            nn += dnum;
          }
        }

        n++;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Fix surface/global neighbor list overflow, "
                 "boost neigh_modify one");

    if (use_history) {
      firsttouch[i] = touchptr;
      firstshear[i] = shearptr;
      ipage_touch->vgot(n);
      dpage_shear->vgot(nn);
    }
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   compute particle/surface interactions
   impart force and torque to spherical particles
------------------------------------------------------------------------- */

void FixSurfaceGlobal::post_force(int vflag)
{
  int i,j,k,ii,jj,inum,jnum,jflag,otherflag;
  double meff,factor_couple;
  double rsq,dr[3],contact[3],ds[3],vs[3],*forces,*torquesi;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch,touch_flag;
  double *history,*allhistory,**firsthistory;

  model->history_update = 1;
  if (update->setupflag) model->history_update = 0;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"surface/global:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    }
  }

  // loop over neighbors of my atoms
  // I is always sphere, J is always line

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (heat_flag) {
    if (tstr)
      Twall = input->variable->compute_equal(tvar);
    model->Tj = Twall;
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
    if (!(mask[i] & groupbit)) continue;
    model->xi = x[i];
    model->radi = radius[i];
    model->vi = v[i];
    model->omegai = omega[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      // Reset model and copy initial geometric data
      model->xj = xsurf[j];
      model->radj = radsurf[j];
      if (use_history) model->touch = touch[jj];

      touch_flag = model->check_contact();

      if (!touch_flag) {
        // unset non-touching neighbors
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      rsq = model->rsq;

      // contact computation for line or tri

      if (dimension == 2) {

        // check for overlap of sphere and line segment
        // jflag = 0 for no overlap, 1 for interior line pt, -1/-2 for end pts
        // if no overlap, just continue
        // for overlap, also return:
        //   contact = nearest point on line to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr
        // NOTE: different for line vs tri

        jflag = overlap_sphere_line(i,j,contact,dr,rsq);

        if (!jflag) {
          if (use_history) {
            touch[jj] = 0;
            history = &allhistory[size_history * jj];
            for (k = 0; k < size_history; k++) history[k] = 0.0;
          }
          continue;
        }

        // if contact = line end pt:
        // check overlap status of line adjacent to the end pt
        // otherflag = 0/1 for this/other line performs calculation
        // NOTE: different for line vs tri

        if (jflag < 0) {
          otherflag = endpt_neigh_check(i,j,jflag);
          if (otherflag) continue;
        }

        // NOTE: add logic to check for coupled contacts and weight them

        factor_couple = 1.0;

      } else {

        // check for overlap of sphere and triangle
        // jflag = 0 for no overlap, 1 for interior line pt,
        //   -1/-2/-3 for 3 edges, -4/-5/-6 for 3 corner pts
        // if no overlap, just continue
        // for overlap, also returns:
        //   contact = nearest point on tri to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr

        jflag = overlap_sphere_tri(i,j,contact,dr,rsq);

        if (!jflag) {
          if (use_history) {
            touch[jj] = 0;
            history = &allhistory[size_history * jj];
            for (k = 0; k < size_history; k++) history[k] = 0.0;
          }
          continue;
        }

        // if contact = tri edge or corner:
        // check status of the tri(s) adjacent to edge or corner
        // otherflag = 0/1 for this/other tri performs calculation

        if (jflag < 0) {
          if (jflag >= -3) {
            otherflag = edge_neigh_check(i,j,jflag);
            if (otherflag) continue;
          } else {
            otherflag = corner_neigh_check(i,j,jflag);
            if (otherflag) continue;
          }
        }

        // NOTE: add logic to check for coupled contacts and weight them

        factor_couple = 1.0;
      }

      // apply new rsq
      
      model->rsq = rsq;

      // meff = effective mass of sphere
      // if I is part of rigid body, use body mass

      meff = rmass[i];
      if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

      // copy additional information and prepare force calculations
      
      model->meff = meff;

      ds[0] = contact[0] - xsurf[j][0];
      ds[1] = contact[1] - xsurf[j][1];
      ds[2] = contact[2] - xsurf[j][2];

      vs[0] = vsurf[j][0] + (omegasurf[j][1] * ds[2] - omegasurf[j][2] * ds[1]);
      vs[1] = vsurf[j][1] + (omegasurf[j][2] * ds[0] - omegasurf[j][0] * ds[2]);
      vs[2] = vsurf[j][2] + (omegasurf[j][0] * ds[1] - omegasurf[j][1] * ds[0]);

      model->vj = vs;
      model->omegaj = omegasurf[j];

      if (heat_flag) model->Ti = temperature[i];

      // pairwise interaction between sphere and surface element

      if (use_history) {
        touch[jj] = 1;
        history = &allhistory[3*jj];
        model->history = history;
      }

      model->dx[0] = dr[0];
      model->dx[1] = dr[1];
      model->dx[2] = dr[2];

      // need to add support coupled contacts
      // is this just multiplying forces (+torques?) by factor_couple?
      
      model->calculate_forces();

      forces = model->forces;
      torquesi = model->torquesi;

      // apply forces & torques
      
      add3(f[i], forces, f[i]);
      add3(torque[i], torquesi, torque[i]);
      if (heat_flag) heatflow[i] += model->dq;
    }
  }
}

/* ----------------------------------------------------------------------
   turn on/off surface motion via fix_modify
   when motion is active, INTITIAL_INTEGRATE is set
------------------------------------------------------------------------- */

int FixSurfaceGlobal::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"move") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    int ifix = modify->find_fix(id);

    if (strcmp(arg[1],"none") == 0) {
      mstyle = NONE;
      modify->fmask[ifix] &= ~INITIAL_INTEGRATE;
      force_reneighbor = 0;
      next_reneighbor = -1;
      move_clear();
      return 2;

    } else if (strcmp(arg[1],"rotate") == 0) {
      if (narg < 9) error->all(FLERR,"Illegal fix_modify command");
      mstyle = ROTATE;
      modify->fmask[ifix] |= INITIAL_INTEGRATE;
      force_reneighbor = 1;
      next_reneighbor = -1;

      rpoint[0] = xscale * utils::numeric(FLERR,arg[2],false,lmp);
      rpoint[1] = yscale * utils::numeric(FLERR,arg[3],false,lmp);
      rpoint[2] = zscale * utils::numeric(FLERR,arg[4],false,lmp);
      raxis[0] = utils::numeric(FLERR,arg[5],false,lmp);
      raxis[1] = utils::numeric(FLERR,arg[6],false,lmp);
      raxis[2] = utils::numeric(FLERR,arg[7],false,lmp);
      rperiod = utils::numeric(FLERR,arg[8],false,lmp);
      if (rperiod <= 0.0) error->all(FLERR,"Illegal fix_modify command");

      if (dimension == 2)
        if (mstyle == ROTATE && (raxis[0] != 0.0 || raxis[1] != 0.0))
          error->all(FLERR,"Fix_modify cannot rotate around "
                     "non z-axis for 2d problem");

      time_origin = update->ntimestep;
      omega_rotate = MY_2PI / rperiod;

      // runit = unit vector along rotation axis

      if (mstyle == ROTATE) {
        double len = MathExtra::len3(raxis);
        if (len == 0.0)
          error->all(FLERR,"Zero length rotation vector with fix_modify");
        MathExtra::normalize3(raxis,runit);
      }

      move_clear();
      move_init();
      return 9;
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::reset_dt()
{
  if (mstyle != NONE)
    error->all(FLERR,"Resetting timestep size is not allowed with "
               "fix surface/global motion");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSurfaceGlobal::memory_usage()
{
  // NOTE: need to include neigh lists

  double bytes = 0.0;
  bytes += npoints*sizeof(Point);
  if (dimension == 2) {
    bytes += nlines*sizeof(Line);
    bytes = nlines*sizeof(Connect2d);
  } else {
    bytes += ntris*sizeof(Tri);
    bytes = ntris*sizeof(Connect3d);
  }
  return bytes;
}

/* ----------------------------------------------------------------------
   extract neighbor lists
------------------------------------------------------------------------- */

void *FixSurfaceGlobal::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"list") == 0) return list;
  else if (strcmp(str,"listhistory") == 0) return listhistory;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceGlobal::image(int *&ivec, double **&darray)
{
  int n;
  double *p1,*p2,*p3;

  if (dimension == 2) {
    n = nlines;

    if (imax == 0) {
      imax = n;
      memory->create(imflag,imax,"surface/global:imflag");
      memory->create(imdata,imax,7,"surface/global:imflag");
    }

    for (int i = 0; i < n; i++) {
      p1 = points[lines[i].p1].x;
      p2 = points[lines[i].p2].x;

      imflag[i] = LINE;
      imdata[i][0] = lines[i].type;
      imdata[i][1] = p1[0];
      imdata[i][2] = p1[1];
      imdata[i][3] = p1[2];
      imdata[i][4] = p2[0];
      imdata[i][5] = p2[1];
      imdata[i][6] = p2[2];
    }

  } else {
    n = ntris;

    if (imax == 0) {
      imax = n;
      memory->create(imflag,imax,"surface/global:imflag");
      memory->create(imdata,imax,10,"surface/global:imflag");
    }

    for (int i = 0; i < n; i++) {
      p1 = points[tris[i].p1].x;
      p2 = points[tris[i].p2].x;
      p3 = points[tris[i].p3].x;

      imflag[i] = TRI;
      imdata[i][0] = tris[i].type;
      imdata[i][1] = p1[0];
      imdata[i][2] = p1[1];
      imdata[i][3] = p1[2];
      imdata[i][4] = p2[0];
      imdata[i][5] = p2[1];
      imdata[i][6] = p2[2];
      imdata[i][7] = p3[0];
      imdata[i][8] = p3[1];
      imdata[i][9] = p3[2];
    }
  }

  ivec = imflag;
  darray = imdata;
  return n;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// 2d geometry methods
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

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

int FixSurfaceGlobal::overlap_sphere_line(int i, int j, double *pt,
                                          double *r, double &rsq)
{
  double a[3],b[3];

  // P1,P2 = end points of line segment

  double *p1 = points[lines[j].p1].x;
  double *p2 = points[lines[j].p2].x;

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

int FixSurfaceGlobal::endpt_neigh_check(int i, int j, int jflag)
{
  // ncheck = # of neighbor lines to check
  // neighs = indices of neighbor lines (including self)

  int ncheck;
  int *neighs;

  if (jflag == -1) {
    if (connect2d[j].np1 == 1) return 0;
    ncheck = connect2d[j].np1;
    neighs = connect2d[j].neigh_p1;
  } else if (jflag == -2) {
    if (connect2d[j].np2 == 1) return 0;
    ncheck = connect2d[j].np2;
    neighs = connect2d[j].neigh_p2;
  }

  // check overlap with each neighbor line
  // if any line has interior overlap, another line computes
  // if all lines have endpt overlap, line with lowest index computes
  // kflag = overlap status with neighbor line
  // kflag = 1, interior overlap
  // kflag = 0, no overlap, should not be possible
  // kflag < 0, overlap at endpt

  int k,kflag;
  double rsq;
  double dr[3],contact[3];

  int linemin = j;

  for (int m = 0; m < ncheck; m++) {
    k = neighs[m];
    if (k == j) continue;   // skip self line
    kflag = overlap_sphere_line(i,k,contact,dr,rsq);
    if (kflag > 0) return 1;
    if (kflag == 0) error->one(FLERR,"Fix surface/global neighbor line overlap is invalid");
    linemin = MIN(linemin,k);
  }

  if (j == linemin) return 0;
  return 1;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// 3d geometry methods
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

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

int FixSurfaceGlobal::overlap_sphere_tri(int i, int j, double *pt,
                                         double *r, double &rsq)
{
  int e12flag,e23flag,e31flag,o12flag,o23flag,o31flag;
  int esum,osum,lineflag;
  double dot;
  double a[3],point[3],edge[3],pvec[3],xproduct[3];

  // P1,P2 = end points of line segment
  // norm = current norm of triangle

  double *p1 = points[tris[j].p1].x;
  double *p2 = points[tris[j].p2].x;
  double *p3 = points[tris[j].p3].x;
  double *norm = tris[j].norm;

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

int FixSurfaceGlobal::nearest_point_line(double *x, double *p1, double *p2,
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

int FixSurfaceGlobal::edge_neigh_check(int i, int j, int jflag)
{
  // ncheck = # of neighbor tris to check
  // neighs = indices of neighbor tris (including self)

  int ncheck;
  int *neighs;

  if (jflag == -1) {
    if (connect3d[j].ne1 == 1) return 0;
    ncheck = connect3d[j].ne1;
    neighs = connect3d[j].neigh_e1;
  } else if (jflag == -2) {
    if (connect3d[j].ne2 == 1) return 0;
    ncheck = connect3d[j].ne2;
    neighs = connect3d[j].neigh_e2;
  } else if (jflag == -3) {
    if (connect3d[j].ne3 == 1) return 0;
    ncheck = connect3d[j].ne3;
    neighs = connect3d[j].neigh_e3;
  }

  // check overlap with each neighbor tri
  // if any tri has interior overlap, another tri computes
  // if all tris have edge overlap, tri with lowest index computes
  // kflag = overlap status with neighbor tri
  // kflag = 1, interior overlap
  // kflag = 0, no overlap, should not be possible
  // kflag < 0, overlap at edge (overlap at corner pt should not be possible)

  int k,kflag;
  double rsq;
  double dr[3],contact[3];

  int trimin = j;

  for (int m = 0; m < ncheck; m++) {
    k = neighs[m];
    if (k == j) continue;   // skip self tri
    kflag = overlap_sphere_tri(i,k,contact,dr,rsq);
    if (kflag > 0) return 1;
    if (kflag == 0) error->one(FLERR,"Fix surface/global neighbor tri overlap is invalid");
    trimin = MIN(trimin,k);
  }

  if (j == trimin) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check overlap status of sphere I with tri J versus any neighbor tris K
   I overlaps J at jflag = -4,-5,-6 for three corners
   return 0 if this tri J performs computation
   return 1 if some other tri K performs computation
------------------------------------------------------------------------- */

int FixSurfaceGlobal::corner_neigh_check(int i, int j, int jflag)
{
  // ncheck = # of neighbor tris to check
  // neighs = indices of neighbor tris (including self)

  int ncheck;
  int *neighs;

  if (jflag == -4) {
    if (connect3d[j].nc1 == 1) return 0;
    ncheck = connect3d[j].nc1;
    neighs = connect3d[j].neigh_c1;
  } else if (jflag == -5) {
    if (connect3d[j].nc2 == 1) return 0;
    ncheck = connect3d[j].nc2;
    neighs = connect3d[j].neigh_c2;
  } else if (jflag == -6) {
    if (connect3d[j].nc3 == 1) return 0;
    ncheck = connect3d[j].nc3;
    neighs = connect3d[j].neigh_c3;
  }

  // check overlap with each neighbor tri
  // if any tri has interior or edge overlap, another tri computes
  // if all tris have corner pt overlap, tri with lowest index computes
  // kflag = overlap status with neighbor tri
  // kflag = 1, interior overlap
  // kflag = 0, no overlap, should not be possible
  // kflag = -1/-2/-3, overlap at edge
  // kflag = -4/-5/-6, overlap at corner pt

  int k,kflag;
  double rsq;
  double dr[3],contact[3];

  int trimin = j;

  for (int m = 0; m < ncheck; m++) {
    k = neighs[m];
    if (k == j) continue;   // skip self tri
    kflag = overlap_sphere_tri(i,k,contact,dr,rsq);
    if (kflag > 0) return 1;
    if (kflag == 0) error->one(FLERR,"Fix surface/global neighbor tri overlap is invalid");
    if (kflag >= -3) return 1;
    trimin = MIN(trimin,k);
  }

  if (j == trimin) return 0;
  return 1;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// initializiation of surfs
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   extract lines or surfs from molecule template ID for one or more molecules
   concatenate into single list of lines and tris
   create list of unique points using hash
------------------------------------------------------------------------- */

void FixSurfaceGlobal::extract_from_molecules(char *molID)
{
  // populate global point/line/tri data structs

  points = nullptr;
  lines = nullptr;
  tris = nullptr;
  npoints = nlines = ntris = 0;
  int maxpoints = 0;

  int imol = atom->find_molecule(molID);
  if (imol == -1)
    error->all(FLERR,"Molecule template ID for fix surface/global does not exist");

  // loop over one or more molecules in molID

  Molecule **onemols = &atom->molecules[imol];
  int nmol = onemols[0]->nset;

  for (int m = 0; m < nmol; m++) {
    if (dimension == 2)
      if (onemols[m]->lineflag == 0)
        error->all(FLERR,"Fix surface/global molecule must have lines");
    if (dimension == 3)
      if (onemols[m]->triflag == 0)
        error->all(FLERR,"Fix surface/global molecule must have triangles");

    int nl = onemols[m]->nlines;
    int nt = onemols[m]->ntris;

    nlines += nl;
    ntris += nt;
    lines = (Line *) memory->srealloc(lines,nlines*sizeof(Line),
                                      "surface/global:lines");
    tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),
                                    "surface/global:tris");

    // create a map
    // key = xyz coords of a point
    // value = index into unique points vector

    std::map<std::tuple<double,double,double>,int> hash;

    // offset line/tri index lists by previous npoints
    // pi,p2,p3 are C-style indices into points vector

    if (dimension == 2) {
      int *molline = onemols[m]->molline;
      int *typeline = onemols[m]->typeline;
      double **epts = onemols[m]->lines;
      int iline = nlines - nl;

      for (int i = 0; i < nl; i++) {
        lines[iline].mol = molline[i];
        lines[iline].type = typeline[i];

        auto key = std::make_tuple(epts[i][0],epts[i][1],0.0);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = epts[i][0];
          points[npoints].x[1] = epts[i][1];
          points[npoints].x[2] = 0.0;
          lines[iline].p1 = npoints;
          npoints++;
        } else lines[iline].p1 = hash[key];

        key = std::make_tuple(epts[i][2],epts[i][3],0.0);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = epts[i][2];
          points[npoints].x[1] = epts[i][3];
          points[npoints].x[2] = 0.0;
          lines[iline].p2 = npoints;
          npoints++;
        } else lines[iline].p2 = hash[key];

        iline++;
      }
    }

    if (dimension == 3) {
      int *moltri = onemols[m]->moltri;
      int *typetri = onemols[m]->typetri;
      double **cpts = onemols[m]->tris;
      int itri = ntris - nt;

      for (int i = 0; i < nt; i++) {
        tris[itri].mol = moltri[i];
        tris[itri].type = typetri[i];

        auto key = std::make_tuple(cpts[i][0],cpts[i][1],cpts[i][2]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][0];
          points[npoints].x[1] = cpts[i][1];
          points[npoints].x[2] = cpts[i][2];
          tris[itri].p1 = npoints;
          npoints++;
        } else tris[itri].p1 = hash[key];

        key = std::make_tuple(cpts[i][3],cpts[i][4],cpts[i][5]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][3];
          points[npoints].x[1] = cpts[i][4];
          points[npoints].x[2] = cpts[i][5];
          tris[itri].p2 = npoints;
          npoints++;
        } else tris[itri].p2 = hash[key];

        key = std::make_tuple(cpts[i][6],cpts[i][7],cpts[i][8]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][6];
          points[npoints].x[1] = cpts[i][7];
          points[npoints].x[2] = cpts[i][8];
          tris[itri].p3 = npoints;
          npoints++;
        } else tris[itri].p3 = hash[key];

        itri++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   extract triangles from an STL file, can be text or binary
   create list of unique points using hash
------------------------------------------------------------------------- */

void FixSurfaceGlobal::extract_from_stlfile(char *filename)
{
  if (dimension == 2)
    error->all(FLERR,"Fix surface/global cannot use an STL file for 2d simulations");

  // read tris from STL file
  // stltris = tri coords internal to STL reader
  
  STLReader *stl = new STLReader(lmp);
  double **stltris;
  ntris = stl->read_file(filename,stltris);

  // create points and tris data structs

  points = nullptr;
  npoints = 0;
  int maxpoints = 0;

  tris = (Tri *) memory->smalloc(ntris*sizeof(Tri),"surface/global:tris");
  
  // create a map
  // key = xyz coords of a point
  // value = index into unique points vector
  
  std::map<std::tuple<double,double,double>,int> hash;

  // loop over STL tris
  // populate points and tris data structs
  // set molecule and type of tri = 1
  
  for (int itri = 0; itri < ntris; itri++) {
    tris[itri].mol = 1;
    tris[itri].type = 1;

    auto key = std::make_tuple(stltris[itri][0],stltris[itri][1],stltris[itri][2]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][0];
      points[npoints].x[1] = stltris[itri][1];
      points[npoints].x[2] = stltris[itri][2];
      tris[itri].p1 = npoints;
      npoints++;
    } else tris[itri].p1 = hash[key];

    key = std::make_tuple(stltris[itri][3],stltris[itri][4],stltris[itri][5]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][3];
      points[npoints].x[1] = stltris[itri][4];
      points[npoints].x[2] = stltris[itri][5];
      tris[itri].p2 = npoints;
      npoints++;
    } else tris[itri].p2 = hash[key];
    
    key = std::make_tuple(stltris[itri][6],stltris[itri][7],stltris[itri][8]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][6];
      points[npoints].x[1] = stltris[itri][7];
      points[npoints].x[2] = stltris[itri][8];
      tris[itri].p3 = npoints;
      npoints++;
    } else tris[itri].p3 = hash[key];
  }
  
  // delete STL reader
  
  delete stl;
}

/* ----------------------------------------------------------------------
   create and initialize Connect2d info for all lines
------------------------------------------------------------------------- */

void FixSurfaceGlobal::connectivity2d_global()
{
  connect2d = (Connect2d *) memory->smalloc(nlines*sizeof(Connect2d),
                                            "surface/global:connect2d");

  // setup end point connectivity lists
  // count # of lines containing each end point
  // create ragged 2d array to contain all point indices, then fill it
  // set neigh_p12 vector ptrs in connect2d to rows of ragged array
  // np12 counts and vectors include self line

  int *counts;
  memory->create(counts,npoints,"surface/global:count");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    counts[lines[i].p1]++;
    counts[lines[i].p2]++;
  }

  memory->create_ragged(plist,npoints,counts,"surface/global:plist");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    plist[lines[i].p1][counts[lines[i].p1]++] = i;
    plist[lines[i].p2][counts[lines[i].p2]++] = i;
  }

  for (int i = 0; i < nlines; i++) {
    connect2d[i].np1 = counts[lines[i].p1];
    if (connect2d[i].np1 == 1) connect2d[i].neigh_p1 = nullptr;
    else connect2d[i].neigh_p1 = plist[lines[i].p1];

    connect2d[i].np2 = counts[lines[i].p2];
    if (connect2d[i].np2 == 1) connect2d[i].neigh_p2 = nullptr;
    else connect2d[i].neigh_p2 = plist[lines[i].p2];
  }

  memory->destroy(counts);
}

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for all triangles
------------------------------------------------------------------------- */

void FixSurfaceGlobal::connectivity3d_global()
{
  int p1,p2,p3;

  connect3d = (Connect3d *) memory->smalloc(ntris*sizeof(Connect3d),
                                            "surface/global:connect3d");
  int **tri2edge;
  memory->create(tri2edge,ntris,3,"surfface/global::tri2edge");

  // create a map
  // key = <p1,p2> indices of 2 points, in either order
  // value = index into count of unique edges

  std::map<std::tuple<int,int>,int> hash;
  int nedges = 0;

  for (int i = 0; i < ntris; i++) {
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;

    auto key1 = std::make_tuple(p1,p2);
    auto key2 = std::make_tuple(p2,p1);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][0] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][0] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][0] = hash[key2];

    key1 = std::make_tuple(p2,p3);
    key2 = std::make_tuple(p3,p2);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][1] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][1] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][1] = hash[key2];

    key1 = std::make_tuple(p3,p1);
    key2 = std::make_tuple(p1,p3);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][2] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][2] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][2] = hash[key2];
  }

  // setup tri edge connectivity lists
  // count # of tris containing each edge
  // create ragged 2d array to contain all edge indices, then fill it
  // set neigh_e123 vector ptrs in connect3d to rows of ragged array
  // ne123 counts and vectors include self tri

  int *counts;
  memory->create(counts,nedges,"surface/global:count");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tri2edge[i][0]]++;
    counts[tri2edge[i][1]]++;
    counts[tri2edge[i][2]]++;
  }

  memory->create_ragged(elist,nedges,counts,"surface/global:elist");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    elist[tri2edge[i][0]][counts[tri2edge[i][0]]++] = i;
    elist[tri2edge[i][1]][counts[tri2edge[i][1]]++] = i;
    elist[tri2edge[i][2]][counts[tri2edge[i][2]]++] = i;
  }

  for (int i = 0; i < ntris; i++) {
    connect3d[i].ne1 = counts[tri2edge[i][0]];
    if (connect3d[i].ne1 == 1) connect3d[i].neigh_e1 = nullptr;
    else connect3d[i].neigh_e1 = elist[tri2edge[i][0]];

    connect3d[i].ne2 = counts[tri2edge[i][1]];
    if (connect3d[i].ne2 == 1) connect3d[i].neigh_e2 = nullptr;
    else connect3d[i].neigh_e2 = elist[tri2edge[i][1]];

    connect3d[i].ne3 = counts[tri2edge[i][2]];
    if (connect3d[i].ne3 == 1) connect3d[i].neigh_e3 = nullptr;
    else connect3d[i].neigh_e3 = elist[tri2edge[i][2]];
  }

  memory->destroy(counts);
  memory->destroy(tri2edge);

  // setup corner point connectivity lists
  // count # of tris containing each point
  // create ragged 2d array to contain all tri indices, then fill it
  // set neigh_c123 vector ptrs in connect3d to rows of ragged array
  // nc123 counts and vectors include self tri

  memory->create(counts,npoints,"surface/global:count");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tris[i].p1]++;
    counts[tris[i].p2]++;
    counts[tris[i].p3]++;
  }

  memory->create_ragged(clist,npoints,counts,"surface/global:clist");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    clist[tris[i].p1][counts[tris[i].p1]++] = i;
    clist[tris[i].p2][counts[tris[i].p2]++] = i;
    clist[tris[i].p3][counts[tris[i].p3]++] = i;
  }

  for (int i = 0; i < ntris; i++) {
    connect3d[i].nc1 = counts[tris[i].p1];
    if (connect3d[i].nc1 == 1) connect3d[i].neigh_c1 = nullptr;
    else connect3d[i].neigh_c1 = clist[tris[i].p1];

    connect3d[i].nc2 = counts[tris[i].p2];
    if (connect3d[i].nc2 == 1) connect3d[i].neigh_c2 = nullptr;
    else connect3d[i].neigh_c2 = clist[tris[i].p2];

    connect3d[i].nc3 = counts[tris[i].p3];
    if (connect3d[i].nc3 == 1) connect3d[i].neigh_c3 = nullptr;
    else connect3d[i].neigh_c3 = clist[tris[i].p3];
  }

  memory->destroy(counts);
}

/* ----------------------------------------------------------------------
   set attributes of all lines or tris
   xsurf,vsurf,omegasurf
   norm of triangles
------------------------------------------------------------------------- */

void FixSurfaceGlobal::surface_attributes()
{
  double delta[3],p12[3],p13[3];
  double *p1,*p2,*p3;

  memory->create(xsurf,nsurf,3,"surface/global:xsurf");
  memory->create(vsurf,nsurf,3,"surface/global:vsurf");
  memory->create(omegasurf,nsurf,3,"surface/global:omegasurf");
  memory->create(radsurf,nsurf,"surface/global:radsurf");

  if (dimension == 2) {
    for (int i = 0; i < nsurf; i++) {
      p1 = points[lines[i].p1].x;
      p2 = points[lines[i].p2].x;
      xsurf[i][0] = 0.5 * (p1[0]+p2[0]);
      xsurf[i][1] = 0.5 * (p1[1]+p2[1]);
      xsurf[i][2] = 0.0;

      MathExtra::sub3(p1,p2,delta);
      radsurf[i] = 0.5 * MathExtra::len3(delta);
    }

  } else {

    for (int i = 0; i < nsurf; i++) {
      p1 = points[tris[i].p1].x;
      p2 = points[tris[i].p2].x;
      p3 = points[tris[i].p3].x;
      xsurf[i][0] = (p1[0]+p2[0]+p3[0]) / 3.0;
      xsurf[i][1] = (p1[1]+p2[1]+p3[1]) / 3.0;
      xsurf[i][2] = (p1[2]+p2[2]+p3[2]) / 3.0;

      MathExtra::sub3(p1,xsurf[i],delta);
      radsurf[i] = MathExtra::lensq3(delta);
      MathExtra::sub3(p2,xsurf[i],delta);
      radsurf[i] = MAX(radsurf[i],MathExtra::lensq3(delta));
      MathExtra::sub3(p3,xsurf[i],delta);
      radsurf[i] = MAX(radsurf[i],MathExtra::lensq3(delta));
      radsurf[i] = sqrt(radsurf[i]);

      MathExtra::sub3(p1,p2,p12);
      MathExtra::sub3(p1,p3,p13);
      MathExtra::cross3(p12,p13,tris[i].norm);
      MathExtra::norm3(tris[i].norm);
    }
  }

  for (int i = 0; i < nsurf; i++) {
    vsurf[i][0] = vsurf[i][1] = vsurf[i][2] = 0.0;
    omegasurf[i][0] = omegasurf[i][1] = omegasurf[i][2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   setup for surface motion
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_init()
{
  memory->create(points_lastneigh,npoints,3,"surface/global:points_lastneigh");
  memory->create(points_original,npoints,3,"surface/global:points_original");
  memory->create(xsurf_original,nsurf,3,"surface/global:xsurf_original");

  for (int i = 0; i < npoints; i++) {
    points_lastneigh[i][0] = points_original[i][0] = points[i].x[0];
    points_lastneigh[i][1] = points_original[i][1] = points[i].x[1];
    points_lastneigh[i][2] = points_original[i][2] = points[i].x[2];
  }

  for (int i = 0; i < nsurf; i++) {
    xsurf_original[i][0] = xsurf[i][0];
    xsurf_original[i][1] = xsurf[i][1];
    xsurf_original[i][2] = xsurf[i][2];
    omegasurf[i][0] = omega_rotate*runit[0];
    omegasurf[i][1] = omega_rotate*runit[1];
    omegasurf[i][2] = omega_rotate*runit[2];
  }
}

/* ----------------------------------------------------------------------
   turn off surface motion and free memory
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_clear()
{
  // reset v,omega to zero

  for (int i = 0; i < nsurf; i++) {
    vsurf[i][0] = vsurf[i][1] = vsurf[i][2] = 0.0;
    omegasurf[i][0] = omegasurf[i][1] = omegasurf[i][2] = 0.0;
  }

  // deallocate memory

  memory->destroy(points_lastneigh);
  memory->destroy(points_original);
  memory->destroy(xsurf_original);
  points_lastneigh = nullptr;
  points_original = nullptr;
  xsurf_original = nullptr;
}

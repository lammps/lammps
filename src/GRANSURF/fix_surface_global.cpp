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
   Contributing authors: Joel Clemmer (SNL), Dan Bolintineanu (SNL)
----------------------------------------------------------------------- */

#include "fix_surface_global.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "graphics.h"
#include "input.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "my_page.h"
#include "nbin_manual.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "nstencil_manual.h"
#include "region.h"
#include "stl_reader.h"
#include "surf_extra.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include <algorithm>
#include <cstring>
#include <map>
#include <unordered_set>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace Granular_NS;
using namespace MathConst;
using namespace MathExtra;
using namespace SurfExtra;

enum{LINEAR,WIGGLE,ROTATE,TRANSROT,VARIABLE};
enum{INTERNAL = 0,EXTERNAL,UNCONNECTED};
enum{NSQ, BIN};

static constexpr double FLATTHRESH = 0.00015230484360876085; // = 1.0-cos(MY_PI/180.0); = 1 degree
static constexpr int DELTAMODEL = 4;
static constexpr int DELTAMOTION = 4;
static constexpr int MAXSURFTYPE = 1024;  // limit for # of surf types
static constexpr int MAXSURFMOL = 1024;   // limit for # of surf mol IDs
static constexpr double BIG = 1.0e20;
static constexpr double EPSILON = 1e-12;

static inline int FLIPSIDE(int nside) {
  if (nside == FixSurface::OPPOSITE_SIDE) return FixSurface::SAME_SIDE;
  else return FixSurface::OPPOSITE_SIDE;
}

/* ---------------------------------------------------------------------- */

FixSurfaceGlobal::FixSurfaceGlobal(LAMMPS *lmp, int narg, char **arg) :
  FixSurface(lmp, narg, arg), tstr(nullptr), nb(nullptr), ns(nullptr)
{
  if (!atom->radius_flag || !atom->omega_flag)
    error->all(FLERR, 2, "Fix surface/global requires atom attributes radius and omega");

  dimension = domain->dimension;

  // process one or more inputs
  // read triangles/lines from molecule template IDs or STL files
  // hash = map to store unique points
  //   key = xyz coords of a point
  //   value = index into unique points vector

  npoints = maxpoints = 0;
  nlines = ntris = 0;
  points = nullptr;
  lines = nullptr;
  tris = nullptr;
  last_setup_bins = -1;

  neigh_style = BIN;

  int ninput = 0;
  std::map<std::tuple<double,double,double,int>,int> hash;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"input") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix surface/global input", error);
      if (strcmp(arg[iarg+1],"mol") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "fix surface/global input mol", error);
        extract_from_molecule(arg[iarg+2],&hash,npoints,maxpoints,points,nlines,lines,ntris,tris);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"stl") == 0) {
        if (iarg+5 > narg) utils::missing_cmd_args(FLERR, "fix surface/global input stl", error);
        int stype = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        int smol = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
        extract_from_stlfile(arg[iarg+4],stype,(tagint) smol,&hash,npoints,maxpoints,points,ntris,tris);
        iarg += 5;
      } else error->all(FLERR, iarg+1, "Unknown fix surface/global input keyword: {}", arg[iarg+1]);
    } else break;

    ninput++;
  }

  if (ninput == 0)
    error->all(FLERR, 2, "Fix surface/global command requires using the input keyword");

  // process one or more granular models
  // disable bonded/history option for now

  class GranularModel* model;
  models = nullptr;
  nmodel = maxmodel = 0;
  heat_flag = 0;
  use_history = 0;
  size_history = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "fix surface/global model", error);

      if (nmodel == maxmodel) {
        maxmodel += DELTAMODEL;
        modeltypes = (ModelTypes *)
          memory->srealloc(models,maxmodel*sizeof(ModelTypes),"surf/global:modeltypes");
        models = (Granular_NS::GranularModel **)
          memory->srealloc(models,maxmodel*sizeof(Granular_NS::GranularModel *),
                           "surf/global:models");
      }

      models[nmodel] = model = new GranularModel(lmp);

      // assign range of particle and surf types to this model
      // use MAXSURFTYPE for now, in case smax keyword extends input surf types

      utils::bounds(FLERR, arg[iarg+1], 1, atom->ntypes,
                    modeltypes[nmodel].plo, modeltypes[nmodel].phi, error);
      utils::bounds(FLERR, arg[iarg+2], 1, MAXSURFTYPE,
                    modeltypes[nmodel].slo, modeltypes[nmodel].shi, error);
      nmodel++;

      model->contact_type = SURFACE;

      int classic_flag = 1;
      iarg += 3;
      if (strcmp(arg[iarg], "granular") == 0) {
        classic_flag = 0;
        iarg += 1;
      }

      if (classic_flag) {
        iarg = model->define_classic_model(arg, iarg, narg);
        if (iarg < narg && strcmp(arg[iarg], "limit_damping") == 0) {
          model->limit_damping = 1;
          iarg++;
        }
      } else {
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
          } else if (strcmp(arg[iarg], "limit_damping") == 0) {
            model->limit_damping = 1;
            iarg++;
          } else break;
        }
      }

      // define default damping sub model
      // if unspecified, takes no args

      if (!model->damping_model) model->construct_sub_model("viscoelastic", DAMPING);

      model->init();

      if (model->beyond_contact) size_history = MAX(size_history + 1, model->size_history);
      else size_history = MAX(size_history, model->size_history);
      if (model->size_history != 0) use_history = 1;
      if (model->beyond_contact)
        error->all(FLERR, "Granular models that extend beyond contact (e.g. JKR) not currenty supported");

    } else break;
  }
  if (nmodel == 0)
    error->all(FLERR, 2, "Fix surface/global command requires model keyword");

  // maxsurftype = max surf type of any input surf (for now)

  maxsurftype = 0;
  if (dimension == 2) {
    for (int i = 0; i < nlines; i++)
      maxsurftype = MAX(maxsurftype,lines[i].type);
  } else {
    for (int i = 0; i < ntris; i++)
      maxsurftype = MAX(maxsurftype,tris[i].type);
  }

  // maxsurfmol = max surf molecule ID of any input surf (for now)

  maxsurfmol = 0;
  if (dimension == 2) {
    for (int i = 0; i < nlines; i++)
      maxsurfmol = MAX(maxsurfmol,lines[i].mol);
  } else {
    for (int i = 0; i < ntris; i++)
      maxsurfmol = MAX(maxsurfmol,tris[i].mol);
  }

  // optional command-line args
  // smaxtype overrides max surf type of input surfs
  // smaxmol overrides max surf mol ID of input surfs
  // flat overrides FLATTHRESH of one degree

  int Twall_defined = 0;
  flatthresh = FLATTHRESH;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"smaxtype") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix surface/global smaxtype", error);
      int smaxtype = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (smaxtype > MAXSURFTYPE)
        error->all(FLERR, iarg+1, "Fix surface/global smaxtype > MAXSURFTYPE");
      if (smaxtype < maxsurftype)
        error->all(FLERR, iarg+1, "Fix surface/global smaxtype < input surf types");
      maxsurftype = smaxtype;
      iarg += 2;
    } else if (strcmp(arg[iarg],"smaxmol") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix surface/global smaxmol", error);
      tagint smaxmol = utils::tnumeric(FLERR,arg[iarg+1],false,lmp);
      if (smaxmol > MAXSURFMOL)
        error->all(FLERR, iarg+1, "Fix surface/global smaxmol > MAXSURFMOL");
      if (smaxmol < maxsurfmol)
        error->all(FLERR, iarg+1, "Fix surface/global smaxmol < input surf types");
      maxsurfmol = smaxmol;
      iarg += 2;
    } else if (strcmp(arg[iarg],"flat") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix surface/global flat", error);
      double flat = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if ((flat < 0.0) || (flat > 90.0))
        error->all(FLERR,"Invalid value {} for fix surface/global flat keyword", flat);
      flatthresh = 1.0 - cos(MY_PI*flat/180.0);
      iarg += 2;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix surface/global temperature", error);
      if (utils::strmatch(arg[iarg+1], "^v_")) {
        tstr = utils::strdup(arg[iarg+1] + 2);
      } else {
        Twall = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      }
      Twall_defined = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"neighbor") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix surface/global neighbor", error);
      if (strcmp(arg[iarg + 1],"nsq") == 0) {
        neigh_style = NSQ;
      } else if (strcmp(arg[iarg + 1],"bin") == 0) {
        neigh_style = BIN;
      } else
        error->all(FLERR, iarg+1, "Unknown fix surface/global neighbor option {}", arg[iarg + 1]);

      iarg += 2;
    } else error->all(FLERR,"Illegal fix surface/global command: {}", arg[iarg]);
  }

  if (heat_flag && !Twall_defined)
    error->all(FLERR, "Must define wall temperature with a heat model");

  // reset modeltypes shi from MAXSURFTYPE to maxsurtype
  // initialize types2model for all particle/surf type pairs
  // check that a model has been assigned to every type pair

  for (int i = 0; i < nmodel; i++)
    if (modeltypes[i].shi == MAXSURFTYPE) modeltypes[i].shi = maxsurftype;

  types2model = new Granular_NS::GranularModel**[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++)
    types2model[i] = new Granular_NS::GranularModel*[maxsurftype+1];
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= maxsurftype; j++)
      types2model[i][j] = nullptr;

  for (int m = 0; m < nmodel; m++)
    for (int i = modeltypes[m].plo; i <= modeltypes[m].phi; i++)
      for (int j = modeltypes[m].slo; j <= modeltypes[m].shi; j++)
        types2model[i][j] = models[m];

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= maxsurftype; j++)
      if (!types2model[i][j])
        error->all(FLERR,"Fix surface/global type pair is missing a granular model");

  // initializations

  if (dimension == 2) nsurf = nlines;
  else nsurf = ntris;

  nmotion = maxmotion = 0;
  motions = nullptr;
  anymove = anymove_variable = 0;

  points_lastneigh = nullptr;
  points_original = nullptr;
  xsurf_original = nullptr;
  pointmove = nullptr;

  neigh_p1 = neigh_p2 = nullptr;
  pwhich_p1 = pwhich_p2 = nullptr;
  nside_p1 = nside_p2 = nullptr;
  aflag_p1 = aflag_p2 = nullptr;
  fflag_p1 = fflag_p2 = nullptr;

  neigh_e1 = neigh_e2 = neigh_e3 = nullptr;
  ewhich_e1 = ewhich_e2 = ewhich_e3 = nullptr;
  nside_e1 = nside_e2 = nside_e3 =nullptr;
  aflag_e1 = aflag_e2 = aflag_e3 = nullptr;
  fflag_e1 = fflag_e2 = fflag_e3 = nullptr;
  neigh_c1 = neigh_c2 = neigh_c3 = nullptr;
  cwhich_c1 = cwhich_c2 = cwhich_c3 = nullptr;
  fflag_c1 = fflag_c2 = fflag_c3 = nullptr;

  connect2d = nullptr;
  connect3d = nullptr;

  xsurf = vsurf = omegasurf = nullptr;
  radsurf = nullptr;
  idsurf = nullptr;

  nsurf_ghost = -1;

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

  mol2motion = new int[maxsurfmol+1];
  for (int i = 0; i <= maxsurfmol; i++)
    mol2motion[i] = -1;

  firsttime = 1;

  // initialize surface attributes

  surface_attributes();

  // error checks on duplicate surfs or zero-size surfs

  if (dimension == 2) check2d();
  else check3d();

  // compute connectivity of triangles/lines
  // create Connect3d or Connect2d data structs

  if (dimension == 2) {
    connectivity2d_global(npoints,nlines,lines,connect2d,neigh_p1,neigh_p2);
    connectivity2d_complete();
  } else {
    nedges = connectivity3d_global(npoints,ntris,tris,connect3d,
                                   neigh_e1,neigh_e2,neigh_e3,
                                   neigh_c1,neigh_c2,neigh_c3);
    connectivity3d_complete();
  }

  // warn if any connections between surfs with different molIDs

  check_molecules();

  // print stats on surfs and their connectivity

  if (dimension == 2) stats2d();
  else stats3d();
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
  memory->destroy(pointmove);

  memory->sfree(modeltypes);
  for (int i = 0; i < nmodel; i++) delete models[i];
  memory->sfree(models);

  for (int i = 1; i <= atom->ntypes; i++) delete [] types2model[i];
  delete [] types2model;

  memory->destroy(neigh_p1);
  memory->destroy(neigh_p2);
  memory->destroy(pwhich_p1);
  memory->destroy(pwhich_p2);
  memory->destroy(nside_p1);
  memory->destroy(nside_p2);
  memory->destroy(aflag_p1);
  memory->destroy(aflag_p2);
  memory->destroy(fflag_p1);
  memory->destroy(fflag_p2);

  memory->destroy(neigh_e1);
  memory->destroy(neigh_e2);
  memory->destroy(neigh_e3);
  memory->destroy(ewhich_e1);
  memory->destroy(ewhich_e2);
  memory->destroy(ewhich_e3);
  memory->destroy(nside_e1);
  memory->destroy(nside_e2);
  memory->destroy(nside_e3);
  memory->destroy(aflag_e1);
  memory->destroy(aflag_e2);
  memory->destroy(aflag_e3);
  memory->destroy(fflag_e1);
  memory->destroy(fflag_e2);
  memory->destroy(fflag_e3);

  memory->destroy(neigh_c1);
  memory->destroy(neigh_c2);
  memory->destroy(neigh_c3);
  memory->destroy(cwhich_c1);
  memory->destroy(cwhich_c2);
  memory->destroy(cwhich_c3);
  memory->destroy(fflag_c1);
  memory->destroy(fflag_c2);
  memory->destroy(fflag_c3);

  memory->sfree(connect2d);
  memory->sfree(connect3d);

  memory->destroy(xsurf);
  memory->destroy(vsurf);
  memory->destroy(omegasurf);
  memory->destroy(radsurf);
  memory->destroy(idsurf);

  memory->destroy(mass_rigid);

  for (int i = 0; i < nmotion; i++) {
    if (motions[i].mstyle == VARIABLE) {
      delete [] motions[i].xvarstr;
      delete [] motions[i].yvarstr;
      delete [] motions[i].zvarstr;
      delete [] motions[i].vxvarstr;
      delete [] motions[i].vyvarstr;
      delete [] motions[i].vzvarstr;
    }
  }

  memory->sfree(motions);
  delete [] mol2motion;

  delete list;
  delete listhistory;
  delete [] zeroes;
  delete [] tstr;

  delete nb;
  delete ns;

  if (use_history)
    modify->delete_fix(id_fix);

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
    id_fix = std::string("NEIGH_HISTORY_SURFACE_GLOBAL_") + id;
    auto cmd = fmt::format("{} all NEIGH_HISTORY {} onesided surface/global", id_fix, size_history);
    fix_history = dynamic_cast<FixNeighHistory *>(modify->add_fix(cmd));
  } else
    fix_history = nullptr;
}

/* ----------------------------------------------------------------------
   mask for INITIAL_INTEGRATE will be set by fix_modify move
---------------------------------------------------------------------- */

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
  // global surface ghosts are created with orthogonal box offsets (see
  // pre_neighbor()), so triclinic boxes are not yet supported

  if (domain->triclinic)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix surface/global does not yet support triclinic simulation boxes");

  dt = update->dt;
  triggersq = 0.25 * neighbor->skin * neighbor->skin;

  // check for compatible heat conduction atom style

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR, Error::NOLASTLINE,
                 "Heat conduction in fix surface/global requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR, Error::NOLASTLINE,
                 "Heat conduction in fix surface/global requires atom style with heatflow property");
  }

  class GranularModel* model;
  int next_index = 0;
  for (int n = 0; n < nmodel; n++) {
    model = models[n];
    for (int i = 0; i < NSUBMODELS; i++) {
      model->sub_models[i]->history_index = next_index;
      next_index += model->sub_models[i]->size_history;
    }
  }
  model->dt = update->dt;

  // one-time setup and allocation of neighbor list
  // wait until now, so neighbor settings have been made
  // normally done in Neighbor->init_pair(), but this list is not registered by the Neighbor class

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
    if (tvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable {} for fix surface/global does not exist", tstr);
    if (!input->variable->equalstyle(tvar))
      error->all(FLERR, "Variable {} for fix surface/global must be an equal style variable", tstr);
  }

  // check on motion variables

  for (int i = 0; i < nmotion; i++) {
    Motion *motion = &motions[i];
    if (motion->mstyle != VARIABLE) continue;

    if (motion->xvarstr) {
      motion->xvar = input->variable->find(motion->xvarstr);
      if (motion->xvar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix_modify move does not exist", motion->xvarstr);
      if (!input->variable->equalstyle(motion->xvar))
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable {} for fix_modify move must be equal style", motion->xvarstr);
    }
    if (motion->yvarstr) {
      motion->yvar = input->variable->find(motion->yvarstr);
      if (motion->yvar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix_modify move does not exist", motion->yvarstr);
      if (!input->variable->equalstyle(motion->yvar))
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable {} for fix_modify move must be equal style", motion->yvarstr);
    }
    if (motion->zvarstr) {
      motion->zvar = input->variable->find(motion->zvarstr);
      if (motion->zvar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix_modify move does not exist", motion->zvarstr);
      if (!input->variable->equalstyle(motion->zvar))
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable {} for fix_modify move must be equal style", motion->zvarstr);
    }

    if (motion->vxvarstr) {
      motion->vxvar = input->variable->find(motion->vxvarstr);
      if (motion->vxvar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix_modify move does not exist", motion->vxvarstr);
      if (!input->variable->equalstyle(motion->vxvar))
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable {} for fix_modify move must be equal style", motion->vxvarstr);
    }
    if (motion->vyvarstr) {
      motion->vyvar = input->variable->find(motion->vyvarstr);
      if (motion->vyvar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix_modify move does not exist", motion->vyvarstr);
      if (!input->variable->equalstyle(motion->vyvar))
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable {} for fix_modify move must be equal style", motion->vyvarstr);
    }
    if (motion->vzvarstr) {
      motion->vzvar = input->variable->find(motion->vzvarstr);
      if (motion->vzvar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix_modify move does not exist", motion->vzvarstr);
      if (!input->variable->equalstyle(motion->vzvar))
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable {} for fix_modify move must be equal style", motion->vzvarstr);
    }
  }

  // initialize pointmove settings
  // done before every run, b/c fix_modify move can be used between runs

  if (anymove) {
    for (int i = 0; i < npoints; i++) pointmove[i] = 0;
  }

  // check for inconsistent surf motion
  // error if two surfs are connected but are assigned to different motions
  // check before every run, b/c fix_modify move can be used between runs

  int ecount = 0;

  if (anymove && dimension == 2) {
    int iline,jline,np1,np2,j,imotion,jmotion,imol,jmol;

    // check connections to either endpoint of each line
    // avoid double counting by jline < iline check

    for (iline = 0; iline < nlines; iline++) {
      imol = lines[iline].mol;
      imotion = mol2motion[imol];

      np1 = connect2d[iline].np1;
      np2 = connect2d[iline].np2;

      for (j = 0; j < np1; j++) {
        jline = connect2d[iline].neigh_p1[j];
        if (jline < iline) continue;
        jmol = lines[jline].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }
      for (j = 0; j < np2; j++) {
        jline = connect2d[iline].neigh_p2[j];
        if (jline < iline) continue;
        jmol = lines[jline].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }
    }
  }

  if (anymove && dimension == 3) {
    int itri,jtri,ne1,ne2,ne3,nc1,nc2,nc3,j,imotion,jmotion,imol,jmol;

    // check connections to each tri's 3 edges and 3 corner pts
    // avoid double counting by jtri < itri check
    // but edge and corner connections both contribute to error count

    for (itri = 0; itri < ntris; itri++) {
      imol = tris[itri].mol;
      imotion = mol2motion[imol];

      ne1 = connect3d[itri].ne1;
      ne2 = connect3d[itri].ne2;
      ne3 = connect3d[itri].ne3;

      for (j = 0; j < ne1; j++) {
        jtri = connect3d[itri].neigh_e1[j];
        if (jtri < itri) continue;
        jmol = tris[jtri].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }
      for (j = 0; j < ne2; j++) {
        jtri = connect3d[itri].neigh_e2[j];
        if (jtri < itri) continue;
        jmol = tris[jtri].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }
      for (j = 0; j < ne3; j++) {
        jtri = connect3d[itri].neigh_e3[j];
        if (jtri < itri) continue;
        jmol = tris[jtri].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }

      nc1 = connect3d[itri].nc1;
      nc2 = connect3d[itri].nc2;
      nc3 = connect3d[itri].nc3;

      for (j = 0; j < nc1; j++) {
        jtri = connect3d[itri].neigh_c1[j];
        if (jtri < itri) continue;
        jmol = tris[jtri].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }
      for (j = 0; j < nc1; j++) {
        jtri = connect3d[itri].neigh_c2[j];
        if (jtri < itri) continue;
        jmol = tris[jtri].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }
      for (j = 0; j < nc3; j++) {
        jtri = connect3d[itri].neigh_c3[j];
        if (jtri < itri) continue;
        jmol = tris[jtri].mol;
        jmotion = mol2motion[jmol];
        if (imotion != jmotion) ecount++;
      }
    }
  }

  if (ecount) error->all(FLERR, Error::NOLASTLINE,
                         "Fix surface/global motions assigned to "
                         "{} connected surf pairs are inconsistent",ecount);
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

void FixSurfaceGlobal::initial_integrate(int /*vflag*/)
{
  int imotion,mstyle;

  // invoke variables for any VARIABLE style motion

  if (anymove_variable) {
    for (int i = 0; i < nmotion; i++) {
      Motion *motion = &motions[i];
      if (motion->mstyle != VARIABLE) continue;
      if (motion->xvarstr) motion->dx = input->variable->compute_equal(motion->xvar);
      if (motion->yvarstr) motion->dy = input->variable->compute_equal(motion->yvar);
      if (motion->zvarstr) motion->dz = input->variable->compute_equal(motion->zvar);
      if (motion->vxvarstr) motion->vx = input->variable->compute_equal(motion->vxvar);
      if (motion->vyvarstr) motion->vy = input->variable->compute_equal(motion->vyvar);
      if (motion->vzvarstr) motion->vz = input->variable->compute_equal(motion->vzvar);
    }
  }

  // invoke appropriate move option for each surf

  for (int i = 0; i < nsurf; i++) {
    if (dimension == 2) imotion = mol2motion[lines[i].mol];
    else imotion = mol2motion[tris[i].mol];
    if (imotion < 0) continue;

    mstyle = motions[imotion].mstyle;
    if (mstyle == LINEAR) move_linear(imotion,i);
    else if (mstyle == WIGGLE) move_wiggle(imotion,i);
    else if (mstyle == ROTATE) move_rotate(imotion,i);
    else if (mstyle == TRANSROT) move_transrotate(imotion,i);
    else if (mstyle == VARIABLE) move_variable(imotion,i);
  }

  // clear pointmove settings

  if (dimension == 2) {
    for (int i = 0; i < nlines; i++) {
      if (mol2motion[lines[i].mol] < 0) continue;
      pointmove[lines[i].p1] = 0;
      pointmove[lines[i].p2] = 0;
    }
  } else {
    for (int i = 0; i < ntris; i++) {
      if (mol2motion[tris[i].mol] < 0) continue;
      pointmove[tris[i].p1] = 0;
      pointmove[tris[i].p2] = 0;
      pointmove[tris[i].p3] = 0;
    }
  }

  // trigger reneighbor if any point has moved skin/2 distance

  double dx,dy,dz,rsq;
  double *pt;

  int triggerflag = 0;

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

  if (triggerflag) next_reneighbor = update->ntimestep;
}

/* ----------------------------------------------------------------------
   build neighbor list for sphere/surf interactions
   I = sphere, J = surf
   similar to methods in neigh_gran.cpp
------------------------------------------------------------------------- */

void FixSurfaceGlobal::pre_neighbor()
{
  int i,j,j2,k,n,nn,dnum,dnumbytes;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,rsq,radsum,cutsq;
  int *neighptr;
  double *valueptr;

  int *npartner;
  tagint **partner;
  double **valuepartner;
  int **firstflag;
  double **firstvalue;
  MyPage<tagint> *ipage_atom;
  MyPage<double> *dpage_atom;

  double **x = atom->x;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double skin = neighbor->skin;
  double dot, dx_proj[3], cutsq_proj;

  list->grow(nlocal,nall);
  if (use_history) listhistory->grow(nlocal,nall);

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  if (use_history) {
    fix_history->otherlist = list;
    fix_history->nlocal_neigh = nlocal;
    npartner = fix_history->get_npartner();         // # of touching partners of each atom
    partner = fix_history->get_partner();           // global atom IDs for the partners
    valuepartner = fix_history->get_valuepartner(); // values for the partners
    ipage_atom = fix_history->get_ipage_atom();     // pages of partner atom IDs
    dpage_atom = fix_history->get_dpage_atom();     // pages of partner values
    firstflag = fix_history->firstflag;       // ptr to each atom's neighbor flag
    firstvalue = fix_history->firstvalue;     // ptr to each atom's values
    dnum = fix_history->get_dnum();
    dnumbytes = dnum * sizeof(double);
  }

  // store current point positions for future neighbor trigger check
  // check is performed in intitial_integrate()

  if (anymove) {
    for (i = 0; i < npoints; i++) {
      points_lastneigh[i][0] = points[i].x[0];
      points_lastneigh[i][1] = points[i].x[1];
      points_lastneigh[i][2] = points[i].x[2];
    }
  }

  int inum = 0;
  ipage->reset();
  if (use_history) {
    ipage_atom->reset();
    dpage_atom->reset();
  }

  if (neigh_style == BIN) {
    if (nb == nullptr) {
      double rmax_surf = 0.0, rmax_surf_all;
      for (j = 0; j < nsurf; j++)
        rmax_surf = MAX(rmax_surf, radsurf[j]);
      MPI_Allreduce(&rmax_surf, &rmax_surf_all, 1, MPI_DOUBLE, MPI_MAX, world);

      // Does NOT yet account for pour/deposit
      double rmax_atom = 0.0, rmax_atom_all;
      for (i = 0; i < nlocal; i++)
        rmax_atom = MAX(rmax_atom, radius[i]);
      MPI_Allreduce(&rmax_atom, &rmax_atom_all, 1, MPI_DOUBLE, MPI_MAX, world);

      //cutneighmax equiv
      double cutoff = skin + rmax_atom_all + rmax_surf_all;

      if (nsurf_ghost == -1) {
        nsurf_ghost = 0;

        // get limits for creating ghosts
        int a;
        double *boxlo = domain->boxlo;
        double *boxhi = domain->boxhi;
        double limitlo[3], limithi[3];
        for (a = 0; a < 3; a++) {
          limitlo[a] = boxlo[a] + cutoff;
          limithi[a] = boxhi[a] - cutoff;
        }

        // count number of ghosts and allocate room
        int ix, iy, iz;
        int imin[3], imax[3];
        int dimension = domain->dimension;
        int *periodicity = domain->periodicity;
        double *prd = domain->prd;
        for (i = 0; i < nsurf; i++) {
          for (a = 0; a < 3; a++) {
            imin[a] = imax[a] = 0;
            if (periodicity[a] && xsurf[i][a] < limitlo[a]) imax[a] = +1;
            if (periodicity[a] && xsurf[i][a] > limithi[a]) imin[a] = -1;
          }

          if (dimension == 2) imin[2] = imax[2] = 0;

          for (ix = imin[0]; ix <= imax[0]; ix++)
            for (iy = imin[1]; iy <= imax[1]; iy++)
              for (iz = imin[2]; iz <= imax[2]; iz++)
                nsurf_ghost += 1;
        }

        memory->grow(xsurf, nsurf + nsurf_ghost, 3, "surface/global:xsurf");
        memory->create(idsurf, nsurf + nsurf_ghost, "surface/global:idsurf");

        // add ghost coordinates
        j = nsurf;
        for (i = 0; i < nsurf; i++) {
          idsurf[i] = i;

          for (a = 0; a < 3; a++) {
            imin[a] = imax[a] = 0;
            if (periodicity[a] && xsurf[i][a] < limitlo[a]) imax[a] = +1;
            if (periodicity[a] && xsurf[i][a] > limithi[a]) imin[a] = -1;
          }

          if (dimension == 2) imin[2] = imax[2] = 0;

          for (ix = imin[0]; ix <= imax[0]; ix++) {
            for (iy = imin[1]; iy <= imax[1]; iy++) {
              for (iz = imin[2]; iz <= imax[2]; iz++) {
                xsurf[j][0] = xsurf[i][0] + ix * prd[0];
                xsurf[j][1] = xsurf[i][1] + iy * prd[1];
                xsurf[j][2] = xsurf[i][2] + iz * prd[2];
                idsurf[j] = i;
                // Todo: edit for triclinic
                j += 1;
              }
            }
          }
        }
      }

      nb = new NBinManual(lmp);
      ns = new NStencilManual(lmp);
      ns->nb = nb;

      nb->assign_neighbor_info(cutoff);
      ns->assign_neighbor_info(cutoff);
    }

    if (last_setup_bins < 0 || domain->box_change) {
      nb->setup_bins(0);
      ns->create_setup();
      ns->create();
      last_setup_bins = update->ntimestep;

      nb->bin_custom_setup(xsurf, nsurf + nsurf_ghost);
      nb->bin_custom(xsurf, nsurf + nsurf_ghost);
    } else if (anymove) {
      nb->bin_custom_setup(xsurf, nsurf + nsurf_ghost);
      nb->bin_custom(xsurf, nsurf + nsurf_ghost);
    }

    nb->bin_atoms_setup(nlocal);
    nb->bin_atoms();

    int *atom2bin = nb->atom2bin;
    int *binhead_surf = nb->binhead_custom;
    int *bins_surf = nb->bins_custom;

    int nstencil = ns->nstencil;
    int *stencil = ns->stencil;

    int ibin, bin_start;
    std::set<int> jadded;

    for (i = 0; i < nlocal; i++) {
      n = 0;
      neighptr = ipage->vget();
      if (use_history) {
        nn = 0;
        valueptr = dpage_atom->vget();
      }

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      radi = radius[i];
      cutsq_proj = radi + skin;
      cutsq_proj *= cutsq_proj;

      ibin = atom2bin[i];
      jadded.clear();

      for (k = 0; k < nstencil; k++) {
        bin_start = binhead_surf[ibin + stencil[k]];
        for (j = bin_start; j >= 0; j = bins_surf[j]) {
          j2 = idsurf[j]; // maps ghosts back to local

          // Skip if already added surface, can happen due to minimum image
          if (jadded.find(j2) != jadded.end()) continue;

          delx = xtmp - xsurf[j2][0];
          dely = ytmp - xsurf[j2][1];
          delz = ztmp - xsurf[j2][2];
          domain->minimum_image(FLERR, delx, dely, delz);
          rsq = delx * delx + dely * dely + delz * delz;
          radsum = radi + radsurf[j2] + skin;
          cutsq = radsum * radsum;

          if (rsq <= cutsq) {
            // since lines/tris flat, separately check distance along normal
            //   should reduce a lot of potential neighbors
            if (dimension == 2) {
              dot = lines[j2].norm[0] * delx + lines[j2].norm[1] * dely + lines[j2].  norm[2] * delz;
              MathExtra::scale3(dot, lines[j2].norm, dx_proj);
            } else {
              dot = tris[j2].norm[0] * delx + tris[j2].norm[1] * dely + tris[j2].norm  [2] * delz;
              MathExtra::scale3(dot, tris[j2].norm, dx_proj);
            }

            rsq = MathExtra::lensq3(dx_proj);
            if (rsq < cutsq_proj) {
              // Note: saves index of surf (like its tag) so will not work
              //       with default FixNeighHist methods that grab partner tags
              neighptr[n] = idsurf[j2];
              jadded.insert(idsurf[j2]);
              n++;
            }
          }
        }
      }

      ilist[inum++] = i;
      firstneigh[i] = neighptr;
      numneigh[i] = n;
      ipage->vgot(n);
      if (ipage->status())
        error->one(FLERR, Error::NOLASTLINE,
                   "Fix surface/global neighbor list overflow, boost neigh_modify one");
    }

    list->inum = inum;
  } else {
    // Nsq algorithm
    for (i = 0; i < nlocal; i++) {
      n = 0;
      neighptr = ipage->vget();
      if (use_history) {
        nn = 0;
        valueptr = dpage_atom->vget();
      }

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      radi = radius[i];
      cutsq_proj = radi + skin;
      cutsq_proj *= cutsq_proj;

      // for now, loop over all surfs
      for (j = 0; j < nsurf; j++) {
        delx = xtmp - xsurf[j][0];
        dely = ytmp - xsurf[j][1];
        delz = ztmp - xsurf[j][2];
        domain->minimum_image(FLERR, delx, dely, delz);
        rsq = delx * delx + dely * dely + delz * delz;
        radsum = radi + radsurf[j] + skin;
        cutsq = radsum * radsum;

        if (rsq <= cutsq) {
          // since lines/tris flat, separately check distance along normal
          //   should reduce a lot of potential neighbors
          if (dimension == 2) {
            dot = lines[j].norm[0] * delx + lines[j].norm[1] * dely + lines[j].norm[2] * delz;
            MathExtra::scale3(dot, lines[j].norm, dx_proj);
          } else {
            dot = tris[j].norm[0] * delx + tris[j].norm[1] * dely + tris[j].norm[2] * delz;
            MathExtra::scale3(dot, tris[j].norm, dx_proj);
          }

          rsq = MathExtra::lensq3(dx_proj);

          if (rsq < cutsq_proj) {
            // Note: saves index of surf (like its tag) so will not work
            //       with default FixNeighHist methods that grab partner tags
            neighptr[n] = j;
            n++;
          }
        }
      }

      ilist[inum++] = i;
      firstneigh[i] = neighptr;
      numneigh[i] = n;
      ipage->vgot(n);
      if (ipage->status())
        error->one(FLERR, Error::NOLASTLINE,
                   "Fix surface/global neighbor list overflow, boost neigh_modify one");
    }

    list->inum = inum;

  }
}

/* ----------------------------------------------------------------------
   compute particle/surface interactions
   impart force and torque to spherical particles
------------------------------------------------------------------------- */

void FixSurfaceGlobal::post_force(int /*vflag*/)
{
  int i, j, k, a, m, ii, jj, inum, jnum, jflag;
  int itype, jtype, external_flag, priority;
  double xtmp, ytmp, ztmp, radi, delx, dely, delz, meff;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firstflag, touch_flag;
  double rsq, rsq_com, rmag, radsum, max_overlap, dot;
  double x_min_image[3], norm[3], dr[3], contact[3], ds[3], xc[3], vc[3], omegac[3];
  double *forces, *torquesi, *history, *allhistory, **firsthistory;

  int jjtmp, nsidej;
  std::vector<int> composite_surfs;
  std::unordered_set<int> processed_contacts;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body", tmp);
    auto *mass_body = (double *) fix_rigid->extract("masstotal", tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid, nmax, "surface/global:mass_rigid");
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
  int *type = atom->type;
  int *mask = atom->mask;

  class GranularModel* model;
  for (auto n = 0; n < nmodel; n++) {
    model = models[n];
    model->history_update = 1;
    model->radj = 0.0;
    if (update->setupflag) model->history_update = 0;
    if (heat_flag) {
      if (tstr)
        Twall = input->variable->compute_equal(tvar);
      model->Tj = Twall;
    }
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  if (use_history) {
    firstflag = fix_history->firstflag;
    firsthistory = fix_history->firstvalue;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    itype = type[i];

    model->xi = x[i];
    model->radi = radius[i];
    model->vi = v[i];
    model->omegai = omega[i];

    // if I is part of rigid body, use body mass
    meff = rmass[i];
    if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    if (use_history) {
      touch = firstflag[i];
      allhistory = firsthistory[i];
    }

    contact_surfs.clear();

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      // Ignore special bonds
      j &= NEIGHMASK;

      delx = xtmp - xsurf[j][0];
      dely = ytmp - xsurf[j][1];
      delz = ztmp - xsurf[j][2];
      domain->minimum_image(FLERR, delx, dely, delz);

      rsq_com = delx * delx + dely * dely + delz * delz;

      // skip contact check if particle/surf are too far apart

      radsum = radi + radsurf[j];
      if (rsq_com > radsum * radsum) {
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // check for contact between particle and line/tri

      x_min_image[0] = delx + xsurf[j][0];
      x_min_image[1] = dely + xsurf[j][1];
      x_min_image[2] = delz + xsurf[j][2];

      if (dimension == 2) {

        // check for overlap of sphere and line segment
        // jflag = 0 for no overlap, 1 for interior line pt, -1/-2 for end pts
        // if no overlap, just continue
        // for overlap, also return:
        //   contact = nearest point on line to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr

        jflag = SurfExtra::
          overlap_sphere_line(x_min_image, radius[i],
                              points[lines[j].p1].x, points[lines[j].p2].x,
                              contact, dr, rsq);
      } else {

        // check for overlap of sphere and triangle
        // jflag = 0 for no overlap, 1 for interior line pt,
        //   -1/-2/-3 for 3 edges, -4/-5/-6 for 3 corner pts
        // if no overlap, just continue
        // for overlap, also returns:
        //   contact = nearest point on tri to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr

        jflag = SurfExtra::
          overlap_sphere_tri(x_min_image, radius[i],
                             points[tris[j].p1].x, points[tris[j].p2].x,
                             points[tris[j].p3].x, tris[j].norm,
                             contact, dr, rsq);
      }

      // unset non-touching neighbors

      if (!jflag) {
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // Find out if contact is on an external edge/corner
      external_flag = INTERNAL;
      if (dimension == 2) {
        MathExtra::copy3(lines[j].norm, norm);
        dot = MathExtra::dot3(norm, dr);
        jtype = lines[j].type;
        if (jflag == -1) external_flag = connect2d[j].external_pt[0];
        if (jflag == -2) external_flag = connect2d[j].external_pt[1];
      } else {
        MathExtra::copy3(tris[j].norm, norm);
        dot = MathExtra::dot3(norm, dr);
        jtype = tris[j].type;
        if (jflag == -1) external_flag = connect3d[j].external_edge[0];
        if (jflag == -2) external_flag = connect3d[j].external_edge[1];
        if (jflag == -3) external_flag = connect3d[j].external_edge[2];
        if (jflag == -4) external_flag = connect3d[j].external_cor[0];
        if (jflag == -5) external_flag = connect3d[j].external_cor[1];
        if (jflag == -6) external_flag = connect3d[j].external_cor[2];
      }

      // Store which side is in contact relative to normal vector
      if (dot >= 0) nsidej = SAME_SIDE;
      else nsidej = OPPOSITE_SIDE;

      rmag = sqrt(rsq);
      MathExtra::scale3(1.0 / rmag, dr, dr);

      priority = 2;
      if (jflag != 1) priority = 1;
      if (dimension == 3 && jflag < -3) priority = 0;

      struct ContactSurf mycontact;
      mycontact.index = j;
      mycontact.neigh_index = jj;
      mycontact.type = jtype;
      mycontact.flag = jflag;
      mycontact.external = external_flag;
      mycontact.nside = nsidej;
      mycontact.overlap = radi - rmag;
      mycontact.weight_contribution = 1.0;
      mycontact.convex_index = -1;
      mycontact.rmag = rmag;
      mycontact.rsq_com = rsq_com;
      mycontact.priority = priority;

      MathExtra::copy3(norm, mycontact.surf_norm);
      MathExtra::copy3(dr, mycontact.dr);
      MathExtra::copy3(dr, mycontact.dr_force);

      contact_surfs.push_back(mycontact);
    }


    if (contact_surfs.size() == 0)
      continue;

    // Sort contacts by overlap and create a map
    std::sort(contact_surfs.begin(), contact_surfs.end(), FixSurface::contact_presort);

    contacts_map.clear();
    for (auto n = 0; n < contact_surfs.size(); n++)
      contacts_map[contact_surfs[n].index] = n;

    // Initial walk to assign consistent sides of surfaces
    //   Won't guarantee will work for v. complex geometries (e.g. Mobius)

    if (dimension == 2) prewalk_connections2d();
    else prewalk_connections3d();

    // Given corrected surface norms, resort contacts
    std::sort(contact_surfs.begin(), contact_surfs.end(), FixSurface::contact_sort);

    for (auto n = 0; n < contact_surfs.size(); n++)
      contacts_map[contact_surfs[n].index] = n;

    processed_contacts.clear();
    for (auto n = 0; n < contact_surfs.size(); n++) {

      j = contact_surfs[n].index;
      if (processed_contacts.find(j) != processed_contacts.end()) continue;

      composite_surfs.clear();
      if (dimension == 2) {
        walk_connections2d(composite_surfs, processed_contacts);
        max_overlap = calculate_2d_forces(composite_surfs);
      } else {
        walk_connections3d(composite_surfs, processed_contacts);
        max_overlap = calculate_3d_forces(composite_surfs);
      }

      if (max_overlap < EPSILON)
        continue;

      // Calculate geometry of contact
      if (composite_surfs.size() > 1) {

        // Calculate overlap-weighted average normal vector
        MathExtra::zero3(dr);
        for (auto it = 0; it < composite_surfs.size(); it++) {
          m = composite_surfs[it];
          if (contact_surfs[m].overlap < EPSILON) continue;
          MathExtra::scaleadd3(contact_surfs[m].overlap * contact_surfs[m].weight_contribution, contact_surfs[m].dr_force, dr, dr);
        }
        MathExtra::norm3(dr);
        MathExtra::scale3(radi - max_overlap, dr);
      } else {
        MathExtra::scale3(radi - max_overlap, contact_surfs[n].dr_force, dr);
      }

      for (a = 0; a < 3; a++)
        xc[a] = x[i][a] - dr[a];

      MathExtra::zero3(vc);
      MathExtra::zero3(omegac);

      jtype = contact_surfs[n].type;
      model = types2model[itype][jtype];
      model->xi = x[i];
      model->radi = radi;
      model->vi = v[i];
      model->omegai = omega[i];
      if (heat_flag) model->Ti = temperature[i];
      model->meff = meff;

      // Correct velocity at contact point, extending from closest surf j
      ds[0] = xc[0] - xsurf[j][0];
      ds[1] = xc[1] - xsurf[j][1];
      ds[2] = xc[2] - xsurf[j][2];
      vc[0] = vsurf[j][0] + (omegasurf[j][1] * ds[2] - omegasurf[j][2] * ds[1]);
      vc[1] = vsurf[j][1] + (omegasurf[j][2] * ds[0] - omegasurf[j][0] * ds[2]);
      vc[2] = vsurf[j][2] + (omegasurf[j][0] * ds[1] - omegasurf[j][1] * ds[0]);

      model->xj = xc;
      model->vj = vc;
      model->omegaj = omegac;

      if (use_history) {
        jj = contact_surfs[n].neigh_index;
        model->touch = touch[jj];
      }

      // guaranteed in contact, but need to calculate intermediate variables
      touch_flag = model->check_contact();

      if (use_history) {
        // Check if another flat contact has a stored history
        if (touch[jj] != 1) {
          for (auto it = 0; it < composite_surfs.size(); it++) {
            m = composite_surfs[it];
            jjtmp = contact_surfs[m].neigh_index;
            if (touch[jjtmp] == 1)
              jj = jjtmp;
          }
        }

        touch[jj] = 1;
        history = &allhistory[size_history * jj];
        model->history = history;
      }

      model->calculate_forces();

      // Sychronize history across flat contacts
      //   can be arbitrary if not all connected flat surfaces are mutually flat
      //   e.g. a hair pin turn where surfs on either end of the 'U' are not flat
      if (use_history) {
        for (auto it = 0; it < composite_surfs.size(); it++) {
          m = composite_surfs[it];
          jjtmp = contact_surfs[m].neigh_index;
          if (jj != jjtmp) {
            touch[jjtmp] = 1;
            for (k = 0; k < size_history; k++)
              allhistory[size_history * jjtmp + k] = history[k];
          }
        }
      }

      forces = model->forces;
      torquesi = model->torquesi;

      add3(f[i], forces, f[i]);
      add3(torque[i], torquesi, torque[i]);
      if (heat_flag) heatflow[i] += model->dq;
    }
  }
}

/* ----------------------------------------------------------------------
   process fix_modify commands specific to fix surface/global
   move and type/region
------------------------------------------------------------------------- */

int FixSurfaceGlobal::modify_param(int narg, char **arg)
{
  // move keyword

  if (strcmp(arg[0],"move") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "fix surface/global modify move", error);

    int lo,hi;
    int *smols = new int[maxsurfmol+1];
    for (int i = 1; i <= maxsurfmol; i++) smols[i] = 0;

    auto fields = Tokenizer(arg[1], ",").as_vector();
    for (int ifield = 0; ifield < fields.size(); ifield++) {
      utils::bounds(FLERR, fields[ifield], 1, maxsurfmol, lo, hi, error);
      for (int i = lo; i <= hi; i++) smols[i] = 1;
    }

    if (strcmp(arg[2],"none") == 0) {
      int count = 0;
      for (int imol = 1; imol <= maxsurfmol; imol++) {
        if (mol2motion[imol] < 0) continue;
        int imotion = mol2motion[imol];
        motions[imotion].active = 0;
        mol2motion[imol] = -1;

        // set vsurf and omegasurf to zero for imol surfs

        int smol;
        for (int i = 0; i < nsurf; i++) {
          if (dimension == 2) smol = lines[i].mol;
          else smol = tris[i].mol;
          if (smol == imol) {
            vsurf[i][0] = vsurf[i][1] = vsurf[i][2] = 0.0;
            omegasurf[i][0] = omegasurf[i][1] = omegasurf[i][2] = 0.0;
            count++;
          }
        }
      }

      // error check if any surf now assigned to inactive motion
      // b/c specification of mols for "none" was incomplete

      int imotion;

      count = 0;
      if (dimension == 2) {
        for (int i = 0; i < nlines; i++) {
          imotion = mol2motion[lines[i].mol];
          if (imotion < 0) continue;
          if (!motions[imotion].active) count++;
        }
      } else {
        for (int i = 0; i < ntris; i++) {
          imotion = mol2motion[tris[i].mol];
          if (imotion < 0) continue;
          if (!motions[imotion].active) count++;
        }
      }

      if (count)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix_modify move none left {} surfs assigned inactive motion", count);

      // stats

      utils::logmesg(lmp,"Fix_modify move:\n");
      utils::logmesg(lmp,"  turned off motion for {} surfs\n",count);

      // reset anymove and anymove_variable
      // if no anymove, deallocate memory and turn off INITIAL_INTEGRATE

      anymove = 0;
      for (int i = 1; i <= maxsurfmol; i++)
        if (mol2motion[i] >= 0) anymove = 1;

      anymove_variable = 0;
      for (int i = 0; i < nmotion; i++)
        if (motions[i].active && motions[i].mstyle == VARIABLE)
          anymove_variable = 1;

      if (!anymove) {
        memory->destroy(points_lastneigh);
        memory->destroy(points_original);
        memory->destroy(xsurf_original);
        memory->destroy(pointmove);
        points_lastneigh = nullptr;
        points_original = nullptr;
        xsurf_original = nullptr;
        pointmove = nullptr;

        int ifix = modify->find_fix(id);
        modify->fmask[ifix] &= ~INITIAL_INTEGRATE;
        force_reneighbor = 0;
        next_reneighbor = -1;
      }

      delete [] smols;
      return 3;
    }

    // new motion operation
    // re-use an inactive motion or add a new motion to list

    int imotion;
    for (imotion = 0; imotion < nmotion; imotion++)
      if (!motions[imotion].active) break;

    if (imotion == nmotion) {
      if (nmotion == maxmotion) {
        maxmotion += DELTAMOTION;
        motions = (Motion *)
          memory->srealloc(motions,maxmotion*sizeof(Motion),"surface/global:motion");
      }
      nmotion++;
    }

    motions[imotion].active = 0;

    // use smols to set mol2motion

    for (int imol = 1; imol <= maxsurfmol; imol++)
      if (smols[imol]) mol2motion[imol] = imotion;

    // if first motion, allocate points and surf memory

    if (!anymove) {
      memory->create(points_lastneigh,npoints,3,"surface/global:points_lastneigh");
      memory->create(points_original,npoints,3,"surface/global:points_original");
      memory->create(xsurf_original,nsurf,3,"surface/global:xsurf_original");
      memory->create(pointmove,npoints,"surface/global:pointmove");
    }

    anymove = 1;

    // parse additional move style arguments

    int styleargs = modify_param_move(&motions[imotion],narg-2,&arg[2]);
    int mstyle = motions[imotion].mstyle;
    if (mstyle == VARIABLE) anymove_variable = 1;
    motions[imotion].time_origin = update->ntimestep;

    int ifix = modify->find_fix(id);
    modify->fmask[ifix] |= INITIAL_INTEGRATE;

    force_reneighbor = 1;
    next_reneighbor = -1;

    // intialize points and surfs in smols

    int imol,p1,p2,p3;
    double omega;
    double *runit;
    int count = 0;

    for (int i = 0; i < nsurf; i++) {
      if (dimension == 2) imol = lines[i].mol;
      else imol = tris[i].mol;
      if (!smols[imol]) continue;
      count++;

      if (dimension == 2) {
        p1 = lines[i].p1;
        p2 = lines[i].p2;
        points_lastneigh[p1][0] = points_original[p1][0] = points[p1].x[0];
        points_lastneigh[p1][1] = points_original[p1][1] = points[p1].x[1];
        points_lastneigh[p1][2] = points_original[p1][2] = points[p1].x[2];
        points_lastneigh[p2][0] = points_original[p2][0] = points[p2].x[0];
        points_lastneigh[p2][1] = points_original[p2][1] = points[p2].x[1];
        points_lastneigh[p2][2] = points_original[p2][2] = points[p2].x[2];
      } else {
        p1 = tris[i].p1;
        p2 = tris[i].p2;
        p3 = tris[i].p3;
        points_lastneigh[p1][0] = points_original[p1][0] = points[p1].x[0];
        points_lastneigh[p1][1] = points_original[p1][1] = points[p1].x[1];
        points_lastneigh[p1][2] = points_original[p1][2] = points[p1].x[2];
        points_lastneigh[p2][0] = points_original[p2][0] = points[p2].x[0];
        points_lastneigh[p2][1] = points_original[p2][1] = points[p2].x[1];
        points_lastneigh[p2][2] = points_original[p2][2] = points[p2].x[2];
        points_lastneigh[p3][0] = points_original[p3][0] = points[p3].x[0];
        points_lastneigh[p3][1] = points_original[p3][1] = points[p3].x[1];
        points_lastneigh[p3][2] = points_original[p3][2] = points[p3].x[2];
      }

      xsurf_original[i][0] = xsurf[i][0];
      xsurf_original[i][1] = xsurf[i][1];
      xsurf_original[i][2] = xsurf[i][2];

      if (mstyle == ROTATE || mstyle == TRANSROT) {
        omega = motions[imotion].omega;
        runit = motions[imotion].unit;
        omegasurf[i][0] = omega*runit[0];
        omegasurf[i][1] = omega*runit[1];
        omegasurf[i][2] = omega*runit[2];
      }
    }

    utils::logmesg(lmp,"Fix_modify move:\n");
    utils::logmesg(lmp,"  turned on motion for {} surfs\n", count);

    delete [] smols;
    return 2 + styleargs;
  }

  // mol/region keyword

  if (strcmp(arg[0],"mol/region") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "fix surface/global modify mol/region", error);

    int smol = utils::inumeric(FLERR,arg[1],false,lmp);
    if ((smol <= 0) || (smol > maxsurfmol))
      error->all(FLERR,"Invalid fix_modify mol/region surf mol");

    auto *region = domain->get_region_by_id(arg[2]);
    if (!region) error->all(FLERR,"Fix_modify mol/region region {} does not exist", arg[2]);

    int count = 0;
    if (dimension == 2) {
      for (int i = 0; i < nlines; i++)
        if (region->match(xsurf[i][0],xsurf[i][1],xsurf[i][2])) {
          lines[i].mol = smol;
          count++;
        }
    } else {
      for (int i = 0; i < ntris; i++)
        if (region->match(xsurf[i][0],xsurf[i][1],xsurf[i][2])) {
          tris[i].mol = smol;
          count++;
        }
    }

    utils::logmesg(lmp,"Fix_modify mol/region:\n");
    utils::logmesg(lmp,"  {} surfs assigned to mol {}\n",count,smol);

    return 3;
  }

  // type/region keyword

  if (strcmp(arg[0],"type/region") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "fix surface/global modify type/region", error);

    int stype = utils::inumeric(FLERR,arg[1],false,lmp);
    if ((stype <= 0) || (stype > maxsurftype))
      error->all(FLERR,"Invalid fix_modify type/region surf type");

    auto *region = domain->get_region_by_id(arg[2]);
    if (!region) error->all(FLERR,"Fix_modify type/region region {} does not exist", arg[2]);

    int count = 0;
    if (dimension == 2) {
      for (int i = 0; i < nlines; i++)
        if (region->match(xsurf[i][0],xsurf[i][1],xsurf[i][2])) {
          lines[i].type = stype;
          count++;
        }
    } else {
      for (int i = 0; i < ntris; i++)
        if (region->match(xsurf[i][0],xsurf[i][1],xsurf[i][2])) {
          tris[i].type = stype;
          count++;
        }
    }

    utils::logmesg(lmp,"Fix_modify type/region:\n");
    utils::logmesg(lmp,"  {} surfs assigned to type {}\n",count,stype);

    return 3;
  }

  // keyword not recognized

  return 0;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceGlobal::modify_param_move(Motion *motion, int narg, char **arg)
{
  if (strcmp(arg[0],"linear") == 0) {
    if (narg < 4) utils::missing_cmd_args(FLERR, "fix surface/global modify move linear", error);
    motion->mstyle = LINEAR;

    if (strcmp(arg[1], "NULL") == 0) motion->vxflag = 0;
    else {
      motion->vxflag = 1;
      motion->vx = utils::numeric(FLERR, arg[4], false, lmp);
    }
    if (strcmp(arg[2], "NULL") == 0) motion->vyflag = 0;
    else {
      motion->vyflag = 1;
      motion->vy = utils::numeric(FLERR, arg[2], false, lmp);
    }
    if (strcmp(arg[3], "NULL") == 0) motion->vzflag = 0;
    else {
      motion->vzflag = 1;
      motion->vz = utils::numeric(FLERR, arg[3], false, lmp);
    }

    if (dimension == 2)
      if (motion->vzflag && (motion->vz != 0.0))
        error->all(FLERR,"Fix_modify move cannot set linear z motion for 2d problem");

    return 4;
  }

  if (strcmp(arg[0],"wiggle") == 0) {
    if (narg < 5) utils::missing_cmd_args(FLERR, "fix surface/global modify move wiggle", error);
    motion->mstyle = WIGGLE;

    if (strcmp(arg[1], "NULL") == 0) motion->axflag = 0;
    else {
      motion->axflag = 1;
      motion->ax = utils::numeric(FLERR, arg[4], false, lmp);
    }
    if (strcmp(arg[2], "NULL") == 0) motion->ayflag = 0;
    else {
      motion->ayflag = 1;
      motion->ay = utils::numeric(FLERR, arg[2], false, lmp);
    }
    if (strcmp(arg[3], "NULL") == 0) motion->azflag = 0;
    else {
      motion->azflag = 1;
      motion->az = utils::numeric(FLERR, arg[3], false, lmp);
    }

    if (dimension == 2)
      if (motion->azflag && (motion->az != 0.0))
        error->all(FLERR,"Fix_modify move cannot set wiggle z motion for 2d problem");

    motion->period = utils::numeric(FLERR, arg[7], false, lmp);
    if (motion->period <= 0.0)
      error->all(FLERR, "Illegal fix_modify move period {}", motion->period);
    motion->omega = MY_2PI / motion->period;

    return 5;
  }

  if (strcmp(arg[0],"rotate") == 0) {
    if (narg < 8) utils::missing_cmd_args(FLERR, "fix surface/global modify move rotate", error);
    motion->mstyle = ROTATE;

    motion->point[0] = utils::numeric(FLERR,arg[1],false,lmp);
    motion->point[1] = utils::numeric(FLERR,arg[2],false,lmp);
    motion->point[2] = utils::numeric(FLERR,arg[3],false,lmp);

    motion->axis[0] = utils::numeric(FLERR,arg[4],false,lmp);
    motion->axis[1] = utils::numeric(FLERR,arg[5],false,lmp);
    motion->axis[2] = utils::numeric(FLERR,arg[6],false,lmp);

    if (dimension == 2)
      if (motion->axis[0] != 0.0 || motion->axis[1] != 0.0)
        error->all(FLERR,"Fix_modify move cannot rotate around non z-axis for 2d problem");

    motion->period = utils::numeric(FLERR,arg[7],false,lmp);
    if (motion->period <= 0.0)
      error->all(FLERR,"Illegal fix_modify move period {}", motion->period);

    motion->omega = MY_2PI / motion->period;

    // runit = unit vector along rotation axis

    double len = MathExtra::len3(motion->axis);
    if (len == 0.0)
      error->all(FLERR,"Fix_modify move zero length rotation vector");
    MathExtra::normalize3(motion->axis,motion->unit);

    return 8;
  }

  if (strcmp(arg[0],"transrot") == 0) {
    if (narg < 11) utils::missing_cmd_args(FLERR, "fix surface/global modify move transrot", error);
    motion->mstyle = TRANSROT;

    error->all(FLERR, "Fix_modify move transrot not yet supported for fix surface/global");

    motion->vxflag = motion->vyflag = motion->vzflag = 1;
    motion->vx = utils::numeric(FLERR, arg[1], false, lmp);
    motion->vy = utils::numeric(FLERR, arg[2], false, lmp);
    motion->vz = utils::numeric(FLERR, arg[3], false, lmp);

    motion->point[0] = utils::numeric(FLERR,arg[4],false,lmp);
    motion->point[1] = utils::numeric(FLERR,arg[5],false,lmp);
    motion->point[2] = utils::numeric(FLERR,arg[6],false,lmp);

    motion->axis[0] = utils::numeric(FLERR,arg[7],false,lmp);
    motion->axis[1] = utils::numeric(FLERR,arg[8],false,lmp);
    motion->axis[2] = utils::numeric(FLERR,arg[9],false,lmp);

    if (dimension == 2) {
      if (motion->vzflag && (motion->vz != 0.0))
        error->all(FLERR,"Fix_modify move cannot set linear z motion for 2d problem");
      if (motion->axis[0] != 0.0 || motion->axis[1] != 0.0)
        error->all(FLERR,"Fix_modify move cannot rotate around non z-axis for 2d problem");
    }

    motion->period = utils::numeric(FLERR,arg[10],false,lmp);
    if (motion->period <= 0.0)
      error->all(FLERR,"Illegal fix_modify move command, period {} must be positive", motion->period);

    motion->omega = MY_2PI / motion->period;

    // runit = unit vector along rotation axis

    double len = MathExtra::len3(motion->axis);
    if (len == 0.0)
      error->all(FLERR,"Fix_modify move zero length rotation vector");
    MathExtra::normalize3(motion->axis,motion->unit);

    return 11;
  }

  if (strcmp(arg[0],"variable") == 0) {
    if (narg < 7) utils::missing_cmd_args(FLERR, "fix surface/global modify move variable", error);
    motion->mstyle = VARIABLE;

    if (strcmp(arg[1], "NULL") == 0)
      motion->xvarstr = nullptr;
    else if (utils::strmatch(arg[1], "^v_")) {
      motion->xvarstr = utils::strdup(arg[1] + 2);
    } else
      error->all(FLERR, "Illegal fix_modify move command");
    if (strcmp(arg[2], "NULL") == 0)
      motion->yvarstr = nullptr;
    else if (utils::strmatch(arg[2], "^v_")) {
      motion->yvarstr = utils::strdup(arg[2] + 2);
    } else
      error->all(FLERR, "Illegal fix_modify move command");
    if (strcmp(arg[3], "NULL") == 0)
      motion->zvarstr = nullptr;
    else if (utils::strmatch(arg[3], "^v_")) {
      motion->zvarstr = utils::strdup(arg[3] + 2);
    } else
      error->all(FLERR, "Illegal fix_modify move command");


    if (strcmp(arg[4], "NULL") == 0)
      motion->vxvarstr = nullptr;
    else if (utils::strmatch(arg[4], "^v_")) {
      motion->vxvarstr = utils::strdup(arg[4] + 2);
    } else
      error->all(FLERR, "Illegal fix_modify move command");
    if (strcmp(arg[5], "NULL") == 0)
      motion->vyvarstr = nullptr;
    else if (utils::strmatch(arg[5], "^v_")) {
      motion->vyvarstr = utils::strdup(arg[5] + 2);
    } else
      error->all(FLERR, "Illegal fix_modify move command");
    if (strcmp(arg[6], "NULL") == 0)
      motion->vzvarstr = nullptr;
    else if (utils::strmatch(arg[6], "^v_")) {
      motion->vzvarstr = utils::strdup(arg[6] + 2);
    } else
      error->all(FLERR, "Illegal fix_modify move command");

    if (dimension == 2 && (motion->zvarstr || motion->vzvarstr))
      error->all(FLERR, "Fix_modify move cannot define z or vz variable for 2d problem");

    return 7;
  }

  error->all(FLERR,"Fix_modify move style {} not recognized", arg[0]);

  return 0;
}

/* ----------------------------------------------------------------------
   memory usage per-surf data and neighbor list
------------------------------------------------------------------------- */

double FixSurfaceGlobal::memory_usage()
{
  double bytes = 0.0;

  // points, lines, tris and connect2d/3d

  bytes += npoints*sizeof(Point);
  if (dimension == 2) {
    bytes += nlines*sizeof(Line);
    bytes += nlines*sizeof(Connect2d);
  } else if (dimension == 3) {
    bytes += ntris*sizeof(Tri);
    bytes += ntris*sizeof(Connect3d);
  }

  bytes += memory->usage(xsurf,nsurf,3);
  bytes += memory->usage(vsurf,nsurf,3);
  bytes += memory->usage(omegasurf,nsurf,3);
  bytes += memory->usage(radsurf,nsurf);

  if (anymove) {
    bytes += memory->usage(points_lastneigh,npoints,3);
    bytes += memory->usage(points_original,npoints,3);
    bytes += memory->usage(xsurf_original,nsurf,3);
    bytes += memory->usage(pointmove,nsurf);
  }

  if (mass_rigid) bytes += memory->usage(mass_rigid,atom->nmax);

  if (imax) {
    bytes += memory->usage(imflag,nsurf);
    if (dimension == 2) bytes += memory->usage(imdata,nsurf,7);
    else if (dimension == 3) bytes += memory->usage(imdata,nsurf,10);
  }

  // ragged connectivity arrays

  if (dimension == 2) {
    int np = 0;
    for (int i = 0; i < nlines; i++) np += connect2d[i].np1 + connect2d[i].np2;
    bytes += 4*np * sizeof(int);
  } else if (dimension == 3) {
    int ne = 0;
    for (int i = 0; i < ntris; i++) ne += connect3d[i].ne1 + connect3d[i].ne2 + connect3d[i].ne3;
    bytes += 4*ne * sizeof(int);
    int nc = 0;
    for (int i = 0; i < ntris; i++) nc += connect3d[i].nc1 + connect3d[i].nc2 + connect3d[i].nc3;
    bytes += 2*nc * sizeof(int);
  }

  // neighbor list

  bytes += list->memory_usage();

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

      imflag[i] = Graphics::LINE;
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

      imflag[i] = Graphics::TRI;
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
// initializiation of surfs
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   error checks on lines
   no zero-length lines and no duplicate lines
------------------------------------------------------------------------- */

void FixSurfaceGlobal::check2d()
{
  // check for zero length lines
  // use coords of points, not indices
  // since (p1,p2) can be zero-length even if p1 != p2

  double *pt1,*pt2;

  int flag = 0;
  for (int i = 0; i < nlines; i++) {
    pt1 = points[lines[i].p1].x;
    pt2 = points[lines[i].p2].x;
    if (pt1[0] == pt2[0] && pt1[1] == pt2[1]) flag++;
  }

  if (flag)
    error->all(FLERR, "Fix surface/global defines {} zero-length lines",flag);

  // check for duplicate lines
  // (p1,p2) is duplicate of any line with same 2 endpoints
  // hash = map of lines
  //   key = <p1,p2> indices of 2 points
  //   value = 1 if already appeared

  int p1,p2;
  std::map<std::tuple<int,int>,int> hash;

  flag = 0;
  for (int i = 0; i < nlines; i++) {
    p1 = lines[i].p1;
    p2 = lines[i].p2;

    auto key1 = std::make_tuple(p1,p2);
    auto key2 = std::make_tuple(p2,p1);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) hash[key1] = 1;
    else flag++;
  }

  if (flag)
    error->all(FLERR, "Fix surface/global defines {} duplicate lines",flag);
}


/* ----------------------------------------------------------------------
   error checks on tris
   no tris with a zero-length edge
   no duplicate tris
------------------------------------------------------------------------- */

void FixSurfaceGlobal::check3d()
{
  // check for zero length tri edges
  // use coords of points, not indices
  // since (p1,p2) can be zero-length even if p1 != p2

  double *pt1,*pt2,*pt3;

  int flag = 0;
  for (int i = 0; i < ntris; i++) {
    pt1 = points[tris[i].p1].x;
    pt2 = points[tris[i].p2].x;
    pt3 = points[tris[i].p3].x;

    if (pt1[0] == pt2[0] && pt1[1] == pt2[1] && pt1[2] == pt2[2]) flag++;
    if (pt2[0] == pt3[0] && pt2[1] == pt3[1] && pt2[2] == pt3[2]) flag++;
    if (pt3[0] == pt1[0] && pt3[1] == pt1[1] && pt3[2] == pt1[2]) flag++;
  }

  if (flag)
    error->all(FLERR, "Fix surface/global defines {} zero-length triangle edges",flag);

  // check for duplicate tris
  // (p1,p2,p3) of any tri with same 3 corner points
  // hash = map of tris
  //   key = <p1,p2,p3> indices of 3 points
  //   value = 1 if already appeared

  int p1,p2,p3;
  std::map<std::tuple<int,int,int>,int> hash;

  flag = 0;
  for (int i = 0; i < ntris; i++) {
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;

    auto key1 = std::make_tuple(p1,p2,p3);
    auto key2 = std::make_tuple(p1,p3,p2);
    auto key3 = std::make_tuple(p2,p1,p3);
    auto key4 = std::make_tuple(p2,p3,p1);
    auto key5 = std::make_tuple(p3,p1,p2);
    auto key6 = std::make_tuple(p3,p2,p1);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end() &&
        hash.find(key3) == hash.end() && hash.find(key4) == hash.end() &&
        hash.find(key5) == hash.end() && hash.find(key6) == hash.end()) hash[key1] = 1;
    else flag++;

  }

  if (flag)
    error->all(FLERR, "Fix surface/global defines {} duplicate triangles",flag);
}

/* ----------------------------------------------------------------------
   issue warning if any connected lines/tris have different molIDs
   molIDs are not used now, but may be in future
   an inter-connected set of surfs should have a single molID
------------------------------------------------------------------------- */

void FixSurfaceGlobal::check_molecules()
{
  int i,m,imol,flag;

  if (dimension == 2) {
    int *neigh_p1, *neigh_p2;
    flag = 0;
    for (i = 0; i < nlines; i++) {
      imol = lines[i].mol;
      neigh_p1 = connect2d[i].neigh_p1;
      neigh_p2 = connect2d[i].neigh_p2;
      for (m = 0; m < connect2d[i].np1; m++)
        if (imol != lines[neigh_p1[m]].mol) flag++;
      for (m = 0; m < connect2d[i].np2; m++)
        if (imol != lines[neigh_p2[m]].mol) flag++;
    }
    flag /= 2;
    if (flag && comm->me == 0)
      error->warning(FLERR,
                     "Fix surface/global endpoint connections between "
                     "different molecule IDs =",flag);

  } else {
    int *neigh_e1, *neigh_e2, *neigh_e3;
    int *neigh_c1, *neigh_c2, *neigh_c3;
    flag = 0;
    for (i = 0; i < ntris; i++) {
      imol = tris[i].mol;
      neigh_e1 = connect3d[i].neigh_e1;
      neigh_e2 = connect3d[i].neigh_e2;
      neigh_e3 = connect3d[i].neigh_e3;
      for (m = 0; m < connect3d[i].ne1; m++)
        if (imol != tris[neigh_e1[m]].mol) flag++;
      for (m = 0; m < connect3d[i].ne2; m++)
        if (imol != tris[neigh_e2[m]].mol) flag++;
      for (m = 0; m < connect3d[i].ne3; m++)
        if (imol != tris[neigh_e3[m]].mol) flag++;
      neigh_c1 = connect3d[i].neigh_c1;
      neigh_c2 = connect3d[i].neigh_c2;
      neigh_c3 = connect3d[i].neigh_c3;
      for (m = 0; m < connect3d[i].nc1; m++)
        if (imol != tris[neigh_c1[m]].mol) flag++;
      for (m = 0; m < connect3d[i].nc2; m++)
        if (imol != tris[neigh_c2[m]].mol) flag++;
      for (m = 0; m < connect3d[i].nc3; m++)
        if (imol != tris[neigh_c3[m]].mol) flag++;
    }
    flag /= 2;
    if (flag && comm->me == 0)
      error->warning(FLERR,
                     "Fix surface/global edge/corner connections between "
                     "different molecule IDs =",flag);
  }
}

/* ----------------------------------------------------------------------
   complete initialization of Connect2d info for all lines
------------------------------------------------------------------------- */

void FixSurfaceGlobal::connectivity2d_complete()
{
  // allocate ragged arrays for other Connect2d fields
  // p12_counts = # of lines connecting to endpoints p12 of each line

  int *p1_counts,*p2_counts;
  memory->create(p1_counts,nlines,"surface/global:p1_counts");
  memory->create(p2_counts,nlines,"surface/global:p2_counts");

  for (int i = 0; i < nlines; i++) {
    p1_counts[i] = connect2d[i].np1;
    p2_counts[i] = connect2d[i].np2;
  }

  memory->create_ragged(pwhich_p1,nlines,p1_counts,"surface/global:pwhich_p1");
  memory->create_ragged(pwhich_p2,nlines,p2_counts,"surface/global:pwhich_p2");
  memory->create_ragged(nside_p1,nlines,p1_counts,"surface/global:nside_p1");
  memory->create_ragged(nside_p2,nlines,p2_counts,"surface/global:nside_p2");
  memory->create_ragged(aflag_p1,nlines,p1_counts,"surface/global:aflag_p1");
  memory->create_ragged(aflag_p2,nlines,p2_counts,"surface/global:aflag_p2");
  memory->create_ragged(fflag_p1,nlines,p1_counts,"surface/global:fflag_p1");
  memory->create_ragged(fflag_p2,nlines,p2_counts,"surface/global:fflag_p2");

  memory->destroy(p1_counts);
  memory->destroy(p2_counts);

  // set connect2d vector ptrs to rows of corresponding ragged arrays

  for (int i = 0; i < nlines; i++) {
    if (connect2d[i].np1 == 0) {
      connect2d[i].pwhich_p1 = nullptr;
      connect2d[i].nside_p1 = nullptr;
      connect2d[i].aflag_p1 = nullptr;
      connect2d[i].fflag_p1 = nullptr;
    } else {
      connect2d[i].pwhich_p1 = pwhich_p1[i];
      connect2d[i].nside_p1 = nside_p1[i];
      connect2d[i].aflag_p1 = aflag_p1[i];
      connect2d[i].fflag_p1 = fflag_p1[i];
    }

    if (connect2d[i].np2 == 0) {
      connect2d[i].pwhich_p2 = nullptr;
      connect2d[i].nside_p2 = nullptr;
      connect2d[i].aflag_p2 = nullptr;
      connect2d[i].fflag_p2 = nullptr;
    } else {
      connect2d[i].pwhich_p2 = pwhich_p2[i];
      connect2d[i].nside_p2 = nside_p2[i];
      connect2d[i].aflag_p2 = aflag_p2[i];
      connect2d[i].fflag_p2 = fflag_p2[i];
    }
  }

  // set connect2d pwhich/nside/aflag for each end point of each line
  // see fsg.h file for an explanation of each vector in Connect2d
  // aflag is based on dot and cross product of 2 connected line normals
  //   cross product is either along +z or -z direction

  int j,m;

  for (int i = 0; i < nlines; i++) {
    for (m = 0; m < connect2d[i].np1; m++) {
      j = connect2d[i].neigh_p1[m];

      int jwhich = -1;
      if (lines[i].p1 == lines[j].p1) jwhich = 0;
      else if (lines[i].p1 == lines[j].p2) jwhich = 1;
      point_connection2d(lines[i].norm,lines[j].norm,0,jwhich,flatthresh,
                         connect2d[i].fflag_p1[m],connect2d[i].nside_p1[m],
                         connect2d[i].aflag_p1[m]);
      connect2d[i].pwhich_p1[m] = jwhich;
    }

    for (m = 0; m < connect2d[i].np2; m++) {
      j = connect2d[i].neigh_p2[m];

      int jwhich = -1;
      if (lines[i].p2 == lines[j].p1) jwhich = 0;
      else if (lines[i].p2 == lines[j].p2) jwhich = 1;
      point_connection2d(lines[i].norm,lines[j].norm,1,jwhich,flatthresh,
                         connect2d[i].fflag_p2[m],connect2d[i].nside_p2[m],
                         connect2d[i].aflag_p2[m]);
      connect2d[i].pwhich_p2[m] = jwhich;
    }
  }

  // deallocate p12_counts

  memory->destroy(p1_counts);
  memory->destroy(p2_counts);

  // determine whether pts are external based on connectivity

  for (int i = 0; i < nsurf; i++) {
    connect2d[i].external_pt[0] = INTERNAL;
    connect2d[i].external_pt[1] = INTERNAL;
  }

  for (int i = 0; i < nsurf; i++) {
    // external if there's a nonflat connection
    for (int n = 0; n < connect2d[i].np1; n++)
      if (connect2d[i].fflag_p1[n] != FLAT)
        connect2d[i].external_pt[0] = EXTERNAL;
    for (int n = 0; n < connect2d[i].np2; n++)
      if (connect2d[i].fflag_p2[n] != FLAT)
        connect2d[i].external_pt[1] = EXTERNAL;

    // or if unconnected on border
    if (connect2d[i].np1 == 0)
      connect2d[i].external_pt[0] = UNCONNECTED;
    if (connect2d[i].np2 == 0)
      connect2d[i].external_pt[1] = UNCONNECTED;
  }
}

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for all triangles
   creates connect3d data structs
------------------------------------------------------------------------- */

void FixSurfaceGlobal::connectivity3d_complete()
{
  // allocate ragged arrays for other Connect3d edge fields
  // e123_counts = # of edges connecting to edges 123 of each tri

  int *e1_counts,*e2_counts,*e3_counts;
  memory->create(e1_counts,ntris,"surface/global:e1_counts");
  memory->create(e2_counts,ntris,"surface/global:e2_counts");
  memory->create(e3_counts,ntris,"surface/global:e3_counts");

  for (int i = 0; i < ntris; i++) {
    e1_counts[i] = connect3d[i].ne1;
    e2_counts[i] = connect3d[i].ne2;
    e3_counts[i] = connect3d[i].ne3;
  }

  memory->create_ragged(ewhich_e1,ntris,e1_counts,"surface:ewhich_e1");
  memory->create_ragged(ewhich_e2,ntris,e2_counts,"surface:ewhich_e2");
  memory->create_ragged(ewhich_e3,ntris,e3_counts,"surface:ewhich_e3");
  memory->create_ragged(nside_e1,ntris,e1_counts,"surface:nside_e1");
  memory->create_ragged(nside_e2,ntris,e2_counts,"surface:nside_e2");
  memory->create_ragged(nside_e3,ntris,e3_counts,"surface:nside_e3");
  memory->create_ragged(aflag_e1,ntris,e1_counts,"surface:aflag_e1");
  memory->create_ragged(aflag_e2,ntris,e2_counts,"surface:aflag_e2");
  memory->create_ragged(aflag_e3,ntris,e3_counts,"surface:aflag_e3");
  memory->create_ragged(fflag_e1,ntris,e1_counts,"surface:fflag_e1");
  memory->create_ragged(fflag_e2,ntris,e2_counts,"surface:fflag_e2");
  memory->create_ragged(fflag_e3,ntris,e3_counts,"surface:fflag_e3");

  memory->destroy(e1_counts);
  memory->destroy(e2_counts);
  memory->destroy(e3_counts);

  // set connect3d vector ptrs to rows of corresponding ragged arrays

  for (int i = 0; i < ntris; i++) {
    if (connect3d[i].ne1 == 0) {
      connect3d[i].ewhich_e1 = nullptr;
      connect3d[i].nside_e1 = nullptr;
      connect3d[i].aflag_e1 = nullptr;
      connect3d[i].fflag_e1 = nullptr;
    } else {
      connect3d[i].ewhich_e1 = ewhich_e1[i];
      connect3d[i].nside_e1 = nside_e1[i];
      connect3d[i].aflag_e1 = aflag_e1[i];
      connect3d[i].fflag_e1 = fflag_e1[i];
    }

    if (connect3d[i].ne2 == 0) {
      connect3d[i].ewhich_e2 = nullptr;
      connect3d[i].nside_e2 = nullptr;
      connect3d[i].aflag_e2 = nullptr;
      connect3d[i].fflag_e2 = nullptr;
    } else {
      connect3d[i].ewhich_e2 = ewhich_e2[i];
      connect3d[i].nside_e2 = nside_e2[i];
      connect3d[i].aflag_e2 = aflag_e2[i];
      connect3d[i].fflag_e2 = fflag_e2[i];
    }

    if (connect3d[i].ne3 == 0) {
      connect3d[i].ewhich_e3 = nullptr;
      connect3d[i].nside_e3 = nullptr;
      connect3d[i].aflag_e3 = nullptr;
      connect3d[i].fflag_e3 = nullptr;
    } else {
      connect3d[i].ewhich_e3 = ewhich_e3[i];
      connect3d[i].nside_e3 = nside_e3[i];
      connect3d[i].aflag_e3 = aflag_e3[i];
      connect3d[i].fflag_e3 = fflag_e3[i];
    }
  }

  // set connect3d edge ewhich/nside/aflag for each edge of each tri
  // see fsg.h file for an explanation of each edge vector in Connect3d
  // aflag is based on dot and cross product of 2 connected tri normals
  //   cross product is either along itri edge or in opposite dir

  int jpfirst,jpsecond;
  double iedge[3];

  int j,m;

  for (int i = 0; i < ntris; i++) {
    for (m = 0; m < connect3d[i].ne1; m++) {
      j = connect3d[i].neigh_e1[m];

      jpfirst = jpsecond = -1;
      if (tris[i].p1 == tris[j].p1) jpfirst = 1;
      else if (tris[i].p1 == tris[j].p2) jpfirst = 2;
      else if (tris[i].p1 == tris[j].p3) jpfirst = 3;

      if (tris[i].p2 == tris[j].p1) jpsecond = 1;
      else if (tris[i].p2 == tris[j].p2) jpsecond = 2;
      else if (tris[i].p2 == tris[j].p3) jpsecond = 3;

      MathExtra::sub3(points[tris[i].p2].x,points[tris[i].p1].x,iedge);
      edge_connection3d(tris[i].norm,tris[j].norm,iedge,jpfirst,jpsecond,flatthresh,
                        connect3d[i].fflag_e1[m],connect3d[i].ewhich_e1[m],
                        connect3d[i].nside_e1[m],connect3d[i].aflag_e1[m]);
    }

    for (m = 0; m < connect3d[i].ne2; m++) {
      j = connect3d[i].neigh_e2[m];

      jpfirst = jpsecond = -1;
      if (tris[i].p2 == tris[j].p1) jpfirst = 1;
      else if (tris[i].p2 == tris[j].p2) jpfirst = 2;
      else if (tris[i].p2 == tris[j].p3) jpfirst = 3;

      if (tris[i].p3 == tris[j].p1) jpsecond = 1;
      else if (tris[i].p3 == tris[j].p2) jpsecond = 2;
      else if (tris[i].p3 == tris[j].p3) jpsecond = 3;

      MathExtra::sub3(points[tris[i].p3].x,points[tris[i].p2].x,iedge);
      edge_connection3d(tris[i].norm,tris[j].norm,iedge,jpfirst,jpsecond,flatthresh,
                        connect3d[i].fflag_e2[m],connect3d[i].ewhich_e2[m],
                        connect3d[i].nside_e2[m],connect3d[i].aflag_e2[m]);
    }

    for (m = 0; m < connect3d[i].ne3; m++) {
      j = connect3d[i].neigh_e3[m];

      jpfirst = jpsecond = -1;
      if (tris[i].p3 == tris[j].p1) jpfirst = 1;
      else if (tris[i].p3 == tris[j].p2) jpfirst = 2;
      else if (tris[i].p3 == tris[j].p3) jpfirst = 3;

      if (tris[i].p1 == tris[j].p1) jpsecond = 1;
      else if (tris[i].p1 == tris[j].p2) jpsecond = 2;
      else if (tris[i].p1 == tris[j].p3) jpsecond = 3;

      MathExtra::sub3(points[tris[i].p1].x,points[tris[i].p3].x,iedge);
      edge_connection3d(tris[i].norm,tris[j].norm,iedge,jpfirst,jpsecond,flatthresh,
                        connect3d[i].fflag_e3[m],connect3d[i].ewhich_e3[m],
                        connect3d[i].nside_e3[m],connect3d[i].aflag_e3[m]);
    }
  }

  // allocate ragged arrays for other Connect3d corner fields
  // c123_counts = # of tris connecting to corners 123 of each tri

  int *c1_counts,*c2_counts,*c3_counts;
  memory->create(c1_counts,ntris,"surface/global:e1_counts");
  memory->create(c2_counts,ntris,"surface/global:e2_counts");
  memory->create(c3_counts,ntris,"surface/global:e3_counts");

  for (int i = 0; i < ntris; i++) {
    c1_counts[i] = connect3d[i].nc1;
    c2_counts[i] = connect3d[i].nc2;
    c3_counts[i] = connect3d[i].nc3;
  }

  memory->create_ragged(cwhich_c1,ntris,c1_counts,"surface:cwhich_c1");
  memory->create_ragged(cwhich_c2,ntris,c2_counts,"surface:cwhich_c2");
  memory->create_ragged(cwhich_c3,ntris,c3_counts,"surface:cwhich_c3");
  memory->create_ragged(nside_c1,ntris,c1_counts,"surface:nside_c1");
  memory->create_ragged(nside_c2,ntris,c2_counts,"surface:nside_c2");
  memory->create_ragged(nside_c3,ntris,c3_counts,"surface:nside_c3");
  memory->create_ragged(fflag_c1,ntris,c1_counts,"surface:fflag_c1");
  memory->create_ragged(fflag_c2,ntris,c2_counts,"surface:fflag_c2");
  memory->create_ragged(fflag_c3,ntris,c3_counts,"surface:fflag_c3");

  memory->destroy(c1_counts);
  memory->destroy(c2_counts);
  memory->destroy(c3_counts);

  // set connect3d vector ptrs to rows of corresponding ragged arrays

  for (int i = 0; i < ntris; i++) {
    if (connect3d[i].nc1 == 0) {
      connect3d[i].cwhich_c1 = nullptr;
      connect3d[i].nside_c1 = nullptr;
      connect3d[i].fflag_c1 = nullptr;
    } else {
      connect3d[i].cwhich_c1 = cwhich_c1[i];
      connect3d[i].nside_c1 = nside_c1[i];
      connect3d[i].fflag_c1 = fflag_c1[i];
    }

    if (connect3d[i].nc2 == 0) {
      connect3d[i].cwhich_c2 = nullptr;
      connect3d[i].nside_c2 = nullptr;
      connect3d[i].fflag_c2 = nullptr;
    } else {
      connect3d[i].cwhich_c2 = cwhich_c2[i];
      connect3d[i].nside_c2 = nside_c2[i];
      connect3d[i].fflag_c2 = fflag_c2[i];
    }

    if (connect3d[i].nc3 == 0) {
      connect3d[i].cwhich_c3 = nullptr;
      connect3d[i].nside_c3 = nullptr;
      connect3d[i].fflag_c3 = nullptr;
    } else {
      connect3d[i].cwhich_c3 = cwhich_c3[i];
      connect3d[i].nside_c3 = nside_c3[i];
      connect3d[i].fflag_c3 = fflag_c3[i];
    }
  }

  // set connect3d cwhich/nside/aflag for each end point of each line
  // see fsg.h file for an explanation of each corner vector in Connect3d
  // aflag is based on dot of 2 connected tri normals to test if flat
  // if not flat, walk around corner via edge connections
  // if all paths have only convex connections, then mark as convex
  // otherwise, concave (may be hard to say if it flip flops)

  for (int i = 0; i < ntris; i++) {
    for (m = 0; m < connect3d[i].nc1; m++) {
      j = connect3d[i].neigh_c1[m];
      int cwhich = -1;
      if (tris[i].p1 == tris[j].p1) cwhich = 0;
      else if (tris[i].p1 == tris[j].p2) cwhich = 1;
      else if (tris[i].p1 == tris[j].p3) cwhich = 2;
      corner_connection3d(tris[i].norm,tris[j].norm,cwhich,flatthresh,
                          connect3d[i].fflag_c1[m],connect3d[i].nside_c1[m]);
      connect3d[i].cwhich_c1[m] = cwhich;
    }

    for (m = 0; m < connect3d[i].nc2; m++) {
      j = connect3d[i].neigh_c2[m];
      int cwhich = -1;
      if (tris[i].p2 == tris[j].p1) cwhich = 0;
      else if (tris[i].p2 == tris[j].p2) cwhich = 1;
      else if (tris[i].p2 == tris[j].p3) cwhich = 2;
      corner_connection3d(tris[i].norm,tris[j].norm,cwhich,flatthresh,
                          connect3d[i].fflag_c2[m],connect3d[i].nside_c2[m]);
      connect3d[i].cwhich_c2[m] = cwhich;
    }

    for (m = 0; m < connect3d[i].nc3; m++) {
      j = connect3d[i].neigh_c3[m];
      int cwhich = -1;
      if (tris[i].p3 == tris[j].p1) cwhich = 0;
      else if (tris[i].p3 == tris[j].p2) cwhich = 1;
      else if (tris[i].p3 == tris[j].p3) cwhich = 2;
      corner_connection3d(tris[i].norm,tris[j].norm,cwhich,flatthresh,
                          connect3d[i].fflag_c3[m],connect3d[i].nside_c3[m]);
      connect3d[i].cwhich_c3[m] = cwhich;
    }
  }

  // determine whether edges/pts are external based on connectivity

  for (int i = 0; i < nsurf; i++)
    for (int a = 0; a < 3; a++)
      connect3d[i].external_edge[a] = INTERNAL;

  for (int i = 0; i < nsurf; i++) {
    // external edge if there's a nonflat connection
    for (int n = 0; n < connect3d[i].ne1; n++)
      if (connect3d[i].fflag_e1[n] != FLAT)
        connect3d[i].external_edge[0] = EXTERNAL;
    for (int n = 0; n < connect3d[i].ne2; n++)
      if (connect3d[i].fflag_e2[n] != FLAT)
        connect3d[i].external_edge[1] = EXTERNAL;
    for (int n = 0; n < connect3d[i].ne3; n++)
      if (connect3d[i].fflag_e3[n] != FLAT)
        connect3d[i].external_edge[2] = EXTERNAL;

    // or unconnected on border
    if (connect3d[i].ne1 == 0)
      connect3d[i].external_edge[0] = UNCONNECTED;
    if (connect3d[i].ne2 == 0)
      connect3d[i].external_edge[1] = UNCONNECTED;
    if (connect3d[i].ne3 == 0)
      connect3d[i].external_edge[2] = UNCONNECTED;

    // corners basically inherit from connected edges
    connect3d[i].external_cor[0] = MAX(connect3d[i].external_edge[0], connect3d[i].external_edge[2]);
    connect3d[i].external_cor[1] = MAX(connect3d[i].external_edge[0], connect3d[i].external_edge[1]);
    connect3d[i].external_cor[2] = MAX(connect3d[i].external_edge[1], connect3d[i].external_edge[2]);
  }
}

/* ----------------------------------------------------------------------
   print stats on all lines and their connections
------------------------------------------------------------------------- */

void FixSurfaceGlobal::stats2d()
{
  double size;
  double delta[3];

  int nconnect = 0;
  int nfree = 0;
  double minsize = BIG;
  double maxsize = 0.0;

  for (int i = 0; i < nlines; i++) {
    nconnect += connect2d[i].np1 + connect2d[i].np2;
    if (connect2d[i].np1 == 0) nfree++;
    if (connect2d[i].np2 == 0) nfree++;
    MathExtra::sub3(points[lines[i].p1].x,points[lines[i].p2].x,delta);
    size = MathExtra::len3(delta);
    minsize = MIN(minsize,size);
    maxsize = MAX(maxsize,size);
  }

  nconnect /= 2;

  if (comm->me == 0) {
    utils::logmesg(lmp,"Fix surface/global line segment creation:\n");
    utils::logmesg(lmp,"  {} lines\n",nlines);
    utils::logmesg(lmp,"  {} line end points\n",npoints);
    utils::logmesg(lmp,"  {} end point connections\n",nconnect);
    utils::logmesg(lmp,"  {} free end points\n",nfree);
    utils::logmesg(lmp,"  {} min line length\n",minsize);
    utils::logmesg(lmp,"  {} max line length\n",maxsize);
  }
}

/* ----------------------------------------------------------------------
   print stats on all tris and their connections
------------------------------------------------------------------------- */

void FixSurfaceGlobal::stats3d()
{
  double size,area;
  double delta[3],edge12[3],edge13[3],cross[3];

  int nconnect_edge = 0;
  int nconnect_corner = 0;
  int nfree_edge = 0;
  int nfree_corner = 0;
  double minedge = BIG;
  double maxedge = 0.0;
  double minarea = BIG;
  double maxarea = 0.0;

  for (int i = 0; i < ntris; i++) {
    nconnect_edge += connect3d[i].ne1 + connect3d[i].ne2 + connect3d[i].ne3;
    nconnect_corner += connect3d[i].nc1 + connect3d[i].nc2 + connect3d[i].nc3;

    if (connect3d[i].ne1 == 0) nfree_edge++;
    if (connect3d[i].ne2 == 0) nfree_edge++;
    if (connect3d[i].ne3 == 0) nfree_edge++;

    // a free corner point requires 2 adjacent edges also have no connections

    if (connect3d[i].nc1 == 0 && (connect3d[i].ne3 == 0 && connect3d[i].ne1 == 0)) nfree_corner++;
    if (connect3d[i].nc2 == 0 && (connect3d[i].ne1 == 0 && connect3d[i].ne2 == 0)) nfree_corner++;
    if (connect3d[i].nc3 == 0 && (connect3d[i].ne2 == 0 && connect3d[i].ne3 == 0)) nfree_corner++;

    MathExtra::sub3(points[tris[i].p1].x,points[tris[i].p2].x,delta);
    size = MathExtra::len3(delta);
    minedge = MIN(minedge,size);
    maxedge = MAX(maxedge,size);
    MathExtra::sub3(points[tris[i].p2].x,points[tris[i].p3].x,delta);
    size = MathExtra::len3(delta);
    minedge = MIN(minedge,size);
    maxedge = MAX(maxedge,size);
    MathExtra::sub3(points[tris[i].p3].x,points[tris[i].p1].x,delta);
    size = MathExtra::len3(delta);
    minedge = MIN(minedge,size);
    maxedge = MAX(maxedge,size);

    MathExtra::sub3(points[tris[i].p2].x,points[tris[i].p1].x,edge12);
    MathExtra::sub3(points[tris[i].p3].x,points[tris[i].p1].x,edge13);
    MathExtra::cross3(edge12,edge13,cross);
    area = 0.5 * MathExtra::len3(cross);
    minarea = MIN(minarea,area);
    maxarea = MAX(maxarea,area);
  }

  nconnect_edge /= 2;
  nconnect_corner /= 2;

  if (comm->me == 0) {
    utils::logmesg(lmp,"Fix surface/global triangle creation:\n");
    utils::logmesg(lmp,"  {} tris\n",ntris);
    utils::logmesg(lmp,"  {} tri edges\n",nedges);
    utils::logmesg(lmp,"  {} tri corner points\n",npoints);
    utils::logmesg(lmp,"  {} edge connections\n",nconnect_edge);
    utils::logmesg(lmp,"  {} corner point connections\n",nconnect_corner);
    utils::logmesg(lmp,"  {} free edges\n",nfree_edge);
    utils::logmesg(lmp,"  {} free corner points\n",nfree_corner);
    utils::logmesg(lmp,"  {} min edge length\n",minedge);
    utils::logmesg(lmp,"  {} max edge length\n",maxedge);
    utils::logmesg(lmp,"  {} min tri area\n",minarea);
    utils::logmesg(lmp,"  {} max tri area\n",maxarea);
  }
}

/* ----------------------------------------------------------------------
   set attributes of all lines or tris
   xsurf,vsurf,omegasurf,norm
------------------------------------------------------------------------- */

void FixSurfaceGlobal::surface_attributes()
{
  double delta[3],p12[3],p13[3];
  double *p1,*p2,*p3;
  double zunit[3] = {0.0,0.0,1.0};

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

      MathExtra::sub3(p2,p1,p12);
      radsurf[i] = 0.5 * MathExtra::len3(p12);

      MathExtra::cross3(zunit,p12,lines[i].norm);
      MathExtra::norm3(lines[i].norm);
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

      MathExtra::sub3(p2,p1,p12);
      MathExtra::sub3(p3,p1,p13);
      MathExtra::cross3(p12,p13,tris[i].norm);
      MathExtra::norm3(tris[i].norm);
    }
  }

  for (int i = 0; i < nsurf; i++) {
    vsurf[i][0] = vsurf[i][1] = vsurf[i][2] = 0.0;
    omegasurf[i][0] = omegasurf[i][1] = omegasurf[i][2] = 0.0;
  }
}

/* -------------------------------------------------------------------------
   linear move: X = X0 + V*dt
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_linear(int imotion, int i)
{
  Motion *motion = &motions[imotion];
  double time_origin = motion->time_origin;
  double delta = (update->ntimestep - time_origin) * dt;

  int vxflag = motion->vxflag;
  int vyflag = motion->vyflag;
  int vzflag = motion->vzflag;
  double vx = motion->vx;
  double vy = motion->vy;
  double vz = motion->vz;

  // points - use of pointmove only moves a point once

  int pindex;
  double *pt;

  if (dimension == 2) {
    pindex = lines[i].p1;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxflag) pt[0] = points_original[i][0] + vx * delta;
      if (vyflag) pt[1] = points_original[i][1] + vy * delta;
      pointmove[pindex] = 1;
    }
    pindex = lines[i].p2;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxflag) pt[0] = points_original[i][0] + vx * delta;
      if (vyflag) pt[1] = points_original[i][1] + vy * delta;
      pointmove[pindex] = 1;
    }

  } else {
    pindex = tris[i].p1;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxflag) pt[0] = points_original[i][0] + vx * delta;
      if (vyflag) pt[1] = points_original[i][1] + vy * delta;
      if (vzflag) pt[2] = points_original[i][2] + vz * delta;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p2;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxflag) pt[0] = points_original[i][0] + vx * delta;
      if (vyflag) pt[1] = points_original[i][1] + vy * delta;
      if (vzflag) pt[2] = points_original[i][2] + vz * delta;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p3;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxflag) pt[0] = points_original[i][0] + vx * delta;
      if (vyflag) pt[1] = points_original[i][1] + vy * delta;
      if (vzflag) pt[2] = points_original[i][2] + vz * delta;
      pointmove[pindex] = 1;
    }
  }

  // xsurf and vsurf

  if (vxflag) {
    vsurf[i][0] = vx;
    xsurf[i][0] = xsurf_original[i][0] + vx * delta;
  }
  if (vyflag) {
    vsurf[i][1] = vy;
    xsurf[i][1] = xsurf_original[i][1] + vy * delta;
  }
  if (vzflag) {
    vsurf[i][2] = vz;
    xsurf[i][2] = xsurf_original[i][2] + vz * delta;
  }
}

/* -------------------------------------------------------------------------
   wiggle move: X = X0 + A sin(w*dt)
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_wiggle(int imotion, int i)
{
  Motion *motion = &motions[imotion];
  double time_origin = motion->time_origin;
  double omega = motion->omega;
  double delta = (update->ntimestep - time_origin) * dt;
  double arg = omega * delta;
  double sine = sin(arg);
  double cosine = cos(arg);

  int axflag = motion->axflag;
  int ayflag = motion->ayflag;
  int azflag = motion->azflag;
  double ax = motion->ax;
  double ay = motion->ay;
  double az = motion->az;

  // points - use of pointmove only moves a point once

  int pindex;
  double *pt;

  if (dimension == 2) {
    pindex = lines[i].p1;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (axflag) pt[0] = points_original[i][0] + ax * sine;
      if (ayflag) pt[1] = points_original[i][1] + ay * sine;
      pointmove[pindex] = 1;
    }
    pindex = lines[i].p2;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (axflag) pt[0] = points_original[i][0] + ax * sine;
      if (ayflag) pt[1] = points_original[i][1] + ay * sine;
      pointmove[pindex] = 1;
    }

  } else {
    pindex = tris[i].p1;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (axflag) pt[0] = points_original[i][0] + ax * sine;
      if (ayflag) pt[1] = points_original[i][1] + ay * sine;
      if (azflag) pt[2] = points_original[i][2] + az * sine;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p2;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (axflag) pt[0] = points_original[i][0] + ax * sine;
      if (ayflag) pt[1] = points_original[i][1] + ay * sine;
      if (azflag) pt[2] = points_original[i][2] + az * sine;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p3;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (axflag) pt[0] = points_original[i][0] + ax * sine;
      if (ayflag) pt[1] = points_original[i][1] + ay * sine;
      if (azflag) pt[2] = points_original[i][2] + az * sine;
      pointmove[pindex] = 1;
    }
  }

  // xsurf and vsurf

  if (axflag) {
    vsurf[i][0] = ax * omega * cosine;
    xsurf[i][0] = xsurf_original[i][0] + ax * sine;
  }
  if (ayflag) {
    vsurf[i][1] = ay * omega * cosine;
    xsurf[i][1] = xsurf_original[i][1] + ay * sine;
  }
  if (azflag) {
    vsurf[i][2] = az * omega * cosine;
    xsurf[i][2] = xsurf_original[i][2] + az * sine;
  }
}

/* -------------------------------------------------------------------------
   rotate move: rotate by right-hand rule around omega
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_rotate(int imotion, int i)
{
  Motion *motion = &motions[imotion];

  double time_origin = motion->time_origin;
  double omega = motion->omega;
  double *rpoint = motion->point;
  double *runit = motion->unit;

  double delta = (update->ntimestep - time_origin) * dt;
  double arg = omega * delta;
  double cosine = cos(arg);
  double sine = sin(arg);

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

  // points - use of pointmove only moves a point once

  int pindex;

  if (dimension == 2) {
    pindex = lines[i].p1;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      pointmove[pindex] = 1;
    }
    pindex = lines[i].p2;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      pointmove[pindex] = 1;
    }

  } else {
    pindex = tris[i].p1;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p2;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p3;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      pointmove[pindex] = 1;
    }
  }

  // xsurf and vsurf

  double ddotr;
  double a[3],b[3],c[3],d[3],disp[3];

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
  vsurf[i][0] = omega * (runit[1]*disp[2] - runit[2]*disp[1]);
  vsurf[i][1] = omega * (runit[2]*disp[0] - runit[0]*disp[2]);
  vsurf[i][2] = omega * (runit[0]*disp[1] - runit[1]*disp[0]);

  // normals

  double p12[3],p13[3];
  double *p1,*p2,*p3;

  if (dimension == 2) {
    double zunit[3] = {0.0,0.0,1.0};
    p1 = points[lines[i].p1].x;
    p2 = points[lines[i].p2].x;
    MathExtra::sub3(p2,p1,p12);
    MathExtra::cross3(zunit,p12,lines[i].norm);
    MathExtra::norm3(lines[i].norm);

  } else {
    p1 = points[tris[i].p1].x;
    p2 = points[tris[i].p2].x;
    p3 = points[tris[i].p3].x;
    MathExtra::sub3(p1,p2,p12);
    MathExtra::sub3(p1,p3,p13);
    MathExtra::cross3(p12,p13,tris[i].norm);
    MathExtra::norm3(tris[i].norm);
  }
}

/* -------------------------------------------------------------------------
   transrotate move:
   rotate by right-hand rule around omega
   add translation after rotation
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_transrotate(int imotion, int i)
{
  Motion *motion = &motions[imotion];

  double time_origin = motion->time_origin;
  double omega = motion->omega;
  double *rpoint = motion->point;
  double *runit = motion->unit;

  double delta = (update->ntimestep - time_origin) * dt;
  double arg = omega * delta;
  double cosine = cos(arg);
  double sine = sin(arg);

  int vxflag = motion->vxflag;
  int vyflag = motion->vyflag;
  int vzflag = motion->vzflag;
  double vx = motion->vx;
  double vy = motion->vy;
  double vz = motion->vz;

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

  // points - use of pointmove only moves a point once

  int pindex;
  double *pt;

  if (dimension == 2) {
    pindex = lines[i].p1;
    pt = points[pindex].x;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      if (vxflag) pt[0] += vx * delta;
      if (vyflag) pt[1] += vy * delta;
      pointmove[pindex] = 1;
    }
    pindex = lines[i].p2;
    pt = points[pindex].x;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      if (vxflag) pt[0] += vx * delta;
      if (vyflag) pt[1] += vy * delta;
      pointmove[pindex] = 1;
    }

  } else {
    pindex = tris[i].p1;
    pt = points[pindex].x;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      if (vxflag) pt[0] += vx * delta;
      if (vyflag) pt[1] += vy * delta;
      if (vzflag) pt[2] += vz * delta;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p2;
    pt = points[pindex].x;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      if (vxflag) pt[0] += vx * delta;
      if (vyflag) pt[1] += vy * delta;
      if (vzflag) pt[2] += vz * delta;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p3;
    pt = points[pindex].x;
    if (!pointmove[pindex]) {
      move_rotate_point(pindex,rpoint,runit,cosine,sine);
      if (vxflag) pt[0] += vx * delta;
      if (vyflag) pt[1] += vy * delta;
      if (vzflag) pt[2] += vz * delta;
      pointmove[pindex] = 1;
    }
  }

  // xsurf and vsurf

  double ddotr;
  double a[3],b[3],c[3],d[3],disp[3];

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
  vsurf[i][0] = omega * (runit[1]*disp[2] - runit[2]*disp[1]);
  vsurf[i][1] = omega * (runit[2]*disp[0] - runit[0]*disp[2]);
  vsurf[i][2] = omega * (runit[0]*disp[1] - runit[1]*disp[0]);

  if (vxflag) xsurf[i][0] += vx*delta;
  if (vyflag) xsurf[i][1] += vy*delta;
  if (vzflag) xsurf[i][2] += vz*delta;
  if (vxflag) vsurf[i][0] += vx;
  if (vyflag) vsurf[i][1] += vy;
  if (vzflag) vsurf[i][2] += vz;

  // normals

  double p12[3],p13[3];
  double *p1,*p2,*p3;

  if (dimension == 2) {
    double zunit[3] = {0.0,0.0,1.0};
    p1 = points[lines[i].p1].x;
    p2 = points[lines[i].p2].x;
    MathExtra::sub3(p2,p1,p12);
    MathExtra::cross3(zunit,p12,lines[i].norm);
    MathExtra::norm3(lines[i].norm);

  } else {
    p1 = points[tris[i].p1].x;
    p2 = points[tris[i].p2].x;
    p3 = points[tris[i].p3].x;
    MathExtra::sub3(p1,p2,p12);
    MathExtra::sub3(p1,p3,p13);
    MathExtra::cross3(p12,p13,tris[i].norm);
    MathExtra::norm3(tris[i].norm);
  }
}

/* -------------------------------------------------------------------------
   rotate point I by right-hand rule around omega
   ------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_rotate_point(int i, double *rpoint, double *runit,
                                         double cosine, double sine)
{
  double a[3],b[3],c[3],d[3],disp[3];

  d[0] = points_original[i][0] - rpoint[0];
  d[1] = points_original[i][1] - rpoint[1];
  d[2] = points_original[i][2] - rpoint[2];

  double ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
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

  double *pt = points[i].x;
  pt[0] = rpoint[0] + c[0] + disp[0];
  pt[1] = rpoint[1] + c[1] + disp[1];
  pt[2] = rpoint[2] + c[2] + disp[2];
}

/* -------------------------------------------------------------------------
   variable move:
   apply displacement variables if not NULL
   apply velocity variables if not NULL
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_variable(int imotion, int i)
{
  Motion *motion = &motions[imotion];

  char *xvarstr = motion->xvarstr;
  char *yvarstr = motion->yvarstr;
  char *zvarstr = motion->zvarstr;
  char *vxvarstr = motion->vxvarstr;
  char *vyvarstr = motion->vyvarstr;
  char *vzvarstr = motion->vzvarstr;
  double dx = motion->dx;
  double dy = motion->dy;
  double dz = motion->dz;
  double vx = motion->vx;
  double vy = motion->vy;
  double vz = motion->vz;

  double dt = update->dt;

  // points - use of pointmove only moves a point once
  // if displacement is variable, set x
  // if only velocity is variable, time-integrate x using v
  // if neither is variable, no change to x

  int pindex;
  double *pt;

  if (dimension == 2) {
    pindex = lines[i].p1;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxvarstr && !xvarstr) pt[0] += dt*dx;
      else if (xvarstr) pt[0] = points_original[i][0] + dx;
      if (vyvarstr && !yvarstr) pt[1] += dt*dy;
      else if (yvarstr) pt[1] = points_original[i][1] + dy;
      pointmove[pindex] = 1;
    }
    pindex = lines[i].p2;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxvarstr && !xvarstr) pt[0] += dt*dx;
      else if (xvarstr) pt[0] = points_original[i][0] + dx;
      if (vyvarstr && !yvarstr) pt[1] += dt*dy;
      else if (yvarstr) pt[1] = points_original[i][1] + dy;
      pointmove[pindex] = 1;
    }

  } else {
    pindex = tris[i].p1;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxvarstr && !xvarstr) pt[0] += dt*dx;
      else if (xvarstr) pt[0] = points_original[i][0] + dx;
      if (vyvarstr && !yvarstr) pt[1] += dt*dy;
      else if (yvarstr) pt[1] = points_original[i][1] + dy;
      if (vzvarstr && !zvarstr) pt[2] += dt*dz;
      else if (zvarstr) pt[2] = points_original[i][2] + dz;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p2;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxvarstr && !xvarstr) pt[0] += dt*dx;
      else if (xvarstr) pt[0] = points_original[i][0] + dx;
      if (vyvarstr && !yvarstr) pt[1] += dt*dy;
      else if (yvarstr) pt[1] = points_original[i][1] + dy;
      if (vzvarstr && !zvarstr) pt[2] += dt*dz;
      else if (zvarstr) pt[2] = points_original[i][2] + dz;
      pointmove[pindex] = 1;
    }
    pindex = tris[i].p3;
    if (!pointmove[pindex]) {
      pt = points[pindex].x;
      if (vxvarstr && !xvarstr) pt[0] += dt*dx;
      else if (xvarstr) pt[0] = points_original[i][0] + dx;
      if (vyvarstr && !yvarstr) pt[1] += dt*dy;
      else if (yvarstr) pt[1] = points_original[i][1] + dy;
      if (vzvarstr && !zvarstr) pt[2] += dt*dz;
      else if (zvarstr) pt[2] = points_original[i][2] + dz;
      pointmove[pindex] = 1;
    }
  }

  // xsurf and vsurf
  // if both displacement and velocity are variable, use them to set both x and v
  // if only displacement is variable, only set x
  // if only velocity is variable, set v and time-integrate x using v
  // if neither is variable, no change to x and v

  if (xvarstr && vxvarstr) {
    vsurf[i][0] = vx;
    xsurf[i][0] = xsurf_original[i][0] + dx;
  } else if (xvarstr) {
    xsurf[i][0] = xsurf_original[i][0] + dx;
  } else if (vxvarstr) {
    vsurf[i][0] = vx;
    xsurf[i][0] += dt*dx;
  }

  if (yvarstr && vyvarstr) {
    vsurf[i][1] = vy;
    xsurf[i][1] = xsurf_original[i][1] + dy;
  } else if (yvarstr) {
    xsurf[i][1] = xsurf_original[i][1] + dy;
  } else if (vyvarstr) {
    vsurf[i][1] = vy;
    xsurf[i][1] += dt*dy;
  }

  if (zvarstr && vzvarstr) {
    vsurf[i][2] = vx;
    xsurf[i][2] = xsurf_original[i][2] + dz;
  } else if (zvarstr) {
    xsurf[i][2] = xsurf_original[i][2] + dz;
  } else if (vzvarstr) {
    vsurf[i][2] = vz;
    xsurf[i][2] += dt*dz;
  }
}

/* ----------------------------------------------------------------------
  walk through contacting connections and determine side of contact
  start with "closest" surf and then walk all 1st, 2nd, ... connections
------------------------------------------------------------------------- */

void FixSurfaceGlobal::prewalk_connections2d()
{
  std::map<int, int> to_walk;
  std::unordered_set<int> walked;

  int j = contact_surfs[0].index;
  to_walk[j] = contact_surfs[0].nside;

  int k, n, m, nsidej, nsidek, nconnect, nc;
  std::tuple<int, int> element;
  while (!to_walk.empty()) {
    auto it = to_walk.begin();
    element = *it;
    j = it->first;
    nsidej = it->second;
    to_walk.erase(it);
    walked.insert(j);

    n = contacts_map[j];

    if (nsidej == OPPOSITE_SIDE)
      MathExtra::negate3(contact_surfs[n].surf_norm);

    for (nconnect = 0; nconnect < (connect2d[j].np1 + connect2d[j].np2); nconnect++) {
      if (nconnect < connect2d[j].np1) {
        k = connect2d[j].neigh_p1[nconnect];
        nsidek = connect2d[j].nside_p1[nconnect];
      } else {
        k = connect2d[j].neigh_p2[nconnect - connect2d[j].np1];
        nsidek = connect2d[j].nside_p2[nconnect - connect2d[j].np1];
      }

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      if (walked.find(k) == walked.end() && to_walk.find(k) == to_walk.end()) {
        // which side is associated with the initial closest surf
        m = contacts_map[k];

        if (nsidej == OPPOSITE_SIDE)
          nsidek = FLIPSIDE(nsidek);

        contact_surfs[m].nside = nsidek;
        to_walk[k] = nsidek;
      }
    }

    // Check if there is another disconnected surf
    if (to_walk.empty()) {
      for (nc = 0; nc < contact_surfs.size(); nc++) {
        j = contact_surfs[nc].index;
        if (walked.find(j) == walked.end())
          to_walk[j] = contact_surfs[nc].nside;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::prewalk_connections3d()
{
  std::map<int, int> to_walk;
  std::unordered_set<int> walked;

  int j = contact_surfs[0].index;
  to_walk[j] = contact_surfs[0].nside;

  int k, n, m, nsidej, nsidek, nconnect, nc, ntotal;
  std::tuple<int, int> element;
  while (!to_walk.empty()) {
    auto it = to_walk.begin();
    element = *it;
    j = it->first;
    nsidej = it->second;
    to_walk.erase(it);
    walked.insert(j);

    n = contacts_map[j];

    if (nsidej == OPPOSITE_SIDE)
      MathExtra::negate3(contact_surfs[n].surf_norm);

    // Loop through edge-connected surfs
    ntotal = connect3d[j].ne1 + connect3d[j].ne2 + connect3d[j].ne3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      if (nconnect < connect3d[j].ne1) {
        nc = nconnect;
        k = connect3d[j].neigh_e1[nc];
        nsidek = connect3d[j].nside_e1[nc];
      } else if (nconnect < connect3d[j].ne1 + connect3d[j].ne2) {
        nc = nconnect - connect3d[j].ne1;
        k = connect3d[j].neigh_e2[nc];
        nsidek = connect3d[j].nside_e2[nc];
      } else {
        nc = nconnect - connect3d[j].ne1 - connect3d[j].ne2;
        k = connect3d[j].neigh_e3[nc];
        nsidek = connect3d[j].nside_e3[nc];
      }

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      if (walked.find(k) == walked.end() && to_walk.find(k) == to_walk.end()) {
        // which side is associated with the initial closest surf
        m = contacts_map[k];
        if (nsidej == OPPOSITE_SIDE)
          nsidek = FLIPSIDE(nsidek);

        contact_surfs[m].nside = nsidek;
        to_walk[k] = nsidek;
      }
    }

    // Loop through corner-connected surfs
    ntotal = connect3d[j].nc1 + connect3d[j].nc2 + connect3d[j].nc3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      if (nconnect < connect3d[j].nc1) {
        nc = nconnect;
        k = connect3d[j].neigh_c1[nc];
        nsidek = connect3d[j].nside_c1[nc];
      } else if (nconnect < connect3d[j].nc1 + connect3d[j].nc2) {
        nc = nconnect - connect3d[j].nc1;
        k = connect3d[j].neigh_c2[nc];
        nsidek = connect3d[j].nside_c2[nc];
      } else {
        nc = nconnect - connect3d[j].nc1 - connect3d[j].nc2;
        k = connect3d[j].neigh_c3[nc];
        nsidek = connect3d[j].nside_c3[nc];
      }

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      if (walked.find(k) == walked.end() && to_walk.find(k) == to_walk.end()) {
        // which side is associated with the initial closest surf
        m = contacts_map[k];
        if (nsidej == OPPOSITE_SIDE)
          nsidek = FLIPSIDE(nsidek);
        contact_surfs[m].nside = nsidek;
        to_walk[k] = nsidek;
      }
    }

    // Check if there is another disconnected surf
    if (to_walk.empty()) {
      for (nc = 0; nc < contact_surfs.size(); nc++) {
        j = contact_surfs[nc].index;
        if (walked.find(j) == walked.end())
          to_walk[j] = contact_surfs[nc].nside;
      }
    }
  }
}

/* ----------------------------------------------------------------------
  walk through flat connections and process any contacts
  start with "closest" surf and then walk all 1st, 2nd, ... connections
------------------------------------------------------------------------- */

void FixSurfaceGlobal::walk_connections2d(std::vector<int> &composite_surfs, std::unordered_set<int> &processed_contacts)
{
  std::set<int> to_walk;
  std::set<int> to_add;

  // Find next closest surface
  int j, n;
  for (n = 0; n < contact_surfs.size(); n++) {
    j = contact_surfs[n].index;

    if (processed_contacts.find(j) == processed_contacts.end()) {
      to_walk.insert(j);
      break;
    }
  }

  int k, m, jflag, aflag, fflag, nconnect, nc, contact_at_joint;
  while (!to_walk.empty()) {
    auto it = to_walk.begin();
    j = *it;
    to_walk.erase(it);

    n = contacts_map[j];
    processed_contacts.insert(j);
    composite_surfs.push_back(n);
    jflag = contact_surfs[n].flag;

    for (nconnect = 0; nconnect < (connect2d[j].np1 + connect2d[j].np2); nconnect++) {
      contact_at_joint = 0; // If j's contact is at j-k joint
      if (nconnect < connect2d[j].np1) {
        nc = nconnect;
        k = connect2d[j].neigh_p1[nc];
        aflag = connect2d[j].aflag_p1[nc];
        fflag = connect2d[j].fflag_p1[nc];
        if (jflag == -1)
          contact_at_joint = 1;
      } else {
        nc = nconnect - connect2d[j].np1;
        k = connect2d[j].neigh_p2[nc];
        aflag = connect2d[j].aflag_p2[nc];
        fflag = connect2d[j].fflag_p2[nc];
        if (jflag == -2)
          contact_at_joint = 2;
      }

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      m = contacts_map[k];

      if (contact_surfs[n].nside == OPPOSITE_SIDE) {
        if (aflag == CONVEX) aflag = CONCAVE;
        else aflag = CONVEX;
      }

      if (fflag == FLAT) {
        // flat, same-type: walk
        if (contact_surfs[n].type == contact_surfs[m].type &&
          processed_contacts.find(k) == processed_contacts.end())
          to_add.insert(k);
      } else if (aflag == CONVEX) {
        // must overlap w/in epsilon or higher priority to hide (can't walk around a composite surf to hide)
        if (contact_surfs[n].overlap > contact_surfs[m].overlap - EPSILON)
          contact_surfs[m].convex_index = j;
      }

      if (contact_at_joint) {
        contact_surfs[n].cindex.push_back(k);
        contact_surfs[n].caflag.push_back(aflag);
      }
    }

    // Add flat surfs to walk list, hiding if appropriate
    for (const int k : to_add) {
      if (contact_surfs[n].convex_index != -1) {
        m = contacts_map[k];
        contact_surfs[m].convex_index = j;
      }
      to_walk.insert(k);
    }
    to_add.clear();
  }
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::walk_connections3d(std::vector<int> &composite_surfs, std::unordered_set<int> &processed_contacts)
{
  std::set<int> to_walk;
  std::set<int> to_add;

  // Find next closest surface
  int j, n;
  for (n = 0; n < contact_surfs.size(); n++) {
    j = contact_surfs[n].index;

    if (processed_contacts.find(j) == processed_contacts.end()) {
      to_walk.insert(j);
      break;
    }
  }

  int k, m, jflag, aflag, fflag, which, nconnect, nc, ntotal, contact_at_joint;
  while (!to_walk.empty()) {
    auto it = to_walk.begin();
    j = *it;
    to_walk.erase(it);

    n = contacts_map[j];
    processed_contacts.insert(j);
    composite_surfs.push_back(n);
    jflag = contact_surfs[n].flag;

    // Loop through edge-connected surfs
    ntotal = connect3d[j].ne1 + connect3d[j].ne2 + connect3d[j].ne3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      contact_at_joint = 0; // If j's contact is at j-k joint
      if (nconnect < connect3d[j].ne1) {
        // e1 = p1+p2
        nc = nconnect;
        k = connect3d[j].neigh_e1[nc];
        aflag = connect3d[j].aflag_e1[nc];
        fflag = connect3d[j].fflag_e1[nc];
        which = 0;
        if (jflag == -1 || jflag == -4 || jflag == -5)
          contact_at_joint = 1;
      } else if (nconnect < connect3d[j].ne1 + connect3d[j].ne2) {
        // e2 = p2+p3
        nc = nconnect - connect3d[j].ne1;
        k = connect3d[j].neigh_e2[nc];
        aflag = connect3d[j].aflag_e2[nc];
        fflag = connect3d[j].fflag_e2[nc];
        which = 1;
        if (jflag == -2 || jflag == -5 || jflag == -6)
          contact_at_joint = 1;
      } else {
        // e3 = p1+p3
        nc = nconnect - connect3d[j].ne1 - connect3d[j].ne2;
        k = connect3d[j].neigh_e3[nc];
        aflag = connect3d[j].aflag_e3[nc];
        fflag = connect3d[j].fflag_e3[nc];
        which = 2;
        if (jflag == -3 || jflag == -4 || jflag == -6)
          contact_at_joint = 1;
      }

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      m = contacts_map[k];

      if (contact_surfs[n].nside == OPPOSITE_SIDE) {
        if (aflag == CONVEX) aflag = CONCAVE;
        else aflag = CONVEX;
      }

      if (fflag == FLAT) {
        // flat, same-type: walk
        if (contact_surfs[n].type == contact_surfs[m].type &&
          processed_contacts.find(k) == processed_contacts.end())
          to_add.insert(k);
      } else if (aflag == CONVEX) {
        // must overlap w/in epsilon or higher priority to hide (can't walk around a composite surf to hide)
        if (contact_surfs[n].overlap > contact_surfs[m].overlap - EPSILON)
          contact_surfs[m].convex_index = j;
      }

      if (contact_at_joint) {
        contact_surfs[n].cindex.push_back(k);
        contact_surfs[n].cwhich.push_back(which);
        contact_surfs[n].caflag.push_back(aflag);
      }
    }

    // Loop through corner-connected surfs to find any other flat connections
    ntotal = connect3d[j].nc1 + connect3d[j].nc2 + connect3d[j].nc3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      contact_at_joint = 0;
      if (nconnect < connect3d[j].nc1) {
        nc = nconnect;
        k = connect3d[j].neigh_c1[nc];
        fflag = connect3d[j].fflag_c1[nc];
        if (jflag == -4)
          contact_at_joint = 1;
      } else if (nconnect < connect3d[j].nc1 + connect3d[j].nc2) {
        nc = nconnect - connect3d[j].nc1;
        k = connect3d[j].neigh_c2[nc];
        fflag = connect3d[j].fflag_c2[nc];
        if (jflag == -5)
          contact_at_joint = 1;
      } else {
        nc = nconnect - connect3d[j].nc1 - connect3d[j].nc2;
        k = connect3d[j].neigh_c3[nc];
        fflag = connect3d[j].fflag_c3[nc];
        if (jflag == -6)
          contact_at_joint = 1;
      }

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      m = contacts_map[k];

      if (fflag == FLAT) {
        // flat, same-type: walk
        if (contact_surfs[n].type == contact_surfs[m].type &&
          processed_contacts.find(k) == processed_contacts.end())
          to_add.insert(k);
      }

      if (contact_at_joint) {
        contact_surfs[n].cindex.push_back(k);
        contact_surfs[n].cwhich.push_back(-1);
        contact_surfs[n].caflag.push_back(-1);
      }
    }

    // Add flat surfs to walk list, hiding if appropriate
    for (const int k : to_add) {
      if (contact_surfs[n].convex_index != -1) {
        m = contacts_map[k];
        contact_surfs[m].convex_index = j;
      }
      to_walk.insert(k);
    }
    to_add.clear();
  }
}

/* ----------------------------------------------------------------------
   Calculate key variables needed for forces (overlap and direction)
------------------------------------------------------------------------- */

double FixSurfaceGlobal::calculate_2d_forces(std::vector<int> &composite_surfs)
{
  int n, m, j, k, external, flag;
  double dot, overlap;
  double jnorm[3], dr[3], kline[3];

  // Check if composite is hidden (convex) and calc max overlaps

  double max_overlap = -BIG;
  double max_overlap_ext = -BIG;
  for (auto it = 0; it < composite_surfs.size(); it++) {
    n = composite_surfs[it];
    j = contact_surfs[n].index;

    if (contact_surfs[n].convex_index != -1)
      contact_surfs[n].overlap = 0.0;

    overlap = contact_surfs[n].overlap;
    external = contact_surfs[n].external;

    max_overlap = MAX(max_overlap, overlap);
    if (external)
      max_overlap_ext = MAX(max_overlap_ext, overlap);
  }

  if (max_overlap < EPSILON)
    return max_overlap;

  // Smooth int/ext based on arbitrary ratio of overlaps
  double w_ext = 0.0;
  if (max_overlap_ext != -BIG) {
    w_ext = max_overlap_ext / max_overlap;
    w_ext *= w_ext;
  }
  double w_int = 1.0 - w_ext;

  // Calculate constraints on force norm
  int i, ck, pt, ptk, caflag;
  double max_dot;
  for (auto it = 0; it < composite_surfs.size(); it++) {
    n = composite_surfs[it];
    j = contact_surfs[n].index;
    flag = contact_surfs[n].flag;
    external = contact_surfs[n].external;
    overlap = contact_surfs[n].overlap;
    MathExtra::copy3(contact_surfs[n].surf_norm, jnorm);
    MathExtra::copy3(contact_surfs[n].dr, dr);

    // To focus smooth contributions towards edge when about to move off
    if (external) {
      contact_surfs[n].weight_contribution = w_ext;
    } else {
      contact_surfs[n].weight_contribution = w_int;
    }

    if (overlap < EPSILON) continue;

    // Interior of line, use default of dr (equiv to surf norm)
    if (flag == 1) continue;

    max_dot = -2;
    ck = -1;
    for (i = 0; i < contact_surfs[n].cindex.size(); i++) {
      k = contact_surfs[n].cindex[i];
      m = contacts_map[k];
      dot = MathExtra::dot3(jnorm, contact_surfs[m].surf_norm);

      if (dot > max_dot) {
        ck = k;
        caflag = contact_surfs[n].caflag[i];
        max_dot = dot;
      }
    }

    if (ck != -1) {
      if (caflag == CONCAVE) {
        MathExtra::copy3(jnorm, contact_surfs[n].dr_force);
      } else {

        // See if dr has component pointing into other (k) line
        if (flag == -1) {
          pt = lines[j].p1;
        } else {
          pt = lines[j].p2;
        }

        if (pt == lines[ck].p1) {
          ptk = lines[ck].p2;
        } else {
          ptk = lines[ck].p1;
        }

        MathExtra::sub3(points[ptk].x, points[pt].x, kline);
        dot = MathExtra::dot3(dr, kline);

        // If so, just use surface norm
        if (dot > 0) {
          m = contacts_map[ck];
          MathExtra::copy3(contact_surfs[m].surf_norm, contact_surfs[n].dr_force);
        }
      }
    }
  }

  return max_overlap;
}

/* ----------------------------------------------------------------------
   Calculate key variables needed for forces (overlap and direction)
    Note: these rules can be quite complicated (particularly for unconnected
          edges). Better formulations surely exist and may be a topic for
          revisiting in the future.
          Currently it relies on a lot of smoothing/interpolation factors
          to ensure the direction of forces vary continuously.
------------------------------------------------------------------------- */

double FixSurfaceGlobal::calculate_3d_forces(std::vector<int> &composite_surfs)
{
  int n, m, i, j, k, external, flag;
  double dot, dist, rmag, overlap;
  double jnorm[3], knorm[3], dr[3];

  // -----------------------------------
  //   total smoothing/hidden status
  // -----------------------------------

  double max_overlap = -BIG;
  double max_overlap_ext = -BIG;
  int uc_flag = 0;

  // Find if surface is hidden and/or whether it's unconnected
  for (auto it = 0; it < composite_surfs.size(); it++) {
    n = composite_surfs[it];
    j = contact_surfs[n].index;

    if (contact_surfs[n].convex_index != -1)
      contact_surfs[n].overlap = 0.0;

    overlap = contact_surfs[n].overlap;
    external = contact_surfs[n].external;

    max_overlap = MAX(max_overlap, overlap);
    if (external == EXTERNAL)
      max_overlap_ext = MAX(max_overlap_ext, overlap);

    if (external == UNCONNECTED)
      uc_flag = 1;
  }

  if (max_overlap < EPSILON)
    return max_overlap;

  // Smooth int/ext based on arbitrary ratio of overlaps
  double w_ext = 0.0;
  if (max_overlap_ext != -BIG) {
    w_ext = max_overlap_ext / max_overlap;
    w_ext *= w_ext;
  }
  double w_int = 1.0 - w_ext;

  // Find primary constraint for all corner/edge connections
  int which1, which2;
  double max_dot1, max_dot2;
  for (auto it = 0; it < composite_surfs.size(); it++) {
    n = composite_surfs[it];
    j = contact_surfs[n].index;
    flag = contact_surfs[n].flag;
    MathExtra::copy3(contact_surfs[n].surf_norm, jnorm);

    if (flag == -4) {
      which1 = 0;
      which2 = 2;
    } else if (flag == -5) {
      which1 = 0;
      which2 = 1;
    } else if (flag == -6) {
      which1 = 2;
      which2 = 1;
    }

    // If multiple constraints (e.g. a T), find which is closest aligned
    max_dot1 = max_dot2 = -2.0;
    contact_surfs[n].ck1 = contact_surfs[n].ck2 = -1;
    for (i = 0; i < contact_surfs[n].cindex.size(); i++) {
      k = contact_surfs[n].cindex[i];
      m = contacts_map[k];
      dot = MathExtra::dot3(jnorm, contact_surfs[m].surf_norm);

      if (flag > -4) {
        if (dot > max_dot1) {
          contact_surfs[n].ck1 = k;
          contact_surfs[n].caflag1 = contact_surfs[n].caflag[i];
          max_dot1 = dot;
        }
      } else {
        if (contact_surfs[n].cwhich[i] == which1 && dot > max_dot1) {
          contact_surfs[n].ck1 = k;
          contact_surfs[n].caflag1 = contact_surfs[n].caflag[i];
          max_dot1 = dot;
        }
        if (contact_surfs[n].cwhich[i] == which2 && dot > max_dot2) {
          contact_surfs[n].ck2 = k;
          contact_surfs[n].caflag2 = contact_surfs[n].caflag[i];
          max_dot2 = dot;
        }
      }
    }
  }

  // If unconnected, calculate distances for smoothing factor
  double w_connect = 1.0;
  if (uc_flag) {
    double max_dist_uc = 0.0;
    for (auto it = 0; it < composite_surfs.size(); it++) {
      n = composite_surfs[it];
      j = contact_surfs[n].index;
      flag = contact_surfs[n].flag;

      if (flag == 1) {
        max_dist_uc = 0.0;
        break;
      }

      which1 = which2 = -1;
      if (flag == -4) {
        which1 = 0;
        which2 = 2;
      } else if (flag == -5) {
        which1 = 0;
        which2 = 1;
      } else if (flag == -6) {
        which1 = 2;
        which2 = 1;
      } else {
        which1 = -(1 + flag); // -1 -> 0, -2 -> 1, -3 ->
        which2 = -1;
      }

      // calculate max in-plane dist from unconnected (UC) edges
      //   Note: ideally would check min, but this doesn't work
      //         along convex turns between UC edges, it'll go to
      //         zero when it aligns with one of the UC edges.
      //         Could revise if we had UC connectivity structures

      rmag = contact_surfs[n].rmag;
      MathExtra::copy3(contact_surfs[n].dr, dr);
      MathExtra::copy3(contact_surfs[n].surf_norm, jnorm);

      if (contact_surfs[n].external == UNCONNECTED) {

        if (flag < -3 && which1 >= 0 && which2 >= 0 && connect3d[j].external_edge[which1] == UNCONNECTED && connect3d[j].external_edge[which2] == UNCONNECTED) {
          // if a corner with 2 unconnected edges, just use distance
          max_dist_uc = MAX(max_dist_uc, rmag);
        } else {
          // else get distance along edge normal pointing outward
          if (connect3d[j].external_edge[0] == UNCONNECTED && (flag == -1 || flag == -4 || flag == -5)) {
            dist = dist_away_from_edge(points[tris[j].p1].x, points[tris[j].p2].x, points[tris[j].p3].x, jnorm, dr);
            if (dist != -1) max_dist_uc = MAX(max_dist_uc, dist * rmag);
          }

          if (connect3d[j].external_edge[1] == UNCONNECTED && (flag == -2 || flag == -5 || flag == -6)) {
            dist = dist_away_from_edge(points[tris[j].p2].x, points[tris[j].p3].x, points[tris[j].p1].x, jnorm, dr);
            if (dist != -1) max_dist_uc = MAX(max_dist_uc, dist * rmag);
          }

          if (connect3d[j].external_edge[2] == UNCONNECTED && (flag == -3 || flag == -4 || flag == -6)) {
            dist = dist_away_from_edge(points[tris[j].p1].x, points[tris[j].p3].x, points[tris[j].p2].x, jnorm, dr);
            if (dist != -1) max_dist_uc = MAX(max_dist_uc, dist * rmag);
          }
        }
      }
    }

    if (max_dist_uc > 0.0)
      w_connect = MAX(0.0, MIN(1.0, 1.0 - max_dist_uc / max_overlap));
  }

  // -----------------------------------
  //    per-surf calculations
  // -----------------------------------

  int pt, pt1, pt2, external1, external2, edge1_uc, edge2_uc;
  double w_in_plane, dot1a, dot2a, dot1xp, dot2xp, dot1ip, dot2ip, w1_in_plane, w2_in_plane, w1, w2, wtmp;
  double line1[3], line2[3], dr_in_plane[3];
  double dr1[3], dr2[3], fn1[3], fn2[3], fntot[3], normave[3];

  for (auto it = 0; it < composite_surfs.size(); it++) {
    n = composite_surfs[it];
    j = contact_surfs[n].index;

    MathExtra::copy3(contact_surfs[n].surf_norm, jnorm);
    MathExtra::copy3(contact_surfs[n].dr, dr);
    overlap = contact_surfs[n].overlap;
    rmag = contact_surfs[n].rmag;
    flag = contact_surfs[n].flag;
    external = contact_surfs[n].external;

    // Weight towards edges when about to move away from composite surf
    //   note this does not weight unconnected
    if (external == EXTERNAL) {
      contact_surfs[n].weight_contribution = w_ext;
    } else if (external == INTERNAL) {
      contact_surfs[n].weight_contribution = w_int;
    }

    // Weight internal contribution to zero when off surf and/or when in plane
    if (external != UNCONNECTED) {
      w_in_plane = 1.0;
      if (w_connect < 1.0) {
        dot = MathExtra::dot3(jnorm, dr);
        dist = rmag * dot;

        // Interpolation function with limits 1 when wc = 1 and 0 when dist = 0
        //   should never see scenario where w_connect = 1.0 and dist = 0...
        if (dist < rmag)
          w_in_plane = MAX(0.0, MIN(1.0, dist / (rmag * (1.0 - w_connect))));
      }

      contact_surfs[n].weight_contribution *= w_connect * w_in_plane;
    }
    if (overlap < EPSILON) continue;

    // Interior of tri, use default of dr (equiv to surf norm)
    if (flag == 1) continue;

    // ---------------------
    //     edge of tri
    // ---------------------

    if (flag < 0 && flag > -4) {
      if (contact_surfs[n].ck1 != -1) {
        m = contacts_map[contact_surfs[n].ck1];
        calculate_3d_edge_force(contact_surfs[n].caflag1, jnorm, contact_surfs[m].surf_norm, dr, contact_surfs[n].dr_force);
      }
    }

    // ---------------------
    //     corner of tri
    // ---------------------

    if (flag < -3) {
      if (flag == -4) {
        pt = tris[j].p1;
        pt1 = tris[j].p2;
        pt2 = tris[j].p3;
        which1 = 0;
        which2 = 2;
      } else if (flag == -5) {
        pt = tris[j].p2;
        pt1 = tris[j].p1;
        pt2 = tris[j].p3;
        which1 = 0;
        which2 = 1;
      } else {    // flag == -6 (overlap classification is bounded to -1..-6)
        pt = tris[j].p3;
        pt1 = tris[j].p1;
        pt2 = tris[j].p2;
        which1 = 2;
        which2 = 1;
      }

      external1 = connect3d[j].external_edge[which1];
      external2 = connect3d[j].external_edge[which2];

      MathExtra::sub3(points[pt].x, points[pt1].x, line1);
      MathExtra::sub3(points[pt].x, points[pt2].x, line2);
      MathExtra::norm3(line1);
      MathExtra::norm3(line2);

      dot = MathExtra::dot3(dr, jnorm);
      MathExtra::scaleadd3(-dot, jnorm, dr, dr_in_plane);
      dist = MathExtra::len3(dr_in_plane);

      // On top of tri, just use surf norm
      if (dist < EPSILON) {
        MathExtra::copy3(jnorm, contact_surfs[n].dr_force);
        return max_overlap;
      }

      MathExtra::scale3(1.0 / dist, dr_in_plane, dr_in_plane);

      // Prevent round off to negative (try removing after surf extra patched?)
      dot1a = MathExtra::dot3(line1, dr);
      dot2a = MathExtra::dot3(line2, dr);
      dot1ip = MAX(0.0, MathExtra::dot3(line1, dr_in_plane));
      dot2ip = MAX(0.0, MathExtra::dot3(line2, dr_in_plane));

      // Perpendicular to both lines, just use dr
      if (dot1ip < EPSILON && dot2ip < EPSILON) {
        MathExtra::copy3(dr, contact_surfs[n].dr_force);
        return max_overlap;
      }

      MathExtra::zero3(fntot);
      w1_in_plane = 1.0;
      w2_in_plane = 1.0;
      edge1_uc = 0;
      edge2_uc = 0;

      // default, use dr w/o component along edge
      dot = MathExtra::dot3(dr, line1);
      MathExtra::scaleadd3(-dot, line1, dr, dr1);
      MathExtra::norm3(dr1);
      MathExtra::copy3(dr1, fn1);

      dot = MathExtra::dot3(dr, line2);
      MathExtra::scaleadd3(-dot, line2, dr, dr2);
      MathExtra::norm3(dr2);
      MathExtra::copy3(dr2, fn2);

      // ---------- Edge 1 ----------

      // check constraints
      dot1xp = -1;
      if (contact_surfs[n].ck1 != -1) {
        m = contacts_map[contact_surfs[n].ck1];

        MathExtra::copy3(contact_surfs[m].surf_norm, knorm);
        calculate_3d_edge_force(contact_surfs[n].caflag1, jnorm, knorm, dr1, fn1);

        if (contact_surfs[n].caflag1 != CONCAVE && contact_surfs[m].convex_index != -1) {
          MathExtra::add3(jnorm, knorm, normave);
          MathExtra::norm3(normave);
          dot1xp = MathExtra::dot3(dr1, normave);
        }
      }

      if (uc_flag && external1 != UNCONNECTED) {
        if (w_connect < 1.0) {
          dist = rmag * MathExtra::dot3(jnorm, dr);
          if (dist < rmag)
            w1_in_plane = MAX(0.0, MIN(1.0, dist / (rmag * (1.0 - w_connect))));
        }
        edge1_uc = 1;
      }

      // ---------- Edge 2 ----------

      dot2xp = -1;
      if (contact_surfs[n].ck2 != -1) {
        m = contacts_map[contact_surfs[n].ck2];

        MathExtra::copy3(contact_surfs[m].surf_norm, knorm);
        calculate_3d_edge_force(contact_surfs[n].caflag2, jnorm, knorm, dr2, fn2);

        if (contact_surfs[n].caflag2 != CONCAVE && contact_surfs[m].convex_index != -1) {
          MathExtra::add3(jnorm, knorm, normave);
          MathExtra::norm3(normave);
          dot2xp = MathExtra::dot3(dr2, normave);
        }
      }

      if (uc_flag && external2 != UNCONNECTED) {
        if (w_connect < 1.0) {
          dist = rmag * MathExtra::dot3(jnorm, dr);
          if (dist < rmag)
            w2_in_plane = MAX(0.0, MIN(1.0, dist / (rmag * (1.0 - w_connect))));
        }
        edge2_uc = 1;
      }

      // ---------- Interpolation ----------

      // Turn off contribution zero when fully aligned (proj ill defined), unless unconnected
      w1 = 1 - dot1a;
      w2 = 1 - dot2a;

      // Turn off contribution when perpendicular to other edge
      w1 *= dot2ip;
      w2 *= dot1ip;

      // When being hidden, turn off when aligned in plane with other edge
      if (dot1xp != -1) w2 *= (1 - dot1xp);
      if (dot2xp != -1) w1 *= (1 - dot2xp);

      // If any component points into other line, cap at surf norm
      //   this can happen b/c default is dr-along-line-component
      //   exception is when outside & unconnected

      dot = MathExtra::dot3(fn1, line2);
      if (dot < 0.0 && (external1 != UNCONNECTED || w_connect == 1.0))
        MathExtra::copy3(jnorm, fn1);

      dot = MathExtra::dot3(fn2, line1);
      if (dot < 0.0 && (external2 != UNCONNECTED || w_connect == 1.0))
        MathExtra::copy3(jnorm, fn2);

      // if it's an unconnected corner w/ a connected edge,
      //   blend in- and out-of-plane solutions

      if (external == UNCONNECTED) {
        if (external1 != UNCONNECTED) {
          wtmp = w_connect * w1_in_plane;
          MathExtra::scaleadd3(wtmp, fn1, 1.0 - wtmp, dr1, fn1);
        }
        if (external2 != UNCONNECTED) {
          wtmp = w_connect * w2_in_plane;
          MathExtra::scaleadd3(wtmp, fn2, 1.0 - wtmp, dr2, fn2);
        }
      }

      MathExtra::scaleadd3(w1, fn1, fntot, fntot);
      MathExtra::scaleadd3(w2, fn2, fntot, fntot);

      // If weights both zero (e.g. if both edges convex) then use dr or norm
      if (w1 < EPSILON && w2 < EPSILON) {
        if (dot1xp != -1 && dot2xp != -1) {
          // corner of some polygon
          MathExtra::copy3(dr, fntot);
        } else if (dot1xp != -1) {
          // if 2 is a flat surf
          MathExtra::copy3(fn1, fntot); // not sure about these 3 options
        } else if (dot2xp != -1) {
          // if 1 is a flat surf
          MathExtra::copy3(fn2, fntot);
        } else {
          // if both are flat surfs
          MathExtra::add3(fn1, fn2, fntot);
        }
      }

      MathExtra::normalize3(fntot, contact_surfs[n].dr_force);
    }
  }

  return max_overlap;
}

/* ----------------------------------------------------------------------
   Calculate force contribution from a triangle edge
------------------------------------------------------------------------- */

void FixSurfaceGlobal::calculate_3d_edge_force(int aflag, double jnorm[3], double knorm[3], double drperp[3], double fn[3])
{
  if (aflag == CONCAVE) {
    MathExtra::copy3(jnorm, fn);
  } else {
    // cannot point beyond knorm
    double dotjk = MathExtra::dot3(jnorm, knorm);
    double dotjr = MathExtra::dot3(jnorm, drperp);
    double dotkr = MathExtra::dot3(knorm, drperp);

    double dotmin = MIN(dotjk, dotjr);
    dotmin = MIN(dotmin, dotkr);

    if (dotmin == dotjk) {
      MathExtra::copy3(drperp, fn);
    } else if (dotmin == dotjr) {
      MathExtra::copy3(knorm, fn);
    } else {
      MathExtra::copy3(jnorm, fn);
    }
  }
}

/* ----------------------------------------------------------------------
   Calculate planar distance of point from a tri edge
     the edge consists of pt1 & pt2 while pt3 is the 3rd point in the tri
     return -1 if dr antialigned with edge normal
------------------------------------------------------------------------- */

double FixSurfaceGlobal::dist_away_from_edge(double* pt1x, double* pt2x, double* pt3x, double *jnorm, double *dr)
{
  double edge[3], edge_norm[3], p3_line[3];
  MathExtra::sub3(pt1x, pt2x, edge);
  MathExtra::cross3(jnorm, edge, edge_norm);
  MathExtra::norm3(edge_norm);

  // Check if dr points into tri
  MathExtra::sub3(pt3x, pt1x, p3_line);
  double dot = MathExtra::dot3(edge_norm, p3_line);
  double sign_flip = 1.0;
  if (dot > 0.0) sign_flip = -1.0;

  dot = sign_flip * MathExtra::dot3(dr, edge_norm);

  if (dot > -EPSILON)
    return dot;

  return -1;
}

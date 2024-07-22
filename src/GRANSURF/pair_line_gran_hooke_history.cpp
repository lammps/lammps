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

#include "pair_line_gran_hooke_history.h"

#include "atom.h"
#include "atom_vec_line.h"
#include "comm.h"
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

/* ---------------------------------------------------------------------- */

PairLineGranHookeHistory::PairLineGranHookeHistory(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  history = 1;
  size_history = 3;

  nmax = 0;
  mass_rigid = nullptr;

  emax = 0;
  endpts = nullptr;

  // set comm size needed by this Pair if used with fix rigid

  comm_forward = 1;

  // create dummy fix as placeholder for FixNeighHistory
  // this is so final order of Modify:fix will conform to input script

  fix_history = nullptr;
  fix_dummy = dynamic_cast<FixDummy *>(
      modify->add_fix("NEIGH_HISTORY_HHLINE_DUMMY" + std::to_string(instance_me) + " all DUMMY"));
}

/* ---------------------------------------------------------------------- */

PairLineGranHookeHistory::~PairLineGranHookeHistory()
{
  delete [] svector;

  if (!fix_history)
    modify->delete_fix("NEIGH_HISTORY_HHLINE_DUMMY" + std::to_string(instance_me));
  else
    modify->delete_fix("NEIGH_HISTORY_HHLINE" + std::to_string(instance_me));

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    delete [] onerad_dynamic;
    delete [] onerad_frozen;
    delete [] maxrad_dynamic;
    delete [] maxrad_frozen;
  }

  memory->destroy(mass_rigid);
  memory->destroy(endpts);
}

/* ---------------------------------------------------------------------- */

void PairLineGranHookeHistory::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,inum,jnum,jflag,kflag,otherflag;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double factor_couple;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,mline;
  double damp,ccel,tor1,tor2,tor3;
  double fn,fs,ft,fs1,fs2,fs3;
  double shrmag,rsht;
  double dr[3],ds[3],vline[3],contact[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

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
        memory->create(mass_rigid,nmax,"line/gran:mass_rigid");
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
  }

  // pre-calculate current end pts of owned+ghost lines
  // only once per reneighbor if surfs not moving

  if (surfmoveflag || neighbor->ago == 0) calculate_endpts();

  // loop over neighbors of my atoms
  // I is always sphere, J is always line

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *line = atom->line;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firstshear = fix_history->firstvalue;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      /*
      if (touch[jj]) {
        shear = &allshear[3*jj];
        printf("PRE step %ld ij %d %d shear %g %g %g\n",
               update->ntimestep,atom->tag[i],tag[j],
               shear[0],shear[1],shear[2]);
      }
      */

      if (rsq >= radsum*radsum) {

        // unset non-touching neighbors
        // do same below if sphere/line are not overlapping

        touch[jj] = 0;
        shear = &allshear[3*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {

        // sanity check that neighbor list is built correctly

        if (line[i] >= 0 || line[j] < 0)
          error->one(FLERR,"Pair line/gran iteraction is invalid");

        // check for overlap of sphere and line segment
        // jflag = 0 for no overlap, 1 for interior line pt, -1/-2 for end pts
        // if no overlap, just continue
        // for overlap, also return:
        //   contact = nearest point on line to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr

        jflag = overlap_sphere_line(i,j,contact,dr,rsq);

        if (!jflag) {
          touch[jj] = 0;
          shear = &allshear[3*jj];
          shear[0] = 0.0;
          shear[1] = 0.0;
          shear[2] = 0.0;
          continue;
        }

        // if contact = line end pt:
        // check overlap status of line adjacent to the end pt
        // otherflag = 0/1 for this/other line performs calculation

        if (jflag < 0) {
          otherflag = endpt_neigh_check(i,j,jflag);
          if (otherflag) continue;
        }

        // NOTE: add logic to persist shear history if contact has changed

        // NOTE: add logic to check for coupled contacts and weight them

        factor_couple = 1.0;

        // ds = vector from line center to contact pt

        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        ds[0] = contact[0] - x[j][0];
        ds[1] = contact[1] - x[j][1];
        ds[2] = contact[2] - x[j][2];

        // vline = velocity of contact pt on line, translation + rotation

        vline[0] = v[j][0] + (omega[j][1]*ds[2] - omega[j][2]*ds[1]);
        vline[1] = v[j][1] + (omega[j][2]*ds[0] - omega[j][0]*ds[2]);
        vline[2] = v[j][2] + (omega[j][0]*ds[1] - omega[j][1]*ds[0]);

        // relative translational velocity

        vr1 = v[i][0] - vline[0];
        vr2 = v[i][1] - vline[1];
        vr3 = v[i][2] - vline[2];

        // normal component

        vnnr = vr1*dr[0] + vr2*dr[1] + vr3*dr[2];
        vn1 = dr[0]*vnnr * rsqinv;
        vn2 = dr[1]*vnnr * rsqinv;
        vn3 = dr[2]*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        wr1 = (radi*omega[i][0]) * rinv;
        wr2 = (radi*omega[i][1]) * rinv;
        wr3 = (radi*omega[i][2]) * rinv;

        // meff = effective mass of sphere and line
        // if I or J is part of rigid body, use body mass
        // if line is not part of rigid body assume infinite mass

        meff = rmass[i];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) meff = mass_rigid[i];
          if (mass_rigid[j] > 0.0) {
            mline = mass_rigid[j];
            meff = meff*mline / (meff+mline);
          }
        }

        // normal forces = Hookian contact + normal velocity damping

        damp = meff*gamman*vnnr*rsqinv;
        ccel = kn*(radi-r)*rinv - damp;
        ccel *= factor_couple;

        // relative velocities

        vtr1 = vt1 - (dr[2]*wr2-dr[1]*wr3);
        vtr2 = vt2 - (dr[0]*wr3-dr[2]*wr1);
        vtr3 = vt3 - (dr[1]*wr1-dr[0]*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;
        shear = &allshear[3*jj];

        if (shearupdate) {
          shear[0] += vtr1*dt;
          shear[1] += vtr2*dt;
          shear[2] += vtr3*dt;
        }
        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
                      shear[2]*shear[2]);

        // rotate shear displacements

        rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
        rsht *= rsqinv;
        if (shearupdate) {
          shear[0] -= rsht*delx;
          shear[1] -= rsht*dely;
          shear[2] -= rsht*delz;
        }

        // tangential force due to tangential velocity damping

        fs1 = - (kt*shear[0] + meff*gammat*vtr1) * factor_couple;;
        fs2 = - (kt*shear[1] + meff*gammat*vtr2) * factor_couple;;
        fs3 = - (kt*shear[2] + meff*gammat*vtr3) * factor_couple;;

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = xmu * fabs(ccel*r);

        if (fs > fn) {
          if (shrmag != 0.0) {
            shear[0] = (fn/fs) * (shear[0] + meff*gammat*vtr1/kt) -
              meff*gammat*vtr1/kt;
            shear[1] = (fn/fs) * (shear[1] + meff*gammat*vtr2/kt) -
              meff*gammat*vtr2/kt;
            shear[2] = (fn/fs) * (shear[2] + meff*gammat*vtr3/kt) -
              meff*gammat*vtr3/kt;
            fs1 *= fn/fs;
            fs2 *= fn/fs;
            fs3 *= fn/fs;
          } else fs1 = fs2 = fs3 = 0.0;
        }

        // total force on sphere

        fx = dr[0]*ccel + fs1;
        fy = dr[1]*ccel + fs2;
        fz = dr[2]*ccel + fs3;

        // sphere force & torque

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        tor1 = rinv * (dr[1]*fs3 - dr[2]*fs2);
        tor2 = rinv * (dr[2]*fs1 - dr[0]*fs3);
        tor3 = rinv * (dr[0]*fs2 - dr[1]*fs1);
        torque[i][0] -= radi*tor1;
        torque[i][1] -= radi*tor2;
        torque[i][2] -= radi*tor3;

        // line force & torque
        // torque applied at contact pt
        // use total force for torque
        //   since opposite force is not norm/tang to line at its end pt

        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;

        torque[j][0] -= ds[1]*fz - ds[2]*fy;
        torque[j][1] -= ds[2]*fx - ds[0]*fz;
        torque[j][2] -= ds[0]*fy - ds[1]*fx;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 0.0,0.0,fx,fy,fz,dr[0],dr[1],dr[2]);
      }

      /*
      if (touch[jj]) {
        shear = &allshear[3*jj];
        printf("POST step %ld ij %d %d shear %g %g %g\n",
               update->ntimestep,atom->tag[i],tag[j],
               shear[0],shear[1],shear[2]);
      }
      */

    }
  }

  // NOTE: should there be virial contributions from boundary tris?

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::settings(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal pair_style command");

  kn = utils::numeric(FLERR,arg[0],false,lmp);
  if (strcmp(arg[1],"NULL") == 0)
    kt = kn * 2.0/7.0;
  else
    kt = utils::numeric(FLERR,arg[1],false,lmp);

  gamman = utils::numeric(FLERR,arg[2],false,lmp);
  if (strcmp(arg[3],"NULL") == 0)
    gammat = 0.5 * gamman;
  else
    gammat = utils::numeric(FLERR,arg[3],false,lmp);

  xmu = utils::numeric(FLERR,arg[4],false,lmp);
  dampflag = utils::inumeric(FLERR,arg[5],false,lmp);
  if (dampflag == 0) gammat = 0.0;

  limit_damping = 0;
  if (narg == 7) {
    if (strcmp(arg[6], "limit_damping") == 0)
      limit_damping = 1;
    else
      error->all(FLERR, "Illegal pair_style command");
  }

  if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 ||
      xmu < 0.0 || xmu > 10000.0 || dampflag < 0 || dampflag > 1)
    error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::coeff(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::init_style()
{
  int i;

  dt = update->dt;

  // error and warning checks

  avec = (AtomVecLine *) atom->style_match("line");
  if (!avec) error->all(FLERR,"Pair line/gran requires atom style line");
  if (!atom->radius_flag || !atom->rmass_flag)
    error->all(FLERR, "Pair line/gran requires atom attributes radius, rmass");
  if (!force->newton_pair)
    error->all(FLERR,"Pair style line/gran requires newton pair on");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair line/gran requires ghost atoms store velocity");

  // need a granular neighbor list

  if (history)
    neighbor->add_request(this,NeighConst::REQ_SIZE | NeighConst::REQ_ONESIDED |
                          NeighConst::REQ_HISTORY);
  else
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_ONESIDED);

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (history && (fix_history == nullptr)) {
    auto cmd = fmt::format("NEIGH_HISTORY_HHLINE{} all NEIGH_HISTORY {} onesided", instance_me, size_history);
    fix_history = dynamic_cast<FixNeighHistory *>(
        modify->replace_fix("NEIGH_HISTORY_HHLINE_DUMMY" + std::to_string(instance_me), cmd, 1));
    fix_history->pair = this;
  }

  // set ptr to FixSurfaceLocal for surf connectivity info

  fsl = nullptr;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style,"surface/local") == 0) {
      if (fsl)
        error->all(FLERR,"Pair line/gran requires single fix surface/local");
      fsl = (FixSurfaceLocal *) modify->fix[m];
    }
  }
  if (!fsl) error->all(FLERR,"Pair line/gran requires a fix surface/local");

  // surfmoveflag = 1 if surfs may move at every step
  // yes if fix move exists and its group includes lines
  // NOTE: are there other conditions, like fix deform or fix npt?

  surfmoveflag = 0;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style,"move") == 0) {
      int groupbit = modify->fix[m]->groupbit;
      int *line = atom->line;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;
      int flag = 0;
      for (int i = 0; i < nlocal; i++) {
        if (line[i] < 0) continue;
        if (mask[i] & groupbit) flag = 1;
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
  // lines cannot be frozen

  int itype;
  for (i = 1; i <= atom->ntypes; i++) {
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
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (line[i] >= 0)
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

  // set fix which stores history info

  if (history) {
    fix_history = dynamic_cast<FixNeighHistory *>(
        modify->get_fix_by_id("NEIGH_HISTORY_HHLINE" + std::to_string(instance_me)));
    if (!fix_history) error->all(FLERR, "Could not find pair fix neigh history ID");
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLineGranHookeHistory::init_one(int i, int j)
{
  if (!allocated) allocate();

  // cutoff = sum of max I,J radii for
  // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

  double cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
  cutoff = MAX(cutoff,maxrad_frozen[i]+maxrad_dynamic[j]);
  cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]);
  return cutoff;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fwrite(&setflag[i][j],sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0)
        utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::write_restart_settings(FILE *fp)
{
  fwrite(&kn,sizeof(double),1,fp);
  fwrite(&kt,sizeof(double),1,fp);
  fwrite(&gamman,sizeof(double),1,fp);
  fwrite(&gammat,sizeof(double),1,fp);
  fwrite(&xmu,sizeof(double),1,fp);
  fwrite(&dampflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&kn,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&kt,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&gamman,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&gammat,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&xmu,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&dampflag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&kn,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kt,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamman,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&gammat,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&xmu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&dampflag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void PairLineGranHookeHistory::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

int PairLineGranHookeHistory::pack_forward_comm(int n, int *list, double *buf,
                                                int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mass_rigid[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairLineGranHookeHistory::unpack_forward_comm(int n, int first,
                                                   double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  for (i = first; i < last; i++)
    mass_rigid[i] = buf[m++];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairLineGranHookeHistory::memory_usage()
{
  double bytes = nmax * sizeof(double);
  bytes = emax*4 * sizeof(double);        // endpts array for line particles
  return bytes;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// private methods specific to PairLineGranHookeHistory
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   compute current end points of owned and ghost lines
   nothing computed for particles that are not lines
------------------------------------------------------------------------- */

void PairLineGranHookeHistory::calculate_endpts()
{
  int i,m;
  double length,theta,dx,dy;
  double *endpt;

  // realloc endpts array if necssary

  if (fsl->nmax_connect > emax) {
    memory->destroy(endpts);
    emax = fsl->nmax_connect;
    memory->create(endpts,emax,4,"line/gran:endpts");
  }

  AtomVecLine::Bonus *bonus = avec->bonus;
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

int PairLineGranHookeHistory::overlap_sphere_line(int i, int j, double *pt,
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

int PairLineGranHookeHistory::endpt_neigh_check(int i, int j, int jflag)
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
  // if all lines have endpt overlap, line with lowest ID computes
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
    if (k < 0) error->one(FLERR,"Pair line/gran neighbor line is missing");
    kflag = overlap_sphere_line(i,k,contact,dr,rsq);
    if (kflag > 0) return 1;
    if (kflag == 0) error->one(FLERR,"Fix surface/global neighbor line overlap is invalid");
    linemin = MIN(linemin,tag[k]);
  }

  if (tag[j] == linemin) return 0;
  return 1;
}

// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#include "fix_nh_uef.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "compute_pressure_uef.h"
#include "compute_temp_uef.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "irregular.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "timer.h"
#include "uef_utils.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{ISO,ANISO,TRICLINIC};

// citation info

static const char cite_user_uef_package[] =
  "UEF package: doi:10.1063/1.4972894\n\n"
  "@Article{NicholsonRutledge16,\n"
  "author = {David A. Nicholson and Gregory C. Rutledge},\n"
  "title = {Molecular Simulation of Flow-Enhanced Nucleation in\n"
  "   {$n$}-Eicosane Melts Under Steady Shear and Uniaxial Extension},\n"
  "journal = {The Journal of Chemical Physics},\n"
  "volume = {145},\n"
  "number = {24},\n"
  "pages = {244903},\n"
  "doi = {10.1063/1.4972894},\n"
  "year = {2016}\n"
  "}\n\n";

/* ----------------------------------------------------------------------
 * Parse fix specific keywords, do some error checking, and initialize
 * temp/pressure fixes
 ---------------------------------------------------------------------- */
FixNHUef::FixNHUef(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg), uefbox(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_user_uef_package);

  //initialization

  erate[0] = erate[1] = 0;

  // default values

  strain[0]=strain[1]= 0;
  ext_flags[0]=ext_flags[1]=ext_flags[2] = true;

  // need to initialize these

  omega_dot[0]=omega_dot[1]=omega_dot[2]=0;

  // parse fix nh/uef specific options

  bool erate_flag = false;
  int iarg = 3;

  while (iarg <narg) {
    if (strcmp(arg[iarg],"erate")==0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix nvt/npt/uef command");
      erate[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      erate[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      erate_flag = true;
      iarg += 3;
    } else if (strcmp(arg[iarg],"strain")==0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix nvt/npt/uef command");
      strain[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      strain[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg],"ext")==0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nvt/npt/uef command");
      if (strcmp(arg[iarg+1],"x")==0)
        ext_flags[1] = ext_flags[2] =  false;
      else if (strcmp(arg[iarg+1],"y")==0)
        ext_flags[0] = ext_flags[2] =  false;
      else if (strcmp(arg[iarg+1],"z")==0)
        ext_flags[0] = ext_flags[1] =  false;
      else if (strcmp(arg[iarg+1],"xy")==0)
        ext_flags[2] = false;
      else if (strcmp(arg[iarg+1],"xz")==0)
        ext_flags[1] = false;
      else if (strcmp(arg[iarg+1],"yz")==0)
        ext_flags[0] = false;
      else if (strcmp(arg[iarg+1],"xyz")!=0)
        error->all(FLERR,"Illegal fix nvt/npt/uef command");

      iarg += 2;
    } else {

      // skip to next argument; argument check for unknown keywords is done in FixNH

      ++iarg;
    }
  }

  if (!erate_flag)
    error->all(FLERR,"Keyword erate must be set for fix nvt/npt/uef command");

  if (mtchain_default_flag) mtchain=1;

  if (!domain->triclinic)
    error->all(FLERR,"Simulation box must be triclinic for fix/nvt/npt/uef");

  // check for conditions that impose a deviatoric stress

  if (pstyle == TRICLINIC)
    error->all(FLERR,"Only normal stresses can be controlled with fix/nvt/npt/uef");
  double erate_tmp[3];
  erate_tmp[0]=erate[0];
  erate_tmp[1]=erate[1];
  erate_tmp[2]=-erate[0]-erate[1];

  if (pstyle == ANISO) {
    if (!(ext_flags[0] & ext_flags[1] & ext_flags[2]))
      error->all(FLERR,"The ext keyword may only be used with iso pressure control");
    for (int k=0;k<3;k++)
      for (int j=0;j<3;j++)
        if (p_flag[k] && p_flag[j]) {
          double tol = 1e-6;
          if ( !nearly_equal(p_start[k],p_start[j],tol)
               || !nearly_equal(p_stop[k],p_stop[j],tol))
            error->all(FLERR,"All controlled stresses must have the same "
                       "value in fix/nvt/npt/uef");
          if ( !nearly_equal(erate_tmp[k],erate_tmp[j],tol)
               || !nearly_equal(erate_tmp[k],erate_tmp[j],tol))
            error->all(FLERR,"Dimensions with controlled stresses must have"\
                       " same strain rate in fix/nvt/npt/uef");
        }
  }

  // conditions that produce a deviatoric stress have already been eliminated.

  deviatoric_flag=0;

  // need pre_exchange and irregular migration

  pre_exchange_flag = 1;
  irregular = new Irregular(lmp);

  // flag that I change the box here (in case of nvt)

  box_change |= BOX_CHANGE_SHAPE;

  // initialize the UEFBox class which computes the box at each step

  uefbox = new UEF_utils::UEFBox();
  uefbox->set_strain(strain[0],strain[1]);

  // reset fixedpoint to the stagnation point. I don't allow fixedpoint
  // to be set by the user.

  fixedpoint[0] = domain->boxlo[0];
  fixedpoint[1] = domain->boxlo[1];
  fixedpoint[2] = domain->boxlo[2];

  // Create temp and pressure computes for nh/uef

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp/uef",id_temp));
  tcomputeflag = 1;

  id_press = utils::strdup(std::string(id) + "_press");
  modify->add_compute(fmt::format("{} all pressure/uef {}",id_press, id_temp));
  pcomputeflag = 1;

  nevery = 1;
}

/* ----------------------------------------------------------------------
 * Erase the UEFBox object and get rid of the pressure compute if the nvt
 * version is being used. Everything else will be done in base destructor
 * ---------------------------------------------------------------------- */
FixNHUef::~FixNHUef()
{
  delete uefbox;
  if (pcomputeflag && !pstat_flag)  {
    modify->delete_compute(id_press);
    delete [] id_press;
  }
}

/* ----------------------------------------------------------------------
 * Make the end_of_step() routine callable
 * ---------------------------------------------------------------------- */
int FixNHUef::setmask()
{
  int mask = FixNH::setmask();
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
 * Run FixNH::init() and do more error checking. Set the pressure
 * pointer in the case that the nvt version is used
 * ---------------------------------------------------------------------- */
void FixNHUef::init()
{
  FixNH::init();


  // find conflict with fix/deform or other box chaging fixes
  for (int i=0; i < modify->nfix; i++)
  {
    if (strcmp(modify->fix[i]->id,id) != 0)
      if ((modify->fix[i]->box_change & BOX_CHANGE_SHAPE) != 0)
        error->all(FLERR,"Can't use another fix which changes box shape with fix/nvt/npt/uef");
  }


  // this will make the pressure compute for nvt
  if (!pstat_flag)
    if (pcomputeflag) {
      int icomp = modify->find_compute(id_press);
      if (icomp<0)
        error->all(FLERR,"Pressure ID for fix/nvt/uef doesn't exist");
      pressure = modify->compute[icomp];

      if (strcmp(pressure->style,"pressure/uef") != 0)
        error->all(FLERR,"Using fix nvt/npt/uef without a compute pressure/uef");
    }

  if (strcmp(temperature->style,"temp/uef") != 0)
    error->all(FLERR,"Using fix nvt/npt/uef without a compute temp/uef");
}

/* ----------------------------------------------------------------------
 * Run FixNH::setup() make sure the box is OK and set the rotation matrix
 * for the first step
 * ---------------------------------------------------------------------- */
void FixNHUef::setup(int j)
{
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  double tol = 1e-4;
  // ensure the box is ok for uef
  bool isok = true;
  isok &= nearly_equal(domain->h[0],box[0][0],tol);
  isok &= nearly_equal(domain->h[1],box[1][1],tol);
  isok &= nearly_equal(domain->h[2],box[2][2],tol);
  isok &= nearly_equal(domain->xy,box[0][1],tol);
  isok &= nearly_equal(domain->yz,box[1][2],tol);
  isok &= nearly_equal(domain->xz,box[0][2],tol);
  if (!isok)
    error->all(FLERR,"Initial box is not close enough to the expected uef box");

  uefbox->get_rot(rot);
  (dynamic_cast<ComputeTempUef*>(temperature))->yes_rot();
  (dynamic_cast<ComputePressureUef*>(pressure))->in_fix = true;
  (dynamic_cast<ComputePressureUef*>(pressure))->update_rot();
  FixNH::setup(j);
}

/* ----------------------------------------------------------------------
 * rotate -> initial integration step -> rotate back
 * ---------------------------------------------------------------------- */
void FixNHUef::initial_integrate(int vflag)
{
  inv_rotate_x(rot);
  inv_rotate_v(rot);
  inv_rotate_f(rot);
  (dynamic_cast<ComputeTempUef*>(temperature))->no_rot();
  FixNH::initial_integrate(vflag);
  rotate_x(rot);
  rotate_v(rot);
  rotate_f(rot);
  (dynamic_cast<ComputeTempUef*>(temperature))->yes_rot();
}

/* ----------------------------------------------------------------------
 * rotate -> initial integration step -> rotate back (RESPA)
 * ---------------------------------------------------------------------- */
void FixNHUef::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  inv_rotate_x(rot);
  inv_rotate_v(rot);
  inv_rotate_f(rot);
  (dynamic_cast<ComputeTempUef*>(temperature))->no_rot();
  FixNH::initial_integrate_respa(vflag,ilevel,iloop);
  rotate_x(rot);
  rotate_v(rot);
  rotate_f(rot);
  (dynamic_cast<ComputeTempUef*>(temperature))->yes_rot();
}

/* ----------------------------------------------------------------------
 * rotate -> final integration step -> rotate back
 * ---------------------------------------------------------------------- */
void FixNHUef::final_integrate()
{
  // update rot here since it must directly follow the virial calculation
  (dynamic_cast<ComputePressureUef*>(pressure))->update_rot();
  inv_rotate_v(rot);
  inv_rotate_f(rot);
  (dynamic_cast<ComputeTempUef*>(temperature))->no_rot();
  FixNH::final_integrate();
  rotate_v(rot);
  rotate_f(rot);
  (dynamic_cast<ComputeTempUef*>(temperature))->yes_rot();
}

/* ----------------------------------------------------------------------
 * at outer level: call this->final_integrate()
 * at other levels: rotate -> 2nd verlet step -> rotate back
 * ---------------------------------------------------------------------- */
void FixNHUef::final_integrate_respa(int ilevel, int /*iloop*/)
{
  // set timesteps by level
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];
  // outermost level - update eta_dot and omega_dot, apply via final_integrate
  // all other levels - NVE update of v
  if (ilevel == nlevels_respa-1) final_integrate();
  else
  {
    inv_rotate_v(rot);
    inv_rotate_f(rot);
    nve_v();
    rotate_v(rot);
    rotate_f(rot);
  }
}

/* ----------------------------------------------------------------------
   SLLOD velocity update in time-reversible (i think) increments
   v -> exp(-edot*dt/2)*v
   v -> v +f/m*dt
   v -> exp(-edot*dt/2)*v
-----------------------------------------------------------------------*/
void FixNHUef::nve_v()
{
  double dtfm;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ex = erate[0]*dtf/2;
  double ey = erate[1]*dtf/2;
  double ez = -ex-ey;
  double e0 = exp(-ex);
  double e1 = exp(-ey);
  double e2 = exp(-ez);
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
        v[i][0] += dtfm*f[i][0];
        v[i][1] += dtfm*f[i][1];
        v[i][2] += dtfm*f[i][2];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
        v[i][0] += dtfm*f[i][0];
        v[i][1] += dtfm*f[i][1];
        v[i][2] += dtfm*f[i][2];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Don't actually move atoms in remap(), just change the box
-----------------------------------------------------------------------*/
void FixNHUef::remap()
{
  double vol = domain->xprd * domain->yprd * domain->zprd;
  double domega = dto*(omega_dot[0]+omega_dot[1]+omega_dot[2])/3.;

  // constant volume strain associated with barostat
  // box scaling
  double ex = dto*omega_dot[0]-domega;
  double ey = dto*omega_dot[1]-domega;
  uefbox->step_deform(ex,ey);
  strain[0] += ex;
  strain[1] += ey;

  // volume change
  vol = vol*exp(3*domega);
  double box[3][3];
  uefbox->get_box(box,vol);
  domain->boxhi[0] = domain->boxlo[0]+box[0][0];
  domain->boxhi[1] = domain->boxlo[1]+box[1][1];
  domain->boxhi[2] = domain->boxlo[2]+box[2][2];
  domain->xy = box[0][1];
  domain->xz = box[0][2];
  domain->yz = box[1][2];
  domain->set_global_box();
  domain->set_local_box();
  uefbox->get_rot(rot);
}

/* ----------------------------------------------------------------------
   SLLOD position update in time-reversible (i think) increments
   x -> exp(edot*dt/2)*x
   x -> x + v*dt
   x -> exp(edot*dt/2)*x
-----------------------------------------------------------------------*/
void FixNHUef::nve_x()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ex = erate[0]*dtv;
  strain[0] += ex;
  double e0 = exp((ex+omega_dot[0]*dtv)/2);
  double ey = erate[1]*dtv;
  strain[1] += ey;
  double e1 = exp((ey+omega_dot[1]*dtv)/2.);
  double ez = -ex -ey;
  double e2 = exp((ez+omega_dot[2]*dtv)/2.);
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x update by full step only for atoms in group
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] *= e0;
      x[i][1] *= e1;
      x[i][2] *= e2;
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
      x[i][0] *= e0;
      x[i][1] *= e1;
      x[i][2] *= e2;
    }
  }
  uefbox->step_deform(ex,ey);
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  domain->boxhi[0] = domain->boxlo[0]+box[0][0];
  domain->boxhi[1] = domain->boxlo[1]+box[1][1];
  domain->boxhi[2] = domain->boxlo[2]+box[2][2];
  domain->xy = box[0][1];
  domain->xz = box[0][2];
  domain->yz = box[1][2];
  domain->set_global_box();
  domain->set_local_box();
  uefbox->get_rot(rot);
}

/* ----------------------------------------------------------------------
 * Do the lattice reduction if necessary.
-----------------------------------------------------------------------*/
void FixNHUef::pre_exchange()
{
  // only need to reset things if the lattice needs to be reduced
  if (uefbox->reduce())
  {
    // go to lab frame
    inv_rotate_x(rot);
    inv_rotate_v(rot);
    inv_rotate_f(rot);
    // get & set the new box and rotation matrix
    double vol = domain->xprd * domain->yprd * domain->zprd;
    double box[3][3];
    uefbox->get_box(box,vol);
    domain->boxhi[0] = domain->boxlo[0]+box[0][0];
    domain->boxhi[1] = domain->boxlo[1]+box[1][1];
    domain->boxhi[2] = domain->boxlo[2]+box[2][2];
    domain->xy = box[0][1];
    domain->xz = box[0][2];
    domain->yz = box[1][2];
    domain->set_global_box();
    domain->set_local_box();
    uefbox->get_rot(rot);

    // rotate to the new upper triangular frame
    rotate_v(rot);
    rotate_x(rot);
    rotate_f(rot);

    // this is a generalization of what is done in domain->image_flip(...)
    int ri[3][3];
    uefbox->get_inverse_cob(ri);
    imageint *image = atom->image;
    int nlocal = atom->nlocal;
    for (int i=0; i<nlocal; i++) {
      int iold[3],inew[3];
      iold[0] = (image[i] & IMGMASK) - IMGMAX;
      iold[1] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      iold[2] = (image[i] >> IMG2BITS) - IMGMAX;
      inew[0] = ri[0][0]*iold[0] + ri[0][1]*iold[1] + ri[0][2]*iold[2];
      inew[1] = ri[1][0]*iold[0] + ri[1][1]*iold[1] + ri[1][2]*iold[2];
      inew[2] = ri[2][0]*iold[0] + ri[2][1]*iold[1] + ri[2][2]*iold[2];
      image[i] = ((imageint) (inew[0] + IMGMAX) & IMGMASK) |
        (((imageint) (inew[1] + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (inew[2] + IMGMAX) & IMGMASK) << IMG2BITS);
    }

    // put all atoms in the new box
    double **x = atom->x;
    for (int i=0; i<nlocal; i++) domain->remap(x[i],image[i]);

    // move atoms to the right processors
    domain->x2lamda(atom->nlocal);
    irregular->migrate_atoms();
    domain->lamda2x(atom->nlocal);
  }
}

/* ----------------------------------------------------------------------
 * The following are routines to rotate between the lab and upper triangular
 * (UT) frames. For most of the time the simulation is in the UT frame.
 * To get to the lab frame, apply the inv_rotate_[..](rot) and to
 * get back to the UT frame apply rotate_[..](rot).
 *
 * Note: the rotate_x() functions also apply a shift to/from the fixedpoint
 * to make the integration a little simpler.
 * ---------------------------------------------------------------------- */
void FixNHUef::rotate_x(double r[3][3])
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double xn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      xn[0]=r[0][0]*x[i][0]+r[0][1]*x[i][1]+r[0][2]*x[i][2];
      xn[1]=r[1][0]*x[i][0]+r[1][1]*x[i][1]+r[1][2]*x[i][2];
      xn[2]=r[2][0]*x[i][0]+r[2][1]*x[i][1]+r[2][2]*x[i][2];
      x[i][0]=xn[0]+domain->boxlo[0];
      x[i][1]=xn[1]+domain->boxlo[1];
      x[i][2]=xn[2]+domain->boxlo[2];
    }
  }
}

void FixNHUef::inv_rotate_x(double r[3][3])
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double xn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      x[i][0] -= domain->boxlo[0];
      x[i][1] -= domain->boxlo[1];
      x[i][2] -= domain->boxlo[2];
      xn[0]=r[0][0]*x[i][0]+r[1][0]*x[i][1]+r[2][0]*x[i][2];
      xn[1]=r[0][1]*x[i][0]+r[1][1]*x[i][1]+r[2][1]*x[i][2];
      xn[2]=r[0][2]*x[i][0]+r[1][2]*x[i][1]+r[2][2]*x[i][2];
      x[i][0]=xn[0];
      x[i][1]=xn[1];
      x[i][2]=xn[2];
    }
  }
}

void FixNHUef::rotate_v(double r[3][3])
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double vn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      vn[0]=r[0][0]*v[i][0]+r[0][1]*v[i][1]+r[0][2]*v[i][2];
      vn[1]=r[1][0]*v[i][0]+r[1][1]*v[i][1]+r[1][2]*v[i][2];
      vn[2]=r[2][0]*v[i][0]+r[2][1]*v[i][1]+r[2][2]*v[i][2];
      v[i][0]=vn[0]; v[i][1]=vn[1]; v[i][2]=vn[2];
    }
  }
}

void FixNHUef::inv_rotate_v(double r[3][3])
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double vn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      vn[0]=r[0][0]*v[i][0]+r[1][0]*v[i][1]+r[2][0]*v[i][2];
      vn[1]=r[0][1]*v[i][0]+r[1][1]*v[i][1]+r[2][1]*v[i][2];
      vn[2]=r[0][2]*v[i][0]+r[1][2]*v[i][1]+r[2][2]*v[i][2];
      v[i][0]=vn[0]; v[i][1]=vn[1]; v[i][2]=vn[2];
    }
  }
}

void FixNHUef::rotate_f(double r[3][3])
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double fn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      fn[0]=r[0][0]*f[i][0]+r[0][1]*f[i][1]+r[0][2]*f[i][2];
      fn[1]=r[1][0]*f[i][0]+r[1][1]*f[i][1]+r[1][2]*f[i][2];
      fn[2]=r[2][0]*f[i][0]+r[2][1]*f[i][1]+r[2][2]*f[i][2];
      f[i][0]=fn[0]; f[i][1]=fn[1]; f[i][2]=fn[2];
    }
  }
}

void FixNHUef::inv_rotate_f(double r[3][3])
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  double fn[3];
  for (int i=0;i<nlocal;i++)
  {
    if (mask[i] & groupbit)
    {
      fn[0]=r[0][0]*f[i][0]+r[1][0]*f[i][1]+r[2][0]*f[i][2];
      fn[1]=r[0][1]*f[i][0]+r[1][1]*f[i][1]+r[2][1]*f[i][2];
      fn[2]=r[0][2]*f[i][0]+r[1][2]*f[i][1]+r[2][2]*f[i][2];
      f[i][0]=fn[0]; f[i][1]=fn[1]; f[i][2]=fn[2];
    }
  }
}

/* ----------------------------------------------------------------------
 * Increase the size of the restart list to add in the strains
 * ---------------------------------------------------------------------- */
int FixNHUef::size_restart_global()
{
  return FixNH::size_restart_global() +2;
}

/* ----------------------------------------------------------------------
 * Pack the strains after packing the default FixNH values
 * ---------------------------------------------------------------------- */
int FixNHUef::pack_restart_data(double *list)
{
  int n = FixNH::pack_restart_data(list);
  list[n++] = strain[0];
  list[n++] = strain[1];
  return n;
}

/* ----------------------------------------------------------------------
 * read and set the strains after the default FixNH values
 * ---------------------------------------------------------------------- */
void FixNHUef::restart(char *buf)
{
  int n = size_restart_global();
  FixNH::restart(buf);
  auto list = (double *) buf;
  strain[0] = list[n-2];
  strain[1] = list[n-1];
  uefbox->set_strain(strain[0],strain[1]);
}

/* ----------------------------------------------------------------------
 * If the step writes a restart, reduce the box beforehand. This makes sure
 * the unique box shape can be found once the restart is read and that
 * all of the atoms lie within the box.
 * This may only be necessary for RESPA runs, but I'm leaving it in anyway.
 * ---------------------------------------------------------------------- */
void FixNHUef::end_of_step()
{
  if (update->ntimestep==output->next_restart)
  {
    pre_exchange();
    domain->x2lamda(atom->nlocal);
    domain->pbc();
    timer->stamp();
    comm->exchange();
    comm->borders();
    domain->lamda2x(atom->nlocal+atom->nghost);
    timer->stamp(Timer::COMM);
    neighbor->build(1);
    timer->stamp(Timer::NEIGH);
  }
}

/* ----------------------------------------------------------------------
 * reduce the simulation box after a run is complete. otherwise it won't
 * be possible to resume from a write_restart since the initialization of
 * the simulation box requires reduced simulation box
 * ---------------------------------------------------------------------- */
void FixNHUef::post_run()
{
  pre_exchange();
  domain->x2lamda(atom->nlocal);
  domain->pbc();
  timer->stamp();
  comm->exchange();
  comm->borders();
  domain->lamda2x(atom->nlocal+atom->nghost);
  timer->stamp(Timer::COMM);
  neighbor->build(1);
  timer->stamp(Timer::NEIGH);
}

/* ----------------------------------------------------------------------
 * public read for rotation matrix
 * ---------------------------------------------------------------------- */
void FixNHUef::get_rot(double r[3][3])
{
  r[0][0] = rot[0][0];
  r[0][1] = rot[0][1];
  r[0][2] = rot[0][2];
  r[1][0] = rot[1][0];
  r[1][1] = rot[1][1];
  r[1][2] = rot[1][2];
  r[2][0] = rot[2][0];
  r[2][1] = rot[2][1];
  r[2][2] = rot[2][2];
}

/* ----------------------------------------------------------------------
 * public read for ext flags
 * ---------------------------------------------------------------------- */
void FixNHUef::get_ext_flags(bool* e)
{
  e[0] = ext_flags[0];
  e[1] = ext_flags[1];
  e[2] = ext_flags[2];
}

/* ----------------------------------------------------------------------
 * public read for simulation box
 * ---------------------------------------------------------------------- */
void FixNHUef::get_box(double b[3][3])
{
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  b[0][0] = box[0][0];
  b[0][1] = box[0][1];
  b[0][2] = box[0][2];
  b[1][0] = box[1][0];
  b[1][1] = box[1][1];
  b[1][2] = box[1][2];
  b[2][0] = box[2][0];
  b[2][1] = box[2][1];
  b[2][2] = box[2][2];
}

/* ----------------------------------------------------------------------
 * comparing floats
 * it's imperfect, but should work provided no infinities
 * ---------------------------------------------------------------------- */
bool FixNHUef::nearly_equal(double a, double b, double epsilon)
{
  double absa = fabs(a);
  double absb = fabs(b);
  double diff = fabs(a-b);
  if (a == b) return true;
  else if ( (absa+absb) < epsilon)
    return diff < epsilon*epsilon;
  else
    return diff/(absa+absb) < epsilon;
}

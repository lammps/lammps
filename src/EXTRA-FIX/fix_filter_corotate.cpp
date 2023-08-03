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
     Contributing author: Lukas Fath (KIT)
     some subroutines are from fix_shake.cpp
   ------------------------------------------------------------------------- */

#include "fix_filter_corotate.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

#include <cctype>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

#define BIG 1.0e20
#define MASSDELTA 0.1

static const char cite_filter_corotate[] =
  "Mollified Impulse Method with Corotational Filter: doi:10.1016/j.jcp.2016.12.024\n\n"
  "@Article{Fath2017,\n"
  " Title ="
  "{A Fast Mollified Impulse Method for Biomolecular Atomistic Simulations},\n"
  " Author = {L. Fath and M. Hochbruck and C. V. Singh},\n"
  " Journal = {Journal of Computational Physics},\n"
  " Year = {2017},\n"
  " Pages = {180--198},\n"
  " Volume = {333},\n\n"
  " Doi = {https://doi.org/10.1016/j.jcp.2016.12.024},\n"
  " ISSN = {0021-9991},\n"
  " Keywords = {Mollified impulse method},\n"
  " Url = {https://www.sciencedirect.com/science/article/pii/S0021999116306787}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixFilterCorotate::FixFilterCorotate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp,narg,arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_filter_corotate);

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  molecular = atom->molecular;
  if (molecular == Atom::ATOMIC)
    error->all(FLERR,"Cannot use fix filter/corotate "
      "with non-molecular system");
  if (molecular == Atom::TEMPLATE)
    error->all(FLERR,"Cannot use fix filter/corotate "
      "with molecular template system");

  // parse args for bond and angle types
  // will be used by find_clusters
  // store args for "b" "a" "t" as flags in (1:n) list for fast access
  // store args for "m" in list of length nmass for looping over
  // for "m" verify that atom masses have been set

  bond_flag = new int[atom->nbondtypes+1];
  for (int i = 1; i <= atom->nbondtypes; i++) bond_flag[i] = 0;
  angle_flag = new int[atom->nangletypes+1];
  for (int i = 1; i <= atom->nangletypes; i++) angle_flag[i] = 0;
  type_flag = new int[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++) type_flag[i] = 0;
  mass_list = new double[atom->ntypes];
  nmass = 0;

  char mode = '\0';
  int next = 3;
  while (next < narg) {
    if (strcmp(arg[next],"b") == 0) mode = 'b';
    else if (strcmp(arg[next],"a") == 0) mode = 'a';
    else if (strcmp(arg[next],"t") == 0) mode = 't';
    else if (strcmp(arg[next],"m") == 0) {
      mode = 'm';
      atom->check_mass(FLERR);

      // break if keyword that is not b,a,t,m

    } else if (isalpha(arg[next][0])) break;

    // read numeric args of b,a,t,m

    else if (mode == 'b') {
      int i = utils::inumeric(FLERR,arg[next],false,lmp);
      if (i < 1 || i > atom->nbondtypes)
        error->all(FLERR,"Invalid bond type index for fix filter/corotate");
      bond_flag[i] = 1;

    } else if (mode == 'a') {
      int i = utils::inumeric(FLERR,arg[next],false,lmp);
      if (i < 1 || i > atom->nangletypes)
        error->all(FLERR,"Invalid angle type index for fix filter/corotate");
      angle_flag[i] = 1;

    } else if (mode == 't') {
      int i = utils::inumeric(FLERR,arg[next],false,lmp);
      if (i < 1 || i > atom->ntypes)
        error->all(FLERR,"Invalid atom type index for fix filter/corotate");
      type_flag[i] = 1;

    } else if (mode == 'm') {
      double massone = utils::numeric(FLERR,arg[next],false,lmp);
      if (massone == 0.0)
        error->all(FLERR,"Invalid atom mass for fix filter/corotate");
      if (nmass == atom->ntypes)
        error->all(FLERR,"Too many masses for fix filter/corotate");
      mass_list[nmass++] = massone;

    } else error->all(FLERR,"Illegal fix filter/corotate command");
    next++;
  }


  // parse optional args
  nlevels_respa = 1; //initialize

  // allocate bond and angle distance arrays, indexed from 1 to n

  bond_distance = new double[atom->nbondtypes+1];
  angle_distance = new double[atom->nangletypes+1];

  //grow_arrays
  array_atom = nullptr;
  shake_flag = nullptr;
  shake_atom = nullptr;
  shake_type = nullptr;

  FixFilterCorotate::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);    //calls grow_arrays

  x_store = nullptr;

  //STUFF
  g = nullptr;
  help2 = nullptr;

  dn1dx = dn2dx = dn3dx = nullptr;
  n1 = n2 = n3 = del1 = del2 = del3 = nullptr;

  memory->grow(help2,15,15,"FilterCorotate:help2");
  memory->grow(n1,3,"FilterCorotate:n1");
  memory->grow(n2,3,"FilterCorotate:n2");
  memory->grow(n3,3,"FilterCorotate:n3");
  memory->grow(del1,3,"FilterCorotate:del1");
  memory->grow(del2,3,"FilterCorotate:del2");
  memory->grow(del3,3,"FilterCorotate:del3");
  memory->grow(dn1dx,3,15,"FilterCorotate:dn1dx");
  memory->grow(dn2dx,3,15,"FilterCorotate:dn2dx");
  memory->grow(dn3dx,3,15,"FilterCorotate:dn3dx");

  memory->grow(g,15,"FilterCorotate:g");

  comm_forward = 3;

  // identify all clusters

  find_clusters();

  // initialize list of clusters to constrain

  maxlist = 0;
  list = nullptr;
  clist_derv = nullptr;
  clist_q0 = nullptr;    //list for derivative and ref. config
  clist_nselect1 = clist_nselect2 = nullptr;
  clist_select1 = clist_select2 = nullptr;

}

/* ---------------------------------------------------------------------- */

FixFilterCorotate::~FixFilterCorotate()
{
  memory->destroy(g);
  memory->destroy(help2);
  memory->destroy(n1);
  memory->destroy(n2);
  memory->destroy(n3);
  memory->destroy(del1);
  memory->destroy(del2);
  memory->destroy(del3);
  memory->destroy(dn1dx);
  memory->destroy(dn2dx);
  memory->destroy(dn3dx);

  atom->delete_callback(id,Atom::GROW);

  // delete locally stored arrays

  memory->destroy(shake_flag);
  memory->destroy(shake_atom);
  memory->destroy(shake_type);
  memory->destroy(array_atom);

  delete [] bond_flag;
  delete [] angle_flag;
  delete [] type_flag;
  delete [] mass_list;

  delete [] bond_distance;
  delete [] angle_distance;

  memory->destroy(list);
  memory->destroy(clist_derv);
  memory->destroy(clist_q0);
  memory->destroy(clist_nselect1);
  memory->destroy(clist_nselect2);
  memory->destroy(clist_select1);
  memory->destroy(clist_select2);
}


/* ---------------------------------------------------------------------- */

int FixFilterCorotate::setmask()
{
  int mask = 0;

  mask |= PRE_NEIGHBOR;
  mask |= PRE_FORCE_RESPA;
  mask |= POST_FORCE_RESPA;

  return mask;
}


/* ---------------------------------------------------------------------- */

void FixFilterCorotate::init()
{
  int i;
  // error if more than one filter
  int count = 0;
  for (i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"filter/corotate") == 0) count++;
  }
  if (count > 1) error->all(FLERR,"More than one fix filter/corotate");

  // check for fix shake:
  count = 0;
  for (i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"shake") == 0) count++;
  }
  if (count > 1)
    error->one(FLERR,"Both fix shake and fix filter/corotate detected.");

  // if rRESPA, find associated fix that must exist
  // could have changed locations in fix list since created
  // set ptrs to rRESPA variables

  if (utils::strmatch(update->integrate_style,"^respa")) {
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
  }
  else error->all(FLERR,"Fix filter/corotate requires rRESPA!");

  // set equilibrium bond distances

  if (force->bond == nullptr)
    error->all(FLERR,"Bond potential must be defined for fix filter/corotate");
  for (i = 1; i <= atom->nbondtypes; i++)
    bond_distance[i] = force->bond->equilibrium_distance(i);

  // set equilibrium angle value

  for (i = 1; i <= atom->nangletypes; i++) {
    angle_distance[i] = force->angle->equilibrium_angle(i);
  }
}

/* ----------------------------------------------------------------------
 *  pre-integrator setup
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::setup_pre_neighbor()
{
  pre_neighbor();
}


/* ----------------------------------------------------------------------
 *  build list of clusters to constrain
 *  if one or more atoms in cluster are on this proc,
 *  this proc lists the cluster exactly once
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::pre_neighbor()
{
  int atom1,atom2,atom3,atom4,atom5;

  // local copies of atom quantities
  // used until next re-neighboring

  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  rmass = atom->rmass;
  type = atom->type;
  nlocal = atom->nlocal;

  // extend size of list if necessary

  if (nlocal > maxlist) {
    maxlist = nlocal;

    memory->destroy(list);
    memory->destroy(clist_derv);
    memory->destroy(clist_q0);
    memory->destroy(clist_nselect1);
    memory->destroy(clist_nselect2);
    memory->destroy(clist_select1);
    memory->destroy(clist_select2);

    memory->create(list,maxlist,"FilterCorotate:list");
    memory->create(clist_derv,maxlist,15,15,"FilterCorotate:derivative");
    memory->create(clist_q0,maxlist,15,"FilterCorotate:q0");
    memory->create(clist_nselect1,maxlist,"FilterCorotate:nselect1");
    memory->create(clist_nselect2,maxlist,"FilterCorotate:nselect2");
    memory->create(clist_select1,maxlist,5,"FilterCorotate:select1");
    memory->create(clist_select2,maxlist,5,"FilterCorotate:select2");
  }

  // build list of clusters I compute

  nlist = 0;

  for (int i = 0; i < nlocal; i++) {
    if (shake_flag[i]) {
      if (shake_flag[i] == 2) {
        atom1 = atom->map(shake_atom[i][0]);
        atom2 = atom->map(shake_atom[i][1]);
        if (atom1 == -1 || atom2 == -1)
          error->one(FLERR,"Cluster atoms {} {} missing on proc {} at step {}",
                     shake_atom[i][0],shake_atom[i][1],me,update->ntimestep);
        if (i <= atom1 && i <= atom2) list[nlist++] = i;

      } else if (shake_flag[i]  == 1 || shake_flag[i]  == 3) {
        atom1 = atom->map(shake_atom[i][0]);
        atom2 = atom->map(shake_atom[i][1]);
        atom3 = atom->map(shake_atom[i][2]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1)
          error->one(FLERR,"Cluster atoms {} {} {} missing on proc {} at step {}",
                     shake_atom[i][0],shake_atom[i][1],shake_atom[i][2],me,update->ntimestep);
        if (i <= atom1 && i <= atom2 && i <= atom3) list[nlist++] = i;

      } else if (shake_flag[i]  == 4) {
        atom1 = atom->map(shake_atom[i][0]);
        atom2 = atom->map(shake_atom[i][1]);
        atom3 = atom->map(shake_atom[i][2]);
        atom4 = atom->map(shake_atom[i][3]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1)
          error->one(FLERR,"Cluster atoms {} {} {} {} missing on proc {} at step {}",
                     shake_atom[i][0],shake_atom[i][1],shake_atom[i][2],shake_atom[i][3],
                     me,update->ntimestep);
        if (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4)
          list[nlist++] = i;

      } else if (shake_flag[i]  == 5) {
        atom1 = atom->map(shake_atom[i][0]);
        atom2 = atom->map(shake_atom[i][1]);
        atom3 = atom->map(shake_atom[i][2]);
        atom4 = atom->map(shake_atom[i][3]);
        atom5 = atom->map(shake_atom[i][4]);
        if (atom1 == -1 || atom2 == -1 || atom3 == -1 || atom4 == -1 || atom5 == -1)
          error->one(FLERR,"Cluster atoms {} {} {} {} {} missing on proc {} at step {}",
                     shake_atom[i][0],shake_atom[i][1], shake_atom[i][2],shake_atom[i][3],
                     shake_atom[i][4],me,update->ntimestep);
        if (i <= atom1 && i <= atom2 && i <= atom3 && i <= atom4 && i <= atom5)
          list[nlist++] = i;
      }
    }
  }

  //also build reference config lists
  int m,N,oxy;
  double r0,r1,r2,theta0;
  double m1,m2,m3,m4,m5,m_all;

  double *mass = atom->mass;
  int *type = atom->type;

  for (int i=0; i<nlist; i++)
    {
      m = list[i];
      N = shake_flag[m];

      //switch cluster type 3 to angle cluster:
      if (N == 3)
        {
          //make it an angle cluster:
          if (shake_type[m][2] == 0)
            shake_type[m][2] = angletype_findset(atom->map(shake_atom[m][0]),
                                                 shake_atom[m][1],shake_atom[m][2],0);
        }

      if (N == 1) N = 3;  //angle cluster

      if (N == 2)    //cluster of size 2:
        {
          atom1 = atom->map(shake_atom[m][0]);
          atom2 = atom->map(shake_atom[m][1]);

          r0 = bond_distance[shake_type[m][0]];

          m1 = mass[type[atom1]];
          m2 = mass[type[atom2]];
          m_all = m1 + m2;

          double a1 = -m2/(m1+m2)*r0;
          double a2 =  m1/(m1+m2)*r0;

          clist_q0[i][0] = a1;
          clist_q0[i][1] = 0;
          clist_q0[i][2] = 0;

          clist_q0[i][3] = a2;
          clist_q0[i][4] = 0;
          clist_q0[i][5] = 0;

          clist_nselect1[i] = 1;
          clist_select1[i][0] = 1;

          clist_nselect2[i] = 0;
          clist_select2[i][0] = 0;
        } else if (N == 3)  //angle cluster
        {
          oxy = atom->map(shake_atom[m][0]);
          atom1 = atom->map(shake_atom[m][1]);
          atom2 = atom->map(shake_atom[m][2]);

          r0 = bond_distance[shake_type[m][0]];
          r1 = bond_distance[shake_type[m][1]];


          theta0 = angle_distance[shake_type[m][2]];

          m1 = mass[type[oxy]];
          m2 = mass[type[atom1]];
          m3 = mass[type[atom2]];
          m_all = m1 + m2 + m3;

          double alpha1 = 0.5*(MY_PI - theta0);

          double xcenter =  m2/m_all*r0*sin(alpha1) + m3/m_all*r1*sin(alpha1);
          double ycenter =  -m2/m_all*r0*cos(alpha1)+ m3/m_all*r1*cos(alpha1);

          double q1 = -xcenter;
          double q2 = -ycenter;

          double q4 = r0*sin(alpha1)-xcenter;
          double q5 = r0*cos(alpha1)-ycenter;

          double q7 = r1*sin(alpha1)-xcenter;
          double q8 = -r1*cos(alpha1)-ycenter;

          clist_q0[i][0] = q1;
          clist_q0[i][1] = q2;
          clist_q0[i][2] = 0;

          clist_q0[i][3] = q4;
          clist_q0[i][4] = q5;
          clist_q0[i][5] = 0;

          clist_q0[i][6] = q7;
          clist_q0[i][7] = q8;
          clist_q0[i][8] = 0;

          clist_nselect1[i] = 2;
          clist_select1[i][0] = 1; clist_select1[i][1] = 2;
          clist_nselect2[i] = 1;
          clist_select2[i][0] = 1;
        }  else if (N == 4)
        {
          oxy = atom->map(shake_atom[m][0]);
          atom1 = atom->map(shake_atom[m][1]);
          atom2 = atom->map(shake_atom[m][2]);
          atom3 = atom->map(shake_atom[m][3]);

          r0 = bond_distance[shake_type[m][0]];
          r1 = bond_distance[shake_type[m][1]];
          r2 = bond_distance[shake_type[m][2]];

          m1 = atom->mass[atom->type[oxy]];
          m2 = atom->mass[atom->type[atom1]];
          m3 = atom->mass[atom->type[atom2]];
          m4 = atom->mass[atom->type[atom3]];
          m_all = m1 + m2 + m3 + m4;

          //how to get these angles?
          double alpha1 = MY_PI/2.57; //roughly 70 degrees
          double alpha2 = MY_PI/3;
          //ENSURE ycenter,zcenter = 0!
          //approximate xcenter, if r0 !=r1 != r2, exact if r0 = r1 = r2
          double xcenter =  (m2*r0+m3*r1+m4*r2)/m_all*cos(alpha1);

          clist_q0[i][0] = -xcenter;
          clist_q0[i][1] = 0.0;
          clist_q0[i][2] = 0.0;
          clist_q0[i][3] = r0*cos(alpha1)-xcenter;
          clist_q0[i][4] = -r0*sin(alpha1);
          clist_q0[i][5] = 0.0;
          clist_q0[i][6] = r1*cos(alpha1)-xcenter;
          clist_q0[i][7] = r1*sin(alpha1)*cos(alpha2);
          clist_q0[i][8] = -r1*sin(alpha1)*sin(alpha2);
          clist_q0[i][9]  = r2*cos(alpha1)-xcenter;
          clist_q0[i][10] = r2*sin(alpha1)*cos(alpha2);
          clist_q0[i][11] = r2*sin(alpha1)*sin(alpha2);

          clist_nselect1[i] = 3;
          clist_nselect2[i] = 2;

          clist_select1[i][0] = 1;clist_select1[i][1] = 2;clist_select1[i][2] = 3;
          clist_select2[i][0] = 2;clist_select2[i][1] = 3;

          //signum ensures correct ordering of three satellites
          //signum = sign(cross(x2-x1,x3-x1))T*(x1-x0))
          double del1[3], del2[3], del3[3];
          del1[0] = x[atom1][0]-x[oxy][0];
          del1[1] = x[atom1][1]-x[oxy][1];
          del1[2] = x[atom1][2]-x[oxy][2];
          domain->minimum_image(del1);

          del2[0] = x[atom2][0]-x[atom1][0];
          del2[1] = x[atom2][1]-x[atom1][1];
          del2[2] = x[atom2][2]-x[atom1][2];
          domain->minimum_image(del2);

          del3[0] = x[atom3][0]-x[atom1][0];
          del3[1] = x[atom3][1]-x[atom1][1];
          del3[2] = x[atom3][2]-x[atom1][2];
          domain->minimum_image(del3);

          double a = (del2[1])*(del3[2]) - (del2[2])*(del3[1]);
          double b = (del2[2])*(del3[0]) - (del2[0])*(del3[2]);
          double c = (del2[0])*(del3[1]) - (del2[1])*(del3[0]);
          int signum = sgn(a*(del1[0]) + b*(del1[1]) + c*(del1[2]));

          if (abs(signum) != 1)
            error->all(FLERR,"Wrong orientation in cluster of size 4"
                       "in fix filter/corotate!");
          clist_q0[i][8] *= signum;
          clist_q0[i][11] *= signum;

        } else if (N == 5) {
        oxy = atom->map(shake_atom[m][0]);
        atom1 = atom->map(shake_atom[m][1]);
        atom2 = atom->map(shake_atom[m][2]);
        atom3 = atom->map(shake_atom[m][3]);
        int c1 = atom->map(shake_atom[m][4]);

        r1 = bond_distance[shake_type[m][3]];
        r0 = bond_distance[shake_type[m][0]];
        theta0 = angle_distance[shake_type[m][2]];

        m1 = atom->mass[atom->type[oxy]];
        m2 = atom->mass[atom->type[atom1]];
        m3 = atom->mass[atom->type[atom2]];
        m4 = atom->mass[atom->type[atom3]];
        m5 = atom->mass[atom->type[c1]];
        m_all = m1 + m2 + m3 + m4 + m5;

        double alpha1 = MY_PI/2.57; //roughly 70 degrees
        double alpha2 = MY_PI/3;
        //ENSURE ycenter,zcenter = 0!
        double xcenter =  -(m2+m3+m4)/m_all*r0*cos(alpha1) +r1*m5/m_all;

        clist_q0[i][0] = -xcenter;
        clist_q0[i][1] = 0.0;
        clist_q0[i][2] = 0.0;
        clist_q0[i][3] = -r0*cos(alpha1)-xcenter;
        clist_q0[i][4] = -r0*sin(alpha1);
        clist_q0[i][5] = 0.0;
        clist_q0[i][6] = -r0*cos(alpha1)-xcenter;
        clist_q0[i][7] = r0*sin(alpha1)*cos(alpha2);
        clist_q0[i][8] = r0*sin(alpha1)*sin(alpha2);
        clist_q0[i][9]  = -r0*cos(alpha1)-xcenter;
        clist_q0[i][10] = r0*sin(alpha1)*cos(alpha2);
        clist_q0[i][11] = -r0*sin(alpha1)*sin(alpha2);
        clist_q0[i][12] = r1-xcenter;
        clist_q0[i][13] = 0.0;
        clist_q0[i][14] = 0.0;
        clist_nselect1[i] = 1;
        clist_nselect2[i] = 2;

        clist_select1[i][0] = 4;
        clist_select2[i][0] = 2;clist_select2[i][1] = 3;

        //signum ensures correct ordering of three satellites
        //signum = sign(cross(x2-x1,x3-x1))T*(x1-x0))
        double del1[3], del2[3], del3[3];
        del1[0] = x[atom1][0]-x[oxy][0];
        del1[1] = x[atom1][1]-x[oxy][1];
        del1[2] = x[atom1][2]-x[oxy][2];
        domain->minimum_image(del1);

        del2[0] = x[atom2][0]-x[atom1][0];
        del2[1] = x[atom2][1]-x[atom1][1];
        del2[2] = x[atom2][2]-x[atom1][2];
        domain->minimum_image(del2);

        del3[0] = x[atom3][0]-x[atom1][0];
        del3[1] = x[atom3][1]-x[atom1][1];
        del3[2] = x[atom3][2]-x[atom1][2];
        domain->minimum_image(del3);

        double a = (del2[1])*(del3[2]) - (del2[2])*(del3[1]);
        double b = (del2[2])*(del3[0]) - (del2[0])*(del3[2]);
        double c = (del2[0])*(del3[1]) - (del2[1])*(del3[0]);
        int signum = sgn(a*(del1[0]) + b*(del1[1]) + c*(del1[2]));

        if (abs(signum)!= 1)
          error->all(FLERR,"Wrong orientation in cluster of size 5"
                     "in fix filter/corotate!");
        clist_q0[i][8] *= signum;
        clist_q0[i][11] *= signum;
      } else {
        error->all(FLERR,"Fix filter/corotate cluster with size > 5"
                   "not yet configured...");
      }
    }
}

/* ----------------------------------------------------------------------
 *  setup, needs to patch missing setup_post_force_respa() in respa...
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::setup(int vflag)
{
  (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
  post_force_respa(vflag,nlevels_respa-1,0);
  (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
}

void FixFilterCorotate::setup_pre_force_respa(int vflag,int ilevel) {
  pre_force_respa(vflag,ilevel,0);
}

//void FixFilterCorotate::setup_post_force_respa(int vflag,int ilevel) {
//  post_force_respa(vflag,ilevel,0);
//}

double FixFilterCorotate::compute_array(int,int)
{
  filter_inner();
  return 1;
}

void FixFilterCorotate::pre_force_respa(int /*vflag*/, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1)
  {
    filter_inner();

    //%Switch pointers x<-->x_filter
    x_store = atom->x;
    atom->x = array_atom;
  }
}

void FixFilterCorotate::post_force_respa(int /*vflag*/, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1)
  {
    atom->x = x_store;

    comm->forward_comm(this); //comm forward of forces for outer filter

    filter_outer();
  }
}

/* ---------------------------------------------------------------------- */

void FixFilterCorotate::filter_inner()
{
  int nall = atom->nlocal + atom->nghost;
  double **x = atom->x;

  //set array_atom to atom->x, in case some atoms are not part of any cluster:
  for (int i=0; i<nall; i++)
  {
    array_atom[i][0] = x[i][0];
    array_atom[i][1] = x[i][1];
    array_atom[i][2] = x[i][2];
  }

  int m;
  //iterate through list:
  for (int i=0; i<nlist; i++)
  {
    m = list[i];

    //filter:
    general_cluster(m,i);
  }
}

/* ---------------------------------------------------------------------- */

void FixFilterCorotate::filter_outer()
{
  double sum1,sum2,sum3;
  double** f = atom->f;
  int m;

  for (int i=0; i<nlist; i++)
  {
    m = list[i];
    int dim = shake_flag[m];
    if (dim == 1) //case: angle
      dim = 3;

    for (int j = 0; j < dim; j++)
    {
      sum1 = sum2 = sum3 = 0;
      for (int k = 0; k < dim; k++)
      {
        //get local tag of k-th atom in cluster i:
        int kk = atom->map(shake_atom[m][k]);
        //apply derivative times force
        sum1 += clist_derv[i][3*j][3*k]*f[kk][0]+clist_derv[i][3*j][3*k+1]*
          f[kk][1]+clist_derv[i][3*j][3*k+2]*f[kk][2];
        sum2 += clist_derv[i][3*j+1][3*k]*f[kk][0]+clist_derv[i][3*j+1][3*k+1]*
          f[kk][1]+clist_derv[i][3*j+1][3*k+2]*f[kk][2];
        sum3 += clist_derv[i][3*j+2][3*k]*f[kk][0]+clist_derv[i][3*j+2][3*k+1]*
          f[kk][1]+clist_derv[i][3*j+2][3*k+2]*f[kk][2];
      }

      g[3*j] = sum1;
      g[3*j+1] = sum2;
      g[3*j+2] = sum3;
    }

    //replace force f with filtered value:
    for (int j = 0; j < dim; j++)
    {
      int kk = atom->map(shake_atom[m][j]);

      f[kk][0] = g[3*j];
      f[kk][1] = g[3*j+1];
      f[kk][2] = g[3*j+2];
    }
  }
}

/* ----------------------------------------------------------------------
 *  identify whether each atom is in a SHAKE cluster
 *  only include atoms in fix group and those bonds/angles specified in input
 *  test whether all clusters are valid
 *  set shake_flag, shake_atom, shake_type values
 *  set bond,angle types negative so will be ignored in neighbor lists
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::find_clusters()
{
  int i,j,m,n;
  int flag,flag_all,nbuf,size;
  double massone;
  tagint *buf;

  if (me == 0 && screen)
    fprintf(screen,"Finding filter clusters ...\n");

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  int nlocal = atom->nlocal;
  int angles_allow = atom->avec->angles_allow;

  // setup ring of procs

  int next = me + 1;
  int prev = me -1;
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  // -----------------------------------------------------
  // allocate arrays for self (1d) and bond partners (2d)
  // max = max # of bond partners for owned atoms = 2nd dim of partner arrays
  // npartner[i] = # of bonds attached to atom i
  // nshake[i] = # of SHAKE bonds attached to atom i
  // partner_tag[i][] = global IDs of each partner
  // partner_mask[i][] = mask of each partner
  // partner_type[i][] = type of each partner
  // partner_massflag[i][] = 1 if partner meets mass criterion, 0 if not
  // partner_bondtype[i][] = type of bond attached to each partner
  // partner_shake[i][] = 1 if SHAKE bonded to partner, 0 if not
  // partner_nshake[i][] = nshake value for each partner
  // -----------------------------------------------------

  int max = 0;

  for (i = 0; i < nlocal; i++) max = MAX(max,nspecial[i][0]);

  int *npartner;
  memory->create(npartner,nlocal,"shake:npartner");
  memory->create(nshake,nlocal,"shake:nshake");

  tagint **partner_tag;
  int **partner_mask,**partner_type,**partner_massflag;
  int **partner_bondtype,**partner_shake,**partner_nshake;
  memory->create(partner_tag,nlocal,max,"shake:partner_tag");
  memory->create(partner_mask,nlocal,max,"shake:partner_mask");
  memory->create(partner_type,nlocal,max,"shake:partner_type");
  memory->create(partner_massflag,nlocal,max,"shake:partner_massflag");
  memory->create(partner_bondtype,nlocal,max,"shake:partner_bondtype");
  memory->create(partner_shake,nlocal,max,"shake:partner_shake");
  memory->create(partner_nshake,nlocal,max,"shake:partner_nshake");

  // -----------------------------------------------------
  // set npartner and partner_tag from special arrays
  // -----------------------------------------------------

  for (i = 0; i < nlocal; i++) {
    npartner[i] = nspecial[i][0];
    for (j = 0; j < npartner[i]; j++)
      partner_tag[i][j] = special[i][j];
  }

  // -----------------------------------------------------
  // set partner_mask, partner_type, partner_massflag, partner_bondtype
  //   for bonded partners
  // requires communication for off-proc partners
  // -----------------------------------------------------

  // fill in mask, type, massflag, bondtype if own bond partner
  // info to store in buf for each off-proc bond = nper = 6
  //   2 atoms IDs in bond, space for mask, type, massflag, bondtype
  // nbufmax = largest buffer needed to hold info from any proc

  int nper = 6;

  nbuf = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      partner_mask[i][j] = 0;
      partner_type[i][j] = 0;
      partner_massflag[i][j] = 0;
      partner_bondtype[i][j] = 0;

      m = atom->map(partner_tag[i][j]);
      if (m >= 0 && m < nlocal) {
        partner_mask[i][j] = mask[m];
        partner_type[i][j] = type[m];
        if (nmass) {
          if (rmass) massone = rmass[m];
          else massone = mass[type[m]];
          partner_massflag[i][j] = masscheck(massone);
        }
        n = bondtype_findset(i,tag[i],partner_tag[i][j],0);
        if (n) partner_bondtype[i][j] = n;
        else {
          n = bondtype_findset(m,tag[i],partner_tag[i][j],0);
          if (n) partner_bondtype[i][j] = n;
        }
      } else nbuf += nper;
    }
  }

  memory->create(buf,nbuf,"filter/corotate:buf");

  // fill buffer with info

  size = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      m = atom->map(partner_tag[i][j]);
      if (m < 0 || m >= nlocal) {
        buf[size] = tag[i];
        buf[size+1] = partner_tag[i][j];
        buf[size+2] = 0;
        buf[size+3] = 0;
        buf[size+4] = 0;
        n = bondtype_findset(i,tag[i],partner_tag[i][j],0);
        if (n) buf[size+5] = n;
        else buf[size+5] = 0;
        size += nper;
      }
    }
  }

  // cycle buffer around ring of procs back to self

  comm->ring(size,sizeof(tagint),buf,1,ring_bonds,buf,(void *)this);

  // store partner info returned to me

  m = 0;
  while (m < size) {
    i = atom->map(buf[m]);
    for (j = 0; j < npartner[i]; j++) {
      if (buf[m+1] == partner_tag[i][j]) break;
    }
    partner_mask[i][j] = buf[m+2];
    partner_type[i][j] = buf[m+3];
    partner_massflag[i][j] = buf[m+4];
    partner_bondtype[i][j] = buf[m+5];
    m += nper;
  }

  memory->destroy(buf);

  // error check for unfilled partner info
  // if partner_type not set, is an error
  // partner_bondtype may not be set if special list is not consistent
  //   with bondatom (e.g. due to delete_bonds command)
  // this is OK if one or both atoms are not in fix group, since
  //   bond won't be SHAKEn anyway
  // else it's an error

  flag = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      if (partner_type[i][j] == 0) flag = 1;
      if (!(mask[i] & groupbit)) continue;
      if (!(partner_mask[i][j] & groupbit)) continue;
      if (partner_bondtype[i][j] == 0) flag = 1;
    }
  }

  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all)
    error->all(FLERR,"Did not find fix filter/corotate partner info");

  // -----------------------------------------------------
  // identify SHAKEable bonds
  // set nshake[i] = # of SHAKE bonds attached to atom i
  // set partner_shake[i][] = 1 if SHAKE bonded to partner, 0 if not
  // both atoms must be in group, bondtype must be > 0
  // check if bondtype is in input bond_flag
  // check if type of either atom is in input type_flag
  // check if mass of either atom is in input mass_list
  // -----------------------------------------------------

  int np;

  for (i = 0; i < nlocal; i++) {
    nshake[i] = 0;
    np = npartner[i];
    for (j = 0; j < np; j++) {
      partner_shake[i][j] = 0;

      if (!(mask[i] & groupbit)) continue;
      if (!(partner_mask[i][j] & groupbit)) continue;
      if (partner_bondtype[i][j] <= 0) continue;

      if (bond_flag[partner_bondtype[i][j]]) {
        partner_shake[i][j] = 1;
        nshake[i]++;
        continue;
      }
      if (type_flag[type[i]] || type_flag[partner_type[i][j]]) {
        partner_shake[i][j] = 1;
        nshake[i]++;
        continue;
      }
      if (nmass) {
        if (partner_massflag[i][j]) {
          partner_shake[i][j] = 1;
          nshake[i]++;
          continue;
        } else {
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          if (masscheck(massone)) {
            partner_shake[i][j] = 1;
            nshake[i]++;
            continue;
          }
        }
      }
    }
  }

  // -----------------------------------------------------
  // set partner_nshake for bonded partners
  // requires communication for off-proc partners
  // -----------------------------------------------------

  // fill in partner_nshake if own bond partner
  // info to store in buf for each off-proc bond =
  //   2 atoms IDs in bond, space for nshake value
  // nbufmax = largest buffer needed to hold info from any proc

  nbuf = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      m = atom->map(partner_tag[i][j]);
      if (m >= 0 && m < nlocal) partner_nshake[i][j] = nshake[m];
      else nbuf += 3;
    }
  }

  memory->create(buf,nbuf,"filter/corotate:buf");

  // fill buffer with info

  size = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < npartner[i]; j++) {
      m = atom->map(partner_tag[i][j]);
      if (m < 0 || m >= nlocal) {
        buf[size] = tag[i];
        buf[size+1] = partner_tag[i][j];
        size += 3;
      }
    }
  }

  // cycle buffer around ring of procs back to self

  comm->ring(size,sizeof(tagint),buf,2,ring_nshake,buf,(void *)this);

  // store partner info returned to me

  m = 0;
  while (m < size) {
    i = atom->map(buf[m]);
    for (j = 0; j < npartner[i]; j++) {
      if (buf[m+1] == partner_tag[i][j]) break;
    }
    partner_nshake[i][j] = buf[m+2];
    m += 3;
  }

  memory->destroy(buf);

  // -----------------------------------------------------
  // error checks
  // no atom with nshake > 3
  // no connected atoms which both have nshake > 1
  // -----------------------------------------------------

  flag = 0;
  for (i = 0; i < nlocal; i++) if (nshake[i] > 4) flag = 1;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Cluster of more than 5 atoms");

  flag = 0;
  for (i = 0; i < nlocal; i++) {
    if (nshake[i] <= 1) continue;
    for (j = 0; j < npartner[i]; j++)
      if (partner_shake[i][j] && partner_nshake[i][j] > 1) flag = 1;
  }
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Clusters are connected");

  // -----------------------------------------------------
  // set SHAKE arrays that are stored with atoms & add angle constraints
  // zero shake arrays for all owned atoms
  // if I am central atom set shake_flag & shake_atom & shake_type
  // for 2-atom clusters, I am central atom if my atom ID < partner ID
  // for 3-atom clusters, test for angle constraint
  //   angle will be stored by this atom if it exists
  //   if angle type matches angle_flag, then it is angle-constrained
  // shake_flag[] = 0 if atom not in SHAKE cluster
  //                2,3,4 = size of bond-only cluster
  //                1 = 3-atom angle cluster
  // shake_atom[][] = global IDs of 2,3,4 atoms in cluster
  //                  central atom is 1st
  //                  for 2-atom cluster, lowest ID is 1st
  // shake_type[][] = bondtype of each bond in cluster
  //                  for 3-atom angle cluster, 3rd value is angletype
  // -----------------------------------------------------

  for (i = 0; i < nlocal; i++) {
    shake_flag[i] = 0;
    shake_atom[i][0] = 0;
    shake_atom[i][1] = 0;
    shake_atom[i][2] = 0;
    shake_atom[i][3] = 0;
    shake_atom[i][4] = 0;
    shake_type[i][0] = 0;
    shake_type[i][1] = 0;
    shake_type[i][2] = 0;
    shake_type[i][3] = 0;

    if (nshake[i] == 1) {
      for (j = 0; j < npartner[i]; j++) {
        if (partner_shake[i][j]) break;
      }
      if (partner_nshake[i][j] == 1 && tag[i] < partner_tag[i][j]) {
        shake_flag[i] = 2;
        shake_atom[i][0] = tag[i];
        shake_atom[i][1] = partner_tag[i][j];
        shake_type[i][0] = partner_bondtype[i][j];
      }
    }

    if (nshake[i] > 1) {
      shake_flag[i] = 1;
      shake_atom[i][0] = tag[i];
      for (j = 0; j < npartner[i]; j++)
        if (partner_shake[i][j]) {
          m = shake_flag[i];
          shake_atom[i][m] = partner_tag[i][j];
          shake_type[i][m-1] = partner_bondtype[i][j];
          shake_flag[i]++;
        }
    }

    if (nshake[i] == 2 && angles_allow) {
      n = angletype_findset(i,shake_atom[i][1],shake_atom[i][2],0);
      if (n <= 0) continue;
      //if (angle_flag[n]) {    Take all angles!
      shake_flag[i] = 1;
      shake_type[i][2] = n;
      //}
    }
  }

  // -----------------------------------------------------
  // set shake_flag,shake_atom,shake_type for non-central atoms
  // requires communication for off-proc atoms
  // -----------------------------------------------------

  // fill in shake arrays for each bond partner I own
  // info to store in buf for each off-proc bond =
  //   all values from shake_flag, shake_atom, shake_type
  // nbufmax = largest buffer needed to hold info from any proc

  nbuf = 0;
  for (i = 0; i < nlocal; i++) {
    if (shake_flag[i] == 0) continue;
    for (j = 0; j < npartner[i]; j++) {
      if (partner_shake[i][j] == 0) continue;
      m = atom->map(partner_tag[i][j]);
      if (m >= 0 && m < nlocal) {
        shake_flag[m] = shake_flag[i];
        shake_atom[m][0] = shake_atom[i][0];
        shake_atom[m][1] = shake_atom[i][1];
        shake_atom[m][2] = shake_atom[i][2];
        shake_atom[m][3] = shake_atom[i][3];
        shake_atom[m][4] = shake_atom[i][4];
        shake_type[m][0] = shake_type[i][0];
        shake_type[m][1] = shake_type[i][1];
        shake_type[m][2] = shake_type[i][2];
        shake_type[m][3] = shake_type[i][3];
      } else nbuf += 11;
    }
  }

  memory->create(buf,nbuf,"filter/corotate:buf");

  // fill buffer with info

  size = 0;
  for (i = 0; i < nlocal; i++) {
    if (shake_flag[i] == 0) continue;
    for (j = 0; j < npartner[i]; j++) {
      if (partner_shake[i][j] == 0) continue;
      m = atom->map(partner_tag[i][j]);
      if (m < 0 || m >= nlocal) {
        buf[size] = partner_tag[i][j];
        buf[size+1] = shake_flag[i];
        buf[size+2] = shake_atom[i][0];
        buf[size+3] = shake_atom[i][1];
        buf[size+4] = shake_atom[i][2];
        buf[size+5] = shake_atom[i][3];
        buf[size+6] = shake_atom[i][4];
        buf[size+7] = shake_type[i][0];
        buf[size+8] = shake_type[i][1];
        buf[size+9] = shake_type[i][2];
        buf[size+10] = shake_type[i][3];
        size += 11;
      }
    }
  }

  // cycle buffer around ring of procs back to self

  comm->ring(size,sizeof(tagint),buf,3,ring_shake,nullptr,(void *)this);

  memory->destroy(buf);

  // -----------------------------------------------------
  // free local memory
  // -----------------------------------------------------

  memory->destroy(npartner);
  memory->destroy(nshake);
  memory->destroy(partner_tag);
  memory->destroy(partner_mask);
  memory->destroy(partner_type);
  memory->destroy(partner_massflag);
  memory->destroy(partner_bondtype);
  memory->destroy(partner_shake);
  memory->destroy(partner_nshake);

  // -----------------------------------------------------
  // print info on SHAKE clusters
  // -----------------------------------------------------

  int count1,count2,count3,count4,count5;
  count1 = count2 = count3 = count4 = count5 = 0;
  for (i = 0; i < nlocal; i++) {
    if (shake_flag[i] == 1) count1++;
    else if (shake_flag[i] == 2) count2++;
    else if (shake_flag[i] == 3) count3++;
    else if (shake_flag[i] == 4) count4++;
    else if (shake_flag[i] == 5) count5++;
  }

  int tmp;
  tmp = count1;
  MPI_Allreduce(&tmp,&count1,1,MPI_INT,MPI_SUM,world);
  tmp = count2;
  MPI_Allreduce(&tmp,&count2,1,MPI_INT,MPI_SUM,world);
  tmp = count3;
  MPI_Allreduce(&tmp,&count3,1,MPI_INT,MPI_SUM,world);
  tmp = count4;
  MPI_Allreduce(&tmp,&count4,1,MPI_INT,MPI_SUM,world);
  tmp = count5;
  MPI_Allreduce(&tmp,&count5,1,MPI_INT,MPI_SUM,world);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  %d = # of size 2 clusters\n",count2/2);
      fprintf(screen,"  %d = # of size 3 clusters\n",count3/3);
      fprintf(screen,"  %d = # of size 4 clusters\n",count4/4);
      fprintf(screen,"  %d = # of size 5 clusters\n",count5/5);
      fprintf(screen,"  %d = # of frozen angles\n",count1/3);
    }
    if (logfile) {
      fprintf(logfile,"  %d = # of size 2 clusters\n",count2/2);
      fprintf(logfile,"  %d = # of size 3 clusters\n",count3/3);
      fprintf(logfile,"  %d = # of size 4 clusters\n",count4/4);
      fprintf(logfile,"  %d = # of size 5 clusters\n",count5/5);
      fprintf(logfile,"  %d = # of frozen angles\n",count1/3);
    }
  }
}

/* ----------------------------------------------------------------------
 *  when receive buffer, scan bond partner IDs for atoms I own
 *  if I own partner:
 *    fill in mask and type and massflag
 *    search for bond with 1st atom and fill in bondtype
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::ring_bonds(int ndatum, char *cbuf, void *ptr)
{
  auto ffptr = (FixFilterCorotate *) ptr;
  Atom *atom = ffptr->atom;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nmass = ffptr->nmass;

  auto buf = (tagint *) cbuf;
  int m,n;
  double massone;

  for (int i = 0; i < ndatum; i += 6) {
    m = atom->map(buf[i+1]);
    if (m >= 0 && m < nlocal) {
      buf[i+2] = mask[m];
      buf[i+3] = type[m];
      if (nmass) {
        if (rmass) massone = rmass[m];
        else massone = mass[type[m]];
        buf[i+4] = ffptr->masscheck(massone);
      }
      if (buf[i+5] == 0) {
        n = ffptr->bondtype_findset(m,buf[i],buf[i+1],0);
        if (n) buf[i+5] = n;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 *  when receive buffer, scan bond partner IDs for atoms I own
 *  if I own partner, fill in nshake value
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::ring_nshake(int ndatum, char *cbuf, void *ptr)
{
  auto ffptr = (FixFilterCorotate *) ptr;
  Atom *atom = ffptr->atom;
  int nlocal = atom->nlocal;

  int *nshake = ffptr->nshake;

  auto buf = (tagint *) cbuf;
  int m;

  for (int i = 0; i < ndatum; i += 3) {
    m = atom->map(buf[i+1]);
    if (m >= 0 && m < nlocal) buf[i+2] = nshake[m];
  }
}

/* ----------------------------------------------------------------------
 *  when receive buffer, scan bond partner IDs for atoms I own
 *  if I own partner, fill in nshake value
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::ring_shake(int ndatum, char *cbuf, void *ptr)
{
  auto ffptr = (FixFilterCorotate *) ptr;
  Atom *atom = ffptr->atom;
  int nlocal = atom->nlocal;

  int *shake_flag = ffptr->shake_flag;
  tagint **shake_atom = ffptr->shake_atom;
  int **shake_type = ffptr->shake_type;

  auto buf = (tagint *) cbuf;
  int m;

  for (int i = 0; i < ndatum; i += 11) {
    m = atom->map(buf[i]);
    if (m >= 0 && m < nlocal) {
      shake_flag[m] = buf[i+1];
      shake_atom[m][0] = buf[i+2];
      shake_atom[m][1] = buf[i+3];
      shake_atom[m][2] = buf[i+4];
      shake_atom[m][3] = buf[i+5];
      shake_atom[m][4] = buf[i+6];
      shake_type[m][0] = buf[i+7];
      shake_type[m][1] = buf[i+8];
      shake_type[m][2] = buf[i+9];
      shake_type[m][3] = buf[i+10];
    }
  }
}

/* ----------------------------------------------------------------------
 *  check if massone is within MASSDELTA of any mass in mass_list
 *  return 1 if yes, 0 if not
 * ------------------------------------------------------------------------- */

int FixFilterCorotate::masscheck(double massone)
{
  for (int i = 0; i < nmass; i++) {
    if (fabs(mass_list[i]-massone) <= MASSDELTA) return 1;
  }
  return 0;
}

/* --------------------------------------------------------------------------
 * filter one cluster of arbitrary size
 * ---------------------------------------------------------------------------*/

void FixFilterCorotate::general_cluster(int index, int index_in_list)
{
  //index: shake index, index_in_list: corresponding index in list
  //get q0, nselect:

  double *q0  = clist_q0[index_in_list];
  int nselect1  = clist_nselect1[index_in_list];
  int nselect2  = clist_nselect2[index_in_list];
  int *select1  = clist_select1[index_in_list];
  int *select2  = clist_select2[index_in_list];

  int N = shake_flag[index];      //number of atoms in cluster
  if (N == 1) N = 3;              //angle cluster

  double**x = atom->x;
  double norm1, norm2, norm3;

  int* list_cluster = new int[N]; // contains local IDs of cluster atoms,
                                  // 0 = center
  auto  m = new double[N];      //contains local mass
  auto r = new double[N];      //contains r[i] = 1/||del[i]||
  auto  del = new double*[N];  //contains del[i] = x_i-x_0
  for (int i = 0; i<N; i++)
    del[i] = new double[3];

  for (int i = 0; i < N; i++) {
    list_cluster[i] = atom->map(shake_atom[index][i]);
    m[i] = atom->mass[atom->type[list_cluster[i]]];
  }

  //%CALC r_i:
  for (int i = 1; i < N; i++) {
    del[i][0] = x[list_cluster[i]][0] - x[list_cluster[0]][0];
    del[i][1] = x[list_cluster[i]][1] - x[list_cluster[0]][1];
    del[i][2] = x[list_cluster[i]][2] - x[list_cluster[0]][2];
    domain->minimum_image(del[i]);
    r[i] = 1.0/sqrt(del[i][0]*del[i][0]+del[i][1]*del[i][1]+
      del[i][2]*del[i][2]);
  }
  //special case i=0: need del[0] for compact notation of array_atom later...
  del[0][0] = del[0][1] = del[0][2] = r[0] = 0.0;

  //calc n1:
  n1[0] = n1[1] = n1[2] = 0.0;

  int k;
  for (int i = 0; i < nselect1; i++) {
    k = select1[i];
    n1[0] += del[k][0]*r[k];
    n1[1] += del[k][1]*r[k];
    n1[2] += del[k][2]*r[k];
  }
  norm1 = 1.0/sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);

  n1[0] *= norm1;
  n1[1] *= norm1;
  n1[2] *= norm1;

  //calc n2:
  n2[0] = n2[1] = n2[2] = 0.0;

  for (int i = 0; i < nselect2; i++) {
    k = select2[i];
    n2[0] += del[k][0]*r[k];
    n2[1] += del[k][1]*r[k];
    n2[2] += del[k][2]*r[k];
  }
  if (!nselect2) //cluster of size 2, has only one direction
    n2[0] = n2[1] = n2[2] = 1.0;

  double alpha = n2[0]*n1[0] + n2[1]*n1[1] + n2[2]*n1[2];

  n2[0] -= alpha*n1[0];
  n2[1] -= alpha*n1[1];
  n2[2] -= alpha*n1[2];

  norm2 = 1.0/sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);

  n2[0] *= norm2;
  n2[1] *= norm2;
  n2[2] *= norm2;

  n3[0] = n1[1]*n2[2] - n1[2]*n2[1];
  n3[1] = n1[2]*n2[0] - n1[0]*n2[2];
  n3[2] = n1[0]*n2[1] - n1[1]*n2[0];

  norm3 = 1.0/sqrt(n3[0]*n3[0]+n3[1]*n3[1]+n3[2]*n3[2]);

  n3[0] *= norm3;
  n3[1] *= norm3;
  n3[2] *= norm3;

  //%x_filter:
  for (int i = 0; i < N; i++) {
    k = list_cluster[i];
    array_atom[k][0] = x[k][0]+q0[3*i]*n1[0]+q0[3*i+1]*n2[0]+q0[3*i+2]*n3[0];
    array_atom[k][1] = x[k][1]+q0[3*i]*n1[1]+q0[3*i+1]*n2[1]+q0[3*i+2]*n3[1];
    array_atom[k][2] = x[k][2]+q0[3*i]*n1[2]+q0[3*i+1]*n2[2]+q0[3*i+2]*n3[2];
  }

  double m_all = 0.0;
  for (int i = 0; i<N; i++)
    m_all += m[i];

  //%x+X_center

  for (int i = 0; i<N; i++) {
    k = list_cluster[i];
    for (int j = 0; j < N; j++) {
      array_atom[k][0] += m[j]/m_all*(del[j][0]-del[i][0]);
      array_atom[k][1] += m[j]/m_all*(del[j][1]-del[i][1]);
      array_atom[k][2] += m[j]/m_all*(del[j][2]-del[i][2]);
    }
  }

  //derivative:
  //dn1dx:

  double **sum1;
  memory->create(sum1,3,3*N,"filter_corotate:sum1");

  for (int i=0; i<3; i++)
    for (int j=0; j<3*N; j++)
      sum1[i][j] = 0;

  double I3mn1n1T[3][3];   //(I_3 - n1n1T)
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++)
      I3mn1n1T[i][j] = -n1[i]*n1[j];
    I3mn1n1T[i][i] += 1.0;
  }

  // sum1 part of dn1dx:

  for (int l = 0; l < nselect1; l++) {
    k = select1[l];
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        double help = del[k][i]*del[k][j]*r[k]*r[k]*r[k];
        sum1[i][j+3*k] -= help;
        sum1[i][j] += help;
      }
      sum1[i][i+3*k] += r[k];
      sum1[i][i] -= r[k];
    }
  }

  //dn1dx = norm1 * I3mn1n1T * sum1

  for (int i=0; i<3; i++) {
    for (int j=0; j<3*N; j++) {
      double sum = 0;
      for (int l = 0; l<3; l++)
        sum += I3mn1n1T[i][l]*sum1[l][j];
      dn1dx[i][j] = norm1*sum;
    }
  }
  memory->destroy(sum1);

  //dn2dx: norm2 * I3mn2n2T * (I3mn1n1T*sum2 - rkn1pn1rk*dn1dx)

  double **sum2;
  memory->create(sum2,3,3*N,"filter_corotate:sum2");
  for (int i=0; i<3; i++)
    for (int j=0; j<3*N; j++)
      sum2[i][j] = 0;

  double I3mn2n2T[3][3];   //(I_3 - n2n2T)
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++)
      I3mn2n2T[i][j] = -n2[i]*n2[j];
    I3mn2n2T[i][i] += 1.0;
  }

  // sum2 part of dn1dx:

  for (int l = 0; l < nselect2; l++) {
    k = select2[l];
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        double help = del[k][i]*del[k][j]*r[k]*r[k]*r[k];
        sum2[i][j+3*k] -= help;
        sum2[i][j] += help;
      }
      sum2[i][i+3*k] += r[k];
      sum2[i][i] -= r[k];
    }
  }

  //prefactor:

  double rkn1pn1rk[3][3];
  double rk[3]; rk[0] = rk[1] = rk[2] = 0.0;

  for (int i = 0; i < nselect2; i++) {
    k = select2[i];
    rk[0] += del[k][0]*r[k];
    rk[1] += del[k][1]*r[k];
    rk[2] += del[k][2]*r[k];
  }

  //rkn1pn1rk = rkT*n1*I3 + n1*rkT

  double scalar = rk[0]*n1[0]+rk[1]*n1[1]+rk[2]*n1[2];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++)
      rkn1pn1rk[i][j] = n1[i]*rk[j];
    rkn1pn1rk[i][i] += scalar;
  }

  //dn2dx: norm2 * I3mn2n2T * (I3mn1n1T*sum2 - rkn1pn1rk*dn1dx)
  //sum3 = (I3mn1n1T*sum2 - rkn1pn1rk*dn1dx)

  double **sum3;
  memory->create(sum3,3,3*N,"filter_corotate:sum3");
  for (int i=0; i<3; i++)
    for (int j=0; j<3*N; j++) {
      double sum = 0;
      for (int l = 0; l<3; l++)
        sum += I3mn1n1T[i][l]*sum2[l][j] - rkn1pn1rk[i][l]*dn1dx[l][j];
      sum3[i][j] = sum;
    }

  memory->destroy(sum2);
  //dn2dx = norm2 * I3mn2n2T * sum3
  for (int i=0; i<3; i++)
    for (int j=0; j<3*N; j++) {
      double sum = 0;
      for (int l = 0; l<3; l++)
        sum += I3mn2n2T[i][l]*sum3[l][j];
      dn2dx[i][j] = norm2*sum;
    }

  memory->destroy(sum3);
  //dn3dx = norm3 * I3mn3n3T * cross
  double I3mn3n3T[3][3];   //(I_3 - n3n3T)
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++)
      I3mn3n3T[i][j] = -n3[i]*n3[j];
    I3mn3n3T[i][i] += 1.0;
  }

  double **cross;
  memory->create(cross,3,3*N,"filter_corotate:cross");

  for (int j=0; j<3*N; j++) {
    cross[0][j] = dn1dx[1][j]*n2[2] -dn1dx[2][j]*n2[1] +
      n1[1]*dn2dx[2][j]-n1[2]*dn2dx[1][j];
    cross[1][j] = dn1dx[2][j]*n2[0] -dn1dx[0][j]*n2[2] +
      n1[2]*dn2dx[0][j]-n1[0]*dn2dx[2][j];
    cross[2][j] = dn1dx[0][j]*n2[1] -dn1dx[1][j]*n2[0] +
      n1[0]*dn2dx[1][j]-n1[1]*dn2dx[0][j];
  }

  for (int i=0; i<3; i++)
    for (int j=0; j<3*N; j++) {
      double sum = 0;
      for (int l = 0; l<3; l++)
        sum += I3mn3n3T[i][l]*cross[l][j];
      dn3dx[i][j] = norm3*sum;
    }

  memory->destroy(cross);
  for (int l=0; l<N; l++)
    for (int i=0; i<3; i++)
      for (int j=0; j<3*N; j++)
        help2[i+3*l][j] = q0[3*l]*dn1dx[i][j] + q0[3*l+1]*dn2dx[i][j] +
          q0[3*l+2]*dn3dx[i][j];

  for (int j=0; j<N; j++)
    for (int l=0; l<N; l++)
      for (int i=0; i<3; i++)
        help2[i+3*l][i+3*j] += m[j]/m_all;

  //need transposed of help2:
  for (int ii = 0; ii < 3*N; ii++)
    for (int jj = 0; jj < 3*N; jj++)
      clist_derv[index_in_list][ii][jj] = help2[jj][ii];

  //free memory
  delete [] list_cluster;
  delete [] m;
  delete [] r;
  for (int i = 0; i<N; i++)
    delete [] del[i];
  delete [] del;
}

int FixFilterCorotate::pack_forward_comm(int n, int *list, double *buf,
                                         int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;
  double**f = atom->f;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    buf[m++] = f[j][0];
    buf[m++] = f[j][1];
    buf[m++] = f[j][2];
  }
  return m;
}

void FixFilterCorotate::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  double **f = atom->f;
  m = 0;
  last = first + n;

  for (i = first; i < last; i++)
  {
    f[i][0] = buf[m++];
    f[i][1] = buf[m++];
    f[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
 *  find a bond between global atom IDs n1 and n2 stored with local atom i
 *  if find it:
 *    if setflag = 0, return bond type
 *    if setflag = -1/1, set bond type to negative/positive and return 0
 *  if do not find it, return 0
 * ------------------------------------------------------------------------- */

int FixFilterCorotate::bondtype_findset(int i, tagint n1, tagint n2,
                                        int setflag)
{
  int m,nbonds;

  tagint *tag = atom->tag;
  tagint **bond_atom = atom->bond_atom;
  nbonds = atom->num_bond[i];

  for (m = 0; m < nbonds; m++) {
    if (n1 == tag[i] && n2 == bond_atom[i][m]) break;
    if (n1 == bond_atom[i][m] && n2 == tag[i]) break;
  }

  if (m < nbonds) {
    if (setflag == 0)
      return atom->bond_type[i][m];

    if ((setflag < 0 && atom->bond_type[i][m] > 0) ||
      (setflag > 0 && atom->bond_type[i][m] < 0))
      atom->bond_type[i][m] = -atom->bond_type[i][m];

  }

  return 0;
}

/* ----------------------------------------------------------------------
 *  find an angle with global end atom IDs n1 and n2 stored with local atom i
 *  if find it:
 *    if setflag = 0, return angle type
 *    if setflag = -1/1, set angle type to negative/positive and return 0
 *  if do not find it, return 0
 * ------------------------------------------------------------------------- */

int FixFilterCorotate::angletype_findset(int i, tagint n1, tagint n2,
                                         int setflag)
{
  int m,nangles;

  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom3 = atom->angle_atom3;
  nangles = atom->num_angle[i];

  for (m = 0; m < nangles; m++) {
    if (n1 == angle_atom1[i][m] && n2 == angle_atom3[i][m]) break;
    if (n1 == angle_atom3[i][m] && n2 == angle_atom1[i][m]) break;
  }

  if (m < nangles) {
    if (setflag == 0) {
      return atom->angle_type[i][m];
    }

    if ((setflag < 0 && atom->angle_type[i][m] > 0) ||
      (setflag > 0 && atom->angle_type[i][m] < 0))
      atom->angle_type[i][m] = -atom->angle_type[i][m];
  }

  return 0;
}

/* ----------------------------------------------------------------------
 *  allocate local atom-based arrays
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::grow_arrays(int nmax)
{
  memory->grow(array_atom,nmax,3,"FilterCorotate:peratomarray");
  memory->grow(shake_flag,nmax,"FilterCorotate::shake_flag");
  memory->grow(shake_atom,nmax,5,"FilterCorotate::shake_atom");
  memory->grow(shake_type,nmax,4,"FilterCorotate::shake_type");
}

/* ----------------------------------------------------------------------
 *  memory usage of local atom-based arrays
 * ------------------------------------------------------------------------- */

double FixFilterCorotate::memory_usage()
{
  double bytes = 0;
  //GROW:
  bytes += (double)3*sizeof(double) + 5*sizeof(tagint) + 5*sizeof(int);
  //clist
  bytes += (double)13*atom->nlocal*sizeof(int);
  bytes += (double)15*16*nlocal*sizeof(double);

  //fixed:
  int nb = atom->nbondtypes+1;
  int na = atom->nangletypes+1;
  int nt = atom->ntypes+1;
  bytes += (double)(nb+na+nt)*sizeof(int);
  bytes += (double)(nt-1+nb+na+15*15+18+10*15)*sizeof(double);

  return bytes;
}


/* ----------------------------------------------------------------------
 *  copy values within local atom-based arrays
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::copy_arrays(int i, int j, int /*delflag*/)
{
  int flag = shake_flag[j] = shake_flag[i];
  if (flag == 1) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_atom[j][2] = shake_atom[i][2];
    shake_type[j][0] = shake_type[i][0];
    shake_type[j][1] = shake_type[i][1];
    shake_type[j][2] = shake_type[i][2];
  } else if (flag == 2) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_type[j][0] = shake_type[i][0];
  } else if (flag == 3) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_atom[j][2] = shake_atom[i][2];
    shake_type[j][0] = shake_type[i][0];
    shake_type[j][1] = shake_type[i][1];
    shake_type[j][2] = shake_type[i][2];
  } else if (flag == 4) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_atom[j][2] = shake_atom[i][2];
    shake_atom[j][3] = shake_atom[i][3];
    shake_type[j][0] = shake_type[i][0];
    shake_type[j][1] = shake_type[i][1];
    shake_type[j][2] = shake_type[i][2];
  } else if (flag == 5) {
    shake_atom[j][0] = shake_atom[i][0];
    shake_atom[j][1] = shake_atom[i][1];
    shake_atom[j][2] = shake_atom[i][2];
    shake_atom[j][3] = shake_atom[i][3];
    shake_atom[j][4] = shake_atom[i][4];
    shake_type[j][0] = shake_type[i][0];
    shake_type[j][1] = shake_type[i][1];
    shake_type[j][2] = shake_type[i][2];
    shake_type[j][3] = shake_type[i][3];
  }
}

/* ----------------------------------------------------------------------
 *  initialize one atom's array values, called when atom is created
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::set_arrays(int i)
{
  shake_flag[i] = 0;
}

/* ----------------------------------------------------------------------
 *  update one atom's array values
 *  called when molecule is created from fix gcmc
 * ------------------------------------------------------------------------- */

void FixFilterCorotate::update_arrays(int i, int atom_offset)
{
  int flag = shake_flag[i];

  if (flag == 1) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
    shake_atom[i][2] += atom_offset;
  } else if (flag == 2) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
  } else if (flag == 3) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
    shake_atom[i][2] += atom_offset;
  } else if (flag == 4) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
    shake_atom[i][2] += atom_offset;
    shake_atom[i][3] += atom_offset;
  } else if (flag == 5) {
    shake_atom[i][0] += atom_offset;
    shake_atom[i][1] += atom_offset;
    shake_atom[i][2] += atom_offset;
    shake_atom[i][3] += atom_offset;
    shake_atom[i][4] += atom_offset;
  }
}


/* ----------------------------------------------------------------------
 *  pack values in local atom-based arrays for exchange with another proc
 * ------------------------------------------------------------------------- */

int FixFilterCorotate::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = shake_flag[i];
  int flag = shake_flag[i];
  if (flag == 1) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_atom[i][2];
    buf[m++] = shake_type[i][0];
    buf[m++] = shake_type[i][1];
    buf[m++] = shake_type[i][2];
  } else if (flag == 2) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_type[i][0];
  } else if (flag == 3) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_atom[i][2];
    buf[m++] = shake_type[i][0];
    buf[m++] = shake_type[i][1];
    buf[m++] = shake_type[i][2];
  } else if (flag == 4) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_atom[i][2];
    buf[m++] = shake_atom[i][3];
    buf[m++] = shake_type[i][0];
    buf[m++] = shake_type[i][1];
    buf[m++] = shake_type[i][2];
  } else if (flag == 5) {
    buf[m++] = shake_atom[i][0];
    buf[m++] = shake_atom[i][1];
    buf[m++] = shake_atom[i][2];
    buf[m++] = shake_atom[i][3];
    buf[m++] = shake_atom[i][4];
    buf[m++] = shake_type[i][0];
    buf[m++] = shake_type[i][1];
    buf[m++] = shake_type[i][2];
    buf[m++] = shake_type[i][3];
  }
  return m;
}

/* ----------------------------------------------------------------------
 *  unpack values in local atom-based arrays from exchange with another proc
 * ------------------------------------------------------------------------- */

int FixFilterCorotate::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  int flag = shake_flag[nlocal] = static_cast<int> (buf[m++]);
  if (flag == 1) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][2] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
    shake_type[nlocal][1] = static_cast<int> (buf[m++]);
    shake_type[nlocal][2] = static_cast<int> (buf[m++]);
  } else if (flag == 2) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
  } else if (flag == 3) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][2] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
    shake_type[nlocal][1] = static_cast<int> (buf[m++]);
    shake_type[nlocal][2] = static_cast<int> (buf[m++]);
  } else if (flag == 4) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][2] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][3] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
    shake_type[nlocal][1] = static_cast<int> (buf[m++]);
    shake_type[nlocal][2] = static_cast<int> (buf[m++]);
  } else if (flag == 5) {
    shake_atom[nlocal][0] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][1] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][2] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][3] = static_cast<tagint> (buf[m++]);
    shake_atom[nlocal][4] = static_cast<tagint> (buf[m++]);
    shake_type[nlocal][0] = static_cast<int> (buf[m++]);
    shake_type[nlocal][1] = static_cast<int> (buf[m++]);
    shake_type[nlocal][2] = static_cast<int> (buf[m++]);
    shake_type[nlocal][3] = static_cast<int> (buf[m++]);
  }
  return m;
}

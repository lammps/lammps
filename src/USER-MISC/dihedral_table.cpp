/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Andrew Jewett (jewett.aij  g m ail)
                        The cyclic tridiagonal matrix solver was borrowed from
                          the "tridiag.c" written by Gerard Jungman for GSL
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#define NDEBUG //(disable assert())
#include <cassert>
using namespace std;

#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "dihedral_table.h"

// Additional header files available for numerical debugging:
//#define DIH_DEBUG_NUM
//#ifdef DIH_DEBUG_NUM
//#include "dihedral_table_dbg.h"    //in "supporting_files/debug_numerical/"
//#endif
// and pointless posturing:
//#include "dihedral_table_nd_mod.h" //in "supporting_files/nd/"

using namespace LAMMPS_NS;
using namespace DIHEDRAL_TABLE_NS;

/* ---------------------------------------------------------------------- */

DihedralTable::DihedralTable(LAMMPS *lmp) : Dihedral(lmp)
{
  ntables = 0;
  tables = NULL;
}

/* ---------------------------------------------------------------------- */

DihedralTable::~DihedralTable()
{
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    //memory->destroy(phi0); <- equilibrium angles not supported
    memory->destroy(tabindex);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralTable::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double edihedral,f1[3],f2[3],f3[3],f4[3];

  double **x = atom->x;
  double **f = atom->f;

  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  // The dihedral angle "phi" is the angle between n123 and n234
  // the planes defined by atoms i1,i2,i3, and i2,i3,i4.
  //
  // Definitions of vectors: vb12, vb23, vb34, perp12on23 
  //                         proj12on23, perp43on32, proj43on32
  //
  //  Note: The positions of the 4 atoms are labeled x[i1], x[i2], x[i3], x[i4]
  //        (which are also vectors)
  //
  //             proj12on23                          proj34on23          
  //             --------->                         ----------->         
  //                           .                                         
  //                          .                                          
  //                         .                                  
  //                  x[i2] .                       x[i3]                 
  //    .                __@----------vb23-------->@ . . . .           . 
  //   /|\                /|                        \                  |
  //    |                /                           \                 |
  //    |               /                             \                |
  // perp12vs23        /                               \               |
  //    |             /                                 \          perp34vs23
  //    |          vb12                                  \             |
  //    |           /                                   vb34           |
  //    |          /                                       \           |
  //    |         /                                         \          |
  //    |        /                                           \         |
  //            @                                             \        |
  //                                                          _\|     \|/
  //         x[i1]                                              @     
  //                                                           
  //                                                           x[i4]   
  //

  double vb12[g_dim]; // displacement vector from atom i1 towards atom i2
  //     vb12[d]       = x[i2][d] - x[i1][d]      (for d=0,1,2)
  double vb23[g_dim]; // displacement vector from atom i2 towards atom i3
  //     vb23[d]       = x[i3][d] - x[i2][d]      (for d=0,1,2)
  double vb34[g_dim]; // displacement vector from atom i3 towards atom i4
  //     vb34[d]       = x[i4][d] - x[i3][d]      (for d=0,1,2)

  //  n123 & n234: These two unit vectors are normal to the planes
  //               defined by atoms 1,2,3 and 2,3,4.
  double n123[g_dim]; //n123=vb12 x vb23 / |vb12 x vb23|  ("x" is cross product)
  double n234[g_dim]; //n234=vb34 x vb23 / |vb34 x vb23|  ("x" is cross product)

  double proj12on23[g_dim];
  //    proj12on23[d] = (vb23[d]/|vb23|) * DotProduct(vb12,vb23)/|vb12|*|vb23|
  double proj34on23[g_dim];
  //    proj34on23[d] = (vb34[d]/|vb23|) * DotProduct(vb34,vb23)/|vb34|*|vb23|
  double perp12on23[g_dim];
  //    perp12on23[d] = v12[d] - proj12on23[d]
  double perp34on23[g_dim];
  //    perp34on23[d] = v34[d] - proj34on23[d]


  edihedral = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;


  for (n = 0; n < ndihedrallist; n++) {

    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // ------ Step 1: Compute the dihedral angle "phi" ------
    //

    // Phi() calculates the dihedral angle.
    // This function also calculates the vectors: 
    // vb12, vb23, vb34, n123, and n234, which we will need later.


    double phi = Phi(x[i1], x[i2], x[i3], x[i4], domain, 
                     vb12, vb23, vb34, n123, n234);


    // ------ Step 2: Compute the gradient of phi with atomic position: ------
    //
    // Gradient variables:
    //
    // dphi_dx1, dphi_dx2, dphi_dx3, dphi_dx4 are the gradients of phi with 
    // respect to the atomic positions of atoms i1, i2, i3, i4, respectively.
    // As an example, consider dphi_dx1.  The d'th element is:
    double dphi_dx1[g_dim]; //                 d phi
    double dphi_dx2[g_dim]; // dphi_dx1[d] = ----------    (partial derivatives)
    double dphi_dx3[g_dim]; //               d x[i1][d]
    double dphi_dx4[g_dim]; //where d=0,1,2 corresponds to x,y,z  (if g_dim==3)

    double dot123             = DotProduct(vb12, vb23);
    double dot234             = DotProduct(vb23, vb34);
    double L23sqr             = DotProduct(vb23, vb23);
    double L23                = sqrt(L23sqr);   // (central bond length)
    double inv_L23sqr = 0.0;
    double inv_L23    = 0.0;
    if (L23sqr != 0.0) {
      inv_L23sqr = 1.0 / L23sqr;
      inv_L23 = 1.0 / L23;
    }
    double neg_inv_L23        = -inv_L23;
    double dot123_over_L23sqr = dot123 * inv_L23sqr;
    double dot234_over_L23sqr = dot234 * inv_L23sqr;

    for (int d=0; d < g_dim; ++d) {
      // See figure above for a visual definitions of these vectors:
      proj12on23[d] = vb23[d] * dot123_over_L23sqr;
      proj34on23[d] = vb23[d] * dot234_over_L23sqr;
      perp12on23[d] = vb12[d] - proj12on23[d];
      perp34on23[d] = vb34[d] - proj34on23[d];
    }


    // --- Compute the gradient vectors dphi/dx1 and dphi/dx4: ---

    // These two gradients point in the direction of n123 and n234,
    // and are scaled by the distances of atoms 1 and 4 from the central axis.
    // Distance of atom 1 to central axis:
    double perp12on23_len = sqrt(DotProduct(perp12on23, perp12on23));
    // Distance of atom 4 to central axis:
    double perp34on23_len = sqrt(DotProduct(perp34on23, perp34on23));

    double inv_perp12on23 = 0.0;
    if (perp12on23_len != 0.0) inv_perp12on23 = 1.0 / perp12on23_len;
    double inv_perp34on23 = 0.0;
    if (perp34on23_len != 0.0) inv_perp34on23 = 1.0 / perp34on23_len;

    for (int d=0; d < g_dim; ++d) {
      dphi_dx1[d] =  n123[d] * inv_perp12on23;
      dphi_dx4[d] =  n234[d] * inv_perp34on23;
    }

    // --- Compute the gradient vectors dphi/dx2 and dphi/dx3: ---
    //
    // This is more tricky because atoms 2 and 3 are shared by both planes 
    // 123 and 234 (the angle between which defines "phi").  Moving either 
    // one of these atoms effects both the 123 and 234 planes
    // Both the 123 and 234 planes intersect with the plane perpendicular to the
    // central bond axis (vb23).  The two lines where these intersections occur
    // will shift when you move either atom 2 or atom 3.  The angle between
    // these lines is the dihedral angle, phi.  We can define four quantities:
    // dphi123_dx2 is the change in "phi" due to the movement of the 123 plane
    //             ...as a result of moving atom 2.
    // dphi234_dx2 is the change in "phi" due to the movement of the 234 plane
    //             ...as a result of moving atom 2.
    // dphi123_dx3 is the change in "phi" due to the movement of the 123 plane
    //             ...as a result of moving atom 3.
    // dphi234_dx3 is the change in "phi" due to the movement of the 234 plane
    //             ...as a result of moving atom 3.

    double proj12on23_len = dot123 * inv_L23;
    double proj34on23_len = dot234 * inv_L23;
    // Interpretation:
    //The magnitude of "proj12on23_len" is the length of the proj12on23 vector.
    //The sign is positive if it points in the same direction as the central
    //bond (vb23).  Otherwise it is negative.  The same goes for "proj34on23".
    //(In the example figure in the comment above, both variables are positive.)

    // The forumula used in the 8 lines below explained here:
    //   "supporting_information/doc/gradient_formula_explanation/"
    double dphi123_dx2_coef = neg_inv_L23 * (L23 + proj12on23_len);
    double dphi234_dx2_coef = inv_L23 * proj34on23_len;

    double dphi234_dx3_coef = neg_inv_L23 * (L23 + proj34on23_len);
    double dphi123_dx3_coef = inv_L23 * proj12on23_len;

    for (int d=0; d < g_dim; ++d) {
      // Recall that the n123 and n234 plane normal vectors are proportional to
      // the dphi/dx1 and dphi/dx2 gradients vectors 
      // It turns out we can save slightly more CPU cycles by expressing
      // dphi/dx2 and dphi/dx3 as linear combinations of dphi/dx1 and dphi/dx2 
      // which we computed already (instead of n123 & n234).
      dphi_dx2[d] = dphi123_dx2_coef*dphi_dx1[d] + dphi234_dx2_coef*dphi_dx4[d];
      dphi_dx3[d] = dphi123_dx3_coef*dphi_dx1[d] + dphi234_dx3_coef*dphi_dx4[d];
    }




    #ifdef DIH_DEBUG_NUM
    // ----- Numerical test? -----

    cerr << "  -- testing gradient for dihedral (n="<<n<<") for atoms ("
         << i1 << "," << i2 << "," << i3 << "," << i4 << ") --" << endl;

    PrintGradientComparison(*this, dphi_dx1, dphi_dx2, dphi_dx3, dphi_dx4, 
                            domain, x[i1], x[i2], x[i3], x[i4]);

    for (int d=0; d < g_dim; ++d) {
      // The sum of all the gradients should be near 0. (translational symmetry)
      cerr <<"sum_gradients["<<d<<"]="<<dphi_dx1[d]<<"+"<<dphi_dx2[d]<<"+"<<dphi_dx3[d]<<"+"<<dphi_dx4[d]<<"="<<dphi_dx1[d]+dphi_dx2[d]+dphi_dx3[d]+dphi_dx4[d]<<endl;
      // These should sum to zero
      assert(abs(dphi_dx1[d]+dphi_dx2[d]+dphi_dx3[d]+dphi_dx4[d]) < 0.0002/L23);
    }
    #endif // #ifdef DIH_DEBUG_NUM




    // ----- Step 3: Calculate the energy and force in the phi direction -----

    // tabulated force & energy
    double u, m_du_dphi; //u = energy.   m_du_dphi = "minus" du/dphi
    assert((0.0 <= phi) && (phi <= TWOPI));

    uf_lookup(type, phi, u, m_du_dphi);


    if (eflag) edihedral = u;

    // ----- Step 4: Calculate the force direction in real space -----

    // chain rule:
    //          d U          d U      d phi
    // -f  =   -----   =    -----  *  -----
    //          d x         d phi      d x
    for(int d=0; d < g_dim; ++d) {
      f1[d] = m_du_dphi * dphi_dx1[d];
      f2[d] = m_du_dphi * dphi_dx2[d];
      f3[d] = m_du_dphi * dphi_dx3[d];
      f4[d] = m_du_dphi * dphi_dx4[d];
    }

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,
               nlocal,newton_bond,edihedral,
               f1,f3,f4,
               vb12[0],vb12[1],vb12[2],
               vb23[0],vb23[1],vb23[2],
               vb34[0],vb34[1],vb34[2]);
  }
} // void DihedralTable::compute()







// single() calculates the dihedral-angle energy of atoms i1, i2, i3, i4.
double DihedralTable::single(int type, int i1, int i2, int i3, int i4)
{
  double vb12[g_dim];
  double vb23[g_dim];
  double vb34[g_dim];
  double n123[g_dim];
  double n234[g_dim];

  double **x = atom->x;

  double phi = Phi(x[i1], x[i2], x[i3], x[i4], domain, 
                   vb12, vb23, vb34, n123, n234);

  double u;
  u_lookup(type, phi, u); //Calculate the energy, and store it in "u"

  return u;
}


/* ---------------------------------------------------------------------- */



void DihedralTable::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(tabindex,n+1,"dihedral:tabindex");
  //memory->create(phi0,n+1,"dihedral:phi0"); <-equilibrium angles not supported
  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}


/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void DihedralTable::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal dihedral_style command");

  if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else error->all(FLERR,"Unknown table style in dihedral style table");

  tablength = force->inumeric(arg[1]);
  if (tablength < 3) 
    error->all(FLERR,"Illegal number of dihedral table entries");
  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}



/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */


void DihedralTable::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal dihedral_coeff command");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);
  
  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *) 
    memory->srealloc(tables,(ntables+1)*sizeof(Table), "dihedral:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[1],arg[2]);
  bcast_table(tb);


  // --- check the angle data for range errors ---
  // ---  and resolve issues with periodicity  ---

  if (tb->ninput < 3) {
    string err_msg;
    err_msg = string("Invalid dihedral table length (")
      + string(arg[2]) + string(").");
    error->one(FLERR,err_msg.c_str());
  }

  // check for monotonicity
  for (int i=0; i < tb->ninput-1; i++) {
    if (tb->phifile[i] >= tb->phifile[i+1]) {
      string err_msg = 
        string("Dihedral table values are not increasing (") + 
        string(arg[2]) + string(").");
      error->all(FLERR,err_msg.c_str());
    }
  }

  // check the range of angles
  double philo = tb->phifile[0];
  double phihi = tb->phifile[tb->ninput-1];
  if (tb->use_degrees) {
    if ((phihi - philo) >= 360) {
      string err_msg;
      err_msg = string("Dihedral table angle range must be < 360 degrees (")
	+string(arg[2]) + string(").");
      error->all(FLERR,err_msg.c_str());
    }
  }
  else {
    if ((phihi - philo) >= TWOPI) {
      string err_msg;
      err_msg = string("Dihedral table angle range must be < 2*PI radians (")
	+ string(arg[2]) + string(").");
      error->all(FLERR,err_msg.c_str());
    }
  }

  // convert phi from degrees to radians
  if (tb->use_degrees) {
    for (int i=0; i < tb->ninput; i++) {
      tb->phifile[i] *= PI/180.0;
      // I assume that if angles are in degrees, then the forces (f=dU/dphi)
      // are specified with "phi" in radians as well.
      tb->ffile[i] *= 180.0/PI;
    }
  }

  // We want all the phi dihedral angles to lie in the range from 0 to 2*PI.
  // But I don't want to restrict users to input their data in this range.
  // We also want the angles to be sorted in increasing order.
  // This messy code fixes these problems with the user's data:
  {
    double *phifile_tmp = new double [tb->ninput];  //temporary arrays 
    double *ffile_tmp = new double [tb->ninput];  //used for sorting
    double *efile_tmp = new double [tb->ninput];

    // After re-imaging, does the range of angles cross the 0 or 2*PI boundary?
    // If so, find the discontinuity:
    int i_discontinuity = tb->ninput; 
    for (int i=0; i < tb->ninput; i++) {
      double phi   = tb->phifile[i];
      // Add a multiple of 2*PI to phi until it lies in the range [0, 2*PI).
      phi -= TWOPI * floor(phi/TWOPI);
      phifile_tmp[i] = phi;
      efile_tmp[i] = tb->efile[i];
      ffile_tmp[i] = tb->ffile[i];
      if ((i>0) && (phifile_tmp[i] < phifile_tmp[i-1])) {
        //There should only be at most one discontinuity, because we have
        //insured that the data was sorted before imaging, and because the
        //range of angle values does not exceed 2*PI.
        assert(i_discontinuity == tb->ninput); //<-should only happen once
        i_discontinuity = i;
      }
    }

    int I = 0;
    for (int i = i_discontinuity; i < tb->ninput; i++) {
      tb->phifile[I] = phifile_tmp[i];
      tb->efile[I] = efile_tmp[i];
      tb->ffile[I] = ffile_tmp[i];
      assert((I==0) || (tb->phifile[I] > tb->phifile[I-1]));
      I++;
    }
    for (int i = 0; i < i_discontinuity; i++) {
      tb->phifile[I] = phifile_tmp[i];
      tb->efile[I] = efile_tmp[i];
      tb->ffile[I] = ffile_tmp[i];
      assert((I==0) || (tb->phifile[I] > tb->phifile[I-1]));
      I++;
    }
    assert(I == tb->ninput);
  }

  // spline read-in and compute r,e,f vectors within table

  spline_table(tb);
  compute_table(tb);

  // Optional: allow the user to print out the interpolated spline tables

  if (me == 0) {
    if (checkU_fname && (strlen(checkU_fname) != 0))
    {
      ofstream checkU_file;
      checkU_file.open(checkU_fname, ios::out);
      for (int i=0; i < tablength; i++) {
        double phi = i*TWOPI/tablength;
        double u = tb->e[i];
        if (tb->use_degrees)
          phi *= 180.0/PI;
        checkU_file << phi << " " << u << "\n";
      }
      checkU_file.close();
    }
    if (checkF_fname && (strlen(checkF_fname) != 0))
    {
      ofstream checkF_file;
      checkF_file.open(checkF_fname, ios::out);
      for (int i=0; i < tablength; i++)
      {
        double phi = i*TWOPI/tablength;
        double f;
        if ((tabstyle == SPLINE) && (tb->f_unspecified)) {
          double dU_dphi = 
            // (If the user did not specify the forces now, AND the user
	    //  selected the "spline" option, (as opposed to "linear")
	    //  THEN the tb->f array is uninitialized, so there's 
	    //  no point to print out the contents of the tb->f[] array. 
	    //  Instead, later on, we will calculate the force using the 
	    //  -cyc_splintD() routine to calculate the derivative of the
	    //  energy spline, using the energy data (tb->e[]).
            //  To be nice and report something, I do the same thing here.)
            cyc_splintD(tb->phi, tb->e, tb->e2, tablength, TWOPI,phi);
          f = -dU_dphi;
        }
        else
	  // Otherwise we calculated the tb->f[] array.  Report its contents.
          f = tb->f[i];
        if (tb->use_degrees) {
          phi *= 180.0/PI;
	  // If the user wants degree angle units, we should convert our 
	  // internal force tables (in energy/radians) to (energy/degrees)
	  f *= PI/180.0;
	}
        checkF_file << phi << " " << f << "\n";
      }
      checkF_file.close();
    } // if (checkF_fname && (strlen(checkF_fname) != 0))
  } // if (me == 0)


  // store ptr to table in tabindex
  int count = 0;
  for (int i = ilo; i <= ihi; i++)
  {
    tabindex[i] = ntables;
    //phi0[i] = tb->phi0; <- equilibrium dihedral angles not supported
    setflag[i] = 1;
    count++;
  }
  ntables++;

  if (count == 0) 
    error->all(FLERR,"Illegal dihedral_coeff command");
} //DihedralTable::coeff()


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void DihedralTable::write_restart(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&tablength,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void DihedralTable::read_restart(FILE *fp)
{
  if (comm->me == 0) {
    fread(&tabstyle,sizeof(int),1,fp);
    fread(&tablength,sizeof(int),1,fp);
  }
  MPI_Bcast(&tabstyle,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tablength,1,MPI_INT,0,world);

  allocate();
}


/* ---------------------------------------------------------------------- */

void DihedralTable::null_table(Table *tb)
{
  tb->phifile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->phi = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
}

/* ---------------------------------------------------------------------- */

void DihedralTable::free_table(Table *tb)
{
  memory->destroy(tb->phifile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);
  
  memory->destroy(tb->phi);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->e2);
  memory->destroy(tb->f2);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void DihedralTable::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];

  // open file

  FILE *fp = fopen(file,"r");
  if (fp == NULL) {
    string err_msg = string("Cannot open file ") + string(file);
    error->one(FLERR,err_msg.c_str());
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL)
      error->one(FLERR,"Did not find keyword in dihedral table file");
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    if (strstr(line,keyword) == line) break;               // matching keyword
    fgets(line,MAXLINE,fp);                         // no match, skip section
    param_extract(tb,line);
    fgets(line,MAXLINE,fp);
    for (int i = 0; i < tb->ninput; i++)
      fgets(line,MAXLINE,fp);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  fgets(line,MAXLINE,fp);
  param_extract(tb,line);
  memory->create(tb->phifile,tb->ninput,"dihedral:phifile");
  memory->create(tb->efile,tb->ninput,"dihedral:efile");
  memory->create(tb->ffile,tb->ninput,"dihedral:ffile");

  // read a,e,f table values from file

  int itmp;
  for (int i = 0; i < tb->ninput; i++) {
    fgets(line,MAXLINE,fp);
    if (tb->f_unspecified)
      sscanf(line,"%d %lg %lg",
             &itmp,&tb->phifile[i],&tb->efile[i]);
    else
      sscanf(line,"%d %lg %lg %lg",
             &itmp,&tb->phifile[i],&tb->efile[i],&tb->ffile[i]);
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file,f2file.  
   I also perform a crude check for force & energy consistency.
------------------------------------------------------------------------- */

void DihedralTable::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"dihedral:e2file");
  memory->create(tb->f2file,tb->ninput,"dihedral:f2file");

  cyc_spline(tb->phifile, tb->efile, tb->ninput, TWOPI, tb->e2file);

  if (! tb->f_unspecified) {
    cyc_spline(tb->phifile, tb->ffile, tb->ninput, TWOPI, tb->f2file);
  }

  // CHECK to help make sure the user calculated forces in a way 
  // which is grossly numerically consistent with the energy table.
  //  --------------------------------------------------
  //             This is an ugly piece of code
  //  --------------------------------------------------
  if (! tb->f_unspecified) {
    int num_disagreements = 0;
    for (int i=0; i<tb->ninput; i++) {

      // Calculate what the force should be at the control points 
      // by using linear interpolation of the derivatives of the energy:

      double phi_i = tb->phifile[i];

      // First deal with periodicity (I hate this part)
      double phi_im1, phi_ip1;
      int im1 = i-1;
      if (im1 < 0) {
        im1 += tb->ninput;
        phi_im1 = tb->phifile[im1] - TWOPI;
      }
      else 
        phi_im1 = tb->phifile[im1];
      int ip1 = i+1;
      if (ip1 >= tb->ninput) {
        ip1 -= tb->ninput;
        phi_ip1 = tb->phifile[ip1] + TWOPI;
      }
      else 
        phi_ip1 = tb->phifile[ip1];

      // Now calculate the midpoints above and below phi_i = tb->phifile[i]
      double phi_lo= 0.5*(phi_im1 + phi_i); //midpoint between phi_im1 and phi_i
      double phi_hi= 0.5*(phi_i + phi_ip1); //midpoint between phi_i and phi_ip1

      // Use a linear approximation to the derivative at these two midpoints
      double dU_dphi_lo = 
        (tb->efile[i] - tb->efile[im1])
        / 
	(phi_i - phi_im1);
      double dU_dphi_hi = 
        (tb->efile[ip1] - tb->efile[i])
        / 
	(phi_ip1 - phi_i);

      // Now calculate the derivative at position 
      // phi_i (=tb->phifile[i]) using linear interpolation

      double a = (phi_i - phi_lo) / (phi_hi - phi_lo);
      double b = (phi_hi - phi_i) / (phi_hi - phi_lo);
      double dU_dphi = a*dU_dphi_lo  +  b*dU_dphi_hi;
      double f = -dU_dphi;
      // Alternately, we could use spline interpolation instead:
      // double f = - splintD(tb->phifile, tb->efile, tb->e2file,
      //                      tb->ninput, TWOPI, tb->phifile[i]);
      // This is the way I originally did it, but I trust
      // the ugly simple linear way above better.
      // Recall this entire block of code doess not calculate
      // anything important.  It does not have to be perfect.
      // We are only checking for stupid user errors here.

      if ((f != 0.0) && 
          (tb->ffile[i] != 0.0) &&
          ((f/tb->ffile[i] < 0.5) || (f/tb->ffile[i] > 2.0))) {
        num_disagreements++;
      }
    } // for (int i=0; i<tb->ninput; i++)

    if (num_disagreements > tb->ninput/2) {
      string msg("Dihedral table has inconsistent forces and energies. (Try \"NOF\".)\n");
      error->all(FLERR,msg.c_str());
    }

  } // check for consistency if (! tb->f_unspecified)

} // DihedralTable::spline_table()


/* ----------------------------------------------------------------------
   compute a,e,f vectors from splined values
------------------------------------------------------------------------- */

void DihedralTable::compute_table(Table *tb)
{
  //delta = table spacing in dihedral angle for tablength (cyclic) bins
  tb->delta = TWOPI / tablength;
  tb->invdelta = 1.0/tb->delta;
  tb->deltasq6 = tb->delta*tb->delta / 6.0;

  // N evenly spaced bins in dihedral angle from 0 to 2*PI
  // phi,e,f = value at lower edge of bin
  // de,df values = delta values of e,f (cyclic, in this case)
  // phi,e,f,de,df are arrays containing "tablength" number of entries

  memory->create(tb->phi,tablength,"dihedral:phi");
  memory->create(tb->e,tablength,"dihedral:e");
  memory->create(tb->de,tablength,"dihedral:de");
  memory->create(tb->f,tablength,"dihedral:f");
  memory->create(tb->df,tablength,"dihedral:df");
  memory->create(tb->e2,tablength,"dihedral:e2");
  memory->create(tb->f2,tablength,"dihedral:f2");

  // Use cubic spline interpolation to calculate the entries in the 
  // internal table. (This is true regardless...even if tabstyle!=SPLINE.)
  for (int i = 0; i < tablength; i++) {
    double phi = i*tb->delta;
    tb->phi[i] = phi;
    tb->e[i]= cyc_splint(tb->phifile,tb->efile,tb->e2file,tb->ninput,TWOPI,phi);
    if (! tb->f_unspecified)
      tb->f[i] = cyc_splint(tb->phifile,tb->ffile,tb->f2file,tb->ninput,TWOPI,phi);
  }

  if (tabstyle == LINEAR) {
    if (tb->f_unspecified)
    {
      // In the linear case, if the user did not specify the forces, then we
      // must generate the "f" array. Do this using linear interpolation of the
      // e array (which itself was generated using spline interpolation above).
      for (int i = 0; i < tablength; i++) {
        int im1 = i-1; if (im1 < 0) im1 += tablength;
        int ip1 = i+1; if (ip1 >= tablength) ip1 -= tablength;
        double dedx = (tb->e[ip1] - tb->e[im1]) / (2.0 * tb->delta);
        // (This is the average of the linear slopes on either side of the node.
        //  Note that the nodes in the internal table are evenly spaced.)
        tb->f[i] = -dedx;
      }
    }

    // Fill the linear interpolation tables (de, df)
    for (int i = 0; i < tablength; i++) {
      int ip1 = i+1; if (ip1 >= tablength) ip1 -= tablength;
      tb->de[i] = tb->e[ip1] - tb->e[i];
      tb->df[i] = tb->f[ip1] - tb->f[i];
    }
  } // if (tabstyle == LINEAR)

  cyc_spline(tb->phi, tb->e, tablength, TWOPI, tb->e2);
  if (! tb->f_unspecified)
    cyc_spline(tb->phi, tb->f, tablength, TWOPI, tb->f2);
}


/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value NOF DEGREES RADIANS
   N is required, other params are optional
------------------------------------------------------------------------- */

void DihedralTable::param_extract(Table *tb, char *line)
{
  //tb->theta0 = 180.0; <- equilibrium angles not supported
  tb->ninput = 0;
  tb->f_unspecified = false; //default
  tb->use_degrees   = true;  //default
  strcpy(checkU_fname, "");
  strcpy(checkF_fname, "");

  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
    }
    else if (strcmp(word,"NOF") == 0) {
      tb->f_unspecified = true;
    }
    else if ((strcmp(word,"DEGREES") == 0) || (strcmp(word,"degrees") == 0)) {
      tb->use_degrees = true;
    }
    else if ((strcmp(word,"RADIANS") == 0) || (strcmp(word,"radians") == 0)) {
      tb->use_degrees = false;
    }
    else if (strcmp(word,"CHECKU") == 0) {
      word = strtok(NULL," \t\n\r\f");
      strncpy(checkU_fname, word, MAXLINE);
    }
    else if (strcmp(word,"CHECKF") == 0) {
      word = strtok(NULL," \t\n\r\f");
      strncpy(checkF_fname, word, MAXLINE);
    }
    // COMMENTING OUT:  equilibrium angles are not supported
    //else if (strcmp(word,"EQ") == 0) {
    //  word = strtok(NULL," \t\n\r\f");
    //  tb->theta0 = atof(word);
    //} 
    else {
      string err_msg("Invalid keyword in dihedral angle table parameters");
      err_msg += string(" (") + string(word) + string(")");
      error->one(FLERR,err_msg.c_str());
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0) 
    error->one(FLERR,"Dihedral table parameters did not set N");

} // DihedralTable::param_extract()


/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,phifile,efile,ffile,
     f_unspecified,use_degrees
------------------------------------------------------------------------- */

void DihedralTable::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->phifile,tb->ninput,"dihedral:phifile");
    memory->create(tb->efile,tb->ninput,"dihedral:efile");
    memory->create(tb->ffile,tb->ninput,"dihedral:ffile");
  }

  MPI_Bcast(tb->phifile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);

  MPI_Bcast(&tb->f_unspecified,1,MPI_INT,0,world);
  MPI_Bcast(&tb->use_degrees,1,MPI_INT,0,world);

  // COMMENTING OUT:  equilibrium angles are not supported
  //MPI_Bcast(&tb->theta0,1,MPI_DOUBLE,0,world);
}

 
namespace LAMMPS_NS {
namespace DIHEDRAL_TABLE_NS {


// -------------------------------------------------------------------
// ---------    The function was stolen verbatim from the    ---------
// ---------    GNU Scientific Library (GSL, version 1.15)   ---------
// -------------------------------------------------------------------

/* Author: Gerard Jungman */
/* for description of method see [Engeln-Mullges + Uhlig, p. 96]
 *
 *      diag[0]  offdiag[0]             0   .....  offdiag[N-1]
 *   offdiag[0]     diag[1]    offdiag[1]   .....
 *            0  offdiag[1]       diag[2]
 *            0           0    offdiag[2]   .....
 *          ...         ...
 * offdiag[N-1]         ...
 *
 */
// -- (A non-symmetric version of this function is also available.) --

enum { //GSL status return codes.
  GSL_FAILURE  = -1,
  GSL_SUCCESS  = 0,
  GSL_ENOMEM   = 8,
  GSL_EZERODIV = 12,
  GSL_EBADLEN  = 19
};


int 
solve_cyc_tridiag(
  const double diag[], size_t d_stride,
  const double offdiag[], size_t o_stride,
  const double b[], size_t b_stride,
  double x[], size_t x_stride,
  size_t N)
{
  int status = GSL_SUCCESS;
  double * delta = (double *) malloc (N * sizeof (double));
  double * gamma = (double *) malloc (N * sizeof (double));
  double * alpha = (double *) malloc (N * sizeof (double));
  double * c = (double *) malloc (N * sizeof (double));
  double * z = (double *) malloc (N * sizeof (double));

  if (delta == 0 || gamma == 0 || alpha == 0 || c == 0 || z == 0) {
    cerr << "Internal Cyclic Spline Error: failed to allocate working space\n";
    exit(GSL_ENOMEM);
  }
  else
    {
      size_t i, j;
      double sum = 0.0;

      /* factor */

      if (N == 1) 
        {
          x[0] = b[0] / diag[0];
          return GSL_SUCCESS;
        }

      alpha[0] = diag[0];
      gamma[0] = offdiag[0] / alpha[0];
      delta[0] = offdiag[o_stride * (N-1)] / alpha[0];

      if (alpha[0] == 0) {
        status = GSL_EZERODIV;
      }

      for (i = 1; i < N - 2; i++)
        {
          alpha[i] = diag[d_stride * i] - offdiag[o_stride * (i-1)] * gamma[i - 1];
          gamma[i] = offdiag[o_stride * i] / alpha[i];
          delta[i] = -delta[i - 1] * offdiag[o_stride * (i-1)] / alpha[i];
          if (alpha[i] == 0) {
            status = GSL_EZERODIV;
          }
        }

      for (i = 0; i < N - 2; i++)
        {
          sum += alpha[i] * delta[i] * delta[i];
        }

      alpha[N - 2] = diag[d_stride * (N - 2)] - offdiag[o_stride * (N - 3)] * gamma[N - 3];

      gamma[N - 2] = (offdiag[o_stride * (N - 2)] - offdiag[o_stride * (N - 3)] * delta[N - 3]) / alpha[N - 2];

      alpha[N - 1] = diag[d_stride * (N - 1)] - sum - alpha[(N - 2)] * gamma[N - 2] * gamma[N - 2];

      /* update */
      z[0] = b[0];
      for (i = 1; i < N - 1; i++)
        {
          z[i] = b[b_stride * i] - z[i - 1] * gamma[i - 1];
        }
      sum = 0.0;
      for (i = 0; i < N - 2; i++)
        {
          sum += delta[i] * z[i];
        }
      z[N - 1] = b[b_stride * (N - 1)] - sum - gamma[N - 2] * z[N - 2];
      for (i = 0; i < N; i++)
        {
          c[i] = z[i] / alpha[i];
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = c[N - 1];
      x[x_stride * (N - 2)] = c[N - 2] - gamma[N - 2] * x[x_stride * (N - 1)];
      if (N >= 3)
        {
          for (i = N - 3, j = 0; j <= N - 3; j++, i--)
            {
              x[x_stride * i] = c[i] - gamma[i] * x[x_stride * (i + 1)] - delta[i] * x[x_stride * (N - 1)];
            }
        }
    }

  if (z != 0)
    free (z);
  if (c != 0)
    free (c);
  if (alpha != 0)
    free (alpha);
  if (gamma != 0)
    free (gamma);
  if (delta != 0)
    free (delta);

  if (status == GSL_EZERODIV) {
    cerr <<"Internal Cyclic Spline Error: Matrix must be positive definite.\n";
    exit(GSL_EZERODIV);
  }

  return status;
} //solve_cyc_tridiag()



/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void cyc_spline(double const *xa, 
                double const *ya, 
                int n,
                double period, 
                double *y2a)
{

  double *diag    = new double[n];   
  double *offdiag = new double[n];   
  double *rhs     = new double[n];
  double xa_im1, xa_ip1;

  // In the cyclic case, there are n equations with n unknows.
  // The for loop sets up the equations we need to solve.
  // Later we invoke the GSL tridiagonal matrix solver to solve them.

  for(int i=0; i < n; i++) {

    // I have to lookup xa[i+1] and xa[i-1].  This gets tricky because of 
    // periodic boundary conditions.  We handle that now.
    int im1 = i-1;
    if (im1<0) {
      im1 += n;
      xa_im1 = xa[im1] - period;
    }
    else
      xa_im1 = xa[im1];

    int ip1 = i+1;
    if (ip1>=n) {
      ip1 -= n;
      xa_ip1 = xa[ip1] + period;
    }
    else
      xa_ip1 = xa[ip1];

    // Recall that we want to find the y2a[] parameters (there are n of them).
    // To solve for them, we have a linear equation with n unknowns
    // (in the cyclic case that is).  For details, the non-cyclic case is 
    // explained in equation 3.3.7 in Numerical Recipes in C, p. 115.
    diag[i]    = (xa_ip1 - xa_im1) / 3.0;
    offdiag[i] = (xa_ip1 - xa[i]) / 6.0;
    rhs[i]     = ((ya[ip1] - ya[i]) / (xa_ip1 - xa[i])) -
                 ((ya[i] - ya[im1]) / (xa[i] - xa_im1));
  }

  // Ordinarily we would have to invert a large matrix (with potentially
  // thousands of rows and columns).  However because this matix happens
  // to be tridiagonal (and cyclic), we can use the following cheap method:
  solve_cyc_tridiag(diag, 1,
		    offdiag, 1,
		    rhs, 1,
		    y2a, 1,
		    n);

  delete [] diag;
  delete [] offdiag;
  delete [] rhs;

} // cyc_spline()



/* ---------------------------------------------------------------------- */

// cyc_splint(): Evaluates a spline at position x, with n control
//           points located at xa[], ya[], and with parameters y2a[]
//           The xa[] must be monotonically increasing and their 
//           range should not exceed period (ie xa[n-1] < xa[0] + period).  
//           x must lie in the range:  [(xa[n-1]-period), (xa[0]+period)]
//           "period" is typically 2*PI.
double cyc_splint(double const *xa, 
                  double const *ya, 
                  double const *y2a, 
                  int n, 
                  double period,
                  double x)
{
  int klo = -1;
  int khi = n;
  int k;
  double xlo = xa[n-1] - period;
  double xhi = xa[0] + period;
  assert((xlo <= x) && (x <= xhi));
  while (khi-klo > 1) {
    k = (khi+klo) >> 1; //(k=(khi+klo)/2)
    if (xa[k] > x) {
      khi = k;
      xhi = xa[k];
    }
    else {
      klo = k;
      xlo = xa[k];
    }
  }
  assert((xlo <= x) && (x <= xhi));

  if (khi == n) khi = 0;
  if (klo ==-1) klo = n-1;

  double h = xhi-xlo;
  double a = (xhi-x) / h;
  double b = (x-xlo) / h;
  double y = a*ya[klo] + b*ya[khi] + 
    ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;

  return y;

} // cyc_splint()



// cyc_splintD(): Evaluate the deriviative of a cyclic spline at position x, 
//           with n control points at xa[], ya[], with parameters y2a[].
//           The xa[] must be monotonically increasing and their 
//           range should not exceed period (ie xa[n-1] < xa[0] + period).  
//           x must lie in the range:  [(xa[n-1]-period), (xa[0]+period)]
//           "period" is typically 2*PI.
double cyc_splintD(double const *xa, 
                   double const *ya, 
                   double const *y2a, 
                   int n, 
                   double period,
                   double x)
{
  int klo = -1;
  int khi = n; // (not n-1)
  int k;
  double xlo = xa[n-1] - period;
  double xhi = xa[0] + period;
  assert((xlo <= x) && (x <= xhi));
  while (khi-klo > 1) {
    k = (khi+klo) >> 1; //(k=(khi+klo)/2)
    if (xa[k] > x) {
      khi = k;
      xhi = xa[k];
    }
    else {
      klo = k;
      xlo = xa[k];
    }
  }
  assert((xlo <= x) && (x <= xhi));

  if (khi == n) khi = 0;
  if (klo ==-1) klo = n-1;

  double yhi = ya[khi];
  double ylo = ya[klo];
  double h = xhi-xlo;
  double g = yhi-ylo;
  double a = (xhi-x) / h;
  double b = (x-xlo) / h;
  // Formula below taken from equation 3.3.5 of "numerical recipes in c"
  // "yD" = the derivative of y
  double yD = g/h - ( (3.0*a*a-1.0)*y2a[klo] - (3.0*b*b-1.0)*y2a[khi] ) * h/6.0;
  // For rerefence: y = a*ylo + b*yhi + 
  //                  ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;

  return yD;

} // cyc_splintD()



} //namespace DIHEDRAL_TABLE_NS
} //namespace LAMMPS_NS

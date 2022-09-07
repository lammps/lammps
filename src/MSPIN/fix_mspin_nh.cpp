/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Akhlak Mahmood

   Contact:
     Department of Materials Science and Engineering,
     North Carolina State University,
     Raleigh, NC, USA

     amahmoo3@ncsu.edu; mahmoodakhlak@gmail.com
------------------------------------------------------------------------- */

#include "fix_mspin_nh.h"

#include "atom.h"
#include "citeme.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace std;
using namespace LAMMPS_NS;

static const char cite_fix_mspin_c[] = "mspin package: doi:10.1021/acs.jctc.1c01253\n\n"
                                       "@Article{Mahmood22,\n"
                                       " author = {A.U. Mahmood and Y.G. Yingling},\n"
                                       " title = {All-Atom Simulation Method for Zeeman Alignment "
                                       "and Dipolar Assembly of Magnetic Nanoparticles},\n"
                                       " journal = {Journal of Chemical Theory and Computation},\n"
                                       " year =    2022,\n"
                                       " volume =  18(5),\n"
                                       " pages =   {3122-3135}\n"
                                       "}\n\n";

/* ---------------------------------------------------------------------- */

FixMspinNH::FixMspinNH(LAMMPS *lmp, int narg, char **arg) :
    FixRigidNH(lmp, narg, arg), mu(nullptr), dq(nullptr), qmcount(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_mspin_c);

  if (rstyle != 1)    // MOLECULE should be 1
    error->all(FLERR, "Fix mspin requires molecule style");

  // set default flag values
  // by default no interactions are on
  // the resulting dynamics should be exactly same as rigid/nh dynamics
  zeeman_flag = 0;
  dipolar_flag = 0;
  uniform_field = 0;

  // this fix contributes to the global energy
  energy_global_flag = 1;

  // disable langevin thermostating
  langflag = 0;

  // demagnetizing/magnetizing factor
  // Sánchez and Raap et. al. Physical Review B 2017, 95 (13), 134421
  alpha = 1.0;

  // a second scaling factor for the dipole moment itself
  beta = 1.0;

  // do not normalize the potential energy by natoms
  // 0/1 if global scalar is intensive/extensive
  extscalar = 0;

  // we need a unit multiplier for Tesla
  qb2f = force->qe2f * 1.0E-5;    // qe2f units * fs/e/A

  mu_0 = 4.6434E-4;    // force/Ampere^2 in real

  int iarg = 2;
  while (iarg < narg) {
    // external B field, Zeeman calculations ON
    if (strcmp(arg[iarg], "bfield") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix mspin bfield keyword");

      // Parse the jacobian components from the args
      // If it's an uniform field, we will treat these as xyz components
      bxdx = qb2f * utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      bydy = qb2f * utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      bzdz = qb2f * utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 3;

      zeeman_flag = 1;

      if (iarg + 1 < narg && strcmp(arg[iarg + 1], "uniform") == 0) {
        // optional 'uniform' specified?
        uniform_field = 1;
        iarg++;
      }
    }

    // turn on dipole dipole interaction
    if (strcmp(arg[iarg], "dpcut") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mspin dpcut keyword");
      dipolar_flag = 1;
      dipole_cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg++;
    }

    // dipolar scaling factor
    if (strcmp(arg[iarg], "alpha") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mspin alpha keyword");
      alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg++;
    }

    // zeeman+dipolar scaling factor
    if (strcmp(arg[iarg], "beta") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix mspin beta keyword");
      beta = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg++;
    }

    iarg++;
  }

  // alpha is only relevant for dipolar interaction
  // check if dipolar interaction is turned on if alpha is not the default
  if (!dipolar_flag && alpha != 1.0)
    error->all(FLERR,
               "Scaling factor alpha applies to dipolar interaction only. Did you forget the dpcut "
               "keyword?");

  // we assume the B field changes only in x or y or z directions
  // dBz/dx = 0, for example, and so on (i.e. a diagonal Jacobian).
  bxdy = bxdz = 0.0;
  bydx = bydz = 0.0;
  bzdx = bzdy = 0.0;

  memory->create(mu, nbody, 3, "mspin:mu");
  memory->create(dq, nbody, 3, "mspin:dq");
  memory->create(qm, nbody, "mspin:qm");
  memory->create(qmcount, nbody, "mspin:mspin_count");

  // initialize constants/counts
  nsum = 0;
  for (int ibody = 0; ibody < nbody; ibody++) nsum += nrigid[ibody];

  if (me == 0) {

    if (screen) fprintf(screen, "\nMSPIN initialization ...\n");
    if (logfile) fprintf(logfile, "\nMSPIN initialization ...\n");

    if (nbody == 1 && dipolar_flag == 1) {
      if (screen) fprintf(screen, "  only 1 rigid body found, turning off dipolar interaction\n");
      if (logfile) fprintf(logfile, "  only 1 rigid body found, turning off dipolar interaction\n");

      dipolar_flag = 0;
    }

    if (dipolar_flag == 1) {
      if (screen)
        fprintf(screen, "  implementing magnetic dipolar interactions with cutoff %f A\n",
                dipole_cutoff);
      if (logfile)
        fprintf(logfile, "  implementing magnetic dipolar interactions with cutoff %f A\n",
                dipole_cutoff);

      if (screen) fprintf(screen, "  dipolar interaction scaling factor alpha %f\n", alpha);
      if (logfile) fprintf(logfile, "  dipolar interaction scaling factor alpha %f\n", alpha);
    }

    if (zeeman_flag == 1) {
      if (uniform_field == 0) {
        if (screen)
          fprintf(screen, "  non-uniform external B field %lf %lf %lf applied\n", bxdx / qb2f,
                  bydy / qb2f, bzdz / qb2f);
        if (logfile)
          fprintf(logfile, "  non-uniform external B field %lf %lf %lf applied\n", bxdx / qb2f,
                  bydy / qb2f, bzdz / qb2f);
      } else {
        if (screen)
          fprintf(screen, "  uniform external B field %lf %lf %lf applied\n", bxdx / qb2f,
                  bydy / qb2f, bzdz / qb2f);
        if (logfile)
          fprintf(logfile, "  uniform external B field %lf %lf %lf applied\n", bxdx / qb2f,
                  bydy / qb2f, bzdz / qb2f);
      }
    }

    if (!zeeman_flag && !dipolar_flag) {
      if (screen)
        fprintf(
            screen,
            "  both zeeman and dipolar interaction is off, dynamics will be same as rigid/nh\n");
      if (logfile)
        fprintf(
            logfile,
            "  both zeeman and dipolar interaction is off, dynamics will be same as rigid/nh\n");
    }

    if (screen) fprintf(screen, "  dipole moment scaling factor beta %f\n", beta);
    if (logfile) fprintf(logfile, "  dipole moment scaling factor beta %f\n", beta);
  }
}

FixMspinNH::~FixMspinNH()
{
  // base destructors are called automatically
  memory->destroy(mu);
  memory->destroy(dq);
  memory->destroy(qm);
  memory->destroy(qmcount);
}

int FixMspinNH::setmask()
{
  int mask = 0;

  // Get the base masks
  mask = FixRigidNH::setmask();

  // Add additional force calculation procedures
  if (zeeman_flag || dipolar_flag) mask |= FixConst::POST_FORCE;
  return mask;
}

void FixMspinNH::init()
{
  // Call base
  FixRigidNH::init();
  calculate_dipoles(1);
}

// called immediately before the first timestep and
// after forces are computed (optional)
void FixMspinNH::setup(int vflag)
{
  // Call base
  FixRigidNH::setup(vflag);
  post_force(vflag);
}

void FixMspinNH::post_force(int vflag)
{
  // Call base
  FixRigidNH::post_force(vflag);

  // Late force computation is not possible
  // since we need to add additional magnetic forces
  if (!earlyflag) compute_forces_and_torques();

  calculate_dipoles(0);

  // compute and add the extra magnetic forces
  if (zeeman_flag == 1) compute_zeeman();
  if (dipolar_flag == 1) compute_dipolar();

  earlyflag = 1;
}

void FixMspinNH::calculate_dipoles(int initialize)
{
  int ibody;
  double *q = atom->qm;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  imageint *image = atom->image;
  double unwrap[3];

  // reset global arrays
  for (int i = 0; i < nbody; i++) {
    mu[i][0] = mu[i][1] = mu[i][2] = 0.0;
    dq[i][0] = dq[i][1] = dq[i][2] = 0.0;
    qm[i] = 0.0;
    if (initialize) qmcount[i] = 0;
  }

  // loop over all proc atoms
  // @todo: instead of looping over all the atoms at each step
  //        is there a way to keep track of the MS atoms?
  for (int i = 0; i < nlocal; i++) {
    // if atom is a part of a rigid body, get the body's id
    if (body[i] < 0)
      continue;
    else
      ibody = body[i];

    // get unwrapped coordinates
    // we are interested in measuring the distance only
    domain->unmap(x[i], image[i], unwrap);

    // atom has qm defined and +ve
    if (q[i] > 0) {
      qm[ibody] = q[i];    // we take the +ve one as qm value
      dq[ibody][0] += unwrap[0];
      dq[ibody][1] += unwrap[1];
      dq[ibody][2] += unwrap[2];
      if (initialize) qmcount[ibody] += 1;
    }
    if (q[i] < 0) {
      dq[ibody][0] -= unwrap[0];
      dq[ibody][1] -= unwrap[1];
      dq[ibody][2] -= unwrap[2];
      if (initialize) qmcount[ibody] += 1;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, dq[0], 3 * nbody, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &qm[0], nbody, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (initialize) MPI_Allreduce(MPI_IN_PLACE, &qmcount[0], nbody, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for (int i = 0; i < nbody; i++) {
    if (initialize && qmcount[i] != 2)
      error->all(FLERR, "Fix mspin requires exactly 2 non-zero qm atoms per particle, found: {}",
                 qmcount[i]);

    // qmag charge in kilo-e/fs units, multiply by 1000
    mu[i][0] = 1000 * beta * qm[i] * dq[i][0];
    mu[i][1] = 1000 * beta * qm[i] * dq[i][1];
    mu[i][2] = 1000 * beta * qm[i] * dq[i][2];

    // printf("Proc [%d]: qm[%d] %f\n", me, i, qm[i]);

    if (me == 0 && initialize) {
      double d = sqrt(dq[i][0] * dq[i][0] + dq[i][1] * dq[i][1] + dq[i][2] * dq[i][2]);
      // print as kilo-e units
      double m = sqrt(mu[i][0] * mu[i][0] + mu[i][1] * mu[i][1] + mu[i][2] * mu[i][2]) / 1000.0;

      if (i == 0 && screen)   fprintf(screen, "\nMSPIN magnetic dipoles ...\n");
      if (i == 0 && logfile)  fprintf(logfile, "\nMSPIN magnetic dipoles ...\n");

      if (screen)
        fprintf(screen, "  Dipole %d: qm = %lf Ke/fs A\td = %lf A\tmu = %lf Ke/fs A^2\n", i + 1,
                qm[i], d, m);
      if (logfile)
        fprintf(logfile, "  Dipole %d: qm = %lf Ke/fs A\td = %lf A\tmu = %lf Ke/fs A^2\n", i + 1,
                qm[i], d, m);
      if (screen)
        fprintf(screen, "  Effective dipole moment [%d]:\tbeta * mu = %lf x 10^-21 A m^2\n", i + 1,
                m * beta * 1.6022);
      if (logfile)
        fprintf(logfile, "  Effective dipole moment [%d]:\tbeta * mu = %lf x 10^-21 A m^2\n", i + 1,
                m * beta * 1.6022);
    }
  }
}

void FixMspinNH::compute_zeeman()
{
  double fx, fy, fz;
  double tx, ty, tz;

  // printf("Proc %d: %s, %d \n", me, FLERR);

  zeeman_pe = 0.0;

  for (int ibody = 0; ibody < nbody; ibody++) {
    // External field will always align the particles.
    // T = mu cross B
    tx = bzdz * mu[ibody][1] - bydy * mu[ibody][2];
    ty = bxdx * mu[ibody][2] - bzdz * mu[ibody][0];
    tz = bydy * mu[ibody][0] - bxdx * mu[ibody][1];

    torque[ibody][0] += tx;
    torque[ibody][1] += ty;
    torque[ibody][2] += tz;

    // printf("Proc [%d]: dT[%d] %f, %f, %f \n", me,
    //             ibody, tx, ty, tz);

    // Non uniform field will also create a force.
    if (uniform_field == 0) {
      // F = mu dot grad B
      fx = mu[ibody][0] * bxdx + mu[ibody][1] * bxdy + mu[ibody][2] * bxdz;
      fy = mu[ibody][0] * bydx + mu[ibody][1] * bydy + mu[ibody][2] * bydz;
      fz = mu[ibody][0] * bzdx + mu[ibody][1] * bzdy + mu[ibody][2] * bzdz;

      fcm[ibody][0] += fx;
      fcm[ibody][1] += fy;
      fcm[ibody][2] += fz;

      // printf("Proc [%d]: dF[%d] %f, %f, %f \n", me,
      //             ibody, fx, fy, fz);
    }

    zeeman_pe -= (mu[ibody][0] * bxdx + mu[ibody][1] * bydy + mu[ibody][2] * bzdz);
  }
}

void FixMspinNH::compute_dipolar()
{
  double delx, dely, delz, fx, fy, fz;
  double rsq, r2inv, rinv, r3inv, r5inv, r7inv;
  double forcecoulx, forcecouly, forcecoulz, crossx, crossy, crossz;
  double tixcoul, tiycoul, tizcoul, tjxcoul, tjycoul, tjzcoul;
  double fq, pdotp, pidotr, pjdotr, pre1, pre2, pre3, pre4;

  int natoms;

  dipolar_pe = 0.0;

  // scaling and prefactors
  fq = alpha * mu_0 / MathConst::MY_4PI;

  // since no of rigid bodies, nbody is usually small
  // we calculate manybody interaction seperately in each proc
  // which scales as nC2, MPI communication can take longer time than this.
  // @todo: there has to be a better way to do it, but it works for now. :)
  for (int ibody = 0; ibody < nbody - 1; ibody++) {
    double iwrap[3];
    domain->unmap(xcm[ibody], imagebody[ibody], iwrap);

    for (int jbody = ibody + 1; jbody < nbody; jbody++) {
      // unwrap the rigid COM
      double jwrap[3];
      domain->unmap(xcm[jbody], imagebody[jbody], jwrap);

      // cm to cm distance
      delx = iwrap[0] - jwrap[0];
      dely = iwrap[1] - jwrap[1];
      delz = iwrap[2] - jwrap[2];

      rsq = delx * delx + dely * dely + delz * delz;

      // @todo: implement modified potential version if cutoff is used
      if (rsq < dipole_cutoff * dipole_cutoff) {
        r2inv = 1.0 / rsq;
        rinv = sqrt(r2inv);
        r3inv = r2inv * rinv;
        r5inv = r3inv * r2inv;
        r7inv = r5inv * r2inv;

        // force dot terms
        pdotp =
            mu[ibody][0] * mu[jbody][0] + mu[ibody][1] * mu[jbody][1] + mu[ibody][2] * mu[jbody][2];
        pidotr = mu[ibody][0] * delx + mu[ibody][1] * dely + mu[ibody][2] * delz;
        pjdotr = mu[jbody][0] * delx + mu[jbody][1] * dely + mu[jbody][2] * delz;

        pre1 = 3.0 * r5inv * pdotp - 15.0 * r7inv * pidotr * pjdotr;
        pre2 = 3.0 * r5inv * pjdotr;
        pre3 = 3.0 * r5inv * pidotr;
        pre4 = -1.0 * r3inv;

        forcecoulx = pre1 * delx + pre2 * mu[ibody][0] + pre3 * mu[jbody][0];
        forcecouly = pre1 * dely + pre2 * mu[ibody][1] + pre3 * mu[jbody][1];
        forcecoulz = pre1 * delz + pre2 * mu[ibody][2] + pre3 * mu[jbody][2];

        crossx = pre4 * (mu[ibody][1] * mu[jbody][2] - mu[ibody][2] * mu[jbody][1]);
        crossy = pre4 * (mu[ibody][2] * mu[jbody][0] - mu[ibody][0] * mu[jbody][2]);
        crossz = pre4 * (mu[ibody][0] * mu[jbody][1] - mu[ibody][1] * mu[jbody][0]);

        tixcoul = crossx + pre2 * (mu[ibody][1] * delz - mu[ibody][2] * dely);
        tiycoul = crossy + pre2 * (mu[ibody][2] * delx - mu[ibody][0] * delz);
        tizcoul = crossz + pre2 * (mu[ibody][0] * dely - mu[ibody][1] * delx);
        tjxcoul = -crossx + pre3 * (mu[jbody][1] * delz - mu[jbody][2] * dely);
        tjycoul = -crossy + pre3 * (mu[jbody][2] * delx - mu[jbody][0] * delz);
        tjzcoul = -crossz + pre3 * (mu[jbody][0] * dely - mu[jbody][1] * delx);

        fx = fq * forcecoulx;
        fy = fq * forcecouly;
        fz = fq * forcecoulz;

        // add dipolar force on the cm
        fcm[ibody][0] += fx;
        fcm[ibody][1] += fy;
        fcm[ibody][2] += fz;
        // Newton's 3rd law
        fcm[jbody][0] -= fx;
        fcm[jbody][1] -= fy;
        fcm[jbody][2] -= fz;

        // add dipolar torque on the cm
        torque[ibody][0] += fq * tixcoul;
        torque[ibody][1] += fq * tiycoul;
        torque[ibody][2] += fq * tizcoul;
        // Newton's 3rd law
        torque[jbody][0] += fq * tjxcoul;
        torque[jbody][1] += fq * tjycoul;
        torque[jbody][2] += fq * tjzcoul;

        // printf("Proc [%d]: dF[%d][%d] %f, %f, %f \n", me,
        //             ibody, jbody, fx, fy, fz);

        // calculate dipolar interaction energy
        dipolar_pe += fq * (r3inv * pdotp - 3.0 * r5inv * pidotr * pjdotr);
      }
    }
  }
}

// return total potential energy
double FixMspinNH::compute_scalar()
{
  // Call base
  double energy = FixRigidNH::compute_scalar();

  // additional energy due to magnetic interactions
  energy += extract_zeeman_pe();
  energy += extract_dipolar_pe();

  return energy;
}

// Return Zeeman PE to output as thermo
double FixMspinNH::extract_zeeman_pe()
{
  return zeeman_flag * zeeman_pe;
}

// Return dipolar PE to output as thermo
double FixMspinNH::extract_dipolar_pe()
{
  return dipolar_flag * dipolar_pe;
}

// COM distance between two rigid bodies
double FixMspinNH::extract_distance(int ibody, int jbody)
{
  double delx, dely, delz, d;

  double iwrap[3];
  double jwrap[3];

  if (ibody > nbody or jbody > nbody)
    error->all(FLERR, "Invalid molecule id for mspin/distance compute");
  if (ibody == 0 or jbody == 0) error->all(FLERR, "Invalid molecule id for mspin/distance compute");

  domain->unmap(xcm[ibody - 1], imagebody[ibody - 1], iwrap);
  domain->unmap(xcm[jbody - 1], imagebody[jbody - 1], jwrap);

  // unwrapped cm to cm distance
  delx = iwrap[0] - jwrap[0];
  dely = iwrap[1] - jwrap[1];
  delz = iwrap[2] - jwrap[2];

  return sqrt(delx * delx + dely * dely + delz * delz);
}

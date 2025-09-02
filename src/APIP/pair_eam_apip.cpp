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
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL) (EAM)
     David Immel (d.immel@fz-juelich.de, FZJ, Germany) for APIP
------------------------------------------------------------------------- */

#include "pair_eam_apip.h"

#include "atom.h"
#include "atom_vec_apip.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEAMAPIP::PairEAMAPIP(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  nmax = 0;
  rho = nullptr;
  fp = nullptr;
  numforce = nullptr;
  type2frho = nullptr;

  nfuncfl = 0;
  funcfl = nullptr;

  setfl = nullptr;
  fs = nullptr;

  frho = nullptr;
  rhor = nullptr;
  z2r = nullptr;
  scale = nullptr;

  rhomax = rhomin = 0.0;

  frho_spline = nullptr;
  rhor_spline = nullptr;
  z2r_spline = nullptr;

  n_non_complex_accumulated = 0;
  time_per_atom = -1;
  time_wall_accumulated = 0;

  lambda_thermostat = true;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEAMAPIP::~PairEAMAPIP()
{
  if (copymode) return;

  memory->destroy(rho);
  memory->destroy(fp);
  memory->destroy(numforce);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete[] type2frho;
    type2frho = nullptr;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);
    memory->destroy(scale);
  }

  if (funcfl) {
    for (int i = 0; i < nfuncfl; i++) {
      delete[] funcfl[i].file;
      memory->destroy(funcfl[i].frho);
      memory->destroy(funcfl[i].rhor);
      memory->destroy(funcfl[i].zr);
    }
    memory->sfree(funcfl);
    funcfl = nullptr;
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete[] setfl->elements[i];
    delete[] setfl->elements;
    memory->destroy(setfl->mass);
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    delete setfl;
    setfl = nullptr;
  }

  if (fs) {
    for (int i = 0; i < fs->nelements; i++) delete[] fs->elements[i];
    delete[] fs->elements;
    memory->destroy(fs->mass);
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
    fs = nullptr;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
}

/* ---------------------------------------------------------------------- */

void PairEAMAPIP::compute(int eflag, int vflag)
{
  // start timers
  double time_wall_start = platform::walltime();
  int n_non_complex = 0;

  int i, j, ii, jj, m, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair, fpair_cl;
  double rsq, r, p, rhoip, rhojp, z2, z2p, recip, phip, psip, phi, psip_cl;
  double *coeff;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    memory->destroy(numforce);
    nmax = atom->nmax;
    memory->create(rho, nmax, "pair:rho");
    memory->create(fp, nmax, "pair:fp");
    memory->create(numforce, nmax, "pair:numforce");
  }

  double **x = atom->x;
  double **f = atom->f;
  double *lambda = atom->apip_lambda;
  int *lambda_required = atom->apip_lambda_required;

  double **f_const_lambda = nullptr;
  double **f_dyn_lambda = nullptr;
  double *e_simple = nullptr;
  double *lambda_const = nullptr;
  if (lambda_thermostat) {
    f_const_lambda = atom->apip_f_const_lambda;
    f_dyn_lambda = atom->apip_f_dyn_lambda;
    e_simple = atom->apip_e_fast;
    lambda_const = atom->apip_lambda_const;
  }
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // avoid double counting due to full neighbour list for two local particles
      // j < i implies j < nlocal
      if (j < i && i < nlocal) { continue; }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq) * rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, nr - 1);
        p -= m;
        p = MIN(p, 1.0);
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

        // do not calculate rho for ghost atoms
        if (j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rho[j] += ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
        }
      }
    }
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i] * rdrho + 1.0;
    m = static_cast<int>(p);
    m = MAX(1, MIN(m, nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0] * p + coeff[1]) * p + coeff[2];
    if (eflag || e_simple) {
      phi = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
      if (rho[i] > rhomax) phi += fp[i] * (rho[i] - rhomax);
      phi *= scale[type[i]][type[i]];
      ev_tally_full(i, 2.0 * lambda[i] * phi, 0.0, 0.0, 0.0, 0.0, 0.0);
      if (e_simple) { e_simple[i] = phi; }
    }
  }

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    numforce[i] = 0;

    // The distances between atoms are not calculated.
    // The neighbour list contains more than the cutoff atoms.
    // The calculation of distances is probably not worth the compute time.
    // lambda_required is used to compute the weight of atoms for load balancing.
    // -> store information about required calculations
    if (lambda_required[i] & ApipLambdaRequired::NO_SIMPLE)
      continue;
    else if (!(lambda_required[i] & ApipLambdaRequired::SIMPLE)) {
      // neither SIMPLE nor NO_SIMPLE set

      // check own atom
      if (lambda[i] != 0 || (lambda_thermostat && lambda_const[i] != 0)) {
        // set own atom
        lambda_required[i] |= ApipLambdaRequired::SIMPLE;
        // set neighbour list
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          if (j < nlocal) lambda_required[j] |= ApipLambdaRequired::SIMPLE;
        }
      } else {
        // check neighbour list
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          if (lambda[j] != 0 || (lambda_thermostat && lambda_const[j] != 0)) {
            lambda_required[i] |= ApipLambdaRequired::SIMPLE;
            // set lambda also for non-ghost j
            if (j < nlocal) lambda_required[j] |= ApipLambdaRequired::SIMPLE;
            break;
          }
        }
      }
      // SIMPLE not set -> set NO_SIMPLE
      if (!(lambda_required[i] & ApipLambdaRequired::SIMPLE)) {
        // go to next atom
        lambda_required[i] |= ApipLambdaRequired::NO_SIMPLE;
        continue;
      }
    }
    n_non_complex++;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // avoid double counting due to full neighbour list for two local particles
      // j < i implies j < nlocal
      if (j < i && i < nlocal) { continue; }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutforcesq) {
        ++numforce[i];
        jtype = type[j];
        r = sqrt(rsq);
        p = r * rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, nr - 1);
        p -= m;
        p = MIN(p, 1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        // scale factor can be applied by thermodynamic integration

        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0] * p + coeff[1]) * p + coeff[2];
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
        z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

        recip = 1.0 / r;
        phi = z2 * recip;
        phip = z2p * recip - phi * recip;
        if (j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rhoip = (coeff[0] * p + coeff[1]) * p + coeff[2];
          psip = lambda[i] * fp[i] * rhojp + lambda[j] * fp[j] * rhoip +
              phip * (lambda[i] + lambda[j]) / 2;
        } else {
          // processor of j calculates the remaining terms (compared to psip in if case)
          psip = lambda[i] * fp[i] * rhojp + phip * 0.5 * lambda[i];
        }

        fpair = -scale[itype][jtype] * psip * recip;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        f[j][0] -= delx * fpair;
        f[j][1] -= dely * fpair;
        f[j][2] -= delz * fpair;
        if (lambda_thermostat) {

          f_dyn_lambda[i][0] += delx * fpair;
          f_dyn_lambda[i][1] += dely * fpair;
          f_dyn_lambda[i][2] += delz * fpair;
          f_dyn_lambda[j][0] -= delx * fpair;
          f_dyn_lambda[j][1] -= dely * fpair;
          f_dyn_lambda[j][2] -= delz * fpair;

          // psip_const
          if (j < nlocal) {
            psip_cl = lambda_const[i] * fp[i] * rhojp + lambda_const[j] * fp[j] * rhoip +
                phip * (lambda_const[i] + lambda_const[j]) / 2;
          } else {
            psip_cl = lambda_const[i] * fp[i] * rhojp + phip * 0.5 * lambda_const[i];
          }
          // calculate fpair_const
          fpair_cl = -scale[itype][jtype] * psip_cl * recip;
          // update f_const_lambda with fpair_const
          f_const_lambda[i][0] += delx * fpair_cl;
          f_const_lambda[i][1] += dely * fpair_cl;
          f_const_lambda[i][2] += delz * fpair_cl;
          f_const_lambda[j][0] -= delx * fpair_cl;
          f_const_lambda[j][1] -= dely * fpair_cl;
          f_const_lambda[j][2] -= delz * fpair_cl;
        }

        if (eflag || e_simple) {
          evdwl = scale[itype][jtype] * phi;
          if (e_simple) {
            e_simple[i] += 0.5 * evdwl;
            if (j < nlocal) e_simple[j] += 0.5 * evdwl;
          }
          ev_tally_full(i, lambda[i] * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
          if (j < nlocal) ev_tally_full(j, lambda[j] * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        if (vflag) ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  // stop timers
  time_wall_accumulated += platform::walltime() - time_wall_start;
  n_non_complex_accumulated += n_non_complex;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEAMAPIP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  delete[] map;
  map = new int[n + 1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n + 1];
  memory->create(type2rhor, n + 1, n + 1, "pair:type2rhor");
  memory->create(type2z2r, n + 1, n + 1, "pair:type2z2r");
  memory->create(scale, n + 1, n + 1, "pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEAMAPIP::settings(int narg, char ** /*arg*/)
{
  if (narg > 0) error->all(FLERR, "Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void PairEAMAPIP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR, "Incorrect args for pair coefficients");

  // parse pair of atom types

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  // read funcfl file if hasn't already been read
  // store filename in Funcfl data struct

  int ifuncfl;
  for (ifuncfl = 0; ifuncfl < nfuncfl; ifuncfl++)
    if (strcmp(arg[2], funcfl[ifuncfl].file) == 0) break;

  if (ifuncfl == nfuncfl) {
    nfuncfl++;
    funcfl = (Funcfl *) memory->srealloc(funcfl, nfuncfl * sizeof(Funcfl), "pair:funcfl");
    read_file(arg[2]);
    funcfl[ifuncfl].file = utils::strdup(arg[2]);
  }

  // set setflag and map only for i,i type pairs
  // set mass of atom type if i = j

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      if (i == j) {
        setflag[i][i] = 1;
        map[i] = ifuncfl;
        atom->set_mass(FLERR, i, funcfl[ifuncfl].mass);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMAPIP::init_style()
{
  // convert read-in file(s) to arrays and spline them
  if (!atom->apip_lambda_flag)
    error->all(FLERR, "Pair style eam/apip requires an atom style with lambda");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style eam/apip requires newton pair on");
  if (!atom->apip_lambda_required_flag)
    error->all(FLERR, "pair style eam/apip requires an atom style with lambda_required.");

  file2array();
  array2spline();

  // communication during computation should be avoided
  // -> do not exchange the derivative of the embedding function by default
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEAMAPIP::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (setflag[i][j] == 0) scale[i][j] = 1.0;
  scale[j][i] = scale[i][j];

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++) cutmax = MAX(cutmax, funcfl[m].cut);
  } else if (setfl)
    cutmax = setfl->cut;
  else if (fs)
    cutmax = fs->cut;

  cutforcesq = cutmax * cutmax;

  return cutmax;
}

/**
  * setup specific to this pair style
  * Determine whether there is a fix lambda_thermostat/apip or not and set
  * lambda_thermostat.
  */

void PairEAMAPIP::setup()
{
  if (modify->get_fix_by_style("^lambda_thermostat/apip$").size() == 0) {
    lambda_thermostat = false;
  } else {
    lambda_thermostat = true;
    if (!atom->apip_lambda_const_flag)
      error->all(
          FLERR,
          "Pair style pace/apip requires an atom style with lambda_const for a local thermostat.");
    if (!atom->apip_e_fast_flag)
      error->all(
          FLERR,
          "Pair style pace/apip requires an atom style with e_simple for a local thermostat.");
    if (!atom->apip_f_const_lambda_flag)
      error->all(FLERR,
                 "Pair style pace/apip requires an atom style with f_const_lambda for a local "
                 "thermostat.");
    if (!atom->apip_f_dyn_lambda_flag)
      error->all(FLERR,
                 "Pair style pace/apip requires an atom style with f_const_lambda for a local "
                 "thermostat.");
  }
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void PairEAMAPIP::read_file(char *filename)
{
  Funcfl *file = &funcfl[nfuncfl - 1];

  // read potential file
  if (comm->me == 0) {
    PotentialFileReader reader(lmp, filename, "eam", unit_convert_flag);

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY, unit_convert);
    try {
      reader.skip_line();

      ValueTokenizer values = reader.next_values(2);
      values.next_int();    // ignore
      file->mass = values.next_double();

      values = reader.next_values(5);
      file->nrho = values.next_int();
      file->drho = values.next_double();
      file->nr = values.next_int();
      file->dr = values.next_double();
      file->cut = values.next_double();

      if ((file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
        error->one(FLERR, "Invalid EAM potential file");

      memory->create(file->frho, (file->nrho + 1), "pair:frho");
      memory->create(file->rhor, (file->nr + 1), "pair:rhor");
      memory->create(file->zr, (file->nr + 1), "pair:zr");

      reader.next_dvector(&file->frho[1], file->nrho);
      reader.next_dvector(&file->zr[1], file->nr);
      reader.next_dvector(&file->rhor[1], file->nr);

      if (unit_convert) {
        const double sqrt_conv = sqrt(conversion_factor);
        for (int i = 1; i <= file->nrho; ++i) file->frho[i] *= conversion_factor;
        for (int j = 1; j <= file->nr; ++j) file->zr[j] *= sqrt_conv;
      }
    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  MPI_Bcast(&file->mass, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nrho, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->drho, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->nr, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->dr, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->cut, 1, MPI_DOUBLE, 0, world);

  if (comm->me != 0) {
    memory->create(file->frho, (file->nrho + 1), "pair:frho");
    memory->create(file->rhor, (file->nr + 1), "pair:rhor");
    memory->create(file->zr, (file->nr + 1), "pair:zr");
  }

  MPI_Bcast(&file->frho[1], file->nrho, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->zr[1], file->nr, MPI_DOUBLE, 0, world);
  MPI_Bcast(&file->rhor[1], file->nr, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairEAMAPIP::file2array()
{
  int i, j, k, m, n;
  int ntypes = atom->ntypes;
  double sixth = 1.0 / 6.0;

  // determine max function params from all active funcfl files
  // active means some element is pointing at it via map

  int active;
  double rmax;
  dr = drho = rmax = rhomax = 0.0;

  for (int i = 0; i < nfuncfl; i++) {
    active = 0;
    for (j = 1; j <= ntypes; j++)
      if (map[j] == i) active = 1;
    if (active == 0) continue;
    Funcfl *file = &funcfl[i];
    dr = MAX(dr, file->dr);
    drho = MAX(drho, file->drho);
    rmax = MAX(rmax, (file->nr - 1) * file->dr);
    rhomax = MAX(rhomax, (file->nrho - 1) * file->drho);
  }

  // set nr,nrho from cutoff and spacings

  nr = std::lround(rmax / dr);
  nrho = std::lround(rhomax / drho);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of funcfl files + 1 for zero array

  nfrho = nfuncfl + 1;
  memory->destroy(frho);
  memory->create(frho, nfrho, nrho + 1, "pair:frho");

  // interpolate each file's frho to a single grid and cutoff

  double r, p, cof1, cof2, cof3, cof4;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nrho; m++) {
      r = (m - 1) * drho;
      p = r / file->drho + 1.0;
      k = static_cast<int>(p);
      k = MIN(k, file->nrho - 2);
      k = MAX(k, 2);
      p -= k;
      p = MIN(p, 2.0);
      cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
      cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
      cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
      cof4 = sixth * p * (p * p - 1.0);
      frho[n][m] = cof1 * file->frho[k - 1] + cof2 * file->frho[k] + cof3 * file->frho[k + 1] +
          cof4 * file->frho[k + 2];
    }
    n++;
  }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho - 1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to file (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0)
      type2frho[i] = map[i];
    else
      type2frho[i] = nfrho - 1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of funcfl files

  nrhor = nfuncfl;
  memory->destroy(rhor);
  memory->create(rhor, nrhor, nr + 1, "pair:rhor");

  // interpolate each file's rhor to a single grid and cutoff

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *file = &funcfl[i];
    for (m = 1; m <= nr; m++) {
      r = (m - 1) * dr;
      p = r / file->dr + 1.0;
      k = static_cast<int>(p);
      k = MIN(k, file->nr - 2);
      k = MAX(k, 2);
      p -= k;
      p = MIN(p, 2.0);
      cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
      cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
      cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
      cof4 = sixth * p * (p * p - 1.0);
      rhor[n][m] = cof1 * file->rhor[k - 1] + cof2 * file->rhor[k] + cof3 * file->rhor[k + 1] +
          cof4 * file->rhor[k + 2];
    }
    n++;
  }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for funcfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++) type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of funcfl files

  nz2r = nfuncfl * (nfuncfl + 1) / 2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r, nr + 1, "pair:z2r");

  // create a z2r array for each file against other files, only for I >= J
  // interpolate zri and zrj to a single grid and cutoff
  // final z2r includes unit conversion of 27.2 eV/Hartree and 0.529 Ang/Bohr

  double zri, zrj;

  n = 0;
  for (i = 0; i < nfuncfl; i++) {
    Funcfl *ifile = &funcfl[i];
    for (j = 0; j <= i; j++) {
      Funcfl *jfile = &funcfl[j];

      for (m = 1; m <= nr; m++) {
        r = (m - 1) * dr;

        p = r / ifile->dr + 1.0;
        k = static_cast<int>(p);
        k = MIN(k, ifile->nr - 2);
        k = MAX(k, 2);
        p -= k;
        p = MIN(p, 2.0);
        cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
        cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
        cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
        cof4 = sixth * p * (p * p - 1.0);
        zri = cof1 * ifile->zr[k - 1] + cof2 * ifile->zr[k] + cof3 * ifile->zr[k + 1] +
            cof4 * ifile->zr[k + 2];

        p = r / jfile->dr + 1.0;
        k = static_cast<int>(p);
        k = MIN(k, jfile->nr - 2);
        k = MAX(k, 2);
        p -= k;
        p = MIN(p, 2.0);
        cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
        cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
        cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
        cof4 = sixth * p * (p * p - 1.0);
        zrj = cof1 * jfile->zr[k - 1] + cof2 * jfile->zr[k] + cof3 * jfile->zr[k + 1] +
            cof4 * jfile->zr[k + 2];

        z2r[n][m] = 27.2 * 0.529 * zri * zrj;
      }
      n++;
    }
  }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairEAMAPIP::array2spline()
{
  rdr = 1.0 / dr;
  rdrho = 1.0 / drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline, nfrho, nrho + 1, 7, "pair:frho");
  memory->create(rhor_spline, nrhor, nr + 1, 7, "pair:rhor");
  memory->create(z2r_spline, nz2r, nr + 1, 7, "pair:z2r");

  for (int i = 0; i < nfrho; i++) interpolate(nrho, drho, frho[i], frho_spline[i]);

  for (int i = 0; i < nrhor; i++) interpolate(nr, dr, rhor[i], rhor_spline[i]);

  for (int i = 0; i < nz2r; i++) interpolate(nr, dr, z2r[i], z2r_spline[i]);
}

/* ---------------------------------------------------------------------- */

void PairEAMAPIP::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5] = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
    spline[m][5] =
        ((spline[m - 2][6] - spline[m + 2][6]) + 8.0 * (spline[m + 1][6] - spline[m - 1][6])) /
        12.0;

  for (int m = 1; m <= n - 1; m++) {
    spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) - 2.0 * spline[m][5] - spline[m + 1][5];
    spline[m][3] = spline[m][5] + spline[m + 1][5] - 2.0 * (spline[m + 1][6] - spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5] / delta;
    spline[m][1] = 2.0 * spline[m][4] / delta;
    spline[m][0] = 3.0 * spline[m][3] / delta;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairEAMAPIP::memory_usage()
{
  double bytes = (double) maxeatom * sizeof(double);
  bytes += (double) maxvatom * 6 * sizeof(double);
  bytes += (double) 2 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   swap fp array with one passed in by caller
------------------------------------------------------------------------- */

void PairEAMAPIP::swap_eam(double *fp_caller, double **fp_caller_hold)
{
  double *tmp = fp;
  fp = fp_caller;
  *fp_caller_hold = tmp;
}

/* ----------------------------------------------------------------------
   set return values for timers and counted particles
------------------------------------------------------------------------- */

void PairEAMAPIP::calculate_time_per_atom()
{
  if (n_non_complex_accumulated > 0)
    time_per_atom = time_wall_accumulated / n_non_complex_accumulated;
  else
    time_per_atom = -1;

  // reset
  time_wall_accumulated = 0;
  n_non_complex_accumulated = 0;
}

/* ---------------------------------------------------------------------- */

void *PairEAMAPIP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;
  dim = 0;
  if (strcmp(str, "eam/apip:time_per_atom") == 0) {
    calculate_time_per_atom();
    return (void *) &time_per_atom;
  }
  return nullptr;
}

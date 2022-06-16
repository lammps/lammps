// clang-format off
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
   Contributing author: Alexander Stukowski
                        Technical University of Darmstadt,
                        Germany Department of Materials Science
------------------------------------------------------------------------- */

#include "pair_eam_cd.h"

#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "tokenizer.h"

#include <cmath>

using namespace LAMMPS_NS;

#define MAXLINE 1024        // This sets the maximum line length in EAM input files.

PairEAMCD::PairEAMCD(LAMMPS *lmp, int _cdeamVersion)
  : PairEAM(lmp), PairEAMAlloy(lmp), cdeamVersion(_cdeamVersion)
{
  single_enable = 0;
  restartinfo = 0;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  rhoB = nullptr;
  D_values = nullptr;
  hcoeff = nullptr;

  // Set communication buffer sizes needed by this pair style.

  if (cdeamVersion == 1) {
    comm_forward = 4;
    comm_reverse = 3;
  } else if (cdeamVersion == 2) {
    comm_forward = 3;
    comm_reverse = 2;
  } else {
    error->all(FLERR,"Invalid eam/cd potential version.");
  }
}

PairEAMCD::~PairEAMCD()
{
  memory->destroy(rhoB);
  memory->destroy(D_values);
  delete[] hcoeff;
}

void PairEAMCD::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rhoip,rhojp,recip,phi;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  // Grow per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    memory->destroy(rhoB);
    memory->destroy(D_values);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(rhoB,nmax,"pair:rhoB");
    memory->create(fp,nmax,"pair:fp");
    memory->create(D_values,nmax,"pair:D_values");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Zero out per-atom arrays.

  int m = nlocal + atom->nghost;
  for (i = 0; i < m; i++) {
    rho[i] = 0.0;
    rhoB[i] = 0.0;
    D_values[i] = 0.0;
  }

  // Stage I

  // Compute rho and rhoB at each local atom site.

  // Additionally calculate the D_i values here if we are using the
  // one-site formulation.  For the two-site formulation we have to
  // calculate the D values in an extra loop (Stage II).

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

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        double r = sqrt(rsq);
        const EAMTableIndex index = radiusToTableIndex(r);
        double localrho = RhoOfR(index, jtype, itype);
        rho[i] += localrho;
        if (jtype == speciesB) rhoB[i] += localrho;
        if (newton_pair || j < nlocal) {
          localrho = RhoOfR(index, itype, jtype);
          rho[j] += localrho;
          if (itype == speciesB) rhoB[j] += localrho;
        }

        if (cdeamVersion == 1 && itype != jtype) {

          // Note: if the i-j interaction is not concentration dependent (because either
          // i or j are not species A or B) then its contribution to D_i and D_j should
          // be ignored.
          // This if-clause is only required for a ternary.

          if ((itype == speciesA && jtype == speciesB)
              || (jtype == speciesA && itype == speciesB)) {
            double Phi_AB = PhiOfR(index, itype, jtype, 1.0 / r);
            D_values[i] += Phi_AB;
            if (newton_pair || j < nlocal)
              D_values[j] += Phi_AB;
          }
        }
      }
    }
  }

  // Communicate and sum densities.

  if (newton_pair) {
    communicationStage = 1;
    comm->reverse_comm(this);
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    EAMTableIndex index = rhoToTableIndex(rho[i]);
    fp[i] = FPrimeOfRho(index, type[i]);
    if (eflag) {
      phi = FofRho(index, type[i]);
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // Communicate derivative of embedding function and densities
  // and D_values (this for one-site formulation only).

  communicationStage = 2;
  comm->forward_comm(this);

  // The electron densities may not drop to zero because then the
  // concentration would no longer be defined.  But the concentration
  // is not needed anyway if there is no interaction with another atom,
  // which is the case if the electron density is exactly zero.
  // That's why the following lines have been commented out.
  //
  //for (i = 0; i < nlocal + atom->nghost; i++) {
  //        if (rho[i] == 0 && (type[i] == speciesA || type[i] == speciesB))
  //                error->one(FLERR,"CD-EAM potential routine: Detected atom with zero electron density.");
  //}

  // Stage II
  // This is only required for the original two-site formulation of the CD-EAM potential.

  if (cdeamVersion == 2) {

    // Compute intermediate value D_i for each atom.

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // This code line is required for ternary alloys.

      if (itype != speciesA && itype != speciesB) continue;

      double x_i = rhoB[i] / rho[i];        // Concentration at atom i.

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];
        if (itype == jtype) continue;

        // This code line is required for ternary alloys.

        if (jtype != speciesA && jtype != speciesB) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < cutforcesq) {
          double r = sqrt(rsq);
          const EAMTableIndex index = radiusToTableIndex(r);

          // The concentration independent part of the cross pair potential.

          double Phi_AB = PhiOfR(index, itype, jtype, 1.0 / r);

          // Average concentration of two sites

          double x_ij = 0.5 * (x_i + rhoB[j]/rho[j]);

          // Calculate derivative of h(x_ij) polynomial function.

          double h_prime = evalHprime(x_ij);

          D_values[i] += h_prime * Phi_AB / (2.0 * rho[i] * rho[i]);
          if (newton_pair || j < nlocal)
            D_values[j] += h_prime * Phi_AB / (2.0 * rho[j] * rho[j]);
        }
      }
    }

    // Communicate and sum D values.

    if (newton_pair) {
      communicationStage = 3;
      comm->reverse_comm(this);
    }
    communicationStage = 4;
    comm->forward_comm(this);
  }

  // Stage III

  // Compute force acting on each atom.

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Concentration at site i
    // The value -1 indicates: no concentration dependence for all interactions of atom i.
    // It will be replaced by the concentration at site i if atom i is either A or B.

    double x_i = -1.0;
    double D_i = 0.0, h_prime_i;

    // This if-clause is only required for ternary alloys.

    if ((itype == speciesA || itype == speciesB) && rho[i] != 0.0) {

      // Compute local concentration at site i.

      x_i = rhoB[i]/rho[i];

      if (cdeamVersion == 1) {

        // Calculate derivative of h(x_i) polynomial function.

        h_prime_i = evalHprime(x_i);
        D_i = D_values[i] * h_prime_i / (2.0 * rho[i] * rho[i]);
      } else if (cdeamVersion == 2) {
        D_i = D_values[i];
      }
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        double r = sqrt(rsq);
        const EAMTableIndex index = radiusToTableIndex(r);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

        rhoip = RhoPrimeOfR(index, itype, jtype);
        rhojp = RhoPrimeOfR(index, jtype, itype);
        fpair = fp[i]*rhojp + fp[j]*rhoip;
        recip = 1.0/r;

        // The value -1 indicates: no concentration dependence for this
        // i-j pair because atom j is not of species A nor B.

        double x_j = -1;

        // This code line is required for ternary alloy.

        if ((jtype == speciesA || jtype == speciesB) && rho[j] != 0.0) {

          // Compute local concentration at site j.

          x_j = rhoB[j]/rho[j];

          double D_j=0.0;
          if (cdeamVersion == 1) {

            // Calculate derivative of h(x_j) polynomial function.

            double h_prime_j = evalHprime(x_j);
            D_j = D_values[j] * h_prime_j / (2.0 * rho[j] * rho[j]);
          } else if (cdeamVersion == 2) {
            D_j = D_values[j];
          }
          double t2 = -rhoB[j];
          if (itype == speciesB) t2 += rho[j];
          fpair += D_j * rhoip * t2;
        }

        // This if-clause is only required for a ternary alloy.
        // Actually we don't need it at all because D_i should be zero
        // anyway if atom i has no concentration dependent interactions
        // (because it is not species A or B).

        if (x_i != -1.0) {
          double t1 = -rhoB[i];
          if (jtype == speciesB) t1 += rho[i];
          fpair += D_i * rhojp * t1;
        }

        double phip;
        double phi = PhiOfR(index, itype, jtype, recip, phip);
        if (itype == jtype || x_i == -1.0 || x_j == -1.0) {

          // Case of no concentration dependence.

          fpair += phip;
        } else {

          // We have a concentration dependence for the i-j interaction.

          double h=0.0;
          if (cdeamVersion == 1) {

            // Calculate h(x_i) polynomial function.

            double h_i = evalH(x_i);

            // Calculate h(x_j) polynomial function.

            double h_j = evalH(x_j);
            h = 0.5 * (h_i + h_j);
          } else if (cdeamVersion == 2) {

            // Average concentration.

            double x_ij = 0.5 * (x_i + x_j);

            // Calculate h(x_ij) polynomial function.

            h = evalH(x_ij);
          }
          fpair += h * phip;
          phi *= h;
        }

        // Divide by r_ij and negate to get forces from gradient.

        fpair /= -r;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) evdwl = phi;
        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairEAMCD::coeff(int narg, char **arg)
{
  PairEAMAlloy::coeff(narg, arg);

  // Make sure the EAM file is a CD-EAM binary alloy.

  if (setfl->nelements < 2)
    error->all(FLERR,"The EAM file must contain at least 2 elements to be "
                    "used with the eam/cd pair style.");

  // Read in the coefficients of the h polynomial from the end of the EAM file.

  read_h_coeff(arg[2]);

  // Determine which atom type is the A species and which is the B
  // species in the alloy.  By default take the first element (index 0)
  // in the EAM file as the A species and the second element (index 1)
  // in the EAM file as the B species.

  speciesA = -1;
  speciesB = -1;
  for (int i = 1; i <= atom->ntypes; i++) {
    if (map[i] == 0) {
      if (speciesA >= 0)
        error->all(FLERR,"The first element from the EAM file may only be mapped to a single atom type.");
      speciesA = i;
    }
    if (map[i] == 1) {
      if (speciesB >= 0)
        error->all(FLERR,"The second element from the EAM file may only be mapped to a single atom type.");
      speciesB = i;
    }
  }
  if (speciesA < 0)
    error->all(FLERR,"The first element from the EAM file must be mapped to exactly one atom type.");
  if (speciesB < 0)
    error->all(FLERR,"The second element from the EAM file must be mapped to exactly one atom type.");
}

/* ----------------------------------------------------------------------
   Reads in the h(x) polynomial coefficients
------------------------------------------------------------------------- */

void PairEAMCD::read_h_coeff(char *filename)
{
  if (comm->me == 0) {

    // Open potential file

    FILE *fptr;
    int convert_flag = unit_convert_flag;
    fptr = utils::open_potential(filename, lmp, &convert_flag);
    if (fptr == nullptr)
      error->one(FLERR,"Cannot open EAMCD potential file {}", filename);

    // h coefficients are stored at the end of the file.
    // Seek to end of file, read last part into a buffer and
    // then skip over lines in buffer until reaching the end.

    if ( (platform::fseek(fptr, platform::END_OF_FILE) < 0)
         || (platform::fseek(fptr, platform::ftell(fptr) - MAXLINE) < 0))
      error->one(FLERR,"Failure to seek to end-of-file for reading h(x) coeffs: {}",
                 utils::getsyserror());

    auto buf = new char[MAXLINE+1];
    auto rv = fread(buf,1,MAXLINE,fptr);
    if (rv == 0) error->one(FLERR,"Failure to read h(x) coeffs: {}", utils::getsyserror());
    buf[rv] = '\0';        // must 0-terminate buffer for string processing

    Tokenizer lines(buf, "\n");
    delete[] buf;

    std::string lastline;
    while (lines.has_next())
      lastline = lines.next();

    ValueTokenizer values(lastline);
    int degree = values.next_int();
    nhcoeff = degree+1;

    if ((int)values.count() != nhcoeff + 1 || nhcoeff < 1)
      error->one(FLERR, "Failed to read h(x) function coefficients in EAM file.");

    delete[] hcoeff;
    hcoeff = new double[nhcoeff];

    for (int i = 0; i < nhcoeff; ++i)
      hcoeff[i] = values.next_double();

    // Close the potential file.

    fclose(fptr);
  }

  MPI_Bcast(&nhcoeff, 1, MPI_INT, 0, world);
  if (comm->me != 0) {
    delete[] hcoeff;
    hcoeff = new double[nhcoeff];
  }
  MPI_Bcast(hcoeff, nhcoeff, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

int PairEAMCD::pack_forward_comm(int n, int *list, double *buf,
                                 int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  if (communicationStage == 2) {
    if (cdeamVersion == 1) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = fp[j];
        buf[m++] = rho[j];
        buf[m++] = rhoB[j];
        buf[m++] = D_values[j];
      }
      return m;
    } else if (cdeamVersion == 2) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = fp[j];
        buf[m++] = rho[j];
        buf[m++] = rhoB[j];
      }
      return m;
    } else return 0;
  } else if (communicationStage == 4) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = D_values[j];
    }
    return m;
  } else return 0;
}

/* ---------------------------------------------------------------------- */

void PairEAMCD::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (communicationStage == 2) {
    if (cdeamVersion == 1) {
      for (i = first; i < last; i++) {
        fp[i] = buf[m++];
        rho[i] = buf[m++];
        rhoB[i] = buf[m++];
        D_values[i] = buf[m++];
      }
    } else if (cdeamVersion == 2) {
      for (i = first; i < last; i++) {
        fp[i] = buf[m++];
        rho[i] = buf[m++];
        rhoB[i] = buf[m++];
      }
    }
  } else if (communicationStage == 4) {
    for (i = first; i < last; i++) {
      D_values[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */
int PairEAMCD::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (communicationStage == 1) {
    if (cdeamVersion == 1) {
      for (i = first; i < last; i++) {
        buf[m++] = rho[i];
        buf[m++] = rhoB[i];
        buf[m++] = D_values[i];
      }
      return m;
    } else if (cdeamVersion == 2) {
      for (i = first; i < last; i++) {
        buf[m++] = rho[i];
        buf[m++] = rhoB[i];
      }
      return m;
    } else return 0;
  } else if (communicationStage == 3) {
    for (i = first; i < last; i++) {
      buf[m++] = D_values[i];
    }
    return m;
  } else return 0;
}

/* ---------------------------------------------------------------------- */

void PairEAMCD::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (communicationStage == 1) {
    if (cdeamVersion == 1) {
      for (i = 0; i < n; i++) {
        j = list[i];
        rho[j] += buf[m++];
        rhoB[j] += buf[m++];
        D_values[j] += buf[m++];
      }
    } else if (cdeamVersion == 2) {
      for (i = 0; i < n; i++) {
        j = list[i];
        rho[j] += buf[m++];
        rhoB[j] += buf[m++];
      }
    }
  } else if (communicationStage == 3) {
    for (i = 0; i < n; i++) {
      j = list[i];
      D_values[j] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */
double PairEAMCD::memory_usage()
{
  double bytes = 2 * nmax * sizeof(double);
  return PairEAMAlloy::memory_usage() + bytes;
}

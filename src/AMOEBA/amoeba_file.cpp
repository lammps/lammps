/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel ator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_amoeba.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cctype>
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

enum{FFIELD,LITERATURE,ATOMTYPE,VDWL,VDWLPAIR,BSTRETCH,SBEND,ABEND,
     PAULI,DISPERSION,UB,OUTPLANE,TORSION,PITORSION,ATOMMULT,
     QPENETRATION,DIPPOLAR,QTRANSFER,UNKNOWN};
enum{ALLINGER,BUFFERED_14_7};
enum{ARITHMETIC,GEOMETRIC,CUBIC_MEAN,R_MIN,SIGMA,DIAMETER,HARMONIC,HHG,W_H};
enum{MUTUAL,OPT,TCG,DIRECT};
enum{NOFRAME,ZONLY,ZTHENX,BISECTOR,ZBISECT,THREEFOLD};
enum{GEAR,ASPC,LSQR};

#define MAXLINE 65536              // crazy big for TORSION-TORSION section
#define MAX_TYPE_PER_GROUP 6       // max types per AMOEBA group
#define MAX_FRAME_PER_TYPE 32      // max multipole frames for any AMOEBA type

#define DELTA_TYPE_CLASS 32
#define DELTA_VDWL_PAIR 16

#define BOHR 0.52917721067         // Bohr in Angstroms

// methods to read, parse, and store info from force field file

/* ----------------------------------------------------------------------
   set default values for items read from PRM and key files
------------------------------------------------------------------------- */

void PairAmoeba::set_defaults()
{
  for (int i = 0; i <= 4; i++) {
    special_hal[i] = 1.0;
    special_repel[i] = 1.0;
    special_disp[i] = 1.0;
    special_mpole[i] = 1.0;
    special_polar_pscale[i] = 1.0;
    special_polar_piscale[i] = 1.0;
    special_polar_wscale[i] = 1.0;
  }

  polar_dscale = 0.0;
  polar_uscale = 0.0;
}

/* ----------------------------------------------------------------------
   read PRM force field file
------------------------------------------------------------------------- */

void PairAmoeba::read_prmfile(char *filename)
{
  int n,nextflag;

  // open potential file

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE],next[MAXLINE];

  if (me == 0) {
    fptr = utils::open_potential(filename,lmp,nullptr);
    if (fptr == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open AMOEBA PRM file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read sections, one at a time
  // Force Field Definition section must come first
  // skip Literature References section
  // Atom Type Definitions must come before any other section
  // other sections can follow in any order

  // NOTE: don't use tokenize when not needed, doc string methods better
  // NOTE: how to insure each section had enough lines?

  int forcefield_flag = 0;
  int atomtype_flag = 0;
  int section;

  while (1) {
    if (me == 0) n = read_section_name(fptr,line);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n < 0) break;
    MPI_Bcast(line,n+1,MPI_CHAR,0,world);

    //printf("Section: %s\n",line);

    if (strstr(line,"Force Field") == line) section = FFIELD;
    else if (strstr(line,"Literature") == line) section = LITERATURE;
    else if (strstr(line,"Atom Type") == line) section = ATOMTYPE;
    else if (strstr(line,"Van der Waals Param") == line) section = VDWL;
    else if (strstr(line,"Van der Waals Pair") == line) section = VDWLPAIR;
    else if (strstr(line,"Bond Stretching") == line) section = BSTRETCH;
    else if (strstr(line,"Stretch-Bend") == line) section = SBEND;
    else if (strstr(line,"Angle Bending") == line) section = ABEND;
    else if (strstr(line,"Pauli Repulsion") == line) section = PAULI;
    else if (strstr(line,"Dispersion Param") == line) section = DISPERSION;
    else if (strstr(line,"Urey-Bradley") == line) section = UB;
    else if (strstr(line,"Out-of-Plane") == line) section = OUTPLANE;
    else if (strstr(line,"Torsional") == line) section = TORSION;
    else if (strstr(line,"Pi-Torsion") == line) section = PITORSION;
    else if (strstr(line,"Atomic Multipole") == line) section = ATOMMULT;
    else if (strstr(line,"Charge Penetration") == line) section = QPENETRATION;
    else if (strstr(line,"Dipole Polarizability") == line) section = DIPPOLAR;
    else if (strstr(line,"Charge Transfer") == line) section = QTRANSFER;
    else {
      if (me == 0) printf("Skipping section: %s\n",line);
      section = UNKNOWN;
    }

    if (forcefield_flag == 0 && section != FFIELD)
      error->all(FLERR,"Force Field is not first section of "
                 "pair amoeba potential file");
    if (section != FFIELD && section != LITERATURE && section != ATOMTYPE &&
        section != UNKNOWN && atomtype_flag == 0)
      error->all(FLERR,"Atom Type section of pair amoeba potential file "
                 "must come before all but Force Field section");

    if (section == FFIELD) forcefield_flag = 1;
    if (section == ATOMTYPE) atomtype_flag = 1;

    if (section == ATOMMULT) {
      for (int i = 1; i <= n_amtype; i++) nmultiframe[i] = 0;
    }

    nextflag = 0;

    while (1) {
      if (me == 0) n = read_section_line(fptr,line,nextflag,next);
      MPI_Bcast(&n,1,MPI_INT,0,world);
      if (n < 0) break;
      MPI_Bcast(line,n+1,MPI_CHAR,0,world);
      if (n == 0) break;    // line starting with #### = next section line

      // convert all chars in line to lower-case

      for (int i = 0; i < n; i++)
        line[i] = tolower(line[i]);

      char *copy;
      char **words;
      int nwords = tokenize(line,words,copy);

      if (section == FFIELD) file_ffield(nwords,words);
      else if (section == LITERATURE) file_literature(nwords,words);
      else if (section == ATOMTYPE) file_atomtype(nwords,words);
      else if (section == VDWL) file_vdwl(nwords,words);
      else if (section == VDWLPAIR) file_vdwl_pair(nwords,words);
      else if (section == BSTRETCH) file_bstretch(nwords,words);
      else if (section == SBEND) file_sbend(nwords,words);
      else if (section == ABEND) file_abend(nwords,words);
      else if (section == PAULI) file_pauli(nwords,words);
      else if (section == DISPERSION) file_dispersion(nwords,words);
      else if (section == UB) file_ub(nwords,words);
      else if (section == OUTPLANE) file_outplane(nwords,words);
      else if (section == TORSION) file_torsion(nwords,words);
      else if (section == PITORSION) file_pitorsion(nwords,words);
      else if (section == ATOMMULT) file_multipole(nwords,words);
      else if (section == QPENETRATION) file_charge_penetration(nwords,words);
      else if (section == DIPPOLAR) file_dippolar(nwords,words);
      else if (section == QTRANSFER) file_charge_transfer(nwords,words);
      else if (section == UNKNOWN) {}

      delete [] copy;
      delete [] words;
    }

    if (n < 0) break;
  }

  if (me == 0) fclose(fptr);

  if (forcefield_flag == 0 || atomtype_flag == 0)
    error->all(FLERR,"Pair amoeba potential file incomplete");
}

/* ----------------------------------------------------------------------
   read optional KEY file of one-line settings
------------------------------------------------------------------------- */

void PairAmoeba::read_keyfile(char *filename)
{
  double aprd,bprd,cprd;

  // default settings for which there are keyword options

  aprd = bprd = cprd = 0.0;

  vdwcut = 9.0;
  vdwtaper = 0.9 * vdwcut;
  repcut = 6.0;
  reptaper = 0.9 * repcut;
  dispcut = 9.0;
  disptaper = 0.9 * dispcut;
  mpolecut = 9.0;
  mpoletaper = 0.65 * mpolecut;
  ctrncut = 6.0;
  ctrntaper = 0.9 * ctrncut;

  ewaldcut = 7.0;
  dewaldcut = 7.0;
  usolvcut = 4.5;
  udiag = 2.0;

  dhal = 0.07;
  ghal = 0.12;

  use_ewald = use_dewald = 0;

  bseorder = 5;
  bsporder = 5;
  bsdorder = 4;

  aeewald = 0.4;
  apewald = 0.4;
  adewald = 0.4;

  use_pred = 0;
  polpred = LSQR;
  politer = 100;
  poleps = 1.0e-6;

  pcgprec = 1;
  pcgpeek = 1.0;
  pcgguess = 1;

  aeewald_key = apewald_key = adewald_key = 0;
  pmegrid_key = dpmegrid_key = 0;

  // done if keyfile not specified by pair_coeff command

  if (!filename) return;

  // open key file

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE],next[MAXLINE];
  if (me == 0) {
    fptr = utils::open_potential(filename,lmp,nullptr);
    if (fptr == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open AMOEBA key file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read lines, one at a time

  int n;
  char *ptr,*keyword;

  while (1) {
    if (me == 0) {
      ptr = fgets(line,MAXLINE,fptr);
      if (!ptr) n = -1;
      else n = strlen(line);
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n < 0) break;
    MPI_Bcast(line,n+1,MPI_CHAR,0,world);
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // convert all chars in line to lower-case

    for (int i = 0; i < n; i++)
      line[i] = tolower(line[i]);

    char *copy;
    char **words;
    int nwords = tokenize(line,words,copy);
    keyword = words[0];

    if (strstr(keyword,"#") || strstr(keyword,"!!") || !isalpha(keyword[0])) {

    } else if (strcmp(keyword,"a-axis") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      aprd = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"b-axis") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      bprd = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"c-axis") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      cprd = utils::numeric(FLERR,words[1],true,lmp);

    } else if (strcmp(keyword,"cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      double cut = utils::numeric(FLERR,words[1],true,lmp);
      vdwcut = repcut = dispcut = mpolecut = ctrncut = ewaldcut = dewaldcut = cut;
      vdwtaper = 0.9 * vdwcut;
      reptaper = 0.9 * repcut;
      disptaper = 0.9 * dispcut;
      mpoletaper = 0.65 * mpolecut;
      ctrntaper = 0.9 * ctrncut;
    } else if (strcmp(keyword,"taper") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      double taper = utils::numeric(FLERR,words[1],true,lmp);
      if (taper >= 1.0) {
        vdwtaper = reptaper = disptaper = mpoletaper = ctrntaper = taper;
      } else {
        taper = -taper;
        vdwtaper = taper * vdwcut;
        reptaper = taper * repcut;
        disptaper = taper * dispcut;
        mpoletaper = taper * mpolecut;
        ctrntaper = taper * ctrncut;
      }
    } else if (strcmp(keyword,"vdw-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      vdwcut = utils::numeric(FLERR,words[1],true,lmp);
      vdwtaper = 0.9 * vdwcut;
    } else if (strcmp(keyword,"repulsion-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      repcut = utils::numeric(FLERR,words[1],true,lmp);
      reptaper = 0.9 * repcut;
    } else if (strcmp(keyword,"dispersion-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      dispcut = utils::numeric(FLERR,words[1],true,lmp);
      disptaper = 0.9 * dispcut;
    } else if (strcmp(keyword,"mpole-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      mpolecut = utils::numeric(FLERR,words[1],true,lmp);
      mpoletaper = 0.65 * mpolecut;
    } else if (strcmp(keyword,"ctrn-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      ctrncut = utils::numeric(FLERR,words[1],true,lmp);
      ctrntaper = 0.9 * ctrncut;
    } else if (strcmp(keyword,"ewald-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      ewaldcut = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"dewald-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      dewaldcut = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"usolve-cutoff") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      usolvcut = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"usolve-diag") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      udiag = utils::numeric(FLERR,words[1],true,lmp);

    } else if (strcmp(keyword,"vdw-taper") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      vdwtaper = utils::numeric(FLERR,words[1],true,lmp);
      if (vdwtaper < 1.0) vdwtaper = -vdwtaper * vdwcut;
    } else if (strcmp(keyword,"repulsion-taper") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      reptaper = utils::numeric(FLERR,words[1],true,lmp);
      if (reptaper < 1.0) reptaper = -reptaper * repcut;
    } else if (strcmp(keyword,"dispersion-taper") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      disptaper = utils::numeric(FLERR,words[1],true,lmp);
      if (disptaper < 1.0) disptaper = -disptaper * vdwcut;
    } else if (strcmp(keyword,"mpole-taper") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      mpoletaper = utils::numeric(FLERR,words[1],true,lmp);
      if (mpoletaper < 1.0) mpoletaper = -mpoletaper * vdwcut;
    } else if (strcmp(keyword,"ctrn-taper") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      ctrntaper = utils::numeric(FLERR,words[1],true,lmp);
      if (ctrntaper < 1.0) ctrntaper = -ctrntaper * vdwcut;

    } else if (strcmp(keyword,"delta-halgren") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      dhal = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"gamma-halgren") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      ghal = utils::numeric(FLERR,words[1],true,lmp);

    } else if (strcmp(keyword,"ewald") == 0) {
      if (nwords != 1) error->all(FLERR,"AMOEBA keyfile line is invalid");
      use_ewald = 1;
    } else if (strcmp(keyword,"dewald") == 0) {
      if (nwords != 1) error->all(FLERR,"AMOEBA keyfile line is invalid");
      use_dewald = 1;

    } else if (strcmp(keyword,"pme-order") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      bseorder = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"ppme-order") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      bsporder = utils::numeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"dpme-order") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      bsdorder = utils::numeric(FLERR,words[1],true,lmp);

    } else if (strcmp(keyword,"pme-grid") == 0) {
      if (nwords != 2 && nwords != 4)
        error->all(FLERR,"AMOEBA keyfile line is invalid");
      if (nwords == 2)
        nefft1 = nefft2 = nefft3 = utils::numeric(FLERR,words[1],true,lmp);
      else {
        nefft1 = utils::numeric(FLERR,words[1],true,lmp);
        nefft2 = utils::numeric(FLERR,words[2],true,lmp);
        nefft3 = utils::numeric(FLERR,words[3],true,lmp);
      }
      pmegrid_key = 1;
    } else if (strcmp(keyword,"dpme-grid") == 0) {
      if (nwords != 2 && nwords != 4)
        error->all(FLERR,"AMOEBA keyfile line is invalid");
      if (nwords == 2)
        ndfft1 = ndfft2 = ndfft3 = utils::numeric(FLERR,words[1],true,lmp);
      else {
        ndfft1 = utils::numeric(FLERR,words[1],true,lmp);
        ndfft2 = utils::numeric(FLERR,words[2],true,lmp);
        ndfft3 = utils::numeric(FLERR,words[3],true,lmp);
      }
      dpmegrid_key = 1;

    } else if (strcmp(keyword,"ewald-alpha") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      aeewald = utils::numeric(FLERR,words[1],true,lmp);
      aeewald_key = 1;
    } else if (strcmp(keyword,"pewald-alpha") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      apewald = utils::numeric(FLERR,words[1],true,lmp);
      apewald_key = 1;
    } else if (strcmp(keyword,"dewald-alpha") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      adewald = utils::numeric(FLERR,words[1],true,lmp);
      adewald_key = 1;

    // polarization options

    } else if (strcmp(words[0],"polarization") == 0) {
      if (strcmp(words[1],"mutual") == 0) poltyp = MUTUAL;
      else if (strstr(words[1],"opt") == words[1]) {
        poltyp = OPT;
        if (strcmp(words[1],"opt") == 0) optorder = 4;
        else optorder = utils::inumeric(FLERR,&words[1][3],true,lmp);
        if (optorder < 1 || optorder > 6)
          error->all(FLERR,"Unrecognized polarization OPTn in AMOEBA FF file");
      } else if (strcmp(words[1],"tcg") == 0)
        error->all(FLERR,"Polarization TCG not yet supported in AMOEBA/HIPPO");
      else if (strcmp(words[1],"direct") == 0) poltyp = DIRECT;
      else error->all(FLERR,"Unrecognized polarization in AMOEBA FF file");

    } else if (strcmp(keyword,"polar-predict") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      if (strcmp(words[1],"gear") == 0) {
        polpred = GEAR;
        maxualt = 7;
      } else if (strcmp(words[1],"aspc") == 0) {
        polpred = ASPC;
        maxualt = 17;
      } else if (strcmp(words[1],"lsqr") == 0) {
        polpred = LSQR;
        maxualt = 7;
      } else error->all(FLERR,"AMOEBA keyfile line is invalid");
      use_pred = 1;
    } else if (strcmp(keyword,"polar-iter") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      politer = utils::inumeric(FLERR,words[1],true,lmp);
    } else if (strcmp(keyword,"polar-eps") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      poleps = utils::numeric(FLERR,words[1],true,lmp);

    } else if (strcmp(keyword,"pcg-precond") == 0) {
      if (nwords != 1) error->all(FLERR,"AMOEBA keyfile line is invalid");
      pcgprec = 1;
    } else if (strcmp(keyword,"pcg-noprecond") == 0) {
      if (nwords != 1) error->all(FLERR,"AMOEBA keyfile line is invalid");
      pcgprec = 0;
    } else if (strcmp(keyword,"pcg-guess") == 0) {
      if (nwords != 1) error->all(FLERR,"AMOEBA keyfile line is invalid");
      pcgguess = 1;
    } else if (strcmp(keyword,"pcg-noguess") == 0) {
      if (nwords != 1) error->all(FLERR,"AMOEBA keyfile line is invalid");
      pcgguess = 0;
    } else if (strcmp(keyword,"pcg-peek") == 0) {
      if (nwords != 2) error->all(FLERR,"AMOEBA keyfile line is invalid");
      pcgpeek = utils::numeric(FLERR,words[1],true,lmp);

    } else {}

    delete [] copy;
    delete [] words;
  }

  // close key file

  if (me == 0) fclose(fptr);

  // cutoff resets for long-range interactions

  if (use_ewald) mpolecut = ewaldcut;
  if (use_dewald) dispcut = dewaldcut;

  // error checks

  if (use_ewald || use_dewald) {
    if (domain->nonperiodic)
      error->all(FLERR,"AMOEBA KSpace requires fully periodic system");
  }

  if (aprd > 0.0 && (!domain->xperiodic || domain->xprd != aprd))
    error->all(FLERR,"AMOEBA abc prd does not match LAMMPS domain");
  if (bprd > 0.0 && (!domain->yperiodic || domain->yprd != bprd))
    error->all(FLERR,"AMOEBA abc prd does not match LAMMPS domain");
  if (cprd > 0.0 && (!domain->zperiodic || domain->zprd != cprd))
    error->all(FLERR,"AMOEBA abc prd does not match LAMMPS domain");
}

/* ---------------------------------------------------------------------- */

int PairAmoeba::read_section_name(FILE *fp, char *line)
{
  char dummy[MAXLINE];

  // loop on read line
  // return -1 if EOF
  // skip line if blank
  // tokenize
  // if ncount <= 2 skip line
  // if 1st word and last word != ## then error
  // pack line with 2nd thru next-to-last words
  // return length of line

  char *ptr,*copy;
  char **words;
  int nwords;

  while (1) {
    ptr = fgets(line,MAXLINE,fp);
    if (!ptr) return -1;
    if (strspn(line," \t\n\r") == strlen(line)) continue;
    nwords = tokenize(line,words,copy);
    if (nwords <= 2) {
      delete [] copy;
      delete [] words;
      continue;
    }
    break;
  }

  if ((strstr(words[0],"##") != words[0]) ||
      (strstr(words[nwords-1],"##") != words[nwords-1]))
    error->one(FLERR,"Section header of pair amoeba potential file is invalid");

  line[0] = '\0';
  for (int i = 1; i < nwords-1; i++) {
    if (i > 1) strcat(line," ");
    strcat(line,words[i]);
  }
  int n = strlen(line);

  delete [] copy;
  delete [] words;

  // skip next 2 lines of section header

  fgets(dummy,MAXLINE,fp);
  fgets(dummy,MAXLINE,fp);

  return n;
}

/* ---------------------------------------------------------------------- */

int PairAmoeba::read_section_line(FILE *fp, char *line,
                                  int &nextflag, char *next)
{
  // loop on read line
  // if next line defined, use it instead of read
  // return -1 if EOF
  // skip line if blank
  // tokenize
  // if first word starts with ####, return 0
  // if first word starts with # or !!, continue
  // append continuation lines to line
  //   until a line is blank, starts with #, starts with alpha char, or EOF
  // save next line read, and set nextflag for next call
  // return length of line

  char *ptr,*copy,*copy_next;
  char **words,**words_next;
  int nwords,nwords_next;

  copy = copy_next = NULL;
  words = words_next = NULL;

  while (1) {
    if (nextflag) {
      strcpy(line,next);
      nextflag = 0;
    } else {
      ptr = fgets(line,MAXLINE,fp);
      if (!ptr) return -1;
    }

    if (strspn(line," \t\n\r") == strlen(line)) continue;
    nwords = tokenize(line,words,copy);
    if (strstr(words[0],"####")) {
      delete [] words;
      delete [] copy;
      return 0;
    }
    if (strstr(words[0],"#") || strstr(words[0],"!!") || !isalpha(words[0][0])) {
      delete [] words;
      delete [] copy;
      words = NULL;
      copy = NULL;
      continue;
    }
    while (1) {
      ptr = fgets(next,MAXLINE,fp);
      if (!ptr) {
        nextflag = 0;
        delete [] words;
        delete [] copy;
        return strlen(line);
      }
      nwords_next = tokenize(next,words_next,copy_next);
      if (nwords_next == 0) break;
      if (words_next[0][0] == '#') break;
      if (isalpha(words_next[0][0])) break;
      strcat(line,next);
      delete [] words_next;
      delete [] copy_next;
    }
    nextflag = 1;
    break;
  }

  delete [] copy;
  delete [] words;
  delete [] copy_next;
  delete [] words_next;

  int n = strlen(line);
  return n;
}

/* ----------------------------------------------------------------------
   tokenize str into white-space separated words
   return nwords = number of words
   return words = vector of ptrs to each word
   also return copystr since words points into it
   if a word is !!, that word and the rest of the line is discarded
   IMPORTANT: caller must delete words and copystr
------------------------------------------------------------------------- */

int PairAmoeba::tokenize(char *str, char **&words, char *&copystr)
{
  int n = strlen(str) + 1;
  copystr = new char[n];
  strcpy(copystr,str);

  int nword = count_words(copystr);
  words = new char*[nword];

  nword = 0;
  char *word = strtok(copystr," \t\n\r\f");
  while (word) {
    words[nword++] = word;
    word = strtok(NULL," \t\n\r\f");
    if (word && strcmp(word,"!!") == 0) break;
  }

  return nword;
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
------------------------------------------------------------------------- */

int PairAmoeba::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = new char[n];
  strcpy(copy,line);

  if (strtok(copy," \t\n\r\f") == NULL) {
    delete [] copy;
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  delete [] copy;
  return n;
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_ffield(int nwords, char **words)
{
  double tmp;

  if (strcmp(words[0],"forcefield") == 0) {
    int n = strlen(words[1]) + 1;
    forcefield = new char[n];
    strcpy(forcefield,words[1]);

  } else if (strcmp(words[0],"bond-cubic") == 0) {
    bond_cubic = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"bond-quartic") == 0) {
    bond_quartic = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"angle-cubic") == 0) {
    angle_cubic = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"angle-quartic") == 0) {
    angle_quartic = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"angle-pentic") == 0) {
    angle_pentic = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"angle-sextic") == 0) {
    angle_sextic = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"opbendtype") == 0) {
    if (strcmp(words[1],"allinger") == 0) opbendtype = ALLINGER;
    else error->all(FLERR,"Unrecognized opbendtype in AMOEBA FF file");
  } else if (strcmp(words[0],"opbend-cubic") == 0) {
    opbend_cubic = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"opbend-quartic") == 0) {
    opbend_quartic = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"opbend-pentic") == 0) {
    opbend_pentic = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"opbend-sextic") == 0) {
    opbend_sextic = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"torsionunit") == 0) {
    torsion_unit = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"vdwtype") == 0) {
    if (strcmp(words[1],"buffered-14-7") == 0) vdwtype = BUFFERED_14_7;
    else error->all(FLERR,"Unrecognized vdwtype in AMOEBA FF file");
  } else if (strcmp(words[0],"radiusrule") == 0) {
    if (strcmp(words[1],"arithmetic") == 0) radius_rule = ARITHMETIC;
    else if (strcmp(words[1],"geometric") == 0) radius_rule = GEOMETRIC;
    else if (strcmp(words[1],"cubic-mean") == 0) radius_rule = CUBIC_MEAN;
    else error->all(FLERR,"Unrecognized radiusrule in AMOEBA FF file");
  } else if (strcmp(words[0],"radiustype") == 0) {
    if (strcmp(words[1],"r-min") == 0) radius_type = R_MIN;
    else if (strcmp(words[1],"sigma") == 0) radius_type = SIGMA;
    else error->all(FLERR,"Unrecognized radiustype in AMOEBA FF file");
  } else if (strcmp(words[0],"radiussize") == 0) {
    if (strcmp(words[1],"diameter") == 0) radius_size = DIAMETER;
    else error->all(FLERR,"Unrecognized radiussize in AMOEBA FF file");
  } else if (strcmp(words[0],"epsilonrule") == 0) {
    if (strcmp(words[1],"arithmetic") == 0) epsilon_rule = ARITHMETIC;
    else if (strcmp(words[1],"geometric") == 0) epsilon_rule = GEOMETRIC;
    else if (strcmp(words[1],"harmonic") == 0) epsilon_rule = HARMONIC;
    else if (strcmp(words[1],"hhg") == 0) epsilon_rule = HHG;
    else if (strcmp(words[1],"w-h") == 0) epsilon_rule = W_H;
    else error->all(FLERR,"Unrecognized epsilonrule in AMOEBA FF file");

  } else if (strcmp(words[0],"dielectric") == 0) {
    am_dielectric = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"polarization") == 0) {
    if (strcmp(words[1],"mutual") == 0) poltyp = MUTUAL;
    else if (strstr(words[1],"opt") == words[1]) {
      poltyp = OPT;
      if (strcmp(words[1],"opt") == 0) optorder = 4;
      else optorder = utils::inumeric(FLERR,&words[1][3],true,lmp);
      if (optorder < 1 || optorder > 6)
        error->all(FLERR,"Unrecognized polarization OPTn in AMOEBA FF file");
    } else if (strcmp(words[1],"tcg") == 0)
      error->all(FLERR,"Polarization TCG not yet supported in AMOEBA/HIPPO");
    else if (strcmp(words[1],"direct") == 0) poltyp = DIRECT;
    else error->all(FLERR,"Unrecognized polarization in AMOEBA FF file");

  // NOTE: enable all variants of special settings
  //       do these need to be set to defaults if don't appear in file?

  } else if (strcmp(words[0],"vdw-12-scale") == 0) {
    special_hal[1] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"vdw-13-scale") == 0) {
    special_hal[2] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"vdw-14-scale") == 0) {
    special_hal[3] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"vdw-15-scale") == 0) {
    special_hal[4]  = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"rep-12-scale") == 0) {
    special_repel[1] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"rep-13-scale") == 0) {
    special_repel[2] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"rep-14-scale") == 0) {
    special_repel[3] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"rep-15-scale") == 0) {
    special_repel[4]  = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"disp-12-scale") == 0) {
    special_disp[1] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"disp-13-scale") == 0) {
    special_disp[2] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"disp-14-scale") == 0) {
    special_disp[3] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"disp-15-scale") == 0) {
    special_disp[4]  = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"mpole-12-scale") == 0) {
    special_mpole[1] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"mpole-13-scale") == 0) {
    special_mpole[2] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"mpole-14-scale") == 0) {
    special_mpole[3] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"mpole-15-scale") == 0) {
    special_mpole[4] = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"polar-12-scale") == 0) {
    special_polar_pscale[1] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"polar-13-scale") == 0) {
    special_polar_pscale[2] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"polar-14-scale") == 0) {
    special_polar_pscale[3] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"polar-15-scale") == 0) {
    special_polar_pscale[4] = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"polar-12-intra") == 0) {
    special_polar_piscale[1] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"polar-13-intra") == 0) {
    special_polar_piscale[2] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"polar-14-intra") == 0) {
    special_polar_piscale[3] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"polar-15-intra") == 0) {
    special_polar_piscale[4] = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"induce-12-scale") == 0) {
    special_polar_wscale[1] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"induce-13-scale") == 0) {
    special_polar_wscale[2] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"induce-14-scale") == 0) {
    special_polar_wscale[3] = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"induce-15-scale") == 0) {
    special_polar_wscale[4] = utils::numeric(FLERR,words[1],true,lmp);

  } else if (strcmp(words[0],"direct-11-scale") == 0) {
    polar_dscale = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"direct-12-scale") == 0) {
    // NOTE: could error check that value = 1
  } else if (strcmp(words[0],"direct-13-scale") == 0) {
  } else if (strcmp(words[0],"direct-14-scale") == 0) {

  } else if (strcmp(words[0],"mutual-11-scale") == 0) {
    polar_uscale = utils::numeric(FLERR,words[1],true,lmp);
  } else if (strcmp(words[0],"mutual-12-scale") == 0) {
    // NOTE: could error check that value = 1
  } else if (strcmp(words[0],"mutual-13-scale") == 0) {
  } else if (strcmp(words[0],"mutual-14-scale") == 0) {

  } else {
    char str[128];
    sprintf(str,"Unrecognized pair amoeba force field definition: %s",words[0]);
    error->all(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_literature(int nwords, char **words)
{
  // do nothing, this section is skipped
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_atomtype(int nwords, char **words)
{
  if (nwords < 8)
    error->all(FLERR,"AMOEBA atom type line is invalid");
  if (strcmp(words[0],"atom") != 0)
    error->all(FLERR,"AMOEBA atom type line is invalid");

  int itype = utils::inumeric(FLERR,words[1],true,lmp);
  int iclass = utils::inumeric(FLERR,words[2],true,lmp);

  // grow per-type and per-class vecs/arrays as needed

  allocate_type_class(itype,iclass);

  n_amtype = MAX(n_amtype,itype);
  n_amclass = MAX(n_amclass,iclass);

  // store words from line

  amtype_defined[itype] = 1;
  amclass_defined[iclass] = 1;
  amtype2class[itype] = iclass;

  atomic_num[itype] = utils::inumeric(FLERR,words[nwords-3],true,lmp);
  am_mass[itype] = utils::numeric(FLERR,words[nwords-2],true,lmp);
  valence[itype] = utils::inumeric(FLERR,words[nwords-1],true,lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_vdwl(int nwords, char **words)
{
  if (nwords != 4 && nwords != 5)
    error->all(FLERR,"AMOEBA Van der Waals line is invalid");
  if (strcmp(words[0],"vdw") != 0)
    error->all(FLERR,"AMOEBA Van der Waals line is invalid");

  int iclass = utils::inumeric(FLERR,words[1],true,lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR,"AMOEBA Van der Waals type index is invalid");

  vdwl_sigma[iclass] = utils::numeric(FLERR,words[2],true,lmp);
  vdwl_eps[iclass] = utils::numeric(FLERR,words[3],true,lmp);
  if (nwords == 4) kred[iclass] = 0.0;
  else kred[iclass] = utils::numeric(FLERR,words[4],true,lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_vdwl_pair(int nwords, char **words)
{
  if (nwords != 5)
    error->all(FLERR,"AMOEBA Van der Waals pair line is invalid");
  if (strcmp(words[0],"vdwpr") != 0)
    error->all(FLERR,"AMOEBA Van der Waals pair line is invalid");

  if (nvdwl_pair == max_vdwl_pair) {
    max_vdwl_pair += DELTA_VDWL_PAIR;
    memory->grow(vdwl_class_pair,max_vdwl_pair,2,"amoeba:vdwl_class_pair");
    memory->grow(vdwl_sigma_pair,max_vdwl_pair,"amoeba:vdwl_sigma_pair");
    memory->grow(vdwl_eps_pair,max_vdwl_pair,"amoeba:vdwl_eps_pair");
  }

  vdwl_class_pair[nvdwl_pair][0] = utils::inumeric(FLERR,words[1],true,lmp);
  vdwl_class_pair[nvdwl_pair][1] = utils::inumeric(FLERR,words[2],true,lmp);
  vdwl_sigma_pair[nvdwl_pair] = utils::numeric(FLERR,words[3],true,lmp);
  vdwl_eps_pair[nvdwl_pair] = utils::numeric(FLERR,words[4],true,lmp);
  nvdwl_pair++;
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_bstretch(int nwords, char **words)
{
  if (nwords < 5)
    error->all(FLERR,"AMOEBA bond stretch line is invalid");
  if (strcmp(words[0],"bond") != 0)
    error->all(FLERR,"AMOEBA bond stretch line is invalid");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_sbend(int nwords, char **words)
{
  if (nwords != 6)
    error->all(FLERR,"AMOEBA strectch-bend line is invalid");
  if (strstr(words[0],"strbnd") != words[0])
    error->all(FLERR,"AMOEBA strectch-bend line is invalid");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_abend(int nwords, char **words)
{
  if (nwords < 6)
    error->all(FLERR,"AMOEBA angle bending line is invalid");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_pauli(int nwords, char **words)
{
  if (nwords < 5)
    error->all(FLERR,"AMOEBA Pauli repulsion line is invalid");
  if (strstr(words[0],"repulsion") != words[0])
    error->all(FLERR,"AMOEBA Pauli repulsion line is invalid");

  int itype = utils::inumeric(FLERR,words[1],true,lmp);
  if (itype < 1 || itype > n_amtype)
    error->all(FLERR,"AMOEBA Pauli repulsion type index is invalid");

  // negate the elepr setting

  sizpr[itype] = utils::numeric(FLERR,words[2],true,lmp);
  dmppr[itype] = utils::numeric(FLERR,words[3],true,lmp);
  elepr[itype] = - fabs(utils::numeric(FLERR,words[4],true,lmp));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_dispersion(int nwords, char **words)
{
  if (nwords < 4)
    error->all(FLERR,"AMOEBA dispersion line is invalid");
  if (strstr(words[0],"dispersion") != words[0])
    error->all(FLERR,"AMOEBA dipersion line is invalid");

  int iclass = utils::inumeric(FLERR,words[1],true,lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR,"AMOEBA dispersion class index is invalid");

  csix[iclass] = utils::numeric(FLERR,words[2],true,lmp);
  adisp[iclass] = utils::numeric(FLERR,words[3],true,lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_ub(int nwords, char **words)
{
  if (nwords != 6)
    error->all(FLERR,"AMOEBA Urey-Bradley line is invalid");
  if (strstr(words[0],"ureybrad") != words[0])
    error->all(FLERR,"AMOEBA Urey-Bradley line is invalid");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_outplane(int nwords, char **words)
{
  if (nwords != 6)
    error->all(FLERR,"AMOEBA out-of-plane bend line is invalid");
  if (strstr(words[0],"opbend") != words[0])
    error->all(FLERR,"AMOEBA out-of-plane bend line is invalid");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_torsion(int nwords, char **words)
{
  if (nwords != 14)
    error->all(FLERR,"AMOEBA torsional line is invalid");
  if (strstr(words[0],"torsion") != words[0])
    error->all(FLERR,"AMOEBA torsional line is invalid");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_pitorsion(int nwords, char **words)
{
  if (nwords != 4)
    error->all(FLERR,"AMOEBA pi-torsion line is invalid");
  if (strstr(words[0],"pitors") != words[0])
    error->all(FLERR,"AMOEBA pi-torsion line is invalid");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_multipole(int nwords, char **words)
{
  if (nwords < 12 || nwords > 15)
    error->all(FLERR,"AMOEBA atomic multipole line is invalid");
  if (strstr(words[0],"multipole") != words[0])
    error->all(FLERR,"AMOEBA atomic multipole line is invalid");

  int itype = utils::inumeric(FLERR,words[1],true,lmp);
  if (itype < 1 || itype > n_amtype)
    error->all(FLERR,"AMOEBA atomic multipole type index is invalid");

  int iframe = nmultiframe[itype];
  if (iframe == MAX_FRAME_PER_TYPE)
    error->all(FLERR,"AMOEBA MAX_FRAME_PER_TYPE is too small");

  int extra;
  if (nwords == 12) {
    zpole[itype][iframe] = xpole[itype][iframe] = ypole[itype][iframe] = 0;
    extra = 2;
  } else if (nwords == 13) {
    zpole[itype][iframe] = utils::inumeric(FLERR,words[2],true,lmp);
    xpole[itype][iframe] = ypole[itype][iframe] = 0;
    extra = 3;
  } else if (nwords == 14) {
    zpole[itype][iframe] = utils::inumeric(FLERR,words[2],true,lmp);
    xpole[itype][iframe] = utils::inumeric(FLERR,words[3],true,lmp);
    ypole[itype][iframe] = 0;
    extra = 4;
  } else if (nwords == 15) {
    zpole[itype][iframe] = utils::inumeric(FLERR,words[2],true,lmp);
    xpole[itype][iframe] = utils::inumeric(FLERR,words[3],true,lmp);
    ypole[itype][iframe] = utils::inumeric(FLERR,words[4],true,lmp);
    extra = 5;
  }

  for (int i = 0; i < 10; i++)
    fpole[itype][iframe][i] = utils::numeric(FLERR,words[extra+i],true,lmp);

  // convert fpole to be 13 values by symmetrizing quadrupole 3x3 matrix
  // xx yx yy zx zy zz --> xx xy xz yz yy yz zx zy zz

  double xx = fpole[itype][iframe][4];
  double yx = fpole[itype][iframe][5];
  double yy = fpole[itype][iframe][6];
  double zx = fpole[itype][iframe][7];
  double zy = fpole[itype][iframe][8];
  double zz = fpole[itype][iframe][9];

  fpole[itype][iframe][4] = xx;
  fpole[itype][iframe][8] = yy;
  fpole[itype][iframe][12] = zz;
  fpole[itype][iframe][5] = fpole[itype][iframe][7] = yx;
  fpole[itype][iframe][6] = fpole[itype][iframe][10] = zx;
  fpole[itype][iframe][9] = fpole[itype][iframe][11] = zy;

  // rescale pole values to real units
  // convert the dipole and quadrupole moments to Angstroms
  // quadrupole terms divided by 3 for use as traceless values

  for (int i = 1; i < 4; i++)
    fpole[itype][iframe][i] *= BOHR;
  for (int i = 4; i < 13; i++)
    fpole[itype][iframe][i] *= BOHR*BOHR / 3.0;

  // set mpaxis from xyz pole values
  // uses both positive and negative values

  mpaxis[itype][iframe] = ZTHENX;
  if (zpole[itype][iframe] == 0) mpaxis[itype][iframe] = NOFRAME;
  if (zpole[itype][iframe] != 0 && xpole[itype][iframe] == 0)
    mpaxis[itype][iframe] = ZONLY;
  if (zpole[itype][iframe] < 0 || xpole[itype][iframe] < 0)
    mpaxis[itype][iframe] = BISECTOR;
  if (xpole[itype][iframe] < 0 && ypole[itype][iframe] < 0)
    mpaxis[itype][iframe] = ZBISECT;
  int xyzmax = MAX(zpole[itype][iframe],xpole[itype][iframe]);
  xyzmax = MAX(xyzmax,ypole[itype][iframe]);
  if (xyzmax < 0) mpaxis[itype][iframe] = THREEFOLD;
  if (mpaxis[itype][iframe] < 0) error->all(FLERR,"Mpaxis value not set");

  // now reset xyz pole to positive values

  if (xpole[itype][iframe] < 0) xpole[itype][iframe] = -xpole[itype][iframe];
  if (ypole[itype][iframe] < 0) ypole[itype][iframe] = -ypole[itype][iframe];
  if (zpole[itype][iframe] < 0) zpole[itype][iframe] = -zpole[itype][iframe];

  nmultiframe[itype]++;
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_charge_penetration(int nwords, char **words)
{
  if (nwords < 4)
    error->all(FLERR,"AMOEBA charge penetration line is invalid");
  if (strstr(words[0],"chgpen") != words[0])
    error->all(FLERR,"AMOEBA charge penetration line is invalid");

  int iclass = utils::inumeric(FLERR,words[1],true,lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR,"AMOEBA charge penetration class index is invalid");

  pcore[iclass] = fabs(utils::numeric(FLERR,words[2],true,lmp));
  palpha[iclass] = utils::numeric(FLERR,words[3],true,lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_dippolar(int nwords, char **words)
{
  int nparams;
  if (amoeba) nparams = 4;
  if (hippo) nparams = 3;

  if (nwords < nparams)
    error->all(FLERR,"AMOEBA dipole polarizability line is invalid");
  if (strstr(words[0],"polarize") != words[0])
    error->all(FLERR,"AMOEBA dipole polarizability line is invalid");

  int itype = utils::inumeric(FLERR,words[1],true,lmp);
  if (itype < 1 || itype > n_amtype)
    error->all(FLERR,"AMOEBA dipole polarizabilty type index is invalid");

  polarity[itype] = utils::numeric(FLERR,words[2],true,lmp);
  pdamp[itype] = pow(polarity[itype],1.0/6.0);
  if (amoeba) thole[itype] = utils::numeric(FLERR,words[3],true,lmp);

  // eventually AMOEBA+ files will set dirdamp

  dirdamp[itype] = 0.0;

  int ngroup = nwords - nparams;
  if (ngroup > MAX_TYPE_PER_GROUP)
    error->all(FLERR,"AMOEBA MAX_TYPE_PER_GROUP is too small");

  npolgroup[itype] = ngroup;
  for (int igroup = 0; igroup < ngroup; igroup++)
    polgroup[itype][igroup] =
      utils::inumeric(FLERR,words[nparams+igroup],true,lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_charge_transfer(int nwords, char **words)
{
  if (nwords < 4)
    error->all(FLERR,"AMOEBA charge transfer line is invalid");
  if (strstr(words[0],"chgtrn") != words[0])
    error->all(FLERR,"AMOEBA charge transfer line is invalid");

  int iclass = utils::inumeric(FLERR,words[1],true,lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR,"AMOEBA charge transfer class index is invalid");

  chgct[iclass] = utils::numeric(FLERR,words[2],true,lmp);
  dmpct[iclass] = utils::numeric(FLERR,words[3],true,lmp);
}

/* ----------------------------------------------------------------------
   initialize per type and per class data to NULL
------------------------------------------------------------------------- */

void PairAmoeba::initialize_type_class()
{
  n_amtype = n_amclass = 0;
  max_amtype = max_amclass = 0;
  nvdwl_pair = max_vdwl_pair = 0;

  // per type data

  amtype_defined = NULL;
  amtype2class = NULL;
  atomic_num = NULL;
  valence = NULL;
  am_mass = NULL;
  am_q = NULL;
  am_mu = NULL;
  npolgroup = NULL;
  polgroup = NULL;
  polarity = NULL;
  pdamp = NULL;
  thole = NULL;
  dirdamp = NULL;
  sizpr = NULL;
  dmppr = NULL;
  elepr = NULL;

  nmultiframe = NULL;
  mpaxis = NULL;
  xpole = NULL;
  ypole = NULL;
  zpole = NULL;
  fpole = NULL;

  // per class data

  amclass_defined = NULL;
  vdwl_eps = NULL;
  vdwl_sigma = NULL;
  kred = NULL;
  csix = NULL;
  adisp = NULL;
  chgct = NULL;
  dmpct = NULL;
  pcore = NULL;
  palpha = NULL;

  // other

  vdwl_class_pair = NULL;
  vdwl_sigma_pair = NULL;
  vdwl_eps_pair = NULL;
}

/* ----------------------------------------------------------------------
   allocate per type and per class data
   vecs/arrays store info for itype = 1 to N_amtype inclusive
   vecs/arrays store info for iclass = 1 to N_amclass inclusive
   itype,iclass = line just read from AMOEBA force field file
------------------------------------------------------------------------- */

void PairAmoeba::allocate_type_class(int itype, int iclass)
{
  if (itype >= max_amtype) {
    while (itype >= max_amtype) max_amtype += DELTA_TYPE_CLASS;

    memory->grow(amtype_defined,max_amtype,"amoeba:amtype_defined");
    memory->grow(amtype2class,max_amtype,"amoeba:amtype2class");

    memory->grow(atomic_num,max_amtype,"amoeba:atomic_num");
    memory->grow(valence,max_amtype,"amoeba:valence");
    memory->grow(am_mass,max_amtype,"amoeba:am_mass");
    memory->grow(am_q,max_amtype,"amoeba:am_q");
    memory->grow(am_mu,max_amtype,3,"amoeba:am_mu");
    memory->grow(npolgroup,max_amtype,"amoeba:npolgroup");
    memory->grow(polgroup,max_amtype,MAX_TYPE_PER_GROUP,"amoeba:polgroup");
    memory->grow(polarity,max_amtype,"amoeba:polarity");
    memory->grow(pdamp,max_amtype,"amoeba:pdamp");
    memory->grow(thole,max_amtype,"amoeba:thole");
    memory->grow(dirdamp,max_amtype,"amoeba:dirdamp");
    memory->grow(sizpr,max_amtype,"amoeba:sizpr");
    memory->grow(dmppr,max_amtype,"amoeba:dmppr");
    memory->grow(elepr,max_amtype,"amoeba:elepr");

    memory->grow(nmultiframe,max_amtype,"amoeba:nummulti");
    memory->grow(mpaxis,max_amtype,MAX_FRAME_PER_TYPE,"amoeba:mpaxis");
    memory->grow(xpole,max_amtype,MAX_FRAME_PER_TYPE,"amoeba:xpole");
    memory->grow(ypole,max_amtype,MAX_FRAME_PER_TYPE,"amoeba:ypole");
    memory->grow(zpole,max_amtype,MAX_FRAME_PER_TYPE,"amoeba:zpole");
    memory->grow(fpole,max_amtype,MAX_FRAME_PER_TYPE,13,"amoeba:fpole");
  }

  if (iclass >= max_amclass) {
    while (iclass >= max_amclass) max_amclass += DELTA_TYPE_CLASS;

    memory->grow(amclass_defined,max_amtype,"amoeba:amclass_defined");

    memory->grow(vdwl_eps,max_amclass,"amoeba:vdwl_eps");
    memory->grow(vdwl_sigma,max_amclass,"amoeba:vdwl_sigma");
    memory->grow(kred,max_amclass,"amoeba:kred");
    memory->grow(csix,max_amclass,"amoeba:csix");
    memory->grow(adisp,max_amclass,"amoeba:adisp");
    memory->grow(chgct,max_amclass,"amoeba:chgct");
    memory->grow(dmpct,max_amclass,"amoeba:dmpct");
    memory->grow(pcore,max_amtype,"amoeba:pcore");
    memory->grow(palpha,max_amtype,"amoeba:palpha");
  }
}

/* ----------------------------------------------------------------------
   deallocate per type and per class data
------------------------------------------------------------------------- */

void PairAmoeba::deallocate_type_class()
{
  // per type data

  memory->destroy(amtype_defined);
  memory->destroy(amtype2class);
  memory->destroy(atomic_num);
  memory->destroy(valence);
  memory->destroy(am_mass);
  memory->destroy(am_q);
  memory->destroy(am_mu);
  memory->destroy(npolgroup);
  memory->destroy(polgroup);
  memory->destroy(polarity);
  memory->destroy(pdamp);
  memory->destroy(thole);
  memory->destroy(dirdamp);
  memory->destroy(sizpr);
  memory->destroy(dmppr);
  memory->destroy(elepr);

  memory->destroy(nmultiframe);
  memory->destroy(mpaxis);
  memory->destroy(xpole);
  memory->destroy(ypole);
  memory->destroy(zpole);
  memory->destroy(fpole);

  // per class data

  memory->destroy(amclass_defined);
  memory->destroy(vdwl_eps);
  memory->destroy(vdwl_sigma);
  memory->destroy(kred);
  memory->destroy(csix);
  memory->destroy(adisp);
  memory->destroy(chgct);
  memory->destroy(dmpct);
  memory->destroy(pcore);
  memory->destroy(palpha);

  // other

  memory->destroy(vdwl_class_pair);
  memory->destroy(vdwl_sigma_pair);
  memory->destroy(vdwl_eps_pair);
}

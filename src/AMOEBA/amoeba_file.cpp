// clang-format off
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

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "utils.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>
#include <cctype>

using namespace LAMMPS_NS;

enum{UNKNOWN,FFIELD,LITERATURE,ATOMTYPE,VDWL,VDWLPAIR,BSTRETCH,SBEND,ABEND,
     PAULI,DISPERSION,UB,OUTPLANE,TORSION,PITORSION,ATOMMULT,
     QPENETRATION,DIPPOLAR,QTRANSFER,END_OF_FILE};
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
  optorder = 0;
  maxualt = 7;
  tcgnab = 0;

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
  // open potential file

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = utils::open_potential(filename, lmp, nullptr);
    if (fptr == nullptr)
      error->one(FLERR, "Cannot open {} PRM file {}: {}", utils::uppercase(mystyle), filename,
                 utils::getsyserror());
  }

  // read sections, one at a time
  // Force Field Definition section must come first
  // skip Literature References section
  // Atom Type Definitions must come before any other section
  // other sections can follow in any order

  bool forcefield_flag = false;
  bool atomtype_flag = false;
  int section;
  int nline = 0;

  while (true) {
    if (me == 0) {
      clearerr(fptr);
      section = END_OF_FILE;
      while (fgets(line, MAXLINE, fptr)) {
        ++nline;
        if (utils::strmatch(line, "^\\s*##\\s+\\S+.*##\\s*$")) {
          auto trimmed = utils::trim(line);
          if (utils::strmatch(trimmed, "^##\\s*Force Field"))
            section = FFIELD;
          else if (utils::strmatch(trimmed, "^##\\s*Literature"))
            section = LITERATURE;
          else if (utils::strmatch(trimmed, "^##\\s*Atom Type"))
            section = ATOMTYPE;
          else if (utils::strmatch(trimmed, "^##\\s*Van der Waals Param"))
            section = VDWL;
          else if (utils::strmatch(trimmed, "^##\\s*Van der Waals Pair"))
            section = VDWLPAIR;
          else if (utils::strmatch(trimmed, "^##\\s*Bond Stretching"))
            section = BSTRETCH;
          else if (utils::strmatch(trimmed, "^##\\s*Stretch-Bend"))
            section = SBEND;
          else if (utils::strmatch(trimmed, "^##\\s*Angle Bending"))
            section = ABEND;
          else if (utils::strmatch(trimmed, "^##\\s*Pauli Repulsion"))
            section = PAULI;
          else if (utils::strmatch(trimmed, "^##\\s*Dispersion Param"))
            section = DISPERSION;
          else if (utils::strmatch(trimmed, "^##\\s*Urey-Bradley"))
            section = UB;
          else if (utils::strmatch(trimmed, "^##\\s*Out-of-Plane"))
            section = OUTPLANE;
          else if (utils::strmatch(trimmed, "^##\\s*Torsional"))
            section = TORSION;
          else if (utils::strmatch(trimmed, "^##\\s*Pi-Torsion"))
            section = PITORSION;
          else if (utils::strmatch(trimmed, "^##\\s*Atomic Multipole"))
            section = ATOMMULT;
          else if (utils::strmatch(trimmed, "^##\\s*Charge Penetration"))
            section = QPENETRATION;
          else if (utils::strmatch(trimmed, "^##\\s*Dipole Polarizability"))
            section = DIPPOLAR;
          else if (utils::strmatch(trimmed, "^##\\s*Charge Transfer"))
            section = QTRANSFER;
          else {
            section = UNKNOWN;
            utils::logmesg(lmp, "Skipping section: {}\n", trimmed.substr(2, trimmed.size() - 4));
          }

          // skip two lines following section head keyword
          fgets(line, MAXLINE, fptr);
          fgets(line, MAXLINE, fptr);
          nline += 2;
          break;
        }
      }
      if (ferror(fptr))
        error->one(FLERR, "Problem reading {} PRM file {}:{} {}", utils::uppercase(mystyle), nline,
                   filename, utils::getsyserror());
      if (feof(fptr)) section = END_OF_FILE;
    }
    MPI_Bcast(&section, 1, MPI_INT, 0, world);
    if (section == END_OF_FILE) break;

    // sanity checks
    if (!forcefield_flag && (section != FFIELD))
      error->all(FLERR, "Force Field is not first section of pair {} potential file", mystyle);

    if ((section > ATOMTYPE) && !atomtype_flag)
      error->all(FLERR,
                 "Atom Type section of pair {} potential file must "
                 "come before all but the Force Field section",
                 mystyle);

    if (section == FFIELD) forcefield_flag = true;
    if (section == ATOMTYPE) atomtype_flag = true;
    if (section == ATOMMULT) {
      for (int i = 1; i <= n_amtype; i++) nmultiframe[i] = 0;
    }

    char next[MAXLINE];
    next[0] = '\0';
    bool has_next = false;
    int n;
    while (true) {
      if (me == 0) {
        while (true) {
          line[0] = '\0';
          n = -1;
          if (has_next) strcpy(line, next);
          has_next = false;
          clearerr(fptr);
          while (fgets(next, MAXLINE, fptr)) {
            ++nline;
            auto trimmed = utils::trim(next);
            // chop off !! comments
            std::size_t pos = trimmed.find("!!");
            if (pos != std::string::npos) trimmed = trimmed.substr(0, pos);

            // append to line if next line starts with a number
            if (utils::is_double(utils::strfind(trimmed, "^\\S+"))) {
              strcat(line, " ");
              strcat(line, trimmed.c_str());
              has_next = false;
            } else {
              strcpy(next, trimmed.c_str());
              has_next = true;
              break;
            }
          }
          if (ferror(fptr))
            error->one(FLERR, "Problem reading {} PRM file {}:{} {}", utils::uppercase(mystyle),
                       nline, filename, utils::getsyserror());

          auto trimmed = utils::trim(line);

          // start of next section
          if (utils::strmatch(trimmed, "^####+$")) {
            n = 0;
            break;
          }

          // skip concatenated line with commented out keyword
          if (utils::strmatch(trimmed, "^#\\w+")) continue;

          // exit loop if line is not empty
          if (!trimmed.empty()) {
            strcpy(line, trimmed.c_str());
            n = strlen(line) + 1;
            break;
          }
          if (feof(fptr)) {
            n = -1;
            break;
          }
        }
      }
      MPI_Bcast(&n, 1, MPI_INT, 0, world);
      if (n < 0) break;
      MPI_Bcast(line, n, MPI_CHAR, 0, world);

      // next section

      if (n == 0) break;

      // convert line to lowercase and get list of words
      // XXX do we need to use lowercase? Preserving case may be useful later for type lables.
      // XXX We could also use utils::split_words() which slower but can handle quotes.
      auto words = Tokenizer(utils::lowercase(line)).as_vector();

      switch (section) {
        case FFIELD:
          file_ffield(words, nline - 1);
          break;
        case LITERATURE:
          file_literature(words, nline - 1);
          break;
        case ATOMTYPE:
          file_atomtype(words, nline - 1);
          break;
        case VDWL:
          file_vdwl(words, nline - 1);
          break;
        case VDWLPAIR:
          file_vdwl_pair(words, nline - 1);
          break;
        case BSTRETCH:
          file_bstretch(words, nline - 1);
          break;
        case SBEND:
          file_sbend(words, nline - 1);
          break;
        case ABEND:
          file_abend(words, nline - 1);
          break;
        case PAULI:
          file_pauli(words, nline - 1);
          break;
        case DISPERSION:
          file_dispersion(words, nline - 1);
          break;
        case UB:
          file_ub(words, nline - 1);
          break;
        case OUTPLANE:
          file_outplane(words, nline - 1);
          break;
        case TORSION:
          file_torsion(words, nline - 1);
          break;
        case PITORSION:
          file_pitorsion(words, nline - 1);
          break;
        case ATOMMULT:
          file_multipole(words, nline - 1);
          break;
        case QPENETRATION:
          file_charge_penetration(words, nline - 1);
          break;
        case DIPPOLAR:
          file_dippolar(words, nline - 1);
          break;
        case QTRANSFER:
          file_charge_transfer(words, nline - 1);
          break;
        case UNKNOWN:
        case END_OF_FILE:
        default:
          ;    // do nothing
      }
    }

    if (n < 0) break;
  }

  if (me == 0) fclose(fptr);

  if (forcefield_flag == 0 || atomtype_flag == 0)
    error->all(FLERR, "Pair {} potential file {} incomplete", mystyle, filename);
}

/* ----------------------------------------------------------------------
   read optional KEY file of one-line settings
------------------------------------------------------------------------- */

void PairAmoeba::read_keyfile(char *filename)
{
  double aprd, bprd, cprd;

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
  char line[MAXLINE];
  if (me == 0) {
    fptr = utils::open_potential(filename, lmp, nullptr);
    if (fptr == nullptr)
      error->one(FLERR, "Cannot open {} key file {}: {}", utils::uppercase(mystyle), filename,
                 utils::getsyserror());
  }

  // read lines, one at a time

  int n;
  char *ptr;

  while (true) {
    if (me == 0) {
      ptr = fgets(line, MAXLINE, fptr);
      if (!ptr)
        n = -1;
      else
        n = strlen(line) + 1;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (n < 0) break;
    MPI_Bcast(line, n, MPI_CHAR, 0, world);

    // skip over empty or comment lines
    auto trimmed = utils::lowercase(utils::trim(line));
    if (trimmed.empty() || utils::strmatch(trimmed, "^#") || utils::strmatch(trimmed, "!!"))
      continue;

    const auto words = Tokenizer(trimmed).as_vector();
    const int nwords = words.size();
    const auto &keyword = words[0];

    if (utils::strmatch(keyword, "^[^a-z]+")) {
      ;    // ignore keywords that do not start with text
    } else if (keyword == "a-axis") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      aprd = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "b-axis") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      bprd = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "c-axis") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      cprd = utils::numeric(FLERR, words[1], false, lmp);

    } else if (keyword == "cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      double cut = utils::numeric(FLERR, words[1], false, lmp);
      vdwcut = repcut = dispcut = mpolecut = ctrncut = ewaldcut = dewaldcut = cut;
      vdwtaper = 0.9 * vdwcut;
      reptaper = 0.9 * repcut;
      disptaper = 0.9 * dispcut;
      mpoletaper = 0.65 * mpolecut;
      ctrntaper = 0.9 * ctrncut;
    } else if (keyword == "taper") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      double taper = utils::numeric(FLERR, words[1], false, lmp);
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
    } else if (keyword == "vdw-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      vdwcut = utils::numeric(FLERR, words[1], false, lmp);
      vdwtaper = 0.9 * vdwcut;
    } else if (keyword == "repulsion-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      repcut = utils::numeric(FLERR, words[1], false, lmp);
      reptaper = 0.9 * repcut;
    } else if (keyword == "dispersion-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      dispcut = utils::numeric(FLERR, words[1], false, lmp);
      disptaper = 0.9 * dispcut;
    } else if (keyword == "mpole-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      mpolecut = utils::numeric(FLERR, words[1], false, lmp);
      mpoletaper = 0.65 * mpolecut;
    } else if (keyword == "ctrn-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      ctrncut = utils::numeric(FLERR, words[1], false, lmp);
      ctrntaper = 0.9 * ctrncut;
    } else if (keyword == "ewald-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      ewaldcut = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "dewald-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      dewaldcut = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "usolve-cutoff") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      usolvcut = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "usolve-diag") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      udiag = utils::numeric(FLERR, words[1], false, lmp);

    } else if (keyword == "vdw-taper") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      vdwtaper = utils::numeric(FLERR, words[1], false, lmp);
      if (vdwtaper < 1.0) vdwtaper = -vdwtaper * vdwcut;
    } else if (keyword == "repulsion-taper") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      reptaper = utils::numeric(FLERR, words[1], false, lmp);
      if (reptaper < 1.0) reptaper = -reptaper * repcut;
    } else if (keyword == "dispersion-taper") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      disptaper = utils::numeric(FLERR, words[1], false, lmp);
      if (disptaper < 1.0) disptaper = -disptaper * vdwcut;
    } else if (keyword == "mpole-taper") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      mpoletaper = utils::numeric(FLERR, words[1], false, lmp);
      if (mpoletaper < 1.0) mpoletaper = -mpoletaper * vdwcut;
    } else if (keyword == "ctrn-taper") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      ctrntaper = utils::numeric(FLERR, words[1], false, lmp);
      if (ctrntaper < 1.0) ctrntaper = -ctrntaper * vdwcut;

    } else if (keyword == "delta-halgren") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      dhal = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "gamma-halgren") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      ghal = utils::numeric(FLERR, words[1], false, lmp);

    } else if (keyword == "ewald") {
      if (nwords != 1) error->all(FLERR, "AMOEBA keyfile line is invalid");
      use_ewald = 1;
    } else if (keyword == "dewald") {
      if (nwords != 1) error->all(FLERR, "AMOEBA keyfile line is invalid");
      use_dewald = 1;

    } else if (keyword == "pme-order") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      bseorder = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "ppme-order") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      bsporder = utils::numeric(FLERR, words[1], false, lmp);
    } else if (keyword == "dpme-order") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      bsdorder = utils::numeric(FLERR, words[1], false, lmp);

    } else if (keyword == "pme-grid") {
      if (nwords != 2 && nwords != 4) error->all(FLERR, "AMOEBA keyfile line is invalid");
      if (nwords == 2)
        nefft1 = nefft2 = nefft3 = utils::numeric(FLERR, words[1], false, lmp);
      else {
        nefft1 = utils::numeric(FLERR, words[1], false, lmp);
        nefft2 = utils::numeric(FLERR, words[2], false, lmp);
        nefft3 = utils::numeric(FLERR, words[3], false, lmp);
      }
      pmegrid_key = 1;
    } else if (keyword == "dpme-grid") {
      if (nwords != 2 && nwords != 4) error->all(FLERR, "AMOEBA keyfile line is invalid");
      if (nwords == 2)
        ndfft1 = ndfft2 = ndfft3 = utils::numeric(FLERR, words[1], false, lmp);
      else {
        ndfft1 = utils::numeric(FLERR, words[1], false, lmp);
        ndfft2 = utils::numeric(FLERR, words[2], false, lmp);
        ndfft3 = utils::numeric(FLERR, words[3], false, lmp);
      }
      dpmegrid_key = 1;

    } else if (keyword == "ewald-alpha") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      aeewald = utils::numeric(FLERR, words[1], false, lmp);
      aeewald_key = 1;
    } else if (keyword == "pewald-alpha") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      apewald = utils::numeric(FLERR, words[1], false, lmp);
      apewald_key = 1;
    } else if (keyword == "dewald-alpha") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      adewald = utils::numeric(FLERR, words[1], false, lmp);
      adewald_key = 1;

      // polarization options

    } else if (keyword == "polarization") {
      if (words[1] == "mutual")
        poltyp = MUTUAL;
      else if (utils::strmatch(words[1], "^opt")) {
        poltyp = OPT;
        if (words[1] == "opt")
          optorder = 4;
        else
          optorder = utils::inumeric(FLERR, &words[1][3], false, lmp);
        if (optorder < 1 || optorder > 6)
          error->all(FLERR, "Unrecognized polarization OPT{} in AMOEBA FF file", optorder);
      } else if (words[1] == "tcg")
        error->all(FLERR, "Polarization TCG not yet supported in AMOEBA/HIPPO");
      else if (words[1] == "direct")
        poltyp = DIRECT;
      else
        error->all(FLERR, "Unrecognized polarization in AMOEBA FF file");

    } else if (keyword == "polar-predict") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      if (words[1] == "gear") {
        polpred = GEAR;
        maxualt = 7;
      } else if (words[1] == "aspc") {
        polpred = ASPC;
        maxualt = 17;
      } else if (words[1] == "lsqr") {
        polpred = LSQR;
        maxualt = 7;
      } else
        error->all(FLERR, "AMOEBA keyfile line is invalid");
      use_pred = 1;
    } else if (keyword == "polar-iter") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      politer = utils::inumeric(FLERR, words[1], false, lmp);
    } else if (keyword == "polar-eps") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      poleps = utils::numeric(FLERR, words[1], false, lmp);

    } else if (keyword == "pcg-precond") {
      if (nwords != 1) error->all(FLERR, "AMOEBA keyfile line is invalid");
      pcgprec = 1;
    } else if (keyword == "pcg-noprecond") {
      if (nwords != 1) error->all(FLERR, "AMOEBA keyfile line is invalid");
      pcgprec = 0;
    } else if (keyword == "pcg-guess") {
      if (nwords != 1) error->all(FLERR, "AMOEBA keyfile line is invalid");
      pcgguess = 1;
    } else if (keyword == "pcg-noguess") {
      if (nwords != 1) error->all(FLERR, "AMOEBA keyfile line is invalid");
      pcgguess = 0;
    } else if (keyword == "pcg-peek") {
      if (nwords != 2) error->all(FLERR, "AMOEBA keyfile line is invalid");
      pcgpeek = utils::numeric(FLERR, words[1], false, lmp);

      // Tinker keywords that LAMMPS can skip

    } else if (keyword == "parameters") {
    } else if (keyword == "verbose") {
    } else if (keyword == "openmp-threads") {
    } else if (keyword == "digits") {
    } else if (keyword == "neighbor-list") {
    } else if (keyword == "tau-temperature") {
    } else if (keyword == "tau-pressure") {

      // error if LAMMPS does not recognize other keywords

    } else
      error->all(FLERR, "LAMMPS does not recognize AMOEBA keyfile keyword {}", keyword);
  }

  // close key file

  if (me == 0) fclose(fptr);

  // cutoff resets for long-range interactions

  if (use_ewald) mpolecut = ewaldcut;
  if (use_dewald) dispcut = dewaldcut;

  // error checks

  if (use_ewald || use_dewald) {
    if (domain->nonperiodic) error->all(FLERR, "AMOEBA KSpace requires fully periodic system");
  }

  if (aprd > 0.0 && (!domain->xperiodic || domain->xprd != aprd))
    error->all(FLERR, "AMOEBA abc prd does not match LAMMPS domain");
  if (bprd > 0.0 && (!domain->yperiodic || domain->yprd != bprd))
    error->all(FLERR, "AMOEBA abc prd does not match LAMMPS domain");
  if (cprd > 0.0 && (!domain->zperiodic || domain->zprd != cprd))
    error->all(FLERR, "AMOEBA abc prd does not match LAMMPS domain");
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_ffield(const std::vector<std::string> &words, int nline)
{
  if (words.size() < 2)
    error->all(FLERR, "Keyword {} without argument(s) in {} PRM file", words[0], mystyle);

  if (words[0] == "forcefield") {
    forcefield = utils::strdup(words[1]);
  } else if (words[0] == "bond-cubic") {
    bond_cubic = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "bond-quartic") {
    bond_quartic = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "angle-cubic") {
    angle_cubic = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "angle-quartic") {
    angle_quartic = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "angle-pentic") {
    angle_pentic = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "angle-sextic") {
    angle_sextic = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "opbendtype") {
    if (words[1] == "allinger") {
      opbendtype = ALLINGER;
    } else {
      error->all(FLERR, "Unrecognized opbendtype {} in {} PRM file", words[1], mystyle);
    }
  } else if (words[0] == "opbend-cubic") {
    opbend_cubic = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "opbend-quartic") {
    opbend_quartic = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "opbend-pentic") {
    opbend_pentic = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "opbend-sextic") {
    opbend_sextic = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "torsionunit") {
    torsion_unit = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "vdwtype") {
    if (words[1] == "buffered-14-7")
      vdwtype = BUFFERED_14_7;
    else
      error->all(FLERR, "Unrecognized vdwtype {} in {} PRM file", words[1], mystyle);
  } else if (words[0] == "radiusrule") {
    if (words[1] == "arithmetic")
      radius_rule = ARITHMETIC;
    else if (words[1] == "geometric")
      radius_rule = GEOMETRIC;
    else if (words[1] == "cubic-mean")
      radius_rule = CUBIC_MEAN;
    else
      error->all(FLERR, "Unrecognized radiusrule {} in {} PRM file", words[1], mystyle);
  } else if (words[0] == "radiustype") {
    if (words[1] == "r-min")
      radius_type = R_MIN;
    else if (words[1] == "sigma")
      radius_type = SIGMA;
    else
      error->all(FLERR, "Unrecognized radiustype {} in {} PRM file", words[1], mystyle);
  } else if (words[0] == "radiussize") {
    if (words[1] == "diameter")
      radius_size = DIAMETER;
    else
      error->all(FLERR, "Unrecognized radiussize {} in {} PRM file", words[1], mystyle);
  } else if (words[0] == "epsilonrule") {
    if (words[1] == "arithmetic")
      epsilon_rule = ARITHMETIC;
    else if (words[1] == "geometric")
      epsilon_rule = GEOMETRIC;
    else if (words[1] == "harmonic")
      epsilon_rule = HARMONIC;
    else if (words[1] == "hhg")
      epsilon_rule = HHG;
    else if (words[1] == "w-h")
      epsilon_rule = W_H;
    else
      error->all(FLERR, "Unrecognized epsilonrule {} in {} PRM file", words[1], mystyle);

  } else if (words[0] == "dielectric") {
    am_dielectric = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "polarization") {
    if (words[1] == "mutual")
      poltyp = MUTUAL;
    else if (utils::strmatch(words[1], "^opt")) {
      poltyp = OPT;
      if (words[1] == "opt")
        optorder = 4;
      else
        optorder = utils::inumeric(FLERR, words[1].c_str() + 3, false, lmp);
      if (optorder < 1 || optorder > 6)
        error->all(FLERR, "Unrecognized polarization {} in {} PRM file line {}", words[1], mystyle);
    } else if (words[1] == "tcg")
      error->all(FLERR, "Polarization TCG not yet supported in AMOEBA/HIPPO");
    else if (words[1] == "direct")
      poltyp = DIRECT;
    else
      error->all(FLERR, "Unrecognized polarization {} in {} PRM file", words[1], mystyle);

  } else if (words[0] == "vdw-12-scale") {
    special_hal[1] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "vdw-13-scale") {
    special_hal[2] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "vdw-14-scale") {
    special_hal[3] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "vdw-15-scale") {
    special_hal[4] = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "rep-12-scale") {
    special_repel[1] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "rep-13-scale") {
    special_repel[2] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "rep-14-scale") {
    special_repel[3] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "rep-15-scale") {
    special_repel[4] = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "disp-12-scale") {
    special_disp[1] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "disp-13-scale") {
    special_disp[2] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "disp-14-scale") {
    special_disp[3] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "disp-15-scale") {
    special_disp[4] = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "mpole-12-scale") {
    special_mpole[1] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "mpole-13-scale") {
    special_mpole[2] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "mpole-14-scale") {
    special_mpole[3] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "mpole-15-scale") {
    special_mpole[4] = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "polar-12-scale") {
    special_polar_pscale[1] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "polar-13-scale") {
    special_polar_pscale[2] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "polar-14-scale") {
    special_polar_pscale[3] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "polar-15-scale") {
    special_polar_pscale[4] = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "polar-12-intra") {
    special_polar_piscale[1] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "polar-13-intra") {
    special_polar_piscale[2] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "polar-14-intra") {
    special_polar_piscale[3] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "polar-15-intra") {
    special_polar_piscale[4] = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "induce-12-scale") {
    special_polar_wscale[1] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "induce-13-scale") {
    special_polar_wscale[2] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "induce-14-scale") {
    special_polar_wscale[3] = utils::numeric(FLERR, words[1], false, lmp);
  } else if (words[0] == "induce-15-scale") {
    special_polar_wscale[4] = utils::numeric(FLERR, words[1], false, lmp);

  } else if (words[0] == "direct-11-scale") {
    polar_dscale = utils::numeric(FLERR, words[1], false, lmp);
  } else if (utils::strmatch(words[0], "^direct-1[234]-scale$")) {
    double tmp = utils::numeric(FLERR, words[1], false, lmp);
    if (tmp != 1.0)
      error->all(FLERR, "{} FF direct-scale 1-2, 1-3, 1-4 values should be 1.0",
                 utils::uppercase(mystyle));
  } else if (words[0] == "mutual-11-scale") {
    polar_uscale = utils::numeric(FLERR, words[1], false, lmp);
  } else if (utils::strmatch(words[0], "^mutual-1[234]-scale$")) {
    double tmp = utils::numeric(FLERR, words[1], false, lmp);
    if (tmp != 1.0)
      error->all(FLERR, "{} FF mutual-scale 1-2, 1-3, 1-4 values should be 1.0",
                 utils::uppercase(mystyle));
    // error if LAMMPS does not recognize keyword
  } else {
    error->all(FLERR, "LAMMPS does not recognize {} PRM file setting on line {}: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));
  }
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_literature(const std::vector<std::string> & /*words*/, int /*nline*/)
{
  // do nothing, this section is skipped
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_atomtype(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "atom")
    error->all(FLERR, "{} PRM file atom type line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < 8)
    error->all(FLERR, "{} PRM file atom type line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int itype = utils::inumeric(FLERR, words[1], false, lmp);
  int iclass = utils::inumeric(FLERR, words[2], false, lmp);

  // grow per-type and per-class vecs/arrays as needed

  allocate_type_class(itype, iclass);
  n_amtype = MAX(n_amtype, itype);
  n_amclass = MAX(n_amclass, iclass);

  // store words from line

  amtype_defined[itype] = 1;
  amclass_defined[iclass] = 1;
  amtype2class[itype] = iclass;

  atomic_num[itype] = utils::inumeric(FLERR, words[words.size() - 3], false, lmp);
  am_mass[itype] = utils::numeric(FLERR, words[words.size() - 2], false, lmp);
  valence[itype] = utils::inumeric(FLERR, words[words.size() - 1], false, lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_vdwl(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "vdw")
    error->all(FLERR, "{} PRM file Van der Waals line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() != 4 && words.size() != 5)
    error->all(FLERR, "{} PRM file Vand der Walls line {} has incorrect length ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int iclass = utils::inumeric(FLERR, words[1], false, lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR, "{} RPM file Van der Waals type index {} on line {} is invalid: {}",
               utils::uppercase(mystyle), iclass, nline, utils::join_words(words, " "));

  vdwl_sigma[iclass] = utils::numeric(FLERR, words[2], false, lmp);
  vdwl_eps[iclass] = utils::numeric(FLERR, words[3], false, lmp);
  if (words.size() == 4)
    kred[iclass] = 0.0;
  else
    kred[iclass] = utils::numeric(FLERR, words[4], false, lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_vdwl_pair(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "vdwpr")
    error->all(FLERR, "{} PRM file Van der Waals pair line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() != 5)
    error->all(FLERR, "{} PRM file Van der Waals pair line {} has incorret length ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  if (nvdwl_pair == max_vdwl_pair) {
    max_vdwl_pair += DELTA_VDWL_PAIR;
    memory->grow(vdwl_class_pair, max_vdwl_pair, 2, "amoeba:vdwl_class_pair");
    memory->grow(vdwl_sigma_pair, max_vdwl_pair, "amoeba:vdwl_sigma_pair");
    memory->grow(vdwl_eps_pair, max_vdwl_pair, "amoeba:vdwl_eps_pair");
  }

  vdwl_class_pair[nvdwl_pair][0] = utils::inumeric(FLERR, words[1], false, lmp);
  vdwl_class_pair[nvdwl_pair][1] = utils::inumeric(FLERR, words[2], false, lmp);
  vdwl_sigma_pair[nvdwl_pair] = utils::numeric(FLERR, words[3], false, lmp);
  vdwl_eps_pair[nvdwl_pair] = utils::numeric(FLERR, words[4], false, lmp);
  nvdwl_pair++;
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_bstretch(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "bond")
    error->all(FLERR, "{} PRM file bond stretch line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < 5)
    error->all(FLERR, "{} PRM file bond stretch line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_sbend(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "strbnd")
    error->all(FLERR, "{} PRM file stretch-bend line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() != 6)
    error->all(FLERR, "{} PRM file stretch-bend line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_abend(const std::vector<std::string> &words, int nline)
{
  if (words.size() < 6)
    error->all(FLERR, "{} PRM file angle bending line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_pauli(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "repulsion")
    error->all(FLERR, "{} PRM file Pauli repulsion line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < 5)
    error->all(FLERR, "{} PRM file Pauli repulsion line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int itype = utils::inumeric(FLERR, words[1], false, lmp);
  if (itype < 1 || itype > n_amtype)
    error->all(FLERR, "{} PRM file Pauli repulsion type index {} on line {} is invalid: {}",
               utils::uppercase(mystyle), itype, nline, utils::join_words(words, " "));

  // negate the elepr setting

  sizpr[itype] = utils::numeric(FLERR, words[2], false, lmp);
  dmppr[itype] = utils::numeric(FLERR, words[3], false, lmp);
  elepr[itype] = -fabs(utils::numeric(FLERR, words[4], false, lmp));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_dispersion(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "dispersion")
    error->all(FLERR, "{} PRM file dispersion line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < 4)
    error->all(FLERR, "{} PRM file dispersion line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int iclass = utils::inumeric(FLERR, words[1], false, lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR, "{} PRM file dispersion class index {} on line {} is invalid: {}",
               utils::uppercase(mystyle), iclass, nline, utils::join_words(words, " "));

  csix[iclass] = utils::numeric(FLERR, words[2], false, lmp);
  adisp[iclass] = utils::numeric(FLERR, words[3], false, lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_ub(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "ureybrad")
    error->all(FLERR, "{} PRM file Urey-Bradley line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() != 6)
    error->all(FLERR, "{} PRM file Urey-Bradley line {} has incorrect length ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_outplane(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "opbend")
    error->all(FLERR, "{} PRM file out-of-plane bend line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() != 6)
    error->all(FLERR, "{} PRM file out-of-plane bend line {} has incorrect length ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_torsion(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "torsion")
    error->all(FLERR, "{} PRM file torsion line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() != 14)
    error->all(FLERR, "{} PRM file torsion line {} has incorrect length ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_pitorsion(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "pitors")
    error->all(FLERR, "{} PRM file pi-torsion line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() != 4)
    error->all(FLERR, "{} PRM file pi-torsion line {} has incorrect length ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_multipole(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "multipole")
    error->all(FLERR, "{} PRM file atomic multipole line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < 12 || words.size() > 15)
    error->all(FLERR, "{} PRM file atomic multipole line {} has incorrect length ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int itype = utils::inumeric(FLERR, words[1], false, lmp);
  if (itype < 1 || itype > n_amtype)
    error->all(FLERR, "{} PRM file atomic multipole type index {} on line {} is invalid: {}",
               utils::uppercase(mystyle), itype, nline, utils::join_words(words, " "));

  int iframe = nmultiframe[itype];
  if (iframe >= MAX_FRAME_PER_TYPE)
    error->all(FLERR, "{} MAX_FRAME_PER_TYPE is too small: {}", utils::uppercase(mystyle), iframe);

  int extra;
  if (words.size() == 12) {
    zpole[itype][iframe] = xpole[itype][iframe] = ypole[itype][iframe] = 0;
    extra = 2;
  } else if (words.size() == 13) {
    zpole[itype][iframe] = utils::inumeric(FLERR, words[2], false, lmp);
    xpole[itype][iframe] = ypole[itype][iframe] = 0;
    extra = 3;
  } else if (words.size() == 14) {
    zpole[itype][iframe] = utils::inumeric(FLERR, words[2], false, lmp);
    xpole[itype][iframe] = utils::inumeric(FLERR, words[3], false, lmp);
    ypole[itype][iframe] = 0;
    extra = 4;
  } else if (words.size() == 15) {
    zpole[itype][iframe] = utils::inumeric(FLERR, words[2], false, lmp);
    xpole[itype][iframe] = utils::inumeric(FLERR, words[3], false, lmp);
    ypole[itype][iframe] = utils::inumeric(FLERR, words[4], false, lmp);
    extra = 5;
  }

  for (int i = 0; i < 10; i++)
    fpole[itype][iframe][i] = utils::numeric(FLERR, words[extra + i], false, lmp);

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

  for (int i = 1; i < 4; i++) fpole[itype][iframe][i] *= BOHR;
  for (int i = 4; i < 13; i++) fpole[itype][iframe][i] *= BOHR * BOHR / 3.0;

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

void PairAmoeba::file_charge_penetration(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "chgpen")
    error->all(FLERR, "{} PRM file charge penetration line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < 4)
    error->all(FLERR, "{} PRM file charge penetration line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int iclass = utils::inumeric(FLERR, words[1], false, lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR, "{} PRM file charge penetration class index {} on line {} is invalid: {}",
               utils::uppercase(mystyle), iclass, nline, utils::join_words(words, " "));

  pcore[iclass] = fabs(utils::numeric(FLERR, words[2], false, lmp));
  palpha[iclass] = utils::numeric(FLERR, words[3], false, lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_dippolar(const std::vector<std::string> &words, int nline)
{
  const std::size_t ndipparams = amoeba ? 4 : 3;
  if (words[0] != "polarize")
    error->all(FLERR, "{} PRM file dipole polariability line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < ndipparams)
    error->all(FLERR, "{} PRM file dipole polarizability line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int itype = utils::inumeric(FLERR, words[1], false, lmp);
  if (itype < 1 || itype > n_amtype)
    error->all(FLERR, "{} PRM file dipole polarizability type index {} on line {} is invalid: {}",
               utils::uppercase(mystyle), itype, nline, utils::join_words(words, " "));

  polarity[itype] = utils::numeric(FLERR, words[2], false, lmp);
  pdamp[itype] = pow(polarity[itype], 1.0 / 6.0);
  if (amoeba) thole[itype] = utils::numeric(FLERR, words[3], false, lmp);

  // eventually AMOEBA+ files will set dirdamp

  dirdamp[itype] = 0.0;

  int ngroup = words.size() - ndipparams;
  if (ngroup > MAX_TYPE_PER_GROUP)
    error->all(FLERR, "{} MAX_TYPE_PER_GROUP is too small: {} vs {}", utils::uppercase(mystyle),
               MAX_TYPE_PER_GROUP, ngroup);

  npolgroup[itype] = ngroup;
  for (int igroup = 0; igroup < ngroup; igroup++)
    polgroup[itype][igroup] = utils::inumeric(FLERR, words[ndipparams + igroup], false, lmp);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::file_charge_transfer(const std::vector<std::string> &words, int nline)
{
  if (words[0] != "chgtrn")
    error->all(FLERR, "{} PRM file charge transfer line {} has invalid format: {}",
               utils::uppercase(mystyle), nline, utils::join_words(words, " "));

  if (words.size() < 4)
    error->all(FLERR, "{} PRM file charge transfer line {} has too few values ({}): {}",
               utils::uppercase(mystyle), nline, words.size(), utils::join_words(words, " "));

  int iclass = utils::inumeric(FLERR, words[1], false, lmp);
  if (iclass < 1 || iclass > n_amclass)
    error->all(FLERR, "{} PRM file charge transfer class index {} on line {} is invalid: {}",
               utils::uppercase(mystyle), iclass, nline, utils::join_words(words, " "));

  chgct[iclass] = utils::numeric(FLERR, words[2], false, lmp);
  dmpct[iclass] = utils::numeric(FLERR, words[3], false, lmp);
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

  amtype_defined = nullptr;
  amtype2class = nullptr;
  atomic_num = nullptr;
  valence = nullptr;
  am_mass = nullptr;
  am_q = nullptr;
  am_mu = nullptr;
  npolgroup = nullptr;
  polgroup = nullptr;
  polarity = nullptr;
  pdamp = nullptr;
  thole = nullptr;
  dirdamp = nullptr;
  sizpr = nullptr;
  dmppr = nullptr;
  elepr = nullptr;

  nmultiframe = nullptr;
  mpaxis = nullptr;
  xpole = nullptr;
  ypole = nullptr;
  zpole = nullptr;
  fpole = nullptr;

  // per class data

  amclass_defined = nullptr;
  vdwl_eps = nullptr;
  vdwl_sigma = nullptr;
  kred = nullptr;
  csix = nullptr;
  adisp = nullptr;
  chgct = nullptr;
  dmpct = nullptr;
  pcore = nullptr;
  palpha = nullptr;

  // other

  vdwl_class_pair = nullptr;
  vdwl_sigma_pair = nullptr;
  vdwl_eps_pair = nullptr;
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

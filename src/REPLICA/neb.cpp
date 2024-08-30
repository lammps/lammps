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

#include "neb.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "finish.h"
#include "fix.h"
#include "fix_neb.h"
#include "math_const.h"
#include "memory.h"
#include "min.h"
#include "modify.h"
#include "output.h"
#include "thermo.h"
#include "timer.h"
#include "tokenizer.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int MAXLINE = 256;
static constexpr int CHUNK = 1024;
static constexpr int ATTRIBUTE_PERLINE = 4;

enum { DEFAULT, TERSE, VERBOSE };

/* ---------------------------------------------------------------------- */

NEB::NEB(LAMMPS *lmp) : Command(lmp), fp(nullptr), all(nullptr), rdist(nullptr)
{
  print_mode = DEFAULT;

  // replica info

  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  me_universe = universe->me;
  uworld = universe->uworld;
  MPI_Comm_rank(world, &me);
}

/* ----------------------------------------------------------------------
   internal NEB constructor, called from TAD
------------------------------------------------------------------------- */

NEB::NEB(LAMMPS *lmp, double etol_in, double ftol_in, int n1steps_in, int n2steps_in, int nevery_in,
         double *buf_init, double *buf_final) :
    NEB(lmp)
{
  double delx, dely, delz;

  etol = etol_in;
  ftol = ftol_in;
  n1steps = n1steps_in;
  n2steps = n2steps_in;
  nevery = nevery_in;

  // generate linear interpolated replica
  double fraction = ireplica / (nreplica - 1.0);
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int ii = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = buf_final[ii] - buf_init[ii];
    dely = buf_final[ii + 1] - buf_init[ii + 1];
    delz = buf_final[ii + 2] - buf_init[ii + 2];
    domain->minimum_image(delx, dely, delz);
    x[i][0] = buf_init[ii] + fraction * delx;
    x[i][1] = buf_init[ii + 1] + fraction * dely;
    x[i][2] = buf_init[ii + 2] + fraction * delz;
    ii += 3;
  }
}

/* ---------------------------------------------------------------------- */

NEB::~NEB()
{
  MPI_Comm_free(&roots);
  memory->destroy(all);
  delete[] rdist;
  if (fp) {
    if (compressed)
      platform::pclose(fp);
    else
      fclose(fp);
  }
}

/* ----------------------------------------------------------------------
   perform NEB on multiple replicas
------------------------------------------------------------------------- */

void NEB::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->universe_all(FLERR, "NEB command before simulation box is defined");

  if (narg < 6) error->universe_all(FLERR, "Illegal NEB command: missing argument(s)");

  etol = utils::numeric(FLERR, arg[0], false, lmp);
  ftol = utils::numeric(FLERR, arg[1], false, lmp);
  n1steps = utils::inumeric(FLERR, arg[2], false, lmp);
  n2steps = utils::inumeric(FLERR, arg[3], false, lmp);
  nevery = utils::inumeric(FLERR, arg[4], false, lmp);

  // error checks

  if (etol < 0.0) error->universe_all(FLERR, fmt::format("Illegal NEB energy tolerance: {}", etol));
  if (ftol < 0.0) error->universe_all(FLERR, fmt::format("Illegal NEB force tolerance: {}", ftol));
  if (nevery <= 0)
    error->universe_all(FLERR, fmt::format("Illegal NEB command every parameter: {}", nevery));
  if (n1steps % nevery)
    error->universe_all(FLERR,
                        fmt::format("NEB N1 value {} incompatible with every {}", n1steps, nevery));
  if (n2steps % nevery)
    error->universe_all(FLERR,
                        fmt::format("NEB N2 value {} incompatible with every {}", n2steps, nevery));

  // error checks

  if (nreplica == 1) error->universe_all(FLERR, "Cannot use NEB with a single replica");
  if (atom->map_style == Atom::MAP_NONE)
    error->universe_all(FLERR, "Cannot use NEB without an atom map");

  // process file-style setting to setup initial configs for all replicas
  int iarg = 5;
  int filecmd = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "final") == 0) {
      if (iarg + 2 > narg)
        error->universe_all(FLERR, "Illegal NEB final command: missing arguments");
      inpfile = arg[iarg + 1];
      readfile(inpfile, 0);
      filecmd = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "each") == 0) {
      if (iarg + 2 > narg)
        error->universe_all(FLERR, "Illegal NEB each command: missing arguments");
      inpfile = arg[iarg + 1];
      readfile(inpfile, 1);
      filecmd = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "none") == 0) {
      filecmd = 1;
      ++iarg;
    } else if (strcmp(arg[iarg], "verbosity") == 0) {
      if (iarg + 2 > narg)
        error->universe_all(FLERR, "Illegal NEB verbosity command: missing arguments");
      if (strcmp(arg[iarg + 1], "verbose") == 0)
        print_mode = VERBOSE;
      else if (strcmp(arg[iarg + 1], "default") == 0)
        print_mode = DEFAULT;
      else if (strcmp(arg[iarg + 1], "terse") == 0)
        print_mode = TERSE;
      else
        error->universe_all(FLERR, fmt::format("Unknown NEB verbosity option {}", arg[iarg + 1]));
      iarg += 2;
    } else
      error->universe_all(FLERR, fmt::format("Unknown NEB command keyword: {}", arg[iarg]));
  }

  if (!filecmd) error->universe_all(FLERR, "NEB is missing 'final', 'each', or 'none' keyword");

  // run the NEB calculation

  run();
}

/* ----------------------------------------------------------------------
   run NEB on multiple replicas
------------------------------------------------------------------------- */

void NEB::run()
{
  // create MPI communicator for root proc from each world

  int color;
  if (me == 0)
    color = 0;
  else
    color = MPI_UNDEFINED;
  MPI_Comm_split(uworld, color, 0, &roots);

  auto fixes = modify->get_fix_by_style("^neb$");
  if (fixes.size() != 1)
    error->universe_all(FLERR, "NEB requires use of exactly one fix neb instance");

  fneb = dynamic_cast<FixNEB *>(fixes[0]);
  if (print_mode == VERBOSE)
    numall = 7;
  else
    numall = 4;
  memory->create(all, nreplica, numall, "neb:all");
  rdist = new double[nreplica];

  // initialize LAMMPS

  update->whichflag = 2;
  update->etol = etol;
  update->ftol = ftol;
  update->multireplica = 1;

  lmp->init();

  if (update->minimize->searchflag)
    error->universe_all(FLERR, "NEB requires a damped dynamics minimizer");

  // setup regular NEB minimization
  FILE *uscreen = universe->uscreen;
  FILE *ulogfile = universe->ulogfile;

  if (me_universe == 0 && uscreen) fprintf(uscreen, "Setting up regular NEB ...\n");

  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + n1steps;
  update->nsteps = n1steps;
  update->max_eval = n1steps;
  if (update->laststep < 0) error->universe_all(FLERR, "Too many timesteps for NEB");

  update->minimize->setup();

  if (me_universe == 0) {
    if (uscreen) {
      fmt::print(uscreen, "    Step     {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} ",
                 "MaxReplicaForce", "MaxAtomForce", "GradV0", "GradV1", "GradVc", "EBF", "EBR",
                 "RDT");

      if (print_mode != TERSE) {
        for (int i = 1; i <= nreplica; ++i)
          fmt::print(uscreen, "{:^14} {:^14} ", "RD" + std::to_string(i), "PE" + std::to_string(i));
      }

      if (print_mode == VERBOSE) {
        for (int i = 1; i <= nreplica; ++i) {
          auto idx = std::to_string(i);
          fmt::print(uscreen, "{:^12}{:^12}{:^12} {:^12} {:^12}{:^12} ", "pathangle" + idx,
                     "angletangrad" + idx, "anglegrad" + idx, "gradV" + idx, "RepForce" + idx,
                     "MaxAtomForce" + idx);
        }
      }
      fprintf(uscreen, "\n");
    }

    if (ulogfile) {
      fmt::print(ulogfile, "    Step     {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} ",
                 "MaxReplicaForce", "MaxAtomForce", "GradV0", "GradV1", "GradVc", "EBF", "EBR",
                 "RDT");

      if (print_mode != TERSE) {
        for (int i = 1; i <= nreplica; ++i)
          fmt::print(ulogfile, "{:^14} {:^14} ", "RD" + std::to_string(i),
                     "PE" + std::to_string(i));
      }

      if (print_mode == VERBOSE) {
        for (int i = 1; i <= nreplica; ++i) {
          auto idx = std::to_string(i);
          fmt::print(ulogfile, "{:^12}{:^12}{:^12} {:^12} {:^12}{:^12} ", "pathangle" + idx,
                     "angletangrad" + idx, "anglegrad" + idx, "gradV" + idx, "RepForce" + idx,
                     "MaxAtomForce" + idx);
        }
      }

      fprintf(ulogfile, "\n");
    }
  }
  print_status();

  // perform regular NEB for n1steps or until replicas converge
  // retrieve PE values from fix NEB and print every nevery iterations
  // break out of while loop early if converged
  // damped dynamic min styles ensure all replicas converge together

  timer->init();
  timer->barrier_start();

  while (update->minimize->niter < n1steps) {
    update->minimize->run(nevery);
    print_status();
    if (update->minimize->stop_condition) break;
  }

  timer->barrier_stop();

  update->minimize->cleanup();

  Finish finish(lmp);
  finish.end(1);

  // switch fix NEB to climbing mode
  // top = replica that becomes hill climber

  double vmax = all[0][0];
  int top = 0;
  for (int m = 1; m < nreplica; m++)
    if (vmax < all[m][0]) {
      vmax = all[m][0];
      top = m;
    }

  // setup climbing NEB minimization
  // must reinitialize minimizer so it re-creates its fix MINIMIZE

  if (me_universe == 0 && uscreen) fprintf(uscreen, "Setting up climbing ...\n");

  if (me_universe == 0) {
    if (uscreen) fprintf(uscreen, "Climbing replica = %d\n", top + 1);
    if (ulogfile) fprintf(ulogfile, "Climbing replica = %d\n", top + 1);
  }

  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + n2steps;
  update->nsteps = n2steps;
  update->max_eval = n2steps;
  if (update->laststep < 0) error->universe_all(FLERR, "Too many timesteps");

  update->minimize->init();
  fneb->rclimber = top;
  update->minimize->setup();

  if (me_universe == 0) {
    if (uscreen) {
      fmt::print(uscreen, "    Step     {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} ",
                 "MaxReplicaForce", "MaxAtomForce", "GradV0", "GradV1", "GradVc", "EBF", "EBR",
                 "RDT");

      if (print_mode != TERSE) {
        for (int i = 1; i <= nreplica; ++i)
          fmt::print(uscreen, "{:^14} {:^14} ", "RD" + std::to_string(i), "PE" + std::to_string(i));
      }

      if (print_mode == VERBOSE) {
        for (int i = 1; i <= nreplica; ++i) {
          auto idx = std::to_string(i);
          fmt::print(uscreen, "{:^12}{:^12}{:^12} {:^12} {:^12}{:^12} ", "pathangle" + idx,
                     "angletangrad" + idx, "anglegrad" + idx, "gradV" + idx, "RepForce" + idx,
                     "MaxAtomForce" + idx);
        }
      }
      fprintf(uscreen, "\n");
    }

    if (ulogfile) {
      fmt::print(ulogfile, "    Step     {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} ",
                 "MaxReplicaForce", "MaxAtomForce", "GradV0", "GradV1", "GradVc", "EBF", "EBR",
                 "RDT");

      if (print_mode != TERSE) {
        for (int i = 1; i <= nreplica; ++i)
          fmt::print(ulogfile, "{:^14} {:^14} ", "RD" + std::to_string(i),
                     "PE" + std::to_string(i));
      }

      if (print_mode == VERBOSE) {
        for (int i = 1; i <= nreplica; ++i) {
          auto idx = std::to_string(i);
          fmt::print(ulogfile, "{:^12}{:^12}{:^12} {:^12} {:^12}{:^12} ", "pathangle" + idx,
                     "angletangrad" + idx, "anglegrad" + idx, "gradV" + idx, "RepForce" + idx,
                     "MaxAtomForce" + idx);
        }
      }
      fprintf(ulogfile, "\n");
    }
  }
  print_status();

  // perform climbing NEB for n2steps or until replicas converge
  // retrieve PE values from fix NEB and print every nevery iterations
  // break induced if converged
  // damped dynamic min styles ensure all replicas converge together

  timer->init();
  timer->barrier_start();

  while (update->minimize->niter < n2steps) {
    update->minimize->run(nevery);
    print_status();
    if (update->minimize->stop_condition) break;
  }

  timer->barrier_stop();

  update->minimize->cleanup();

  finish.end(1);

  update->whichflag = 0;
  update->multireplica = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   read initial config atom coords from file
   flag = 0
   only first replica opens file and reads it
   first replica bcasts lines to all replicas
   final replica stores coords
   intermediate replicas interpolate from coords
   new coord = replica fraction between current and final state
   initial replica does nothing
   flag = 1
   each replica (except first) opens file and reads it
   each replica stores coords
   initial replica does nothing
------------------------------------------------------------------------- */

void NEB::readfile(char *file, int flag)
{
  int i, nchunk, eofflag, nlines;
  tagint tag;
  char *eof, *start, *next, *buf;
  char line[MAXLINE] = {'\0'};
  double delx, dely, delz;

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen, "Reading NEB coordinate file(s) ...\n");

  // flag = 0, universe root reads header of file, bcast to universe
  // flag = 1, each replica's root reads header of file, bcast to world
  //   but explicitly skip first replica

  if (flag == 0) {
    if (me_universe == 0) {
      open(file);
      while (true) {
        eof = fgets(line, MAXLINE, fp);
        if (eof == nullptr) error->one(FLERR, "Unexpected end of NEB file");
        start = &line[strspn(line, " \t\n\v\f\r")];
        if (*start != '\0' && *start != '#') break;
      }
      int rv = sscanf(line, "%d", &nlines);
      if (rv != 1) nlines = -1;
    }
    MPI_Bcast(&nlines, 1, MPI_INT, 0, uworld);
    if (nlines < 0) error->universe_all(FLERR, "Incorrectly formatted NEB file");
  } else {
    if (me == 0) {
      if (ireplica) {
        open(file);
        while (true) {
          eof = fgets(line, MAXLINE, fp);
          if (eof == nullptr) error->one(FLERR, "Unexpected end of NEB file");
          start = &line[strspn(line, " \t\n\v\f\r")];
          if (*start != '\0' && *start != '#') break;
        }
        int rv = sscanf(line, "%d", &nlines);
        if (rv != 1) nlines = -1;
      } else
        nlines = 0;
    }
    MPI_Bcast(&nlines, 1, MPI_INT, 0, world);
    if (nlines < 0) error->universe_all(FLERR, "Incorrectly formatted NEB file");
  }

  auto buffer = new char[CHUNK * MAXLINE];
  double fraction = ireplica / (nreplica - 1.0);
  double **x = atom->x;
  int nlocal = atom->nlocal;

  // loop over chunks of lines read from file
  // two versions of read_lines_from_file() for world vs universe bcast
  // count # of atom coords changed so can check for invalid atom IDs in file

  int ncount = 0, nread = 0;
  while (nread < nlines) {
    nchunk = MIN(nlines - nread, CHUNK);
    if (flag == 0)
      eofflag =
          utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, universe->me, universe->uworld);
    else
      eofflag = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eofflag) error->all(FLERR, "Unexpected end of NEB file");

    buf = buffer;
    next = strchr(buf, '\n');
    *next = '\0';
    int nwords = utils::count_words(utils::trim_comment(buf));
    *next = '\n';

    if (nwords != ATTRIBUTE_PERLINE) error->all(FLERR, "Incorrect atom format in NEB file");

    // loop over lines of atom coords
    // tokenize the line into values

    for (i = 0; i < nchunk; i++) {
      next = strchr(buf, '\n');
      *next = '\0';

      try {
        ValueTokenizer values(buf, " \t\n\r\f");

        // adjust atom coord based on replica fraction
        // for flag = 0, interpolate for intermediate and final replicas
        // for flag = 1, replace existing coord with new coord
        // ignore image flags of replica x
        // displacement from first replica is via minimum image convention
        // if x of some replica is across periodic boundary:
        //   new x may be outside box
        //   will be remapped back into box when simulation starts
        //   its image flags will then be adjusted

        tag = values.next_tagint();
        int m = atom->map(tag);
        if (m >= 0 && m < nlocal) {
          ncount++;

          delx = values.next_double() - x[m][0];
          dely = values.next_double() - x[m][1];
          delz = values.next_double() - x[m][2];

          domain->minimum_image(delx, dely, delz);

          if (flag == 0) {
            x[m][0] += fraction * delx;
            x[m][1] += fraction * dely;
            x[m][2] += fraction * delz;
          } else {
            x[m][0] += delx;
            x[m][1] += dely;
            x[m][2] += delz;
          }
        }
      } catch (std::exception &e) {
        error->universe_one(FLERR, "Incorrectly formatted NEB file: " + std::string(e.what()));
      }
      buf = next + 1;
    }
    nread += nchunk;
  }

  // check that all atom IDs in file were found by a proc

  if (flag == 0) {
    int ntotal;
    MPI_Allreduce(&ncount, &ntotal, 1, MPI_INT, MPI_SUM, uworld);
    if (ntotal != nreplica * nlines) error->universe_all(FLERR, "Invalid atom IDs in NEB file");
  } else {
    int ntotal;
    MPI_Allreduce(&ncount, &ntotal, 1, MPI_INT, MPI_SUM, world);
    if (ntotal != nlines) error->all(FLERR, "Invalid atom IDs in NEB file");
  }

  // clean up
  delete[] buffer;

  if (flag == 0) {
    if (me_universe == 0) {
      if (compressed)
        platform::pclose(fp);
      else
        fclose(fp);
    }
  } else {
    if (me == 0 && ireplica) {
      if (compressed)
        platform::pclose(fp);
      else
        fclose(fp);
    }
  }
  fp = nullptr;
}

/* ----------------------------------------------------------------------
   universe proc 0 opens NEB data file
   test if compressed
------------------------------------------------------------------------- */

void NEB::open(char *file)
{
  compressed = 0;
  if (platform::has_compress_extension(file)) {
    compressed = 1;
    fp = platform::compressed_read(file);
    if (!fp) error->one(FLERR, "Cannot open compressed file {}: {}", file, utils::getsyserror());
  } else
    fp = fopen(file, "r");

  if (fp == nullptr) error->one(FLERR, "Cannot open file {}: {}", file, utils::getsyserror());
}

/* ----------------------------------------------------------------------
   query fix NEB for info on each replica
   universe proc 0 prints current NEB status
------------------------------------------------------------------------- */

void NEB::print_status()
{
  double fnorm2 = sqrt(update->minimize->fnorm_sqr());
  double fnorminf = update->minimize->fnorm_inf();
  double fmaxreplica = 0.0;
  double fmaxatom = 0.0;

  if (me == 0) {
    MPI_Allreduce(&fnorm2, &fmaxreplica, 1, MPI_DOUBLE, MPI_MAX, roots);
    MPI_Allreduce(&fnorminf, &fmaxatom, 1, MPI_DOUBLE, MPI_MAX, roots);

    if (print_mode == VERBOSE) {
      freplica = new double[nreplica];
      MPI_Allgather(&fnorm2, 1, MPI_DOUBLE, &freplica[0], 1, MPI_DOUBLE, roots);
      fmaxatomInRepl = new double[nreplica];
      MPI_Allgather(&fnorminf, 1, MPI_DOUBLE, &fmaxatomInRepl[0], 1, MPI_DOUBLE, roots);
    }
  }

  double one[7];
  one[0] = fneb->veng;
  one[1] = fneb->plen;
  one[2] = fneb->nlen;
  one[3] = fneb->gradlen;

  if (print_mode == VERBOSE) {
    one[4] = fneb->dotpath;
    one[5] = fneb->dottangrad;
    one[6] = fneb->dotgrad;
  }

  if (output->thermo->normflag) one[0] /= atom->natoms;
  if (me == 0) MPI_Allgather(one, numall, MPI_DOUBLE, &all[0][0], numall, MPI_DOUBLE, roots);
  MPI_Bcast(&all[0][0], numall * nreplica, MPI_DOUBLE, 0, world);

  rdist[0] = 0.0;
  for (int i = 1; i < nreplica; i++) rdist[i] = rdist[i - 1] + all[i][1];
  double endpt = rdist[nreplica - 1] = rdist[nreplica - 2] + all[nreplica - 2][2];
  for (int i = 1; i < nreplica; i++) rdist[i] /= endpt;

  // look up GradV for the initial, final, and climbing replicas
  // these are identical to fnorm2, but to be safe we
  // take them straight from fix_neb

  double gradvnorm0, gradvnorm1, gradvnormc;

  int irep;
  irep = 0;
  gradvnorm0 = all[irep][3];
  irep = nreplica - 1;
  gradvnorm1 = all[irep][3];
  irep = fneb->rclimber;
  if (irep > -1) {
    gradvnormc = all[irep][3];
    ebf = all[irep][0] - all[0][0];
    ebr = all[irep][0] - all[nreplica - 1][0];
  } else {
    double vmax = all[0][0];
    int top = 0;
    for (int m = 1; m < nreplica; m++)
      if (vmax < all[m][0]) {
        vmax = all[m][0];
        top = m;
      }
    irep = top;
    gradvnormc = all[irep][3];
    ebf = all[irep][0] - all[0][0];
    ebr = all[irep][0] - all[nreplica - 1][0];
  }

  if (me_universe == 0) {
    constexpr double todeg = 180.0 / MY_PI;
    std::string mesg =
        fmt::format("{:10}   {:<14.8g}   {:<14.8g} ", update->ntimestep, fmaxreplica, fmaxatom);
    mesg += fmt::format("{:<14.8g} {:<14.8g} {:<14.8g} ", gradvnorm0, gradvnorm1, gradvnormc);
    mesg += fmt::format("{:<14.8g} {:<14.8g} {:<14.8g} ", ebf, ebr, endpt);
    if (print_mode != TERSE) {
      for (int i = 0; i < nreplica; i++)
        mesg += fmt::format("{:<14.8g} {:<14.8g} ", rdist[i], all[i][0]);
    }
    if (print_mode == VERBOSE) {
      mesg += fmt::format("{:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g}", NAN,
                          180 - acos(all[0][5]) * todeg, 180 - acos(all[0][6]) * todeg, all[0][3],
                          freplica[0], fmaxatomInRepl[0]);
      for (int i = 1; i < nreplica - 1; i++)
        mesg +=
            fmt::format("{:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g}",
                        180 - acos(all[i][4]) * todeg, 180 - acos(all[i][5]) * todeg,
                        180 - acos(all[i][6]) * todeg, all[i][3], freplica[i], fmaxatomInRepl[i]);
      mesg += fmt::format("{:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g} {:<12.5g}", NAN,
                          180 - acos(all[nreplica - 1][5]) * todeg, NAN, all[nreplica - 1][3],
                          freplica[nreplica - 1], fmaxatomInRepl[nreplica - 1]);
    }
    mesg += "\n";

    if (universe->uscreen) fputs(mesg.c_str(), universe->uscreen);
    if (universe->ulogfile) {
      fputs(mesg.c_str(), universe->ulogfile);
      fflush(universe->ulogfile);
    }
  }
  if ((me == 0) && (print_mode == VERBOSE)) {
    delete[] freplica;
    delete[] fmaxatomInRepl;
  }
}

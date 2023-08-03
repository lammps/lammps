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

// lmptype.h must be first b/c this file uses MAXBIGINT and includes mpi.h
// due to OpenMPI bug which sets INT64_MAX via its mpi.h
//   before lmptype.h can set flags to ensure it is done correctly

#include "thermo.h"

#include "angle.h"
#include "arg_info.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "input.h"
#include "kspace.h"
#include "lattice.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "timer.h"
#include "tokenizer.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

using namespace LAMMPS_NS;
using namespace MathConst;

// CUSTOMIZATION: add a new keyword by adding it to this list:

// step, elapsed, elaplong, dt, time, cpu, tpcpu, spcpu, cpuremain,
// part, timeremain
// atoms, temp, press, pe, ke, etotal
// evdwl, ecoul, epair, ebond, eangle, edihed, eimp, emol, elong, etail
// enthalpy, ecouple, econserve
// vol, density, lx, ly, lz, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
// xlat, ylat, zlat
// bonds, angles, dihedrals, impropers
// pxx, pyy, pzz, pxy, pxz, pyz
// fmax, fnorm, nbuild, ndanger
// cella, cellb, cellc, cellalpha, cellbeta, cellgamma

// CUSTOMIZATION: add a new thermo style by adding a constant to the enumerator,
// define a new string constant with the keywords and provide default formats.

enum { ONELINE, MULTILINE, YAMLLINE };
// style "one"
static constexpr char ONE[] = "step temp epair emol etotal press";
#define FORMAT_FLOAT_ONE_DEFAULT "% -14.8g"
#define FORMAT_INT_ONE_DEFAULT "%10d "
// style "multi"
static constexpr char MULTI[] =
    "etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press";
#define FORMAT_FLOAT_MULTI_DEFAULT "%14.4f"
#define FORMAT_INT_MULTI_DEFAULT "%14d"
// style "yaml"
static constexpr char YAML[] = "step temp ke pe ebond eangle edihed eimp evdwl ecoul elong press";
#define FORMAT_FLOAT_YAML_DEFAULT "%.15g"
#define FORMAT_INT_YAML_DEFAULT "%d"

#define FORMAT_MULTI_HEADER "------------ Step {:14} ----- CPU = {:12.7g} (sec) -------------"

enum { SCALAR, VECTOR, ARRAY };

static constexpr char id_temp[] = "thermo_temp";
static constexpr char id_press[] = "thermo_press";
static constexpr char id_pe[] = "thermo_pe";

static char fmtbuf[512];
#define DELTA 8

/* ---------------------------------------------------------------------- */

Thermo::Thermo(LAMMPS *_lmp, int narg, char **arg) :
    Pointers(_lmp), style(nullptr), vtype(nullptr), field2index(nullptr), argindex1(nullptr),
    argindex2(nullptr), temperature(nullptr), pressure(nullptr), pe(nullptr)
{
  style = utils::strdup(arg[0]);

  // set thermo_modify defaults

  lineflag = ONELINE;
  modified = 0;
  normuserflag = 0;
  lostflag = lostbond = Thermo::ERROR;
  lostbefore = warnbefore = 0;
  flushflag = 0;
  ntimestep = -1;

  // set style and corresponding lineflag
  // custom style builds its own line of keywords, including wildcard expansion
  // CUSTOMIZATION: add a new thermo style by adding it to the if statement
  // set line string with default keywords if not custom style.

  if (strcmp(style, "one") == 0) {
    line = ONE;
    lineflag = ONELINE;
  } else if (strcmp(style, "multi") == 0) {
    line = MULTI;
    lineflag = MULTILINE;
  } else if (strcmp(style, "yaml") == 0) {
    line = YAML;
    lineflag = YAMLLINE;

  } else if (strcmp(style, "custom") == 0) {
    if (narg == 1) error->all(FLERR, "Illegal thermo style custom command");

    // expand args if any have wildcard character "*"

    int expand = 0;
    char **earg;
    int nvalues = utils::expand_args(FLERR, narg - 1, &arg[1], 0, earg, lmp);
    if (earg != &arg[1]) expand = 1;

    line.clear();
    for (int iarg = 0; iarg < nvalues; iarg++) {
      line += earg[iarg];
      line += ' ';
    }

    // if wildcard expansion occurred, free earg memory from exapnd_args()

    if (expand) {
      for (int i = 0; i < nvalues; i++) delete[] earg[i];
      memory->sfree(earg);
    }

  } else
    error->all(FLERR, "Illegal thermo style command");

  index_temp = index_press_scalar = index_press_vector = index_pe = -1;

  // count fields in line
  // allocate per-field memory
  // process line of keywords

  nfield_initial = utils::trim_and_count_words(line);
  allocate();
  parse_fields(line);
}

/* ---------------------------------------------------------------------- */

Thermo::~Thermo()
{
  delete[] style;
  deallocate();
}

/* ---------------------------------------------------------------------- */

void Thermo::init()
{
  // set normvalue to default setting unless user has specified it

  if (normuserflag)
    normvalue = normuser;
  else if (strcmp(update->unit_style, "lj") == 0)
    normvalue = 1;
  else
    normvalue = 0;

  // add Volume field if volume changes and not style = custom
  // this check must come after domain init, so box_change is set

  nfield = nfield_initial;
  if (domain->box_change && strcmp(style, "custom") != 0)
    addfield("Volume", &Thermo::compute_vol, FLOAT);

  // set format string for each field
  // include keyword if lineflag = MULTILINE
  // add '/n' every 3 values if lineflag = MULTILINE
  // add trailing '/n' to last value
  // add YAML list item prefix for lineflag = YAMLLINE

  ValueTokenizer *format_line = nullptr;
  if (format_line_user.size()) format_line = new ValueTokenizer(format_line_user);

  field_data.clear();
  field_data.resize(nfield);
  std::string format_this, format_line_user_def;
  for (int i = 0; i < nfield; i++) {

    format[i].clear();
    format_this.clear();
    format_line_user_def.clear();

    if (vtype[i] == FLOAT) {
      field_data[i] = (double) 0.0;
    } else if (vtype[i] == INT) {
      field_data[i] = (int) 0;
    } else if (vtype[i] == BIGINT) {
      field_data[i] = (bigint) 0;
    }

    if ((lineflag == MULTILINE) && ((i % 3) == 0)) format[i] += "\n";
    if ((lineflag == YAMLLINE) && (i == 0)) format[i] += "  - [";
    if (format_line) format_line_user_def = format_line->next_string();

    if (format_column_user[i].size())
      format_this = format_column_user[i];
    else if (vtype[i] == FLOAT) {
      if (format_float_user.size())
        format_this = format_float_user;
      else if (format_line_user_def.size())
        format_this = format_line_user_def;
      else if (lineflag == ONELINE)
        format_this = FORMAT_FLOAT_ONE_DEFAULT;
      else if (lineflag == MULTILINE)
        format_this = FORMAT_FLOAT_MULTI_DEFAULT;
      else if (lineflag == YAMLLINE)
        format_this = FORMAT_FLOAT_YAML_DEFAULT;
    } else if ((vtype[i] == INT) || (vtype[i] == BIGINT)) {
      if (format_int_user.size()) {
        if (vtype[i] == INT)
          format_this = format_int_user;
        else
          format_this = format_bigint_user;
      } else if (format_line_user_def.size()) {
        format_this = format_line_user_def;
      } else {
        if (lineflag == ONELINE)
          format_this = FORMAT_INT_ONE_DEFAULT;
        else if (lineflag == MULTILINE)
          format_this = FORMAT_INT_MULTI_DEFAULT;
        else
          format_this = FORMAT_INT_YAML_DEFAULT;
        if (vtype[i] == BIGINT) {
          // replace "d" in int format with bigint format specifier
          auto found = format_this.find('%');
          found = format_this.find('d', found);
          format_this = format_this.replace(found, 1, std::string(BIGINT_FORMAT).substr(1));
        }
      }
    }

    if (lineflag == ONELINE)
      format[i] += format_this + " ";
    else if (lineflag == YAMLLINE)
      format[i] += format_this + ", ";
    else {
      if (keyword_user[i].size())
        format[i] += fmt::format("{:<8} = {} ", keyword_user[i], format_this);
      else
        format[i] += fmt::format("{:<8} = {} ", keyword[i], format_this);
    }
  }

  // chop off trailing blank or add closing bracket if needed and then add newline
  if (lineflag == ONELINE)
    format[nfield - 1].resize(format[nfield - 1].size() - 1);
  else if (lineflag == MULTILINE)
    format[nfield - 1].resize(format[nfield - 1].size() - 1);
  else if (lineflag == YAMLLINE)
    format[nfield - 1] += ']';
  format[nfield - 1] += '\n';

  delete format_line;

  // find current ptr for each Compute ID

  for (int i = 0; i < ncompute; i++) {
    computes[i] = modify->get_compute_by_id(id_compute[i]);
    if (!computes[i]) error->all(FLERR, "Could not find thermo compute with ID {}", id_compute[i]);
  }

  // find current ptr for each Fix ID
  // check that fix frequency is acceptable with thermo output frequency

  for (int i = 0; i < nfix; i++) {
    fixes[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fixes[i]) error->all(FLERR, "Could not find thermo fix ID {}", id_fix[i]);

    if (output->thermo_every % fixes[i]->global_freq)
      error->all(FLERR, "Thermo and fix {} not computed at compatible times", id_fix[i]);
  }

  // find current ptr for each Variable ID

  for (int i = 0; i < nvariable; i++) {
    variables[i] = input->variable->find(id_variable[i]);
    if (variables[i] < 0) error->all(FLERR, "Could not find thermo variable {}", id_variable[i]);
  }

  // set ptrs to keyword-specific Compute objects

  if (index_temp >= 0) temperature = computes[index_temp];
  if (index_press_scalar >= 0) pressure = computes[index_press_scalar];
  if (index_press_vector >= 0) pressure = computes[index_press_vector];
  if (index_pe >= 0) pe = computes[index_pe];
}

/* ----------------------------------------------------------------------
   Print thermo style header text
   if thermo style "multi":
      nothing

   if thermo style "one":
      <keyword1> <keyword2> <keyword3> ...
   each column head is centered with default formatting
   width should match the width of the default format string

   if thermo style "yaml":
      ---
      keywords: [ <keyword1>, <keyword2>, <keyword3> ...]
      data:

   ---------------------------------------------------------------------- */

void Thermo::header()
{
  if (lineflag == MULTILINE) return;

  std::string hdr;
  if (lineflag == YAMLLINE) hdr = "---\nkeywords: [";
  for (int i = 0; i < nfield; i++) {
    auto head = keyword[i];
    if (keyword_user[i].size()) head = keyword_user[i];
    if (lineflag == ONELINE) {
      if (vtype[i] == FLOAT)
        hdr += fmt::format("{:^14} ", head);
      else if ((vtype[i] == INT) || (vtype[i] == BIGINT))
        hdr += fmt::format("{:^11} ", head);
    } else if (lineflag == YAMLLINE)
      hdr += fmt::format("'{}', ", head);
  }
  if (lineflag == YAMLLINE)
    hdr += "]\ndata:";
  else
    hdr.resize(hdr.size() - 1);    // chop off trailing blank

  if (comm->me == 0) utils::logmesg(lmp, hdr + "\n");
}

/* ---------------------------------------------------------------------- */

void Thermo::footer()
{
  if (lineflag == YAMLLINE) utils::logmesg(lmp, "...\n");
}

/* ---------------------------------------------------------------------- */

void Thermo::compute(int flag)
{
  int i;

  firststep = flag;
  ntimestep = update->ntimestep;

  // check for lost atoms
  // turn off normflag if natoms = 0 to avoid divide by 0

  natoms = atom->natoms = lost_check();
  if (natoms == 0)
    normflag = 0;
  else
    normflag = normvalue;

  // invoke Compute methods needed for thermo keywords

  for (i = 0; i < ncompute; i++)
    if (compute_which[i] == SCALAR) {
      if (!(computes[i]->invoked_flag & Compute::INVOKED_SCALAR)) {
        computes[i]->compute_scalar();
        computes[i]->invoked_flag |= Compute::INVOKED_SCALAR;
      }
    } else if (compute_which[i] == VECTOR) {
      if (!(computes[i]->invoked_flag & Compute::INVOKED_VECTOR)) {
        computes[i]->compute_vector();
        computes[i]->invoked_flag |= Compute::INVOKED_VECTOR;
      }
    } else if (compute_which[i] == ARRAY) {
      if (!(computes[i]->invoked_flag & Compute::INVOKED_ARRAY)) {
        computes[i]->compute_array();
        computes[i]->invoked_flag |= Compute::INVOKED_ARRAY;
      }
    }

  // if lineflag = MULTILINE, prepend step/cpu header line

  line.clear();
  if (lineflag == MULTILINE) {
    double cpu;
    if (flag)
      cpu = timer->elapsed(Timer::TOTAL);
    else
      cpu = 0.0;
    line += fmt::format(FORMAT_MULTI_HEADER, ntimestep, cpu);
  }

  // add each thermo value to line with its specific format
  field_data.clear();
  field_data.resize(nfield);

  for (ifield = 0; ifield < nfield; ifield++) {
    (this->*vfunc[ifield])();
    if (vtype[ifield] == FLOAT) {
      snprintf(fmtbuf, sizeof(fmtbuf), format[ifield].c_str(), dvalue);
      line += fmtbuf;
      field_data[ifield] = dvalue;
    } else if (vtype[ifield] == INT) {
      snprintf(fmtbuf, sizeof(fmtbuf), format[ifield].c_str(), ivalue);
      line += fmtbuf;
      field_data[ifield] = ivalue;
    } else if (vtype[ifield] == BIGINT) {
      snprintf(fmtbuf, sizeof(fmtbuf), format[ifield].c_str(), bivalue);
      line += fmtbuf;
      field_data[ifield] = bivalue;
    }
  }

  // print line to screen and logfile

  if (comm->me == 0) {
    utils::logmesg(lmp, line);
    if (flushflag) utils::flush_buffers(lmp);
  }

  // set to 1, so that subsequent invocations of CPU time will be non-zero
  // e.g. via variables in print command

  firststep = 1;
}

/* ----------------------------------------------------------------------
   check for lost atoms, return current number of atoms
   also could number of warnings across MPI ranks and update total
------------------------------------------------------------------------- */

bigint Thermo::lost_check()
{
  // ntotal = current # of atoms, and Error class warnings

  bigint nlocal[2], ntotal[2] = {0, 0};
  nlocal[0] = atom->nlocal;
  nlocal[1] = error->get_numwarn();
  MPI_Allreduce(nlocal, ntotal, 2, MPI_LMP_BIGINT, MPI_SUM, world);
  if (ntotal[0] < 0) error->all(FLERR, "Too many total atoms");

  // print notification, if future warnings will be ignored
  bigint maxwarn = error->get_maxwarn();
  if ((maxwarn > 0) && (warnbefore == 0) && (ntotal[1] > maxwarn)) {
    warnbefore = 1;
    if (comm->me == 0)
      error->message(FLERR,
                     "WARNING: Too many warnings: {} vs {}. All future warnings will be suppressed",
                     ntotal[1], maxwarn);
  }
  error->set_allwarn(MIN(MAXSMALLINT, ntotal[1]));

  // no lost atoms, nothing else to do.
  if (ntotal[0] == atom->natoms) return ntotal[0];

  // if not checking or already warned, just return
  if (lostflag == Thermo::IGNORE) return ntotal[0];
  if (lostflag == Thermo::WARN && lostbefore == 1) { return ntotal[0]; }

  // error message

  if (lostflag == Thermo::ERROR)
    error->all(FLERR, "Lost atoms: original {} current {}", atom->natoms, ntotal[0]);

  // warning message

  if (comm->me == 0)
    error->warning(FLERR, "Lost atoms: original {} current {}", atom->natoms, ntotal[0]);

  // reset total atom count

  atom->natoms = ntotal[0];
  lostbefore = 1;
  return ntotal[0];
}

/* ----------------------------------------------------------------------
   modify thermo parameters
------------------------------------------------------------------------- */

void Thermo::modify_params(int narg, char **arg)
{
  if (narg == 0) utils::missing_cmd_args(FLERR, "thermo_modify", error);

  modified = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "temp") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify temp", error);
      if (index_temp < 0) error->all(FLERR, "Thermo style does not use temp");
      delete[] id_compute[index_temp];
      id_compute[index_temp] = utils::strdup(arg[iarg + 1]);

      temperature = modify->get_compute_by_id(arg[iarg + 1]);
      if (!temperature)
        error->all(FLERR, "Could not find thermo_modify temperature compute {}", arg[iarg + 1]);

      if (temperature->tempflag == 0)
        error->all(FLERR, "Thermo_modify compute {} does not compute temperature", arg[iarg + 1]);
      if (temperature->igroup != 0 && comm->me == 0)
        error->warning(FLERR, "Temperature for thermo pressure is not for group all");

      // reset id_temp of pressure to new temperature ID
      // either pressure currently being used by thermo or "thermo_press"

      Compute *pcompute;
      if (index_press_scalar >= 0) {
        pcompute = modify->get_compute_by_id(id_compute[index_press_scalar]);
        if (!pcompute)
          error->all(FLERR, "Pressure compute {} for thermo output does not exist",
                     id_compute[index_press_scalar]);
      } else if (index_press_vector >= 0) {
        pcompute = modify->get_compute_by_id(id_compute[index_press_vector]);
        if (!pcompute)
          error->all(FLERR, "Pressure compute {} for thermo output does not exist",
                     id_compute[index_press_vector]);
      } else
        pcompute = modify->get_compute_by_id("thermo_press");

      pcompute->reset_extra_compute_fix(arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "press") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify press", error);
      if (index_press_scalar < 0 && index_press_vector < 0)
        error->all(FLERR, "Thermo style does not use press");

      if (index_press_scalar >= 0) {
        delete[] id_compute[index_press_scalar];
        id_compute[index_press_scalar] = utils::strdup(arg[iarg + 1]);
      }
      if (index_press_vector >= 0) {
        delete[] id_compute[index_press_vector];
        id_compute[index_press_vector] = utils::strdup(arg[iarg + 1]);
      }

      pressure = modify->get_compute_by_id(arg[iarg + 1]);
      if (!pressure)
        error->all(FLERR, "Could not find thermo_modify pressure compute {}", arg[iarg + 1]);

      if (pressure->pressflag == 0)
        error->all(FLERR, "Thermo_modify compute {} does not compute pressure", arg[iarg + 1]);

      iarg += 2;

    } else if (strcmp(arg[iarg], "lost") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify lost", error);
      if (strcmp(arg[iarg + 1], "ignore") == 0)
        lostflag = Thermo::IGNORE;
      else if (strcmp(arg[iarg + 1], "warn") == 0)
        lostflag = Thermo::WARN;
      else if (strcmp(arg[iarg + 1], "error") == 0)
        lostflag = Thermo::ERROR;
      else
        error->all(FLERR, "Unknown thermo_modify lost argument: {}", arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "lost/bond") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify lost/bond", error);
      if (strcmp(arg[iarg + 1], "ignore") == 0)
        lostbond = Thermo::IGNORE;
      else if (strcmp(arg[iarg + 1], "warn") == 0)
        lostbond = Thermo::WARN;
      else if (strcmp(arg[iarg + 1], "error") == 0)
        lostbond = Thermo::ERROR;
      else
        error->all(FLERR, "Unknown thermo_modify lost/bond argument: {}", arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "warn") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify warn", error);
      if (strcmp(arg[iarg + 1], "ignore") == 0)
        error->set_maxwarn(-1);
      else if (strcmp(arg[iarg + 1], "always") == 0)
        error->set_maxwarn(0);
      else if (strcmp(arg[iarg + 1], "reset") == 0) {
        error->set_numwarn(0);
        warnbefore = 0;
      } else if (strcmp(arg[iarg + 1], "default") == 0) {
        warnbefore = 0;
        error->set_numwarn(0);
        error->set_maxwarn(100);
      } else
        error->set_maxwarn(utils::inumeric(FLERR, arg[iarg + 1], false, lmp));
      iarg += 2;

    } else if (strcmp(arg[iarg], "norm") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify norm", error);
      normuserflag = 1;
      normuser = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg], "flush") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify flush", error);
      flushflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg], "line") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify line", error);
      if (strcmp(arg[iarg + 1], "one") == 0)
        lineflag = ONELINE;
      else if (strcmp(arg[iarg + 1], "multi") == 0)
        lineflag = MULTILINE;
      else if (strcmp(arg[iarg + 1], "yaml") == 0)
        lineflag = YAMLLINE;
      else
        error->all(FLERR, "Unknown thermo_modify line argument: {}", arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "colname") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify colname", error);
      if (strcmp(arg[iarg + 1], "default") == 0) {
        for (auto &item : keyword_user) item.clear();
        iarg += 2;
      } else {
        if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "thermo_modify colname", error);
        int icol = -1;
        if (utils::is_integer(arg[iarg + 1])) {
          icol = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
          if (icol < 0) icol = nfield_initial + icol + 1;
          icol--;
        } else {
          try {
            icol = key2col.at(arg[iarg + 1]);
          } catch (std::out_of_range &) {
            icol = -1;
          }
        }
        if ((icol < 0) || (icol >= nfield_initial))
          error->all(FLERR, "Invalid thermo_modify colname argument: {}", arg[iarg + 1]);
        keyword_user[icol] = arg[iarg + 2];
        iarg += 3;
      }
    } else if (strcmp(arg[iarg], "format") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "thermo_modify format", error);

      if (strcmp(arg[iarg + 1], "none") == 0) {
        format_line_user.clear();
        format_int_user.clear();
        format_bigint_user.clear();
        format_float_user.clear();
        for (auto &item : format_column_user) item.clear();
        iarg += 2;
        continue;
      }

      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "thermo_modify format", error);

      if (strcmp(arg[iarg + 1], "line") == 0) {
        format_line_user = arg[iarg + 2];
      } else if (strcmp(arg[iarg + 1], "int") == 0) {
        format_int_user = arg[iarg + 2];
        // replace "d" in format_int_user with bigint format specifier
        auto found = format_int_user.find('%');
        found = format_int_user.find('d', found);
        if (found == std::string::npos)
          error->all(FLERR, "Thermo_modify int format does not contain a d conversion character");
        format_bigint_user =
            format_int_user.replace(found, 1, std::string(BIGINT_FORMAT).substr(1));
      } else if (strcmp(arg[iarg + 1], "float") == 0) {
        format_float_user = arg[iarg + 2];
      } else if (utils::strmatch(arg[iarg + 1], "^\\d*\\*\\d*$")) {
        // handles cases such as 2*6; currently doesn't allow negatives
        int nlo, nhi;
        utils::bounds(FLERR, arg[iarg + 1], 1, nfield_initial, nlo, nhi, error);
        int icol = -1;
        for (int i = nlo - 1; i < nhi; i++) {
          if (i < 0) icol = nfield_initial + i + 1; // doesn't happen currently
          else icol = i;
          if (icol < 0 || (icol >= nfield_initial))
            error->all(FLERR, "Invalid thermo_modify format argument: {}",
              arg[iarg + 1]);
          format_column_user[icol] = arg[iarg + 2];
        }
      } else {
        int icol = -1;
        if (utils::is_integer(arg[iarg + 1])) {
          icol = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
          if (icol < 0) icol = nfield_initial + icol + 1;
          icol--;
        } else {
          try {
            icol = key2col.at(arg[iarg + 1]);
          } catch (std::out_of_range &) {
            icol = -1;
          }
        }
        if ((icol < 0) || (icol >= nfield_initial))
          error->all(FLERR, "Invalid thermo_modify format argument: {}", arg[iarg + 1]);
        format_column_user[icol] = arg[iarg + 2];
      }
      iarg += 3;

    } else
      error->all(FLERR, "Unknown thermo_modify keyword: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   allocate all per-field memory
------------------------------------------------------------------------- */

void Thermo::allocate()
{
  // n = specified fields + Volume field (added at run time)

  const int n = nfield_initial + 1;

  keyword.resize(n);
  format.resize(n);
  format_column_user.resize(n);
  keyword_user.resize(n);
  for (int i = 0; i < n; i++) {
    keyword[i].clear();
    format[i].clear();
    format_column_user[i].clear();
    keyword_user[i].clear();
  }

  vfunc = new FnPtr[n];
  vtype = new int[n];

  field2index = new int[n];
  argindex1 = new int[n];
  argindex2 = new int[n];

  // factor of 3 is max number of computes a single field can add

  ncompute = 0;
  id_compute = new char *[3 * n];
  compute_which = new int[3 * n];
  computes = new Compute *[3 * n];

  nfix = 0;
  id_fix = new char *[n];
  fixes = new Fix *[n];

  nvariable = 0;
  id_variable = new char *[n];
  variables = new int[n];

  int i = 0;
  key2col.clear();
  for (const auto &item : utils::split_words(line)) { key2col[item] = i++; }
}

/* ----------------------------------------------------------------------
   deallocate all per-field memory
------------------------------------------------------------------------- */

void Thermo::deallocate()
{
  delete[] vfunc;
  delete[] vtype;

  delete[] field2index;
  delete[] argindex1;
  delete[] argindex2;

  for (int i = 0; i < ncompute; i++) delete[] id_compute[i];
  delete[] id_compute;
  delete[] compute_which;
  delete[] computes;

  for (int i = 0; i < nfix; i++) delete[] id_fix[i];
  delete[] id_fix;
  delete[] fixes;

  for (int i = 0; i < nvariable; i++) delete[] id_variable[i];
  delete[] id_variable;
  delete[] variables;
}

/* ----------------------------------------------------------------------
   parse list of thermo keywords from str
   set compute flags (temp, press, pe, etc)
------------------------------------------------------------------------- */

void Thermo::parse_fields(const std::string &str)
{
  nfield = 0;

  // CUSTOMIZATION: add a new thermo keyword by adding it to the if statements

  ValueTokenizer keywords(str);
  while (keywords.has_next()) {
    std::string word = keywords.next_string();

    if (word == "step") {
      addfield("Step", &Thermo::compute_step, BIGINT);
    } else if (word == "elapsed") {
      addfield("Elapsed", &Thermo::compute_elapsed, BIGINT);
    } else if (word == "elaplong") {
      addfield("Elaplong", &Thermo::compute_elapsed_long, BIGINT);
    } else if (word == "dt") {
      addfield("Dt", &Thermo::compute_dt, FLOAT);
    } else if (word == "time") {
      addfield("Time", &Thermo::compute_time, FLOAT);
    } else if (word == "cpu") {
      addfield("CPU", &Thermo::compute_cpu, FLOAT);
    } else if (word == "tpcpu") {
      addfield("T/CPU", &Thermo::compute_tpcpu, FLOAT);
    } else if (word == "spcpu") {
      addfield("S/CPU", &Thermo::compute_spcpu, FLOAT);
    } else if (word == "cpuremain") {
      addfield("CPULeft", &Thermo::compute_cpuremain, FLOAT);
    } else if (word == "part") {
      addfield("Part", &Thermo::compute_part, INT);
    } else if (word == "timeremain") {
      addfield("TimeoutLeft", &Thermo::compute_timeremain, FLOAT);

    } else if (word == "atoms") {
      addfield("Atoms", &Thermo::compute_atoms, BIGINT);
    } else if (word == "temp") {
      addfield("Temp", &Thermo::compute_temp, FLOAT);
      index_temp = add_compute(id_temp, SCALAR);
    } else if (word == "press") {
      addfield("Press", &Thermo::compute_press, FLOAT);
      index_press_scalar = add_compute(id_press, SCALAR);
    } else if (word == "pe") {
      addfield("PotEng", &Thermo::compute_pe, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "ke") {
      addfield("KinEng", &Thermo::compute_ke, FLOAT);
      index_temp = add_compute(id_temp, SCALAR);
    } else if (word == "etotal") {
      addfield("TotEng", &Thermo::compute_etotal, FLOAT);
      index_temp = add_compute(id_temp, SCALAR);
      index_pe = add_compute(id_pe, SCALAR);

    } else if (word == "evdwl") {
      addfield("E_vdwl", &Thermo::compute_evdwl, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "ecoul") {
      addfield("E_coul", &Thermo::compute_ecoul, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "epair") {
      addfield("E_pair", &Thermo::compute_epair, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "ebond") {
      addfield("E_bond", &Thermo::compute_ebond, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "eangle") {
      addfield("E_angle", &Thermo::compute_eangle, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "edihed") {
      addfield("E_dihed", &Thermo::compute_edihed, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "eimp") {
      addfield("E_impro", &Thermo::compute_eimp, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "emol") {
      addfield("E_mol", &Thermo::compute_emol, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "elong") {
      addfield("E_long", &Thermo::compute_elong, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "etail") {
      addfield("E_tail", &Thermo::compute_etail, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);

    } else if (word == "enthalpy") {
      addfield("Enthalpy", &Thermo::compute_enthalpy, FLOAT);
      index_temp = add_compute(id_temp, SCALAR);
      index_press_scalar = add_compute(id_press, SCALAR);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "ecouple") {
      addfield("Ecouple", &Thermo::compute_ecouple, FLOAT);
      index_pe = add_compute(id_pe, SCALAR);
    } else if (word == "econserve") {
      addfield("Econserve", &Thermo::compute_econserve, FLOAT);
      index_temp = add_compute(id_temp, SCALAR);
      index_pe = add_compute(id_pe, SCALAR);

    } else if (word == "vol") {
      addfield("Volume", &Thermo::compute_vol, FLOAT);
    } else if (word == "density") {
      addfield("Density", &Thermo::compute_density, FLOAT);
    } else if (word == "lx") {
      addfield("Lx", &Thermo::compute_lx, FLOAT);
    } else if (word == "ly") {
      addfield("Ly", &Thermo::compute_ly, FLOAT);
    } else if (word == "lz") {
      addfield("Lz", &Thermo::compute_lz, FLOAT);

    } else if (word == "xlo") {
      addfield("Xlo", &Thermo::compute_xlo, FLOAT);
    } else if (word == "xhi") {
      addfield("Xhi", &Thermo::compute_xhi, FLOAT);
    } else if (word == "ylo") {
      addfield("Ylo", &Thermo::compute_ylo, FLOAT);
    } else if (word == "yhi") {
      addfield("Yhi", &Thermo::compute_yhi, FLOAT);
    } else if (word == "zlo") {
      addfield("Zlo", &Thermo::compute_zlo, FLOAT);
    } else if (word == "zhi") {
      addfield("Zhi", &Thermo::compute_zhi, FLOAT);

    } else if (word == "xy") {
      addfield("Xy", &Thermo::compute_xy, FLOAT);
    } else if (word == "xz") {
      addfield("Xz", &Thermo::compute_xz, FLOAT);
    } else if (word == "yz") {
      addfield("Yz", &Thermo::compute_yz, FLOAT);

    } else if (word == "xlat") {
      addfield("Xlat", &Thermo::compute_xlat, FLOAT);
    } else if (word == "ylat") {
      addfield("Ylat", &Thermo::compute_ylat, FLOAT);
    } else if (word == "zlat") {
      addfield("Zlat", &Thermo::compute_zlat, FLOAT);

    } else if (word == "bonds") {
      addfield("Bonds", &Thermo::compute_bonds, BIGINT);
    } else if (word == "angles") {
      addfield("Angles", &Thermo::compute_angles, BIGINT);
    } else if (word == "dihedrals") {
      addfield("Diheds", &Thermo::compute_dihedrals, BIGINT);
    } else if (word == "impropers") {
      addfield("Impros", &Thermo::compute_impropers, BIGINT);

    } else if (word == "pxx") {
      addfield("Pxx", &Thermo::compute_pxx, FLOAT);
      index_press_vector = add_compute(id_press, VECTOR);
    } else if (word == "pyy") {
      addfield("Pyy", &Thermo::compute_pyy, FLOAT);
      index_press_vector = add_compute(id_press, VECTOR);
    } else if (word == "pzz") {
      addfield("Pzz", &Thermo::compute_pzz, FLOAT);
      index_press_vector = add_compute(id_press, VECTOR);
    } else if (word == "pxy") {
      addfield("Pxy", &Thermo::compute_pxy, FLOAT);
      index_press_vector = add_compute(id_press, VECTOR);
    } else if (word == "pxz") {
      addfield("Pxz", &Thermo::compute_pxz, FLOAT);
      index_press_vector = add_compute(id_press, VECTOR);
    } else if (word == "pyz") {
      addfield("Pyz", &Thermo::compute_pyz, FLOAT);
      index_press_vector = add_compute(id_press, VECTOR);

    } else if (word == "fmax") {
      addfield("Fmax", &Thermo::compute_fmax, FLOAT);
    } else if (word == "fnorm") {
      addfield("Fnorm", &Thermo::compute_fnorm, FLOAT);

    } else if (word == "nbuild") {
      addfield("Nbuild", &Thermo::compute_nbuild, BIGINT);
    } else if (word == "ndanger") {
      addfield("Ndanger", &Thermo::compute_ndanger, BIGINT);

    } else if (word == "cella") {
      addfield("Cella", &Thermo::compute_cella, FLOAT);
    } else if (word == "cellb") {
      addfield("Cellb", &Thermo::compute_cellb, FLOAT);
    } else if (word == "cellc") {
      addfield("Cellc", &Thermo::compute_cellc, FLOAT);
    } else if (word == "cellalpha") {
      addfield("CellAlpha", &Thermo::compute_cellalpha, FLOAT);
    } else if (word == "cellbeta") {
      addfield("CellBeta", &Thermo::compute_cellbeta, FLOAT);
    } else if (word == "cellgamma") {
      addfield("CellGamma", &Thermo::compute_cellgamma, FLOAT);

      // compute value = c_ID, fix value = f_ID, variable value = v_ID
      // count trailing [] and store int arguments

    } else {
      ArgInfo argi(word);

      if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_type() == ArgInfo::NONE) ||
          (argi.get_dim() > 2))
        error->all(FLERR, "Unknown keyword '{}' in thermo_style custom command", word);

      // process zero or one or two trailing brackets
      // argindex1,argindex2 = int inside each bracket pair, 0 if no bracket

      argindex1[nfield] = argi.get_index1();
      argindex2[nfield] = (argi.get_dim() > 1) ? argi.get_index2() : 0;

      if (argi.get_type() == ArgInfo::COMPUTE) {
        auto icompute = modify->get_compute_by_id(argi.get_name());
        if (!icompute)
          error->all(FLERR, "Could not find thermo custom compute ID: {}", argi.get_name());
        if (argindex1[nfield] == 0 && icompute->scalar_flag == 0)
          error->all(FLERR, "Thermo compute does not compute scalar");
        if (argindex1[nfield] > 0 && argindex2[nfield] == 0) {
          if (icompute->vector_flag == 0)
            error->all(FLERR, "Thermo compute does not compute vector");
          if (argindex1[nfield] > icompute->size_vector && icompute->size_vector_variable == 0)
            error->all(FLERR, "Thermo compute vector is accessed out-of-range");
        }
        if (argindex1[nfield] > 0 && argindex2[nfield] > 0) {
          if (icompute->array_flag == 0) error->all(FLERR, "Thermo compute does not compute array");
          if (argindex1[nfield] > icompute->size_array_rows &&
              icompute->size_array_rows_variable == 0)
            error->all(FLERR, "Thermo compute array is accessed out-of-range");
          if (argindex2[nfield] > icompute->size_array_cols)
            error->all(FLERR, "Thermo compute array is accessed out-of-range");
        }

        if (argindex1[nfield] == 0)
          field2index[nfield] = add_compute(argi.get_name(), SCALAR);
        else if (argindex2[nfield] == 0)
          field2index[nfield] = add_compute(argi.get_name(), VECTOR);
        else
          field2index[nfield] = add_compute(argi.get_name(), ARRAY);
        addfield(word.c_str(), &Thermo::compute_compute, FLOAT);

      } else if (argi.get_type() == ArgInfo::FIX) {
        auto ifix = modify->get_fix_by_id(argi.get_name());
        if (!ifix) error->all(FLERR, "Could not find thermo custom fix ID: {}", argi.get_name());
        if (argindex1[nfield] == 0 && ifix->scalar_flag == 0)
          error->all(FLERR, "Thermo fix does not compute scalar");
        if (argindex1[nfield] > 0 && argindex2[nfield] == 0) {
          if (ifix->vector_flag == 0) error->all(FLERR, "Thermo fix does not compute vector");
          if (argindex1[nfield] > ifix->size_vector && ifix->size_vector_variable == 0)
            error->all(FLERR, "Thermo fix vector is accessed out-of-range");
        }
        if (argindex1[nfield] > 0 && argindex2[nfield] > 0) {
          if (ifix->array_flag == 0) error->all(FLERR, "Thermo fix does not compute array");
          if (argindex1[nfield] > ifix->size_array_rows && ifix->size_array_rows_variable == 0)
            error->all(FLERR, "Thermo fix array is accessed out-of-range");
          if (argindex2[nfield] > ifix->size_array_cols)
            error->all(FLERR, "Thermo fix array is accessed out-of-range");
        }

        field2index[nfield] = add_fix(argi.get_name());
        addfield(word.c_str(), &Thermo::compute_fix, FLOAT);

      } else if (argi.get_type() == ArgInfo::VARIABLE) {
        int n = input->variable->find(argi.get_name());
        if (n < 0)
          error->all(FLERR, "Could not find thermo custom variable name: {}", argi.get_name());
        if (argindex1[nfield] == 0 && input->variable->equalstyle(n) == 0)
          error->all(FLERR, "Thermo custom variable is not equal-style variable");
        if (argindex1[nfield] && input->variable->vectorstyle(n) == 0)
          error->all(FLERR, "Thermo custom variable is not vector-style variable");
        if (argindex2[nfield]) error->all(FLERR, "Thermo custom variable cannot have two indices");

        field2index[nfield] = add_variable(argi.get_name());
        addfield(word.c_str(), &Thermo::compute_variable, FLOAT);
      }
    }
  }
  field_data.clear();
  field_data.resize(nfield);
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Thermo::addfield(const char *key, FnPtr func, int typeflag)
{
  keyword[nfield] = key;
  vfunc[nfield] = func;
  vtype[nfield] = typeflag;
  nfield++;
}

/* ----------------------------------------------------------------------
   add compute ID to list of Compute objects to call
   return location of where this Compute is in list
   if already in list with same which, do not add, just return index
------------------------------------------------------------------------- */

int Thermo::add_compute(const char *id, int which)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if ((strcmp(id, id_compute[icompute]) == 0) && which == compute_which[icompute]) break;
  if (icompute < ncompute) return icompute;

  id_compute[ncompute] = utils::strdup(id);
  compute_which[ncompute] = which;
  ncompute++;
  return ncompute - 1;
}

/* ----------------------------------------------------------------------
   add fix ID to list of Fix objects to call
------------------------------------------------------------------------- */

int Thermo::add_fix(const char *id)
{
  id_fix[nfix] = utils::strdup(id);
  nfix++;
  return nfix - 1;
}

/* ----------------------------------------------------------------------
   add variable ID to list of Variables to evaluate
------------------------------------------------------------------------- */

int Thermo::add_variable(const char *id)
{
  id_variable[nvariable] = utils::strdup(id);
  nvariable++;
  return nvariable - 1;
}

/* ----------------------------------------------------------------------
   check whether temperature compute is defined, active, and needs invoking
------------------------------------------------------------------------- */

void Thermo::check_temp(const std::string &keyword)
{
  if (!temperature)
    error->all(FLERR, "Thermo keyword {} in variable requires thermo to use/init temperature",
               keyword);
  if (!temperature->is_initialized())
    error->all(FLERR,"Thermo keyword {} cannot be invoked before initialization by a run",keyword);
  if (!(temperature->invoked_flag & Compute::INVOKED_SCALAR)) {
    temperature->compute_scalar();
    temperature->invoked_flag |= Compute::INVOKED_SCALAR;
  }
}

/* ----------------------------------------------------------------------
   check whether potential energy compute is defined, active, and needs invoking
------------------------------------------------------------------------- */

void Thermo::check_pe(const std::string &keyword)
{
  if (update->eflag_global != update->ntimestep)
    error->all(FLERR, "Energy was not tallied on needed timestep");
  if (!pe)
    error->all(FLERR, "Thermo keyword {} in variable requires thermo to use/init potential energy",
               keyword);
  if (!pe->is_initialized())
    error->all(FLERR,"Thermo keyword {} cannot be invoked before initialization by a run",keyword);
  if (!(pe->invoked_flag & Compute::INVOKED_SCALAR)) {
    pe->compute_scalar();
    pe->invoked_flag |= Compute::INVOKED_SCALAR;
  }
}

/* ----------------------------------------------------------------------
   check whether scalar pressure compute is defined, active, and needs invoking
------------------------------------------------------------------------- */

void Thermo::check_press_scalar(const std::string &keyword)
{
  if (!pressure)
    error->all(FLERR, "Thermo keyword {} in variable requires thermo to use/init press", keyword);
  if (!pressure->is_initialized())
    error->all(FLERR,"Thermo keyword {} cannot be invoked before initialization by a run",keyword);
  if (!(pressure->invoked_flag & Compute::INVOKED_SCALAR)) {
    pressure->compute_scalar();
    pressure->invoked_flag |= Compute::INVOKED_SCALAR;
  }
}

/* ----------------------------------------------------------------------
   check whether tensor pressure compute is defined, active, and needs invoking
------------------------------------------------------------------------- */

void Thermo::check_press_vector(const std::string &keyword)
{
  if (!pressure)
    error->all(FLERR, "Thermo keyword {} in variable requires thermo to use/init press", keyword);
  if (!pressure->is_initialized())
    error->all(FLERR,"Thermo keyword {} cannot be invoked before initialization by a run",keyword);
  if (!(pressure->invoked_flag & Compute::INVOKED_VECTOR)) {
    pressure->compute_vector();
    pressure->invoked_flag |= Compute::INVOKED_VECTOR;
  }
}

/* ----------------------------------------------------------------------
   compute a single thermodynamic value, word is any keyword in custom list
   called when a variable is evaluated by Variable class
   return value as double in answer
   return 0 if str is recognized keyword, 1 if unrecognized
   CUSTOMIZATION: add a new keyword by adding a suitable if statement
------------------------------------------------------------------------- */

int Thermo::evaluate_keyword(const std::string &word, double *answer)
{
  // turn off normflag if natoms = 0 to avoid divide by 0
  // normflag must be set for lo-level thermo routines that may be invoked

  natoms = atom->natoms;
  if (natoms == 0)
    normflag = 0;
  else
    normflag = normvalue;

  // invoke a lo-level thermo routine to compute the variable value
  // if keyword requires a compute, error if thermo doesn't use the compute
  // for keywords that use energy (evdwl, ebond, etc):
  //   check if energy was tallied on this timestep and set pe->invoked_flag
  //   this will trigger next timestep for energy tallying via addstep()
  //   this means keywords that use pe (pe, etotal, enthalpy)
  //     need to always invoke it even if invoked_flag is set,
  //     because evdwl/etc may have set invoked_flag w/out
  //       actually invoking pe->compute_scalar()

  if (word == "step") {
    compute_step();
    dvalue = bivalue;

  } else if (word == "elapsed") {
    if (update->whichflag == 0)
      error->all(FLERR, "The variable thermo keyword elapsed cannot be used between runs");
    compute_elapsed();
    dvalue = bivalue;

  } else if (word == "elaplong") {
    if (update->whichflag == 0)
      error->all(FLERR, "The variable thermo keyword elaplong cannot be used between runs");
    compute_elapsed_long();
    dvalue = bivalue;

  } else if (word == "dt") {
    compute_dt();

  } else if (word == "time") {
    compute_time();

  } else if (word == "cpu") {
    if (update->whichflag == 0)
      error->all(FLERR, "The variable thermo keyword cpu cannot be used between runs");
    compute_cpu();

  } else if (word == "tpcpu") {
    if (update->whichflag == 0)
      error->all(FLERR, "The variable thermo keyword tpcpu cannot be used between runs");
    compute_tpcpu();

  } else if (word == "spcpu") {
    if (update->whichflag == 0)
      error->all(FLERR, "The variable thermo keyword spcpu cannot be used between runs");
    compute_spcpu();

  } else if (word == "cpuremain") {
    if (update->whichflag == 0)
      error->all(FLERR, "The variable thermo keyword cpuremain cannot be used between runs");
    compute_cpuremain();

  } else if (word == "part") {
    compute_part();
    dvalue = ivalue;

  } else if (word == "timeremain") {
    compute_timeremain();

  } else if (word == "atoms") {
    compute_atoms();
    dvalue = bivalue;

  } else if (word == "bonds") {
    compute_bonds();
    dvalue = bivalue;

  } else if (word == "angles") {
    compute_angles();
    dvalue = bivalue;

  } else if (word == "dihedrals") {
    compute_dihedrals();
    dvalue = bivalue;

  } else if (word == "impropers") {
    compute_impropers();
    dvalue = bivalue;

  } else if (word == "temp") {
    check_temp(word);
    compute_temp();

  } else if (word == "press") {
    check_press_scalar(word);
    compute_press();

  } else if (word == "pe") {
    check_pe(word);
    compute_pe();

  } else if (word == "ke") {
    check_temp(word);
    compute_ke();

  } else if (word == "etotal") {
    check_pe(word);
    check_temp(word);
    compute_etotal();

  } else if (word == "evdwl") {
    check_pe(word);
    compute_evdwl();

  } else if (word == "ecoul") {
    check_pe(word);
    compute_ecoul();

  } else if (word == "epair") {
    check_pe(word);
    compute_epair();

  } else if (word == "ebond") {
    check_pe(word);
    compute_ebond();

  } else if (word == "eangle") {
    check_pe(word);
    compute_eangle();

  } else if (word == "edihed") {
    check_pe(word);
    compute_edihed();

  } else if (word == "eimp") {
    check_pe(word);
    compute_eimp();

  } else if (word == "emol") {
    check_pe(word);
    compute_emol();

  } else if (word == "elong") {
    check_pe(word);
    compute_elong();

  } else if (word == "etail") {
    if (update->eflag_global != update->ntimestep)
      error->all(FLERR, "Energy was not tallied on needed timestep for thermo keyword etail");
    compute_etail();

  } else if (word == "enthalpy") {
    check_pe(word);
    check_temp(word);
    check_press_scalar(word);
    compute_enthalpy();

  } else if (word == "ecouple")
    compute_ecouple();

  else if (word == "econserve") {
    check_pe(word);
    check_temp(word);
    compute_econserve();

  } else if (word == "vol")
    compute_vol();
  else if (word == "density")
    compute_density();
  else if (word == "lx")
    compute_lx();
  else if (word == "ly")
    compute_ly();
  else if (word == "lz")
    compute_lz();

  else if (word == "xlo")
    compute_xlo();
  else if (word == "xhi")
    compute_xhi();
  else if (word == "ylo")
    compute_ylo();
  else if (word == "yhi")
    compute_yhi();
  else if (word == "zlo")
    compute_zlo();
  else if (word == "zhi")
    compute_zhi();

  else if (word == "xy")
    compute_xy();
  else if (word == "xz")
    compute_xz();
  else if (word == "yz")
    compute_yz();

  else if (word == "xlat")
    compute_xlat();
  else if (word == "ylat")
    compute_ylat();
  else if (word == "zlat")
    compute_zlat();

  else if (word == "pxx") {
    check_press_vector(word);
    compute_pxx();

  } else if (word == "pyy") {
    check_press_vector(word);
    compute_pyy();

  } else if (word == "pzz") {
    check_press_vector(word);
    compute_pzz();

  } else if (word == "pxy") {
    check_press_vector(word);
    compute_pxy();

  } else if (word == "pxz") {
    check_press_vector(word);
    compute_pxz();

  } else if (word == "pyz") {
    check_press_vector(word);
    compute_pyz();
  }

  else if (word == "fmax")
    compute_fmax();
  else if (word == "fnorm")
    compute_fnorm();

  else if (word == "nbuild") {
    compute_nbuild();
    dvalue = bivalue;
  } else if (word == "ndanger") {
    compute_ndanger();
    dvalue = bivalue;
  }

  else if (word == "cella")
    compute_cella();
  else if (word == "cellb")
    compute_cellb();
  else if (word == "cellc")
    compute_cellc();
  else if (word == "cellalpha")
    compute_cellalpha();
  else if (word == "cellbeta")
    compute_cellbeta();
  else if (word == "cellgamma")
    compute_cellgamma();

  else
    return 1;

  *answer = dvalue;
  return 0;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
   compute/fix are normalized by atoms if returning extensive value
   variable value is not normalized (formula should normalize if desired)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Thermo::compute_compute()
{
  int m = field2index[ifield];
  Compute *compute = computes[m];

  // check for out-of-range access if vector/array is variable length

  if (compute_which[m] == SCALAR) {
    dvalue = compute->scalar;
    if (normflag && compute->extscalar) dvalue /= natoms;
  } else if (compute_which[m] == VECTOR) {
    if (compute->size_vector_variable && argindex1[ifield] > compute->size_vector)
      dvalue = 0.0;
    else
      dvalue = compute->vector[argindex1[ifield] - 1];
    if (normflag) {
      if (compute->extvector == 0)
        return;
      else if (compute->extvector == 1)
        dvalue /= natoms;
      else if (compute->extlist[argindex1[ifield] - 1])
        dvalue /= natoms;
    }
  } else {
    if (compute->size_array_rows_variable && argindex1[ifield] > compute->size_array_rows)
      dvalue = 0.0;
    else
      dvalue = compute->array[argindex1[ifield] - 1][argindex2[ifield] - 1];
    if (normflag && compute->extarray) dvalue /= natoms;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_fix()
{
  int m = field2index[ifield];
  Fix *fix = fixes[m];

  if (argindex1[ifield] == 0) {
    dvalue = fix->compute_scalar();
    if (normflag && fix->extscalar) dvalue /= natoms;
  } else if (argindex2[ifield] == 0) {
    dvalue = fix->compute_vector(argindex1[ifield] - 1);
    if (normflag) {
      if (fix->extvector == 0)
        return;
      else if (fix->extvector == 1)
        dvalue /= natoms;
      else if (fix->extlist[argindex1[ifield] - 1])
        dvalue /= natoms;
    }
  } else {
    dvalue = fix->compute_array(argindex1[ifield] - 1, argindex2[ifield] - 1);
    if (normflag && fix->extarray) dvalue /= natoms;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_variable()
{
  int iarg = argindex1[ifield];

  if (iarg == 0)
    dvalue = input->variable->compute_equal(variables[field2index[ifield]]);
  else {
    double *varvec;
    int nvec = input->variable->compute_vector(variables[field2index[ifield]], &varvec);
    if (nvec < iarg)
      dvalue = 0.0;
    else
      dvalue = varvec[iarg - 1];
  }
}

/* ----------------------------------------------------------------------
   one method for every keyword thermo can output
   called by compute() or evaluate_keyword()
   compute will have already been called
   set ivalue/dvalue/bivalue if value is int/double/bigint
   CUSTOMIZATION: add a new thermo keyword by adding a new method
------------------------------------------------------------------------- */

void Thermo::compute_step()
{
  bivalue = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_elapsed()
{
  bivalue = update->ntimestep - update->firststep;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_elapsed_long()
{
  bivalue = update->ntimestep - update->beginstep;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_dt()
{
  dvalue = update->dt;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_time()
{
  dvalue = update->atime + (update->ntimestep - update->atimestep) * update->dt;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cpu()
{
  if (firststep == 0)
    dvalue = 0.0;
  else
    dvalue = timer->elapsed(Timer::TOTAL);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_tpcpu()
{
  double new_cpu;
  double new_time = update->ntimestep * update->dt;

  if (firststep == 0) {
    new_cpu = 0.0;
    dvalue = 0.0;
  } else {
    new_cpu = timer->elapsed(Timer::TOTAL);
    double cpu_diff = new_cpu - last_tpcpu;
    double time_diff = new_time - last_time;
    if (time_diff > 0.0 && cpu_diff > 0.0)
      dvalue = time_diff / cpu_diff;
    else
      dvalue = 0.0;
  }

  last_time = new_time;
  last_tpcpu = new_cpu;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_spcpu()
{
  double new_cpu;
  int new_step = update->ntimestep;

  if (firststep == 0) {
    new_cpu = 0.0;
    dvalue = 0.0;
  } else {
    new_cpu = timer->elapsed(Timer::TOTAL);
    double cpu_diff = new_cpu - last_spcpu;
    int step_diff = new_step - last_step;
    if (cpu_diff > 0.0)
      dvalue = step_diff / cpu_diff;
    else
      dvalue = 0.0;
  }

  last_step = new_step;
  last_spcpu = new_cpu;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cpuremain()
{
  if (firststep == 0)
    dvalue = 0.0;
  else
    dvalue = timer->elapsed(Timer::TOTAL) * (update->laststep - update->ntimestep) /
        (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_part()
{
  ivalue = universe->iworld;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_timeremain()
{
  dvalue = timer->get_timeout_remain();
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_atoms()
{
  bivalue = group->count_all();
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_temp()
{
  dvalue = temperature->scalar;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_press()
{
  dvalue = pressure->scalar;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pe()
{
  dvalue = pe->scalar;
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ke()
{
  dvalue = temperature->scalar;
  dvalue *= 0.5 * temperature->dof * force->boltz;
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_etotal()
{
  compute_pe();
  double dvalue_pe = dvalue;
  compute_ke();
  dvalue += dvalue_pe;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ecouple()
{
  dvalue = modify->energy_couple();
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_econserve()
{
  compute_etotal();
  double dvalue_etotal = dvalue;
  compute_ecouple();
  dvalue += dvalue_etotal;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_evdwl()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl;
  MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);

  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }

  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ecoul()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_coul;
  MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);
  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_epair()
{
  double tmp = 0.0;
  if (force->pair) tmp += force->pair->eng_vdwl + force->pair->eng_coul;
  MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);

  if (force->kspace) dvalue += force->kspace->energy;
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue += force->pair->etail / volume;
  }

  if (normflag) dvalue /= natoms;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ebond()
{
  if (force->bond) {
    double tmp = force->bond->energy;
    MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);
    if (normflag) dvalue /= natoms;
  } else
    dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eangle()
{
  if (force->angle) {
    double tmp = force->angle->energy;
    MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);
    if (normflag) dvalue /= natoms;
  } else
    dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_edihed()
{
  if (force->dihedral) {
    double tmp = force->dihedral->energy;
    MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);
    if (normflag) dvalue /= natoms;
  } else
    dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_eimp()
{
  if (force->improper) {
    double tmp = force->improper->energy;
    MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);
    if (normflag) dvalue /= natoms;
  } else
    dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_emol()
{
  double tmp = 0.0;
  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) tmp += force->bond->energy;
    if (force->angle) tmp += force->angle->energy;
    if (force->dihedral) tmp += force->dihedral->energy;
    if (force->improper) tmp += force->improper->energy;
    MPI_Allreduce(&tmp, &dvalue, 1, MPI_DOUBLE, MPI_SUM, world);
    if (normflag) dvalue /= natoms;
  } else
    dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_elong()
{
  if (force->kspace) {
    dvalue = force->kspace->energy;
    if (normflag) dvalue /= natoms;
  } else
    dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_etail()
{
  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    dvalue = force->pair->etail / volume;
    if (normflag) dvalue /= natoms;
  } else
    dvalue = 0.0;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_enthalpy()
{
  compute_etotal();
  double etmp = dvalue;

  compute_vol();
  double vtmp = dvalue;
  if (normflag) vtmp /= natoms;

  compute_press();
  double ptmp = dvalue;

  dvalue = etmp + ptmp * vtmp / (force->nktv2p);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_vol()
{
  if (domain->dimension == 3)
    dvalue = domain->xprd * domain->yprd * domain->zprd;
  else
    dvalue = domain->xprd * domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_density()
{
  double mass = group->mass(0);
  compute_vol();
  dvalue = force->mv2d * mass / dvalue;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_lx()
{
  dvalue = domain->xprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ly()
{
  dvalue = domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_lz()
{
  dvalue = domain->zprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xlo()
{
  dvalue = domain->boxlo[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xhi()
{
  dvalue = domain->boxhi[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ylo()
{
  dvalue = domain->boxlo[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_yhi()
{
  dvalue = domain->boxhi[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_zlo()
{
  dvalue = domain->boxlo[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_zhi()
{
  dvalue = domain->boxhi[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xy()
{
  dvalue = domain->xy;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xz()
{
  dvalue = domain->xz;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_yz()
{
  dvalue = domain->yz;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_xlat()
{
  dvalue = domain->lattice->xlattice;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ylat()
{
  dvalue = domain->lattice->ylattice;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_zlat()
{
  dvalue = domain->lattice->zlattice;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_bonds()
{
  bivalue = atom->nbonds;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_angles()
{
  bivalue = atom->nangles;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_dihedrals()
{
  bivalue = atom->ndihedrals;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_impropers()
{
  bivalue = atom->nimpropers;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxx()
{
  dvalue = pressure->vector[0];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyy()
{
  dvalue = pressure->vector[1];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pzz()
{
  dvalue = pressure->vector[2];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxy()
{
  dvalue = pressure->vector[3];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pxz()
{
  dvalue = pressure->vector[4];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_pyz()
{
  dvalue = pressure->vector[5];
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_fmax()
{
  double **f = atom->f;
  int nlocal = atom->nlocal;

  double max = 0.0;
  for (int i = 0; i < nlocal; i++) {
    max = MAX(max, fabs(f[i][0]));
    max = MAX(max, fabs(f[i][1]));
    max = MAX(max, fabs(f[i][2]));
  }
  double maxall;
  MPI_Allreduce(&max, &maxall, 1, MPI_DOUBLE, MPI_MAX, world);
  dvalue = maxall;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_fnorm()
{
  double **f = atom->f;
  int nlocal = atom->nlocal;

  double dot = 0.0;
  for (int i = 0; i < nlocal; i++) dot += f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
  double dotall;
  MPI_Allreduce(&dot, &dotall, 1, MPI_DOUBLE, MPI_SUM, world);
  dvalue = sqrt(dotall);
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_nbuild()
{
  bivalue = neighbor->ncalls;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_ndanger()
{
  bivalue = neighbor->ndanger;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cella()
{
  dvalue = domain->xprd;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cellb()
{
  if (!domain->triclinic)
    dvalue = domain->yprd;
  else {
    double *h = domain->h;
    dvalue = sqrt(h[1] * h[1] + h[5] * h[5]);
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cellc()
{
  if (!domain->triclinic)
    dvalue = domain->zprd;
  else {
    double *h = domain->h;
    dvalue = sqrt(h[2] * h[2] + h[3] * h[3] + h[4] * h[4]);
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cellalpha()
{
  if (!domain->triclinic)
    dvalue = 90.0;
  else {

    // Cos(alpha) = (xy.xz + ly.yz)/(b.c)

    double *h = domain->h;
    double cosalpha = (h[5] * h[4] + h[1] * h[3]) /
        sqrt((h[1] * h[1] + h[5] * h[5]) * (h[2] * h[2] + h[3] * h[3] + h[4] * h[4]));
    dvalue = acos(cosalpha) * 180.0 / MY_PI;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cellbeta()
{
  if (!domain->triclinic)
    dvalue = 90.0;
  else {

    // Cos(beta) = xz/c

    double *h = domain->h;
    double cosbeta = h[4] / sqrt(h[2] * h[2] + h[3] * h[3] + h[4] * h[4]);
    dvalue = acos(cosbeta) * 180.0 / MY_PI;
  }
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_cellgamma()
{
  if (!domain->triclinic)
    dvalue = 90.0;
  else {

    // Cos(gamma) = xy/b

    double *h = domain->h;
    double cosgamma = h[5] / sqrt(h[1] * h[1] + h[5] * h[5]);
    dvalue = acos(cosgamma) * 180.0 / MY_PI;
  }
}

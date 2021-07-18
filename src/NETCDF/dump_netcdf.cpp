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
   Contributing author: Lars Pastewka (University of Freiburg)
------------------------------------------------------------------------- */

#if defined(LMP_HAS_NETCDF)

#include "dump_netcdf.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "thermo.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <netcdf.h>

using namespace LAMMPS_NS;
using namespace MathConst;

static const char NC_FRAME_STR[]         = "frame";
static const char NC_SPATIAL_STR[]       = "spatial";
static const char NC_VOIGT_STR[]         = "Voigt";
static const char NC_ATOM_STR[]          = "atom";
static const char NC_CELL_SPATIAL_STR[]  = "cell_spatial";
static const char NC_CELL_ANGULAR_STR[]  = "cell_angular";
static const char NC_LABEL_STR[]         = "label";

static const char NC_TIME_STR[]          = "time";
static const char NC_CELL_ORIGIN_STR[]   = "cell_origin";
static const char NC_CELL_LENGTHS_STR[]  = "cell_lengths";
static const char NC_CELL_ANGLES_STR[]   = "cell_angles";

static const char NC_UNITS_STR[]         = "units";
static const char NC_SCALE_FACTOR_STR[]  = "scale_factor";

static constexpr int THIS_IS_A_FIX      = -1;
static constexpr int THIS_IS_A_COMPUTE  = -2;
static constexpr int THIS_IS_A_VARIABLE = -3;
static constexpr int THIS_IS_A_BIGINT   = -4;

/* ---------------------------------------------------------------------- */

#define NCERR(x) ncerr(x, nullptr, __LINE__)
#define NCERRX(x, descr) ncerr(x, descr, __LINE__)
#if !defined(NC_64BIT_DATA)
#define NC_64BIT_DATA NC_64BIT_OFFSET
#endif

/* ---------------------------------------------------------------------- */

DumpNetCDF::DumpNetCDF(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  // arrays for data rearrangement

  sort_flag = 1;
  sortcol = 0;
  binary = 1;
  flush_flag = 0;

  if (multiproc)
    error->all(FLERR,"Multi-processor writes are not supported.");
  if (append_flag && multifile)
    error->all(FLERR,"Cannot append when writing to multiple files.");

  perat = new nc_perat_t[nfield];

  for (int i = 0; i < nfield; i++) {
    perat[i].dims = 0;
  }

  n_perat = 0;
  for (int i = 0; i < nfield; i++) {
    int idim = 0;
    int ndims = 1;
    std::string mangled = earg[i];
    bool constant = false;

    // name mangling
    // in the AMBER specification
    if ((mangled == "x") || (mangled == "y") || (mangled == "z")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      mangled = "coordinates";
    } else if ((mangled == "vx") || (mangled == "vy") || (mangled == "vz")) {
      idim = mangled[1] - 'x';
      ndims = 3;
      mangled = "velocities";
    } else if ((mangled == "xs") || (mangled == "ys") || (mangled == "zs")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      mangled = "scaled_coordinates";
    } else if ((mangled == "xu") || (mangled == "yu") || (mangled == "zu")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      mangled = "unwrapped_coordinates";
    } else if ((mangled == "fx") || (mangled == "fy") || (mangled == "fz")) {
      idim = mangled[1] - 'x';
      ndims = 3;
      mangled = "forces";
    } else if ((mangled == "mux") || (mangled == "muy") || (mangled == "muz")) {
      idim = mangled[2] - 'x';
      ndims = 3;
      mangled = "mu";
    } else if (utils::strmatch(mangled, "^c_")) {
      std::size_t found = mangled.find('[');
      if (found != std::string::npos) {
        if (mangled.find(']',found) == std::string::npos)
          error->all(FLERR,"Missing ']' in dump command");
        idim = mangled[found+1] - '1';
        mangled = mangled.substr(0,found);
        ndims = THIS_IS_A_COMPUTE;
      }
    } else if (utils::strmatch(mangled, "^f_")) {
      std::size_t found = mangled.find('[');
      if (found != std::string::npos) {
        if (mangled.find(']',found) == std::string::npos)
          error->all(FLERR,"Missing ']' in dump command");
        idim = mangled[found+1] - '1';
        mangled = mangled.substr(0,found);
        ndims = THIS_IS_A_FIX;
      }
    } else if (utils::strmatch(mangled, "^v_")) {
      idim = 0;
      ndims = THIS_IS_A_VARIABLE;
    }

    // find mangled name
    int inc = -1;
    for (int j = 0; j < n_perat && inc < 0; j++) {
      if (mangled == perat[j].name) {
        inc = j;
      }
    }

    if (inc < 0) {
      // this has not yet been defined
      inc = n_perat;
      perat[inc].dims = ndims;
      if (ndims < 0) ndims = DUMP_NC_MAX_DIMS;
      for (int j = 0; j < DUMP_NC_MAX_DIMS; j++) {
        perat[inc].field[j] = -1;
      }
      strncpy(perat[inc].name, mangled.c_str(), NC_FIELD_NAME_MAX);
      n_perat++;
    }

    perat[inc].constant = constant;
    perat[inc].ndumped = 0;
    perat[inc].field[idim] = i;
  }

  n_buffer = 0;
  int_buffer = nullptr;
  double_buffer = nullptr;

  double_precision = false;

  thermo = false;
  thermovar = nullptr;

  framei = 0;
}

/* ---------------------------------------------------------------------- */

DumpNetCDF::~DumpNetCDF()
{
  closefile();

  delete[] perat;
  if (thermovar) delete[] thermovar;

  if (int_buffer) memory->sfree(int_buffer);
  if (double_buffer) memory->sfree(double_buffer);
}

/* ---------------------------------------------------------------------- */

void DumpNetCDF::openfile()
{
  char *filecurrent = filename;
  if (multifile && !singlefile_opened) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s", filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }

  if (thermo && !singlefile_opened) {
    if (thermovar)  delete[] thermovar;
    thermovar = new int[output->thermo->nfield];
  }

  // now the computes and fixes have been initialized, so we can query
  // for the size of vector quantities
  for (int i = 0; i < n_perat; i++) {
    if (perat[i].dims == THIS_IS_A_COMPUTE) {
      int j = -1;
      for (int k = 0; k < DUMP_NC_MAX_DIMS; k++) {
        if (perat[i].field[k] >= 0) {
          j = field2index[perat[i].field[0]];
        }
      }
      if (j < 0)
        error->all(FLERR,"Internal error.");
      if (!compute[j]->peratom_flag)
        error->all(FLERR,"compute does not provide per atom data");
      perat[i].dims = compute[j]->size_peratom_cols;
      if (perat[i].dims > DUMP_NC_MAX_DIMS)
        error->all(FLERR,"perat[i].dims > DUMP_NC_MAX_DIMS");
    } else if (perat[i].dims == THIS_IS_A_FIX) {
      int j = -1;
      for (int k = 0; k < DUMP_NC_MAX_DIMS; k++) {
        if (perat[i].field[k] >= 0) {
          j = field2index[perat[i].field[0]];
        }
      }
      if (j < 0)
        error->all(FLERR,"Internal error.");
      if (!fix[j]->peratom_flag)
        error->all(FLERR,"fix does not provide per atom data");
      perat[i].dims = fix[j]->size_peratom_cols;
      if (perat[i].dims > DUMP_NC_MAX_DIMS)
        error->all(FLERR,"perat[i].dims > DUMP_NC_MAX_DIMS");
    } else if (perat[i].dims == THIS_IS_A_VARIABLE) {
      error->all(FLERR,"Dump netcdf currently does not support dumping variables");
    }
  }

  // get total number of atoms
  ntotalgr = group->count(igroup);
  for (int i = 0; i < DUMP_NC_MAX_DIMS; i++) {
    vector_dim[i] = -1;
  }

  if (filewriter) {
    if (append_flag && !multifile) {
      // Fixme! Perform checks if dimensions and variables conform with
      // data structure standard.
      if (not utils::file_is_readable(filecurrent))
        error->all(FLERR, "cannot append to non-existent file {}",filecurrent);

      if (singlefile_opened) return;
      singlefile_opened = 1;

      NCERRX( nc_open(filecurrent, NC_WRITE, &ncid), filecurrent );

      // dimensions
      NCERRX( nc_inq_dimid(ncid, NC_FRAME_STR, &frame_dim), NC_FRAME_STR );
      NCERRX( nc_inq_dimid(ncid, NC_ATOM_STR, &atom_dim), NC_ATOM_STR );
      NCERRX( nc_inq_dimid(ncid, NC_CELL_SPATIAL_STR, &cell_spatial_dim), NC_CELL_SPATIAL_STR );
      NCERRX( nc_inq_dimid(ncid, NC_CELL_ANGULAR_STR, &cell_angular_dim), NC_CELL_ANGULAR_STR );
      NCERRX( nc_inq_dimid(ncid, NC_LABEL_STR, &label_dim), NC_LABEL_STR );

      for (int i = 0; i < n_perat; i++) {
        int dims = perat[i].dims;
        if (vector_dim[dims] < 0) {
          char dimstr[1024];
          if (dims == 3) {
            strcpy(dimstr, NC_SPATIAL_STR);
          } else if (dims == 6) {
            strcpy(dimstr, NC_VOIGT_STR);
          } else {
            sprintf(dimstr, "vec%i", dims);
          }
          if (dims != 1) {
            NCERRX( nc_inq_dimid(ncid, dimstr, &vector_dim[dims]), dimstr );
          }
        }
      }

      // default variables
      NCERRX( nc_inq_varid(ncid, NC_SPATIAL_STR, &spatial_var), NC_SPATIAL_STR );
      NCERRX( nc_inq_varid(ncid, NC_CELL_SPATIAL_STR, &cell_spatial_var), NC_CELL_SPATIAL_STR);
      NCERRX( nc_inq_varid(ncid, NC_CELL_ANGULAR_STR, &cell_angular_var), NC_CELL_ANGULAR_STR);

      NCERRX( nc_inq_varid(ncid, NC_TIME_STR, &time_var), NC_TIME_STR );
      NCERRX( nc_inq_varid(ncid, NC_CELL_ORIGIN_STR, &cell_origin_var), NC_CELL_ORIGIN_STR );
      NCERRX( nc_inq_varid(ncid, NC_CELL_LENGTHS_STR, &cell_lengths_var), NC_CELL_LENGTHS_STR);
      NCERRX( nc_inq_varid(ncid, NC_CELL_ANGLES_STR, &cell_angles_var), NC_CELL_ANGLES_STR);

      for (int i = 0; i < n_perat; i++) {
        NCERRX( nc_inq_varid(ncid, perat[i].name, &perat[i].var), perat[i].name );
      }

      // perframe variables
      if (thermo) {
        Thermo *th = output->thermo;
        for (int i = 0; i < th->nfield; i++) {
          NCERRX( nc_inq_varid(ncid, th->keyword[i], &thermovar[i]), th->keyword[i] );
        }
      }

      size_t nframes;
      NCERR( nc_inq_dimlen(ncid, frame_dim, &nframes) );
      // framei == -1 means append to file, == -2 means override last frame
      // Note that in the input file this translates to 'yes', '-1', etc.

      if (framei <= 0) framei = nframes+framei+1;
      if (framei < 1)  framei = 1;
    } else {
      if (framei != 0 && !multifile)
        error->all(FLERR,"at keyword requires use of 'append yes'");

      int dims[NC_MAX_VAR_DIMS];
      size_t index[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
      double d[1];

      if (singlefile_opened) return;
      singlefile_opened = 1;

      NCERRX( nc_create(filecurrent, NC_64BIT_DATA, &ncid), filecurrent );

      // dimensions
      NCERRX( nc_def_dim(ncid, NC_FRAME_STR, NC_UNLIMITED, &frame_dim), NC_FRAME_STR );
      NCERRX( nc_def_dim(ncid, NC_ATOM_STR, ntotalgr, &atom_dim), NC_ATOM_STR );
      NCERRX( nc_def_dim(ncid, NC_CELL_SPATIAL_STR, 3, &cell_spatial_dim), NC_CELL_SPATIAL_STR );
      NCERRX( nc_def_dim(ncid, NC_CELL_ANGULAR_STR, 3, &cell_angular_dim), NC_CELL_ANGULAR_STR );
      NCERRX( nc_def_dim(ncid, NC_LABEL_STR, 10, &label_dim), NC_LABEL_STR );

      for (int i = 0; i < n_perat; i++) {
        int dim = perat[i].dims;
        if (vector_dim[dim] < 0) {
          char dimstr[1024];
          if (dim == 3) {
            strcpy(dimstr, NC_SPATIAL_STR);
          } else if (dim == 6) {
            strcpy(dimstr, NC_VOIGT_STR);
          } else {
            sprintf(dimstr, "vec%i", dim);
          }
          if (dim != 1) {
            NCERRX( nc_def_dim(ncid, dimstr, dim, &vector_dim[dim]), dimstr );
          }
        }
      }

      // default variables
      dims[0] = 0;
      NCERRX( nc_def_var(ncid, NC_SPATIAL_STR, NC_CHAR, 1, dims, &spatial_var), NC_SPATIAL_STR );
      NCERRX( nc_def_var(ncid, NC_CELL_SPATIAL_STR, NC_CHAR, 1, dims, &cell_spatial_var), NC_CELL_SPATIAL_STR );
      dims[0] = 0;
      dims[1] = label_dim;
      NCERRX( nc_def_var(ncid, NC_CELL_ANGULAR_STR, NC_CHAR, 2, dims, &cell_angular_var), NC_CELL_ANGULAR_STR );

      dims[0] = frame_dim;
      NCERRX( nc_def_var(ncid, NC_TIME_STR, NC_DOUBLE, 1, dims, &time_var), NC_TIME_STR);
      dims[0] = frame_dim;
      dims[1] = cell_spatial_dim;
      NCERRX( nc_def_var(ncid, NC_CELL_ORIGIN_STR, NC_DOUBLE, 2, dims, &cell_origin_var), NC_CELL_ORIGIN_STR );
      NCERRX( nc_def_var(ncid, NC_CELL_LENGTHS_STR, NC_DOUBLE, 2, dims, &cell_lengths_var), NC_CELL_LENGTHS_STR );
      dims[0] = frame_dim;
      dims[1] = cell_angular_dim;
      NCERRX( nc_def_var(ncid, NC_CELL_ANGLES_STR, NC_DOUBLE, 2, dims, &cell_angles_var), NC_CELL_ANGLES_STR );

      // variables specified in the input file
      dims[0] = frame_dim;
      dims[1] = atom_dim;
      dims[2] = 0;

      for (int i = 0; i < n_perat; i++) {
        nc_type xtype;

        // Type mangling
        if (vtype[perat[i].field[0]] == Dump::INT) {
          xtype = NC_INT;
        } else if (vtype[perat[i].field[0]] == Dump::BIGINT) {
          xtype = NC_INT64;
        } else if (vtype[perat[i].field[0]] == Dump::STRING) {
          error->all(FLERR,"Dump netcdf currently does not support dumping string properties");
        } else {
          if (double_precision)
            xtype = NC_DOUBLE;
          else
            xtype = NC_FLOAT;
        }

        if (perat[i].constant) {
          // this quantity will only be written once
          if (perat[i].dims == 1) {
            NCERRX( nc_def_var(ncid, perat[i].name, xtype, 1, dims+1, &perat[i].var), perat[i].name );
          } else {
            // this is a vector
            dims[1] = vector_dim[perat[i].dims];
            NCERRX( nc_def_var(ncid, perat[i].name, xtype, 2, dims+1, &perat[i].var), perat[i].name );
          }
        } else {
          if (perat[i].dims == 1) {
            NCERRX( nc_def_var(ncid, perat[i].name, xtype, 2, dims, &perat[i].var), perat[i].name );
          } else {
            // this is a vector
            dims[2] = vector_dim[perat[i].dims];
            NCERRX( nc_def_var(ncid, perat[i].name, xtype, 3, dims, &perat[i].var), perat[i].name );
          }
        }
      }

      // perframe variables
      if (thermo) {
        Thermo *th = output->thermo;
        for (int i = 0; i < th->nfield; i++) {
          if (th->vtype[i] == Thermo::FLOAT) {
            NCERRX( nc_def_var(ncid, th->keyword[i], NC_DOUBLE, 1, dims,
                               &thermovar[i]), th->keyword[i] );
          } else if (th->vtype[i] == Thermo::INT) {
            NCERRX( nc_def_var(ncid, th->keyword[i], NC_INT, 1, dims,
                               &thermovar[i]), th->keyword[i] );
          } else if (th->vtype[i] == Thermo::BIGINT) {
#if defined(LAMMPS_SMALLBIG) || defined(LAMMPS_BIGBIG)
            NCERRX( nc_def_var(ncid, th->keyword[i], NC_INT64, 1, dims,
                               &thermovar[i]), th->keyword[i] );
#else
            NCERRX( nc_def_var(ncid, th->keyword[i], NC_LONG, 1, dims,
                               &thermovar[i]), th->keyword[i] );
#endif
          }
        }
      }

      // attributes
      NCERR( nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 5, "AMBER") );
      NCERR( nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion", 3, "1.0") );

      NCERR( nc_put_att_text(ncid, NC_GLOBAL, "program", 6, "LAMMPS") );
      NCERR( nc_put_att_text(ncid, NC_GLOBAL, "programVersion",strlen(lmp->version), lmp->version) );

      // units
      if (!strcmp(update->unit_style, "lj")) {
        NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR, 2, "lj") );
        NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR, 2, "lj") );
        NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR, 2, "lj") );
      } else if (!strcmp(update->unit_style, "real")) {
        NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR, 11, "femtosecond") );
        NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR, 8, "Angstrom") );
        NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR, 8, "Angstrom") );
      } else if (!strcmp(update->unit_style, "metal")) {
        NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR, 10, "picosecond") );
        NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR, 8, "Angstrom") );
        NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR, 8, "Angstrom") );
      } else if (!strcmp(update->unit_style, "si")) {
        NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR, 6, "second") );
        NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR, 5, "meter") );
        NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR, 5, "meter") );
      } else if (!strcmp(update->unit_style, "cgs")) {
        NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR, 6, "second") );
        NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR, 10, "centimeter") );
        NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR, 10, "centimeter") );
      } else if (!strcmp(update->unit_style, "electron")) {
        NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR, 11, "femtosecond") );
        NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR, 4, "Bohr") );
        NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR, 4, "Bohr") );
      } else {
        error->all(FLERR,"Unsupported unit style: {}", update->unit_style);
      }

      NCERR( nc_put_att_text(ncid, cell_angles_var, NC_UNITS_STR,6, "degree") );

      d[0] = update->dt;
      NCERR( nc_put_att_double(ncid, time_var, NC_SCALE_FACTOR_STR,NC_DOUBLE, 1, d) );
      d[0] = 1.0;
      NCERR( nc_put_att_double(ncid, cell_origin_var, NC_SCALE_FACTOR_STR,NC_DOUBLE, 1, d) );
      d[0] = 1.0;
      NCERR( nc_put_att_double(ncid, cell_lengths_var, NC_SCALE_FACTOR_STR,NC_DOUBLE, 1, d) );

      /*
       * Finished with definition
       */

      NCERR( nc_enddef(ncid) );

      /*
       * Write label variables
       */

      NCERR( nc_put_var_text(ncid, spatial_var, "xyz") );
      NCERR( nc_put_var_text(ncid, cell_spatial_var, "abc") );
      index[0] = 0;
      index[1] = 0;
      count[0] = 1;
      count[1] = 5;
      NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "alpha") );
      index[0] = 1;
      count[1] = 4;
      NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "beta") );
      index[0] = 2;
      count[1] = 5;
      NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "gamma") );

      append_flag = 1;
      framei = 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpNetCDF::closefile()
{
  if (filewriter && singlefile_opened) {
    NCERR( nc_close(ncid) );
    singlefile_opened = 0;
    // write to next frame upon next open
    if (multifile)
      framei = 1;
    else {
      // append next time DumpNetCDF::openfile is called
      append_flag = 1;
      framei++;
    }
  }
}

/* ---------------------------------------------------------------------- */

template <typename T>
int nc_put_var1_bigint(int ncid, int varid, const size_t index[], const T* tp)
{
  return nc_put_var1_int(ncid, varid, index, tp);
}

template <>
int nc_put_var1_bigint<long>(int ncid, int varid, const size_t index[],
                        const long* tp)
{
  return nc_put_var1_long(ncid, varid, index, tp);
}

template <>
int nc_put_var1_bigint<long long>(int ncid, int varid, const size_t index[],
                             const long long* tp)
{
  return nc_put_var1_longlong(ncid, varid, index, tp);
}

template <typename T>
int nc_put_vara_bigint(int ncid, int varid, const size_t start[],
                       const size_t count[], const T* tp)
{
  return nc_put_vara_int(ncid, varid, start, count, tp);
}

template <>
int nc_put_vara_bigint<long>(int ncid, int varid, const size_t start[],
                             const size_t count[], const long* tp)
{
  return nc_put_vara_long(ncid, varid, start, count, tp);
}

template <>
int nc_put_vara_bigint<long long>(int ncid, int varid, const size_t start[],
                                  const size_t count[], const long long* tp)
{
  return nc_put_vara_longlong(ncid, varid, start, count, tp);
}

template <typename T>
int nc_put_vars_bigint(int ncid, int varid, const size_t start[],
                       const size_t count[], const ptrdiff_t stride[],
                       const T* tp)
{
  return nc_put_vars_int(ncid, varid, start, count, stride, tp);
}

template <>
int nc_put_vars_bigint<long>(int ncid, int varid, const size_t start[],
                             const size_t count[], const ptrdiff_t stride[],
                             const long* tp)
{
  return nc_put_vars_long(ncid, varid, start, count, stride, tp);
}

template <>
int nc_put_vars_bigint<long long>(int ncid, int varid, const size_t start[],
                                  const size_t count[],
                                  const ptrdiff_t stride[],
                                  const long long* tp)
{
  return nc_put_vars_longlong(ncid, varid, start, count, stride, tp);
}

void DumpNetCDF::write()
{
  // open file

  openfile();

  // need to write per-frame (global) properties here since they may come
  // from computes. write_header below is only called from the writing
  // processes, but modify->compute[j]->compute_* must be called from all
  // processes.

  size_t start[2];

  start[0] = framei-1;
  start[1] = 0;

  if (thermo) {
    Thermo *th = output->thermo;
    for (int i = 0; i < th->nfield; i++) {
      th->call_vfunc(i);
      if (filewriter) {
        if (th->vtype[i] == Thermo::FLOAT) {
          NCERRX( nc_put_var1_double(ncid, thermovar[i], start,
                                     &th->dvalue),
                  th->keyword[i] );
        } else if (th->vtype[i] == Thermo::INT) {
          NCERRX( nc_put_var1_int(ncid, thermovar[i], start, &th->ivalue),
                  th->keyword[i] );
        } else if (th->vtype[i] == Thermo::BIGINT) {
          NCERRX( nc_put_var1_bigint(ncid, thermovar[i], start, &th->bivalue),
                  th->keyword[i] );
        }
      }
    }
  }

  // call write of superclass

  Dump::write();

  // close file. this ensures data is flushed and mimized data corruption

  closefile();
}

/* ---------------------------------------------------------------------- */

void DumpNetCDF::write_header(bigint n)
{
  size_t start[2];

  start[0] = framei-1;
  start[1] = 0;

  if (filewriter) {
    size_t count[2];
    double time, cell_origin[3], cell_lengths[3], cell_angles[3];

    time = update->ntimestep;
    if (domain->triclinic == 0) {
      cell_origin[0] = domain->boxlo[0];
      cell_origin[1] = domain->boxlo[1];
      cell_origin[2] = domain->boxlo[2];

      cell_lengths[0] = domain->xprd;
      cell_lengths[1] = domain->yprd;
      cell_lengths[2] = domain->zprd;

      cell_angles[0] = 90;
      cell_angles[1] = 90;
      cell_angles[2] = 90;
    } else {
      double cosalpha, cosbeta, cosgamma;
      double *h = domain->h;

      cell_origin[0] = domain->boxlo[0];
      cell_origin[1] = domain->boxlo[1];
      cell_origin[2] = domain->boxlo[2];

      cell_lengths[0] = domain->xprd;
      cell_lengths[1] = sqrt(h[1]*h[1]+h[5]*h[5]);
      cell_lengths[2] = sqrt(h[2]*h[2]+h[3]*h[3]+h[4]*h[4]);

      cosalpha = (h[5]*h[4]+h[1]*h[3])/
          sqrt((h[1]*h[1]+h[5]*h[5])*(h[2]*h[2]+h[3]*h[3]+h[4]*h[4]));
      cosbeta = h[4]/sqrt(h[2]*h[2]+h[3]*h[3]+h[4]*h[4]);
      cosgamma = h[5]/sqrt(h[1]*h[1]+h[5]*h[5]);

      cell_angles[0] = acos(cosalpha)*180.0/MY_PI;
      cell_angles[1] = acos(cosbeta)*180.0/MY_PI;
      cell_angles[2] = acos(cosgamma)*180.0/MY_PI;
    }

    // Recent AMBER conventions say that non-periodic boundaries should have
    // 'cell_lengths' set to zero.
    for (int dim = 0; dim < 3; dim++) {
        if (!domain->periodicity[dim])
            cell_lengths[dim] = 0.0;
    }

    count[0] = 1;
    count[1] = 3;
    NCERR( nc_put_var1_double(ncid, time_var, start, &time) );
    NCERR( nc_put_vara_double(ncid, cell_origin_var, start, count, cell_origin) );
    NCERR( nc_put_vara_double(ncid, cell_lengths_var, start, count, cell_lengths) );
    NCERR( nc_put_vara_double(ncid, cell_angles_var, start, count, cell_angles) );
  }

  ndata = n;
  blocki = 0;
}


/* ----------------------------------------------------------------------
   write data lines to file in a block-by-block style
   write head of block (mass & element name) only if has atoms of the type
------------------------------------------------------------------------- */

void DumpNetCDF::write_data(int n, double *mybuf)
{
  size_t start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  ptrdiff_t stride[NC_MAX_VAR_DIMS];

  if (!int_buffer) {
    n_buffer = n;
    int_buffer = (bigint *)
      memory->smalloc(n*sizeof(bigint),"dump::int_buffer");
    double_buffer = (double *)
      memory->smalloc(n*sizeof(double),"dump::double_buffer");
  }

  if (n > n_buffer) {
    n_buffer = n;
    int_buffer = (bigint *)
      memory->srealloc(int_buffer, n*sizeof(bigint),"dump::int_buffer");
    double_buffer = (double *)
      memory->srealloc(double_buffer, n*sizeof(double),"dump::double_buffer");
  }

  start[0] = framei-1;
  start[1] = blocki;
  start[2] = 0;

  count[0] = 1;
  count[1] = n;
  count[2] = 1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 3;

  for (int i = 0; i < n_perat; i++) {
    int iaux = perat[i].field[0];

    if (vtype[iaux] == Dump::INT || vtype[iaux] == Dump::BIGINT) {
      // integers
      if (perat[i].dims > 1) {

        for (int idim = 0; idim < perat[i].dims; idim++) {
          iaux = perat[i].field[idim];

          if (iaux >= 0) {
            if (vtype[iaux] == Dump::INT) {
              for (int j = 0; j < n; j++, iaux+=size_one) {
                int_buffer[j] = static_cast<int>(mybuf[iaux]);
              }
            } else { // Dump::BIGINT
              for (int j = 0; j < n; j++, iaux+=size_one) {
                int_buffer[j] = static_cast<bigint>(mybuf[iaux]);
              }
            }

            start[2] = idim;

            if (perat[i].constant) {
              if (perat[i].ndumped < ntotalgr) {
                NCERR( nc_put_vars_bigint(ncid, perat[i].var,start+1, count+1, stride+1,int_buffer) );
                perat[i].ndumped += n;
              }
            } else
              NCERR( nc_put_vars_bigint(ncid, perat[i].var, start, count, stride,int_buffer) );
          }
        }
      } else {
        if (vtype[iaux] == Dump::INT) {
          for (int j = 0; j < n; j++, iaux+=size_one) {
              int_buffer[j] = static_cast<int>(mybuf[iaux]);
          }
        } else { // Dump::BIGINT
          for (int j = 0; j < n; j++, iaux+=size_one) {
              int_buffer[j] = static_cast<bigint>(mybuf[iaux]);
          }
        }

        if (perat[i].constant) {
          if (perat[i].ndumped < ntotalgr) {
            NCERR( nc_put_vara_bigint(ncid, perat[i].var, start+1, count+1,int_buffer) );
            perat[i].ndumped += n;
          }
        } else
          NCERR( nc_put_vara_bigint(ncid, perat[i].var, start, count,int_buffer) );
      }
    } else {
      // doubles
      if (perat[i].dims > 1) {

        for (int idim = 0; idim < perat[i].dims; idim++) {
          iaux = perat[i].field[idim];

          if (iaux >= 0) {
            for (int j = 0; j < n; j++, iaux+=size_one) {
                double_buffer[j] = mybuf[iaux];
            }

            start[2] = idim;

            if (perat[i].constant) {
              if (perat[i].ndumped < ntotalgr) {
                  NCERR( nc_put_vars_double(ncid, perat[i].var,
                                            start+1, count+1, stride+1,
                                            double_buffer) );
                  perat[i].ndumped += n;
              }
            } else
              NCERR( nc_put_vars_double(ncid, perat[i].var, start, count,
                                        stride, double_buffer) );
          }
        }
      } else {
        for (int j = 0; j < n; j++, iaux+=size_one) {
            double_buffer[j] = mybuf[iaux];
        }

        if (perat[i].constant) {
          if (perat[i].ndumped < ntotalgr) {
            NCERR( nc_put_vara_double(ncid, perat[i].var, start+1, count+1,
                                      double_buffer) );
            perat[i].ndumped += n;
          }
        } else
          NCERR( nc_put_vara_double(ncid, perat[i].var, start, count,
                                    double_buffer) );
      }
    }
  }

  blocki += n;
}

/* ---------------------------------------------------------------------- */

int DumpNetCDF::modify_param(int narg, char **arg)
{
  int iarg = 0;
  if (strcmp(arg[iarg],"double") == 0) {
    iarg++;
    if (iarg >= narg)
      error->all(FLERR,"expected 'yes' or 'no' after 'double' keyword.");
    if (strcmp(arg[iarg],"yes") == 0) {
      double_precision = true;
    } else if (strcmp(arg[iarg],"no") == 0) {
      double_precision = false;
    } else error->all(FLERR,"expected 'yes' or 'no' after 'double' keyword.");
    iarg++;
    return 2;
  } else if (strcmp(arg[iarg],"at") == 0) {
    iarg++;
    if (iarg >= narg)
      error->all(FLERR,"expected additional arg after 'at' keyword.");
    framei = utils::inumeric(FLERR,arg[iarg],false,lmp);
    if (framei == 0) error->all(FLERR,"frame 0 not allowed for 'at' keyword.");
    else if (framei < 0) framei--;
    iarg++;
    return 2;
  } else if (strcmp(arg[iarg],"thermo") == 0) {
    iarg++;
    if (iarg >= narg)
      error->all(FLERR,"expected 'yes' or 'no' after 'thermo' keyword.");
    if (strcmp(arg[iarg],"yes") == 0) {
      thermo = true;
    } else if (strcmp(arg[iarg],"no") == 0) {
      thermo = false;
    } else error->all(FLERR,"expected 'yes' or 'no' after 'thermo' keyword.");
    iarg++;
    return 2;
  } else return 0;
}

/* ---------------------------------------------------------------------- */

void DumpNetCDF::ncerr(int err, const char *descr, int line)
{
  if (err != NC_NOERR) {
    if (descr) error->one(FLERR,"NetCDF failed with error '{}' (while accessing '{}') "
                          " in line {} of {}.", nc_strerror(err), descr, line, __FILE__);
    else error->one(FLERR,"NetCDF failed with error '{}' in line {} of {}.",
                    nc_strerror(err), line, __FILE__);
  }
}

#endif /* defined(LMP_HAS_NETCDF) */

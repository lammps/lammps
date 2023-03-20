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
   Contributing author: Lars Pastewka (University of Freiburg)
------------------------------------------------------------------------- */

#if defined(LMP_HAS_PNETCDF)

#include "dump_netcdf_mpiio.h"
#include "netcdf_units.h"

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
#include <pnetcdf.h>

using namespace LAMMPS_NS;
using namespace MathConst;
using NetCDFUnits::Quantity;
using NetCDFUnits::get_unit_for;
using NetCDFUnits::LMP_MAX_VAR_DIMS;

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

/* ---------------------------------------------------------------------- */

#define NCERR(x) ncerr(x, nullptr, __LINE__)
#define NCERRX(x, descr) ncerr(x, descr, __LINE__)
#if !defined(NC_64BIT_DATA)
#define NC_64BIT_DATA NC_64BIT_OFFSET
#endif

/* ---------------------------------------------------------------------- */

DumpNetCDFMPIIO::DumpNetCDFMPIIO(LAMMPS *lmp, int narg, char **arg) :
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
    int quantity = Quantity::UNKNOWN;

    // name mangling
    // in the AMBER specification
    if ((mangled == "x") || (mangled == "y") || (mangled == "z")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      mangled = "coordinates";
      quantity = Quantity::DISTANCE;
    } else if ((mangled == "vx") || (mangled == "vy") || (mangled == "vz")) {
      idim = mangled[1] - 'x';
      ndims = 3;
      mangled = "velocities";
      quantity = Quantity::VELOCITY;
    } else if ((mangled == "xs") || (mangled == "ys") || (mangled == "zs")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      mangled = "scaled_coordinates";
      // no unit for scaled coordinates
    } else if ((mangled == "xu") || (mangled == "yu") || (mangled == "zu")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      mangled = "unwrapped_coordinates";
      quantity = Quantity::DISTANCE;
    } else if ((mangled == "fx") || (mangled == "fy") || (mangled == "fz")) {
      idim = mangled[1] - 'x';
      ndims = 3;
      mangled = "forces";
      quantity = Quantity::FORCE;
    } else if ((mangled == "mux") || (mangled == "muy") || (mangled == "muz")) {
      idim = mangled[2] - 'x';
      ndims = 3;
      mangled = "mu";
      quantity = Quantity::DIPOLE_MOMENT;
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
      if (ndims < 0) ndims = DUMP_NC_MPIIO_MAX_DIMS;
      for (int j = 0; j < DUMP_NC_MPIIO_MAX_DIMS; j++) {
        perat[inc].field[j] = -1;
      }
      strncpy(perat[inc].name, mangled.c_str(), NC_MPIIO_FIELD_NAME_MAX-1);
      n_perat++;
    }

    perat[inc].field[idim] = i;
    perat[inc].quantity = quantity;
  }

  n_buffer = 0;
  int_buffer = nullptr;
  double_buffer = nullptr;

  type_nc_real = NC_FLOAT;

  thermo = false;
  thermovar = nullptr;

  framei = 0;
}

/* ---------------------------------------------------------------------- */

DumpNetCDFMPIIO::~DumpNetCDFMPIIO()
{
  closefile();

  delete[] perat;
  delete[] thermovar;

  if (int_buffer) memory->sfree(int_buffer);
  if (double_buffer) memory->sfree(double_buffer);
}

/* ---------------------------------------------------------------------- */

void DumpNetCDFMPIIO::openfile()
{
  std::string filecurrent = filename;
  if (multifile && !singlefile_opened)
    filecurrent = utils::star_subst(filename, update->ntimestep, padflag);

  if (thermo && !singlefile_opened) {
    delete[] thermovar;
    thermovar = new int[output->thermo->nfield];
  }

  // now the computes and fixes have been initialized, so we can query
  // for the size of vector quantities
  for (int i = 0; i < n_perat; i++) {
    if (perat[i].dims == THIS_IS_A_COMPUTE) {
      int j = -1;
      for (int k = 0; k < DUMP_NC_MPIIO_MAX_DIMS; k++) {
        if (perat[i].field[k] >= 0) {
          j = field2index[perat[i].field[0]];
        }
      }
      if (j < 0)
        error->all(FLERR,"Internal error.");
      if (!compute[j]->peratom_flag)
        error->all(FLERR,"compute does not provide per atom data");
      perat[i].dims = compute[j]->size_peratom_cols;
      if (perat[i].dims > DUMP_NC_MPIIO_MAX_DIMS)
        error->all(FLERR,"perat[i].dims > DUMP_NC_MPIIO_MAX_DIMS");
    } else if (perat[i].dims == THIS_IS_A_FIX) {
      int j = -1;
      for (int k = 0; k < DUMP_NC_MPIIO_MAX_DIMS; k++) {
        if (perat[i].field[k] >= 0) {
          j = field2index[perat[i].field[0]];
        }
      }
      if (j < 0)
        error->all(FLERR,"Internal error.");
      if (!fix[j]->peratom_flag)
        error->all(FLERR,"fix does not provide per atom data");
      perat[i].dims = fix[j]->size_peratom_cols;
      if (perat[i].dims > DUMP_NC_MPIIO_MAX_DIMS)
        error->all(FLERR,"perat[i].dims > DUMP_NC_MPIIO_MAX_DIMS");
    } else if (perat[i].dims == THIS_IS_A_VARIABLE) {
      error->all(FLERR,"Dump netcdf/mpiio currently does not support dumping variables");
    }
  }

  // get total number of atoms
  ntotalgr = group->count(igroup);
  for (int i = 0; i < DUMP_NC_MPIIO_MAX_DIMS; i++) {
    vector_dim[i] = -1;
  }

  if (append_flag && !multifile) {
    // Fixme! Perform checks if dimensions and variables conform with
    // data structure standard.
    if (!platform::file_is_readable(filecurrent))
      error->all(FLERR, "cannot append to non-existent file {}", filecurrent);

    if (singlefile_opened) return;
    singlefile_opened = 1;

    NCERRX( ncmpi_open(world, filecurrent.c_str(), NC_WRITE, MPI_INFO_NULL, &ncid),
            filecurrent.c_str() );

    // dimensions
    NCERRX( ncmpi_inq_dimid(ncid, NC_FRAME_STR, &frame_dim), NC_FRAME_STR );
    NCERRX( ncmpi_inq_dimid(ncid, NC_ATOM_STR, &atom_dim), NC_ATOM_STR );
    NCERRX( ncmpi_inq_dimid(ncid, NC_CELL_SPATIAL_STR, &cell_spatial_dim), NC_CELL_SPATIAL_STR );
    NCERRX( ncmpi_inq_dimid(ncid, NC_CELL_ANGULAR_STR, &cell_angular_dim), NC_CELL_ANGULAR_STR );
    NCERRX( ncmpi_inq_dimid(ncid, NC_LABEL_STR, &label_dim), NC_LABEL_STR );

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
          NCERRX( ncmpi_inq_dimid(ncid, dimstr, &vector_dim[dim]), dimstr );
        }
      }
    }

    // default variables
    NCERRX( ncmpi_inq_varid(ncid, NC_SPATIAL_STR, &spatial_var), NC_SPATIAL_STR );
    NCERRX( ncmpi_inq_varid(ncid, NC_CELL_SPATIAL_STR, &cell_spatial_var), NC_CELL_SPATIAL_STR);
    NCERRX( ncmpi_inq_varid(ncid, NC_CELL_ANGULAR_STR, &cell_angular_var), NC_CELL_ANGULAR_STR);

    NCERRX( ncmpi_inq_varid(ncid, NC_TIME_STR, &time_var), NC_TIME_STR );
    NCERRX( ncmpi_inq_varid(ncid, NC_CELL_ORIGIN_STR, &cell_origin_var), NC_CELL_ORIGIN_STR );
    NCERRX( ncmpi_inq_varid(ncid, NC_CELL_LENGTHS_STR, &cell_lengths_var), NC_CELL_LENGTHS_STR);
    NCERRX( ncmpi_inq_varid(ncid, NC_CELL_ANGLES_STR, &cell_angles_var), NC_CELL_ANGLES_STR);

    // variables specified in the input file
    for (int i = 0; i < n_perat; i++) {
      NCERRX( ncmpi_inq_varid(ncid, perat[i].name, &perat[i].var), perat[i].name );
    }

    // perframe variables
    if (thermo) {
      Thermo *th = output->thermo;
      for (int i = 0; i < th->nfield; i++) {
        NCERRX( ncmpi_inq_varid(ncid, th->keyword[i].c_str(), &thermovar[i]), th->keyword[i].c_str() );
      }
    }

    MPI_Offset nframes;
    NCERR( ncmpi_inq_dimlen(ncid, frame_dim, &nframes) );
    // framei == -1 means append to file, == -2 means override last frame
    // Note that in the input file this translates to 'yes', '-1', etc.
    if (framei <= 0) framei = nframes+framei+1;
    if (framei < 1)  framei = 1;
  } else {
    if (framei != 0 && !multifile)
      error->all(FLERR,"at keyword requires use of 'append yes'");

    int dims[LMP_MAX_VAR_DIMS];
    MPI_Offset index[LMP_MAX_VAR_DIMS], count[LMP_MAX_VAR_DIMS];

    if (singlefile_opened) return;
    singlefile_opened = 1;

    NCERRX( ncmpi_create(world, filecurrent.c_str(), NC_64BIT_DATA, MPI_INFO_NULL, &ncid),
            filecurrent.c_str() );

    // dimensions
    NCERRX( ncmpi_def_dim(ncid, NC_FRAME_STR, NC_UNLIMITED, &frame_dim), NC_FRAME_STR );
    NCERRX( ncmpi_def_dim(ncid, NC_ATOM_STR, ntotalgr, &atom_dim), NC_ATOM_STR );
    NCERRX( ncmpi_def_dim(ncid, NC_CELL_SPATIAL_STR, 3, &cell_spatial_dim), NC_CELL_SPATIAL_STR );
    NCERRX( ncmpi_def_dim(ncid, NC_CELL_ANGULAR_STR, 3, &cell_angular_dim), NC_CELL_ANGULAR_STR );
    NCERRX( ncmpi_def_dim(ncid, NC_LABEL_STR, 10, &label_dim), NC_LABEL_STR );

    if (vector_dim[3] < 0)
      NCERRX( ncmpi_def_dim(ncid, NC_SPATIAL_STR, 3, &vector_dim[3]), NC_SPATIAL_STR );
    if (vector_dim[6] < 0)
      NCERRX( ncmpi_def_dim(ncid, NC_VOIGT_STR, 6, &vector_dim[6]), NC_VOIGT_STR );

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
          NCERRX( ncmpi_def_dim(ncid, dimstr, dim, &vector_dim[dim]), dimstr );
        }
      }
    }

    // default variables
    dims[0] = vector_dim[3];
    NCERRX( ncmpi_def_var(ncid, NC_SPATIAL_STR, NC_CHAR, 1, dims, &spatial_var), NC_SPATIAL_STR );
    NCERRX( ncmpi_def_var(ncid, NC_CELL_SPATIAL_STR, NC_CHAR, 1, dims, &cell_spatial_var), NC_CELL_SPATIAL_STR );
    dims[0] = vector_dim[3];
    dims[1] = label_dim;
    NCERRX( ncmpi_def_var(ncid, NC_CELL_ANGULAR_STR, NC_CHAR, 2, dims, &cell_angular_var), NC_CELL_ANGULAR_STR );
    dims[0] = frame_dim;
    NCERRX( ncmpi_def_var(ncid, NC_TIME_STR, type_nc_real, 1, dims, &time_var), NC_TIME_STR);
    dims[0] = frame_dim;
    dims[1] = cell_spatial_dim;
    NCERRX( ncmpi_def_var(ncid, NC_CELL_ORIGIN_STR, type_nc_real, 2, dims, &cell_origin_var), NC_CELL_ORIGIN_STR );
    NCERRX( ncmpi_def_var(ncid, NC_CELL_LENGTHS_STR, type_nc_real, 2, dims, &cell_lengths_var), NC_CELL_LENGTHS_STR );
    dims[0] = frame_dim;
    dims[1] = cell_angular_dim;
    NCERRX( ncmpi_def_var(ncid, NC_CELL_ANGLES_STR, type_nc_real, 2, dims, &cell_angles_var), NC_CELL_ANGLES_STR );

    // variables specified in the input file
    dims[0] = frame_dim;
    dims[1] = atom_dim;
    dims[2] = vector_dim[3];

    for (int i = 0; i < n_perat; i++) {
      nc_type xtype;

      // Type mangling
      if (vtype[perat[i].field[0]] == Dump::INT) {
        xtype = NC_INT;
      } else if (vtype[perat[i].field[0]] == Dump::BIGINT) {
        xtype = NC_INT64;
      } else {
        xtype = type_nc_real;
      }

      if (perat[i].dims == 1) {
        NCERRX( ncmpi_def_var(ncid, perat[i].name, xtype, 2, dims, &perat[i].var), perat[i].name );
      } else {
        // this is a vector
        dims[2] = vector_dim[perat[i].dims];
        NCERRX( ncmpi_def_var(ncid, perat[i].name, xtype, 3, dims, &perat[i].var), perat[i].name );
      }

      std::string unit = get_unit_for(update->unit_style, perat[i].quantity, error);
      if (!unit.empty()) {
        NCERR( ncmpi_put_att_text(ncid, perat[i].var, NC_UNITS_STR, unit.size(), unit.c_str()) );
      }
    }

    // perframe variables
    if (thermo) {
      Thermo *th = output->thermo;
      for (int i = 0; i < th->nfield; i++) {
        if (th->vtype[i] == Thermo::FLOAT) {
          NCERRX( ncmpi_def_var(ncid, th->keyword[i].c_str(), type_nc_real, 1, dims, &thermovar[i]), th->keyword[i].c_str() );
        } else if (th->vtype[i] == Thermo::INT) {
          NCERRX( ncmpi_def_var(ncid, th->keyword[i].c_str(), NC_INT, 1, dims, &thermovar[i]), th->keyword[i].c_str() );
        } else if (th->vtype[i] == Thermo::BIGINT) {
#if defined(LAMMPS_SMALLBIG) || defined(LAMMPS_BIGBIG)
          NCERRX( ncmpi_def_var(ncid, th->keyword[i].c_str(), NC_INT64, 1, dims, &thermovar[i]), th->keyword[i].c_str() );
#else
          NCERRX( ncmpi_def_var(ncid, th->keyword[i].c_str(), NC_LONG, 1, dims, &thermovar[i]), th->keyword[i].c_str() );
#endif
        }
      }
    }

    // attributes
    NCERR( ncmpi_put_att_text(ncid, NC_GLOBAL, "Conventions", 5, "AMBER") );
    NCERR( ncmpi_put_att_text(ncid, NC_GLOBAL, "ConventionVersion", 3, "1.0") );

    NCERR( ncmpi_put_att_text(ncid, NC_GLOBAL, "program", 6, "LAMMPS") );
    NCERR( ncmpi_put_att_text(ncid, NC_GLOBAL, "programVersion", strlen(lmp->version), lmp->version) );

    // units & scale
    std::string unit = get_unit_for(update->unit_style, Quantity::TIME, error);
    NCERR( ncmpi_put_att_text(ncid, time_var, NC_UNITS_STR, unit.size(), unit.c_str()) );

    unit = get_unit_for(update->unit_style, Quantity::DISTANCE, error);
    NCERR( ncmpi_put_att_text(ncid, cell_origin_var, NC_UNITS_STR, unit.size(), unit.c_str()) );
    NCERR( ncmpi_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR, unit.size(), unit.c_str()) );

    NCERR( ncmpi_put_att_text(ncid, cell_angles_var, NC_UNITS_STR, 6, "degree") );

    float scale[1] = {static_cast<float>(update->dt)};
    NCERR( ncmpi_put_att_float(ncid, time_var, NC_SCALE_FACTOR_STR, NC_FLOAT, 1, scale) );

    /*
     * Finished with definition
     */

    NCERR( ncmpi_enddef(ncid) );

    /*
     * Write label variables
     */

    NCERR( ncmpi_begin_indep_data(ncid) );

    if (filewriter) {
      NCERR( ncmpi_put_var_text(ncid, spatial_var, "xyz") );
      NCERR( ncmpi_put_var_text(ncid, cell_spatial_var, "abc") );
      index[0] = 0;
      index[1] = 0;
      count[0] = 1;
      count[1] = 5;
      NCERR( ncmpi_put_vara_text(ncid, cell_angular_var, index, count, "alpha") );
      index[0] = 1;
      count[1] = 4;
      NCERR( ncmpi_put_vara_text(ncid, cell_angular_var, index, count, "beta") );
      index[0] = 2;
      count[1] = 5;
      NCERR( ncmpi_put_vara_text(ncid, cell_angular_var, index, count, "gamma") );
    }

    NCERR( ncmpi_end_indep_data(ncid) );

    append_flag = 1;
    framei = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DumpNetCDFMPIIO::closefile()
{
  if (singlefile_opened) {
    NCERR( ncmpi_close(ncid) );
    singlefile_opened = 0;
    // write to next frame upon next open
    if (multifile)
      framei = 1;
    else {
      // append next time DumpNetCDFMPIIO::openfile is called
      append_flag = 1;
      framei++;
    }
  }
}

/* ---------------------------------------------------------------------- */

template <typename T>
int ncmpi_put_var1_bigint(int ncid, int varid, const MPI_Offset index[], const T* tp)
{
  return ncmpi_put_var1_int(ncid, varid, index, tp);
}

template <>
int ncmpi_put_var1_bigint<long>(int ncid, int varid,
                                const MPI_Offset index[], const long* tp)
{
  return ncmpi_put_var1_long(ncid, varid, index, tp);
}

template <>
int ncmpi_put_var1_bigint<long long>(int ncid, int varid,
                                     const MPI_Offset index[], const long long* tp)
{
  return ncmpi_put_var1_longlong(ncid, varid, index, tp);
}

template <typename T>
int ncmpi_put_vara_bigint_all(int ncid, int varid, const MPI_Offset start[],
                              const MPI_Offset count[], const T* tp)
{
  return ncmpi_put_vara_int_all(ncid, varid, start, count, tp);
}

template <>
int ncmpi_put_vara_bigint_all<long>(int ncid, int varid, const MPI_Offset start[],
                                    const MPI_Offset count[], const long* tp)
{
  return ncmpi_put_vara_long_all(ncid, varid, start, count, tp);
}

template <>
int ncmpi_put_vara_bigint_all<long long>(int ncid, int varid, const MPI_Offset start[],
                                         const MPI_Offset count[], const long long* tp)
{
  return ncmpi_put_vara_longlong_all(ncid, varid, start, count, tp);
}

template <typename T>
int ncmpi_put_vars_bigint_all(int ncid, int varid, const MPI_Offset start[],
                              const MPI_Offset count[], const MPI_Offset stride[], const T* tp)
{
  return ncmpi_put_vars_int_all(ncid, varid, start, count, stride, tp);
}

template <>
int ncmpi_put_vars_bigint_all<long>(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                    const MPI_Offset stride[], const long* tp)
{
  return ncmpi_put_vars_long_all(ncid, varid, start, count, stride, tp);
}

template <>
int ncmpi_put_vars_bigint_all<long long>(int ncid, int varid, const MPI_Offset start[], const MPI_Offset count[],
                                         const MPI_Offset stride[], const long long* tp)
{
  return ncmpi_put_vars_longlong_all(ncid, varid, start, count, stride, tp);
}

void DumpNetCDFMPIIO::write()
{
  // open file

  openfile();

  // need to write per-frame (global) properties here since they may come
  // from computes. write_header below is only called from the writing
  // processes, but modify->compute[j]->compute_* must be called from all
  // processes.

  MPI_Offset start[2];

  start[0] = framei-1;
  start[1] = 0;

  NCERR( ncmpi_begin_indep_data(ncid) );

  if (thermo) {
    Thermo *th = output->thermo;
    for (int i = 0; i < th->nfield; i++) {
      th->call_vfunc(i);
      if (filewriter) {
        if (th->vtype[i] == Thermo::FLOAT) {
          NCERRX( ncmpi_put_var1_double(ncid, thermovar[i], start,
                                        &th->dvalue),
                  th->keyword[i].c_str() );
        } else if (th->vtype[i] == Thermo::INT) {
          NCERRX( ncmpi_put_var1_int(ncid, thermovar[i], start, &th->ivalue),
                  th->keyword[i].c_str() );
        } else if (th->vtype[i] == Thermo::BIGINT) {
          NCERRX( ncmpi_put_var1_bigint(ncid, thermovar[i], start, &th->bivalue),
                  th->keyword[i].c_str() );
        }
      }
    }
  }

 // write timestep header

  write_time_and_cell();

  NCERR( ncmpi_end_indep_data(ncid) );

  // nme = # of dump lines this proc contributes to dump

  nme = count();
  int *block_sizes = new int[comm->nprocs];
  MPI_Allgather(&nme, 1, MPI_INT, block_sizes, 1, MPI_INT, world);
  blocki = 0;
  for (int i = 0; i < comm->me; i++)  blocki += block_sizes[i];
  delete[] block_sizes;

  // ensure buf is sized for packing and communicating
  // use nme to ensure filewriter proc can receive info from others
  // limit nme*size_one to int since used as arg in MPI calls

  if (nme > maxbuf) {
    if ((bigint) nme * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // pack my data into buf

  pack(nullptr);

  // each process writes its data

  write_data(nme, buf);

  // close file. this ensures data is flushed and minimizes data corruption

  closefile();
}

/* ---------------------------------------------------------------------- */

void DumpNetCDFMPIIO::write_time_and_cell()
{
  MPI_Offset start[2];

  start[0] = framei-1;
  start[1] = 0;

  MPI_Offset count[2];
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
  if (filewriter) {
    NCERR( ncmpi_put_var1_double(ncid, time_var, start, &time) );
    NCERR( ncmpi_put_vara_double(ncid, cell_origin_var, start, count, cell_origin) );
    NCERR( ncmpi_put_vara_double(ncid, cell_lengths_var, start, count, cell_lengths) );
    NCERR( ncmpi_put_vara_double(ncid, cell_angles_var, start, count, cell_angles) );
  }
}


/* ----------------------------------------------------------------------
   write data lines to file in a block-by-block style
   write head of block (mass & element name) only if has atoms of the type
------------------------------------------------------------------------- */

void DumpNetCDFMPIIO::write_data(int n, double *mybuf)
{
  MPI_Offset start[LMP_MAX_VAR_DIMS], count[LMP_MAX_VAR_DIMS], stride[LMP_MAX_VAR_DIMS];

  if (!int_buffer) {
    n_buffer = std::max(1, n);
    int_buffer = (bigint *)
      memory->smalloc(n_buffer*sizeof(bigint),"dump::int_buffer");
    double_buffer = (double *)
      memory->smalloc(n_buffer*sizeof(double),"dump::double_buffer");
  }

  if (n > n_buffer) {
    n_buffer = std::max(1, n);
    int_buffer = (bigint *)
      memory->srealloc(int_buffer, n_buffer*sizeof(bigint),"dump::int_buffer");
    double_buffer = (double *)
      memory->srealloc(double_buffer, n_buffer*sizeof(double), "dump::double_buffer");
  }

  start[0] = framei-1;
  start[1] = blocki;
  start[2] = 0;

  if (n == 0) {
    /* If there is no data, we need to make sure the start values don't exceed
       dimension bounds. Just set them to zero. */
    start[1] = 0;
  }

  count[0] = 1;
  count[1] = n;
  count[2] = 1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 3;

  for (int i = 0; i < n_perat; i++) {
    int iaux = perat[i].field[0];
    if (iaux < 0 || iaux >= size_one)
      error->one(FLERR, "Internal error: name = {}, iaux = {}, size_one = {}", perat[i].name, iaux, size_one);

    if (vtype[iaux] == Dump::INT || vtype[iaux] == Dump::BIGINT) {
      // integers
      if (perat[i].dims > 1) {

        for (int idim = 0; idim < perat[i].dims; idim++) {
          iaux = perat[i].field[idim];

          if (iaux >= 0) {
            if (iaux >= size_one)
              error->one(FLERR, "Internal error: name = {}, iaux = {}, size_one = {}", perat[i].name, iaux, size_one);

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
            NCERRX( ncmpi_put_vars_bigint_all(ncid, perat[i].var, start, count,stride, int_buffer), perat[i].name );
          }
        }
      } else {
        for (int j = 0; j < n; j++, iaux+=size_one) {
            int_buffer[j] = mybuf[iaux];
        }

        NCERRX( ncmpi_put_vara_bigint_all(ncid, perat[i].var, start, count, int_buffer), perat[i].name );
      }
    } else {
      // doubles
      if (perat[i].dims > 1) {

        for (int idim = 0; idim < perat[i].dims; idim++) {
          iaux = perat[i].field[idim];

          if (iaux >= 0) {
            if (iaux >= size_one)
              error->one(FLERR, "Internal error: name = {}, iaux = {}, size_one = {}", perat[i].name, iaux, size_one);

            for (int j = 0; j < n; j++, iaux+=size_one) {
                double_buffer[j] = mybuf[iaux];
            }

            start[2] = idim;
            NCERRX( ncmpi_put_vars_double_all(ncid, perat[i].var, start, count,
                                              stride, double_buffer), perat[i].name );
          }
        }
      } else {
        for (int j = 0; j < n; j++, iaux+=size_one) {
            double_buffer[j] = mybuf[iaux];
        }

        NCERRX( ncmpi_put_vara_double_all(ncid, perat[i].var, start, count,
                                          double_buffer), perat[i].name );
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpNetCDFMPIIO::modify_param(int narg, char **arg)
{
  int iarg = 0;
  if (strcmp(arg[iarg],"double") == 0) {
    iarg++;
    if (iarg >= narg) error->all(FLERR,"expected 'yes' or 'no' after 'double' keyword.");

    if (utils::logical(FLERR,arg[iarg],false,lmp) == 1)
      type_nc_real = NC_DOUBLE;
    else
      type_nc_real = NC_FLOAT;

    iarg++;
    return 2;
  } else if (strcmp(arg[iarg],"at") == 0) {
    iarg++;
    if (iarg >= narg) error->all(FLERR,"expected additional arg after 'at' keyword.");
    framei = utils::inumeric(FLERR,arg[iarg],false,lmp);
    if (framei == 0) error->all(FLERR,"frame 0 not allowed for 'at' keyword.");
    else if (framei < 0) framei--;
    iarg++;
    return 2;
  } else if (strcmp(arg[iarg],"thermo") == 0) {
    iarg++;
    if (iarg >= narg) error->all(FLERR,"expected 'yes' or 'no' after 'thermo' keyword.");
    thermo = utils::logical(FLERR,arg[iarg],false,lmp) == 1;
    iarg++;
    return 2;
  } else return 0;
}

/* ---------------------------------------------------------------------- */

void DumpNetCDFMPIIO::ncerr(int err, const char *descr, int line)
{
  if (err != NC_NOERR) {
    if (descr) error->one(__FILE__, line, "NetCDF failed with error '{}' (while accessing '{}') ",
                          ncmpi_strerror(err), descr);
    else error->one(__FILE__, line,"NetCDF failed with error '{}'.", ncmpi_strerror(err));
  }
}

#endif /* defined(LMP_HAS_PNETCDF) */

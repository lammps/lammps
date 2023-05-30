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
   Contributing author: Norbert Podhorszki (Oak Ridge National Laboratory)
------------------------------------------------------------------------- */

#include "reader_adios.h"

#include "comm.h"
#include "error.h"
#include "memory.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "adios2.h"
#include "adios_common.h"

using namespace LAMMPS_NS;

// also in read_dump.cpp

enum { ID, TYPE, X, Y, Z, VX, VY, VZ, Q, IX, IY, IZ, FX, FY, FZ };
enum { UNSET, NOSCALE_NOWRAP, NOSCALE_WRAP, SCALE_NOWRAP, SCALE_WRAP };

#define SMALL 1.0e-6

namespace LAMMPS_NS {
class ReadADIOSInternal {

 public:
  ReadADIOSInternal() = default;
  ~ReadADIOSInternal() = default;

  // name of adios group, referrable in adios2_config.xml
  const std::string ioName = "read_dump";
  adios2::ADIOS *ad = nullptr;    // adios object
  adios2::IO io;                  // adios group of variables and attributes in this dump
  adios2::Engine fh;              // adios file/stream handle object
  // ADIOS input variables we need to change every step
  adios2::Variable<uint64_t> varNtimestep;
  adios2::Variable<uint64_t> varNatoms;
  adios2::Variable<double> varAtoms;
  // list of column names for the atom table
  // (individual list of 'columns' string)
  std::vector<std::string> columnNames;
  float timeout = 0.0;
};
}    // namespace LAMMPS_NS

/* ---------------------------------------------------------------------- */

ReaderADIOS::ReaderADIOS(LAMMPS *lmp) : Reader(lmp)
{
  fieldindex = nullptr;
  nAtoms = 0;
  nAtomsTotal = 0;
  atomOffset = 0;
  nstep = 0;
  nid = 0;
  me = comm->me;

  // create a default adios2_config.xml if it doesn't exist yet.
  FILE *cfgfp = fopen("adios2_config.xml", "r");
  if (!cfgfp) {
    cfgfp = fopen("adios2_config.xml", "w");
    if (cfgfp) fputs(default_config, cfgfp);
  }
  if (cfgfp) fclose(cfgfp);

  internal = new ReadADIOSInternal();
  try {
#if defined(MPI_STUBS)
    internal->ad = new adios2::ADIOS("adios2_config.xml");
#else
    internal->ad = new adios2::ADIOS("adios2_config.xml", world);
#endif
  } catch (std::ios_base::failure &e) {
    error->one(FLERR, "ADIOS initialization failed with error: {}", e.what());
  }

  /* Define the group holding all variables and attributes  */
  internal->io = internal->ad->DeclareIO(internal->ioName);
}

/* ---------------------------------------------------------------------- */

ReaderADIOS::~ReaderADIOS()
{
  if (me == 0) memory->destroy(fieldindex);
  internal->columnNames.clear();
  if (internal->fh) internal->fh.Close();
  delete internal->ad;
  delete internal;
}

/* ----------------------------------------------------------------------
   pass on settings to find and load the proper plugin
   Called by all processors.
------------------------------------------------------------------------- */
void ReaderADIOS::settings(int narg, char **arg)
{
  int idx = 0;
  while (idx < narg) {
    if (!strcmp(arg[idx], "timeout")) {
      if (idx + 1 < narg) {
        internal->timeout = std::stof(arg[idx + 1]);
        internal->io.SetParameter("OpenTimeoutSecs", arg[idx + 1]);
        ++idx;
      } else {
        error->one(FLERR, "Missing value for 'timeout' option for ADIOS read_dump command");
      }
    }
    ++idx;
  }
}

/* ----------------------------------------------------------------------
   try to open given file
   Every process must call this Collective operation
------------------------------------------------------------------------- */

void ReaderADIOS::open_file(const std::string &file)
{
  // close open file, if needed.
  if (internal->fh) internal->fh.Close();

  try {
#if defined(MPI_STUBS)
    internal->fh = internal->io.Open(file, adios2::Mode::Read);
#else
    internal->fh = internal->io.Open(file, adios2::Mode::Read, world);
#endif
  } catch (std::ios_base::failure &e) {
    error->one(FLERR, "Error opening file {}: {}", file, e.what());
  }
  if (!internal->fh) error->one(FLERR, "Cannot open file {} using ADIOS", file);
}

/* ----------------------------------------------------------------------
   close current file
   Every process must call this Collective operation
------------------------------------------------------------------------- */

void ReaderADIOS::close_file()
{
  // close open file, if needed.
  if (internal->fh) { internal->fh.Close(); }
}

/* ----------------------------------------------------------------------
   read and return time stamp from dump file
   if first read reaches end-of-file, return 1 so caller can open next file
   Called by all processors.
------------------------------------------------------------------------- */

int ReaderADIOS::read_time(bigint &ntimestep)
{
  adios2::StepStatus status = internal->fh.BeginStep(adios2::StepMode::Read, internal->timeout);

  switch (status) {
    case adios2::StepStatus::EndOfStream:
    case adios2::StepStatus::NotReady:
    case adios2::StepStatus::OtherError:
      return 1;
    default:
      break;
  }

  internal->varNtimestep = internal->io.InquireVariable<uint64_t>("ntimestep");

  if (!internal->varNtimestep)
    error->one(FLERR, "Did not find 'ntimestep' variable in ADIOS file {}", internal->fh.Name());

  ntimestep = static_cast<bigint>(internal->varNtimestep.Max());
  return 0;
}

/* ----------------------------------------------------------------------
   skip snapshot for given timestamp
    Called by all processors.
------------------------------------------------------------------------- */

void ReaderADIOS::skip()
{
  internal->fh.EndStep();
}

/* ----------------------------------------------------------------------
   read remaining header info:
     return natoms
     box bounds, triclinic (inferred), fieldflag (1 if any fields not found),
     xyz flags = from input scaleflag & wrapflag
   if fieldflag set:
     match Nfield fields to per-atom column labels
     allocate and set fieldindex = which column each field maps to
     fieldtype = X,VX,IZ etc
     fieldlabel = user-specified label or nullptr if use fieldtype default
   xyz flag = scaledflag if has fieldlabel name, else set by x,xs,xu,xsu
   only called by proc 0
------------------------------------------------------------------------- */

bigint ReaderADIOS::read_header(double box[3][3], int &boxinfo, int &triclinic, int fieldinfo,
                                int nfield, int *fieldtype, char **fieldlabel, int scaleflag,
                                int wrapflag, int &fieldflag, int &xflag, int &yflag, int &zflag)
{
  nid = 0;

  // signal that we have no box info at all so far.

  internal->varNatoms = internal->io.InquireVariable<uint64_t>("natoms");
  if (!internal->varNatoms)
    error->one(FLERR, "Did not find 'natoms' variable in ADIOS file {}", internal->fh.Name());

  /* nAtoms */
  nAtomsTotal = internal->varNatoms.Max();
  uint64_t rem = nAtomsTotal % comm->nprocs;
  nAtoms = nAtomsTotal / comm->nprocs;
  atomOffset = comm->me * nAtoms;
  if (comm->me < (int) rem) {
    ++nAtoms;
    atomOffset += comm->me;
  } else {
    atomOffset += rem;
  }

  /* triclinic */
  adios2::Attribute<int32_t> attTriclinic = internal->io.InquireAttribute<int32_t>("triclinic");
  if (!attTriclinic)
    error->one(FLERR, "Did not find 'triclinic' attribute in ADIOS file {}", internal->fh.Name());

  triclinic = attTriclinic.Data()[0];

  /* read Box */
  adios2::Variable<double> varBoxxlo = internal->io.InquireVariable<double>("boxxlo");
  adios2::Variable<double> varBoxxhi = internal->io.InquireVariable<double>("boxxhi");
  adios2::Variable<double> varBoxylo = internal->io.InquireVariable<double>("boxylo");
  adios2::Variable<double> varBoxyhi = internal->io.InquireVariable<double>("boxyhi");
  adios2::Variable<double> varBoxzlo = internal->io.InquireVariable<double>("boxzlo");
  adios2::Variable<double> varBoxzhi = internal->io.InquireVariable<double>("boxzhi");

  box[0][0] = varBoxxlo.Max();
  box[0][1] = varBoxxhi.Max();
  box[0][2] = 0.0;
  box[1][0] = varBoxylo.Max();
  box[1][1] = varBoxyhi.Max();
  box[1][2] = 0.0;
  box[2][0] = varBoxzlo.Max();
  box[2][1] = varBoxzhi.Max();
  box[2][2] = 0.0;

  if (triclinic) {
    adios2::Variable<double> varBoxxy = internal->io.InquireVariable<double>("boxxy");
    adios2::Variable<double> varBoxxz = internal->io.InquireVariable<double>("boxxz");
    adios2::Variable<double> varBoxyz = internal->io.InquireVariable<double>("boxyz");

    box[0][2] = varBoxxy.Max();
    box[1][2] = varBoxxz.Max();
    box[2][2] = varBoxyz.Max();
  }

  boxinfo = 1;

  // if no field info requested, just return

  if (!fieldinfo) return nAtoms;

  memory->create(fieldindex, nfield, "read_dump:fieldindex");

  /* Columns */
  adios2::Attribute<std::string> attColumns = internal->io.InquireAttribute<std::string>("columns");

  std::vector<std::string> labelVector = attColumns.Data();
  int nwords = labelVector.size();
  std::map<std::string, int> labels;
  for (int i = 0; i < nwords; ++i) { labels.emplace(labelVector[i], i); }

  int s_index, u_index, su_index;
  xflag = UNSET;
  yflag = UNSET;
  zflag = UNSET;

  // copy fieldtype list for supported fields

  for (int i = 0; i < nfield; i++) {
    if (fieldlabel[i]) {
      fieldindex[i] = find_label(fieldlabel[i], labels);
      if (fieldtype[i] == X)
        xflag = 2 * scaleflag + wrapflag + 1;
      else if (fieldtype[i] == Y)
        yflag = 2 * scaleflag + wrapflag + 1;
      else if (fieldtype[i] == Z)
        zflag = 2 * scaleflag + wrapflag + 1;
    }

    else if (fieldtype[i] == ID)
      fieldindex[i] = find_label("id", labels);
    else if (fieldtype[i] == TYPE)
      fieldindex[i] = find_label("type", labels);

    else if (fieldtype[i] == X) {
      fieldindex[i] = find_label("x", labels);
      xflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("xs", labels);
        u_index = find_label("xu", labels);
        su_index = find_label("xsu", labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          xflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          xflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          xflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == Y) {
      fieldindex[i] = find_label("y", labels);
      yflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("ys", labels);
        u_index = find_label("yu", labels);
        su_index = find_label("ysu", labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          yflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          yflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          yflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == Z) {
      fieldindex[i] = find_label("z", labels);
      zflag = NOSCALE_WRAP;
      if (fieldindex[i] < 0) {
        fieldindex[i] = nwords;
        s_index = find_label("zs", labels);
        u_index = find_label("zu", labels);
        su_index = find_label("zsu", labels);
        if (s_index >= 0 && s_index < fieldindex[i]) {
          fieldindex[i] = s_index;
          zflag = SCALE_WRAP;
        }
        if (u_index >= 0 && u_index < fieldindex[i]) {
          fieldindex[i] = u_index;
          zflag = NOSCALE_NOWRAP;
        }
        if (su_index >= 0 && su_index < fieldindex[i]) {
          fieldindex[i] = su_index;
          zflag = SCALE_NOWRAP;
        }
      }
      if (fieldindex[i] == nwords) fieldindex[i] = -1;

    } else if (fieldtype[i] == VX)
      fieldindex[i] = find_label("vx", labels);
    else if (fieldtype[i] == VY)
      fieldindex[i] = find_label("vy", labels);
    else if (fieldtype[i] == VZ)
      fieldindex[i] = find_label("vz", labels);

    else if (fieldtype[i] == FX)
      fieldindex[i] = find_label("fx", labels);
    else if (fieldtype[i] == FY)
      fieldindex[i] = find_label("fy", labels);
    else if (fieldtype[i] == FZ)
      fieldindex[i] = find_label("fz", labels);

    else if (fieldtype[i] == Q)
      fieldindex[i] = find_label("q", labels);

    else if (fieldtype[i] == IX)
      fieldindex[i] = find_label("ix", labels);
    else if (fieldtype[i] == IY)
      fieldindex[i] = find_label("iy", labels);
    else if (fieldtype[i] == IZ)
      fieldindex[i] = find_label("iz", labels);
  }

  // set fieldflag = -1 if any unfound fields

  fieldflag = 0;
  for (int i = 0; i < nfield; i++)
    if (fieldindex[i] < 0) fieldflag = -1;

  return nAtoms;
}

/* ----------------------------------------------------------------------
   read N atom lines from dump file
   stores appropriate values in fields array
   return 0 if success, 1 if error
   only called by proc 0
------------------------------------------------------------------------- */

void ReaderADIOS::read_atoms(int n, int nfield, double **fields)
{
  /* Read Atoms */
  /* This is the firsts (and last) read of array data, so we
   * call EndStep() here instead of PerformGets()
   */

  adios2::Variable<double> varAtoms = internal->io.InquireVariable<double>("atoms");

  if ((uint64_t) n != nAtoms)
    error->one(FLERR,
               "ReaderADIOS::read_atoms() expects 'n={}' equal to the number of "
               "atoms (={}) for process {} in ADIOS file {}.",
               n, nAtoms, comm->me, internal->fh.Name());

  size_t ncols = varAtoms.Count()[1];
  varAtoms.SetSelection({{atomOffset, 0}, {nAtoms, ncols}});

  std::vector<double> table;

  internal->fh.Get<double>(varAtoms, table);
  // EndStep or PerformGets required to make the read happen
  internal->fh.EndStep();

  size_t idx;
  for (uint64_t i = 0; i < nAtoms; i++) {
    idx = i * ncols;
    for (int m = 0; m < nfield; m++) { fields[i][m] = table[idx + fieldindex[m]]; }
  }
}

/* ----------------------------------------------------------------------
   match label to any of N labels
   return index of match or -1 if no match
------------------------------------------------------------------------- */

int ReaderADIOS::find_label(const std::string &label, const std::map<std::string, int> &labels)
{
  auto it = labels.find(label);
  if (it != labels.end()) { return it->second; }
  return -1;
}

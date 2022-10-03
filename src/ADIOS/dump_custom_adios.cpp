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
   Contributing author: Norbert Podhorszki (ORNL)
------------------------------------------------------------------------- */

#include "dump_custom_adios.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

#include "adios2.h"
#include "adios_common.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
class DumpCustomADIOSInternal {

 public:
  DumpCustomADIOSInternal() = default;
  ~DumpCustomADIOSInternal() = default;

  // name of adios group, referrable in adios2_config.xml
  const std::string ioName = "custom";
  adios2::ADIOS *ad = nullptr;    // adios object
  adios2::IO io;                  // adios group of variables and attributes in this dump
  adios2::Engine fh;              // adios file/stream handle object
  // one ADIOS output variable we need to change every step
  adios2::Variable<double> varAtoms;
  // list of column names for the atom table
  // (individual list of 'columns' string)
  std::vector<std::string> columnNames;
};
}    // namespace LAMMPS_NS

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::DumpCustomADIOS(LAMMPS *lmp, int narg, char **arg) : DumpCustom(lmp, narg, arg)
{
  // create a default adios2_config.xml if it doesn't exist yet.
  FILE *cfgfp = fopen("adios2_config.xml", "r");
  if (!cfgfp) {
    cfgfp = fopen("adios2_config.xml", "w");
    if (cfgfp) fputs(default_config, cfgfp);
  }
  if (cfgfp) fclose(cfgfp);

  internal = new DumpCustomADIOSInternal();
  try {
    internal->ad = new adios2::ADIOS("adios2_config.xml", world, adios2::DebugON);
  } catch (std::ios_base::failure &e) {
    error->all(FLERR, "ADIOS initialization failed with error: {}", e.what());
  }

  internal->columnNames.reserve(nfield);
  for (int i = 0; i < nfield; ++i) { internal->columnNames.emplace_back(earg[i]); }
}

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::~DumpCustomADIOS()
{
  internal->columnNames.clear();
  if (internal->fh) { internal->fh.Close(); }
  delete internal->ad;
  delete internal;
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::openfile()
{
  if (multifile) {
    // if one file per timestep, replace '*' with current timestep
    auto filecurrent = utils::star_subst(filename, update->ntimestep, padflag);
    internal->fh = internal->io.Open(filecurrent, adios2::Mode::Write, world);
    if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filecurrent);
  } else {
    if (!singlefile_opened) {
      internal->fh = internal->io.Open(filename, adios2::Mode::Write, world);
      if (!internal->fh) error->one(FLERR, "Cannot open dump file {}", filename);
      singlefile_opened = 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::write()
{
  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();

  // ntotal = total # of atoms in snapshot
  // atomOffset = sum of # of atoms up to this proc (exclusive prefix sum)

  bigint bnme = nme;
  MPI_Allreduce(&bnme, &ntotal, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  bigint atomOffset;    // sum of all atoms on processes 0..me-1
  MPI_Scan(&bnme, &atomOffset, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  atomOffset -= nme;    // exclusive prefix sum needed

  // Now we know the global size and the local subset size and offset
  // of the atoms table
  auto nAtomsGlobal = static_cast<size_t>(ntotal);
  auto startRow = static_cast<size_t>(atomOffset);
  auto nAtomsLocal = static_cast<size_t>(nme);
  auto nColumns = static_cast<size_t>(size_one);
  internal->varAtoms.SetShape({nAtomsGlobal, nColumns});
  internal->varAtoms.SetSelection({{startRow, 0}, {nAtomsLocal, nColumns}});

  // insure filewriter proc can receive everyone's info
  // limit nmax*size_one to int since used as arg in MPI_Rsend() below
  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nme > maxbuf) {
    if ((bigint) nme * size_one > MAXSMALLINT) error->all(FLERR, "Too much per-proc info for dump");
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf, (maxbuf * size_one), "dump:buf");
  }
  if (sort_flag && sortcol == 0 && nme > maxids) {
    maxids = nme;
    memory->destroy(ids);
    memory->create(ids, maxids, "dump:ids");
  }

  if (sort_flag && sortcol == 0)
    pack(ids);
  else
    pack(nullptr);
  if (sort_flag) sort();

  openfile();
  internal->fh.BeginStep();
  // write info on data as scalars (by me==0)
  if (me == 0) {
    internal->fh.Put<uint64_t>("ntimestep", update->ntimestep);
    internal->fh.Put<int>("nprocs", nprocs);

    internal->fh.Put<double>("boxxlo", boxxlo);
    internal->fh.Put<double>("boxxhi", boxxhi);
    internal->fh.Put<double>("boxylo", boxylo);
    internal->fh.Put<double>("boxyhi", boxyhi);
    internal->fh.Put<double>("boxzlo", boxzlo);
    internal->fh.Put<double>("boxzhi", boxzhi);

    if (domain->triclinic) {
      internal->fh.Put<double>("boxxy", boxxy);
      internal->fh.Put<double>("boxxz", boxxz);
      internal->fh.Put<double>("boxyz", boxyz);
    }
  }
  // Everyone needs to write scalar variables that are used as dimensions and
  // offsets of arrays
  internal->fh.Put<uint64_t>("natoms", ntotal);
  internal->fh.Put<int>("ncolumns", size_one);
  internal->fh.Put<uint64_t>("nme", bnme);
  internal->fh.Put<uint64_t>("offset", atomOffset);
  // now write the atoms
  internal->fh.Put<double>("atoms", buf);
  internal->fh.EndStep();    // I/O will happen now...

  if (multifile) { internal->fh.Close(); }
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::init_style()
{
  // assemble column string from defaults and user values

  delete[] columns;
  std::string combined;
  int icol = 0;
  for (const auto &item : utils::split_words(columns_default)) {
    if (combined.size()) combined += " ";
    if (keyword_user[icol].size())
      combined += keyword_user[icol];
    else
      combined += item;
    ++icol;
  }
  columns = utils::strdup(combined);

  // setup boundary string

  domain->boundary_string(boundstr);

  // remove % from filename since ADIOS always writes a global file with
  // data/metadata
  char *ptr = strchr(filename, '%');
  if (ptr) {
    while (*ptr) {
      ptr[0] = ptr[1];
      ++ptr;
    }
  }

  /* The next four loops are copied from dump_custom_mpiio, but nothing is
   * done with them.
   * It is unclear why we need them here.
   * For metadata, variable[] will be written out as an ADIOS attribute if
   * nvariable>0
   */
  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable
  for (int i = 0; i < ncompute; i++) {
    compute[i] = modify->get_compute_by_id(id_compute[i]);
    if (!compute[i])
      error->all(FLERR, "Could not find dump custom/adios compute ID {}", id_compute[i]);
  }

  for (int i = 0; i < nfix; i++) {
    fix[i] = modify->get_fix_by_id(id_fix[i]);
    if (!fix[i]) error->all(FLERR, "Could not find dump custom/adios fix ID {}", id_fix[i]);
    if (nevery % fix[i]->peratom_freq)
      error->all(FLERR, "dump custom/adios and fix {} with ID {} not computed at compatible times",
                 fix[i]->style, id_fix[i]);
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) error->all(FLERR, "Could not find dump custom/adios variable name");
    variable[i] = ivariable;
  }

  // set index and check validity of region
  if (idregion && !domain->get_region_by_id(idregion))
    error->all(FLERR, "Region {} for dump custom/adios does not exist", idregion);

  /* Define the group of variables for the atom style here since it's a fixed
   * set */
  internal->io = internal->ad->DeclareIO(internal->ioName);
  if (!internal->io.InConfigFile()) {
    // if not defined by user, we can change the default settings
    // BPFile is the default writer
    internal->io.SetEngine("BPFile");
    int num_aggregators = multiproc;
    if (num_aggregators == 0) num_aggregators = 1;
    auto nstreams = std::to_string(num_aggregators);
    internal->io.SetParameters({{"substreams", nstreams}});
    if (me == 0)
      utils::logmesg(lmp, "ADIOS method for {} is n-to-m (aggregation with {} writers)\n", filename,
                     nstreams);
  }

  internal->io.DefineVariable<uint64_t>("ntimestep");
  internal->io.DefineVariable<uint64_t>("natoms");

  internal->io.DefineVariable<int>("nprocs");
  internal->io.DefineVariable<int>("ncolumns");

  internal->io.DefineVariable<double>("boxxlo");
  internal->io.DefineVariable<double>("boxxhi");
  internal->io.DefineVariable<double>("boxylo");
  internal->io.DefineVariable<double>("boxyhi");
  internal->io.DefineVariable<double>("boxzlo");
  internal->io.DefineVariable<double>("boxzhi");

  internal->io.DefineVariable<double>("boxxy");
  internal->io.DefineVariable<double>("boxxz");
  internal->io.DefineVariable<double>("boxyz");

  internal->io.DefineAttribute<int>("triclinic", domain->triclinic);

  int *boundaryptr = reinterpret_cast<int *>(domain->boundary);
  internal->io.DefineAttribute<int>("boundary", boundaryptr, 6);

  auto nColumns = static_cast<size_t>(size_one);
  internal->io.DefineAttribute<std::string>("columns", internal->columnNames.data(), nColumns);
  internal->io.DefineAttribute<std::string>("columnstr", columns);
  internal->io.DefineAttribute<std::string>("boundarystr", boundstr);
  internal->io.DefineAttribute<std::string>("LAMMPS/dump_style", "custom");
  internal->io.DefineAttribute<std::string>("LAMMPS/version", lmp->version);
  internal->io.DefineAttribute<std::string>("LAMMPS/num_ver", std::to_string(lmp->num_ver));

  internal->io.DefineVariable<uint64_t>("nme",
                                        {adios2::LocalValueDim});    // local dimension variable
  internal->io.DefineVariable<uint64_t>("offset",
                                        {adios2::LocalValueDim});    // local dimension variable

  // atom table size is not known at the moment
  // it will be correctly defined at the moment of write
  size_t UnknownSizeYet = 1;
  internal->varAtoms = internal->io.DefineVariable<double>(
      "atoms", {UnknownSizeYet, nColumns}, {UnknownSizeYet, 0}, {UnknownSizeYet, nColumns});
}

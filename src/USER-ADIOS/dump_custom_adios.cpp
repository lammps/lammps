/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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
#include <cmath>
#include <cstring>
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include "adios2.h"

using namespace LAMMPS_NS;

#define MAX_TEXT_HEADER_SIZE 4096
#define DUMP_BUF_CHUNK_SIZE 16384
#define DUMP_BUF_INCREMENT_SIZE 4096

enum {
    ID,
    MOL,
    TYPE,
    ELEMENT,
    MASS,
    X,
    Y,
    Z,
    XS,
    YS,
    ZS,
    XSTRI,
    YSTRI,
    ZSTRI,
    XU,
    YU,
    ZU,
    XUTRI,
    YUTRI,
    ZUTRI,
    XSU,
    YSU,
    ZSU,
    XSUTRI,
    YSUTRI,
    ZSUTRI,
    IX,
    IY,
    IZ,
    VX,
    VY,
    VZ,
    FX,
    FY,
    FZ,
    Q,
    MUX,
    MUY,
    MUZ,
    MU,
    RADIUS,
    DIAMETER,
    OMEGAX,
    OMEGAY,
    OMEGAZ,
    ANGMOMX,
    ANGMOMY,
    ANGMOMZ,
    TQX,
    TQY,
    TQZ,
    SPIN,
    ERADIUS,
    ERVEL,
    ERFORCE,
    COMPUTE,
    FIX,
    VARIABLE
};
enum { LT, LE, GT, GE, EQ, NEQ };
enum { INT, DOUBLE, STRING, BIGINT }; // same as in DumpCustom

namespace LAMMPS_NS
{
class DumpCustomADIOSInternal
{

public:
    DumpCustomADIOSInternal(){};
    ~DumpCustomADIOSInternal() = default;

    // name of adios group, referrable in adios2_config.xml
    const std::string ioName = "custom";
    adios2::ADIOS *ad = nullptr; // adios object
    adios2::IO io;     // adios group of variables and attributes in this dump
    adios2::Engine fh; // adios file/stream handle object
    // one ADIOS output variable we need to change every step
    adios2::Variable<double> varAtoms;
    // list of column names for the atom table
    // (individual list of 'columns' string)
    std::vector<std::string> columnNames;
};
}

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::DumpCustomADIOS(LAMMPS *lmp, int narg, char **arg)
: DumpCustom(lmp, narg, arg)
{
    internal = new DumpCustomADIOSInternal();
    internal->ad =
        new adios2::ADIOS("adios2_config.xml", world, adios2::DebugON);

    // if (screen) fprintf(screen, "DumpCustomADIOS constructor: nvariable=%d
    // id_variable=%p, variables=%p, nfield=%d, earg=%p\n", nvariable,
    // id_variable, variable, nfield, earg);
    internal->columnNames.reserve(nfield);
    for (int i = 0; i < nfield; ++i) {
        internal->columnNames.push_back(earg[i]);
        // if (screen) fprintf(screen, "earg[%d] = '%s'\n", i, earg[i]);
    }
}

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::~DumpCustomADIOS()
{
    internal->columnNames.clear();
    if (internal->fh) {
        internal->fh.Close();
    }
    delete internal->ad;
    delete internal;
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::openfile()
{
    if (multifile) {
        // if one file per timestep, replace '*' with current timestep
        char *filestar = strdup(filename);
        char *filecurrent = new char[strlen(filestar) + 16];
        char *ptr = strchr(filestar, '*');
        *ptr = '\0';
        if (padflag == 0)
            sprintf(filecurrent, "%s" BIGINT_FORMAT "%s", filestar,
                    update->ntimestep, ptr + 1);
        else {
            char bif[8], pad[16];
            strcpy(bif, BIGINT_FORMAT);
            sprintf(pad, "%%s%%0%d%s%%s", padflag, &bif[1]);
            sprintf(filecurrent, pad, filestar, update->ntimestep, ptr + 1);
        }
        internal->fh =
            internal->io.Open(filecurrent, adios2::Mode::Write, world);
        if (!internal->fh) {
            char str[128];
            sprintf(str, "Cannot open dump file %s", filecurrent);
            error->one(FLERR, str);
        }
        free(filestar);
        delete[] filecurrent;
    } else {
        if (!singlefile_opened) {
            internal->fh =
                internal->io.Open(filename, adios2::Mode::Write, world);
            if (!internal->fh) {
                char str[128];
                sprintf(str, "Cannot open dump file %s", filename);
                error->one(FLERR, str);
            }
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

    bigint atomOffset; // sum of all atoms on processes 0..me-1
    MPI_Scan(&bnme, &atomOffset, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    atomOffset -= nme; // exclusive prefix sum needed

    // Now we know the global size and the local subset size and offset
    // of the atoms table
    size_t nAtomsGlobal = static_cast<size_t>(ntotal);
    size_t startRow = static_cast<size_t>(atomOffset);
    size_t nAtomsLocal = static_cast<size_t>(nme);
    size_t nColumns = static_cast<size_t>(size_one);
    internal->varAtoms.SetShape({nAtomsGlobal, nColumns});
    internal->varAtoms.SetSelection({{startRow, 0}, {nAtomsLocal, nColumns}});

    // insure filewriter proc can receive everyone's info
    // limit nmax*size_one to int since used as arg in MPI_Rsend() below
    // pack my data into buf
    // if sorting on IDs also request ID list from pack()
    // sort buf as needed

    if (nme > maxbuf) {
        if ((bigint)nme * size_one > MAXSMALLINT)
            error->all(FLERR, "Too much per-proc info for dump");
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
        pack(NULL);
    if (sort_flag)
        sort();

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
    internal->fh.EndStep(); // I/O will happen now...

    if (multifile) {
        internal->fh.Close();
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::init_style()
{

    // setup boundary string

    domain->boundary_string(boundstr);

    // remove % from filename since ADIOS always writes a global file with
    // data/metadata
    int len = strlen(filename);
    char *ptr = strchr(filename, '%');
    if (ptr) {
        *ptr = '\0';
        char *s = new char[len - 1];
        sprintf(s, "%s%s", filename, ptr + 1);
        strncpy(filename, s, len);
    }

    /* The next four loops are copied from dump_custom_mpiio, but nothing is
     * done with them.
     * It is unclear why we need them here.
     * For metadata, variable[] will be written out as an ADIOS attribute if
     * nvariable>0
     */
    // find current ptr for each compute,fix,variable
    // check that fix frequency is acceptable
    int icompute;
    for (int i = 0; i < ncompute; i++) {
        icompute = modify->find_compute(id_compute[i]);
        if (icompute < 0)
            error->all(FLERR, "Could not find dump custom compute ID");
        compute[i] = modify->compute[icompute];
    }

    int ifix;
    for (int i = 0; i < nfix; i++) {
        ifix = modify->find_fix(id_fix[i]);
        if (ifix < 0)
            error->all(FLERR, "Could not find dump custom fix ID");
        fix[i] = modify->fix[ifix];
        if (nevery % modify->fix[ifix]->peratom_freq)
            error->all(FLERR,
                       "Dump custom and fix not computed at compatible times");
    }

    int ivariable;
    for (int i = 0; i < nvariable; i++) {
        ivariable = input->variable->find(id_variable[i]);
        if (ivariable < 0)
            error->all(FLERR, "Could not find dump custom variable name");
        variable[i] = ivariable;
    }

    // set index and check validity of region
    if (iregion >= 0) {
        iregion = domain->find_region(idregion);
        if (iregion == -1)
            error->all(FLERR, "Region ID for dump custom does not exist");
    }

    /* Define the group of variables for the atom style here since it's a fixed
     * set */
    internal->io = internal->ad->DeclareIO(internal->ioName);
    if (!internal->io.InConfigFile()) {
        // if not defined by user, we can change the default settings
        // BPFile is the default writer
        internal->io.SetEngine("BPFile");
        int num_aggregators = multiproc;
        if (num_aggregators == 0)
            num_aggregators = 1;
        char nstreams[128];
        sprintf(nstreams, "%d", num_aggregators);
        internal->io.SetParameters({{"substreams", nstreams}});
        if (me == 0 && screen)
            fprintf(
                screen,
                "ADIOS method for %s is n-to-m (aggregation with %s writers)\n",
                filename, nstreams);
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

    size_t nColumns = static_cast<size_t>(size_one);
    internal->io.DefineAttribute<std::string>(
        "columns", internal->columnNames.data(), nColumns);
    internal->io.DefineAttribute<std::string>("columnstr", columns);
    internal->io.DefineAttribute<std::string>("boundarystr", boundstr);
    internal->io.DefineAttribute<std::string>("LAMMPS/dump_style", "atom");
    internal->io.DefineAttribute<std::string>("LAMMPS/version",
                                              universe->version);
    internal->io.DefineAttribute<std::string>("LAMMPS/num_ver",
                                              universe->num_ver);

    internal->io.DefineVariable<uint64_t>(
        "nme", {adios2::LocalValueDim}); // local dimension variable
    internal->io.DefineVariable<uint64_t>(
        "offset", {adios2::LocalValueDim}); // local dimension variable

    // atom table size is not known at the moment
    // it will be correctly defined at the moment of write
    size_t UnknownSizeYet = 1;
    internal->varAtoms = internal->io.DefineVariable<double>(
        "atoms", {UnknownSizeYet, nColumns}, {UnknownSizeYet, 0},
        {UnknownSizeYet, nColumns});
}

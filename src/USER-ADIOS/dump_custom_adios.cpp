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
   Contributing author: Paul Coffman (IBM)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "dump_custom_adios.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "universe.h"
#include "memory.h"
#include "error.h"
#include <stdlib.h>

using namespace LAMMPS_NS;

#define MAX_TEXT_HEADER_SIZE 4096
#define DUMP_BUF_CHUNK_SIZE 16384
#define DUMP_BUF_INCREMENT_SIZE 4096

enum{ID,MOL,TYPE,ELEMENT,MASS,
     X,Y,Z,XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q,MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,SPIN,ERADIUS,ERVEL,ERFORCE,
     COMPUTE,FIX,VARIABLE};
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,STRING,BIGINT};    // same as in DumpCustom

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::DumpCustomADIOS(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
    ad = new adios2::ADIOS("adios2_config.xml", world, adios2::DebugON);
    groupsize = 0;
    //if (screen) fprintf(screen, "DumpCustomADIOS constructor: nvariable=%d id_variable=%p, variables=%p, nfield=%d, earg=%p\n", nvariable, id_variable, variable, nfield, earg);
    columnNames.reserve(nfield);
    for (int i = 0; i < nfield; ++i) {
        columnNames.push_back(earg[i]);
        //if (screen) fprintf(screen, "earg[%d] = '%s'\n", i, earg[i]);
    }
}

/* ---------------------------------------------------------------------- */

DumpCustomADIOS::~DumpCustomADIOS()
{
    columnNames.clear();
    if (fh)
    {
        fh.Close();
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::openfile()
{
    if (multifile) {
        // if one file per timestep, replace '*' with current timestep
        char *filestar = strdup(filename);
        char *filecurrent = new char[strlen(filestar) + 16];
        char *ptr = strchr(filestar,'*');
        *ptr = '\0';
        if (padflag == 0)
            sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
                    filestar,update->ntimestep,ptr+1);
        else {
            char bif[8],pad[16];
            strcpy(bif,BIGINT_FORMAT);
            sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
            sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
        }
        fh = io.Open(filecurrent, adios2::Mode::Write, world);
        if (!fh) {
            char str[128];
            sprintf(str,"Cannot open dump file %s",filecurrent);
            error->one(FLERR,str);
        }
        free(filestar);
        delete [] filecurrent;
    }
    else
    {
        if (!singlefile_opened)
        {
            fh = io.Open(filename, adios2::Mode::Write, world);
            if (!fh) {
                char str[128];
                sprintf(str,"Cannot open dump file %s",filename);
                error->one(FLERR,str);
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
    MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

    bigint atomOffset; // sum of all atoms on processes 0..me-1
    MPI_Scan (&bnme, &atomOffset, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    atomOffset -= nme; // exclusive prefix sum needed

    // Now we know the global size and the local subset size and offset
    // of the atoms table 
    size_t nAtomsGlobal = static_cast<size_t>(ntotal);
    size_t startRow = static_cast<size_t>(atomOffset);
    size_t nAtomsLocal = static_cast<size_t>(nme);
    size_t nColumns = static_cast<size_t>(size_one);
    varAtoms.SetShape({nAtomsGlobal,nColumns});
    varAtoms.SetSelection({{startRow, 0}, {nAtomsLocal,nColumns}});

    // insure filewriter proc can receive everyone's info
    // limit nmax*size_one to int since used as arg in MPI_Rsend() below
    // pack my data into buf
    // if sorting on IDs also request ID list from pack()
    // sort buf as needed

    if (nme > maxbuf) {
        if ((bigint) nme * size_one > MAXSMALLINT)
            error->all(FLERR,"Too much per-proc info for dump");
        maxbuf = nme;
        memory->destroy(buf);
        memory->create(buf,(maxbuf*size_one),"dump:buf");
    }
    if (sort_flag && sortcol == 0 && nme > maxids) {
        maxids = nme;
        memory->destroy(ids);
        memory->create(ids,maxids,"dump:ids");
    }

    if (sort_flag && sortcol == 0) pack(ids);
    else pack(NULL);
    if (sort_flag) sort();

    // Calculate data size written by this process
    groupsize = nme * size_one * sizeof(double); // size of atoms data on this process
    groupsize += 3*sizeof(uint64_t) + 1*sizeof(int); // scalars written by each process
    if (me == 0) {
        groupsize += 1*sizeof(uint64_t) + 1*sizeof(int) + 6*sizeof(double); // scalars
        if (domain->triclinic) {
            groupsize += 3*sizeof(double); // boxxy, boxxz, boxyz
        }
    }

    openfile();
    fh.BeginStep();
    // write info on data as scalars (by me==0)
    if (me == 0) {
        fh.Put<uint64_t>("ntimestep",   update->ntimestep);
        fh.Put<int>("nprocs",      nprocs);

        fh.Put<double>("boxxlo", boxxlo);
        fh.Put<double>("boxxhi", boxxhi);
        fh.Put<double>("boxylo", boxylo);
        fh.Put<double>("boxyhi", boxyhi);
        fh.Put<double>("boxzlo", boxzlo);
        fh.Put<double>("boxzhi", boxzhi);

        if (domain->triclinic) {
            fh.Put<double>("boxxy", boxxy);
            fh.Put<double>("boxxz", boxxz);
            fh.Put<double>("boxyz", boxyz);
        }
    }
    // Everyone needs to write scalar variables that are used as dimensions and offsets of arrays
    fh.Put<uint64_t>("natoms",   ntotal);
    fh.Put<int>("ncolumns", size_one);
    fh.Put<uint64_t>("nme",      bnme);
    fh.Put<uint64_t>("offset",   atomOffset);
    // now write the atoms
    fh.Put<double>("atoms",    buf);
    fh.EndStep();// I/O will happen now...

    if (multifile)
    {
        fh.Close();
    }
}

/* ---------------------------------------------------------------------- */

void DumpCustomADIOS::init_style()
{

    // setup boundary string

    domain->boundary_string(boundstr);

    // remove % from filename since ADIOS always writes a global file with data/metadata
    int len = strlen(filename);
    char *ptr = strchr(filename,'%');
    if (ptr) {
        *ptr = '\0';
        char *s = new char[len-1];
        sprintf(s,"%s%s",filename,ptr+1);
        strncpy(filename,s,len);
    }

    /* The next four loops are copied from dump_custom_mpiio, but nothing is done with them.
     * It is unclear why we need them here.
     * For metadata, variable[] will be written out as an ADIOS attribute if nvariable>0
     */
    // find current ptr for each compute,fix,variable
    // check that fix frequency is acceptable
    int icompute;
    for (int i = 0; i < ncompute; i++) {
        icompute = modify->find_compute(id_compute[i]);
        if (icompute < 0) error->all(FLERR,"Could not find dump custom compute ID");
        compute[i] = modify->compute[icompute];
    }

    int ifix;
    for (int i = 0; i < nfix; i++) {
        ifix = modify->find_fix(id_fix[i]);
        if (ifix < 0) error->all(FLERR,"Could not find dump custom fix ID");
        fix[i] = modify->fix[ifix];
        if (nevery % modify->fix[ifix]->peratom_freq)
            error->all(FLERR,"Dump custom and fix not computed at compatible times");
    }

    int ivariable;
    for (int i = 0; i < nvariable; i++) {
        ivariable = input->variable->find(id_variable[i]);
        if (ivariable < 0)
            error->all(FLERR,"Could not find dump custom variable name");
        variable[i] = ivariable;
    }

    // set index and check validity of region
    if (iregion >= 0) {
        iregion = domain->find_region(idregion);
        if (iregion == -1)
            error->all(FLERR,"Region ID for dump custom does not exist");
    }

    /* Define the group of variables for the atom style here since it's a fixed set */
    io = ad->DeclareIO(ioName);
    if (!io.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default writer
        io.SetEngine("BPFile");
        int num_aggregators = multiproc;
        if (num_aggregators == 0)
            num_aggregators = 1;
        char nstreams[128];
        sprintf (nstreams, "%d", num_aggregators);
        io.SetParameters({{"substreams", nstreams}});
        if (me==0 && screen) fprintf(screen, "ADIOS method for %s is n-to-m (aggregation with %s writers)\n", filename, nstreams);
    }


    io.DefineVariable<uint64_t>("ntimestep");
    io.DefineVariable<uint64_t>("natoms");

    io.DefineVariable<int>("nprocs");
    io.DefineVariable<int>("ncolumns");

    io.DefineVariable<double>("boxxlo");
    io.DefineVariable<double>("boxxhi");
    io.DefineVariable<double>("boxylo");
    io.DefineVariable<double>("boxyhi");
    io.DefineVariable<double>("boxzlo");
    io.DefineVariable<double>("boxzhi");

    io.DefineVariable<double>("boxxy");
    io.DefineVariable<double>("boxxz");
    io.DefineVariable<double>("boxyz");

    io.DefineAttribute<int>("triclinic",  domain->triclinic);

    int *boundaryptr = reinterpret_cast<int*>(domain->boundary);
    io.DefineAttribute<int>("boundary", boundaryptr, 6);

    size_t nColumns = static_cast<size_t>(size_one);
    io.DefineAttribute<std::string>("columns",  columnNames.data(), nColumns);
    io.DefineAttribute<std::string>("columnstr", columns);
    io.DefineAttribute<std::string>("boundarystr",  boundstr);
    io.DefineAttribute<std::string>("LAMMPS/dump_style",  "atom");
    io.DefineAttribute<std::string>("LAMMPS/version",  universe->version);
    io.DefineAttribute<std::string>("LAMMPS/num_ver",  universe->num_ver);

    io.DefineVariable<uint64_t>("nme", {adios2::LocalValueDim}); // local dimension variable
    io.DefineVariable<uint64_t>("offset", {adios2::LocalValueDim}); // local dimension variable

    // atom table size is not known at the moment
    // it will be correctly defined at the moment of write
    size_t UnknownSizeYet = 1;
    varAtoms = io.DefineVariable<double>("atoms", 
                                         {UnknownSizeYet,nColumns}, 
                                         {UnknownSizeYet, 0}, 
                                         {UnknownSizeYet,nColumns});

}

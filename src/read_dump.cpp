/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Timothy Sirk
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "read_dump.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "irregular.h"

using namespace LAMMPS_NS;

#define MAXATOM 2000   //max atoms for one comm
#define MAXLINE 1024    //max char in a line

#define XYZ  1
#define XYZV 2
#define XYZI 4

/* ---------------------------------------------------------------------- */

ReadDump::ReadDump(LAMMPS *lmp) : Pointers(lmp)
{
    line = new char[MAXLINE];
    keyword = new char[MAXLINE];
    narg = maxarg = 0;
    arg = NULL;
    boxsize= new double[9]; // orth or tri
}

/* ---------------------------------------------------------------------- */

ReadDump::~ReadDump()
{
    delete [] line;
    delete [] keyword;
    delete [] boxsize;
    if(me==0) memory->destroy(buf);
}

/* ---------------------------------------------------------------------- */

void ReadDump::command(int narg, char **arg)
{
    int ttime, assigned=0;
    if (narg != 2) error->all(FLERR,"Illegal read_dump command");

    // set some basic output info
    // important for use with rerun
    setup(false);

    // clear old per-atom properties
    if(clear) clearAtom();

    // read a frame for given time
    // pack whole frame into buffer
    // optionally close file
    if (me == 0)
    {
        open(arg[0]);
        // set target step
        ttime=atoi(arg[1]);
        // read until requested step
        findFrame(ttime);
        // check header and load the frame
        getHeader();
        packFrame(true);
    }

    // communication to size buffer, etc
    commBuffInfo();

    // broadcast dump in chunks
    // ratoms = remaining atoms to send
    // nchunk = atoms to send this time

    while(ratoms>0)
    {
        nchunk=MIN(ratoms,MAXATOM);

        // nchunk to all procs
        MPI_Bcast(&nchunk,1,MPI_INT,0,world);
        // let procs find owned atoms
        // pass in how many atoms already sent
        sendCoord(ndatoms-ratoms);
        // copy in new data to owned atoms
        // keep track of how many atoms assigned
        assigned += updateCoord();
        ratoms=ratoms-MAXATOM;
    }

    // move atoms back inside simulation box and to new processors
    migrateAtoms();


    // summarize
    int ntotal=0;
    char str[128];
    MPI_Allreduce(&assigned,&ntotal,1,MPI_INT,MPI_SUM,world);
    sprintf(str,"Assigned %d of %d atoms from dump", ntotal,ndatoms);
    if(me==0) error->message(FLERR,str);

    // update the image flags
    // don't bother to chunk
    if ((exist & XYZI) == XYZI){
       updateImages();
       memory->destroy(imagebuf);    
    }

}

/* ----------------------------------------------------------------------
settings for read_dump
------------------------------------------------------------------------- */

void ReadDump::setup(bool rerun)
{
    MPI_Comm_rank(world,&me);

    // binary existance flags, add more as needed
    // unset later if not present
    exist = XYZ | XYZV | XYZI;

    //set flags for default and rerun
    if(rerun)
    {
        quiet=true; // not verbose
        fileclose=false; // leave dump file open
        clear=true; // clear per-atom properties before update
    }
    else
    {
        quiet=false;
        fileclose=true;
        clear=false;
    }
}

/* ----------------------------------------------------------------------
   clear arrays
------------------------------------------------------------------------- */

void ReadDump::clearAtom()
{
    // per atom info to clear

    double **v = atom->v;
    for (int i = 0; i < atom->nlocal+atom->nghost; i++)
    {
        v[i][0] = 0.0;
        v[i][1] = 0.0;
        v[i][2] = 0.0;
    }
}


/* ----------------------------------------------------------------------
   proc 0 opens dump file
   test if gzipped
------------------------------------------------------------------------- */

void ReadDump::open(char *file)
{
    compressed = 0;
    char *suffix = file + strlen(file) - 3;
    if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
    if (!compressed) fp = fopen(file,"r");
    else
    {
        error->one(FLERR,"Cannot open gzipped file");
    }

    if (fp == NULL)
    {
        char str[128];
        sprintf(str,"Cannot open file %s",file);
        error->one(FLERR,str);
    }
}

/* ----------------------------------------------------------------------
   find right frame
------------------------------------------------------------------------- */

void ReadDump::findFrame(int ttime)
{

    // loop until section found with matching keyword
    char* keyword="TIMESTEP";
    long int offset=0;

    while (1)
    {

        if (fgets(line,MAXLINE,fp) == NULL) error->one(FLERR,"Unexpected end of file");
        if (!strstr(line,keyword)) error->one(FLERR,"Bad dump format");

        // found TIMESTEP
        time=atoi(fgets(line,MAXLINE,fp));
        line=fgets(line,MAXLINE,fp);

        // position for NUMBER OF ATOMS
        offset=ftell(fp);

        line=fgets(line,MAXLINE,fp);
        ndatoms=atoi(line);

        if(time==ttime) break;

        // no match, drop remaining header and atoms
        skip_lines(5);
        skip_lines(ndatoms);
    }

    // found frame
    // reset file position to top of frame
    fseek(fp, offset, SEEK_SET);
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
------------------------------------------------------------------------- */

void ReadDump::skip_lines(int n)
{
    char *eof;
    for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of dump file");
}

/* ----------------------------------------------------------------------
   read header
------------------------------------------------------------------------- */

void ReadDump::getHeader()
{
    bstride=1; // buffer stride must be at least 1
    int count=0,xs=1;

    // get num dump atoms
    line=fgets(line,MAXLINE,fp);
    ndatoms=atoi(line);

    if (ndatoms != atom->natoms && !quiet)
        error->warning(FLERR,"Different Natoms in dump snapshot");

    fgets(line,MAXLINE,fp);
    if (!strstr(line,"ITEM: BOX BOUNDS")) error->one(FLERR,"No box bounds"); // format check

    //box sizes
    if (domain->triclinic == 0)
    {
        fgets(line,MAXLINE,fp);
        sscanf (line,"%lf %lf",&boxsize[0],&boxsize[1]); // xlo xhi
        fgets(line,MAXLINE,fp);
        sscanf (line,"%lf %lf",&boxsize[3],&boxsize[4]);
        fgets(line,MAXLINE,fp);
        sscanf (line,"%lf %lf",&boxsize[6],&boxsize[7]);
    }
    else
    {
        fgets(line,MAXLINE,fp);
        sscanf (line,"%lf %lf %lf",&boxsize[0],&boxsize[1],&boxsize[2]);     // xlo xhi xy
        fgets(line,MAXLINE,fp);
        sscanf (line,"%lf %lf %lf",&boxsize[3],&boxsize[4],&boxsize[5]);
        fgets(line,MAXLINE,fp);
        sscanf (line,"%lf %lf %lf",&boxsize[6],&boxsize[7],&boxsize[8]);
    }

    // tokenize custom headers, can be any order
    // record position of x,v,etc in dump header

    fgets(line,MAXLINE,fp);
    if (!strstr(line,"ATOMS")) error->one(FLERR,"Misplaced or Missing ATOMS section");

    // position header, offset by two for zero index
    int nheader=-2;
    for (int i=0; i<10; i++) header[i]=-1;
    char *pch = strtok(line," \n");

    while (pch != NULL)
    {

        // fields to read
        if(!strcmp(pch,"id")) header[0]=nheader;
        if(!strcmp(pch,"x"))  header[1]=nheader;
        if(!strcmp(pch,"y"))  header[2]=nheader;
        if(!strcmp(pch,"z"))  header[3]=nheader;
        if(!strcmp(pch,"xu")) header[1]=nheader;
        if(!strcmp(pch,"yu")) header[2]=nheader;
        if(!strcmp(pch,"zu")) header[3]=nheader;
        if(!strcmp(pch,"vx")) header[4]=nheader;
        if(!strcmp(pch,"vy")) header[5]=nheader;
        if(!strcmp(pch,"vz")) header[6]=nheader;
        if(!strcmp(pch,"ix")) header[7]=nheader;
        if(!strcmp(pch,"iy")) header[8]=nheader;
        if(!strcmp(pch,"iz")) header[9]=nheader;

        // scaled coordinates
        if(!strcmp(pch,"xs") || !strcmp(pch,"xsu"))
        {
            header[1]=nheader;
            xs=xs<<1;
        }
        if(!strcmp(pch,"ys") || !strcmp(pch,"ysu"))
        {
            header[2]=nheader;
            xs=xs<<1;
        }
        if(!strcmp(pch,"zs") || !strcmp(pch,"zsu"))
        {
            header[3]=nheader;
            xs=xs<<1;
        }

        pch = strtok(NULL, " \n");
        nheader++;
    }

    // error if no id, warn for others

    // check for id (required)
    if(header[0]<0)
    {
        error->one(FLERR,"Missing ID in dump header");
    }

    // check for coordinates (optional)
    for(int i=1; i<4; i++)
        if(header[i]<0)
        {
            if(!quiet) error->warning(FLERR,"Missing coordinates in dump header");
            exist=exist^XYZ;
            break;
        }

    // check for velocities (optional)
    for(int i=4; i<7; i++)
        if(header[i]<0)
        {
            if(!quiet) error->warning(FLERR,"No velocities in dump header");
            exist=exist^XYZV;
            break;
        }

    // check for image flags (optional)
    for(int i=7; i<10; i++)
        if(header[i]<0)
        {
            if(!quiet) error->warning(FLERR,"No image flags in dump header");
            exist=exist^XYZI;
            break;
        }

    // map headers position to buffer position
    // must have at least the ID

    hndx[0]=count++;

    if ((exist & XYZ) == XYZ)
    {
        if(!quiet) error->message(FLERR,"Reading coordinates");
        bstride=bstride+3;
        hndx[1]=count++;
        hndx[2]=count++;
        hndx[3]=count++;

        // check scale coordinates
        // all xs,ys,zs must be present
        if (xs==8)
        {
            scale[0]=boxsize[1]-boxsize[0];
            scale[1]=boxsize[4]-boxsize[3];
            scale[2]=boxsize[7]-boxsize[6];
            if(!quiet) error->message(FLERR,"Found scaled coordinates");
        }
        else scale[0]=scale[1]=scale[2]=1;
    }

    if ((exist & XYZV) == XYZV)
    {
        if(!quiet) error->message(FLERR,"Reading velocities");
        bstride=bstride+3;
        hndx[4]=count++;
        hndx[5]=count++;
        hndx[6]=count++;
    }

    if ((exist & XYZI) == XYZI)
    {
        if(!quiet) error->message(FLERR,"Reading image flags");
        hndx[7]=count++;
        hndx[8]=count++;
        hndx[9]=count++;
    }
}

/* ----------------------------------------------------------------------
   read header
------------------------------------------------------------------------- */

void ReadDump::packFrame(bool fileclose)
{
    char *pch;
    double tmp;
    int nheader;

    //allocate
    memory->create(buf,ndatoms,bstride,"command:buf");

    // pack entire frame into a buffer

    if((exist & XYZI) == XYZI)
    {
        memory->create(imagebuf,4*ndatoms,"command:imagebuf");
    }


    for(int i=0; i<ndatoms; i++)
    {
        fgets(line,MAXLINE,fp);
        pch = strtok(line, " ");
        nheader=0;

        while (pch != NULL)
        {
            tmp=atof(pch);

            if(nheader==header[0])
            {
                buf[i][hndx[0]]=tmp; // id
                if((exist & XYZI) == XYZI) imagebuf[i*4]=(int)tmp; // id for images (optional)
            }

            if((exist & XYZ) == XYZ)
            {
                if(nheader==header[1])  buf[i][hndx[1]]=tmp*scale[0];
                if(nheader==header[2])  buf[i][hndx[2]]=tmp*scale[1];
                if(nheader==header[3])  buf[i][hndx[3]]=tmp*scale[2];
            }

            if((exist & XYZV) == XYZV)
            {
                if(nheader==header[4])  buf[i][hndx[4]]=tmp;
                if(nheader==header[5])  buf[i][hndx[5]]=tmp;
                if(nheader==header[6])  buf[i][hndx[6]]=tmp;
            }

            if((exist & XYZI) == XYZI)
            {

                // save image flags into seperate array
                // can't compress to binary yet b/c ix iy iz can be in any order

                if(nheader==header[7])  imagebuf[i*4+1]=(int)tmp;
                if(nheader==header[8])  imagebuf[i*4+2]=(int)tmp;
                if(nheader==header[9])  imagebuf[i*4+3]=(int)tmp;
            }

            nheader++;
            pch = strtok (NULL, " ");
        }
    }

    if(fileclose) fclose(fp);
}

/* ----------------------------------------------------------------------
Communicate some info about this frame
------------------------------------------------------------------------- */

void ReadDump::commBuffInfo()
{

    // total number of atoms
    MPI_Bcast(&ndatoms,1,MPI_INT,0,world);
    // remaining atoms to comm
    // begin with all atoms
    ratoms=ndatoms;

    // width of buffer
    MPI_Bcast(&bstride,1,MPI_INT,0,world);

    // update time
    MPI_Bcast(&time,1,MPI_INT,0,world);
    update->ntimestep=(bigint)time;

    // existance bitflags
    MPI_Bcast(&exist,1,MPI_INT,0,world);

    // box size
    MPI_Bcast(boxsize,9,MPI_DOUBLE,0,world);
}


/* ----------------------------------------------------------------------
broadcast all the dump atom info from a buffer. Could improve by
sending just sending one message to each proc that contains owned atoms
or just have each proc read the file?
------------------------------------------------------------------------- */

void ReadDump::sendCoord(int last)
{

    if (domain->triclinic == 0)
    {
        domain->boxlo[0]=boxsize[0];
        domain->boxhi[0]=boxsize[1];
        domain->boxlo[1]=boxsize[3];
        domain->boxhi[1]=boxsize[4];
        domain->boxlo[2]=boxsize[6];
        domain->boxhi[2]=boxsize[7];
    }
    else
    {
        domain->boxlo[0]=boxsize[0];
        domain->boxhi[0]=boxsize[1];
        domain->xy=boxsize[2];
        domain->boxlo[1]=boxsize[3];
        domain->boxhi[1]=boxsize[4];
        domain->xy=boxsize[5];
        domain->boxlo[2]=boxsize[6];
        domain->boxhi[2]=boxsize[7];
        domain->xy=boxsize[8];
    }

    int bsize=nchunk*bstride;
    xbuffer = new double[bsize];

    //fill buffer on proc 0

    if(me==0)
    {
        int z=0;
        for(int y=0; y<nchunk; y++)
        {
            for(int j=0; j<bstride; j++)
                xbuffer[z+j]=buf[y+last][j];
            z=z+bstride;
        }
    }

    MPI_Bcast(xbuffer,bsize,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
update atoms using dump info
------------------------------------------------------------------------- */

int ReadDump::updateCoord()
{
    double** x=atom->x;
    double** v=atom->v;
    int nlocal=atom->nlocal;
    int assigned=0;
    int j=0;
    int id=0,gid=0;

    // customize with more per-atom data
    // nchunk atoms per comm

    // need a global->local map
    // make an array if no mapping exists

    if(atom->map_style==0)
    {
        atom->map_style=1;
        atom->map_init();
        atom->map_set();
    }

    for(int n=0; n<nchunk; n++)
    {
        j=n*bstride;
        gid=(int)xbuffer[j++];
        id=atom->map( gid );

        if(id >= 0 && id < nlocal)
        {
            // add more fields as needed
            if((exist & XYZ) == XYZ)
            {
                x[id][0]=xbuffer[j++];
                x[id][1]=xbuffer[j++];
                x[id][2]=xbuffer[j++];
            }

            if((exist & XYZV) == XYZV)
            {
                v[id][0]=xbuffer[j++];
                v[id][1]=xbuffer[j++];
                v[id][2]=xbuffer[j++];
            }

            assigned++;
        }
    }

    delete [] xbuffer;

    // return the number of atoms updated
    return assigned;
}

/* ----------------------------------------------------------------------
  move atoms back inside simulation box and to new processors
  use remap() instead of pbc() in case atoms moved a long distance
  use irregular() in case atoms moved a long distance
------------------------------------------------------------------------- */

void ReadDump::migrateAtoms()
{
    double **x = atom->x;
    int nlocal = atom->nlocal;

    // convert coord if tric box
    if (domain->triclinic) domain->x2lamda(atom->nlocal);

    // reset box and remap before communicating
    // remap will modify image flags
    domain->reset_box();
    for (int i = 0; i < nlocal; i++) domain->remap(x[i],atom->image[i]);

    // irregular comm
    Irregular *irregular = new Irregular(lmp);
    irregular->migrate_atoms();
    delete irregular;

    // convert coord if tric box
    if (domain->triclinic) domain->lamda2x(atom->nlocal);
}


/* ----------------------------------------------------------------------
update the atoms image flags. must be done after remap.
------------------------------------------------------------------------- */

void ReadDump::updateImages()
{
    int gid, lid; // global and local id

    // allocate for images, except on proc 0
    if(me != 0)  memory->create(imagebuf,4*ndatoms,"command:imagebuf");

    MPI_Bcast(imagebuf,ndatoms*4,MPI_INT,0,world);

    for(int i=0; i<ndatoms; i++)
    {

        gid=(int)imagebuf[i*4];
        lid=atom->map(gid);

        if(lid >= 0 && lid < atom->nlocal)
        {
            atom->image[lid] = ((imagebuf[i*4+3] + 512 & 1023) << 20) | ((imagebuf[i*4+2] + 512 & 1023)  << 10) | (imagebuf[i*4+1]+ 512 & 1023);
        }
    }

}

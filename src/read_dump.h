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

#ifdef COMMAND_CLASS

CommandStyle(read_dump,ReadDump)

#else

#ifndef LMP_READ_DUMP_H
#define LMP_READ_DUMP_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS
{

class ReadDump : protected Pointers
{

public:

    ReadDump(class LAMMPS *);
    ~ReadDump();
    void command(int, char **);
    void setup(bool); // bool for rerun
    void clearAtom();
    void findFrame(int); // which frame
    void getHeader();
    void packFrame(bool);  // bool to close file
    void commBuffInfo();
    void sendCoord(int);
    int updateCoord();
    void migrateAtoms();
    void updateImages();

    // file methods
    void open(char *);
    void skip_lines(int);

    int me,ndatoms,ratoms,nchunk,last,exist,bstride,time;
    bool fileclose,quiet,clear;

private:

    char *line,*keyword;
    FILE *fp;
    int narg,maxarg,compressed;
    char **arg;
    double **buf;
    double *xbuffer;
    double *boxsize;
    int hndx[10];
    int header[10];
    double scale[3];
    int* imagebuf;


};

}

#endif
#endif


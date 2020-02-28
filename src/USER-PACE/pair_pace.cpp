//
// Created by lysogy36 on 27.02.20.
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_pace.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"


#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4

//added YL
int elements_num = 104;
char const *const elements[104] = {"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                                   "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
                                   "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
                                   "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                                   "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                                   "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                                   "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                                   "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
};

int AtomicNumberByName(char *elname) {
    int i = 0;
    for (i = 1; i < elements_num; i++)
        if (strcmp(elname, elements[i]) == 0)
            return i;
    return -1;
}

/* ---------------------------------------------------------------------- */
PairPACE::PairPACE(LAMMPS *lmp) : Pair(lmp) {
    //single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    nelements = 0;

    aceptr = NULL;
    potential_file_name = NULL;
    elements = NULL;
    map = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPACE::~PairPACE() {
    if (copymode) return;

    if (elements)
        for (int i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;


    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        delete[] map;
    }

    delete[] potential_file_name;
    delete aceptr;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(map);
    }
}

/* ---------------------------------------------------------------------- */

void PairPACE::compute(int eflag, int vflag) {
    //  printf("--> PairPACE::compute\n");
    int i, j, k, ii, jj, kk, inum, jnum;
    int itype, jtype, ktype, iparam_ij, iparam_ijk;
    tagint itag, jtag;
    double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
    double rsq, rsq1, rsq2;
    double delr1[3], delr2[3], fi[3], fj[3], fk[3];
    double zeta_ij, prefactor;
    int *ilist, *jlist, *numneigh, **firstneigh;

    evdwl = 0.0;
    ev_init(eflag, vflag);

    // downwards modified by RD

    double **x = atom->x;
    double **f = atom->f;
    tagint *tag = atom->tag;
    int *type = atom->type;
    // number of atoms in cell
    int nlocal = atom->nlocal;
    //printf("nlocal=%d\n",nlocal);
    int newton_pair = force->newton_pair;
    //printf("newton_pair=%d\n",newton_pair);
    const double cutshortsq = cutmax * cutmax;
    // number of atoms including ghost atoms
    int nall = nlocal + atom->nghost;
    //printf("nall=%d\n",nall);

    // inum: length of the neighborlists list
    inum = list->inum;
    //printf("inum=%d\n",inum);
    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;
    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;
    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;

    // RD new stuff
    int n;
    int nele;
    int neicount;
    int nei;
    int natom;
    int *nlist;
    int *nstart;
    int *nstop;
    int *occlist;
    int *occ;
    double *dlist;
    double *elist;
    double *eout;
    double *flist;

    //parameters read from pair_style:
    // pair_style ace nelements [1] rcut[8.7] xcfunctional[pbe/lda] smear[0.1] nele[2] elelist[Al Co]
    // pair_style ace 8.7 pbe 0.1 1 Al

    //double rcut = 8.7;
    //double smear = 0.1;
    //char elelist[2];
    //char xcfunctional[]="pbe";

    double fxtmp, fytmp, fztmp;

    // Aidan Thompson told RD (26 July 2019) that practically always holds:
    // inum = nlocal
    // i = ilist(ii) < inum
    // j = jlist(jj) < nall
    // neighborlist contains meighbor atoms plus skin atoms,
    //       skin atoms can be removed by setting skin to zero but here
    //       they are disregarded anyway

    // TODO: need to get elements properly from LAMMPS
    // bad: need proper string assignment
    // - make memory allocation more efficient and revoid regular malloc/free operations
    // - allocate, populate and distribute three vectors gvec1, gvec2, gvec3


    //nele = 1;
    //  strncpy(elelist, "Al", 2);
    //  strncpy(xcfunctional, "pbe", 3);

    if (inum != nlocal) {
        printf("inum: %d nlocal: %d are different.\n", inum, nlocal);
        exit(0);
    }

    // find number of neighbors and total number of atoms, including ghost atoms
    natom = nlocal;
    // printf("natom=%d\n",natom);
    nei = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        //        printf("ilist[%d]=%d\n",ii,i);
        jnum = numneigh[i];
        //        printf("numneigh[%d]=%d\n",i,jnum);
        nei = nei + jnum;
        //        printf("firstneigh[%d]/jlist=\n",i);
        //        jlist = firstneigh[i];
        // loop all neighbors of atom i
        //        for (jj = 0; jj < jnum; jj++) {
        //            if (jlist[jj]<natom)
        //                printf("\n%d\n", jlist[jj]);
        //            else
        //                printf("%d ", jlist[jj]);
        //        }
        //        printf("\n");
    }

    // allocate space for my neighborlist
    // (TODO: need other variable to tell if array sizes have changed since last call to avoid frequent allocation)
    nstart = (int *) malloc(nlocal * sizeof(int *));
    nstop = (int *) malloc(nlocal * sizeof(int *));
    occ = (int *) malloc(nlocal * sizeof(int *));
    nlist = (int *) malloc(nei * sizeof(int *));
    occlist = (int *) malloc(nei * sizeof(int *));
    dlist = (double *) malloc(nei * sizeof(double *));
    elist = (double *) malloc(3 * nei * sizeof(double *));
    flist = (double *) malloc(3 * nei * sizeof(double *));
    eout = (double *) malloc(nlocal * sizeof(double *));

    // loop all atoms and neighbors
    neicount = 0;
    for (ii = 0; ii < inum; ii++) { //4 times
        i = ilist[ii];
        k = map[type[i]];
        occ[ii] = k;
        // this to be removed once fully integrated in LAMMPS and Fortran call fixed
        if (occ[ii] != 0) {
            printf("Error: Cannot handle multiple species.\n");
            printf("Stopping.\n");
            exit(0);
        }
        nstart[ii] = neicount;
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        //	printf("Position of at.%d=(%f,%f,%f)\n",i,x[i][0],x[i][1],x[i][2]);
        jlist = firstneigh[i];
        jnum = numneigh[i];

        // loop all neighbors of atom i
        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;

            if (j >= nall) {
                printf("Error: Atoms index too large.\n");
                printf("Stopping.\n");
                exit(0);
            }

            nlist[neicount] = j;
            k = map[type[j]];
            occlist[neicount] = k;

            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx * delx + dely * dely + delz * delz;
            rsq = sqrt(rsq);
            // store distances
            dlist[neicount] = rsq;
            if (rsq < 1.e-6) {
                printf("Error: Atoms with identical positions.\n");
                printf("Stopping.\n");
                exit(0);
            }
            // store direction vectors
            elist[3 * neicount + 0] = delx / rsq;
            elist[3 * neicount + 1] = dely / rsq;
            elist[3 * neicount + 2] = delz / rsq;
            neicount++;
        }
        nstop[ii] = neicount - 1;
    }

    if (nei != neicount) {
        printf("Something wrong with neighbor list translation: neicount\n");
        printf("Stopping.\n");
        exit(0);
    }
//    printf("neicount=%d\n",neicount);

    // up by +1 to run in Fortran style
    for (ii = 0; ii < neicount; ii++) {
        nlist[ii] = nlist[ii] + 1;
        occlist[ii] = occlist[ii] + 1;
    }
    for (ii = 0; ii < natom; ii++) {
        nstart[ii] = nstart[ii] + 1;
        nstop[ii] = nstop[ii] + 1;
        occ[ii] = occ[ii] + 1;
    }






    // check if this is necessary
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        f[i][0] = 0.0;
        f[i][1] = 0.0;
        f[i][2] = 0.0;
        for (jj = 0; jj < jnum; jj++) {
            jlist = firstneigh[i];
            j = jlist[jj];
            j &= NEIGHMASK;
            f[j][0] = 0.0;
            f[j][1] = 0.0;
            f[j][2] = 0.0;
        }
    }


    //TODO: check the force summation
    k = 0;
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        // xtmp = x[i][0];
        // ytmp = x[i][1];
        //ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        //fxtmp = fytmp = fztmp = 0.0;

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
            //  delx = xtmp - x[j][0];
            //dely = ytmp - x[j][1];
            //delz = ztmp - x[j][2];
            //rsq = delx*delx + dely*dely + delz*delz;
            fxtmp = flist[k + 0];
            fytmp = flist[k + 1];
            fztmp = flist[k + 2];

            //if(fxtmp!=0.0 || fytmp!=0.0 || fztmp!=0.0){
            //                printf("Pair force (%d-%d) = (%f, %f, %f)\n",i,j,fxtmp,fytmp,fztmp);
            //}
            f[i][0] -= fxtmp;
            f[i][1] -= fytmp;
            f[i][2] -= fztmp;
            //            if (newton_pair) {
            f[j][0] += fxtmp;
            f[j][1] += fytmp;
            f[j][2] += fztmp;
            //}
            k = k + 3;

            if (vflag)
                ev_tally_xyz(i, j, nlocal, newton_pair,
                             0.0, 0.0, fxtmp, fytmp, fztmp, delx, dely, delz);
        }
        if (eflag) {
            ev_tally_full(i, 2.0 * eout[i], 0.0, 0.0, 0.0, 0.0, 0.0);
        }

        // note: incorrect force print out as still inside of loop
        //        printf("Total force II on at.%d=(%f,%f,%f)\n",i,f[i][0],f[i][1],f[i][2]);
    }

    //printf("LAMMPS force=\n");
    //for (ii = 0; ii < inum; ii++) {
    //    i = ilist[ii];
    //    printf("(%f, %f, %f)\n",f[i][0], f[i][1], f[i][2]);
    //
    //}


    free(nlist);
    free(occlist);
    free(dlist);
    free(elist);
    free(nstart);;
    free(nstop);
    free(occ);
    free(flist);
    free(eout);

    printf("vflag_fdotr=%d\n", vflag_fdotr);
    if (vflag_fdotr) virial_fdotr_compute();


    // end modifications RD
    //-------------------------------------------------------


}

/* ---------------------------------------------------------------------- */

void PairPACE::allocate() {
    //   printf("--> PairPACE::allocate\n");
    allocated = 1;
    int n = atom->ntypes;

    memory->create(setflag, n + 1, n + 1, "pair:setflag");
    memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
    memory->create(map, n + 1, "pair:map");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPACE::settings(int narg, char **arg) {
    //   printf("--> PairPACE::settings\n");
    //pair_style 0:rcut 1:xcfunctional 2:smear 3..:elelist[1..nele]
    if (narg > 0)
        error->all(FLERR,
                   "Illegal pair_style command. Correct form:\npair_style pace");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACE::coeff(int narg, char **arg) {
    for (int i = 0; i < narg; i++) {
        printf("PairPACE::coeff: arg[%d] = %s\n", i, arg[i]);
    }

    if (narg < 4)
        error->all(FLERR,
                   "Incorrect args for pair coefficients. Correct form:\npair_coeff * * <potential.ace> elem1 elem2 ...");

    if (!allocated) allocate();

    //number of provided elements in pair_coeff line
    int ntypes_coeff = narg - 3;

    if (ntypes_coeff != atom->ntypes) {
        char error_message[1024];
        sprintf(error_message,
                "Incorrect args for pair coefficients. You provided %d elements in pair_coeff, but structure has %d atom types",
                ntypes_coeff, atom->ntypes);
        error->all(FLERR, error_message);
    }

    char *type1 = arg[0];
    char *type2 = arg[1];
    char *potential_file_name = arg[2];

    // insure I,J args are * *

    if (strcmp(type1, "*") != 0 || strcmp(type2, "*") != 0)
        error->all(FLERR, "Incorrect args for pair coefficients");

    elelist = new char *[ntypes_coeff];
    elenumbers = new int[ntypes_coeff];
    for (int i = 0; i < ntypes_coeff; i++) {
        n = strlen(arg[shift + i]) + 1;
        elelist[i] = new char[n];
        strcpy(elelist[i], arg[i + shift]);

        int atomicnumber = AtomicNumberByName(arg[i + shift]);
        if (atomicnumber == -1) {
            sprintf(buf, "%s is not a valid element name", arg[i + shift]);
            error->all(FLERR, buf);
        }
        elenumbers[i] = atomicnumber;
        printf("element: %s -> %d\n", arg[i + shift], elenumbers[i]);
    }


    // insure I,J args are * *
    if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
        error->all(FLERR, "Incorrect args for pair coefficients");

    // read args that map atom types to elements in potential file
    // map[i] = which element the Ith atom type is, -1 if NULL
    // nelements = # of unique elements
    // elements = list of element names

    if (elelist) {
        for (i = 0; i < nelements; i++) delete[] elelist[i];
        delete[] elelist;
    }
    elelist = new char *[atom->ntypes];
    for (i = 0; i < atom->ntypes; i++) elelist[i] = NULL;

    nelements = 0;
    for (i = shift; i < narg; i++) {
        if (strcmp(arg[i], "NULL") == 0) {
            map[i - shift + 1] = -1;
            continue;
        }
        for (j = 0; j < nelements; j++)
            if (strcmp(arg[i], elelist[j]) == 0) break;
        map[i - shift + 1] = j;
        if (j == nelements) {
            n = strlen(arg[i]) + 1;
            elelist[j] = new char[n];
            strcpy(elelist[j], arg[i]);
            nelements++;
        }
    }

    for (int i = 0; i < nele; i++) {
        printf("elelist[%d]=%s\n", i, elelist[i]);
    }

    // clear setflag since coeff() called once with I,J = * *

    n = atom->ntypes;
    for (i = 1; i <= n; i++)
        for (j = i; j <= n; j++)
            setflag[i][j] = 0;

    // set setflag i,j for type pairs where both are mapped to elements

    int count = 0;
    for (i = 1; i <= n; i++)
        for (j = i; j <= n; j++)
            if (map[i] >= 0 && map[j] >= 0) {
                setflag[i][j] = 1;
                count++;
            }

    if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

    printf("nelements = %d\n", nelements);


    printf("xcfunctional=%s\n", xcfunctional);
    nxc = -1;
    if (strcmp(xcfunctional, "new") == 0) {
        nxc = 0;//generate reference input file only
    }
    if (strcmp(xcfunctional, "pbe") == 0) {
        nxc = 1;//PBE
    }
    if (strcmp(xcfunctional, "lda") == 0) {
        nxc = 2; //LDA
    }
    if (nxc == -1) error->all(FLERR, "Exchange functional not recognized");
    else printf("Reference exchange functional: %s (=%d)\n", xcfunctional, nxc);

    // up 1 for index in Fortran style
    for (i = 0; i < nele; i++) {
        elenumbers[i] = elenumbers[i] + 1;
    }

    // here everything is in place for ace call I: ini
    //acelibini_(&nele, &rcut, &smear, &nxc, elenumbers);
    // printf("acelibini_done\n");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPACE::init_style() {
    if (atom->tag_enable == 0)
        error->all(FLERR, "Pair style ACE requires atom IDs");
    if (force->newton_pair == 0)
        error->all(FLERR, "Pair style ACE requires newton pair on");

    // need a full neighbor list

    int irequest = neighbor->request(this, instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPACE::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

    return cutmax;
}

/* ---------------------------------------------------------------------- */


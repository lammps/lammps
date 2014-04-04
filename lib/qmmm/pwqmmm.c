/* Copyright (C) 2013 Quantum ESPRESSO group
 * This file is distributed under the terms of the
 * GNU General Public License. See the file `License'
 * in the root directory of the present distribution,
 * or http://www.gnu.org/copyleft/gpl.txt . */

/* Toplevel wrapper code to launch QM/MM calculations via MPI */

#include <mpi.h>

#include <unistd.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "libqecouple.h"
#include "libqmmm.h"

#include "library.h"

static const char delim[] = " \t\n\r";

int main(int argc, char **argv)
{
    MPI_Comm intra_comm, qm_comm, mm_comm;
    int me, ncpu, nqm, key, retval;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&ncpu);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);

    /* we accept just the qm/mm file as the one argument */
    if ((argc < 2) || (argc > 3)) {
        if (me == 0) fprintf(stderr,"\n usage: %s <qmmm input file> "
                             "[<number of mm procs>]\n\n",argv[0]);
        MPI_Finalize();
        return -1;
    }

    /* parse the qm/mm file */
    if (read_qmmm_config(argv[1],&qmmmcfg)) {
        if (me == 0) fputs("\n Error parsing QM/MM config file\n\n",stderr);
        MPI_Finalize();
        return -2;
    }

    key = 2;
    if (argc == 3) {
        key = atoi(argv[2]);
    }
    qmmmcfg.nmm = key;

    /* sanity checks */
    qmmmcfg.comm_mode = QMMM_COMM_MPI;
    nqm = ncpu - qmmmcfg.nmm;
    retval = 0;

    if (me == 0) {
        const char *msg;

        msg = check_qmmm_config(&qmmmcfg);
        
        if ((nqm < 1) || (qmmmcfg.nmm < 2)) {
            msg = "Need at least 2 MM and 1 QM processes";
        }
        if (msg != NULL) {
            retval = 1;
            fprintf(stderr,"\n%s\n\n",msg);
        }
    }
    MPI_Bcast(&retval,1,MPI_INT,0,MPI_COMM_WORLD);
    if (retval != 0) {
        MPI_Finalize();
        return -3;
    }

    if (me == 0) {
        write_qmmm_config(NULL,&qmmmcfg);

        printf("\nRunning QM/MM calculation for %d steps on %d procs\n\n",
               qmmmcfg.steps, ncpu);
        printf("QM system: procs: % 4d | directory: %-16s | input: %-16s|\n",
               nqm,qmmmcfg.qmdir,qmmmcfg.qminp);
        printf("MM master: procs: % 4d | directory: %-16s | input: %-16s|\n",
               qmmmcfg.nmm-1,qmmmcfg.madir,qmmmcfg.mainp);
        printf("MM slave:  procs: % 4d | directory: %-16s | input: %-16s|\n\n",
               1,qmmmcfg.sldir,qmmmcfg.slinp);
    }

    /* set up MPI communication */

    /* process partitioning: QM | MM master | MM slave */
    qmmmcfg.role = QMMM_ROLE_QM;
    if (me >= nqm) qmmmcfg.role = QMMM_ROLE_MASTER;
    if (me == (ncpu-1)) qmmmcfg.role = QMMM_ROLE_SLAVE;

    MPI_Comm_split(MPI_COMM_WORLD, qmmmcfg.role, me, &intra_comm);
    qmmmcfg.my_comm = MPI_Comm_c2f(intra_comm);

    /* qm to mm-master inter communicator */
    key = MPI_UNDEFINED;
    if ((me == 0) || (me == nqm)) key = 1;

    MPI_Comm_split(MPI_COMM_WORLD, key, ncpu-me, &qm_comm);
    qmmmcfg.qm_comm = MPI_Comm_c2f(qm_comm);

    /* mm-slave to mm-master inter communicator */
    key = MPI_UNDEFINED;
    if ((me == (ncpu-1)) || (me == nqm)) key = 1;

    MPI_Comm_split(MPI_COMM_WORLD, key, me, &mm_comm);
    qmmmcfg.mm_comm = MPI_Comm_c2f(mm_comm);

    if (qmmmcfg.role == QMMM_ROLE_QM) {
        FILE *fp;
        int nimage,npots,npool,ntg,nband,ndiag;

        MPI_Comm_rank(intra_comm,&me);

        if (chdir(qmmmcfg.qmdir)) {
            if (me == 0) fputs("failure to change into QM directory\n",stderr);
            MPI_Abort(MPI_COMM_WORLD, 10);
        }

        /* check if input file is available and readable */
        fp = fopen(qmmmcfg.qminp,"r");
        if (fp == NULL) {
            if (me == 0) fprintf(stderr,"failure to open QM input for "
                                 "reading: %s\n", strerror(errno));
            MPI_Abort(MPI_COMM_WORLD, 20);
        } else fclose(fp);

        /* redirect output to file, if requested */
        if (strcmp(qmmmcfg.qmout,"NULL") != 0) {
            fp = freopen(qmmmcfg.qmout,"w",stdout);
            if (fp == NULL) {
                if (me == 0) fprintf(stderr,"failure to open QM output for "
                                     "writing: %s\n", strerror(errno));
                MPI_Abort(MPI_COMM_WORLD, 30);
            }
        }

        /* parse additional command line flags for pw.x */
        nimage=npots=npool=ntg=nband=ndiag=1;

        if (qmmmcfg.qmarg != NULL) {
            char *ptr = strtok(qmmmcfg.qmarg,delim);
            do {
                /* -nimage is not supported */
                if (strncmp("-npot",ptr,5) == 0) {
                    ptr=strtok(NULL,delim);
                    npots=atoi(ptr);
                } else if ((strncmp("-nk",ptr,3) == 0)
                           || (strncmp("-npoo",ptr,5) == 0)) {
                    ptr=strtok(NULL,delim);
                    npool=atoi(ptr);
                } else if (strncmp("-nt",ptr,3) == 0) {
                    ptr=strtok(NULL,delim);
                    ntg=atoi(ptr);
                } else if (strncmp("-nb",ptr,3) == 0) {
                    ptr=strtok(NULL,delim);
                    nband=atoi(ptr);
                } else if ((strncmp("-nd",ptr,3) == 0)
                           || (strncmp("-no",ptr,3) == 0)
                           || (strcmp("-nproc_diag",ptr) == 0)
                           || (strcmp("-nproc_ortho",ptr) == 0)) {
                    ptr=strtok(NULL,delim);
                    ndiag=atoi(ptr);
                } else {
                    if (me == 0) {
                        fprintf(stderr,"QM: ignoring flag '%s", ptr);
                        ptr=strtok(NULL,delim);
                        fprintf(stderr," %s'\n",ptr);
                    } else strtok(NULL,delim);
                }
            } while ((ptr=strtok(NULL,delim)));
        }

        retval = 0;
        if (me == 0) fprintf(stderr,"QM: nimage: %d  npots: %d  npools: %d  "
                             "ntg: %d  nband: %d  ndiag: %d\n",
                             nimage,npots,npool,ntg,nband,ndiag);

        /* setup and call Q-E. */
        c2qmmm_mpi_config(qmmmcfg.qmmm_mode, qmmmcfg.qm_comm,
                          qmmmcfg.verbose, qmmmcfg.steps);
        c2libpwscf(qmmmcfg.my_comm, nimage, npots, npool, ntg, nband, ndiag,
                   &retval, qmmmcfg.qminp);

        if (strcmp(qmmmcfg.qmout,"NULL") != 0)
            fclose(fp);

    } else if (qmmmcfg.role == QMMM_ROLE_MASTER) {
        FILE *fp;
        char *cuda, *echo, *suffix;
        void *lmp;

        MPI_Comm_rank(intra_comm,&me);

        if (chdir(qmmmcfg.madir)) {
            if (me == 0)
                fputs("failure to change into MM master directory\n",stderr);
            MPI_Abort(MPI_COMM_WORLD, 10);
        }

        /* check if input file is available and readable */
        fp = fopen(qmmmcfg.mainp,"r");
        if (fp == NULL) {
            if (me == 0) fprintf(stderr,"failure to open MM master input "
                                 "for reading: %s\n", strerror(errno));
            MPI_Abort(MPI_COMM_WORLD, 20);
        } else fclose(fp);

        /* redirect output to file, if requested */
        if (strcmp(qmmmcfg.maout,"NULL") != 0) {
            fp = freopen(qmmmcfg.maout,"w",stdout);
            if (fp == NULL) {
                if (me == 0) fprintf(stderr,"failure to open MM master output "
                                     "for writing: %s\n", strerror(errno));
                MPI_Abort(MPI_COMM_WORLD, 30);
            }
        }

        /* parse additional support command line flags for LAMMPS */  

        if (qmmmcfg.maarg != NULL) {
            char *ptr = strtok(qmmmcfg.maarg,delim);
            do {
                if ((strncmp("-c",ptr,2) == 0)
                    || (strncmp("-cuda",ptr,5) == 0)) {
                    ptr=strtok(NULL,delim);
                    cuda=strdup(ptr);
                } else if ((strncmp("-e",ptr,2) == 0)
                           || (strncmp("-echo",ptr,5) == 0)) {
                    ptr=strtok(NULL,delim);
                    echo=strdup(ptr);
                } else if ((strncmp("-sf",ptr,3) == 0)
                           || (strncmp("-suffix",ptr,7) == 0)) {
                    ptr=strtok(NULL,delim);
                    suffix=strdup(ptr);
                } else {
                    if (me == 0)
                        fprintf(stderr,"unsupported LAMMPS flag: %s\n",ptr);
                    MPI_Abort(MPI_COMM_WORLD, 40);
                }
            } while ((ptr=strtok(NULL,delim)));
        }
        char *lmpargs[5];
        lmpargs[0] = strdup("MM Master");
        lmpargs[1] = strdup("-log");
        lmpargs[2] = strdup("none");
        lmpargs[3] = strdup("-echo");
        lmpargs[4] = strdup("screen");

        lammps_open(5, lmpargs, intra_comm, &lmp);
        lammps_file(lmp,qmmmcfg.mainp);

        char runcmd[1024];
        sprintf(runcmd,"run %d\n",qmmmcfg.steps);
        lammps_command(lmp,runcmd);
        if (qmmmcfg.restart != NULL) {
          sprintf(runcmd,"write_restart %s\n",qmmmcfg.restart);
          lammps_command(lmp,runcmd);
        }
        lammps_close(lmp);

        fputs("MM: master done\n",stderr);
        if (strcmp(qmmmcfg.maout,"NULL") != 0)
            fclose(fp);

        retval = 0;

    } else if (qmmmcfg.role == QMMM_ROLE_SLAVE) {
        FILE *fp;
        char *cuda, *echo, *suffix;
        void *lmp;

        MPI_Comm_rank(intra_comm,&me);

        if (chdir(qmmmcfg.sldir)) {
            if (me == 0)
                fputs("failure to change into MM slave directory\n",stderr);
            MPI_Abort(MPI_COMM_WORLD, 10);
        }

        /* check if input file is available and readable */
        fp = fopen(qmmmcfg.slinp,"r");
        if (fp == NULL) {
            if (me == 0) fprintf(stderr,"failure to open MM slave input "
                                 "for reading: %s\n", strerror(errno));
            MPI_Abort(MPI_COMM_WORLD, 20);
        } else fclose(fp);

        /* redirect output to file, if requested */
        if (strcmp(qmmmcfg.slout,"NULL") != 0) {
            fp = freopen(qmmmcfg.slout,"w",stdout);
            if (fp == NULL) {
                if (me == 0) fprintf(stderr,"failure to open MM slave output "
                                     "for writing: %s\n", strerror(errno));
                MPI_Abort(MPI_COMM_WORLD, 30);
            }
        }

        /* parse additional support command line flags for LAMMPS */  

        if (qmmmcfg.slarg != NULL) {
            char *ptr = strtok(qmmmcfg.maarg,delim);
            do {
                if ((strncmp("-c",ptr,2) == 0)
                    || (strncmp("-cuda",ptr,5) == 0)) {
                    ptr=strtok(NULL,delim);
                    cuda=strdup(ptr);
                } else if ((strncmp("-e",ptr,2) == 0)
                           || (strncmp("-echo",ptr,5) == 0)) {
                    ptr=strtok(NULL,delim);
                    echo=strdup(ptr);
                } else if ((strncmp("-sf",ptr,3) == 0)
                           || (strncmp("-suffix",ptr,7) == 0)) {
                    ptr=strtok(NULL,delim);
                    suffix=strdup(ptr);
                } else {
                    if (me == 0)
                        fprintf(stderr,"unsupported LAMMPS flag: %s\n",ptr);
                    MPI_Abort(MPI_COMM_WORLD, 40);
                }
            } while ((ptr=strtok(NULL,delim)));
        }

        char *lmpargs[5];
        lmpargs[0] = strdup("MM slave");
        lmpargs[1] = strdup("-log");
        lmpargs[2] = strdup("none");
        lmpargs[3] = strdup("-echo");
        lmpargs[4] = strdup("screen");

        lammps_open(5, lmpargs, intra_comm, &lmp);
        lammps_file(lmp,qmmmcfg.slinp);

        char runcmd[64];
        sprintf(runcmd,"run %d\n",qmmmcfg.steps);
        lammps_command(lmp,runcmd);
        lammps_close(lmp);

        fputs("MM: slave done\n",stderr);
        if (strcmp(qmmmcfg.slout,"NULL") != 0)
            fclose(fp);
        retval = 0;

    } else {
        fputs("\n how on earth did you end up here?\n\n",stderr);
        retval = 0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);

    if (me == 0) fputs("\nQM/MM run complete.\n",stderr);
    MPI_Finalize();
    free_qmmm_config(&qmmmcfg);

    return retval;
}


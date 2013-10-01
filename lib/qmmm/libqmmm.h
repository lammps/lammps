/*
 * This file is distributed under the terms of the
 * GNU General Public License. See the file `License'
 * in the root directory of the present distribution,
 * or http://www.gnu.org/copyleft/gpl.txt .
 *
 * common definitions, APIs and global data for QM/MM interface
 */

#ifndef QE_LIBQMMM_H
#define QE_LIBQMMM_H

#ifdef __cplusplus
extern "C" {
#endif

#define QMMM_COMM_NONE  0
#define QMMM_COMM_MPI   1
#define QMMM_COMM_SHM   2

#define QMMM_MODE_OFF  -1
#define QMMM_MODE_NONE  0
#define QMMM_MODE_MECH  1
#define QMMM_MODE_ELEC  2

#define QMMM_ROLE_QM     1
#define QMMM_ROLE_MASTER 2
#define QMMM_ROLE_SLAVE  3

#define QMMM_OK     0
#define QMMM_ERROR -1

typedef struct {
    int comm_mode, qmmm_mode;    /* communication and coupling mode */
    char *qmdir, *madir, *sldir; /* directories to run codes in */
    char *qminp, *mainp, *slinp; /* input files for codes */
    char *qmout, *maout, *slout; /* stdout files for codes */
    char *qmcmd, *macmd, *slcmd; /* command to run codes (SHMEM only) */
    char *qmarg, *maarg, *slarg; /* extra flags to pass to code */
    int role;                    /* role of this rank */
    int steps;                   /* number of MD steps */
    int nmm;                     /* tasks reserved for MD (master and slave) */
    char *handle;                /* handle for SHEMEM communication */
    int my_comm, qm_comm, mm_comm; /* Fortran representation of MPI_Comms */
} qmmm_config_t;

extern qmmm_config_t qmmmcfg;    

int read_qmmm_config(const char *, qmmm_config_t *);
const char *check_qmmm_config(qmmm_config_t *);
void free_qmmm_config(qmmm_config_t *);

#ifdef __cplusplus
}
#endif

#endif /* QE_LIBQMMM_H */

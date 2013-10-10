/*
 * This file is distributed under the terms of the
 * GNU General Public License. See the file `License'
 * in the root directory of the present distribution,
 * or http://www.gnu.org/copyleft/gpl.txt .
 *
 * generic QM/MM support library
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "libqmmm.h"

#define BUF_SIZE 1024

/* global variable for global QM/MM configuration */
qmmm_config_t qmmmcfg;

/* local helper function: advance char pointer beyond leading whitespace */
static char *skip_whitespace(char *ptr)
{
    while ((*ptr == ' ') || (*ptr == '\t')
           || (*ptr == '\r') || (*ptr == '\n')) ++ptr;
    return ptr;
}

/* local helper function: trim string to remove trailing whitespace */
static void trim_whitespace(char *ptr)
{
    while ((*ptr != ' ') && (*ptr != '\t')
           && (*ptr != '\r') && (*ptr != '\n')) ++ptr;
    *ptr = '\0';
}


/* read and parse global QM/MM configuration file and
   store the result in a qmmm_config_t struct */
int read_qmmm_config(const char *file, qmmm_config_t *cfg)
{
    FILE *fp;
    char *ptr;
    int i, lineno, len;
    char buf[BUF_SIZE];

    /* need to have config file and storage for config provided */
    if ((file == NULL) || (cfg == NULL))
        return QMMM_ERROR;

    /* clear config */
    memset(cfg, 0, sizeof(qmmm_config_t));

    fp = fopen(file,"r");
    if (fp == NULL) return QMMM_ERROR;

    lineno = 0;
    while (fgets(buf,BUF_SIZE,fp)) {
        ++lineno;

        /* remove comments */
        ptr = strchr(buf,'#');
        if (ptr != NULL) *ptr = '\0';

        /* skip over leading whitespace */
        ptr = skip_whitespace(buf);
        len = strlen(ptr);
        if (len == 0) continue;

        /* convert keyword to lower case */
        for (i=0; i < len; ++i) {
            if ((ptr[i]=='\0') || (ptr[i]==' ') || (ptr[i]=='\t'))
                break;
            ptr[i] = tolower(ptr[i]);
        }   

        /* handle keywords */
        if (strncmp(ptr,"mode",4) == 0) {
            ptr = skip_whitespace(ptr+4);
            if ((*ptr=='m') || (*ptr=='M'))
                cfg->qmmm_mode = QMMM_MODE_MECH;
            else if ((*ptr=='e') || (*ptr=='E'))
                cfg->qmmm_mode = QMMM_MODE_ELEC;
            else if ((*ptr=='o') || (*ptr=='O'))
                cfg->qmmm_mode = QMMM_MODE_OFF;
            else cfg->qmmm_mode = QMMM_ERROR;

        } else if (strncmp(ptr,"comm",4) == 0) {
            ptr = skip_whitespace(ptr+4);
            if ((*ptr=='m') || (*ptr=='M'))
                cfg->comm_mode = QMMM_COMM_MPI;
            else if ((*ptr=='s') || (*ptr=='S'))
                cfg->comm_mode = QMMM_COMM_SHM;
            else cfg->comm_mode = QMMM_ERROR;

        } else if (strncmp(ptr,"steps",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            cfg->steps = atoi(ptr);

        } else if (strncmp(ptr,"verbose",7) == 0) {
            ptr = skip_whitespace(ptr+7);
            trim_whitespace(ptr);
            cfg->verbose = atoi(ptr);

        } else if (strncmp(ptr,"handle",6) == 0) {
            ptr = skip_whitespace(ptr+6);
            trim_whitespace(ptr);
            if (cfg->handle) free(cfg->handle);
            cfg->handle = strdup(ptr);

        } else if (strncmp(ptr,"restart",7) == 0) {
            ptr = skip_whitespace(ptr+7);
            trim_whitespace(ptr);
            if (cfg->restart) free(cfg->restart);
            cfg->restart = strdup(ptr);

        } else if (strncmp(ptr,"qmdir",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->qmdir) free(cfg->qmdir);
            cfg->qmdir = strdup(ptr);

        } else if (strncmp(ptr,"madir",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->madir) free(cfg->madir);
            cfg->madir = strdup(ptr);

        } else if (strncmp(ptr,"sldir",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->sldir) free(cfg->sldir);
            cfg->sldir = strdup(ptr);

        } else if (strncmp(ptr,"qminp",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->qminp) free(cfg->qminp);
            cfg->qminp = strdup(ptr);

        } else if (strncmp(ptr,"mainp",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->mainp) free(cfg->mainp);
            cfg->mainp = strdup(ptr);

        } else if (strncmp(ptr,"slinp",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->slinp) free(cfg->slinp);
            cfg->slinp = strdup(ptr);

        } else if (strncmp(ptr,"qmout",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->qmout) free(cfg->qmout);
            cfg->qmout = strdup(ptr);

        } else if (strncmp(ptr,"maout",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->maout) free(cfg->maout);
            cfg->maout = strdup(ptr);

        } else if (strncmp(ptr,"slout",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            trim_whitespace(ptr);
            if (cfg->slout) free(cfg->slout);
            cfg->slout = strdup(ptr);

        } else if (strncmp(ptr,"qmcmd",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            if (cfg->qmcmd) free(cfg->qmcmd);
            cfg->qmcmd = strdup(ptr);

        } else if (strncmp(ptr,"macmd",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            if (cfg->macmd) free(cfg->macmd);
            cfg->macmd = strdup(ptr);

        } else if (strncmp(ptr,"slcmd",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            if (cfg->slcmd) free(cfg->slcmd);
            cfg->slcmd = strdup(ptr);

        } else if (strncmp(ptr,"qmarg",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            if (cfg->qmarg) free(cfg->qmarg);
            cfg->qmarg = strdup(ptr);

        } else if (strncmp(ptr,"maarg",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            if (cfg->maarg) free(cfg->maarg);
            cfg->maarg = strdup(ptr);

        } else if (strncmp(ptr,"slarg",5) == 0) {
            ptr = skip_whitespace(ptr+5);
            if (cfg->slarg) free(cfg->slarg);
            cfg->slarg = strdup(ptr);

        } else {
            fprintf(stderr,"\nParse error in line %d of file %s\n>>%s<<\n\n",
                    lineno, file, buf);
            return QMMM_ERROR;
        }
        
        if (feof(fp)) break;
    }

    fclose(fp);
    return QMMM_OK;
}


/* perform consistency checks on a qmmm_config_t struct */
const char *check_qmmm_config(qmmm_config_t *cfg) {
    const char *msg = NULL;

    if (cfg->qmmm_mode == QMMM_MODE_NONE) {
        msg = "QM/MM coupling mode not set";
    } else if (cfg->steps < 1) {
        msg = "Number of QM/MM steps must be > 0";
    } else if (cfg->qmdir == NULL) {
        msg = "QM directory not set";
    } else if (cfg->madir == NULL) {
        msg = "MM master directory not set";
    } else if (cfg->sldir == NULL) {
        msg = "MM slave directory not set";
    } else if (cfg->qminp == NULL) {
        msg = "QM input file not set";
    } else if (cfg->mainp == NULL) {
        msg = "MM master input file not set";
    } else if (cfg->slinp == NULL) {
        msg = "MM slave input file not set";
    } else if (cfg->qmout == NULL) {
        msg = "QM output file not set";
    } else if (cfg->maout == NULL) {
        msg = "MM master output file not set";
    } else if (cfg->slout == NULL) {
        msg = "MM slave output file not set";    
    }
    return msg;
}

/* take the existing global QM/MM configuration and write it to a file */
int write_qmmm_config(const char *file, qmmm_config_t *cfg)
{
    FILE *fp;
    time_t now;

    if (cfg == NULL)
        return QMMM_ERROR;

    if (file == NULL)
        fp = stderr;
    else {
        fp = fopen(file,"w");
        if (fp == NULL) return QMMM_ERROR;
    }

    now = time(NULL);
    fprintf(fp,"# QM/MM global configuration. %s\n",ctime(&now));
    fprintf(fp,"# comm_mode: %d\n", cfg->comm_mode);
    if (cfg->qmmm_mode == QMMM_MODE_OFF) {
        fputs("mode off\n",fp);
    } else if (cfg->qmmm_mode == QMMM_MODE_NONE) {
        fputs("mode none\n",fp);
    } else if (cfg->qmmm_mode == QMMM_MODE_MECH) {
        fputs("mode mechanical\n",fp);
    } else if (cfg->qmmm_mode == QMMM_MODE_ELEC) {
        fputs("mode electrostatic\n",fp);
    } else fputs("mode unknown\n",fp);

    if (cfg->restart) fprintf(fp,"restart %s\n",cfg->restart);
    if (cfg->handle)  fprintf(fp,"handle %s\n", cfg->handle);
    fprintf(fp,"steps %d\n",cfg->steps);
    fprintf(fp,"verbose %d\n",cfg->verbose);

    fputs("\n# QM system configuration:\n",fp);
    if (cfg->qmdir) fprintf(fp,"qmdir %s\n", cfg->qmdir);
    else fputs("qmdir (required input not set)\n",fp);

    if (cfg->qminp) fprintf(fp,"qminp %s\n", cfg->qminp);
    else fputs("qminp (required input not set)\n",fp);

    if (cfg->qmout) fprintf(fp,"qmout %s\n", cfg->qmout);
    else fputs("qmout NULL\n",fp);

    if (cfg->qmcmd) fprintf(fp,"qmcmd %s\n", cfg->qmcmd);
    if (cfg->qmarg) fprintf(fp,"qmarg %s\n", cfg->qmarg);

    fputs("\n# MM master system configuration:\n",fp);
    if (cfg->madir) fprintf(fp,"madir %s\n", cfg->madir);
    else fputs("madir (required input not set)\n",fp);

    if (cfg->mainp) fprintf(fp,"mainp %s\n", cfg->mainp);
    else fputs("mainp (required input not set)\n",fp);

    if (cfg->maout) fprintf(fp,"maout %s\n", cfg->maout);
    else fputs("maout NULL\n",fp);

    if (cfg->macmd) fprintf(fp,"macmd %s\n", cfg->macmd);
    if (cfg->maarg) fprintf(fp,"maarg %s\n", cfg->maarg);

    fputs("\n# MM slave system configuration:\n",fp);
    if (cfg->sldir) fprintf(fp,"sldir %s\n", cfg->sldir);
    else fputs("sldir (required input not set)\n",fp);

    if (cfg->slinp) fprintf(fp,"slinp %s\n", cfg->slinp);
    else fputs("slinp (required input not set)\n",fp);

    if (cfg->slout) fprintf(fp,"slout %s\n", cfg->slout);
    else fputs("slout NULL\n",fp);

    if (cfg->slcmd) fprintf(fp,"slcmd %s\n", cfg->slcmd);
    if (cfg->slarg) fprintf(fp,"slarg %s\n", cfg->slarg);

    if (file != NULL) fclose(fp);
    return QMMM_OK;
}

/* free storage associated with qmmm_config_t struct */
void free_qmmm_config(qmmm_config_t *cfg) {

    if (cfg->qmdir != NULL) free(cfg->qmdir);
    if (cfg->madir != NULL) free(cfg->madir);
    if (cfg->sldir != NULL) free(cfg->sldir);

    if (cfg->qminp != NULL) free(cfg->qminp);
    if (cfg->mainp != NULL) free(cfg->mainp);
    if (cfg->slinp != NULL) free(cfg->slinp);

    if (cfg->qmout != NULL) free(cfg->qmout);
    if (cfg->maout != NULL) free(cfg->maout);
    if (cfg->slout != NULL) free(cfg->slout);

    if (cfg->qmcmd != NULL) free(cfg->qmcmd);
    if (cfg->macmd != NULL) free(cfg->macmd);
    if (cfg->slcmd != NULL) free(cfg->slcmd);

    if (cfg->qmarg != NULL) free(cfg->qmarg);
    if (cfg->maarg != NULL) free(cfg->maarg);
    if (cfg->slarg != NULL) free(cfg->slarg);
}

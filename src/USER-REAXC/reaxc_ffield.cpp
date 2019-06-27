/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include "pair_reaxc.h"
#include "error.h"
#include "reaxc_ffield.h"
#include "reaxc_tool_box.h"


char Read_Force_Field( FILE *fp, reax_interaction *reax,
                       control_params *control )
{
  char    *s;
  char   **tmp;
  char ****tor_flag;
  int      c, i, j, k, l, m, n, o, p, cnt;
  int lgflag = control->lgflag;
  int errorflag = 1;
  double     val;
  MPI_Comm comm;
  int me;
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &me);

  s = (char*) malloc(sizeof(char)*MAX_LINE);
  tmp = (char**) malloc(sizeof(char*)*MAX_TOKENS);
  for (i=0; i < MAX_TOKENS; i++)
    tmp[i] = (char*) malloc(sizeof(char)*MAX_TOKEN_LEN);

  /* reading first header comment */
  fgets( s, MAX_LINE, fp );

  /* line 2 is number of global parameters */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );

  /* reading the number of global parameters */
  n = atoi(tmp[0]);
  if (n < 1) {
    if (me == 0)
      control->error_ptr->warning( FLERR, "Number of globals in ffield file is 0. The file will not be read." );
    fclose(fp);
    free(s);
    free(tmp);
    return 1;
  }

  reax->gp.n_global = n;
  reax->gp.l = (double*) malloc(sizeof(double)*n);

  /* see reax_types.h for mapping between l[i] and the lambdas used in ff */
  for (i=0; i < n; i++) {
    fgets(s,MAX_LINE,fp);
    c = Tokenize(s,&tmp);

    val = (double) atof(tmp[0]);
    reax->gp.l[i] = val;
  }

  control->bo_cut    = 0.01 * reax->gp.l[29];
  control->nonb_low  = reax->gp.l[11];
  control->nonb_cut  = reax->gp.l[12];

  /* next line is number of atom types and some comments */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  reax->num_atom_types = atoi(tmp[0]);

  /* 3 lines of comments */
  fgets(s,MAX_LINE,fp);
  fgets(s,MAX_LINE,fp);
  fgets(s,MAX_LINE,fp);

  /* Allocating structures in reax_interaction */
  reax->sbp = (single_body_parameters*)
    scalloc(control->error_ptr,  reax->num_atom_types, sizeof(single_body_parameters), "sbp");
  reax->tbp = (two_body_parameters**)
    scalloc(control->error_ptr,  reax->num_atom_types, sizeof(two_body_parameters*), "tbp");
  reax->thbp= (three_body_header***)
    scalloc(control->error_ptr,  reax->num_atom_types, sizeof(three_body_header**), "thbp");
  reax->hbp = (hbond_parameters***)
    scalloc(control->error_ptr,  reax->num_atom_types, sizeof(hbond_parameters**), "hbp");
  reax->fbp = (four_body_header****)
    scalloc(control->error_ptr,  reax->num_atom_types, sizeof(four_body_header***), "fbp");
  tor_flag  = (char****)
    scalloc(control->error_ptr,  reax->num_atom_types, sizeof(char***), "tor_flag");

  for( i = 0; i < reax->num_atom_types; i++ ) {
    reax->tbp[i] = (two_body_parameters*)
      scalloc(control->error_ptr,  reax->num_atom_types, sizeof(two_body_parameters), "tbp[i]");
    reax->thbp[i]= (three_body_header**)
      scalloc(control->error_ptr,  reax->num_atom_types, sizeof(three_body_header*), "thbp[i]");
    reax->hbp[i] = (hbond_parameters**)
      scalloc(control->error_ptr,  reax->num_atom_types, sizeof(hbond_parameters*), "hbp[i]");
    reax->fbp[i] = (four_body_header***)
      scalloc(control->error_ptr,  reax->num_atom_types, sizeof(four_body_header**), "fbp[i]");
    tor_flag[i]  = (char***)
      scalloc(control->error_ptr,  reax->num_atom_types, sizeof(char**), "tor_flag[i]");

    for( j = 0; j < reax->num_atom_types; j++ ) {
      reax->thbp[i][j]= (three_body_header*)
        scalloc(control->error_ptr,  reax->num_atom_types, sizeof(three_body_header), "thbp[i,j]");
      reax->hbp[i][j] = (hbond_parameters*)
        scalloc(control->error_ptr,  reax->num_atom_types, sizeof(hbond_parameters), "hbp[i,j]");
      reax->fbp[i][j] = (four_body_header**)
        scalloc(control->error_ptr,  reax->num_atom_types, sizeof(four_body_header*), "fbp[i,j]");
      tor_flag[i][j]  = (char**)
        scalloc(control->error_ptr,  reax->num_atom_types, sizeof(char*), "tor_flag[i,j]");

      for (k=0; k < reax->num_atom_types; k++) {
        reax->fbp[i][j][k] = (four_body_header*)
          scalloc(control->error_ptr,  reax->num_atom_types, sizeof(four_body_header), "fbp[i,j,k]");
        tor_flag[i][j][k]  = (char*)
          scalloc(control->error_ptr,  reax->num_atom_types, sizeof(char), "tor_flag[i,j,k]");
      }
    }
  }

  reax->gp.vdw_type = 0;

  char errmsg[1024];

  for( i = 0; i < reax->num_atom_types; i++ ) {
    /* line one */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    /* Sanity checks */
    if (c == 2 && !lgflag)
        control->error_ptr->all(FLERR, "Force field file requires using 'lgvdw yes'");	

    if (c < 9) {
      snprintf (errmsg, 1024, "Missing parameter(s) in line %s", s);
      control->error_ptr->all(FLERR, errmsg);
    }

    for( j = 0; j < (int)(strlen(tmp[0])); ++j )
      reax->sbp[i].name[j] = toupper( tmp[0][j] );

    val = atof(tmp[1]); reax->sbp[i].r_s        = val;
    val = atof(tmp[2]); reax->sbp[i].valency    = val;
    val = atof(tmp[3]); reax->sbp[i].mass       = val;
    val = atof(tmp[4]); reax->sbp[i].r_vdw      = val;
    val = atof(tmp[5]); reax->sbp[i].epsilon    = val;
    val = atof(tmp[6]); reax->sbp[i].gamma      = val;
    val = atof(tmp[7]); reax->sbp[i].r_pi       = val;
    val = atof(tmp[8]); reax->sbp[i].valency_e  = val;
    reax->sbp[i].nlp_opt = 0.5 * (reax->sbp[i].valency_e-reax->sbp[i].valency);

    /* line two */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    /* Sanity check */
    if (c < 8) {
      snprintf (errmsg, 1024, "Missing parameter(s) in line %s", s);
      control->error_ptr->all(FLERR, errmsg);
    }

    val = atof(tmp[0]); reax->sbp[i].alpha      = val;
    val = atof(tmp[1]); reax->sbp[i].gamma_w    = val;
    val = atof(tmp[2]); reax->sbp[i].valency_boc= val;
    val = atof(tmp[3]); reax->sbp[i].p_ovun5    = val;
    val = atof(tmp[4]);
    val = atof(tmp[5]); reax->sbp[i].chi        = val;
    val = atof(tmp[6]); reax->sbp[i].eta        = 2.0 * val;
    val = atof(tmp[7]); reax->sbp[i].p_hbond = (int) val;

    /* line 3 */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    /* Sanity check */
    if (c < 8) {
      snprintf (errmsg, 1024, "Missing parameter(s) in line %s", s);
      control->error_ptr->all(FLERR, errmsg);
    }

    val = atof(tmp[0]); reax->sbp[i].r_pi_pi    = val;
    val = atof(tmp[1]); reax->sbp[i].p_lp2      = val;
    val = atof(tmp[2]);
    val = atof(tmp[3]); reax->sbp[i].b_o_131    = val;
    val = atof(tmp[4]); reax->sbp[i].b_o_132    = val;
    val = atof(tmp[5]); reax->sbp[i].b_o_133    = val;
    val = atof(tmp[6]);
    val = atof(tmp[7]);

    /* line 4  */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    /* Sanity check */
    if (c < 8) {
      snprintf (errmsg, 1024, "Missing parameter(s) in line %s", s);
      control->error_ptr->all(FLERR, errmsg);
    }

    val = atof(tmp[0]); reax->sbp[i].p_ovun2    = val;
    val = atof(tmp[1]); reax->sbp[i].p_val3     = val;
    val = atof(tmp[2]);
    val = atof(tmp[3]); reax->sbp[i].valency_val= val;
    val = atof(tmp[4]); reax->sbp[i].p_val5     = val;
    val = atof(tmp[5]); reax->sbp[i].rcore2     = val;
    val = atof(tmp[6]); reax->sbp[i].ecore2     = val;
    val = atof(tmp[7]); reax->sbp[i].acore2     = val;

    /* line 5, only if lgvdw is yes */
    if (lgflag) {
      fgets( s, MAX_LINE, fp );
      c = Tokenize( s, &tmp );

      /* Sanity check */
      if (c > 2) {
        control->error_ptr->all(FLERR,"Force field file incompatible with 'lgvdw yes'");
      }

      val = atof(tmp[0]); reax->sbp[i].lgcij           = val;
      val = atof(tmp[1]); reax->sbp[i].lgre           = val;
    }

    if (reax->sbp[i].rcore2>0.01 && reax->sbp[i].acore2>0.01) { // Inner-wall
      if (reax->sbp[i].gamma_w>0.5) { // Shielding vdWaals
        if (reax->gp.vdw_type != 0 && reax->gp.vdw_type != 3) {
          if (errorflag && (me == 0)) {
            char errmsg[512];
            snprintf(errmsg, 512, "VdWaals-parameters for element %s "
                   "indicate inner wall+shielding, but earlier "
                   "atoms indicate different vdWaals-method. "
                   "This may cause division-by-zero errors. "
                   "Keeping vdWaals-setting for earlier atoms.",
                   reax->sbp[i].name);
            control->error_ptr->warning(FLERR,errmsg);
          }
          errorflag = 0;
        } else {
          reax->gp.vdw_type = 3;
        }
      } else {  // No shielding vdWaals parameters present
        if (reax->gp.vdw_type != 0 && reax->gp.vdw_type != 2) {
          if (me == 0) {
            char errmsg[512];
            snprintf(errmsg, 512, "VdWaals-parameters for element %s "
                   "indicate inner wall without shielding, but earlier "
                   "atoms indicate different vdWaals-method. "
                   "This may cause division-by-zero errors. "
                   "Keeping vdWaals-setting for earlier atoms.",
                   reax->sbp[i].name);
            control->error_ptr->warning(FLERR,errmsg);
          }
        } else {
          reax->gp.vdw_type = 2;
        }
      }
    } else { // No Inner wall parameters present
      if (reax->sbp[i].gamma_w>0.5) { // Shielding vdWaals
        if (reax->gp.vdw_type != 0 && reax->gp.vdw_type != 1) {
          if (me == 0) {
            char errmsg[512];
            snprintf(errmsg, 512, "VdWaals parameters for element %s "
                   "indicate shielding without inner wall, but earlier "
                   "elements indicate different vdWaals-method. "
                   "This may cause division-by-zero errors. "
                   "Keeping vdWaals-setting for earlier atoms.",
                   reax->sbp[i].name);
            control->error_ptr->warning(FLERR,errmsg);
          }
        } else {
          reax->gp.vdw_type = 1;
        }
      } else {
        char errmsg[256];
        snprintf(errmsg, 256, "Inconsistent vdWaals-parameters: "
                 "No shielding or inner-wall set for element %s",
                 reax->sbp[i].name);
        control->error_ptr->all(FLERR, errmsg);
      }
    }
  }

  /* Equate vval3 to valf for first-row elements (25/10/2004) */
  for( i = 0; i < reax->num_atom_types; i++ )
    if( reax->sbp[i].mass < 21 &&
        reax->sbp[i].valency_val != reax->sbp[i].valency_boc ) {
      if (me == 0) {
        char errmsg[256];
        snprintf(errmsg, 256, "Changed valency_val to valency_boc for %s",
               reax->sbp[i].name);
        control->error_ptr->warning(FLERR,errmsg);
      }
      reax->sbp[i].valency_val = reax->sbp[i].valency_boc;
    }

  /* next line is number of two body combination and some comments */
  fgets(s,MAX_LINE,fp);
  c=Tokenize(s,&tmp);

  if (c == 2 && !lgflag) {
      control->error_ptr->all(FLERR, "Force field file requires using 'lgvdw yes'");
    }

  l = atoi(tmp[0]);

  /* a line of comments */
  fgets(s,MAX_LINE,fp);

  for (i=0; i < l; i++) {
    /* line 1 */
    fgets(s,MAX_LINE,fp);
    c=Tokenize(s,&tmp);

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;

    if (j < reax->num_atom_types && k < reax->num_atom_types) {

      val = atof(tmp[2]); reax->tbp[j][k].De_s      = val;
      reax->tbp[k][j].De_s      = val;
      val = atof(tmp[3]); reax->tbp[j][k].De_p      = val;
      reax->tbp[k][j].De_p      = val;
      val = atof(tmp[4]); reax->tbp[j][k].De_pp     = val;
      reax->tbp[k][j].De_pp     = val;
      val = atof(tmp[5]); reax->tbp[j][k].p_be1     = val;
      reax->tbp[k][j].p_be1     = val;
      val = atof(tmp[6]); reax->tbp[j][k].p_bo5     = val;
      reax->tbp[k][j].p_bo5     = val;
      val = atof(tmp[7]); reax->tbp[j][k].v13cor    = val;
      reax->tbp[k][j].v13cor    = val;

      val = atof(tmp[8]); reax->tbp[j][k].p_bo6     = val;
      reax->tbp[k][j].p_bo6     = val;
      val = atof(tmp[9]); reax->tbp[j][k].p_ovun1 = val;
      reax->tbp[k][j].p_ovun1 = val;

      /* line 2 */
      fgets(s,MAX_LINE,fp);
      c=Tokenize(s,&tmp);

      val = atof(tmp[0]); reax->tbp[j][k].p_be2     = val;
      reax->tbp[k][j].p_be2     = val;
      val = atof(tmp[1]); reax->tbp[j][k].p_bo3     = val;
      reax->tbp[k][j].p_bo3     = val;
      val = atof(tmp[2]); reax->tbp[j][k].p_bo4     = val;
      reax->tbp[k][j].p_bo4     = val;
      val = atof(tmp[3]);

      val = atof(tmp[4]); reax->tbp[j][k].p_bo1     = val;
      reax->tbp[k][j].p_bo1     = val;
      val = atof(tmp[5]); reax->tbp[j][k].p_bo2     = val;
      reax->tbp[k][j].p_bo2     = val;
      val = atof(tmp[6]); reax->tbp[j][k].ovc       = val;
      reax->tbp[k][j].ovc       = val;

      val = atof(tmp[7]);
    }
  }

  for (i=0; i < reax->num_atom_types; i++)
    for (j=i; j < reax->num_atom_types; j++) {
      reax->tbp[i][j].r_s = 0.5 *
        (reax->sbp[i].r_s + reax->sbp[j].r_s);
      reax->tbp[j][i].r_s = 0.5 *
        (reax->sbp[j].r_s + reax->sbp[i].r_s);

      reax->tbp[i][j].r_p = 0.5 *
        (reax->sbp[i].r_pi + reax->sbp[j].r_pi);
      reax->tbp[j][i].r_p = 0.5 *
        (reax->sbp[j].r_pi + reax->sbp[i].r_pi);

      reax->tbp[i][j].r_pp = 0.5 *
        (reax->sbp[i].r_pi_pi + reax->sbp[j].r_pi_pi);
      reax->tbp[j][i].r_pp = 0.5 *
        (reax->sbp[j].r_pi_pi + reax->sbp[i].r_pi_pi);


      reax->tbp[i][j].p_boc3 =
        sqrt(reax->sbp[i].b_o_132 *
             reax->sbp[j].b_o_132);
      reax->tbp[j][i].p_boc3 =
        sqrt(reax->sbp[j].b_o_132 *
             reax->sbp[i].b_o_132);

      reax->tbp[i][j].p_boc4 =
        sqrt(reax->sbp[i].b_o_131 *
             reax->sbp[j].b_o_131);
      reax->tbp[j][i].p_boc4 =
        sqrt(reax->sbp[j].b_o_131 *
             reax->sbp[i].b_o_131);

      reax->tbp[i][j].p_boc5 =
        sqrt(reax->sbp[i].b_o_133 *
             reax->sbp[j].b_o_133);
      reax->tbp[j][i].p_boc5 =
        sqrt(reax->sbp[j].b_o_133 *
             reax->sbp[i].b_o_133);


      reax->tbp[i][j].D =
        sqrt(reax->sbp[i].epsilon *
             reax->sbp[j].epsilon);

      reax->tbp[j][i].D =
        sqrt(reax->sbp[j].epsilon *
             reax->sbp[i].epsilon);

      reax->tbp[i][j].alpha =
        sqrt(reax->sbp[i].alpha *
             reax->sbp[j].alpha);

      reax->tbp[j][i].alpha =
        sqrt(reax->sbp[j].alpha *
             reax->sbp[i].alpha);

      reax->tbp[i][j].r_vdW =
        2.0 * sqrt(reax->sbp[i].r_vdw * reax->sbp[j].r_vdw);

      reax->tbp[j][i].r_vdW =
        2.0 * sqrt(reax->sbp[j].r_vdw * reax->sbp[i].r_vdw);

      reax->tbp[i][j].gamma_w =
        sqrt(reax->sbp[i].gamma_w *
             reax->sbp[j].gamma_w);

      reax->tbp[j][i].gamma_w =
        sqrt(reax->sbp[j].gamma_w *
             reax->sbp[i].gamma_w);

      reax->tbp[i][j].gamma =
        pow(reax->sbp[i].gamma *
            reax->sbp[j].gamma,-1.5);

      reax->tbp[j][i].gamma =
        pow(reax->sbp[j].gamma *
            reax->sbp[i].gamma,-1.5);

      // additions for additional vdWaals interaction types - inner core

      reax->tbp[i][j].rcore = reax->tbp[j][i].rcore =
        sqrt( reax->sbp[i].rcore2 * reax->sbp[j].rcore2 );

      reax->tbp[i][j].ecore = reax->tbp[j][i].ecore =
        sqrt( reax->sbp[i].ecore2 * reax->sbp[j].ecore2 );

      reax->tbp[i][j].acore = reax->tbp[j][i].acore =
        sqrt( reax->sbp[i].acore2 * reax->sbp[j].acore2 );

      // additions for additional vdWalls interaction types lg correction

      reax->tbp[i][j].lgcij = reax->tbp[j][i].lgcij =
        sqrt( reax->sbp[i].lgcij * reax->sbp[j].lgcij );

      reax->tbp[i][j].lgre = reax->tbp[j][i].lgre = 2.0 * reax->gp.l[35] *
        sqrt( reax->sbp[i].lgre*reax->sbp[j].lgre );

    }

  fgets(s,MAX_LINE,fp);
  c=Tokenize(s,&tmp);
  l = atoi(tmp[0]);

  for (i=0; i < l; i++) {
    fgets(s,MAX_LINE,fp);
    c=Tokenize(s,&tmp);

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;

    if (j < reax->num_atom_types && k < reax->num_atom_types)        {
      val = atof(tmp[2]);
      if (val > 0.0) {
        reax->tbp[j][k].D = val;
        reax->tbp[k][j].D = val;
      }

      val = atof(tmp[3]);
      if (val > 0.0) {
        reax->tbp[j][k].r_vdW = 2 * val;
        reax->tbp[k][j].r_vdW = 2 * val;
      }

      val = atof(tmp[4]);
      if (val > 0.0) {
        reax->tbp[j][k].alpha = val;
        reax->tbp[k][j].alpha = val;
      }

      val = atof(tmp[5]);
      if (val > 0.0) {
        reax->tbp[j][k].r_s = val;
        reax->tbp[k][j].r_s = val;
      }

      val = atof(tmp[6]);
      if (val > 0.0) {
        reax->tbp[j][k].r_p = val;
        reax->tbp[k][j].r_p = val;
      }

      val = atof(tmp[7]);
      if (val > 0.0) {
        reax->tbp[j][k].r_pp = val;
        reax->tbp[k][j].r_pp = val;
      }

      val = atof(tmp[8]);
      if (val >= 0.0) {
        reax->tbp[j][k].lgcij = val;
        reax->tbp[k][j].lgcij = val;
      }
    }
  }

  for( i = 0; i < reax->num_atom_types; ++i )
    for( j = 0; j < reax->num_atom_types; ++j )
      for( k = 0; k < reax->num_atom_types; ++k )
        reax->thbp[i][j][k].cnt = 0;

  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  l = atoi( tmp[0] );

  for( i = 0; i < l; i++ ) {
    fgets(s,MAX_LINE,fp);
    c=Tokenize(s,&tmp);

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;
    m = atoi(tmp[2]) - 1;

    if (j < reax->num_atom_types && k < reax->num_atom_types &&
        m < reax->num_atom_types) {
      cnt = reax->thbp[j][k][m].cnt;
      reax->thbp[j][k][m].cnt++;
      reax->thbp[m][k][j].cnt++;

      val = atof(tmp[3]);
      reax->thbp[j][k][m].prm[cnt].theta_00 = val;
      reax->thbp[m][k][j].prm[cnt].theta_00 = val;

      val = atof(tmp[4]);
      reax->thbp[j][k][m].prm[cnt].p_val1 = val;
      reax->thbp[m][k][j].prm[cnt].p_val1 = val;

      val = atof(tmp[5]);
      reax->thbp[j][k][m].prm[cnt].p_val2 = val;
      reax->thbp[m][k][j].prm[cnt].p_val2 = val;

      val = atof(tmp[6]);
      reax->thbp[j][k][m].prm[cnt].p_coa1 = val;
      reax->thbp[m][k][j].prm[cnt].p_coa1 = val;

      val = atof(tmp[7]);
      reax->thbp[j][k][m].prm[cnt].p_val7 = val;
      reax->thbp[m][k][j].prm[cnt].p_val7 = val;

      val = atof(tmp[8]);
      reax->thbp[j][k][m].prm[cnt].p_pen1 = val;
      reax->thbp[m][k][j].prm[cnt].p_pen1 = val;

      val = atof(tmp[9]);
      reax->thbp[j][k][m].prm[cnt].p_val4 = val;
      reax->thbp[m][k][j].prm[cnt].p_val4 = val;
    }
  }

  /* clear all entries first */
  for( i = 0; i < reax->num_atom_types; ++i )
    for( j = 0; j < reax->num_atom_types; ++j )
      for( k = 0; k < reax->num_atom_types; ++k )
        for( m = 0; m < reax->num_atom_types; ++m ) {
          reax->fbp[i][j][k][m].cnt = 0;
          tor_flag[i][j][k][m] = 0;
        }

  /* next line is number of 4-body params and some comments */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  l = atoi( tmp[0] );

  for( i = 0; i < l; i++ ) {
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;
    m = atoi(tmp[2]) - 1;
    n = atoi(tmp[3]) - 1;

    if (j >= 0 && n >= 0) { // this means the entry is not in compact form
      if (j < reax->num_atom_types && k < reax->num_atom_types &&
          m < reax->num_atom_types && n < reax->num_atom_types) {
        tor_flag[j][k][m][n] = 1;
        tor_flag[n][m][k][j] = 1;

        reax->fbp[j][k][m][n].cnt = 1;
        reax->fbp[n][m][k][j].cnt = 1;

        val = atof(tmp[4]);
        reax->fbp[j][k][m][n].prm[0].V1 = val;
        reax->fbp[n][m][k][j].prm[0].V1 = val;

        val = atof(tmp[5]);
        reax->fbp[j][k][m][n].prm[0].V2 = val;
        reax->fbp[n][m][k][j].prm[0].V2 = val;

        val = atof(tmp[6]);
        reax->fbp[j][k][m][n].prm[0].V3 = val;
        reax->fbp[n][m][k][j].prm[0].V3 = val;

        val = atof(tmp[7]);
        reax->fbp[j][k][m][n].prm[0].p_tor1 = val;
        reax->fbp[n][m][k][j].prm[0].p_tor1 = val;

        val = atof(tmp[8]);
        reax->fbp[j][k][m][n].prm[0].p_cot1 = val;
        reax->fbp[n][m][k][j].prm[0].p_cot1 = val;
      }
    } else { /* This means the entry is of the form 0-X-Y-0 */
      if (k < reax->num_atom_types && m < reax->num_atom_types)
        for( p = 0; p < reax->num_atom_types; p++ )
          for( o = 0; o < reax->num_atom_types; o++ ) {
            reax->fbp[p][k][m][o].cnt = 1;
            reax->fbp[o][m][k][p].cnt = 1;

            if (tor_flag[p][k][m][o] == 0) {
              reax->fbp[p][k][m][o].prm[0].V1 = atof(tmp[4]);
              reax->fbp[p][k][m][o].prm[0].V2 = atof(tmp[5]);
              reax->fbp[p][k][m][o].prm[0].V3 = atof(tmp[6]);
              reax->fbp[p][k][m][o].prm[0].p_tor1 = atof(tmp[7]);
              reax->fbp[p][k][m][o].prm[0].p_cot1 = atof(tmp[8]);
            }

            if (tor_flag[o][m][k][p] == 0) {
              reax->fbp[o][m][k][p].prm[0].V1 = atof(tmp[4]);
              reax->fbp[o][m][k][p].prm[0].V2 = atof(tmp[5]);
              reax->fbp[o][m][k][p].prm[0].V3 = atof(tmp[6]);
              reax->fbp[o][m][k][p].prm[0].p_tor1 = atof(tmp[7]);
              reax->fbp[o][m][k][p].prm[0].p_cot1 = atof(tmp[8]);
            }
          }
    }
  }



  /* next line is number of hydrogen bond params and some comments */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  l = atoi( tmp[0] );

  for( i = 0; i < reax->num_atom_types; ++i )
    for( j = 0; j < reax->num_atom_types; ++j )
      for( k = 0; k < reax->num_atom_types; ++k )
        reax->hbp[i][j][k].r0_hb = -1.0;

  for( i = 0; i < l; i++ ) {
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;
    m = atoi(tmp[2]) - 1;


    if (j < reax->num_atom_types && m < reax->num_atom_types) {
      val = atof(tmp[3]);
      reax->hbp[j][k][m].r0_hb = val;

      val = atof(tmp[4]);
      reax->hbp[j][k][m].p_hb1 = val;

      val = atof(tmp[5]);
      reax->hbp[j][k][m].p_hb2 = val;

      val = atof(tmp[6]);
      reax->hbp[j][k][m].p_hb3 = val;
    }
  }

  /* deallocate helper storage */
  for( i = 0; i < MAX_TOKENS; i++ )
    free( tmp[i] );
  free( tmp );
  free( s );


  /* deallocate tor_flag */
  for( i = 0; i < reax->num_atom_types; i++ ) {
    for( j = 0; j < reax->num_atom_types; j++ ) {
      for( k = 0; k < reax->num_atom_types; k++ ) {
        free( tor_flag[i][j][k] );
      }
      free( tor_flag[i][j] );
    }
    free( tor_flag[i] );
  }
  free( tor_flag );

  // close file

  fclose(fp);

  return SUCCESS;
}

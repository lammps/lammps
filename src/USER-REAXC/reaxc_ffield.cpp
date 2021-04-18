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
  <https://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxff_api.h"

#include "error.h"
#include "memory.h"
#include "text_file_reader.h"
#include "utils.h"

#include <cctype>
#include <cmath>
#include <string>

using LAMMPS_NS::utils::open_potential;
using LAMMPS_NS::utils::getsyserror;

namespace ReaxFF {

  class parser_error : public std::exception {
    std::string message;
  public:
    parser_error(const std::string &mesg) { message = mesg; }
    const char *what() const noexcept { return message.c_str(); }
  };

  void Read_Force_Field(const char *filename, reax_interaction *reax,
                        control_params *control, MPI_Comm world)
  {
    char ****tor_flag;
    auto error = control->error_ptr;
    auto lmp = control->lmp_ptr;
    auto memory = control->lmp_ptr->memory;

    // read and parse the force field only on rank 0

#define THROW_ERROR(txt)                                                \
    throw parser_error(fmt::format("{}:{}: {}",filename,lineno,txt))

    if (control->me == 0) {
      FILE *fp = LAMMPS_NS::utils::open_potential(filename, lmp, nullptr);
      if (!fp)
        error->one(FLERR,fmt::format("The ReaxFF parameter file {} cannot be opened: {}",
                                     filename, getsyserror()));
      LAMMPS_NS::TextFileReader reader(fp, "ReaxFF parameter");
      reader.ignore_comments = false;

      try {
        int i,j,k,l,m,n,lineno = 0;

        // skip header comment line

        reader.skip_line();
        ++lineno;

        // set some defaults

        reax->gp.vdw_type = 0;

        // get number of global parameters

        auto values = reader.next_values(0);
        n = values.next_int();
        reax->gp.n_global = n;
        ++lineno;

        if (n < 1)
          THROW_ERROR("Invalid number of global parameters");

        memory->create(reax->gp.l,n,"reaxff:gp.l");

        // see reaxff_types.h for mapping between l[i] and the lambdas used in ff

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          reax->gp.l[i] = values.next_double();
        }

        // next line is number of atom types followed by 3 lines of comments

        values = reader.next_values(0);
        n = values.next_int();
        reax->num_atom_types = n;
        reader.skip_line();
        reader.skip_line();
        reader.skip_line();
        lineno += 4;

        // allocate and clear storage for ffield data

        memory->create(reax->sbp,n,"reaxff:sbp");
        memory->create(reax->tbp,n,n,"reaxff:tbp");
        memory->create(reax->thbp,n,n,n,"reaxff:thbp");
        memory->create(reax->hbp,n,n,n,"reaxff:hbp");
        memory->create(reax->fbp,n,n,n,n,"reaxff:fbp");
        memory->create(tor_flag,n,n,n,n,"reaxff:tor_flag");
        memset(&(reax->sbp[0]),0,sizeof(single_body_parameters)*n);
        memset(&(reax->tbp[0][0]),0,sizeof(two_body_parameters)*n*n);
        memset(&(reax->thbp[0][0][0]),0,sizeof(three_body_header)*n*n*n);
        memset(&(reax->hbp[0][0][0]),0,sizeof(hbond_parameters)*n*n*n);
        memset(&(reax->fbp[0][0][0][0]),0,sizeof(four_body_header)*n*n*n*n);
        memset(&tor_flag[0][0][0][0],0,sizeof(char)*n*n*n*n);

        // atomic parameters
        // four lines per atom type, or 5 if lgvdw != 0
        // the first starts with the symbol and has 9 words
        // the next three have 8 words
        // the fifth will have 2 words, if present

        const int lgflag = control->lgflag;
        const int ntypes = n;
        for (i = 0; i < ntypes; ++i) {

          // line one

          values = reader.next_values(0);
          ++lineno;

          if ((values.count() < 8) && !lgflag)
            THROW_ERROR("This force field file requires using 'lgvdw yes'");
          if (values.count() < 9)
            THROW_ERROR("Invalid force field file format");

          auto element = values.next_string();
          int len = MIN(element.size(),3);
          for (j = 0; j < len; ++j)
            reax->sbp[i].name[j] = toupper(element[j]);
          reax->sbp[i].name[j] = '\0';

          reax->sbp[i].r_s        = values.next_double();
          reax->sbp[i].valency    = values.next_double();
          reax->sbp[i].mass       = values.next_double();
          reax->sbp[i].r_vdw      = values.next_double();
          reax->sbp[i].epsilon    = values.next_double();
          reax->sbp[i].gamma      = values.next_double();
          reax->sbp[i].r_pi       = values.next_double();
          reax->sbp[i].valency_e  = values.next_double();
          reax->sbp[i].nlp_opt = 0.5 * (reax->sbp[i].valency_e-reax->sbp[i].valency);

          // line two

          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 8)
            THROW_ERROR("Invalid force field file format");

          reax->sbp[i].alpha      = values.next_double();
          reax->sbp[i].gamma_w    = values.next_double();
          reax->sbp[i].valency_boc= values.next_double();
          reax->sbp[i].p_ovun5    = values.next_double();
          values.skip();
          reax->sbp[i].chi        = values.next_double();
          reax->sbp[i].eta        = 2.0*values.next_double();
          reax->sbp[i].p_hbond = (int) values.next_double();

          // line three

          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 8)
            THROW_ERROR("Invalid force field file format");

          reax->sbp[i].r_pi_pi    = values.next_double();
          reax->sbp[i].p_lp2      = values.next_double();
          values.skip();
          reax->sbp[i].b_o_131    = values.next_double();
          reax->sbp[i].b_o_132    = values.next_double();
          reax->sbp[i].b_o_133    = values.next_double();

          // line four

          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 8)
            THROW_ERROR("Invalid force field file format");

          reax->sbp[i].p_ovun2    = values.next_double();
          reax->sbp[i].p_val3     = values.next_double();
          values.skip();
          reax->sbp[i].valency_val= values.next_double();
          reax->sbp[i].p_val5     = values.next_double();
          reax->sbp[i].rcore2     = values.next_double();
          reax->sbp[i].ecore2     = values.next_double();
          reax->sbp[i].acore2     = values.next_double();

          // read line five only when lgflag != 0

          if (lgflag) {
            values = reader.next_values(0);
            ++lineno;
            if (values.count() < 2)
              THROW_ERROR("Invalid force field file format");
            reax->sbp[i].lgcij    = values.next_double();
            reax->sbp[i].lgre     = values.next_double();
          } else reax->sbp[i].lgcij = reax->sbp[i].lgre = 0.0;

          // van der Waals settings check:

          // Inner-wall?
          if ((reax->sbp[i].rcore2 > 0.01) && (reax->sbp[i].acore2 > 0.01)) {
            // Shielding van der Waals?
            if (reax->sbp[i].gamma_w > 0.5) {
              if (reax->gp.vdw_type != 0 && reax->gp.vdw_type != 3) {
                const auto errmsg
                  = fmt::format("Van der Waals parameters for element {} "
                                "indicate inner wall+shielding, but earlier "
                                "atoms indicate a different van der Waals "
                                "method. This may cause division-by-zero "
                                "errors. Keeping van der Waals setting for "
                                "earlier atoms.",reax->sbp[i].name);
                error->warning(FLERR,errmsg);
              } else {
                reax->gp.vdw_type = 3;
              }
            } else {  // No shielding van der Waals parameters present
              if ((reax->gp.vdw_type != 0) && (reax->gp.vdw_type != 2)) {
                const auto errmsg
                  = fmt::format("Van der Waals parameters for element {} "
                                "indicate inner wall withou shielding, but "
                                "earlier atoms indicate a different van der "
                                "Waals-method. This may cause division-by-"
                                "zero errors. Keeping van der Waals setting "
                                "for earlier atoms.", reax->sbp[i].name);
                error->warning(FLERR,errmsg);
              } else {
                reax->gp.vdw_type = 2;
              }
            }
          } else { // No Inner wall parameters present
            if (reax->sbp[i].gamma_w > 0.5) { // Shielding vdWaals
              if ((reax->gp.vdw_type != 0) && (reax->gp.vdw_type != 1)) {
                const auto errmsg
                  = fmt::format("Van der Waals parameters for element {} "
                                "indicate shielding without inner wall, but "
                                "earlier elements indicate a different van der "
                                "Waals method. This may cause division-by-zero "
                                "errors. Keeping van der Waals setting for "
                                "earlier atoms.", reax->sbp[i].name);
                error->warning(FLERR,errmsg);
              } else {
                reax->gp.vdw_type = 1;
              }
            } else {
              error->one(FLERR,fmt::format("Inconsistent van der Waals "
                                           "parameters: No shielding or inner "
                                           "wall set for element {}",
                                           reax->sbp[i].name));
            }
          }
        }

        /* Equate vval3 to valf for first-row elements (25/10/2004) */
        for (i = 0; i < ntypes; i++) {
          if ((reax->sbp[i].mass < 21) &&
              (reax->sbp[i].valency_val != reax->sbp[i].valency_boc)) {
            error->warning(FLERR,fmt::format("Changed valency_val to valency"
                                             "_boc for {}", reax->sbp[i].name));
            reax->sbp[i].valency_val = reax->sbp[i].valency_boc;
          }
        }

        // next line is number of two body parameters followed by 1 comment line

        values = reader.next_values(0);
        n = values.next_int();
        reader.skip_line();
        lineno += 2;

        for (i = 0; i < n; ++i) {

          // first line

          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 10)
            THROW_ERROR("Invalid force field file format");

          j = values.next_int() - 1;
          k = values.next_int() - 1;

          if ((j < 0) || (k < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes)) {
            reax->tbp[j][k].De_s    = reax->tbp[k][j].De_s    = values.next_double();
            reax->tbp[j][k].De_p    = reax->tbp[k][j].De_p    = values.next_double();
            reax->tbp[j][k].De_pp   = reax->tbp[k][j].De_pp   = values.next_double();
            reax->tbp[j][k].p_be1   = reax->tbp[k][j].p_be1   = values.next_double();
            reax->tbp[j][k].p_bo5   = reax->tbp[k][j].p_bo5   = values.next_double();
            reax->tbp[j][k].v13cor  = reax->tbp[k][j].v13cor  = values.next_double();
            reax->tbp[j][k].p_bo6   = reax->tbp[k][j].p_bo6   = values.next_double();
            reax->tbp[j][k].p_ovun1 = reax->tbp[k][j].p_ovun1 = values.next_double();
          }

          // second line

          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 8)
            THROW_ERROR("Invalid force field file format");

          if ((j < ntypes) && (k < ntypes)) {
            reax->tbp[j][k].p_be2 = reax->tbp[k][j].p_be2 = values.next_double();
            reax->tbp[j][k].p_bo3 = reax->tbp[k][j].p_bo3 = values.next_double();
            reax->tbp[j][k].p_bo4 = reax->tbp[k][j].p_bo4 = values.next_double();
            values.skip();
            reax->tbp[j][k].p_bo1 = reax->tbp[k][j].p_bo1 = values.next_double();
            reax->tbp[j][k].p_bo2 = reax->tbp[k][j].p_bo2 = values.next_double();
            reax->tbp[j][k].ovc   = reax->tbp[k][j].ovc   = values.next_double();
          }
        }

        for (i=0; i < ntypes; ++i) {
          for (j=i; j < ntypes; ++j) {
            reax->tbp[i][j].r_s = reax->tbp[j][i].r_s =
              0.5*(reax->sbp[j].r_s + reax->sbp[i].r_s);

            reax->tbp[i][j].r_p = reax->tbp[j][i].r_p =
              0.5*(reax->sbp[j].r_pi + reax->sbp[i].r_pi);

            reax->tbp[i][j].r_pp = reax->tbp[j][i].r_pp =
              0.5*(reax->sbp[j].r_pi_pi + reax->sbp[i].r_pi_pi);

            reax->tbp[i][j].p_boc3 = reax->tbp[j][i].p_boc3 =
              sqrt(reax->sbp[j].b_o_132 * reax->sbp[i].b_o_132);

            reax->tbp[i][j].p_boc4 = reax->tbp[j][i].p_boc4 =
              sqrt(reax->sbp[j].b_o_131 * reax->sbp[i].b_o_131);

            reax->tbp[i][j].p_boc5 = reax->tbp[j][i].p_boc5 =
              sqrt(reax->sbp[j].b_o_133 * reax->sbp[i].b_o_133);

            reax->tbp[i][j].D = reax->tbp[j][i].D =
              sqrt(reax->sbp[j].epsilon * reax->sbp[i].epsilon);

            reax->tbp[i][j].alpha = reax->tbp[j][i].alpha =
              sqrt(reax->sbp[j].alpha * reax->sbp[i].alpha);

            reax->tbp[i][j].r_vdW = reax->tbp[j][i].r_vdW =
              2.0*sqrt(reax->sbp[j].r_vdw * reax->sbp[i].r_vdw);

            reax->tbp[i][j].gamma_w = reax->tbp[j][i].gamma_w =
              sqrt(reax->sbp[j].gamma_w * reax->sbp[i].gamma_w);

            reax->tbp[i][j].gamma = reax->tbp[j][i].gamma =
              pow(reax->sbp[j].gamma * reax->sbp[i].gamma,-1.5);

            // additions for additional vdWaals interaction types - inner core

            reax->tbp[i][j].rcore = reax->tbp[j][i].rcore =
              sqrt(reax->sbp[i].rcore2 * reax->sbp[j].rcore2);

            reax->tbp[i][j].ecore = reax->tbp[j][i].ecore =
              sqrt(reax->sbp[i].ecore2 * reax->sbp[j].ecore2);

            reax->tbp[i][j].acore = reax->tbp[j][i].acore =
              sqrt(reax->sbp[i].acore2 * reax->sbp[j].acore2);

            // additions for additional vdWalls interaction types lg correction

            reax->tbp[i][j].lgcij = reax->tbp[j][i].lgcij =
              sqrt(reax->sbp[i].lgcij * reax->sbp[j].lgcij);

            reax->tbp[i][j].lgre = reax->tbp[j][i].lgre
              = 2.0 * reax->gp.l[35] * sqrt(reax->sbp[i].lgre*reax->sbp[j].lgre);
          }
        }

        // next line is number of two body off-diagonal parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        double val;
        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          if ((int)values.count() < 8 + lgflag)
            THROW_ERROR("Invalid force field file format");

          j = values.next_int() - 1;
          k = values.next_int() - 1;

          if ((j < 0) || (k < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes)) {
            val = values.next_double();
            if (val > 0.0) reax->tbp[j][k].D = reax->tbp[k][j].D = val;

            val = values.next_double();
            if (val > 0.0) reax->tbp[j][k].r_vdW = reax->tbp[k][j].r_vdW = 2*val;

            val = values.next_double();
            if (val > 0.0) reax->tbp[j][k].alpha = reax->tbp[k][j].alpha = val;

            val = values.next_double();
            if (val > 0.0) reax->tbp[j][k].r_s = reax->tbp[k][j].r_s = val;

            val = values.next_double();
            if (val > 0.0) reax->tbp[j][k].r_p = reax->tbp[k][j].r_p = val;

            val = values.next_double();
            if (val > 0.0) reax->tbp[j][k].r_pp = reax->tbp[k][j].r_pp = val;

            if (lgflag) {
              val = values.next_double();
              if (val >= 0.0) reax->tbp[j][k].lgcij = reax->tbp[k][j].lgcij = val;
            }
          }
        }

        // next line is number of three body parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        int cnt;
        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 10)
            THROW_ERROR("Invalid force field file format");

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;

          if ((j < 0) || (k < 0) || (l < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes) && (l < ntypes)) {

            cnt = reax->thbp[j][k][l].cnt;
            reax->thbp[j][k][l].cnt++;
            reax->thbp[l][k][j].cnt++;

            val = values.next_double();
            reax->thbp[j][k][l].prm[cnt].theta_00 = val;
            reax->thbp[l][k][j].prm[cnt].theta_00 = val;

            val = values.next_double();
            reax->thbp[j][k][l].prm[cnt].p_val1 = val;
            reax->thbp[l][k][j].prm[cnt].p_val1 = val;

            val = values.next_double();
            reax->thbp[j][k][l].prm[cnt].p_val2 = val;
            reax->thbp[l][k][j].prm[cnt].p_val2 = val;

            val = values.next_double();
            reax->thbp[j][k][l].prm[cnt].p_coa1 = val;
            reax->thbp[l][k][j].prm[cnt].p_coa1 = val;

            val = values.next_double();
            reax->thbp[j][k][l].prm[cnt].p_val7 = val;
            reax->thbp[l][k][j].prm[cnt].p_val7 = val;

            val = values.next_double();
            reax->thbp[j][k][l].prm[cnt].p_pen1 = val;
            reax->thbp[l][k][j].prm[cnt].p_pen1 = val;

            val = values.next_double();
            reax->thbp[j][k][l].prm[cnt].p_val4 = val;
            reax->thbp[l][k][j].prm[cnt].p_val4 = val;
          }
        }

        // next line is number of four body parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 9)
            THROW_ERROR("Invalid force field file format");

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;
          m = values.next_int() - 1;

          if ((j < -1) || (k < 0) || (l < 0) || (m < -1))
            THROW_ERROR("Inconsistent force field file");

          if ((j >= 0) && (m >= 0)) { // this means the entry is not in compact form

            if ((j < ntypes) && (k < ntypes) && (l < ntypes) && (m < ntypes)) {

              tor_flag[j][k][l][m] = 1;
              tor_flag[m][l][k][j] = 1;

              reax->fbp[j][k][l][m].cnt = 1;
              reax->fbp[m][l][k][j].cnt = 1;

              val = values.next_double();
              reax->fbp[j][k][l][m].prm[0].V1 = val;
              reax->fbp[m][l][k][j].prm[0].V1 = val;

              val = values.next_double();
              reax->fbp[j][k][l][m].prm[0].V2 = val;
              reax->fbp[m][l][k][j].prm[0].V2 = val;

              val = values.next_double();
              reax->fbp[j][k][l][m].prm[0].V3 = val;
              reax->fbp[m][l][k][j].prm[0].V3 = val;

              val = values.next_double();
              reax->fbp[j][k][l][m].prm[0].p_tor1 = val;
              reax->fbp[m][l][k][j].prm[0].p_tor1 = val;

              val = values.next_double();
              reax->fbp[j][k][l][m].prm[0].p_cot1 = val;
              reax->fbp[m][l][k][j].prm[0].p_cot1 = val;
            }

          } else { /* This means the entry is of the form 0-X-Y-0 */

            if ((k < ntypes) && (l < ntypes)) {
              const double val1 = values.next_double();
              const double val2 = values.next_double();
              const double val3 = values.next_double();
              const double val4 = values.next_double();
              const double val5 = values.next_double();

              for (int p = 0; p < ntypes; ++p) {
                for (int o = 0; o < ntypes; ++o) {
                  reax->fbp[p][k][l][o].cnt = 1;
                  reax->fbp[o][l][k][p].cnt = 1;

                  if (tor_flag[p][k][l][o] == 0) {
                    reax->fbp[p][k][l][o].prm[0].V1 = val1;
                    reax->fbp[p][k][l][o].prm[0].V2 = val2;
                    reax->fbp[p][k][l][o].prm[0].V3 = val3;
                    reax->fbp[p][k][l][o].prm[0].p_tor1 = val4;
                    reax->fbp[p][k][l][o].prm[0].p_cot1 = val5;
                  }

                  if (tor_flag[o][l][k][p] == 0) {
                    reax->fbp[o][l][k][p].prm[0].V1 = val1;
                    reax->fbp[o][l][k][p].prm[0].V2 = val2;
                    reax->fbp[o][l][k][p].prm[0].V3 = val3;
                    reax->fbp[o][l][k][p].prm[0].p_tor1 = val4;
                    reax->fbp[o][l][k][p].prm[0].p_cot1 = val5;
                  }
                }
              }
            }
          }
        }

        // next line is number of hydrogen bond parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        for (i = 0; i < ntypes; ++i)
          for (j = 0; j < ntypes; ++j)
            for (k = 0; k < ntypes; ++k)
              reax->hbp[i][j][k].r0_hb = -1.0;

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          if (values.count() < 7)
            THROW_ERROR("Invalid force field file format");

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;

          if ((j < 0) || (k < 0) || (l < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes) && (l < ntypes)) {
            reax->hbp[j][k][l].r0_hb = values.next_double();
            reax->hbp[j][k][l].p_hb1 = values.next_double();
            reax->hbp[j][k][l].p_hb2 = values.next_double();
            reax->hbp[j][k][l].p_hb3 = values.next_double();
          }
        }

        memory->destroy(tor_flag);
      } catch (std::exception &e) {
        error->one(FLERR,e.what());
      }
    }

    // broadcast global parameters and allocate list on ranks != 0
    MPI_Bcast(&reax->gp,sizeof(global_parameters),MPI_CHAR,0,world);
    if (control->me != 0) memory->create(reax->gp.l,reax->gp.n_global,"reaxff:gp.l");
    MPI_Bcast(reax->gp.l,reax->gp.n_global,MPI_DOUBLE,0,world);

    // allocate storage for atom type based data
    MPI_Bcast(&reax->num_atom_types,1,MPI_INT,0,world);
    if (control->me != 0) {
      const int n = reax->num_atom_types;
      memory->create(reax->sbp,n,"reaxff:sbp");
      memory->create(reax->tbp,n,n,"reaxff:tbp");
      memory->create(reax->thbp,n,n,n,"reaxff:thbp");
      memory->create(reax->hbp,n,n,n,"reaxff:hbp");
      memory->create(reax->fbp,n,n,n,n,"reaxff:fbp");
    }

    // broadcast type specific force field data
    const int n = reax->num_atom_types;
    MPI_Bcast(&(reax->sbp[0]),sizeof(single_body_parameters)*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->tbp[0][0]),sizeof(two_body_parameters)*n*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->thbp[0][0][0]),sizeof(three_body_header)*n*n*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->hbp[0][0][0]),sizeof(hbond_parameters)*n*n*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->fbp[0][0][0][0]),sizeof(four_body_header)*n*n*n*n,MPI_CHAR,0,world);

    // apply parameters to various settings

    control->bo_cut    = 0.01 * reax->gp.l[29];
    control->nonb_low  = reax->gp.l[11];
    control->nonb_cut  = reax->gp.l[12];
  }
#undef THROW_ERROR
}

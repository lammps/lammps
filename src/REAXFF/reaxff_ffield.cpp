// clang-format off
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
  Algorithmic Techniques", Parallel Computing,  38, 245-259 (2012).

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

#include <cmath>
#include <cstring>
#include <exception>
#include <string>

using LAMMPS_NS::utils::open_potential;
using LAMMPS_NS::utils::getsyserror;
using LAMMPS_NS::utils::strmatch;
using LAMMPS_NS::utils::uppercase;
using LAMMPS_NS::EOFException;
using LAMMPS_NS::ValueTokenizer;

namespace ReaxFF {

  class ffield_parser_error : public std::exception {
    std::string message;
  public:
    explicit ffield_parser_error(const std::string &mesg) { message = mesg; }
    const char *what() const noexcept override { return message.c_str(); }
  };

  void Read_Force_Field(const char *filename, reax_interaction *reax,
                        control_params *control, MPI_Comm world)
  {
    char ****tor_flag;
    auto error = control->error_ptr;
    auto lmp = control->lmp_ptr;
    auto memory = control->lmp_ptr->memory;

    // read and parse the force field only on rank 0

#define THROW_ERROR(txt) throw ffield_parser_error(fmt::format("{}:{}: {}", filename, lineno, txt))

#define CHECK_COLUMNS(want)                                                          \
  if (values.count() < static_cast<std::size_t>(want))                               \
    throw ffield_parser_error(fmt::format("{}:{}: Invalid force field file format: " \
                                          " expected {} columns but found {}",       \
                                          filename, lineno, want, values.count()))

    if (control->me == 0) {
      FILE *fp = open_potential(filename, lmp, nullptr);
      if (!fp)
        error->one(FLERR,"The ReaxFF parameter file {} cannot be opened: {}",
                   filename, getsyserror());
      LAMMPS_NS::TextFileReader reader(fp, "ReaxFF parameter");
      reader.ignore_comments = false;

      try {
        int i,j,k,l,m,n,lineno = 0;

        // check if header comment line is present

        auto line = reader.next_line();
        if (strmatch(line, "^\\s*[0-9]+\\s+!.*general parameters.*"))
          THROW_ERROR("First line of ReaxFF potential file must be a comment or empty");
        ++lineno;

        // set some defaults
        auto &gp = reax->gp;

        gp.vdw_type = 0;

        // get number of global parameters

        auto values = reader.next_values(0);
        n = values.next_int();
        gp.n_global = n;
        ++lineno;

        if (n < 1)
          THROW_ERROR("Invalid number of global parameters");

        memory->destroy(gp.l);
        memory->create(gp.l,n,"reaxff:gp.l");

        // see reaxff_types.h for mapping between l[i] and the lambdas used in ff

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          gp.l[i] = values.next_double();
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
        auto &sbp = reax->sbp;
        auto &tbp = reax->tbp;
        auto &thbp = reax->thbp;
        auto &hbp = reax->hbp;
        auto &fbp = reax->fbp;

        memory->destroy(sbp);
        memory->destroy(tbp);
        memory->destroy(thbp);
        memory->destroy(hbp);
        memory->destroy(fbp);
        memory->create(sbp,n,"reaxff:sbp");
        memory->create(tbp,n,n,"reaxff:tbp");
        memory->create(thbp,n,n,n,"reaxff:thbp");
        memory->create(hbp,n,n,n,"reaxff:hbp");
        memory->create(fbp,n,n,n,n,"reaxff:fbp");
        memory->create(tor_flag,n,n,n,n,"reaxff:tor_flag");
        memset(&(sbp[0]),0,sizeof(single_body_parameters)*n);
        memset(&(tbp[0][0]),0,sizeof(two_body_parameters)*n*n);
        memset(&(thbp[0][0][0]),0,sizeof(three_body_header)*n*n*n);
        memset(&(hbp[0][0][0]),0,sizeof(hbond_parameters)*n*n*n);
        memset(&(fbp[0][0][0][0]),0,sizeof(four_body_header)*n*n*n*n);
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
          CHECK_COLUMNS(9);

          // copy element symbol in uppercase and truncate stored element symbol if necessary
          auto element = uppercase(values.next_string());
          strncpy(sbp[i].name,element.c_str(),3);
          sbp[i].name[3] = '\0';

          sbp[i].r_s        = values.next_double();
          sbp[i].valency    = values.next_double();
          sbp[i].mass       = values.next_double();
          sbp[i].r_vdw      = values.next_double();
          sbp[i].epsilon    = values.next_double();
          sbp[i].gamma      = values.next_double();
          sbp[i].r_pi       = values.next_double();
          sbp[i].valency_e  = values.next_double();
          sbp[i].nlp_opt = 0.5 * (sbp[i].valency_e-sbp[i].valency);

          // line two

          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(8);

          sbp[i].alpha      = values.next_double();
          sbp[i].gamma_w    = values.next_double();
          sbp[i].valency_boc= values.next_double();
          sbp[i].p_ovun5    = values.next_double();
          values.skip();
          sbp[i].chi        = values.next_double();
          sbp[i].eta        = 2.0*values.next_double();
          sbp[i].p_hbond = (int) values.next_double();

          // line three

          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(8);

          sbp[i].r_pi_pi    = values.next_double();
          sbp[i].p_lp2      = values.next_double();
          values.skip();
          sbp[i].b_o_131    = values.next_double();
          sbp[i].b_o_132    = values.next_double();
          sbp[i].b_o_133    = values.next_double();
          sbp[i].bcut_acks2  = values.next_double();

          // line four

          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(8);

          sbp[i].p_ovun2    = values.next_double();
          sbp[i].p_val3     = values.next_double();
          values.skip();
          sbp[i].valency_val= values.next_double();
          sbp[i].p_val5     = values.next_double();
          sbp[i].rcore2     = values.next_double();
          sbp[i].ecore2     = values.next_double();
          sbp[i].acore2     = values.next_double();

          // read line five only when lgflag != 0

          if (lgflag) {
            values = reader.next_values(0);
            ++lineno;
            CHECK_COLUMNS(2);
            sbp[i].lgcij    = values.next_double();
            sbp[i].lgre     = values.next_double();
          } else sbp[i].lgcij = sbp[i].lgre = 0.0;

          // van der Waals settings check:

          // Inner-wall?
          if ((sbp[i].rcore2 > 0.01) && (sbp[i].acore2 > 0.01)) {
            // Shielding van der Waals?
            if (sbp[i].gamma_w > 0.5) {
              if (gp.vdw_type != 0 && gp.vdw_type != 3) {
                error->warning(FLERR, "Van der Waals parameters for element {} "
                               "indicate inner wall+shielding, but earlier "
                               "atoms indicate a different van der Waals "
                               "method. This may cause division-by-zero "
                               "errors. Keeping van der Waals setting for "
                               "earlier atoms.",sbp[i].name);

              } else {
                gp.vdw_type = 3;
              }
            } else {  // No shielding van der Waals parameters present
              if ((gp.vdw_type != 0) && (gp.vdw_type != 2)) {
                error->warning(FLERR, "Van der Waals parameters for element {} "
                               "indicate inner wall withou shielding, but "
                               "earlier atoms indicate a different van der "
                               "Waals-method. This may cause division-by-"
                               "zero errors. Keeping van der Waals setting "
                               "for earlier atoms.", sbp[i].name);
              } else {
                gp.vdw_type = 2;
              }
            }
          } else { // No Inner wall parameters present
            if (sbp[i].gamma_w > 0.5) { // Shielding vdWaals
              if ((gp.vdw_type != 0) && (gp.vdw_type != 1)) {
                error->warning(FLERR, "Van der Waals parameters for element {} "
                               "indicate shielding without inner wall, but "
                               "earlier elements indicate a different van der "
                               "Waals method. This may cause division-by-zero "
                               "errors. Keeping van der Waals setting for "
                               "earlier atoms.", sbp[i].name);
              } else {
                gp.vdw_type = 1;
              }
            } else {
              error->one(FLERR, "Inconsistent van der Waals parameters: "
                         "No shielding or inner wall set for element {}", sbp[i].name);
            }
          }
        }

        /* Equate vval3 to valf for first-row elements (25/10/2004) */
        for (i = 0; i < ntypes; i++) {
          if ((sbp[i].mass < 21) &&
              (sbp[i].valency_val != sbp[i].valency_boc)) {
            error->warning(FLERR, "Changed valency_val to valency_boc for {}", sbp[i].name);
            sbp[i].valency_val = sbp[i].valency_boc;
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
          CHECK_COLUMNS(10);

          j = values.next_int() - 1;
          k = values.next_int() - 1;

          if ((j < 0) || (k < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes)) {
            tbp[j][k].De_s    = tbp[k][j].De_s    = values.next_double();
            tbp[j][k].De_p    = tbp[k][j].De_p    = values.next_double();
            tbp[j][k].De_pp   = tbp[k][j].De_pp   = values.next_double();
            tbp[j][k].p_be1   = tbp[k][j].p_be1   = values.next_double();
            tbp[j][k].p_bo5   = tbp[k][j].p_bo5   = values.next_double();
            tbp[j][k].v13cor  = tbp[k][j].v13cor  = values.next_double();
            tbp[j][k].p_bo6   = tbp[k][j].p_bo6   = values.next_double();
            tbp[j][k].p_ovun1 = tbp[k][j].p_ovun1 = values.next_double();
          }

          // second line

          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(7);

          if ((j < ntypes) && (k < ntypes)) {
            tbp[j][k].p_be2 = tbp[k][j].p_be2 = values.next_double();
            tbp[j][k].p_bo3 = tbp[k][j].p_bo3 = values.next_double();
            tbp[j][k].p_bo4 = tbp[k][j].p_bo4 = values.next_double();
            values.skip();
            tbp[j][k].p_bo1 = tbp[k][j].p_bo1 = values.next_double();
            tbp[j][k].p_bo2 = tbp[k][j].p_bo2 = values.next_double();
            // if the 8th value is missing use 0.0
            if (values.has_next())
              tbp[j][k].ovc   = tbp[k][j].ovc   = values.next_double();
            else
              tbp[j][k].ovc   = tbp[k][j].ovc   = 0.0;
          }
        }

        for (i=0; i < ntypes; ++i) {
          for (j=i; j < ntypes; ++j) {
            tbp[i][j].r_s     = tbp[j][i].r_s     = 0.5*(sbp[j].r_s + sbp[i].r_s);
            tbp[i][j].r_p     = tbp[j][i].r_p     = 0.5*(sbp[j].r_pi + sbp[i].r_pi);
            tbp[i][j].r_pp    = tbp[j][i].r_pp    = 0.5*(sbp[j].r_pi_pi + sbp[i].r_pi_pi);
            tbp[i][j].p_boc3  = tbp[j][i].p_boc3  = sqrt(sbp[j].b_o_132 * sbp[i].b_o_132);
            tbp[i][j].p_boc4  = tbp[j][i].p_boc4  = sqrt(sbp[j].b_o_131 * sbp[i].b_o_131);
            tbp[i][j].p_boc5  = tbp[j][i].p_boc5  = sqrt(sbp[j].b_o_133 * sbp[i].b_o_133);
            tbp[i][j].D       = tbp[j][i].D       = sqrt(sbp[j].epsilon * sbp[i].epsilon);
            tbp[i][j].alpha   = tbp[j][i].alpha   = sqrt(sbp[j].alpha * sbp[i].alpha);
            tbp[i][j].r_vdW   = tbp[j][i].r_vdW   = 2.0*sqrt(sbp[j].r_vdw * sbp[i].r_vdw);
            tbp[i][j].gamma_w = tbp[j][i].gamma_w = sqrt(sbp[j].gamma_w * sbp[i].gamma_w);
            tbp[i][j].gamma   = tbp[j][i].gamma   = pow(sbp[j].gamma * sbp[i].gamma,-1.5);

            // additions for additional vdWaals interaction types - inner core

            tbp[i][j].rcore   = tbp[j][i].rcore   = sqrt(sbp[i].rcore2 * sbp[j].rcore2);
            tbp[i][j].ecore   = tbp[j][i].ecore   = sqrt(sbp[i].ecore2 * sbp[j].ecore2);
            tbp[i][j].acore   = tbp[j][i].acore   = sqrt(sbp[i].acore2 * sbp[j].acore2);

            // additions for additional vdWalls interaction types lg correction

            tbp[i][j].lgcij   = tbp[j][i].lgcij   = sqrt(sbp[i].lgcij * sbp[j].lgcij);
            tbp[i][j].lgre    = tbp[j][i].lgre    = 2.0*gp.l[35]*sqrt(sbp[i].lgre*sbp[j].lgre);
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
          CHECK_COLUMNS(8 + lgflag);

          j = values.next_int() - 1;
          k = values.next_int() - 1;

          if ((j < 0) || (k < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes)) {
            val = values.next_double();
            if (val > 0.0) tbp[j][k].D = tbp[k][j].D = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_vdW = tbp[k][j].r_vdW = 2*val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].alpha = tbp[k][j].alpha = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_s = tbp[k][j].r_s = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_p = tbp[k][j].r_p = val;

            val = values.next_double();
            if (val > 0.0) tbp[j][k].r_pp = tbp[k][j].r_pp = val;

            if (lgflag) {
              val = values.next_double();
              if (val >= 0.0) tbp[j][k].lgcij = tbp[k][j].lgcij = val;
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
          CHECK_COLUMNS(10);

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;

          if ((j < 0) || (k < 0) || (l < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes) && (l < ntypes)) {

            cnt = thbp[j][k][l].cnt;
            thbp[j][k][l].cnt++;
            thbp[l][k][j].cnt++;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].theta_00 = val;
            thbp[l][k][j].prm[cnt].theta_00 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val1 = val;
            thbp[l][k][j].prm[cnt].p_val1 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val2 = val;
            thbp[l][k][j].prm[cnt].p_val2 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_coa1 = val;
            thbp[l][k][j].prm[cnt].p_coa1 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val7 = val;
            thbp[l][k][j].prm[cnt].p_val7 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_pen1 = val;
            thbp[l][k][j].prm[cnt].p_pen1 = val;

            val = values.next_double();
            thbp[j][k][l].prm[cnt].p_val4 = val;
            thbp[l][k][j].prm[cnt].p_val4 = val;
          }
        }

        // next line is number of four body parameters

        values = reader.next_values(0);
        n = values.next_int();
        ++lineno;

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(9);

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;
          m = values.next_int() - 1;

          if ((j < -1) || (k < 0) || (l < 0) || (m < -1))
            THROW_ERROR("Inconsistent force field file");

          const double val1 = values.next_double();
          const double val2 = values.next_double();
          const double val3 = values.next_double();
          const double val4 = values.next_double();
          const double val5 = values.next_double();

          if ((j >= 0) && (m >= 0)) { // this means the entry is not in compact form

            if ((j < ntypes) && (k < ntypes) && (l < ntypes) && (m < ntypes)) {

              tor_flag[j][k][l][m] = tor_flag[m][l][k][j] = 1;
              fbp[j][k][l][m].cnt  = fbp[m][l][k][j].cnt  = 1;

              fbp[j][k][l][m].prm[0].V1     = fbp[m][l][k][j].prm[0].V1     = val1;
              fbp[j][k][l][m].prm[0].V2     = fbp[m][l][k][j].prm[0].V2     = val2;
              fbp[j][k][l][m].prm[0].V3     = fbp[m][l][k][j].prm[0].V3     = val3;
              fbp[j][k][l][m].prm[0].p_tor1 = fbp[m][l][k][j].prm[0].p_tor1 = val4;
              fbp[j][k][l][m].prm[0].p_cot1 = fbp[m][l][k][j].prm[0].p_cot1 = val5;
            }

          } else { /* This means the entry is of the form 0-X-Y-0 */

            if ((k < ntypes) && (l < ntypes)) {
              for (int p = 0; p < ntypes; ++p) {
                for (int o = 0; o < ntypes; ++o) {

                  if (tor_flag[p][k][l][o] == 0) {
                    fbp[p][k][l][o].cnt = 1;
                    fbp[p][k][l][o].prm[0].V1 = val1;
                    fbp[p][k][l][o].prm[0].V2 = val2;
                    fbp[p][k][l][o].prm[0].V3 = val3;
                    fbp[p][k][l][o].prm[0].p_tor1 = val4;
                    fbp[p][k][l][o].prm[0].p_cot1 = val5;
                  }

                  if (tor_flag[o][l][k][p] == 0) {
                    fbp[o][l][k][p].cnt = 1;
                    fbp[o][l][k][p].prm[0].V1 = val1;
                    fbp[o][l][k][p].prm[0].V2 = val2;
                    fbp[o][l][k][p].prm[0].V3 = val3;
                    fbp[o][l][k][p].prm[0].p_tor1 = val4;
                    fbp[o][l][k][p].prm[0].p_cot1 = val5;
                  }
                }
              }
            }
          }
        }

        // next line is number of hydrogen bond parameters. that block may be missing

        for (i = 0; i < ntypes; ++i)
          for (j = 0; j < ntypes; ++j)
            for (k = 0; k < ntypes; ++k)
              hbp[i][j][k].r0_hb = -1.0;

        auto thisline = reader.next_line();
        if (!thisline) throw EOFException("ReaxFF parameter file has no hydrogen bond parameters");

        values = ValueTokenizer(thisline);
        n = values.next_int();
        ++lineno;

        for (i = 0; i < n; ++i) {
          values = reader.next_values(0);
          ++lineno;
          CHECK_COLUMNS(7);

          j = values.next_int() - 1;
          k = values.next_int() - 1;
          l = values.next_int() - 1;

          if ((j < 0) || (k < 0) || (l < 0))
            THROW_ERROR("Inconsistent force field file");

          if ((j < ntypes) && (k < ntypes) && (l < ntypes)) {
            hbp[j][k][l].r0_hb = values.next_double();
            hbp[j][k][l].p_hb1 = values.next_double();
            hbp[j][k][l].p_hb2 = values.next_double();
            hbp[j][k][l].p_hb3 = values.next_double();
          }
        }

        memory->destroy(tor_flag);
      } catch (EOFException &e) {
        error->warning(FLERR, e.what());
      } catch (std::exception &e) {
        error->one(FLERR,e.what());
      }
      fclose(fp);
    }

    // broadcast global parameters and allocate list on ranks != 0
    MPI_Bcast(&reax->gp,sizeof(global_parameters),MPI_CHAR,0,world);
    if (control->me != 0) memory->create(reax->gp.l,reax->gp.n_global,"reaxff:gp.l");
    MPI_Bcast(reax->gp.l,reax->gp.n_global,MPI_DOUBLE,0,world);

    // allocate storage for atom type based data
    MPI_Bcast(&reax->num_atom_types,1,MPI_INT,0,world);
    const int n = reax->num_atom_types;
    if (control->me != 0) {
      memory->create(reax->sbp,n,"reaxff:sbp");
      memory->create(reax->tbp,n,n,"reaxff:tbp");
      memory->create(reax->thbp,n,n,n,"reaxff:thbp");
      memory->create(reax->hbp,n,n,n,"reaxff:hbp");
      memory->create(reax->fbp,n,n,n,n,"reaxff:fbp");
    }

    // broadcast type specific force field data
    MPI_Bcast(&(reax->sbp[0]),sizeof(single_body_parameters)*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->tbp[0][0]),sizeof(two_body_parameters)*n*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->thbp[0][0][0]),sizeof(three_body_header)*n*n*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->hbp[0][0][0]),sizeof(hbond_parameters)*n*n*n,MPI_CHAR,0,world);
    MPI_Bcast(&(reax->fbp[0][0][0][0]),sizeof(four_body_header)*n*n*n*n,MPI_CHAR,0,world);

    // apply global parameters to global control settings

    control->bo_cut    = 0.01 * reax->gp.l[29];
    control->nonb_low  = reax->gp.l[11];
    control->nonb_cut  = reax->gp.l[12];
  }
#undef THROW_ERROR
#undef CHECK_COLUMNS
}

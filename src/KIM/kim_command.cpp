/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U),
                         Yaser Afshar (UMN)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-2.1.0 (and newer) package
------------------------------------------------------------------------- */

#include "kim_command.h"

#include "citeme.h"
#include "error.h"

// include KIM sub-command headers here
#include "kim_init.h"
#include "kim_interactions.h"
#include "kim_param.h"
#include "kim_property.h"
#include "kim_query.h"

using namespace LAMMPS_NS;

static constexpr const char *const cite_openkim =
    "OpenKIM Project: doi:10.1007/s11837-011-0102-6\n\n"
    "@Article{tadmor:elliott:2011,\n"
    " author = {E. B. Tadmor and R. S. Elliott and J. P. Sethna and R. E. Miller "
    "and C. A. Becker},\n"
    " title = {The potential of atomistic simulations and the {K}nowledgebase of "
    "{I}nteratomic {M}odels},\n"
    " journal = {{JOM}},\n"
    " year =    2011,\n"
    " volume =  63,\n"
    " number =  17,\n"
    " pages =   {17},\n"
    " doi =     {10.1007/s11837-011-0102-6}\n"
    "}\n\n";

static constexpr const char *const cite_openkim_query =
    "OpenKIM query: doi:10.1063/5.0014267\n\n"
    "@Article{karls:bierbaum:2020,\n"
    " author = {D. S. Karls and M. Bierbaum and A. A. Alemi and R. S. Elliott "
    "and J. P. Sethna and E. B. Tadmor},\n"
    " title = {The {O}pen{KIM} processing pipeline: {A} cloud-based automatic "
    "material property computation engine},\n"
    " journal = {{T}he {J}ournal of {C}hemical {P}hysics},\n"
    " year =    2020,\n"
    " volume =  153,\n"
    " number =  6,\n"
    " pages =   {064104},\n"
    " doi =     {10.1063/5.0014267}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

void KimCommand::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal kim command");

  const std::string subcmd(arg[0]);
  narg--;
  arg++;

  if (lmp->citeme) lmp->citeme->add(cite_openkim);

  if (subcmd == "init") {
    auto cmd = new KimInit(lmp);
    cmd->command(narg, arg);
    delete cmd;
  } else if (subcmd == "interactions") {
    auto cmd = new KimInteractions(lmp);
    cmd->command(narg, arg);
    delete cmd;
  } else if (subcmd == "param") {
    auto cmd = new KimParam(lmp);
    cmd->command(narg, arg);
    delete cmd;
  } else if (subcmd == "property") {
    auto cmd = new KimProperty(lmp);
    cmd->command(narg, arg);
    delete cmd;
  } else if (subcmd == "query") {
    if (lmp->citeme) lmp->citeme->add(cite_openkim_query);
    auto cmd = new KimQuery(lmp);
    cmd->command(narg, arg);
    delete cmd;
  } else
    error->all(FLERR, "Unknown kim subcommand {}", subcmd);
}

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
   Registration table for the granular sub-models of the GRANULAR package.

   This is the list of all granular sub-models that the contact models of
   "pair_style granular" and "fix wall/gran" can select by keyword.  Unlike
   the name-keyed styles (pair, fix, compute, ...), granular sub-models cannot
   be loaded as plugins, so this list is a plain checked-in table instead of an
   auto-generated style header.  It replaces the "style_gran_sub_mod.h" file
   that the build system used to generate.

   To add a new granular sub-model, only this file needs to be edited:

     (1) write the sub-model class in a gran_sub_mod_<group>.{h,cpp} file pair,
         deriving it from the matching base class (e.g. GranSubModNormal),
     (2) add an "#include" for its header in the 1st customization section
         below, and
     (3) add one row to the gran_sub_mod_table[] array in the 2nd customization
         section below: the keyword used in the input script, the class name,
         and the sub-model type (NORMAL, DAMPING, TANGENTIAL, ROLLING, TWISTING,
         or HEAT).

   The same keyword may appear for different sub-model types (e.g. "none" exists
   for every type); the type column disambiguates them.
------------------------------------------------------------------------- */

#include "gran_sub_mod.h"
#include "granular_model.h"

// 1st customization section: add the #include for a new sub-model header here

#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_heat.h"
#include "gran_sub_mod_normal.h"
#include "gran_sub_mod_rolling.h"
#include "gran_sub_mod_tangential.h"
#include "gran_sub_mod_twisting.h"

namespace LAMMPS_NS::Granular_NS {

// factory function: one instance is created for each entry in the table below

template <typename T> static GranSubMod *creator(GranularModel *gm, LAMMPS *lmp)
{
  return new T(gm, lmp);
}

// clang-format off

// 2nd customization section: add a row for a new sub-model here.  The columns
// are the input-script keyword, the sub-model class name, and the sub-model
// type.  The order of the rows does not matter.

LMP_REGISTRY_CONST GranSubModInfo gran_sub_mod_table[] = {

  // normal models
  { "none",                   &creator<GranSubModNormalNone>,                   NORMAL },
  { "hooke",                  &creator<GranSubModNormalHooke>,                  NORMAL },
  { "hertz",                  &creator<GranSubModNormalHertz>,                  NORMAL },
  { "hertz/material",         &creator<GranSubModNormalHertzMaterial>,          NORMAL },
  { "dmt",                    &creator<GranSubModNormalDMT>,                    NORMAL },
  { "jkr",                    &creator<GranSubModNormalJKR>,                    NORMAL },
  { "mdr",                    &creator<GranSubModNormalMDR>,                    NORMAL },

  // damping models
  { "none",                   &creator<GranSubModDampingNone>,                  DAMPING },
  { "velocity",               &creator<GranSubModDampingVelocity>,              DAMPING },
  { "mass_velocity",          &creator<GranSubModDampingMassVelocity>,          DAMPING },
  { "viscoelastic",           &creator<GranSubModDampingViscoelastic>,          DAMPING },
  { "tsuji",                  &creator<GranSubModDampingTsuji>,                  DAMPING },
  { "coeff_restitution",      &creator<GranSubModDampingCoeffRestitution>,      DAMPING },
  { "mdr",                    &creator<GranSubModDampingMDR>,                    DAMPING },

  // tangential models
  { "none",                   &creator<GranSubModTangentialNone>,               TANGENTIAL },
  { "linear_nohistory",       &creator<GranSubModTangentialLinearNoHistory>,    TANGENTIAL },
  { "linear_history",         &creator<GranSubModTangentialLinearHistory>,      TANGENTIAL },
  { "linear_history_classic", &creator<GranSubModTangentialLinearHistoryClassic>, TANGENTIAL },
  { "mindlin_classic",        &creator<GranSubModTangentialMindlinClassic>,     TANGENTIAL },
  { "mindlin",                &creator<GranSubModTangentialMindlin>,            TANGENTIAL },
  { "mindlin/force",          &creator<GranSubModTangentialMindlinForce>,       TANGENTIAL },
  { "mindlin_rescale",        &creator<GranSubModTangentialMindlinRescale>,     TANGENTIAL },
  { "mindlin_rescale/force",  &creator<GranSubModTangentialMindlinRescaleForce>, TANGENTIAL },

  // rolling models
  { "none",                   &creator<GranSubModRollingNone>,                  ROLLING },
  { "sds",                    &creator<GranSubModRollingSDS>,                   ROLLING },

  // twisting models
  { "none",                   &creator<GranSubModTwistingNone>,                 TWISTING },
  { "marshall",               &creator<GranSubModTwistingMarshall>,             TWISTING },
  { "sds",                    &creator<GranSubModTwistingSDS>,                  TWISTING },

  // heat conduction models
  { "none",                   &creator<GranSubModHeatNone>,                     HEAT },
  { "radius",                 &creator<GranSubModHeatRadius>,                   HEAT },
  { "area",                   &creator<GranSubModHeatArea>,                     HEAT }
};
// clang-format on

LMP_REGISTRY_CONST int num_gran_sub_mod = sizeof(gran_sub_mod_table) / sizeof(GranSubModInfo);

}    // namespace LAMMPS_NS::Granular_NS

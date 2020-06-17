/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
#include "input.h"
#include "lammps.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::Eq;

class AtomStyleTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"SimpleCommandsTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }
};

TEST_F(AtomStyleTest, atomic)
{
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->atom->nlocal, 0);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_EQ(lmp->atom->nmax, 1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->nellipsoids, 0);
    ASSERT_EQ(lmp->atom->nlines, 0);
    ASSERT_EQ(lmp->atom->ntris, 0);
    ASSERT_EQ(lmp->atom->nbodies, 0);
    ASSERT_EQ(lmp->atom->nbonds, 0);
    ASSERT_EQ(lmp->atom->nangles, 0);
    ASSERT_EQ(lmp->atom->ndihedrals, 0);
    ASSERT_EQ(lmp->atom->nimpropers, 0);
    ASSERT_EQ(lmp->atom->ntypes, 0);
    ASSERT_EQ(lmp->atom->nbondtypes, 0);
    ASSERT_EQ(lmp->atom->nangletypes, 0);
    ASSERT_EQ(lmp->atom->ndihedraltypes, 0);
    ASSERT_EQ(lmp->atom->nimpropertypes, 0);
    ASSERT_EQ(lmp->atom->bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->improper_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_improper_per_atom, 0);

    ASSERT_NE(lmp->atom->tag, nullptr);
    ASSERT_NE(lmp->atom->type, nullptr);
    ASSERT_NE(lmp->atom->mask, nullptr);
    ASSERT_NE(lmp->atom->image, nullptr);
    ASSERT_NE(lmp->atom->x, nullptr);
    ASSERT_NE(lmp->atom->v, nullptr);
    ASSERT_NE(lmp->atom->f, nullptr);
    ASSERT_EQ(lmp->atom->q, nullptr);
    ASSERT_EQ(lmp->atom->mu, nullptr);
    ASSERT_EQ(lmp->atom->omega, nullptr);
    ASSERT_EQ(lmp->atom->angmom, nullptr);
    ASSERT_EQ(lmp->atom->torque, nullptr);
    ASSERT_EQ(lmp->atom->radius, nullptr);
    ASSERT_EQ(lmp->atom->rmass, nullptr);
    ASSERT_EQ(lmp->atom->ellipsoid, nullptr);
    ASSERT_EQ(lmp->atom->line, nullptr);
    ASSERT_EQ(lmp->atom->tri, nullptr);
    ASSERT_EQ(lmp->atom->body, nullptr);
    ASSERT_EQ(lmp->atom->molecule, nullptr);
    ASSERT_EQ(lmp->atom->molindex, nullptr);
    ASSERT_EQ(lmp->atom->molatom, nullptr);
    ASSERT_EQ(lmp->atom->num_bond, nullptr);
    ASSERT_EQ(lmp->atom->bond_type, nullptr);
    ASSERT_EQ(lmp->atom->bond_atom, nullptr);
    ASSERT_EQ(lmp->atom->num_angle, nullptr);
    ASSERT_EQ(lmp->atom->angle_type, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom1, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom2, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom3, nullptr);
    ASSERT_EQ(lmp->atom->num_dihedral, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_type, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom1, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom2, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom3, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom4, nullptr);
    ASSERT_EQ(lmp->atom->num_improper, nullptr);
    ASSERT_EQ(lmp->atom->improper_type, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom1, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom2, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom3, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom4, nullptr);
    ASSERT_EQ(lmp->atom->maxspecial, 1);
    ASSERT_EQ(lmp->atom->nspecial, nullptr);
    ASSERT_EQ(lmp->atom->special, nullptr);
    ASSERT_EQ(lmp->atom->vfrac, nullptr);
    ASSERT_EQ(lmp->atom->s0, nullptr);
    ASSERT_EQ(lmp->atom->x0, nullptr);
    ASSERT_EQ(lmp->atom->sp, nullptr);
    ASSERT_EQ(lmp->atom->fm, nullptr);
    ASSERT_EQ(lmp->atom->fm_long, nullptr);
    ASSERT_EQ(lmp->atom->spin, nullptr);
    ASSERT_EQ(lmp->atom->eradius, nullptr);
    ASSERT_EQ(lmp->atom->ervel, nullptr);
    ASSERT_EQ(lmp->atom->erforce, nullptr);
    ASSERT_EQ(lmp->atom->ervelforce, nullptr);
    ASSERT_EQ(lmp->atom->cs, nullptr);
    ASSERT_EQ(lmp->atom->csforce, nullptr);
    ASSERT_EQ(lmp->atom->vforce, nullptr);
    ASSERT_EQ(lmp->atom->etag, nullptr);
    ASSERT_EQ(lmp->atom->uCond, nullptr);
    ASSERT_EQ(lmp->atom->uMech, nullptr);
    ASSERT_EQ(lmp->atom->uChem, nullptr);
    ASSERT_EQ(lmp->atom->uCG, nullptr);
    ASSERT_EQ(lmp->atom->uCGnew, nullptr);
    ASSERT_EQ(lmp->atom->duChem, nullptr);
    ASSERT_EQ(lmp->atom->dpdTheta, nullptr);
    ASSERT_EQ(lmp->atom->cc, nullptr);
    ASSERT_EQ(lmp->atom->cc_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_temp, nullptr);
    ASSERT_EQ(lmp->atom->edpd_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_cv, nullptr);
    ASSERT_EQ(lmp->atom->length, nullptr);
    ASSERT_EQ(lmp->atom->buckling, nullptr);
    ASSERT_EQ(lmp->atom->bond_nt, nullptr);
    ASSERT_EQ(lmp->atom->contact_radius, nullptr);
    ASSERT_EQ(lmp->atom->smd_data_9, nullptr);
    ASSERT_EQ(lmp->atom->smd_stress, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_rate, nullptr);
    ASSERT_EQ(lmp->atom->damage, nullptr);
    ASSERT_EQ(lmp->atom->rho, nullptr);
    ASSERT_EQ(lmp->atom->drho, nullptr);
    ASSERT_EQ(lmp->atom->esph, nullptr);
    ASSERT_EQ(lmp->atom->desph, nullptr);
    ASSERT_EQ(lmp->atom->cv, nullptr);
    ASSERT_EQ(lmp->atom->vest, nullptr);
    ASSERT_EQ(lmp->atom->nmolecule, 0);
    ASSERT_EQ(lmp->atom->molecules, nullptr);
    ASSERT_EQ(lmp->atom->nivector, 0);
    ASSERT_EQ(lmp->atom->ndvector, 0);
    ASSERT_EQ(lmp->atom->iname, nullptr);
    ASSERT_EQ(lmp->atom->dname, nullptr);
    ASSERT_EQ(lmp->atom->mass, nullptr);
    ASSERT_EQ(lmp->atom->mass_setflag, nullptr);
    ASSERT_EQ(lmp->atom->nextra_grow, 0);
    ASSERT_EQ(lmp->atom->nextra_restart, 0);
    ASSERT_EQ(lmp->atom->nextra_border, 0);
    ASSERT_EQ(lmp->atom->nextra_grow_max, 0);
    ASSERT_EQ(lmp->atom->nextra_restart_max, 0);
    ASSERT_EQ(lmp->atom->nextra_border_max, 0);
    ASSERT_EQ(lmp->atom->nextra_store, 0);
    ASSERT_EQ(lmp->atom->extra_grow, nullptr);
    ASSERT_EQ(lmp->atom->extra_restart, nullptr);
    ASSERT_EQ(lmp->atom->extra_border, nullptr);
    ASSERT_EQ(lmp->atom->extra, nullptr);
    ASSERT_EQ(lmp->atom->sametag, nullptr);
    ASSERT_EQ(lmp->atom->map_style, 0);
    ASSERT_EQ(lmp->atom->map_user, 0);
    ASSERT_EQ(lmp->atom->map_tag_max, -1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("atom_style charge");
    lmp->input->one("atom_style atomic");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->atom->nlocal, 0);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_EQ(lmp->atom->nmax, 1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->nellipsoids, 0);
    ASSERT_EQ(lmp->atom->nlines, 0);
    ASSERT_EQ(lmp->atom->ntris, 0);
    ASSERT_EQ(lmp->atom->nbodies, 0);
    ASSERT_EQ(lmp->atom->nbonds, 0);
    ASSERT_EQ(lmp->atom->nangles, 0);
    ASSERT_EQ(lmp->atom->ndihedrals, 0);
    ASSERT_EQ(lmp->atom->nimpropers, 0);
    ASSERT_EQ(lmp->atom->ntypes, 0);
    ASSERT_EQ(lmp->atom->nbondtypes, 0);
    ASSERT_EQ(lmp->atom->nangletypes, 0);
    ASSERT_EQ(lmp->atom->ndihedraltypes, 0);
    ASSERT_EQ(lmp->atom->nimpropertypes, 0);
    ASSERT_EQ(lmp->atom->bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->improper_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_improper_per_atom, 0);

    ASSERT_NE(lmp->atom->tag, nullptr);
    ASSERT_NE(lmp->atom->type, nullptr);
    ASSERT_NE(lmp->atom->mask, nullptr);
    ASSERT_NE(lmp->atom->image, nullptr);
    ASSERT_NE(lmp->atom->x, nullptr);
    ASSERT_NE(lmp->atom->v, nullptr);
    ASSERT_NE(lmp->atom->f, nullptr);
    ASSERT_EQ(lmp->atom->q, nullptr);
    ASSERT_EQ(lmp->atom->mu, nullptr);
    ASSERT_EQ(lmp->atom->omega, nullptr);
    ASSERT_EQ(lmp->atom->angmom, nullptr);
    ASSERT_EQ(lmp->atom->torque, nullptr);
    ASSERT_EQ(lmp->atom->radius, nullptr);
    ASSERT_EQ(lmp->atom->rmass, nullptr);
    ASSERT_EQ(lmp->atom->ellipsoid, nullptr);
    ASSERT_EQ(lmp->atom->line, nullptr);
    ASSERT_EQ(lmp->atom->tri, nullptr);
    ASSERT_EQ(lmp->atom->body, nullptr);
    ASSERT_EQ(lmp->atom->molecule, nullptr);
    ASSERT_EQ(lmp->atom->molindex, nullptr);
    ASSERT_EQ(lmp->atom->molatom, nullptr);
    ASSERT_EQ(lmp->atom->num_bond, nullptr);
    ASSERT_EQ(lmp->atom->bond_type, nullptr);
    ASSERT_EQ(lmp->atom->bond_atom, nullptr);
    ASSERT_EQ(lmp->atom->num_angle, nullptr);
    ASSERT_EQ(lmp->atom->angle_type, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom1, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom2, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom3, nullptr);
    ASSERT_EQ(lmp->atom->num_dihedral, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_type, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom1, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom2, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom3, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom4, nullptr);
    ASSERT_EQ(lmp->atom->num_improper, nullptr);
    ASSERT_EQ(lmp->atom->improper_type, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom1, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom2, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom3, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom4, nullptr);
    ASSERT_EQ(lmp->atom->maxspecial, 1);
    ASSERT_EQ(lmp->atom->nspecial, nullptr);
    ASSERT_EQ(lmp->atom->special, nullptr);
    ASSERT_EQ(lmp->atom->vfrac, nullptr);
    ASSERT_EQ(lmp->atom->s0, nullptr);
    ASSERT_EQ(lmp->atom->x0, nullptr);
    ASSERT_EQ(lmp->atom->sp, nullptr);
    ASSERT_EQ(lmp->atom->fm, nullptr);
    ASSERT_EQ(lmp->atom->fm_long, nullptr);
    ASSERT_EQ(lmp->atom->spin, nullptr);
    ASSERT_EQ(lmp->atom->eradius, nullptr);
    ASSERT_EQ(lmp->atom->ervel, nullptr);
    ASSERT_EQ(lmp->atom->erforce, nullptr);
    ASSERT_EQ(lmp->atom->ervelforce, nullptr);
    ASSERT_EQ(lmp->atom->cs, nullptr);
    ASSERT_EQ(lmp->atom->csforce, nullptr);
    ASSERT_EQ(lmp->atom->vforce, nullptr);
    ASSERT_EQ(lmp->atom->etag, nullptr);
    ASSERT_EQ(lmp->atom->uCond, nullptr);
    ASSERT_EQ(lmp->atom->uMech, nullptr);
    ASSERT_EQ(lmp->atom->uChem, nullptr);
    ASSERT_EQ(lmp->atom->uCG, nullptr);
    ASSERT_EQ(lmp->atom->uCGnew, nullptr);
    ASSERT_EQ(lmp->atom->duChem, nullptr);
    ASSERT_EQ(lmp->atom->dpdTheta, nullptr);
    ASSERT_EQ(lmp->atom->cc, nullptr);
    ASSERT_EQ(lmp->atom->cc_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_temp, nullptr);
    ASSERT_EQ(lmp->atom->edpd_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_cv, nullptr);
    ASSERT_EQ(lmp->atom->length, nullptr);
    ASSERT_EQ(lmp->atom->buckling, nullptr);
    ASSERT_EQ(lmp->atom->bond_nt, nullptr);
    ASSERT_EQ(lmp->atom->contact_radius, nullptr);
    ASSERT_EQ(lmp->atom->smd_data_9, nullptr);
    ASSERT_EQ(lmp->atom->smd_stress, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_rate, nullptr);
    ASSERT_EQ(lmp->atom->damage, nullptr);
    ASSERT_EQ(lmp->atom->rho, nullptr);
    ASSERT_EQ(lmp->atom->drho, nullptr);
    ASSERT_EQ(lmp->atom->esph, nullptr);
    ASSERT_EQ(lmp->atom->desph, nullptr);
    ASSERT_EQ(lmp->atom->cv, nullptr);
    ASSERT_EQ(lmp->atom->vest, nullptr);
    ASSERT_EQ(lmp->atom->nmolecule, 0);
    ASSERT_EQ(lmp->atom->molecules, nullptr);
    ASSERT_EQ(lmp->atom->nivector, 0);
    ASSERT_EQ(lmp->atom->ndvector, 0);
    ASSERT_EQ(lmp->atom->iname, nullptr);
    ASSERT_EQ(lmp->atom->dname, nullptr);
    ASSERT_EQ(lmp->atom->mass, nullptr);
    ASSERT_EQ(lmp->atom->mass_setflag, nullptr);
    ASSERT_EQ(lmp->atom->nextra_grow, 0);
    ASSERT_EQ(lmp->atom->nextra_restart, 0);
    ASSERT_EQ(lmp->atom->nextra_border, 0);
    ASSERT_EQ(lmp->atom->nextra_grow_max, 0);
    ASSERT_EQ(lmp->atom->nextra_restart_max, 0);
    ASSERT_EQ(lmp->atom->nextra_border_max, 0);
    ASSERT_EQ(lmp->atom->nextra_store, 0);
    ASSERT_EQ(lmp->atom->extra_grow, nullptr);
    ASSERT_EQ(lmp->atom->extra_restart, nullptr);
    ASSERT_EQ(lmp->atom->extra_border, nullptr);
    ASSERT_EQ(lmp->atom->extra, nullptr);
    ASSERT_EQ(lmp->atom->sametag, nullptr);
    ASSERT_EQ(lmp->atom->map_style, 0);
    ASSERT_EQ(lmp->atom->map_user, 0);
    ASSERT_EQ(lmp->atom->map_tag_max, -1);
}

TEST_F(AtomStyleTest, charge)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("atom_style charge");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("charge"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->atom->nlocal, 0);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_EQ(lmp->atom->nmax, 1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->nellipsoids, 0);
    ASSERT_EQ(lmp->atom->nlines, 0);
    ASSERT_EQ(lmp->atom->ntris, 0);
    ASSERT_EQ(lmp->atom->nbodies, 0);
    ASSERT_EQ(lmp->atom->nbonds, 0);
    ASSERT_EQ(lmp->atom->nangles, 0);
    ASSERT_EQ(lmp->atom->ndihedrals, 0);
    ASSERT_EQ(lmp->atom->nimpropers, 0);
    ASSERT_EQ(lmp->atom->ntypes, 0);
    ASSERT_EQ(lmp->atom->nbondtypes, 0);
    ASSERT_EQ(lmp->atom->nangletypes, 0);
    ASSERT_EQ(lmp->atom->ndihedraltypes, 0);
    ASSERT_EQ(lmp->atom->nimpropertypes, 0);
    ASSERT_EQ(lmp->atom->bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->improper_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_improper_per_atom, 0);

    ASSERT_NE(lmp->atom->tag, nullptr);
    ASSERT_NE(lmp->atom->type, nullptr);
    ASSERT_NE(lmp->atom->mask, nullptr);
    ASSERT_NE(lmp->atom->image, nullptr);
    ASSERT_NE(lmp->atom->x, nullptr);
    ASSERT_NE(lmp->atom->v, nullptr);
    ASSERT_NE(lmp->atom->f, nullptr);
    ASSERT_NE(lmp->atom->q, nullptr);
    ASSERT_EQ(lmp->atom->mu, nullptr);
    ASSERT_EQ(lmp->atom->omega, nullptr);
    ASSERT_EQ(lmp->atom->angmom, nullptr);
    ASSERT_EQ(lmp->atom->torque, nullptr);
    ASSERT_EQ(lmp->atom->radius, nullptr);
    ASSERT_EQ(lmp->atom->rmass, nullptr);
    ASSERT_EQ(lmp->atom->ellipsoid, nullptr);
    ASSERT_EQ(lmp->atom->line, nullptr);
    ASSERT_EQ(lmp->atom->tri, nullptr);
    ASSERT_EQ(lmp->atom->body, nullptr);
    ASSERT_EQ(lmp->atom->molecule, nullptr);
    ASSERT_EQ(lmp->atom->molindex, nullptr);
    ASSERT_EQ(lmp->atom->molatom, nullptr);
    ASSERT_EQ(lmp->atom->num_bond, nullptr);
    ASSERT_EQ(lmp->atom->bond_type, nullptr);
    ASSERT_EQ(lmp->atom->bond_atom, nullptr);
    ASSERT_EQ(lmp->atom->num_angle, nullptr);
    ASSERT_EQ(lmp->atom->angle_type, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom1, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom2, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom3, nullptr);
    ASSERT_EQ(lmp->atom->num_dihedral, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_type, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom1, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom2, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom3, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom4, nullptr);
    ASSERT_EQ(lmp->atom->num_improper, nullptr);
    ASSERT_EQ(lmp->atom->improper_type, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom1, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom2, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom3, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom4, nullptr);
    ASSERT_EQ(lmp->atom->maxspecial, 1);
    ASSERT_EQ(lmp->atom->nspecial, nullptr);
    ASSERT_EQ(lmp->atom->special, nullptr);
    ASSERT_EQ(lmp->atom->vfrac, nullptr);
    ASSERT_EQ(lmp->atom->s0, nullptr);
    ASSERT_EQ(lmp->atom->x0, nullptr);
    ASSERT_EQ(lmp->atom->sp, nullptr);
    ASSERT_EQ(lmp->atom->fm, nullptr);
    ASSERT_EQ(lmp->atom->fm_long, nullptr);
    ASSERT_EQ(lmp->atom->spin, nullptr);
    ASSERT_EQ(lmp->atom->eradius, nullptr);
    ASSERT_EQ(lmp->atom->ervel, nullptr);
    ASSERT_EQ(lmp->atom->erforce, nullptr);
    ASSERT_EQ(lmp->atom->ervelforce, nullptr);
    ASSERT_EQ(lmp->atom->cs, nullptr);
    ASSERT_EQ(lmp->atom->csforce, nullptr);
    ASSERT_EQ(lmp->atom->vforce, nullptr);
    ASSERT_EQ(lmp->atom->etag, nullptr);
    ASSERT_EQ(lmp->atom->uCond, nullptr);
    ASSERT_EQ(lmp->atom->uMech, nullptr);
    ASSERT_EQ(lmp->atom->uChem, nullptr);
    ASSERT_EQ(lmp->atom->uCG, nullptr);
    ASSERT_EQ(lmp->atom->uCGnew, nullptr);
    ASSERT_EQ(lmp->atom->duChem, nullptr);
    ASSERT_EQ(lmp->atom->dpdTheta, nullptr);
    ASSERT_EQ(lmp->atom->cc, nullptr);
    ASSERT_EQ(lmp->atom->cc_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_temp, nullptr);
    ASSERT_EQ(lmp->atom->edpd_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_cv, nullptr);
    ASSERT_EQ(lmp->atom->length, nullptr);
    ASSERT_EQ(lmp->atom->buckling, nullptr);
    ASSERT_EQ(lmp->atom->bond_nt, nullptr);
    ASSERT_EQ(lmp->atom->contact_radius, nullptr);
    ASSERT_EQ(lmp->atom->smd_data_9, nullptr);
    ASSERT_EQ(lmp->atom->smd_stress, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_rate, nullptr);
    ASSERT_EQ(lmp->atom->damage, nullptr);
    ASSERT_EQ(lmp->atom->rho, nullptr);
    ASSERT_EQ(lmp->atom->drho, nullptr);
    ASSERT_EQ(lmp->atom->esph, nullptr);
    ASSERT_EQ(lmp->atom->desph, nullptr);
    ASSERT_EQ(lmp->atom->cv, nullptr);
    ASSERT_EQ(lmp->atom->vest, nullptr);
    ASSERT_EQ(lmp->atom->nmolecule, 0);
    ASSERT_EQ(lmp->atom->molecules, nullptr);
    ASSERT_EQ(lmp->atom->nivector, 0);
    ASSERT_EQ(lmp->atom->ndvector, 0);
    ASSERT_EQ(lmp->atom->iname, nullptr);
    ASSERT_EQ(lmp->atom->dname, nullptr);
    ASSERT_EQ(lmp->atom->mass, nullptr);
    ASSERT_EQ(lmp->atom->mass_setflag, nullptr);
    ASSERT_EQ(lmp->atom->nextra_grow, 0);
    ASSERT_EQ(lmp->atom->nextra_restart, 0);
    ASSERT_EQ(lmp->atom->nextra_border, 0);
    ASSERT_EQ(lmp->atom->nextra_grow_max, 0);
    ASSERT_EQ(lmp->atom->nextra_restart_max, 0);
    ASSERT_EQ(lmp->atom->nextra_border_max, 0);
    ASSERT_EQ(lmp->atom->nextra_store, 0);
    ASSERT_EQ(lmp->atom->extra_grow, nullptr);
    ASSERT_EQ(lmp->atom->extra_restart, nullptr);
    ASSERT_EQ(lmp->atom->extra_border, nullptr);
    ASSERT_EQ(lmp->atom->extra, nullptr);
    ASSERT_EQ(lmp->atom->sametag, nullptr);
    ASSERT_EQ(lmp->atom->map_style, 0);
    ASSERT_EQ(lmp->atom->map_user, 0);
    ASSERT_EQ(lmp->atom->map_tag_max, -1);
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}

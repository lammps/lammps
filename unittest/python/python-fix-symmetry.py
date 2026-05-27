"""
Integration test for fix symmetry: drives an 8-atom P-1 inversion-symmetric
LJ system from the Python wrapper, runs a few NVE steps, and verifies that
each (asym, image) pair stays related by lattice inversion to ~FP precision.
"""

import math
import os
import sys
import tempfile
import unittest

from lammps import lammps, LAMMPS_DOUBLE_2D


def _has_symmetry():
    try:
        machine = os.environ.get("LAMMPS_MACHINE_NAME")
        lmp = lammps(name=machine)
        ok = lmp.has_style("fix", "symmetry")
        lmp.close()
        return ok
    except Exception:
        return False


HAS_SYMMETRY = _has_symmetry()


@unittest.skipUnless(HAS_SYMMETRY, "SYMMETRY package not built into this LAMMPS")
class FixSymmetryPython(unittest.TestCase):

    L = 8.0

    @classmethod
    def _write_json(cls, content):
        f = tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False, dir=os.getcwd()
        )
        f.write(content)
        f.close()
        return f.name

    def setUp(self):
        machine = os.environ.get("LAMMPS_MACHINE_NAME")
        self.lmp = lammps(
            name=machine,
            cmdargs=["-nocite", "-log", "none", "-echo", "screen"],
        )
        self._tmpfiles = []

    def tearDown(self):
        self.lmp.close()
        for p in self._tmpfiles:
            try:
                os.unlink(p)
            except OSError:
                pass

    def _set_up_pinv_system(self):
        L = self.L
        cmds = [
            "units lj",
            "atom_style atomic",
            "boundary p p p",
            "atom_modify map array",
            f"region box block 0 {L} 0 {L} 0 {L}",
            "create_box 1 box",
            "create_atoms 1 single 1.2 1.5 1.8",
            "create_atoms 1 single 3.0 1.2 0.9",
            "create_atoms 1 single 1.5 3.0 1.1",
            "create_atoms 1 single 0.9 1.1 3.2",
            f"create_atoms 1 single {L-1.2} {L-1.5} {L-1.8}",
            f"create_atoms 1 single {L-3.0} {L-1.2} {L-0.9}",
            f"create_atoms 1 single {L-1.5} {L-3.0} {L-1.1}",
            f"create_atoms 1 single {L-0.9} {L-1.1} {L-3.2}",
            "mass 1 1.0",
            "pair_style lj/cut 2.5",
            "pair_coeff * * 1.0 1.0",
            "velocity all create 0.7 12345 mom yes rot yes",
            "fix 1 all nve",
        ]
        for c in cmds:
            self.lmp.command(c)

    def _gather_positions(self):
        """Return positions as a dict: {tag (1-based int): (x, y, z)}.

        Uses extract_atom() rather than gather_atoms() so the test works
        under -DLAMMPS_BIGBIG (where the C gather_atoms function aborts
        because its tag buffer is sized for 4-byte ints). Tests are run
        serial, so nlocal == natoms and one pass covers every atom.
        """
        n = self.lmp.get_natoms()
        tag = self.lmp.extract_atom("id")
        x = self.lmp.extract_atom("x", LAMMPS_DOUBLE_2D)
        return {tag[i]: (x[i][0], x[i][1], x[i][2]) for i in range(n)}

    def _gather_forces(self):
        """Return forces as a dict: {tag (1-based int): (fx, fy, fz)}."""
        n = self.lmp.get_natoms()
        tag = self.lmp.extract_atom("id")
        f = self.lmp.extract_atom("f", LAMMPS_DOUBLE_2D)
        return {tag[i]: (f[i][0], f[i][1], f[i][2]) for i in range(n)}

    # --------------------------------------------------------------------

    def test_pinv_preserves_inversion_symmetry(self):
        """After 50 NVE steps, asym + image == L (mod L, componentwise)."""
        symfile = self._write_json(
            r"""{
                "group_name": "P-1",
                "lattice":    "triclinic",
                "ops": [
                    { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
                    { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
                ],
                "orbits": [
                    { "tags": [1, 5] }, { "tags": [2, 6] },
                    { "tags": [3, 7] }, { "tags": [4, 8] }
                ]
            }"""
        )
        self._tmpfiles.append(symfile)

        self._set_up_pinv_system()
        self.lmp.command(f"fix sym all symmetry {symfile}")
        self.lmp.command("run 50")

        pos = self._gather_positions()
        for asym in (1, 2, 3, 4):
            image = asym + 4
            xa = pos[asym]
            xb = pos[image]
            for k in range(3):
                # sum/L is some integer N (the lattice translation chosen at
                # init). Subtract its nearest integer to expose any drift.
                s = (xa[k] + xb[k]) / self.L
                residual = s - round(s)
                self.assertAlmostEqual(
                    residual,
                    0.0,
                    delta=1.0e-12,
                    msg=(
                        f"orbit {asym} component {k}: "
                        f"xa={xa[k]:.15g} xb={xb[k]:.15g} sum/L={s:.15g}"
                    ),
                )

    def test_force_is_exactly_antisymmetric(self):
        """In post_force the forces on orbit partners must be exactly opposite."""
        symfile = self._write_json(
            r"""{
                "group_name": "P-1", "lattice": "triclinic",
                "ops": [
                    { "R": [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]], "t": [0,0,0] },
                    { "R": [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]], "t": [0,0,0] }
                ],
                "orbits": [
                    { "tags": [1, 5] }, { "tags": [2, 6] },
                    { "tags": [3, 7] }, { "tags": [4, 8] }
                ]
            }"""
        )
        self._tmpfiles.append(symfile)

        self._set_up_pinv_system()
        self.lmp.command(f"fix sym all symmetry {symfile}")
        # run 0 still executes setup() which calls post_force/end_of_step
        self.lmp.command("run 0")

        f = self._gather_forces()
        for asym in (1, 2, 3, 4):
            image = asym + 4
            for k in range(3):
                self.assertAlmostEqual(
                    f[asym][k] + f[image][k], 0.0,
                    delta=1.0e-12,
                    msg=f"orbit {asym} comp {k}: fa={f[asym][k]} fb={f[image][k]}",
                )

    def test_wyckoff_atom_stays_on_mirror(self):
        """An atom declared on a y-mirror Wyckoff site must stay at y=L/2."""
        L = 8.0
        symfile = self._write_json(
            r"""{
                "group_name": "Pm", "lattice": "triclinic",
                "ops": [
                    { "R": [[1, 0, 0], [0,  1, 0], [0, 0, 1]], "t": [0, 0, 0] },
                    { "R": [[1, 0, 0], [0, -1, 0], [0, 0, 1]], "t": [0, 0, 0] }
                ],
                "orbits": [
                    { "tags": [1, 2] },
                    { "tags": [3, 3], "site_symmetry": [2] }
                ]
            }"""
        )
        self._tmpfiles.append(symfile)

        for c in [
            "units lj", "atom_style atomic", "boundary p p p", "atom_modify map array",
            f"region box block 0 {L} 0 {L} 0 {L}", "create_box 1 box",
            "create_atoms 1 single 1.2 2.5 3.0",
            f"create_atoms 1 single 1.2 {L-2.5} 3.0",
            f"create_atoms 1 single 1.5 {0.5*L} 2.0",
            "mass 1 1.0",
            "pair_style lj/cut 2.5", "pair_coeff * * 1.0 1.0",
            "velocity all create 0.7 12345 mom yes rot yes",
            "fix 1 all nve",
            f"fix sym all symmetry {symfile}",
            "run 50",
        ]:
            self.lmp.command(c)

        pos = self._gather_positions()
        # Tag 3 on the mirror.
        self.assertAlmostEqual(pos[3][1], 0.5 * L, delta=1.0e-12,
                               msg=f"Wyckoff atom drifted: y3={pos[3][1]}")
        # Tags 1 and 2 still y-mirror partners.
        self.assertAlmostEqual(pos[1][1] + pos[2][1], L, delta=1.0e-12)
        self.assertAlmostEqual(pos[1][0], pos[2][0], delta=1.0e-12)
        self.assertAlmostEqual(pos[1][2], pos[2][2], delta=1.0e-12)

    def test_p1_does_not_move_atoms(self):
        """Identity-only group must leave all positions exactly fixed."""
        symfile = self._write_json(
            r"""{
                "group_name": "P1", "lattice": "triclinic",
                "ops":   [ { "R": [[1,0,0],[0,1,0],[0,0,1]], "t": [0,0,0] } ],
                "orbits": []
            }"""
        )
        self._tmpfiles.append(symfile)

        # Static system: zero velocity and an empty fix sym should not move atoms.
        for c in [
            "units lj", "atom_style atomic", "boundary p p p", "atom_modify map array",
            "region box block 0 8 0 8 0 8", "create_box 1 box",
            "create_atoms 1 single 1.2 1.5 1.8",
            "create_atoms 1 single 6.8 6.5 6.2",
            "mass 1 1.0",
            "pair_style lj/cut 2.5", "pair_coeff * * 1.0 1.0",
            "velocity all set 0 0 0",
            "fix 1 all nve",
            f"fix sym all symmetry {symfile}",
        ]:
            self.lmp.command(c)
        before = self._gather_positions()
        self.lmp.command("run 10")
        after = self._gather_positions()
        for t in before:
            for k in range(3):
                self.assertAlmostEqual(
                    before[t][k], after[t][k], delta=1.0e-14,
                    msg=f"P1 fix moved atom {t} component {k}",
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)

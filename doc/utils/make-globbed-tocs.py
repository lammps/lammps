#!/usr/bin/env python3

"""
Script to emulate globbed toctrees with the added feature
to avoid adding files that are already included elsewhere.
"""

import os, re
import filecmp
import tempfile
import shutil
from glob import glob
from argparse import ArgumentParser

LAMMPS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

parser = ArgumentParser(prog='make-globbed-tocs.py',
                        description="Create toctree files from patterns with exclusions")
parser.add_argument("-v", "--verbose", action='store_true', help="Enable verbose output")
parser.add_argument("-d", "--doc", help="Path to LAMMPS documentation sources")

args = parser.parse_args()
verbose = args.verbose
if args.doc:
    docsrc = os.path.realpath(args.doc)
else:
    docsrc = os.path.join(LAMMPS_DIR, 'doc', 'src')

if not os.path.isdir(docsrc):
    sys.exit(f"LAMMPS manual source dir {docsrc} does not exist")

def glob_tocfile(style, name, head, exclude):

    newname = None
    exclude_re = re.compile(exclude)
    if verbose: print("Processing style ", style)
    with tempfile.NamedTemporaryFile(prefix=style + '.', delete=False) as f:
        newname = f.name
        if verbose: print("Temporary file: ", newname)
        f.write(head.encode('utf-8'))
        for doc in sorted(glob(os.path.join(docsrc, style + '_*.rst'))):
            d,e = os.path.splitext(os.path.basename(doc))
            if exclude_re.match(d):
                if verbose: print("skipping file ", d)
                continue
            if verbose: print("processing file ", d)
            f.write(("   " + d + "\n").encode('utf-8'))

    oldname = os.path.join(docsrc, name)
    if os.path.exists(oldname) and filecmp.cmp(newname, oldname, shallow=False):
        print("File ", name, " is unchanged")
        os.remove(newname)
    else:
        print("Overwriting file ", name, " with new version")
        shutil.move(newname, os.path.join(docsrc, name))

##################################

pair_head = """Pair Styles
###########

.. toctree::
   :maxdepth: 1

"""
glob_tocfile('pair', 'pairs.rst', pair_head, r"pair_(coeff|modify|style|write)")

bond_head = """Bond Styles
###########

.. toctree::
   :maxdepth: 1

"""
glob_tocfile('bond', 'bonds.rst', bond_head, r"bond_(coeff|modify|style|write)")

angle_head = """Angle Styles
############

.. toctree::
   :maxdepth: 1

"""
glob_tocfile('angle', 'angles.rst', angle_head, r"angle_(coeff|modify|style|write)")

dihedral_head = """Dihedral Styles
###############

.. toctree::
   :maxdepth: 1

"""
glob_tocfile('dihedral', 'dihedrals.rst', dihedral_head, r"dihedral_(coeff|modify|style|write)")

improper_head = """Improper Styles
###############

.. toctree::
   :maxdepth: 1

"""
glob_tocfile('improper', 'impropers.rst', improper_head, r"improper_(coeff|modify|style|write)")

compute_head = """Compute Styles
###############

.. toctree::
   :maxdepth: 1

"""
glob_tocfile('compute', 'computes.rst', compute_head, r"compute_modify")

fix_head = """Fix Styles
##########

.. toctree::
   :maxdepth: 1

"""
glob_tocfile('fix', 'fixes.rst', fix_head, r"fix_modify(|_atc_commands)")

#!/bin/bash
# Regenerate codemeta.json from CITATION.cff.
#
# CITATION.cff is the maintained source of truth for the software metadata
# (version, release date, authors, license, keywords).  This script validates
# it, converts it to CodeMeta with cffconvert, and merges in the fields that
# the Citation File Format cannot express.  Run it after updating CITATION.cff
# for a release and commit the result; the metadata-check workflow fails when
# codemeta.json is out of sync with CITATION.cff.

set -euo pipefail
cd "$(dirname "$0")"

if ! command -v cffconvert > /dev/null 2>&1; then
    echo "cffconvert is required. Install it with: python3 -m pip install --user cffconvert" >&2
    exit 1
fi

cffconvert --validate -i CITATION.cff

cffconvert -f codemeta -i CITATION.cff | python3 -c '
import json
import sys

meta = json.load(sys.stdin)

# fields with no equivalent in the Citation File Format; development of
# LAMMPS started in the mid-1990s in Fortran. The current source code
# repository for the C++ version was initially managed by subversion
# and later converted to git
meta.update({
    "dateCreated": "2006-09-21",
    "programmingLanguage": ["C++", "Python", "CUDA", "Fortran", "C"],
    "operatingSystem": ["Linux", "macOS", "Windows", "FreeBSD"],
    "issueTracker": "https://github.com/lammps/lammps-gui/issues",
    "developmentStatus": "active",
})

with open("codemeta.json", "w") as out:
    json.dump(meta, out, indent=2, sort_keys=True)
    out.write("\n")
'
echo "codemeta.json updated"

#!/bin/bash
# script to emulate globbed toctrees that allows to skip files that were entered already elsewhere

if [ $# != 1 ]
then
    echo "Usage: $0 <rstdir>"
    exit 1
fi

RSTDIR="$1"
TMPNAM=${RANDOM}
TMPDIR=${TMPDIR-/tmp}

# pairs.rst
cp ${RSTDIR}/pairs.rst.in ${TMPDIR}/${TMPNAM}.pairs.rst
for f in $(echo ${RSTDIR}/pair_*.rst | sed -e "s@${RSTDIR}/@@g" -e 's@\.rst@@g' -e 's@pair_\(coeff\|modify\|style\|write\)@@g' | sort )
do \
    echo "   $f" >>  ${TMPDIR}/${TMPNAM}.pairs.rst
done
cmp -s ${TMPDIR}/${TMPNAM}.pairs.rst ${RSTDIR}/pairs.rst || mv -vf ${TMPDIR}/${TMPNAM}.pairs.rst ${RSTDIR}/pairs.rst
rm -f  ${TMPDIR}/${TMPNAM}.pairs.rst

# bonds.rst
cp ${RSTDIR}/bonds.rst.in ${TMPDIR}/${TMPNAM}.bonds.rst
for f in $(echo ${RSTDIR}/bond_*.rst | sed -e "s@${RSTDIR}/@@g" -e "s@\\.rst@@g" -e 's@bond_\(coeff\|modify\|style\|write\)@@g' | sort )
do \
    echo "   $f" >>  ${TMPDIR}/${TMPNAM}.bonds.rst
done
cmp -s ${TMPDIR}/${TMPNAM}.bonds.rst ${RSTDIR}/bonds.rst || mv -vf ${TMPDIR}/${TMPNAM}.bonds.rst ${RSTDIR}/bonds.rst
rm -f  ${TMPDIR}/${TMPNAM}.bonds.rst

# angles.rst
cp ${RSTDIR}/angles.rst.in ${TMPDIR}/${TMPNAM}.angles.rst
for f in $(echo ${RSTDIR}/angle_*.rst | sed -e "s@${RSTDIR}/@@g" -e "s@\\.rst@@g" -e 's@angle_\(coeff\|modify\|style\|write\)@@g' | sort )
do \
    echo "   $f" >>  ${TMPDIR}/${TMPNAM}.angles.rst
done
cmp -s ${TMPDIR}/${TMPNAM}.angles.rst ${RSTDIR}/angles.rst || mv -vf ${TMPDIR}/${TMPNAM}.angles.rst ${RSTDIR}/angles.rst
rm -f  ${TMPDIR}/${TMPNAM}.angles.rst

# dihedrals.rst
cp ${RSTDIR}/dihedrals.rst.in ${TMPDIR}/${TMPNAM}.dihedrals.rst
for f in $(echo ${RSTDIR}/dihedral_*.rst | sed -e "s@${RSTDIR}/@@g" -e "s@\\.rst@@g" -e 's@dihedral_\(coeff\|modify\|style\|write\)@@g' | sort )
do \
    echo "   $f" >>  ${TMPDIR}/${TMPNAM}.dihedrals.rst
done
cmp -s ${TMPDIR}/${TMPNAM}.dihedrals.rst ${RSTDIR}/dihedrals.rst || mv -vf ${TMPDIR}/${TMPNAM}.dihedrals.rst ${RSTDIR}/dihedrals.rst
rm -f  ${TMPDIR}/${TMPNAM}.dihedrals.rst

# impropers.rst
cp ${RSTDIR}/impropers.rst.in ${TMPDIR}/${TMPNAM}.impropers.rst
for f in $(echo ${RSTDIR}/improper_*.rst | sed -e "s@${RSTDIR}/@@g" -e "s@\\.rst@@g" -e 's@improper_\(coeff\|modify\|style\|write\)@@g' | sort )
do \
    echo "   $f" >>  ${TMPDIR}/${TMPNAM}.impropers.rst
done
cmp -s ${TMPDIR}/${TMPNAM}.impropers.rst ${RSTDIR}/impropers.rst || mv -vf ${TMPDIR}/${TMPNAM}.impropers.rst ${RSTDIR}/impropers.rst
rm -f  ${TMPDIR}/${TMPNAM}.impropers.rst

# computes.rst
cp ${RSTDIR}/computes.rst.in ${TMPDIR}/${TMPNAM}.computes.rst
for f in $(echo ${RSTDIR}/compute_*.rst | sed -e "s@${RSTDIR}/@@g" -e "s@\\.rst@@g" -e 's@compute_modify@@g' | sort )
do \
    echo "   $f" >>  ${TMPDIR}/${TMPNAM}.computes.rst
done
cmp -s ${TMPDIR}/${TMPNAM}.computes.rst ${RSTDIR}/computes.rst || mv -vf ${TMPDIR}/${TMPNAM}.computes.rst ${RSTDIR}/computes.rst
rm -f  ${TMPDIR}/${TMPNAM}.computes.rst

# fixs.rst
cp ${RSTDIR}/fixes.rst.in ${TMPDIR}/${TMPNAM}.fixes.rst
for f in $(echo ${RSTDIR}/fix_*.rst | sed -e "s@${RSTDIR}/@@g" -e "s@\\.rst@@g" -e 's@fix_\(modify\|modify_atc_commands\)@@g' | sort )
do \
    echo "   $f" >>  ${TMPDIR}/${TMPNAM}.fixes.rst
done
cmp -s ${TMPDIR}/${TMPNAM}.fixes.rst ${RSTDIR}/fixes.rst || mv -vf ${TMPDIR}/${TMPNAM}.fixes.rst ${RSTDIR}/fixes.rst
rm -f  ${TMPDIR}/${TMPNAM}.fixes.rst

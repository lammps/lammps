#!/bin/bash
# Create a signed list of SHA-256 checksums for the files of a LAMMPS release
#
# All files that are attached to the (draft) release on GitHub are downloaded
# to a temporary folder, so that the checksums are computed for what was
# actually uploaded and with the file names that are used on the release page.
# The list of checksums is written to the file SHA256SUMS and the signature
# for it to the file SHA256SUMS.asc in the current folder.  Both files must be
# uploaded *before* the release is published, because the files of a release
# can no longer be changed after that.

if [ $# -ne 1 ]
then
    echo "usage: $0 <release tag>"
    exit 1
fi
tag=$1
repo=lammps/lammps

# sign with the same program that git uses for signing the release tags
gpg=$(git config --get gpg.program || echo gpg)

for cmd in gh sha256sum ${gpg}
do
    if ! type ${cmd} > /dev/null 2>&1
    then
        echo "Need the '${cmd}' command to run this script"
        exit 2
    fi
done

tmpdir=$(mktemp -d) || exit 3
trap 'rm -rf "${tmpdir}"' EXIT

# download all files of the release and ignore checksum files from a previous run
echo "Downloading files of release ${tag}"
gh release download ${tag} --repo ${repo} --dir "${tmpdir}" || exit 4
rm -f "${tmpdir}/SHA256SUMS" "${tmpdir}/SHA256SUMS.asc"
if [ -z "$(ls -A "${tmpdir}")" ]
then
    echo "Release ${tag} has no files"
    exit 4
fi

# create list of checksums and sign it
rm -f SHA256SUMS SHA256SUMS.asc
(cd "${tmpdir}" && sha256sum -- *) > SHA256SUMS || exit 5
${gpg} --armor --detach-sign --output SHA256SUMS.asc SHA256SUMS || exit 6

# confirm that the signature and the checksums are valid
${gpg} --verify SHA256SUMS.asc SHA256SUMS || exit 7
sumfile="${PWD}/SHA256SUMS"
(cd "${tmpdir}" && sha256sum -c "${sumfile}") || exit 8

echo
echo "Created SHA256SUMS and SHA256SUMS.asc.  Upload them to the release with:"
echo "gh release upload ${tag} --repo ${repo} --clobber SHA256SUMS SHA256SUMS.asc"

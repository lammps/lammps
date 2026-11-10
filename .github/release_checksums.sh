#!/bin/bash
# Create a signed list of SHA-256 checksums for the files of a LAMMPS release
#
# The checksums are not computed locally, but taken from the SHA-256 digests
# that GitHub computes for each uploaded file, so that they belong to what was
# actually uploaded and use the file names of the release page.
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

for cmd in gh ${gpg}
do
    if ! type ${cmd} > /dev/null 2>&1
    then
        echo "Need the '${cmd}' command to run this script"
        exit 2
    fi
done

tmpfile=$(mktemp) || exit 3
trap 'rm -f "${tmpfile}"' EXIT

# get the digests of all files of the release, but not of the checksum files
# from a previous run.  each line is "<upload state>/<digest>  <file name>".
echo "Getting checksums of the files of release ${tag} from GitHub"
gh release view ${tag} --repo ${repo} --json assets --jq \
   '.assets | sort_by(.name | ascii_downcase) | .[]
    | select(.name != "SHA256SUMS" and .name != "SHA256SUMS.asc")
    | "\(.state)/\(.digest)  \(.name)"' > "${tmpfile}" || exit 4
if [ ! -s "${tmpfile}" ]
then
    echo "Release ${tag} has no files"
    exit 4
fi
if grep -Ev '^uploaded/sha256:[0-9a-f]{64}  ' "${tmpfile}"
then
    echo "No SHA-256 checksum available for the file(s) above"
    exit 5
fi

# create list of checksums and sign it
rm -f SHA256SUMS SHA256SUMS.asc
sed -e 's,^uploaded/sha256:,,' "${tmpfile}" > SHA256SUMS || exit 5
${gpg} --armor --detach-sign --output SHA256SUMS.asc SHA256SUMS || exit 6

# confirm that the signature is valid
${gpg} --verify SHA256SUMS.asc SHA256SUMS || exit 7

echo
echo "Created SHA256SUMS and SHA256SUMS.asc.  Upload them to the release with:"
echo "gh release upload ${tag} --repo ${repo} --clobber SHA256SUMS SHA256SUMS.asc"

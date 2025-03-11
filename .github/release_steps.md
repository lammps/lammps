# LAMMPS Release Steps

The following notes chronicle the current steps for preparing and
publishing LAMMPS releases.  For definitions of LAMMPS versions and
releases, please refer to [the corresponding section in the LAMMPS
manual](https://docs.lammps.org/Manual_version.html).

## LAMMPS Feature Release

A LAMMPS feature release is currently prepared after about 500 to 750
commits to the 'develop' branch or after a period of four weeks up to
two months.  This is not a fixed rule, though, since external
circumstances can cause delays in preparing a release, or pull requests
that are desired to be merged for the release are not yet completed.

### Preparing a 'next\_release' branch

Create a 'next\_release' branch off 'develop' and make the following changes:

- set the LAMMPS\_VERSION define to the planned release date in
  src/version.h in the format "D Mmm YYYY" or "DD Mmm YYYY"
- remove the LAMMPS\_UPDATE define in src/version.h
- update the release date in doc/lammps.1
- update all TBD arguments for ..versionadded::, ..versionchanged::
  ..deprecated:: to the planned release date in the format "DMmmYYYY" or
  "DDMmmYYYY"
- check release notes for merged new features and check if
  ..versionadded:: or ..versionchanged:: are missing and need to be
  added

Submit this pull request.  This is the last pull request merged for the
release and should not contain any other changes. (Exceptions: this
document, last minute trivial(!) changes).

This PR shall not be merged before **all** pending tests have completed
and cleared.  We currently use a mix of automated tests running on
either Temple's Jenkins cluster or GitHub workflows.  Those include time
consuming tests not run on pull requests.  If needed, a bug-fix pull
request should be created and merged to clear all tests.

### Create release on GitHub

When all pending pull requests for the release are merged and have
cleared testing, the 'next\_release' branch is merged into 'develop'.

Check out or update the 'develop' branch locally, pull the latest
changes, merge them into 'release' with a fast forward(!) merge, and
apply a suitable release tag (for historical reasons the tag starts with
"patch_" followed by the date, and finally push everything back to
GitHub.  There should be no commits made to 'release' but only
fast forward merges.  Example:

```
git checkout develop
git pull
git checkout release
git pull
git merge --ff-only develop
git tag -s -m "LAMMPS feature release 4 February 2025" patch_4Feb2025
git push git@github.com:lammps/lammps.git --tags develop release
```

Applying this tag will trigger two actions on the Temple Jenkins cluster:
- The online manual at https://docs.lammps.org/ will be updated to the
  state of the 'release' branch.  Merges to the 'develop' branch will
  trigger updating https://docs.lammps.org/latest/ so by reviewing the
  version of the manual under the "latest" URL, it is possible to preview
  what the updated release documentation will look like.
- A downloadable tar archive of the LAMMPS distribution that includes the
  html format documentation and a PDF of the manual will be created and
  uploaded to the download server at https://download.lammps.org/tars
  Note that the file is added, but the `index.html` file is not updated,
  so it is not yet publicly visible.

Go to https://github.com/lammps/lammps/releases and create a new (draft)
release page with a summary of all the changes included and references
to the pull requests they were merged from or check the existing draft
for any necessary changes from pull requests that were merged but are
not listed.  Then select the applied tag for the release in the "Choose
a tag" drop-down list. Go to the bottom of the list and select the "Set
as pre-release" checkbox.  The "Set as the latest release" button is
reserved for stable releases and updates to them.

If everything is in order, you can click on the "Publish release"
button.  Otherwise, click on "Save draft" and finish pending tasks until
you can return to edit the release page and publish it.

### Prepare pre-compiled packages, update packages to GitHub

A suitable build environment is provided with the
https://download.lammps.org/static/fedora41_musl_mingw.sif container
image.  The corresponding container build definition file is maintained
in the tools/singularity folder of the LAMMPS source distribution.

#### Fully portable static Linux x86_64 non-MPI binaries

The following commands use the Fedora container to build a fully static
LAMMPS installation using a musl-libc cross-compiler, install it into a
`lammps-static` folder, and create a tarball called
`lammps-linux-x86_64-4Feb2025.tar.gz` (or using a corresponding date
with a future release) from the `lammps-static` folder.

``` sh
rm -rf release-packages
mkdir release-packages
cd release-packages
wget https://download.lammps.org/static/fedora41_musl.sif
apptainer shell fedora41_musl.sif
git clone -b release --depth 10 https://github.com/lammps/lammps.git lammps-release
cmake -S lammps-release/cmake -B build-release -G Ninja -D CMAKE_INSTALL_PREFIX=$PWD/lammps-static -D CMAKE_TOOLCHAIN_FILE=/usr/musl/share/cmake/linux-musl.cmake -C lammps-release/cmake/presets/most.cmake -C lammps-release/cmake/presets/kokkos-openmp.cmake -D DOWNLOAD_POTENTIALS=OFF -D BUILD_MPI=OFF -D BUILD_TESTING=OFF -D CMAKE_BUILD_TYPE=Release -D PKG_ATC=ON -D PKG_AWPMD=ON -D PKG_MANIFOLD=ON -D PKG_MESONT=ON -D PKG_MGPT=ON -D PKG_ML-PACE=ON -D PKG_ML-RANN=ON -D PKG_MOLFILE=ON -D PKG_PTM=ON -D PKG_QTB=ON -D PKG_SMTBQ=ON
cmake --build build-release --target all
cmake --build build-release --target install
/usr/musl/bin/x86_64-linux-musl-strip lammps-static/bin/*
tar -czvvf ../lammps-linux-x86_64-4Feb2025.tar.gz lammps-static
exit # fedora 41 container
cd ..
```

The resulting tar archive can be uploaded to the GitHub release page with:

``` sh
gh release upload patch_4Feb2025 lammps-linux-x86_64-4Feb2025.tar.gz
```

#### Linux x86_64 Flatpak bundle with GUI included

Make sure you have the `flatpak` and `flatpak-builder` packages
installed locally (they require binaries that run with elevated
privileges and thus cannot be used from the container) and build a
LAMMPS and LAMMPS-GUI flatpak bundle in the `release-packages` folder
with:

``` sh
cd release-packages
flatpak --user remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
flatpak-builder  --force-clean --verbose --repo=$PWD/flatpak-repo --install-deps-from=flathub --state-dir=$PWD --user --ccache --default-branch=release flatpak-build lammps-release/tools/lammps-gui/org.lammps.lammps-gui.yml
flatpak build-bundle --runtime-repo=https://flathub.org/repo/flathub.flatpakrepo --verbose $PWD/flatpak-repo ../LAMMPS-Linux-x86_64-GUI-4Feb2025.flatpak org.lammps.lammps-gui release
cd ..
```

The resulting flatpak bundle file can be uploaded to the GitHub release page with:

``` sh
gh release upload patch_4Feb2025 LAMMPS-Linux-x86_64-GUI-4Feb2025.flatpak
```

#### LAMMPS Source tarball

The container for the static binary can also be used to prepare the source
tarball including the HTML and PDF manual (this is currently done automatically
when the releases is created and the tarball uploaded to https://download.lammps.org/tars/).
The steps are as follows:

``` sh
cd release-packages
apptainer shell fedora41_musl_mingw.sif
cd lammps-release
rm -f ../release.tar*
git archive --output=../release.tar --prefix=lammps-4Feb2025/ HEAD
cd doc
make clean-all
make html pdf
tar -rf ../../release.tar --transform 's,^,lammps-4Feb2025/doc/,' html Manual.pdf
gzip -9v ../../release.tar
mv ../../release.tar.gz ../../lammps-src-4Feb2025.tar.gz
exit # fedora41 container
cd ..
```

The resulting source tarball can be uploaded to the GitHub release page with:

``` sh
gh release upload patch_4Feb2025 lammps-src-4Feb2025.tar.gz
```

#### Build Windows Installer Packages with MinGW Linux-to-Windows Cross-compiler

The various Windows installer packages can also be built with
apptainer container image.

``` sh
cd release-packages
apptainer shell fedora41_musl_mingw.sif
git clone --depth 10 https://github.com/lammps/lammps-packages.git lammps-packages
cd lammps-packages/mingw-cross
ln -sf ../../lammps-release lammps
./buildall.sh release >& mk.log & less +F mk.log
```

The installer with the GUI included can be uploaded to the GitHub release page with:

``` sh
ln -sf LAMMPS-64bit-GUI-4Feb2025.exe LAMMPS-Win10-64bit-GUI-4Feb2025.exe
gh release upload patch_4Feb2025 LAMMPS-Win10-64bit-GUI-4Feb2025.exe
```

The symbolic link is used to have a consistent naming scheme for the packages
attached to the GitHub release page.

#### Clean up:

``` sh
cd ..
rm -r release-packages
```

#### Build Multi-arch App-bundle for macOS

Building app-bundles for macOS is not as easily automated and portable
as some of the other steps.  It requires a machine actually running
macOS.  In that machine the Xcode compiler package needs to be
installed. This also includes tools for building and manipulating disk
images.  This compiler supports building executables for both, the
x86_64 and the arm64 architectures.  This requires building with CMake
and using the CMake settings:

``` sh
-D CMAKE_OSX_ARCHITECTURES=arm64;x86_64
-D CMAKE_OSX_DEPLOYMENT_TARGER=11.0
```

This will add the compiler flags `-arch arm64 -arch x86_64
-mmacosx-version-min=11.0` and thus produce object for both
architectures and support for macOS versions back to version 11 (aka Big
Sur).  With these settings the following libraries should be compiled
and installed (e.g. to `$HOME/.local`) as static libraries only:
- libomp taken from the LLVM/Clang source distribution (to support OpenMP)
- jpeg
- zlib
- png
- Qt (for LAMMPS-GUI)

When configuring LAMMPS the `cmake/presets/clang.cmake` should be used
and as many packages as possible enabled. For LAMMPS-GUI, MPI should be
disabled with `-D BUILD_MPI=OFF` and LAMMPS-GUI enabled with 
`-D BUILD_LAMMPS_GUI=ON`.  If the CMake configuration is successful,
settings for building a macOS app-bundle are enabled and with `cmake
--build build --target dmg` extra steps will be executed that will build
a macOS application installer image under the name
`LAMMPS_GUI-macOS-multiarch-4Feb2025.dmg`

The application image can be uploaded to the GitHub release page with:

``` sh
ln -sf LAMMPS_GUI-macOS-multiarch-4Feb2025.dmg LAMMPS-macOS-multiarch-GUI-4Feb2025.dmg
gh release upload patch_4Feb2025 LAMMPS-macOS-multiarch-GUI-4Feb2025.dmg
```

The symbolic link is used to have a consistent naming scheme for the packages
attached to the GitHub release page.

We are currently building the application images on macOS 12 (aka Monterey).

#### Build Linux x86_64 binary tarball on Ubuntu 20.04LTS

While the flatpak Linux version uses portable runtime libraries provided
by the flatpak environment, we also build regular Linux executables that
use a wrapper script and matching shared libraries in a tarball.  To be
compatible with many Linux distributions, one has to build this on a
very old Linux distribution, since most Linux system libraries are
usually backward compatible but not forward compatible.  This is
currently done on an Ubuntu 20.04LTS system. Once LAMMPS moves to
require CMake 3.20 and C++17, we will have to move to Ubuntu 22.04LTS.
This installation (either on a real or a virtual machine) should have
the packages installed that are indicated in
`tools/singularity/ubuntu20.04.def` plus Qt version 5.x with development
headers, so that LAMMPS-GUI can be compiled.

Also the building of the binary tarball and setup of the bundled
libraries and wrapper scripts is automated and can executed with `cmake
--build build --target tgz`.  This should produce a file
`LAMMPS_GUI-Linux-amd64-4Feb2025.tar.gz` which can be uploaded to the
GitHub release page with:

``` sh
ln -sf LAMMPS_GUI-Linux-amd64-4Feb2025.tar.gz LAMMPS-Linux-x86_64-GUI-4Feb2025.tar.gz
gh release upload patch_4Feb2025 LAMMPS-Linux-x86_64-GUI-4Feb2025.tar.gz
```

### Update download page on LAMMPS website

Check out the LAMMPS website repo
https://github.com/lammps/lammps-website.git and edit the file
`src/download.txt` for the new release.  Test translation with `make
html` and review `html/download.html` Then add and commit to git and
push the changes to GitHub.  The Temple Jenkis cluster will
automatically update https://www.lammps.org/download.html accordingly.

Also notify Steve of the release so he can update `src/bug.txt` on the
website from the available release notes.

## LAMMPS Stable Release

A LAMMPS stable release is prepared about once per year in the months
July, August, or September.  One (or two, if needed) feature releases
before the stable release shall contain only bug fixes or minor feature
updates in optional packages.  Also substantial changes to the core of
the code shall be applied rather toward the beginning of a development
cycle between two stable releases than toward the end.  The intention is
to stablilize significant change to the core and have outside users and
developers try them out during the development cycle; the sooner the
changes are included, the better chances for spotting peripheral bugs
and issues.

### Prerequesites

Before making a stable release all remaining backported bugfixes shall
be released as a (final) stable update release (see below).

A LAMMPS stable release process starts like a feature release (see
above), only that this feature release is called a "Stable Release
Candidate" and no assets are uploaded to GitHub.

### Synchronize 'maintenance' branch with 'release'

The state of the 'release' branch is then transferred to the
'maintenance' branch (which will have diverged significantly from
'release' due to the selectively backported bug fixes).

### Fast-forward merge of 'maintenance' into 'stable' and apply tag

At this point it should be possible to do a fast-forward merge of
'maintenance' to 'stable' and then apply the stable\_DMmmYYYY tag.

### Push branches and tags



## LAMMPS Stable Update Release

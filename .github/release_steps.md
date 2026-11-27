# LAMMPS Release Steps

The following notes chronicle the current steps for preparing and
publishing LAMMPS and LAMMPS-GUI releases on the build server at Temple
U (a server called `rstudio` for historical reasons).

For definitions of LAMMPS versions and releases, please refer to [the
corresponding section in the LAMMPS
manual](https://docs.lammps.org/Manual_version.html).

## Current setup of build server:

- Dell PowerEdge R640 rack-mounted server
- 2x Intel Xeon Gold 6248 (2.5 GHz base, 3.9 GHz max turbo), 20 cores each with SMT enabled
- 384GB RAM, 2x 450 GiB SATA SSD: 1 GB RAID-1 -> /boot, 880 GiB RAID-0 /
- Rocky Linux 10.2 with EPEL as host operating system
- AlmaLinux 9 container via podman/toolbox, for Linux tar.gz, liblammps.so.0
- Fedora 44 container via podman/toolbox, for MinGW64 cross-compiler,
  flatpak, static binary, source tarball, PDF manual
- macOS 12 Monterey virtual machine via quickemu using qemu/kvm, headless
  must be *launched* from within Fedora44 container but you can enter
  *directly* with ssh to localhost on port 22230 from `rstudio`
  via public key (`ssh hackintosh`), also from ws01 via proxyjump.

## Directory structure

### LAMMPS

Use a checkout of LAMMPS created on Rocky 10 with a 'lammps' folder for
regular work, a worktree 'lammps-release' for feature releases, and a
worktree 'lammps-stable' for stable releases and updates.  A build
without a suitable release tag will build packages using the current
state of the 'lammps' checkout and will only copy results to
`${HOME}/release` but not to the webserver's testing folder.  Checkout
and worktrees are shared across the containers.  As a virtual machine,
macOS requires its own checkout and worktrees but follows the same
directory structure.

Build folders: | Used for:
---|---
*build-pak*    | for Linux flatpak, source tarball, PDF manual
*build-lib*    | for liblammps shared lib for use with LAMMPS-GUI
*build-mingw*  | for liblammps.dll and Windows installers
*build-gui*    | for LAMMPS + LAMMPS-GUI and tar.gz package or dmg installer
*build-static* | for fully static LAMMPS Linux binary package

Install folders: | Used for:
---|---
*lammps-static* | for the directory tree with the static LAMMPS binaries

### LAMMPS-GUI

Use a checkout of LAMMPS-GUI created on Rocky 10 with a 'lammps-gui' folder for
regular work, and a worktree 'lammps-gui-release' for releases.  Git
checkout and worktrees are shared across the containers.  As a virtual
machine, macOS requires its own checkout and worktrees but follows the
same scheme.

Build folders: | Used for:
---|---
*build*       |  for Linux flatpak, source tarball, and PDF manual
*build-mingw* |  for Windows installer package, plugin mode
*build-gui*   |  for tar.gz package or dmg installer, plugin mode

## Collecting release assets

Release assets are collected in $HOME/release before uploading them to
GitHub with `gh release upload`. Release versions are also copied with
`scp` to `publish@www.lammps.org:/var/www/lammps/download/testing` so
they can be downloaded from https://download.lammps.org/testing/ and
tested before finalizing the release and locking the version.

## Scripts for building and uploading release assets

The [lammps-release GitHub repository](https://github.com/lammps/lammps-release)
contains the scripts used for preparing and uploading the tarballs,
binaries and installer packages on the build machine. The
`${HOME}/release/` folder on the build machine is a clone of that
repository.  Thus the scripts are also the authoritative documentation
of the commands that create all release assets.  Together with this
document they should describe the process completely and automate most
of it.  If something is not well explained or details are missing,
please submit an issue in the LAMMPS repository.

Due to the sensitive nature of those scripts revealing details of the
internal setup, the repository is private and access is limited to the
members of the LAMMPS organization.

The scripts for the Linux containers use `/bin/bash` syntax, the
scripts for macOS `/bin/zsh` (which is largely compatible with bash).

### Makefile

Target | Purpose
---|---
check  |  run "shellcheck" on script files
spell  |  run "codespell" on script files
clean  |  delete editor and git backups and temporary files
purge  |  "make clean" and also delete all release asset files

### Creating and maintaining virtual machine and containers

TODO.

### Launching and stopping virtual machine and containers

The `./boot-build.sh` script will activate all required containers and
launch the macOS virtual machine and apply a 2 minute delay to wait
until the virtual machine has booted to the point where one can enter it
with ssh.  As the last step of the boot, all scripts for macOS are
synchronized with the release folder of the virtual machine.  There is
no git checkout since no development of the build scripts is planned to
be done there.

The `./shutdown-build.sh` script performs the inverse process and
first shuts down the macOS virtual machine, waits for 30 seconds
and then stops all active containers.

### macOS virtual machine

``` console
toolbox enter fedora-toolbox-44
./quickemu/macos-monterey.conf
exit
```

This will launch the container and virtual machine within the container.
They will keep running after `exit`.  Wait about 2 minutes until macOS is
fully booted.  Then you can log in password-less with `ssh hackintosh`.

### Fedora 44 container

``` console
toolbox enter fedora-toolbox-44
```

This launches the container and it will keep running.  Launching the
macOS virtual machine or running `./boot-build.sh` also starts it.  The
container has a password-less sudo configuration for software
installation and updates.

### AlmaLinux 9 container

``` console
toolbox enter almalinux-9
```

This launches the container and it will keep running.  Running
`./boot-build.sh` also starts it.  The container has a password-less
sudo configuration for software installation and updates.

## LAMMPS Online Manual updates

The online manual pages are now updated every two hours by a cron job
that runs the scripts in the
[lammps-docs](https://github.com/lammps/lammps-docs/) repository, if
there have been changes to the `stable`, `release`, or `develop`
branches.  The [LAMMPS Testing Dashboard](https://dashboard.lammps.org)
has a panel showing the status of the online manual builds and reports
when the last build succeeded, if there were failed builds, failed spell
checks, or failed synchronizations with the webserver.

## LAMMPS Feature Release

A LAMMPS feature release is currently prepared after about 750 to 1000
commits to the 'develop' branch or after a period of six weeks up to two
months.  This is not a fixed rule, though, since external circumstances
can cause delays in preparing a release, or pull requests that are
desired to be merged for the release are not yet completed.

### Preparing a 'next\_release' branch

Create a 'next\_release' branch off 'develop' and make the following changes:

- set the LAMMPS\_VERSION define to the planned release date in
  `src/version.h` in the format "D Mmm YYYY" or "DD Mmm YYYY"
- remove the LAMMPS\_UPDATE define in src/version.h
- update the release date in `doc/lammps.1`
- check release notes for merged new features and check if
  `.. versionadded::` or `.. versionchanged::` are missing and need to be
  added
- update all TBD arguments for `.. versionadded::`, `.. versionchanged::`,
  `.. versionremoved::`, and `.. deprecated::` to the planned release
  date in the format "DMmmYYYY" or "DDMmmYYYY";
  `make -C src check-version` lists the remaining TBD entries
- check at https://github.com/lammps/lammps-gui/releases whether there
  is a newer release of LAMMPS-GUI.  If there is, update the `GIT_TAG`
  setting for LAMMPS-GUI in `cmake/Modules/Packaging.cmake` to the tag
  of that release, since the pre-compiled packages with LAMMPS-GUI
  included are built from that tag.  Only the tag of a published
  LAMMPS-GUI release may be used and never the name of a branch: the
  releases of LAMMPS-GUI are immutable, so that the tag cannot be
  changed later.  This can be confirmed with
  `gh release verify --repo lammps/lammps-gui <tag>`

Submit this pull request.  This is the last pull request merged for the
release and should not contain any other changes. (Exceptions: internal
documents like this one and last-minute trivial(!) changes).

This pull request shall not be merged before **all** pending tests have
completed and cleared.  We currently use automated tests with GitHub
workflows and automated nightly tests.  There is also a Coverity Scan
run with static code analysis done every Thursday evening or when
triggered manually.  The various test outputs are consolidated on the
[LAMMPS Testing Dashboard](https://dashboard.lammps.org)
website and it should be consulted to decide which issues need to be
addressed before the source code is suitable for a release.

### Create a draft release on GitHub

When all pending pull requests for the release are merged and have
cleared testing, the 'next\_release' branch is merged into 'develop'.

Check out or update the 'develop' branch locally, pull the latest
changes, merge them into 'release' with a fast-forward(!) merge, and
apply a suitable release tag (for historical reasons the tag starts with
"patch_" followed by the version date), and finally push the 'release'
branch and the tag to GitHub.  There should be no commits made to
'release' but only fast-forward merges.  Example:

``` console
git checkout develop
git pull
git checkout release
git pull
git merge --ff-only develop
git tag -s -m "LAMMPS feature release 4 February 2025" patch_4Feb2025
git verify-tag patch_4Feb2025
git push git@github.com:lammps/lammps.git release patch_4Feb2025
```

The `git verify-tag` command must report a good signature.  If it does
not, the tag was created without a signature and must be deleted and
created again.

Go to https://github.com/lammps/lammps/releases and create a new draft
release page with a summary of all the changes included and references
to the pull requests they were merged from or check the existing draft
for any necessary changes from pull requests that were merged but are
not listed.  The list of changes should follow a `## Changelog`
headline, since automated assessments of the project search the release
notes for this word.  Then select the applied tag for the release in the
"Choose a tag" drop-down list.  Go to the bottom of the list and select
the "Set as pre-release" checkbox.  The "Set as the latest release"
button is reserved for stable releases and updates to them.

Releases are now "immutable", so neither the tag hash nor any uploaded
assets can be changed after the release is published.  Only the release
notes text may be changed, e.g. to document issues with the uploaded
assets.  Thus only when *everything* is in order *and* all required
assets (source, binary packages for different platforms, manual PDF,
signed list of checksums) are uploaded in a suitable version, you can
click on the "Publish release" button.  Otherwise, click on "Save draft"
and finish pending tasks until you can return to edit the release page,
update assets, and publish it when it is ready.

### Prepare pre-compiled packages, upload packages to GitHub

Creating pre-compiled packages is now largely scripted using scripts in
the [LAMMPS release scripts
repository](https://github.com/lammps/lammps-release).  The following
paragraphs list the individual steps and highlight important commands.
However, the full build process for a feature release can be
performed by executing `./create-all-lammps.sh release` in the
``${HOME}/release`` folder of the build machine. This essentially executes:

``` bash
branch=release
./boot-build.sh
ssh hackintosh ./release/create-dmg-macos.sh "${branch}"
toolbox run -c almalinux-9 "${PWD}/create-tgz-linux.sh" "${branch}"
toolbox run -c fedora-toolbox-44 "${PWD}/create-static-linux.sh" "${branch}"
toolbox run -c fedora-toolbox-44 "${PWD}/create-pak-other.sh" "${branch}"
toolbox run -c fedora-toolbox-44 "${PWD}/create-exe-windows.sh" "${branch}"
./shutdown-build.sh
```

Here are the individual steps with some additional explanations.

#### Fully portable static Linux x86_64 non-MPI binaries

The `./create-static-linux.sh release` command when run in the Fedora 44
container as described earlier builds a fully static LAMMPS installation
using a musl-libc cross-compiler, installs it into a `lammps-static`
folder, and creates a tarball called
`lammps-linux-x86_64-<version>.tar.gz` (with `<version>` standing in for
the actual release date) from the `lammps-static` folder.

The core steps from the script are the following:

``` bash
make -C doc clean-all
rm -rf build-static lammps-static
cmake -S cmake -B build-static -G Ninja -D FFT=KISS -D BUILD_MPI=OFF \
      -D CMAKE_INSTALL_PREFIX="${PWD}/lammps-static" -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_TOOLCHAIN_FILE=/usr/musl/share/cmake/linux-musl.cmake \
      -C cmake/presets/most.cmake -C cmake/presets/kokkos-openmp.cmake \
      -D USE_INTERNAL_LINALG=ON  -D DOWNLOAD_POTENTIALS=OFF -D DOWNLOAD_VORO=ON \
      -D PKG_APIP=ON -D PKG_MANIFOLD=ON -D PKG_MESONT=ON -D PKG_MGPT=ON \
      -D PKG_ML-PACE=ON -D PKG_PTM=ON -D PKG_QTB=ON -D PKG_ML-RUNNER=ON \
      -D PKG_SMTBQ=ON
cmake --build build-static --target all
cmake --build build-static --target install
/usr/musl/bin/x86_64-linux-musl-strip -g lammps-static/bin/*
tarball="build-static/lammps-linux-x86_64-<version>.tar.gz"
tar -czvvf "${tarball}" lammps-static
```

#### LAMMPS Source tarball and PDF Manual

The scripts to build the LAMMPS source tarball (which includes HTML and
the PDF versions of the manual) and a standalone PDF manual are
integrated into the general CMake build as custom targets.  Thus, for
scripted release packaging we have to set up a build folder, but without
compiling LAMMPS itself.  The `./create-pak-other.sh release` command
when run in the Fedora 44 container creates
`lammps-src-<version>.tar.gz` (with `<version>` standing in for the
actual release date) and `Manual-<version>.pdf`.

The core steps from the script are the following:

``` bash
make -C doc clean-all
rm -rf build-pak
cmake -S cmake -B build-pak -D CMAKE_BUILD_TYPE=Release
cmake --build build-pak --target gitversion
cmake --build build-pak --target tarball
```

#### Linux x86_64 Flatpak bundle with GUI included

Creating the [Flatpak](https://flatpak.org) application bundle follows
the same method as for the source tarball, so its creation is included
when running `./create-pak-other.sh release` mentioned above as
`LAMMPS-Linux-x86_64-GUI-<version>.flatpak`.

The core steps from the script are the following:

``` bash
make -C doc clean-all
rm -rf build-pak
cmake -S cmake -B build-pak -D CMAKE_BUILD_TYPE=Release || exit 5
cmake --build build-pak --target gitversion || exit 6
cp cmake/packaging/org.lammps.lammps-gui.yml cmake/packaging/org.lammps.lammps-gui.yml.bak
sed -i -e "s/branch: .*/branch: ${branch}/"  cmake/packaging/org.lammps.lammps-gui.yml
cmake --build build-pak --target flatpak || exit 9
mv cmake/packaging/org.lammps.lammps-gui.yml.bak cmake/packaging/org.lammps.lammps-gui.yml
```

#### Build Linux x86_64 binary tarball with GUI on AlmaLinux 9

While the LAMMPS flatpak bundle for Linux uses its own runtime provided
by the flatpak environment, we also build regular Linux executables that
use a wrapper script with matching shared libraries in a tarball.  To be
compatible with many Linux distributions, one has to build this on an
older Linux distribution that still complies with the build requirements
of LAMMPS and LAMMPS-GUI and then depend on the backward compatibility
of the core Linux system libraries for portability.  This is currently
RHEL 9, because its GNU libc version is one step older than the next closest
match, Ubuntu 22.04LTS. We use an AlmaLinux 9 container to stand in for
RHEL 9.

The `./create-tgz-linux.sh release` command when run in the AlmaLinux 9
container creates `LAMMPS-Linux-x86_64-GUI-<version>.tar.gz` (with
`<version>` standing in for the actual release date).  The core steps
from the script are the following:

``` bash
make -C doc clean-all
rm -rf build-gui
cmake -S cmake -B build-gui -G Ninja -C cmake/presets/lammps-gui.cmake \
      -D DOWNLOAD_VORO=ON -D ENABLE_TESTING=OFF \
      -D BUILD_LAMMPS_GUI=ON -D LAMMPS_GUI_USE_PLUGIN=OFF
cmake --build build-gui --target tgz
```

Note that we have to force downloading the voro++ library so that
CMake won't pick up a version that may be installed locally to support
some other application.  The detailed steps of determining the necessary
shared libraries for the wrapper script are part of the LAMMPS sources
and executed with the `tgz` custom target.

#### Build multi-arch app-bundle with GUI for macOS

Creating app-bundles for macOS currently requires either running on
Apple hardware or a virtual machine (which will only work with x86_64
CPUs).  The system that Apple uses for applications is a bit of a hybrid
between the flatpak and the wrapper script approach used on Linux.  For
maximum compatibility, we are currently using a virtual machine running
macOS 12 Monterey which has the native Xcode compiler package installed.
This also includes tools for building and manipulating disk images.  On
top of this we use the dmgbuild Python package.  The Xcode compiler is
derived from LLVM/Clang and supports building executables for both the
x86_64 and the arm64 architectures.  This requires building with CMake
and using the CMake settings:

``` sh
-D CMAKE_OSX_ARCHITECTURES="arm64;x86_64"
-D CMAKE_OSX_DEPLOYMENT_TARGET=12.0
```

This will add the compiler flags `-arch arm64 -arch x86_64
-mmacosx-version-min=12.0` and thus produce objects for both
architectures and support for macOS versions back to version 12.  On top
of this the virtual machine has a host of additional libraries and
packages compiled locally with the same settings, most importantly
libomp (taken from the LLVM/Clang distribution) to support OpenMP (the
compiler supports processing OpenMP directives, but the Xcode
distribution is missing the header and runtime library).  The
command `./create-dmg-macos.sh release` performs the necessary steps:

``` bash
rm -rf build-gui
cmake -S cmake -B build-gui -D CMAKE_BUILD_TYPE=Release -G Ninja -D FFT=KISS \
      -C cmake/presets/macos-multiarch.cmake -C cmake/presets/kokkos-openmp.cmake \
      -C cmake/presets/most.cmake  -D BUILD_LAMMPS_GUI=ON -D DOWNLOAD_POTENTIALS=OFF \
      -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=OFF -D USE_INTERNAL_LINALG=ON \
      -D CMAKE_CXX_COMPILER_LAUNCHER=ccache -D CMAKE_C_COMPILER_LAUNCHER=ccache \
      -D PKG_ML-PACE=ON -D PKG_APIP=ON -D PKG_MANIFOLD=ON -D PKG_MGPT=ON \
      -D PKG_MOLFILE=ON -D PKG_PTM=ON -D PKG_QTB=ON -D PKG_SMTBQ=ON \
      -D DOWNLOAD_VORO:BOOL=OFF -D VORO_INCLUDE_DIR:PATH=/Users/akohlmey/.local/include/voro++ \
      -D VORO_LIBRARY:FILEPATH=/Users/akohlmey/.local/lib/libvoro++.a \
      -D OpenMP_libomp_LIBRARY:FILEPATH=/Users/akohlmey/.local/lib/libomp.a
cmake --build build-gui --target dmg
```

It is important to note that in this case, and unlike in the Linux case,
voro++ *must not* be downloaded because the automated compilation will
not apply the multi-arch settings.

#### Build Windows Installer Packages with MinGW Linux-to-Windows Cross-compiler

The various Windows installer packages can also be built with the Fedora 44
container.  You need to run `./create-exe-windows.sh release` to
configure a build folder and then run the
`./cmake/packaging/cmake-win-on-linux.py` script multiple times with the
suitable settings.  Unlike for LAMMPS-GUI, support for creating signed
executables is not yet tested and thus currently disabled.

The core steps from the script are the following:

``` bash
make -C doc clean-all
rm -rf build-mingw
mingw64-cmake -S cmake -B build-mingw -G Ninja -D FFT=KISS \
              -C cmake/presets/lammps-gui.cmake -D CMAKE_BUILD_TYPE=Release
cmake --build build-mingw --target gitversion
pushd build-mingw
env SIGN_DISABLE=1 ../cmake/packaging/cmake-win-on-linux.py -p no -y no -u no
env SIGN_DISABLE=1 ../cmake/packaging/cmake-win-on-linux.py -p no -y yes -u no
env SIGN_DISABLE=1 ../cmake/packaging/cmake-win-on-linux.py -p ms -y no -u no
env SIGN_DISABLE=1 ../cmake/packaging/cmake-win-on-linux.py -p ms -y yes -u no
env SIGN_DISABLE=1 ../cmake/packaging/cmake-win-on-linux.py -p no -y no -u yes
popd
```

### Upload release assets to GitHub and download.lammps.org

If completed successfully, the scripts explained above will have
copied all created release assets to the ``${HOME}/release`` folder
on the build machine and to https://download.lammps.org/testing/.

In the ``${HOME}/release`` folder all Windows installer packages
except for `LAMMPS-Win10-x64-GUI-<version>.exe` can be deleted.
The remaining packages can be uploaded to the draft release using
the GitHub CLI, gh.  You can check for the available draft with:

``` bash
gh -R lammps/lammps release list
```

The upload then is done with, e.g.:

``` bash
gh -R lammps/lammps release upload --clobber patch_<version> \
    lammps-src-<version>.tar.gz \
    lammps-linux-x86_64-<version>.tar.gz \
    LAMMPS-Linux-x86_64-GUI-<version>.flatpak \
    LAMMPS-Linux-x86_64-GUI-<version>.tar.gz \
    LAMMPS-macOS-multiarch-GUI-<version>.dmg \
    LAMMPS-Win10-x64-GUI-<version>.exe \
    Manual-<version>.pdf
```

In addition the `lammps-gui-v<guiversion>.pdf` file from the
corresponding LAMMPS-GUI version needs to be downloaded from the
LAMMPS-GUI release and uploaded here.

On the web server, the `lammps-src-<version>.tar.gz` and
`Manual-<version>.pdf` files need to be moved to the `tars` folder, renamed
as needed and symlinks and SHA256SUMS files updated.  The
`lammps-linux-x86_64-<version>.tar.gz` file goes into the `static`
folder, where its SHA256SUMS file is updated as well.  All package
files (including the ones deleted in ``${HOME}/release``) need to be
moved into their corresponding subfolders of the `packages` folder and
symlinks and SHA256SUMS files updated. Once this is completed to
satisfaction, the `index.html` files need to be updated by running the
`mkhtmlindex.sh` scripts (`./mkhtmlindex.sh tars` and `./mkhtmlindex.sh
linux macos windows flatpak`).

### Create and upload signed list of checksums

This is the last step before publishing the release and requires that
*all* other files have been uploaded to the GitHub release page.  The
script `.github/release_checksums.sh` takes the SHA-256 checksums that
GitHub computed for all files of the release page (no download is
needed), records them in a file `SHA256SUMS`, and signs this file with
the same GPG key that is used for signing the release tag, which
creates the file `SHA256SUMS.asc`.  Users can then confirm with these
two files that their downloads are complete and unmodified (see the
`SECURITY.md` file for the commands).

``` sh
bash .github/release_checksums.sh patch_4Feb2025
gh release upload patch_4Feb2025 --clobber SHA256SUMS SHA256SUMS.asc
```

If a file of the draft release is replaced later, these two commands
must be repeated.  They must be run from a checkout where the signing
key for the git tags is configured to use the default GPG key used for
signing the checksum file.  The script writes both files into the
current folder; delete them there after the upload.

### Publish the release

If all previous steps were completed successfully and all tests of the
draft release worked out, you can now finalize the release by publishing
it.  As explained before, **this step is irreversible**.  After
publishing, a release may be deleted, but the release tag cannot be
deleted or replaced, so all release steps need to be repeated with a new
release date.

After publishing, add `#define LAMMPS_UPDATE "Development"` back to
`src/version.h` on the 'develop' branch to flag it as development
version again.

### Update download page on LAMMPS website

Check out the LAMMPS website repo at https://github.com/lammps/website.git
and edit the file `src/content/download/_index.md` for the new release.
Test the site build with `make check`, `make lint` and `make spell`.
You can check the page locally with `make serve` and opening the
presented localhost URL; any edits will be tracked live.

Also add a release announcement to `src/content/news/` using one of the
existing announcements as a template.  Check and review as before.

Then add and commit to git and push the changes to GitHub.
A cron job will automatically update https://www.lammps.org/ accordingly.

Also notify Steve of the release so he can update
`src/content/about/contributors.md` from the available release notes.

## LAMMPS Stable Release

A LAMMPS stable release is prepared about once per year in the months
July, August, or September.  One (or two, if needed) feature releases
before the stable release shall contain only bug fixes or minor feature
updates in optional packages.  Also substantial changes to the core of
the code shall be applied rather toward the beginning of a development
cycle between two stable releases than toward the end.  The intention is
to stabilize significant changes to the core and have outside users and
developers try them out during the development cycle; the sooner the
changes are included, the better chances for spotting peripheral bugs
and issues.

### Prerequisites

Before making a new stable release all remaining backported bugfixes
shall be released as a (final) stable update for the current stable
release.

A LAMMPS stable release process starts like a feature release (see
above), only that this feature release is called a "Stable Release
Candidate" and no assets are uploaded to GitHub.

The pull request for the release candidate is also the time to review
the file `.github/security-insights.yml`, which describes the project
to automated tools (contacts, documentation links, list of project
repositories, checking tools in use).  Correct what has changed and set
the `last-reviewed` date (and `last-updated`, if there were changes).

### Create feature release first

Create a LAMMPS feature release with the corresponding `patch_<version>`
tag first and follow the usual scheme, but do *not* create and upload
any assets.  The release notes should cover the changes since the
previous feature release, and the release is flagged as pre-release
as usual.

### Synchronize 'maintenance' branch with 'release'

The state of the 'release' branch is then transferred to the
'maintenance' branch (which will have diverged significantly from
'release' due to the selectively backported bug fixes).  Merge 'release'
into 'maintenance', resolve all conflicts in favor of 'release', and
remove any remaining differences with a follow-up commit.  Afterwards,
`git diff release maintenance` must only show the `branch: stable`
setting in `cmake/packaging/org.lammps.lammps-gui.yml`.  Pay attention
to the `src/version.h` file.  It must not have a `LAMMPS_UPDATE` line.

### Merge 'maintenance' into 'stable' and apply tag

Now merge 'maintenance' into 'stable' (a fast-forward merge, if 'stable'
has no commits that are missing in 'maintenance', otherwise a regular
merge) and then apply the stable\_DMmmYYYY tag.  Since 'stable' is
protected, pushing it directly instead of through a pull request
requires administrator rights.

### Push branches and tags and create draft release

The stable release notes provide a summary of the changes over
the last year since the last stable release.  The stable releases
(and stable updates) should be flagged as the "latest release".

### Create and upload release assets and publish release

From here on the steps are pretty much the same as for a feature release,
only that the scripts get the 'stable' argument instead of 'release'.
After the release is published, set `LAMMPS_UPDATE` in `src/version.h`
on 'maintenance' back to "Maintenance".

## LAMMPS Stable Update Release

After making a stable release, bugfixes from the 'develop' branch are
selectively backported to the 'maintenance' branch.  This is done with
`git cherry-pick <commit hash>` wherever possible.  When fixing bugs
with the help of an AI coding assistant, the backport is usually also
requested, if the fix applies.  The LAMMPS\_UPDATE define in
"src/version.h" is usually set to "Maintenance" for the 'maintenance'
branch to signal its provenance, similar to "Development" for 'develop'.

### Prerequisites

To simplify keeping track of the backported changes, a pull request is
created that would merge the `maintenance` branch into `stable` and the
pull request description serves as a draft for the release notes. As
bugfixes are backported and added to the `maintenance` branch, new items
are added to the list of bugfixes, if needed with reference to the pull
requests they were taken from.

When a sufficient number of bugfixes has accumulated or an urgent or
important bugfix needs to be distributed a new stable update release is
made.  Since the 'stable' branch is protected, we need to submit a pull
request merging 'maintenance' into 'stable'.  Before merging, set
LAMMPS\_UPDATE in "src/version.h" to "Update #" with "#" indicating the
update count (1, 2, and so on).  Also draft suitable release notes under
https://github.com/lammps/lammps/releases. These can usually be
transferred from the pull request description.

Before that step, several checks need to be performed and changes applied
to the 'maintenance' branch:
- LAMMPS\_UPDATE define in "src/version.h" is changed to "Update #"
- `make -C src check-version` must not list any TBD entries; remaining
  ones are renamed to <stable version>-backport
  (e.g. `.. versionchanged:: 30Sep2026-backport`)
- ensure that GIT\_TAG for the LAMMPS-GUI external target points to the
  LAMMPS-GUI release used for the stable release or a newer bug-fix
  release of it (e.g. v3.1.2 for LAMMPS version 30Sep2026 Update 1).
- ensure that `cmake/packaging/org.lammps.lammps-gui.yml` uses
  `branch: stable` (the release script sets it anyway, but manual
  flatpak builds use the committed value)

If some last minute changes or fixes turn up before the release is
complete, the release process is aborted and the release tag is deleted (it is
protected only after the release is final).

### Merge pull request of 'maintenance' into 'stable', apply tag

The release process begins with merging the pull request, after the
sufficient approval was obtained.  In the local worktree for the
maintenance branch, update both 'maintenance' and 'stable', apply the
stable\_DMmmYYYY\_update# tag to the head of 'stable', switch to
'maintenance' and merge 'stable' (should be a fast-forward) and revert
the LAMMPS\_UPDATE define in "src/version.h" back to "Maintenance".
Commit and push 'maintenance', 'stable' and the tag to the github remote and
our fork.  Example:

```
git checkout stable
git pull
git tag -s -m "Update 1 for Stable LAMMPS version 30 September 2026" stable_30Sep2026_update1
git verify-tag stable_30Sep2026_update1
git checkout maintenance
git pull
git merge stable
# set LAMMPS_UPDATE in src/version.h back to "Maintenance"
git commit -m 'flag maintenance branch with "Maintenance" again' src/version.h
git push git@github.com:lammps/lammps.git maintenance stable stable_30Sep2026_update1
```

Associate draft release notes with the new tag, select "latest release",
and save the draft.  Saving as a draft is required because we need
to use the tag for uploading release assets with `gh release upload`.
The `git verify-tag` command should show that the tag was properly
GPG signed.

### Build and upload binary packages and source tarball to GitHub

The build procedure is the same as for the feature releases, only
that packages are built from the 'stable' branch.  This includes
creating and uploading the signed list of checksums as the last step.
The download server and the website are updated as for a feature
release.

### Publish the release

Do not publish until everything is complete and the SHA256SUMS with a
corresponding GPG signature using the release GPG key have been
uploaded.  After the release is published, it becomes "immutable" and
the release tag protected.  So it can only be updated by withdrawing it
and creating a new release using a different version tag.

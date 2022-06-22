# Apptainer (aka Singularity) container definitions for compiling/testing LAMMPS

The *.def files in this folder can be used to build container images
for [Apptainer](https://apptainer.org) (previously called
[Singularity](https://sylabs.io)), suitable for compiling and testing
LAMMPS on a variety of OS variants with support for most standard
packages and - for some of them - also building/spellchecking the manual
in all supported formats.  This allows to test and debug LAMMPS code on
different OS variants without doing a full installation on your development
workstation, e.g. when bugs are reported that can only be reproduced on
a specific OS or with specific (mostly older) versions of tools,
compilers, or libraries.

Here is a workflow for testing a compilation of LAMMPS with a locally
built CentOS 7.x Singularity container.  For Apptainer replace the
`singularity` command with `apptainer`.

```
cd some/work/directory
git clone --depth 500  https://github.com/lammps/lammps.git lammps
mkdir build-centos7
cd build-centos7
sudo singularity build centos7.sif ../tools/singularity/centos7.def
singularity exec centos7.sif bash --login
cmake -C ../cmake/presets/most.cmake ../cmake
make
```


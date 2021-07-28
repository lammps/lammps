# Singularity container definitions for compiling/testing LAMMPS

The *.def files in this folder can be used to build container images
for [Singularity](https://sylabs.io), suitable for compiling and testing
LAMMPS on a variety of OS variants with support for most standard
packages and - for some of them - also building/spellchecking the manual
in all supported formats. This allows to test and debug LAMMPS code on
different OS variants without doing a full installation on your development
workstation, e.g. when bugs are reported that can only be reproduced on
a specific OS or with specific (mostly older) versions of tools,
compilers, or libraries.

Ready-to-use container images built from some these definition files are
occasionally uploaded to the container library at sylabs.io. They
can be found here: https://cloud.sylabs.io/library/lammps/default/lammps_development#
and will be signed with a GPG key that has the fingerprint:
EEA103764C6C633EDC8AC428D9B44E93BF0C375A

Here is a workflow for testing a compilation of LAMMPS with a locally
built CentOS 7.x singularity container.

```
cd some/work/directory
git clone --depth 500  git://github.com/lammps/lammps.git lammps
mkdir build-centos7
cd build-centos7
sudo singularity build centos7.sif ../tools/singularity/centos7.def
singularity exec centos7.sif bash --login
cmake -C ../cmake/presets/most.cmake ../cmake
make
```

And here is the equivalent workflow for testing a compilation of LAMMPS
using a pre-built Ubuntu 18.04LTS singularity container.

```
cd some/work/directory
git clone --depth 500  git://github.com/lammps/lammps.git lammps
mkdir build-ubuntu18
cd build-ubuntu18
singularity pull library://lammps/default/lammps_development:ubuntu18.04
singularity exec lammps_development_ubuntu18.04.sif bash --login
cmake -C ../cmake/presets/most.cmake ../cmake
make
```

| Currently available:           | Description                                    |
| ------------------------------ | ---------------------------------------------- |
| centos7.def                    | CentOS 7.x with EPEL enabled, no LaTeX         |
| centos8.def                    | CentOS 8.x with EPEL enabled                   |
| fedora34_mingw.def             | Fedora 34 with MinGW cross-compiler toolchain  |
| ubuntu16.04.def                | Ubuntu 16.04LTS with MPI == OpenMPI, no LaTeX  |
| ubuntu18.04.def                | Ubuntu 18.04LTS with MPI == OpenMPI            |
| ubuntu18.04_amd_rocm.def       | Ubuntu 18.04LTS with AMD ROCm toolkit          |
| ubuntu18.04_gpu.def            | Ubuntu 18.04LTS with -"- plus Nvidia CUDA 11.0 |
| ubuntu18.04_nvidia.def         | Ubuntu 18.04LTS with Nvidia CUDA 11.0 toolkit  |
| ubuntu18.04_intel_opencl.def   | Ubuntu 18.04LTS with Intel OpenCL runtime      |
| ubuntu20.04.def                | Ubuntu 20.04LTS with MPI == OpenMPI            |
| ubuntu20.04_amd_rocm.def       | Ubuntu 20.04LTS with AMD ROCm toolkit          |
| ubuntu20.04_gpu.def            | Ubuntu 20.04LTS with -"- plus Nvidia CUDA 11.0 |
| ubuntu20.04_nvidia.def         | Ubuntu 20.04LTS with Nvidia CUDA 11.0 toolkit  |
| ubuntu20.04_intel_opencl.def   | Ubuntu 20.04LTS with Intel OpenCL runtime      |
| ------------------------------ | ---------------------------------------------- |

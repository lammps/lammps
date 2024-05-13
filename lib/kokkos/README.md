![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Kokkos: Core Libraries

Kokkos Core implements a programming model in C++ for writing performance portable
applications targeting all major HPC platforms. For that purpose it provides
abstractions for both parallel execution of code and data management.
Kokkos is designed to target complex node architectures with N-level memory
hierarchies and multiple types of execution resources. It currently can use
CUDA, HIP, SYCL, HPX, OpenMP and C++ threads as backend programming models with several other
backends in development.

**Kokkos Core is part of the Kokkos C++ Performance Portability Programming EcoSystem.**

For the complete documentation, click below:

# [kokkos.github.io/kokkos-core-wiki](https://kokkos.github.io/kokkos-core-wiki)

# Learning about Kokkos

To start learning about Kokkos:

- [Kokkos Lectures](https://kokkos.github.io/kokkos-core-wiki/videolectures.html): they contain a mix of lecture videos and hands-on exercises covering all the important Kokkos Ecosystem capabilities.

- [Programming guide](https://kokkos.github.io/kokkos-core-wiki/programmingguide.html): contains in "narrative" form a technical description of the programming model, machine model, and the main building blocks like the Views and parallel dispatch.

- [API reference](https://kokkos.github.io/kokkos-core-wiki/): organized by category, i.e., [core](https://kokkos.github.io/kokkos-core-wiki/API/core-index.html), [algorithms](https://kokkos.github.io/kokkos-core-wiki/API/algorithms-index.html) and [containers](https://kokkos.github.io/kokkos-core-wiki/API/containers-index.html) or, if you prefer, in [alphabetical order](https://kokkos.github.io/kokkos-core-wiki/API/alphabetical.html).

- [Use cases and Examples](https://kokkos.github.io/kokkos-core-wiki/usecases.html): a series of examples ranging from how to use Kokkos with MPI to Fortran interoperability.

For questions find us on Slack: https://kokkosteam.slack.com or open a GitHub issue.

For non-public questions send an email to: *crtrott(at)sandia.gov*

# Contributing to Kokkos

Please see [this page](https://kokkos.github.io/kokkos-core-wiki/contributing.html) for details on how to contribute.

# Requirements, Building and Installing

All requirements including minimum and primary tested compiler versions can be found [here](https://kokkos.github.io/kokkos-core-wiki/requirements.html).

Building and installation instructions are described [here](https://kokkos.github.io/kokkos-core-wiki/building.html).

# Citing Kokkos

Please see the [following page](https://kokkos.github.io/kokkos-core-wiki/citation.html).

# License

[![License](https://img.shields.io/badge/License-Apache--2.0_WITH_LLVM--exception-blue)](https://spdx.org/licenses/LLVM-exception.html)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.

The full license statement used in all headers is available [here](https://kokkos.org/kokkos-core-wiki/license.html) or
[here](https://github.com/kokkos/kokkos/blob/develop/LICENSE).

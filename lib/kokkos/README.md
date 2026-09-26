
<div align="center">

[<img src="https://github.com/kokkos/kokkos.github.io/blob/main/assets/img/kokkos-logo.png" width="50%">](https://kokkos.org)

</div>

<div align="center">

[<img src="https://github.com/hpsfoundation/hpsf-logos/blob/main/Logos/PNG/Horizontal/HPSF_horizontal-tagline-color.png"
  width="23%" style="vertical-align: middle;">](https://hpsf.io)&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
[<img src="https://github.com/hpsfoundation/hpsf-logos/blob/main/Badges/HPSF_Project_Badge_Established.png?raw=true"
  width="8%" style="vertical-align: middle;">](https://hpsf.io/projects/#tab-established)&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;
[<img src="https://www.linuxfoundation.org/hubfs/lf-stacked-color.svg"
  width="23%" style="vertical-align: middle;">](https://linuxfoundation.org)

</div>

[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/9344/badge)](https://www.bestpractices.dev/projects/9344)

# Kokkos: Core Libraries

Kokkos Core implements a programming model in C++ for writing performance portable
applications targeting all major HPC platforms. For that purpose it provides
abstractions for both parallel execution of code and data management.
Kokkos is designed to target complex node architectures with N-level memory
hierarchies and multiple types of execution resources. It currently can use
CUDA, HIP, SYCL, HPX, OpenMP and C++ threads as backend programming models with several other
backends in development.

**Kokkos Core is part of the [Kokkos C++ Performance Portability Programming Ecosystem](https://kokkos.org).**

Kokkos is a [Linux Foundation](https://linuxfoundation.org) project.

## Learning about Kokkos

To start learning about Kokkos:

- [Kokkos Lectures](https://kokkos.org/kokkos-core-wiki/tutorials-and-examples/video-lectures.html): they contain a mix of lecture videos and hands-on exercises covering all the important capabilities.

- [Programming guide](https://kokkos.org/kokkos-core-wiki/programmingguide.html): contains in "narrative" form a technical description of the programming model, machine model, and the main building blocks like the Views and parallel dispatch.

- [API reference](https://kokkos.org/kokkos-core-wiki/): organized by category, i.e., [core](https://kokkos.org/kokkos-core-wiki/API/core-index.html), [algorithms](https://kokkos.org/kokkos-core-wiki/API/algorithms-index.html), [containers](https://kokkos.org/kokkos-core-wiki/API/containers-index.html), and [simd](https://kokkos.org/kokkos-core-wiki/API/simd-index.html).

- [Use cases and Examples](https://kokkos.org/kokkos-core-wiki/tutorials-and-examples/use-cases-and-examples.html): a serie of examples ranging from how to use Kokkos with MPI to Fortran interoperability.

## Obtaining Kokkos

The latest release of Kokkos can be obtained from the [GitHub releases page](https://github.com/kokkos/kokkos/releases/latest).

The current release is [5.2.0](https://github.com/kokkos/kokkos/releases/tag/5.2.0).

```bash
curl -OJ -L https://github.com/kokkos/kokkos/releases/download/5.2.0/kokkos-5.2.0.tar.gz
# Or with wget
wget https://github.com/kokkos/kokkos/releases/download/5.2.0/kokkos-5.2.0.tar.gz
# Or with git
git clone --depth=2 --branch 5.2.0 https://github.com/kokkos/kokkos.git
```

To clone the latest development version of Kokkos from GitHub:

```bash
git clone --branch develop  https://github.com/kokkos/kokkos.git
```

### Building Kokkos

To build Kokkos, you will need to have a C++ compiler that supports C++20 or later.
All requirements including minimum and primary tested compiler versions can be found [here](https://kokkos.org/kokkos-core-wiki/get-started/requirements.html).

Building and installation instructions are described [here](https://kokkos.org/kokkos-core-wiki/get-started/building-from-source.html#configuring-and-building-kokkos).

You can also install Kokkos using [Spack](https://spack.io/): `spack install kokkos`. [Available configuration options](https://packages.spack.io/package.html?name=kokkos) can be displayed using `spack info kokkos`.

## For the complete documentation: [kokkos.org/kokkos-core-wiki/](https://kokkos.org/kokkos-core-wiki/)

## Support

For questions find us on Slack: https://kokkosteam.slack.com or open a GitHub issue.

For non-public questions send an email to: *crtrott(at)sandia.gov*

## Contributing

Please see [this page](https://kokkos.org/kokkos-core-wiki/contributing.html) for details on how to contribute.

## Citing Kokkos

Please see the [following page](https://kokkos.org/kokkos-core-wiki/citation.html).

## License

[![License](https://img.shields.io/badge/License-Apache--2.0_WITH_LLVM--exception-blue)](https://spdx.org/licenses/LLVM-exception.html)

The full license statement used in all headers is available [here](https://kokkos.org/kokkos-core-wiki/license.html) or
[here](https://github.com/kokkos/kokkos/blob/develop/LICENSE).

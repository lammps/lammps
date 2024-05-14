#!/usr/bin/env python
"""
Install.py tool to download, compile, and setup the ml-metatensor LAMMPS package.
This automates the steps described in the README file in this dir.
"""

import glob
import os
import shutil
import subprocess
import sys
from argparse import ArgumentParser

sys.path.append("..")
from install_helpers import checkmd5sum, fullpath, geturl  # noqa F401

if sys.version_info < (3, 6):
    sys.exit("this script requires at least Python 3.6")

# settings

HERE = fullpath(".")
METATENSOR_CORE_VERSION_DEFAULT = "0.1.7"
METATENSOR_TORCH_VERSION_DEFAULT = "0.4.0"
LIBTORCH_VERSION_DEFAULT = "2.2.2"

# known checksums for different versions. used to validate the download.
METATENSOR_CORE_CHECKSUMS = {
    "0.1.6": "3edaf5ffd1af8892965f1f618c04083d",
    "0.1.7": "43aeb651e6d040fbc82d4b5191e2f6cb",
    "0.1.8": "22ae27fb3b26f356ca6e477326be6470",
}
METATENSOR_TORCH_CHECKSUMS = {
    "0.4.0": "d10ddfaf0213ec351a4386aae7d89dd0",
}

GITHUB_RELEASES = "https://github.com/lab-cosmo/metatensor/releases/download"

CMAKE_EXE = os.environ.get("CMAKE", "cmake")


def download_unpack(url, unpacked_dir, checksum=None, verbose=False):
    file_name = os.path.basename(url)
    print(f"===> Downloading {file_name}")

    geturl(url, os.path.basename(url))
    if checksum is not None:
        if not checkmd5sum(checksum, file_name):
            sys.exit(f"ERROR: checksum for {file_name} does not match!")

    print(f"===> Unpacking to {unpacked_dir}")
    if verbose:
        stdout = sys.stdout
        stderr = sys.stderr
    else:
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE

    shutil.rmtree(unpacked_dir, ignore_errors=True)
    subprocess.run(
        [CMAKE_EXE, "-E", "tar", "xf", file_name],
        stdout=stdout,
        stderr=stderr,
        check=True,
    )
    extracted = file_name[:-7]
    shutil.move(extracted, unpacked_dir)


def install_with_pip(python, package, verbose=False):
    print(f"===> Installing {package} with pip")
    if verbose:
        stdout = sys.stdout
        stderr = sys.stderr
    else:
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE

    subprocess.run(
        [
            python,
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            package,
        ],
        stdout=stdout,
        stderr=stderr,
        check=True,
    )


def copy_from_pip(python, package, paths, install_dir, verbose=False):
    print(f"===> Copying files for {package} from Python installation")
    prefix = subprocess.check_output(
        [python, "-c", f"import {package}; print({package}.__file__)"],
        encoding="utf8",
    )
    prefix = os.path.dirname(prefix.strip())

    for src, dst in paths:
        shutil.copytree(
            os.path.join(prefix, src),
            os.path.join(install_dir, dst),
            dirs_exist_ok=True,
        )


def build_cmake(source_dir, install_dir, cmake_opts, verbose=False):
    print(f"===> Building {source_dir} with cmake")
    if verbose:
        stdout = sys.stdout
        stderr = sys.stderr
    else:
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE

    build_dir = os.path.join(source_dir, "build")

    subprocess.run(
        [
            CMAKE_EXE,
            "-S",
            source_dir,
            "-B",
            build_dir,
            *cmake_opts,
            f"-DCMAKE_INSTALL_PREFIX={install_dir}",
            "-DCMAKE_BUILD_TYPE=Release",
        ],
        stdout=stdout,
        stderr=stderr,
        check=True,
    )

    subprocess.run(
        [
            CMAKE_EXE,
            "--build",
            build_dir,
            "--config",
            "Release",
            "--target",
            "install",
        ],
        stdout=stdout,
        stderr=stderr,
        check=True,
    )


if __name__ == "__main__":

    parser = ArgumentParser(
        prog="Install.py", description="LAMMPS library build wrapper script"
    )

    # help message

    HELP = """
This script tries to install all the dependencies of the ml-metatensor package. This
includes libtorch, metatensor and metatensor-torch. Different options are available for
the different dependencies:

- libtorch can be downloaded from pip, or taken from another installation on your
  system.
    - If downloaded by pip, we will create a virtual environment and install it there,
      unless `--python` is given as an option, in which case we will try to install it
      using the provided python executable.
    - If you want to use another installation of libtorch, please use the `--no-torch`
      option and export the TORCH_PREFIX environment variable containing the path of the
      installation.
- metatensor and metatensor-torch can be downloaded from pip or built from sources
    - if downloaded by pip, the same options as torch applies. This is triggered by the
      `--metatensor-use-pip` option.
    - if building from sources, you will need to install cmake and a rust compiler on
      your system (we suggest https://rustup.rs/ to install a rust compiler).


Syntax from src dir: make lib-metatensor args="-b ..."
Syntax from lib dir: python Install.py -b ...

Examples:

# install with default versions and settings
make lib-metatensor args="-b"

# install specified version of libtorch
make lib-metatensor args="-b --torch-version <version>"

# install using python3 in the PATH
make lib-metatensor args="-b --python $(which python3)"

# build metatensor from sources
make lib-metatensor args="-b --metatensor-from-sources"
    """

    parser.add_argument(
        "-b",
        "--build",
        action="store_true",
        help="download and build metatensor and metatensor-torch libraries",
    )
    parser.add_argument(
        "--torch-version",
        default=LIBTORCH_VERSION_DEFAULT,
        help=f"version of libtorch to download (default: {LIBTORCH_VERSION_DEFAULT})",
    )
    parser.add_argument(
        "--python",
        help="path to a Python executable to use, this bypass the use of venv",
    )

    parser.add_argument(
        "--no-torch",
        action="store_true",
        default=False,
        help="disabled download of libtorch",
    )
    parser.add_argument(
        "--metatensor-from-sources",
        action="store_true",
        help="build metatensor and metatensor-torch from sources",
    )

    parser.add_argument(
        "--metatensor-version",
        default=METATENSOR_CORE_VERSION_DEFAULT,
        choices=METATENSOR_CORE_CHECKSUMS.keys(),
        help=(
            "version of metatensor-core to download and build "
            f"(default: {METATENSOR_CORE_VERSION_DEFAULT})"
        ),
    )
    parser.add_argument(
        "--metatensor-torch-version",
        default=METATENSOR_TORCH_VERSION_DEFAULT,
        choices=METATENSOR_TORCH_CHECKSUMS.keys(),
        help=(
            "version of metatensor-torch to download and build "
            f"(default: {METATENSOR_TORCH_VERSION_DEFAULT})"
        ),
    )
    parser.add_argument(
        "-vv",
        "--verbose",
        action="store_true",
        help="be more verbose about is happening while this script runs",
    )

    args = parser.parse_args()

    # print help message and exit, if build option is not given
    if not args.build:
        parser.print_help()
        sys.exit(HELP)

    do_build = args.build
    verbose = args.verbose

    metatensor_core_version = args.metatensor_version
    metatensor_torch_version = args.metatensor_torch_version
    metatensor_from_sources = args.metatensor_from_sources

    do_torch = not args.no_torch
    torch_version = args.torch_version

    # If we don't know the version of torch, we need to build metatensor-torch from
    # sources.
    if args.no_torch and not metatensor_from_sources:
        raise Exception("--no-torch requires --metatensor-from-sources")

    python = args.python

    if python is None:
        python = sys.executable
        do_venv = do_torch or not metatensor_from_sources
    else:
        do_venv = False

    if verbose:
        stdout = sys.stdout
        stderr = sys.stderr
    else:
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE

    if do_venv:
        venv_path = os.path.join(HERE, "virtualenv")
        print(f"===> Creating virtual environment at {venv_path}")

        try:
            import venv
        except ImportError:
            sys.exit(
                "could not import 'venv', make sure 'python-venv' package is installed"
            )

        builder = venv.EnvBuilder(
            with_pip=True,
            symlinks=not sys.platform.startswith("win"),
        )
        builder.create(venv_path)

        if sys.platform.startswith("win"):
            python = os.path.join(venv_path, "Scripts", "python.exe")
        else:
            python = os.path.join(venv_path, "bin", "python")

        if verbose:
            print("===> Upgrading pip in virtualenv")

        subprocess.run(
            [python, "-m", "pip", "install", "--upgrade", "pip"],
            stdout=stdout,
            stderr=stderr,
            check=True,
        )

    install_dir = os.path.join(HERE, "usr")
    os.makedirs(os.path.join(install_dir, "lib"), exist_ok=True)
    os.makedirs(os.path.join(install_dir, "include"), exist_ok=True)
    os.makedirs(os.path.join(install_dir, "share"), exist_ok=True)

    cmake_prefix_path = []
    # Torch needs C++17 to compile
    extra_cxx_flags = ["-std=c++17"]

    torch_extra_include = os.path.join("torch", "csrc", "api", "include")
    if do_torch:
        install_with_pip(python, f"torch=={torch_version}", verbose=verbose)
        install_with_pip(python, "numpy", verbose=verbose)

        torch_prefix = subprocess.check_output(
            [python, "-c", "import torch; print(torch.__file__)"],
            encoding="utf8",
        )
        torch_prefix = os.path.dirname(torch_prefix.strip())

        if sys.platform.startswith("linux"):
            # the pip version of torch uses the pre-cxx11 ABI
            extra_cxx_flags.append("-D_GLIBCXX_USE_CXX11_ABI=0")

    else:
        torch_prefix = os.environ.get("TORCH_PREFIX")

    if torch_prefix is not None:
        cmake_prefix_path.append(torch_prefix)

        if os.path.exists(os.path.join(torch_prefix, "include")):
            extra_cxx_flags.append(f"-I{os.path.join(torch_prefix, 'include')}")
            extra_cxx_flags.append(
                f"-I{os.path.join(torch_prefix, 'include', torch_extra_include)}"
            )

    if metatensor_from_sources:
        # Download and build metatensor-core from source
        download_unpack(
            url=f"{GITHUB_RELEASES}/metatensor-core-v{metatensor_core_version}/metatensor-core-cxx-{metatensor_core_version}.tar.gz",  # noqa: E501
            unpacked_dir="metatensor-core",
            checksum=METATENSOR_CORE_CHECKSUMS.get(metatensor_core_version),
            verbose=verbose,
        )

        build_cmake(
            source_dir="metatensor-core",
            install_dir=install_dir,
            cmake_opts=[
                "-DBUILD_SHARED_LIBS=ON",
                "-DMETATENSOR_INSTALL_BOTH_STATIC_SHARED=OFF",
            ],
            verbose=verbose,
        )

        # Download and build metatensor-torch from source
        cmake_prefix_path.append(f"{install_dir}/lib/")
        download_unpack(
            url=f"{GITHUB_RELEASES}/metatensor-torch-v{metatensor_torch_version}/metatensor-torch-cxx-{metatensor_torch_version}.tar.gz",  # noqa: E501
            unpacked_dir="metatensor-torch",
            checksum=METATENSOR_TORCH_CHECKSUMS.get(metatensor_torch_version),
            verbose=verbose,
        )

        build_cmake(
            source_dir="metatensor-torch",
            install_dir=install_dir,
            cmake_opts=[
                "-DBUILD_SHARED_LIBS=ON",
                f"-DCMAKE_PREFIX_PATH={';'.join(cmake_prefix_path)}",
            ],
            verbose=verbose,
        )
    else:
        install_with_pip(
            python,
            f"metatensor-core=={metatensor_core_version}",
            verbose=verbose,
        )
        copy_from_pip(
            python,
            "metatensor",
            [("lib", "lib"), ("include", "include")],
            install_dir=install_dir,
            verbose=verbose,
        )

        install_with_pip(
            python,
            f"metatensor-torch=={metatensor_torch_version}",
            verbose=verbose,
        )
        torch_version_prefix = "torch-" + ".".join(torch_version.split(".")[:2])
        copy_from_pip(
            python,
            "metatensor.torch",
            [
                (os.path.join(torch_version_prefix, "lib"), "lib"),
                (os.path.join(torch_version_prefix, "include"), "include"),
            ],
            install_dir=install_dir,
            verbose=verbose,
        )

    print("===> Creating Makefile.lammps")
    with open("Makefile.lammps", "w") as fd:
        fd.write("# autogenerated file\n\n\n")

        fd.write(f"metatensor_SYSINC = -I{os.path.join(install_dir, 'include')}")
        fd.write(f" {' '.join(extra_cxx_flags)}")
        fd.write("\n\n")

        fd.write("metatensor_SYSLIB = -lmetatensor -lmetatensor_torch")
        fd.write(" -ltorch -lc10 -ltorch_cpu")
        if len(glob.glob(os.path.join(torch_prefix, "lib", "*torch_cuda*"))) != 0:
            fd.write(" -ltorch_cuda -lc10_cuda")
        fd.write("\n\n")

        fd.write(f"metatensor_SYSPATH = -L{os.path.join(install_dir, 'lib')}")
        fd.write(f" -L{os.path.join(torch_prefix, 'lib')}")

        # set the rpath on the final binary
        fd.write(f" -Wl,-rpath,{os.path.join(install_dir, 'lib')}")
        fd.write(f" -Wl,-rpath,{os.path.join(torch_prefix, 'lib')}")

        if sys.platform.startswith("linux"):
            # force the use of rpath instead of runpath
            fd.write(" -Wl,--disable-new-dtags")

        fd.write("\n\n")

    print("===> All done!")

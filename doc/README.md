# Generation of LAMMPS Documentation

The generation of all the documentation is managed by the Makefile inside the
`doc/` folder.

## Usage:

```bash
make html         # generate HTML using Sphinx
make pdf          # generate PDF using htmldoc
make clean        # remove generated RST files
make clean-all    # remove entire build folder and any cached data
```

## Installing prerequisites

To run the documention build toolchain Python 3 and virtualenv have to be
installed. Here are instructions for common setups:

### Ubuntu

```bash
sudo apt-get install python-virtualenv
```

### Fedora

```
sudo yum install python-virtualenv
```

### MacOS X

## Python 3

Download the latest Python 3 MacOS X package from https://www.python.org and install it.
This will install both Python 3 and pip3.

## virtualenv

Once Python 3 is installed, open a Terminal and type `pip3 install virtualenv`. This will
install virtualenv from the Python Package Index.

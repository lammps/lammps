# LAMMPS Documentation Utilities

[![Build Status](https://travis-ci.org/rbberger/lammps-doc-utils.svg?branch=master)](https://travis-ci.org/rbberger/lammps-doc-utils)

[![Coverage Status](https://coveralls.io/repos/github/rbberger/lammps-doc-utils/badge.svg?branch=master)](https://coveralls.io/github/rbberger/lammps-doc-utils?branch=master)

This repository contains a set of utilities to convert existing LAMMPS
documentation text files into ReStructured Text. These files can then be used
to generate documentation using Sphinx (www.sphinx-doc.org).

The goal of these tools is to simplify the transition to this new format by
automatically transforming existing formatting and adding some semantic
replacements.

## Contents

* A Python port of txt2html
* A new txt2rst utility which transforms txt files directly into ReStructured
  Text for Sphinx
* Unit tests

## Prerequisites

* Currently these tools only run on Python 3

## Sphinx requirements

The generated ReStructured Text assumes you've installed both `sphinx` and
`sphinxcontrib-images`. Both can be installed through pip:

```
pip3 install sphinx
pip3 install sphinxcontrib-images
```

Once Sphinx is set up adjust your configuration `conf.py` to use the
`sphinxcontrib-images` extension:

```python
extensions = [
              ...
              'sphinxcontrib.images',
              ...
]
```

## Installation

1. Clone this repository
2. Install using setup.py

   ```bash
   python setup.py install
   ```

## Usage

Convert .txt files using `txt2rst`:

```bash
# single files
txt2rst Manual.txt > Manual.rst

# multiple files
txt2rst *.txt
```

## Backwards compatibility with txt2html

### RST portions

Since the original txt2html passes through HTML tags the following syntax can
be used to invalidate RST code:

```html
<!-- RST

.. toctree::
   :maxdepth: 2
   :numbered:
   
   Section_intro
   Section_start
   Section_commands

END_RST -->
```

### HTML_ONLY portions

At the same time, some contents of the original txt files might no longer be
needed when compiling the RST files using Sphinx. These regions can be marked
using a HTML_ONLY region:

```html
<!-- HTML_ONLY -->
<HEAD>
<TITLE>LAMMPS Users Manual</TITLE>
...
</HEAD>

<BODY>
<!-- END_HTML_ONLY -->
```

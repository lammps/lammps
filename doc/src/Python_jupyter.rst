Using LAMMPS in IPython notebooks and Jupyter
=============================================

If the LAMMPS Python package is installed for the same Python interpreter as
`IPython <ipython>`_, you can use LAMMPS directly inside of an IPython notebook inside of
Jupyter. `Jupyter <juypter>`_ is a powerful integrated development environment (IDE) for
many dynamic languages like Python, Julia and others, which operates inside of
any web browser. Besides auto-completion and syntax highlighting it allows you
to create formatted documents using Markup, mathematical formulas, graphics and
animations intermixed with executable Python code. It is a great format for
tutorials and showcasing your latest research.

The easiest way to install it is via ``pip``:

.. code-block:: bash

   pip install --user jupyter

To launch an instance of Jupyter simply run the following command inside your
Python environment:

.. code-block:: bash

   jupyter notebook

Interactive Python Examples
---------------------------

Examples of IPython notebooks can be found in the ``python/examples/ipython``
subdirectory. They require LAMMPS to be compiled as shared library with PYTHON,
PNG, JPEG and FFMPEG support.

To open these notebooks launch ``jupyter notebook index.ipynb`` inside this
directory. The opened file provides an overview of the available examples.

- Example 1: Using LAMMPS with Python (``simple.ipynb``)
- Example 2: Analyzing LAMMPS thermodynamic data (``thermo.ipynb``)
- Example 3: Working with Per-Atom Data (``atoms.ipynb``)
- Example 4: Working with LAMMPS variables (``variables.ipynb``)
- Example 5: Validating a dihedral potential (``dihedrals/dihedral.ipynb``)
- Example 6: Running a Monte Carlo relaxation (``montecarlo/mc.ipynb``)

.. note::

   Typically clicking a link in Jupyter will open a new tab, which might be blocked by your pop-up blocker.

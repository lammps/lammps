from setuptools import setup

setup(name='LAMMPS Documentation Utilities',
      version='2.0.0',
      description='Utilities to convert existing LAMMPS documentation text files into ReStructured Text',
      url='https://github.com/rbberger/lammps-doc-utils',
      author='Richard Berger',
      author_email='richard.berger@outlook.com',
      license='GPL',
      packages=['lammpsdoc'],
      test_suite='nose.collector',
      tests_require=['nose'],
      entry_points = {
          "console_scripts": ['txt2html = lammpsdoc.txt2html:main',
                              'txt2rst  = lammpsdoc.txt2rst:main',
                              'rst_anchor_check = lammpsdoc.rst_anchor_check:main ']
      },
)

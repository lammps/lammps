from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# To build do: python setup.py build_ext --inplace

ext_modules = [Extension("lmp2radii",["lmp2radii.pyx"])]
setup(
    name = 'lmp2radii',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)

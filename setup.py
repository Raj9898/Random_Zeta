from distutils.core import setup
from Cython.Build import cythonize

setup(
    name='Zeta Project',
    ext_modules=cythonize("zetaFunctions.pyx"),
    zip_safe=False,
)

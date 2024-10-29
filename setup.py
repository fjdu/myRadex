from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize
from Cython.Compiler import Options
Options.docstrings = True

extension1 = Extension(
    name="myRadex",
    sources=["wrapper_for_cython.pyx"],
    libraries=["my_radex", "gfortran"],
    library_dirs=["./", '/usr/local/Cellar/gcc/14.2.0/lib/gcc/current/'],
    include_dirs=["./"],
    depends=['setup.py', 'makefile'],
    extra_compile_args=["-std=c++11"],
    extra_link_args=["-mmacosx-version-min=13.0"]
)
# gfortran --print-file-name libgfortran.a

setup(
    name="myRadex",
    version='0.1',
    description='A radex wrapper',
    author='Fujun Du',
    author_email='fujun.du@gmail.com',
    ext_modules=cythonize([extension1], language_level=3, gdb_debug=True, compiler_directives={'linetrace': True})
)

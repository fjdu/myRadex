from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize
from Cython.Compiler import Options
Options.docstrings = True

import subprocess
import os

def find_gfortran_lib():
    """Find gfortran library directory using gfortran itself"""
    lib_dirs = []
    
    # Method 1: Find libgfortran directly
    try:
        result = subprocess.check_output(
            ['gfortran', '--print-file-name=libgfortran.so'],
            stderr=subprocess.DEVNULL,
            universal_newlines=True
        ).strip()
        
        if result and result != 'libgfortran.so':  # Found actual path, not just filename
            lib_dir = os.path.dirname(result)
            lib_dirs.append(lib_dir)
            print(f"Found gfortran library at: {lib_dir}")
            
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("gfortran not found or failed to query")
    
    # Method 2: Also try libgfortran.a for static linking
    try:
        result = subprocess.check_output(
            ['gfortran', '--print-file-name=libgfortran.a'],
            stderr=subprocess.DEVNULL,
            universal_newlines=True
        ).strip()
        
        if result and result != 'libgfortran.a':
            lib_dir = os.path.dirname(result)
            if lib_dir not in lib_dirs:
                lib_dirs.append(lib_dir)
                print(f"Found gfortran static library at: {lib_dir}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    
    return lib_dirs


extension1 = Extension(
    name="myRadex",
    sources=["wrapper_for_cython.pyx"],
    libraries=["my_radex", "gfortran"],
    library_dirs=["./"] + find_gfortran_lib(),
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

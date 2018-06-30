from __future__ import print_function
import sys
from warnings import warn
from platform import system
from setuptools import setup, Extension
from Cython.Distutils import build_ext

persUF = Extension(name="persUF",
                   sources=["persUF.pyx"],
                   language="c++",
                   extra_compile_args=["-Wmaybe-uninitialized", "-Wunused-but-set-variable", "-lboost_filesystem", "-std=c++11"])
                   extra_link_args=["-std=c++11"])

setup(
    name="persUF",
    author="Mathieu Carriere",
    author_email="mathieu.carriere3@gmail.com",
    ext_modules = [persUF],
    cmdclass = {"build_ext": build_ext},
)

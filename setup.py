from setuptools       import setup, Extension
from Cython.Build     import cythonize
from Cython.Distutils import build_ext

persUF = Extension(name="persUF",
                   sources=["persUF.pyx"],
                   language="c++",
                   extra_compile_args=["-Wmaybe-uninitialized", "-Wunused-but-set-variable", "-lboost_filesystem", "-std=c++11"])

setup(
    name="persUF",
    author="Mathieu Carriere",
    author_email="mathieu.carriere3@gmail.com",
    ext_modules = cythonize(persUF),
    cmdclass = {"build_ext": build_ext},
)

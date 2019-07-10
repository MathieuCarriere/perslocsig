from setuptools       import setup, Extension
from Cython.Build     import cythonize
from Cython.Distutils import build_ext

perslocsig = Extension(name="perslocsig",
                       sources=["perslocsig/perslocsig.pyx"],
                       language="c++",
                       extra_compile_args=["-Wmaybe-uninitialized", "-Wunused-but-set-variable", "-lboost_filesystem", "-std=c++11"])

setup(name="perslocsig",
      author="Mathieu Carriere",
      author_email="mathieu.carriere3@gmail.com",
      ext_modules = cythonize([perslocsig]),
      cmdclass = {"build_ext": build_ext})

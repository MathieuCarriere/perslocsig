This code is an implementation of the algorithm described and used in
Stable topological signatures for points on 3D shapes, Carriere et al, SGP 2015
It requires python, Boost and C++11.

**********

persistence_UF.h contains the C++ code to compute 1D persistence using symmetry. It requires a neighborhood graph and a function defined on the vertices of this graph.
To compile the code, run the following command in a terminal:

python setup.py build_ext

This will produce a .so file that you can import in a python shell. Then, you can call geodesic_diagrams_on_shape(name.off, indices), which computes the persistence diagrams 
of the points given by indices on the shape name.off. For instance, calling geodesic_diagrams_on_shape("61.off", [:10]) will compute the diagrams associated to the first 10 points of shape 61.off.

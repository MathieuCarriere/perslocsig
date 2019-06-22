# Local persistence on 3D shapes with Union Find 

This code is an implementation of the algorithm described and used in
Stable topological signatures for points on 3D shapes, Carriere et al, SGP 2015
It requires python, Boost and C++11. To compile the code, run the following command in a terminal:

python setup.py build_ext

This will produce a .so file that you can import in a python shell. 
Then, you can call geodesic_diagrams_on_shape(name.off, indices), which computes the persistence diagrams 
of the points given by indices on the shape name.off. 
For instance, calling geodesic_diagrams_on_shape("61.off", np.arange(10)) will compute the diagrams associated to the first 10 points of shape 61.off.
Check the notebook 3DSigTDA.ipynb for more detailed examples.

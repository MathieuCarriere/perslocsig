# Local persistence on 3D shapes with Union Find 

This code is an implementation of the algorithm described and used in
_Stable topological signatures for points on 3D shapes, Carriere et al, SGP 2015_.
It requires python, Boost and C++11. To compile the code, run the following command in a terminal:

	$ git clone https://github.com/MathieuCarriere/perslocsig
	$ cd perslocsig
	$ (sudo) pip install .

Then, you can import `perslocsig` in python, an especially the function `compute_geodesic_persistence_diagrams`, which has arguments `name.off` and `indices`, 
and which computes the persistence diagrams of the points given by `indices` on the shape `name.off`. 
For instance, running `compute_geodesic_persistence_diagrams("61.off", np.arange(10))` will compute the persistence diagrams associated to the first 10 points of the shape `61.off`.
You can check the notebook `example/perslocsig.ipynb` for more details.

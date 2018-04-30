This code is an implementation of the algorithm described and used in
Stable topological signatures for points on 3D shapes, Carriere et al, SGP 2015
It requires Boost and C++11.

**********

persistence_UF.cpp contains the C++ code to compute 1D persistence using symmetry. It requires a neighborhood graph and a function defined on the vertices of this graph.
To compile the code, run the following command in a terminal:

g++ persistence_UF.cpp -o persUF -std=c++11

Assuming the function is given in a file F.txt as follows:

line 1 <func_vertex_0> <func_vertex_1> <func_vertex_2> <func_vertex_3> ... <func_vertex_n>
line 2 <empty>

and the neighborhood graph is given in a file N.txt as follows:

line 1 <idx_0th_neighbor_of_0> <idx_1st_neighbor_of_0> <idx_2nd_neighbor_of_0> ... <idx_N0th_neighbor_of_0>
line 2 <idx_0th_neighbor_of_1> <idx_1st_neighbor_of_1> <idx_2nd_neighbor_of_1> ... <idx_N1th_neighbor_of_1>
.
.
.
line n <idx_0th_neighbor_of_n> <idx_1st_neighbor_of_n> <idx_2nd_neighbor_of_n> ... <idx_Nnth_neighbor_of_n>

where point 0 has N0 neighbors, point 1 has N1 neighbors and so on, then the persistence diagram of point idx is computed by running the following command in a terminal:

(cat F.txt && cat N.txt) | ./persUF 0

Note that running:

(cat F.txt && cat N.txt) | ./persUF 1

will give you the indices of the pair of vertices that lead to the corresponding point in the persistence diagram.
In the specific case of 3D shapes given in .off file format, we provide graph_geodesic.cpp that computes the geodesic distances to a given point using Boost's Dijkstra algorithm, 
and graph_neighbors.cpp that retrieves the neighbors from the triangulation.
For instance, computing the persistence diagram of point x of the 3D shape given in file.off is done with:

( ./graph_geodesic file.off x && ./graph_neighbors) | ./persUF 0

We also provide compute_PD.sh, which is a bash script that computes the persistence diagrams for all points of all shapes of a given 3D shape category ("airplane", "hand"...).

 

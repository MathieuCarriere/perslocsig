Warning!!! The code has only been tested with g++ on Linux.
After compiling, running the executables in a terminal with no arguments will show you the arguments the executables expect.

**********

persistence_UF contains the C++ code to compute 1D persistence on 3D shapes using symmetry. It needs the "trimesh" library and the "ANN" one (for point clouds only). If you only want to compute PDs on 3D shapes and do not have access to the "ANN" library, comment the C++ code in persistence.cpp that deals with point clouds. To compile the "trimesh" library type the following command in a terminal:

cd trimesh/
mkdir lib.Linux64
cd libsrc
make
cd ../..

To compile the code, run the following command in a terminal.

mkdir bin
cd persistence_UF/
g++ persistence.cpp -o ../bin/persistence -I ../trimesh/include/ -I ../ann/ann_1.1.1/include/ -L ../trimesh/lib.Linux64/ -L ../ann/ann_1.1.1/lib/ -ltrimesh -lANN -fopenmp
cd ../

Parameters for the "persistence" executable that is created in the bin folder are:

1. Type of the object -- "3DShape" or "PointCloud"

2. Name of the object -- e.g.: "1.off"

3. Type of the function that induces the PD -- "loc" for functions anchored at base points (e.g. distance functions, as in the article) or "glo" for global functions

4. Name of the functions -- e.g. "geo" for distance functions, but also "hks", "cur", "ecc"... Warning!!! The code expects that a binary file containing the values of the functions is present in the directory. You may need to modify the code that reads the functions values ("shape.h") according to the shape of your function bianry file.

5. Name of the data structure used to compute the PD -- always type "UF" for Union-Find if you deal with 3D shapes. "filtration" is when you already have executables that can compute PDs with a give filtration file.

6. (optional) delta parameter for the neighborhood graph if you have point clouds. If you only deal with 3D shapes, leave this field blank.

**********

signature.cpp contains the code to compute the topological signature on a given PD. Type the following in a terminal.

g++ signature.cpp -o bin/signature

Parameters for the "signature" executable are:

1. Name of your PD -- e.g.: "PersistenceDiagram_geo_loc_Value_1"

2. Quantity to compute -- "SC" for the signature vectors or "AGD" for the means of the previous vectors

3. Dimension. Needed only if you compute "AGD". Specify the number of components you want to use to compute the mean.
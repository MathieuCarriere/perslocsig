from cython cimport numeric
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

cdef extern from "persistence_UF.h":
    vector[pair[vector[pair[int,int]], vector[pair[double,double]]]] compute_PDs_from_off(string, vector[int])

def geodesic_diagrams_on_shape(name, idx):
    return compute_PDs_from_off(name, idx)

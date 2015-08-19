/* 
 * File:   main.cpp
 * Author: mathieu
 *
 * Created on January 22, 2015, 2:00 PM
 */

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <cmath>
#include <math.h>
#include <omp.h>
 
#include "objects/shape.h"
#include "objects/point_cloud.h"
#include "persistence_diagram/persistence.h"
  
using namespace std;

int main(int argc, char** argv) {
     
    // Read arguments.
    // ***************

    if (argc <= 5){ cout << "Inputs: <object type> <object_name> <loc/glo> <function> <UF/filtration> (<delta>)" << endl; return 0; }	
    
    string object = argv[1];
    string object_name = argv[2];
    string prop = argv[3];
    string function = argv[4];
    string method = argv[5];
    
    // Process.
    // ********
    
    // 3D Shape.
    // *********
    
    if (object.compare("3DShape") == 0){    
    
        size_t len = object_name.size();
        string id = object_name.substr(0,len-4);
        int val = atoi(id.c_str());
        
        trimesh::TriMesh shape = read_3Dshape(object_name);
        
        vector<vector<double> > locfunc;
        vector<double> globfunc;
        vector<vector<pair<int,int> > > local_PD_index;
        vector<vector<pair<double,double> > > local_PD_value;
        vector<pair<int,int> > global_PD_index;
        vector<pair<double,double> > global_PD_value;
        
        // Scalar function.
        // ****************
        
        if (prop.compare("loc") == 0)
            read_3Dshape_local(val, &locfunc, &shape, function);
        
        if (prop.compare("glo") == 0)
            read_3Dshape_global(val, &globfunc, &shape, function);
       
        // Persistence Diagram.
        // ********************
        
        if (prop.compare("loc") == 0){
            if (method.compare("UF") == 0){
                cout << "Warning: UF is not suitable for essential 1D persistence." << endl;
		vector<vector<int> > neighb = shape.neighbors;
                compute_3Dshape_local_PD_with_UF(val, function, prop, &local_PD_index, &local_PD_value, &locfunc, &neighb);
            }
            if (method.compare("filtration") == 0)
                compute_3Dshape_local_PD_with_filtration(val, &locfunc, &shape);
        }
        
        if (prop.compare("glo") == 0){
            if (method.compare("UF") == 0){
                cout << "Warning: UF is not suitable for essential 1D persistence." << endl;
                vector<vector<int> > neighb = shape.neighbors;
                compute_3Dshape_global_PD_with_UF(val, function, prop, &global_PD_index, &global_PD_value, &globfunc, &neighb);
            }
            
            if (method.compare("filtration") == 0)
                compute_3Dshape_global_PD_with_filtration(val, &globfunc, &shape);
        }
        
    }
    
    // Point Cloud.
    // ************
    
    if (object.compare("PointCloud") == 0){

        if (argc == 6){ cout << "Please, provide <delta>" << endl; return 0; }	
        
        double delta = atof(argv[6]);
        vector<vector<int> > neighb; vector<vector<double> > dist;
        build_neighborhood_graph(object_name,delta,&neighb,&dist);

        size_t len = object_name.size();
        string id = object_name.substr(0,len-4);
        int val = atoi(id.c_str());

        vector<vector<double> > locfunc;
        vector<double> globfunc;
        vector<vector<pair<int,int> > > local_PD_index;
        vector<vector<pair<double,double> > > local_PD_value;
        vector<pair<int,int> > global_PD_index;
        vector<pair<double,double> > global_PD_value;

        // Scalar function.
        // ****************

        if (prop.compare("loc") == 0)
            read_pointcloud_local(val, &locfunc, &neighb, &dist, function);

        if (prop.compare("glo") == 0)
            read_pointcloud_global(val, &globfunc, &neighb, &dist, function);

        if (prop.compare("loc") == 0){
            if (method.compare("UF") == 0){
                cout << "Warning: UF is not suitable for essential 1D persistence." << endl;
                compute_3Dshape_local_PD_with_UF(val, function, prop, &local_PD_index, &local_PD_value, &locfunc, &neighb);
            }

            if (method.compare("filtration") == 0){
                cout << "Warning: computes homology of dim <= 2." << endl;
                compute_PtCloud_local_PD_with_filtration(val, &locfunc, &neighb);
            }
        }
        
        if (prop.compare("glo") == 0){
            if (method.compare("UF") == 0){
                cout << "Warning: UF is not suitable for essential 1D persistence." << endl;
                compute_3Dshape_global_PD_with_UF(val, function, prop, &global_PD_index, &global_PD_value, &globfunc, &neighb);
            }

            if (method.compare("filtration") == 0){
                cout << "Warning: computes homology of dim <= 2." << endl;
                compute_PtCloud_global_PD_with_filtration(val, &globfunc, &neighb);
            }
        }
        
        
    }
    
    
    return 0;
}


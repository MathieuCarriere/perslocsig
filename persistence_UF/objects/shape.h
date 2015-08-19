/* 
 * File:   shape.h
 * Author: mathieu
 *
 * Created on January 22, 2015, 2:35 PM
 */

#ifndef SHAPE_H
#define	SHAPE_H

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

#include "TriMesh.h"

using namespace std;

// Read input file.
// ****************

trimesh::TriMesh read_3Dshape(string object_name){
   
    trimesh::TriMesh *m = trimesh::TriMesh::read((char*) object_name.c_str());
    
    if (!m)
        cout << "Error while reading mesh" << endl;

    m->need_faces();
    m->need_normals();
    m->need_neighbors();
    
    return *m;

}

// Read input local function.
// **************************

void read_3Dshape_local(int val, vector<vector<double> >* loc, trimesh::TriMesh* m, string function){
    
  int num_points = (*m).vertices.size();
    
    // Geodesic distance functions.
    // ****************************

    if (function.compare("geo") == 0){
        
        char dmat[100];
        sprintf(dmat, "dist_matrix_%d", val);
        ifstream ifile(dmat);
    
        if (!ifile){
            cout << "Geodesics not computed." << endl;
            return;
        }
        
        else{
            cout << "Reading geodesics... ";
            vector<double> val;
            double values[num_points];
        
            for (int nb = 0; nb < num_points; nb++){
                int pos = nb*num_points*8;
                ifile.clear();
                ifile.seekg(pos);
                ifile.read((char*) values, 8*num_points);
            
                val.clear();
                for (int k = 0; k < num_points; k++)
                    val.push_back(values[k]);
            
                (*loc).push_back(val);
            }
            cout << "Done." << endl; 
        }
    }

    // HKS distance functions.
    // ***********************
    
    if (function.compare("hks") == 0){
        
        char dmat[100];
        sprintf(dmat, "hks_matrix_%d", val);
        ifstream ifile(dmat);
        
        if (!ifile){
            cout << "HKS not computed." << endl;
            return;
        }
        
        else{
            cout << "Reading hks... ";
            vector<double> val;
            double values[num_points];
        
            for (int nb = 0; nb < num_points; nb++){
                int pos = nb*num_points*8;
                ifile.clear();
                ifile.seekg(pos);
                ifile.read((char*) values, 8*num_points);
            
                val.clear();
                for (int k = 0; k < num_points; k++)
                    val.push_back(values[k]);
            
                (*loc).push_back(val);
            }
            cout << "Done." << endl; 
        }
    }
    
    
}

// Read input global function.
// **************************

void read_3Dshape_global(int val, vector<double>* scal, trimesh::TriMesh* m, string function){

  int num_points = (*m).vertices.size();
    
    // Curvature.
    // **********

    if (function.compare("cur") == 0){
        
        char curv[100];
        sprintf(curv, "curv_points_%d", val);
        ifstream ifile(curv);
        
        ifile.clear();
        int numfeat;
        ifile.read((char*) &numfeat, 4);
        
        vector<vector<double> > scal2(numfeat);
        for (int i = 0; i < numfeat; i++)
            for (int j = 0; j < num_points; j++)
                scal2[i].push_back(0);
    
        if (!ifile){
            cout << "Curvature not computed." << endl;
            return;
        }
        
        else{
            
            cout << "Reading curvature... ";
            double values[num_points];
            
            for (int k = 0; k < num_points; k++){
           
                ifile.read((char*) values, 8*numfeat);
            
                for (int l = 0; l < numfeat; l++)
                    scal2[l][k] = values[l];
            }
            
            (*scal) = scal2[6];
            
            cout << "Done." << endl; 
        }
    }

    // Eccentricity.
    // *************
    
    if (function.compare("ecc") == 0){
        
        char dist[100];
        sprintf(dist, "dist_matrix_%d", val);
        ifstream ifile(dist);
    
        if (!ifile){
            cout << "Geodesics not computed." << endl;
            return;
        }
        
        else{
            
            cout << "Reading geodesics... ";
            double values[num_points];
        
            for (int nb = 0; nb < num_points; nb++){
                
                int pos = nb*num_points*8;
                ifile.clear();
                ifile.seekg(pos);
                ifile.read((char*) values, 8*num_points);
                
                double ecc = values[0];
                for (int nbb = 1; nbb < num_points; nbb++)
                    if (ecc < values[nbb])
                        ecc = values[nbb];
                
                (*scal).push_back(ecc);
                
            }
            cout << "Done." << endl;
        }
    }

    // Height & Minus Height.
    // **********************

    if (function.compare("height") == 0){

        char coord[100]; sprintf(coord, "%dZ", val); ifstream ifile(coord);
    
        if (!ifile){
            cout << "Heights not computed." << endl;
            return;
        }

        else{
            
            cout << "Reading heights... ";
        
            for (int nb = 0; nb < num_points; nb++){
                
                double z; ifile >> z;
                (*scal).push_back(z);
      
            }
            cout << "Done." << endl;
        }
    }

    if (function.compare("minus_height") == 0){

        char coord[100]; sprintf(coord, "%dZ", val); ifstream ifile(coord);
    
        if (!ifile){
            cout << "Heights not computed." << endl;
            return;
        }

        else{
            
            cout << "Reading heights... ";
        
            for (int nb = 0; nb < num_points; nb++){
                
                double z; ifile >> z;
                (*scal).push_back(-z);
      
            }
            cout << "Done." << endl;
        }
    }

}

#endif	/* SHAPE_H */


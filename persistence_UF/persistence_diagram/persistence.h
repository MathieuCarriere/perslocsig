/* 
 * File:   persistence.h
 * Author: mathieu
 *
 * Created on January 22, 2015, 4:05 PM
 */

#ifndef PERSISTENCE_H
#define	PERSISTENCE_H

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
#include <string>

#include "union_find.h"

using namespace std;

struct vertex {
    unsigned int index;
    double value;
};

struct by_value {
    bool operator()(vertex const &a, vertex const &b){ 
        return a.value > b.value;
    }  
};

struct by_persistence_val {
    bool operator()(pair<double,double> const &a, pair<double,double> const &b){ 
        return (a.first-a.second) > (b.first-b.second);
    }  
};

vector<double> f;

struct by_persistence_ind {
    bool operator()(pair<int,int> const &a, pair<int,int> const &b){ 
        return (f[a.first]-f[a.second]) > (f[b.first]-f[b.second]);
    }  
};

void compute_3Dshape_local_PD_with_UF(int val, string function, string prop, vector<vector<pair<int,int> > >* PD_index, \
                                      vector<vector<pair<double,double> > >* PD_value, \
                                      vector<vector<double> >* func, \
                                      vector<vector<int> >* neighb){

    char pdind[36];
    sprintf(pdind, "PersistenceDiagram_%s_%s_Index_%d",(char*) function.c_str(),(char*) prop.c_str(), val);
    char pdval[36];
    sprintf(pdval, "PersistenceDiagram_%s_%s_Value_%d",(char*) function.c_str(),(char*) prop.c_str(), val);

    cout << "Computing persistence..." << endl;
    
    int num_points = (*neighb).size();
    
    vector<int> neighbors; vector<int> neighborsp;
    vector<pair<int,int> > pd_ind;
    vector<pair<double,double> > pd_val;
    vector<vertex> verts(num_points);
    
    vector<pair<int,int> > entries;
    vector<int> v_vect;
    vector<int> parents(num_points);
    
    for (int i = 0; i < num_points; i++){
    
        if (i % 1000 == 0)
            cout << "    source = " << i << "/" << num_points << endl;
        
        entries.clear();
        v_vect.clear();
        pd_ind.clear();
        pd_val.clear();
        neighbors.clear();
        
        f = (*func)[i];
        
        for (int j = 0; j < f.size(); j++){
            verts[j].index = j;
            verts[j].value = f[j];
        }
        
        sort(verts.begin(), verts.end(), by_value());
        
        for (int j = 0; j < f.size(); j++){
        
            int id = verts[j].index;
            int max_neighb = id;
            parents[id] = id;
            double value = verts[j].value;
            neighbors.clear(); neighborsp.clear();
            
            // Find neighbors of current node.
            // *******************************
            
            // int num_neighb = shape.neighbors[id].size();
            int num_neighb = (*neighb)[id].size();
            
            for (int k = 0; k < num_neighb; k++){
                
                // int n = shape.neighbors[id][k];
	        int n = (*neighb)[id][k];
                if (f[n] > f[id])
                    neighbors.push_back(n);

            }
            
            
            
            // If current node is not a peak, attach it to peak with highest value.
            // ********************************************************************
            
            if (neighbors.size() >= 1){
                
                int num_neighbors = neighbors.size();
                
                for (int k = 0; k < num_neighbors; k++){
                
                    int n = UF_find(neighbors[k], &parents);
                    int pn;
                    
                    for (int k = 0; k < entries.size(); k++)
                        if (entries[k].first == n){
                            pn = v_vect[k];
                            neighborsp.push_back(pn);
                            break;
                        }
                    
                    if (  (pn != max_neighb && f[pn] > f[max_neighb])  ){
                        max_neighb = pn;
                        parents[id] = n;
                    }
                    
                }
                
                int v = max_neighb;
                int vp = UF_find(v, &parents);
                
                for (int k = 0; k < entries.size(); k++)
                    if (entries[k].first == vp){
                        entries[k].second += 1;
                        break;
                    }
                
                
                // and merge all other cc.
                // ***********************
                
                for (int k = 0; k < num_neighbors; k++){
                
                    int n = neighbors[k];
                    int indn, indv;
                    int np = UF_find(n, &parents);
                    int vp = UF_find(v, &parents);
                    
                    for (int k = 0; k < entries.size(); k++)
                        if (entries[k].first == vp){
                            indv = k;
                            break;
                        }
                    
                    for (int k = 0; k < entries.size(); k++)
                        if (entries[k].first == np){
                            indn = k;
                            break;
                        }
                    
                    if (np != vp){
                        
                        pair<int,int> pt1(v_vect[indn],id);
                        pair<double,double> pt2(f[v_vect[indn]],value);
                        pd_ind.push_back(pt1);
                        pd_val.push_back(pt2);
                    
                        UF_union(n, v, v, &entries, &parents, &v_vect);
     
                    }
                    
                }
            
            }
            
            // If current node is a peak, create new cc.
            // *****************************************
            
            else{
                
                pair<int,int> new_cc(id,1);
                entries.push_back(new_cc);
                v_vect.push_back(id);
                parents[id] = id;
            
            }
            
            if (j == (f.size()-1)){
                pair<int,int> pt1(v_vect[0],id);
                pair<double,double> pt2(f[v_vect[0]],value);
                pd_ind.push_back(pt1);
                pd_val.push_back(pt2);
            }
            
        }

        sort(pd_ind.begin(), pd_ind.end(), by_persistence_ind());
        sort(pd_val.begin(), pd_val.end(), by_persistence_val());
        
        (*PD_index).push_back(pd_ind);
        (*PD_value).push_back(pd_val);
    
    }
    
    cout << "Persistence computed." << endl;
    
    /*
    ofstream PDIND0(pdind0); ofstream PDVAL0(pdval0);
    ofstream PDIND1(pdind1); ofstream PDVAL1(pdval1);
    ofstream PDIND2(pdind2); ofstream PDVAL2(pdval2);
    
    for (int i = 0; i < num_points; i++){
        vector<pair<int,int> > pdi = (*PD_index)[i];
        vector<pair<double,double> > pdv = (*PD_value)[i];
        int num_h = pdi.size();
        PDIND0 << pdi[0].second << " -" << endl; PDVAL0 << pdv[0].second << " inf" << endl;
        PDIND2 << pdi[0].first << " -" << endl; PDVAL2 << pdv[0].first << " inf" << endl;
        for (int j = 1; j < num_h; j++){
            PDIND1 << pdi[j].second << " " << pdi[j].first << " "; 
            PDVAL1 << pdv[j].second << " " << pdv[j].first << " ";
        }
        PDIND1 << endl; PDVAL1 << endl;
    }
    
    PDIND0.close(); PDVAL0.close();
    PDIND1.close(); PDVAL1.close();
    PDIND2.close(); PDVAL2.close();
    */
    
    ofstream PDIND(pdind); ofstream PDVAL(pdval);
    
    for (int i = 0; i < num_points; i++){
        vector<pair<int,int> > pdi = (*PD_index)[i];
        vector<pair<double,double> > pdv = (*PD_value)[i];
        int num_h = pdi.size();
        for (int j = 0; j < num_h; j++){
            PDIND << pdi[j].second << " " << pdi[j].first << " "; 
            PDVAL << pdv[j].second << " " << pdv[j].first << " ";
        }
        PDIND << endl; PDVAL << endl;
    }
    
    PDIND.close(); PDVAL.close();
    
}

void compute_3Dshape_global_PD_with_UF(int val, string function, string prop, vector<pair<int,int> >* PD_index, \
                                       vector<pair<double,double> >* PD_value, \
                                       vector<double>* globfunc, \
                                       vector<vector<int> >* neighb){
    
    char pdind[100];
    sprintf(pdind, "PersistenceDiagram_%s_%s_Index_%d", (char*) function.c_str(), (char*) prop.c_str(), val);
    char pdval[100];
    sprintf(pdval, "PersistenceDiagram_%s_%s_Value_%d", (char*) function.c_str(), (char*) prop.c_str(), val);
    
    cout << "Computing persistence..." << endl;
    
    int num_points = (*neighb).size();
    
    vector<int> neighbors, neighborsp;
    vector<pair<int,int> > pd_ind;
    vector<pair<double,double> > pd_val;
    vector<vertex> verts(num_points);
    
    vector<pair<int,int> > entries;
    vector<int> v_vect;
    vector<int> parents(num_points);
        
    entries.clear();
    v_vect.clear();
    pd_ind.clear();
    pd_val.clear();
    neighbors.clear();
        
    f = (*globfunc);
    
    for (int j = 0; j < f.size(); j++){
        verts[j].index = j; verts[j].value = f[j];}
        
    sort(verts.begin(), verts.end(), by_value());
        
    for (int j = 0; j < f.size(); j++){
        
            int id = verts[j].index;
            int max_neighb = id;
            parents[id] = id;
            double value = verts[j].value;
            neighbors.clear(); neighborsp.clear();
            
            // Find neighbors of current node.
            // *******************************
            
            int num_neighb = (*neighb)[id].size();
            
            for (int k = 0; k < num_neighb; k++){
                
	        int n = (*neighb)[id][k];
                if (f[n] > f[id])
                    neighbors.push_back(n);

            }
            
            // If current node is not a peak, attach it to peak with highest value.
            // ********************************************************************
            
            if (neighbors.size() >= 1){
                
                int num_neighbors = neighbors.size();
                
                for (int k = 0; k < num_neighbors; k++){
                
                    int n = UF_find(neighbors[k], &parents);
                    int pn;
                    
                    for (int k = 0; k < entries.size(); k++)
                        if (entries[k].first == n){
                            pn = v_vect[k];
                            neighborsp.push_back(pn);
                            break;
                        }
                    
                    if (  (pn != max_neighb && f[pn] > f[max_neighb])  ){
                        max_neighb = pn;
                        parents[id] = n;
                    }
                    
                }
                
                int v = max_neighb;
                int vp = UF_find(v, &parents);
                
                for (int k = 0; k < entries.size(); k++)
                    if (entries[k].first == vp){
                        entries[k].second += 1;
                        break;
                    }
                
                
                // and merge all other cc.
                // ***********************
                
                for (int k = 0; k < num_neighbors; k++){
                
                    int n = neighbors[k];
                    int indn, indv;
                    int np = UF_find(n, &parents);
                    int vp = UF_find(v, &parents);
                    
                    for (int k = 0; k < entries.size(); k++)
                        if (entries[k].first == vp){
                            indv = k;
                            break;
                        }
                    
                    for (int k = 0; k < entries.size(); k++)
                        if (entries[k].first == np){
                            indn = k;
                            break;
                        }
                    
                    if (np != vp){
                        
                        pair<int,int> pt1(v_vect[indn],id);
                        pair<double,double> pt2(f[v_vect[indn]],value);
                        pd_ind.push_back(pt1);
                        pd_val.push_back(pt2);
                    
                        UF_union(n, v, v, &entries, &parents, &v_vect);
     
                    }
                    
                }
            
            }
            
            // If current node is a peak, create new cc.
            // *****************************************
            
            else{
                
                pair<int,int> new_cc(id,1);
                entries.push_back(new_cc);
                v_vect.push_back(id);
                parents[id] = id;
            
            }
            
            if (j == (f.size()-1)){
                pair<int,int> pt1(v_vect[0],id);
                pair<double,double> pt2(f[v_vect[0]],value);
                pd_ind.push_back(pt1);
                pd_val.push_back(pt2);
            }
            
    }

    sort(pd_ind.begin(), pd_ind.end(), by_persistence_ind());
    sort(pd_val.begin(), pd_val.end(), by_persistence_val());
    
    int num_h = pd_ind.size();
    
    (*PD_index) = pd_ind; (*PD_value) = pd_val;
    
    cout << "Persistence computed." << endl;
    
    /*
    ofstream PDIND0(pdind0); ofstream PDVAL0(pdval0);
    ofstream PDIND1(pdind1); ofstream PDVAL1(pdval1);
    ofstream PDIND2(pdind2); ofstream PDVAL2(pdval2);
    
    vector<pair<int,int> > pdi = (*PD_index);
    vector<pair<double,double> > pdv = (*PD_value);

    PDIND0 << pdi[0].second << " -" << endl; PDVAL0 << pdv[0].second << " inf" << endl;
    PDIND2 << pdi[0].first << " -" << endl; PDVAL2 << pdv[0].first << " inf" << endl;
    for (int j = 1; j < num_h; j++){
        PDIND1 << pdi[j].second << " " << pdi[j].first << " "; 
        PDVAL1 << pdv[j].second << " " << pdv[j].first << " ";
    }

    PDIND0.close(); PDVAL0.close();
    PDIND1.close(); PDVAL1.close();
    PDIND2.close(); PDVAL2.close();
    */
    
    vector<pair<int,int> > pdi = (*PD_index);
    vector<pair<double,double> > pdv = (*PD_value);
    
    ofstream PDIND(pdind); ofstream PDVAL(pdval);
    
    for (int j = 0; j < num_h; j++){
        PDIND << pdi[j].second << " " << pdi[j].first << " "; 
        PDVAL << pdv[j].second << " " << pdv[j].first << " ";
    }
    PDIND << endl; PDVAL << endl;
    
    
    PDIND.close(); PDVAL.close();
    
    
}

void compute_3Dshape_local_PD_with_filtration(int val, vector<vector<double> >* func, \
                                              trimesh::TriMesh* shape){
    
    int num_points = (*shape).vertices.size();
    int num_faces = (*shape).faces.size();
    
    vector<vertex> verts(num_points);
    
    cout << "Computing persistence..." << endl;
    
    for (int i = 0; i < num_points; i++){
        
        char strFil[10];
        sprintf(strFil, "Filtration");
        ofstream FIL(strFil);
        
        if (i % 100 == 0)
            cout << "    source = " << i << "/" << num_points << endl;
      
        f = (*func)[i];
        
        for (int j = 0; j < num_points; j++){
            verts[j].index = j;
            verts[j].value = f[j];
        }
        
        sort(verts.begin(), verts.end(), by_value());
        
        // Add simplices of dimension 0.
        // *****************************
        
        for (int j = 0; j < num_points; j++)
            FIL << 0 << " " << j << " " << f[j] << endl;
        
        // Add simplices of dimension 1&2.
        // *******************************
        
        for (int j = 0; j < num_faces; j++){
            
	    trimesh::TriMesh::Face face = (*shape).faces[j];
            int v0 = face.v[0]; int v1 = face.v[1]; int v2 = face.v[2];
            double x; 
            
            x = max(f[v0],max(f[v1],f[v2]));
            FIL << 2 << " " << v0 << " " << v1 << " " << v2 << " " << x << endl; 

            x = max(f[v0],f[v1]);
            FIL << 1 << " " << v0 << " " << v1 << " " << x << endl;
            
            x = max(f[v1],f[v2]);
            FIL << 1 << " " << v1 << " " << v2 << " " << x << endl;
            
            x = max(f[v0],f[v2]);
            FIL << 1 << " " << v0 << " " << v2 << " " << x << endl;
            
        }
        
        char command[100];
        sprintf(command, "./zomorodian_persistence Filtration 3");
        char command0[100];
        sprintf(command0, "cat PersistenceDiagram0_Value_%d PersistenceDiagram0 >> PersistenceDiagram0_Value_%d 2>/dev/null",val,val);
        char command1[100];
        sprintf(command1, "cat PersistenceDiagram1_Value_%d PersistenceDiagram1 >> PersistenceDiagram1_Value_%d 2>/dev/null",val,val);
        char command2[100];
        sprintf(command2, "cat PersistenceDiagram2_Value_%d PersistenceDiagram2 >> PersistenceDiagram2_Value_%d 2>/dev/null",val,val);
                         
        system(command);
        system(command0);
        system(command1);
        system(command2);
        
    }
    
    cout << "Persistence computed." << endl;
    
    char command[100];
    sprintf(command, "rm PersistenceDiagram0 PersistenceDiagram1 PersistenceDiagram2 Filtration");
    system(command);
    
}

void compute_3Dshape_global_PD_with_filtration(int val, vector<double>* func, \
                                               trimesh::TriMesh* shape){
    
    int num_points = (*shape).vertices.size();
    int num_faces = (*shape).faces.size();
    
    vector<vertex> verts(num_points);
        
    char strFil[10];
    sprintf(strFil, "Filtration");
    ofstream FIL(strFil);
      
    f = (*func);
        
    for (int j = 0; j < num_points; j++){
        verts[j].index = j;
        verts[j].value = f[j];
    }
        
    sort(verts.begin(), verts.end(), by_value());
        
    // Add simplices of dimension 0.
    // *****************************
        
    for (int j = 0; j < num_points; j++)
        FIL << 0 << " " << j << " " << f[j] << endl;
        
    // Add simplices of dimension 1&2.
    // *******************************
        
    for (int j = 0; j < num_faces; j++){
            
        trimesh::TriMesh::Face face = (*shape).faces[j];
        int v0 = face.v[0]; int v1 = face.v[1]; int v2 = face.v[2];
        double x; 
            
        x = max(f[v0],max(f[v1],f[v2]));
        FIL << 2 << " " << v0 << " " << v1 << " " << v2 << " " << x << endl; 

        x = max(f[v0],f[v1]);
        FIL << 1 << " " << v0 << " " << v1 << " " << x << endl;
            
        x = max(f[v1],f[v2]);
        FIL << 1 << " " << v1 << " " << v2 << " " << x << endl;
            
        x = max(f[v0],f[v2]);
        FIL << 1 << " " << v0 << " " << v2 << " " << x << endl;
            
    }
    
    cout << "Computing persistence..." << endl;
        
    char command[100];
    sprintf(command, "./zomorodian_persistence Filtration 3");
    char command0[100];
    sprintf(command0, "mv PersistenceDiagram0 PersistenceDiagram0_Value_%d",val);
    char command1[100];
    sprintf(command1, "mv PersistenceDiagram1 PersistenceDiagram1_Value_%d",val);
    char command2[100];
    sprintf(command2, "mv PersistenceDiagram2 PersistenceDiagram2_Value_%d",val);
    
    cout << "Persistence computed." << endl;
                         
    system(command);
    system(command0);
    system(command1);
    system(command2);
    
}

void compute_PtCloud_local_PD_with_filtration(int val, vector<vector<double> >* func, \
                                              vector<vector<int> >* neighb){
    
    int num_points = (*func).size();
    
    vector<vertex> verts(num_points);
    
    cout << "Computing persistence..." << endl;
    
    for (int i = 0; i < num_points; i++){
        
        char strFil[10];
        sprintf(strFil, "Filtration");
        ofstream FIL(strFil);
        
        if (i % 100 == 0)
            cout << "    source = " << i << "/" << num_points << endl;
      
        f = (*func)[i];
        
        for (int j = 0; j < num_points; j++){
            verts[j].index = j;
            verts[j].value = f[j];
        }
        
        sort(verts.begin(), verts.end(), by_value());
        
        // Add simplices of dimension 0.
        // *****************************
        
        for (int j = 0; j < num_points; j++)
            FIL << 0 << " " << j << " " << f[j] << endl;
        
        // Add simplices of dimension 1&2.
        // *******************************
        
        for (int j = 0; j < num_points; j++){            
        int v0 = j; vector<int> nei1 = (*neighb)[j];
        for (int k = 0; k < nei1.size(); k++){
            int v1 = nei1[k]; vector<int> nei2 = (*neighb)[v1];
            for (int l = 0; l < nei2.size(); l++){
                int v2 = nei2[l]; vector<int> nei3 = (*neighb)[v2];
                if (find(nei3.begin(),nei3.end(),v0)!=nei3.end()){
                    
                    double x; 
            
		    x = max(f[v0],max(f[v1],f[v2]));
        	    FIL << 2 << " " << v0 << " " << v1 << " " << v2 << " " << x << endl; 

        	    x = max(f[v0],f[v1]);
    		    FIL << 1 << " " << v0 << " " << v1 << " " << x << endl;
            
        	    x = max(f[v1],f[v2]);
            	    FIL << 1 << " " << v1 << " " << v2 << " " << x << endl;
            
        	    x = max(f[v0],f[v2]);
            	    FIL << 1 << " " << v0 << " " << v2 << " " << x << endl;
		
		}
            }
        }}    
    
        
        char command[100];
        sprintf(command, "./zomorodian_persistence Filtration 3");
        char command0[100];
        sprintf(command0, "cat PersistenceDiagram0_Value_%d PersistenceDiagram0 >> PersistenceDiagram0_Value_%d 2>/dev/null",val,val);
        char command1[100];
        sprintf(command1, "cat PersistenceDiagram1_Value_%d PersistenceDiagram1 >> PersistenceDiagram1_Value_%d 2>/dev/null",val,val);
        char command2[100];
        sprintf(command2, "cat PersistenceDiagram2_Value_%d PersistenceDiagram2 >> PersistenceDiagram2_Value_%d 2>/dev/null",val,val);
                         
        system(command);
        system(command0);
        system(command1);
        system(command2);
        
    }
    
    cout << "Persistence computed." << endl;
    
    char command[100];
    sprintf(command, "rm PersistenceDiagram0 PersistenceDiagram1 PersistenceDiagram2 Filtration");
    system(command);
    
}

void compute_PtCloud_global_PD_with_filtration(int val, vector<double>* func, \
                                               vector<vector<int> >* neighb){
    
    int num_points = (*func).size();
    
    vector<vertex> verts(num_points);
        
    char strFil[10];
    sprintf(strFil, "Filtration");
    ofstream FIL(strFil);
      
    f = (*func);
        
    for (int j = 0; j < num_points; j++){
        verts[j].index = j;
        verts[j].value = f[j];
    }
        
    sort(verts.begin(), verts.end(), by_value());
        
    // Add simplices of dimension 0.
    // *****************************
        
    for (int j = 0; j < num_points; j++)
        FIL << 0 << " " << j << " " << f[j] << endl;
        
    // Add simplices of dimension 1&2.
    // *******************************
        
    for (int j = 0; j < num_points; j++){            
        int v0 = j; vector<int> nei1 = (*neighb)[j];
        for (int k = 0; k < nei1.size(); k++){
            int v1 = nei1[k]; vector<int> nei2 = (*neighb)[v1];
            for (int l = 0; l < nei2.size(); l++){
                int v2 = nei2[l]; vector<int> nei3 = (*neighb)[v2];
                if (find(nei3.begin(),nei3.end(),v0)!=nei3.end()){
                    
                    double x; 
            
		    x = max(f[v0],max(f[v1],f[v2]));
        	    FIL << 2 << " " << v0 << " " << v1 << " " << v2 << " " << x << endl; 

        	    x = max(f[v0],f[v1]);
    		    FIL << 1 << " " << v0 << " " << v1 << " " << x << endl;
            
        	    x = max(f[v1],f[v2]);
            	    FIL << 1 << " " << v1 << " " << v2 << " " << x << endl;
            
        	    x = max(f[v0],f[v2]);
            	    FIL << 1 << " " << v0 << " " << v2 << " " << x << endl;
		
		}
            }
        }    
    }
    
    cout << "Computing persistence..." << endl;
        
    char command[100];
    sprintf(command, "./zomorodian_persistence Filtration 3");
    char command0[100];
    sprintf(command0, "mv PersistenceDiagram0 PersistenceDiagram0_Value_%d",val);
    char command1[100];
    sprintf(command1, "mv PersistenceDiagram1 PersistenceDiagram1_Value_%d",val);
    char command2[100];
    sprintf(command2, "mv PersistenceDiagram2 PersistenceDiagram2_Value_%d",val);
    
    cout << "Persistence computed." << endl;
                         
    system(command);
    system(command0);
    system(command1);
    system(command2);
    
}
 
#endif	/* PERSISTENCE_H */


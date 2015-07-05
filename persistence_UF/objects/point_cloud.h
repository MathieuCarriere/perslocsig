#ifndef POINT_CLOUD_H
#define	POINT_CLOUD_H

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

#include "ANN/ANN.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <list>
#include <set>

// Dijkstra for shortest paths.
// ****************************

typedef int vertex_t;
typedef double weight_t;

struct edge {
    vertex_t target;
    weight_t weight;
    edge(vertex_t arg_target, weight_t arg_weight)
        : target(arg_target), weight(arg_weight) { }
};

typedef std::map<vertex_t, std::list<edge> > adjacency_map_t;

template <typename T1, typename T2>
struct pair_first_less
{
    bool operator()(std::pair<T1,T2> p1, std::pair<T1,T2> p2)
    {
        return p1.first < p2.first;
    }
};


void DijkstraComputePaths(vertex_t source,
                          adjacency_map_t& adjacency_map,
                          std::map<vertex_t, weight_t>& min_distance,
                          std::map<vertex_t, vertex_t>& previous)
{
    for (adjacency_map_t::iterator vertex_iter = adjacency_map.begin();
         vertex_iter != adjacency_map.end();
         vertex_iter++)
    {
        vertex_t v = vertex_iter->first;
        //min_distance[v] = std::numeric_limits< double >::infinity();
        min_distance[v] = 1e100;
    }
    min_distance[source] = 0;
    std::set< std::pair<weight_t, vertex_t>,
              pair_first_less<weight_t, vertex_t> > vertex_queue;
    for (adjacency_map_t::iterator vertex_iter = adjacency_map.begin();
         vertex_iter != adjacency_map.end();
         vertex_iter++)
    {
        vertex_t v = vertex_iter->first;
        vertex_queue.insert(std::pair<weight_t, vertex_t>(min_distance[v], v));
    }

    while (!vertex_queue.empty()) {
        vertex_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());

        // Visit each edge exiting u
        for (std::list<edge>::iterator edge_iter = adjacency_map[u].begin();
             edge_iter != adjacency_map[u].end();
             edge_iter++)
        {
            vertex_t v = edge_iter->target;
            weight_t weight = edge_iter->weight;
            weight_t distance_through_u = min_distance[u] + weight;
	    if (distance_through_u < min_distance[v]) {
	        vertex_queue.erase(std::pair<weight_t, vertex_t>(min_distance[v], v));

	        min_distance[v] = distance_through_u;
	        previous[v] = u;
	        vertex_queue.insert(std::pair<weight_t, vertex_t>(min_distance[v], v));
	    }
        }
    }
}

std::list<vertex_t> DijkstraGetShortestPathTo(
    vertex_t target, std::map<vertex_t, vertex_t>& previous)
{
    std::list<vertex_t> path;
    std::map<vertex_t, vertex_t>::iterator prev;
    vertex_t vertex = target;
    path.push_front(vertex);
    while((prev = previous.find(vertex)) != previous.end())
    {
        vertex = prev->second;
        path.push_front(vertex);
    }
    return path;
}

using namespace std;

// Read input file & compute Neighborhood Graph.
// *********************************************

void build_neighborhood_graph(string name, double delta, vector<vector<int> >* edges, vector<vector<double> >* dist){

    ifstream PC((char*) name.c_str());
    vector<vector<double> > point_cloud;
    string line;
    
    while(getline(PC,line)){
    
        stringstream stream(line);
        vector<double> v; double x;
        while(stream >> x)
            	v.push_back(x);
        point_cloud.push_back(v);
        vector<int> vec; vec.clear(); (*edges).push_back(vec); 
        vector<double> vec2; vec2.clear(); (*dist).push_back(vec2);
    
    }
    
    int num_points = point_cloud.size(); int dim = point_cloud[0].size();
    ANNpointArray dataPts; dataPts = annAllocPts(num_points, dim);
    for(int i = 0; i < num_points; i++)
        for(int j = 0; j < dim; j++)
            dataPts[i][j] = point_cloud[i][j];
    
    ANNkd_tree* kd_tree = new ANNkd_tree(dataPts,num_points,dim);
    
    cout << "Computing neighborhood graph (1-skeleton Rips)..." << endl;
    
    for(int i = 0; i < num_points; i++){
        
        ANNpoint p = annAllocPt(dim);
        p[0] = point_cloud[i][0]; p[1] = point_cloud[i][1]; p[2] = point_cloud[i][2];
        ANNdist d = delta*delta;
        int k = 1000;
        ANNidxArray idxs = new ANNidx[k];
        ANNdistArray distances = new ANNdist[k];
        
        kd_tree->annkFRSearch(p,d,k,idxs,distances);
        int iter = 1;
        while(idxs[iter] != -1){
            (*edges)[i].push_back(idxs[iter]); 
            (*dist)[i].push_back(sqrt(distances[iter]));
            iter++;
        }
        
        if(i%1000 == 0)
            cout << "    i = " << i << "/" << num_points << " -- " << iter-1 << " neighbors" << endl;
        
        delete[] p; delete[] idxs; delete[] distances; 
    }

}

// Read input local function.
// **************************

void read_pointcloud_local(int val, vector<vector<double> >* loc, vector<vector<int> >* neighb, vector<vector<double> >* dist, string function){
    
    int num_points = (*neighb).size();

    // Geodesic distance functions.
    // ****************************
    
    if (function.compare("geo") == 0){
        
        char dmat[100];
        sprintf(dmat, "dist_matrix_%d", val);
        ifstream ifile(dmat);
    
        if (!ifile){
            cout << "Computing geodesics..." << endl;
            
            // VERY SLOW !!! Use Boost instead
            
            ofstream ofile(dmat);         
            adjacency_map_t adjacency_map;

            for(int i = 0; i < num_points; i++)
                for (int j = 0; j < neighb[i].size(); j++)
                    adjacency_map[i].push_back(edge((*neighb)[i][j],  (*dist)[i][j]));
            
            std::map<vertex_t, weight_t> min_distance;
            std::map<vertex_t, vertex_t> previous;
            for(int i = 0; i < num_points; i++){
                if (i%1000 == 0)
                    cout << "    i = " << i << "/" << num_points << endl;
                DijkstraComputePaths(i, adjacency_map, min_distance, previous);
                double values[num_points];
                vector<double> vval; vval.clear();
                for (int j = 0; j < num_points; j++){
                    values[j] = min_distance[j];
                    vval.push_back(min_distance[j]);
                }
                (*loc).push_back(vval);
                ofile.write((char*) values, 8*num_points);
                
            }
            
            ofile.close();
            
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
        sprintf(dmat, "hks_0.05_matrix_%d", val);
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

// Read input global functions.
// ****************************

void read_pointcloud_global(int val, vector<double>* scal, vector<vector<int> >* neighb, vector<vector<double> >* dist, string function){

    int num_points = (*neighb).size();

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

    // Particular geodesic distance function.
    // **************************************
    
    if (function.compare("geo") == 0){
        
        char dist[100]; sprintf(dist, "dist_matrix_%d", val); ifstream ifile(dist);
    
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

#endif	/* POINT_CLOUD_H */


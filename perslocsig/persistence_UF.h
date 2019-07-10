/*    Author:       Mathieu Carriere
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <string>

#include "union_find.h"
#include "graph_geodesic.h"
#include "graph_neighbors.h"

using namespace std;

struct vertex {
    unsigned int index;
    double value;
};

using Persistence_idx_val = pair<vector<pair<int,int> >, vector<pair<double,double> > >;

Persistence_idx_val compute_PD_with_UF(vector<double> & f, vector<vector<int> > & neighb){

    int num_points = neighb.size();
    vector<int> neighbors, neighborsp;
    vector<pair<int,int> > pd_ind;
    vector<pair<double,double> > pd_val;
    vector<vertex> verts(num_points);
    vector<pair<int,int> > entries;
    vector<int> v_vect;
    vector<int> parents(num_points);

    for (int j = 0; j < num_points; j++){
        verts[j].index = j; verts[j].value = f[j];}

    sort(verts.begin(), verts.end(), [=](const vertex & a, const vertex & b){return a.value > b.value;});

    for (int j = 0; j < num_points; j++){

            int id = verts[j].index;
            int max_neighb = id;
            parents[id] = id;
            double value = verts[j].value;
            neighbors.clear(); neighborsp.clear();

            // Find neighbors of current node.
            // *******************************

            int num_neighb = neighb[id].size();

            for (int k = 0; k < num_neighb; k++){

                int n = neighb[id][k];
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

                    for (size_t k = 0; k < entries.size(); k++)
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

                for (size_t k = 0; k < entries.size(); k++)
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

                    for (size_t k = 0; k < entries.size(); k++)
                        if (entries[k].first == vp){
                            indv = k;
                            break;
                        }

                    for (size_t k = 0; k < entries.size(); k++)
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

            if (j == (num_points-1)){
                pair<int,int> pt1(v_vect[0],id);
                pair<double,double> pt2(f[v_vect[0]],value);
                pd_ind.push_back(pt1);
                pd_val.push_back(pt2);
            }

    }

    sort(pd_ind.begin(), pd_ind.end(), [=](const pair<int,int> & a, const pair<int,int> & b){return abs(f[a.first]-f[a.second]) > abs(f[b.first]-f[b.second]);});
    sort(pd_val.begin(), pd_val.end(), [=](const pair<double,double> & a, const pair<double,double> & b){return abs(a.first-a.second) > abs(b.first-b.second);});

    Persistence_idx_val diag; diag.first = pd_ind; diag.second = pd_val;

    return diag;

}

vector<Persistence_idx_val> compute_PDs_from_off(string name, vector<int> idx){

  vector<vector<int> > neighbors = graph_neighbors(name);
  vector<vector<double> > dist   = graph_geodesic(name);
  int numpts = idx.size();
  vector<Persistence_idx_val> PDs(numpts);
  for(int i = 0; i < numpts; i++)  PDs[i] = compute_PD_with_UF(dist[idx[i]], neighbors);

  return PDs;

}

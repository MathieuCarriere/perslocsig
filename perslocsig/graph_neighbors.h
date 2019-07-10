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
#include <vector>
#include <cmath>
#include <string>
#include <iterator>
#include <algorithm>

using namespace std;

vector<vector<int> > graph_neighbors(string name){

    int numpts, numfaces, numedges;

    ifstream input(name); string line; int i;
    getline(input, line); // Read "OFF"
    getline(input, line); stringstream stream(line); stream >> numpts; stream >> numfaces; stream >> numedges; vector<vector<int> > neighbors(numpts);

    i = 0;
    while (i < numpts) { getline(input, line); i++; }

    i = 0;
    while (i < numfaces) {
      getline(input, line); stringstream iss(line); vector<int> simplex;
      simplex.assign(istream_iterator<int>(iss), istream_iterator<int>());
      int dim = simplex[0];
      for (int j = 1; j <= dim; j++){
        for (int k = j + 1; k <= dim; k++){
          neighbors[simplex[j]].push_back(simplex[k]); neighbors[simplex[k]].push_back(simplex[j]);}}
      i++;
    }

    for (int j = 0; j < numpts; j++){
      vector<int> neighbors_j = neighbors[j];
      sort(neighbors[j].begin(), neighbors[j].end()); vector<int>::iterator it = unique(neighbors[j].begin(), neighbors[j].end()); neighbors[j].resize(distance(neighbors[j].begin(), it));
    }

    return neighbors;

}

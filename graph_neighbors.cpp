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

int main (int argc, char** argv) {

    string name = argv[1]; int numpts, numfaces, numedges;

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
      sort(neighbors_j.begin(), neighbors_j.end()); vector<int>::iterator it = unique(neighbors_j.begin(), neighbors_j.end()); neighbors_j.resize(distance(neighbors_j.begin(), it));
      for (int k = 0; k < neighbors_j.size(); k++)  cout << neighbors_j[k] << " ";
      cout << endl;
    }

    return 0;
}

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

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>

using Graph = boost::subgraph<
    boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                          boost::property<boost::edge_index_t, int, boost::property<boost::edge_weight_t, double> > > >;
using Vertex_t = boost::graph_traits<Graph>::vertex_descriptor;
using Index_map = boost::property_map<Graph, boost::vertex_index_t>::type;
using Weight_map = boost::property_map<Graph, boost::edge_weight_t>::type;

using namespace std;

vector<vector<double> > graph_geodesic(string name){

    int numpts, numfaces, numedges;
    Graph graph; vector<Vertex_t> vertices; vector<vector<double> > distances; vector<vector<double> > point_cloud;

    // Read file.

    ifstream input(name); string line; int i;
    getline(input, line); // Read "OFF"
    getline(input, line); stringstream stream(line); stream >> numpts; stream >> numfaces; stream >> numedges;

    i = 0;
    while (i < numpts) {
      getline(input, line); stringstream iss(line); vector<double> point;
      point.assign(istream_iterator<double>(iss), istream_iterator<double>()); point_cloud.emplace_back(point.begin(), point.begin() + 3);
      vertices.push_back(boost::add_vertex(graph)); distances.emplace_back(numpts);
      i++;
    }

    i = 0;
    while (i < numfaces) {
      getline(input, line); stringstream iss(line); vector<int> simplex;
      simplex.assign(istream_iterator<int>(iss), istream_iterator<int>());
      int dim = simplex[0];
      for (int j = 1; j <= dim; j++)
        for (int k = j + 1; k <= dim; k++)
          boost::add_edge(vertices[simplex[j]], vertices[simplex[k]], graph);
      i++;
    }

    // Compute distances.

    Index_map index = boost::get(boost::vertex_index, graph); Weight_map weight = boost::get(boost::edge_weight, graph); boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei){
      int source = index[boost::source(*ei, graph)]; int target = index[boost::target(*ei, graph)];
      double dis = 0; for (int k = 0; k < 3; k++) dis += (point_cloud[source][k] - point_cloud[target][k])*(point_cloud[source][k] - point_cloud[target][k]); dis = sqrt(dis);
      boost::put(weight, *ei, dis);
    }

    // Dijkstra.

    for (int idx = 0; idx < numpts; idx++)
      boost::dijkstra_shortest_paths(graph, vertices[idx], boost::weight_map(weight).distance_map(boost::make_iterator_property_map(distances[idx].begin(), index)));

    return distances;
}

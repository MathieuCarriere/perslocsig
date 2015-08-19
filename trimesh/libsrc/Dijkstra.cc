#include "Dijkstra.h"
#include <cassert>

using namespace std;

namespace trimesh {
GeodesicTraversal::GeodesicTraversal(TriMesh& mesh) {
	_mesh = mesh;
	_visited = new bool[_mesh.faces.size()];
	_currentGeodesicDistance = new float[_mesh.faces.size()];
	_mesh.need_pointareas();
	_mesh.need_neighbors();
	_mesh.need_across_edge();
	_mesh.need_adjacentfaces();
}

GeodesicTraversal::~GeodesicTraversal() {
	delete[] _visited;
	delete[] _currentGeodesicDistance;
}

float* GeodesicTraversal::getGeodDistances() {
	return _currentGeodesicDistance;
}
float GeodesicTraversal::traverse(int sourceVertex) {
	assert( _pq.empty() );

	_pq.push( vertexForDijkstra( sourceVertex, 0.0f) );
	for (int r = 0; r < _mesh.vertices.size(); r++) {
		_currentGeodesicDistance[r] = FLT_MAX;
		_visited[r] = false;
	}
	_currentGeodesicDistance[sourceVertex] = 0.0f;
	float maxGeodesicDistance = 0.0f;

	while (!_pq.empty()) {
		vertexForDijkstra v = _pq.top();
		int vi = v.getVertexIndex();
		_pq.pop();
		if (_visited[vi]) {
			continue;
		}
		_visited[vi] = true;
		if (v.getGeodDistance() > maxGeodesicDistance) {
			maxGeodesicDistance = v.getGeodDistance();
		}

		for (int j = 0; j < _mesh.neighbors[vi].size(); j++) {
			int vj = _mesh.neighbors[vi][j];
			if (!_visited[vj]) {
				float dist_vivj = dist(_mesh.vertices[vi], _mesh.vertices[vj]);
				if ( _currentGeodesicDistance[vj] > _currentGeodesicDistance[vi] + dist_vivj ) {
					_currentGeodesicDistance[vj] = _currentGeodesicDistance[vi] + dist_vivj;
					_pq.push( vertexForDijkstra(vj, _currentGeodesicDistance[vj]) );
				}
			}
		}
	}

	return maxGeodesicDistance;
}


vector<int> GeodesicTraversal::traverse(int sourceVertex, float maxGeodesicDistance) {
	assert( _pq.empty() );

	vector<int> visitedVertices;

	_pq.push( vertexForDijkstra( sourceVertex, 0.0f) );
	for (int r = 0; r < _mesh.vertices.size(); r++) {
		_currentGeodesicDistance[r] = FLT_MAX;
		_visited[r] = false;
	}
	_currentGeodesicDistance[sourceVertex] = 0.0f;

	while (!_pq.empty()) {
		vertexForDijkstra v = _pq.top();
		int vi = v.getVertexIndex();
		_pq.pop();
		if (_visited[vi]) {
			continue;
		}
		_visited[vi] = true;
		visitedVertices.push_back( vi );
		if (v.getGeodDistance() > maxGeodesicDistance) {
			continue;
		}

		for (int j = 0; j < _mesh.neighbors[vi].size(); j++) {
			int vj = _mesh.neighbors[vi][j];
			if (!_visited[vj]) {
				float dist_vivj = dist(_mesh.vertices[vi], _mesh.vertices[vj]);
				if ( _currentGeodesicDistance[vj] > _currentGeodesicDistance[vi] + dist_vivj ) {
					_currentGeodesicDistance[vj] = _currentGeodesicDistance[vi] + dist_vivj;
					_pq.push( vertexForDijkstra(vj, _currentGeodesicDistance[vj]) );
				}
			}
		}
	}

	return visitedVertices;
}


vector<int> GeodesicTraversal::traverseFaces(int sourceFace, float maxGeodesicDistance, int minNumFaces, float checkNormalCompatibility) {
	assert( _pq.empty() );

	vector<int> visitedFaces;

	_pq.push( vertexForDijkstra( sourceFace, 0.0f) );
	for (int r = 0; r < _mesh.faces.size(); r++) {
		_currentGeodesicDistance[r] = FLT_MAX;
		_visited[r] = false;
	}
	_currentGeodesicDistance[sourceFace] = 0.0f;

	while (!_pq.empty()) {
		vertexForDijkstra v = _pq.top();
		int vi = v.getVertexIndex();
		_pq.pop();
		if (_visited[vi]) {
			continue;
		}
		_visited[vi] = true;

		if ( ((_mesh.faces[vi].facenormal DOT _mesh.faces[sourceFace].facenormal ) > checkNormalCompatibility) || (visitedFaces.size() <= minNumFaces) ) { // reject faces with incompatible normals
			visitedFaces.push_back( vi );
		} else {
			continue;
		}
		if ( (v.getGeodDistance() > maxGeodesicDistance) && (visitedFaces.size() > minNumFaces) ) {
			continue;
		}

		for (int j = 0; j < 3; j++) {
			int vj = _mesh.across_edge[vi][j];
			if (!_visited[vj]) {
				float dist_vivj = _mesh.dist_faces_across_edge[vi][j];
//				float dist_vivj = dist(_mesh.faces[vi].faceCenter, _mesh.faces[vj].faceCenter);
				if ( _currentGeodesicDistance[vj] > _currentGeodesicDistance[vi] + dist_vivj ) {
					_currentGeodesicDistance[vj] = _currentGeodesicDistance[vi] + dist_vivj;
					_pq.push( vertexForDijkstra(vj, _currentGeodesicDistance[vj]) );
				}
			}
		}
	}

	return visitedFaces;
}


float GeodesicTraversal::getMaxGeodesicDistance(int subsample) {
	float maxGeodesicDistance = 0.0f; 

	for (int i = 0; i < _mesh.vertices.size(); i+=subsample) {
		//if (i % (subsample*PRINT_EVERY_N) == 0) {
		//	std::cout << 100.0f * (float)i / (float)_mesh.vertices.size() << "% complete\t\t\t\r";
		//	std::cout.flush();
		//}

		float maxGeodesicDistance_i = traverse(i);
		if (maxGeodesicDistance_i > maxGeodesicDistance) {
			maxGeodesicDistance = maxGeodesicDistance_i;
		}
	}
	//std::cout << 100.0f << "% complete\t\t\t\r";

	return maxGeodesicDistance / 2.0f;
}

float GeodesicTraversal::getMeanMaxGeodesicDistance(int subsample) {
	float meanMaxGeodesicDistance = 0.0f; 
	float samples = 0.0f; 

	for (int i = 0; i < _mesh.vertices.size(); i+=subsample) {
		//if (i % (subsample*PRINT_EVERY_N) == 0) {
		//	std::cout << 100.0f * (float)i / (float)_mesh.vertices.size() << "% complete\t\t\t\r";
		//	std::cout.flush();
		//}

		meanMaxGeodesicDistance += traverse(i);
		samples += 1;
	}
	meanMaxGeodesicDistance = meanMaxGeodesicDistance / samples;
	//std::cout << 100.0f << "% complete\t\t\t\r";

	return meanMaxGeodesicDistance / 2.0f;
}


float GeodesicTraversal::getMeanGeodesicDistance(int subsample) {
	float meanGeodesicDistance = 0.0f; 
	float samples = 0.0f;

	for (int i = 0; i < _mesh.vertices.size(); i+=subsample) {
		//if (i % (subsample*PRINT_EVERY_N) == 0) {
		//	std::cout << 100.0f * (float)i / (float)_mesh.vertices.size() << "% complete\t\t\t\r";
		//	std::cout.flush();
		//}

		float meanGeodesicDistance_i = 0.0f; 
		float samples_i = 0.0f;
		traverse(i);
		for (int j = 0; j < _mesh.vertices.size(); j+=subsample) {
			if ( (i != j) && (_currentGeodesicDistance[j] != FLT_MAX) ) {
				meanGeodesicDistance_i += _currentGeodesicDistance[j] * _mesh.pointareas[j]; 
				samples_i += _mesh.pointareas[j];
			}
		}
		meanGeodesicDistance += meanGeodesicDistance_i / samples_i;
		samples += 1;
	}
	meanGeodesicDistance = meanGeodesicDistance / samples;
	//std::cout << 100.0f << "% complete\t\t\t\r";

	return meanGeodesicDistance;
}


float GeodesicTraversal::getMedianGeodesicDistance(int subsample) {
	std::vector<float> medianGeodesicDistances;

	for (int i = 0; i < _mesh.vertices.size(); i+=subsample) {
		//if (i % (subsample*PRINT_EVERY_N) == 0) {
		//	std::cout << 100.0f * (float)i / (float)_mesh.vertices.size() << "% complete\t\t\t\r";
		//	std::cout.flush();
		//}
		vector<int> visitedVertices = traverse(i, FLT_MAX);
		medianGeodesicDistances.push_back( _currentGeodesicDistance[ visitedVertices[ _mesh.vertices.size()/2 ] ] );
	}
	//std::cout << 100.0f << "% complete\t\t\t\r";

	return median(medianGeodesicDistances);
}

float GeodesicTraversal::getPercentileGeodesicDistance(int subsample, float k, bool normalize) {
	std::vector<float> medianGeodesicDistances;

	for (int i = 0; i < _mesh.vertices.size(); i+=subsample) {
		//if (i % (subsample*PRINT_EVERY_N) == 0) {
		//	std::cout << 100.0f * (float)i / (float)_mesh.vertices.size() << "% complete\t\t\t\r";
		//	std::cout.flush();
		//}
		vector<int> visitedVertices = traverse(i, FLT_MAX);
		medianGeodesicDistances.push_back( _currentGeodesicDistance[ visitedVertices[ (int)floor(k * (float)_mesh.vertices.size()) ] ] );
	}
	//std::cout << 100.0f << "% complete\t\t\t\r";

	if (normalize)
		return percentile(medianGeodesicDistances, k) * (0.5f / k);
	else
		return percentile(medianGeodesicDistances, k);
}

};

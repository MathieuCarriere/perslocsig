#ifndef __DIJKSTRA_H
#define __DIJKSTRA_H

#include <queue>
#include <float.h>
#include <assert.h>

#include "TriMesh.h"
#include "feature.h"

#define PRINT_EVERY_N 3000

namespace trimesh{
class vertexForDijkstra {
private:
	float _di;
	int _vi;

public:
  vertexForDijkstra(int vi, float di): _vi(vi), _di(di) { }
  bool operator < (const vertexForDijkstra& v) const {
    return _di > v._di;
  }
  float getGeodDistance() { 
	  return _di;
  }
  int getVertexIndex() { 
	  return _vi;
  }
};

class GeodesicTraversal {
private: 
	TriMesh	_mesh;
	bool* _visited;
	float* _currentGeodesicDistance;
	std::priority_queue<vertexForDijkstra> _pq;
public:
	GeodesicTraversal(TriMesh& mesh);
	float traverse(int sourceVertex);
	std::vector<int>  traverse(int sourceVertex, float maxGeodesicDistance);
	std::vector<int>  traverseFaces(int sourceFace, float maxGeodesicDistance = FLT_MAX, int minFaces = 0, float checkNormalCompatibility = -FLT_MAX);
	float getMaxGeodesicDistance(int subsample = 1);	
	float getMeanMaxGeodesicDistance(int subsample = 1);	
	float getMeanGeodesicDistance(int subsample = 1);
	float getMedianGeodesicDistance(int subsample = 1);
	float getPercentileGeodesicDistance(int subsample = 1, float k = 0.3f, bool normalize = false);
	float* getGeodDistances();
	~GeodesicTraversal();
};

};

# endif

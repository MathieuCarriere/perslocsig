#pragma once

#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "Dijkstra.h"

#define EPSILON	1e-30f
#define SAMPLING_RATE_FOR_MESH_NORMALIZATION 5000

// Curvature parameters. 64 features.
// **********************************

#define NUM_CURVATURE_FEATURES_PER_SCALE 16
#define NUM_SCALES 4
#define SCALE1_DIV 50.0f
#define SCALE2_DIV 25.0f
#define SCALE3_DIV 10.0f
#define SCALE4_DIV 5.0f
const float SCALES[NUM_SCALES] = {(1.0f / SCALE1_DIV), (1.0f / SCALE2_DIV), (1.0f / SCALE3_DIV), (1.0f / SCALE4_DIV)};

// PCA parameters. 60 features.
// ****************************

#define NUM_PCA_FEATURES_PER_SCALE 12
#define NUM_PSCALES 5
#define PSCALE1_DIV 10.0f
#define PSCALE2_DIV 5.0f
#define PSCALE3_DIV 2.5f
#define PSCALE4_DIV 1.5f
#define PSCALE5_DIV 1.0f
const float PSCALES[NUM_PSCALES] = {(1.0f / PSCALE1_DIV), (1.0f / PSCALE2_DIV), (1.0f / PSCALE3_DIV), (1.0f / PSCALE4_DIV), (1.0f / PSCALE5_DIV)};

// Shape Context & Geodesic parameters. 270 + 18 = 288 features.
// *************************************************************

#define NUM_BINNING_TYPES_GD 6
#define NUM_FEATURES_GD 18 // mean, squared mean, tripled mean, quad mean, variance, skewness, kurtosis, percentiles
#define ACCELERATE_N2_COMPUTATIONS true
const int NUM_GD_BINS[ NUM_BINNING_TYPES_GD ] = {5, 5, 5, 10, 10, 10};
const int NUM_ANGLE_BINS[ NUM_BINNING_TYPES_GD ] = {4, 6, 8, 4, 6, 8};
const float GD_LOGP[ NUM_BINNING_TYPES_GD ] = {0.90f, 0.90f, 0.90f, 0.90f, 0.90f, 0.90f};
const float GD_MAX[ NUM_BINNING_TYPES_GD ] = {1.20f, 1.20f, 1.20f, 1.20f, 1.20f, 1.20f};

// SDF & VSI parameters. 72 + 24 = 96 features. 
// ********************************************

#define NUM_SDF_FEATURES_PER_BASE 4
#define NUM_SDF_FEATURES_PER_ANGLE 6
#define NUM_ANGLES_SDF 3
#define ANGLE1_SDF 15.0f
#define ANGLE2_SDF 30.0f
#define ANGLE3_SDF 45.0f
#define NUM_VSI_FEATURES 6
#define ANGLE_VSI 45.0f
const float ANGLES_SDF[NUM_ANGLES_SDF] = {ANGLE1_SDF, ANGLE2_SDF, ANGLE3_SDF};
#define ANGLE_STEP ANGLE2_SDF-ANGLE1_SDF
#define ROT_ANGLE_STEP 45.0f    // was 90
const float BASES_SDF[NUM_SDF_FEATURES_PER_BASE] = {1.0f, 2.0f, 4.0f, 8.0f}; // always start from 1, make sure the fiest element of BASES_SDF is 1 (no use of log then).

// Spin Image parameters. 100 features.
// ************************************

#define SPIN_DISTANCE_SUPPORT 2.0f
#define SPIN_RESOLUTION 10.0f
#define SPIN_BIN_SIZE (SPIN_DISTANCE_SUPPORT / SPIN_RESOLUTION)

// Contextual probability parameters. 35*NbLabels features.
// ********************************************************

#define NUM_GD_BINS_SCCLASSPROB	5
#define GD_LOGP_SCCLASSPROB 0.90f
#define GD_MAX_SCCLASSPROB 1.20f
#define ED_LOGP_SCCLASSPROB 0.90f
#define ED_MAX_SCCLASSPROB 1.20f
#define USE_PCADIR_PART_LABELS_MAX 1

// Edge parameters. 120 features.
// ******************************

#define NUM_EDGE_FEATURES 120
#define NUM_FACEDIHEDRALANGLE_FEATURES_PERSCALE	2
#define NUM_DIHEDRALANGLE_FEATURES_PERSCALE 10
#define NUM_EDGE_SCALES	4

#define NUM_EDGE_BINS_SCCLASSPROB 4
#define EDGE_LOGP_SCCLASSPROB 1.00f
#define EDGE_MAX_SCCLASSPROB 0.08f

using namespace std;

namespace trimesh {

inline bool isemptychar(char *c) {
	if (c == NULL)
		return true;

	int i = 0;
	while( c[i] != '\0' && c[i] != 0) {
		if ( isspace( c[i++] ) == false ) {
			return false;
		}
	}
	return true;
}

inline double* covar(double* D, unsigned int N) {

	double *CD = new double[9];

	CD[0] = 0.0f;
	for (unsigned int r = 0; r < N; r++)
		CD[0]  += D[r*3] * D[r*3];
	for (unsigned int j = 0; j <= 1; j++) {
		CD[1 + j*3] = 0.0f;
		for (unsigned int r = 0; r < N; r++) {
			CD[1 + j*3] +=  D[1 + r*3] * D[j + r*3];
		}
	}
	for (unsigned int j = 0; j <= 2; j++) {
		CD[2 + j*3] = 0.0f;
		for (unsigned int r = 0; r < N; r++) {
			CD[2 + j*3] +=  D[2 + r*3] * D[j + r*3];
		}
	}

	CD[3] = CD[1];
	CD[6] = CD[2];
	CD[7] = CD[5];

	return CD;
}

template <class T>
inline T acosd(T d) {
	d = (d<-1.0f)?-1.0f:d;
	d = (d>1.0f)?1.0f:d;
	return acos(d) * (180.0f / M_PI);
}

inline float getAngleFromUnitVectors(vec& v1, vec &v2 ) {
	float cosphi = v1 DOT v2;
	float sinphi = len(v1 CROSS v2);

	return (float)( atan2( sinphi, cosphi ) * (180.0f / M_PI) );
}

inline float getAngleFromUnitVectors(vec v1, vec v2 ) {
	float cosphi = v1 DOT v2;
	float sinphi = len(v1 CROSS v2);

	return (float)( atan2( sinphi, cosphi ) * (180.0f / M_PI) );
}


template <class T>
inline T logabs(T d) {
	return fabs( log( abs(d) + EPSILON ) );
}

template <class T>
inline T cosd(T d) {
	return cos( (d / 180.0f) * M_PI );
}

template <class T>
inline T sind(T d) {
	return sin( (d / 180.0f) * M_PI );
}


template <class T>
inline vector<T> removePositiveOutliers(vector<T>& v, float percentile) {
	vector<T> result( v.size() );
	for (int i = 0; i < v.size(); i++) {
		result[i] = fabs(v[i]);
	}
	std::sort( result.begin(), result.end() );
	float outlierValue = result[ (unsigned int)floor(percentile*(float)v.size()) ];
	for (int i = 0; i < v.size(); i++) {
		if (fabs(v[i]) < outlierValue) {
			result[i] = v[i];
		} else {
			result[i] = sgn(v[i]) * outlierValue;
		}
	}
	return result;
}

template <class T>
inline vector<T> removeOutliers(vector<T>& v, float percentile) {
	vector<T> result( v.size() );
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i];
	}
	std::sort( result.begin(), result.end() );
	float outlierValueMax = result[ (unsigned int)floor(percentile*(float)v.size()) ];
	float outlierValueMin = result[ (unsigned int)floor((1.0f-percentile)*(float)v.size()) ];
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < outlierValueMin) {
			result[i] = outlierValueMin;
		} else if (v[i] > outlierValueMax) {
			result[i] = outlierValueMax;
		} else {
			result[i] = v[i];
		}
	}
	return result;
}

template <class T>
inline vector<T> removeOutliers(vector<T>& v, float percentilePositive, float percentileNegative) {
	vector<T> result( v.size() );
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i];
	}
	std::sort( result.begin(), result.end() );
	float outlierValueMax = result[ (unsigned int)floor(percentilePositive*(float)v.size()) ];
	float outlierValueMin = result[ (unsigned int)floor(percentileNegative*(float)v.size()) ];
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < outlierValueMin) {
			result[i] = outlierValueMin;
		} else if (v[i] > outlierValueMax) {
			result[i] = outlierValueMax;
		} else {
			result[i] = v[i];
		}
	}
	return result;
}

template<class T> 
struct IndexCmp {
	T V;
	IndexCmp(T _V) : V(_V) {}
	bool operator()(const unsigned int a, const unsigned int b) const { 
		return V[a] < V[b]; 
	}
};


template <class T>
inline unsigned int* sortArrayIndices(T* A, int size) {
	vector<unsigned int> VI(size);
	for (unsigned int i = 0; i < size; i++) {
		VI[i] = i;
	}		
	vector<T> V(A, A+size);
	std::sort( VI.begin(), VI.end(), IndexCmp<vector<T>&>(V) );
	unsigned int *I = new unsigned int[size];
	for (unsigned int i = 0; i < size; i++) {
		I[i] = VI[i];
	}		
	return I;
};

template <class T>
inline std::vector<unsigned int> sortVectorIndices(vector<T>& V) {
	vector<unsigned int> VI(V.size());
	for (unsigned int i = 0; i < V.size(); i++) {
		VI[i] = i;
	}		
	std::sort( VI.begin(), VI.end(), IndexCmp<vector<T>&>(V) );
	return VI;
};

template <class T>
inline T maximum(vector<T>& v) {
	T result = -std::numeric_limits<T>::max();
	for (int i = 0; i < v.size(); i++) {
		if (v[i] > result) {
			result = v[i];
		}
	}
	return result;
}

template <class T>
inline T minimum(vector<T>& v) {
	T result = std::numeric_limits<T>::max(); 
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < result) {
			result = v[i];
		}
	}
	return result;
}


template <class T>
inline T mean(vector<T>& v) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += v[i];
	}
	result /= (T)v.size();
	return result;
}

template <class T>
inline T weightedmean(vector<T>& v, vector<T>& w) {
	T result = 0;
	T sumw = 0;
	for (int i = 0; i < v.size(); i++) {
		result += v[i] * w[i];
		sumw += w[i];
	}

	result /= sumw;
	return result;
}

template <class T>
inline T median(vector<T> v) {
	std::sort(v.begin(), v.end());
	if ( v.size() % 2 == 0 ) {
		return ( v[ v.size() / 2 ] +  v[ (v.size() / 2) - 1 ] ) / 2.0f;
	} else {
		return v[ v.size() / 2 ];
	}
}

template <class T>
inline T weightedmedian(vector<T>& v, vector<T>& w) {
        
	vector<unsigned int> VI = sortVectorIndices(v);
        
	double sumw = 0.0f, csumw = 0.0f;
	int i;
	for (i = 0; i < w.size(); i++) 
		sumw += (double)w[i];
	
	for (i = 0; i < w.size(); i++) {
		csumw += (double)w[ VI[i] ] / (double)sumw;
		if (csumw >= .5f)
			break;
	}
        
	return v[ VI[i] ];
}

template <class T>
inline T percentile(vector<T> v, float k) {
	std::sort(v.begin(), v.end());
	return v[ (int)floor( k * (float)v.size() ) ];
}


template <class T>
inline T mean2(vector<T>& v) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += v[i] * v[i];
	}
	result /= (T)v.size();
	return result;
}

template <class T>
inline T mean4(vector<T>& v) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += v[i] * v[i] * v[i] * v[i];
	}
	result /= (T)v.size();
	return result;
}

template <class T>
inline T meansqrt(vector<T>& v) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += sqrt(v[i]);
	}
	result /= (T)v.size();
	return result;
}

template <class T>
inline T meansqrtsqrt(vector<T>& v) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += sqrt(sqrt(v[i]));
	}
	result /= (T)v.size();
	return result;
}


template <class T>
inline T variance(vector<T>& v, T meanValue) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += (v[i] - meanValue) * (v[i] - meanValue);
	}
	result /= (T)v.size();
	return result;
}

template <class T>
inline T skewness(vector<T>& v, T meanValue, T varianceValue) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += (v[i] - meanValue) * (v[i] - meanValue) * (v[i] - meanValue);
	}
	result /= (T)v.size();
	result /= sqrt( varianceValue*varianceValue*varianceValue + EPSILON);
	return result;
}

template <class T>
inline T kurtosis(vector<T>& v, T meanValue, T varianceValue) {
	T result = 0;
	for (int i = 0; i < v.size(); i++) {
		result += (v[i] - meanValue) * (v[i] - meanValue) * (v[i] - meanValue) * (v[i] - meanValue);
	}
	result /= (T)v.size();
	result /= ( varianceValue*varianceValue  + EPSILON);
	return result;
}

template <class T>
inline std::vector<T> unique(vector<T>& v) {
	std::vector<T> ret;
	ret.push_back( v[0] );
	for (int i = 1; i < v.size(); i++) {
		bool isunique = true;
		for (int j = 0; j < ret.size(); j++) {
			if (v[i] == ret[j]) {
				isunique = false;
				break;
			}
		}
		if (isunique == true) {
			ret.push_back( v[i] );
		}
	}

	return ret;
}





class Bin2D {
public:
	float x1;
	float x2;
	float y1;
	float y2;
	float cx, cy;
	float value;
	Bin2D() { value = 0.0f; }
	Bin2D(float startx, float endx, float starty, float endy): x1(startx), x2(endx), y1(starty), y2(endy) { value = 0.0f; }
	bool contains(float x, float y) { return (x >= x1) & (x <= x2) & (y >= y1) & (y <= y2); }
	void increaseValue(float v) { value += v; }
	void reset() { value = 0.0f; }
	void setCenter() { cx = (x1 + x2) / 2.0f; cy = (y1 + y2) / 2.0f; }
};

static inline std::ostream &operator << (std::ostream &os, const Bin2D &b) { 
	os << "x=[" << b.x1 << ", " << b.x2 << "], y=[" << b.y1 << ", " << b.y2 << "], value=" << b.value << std::endl;
	return os;
}

class FeatureSet {
public:
	float** FEATURES;
	int numFeatures;
	float** FEATURES2;
	int numFeatures2;
	int *LABELS;
	float* WEIGHTS;

	FeatureSet(): FEATURES(NULL), numFeatures(0), FEATURES2(NULL), numFeatures2(0), LABELS(NULL), WEIGHTS(NULL) {};
	~FeatureSet() {
		if (FEATURES != NULL) {
			for (int i = 0; i < numFeatures; i++) {
				delete[] FEATURES[i];
			}
			delete[] FEATURES;
		}
		if (FEATURES2 != NULL) {
			for (int i = 0; i < numFeatures2; i++) {
				delete[] FEATURES2[i];
			}
			delete[] FEATURES2;
		}
		if (LABELS != NULL)
			delete[] LABELS;
		if (WEIGHTS != NULL)
			delete[] WEIGHTS;
	}
};



// this is old code, replace with CGAL's AABB for large mesh datasets
inline float checkClosestLineIntersection(vec& line, point& origin, int f, TriMesh* mesh, bool *done, bool excludeCurrentFace = true) {

	// line is assumed to be normalized vector
	for (int i = 0; i < mesh->faces.size(); i++) {
		done[i] = false;
	}
	done[f] = excludeCurrentFace;
	float bboxdiagonal = len(mesh->bbox.max - mesh->bbox.min);		
	float currentDist = bboxdiagonal;

	for (float checkMaxDist = .2 * bboxdiagonal; checkMaxDist <= bboxdiagonal; checkMaxDist +=  .2 * bboxdiagonal) {
		float checkMaxDist2 = checkMaxDist*checkMaxDist;

		for (int i = 0; i < mesh->faces.size(); i++) {
			if (done[i] == true) {
				continue;
			}
			if ( len2(mesh->faces[i].faceCenter - origin) > checkMaxDist2) {
				continue;
			} else {
				done[i] = true;
			}
			//if ( (line DOT mesh->faces[i].facenormal) < 0) {
			//	continue;
			//}

			vec e1 = mesh->vertices[mesh->faces[i][1]] - mesh->vertices[mesh->faces[i][0]];
			vec e2 = mesh->vertices[mesh->faces[i][2]] - mesh->vertices[mesh->faces[i][0]];
			vec p = line CROSS e2;
			float tmp = p DOT e1;

			if ( tmp < EPSILON && tmp > -EPSILON) 
				continue;

			tmp = 1.0f / tmp;
			vec s = origin - mesh->vertices[mesh->faces[i][0]];
			float u = tmp * (s DOT p);
			if (u < 0.0 || u > 1.0)
				continue;

			vec q = s CROSS e1;
			float v = tmp * (line DOT q);
			if (v < 0.0 || v > 1.0 || u + v > 1.0)
				continue;


			float t = tmp * (e2 DOT q);		
			if (t <= 0) 
				continue;

			if (t < currentDist)
				currentDist = t;
		}

		if (currentDist != bboxdiagonal) {
			return currentDist;
		}
	}

	return currentDist;
}

void writeDebugInfoToFile( TriMesh *m, std::string c, FeatureSet* features, int numFaces, int pos, int fid );
TriMesh* processMesh( TriMesh *m, std::string id, bool writeDebugInfo, std::string USE_NORMALIZATION_FEATURE  );
FeatureSet* exportCurvatureFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numFaces, bool writeDebugInfo, bool returnNumFeaturesOnly );
FeatureSet* exportPCAFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numFaces, bool writeDebugInfo, bool returnNumFeaturesOnly );
FeatureSet* exportSCFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numFaces, float minf, float maxf, bool writeDebugInfo, bool returnNumFeaturesOnly);
FeatureSet* exportSDFFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numFaces, bool writeDebugInfo, bool returnNumFeaturesOnly);
FeatureSet* exportSpinImageFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numFaces, bool writeDebugInfo, bool returnNumFeaturesOnly);
FeatureSet* exportKSFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numFaces, bool writeDebugInfo, bool returnNumFeaturesOnly);
FeatureSet* exportSCClassProbFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numFaces, float** PH, int numlabels, bool writeDebugInfo, bool returnNumFeaturesOnly);

// Warning: m need to have its SDF features computed
FeatureSet* exportEdgeFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int numEdges, std::vector<int> labels, float boundaryWeight, bool writeDebugInfo, bool returnNumFeaturesOnly );
FeatureSet* exportSCEdgeFeatures(TriMesh* m, std::string id, FeatureSet* features, int pos, int fpos, float** PH, bool writeDebugInfo );

};


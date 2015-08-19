#include "TriMesh.h"

#include "../../clapack/INCLUDE/blaswrap.h"
#include "../../clapack/INCLUDE/clapack.h"
#include "../../clapack/INCLUDE/f2c.h"

namespace trimesh{

class waPCA {	
public:
	std::vector<vec> alignedDirections; 
	std::vector<float> maxDistance;

	waPCA( std::vector<point>& points, std::vector<float>& weights, std::vector<point>& allpoints ) {
		if (points.size() < 10) {
			alignedDirections.push_back( vec(0,0,0) );
			alignedDirections.push_back( vec(0,0,0) );
			alignedDirections.push_back( vec(0,0,0) );
			maxDistance.push_back(0);
			maxDistance.push_back(0);
			maxDistance.push_back(0);
			return;
		}
		char JOBZ = 'V';
		char UPLO = 'L';
		integer N = 3; 
		double EIG[3]; // eigenvalues in ascending order
		double WORK[5 * 8];
		integer LWORK = 5 * 8;
		integer INFO;

		double *V = new double[ points.size() * 3 ];
		float totalweight = 0.0f;
	
		float M[3] = {0.0f, 0.0f, 0.0f};
		for (int j = 0; j < points.size(); j++) {
			float x = points[j][0] * weights[j];
			float y = points[j][1] * weights[j];
			float z = points[j][2] * weights[j];
			M[ 0 ] += x;
			M[ 1 ] += y;
			M[ 2 ] += z;
			totalweight += weights[j];
		}
		M[0] /= totalweight;
		M[1] /= totalweight;
		M[2] /= totalweight;
		//float M[3] = {0.0f, 0.0f, 0.0f};
		//std::vector<float> xs, ys, zs;
		//for (int j = 0; j < points.size(); j++) {
		//	xs.push_back(points[j][0]);
		//	ys.push_back(points[j][1]); 
		//	zs.push_back(points[j][2]);  
		//}
		//M[0] = weightedmedian( xs, weights );
		//M[1] = weightedmedian( ys, weights );
		//M[2] = weightedmedian( zs, weights );
		//xs.clear(); ys.clear(); zs.clear();


		for (int j = 0; j < points.size(); j++) {
			float w = sqrt( weights[j] );
			float x = (points[j][0] - M[0]) * w;
			float y = (points[j][1] - M[1]) * w;
			float z = (points[j][2] - M[2]) * w;
			V[ j * 3 ] = x;
			V[ 1 + j * 3 ] = y;
			V[ 2 + j * 3 ] = z;
		}

		double* COV = covar(V, points.size());
		for (int j = 0; j < 9; j++) {
			COV[j] /= totalweight;
		}
		dsyev_( &JOBZ, &UPLO, &N, COV, &N, EIG, WORK, &LWORK, &INFO );
		if (INFO != 0) {
			std::cerr << std::endl << "Part of label possibly did not acquire right pca vales during dsyev_ call! INFO = " << INFO << std::endl;
		}
		vector<vec> pcaDirections(3); 
		pcaDirections[0] = vec( COV[0], COV[1], COV[2] );
		pcaDirections[1] = vec( COV[3], COV[4], COV[5] );
		pcaDirections[2] = vec( COV[6], COV[7], COV[8] );
		float directionVariation[3] = {0.0f, 0.0f, 0.0f};
		delete[] V;
		delete[] COV;

		std::vector<float> p0s, p1s, p2s;
		for (int j = 0; j < points.size(); j++) {
			float x = (points[j][0] - M[0]);
			float y = (points[j][1] - M[1]);
			float z = (points[j][2] - M[2]);
			float p0 = x * pcaDirections[0][0] + y * pcaDirections[0][1] + z * pcaDirections[0][2];
			float p1 = x * pcaDirections[1][0] + y * pcaDirections[1][1] + z * pcaDirections[1][2];
			float p2 = x * pcaDirections[2][0] + y * pcaDirections[2][1] + z * pcaDirections[2][2];

			//directionVariation[0] += fabs(p0) * weights[j];
			//directionVariation[1] += fabs(p1) * weights[j];
			//directionVariation[2] += fabs(p2) * weights[j];
			p0s.push_back(fabs(p0));
			p1s.push_back(fabs(p1)); 
			p2s.push_back(fabs(p2));  
		}

		//directionVariation[0] /= totalweight;
		//directionVariation[1] /= totalweight;
		//directionVariation[2] /= totalweight;
		directionVariation[0] = weightedmedian(p0s, weights);   // CCY median with weightedmedian (weights are face areas)
		directionVariation[1] = weightedmedian(p1s, weights);
		directionVariation[2] = weightedmedian(p2s, weights);
		p0s.clear(); 
		p1s.clear(); 
		p2s.clear(); 

		unsigned int *si = sortArrayIndices(directionVariation, 3);		
		alignedDirections.push_back( pcaDirections[ si[0] ] );
		alignedDirections.push_back( pcaDirections[ si[1] ] );
		alignedDirections.push_back( pcaDirections[ si[2] ] );		
		delete[] si;


		float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX, maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;
		for (int j = 0; j < allpoints.size(); j++) {
			float x = (allpoints[j][0] - M[0]);
			float y = (allpoints[j][1] - M[1]);
			float z = (allpoints[j][2] - M[2]);
			float p0 = x * alignedDirections[0][0] + y * alignedDirections[0][1] + z * alignedDirections[0][2];
			float p1 = x * alignedDirections[1][0] + y * alignedDirections[1][1] + z * alignedDirections[1][2];
			float p2 = x * alignedDirections[2][0] + y * alignedDirections[2][1] + z * alignedDirections[2][2];

			if ( p0 < minx ) {
				minx = p0;
			}
			if ( p1 < miny ) {
				miny = p1;
			}
			if ( p2 < minz ) {
				minz = p2;
			}
			if ( p0 > maxx ) {
				maxx = p0;
			}
			if ( p1 > maxy ) {
				maxy = p1;
			}
			if ( p2 > maxz ) {
				maxz = p2;
			}
		}
		maxDistance.resize(3);
		maxDistance[0] = maxx - minx;
		maxDistance[1] = maxy - miny;
		maxDistance[2] = maxz - minz;	
	}
};

};

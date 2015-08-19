/*
Szymon Rusinkiewicz
Princeton University

TriMesh_normals.cc
Compute per-vertex normals for TriMeshes

For meshes, uses average of per-face normals, weighted according to:
  Max, N.
  "Weights for Computing Vertex Normals from Facet Normals,"
  Journal of Graphics Tools, Vol. 4, No. 2, 1999.

For raw point clouds, fits plane to k nearest neighbors.
*/

#include "TriMesh.h"
#include "KDtree.h"
#include "lineqn.h"
using namespace std;


namespace trimesh {

// Compute per-vertex normals
void TriMesh::need_normals()
{
	// Nothing to do if we already have normals
	int nv = vertices.size();
	if (int(normals.size()) == nv)
		return;

	//dprintf("Computing normals... ");
	normals.clear();
	normals.resize(nv);
        
        int nf = faces.size();

	if (need_faces(), !faces.empty()) {
		// Compute from faces
#pragma omp parallel for
		for (int i = 0; i < nf; i++) {
			const point &p0 = vertices[faces[i][0]];
			const point &p1 = vertices[faces[i][1]];
			const point &p2 = vertices[faces[i][2]];
			vec a = p0-p1, b = p1-p2, c = p2-p0;
			float l2a = len2(a), l2b = len2(b), l2c = len2(c);
			if (!l2a || !l2b || !l2c)
				continue;
			vec facenormal = a CROSS b;
			normals[faces[i][0]] += facenormal * (1.0f / (l2a * l2c));
			normals[faces[i][1]] += facenormal * (1.0f / (l2b * l2a));
			normals[faces[i][2]] += facenormal * (1.0f / (l2c * l2b));
                        vec bias(1e-37f, 1e-37f, 1e-37f);
                        faces[i].facenormal = (a CROSS b) + bias;
			faces[i].faceArea = .5f * len( faces[i].facenormal );
                        faces[i].faceCenter = (p0+p1+p2)/3.0f;
		}
	} else if (!tstrips.empty()) {
		// Compute from tstrips
		const int *t = &tstrips[0], *end = t + tstrips.size();
		while (likely(t < end)) {
			int striplen = *t - 2;
			t += 3;
			bool flip = false;
			for (int i = 0; i < striplen; i++, t++, flip = !flip) {
				const point &p0 = vertices[*(t-2)];
				const point &p1 = vertices[*(t-1)];
				const point &p2 = vertices[* t   ];
				vec a = p0-p1, b = p1-p2, c = p2-p0;
				float l2a = len2(a), l2b = len2(b), l2c = len2(c);
				if (!l2a || !l2b || !l2c)
					continue;
				vec facenormal = flip ? (b CROSS a) : (a CROSS b);
				normals[*(t-2)] += facenormal * (1.0f / (l2a * l2c));
				normals[*(t-1)] += facenormal * (1.0f / (l2b * l2a));
				normals[* t   ] += facenormal * (1.0f / (l2c * l2b));
			}
		}
	} else {
		// Find normals of a point cloud
		const int k = 6;
		const vec ref(0, 0, 1);
		KDtree kd(vertices);
#pragma omp parallel for
		for (int i = 0; i < nv; i++) {
			vector<const float *> knn;
			kd.find_k_closest_to_pt(knn, k, vertices[i]);
			int actual_k = knn.size();
			if (actual_k < 3) {
				dprintf("Warning: not enough points for vertex %d\n", i);
				normals[i] = ref;
				continue;
			}
			// Compute covariance
			float C[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
			// The below loop starts at 1, since element 0     
			// is just vertices[i] itself 
			for (int j = 1; j < actual_k; j++) {
				vec d = point(knn[j]) - vertices[i];
				for (int l = 0; l < 3; l++)
					for (int m = 0; m < 3; m++)
						C[l][m] += d[l] * d[m];
			}
			float e[3];
			eigdc<float,3>(C, e);
			normals[i] = vec(C[0][0], C[1][0], C[2][0]);
			if ((normals[i] DOT ref) < 0.0f)
				normals[i] = -normals[i];
		}
	}

	// Make them all unit-length
#pragma omp parallel for
	for (int i = 0; i < nv; i++)
		normalize(normals[i]);

        totalFaceArea = 0;
#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
                vec vec0 = normals[faces[i][0]];
                vec vec1 = normals[faces[i][1]];
                vec vec2 = normals[faces[i][2]];
                Vec<3,float> sumv = vec0+vec1+vec2;
		faces[i].facenormal = normalize( sumv );
		totalFaceArea += faces[i].faceArea;
	}

	//dprintf("Done.\n");
}

}; // namespace trimesh

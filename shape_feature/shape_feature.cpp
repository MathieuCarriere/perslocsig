#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

#include "TriMesh.h"
#include "Dijkstra.h"

void read_labels(int t_int, std::vector<int>* t_labels_faces);

int main(int argc, char** argv){
    
    for (int i = 1; i < argc; i++){
    
    string name = argv[i];
    size_t len = name.size();
    string id = name.substr(0,len-4);
    int val = atoi(id.c_str());
    
    cout << "*********************************************************" << endl;
    cout << "Mesh" << " : " << id << " " << endl << endl;
    
    char dist[100];
    sprintf(dist, "dist_matrix_%s",(char*) id.c_str());
    ifstream distfile(dist);
    
    trimesh::TriMesh* m = trimesh::TriMesh::read((char*) name.c_str());
    trimesh::GeodesicTraversal Geod(*m);
    
    int num_vertices = m->vertices.size();
    int num_faces = m->faces.size();
    
    if (!distfile){
        ofstream Distfile(dist, ios::out | ios::binary);
        cout << endl << "Computing geodesics..." << endl;
        for (int i = 0; i < num_vertices; i++){
            if (i%1000 == 0)
                cout << "    source point = " << i << endl;
            Geod.traverse(i);
            for (int j = 0; j < num_vertices; j++){
                double d = Geod.getGeodDistances()[j];
                Distfile.write((char*) &d, 8);
            }
        }
        cout << "Geodesics computed." << endl;
    }
    else
        cout << "Geodesics already computed." << endl;
    
    
    trimesh::FeatureSet *features = NULL;
    
    trimesh::processMesh(m,id,true,"PERCENTILE_GEOD_DISTANCE");
    features = trimesh::exportCurvatureFeatures(m,id,features,0,num_faces,true,false);
    features = trimesh::exportPCAFeatures(m,id,features,0,num_faces,true,false);
    features = trimesh::exportSCFeatures(m,id,features,0,num_faces,-90,90,true,false);
    features = trimesh::exportSDFFeatures(m,id,features,0,num_faces,true,false);
    features = trimesh::exportSpinImageFeatures(m,id,features,0,num_faces,true,false);
    
    
    if (val != 0){
       
        std::vector<int> lb(num_faces);
        read_labels(val,&lb);
        
        m->need_across_edge();
    
        float boundariesPercentage = 0;
        for (int j = 0; j < m->faces.size(); j++) {
                for (int k = 0; k < 3; k++) {
                        boundariesPercentage += ( lb[j] != lb[ m->across_edge[j][k] ] );
                }
        }
        boundariesPercentage /= 3.0f * (float) num_faces;
        
        features = trimesh::exportEdgeFeatures(m,id,features,0,3*num_faces,lb,1.0f/boundariesPercentage,true,false);
        
        cout << endl;
    }
    
    
    cout << endl;
    
    }
    
    return 0;
    
}

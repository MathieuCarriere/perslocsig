#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <cmath>
#include <math.h>
#include <omp.h>

using namespace std;

bool cmp(double a, double b){return a > b;}

int main(int argc, char** argv) {

	string PersDiag = argv[1];

	cout << "Computing signature..." << endl;
    
	size_t len = PersDiag.size(); string id = PersDiag.substr(33,len); char* f = &PersDiag[19]; char* p = &PersDiag[23];

	char sign[100]; double x; int max_h = 0;

	
	if ((*p)==*("g"))
		if ((*f)==*("e"))
	    		sprintf(sign,"signature_ecc_points_%s", (char*) id.c_str());
	if ((*p)==*("l"))
		if ((*f)==*("g"))
	    		sprintf(sign,"signature_geo_points_%s", (char*) id.c_str());
        
        
	if ((*f)==*("h"))
		sprintf(sign,"signature_hks_points_%s", (char*) id.c_str());

	ofstream SIGN(sign, ios::out | ios::binary); ifstream PD((char*) PersDiag.c_str()); string line;

	while(getline(PD,line)){

		vector<double> v; v.clear();
		stringstream stream(line);
        
        	while(stream >> x)
            		v.push_back(x);
        
        	int num = v.size()/2;

		if (num*(num-1)/2+1 > max_h)
			max_h = min(num*(num-1)/2+1,1000);

	}

	int k = 0;
    
        SIGN.write((char*) &max_h, 4); PD.clear(); PD.seekg(0);
    
    	while(getline(PD,line)){
    
        	 vector<double> v; v.clear();
        	 vector<double> cm; cm.clear();
        
        	 stringstream stream(line);

		if (k%1000 == 0)
			cout << "    source = " << k << endl;
        
        	 while(stream >> x)
            		 v.push_back(x);
        
        	 int num = v.size()/2;
        
        	 if (num != 0){
        		for (int i = 0; i < num-1; i++){
            			for(int j = i+1; j < num; j++){
                
                		double v11 = v[2*i]; double v12 = v[2*i+1];
                		double v21 = v[2*j]; double v22 = v[2*j+1];
                
                		double p1 = (log(1+v12)-log(1+v11))/2;
                		double p2 = (log(1+v22)-log(1+v21))/2;
                		double d = max(abs(log(1+v11)-log(1+v21)),abs(log(1+v12)-log(1+v22)));

                		cm.push_back(min(d,min(p1,p2)));
                
            			}
        		}
        
        		sort(cm.begin(), cm.end(), cmp);
                        cm.insert(cm.begin(),log(1+v[1]));

			int n = cm.size();
			if(n <= max_h){
			for (int j = 0; j < n; j++){
				x = cm[j];
    				SIGN.write((char*) &x, 8);
			}
			for (unsigned int j = n; j < max_h; j++){
    				x = 0;
    				SIGN.write((char*) &x, 8);
			}}
			else{
			for (int j = 0; j < max_h; j++){
				x = cm[j];
    				SIGN.write((char*) &x, 8);
			}}

        		k += 1;
	        }

	        else{

			cm.push_back(log(1+v[0]));

		        int n = cm.size();
			for (int j = 0; j < n; j++){
				x = cm[j];
				SIGN.write((char*) &x, 8);
			}
			for (unsigned int j = n; j < max_h; j++){
				x = 0;
				SIGN.write((char*) &x, 8);
			}

	       }
	}
    	
    
    return 0;
}

/* 
 * File:   union_find.h
 * Author: mathieu
 *
 * Created on January 22, 2015, 3:37 PM
 */

#ifndef UNION_FIND_H
#define	UNION_FIND_H

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

using namespace std;

int UF_find(int i, vector<int>* parents){

    if ((*parents)[i] != i){
        (*parents)[i] = UF_find((*parents)[i], parents);
        return (*parents)[i];
    }

    return i;
    
}

void UF_union(int i, int j, int v, vector<pair<int,int> >* entries, vector<int>* parents, vector<int>* v_vect){

    int x = UF_find(i, parents);
    int y = UF_find(j, parents);
    
    int ind_x, ind_y;
    
    for (int i = 0; i < (*entries).size(); i++)
        if ((*entries)[i].first == x){
            ind_x = i;
            break;
        }
    
    for (int i = 0; i < (*entries).size(); i++)
        if ((*entries)[i].first == y){
            ind_y = i;
            break;
        }
            
    
    if (x != y){
    
        if ((*entries)[ind_x].second > (*entries)[ind_y].second){
            (*parents)[y] = x;
            (*entries)[ind_x].second += (*entries)[ind_y].second;
            (*v_vect)[ind_x] = v;
            (*entries).erase((*entries).begin()+ind_y);
            (*v_vect).erase((*v_vect).begin()+ind_y);
        }
        else{
            (*parents)[x] = y;
            (*entries)[ind_y].second += (*entries)[ind_x].second;
            (*v_vect)[ind_y] = v;
            (*entries).erase((*entries).begin()+ind_x);
            (*v_vect).erase((*v_vect).begin()+ind_x);
        }
    }
    
}

#endif	/* UNION_FIND_H */


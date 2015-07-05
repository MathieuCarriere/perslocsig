#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

void read_labels(int t_int, std::vector<int>* t_labels_faces){

    char t_lab[100];
    sprintf(t_lab,"%d_labels.txt",t_int);
    //sprintf(t_lab,"%d.lab",t_int);
    
    std::cout << std::endl << "Reading labels... ";
    
    std::ifstream T_LAB(t_lab);
    std::string line;
   
    int k = 0;
    int lab,n;
    
    while(std::getline(T_LAB,line)){
        
        if (k % 2 == 0){
            
            // HUMAN.
            
            if (t_int >= 1 && t_int <= 20){
            
                char* c0 = &line[0];
                char* c1 = &line[1];
                char* c5 = &line[5];
            
                if ((*c0 == *("l")) && (*c5 == *("l")))
                        lab = 0;
                if ((*c0 == *("t")))
                        lab = 1;
                if ((*c0 == *("u")) && (*c5 == *("l")))
                        lab = 2;
                if ((*c0 == *("f")))
                        lab = 3;
                if ((*c0 == *("h")) && (*c1 == *("e")))
                        lab = 4;
                if ((*c0 == *("l")) && (*c5 == *("a")))
                        lab = 5;
                if ((*c0 == *("u")) && (*c5 == *("a")))
                        lab = 6;
                if ((*c0 == *("h")) && (*c1 == *("a")))
                        lab = 7;
            }
            
            // CUP.
            
            if (t_int >= 21 && t_int <= 40){
            
                char* c0 = &line[0];
         
                if ((*c0 == *("b")))
                        lab = 0;
                if ((*c0 == *("h")))
                        lab = 1;
            }
            
            // ANT.
            
            if (t_int >= 41 && t_int <= 60){
            
                char* c0 = &line[0];
            
                if (*c0 == *("l"))
                        lab = 1;
                if (*c0 == *("s"))
                        lab = 2;
                if (*c0 == *("m"))
                        lab = 0;
                
            }
            
            // AIRPLANE.
            
            if (t_int >= 61 && t_int <= 80){
            
                char* c0 = &line[0];
            
                if (*c0 == *("b"))
                        lab = 1;
                if (*c0 == *("s"))
                        lab = 3;
                if (*c0 == *("w"))
                        lab = 0;
                if (*c0 == *("r"))
                        lab = 2;
            }
            
            // ANT.
            
            if (t_int >= 81 && t_int <= 100){
            
                char* c0 = &line[0];
            
                if (*c0 == *("b"))
                        lab = 2;
                if (*c0 == *("t"))
                        lab = 3;
                if (*c0 == *("h"))
                        lab = 0;
                if (*c0 == *("a"))
                        lab = 4;
                if (*c0 == *("l"))
                        lab = 1;
            }
            
            // CHAIR.
            
            if (t_int >= 101 && t_int <= 120){
            
                char* c0 = &line[0];
            
                if (*c0 == *("s"))
                        lab = 2;
                if (*c0 == *("l"))
                        lab = 3;
                if (*c0 == *("a"))
                        lab = 0;
                if (*c0 == *("b"))
                        lab = 1;
            }
            
            // OCTOPUS.
            
            if (t_int >= 121 && t_int <= 140){
            
                char* c0 = &line[0];
            
                if (*c0 == *("b"))
                        lab = 0;
                if (*c0 == *("l"))
                        lab = 1;
            }
            
            // TABLE.
            
            if (t_int >= 141 && t_int <= 160){
            
                char* c0 = &line[0];
            
                if (*c0 == *("s"))
                        lab = 0;
                if (*c0 == *("l"))
                        lab = 1;
            }
            
            // TEDDY.
            
            if (t_int >= 161 && t_int <= 180){
            
                char* c0 = &line[0];
                char* c1 = &line[1];
            
                if (*c0 == *("t"))
                        lab = 0;
                if (*c0 == *("h")){
                    if (*c1 == *("a"))
                        lab = 1;
                    if (*c1 == *("e"))
                        lab = 2;
                }
                if (*c0 == *("l"))
                        lab = 3;
                if (*c0 == *("e"))
                        lab = 4;
            }
            
            // HAND.
            
            if (t_int >= 181 && t_int <= 200){
            
                char* c0 = &line[0];
                char* c6 = &line[6];
            
                if (*c0 == *("t"))
                        lab = 0;
                if (*c0 == *("h"))
                        lab = 1;
                if (*c0 == *("f")){
                if (*c6 == *("1"))
                        lab = 2;
                if (*c6 == *("2"))
                        lab = 3;
                if (*c6 == *("3"))
                        lab = 4;
                if (*c6 == *("4"))
                        lab = 5;
                }
            }
            
            // PLIER.
            
            if (t_int >= 201 && t_int <= 220){
            
                char* c0 = &line[0];             
            
                if (*c0 == *("c"))
                        lab = 0;
                if (*c0 == *("h"))
                        lab = 1;
                if (*c0 == *("e"))
                        lab = 2;
                
            }
            
            // FISH.
            
            if (t_int >= 221 && t_int <= 240){
            
                char* c0 = &line[0];
            
                if (*c0 == *("b"))
                        lab = 0;
                if (*c0 == *("f"))
                        lab = 1;
                if (*c0 == *("t"))
                        lab = 2;
                
            }
            
            // BIRD.
            
            if (t_int >= 241 && t_int <= 260){
            
                char* c0 = &line[0];
            
                if (*c0 == *("b"))
                        lab = 1;
                if (*c0 == *("h"))
                        lab = 3;
                if (*c0 == *("w"))
                        lab = 0;
                if (*c0 == *("t"))
                        lab = 2;
            }
            
            // ARMADILLO.
            
            if (t_int >= 281 && t_int <= 300){
            
                char* c0 = &line[0];
                char* c1 = &line[1];
                char* c5 = &line[5];
            
                if ((*c0 == *("l")) && (*c5 == *("L")))
                        lab = 0;
                if ((*c0 == *("t"))){
                    if ((*c1 == *("o")))
                        lab = 1;
                    if ((*c1 == *("a")))
                        lab = 2;
                }
                if ((*c0 == *("u")) && (*c5 == *("L")))
                        lab = 4;
                if ((*c0 == *("f")))
                        lab = 3;
                if ((*c0 == *("b")))
                        lab = 5;
                if ((*c0 == *("e")))
                        lab = 6;
                if ((*c0 == *("h")) && (*c1 == *("e")))
                        lab = 7;
                if ((*c0 == *("l")) && (*c5 == *("A")))
                        lab = 8;
                if ((*c0 == *("u")) && (*c5 == *("A")))
                        lab = 9;
                if ((*c0 == *("h")) && (*c1 == *("a")))
                        lab = 7;
            }
            
            // BUST.
            
            if (t_int >= 301 && t_int <= 320){
            
                char* c0 = &line[0];
                char* c1 = &line[1];
            
                if (*c0 == *("h"))
                        lab = 0;
                if (*c0 == *("b"))
                        lab = 1;
                if (*c0 == *("e"))
                        lab = 2;
                if (*c0 == *("f"))
                        lab = 3;
                if (*c0 == *("n")){
                    if (*c1 == *("o"))
                        lab = 4;
                    if (*c1 == *("e"))
                        lab = 5;
                    
                }
            }
            
            // MECH.
            
            if (t_int >= 321 && t_int <= 340){
            
                char* c0 = &line[0];
            
                if (*c0 == *("b"))
                        lab = 0;
                if (*c0 == *("s"))
                        lab = 1;

            }
            
            // BEARING.
            
            if (t_int >= 341 && t_int <= 360){
            
                char* c0 = &line[0];
                char* c4 = &line[4];
            
                if (*c0 == *("b")){
                    if (*c4 == *("o"))
                        lab = 0;
                    if (*c4 == *("1"))
                        lab = 1;
                    if (*c4 == *("2"))
                        lab = 2;
                }
                if (*c0 == *("s"))
                        lab = 3;
                if (*c0 == *("t"))
                        lab = 4;
                if (*c0 == *("m"))
                        lab = 5;

            }
            
            // VASE.
            
            if (t_int >= 361 && t_int <= 380){
            
                char* c0 = &line[0];
                char* c1 = &line[1];
            
                if (*c0 == *("b")){
                    if (*c1 == *("o"))
                        lab = 0;
                    if (*c1 == *("a"))
                        lab = 1;
                }
                if (*c0 == *("s"))
                        lab = 3;
                if (*c0 == *("t"))
                        lab = 4;
                if (*c0 == *("h"))
                        lab = 2;

            }
            
            // FOURLEGS.
            
            if (t_int >= 381 && t_int <= 400){
            
                char* c0 = &line[0];
                char* c1 = &line[1];
                
                if (*c0 == *("t")){
                    if (*c1 == *("o"))
                        lab = 0;
                    if (*c1 == *("a"))
                        lab = 2;
                }
                if (*c0 == *("e"))
                        lab = 5;
                if (*c0 == *("l"))
                        lab = 1;
                if (*c0 == *("n"))
                        lab = 3;
                if (*c0 == *("h"))
                        lab = 4;
            }
            
            
        }
        
        std::stringstream stream(line);
        
        if(k % 2 == 1)
            while(stream >> n)
                (*t_labels_faces)[n-1] = lab;
        
        k += 1;
        line.clear();
        
    }
    
    std::cout << "Done." << std::endl;
}

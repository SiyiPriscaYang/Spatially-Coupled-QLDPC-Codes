#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <string>
#include "tools.h"
#include "MD.h"
#include "SC.h"
#include "sc_qldpc_optimization.h"

using namespace std;

int main(int argc, char* argv[]) {

    vector<vector<int> > Code1_3_8{{2,1,3,8,4,8,3,3,2,0,6,1,6,6,2,5,6,8,2,0,4,1,5,7},{2,2,6,5,6,3,1,0,7,6,2,0,0,4,3,8,6,0,0,7,5,8,5,3}};// 0,11,5113
    vector<vector<int> > Code2_3_8{{2, 3, 5, 3, 4, 0, 7, 0, 3, 6, 1, 6, 6, 0, 2, 8, 4, 8, 2, 7, 8, 5, 5, 1},{3, 3, 6, 6, 8, 3, 2, 5, 1, 4, 2, 0, 0, 4, 7, 8, 5, 1, 0, 7, 5, 8, 6, 2}};// 110,264,48142,
    vector<vector<int> > Code3_3_8{{8, 1, 2, 2, 6, 8, 6, 7, 0, 2, 4, 5, 3, 7, 1, 5, 8, 6, 8, 6, 2, 2, 1, 0},{8, 8, 2, 6, 6, 0, 3, 0, 7, 6, 6, 2, 1, 8, 2, 5, 3, 1, 8, 5, 2, 3, 8, 7}};// 0,11,5131, memort 2,2
    vector<vector<int> > Code4_3_8{{0, 1, 7, 2, 5, 8, 6, 8, 0, 2, 4, 8, 3, 6, 0, 5, 3, 6, 8, 6, 2, 2, 1, 0},{0, 5, 1, 6, 3, 3, 2, 0, 6, 6, 8, 0, 1, 8, 2, 5, 0, 4, 8, 8, 2, 6, 2, 7}};// 66,143,23120,
    vector<vector<int> > Code5_3_8{{1,3,0,2,0,0,1,2,0,0,3,1,1,2,2,2,2,0,0,0,2,3,0,1},{1,2,1,3,0,2,1,0,1,1,0,0,1,0,2,2,2,0,3,0,3,1,0,1}};
    vector<vector<int> > Code6_3_8{{2,4,3,1,2,3,2,3,3,3,0,5,4,5,0,1,2,0,2,3,0,0,3,5},{5,1,3,0,2,0,3,5,0,0,1,5,3,2,2,3,1,5,5,1,3,5,3,2}};
    vector<vector<int> > Code7_3_8{{14, 12, 3, 3, 9,  11, 0, 2, 15, 6, 8, 12, 1, 14, 0, 7, 4, 0, 15, 10, 15, 3, 11, 12},{14, 9,  5, 5, 12, 0,  7, 3, 12, 3, 6, 2,  1, 15, 0, 8, 0, 1, 15, 14, 7,  2, 4,  15}};//memory 3,3
    vector<vector<int> > Code1_3_7{{2, 5, 6, 8, 0, 6, 5, 6, 7, 2, 6, 5, 0, 6, 7, 0, 2, 1, 7, 8, 5},{2, 1, 3, 2, 0, 7, 6, 1, 8, 8, 6, 8, 0, 5, 6, 6, 2, 3, 0, 5, 1}};// 0,0,2316,
    vector<vector<int> > Code2_3_7{{8,3,1,7,7,2,1,2,4,8,1,6,0,4,7,5,5,4,3,6,0},{5,0,4,8,1,6,1,6,8,1,4,5,4,3,3,7,2,7,7,0,2}};//30,150,17804
    vector<vector<int> > Code3_3_7{{6, 2, 5, 1, 6, 0, 8, 4, 3, 8, 0, 1, 2, 2, 2, 6, 0, 8, 1, 4, 1},{6, 1, 0, 6, 7, 5, 2, 8, 3, 0, 2, 2, 4, 6, 1, 8, 3, 6, 8, 6, 0}};// 0,0,2466,
    vector<vector<int> > Code4_3_7{{6, 2, 3, 0, 6, 7, 0, 2, 6, 8, 5, 6, 8, 7, 2, 0, 0, 8, 1, 4, 1},{2, 5, 7, 6, 7, 2, 2, 8, 0, 0, 1, 6, 4, 6, 1, 8, 3, 6, 8, 0, 0}};// 40,140,18710,
    vector<vector<int> > Code5_3_7{{0, 3, 1, 3, 2, 1, 2, 1, 0, 2, 3, 1, 3, 0, 2, 3, 0, 0, 3, 2, 1},{1, 3, 3, 0, 0, 2, 1, 2, 1, 0, 2, 3, 3, 1, 3, 0, 2, 1, 3, 1, 2}};
    vector<vector<int> > Code6_3_7{{0,2,0,3,2,5,3,5,0,3,4,1,2,5,0,4,2,2,3,3,1},{2,4,5,2,0,5,0,3,2,2,4,2,3,5,4,3,0,0,5,2,4}};
    vector<vector<int> > Code7_3_7{{4,10,3,12,13,1,5,7,8,8,11,7,11,3,1,6,12,3,2,0,15},{4,11,13,12,7,3,12,3,13,15,10,6,1,0,4,1,5,7,14,15,3}};


    vector<vector<int> > P=Code2_3_7_2;
    string filename = "Code2_3_7.txt";

    int r1=3;
    int r2=3;
    int n1=7;
    int n2=7;
    int m1=2;
    int m2=2;
    int L1=10;
    int L2=10;
    int r=r1*n2+r2*n1;
    int n=n1*n2+r1*r2;

    vector<vector<int> > a= initialize(m1,m2);
    vector<vector<int> > b=a;

    MD md(a,b);
    vector<vector<int> > Ha(r1,vector<int> (n1,1));
    vector<vector<int> > Hb(r2,vector<int> (n2,1));
    comcodeqsc Hqsc(Ha,Hb);

    vector<vector<int> > H=md.SC(P,Hqsc.H[0].elist,Hqsc.H[1].elist,Hqsc.n1,Hqsc.n2,Hqsc.r1,Hqsc.r2,L1,L2);

    ofstream outfile(filename);
    for(int i=0;i<H.size();i++){
        for(int j=0;j<H[0].size();j++)
            outfile<<H[i][j]<<" ";
        outfile<<endl;
    }
    outfile.close();


    return 0;
}

//
// Created by Siyi Yang on 2/15/24.
//

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

//    double p = stod(argv[1]);
//    int num_trials = stoi(argv[2]);
//    int num_iter = stoi(argv[3]);
//    int job_index = stoi(argv[4]);
int main(int argc, char* argv[]) {

    int r1=3;
    int r2=3;
    int n1=7;
    int n2=7;
    int m1=2;
    int m2=2;
    int L1=10;
    int L2=10;

    vector<vector<int> > P_3_8_m2_unf{{2,1,3,8,4,8,3,3,2,0,6,1,6,6,2,5,6,8,2,0,4,1,5,7},{2,2,6,5,6,3,1,0,7,6,2,0,0,4,3,8,6,0,0,7,5,8,5,3}};// 0,11,5113
    vector<vector<int> > P_3_8_m2_gd{{8,1,2,2,6,8,6,7,0,2,4,5,3,7,1,5,8,6,8,6,2,2,1,0},{8,8,2,6,6,0,3,0,7,6,6,2,1,8,2,5,3,1,8,5,2,3,8,7}};// 0,11,5131, memort 2,2
    vector<vector<int> > P_3_8_m3{{14,12,3,3,9,11,0,2,15,6,8,12,1,14,0,7,4,0,15,10,15,3,11,12},{14,9,5,5,12,0,7,3,12,3,6,2,1,15,0,8,0,1,15,14,7,2,4,15}};//memory 3,3
    vector<vector<int> > P_3_7_m2_unf{{2,5,6,8,0,6,5,6,7,2,6,5,0,6,7,0,2,1,7,8,5},{2,1,3,2,0,7,6,1,8,8,6,8,0,5,6,6,2,3,0,5,1}};// 0,0,2316,
    vector<vector<int> > P_3_7_m2_gd{{6,2,5,1,6,0,8,4,3,8,0,1,2,2,2,6,0,8,1,4,1},{6,1,0,6,7,5,2,8,3,0,2,2,4,6,1,8,3,6,8,6,0}};// 0,0,2466,
    vector<vector<int> > P_3_7_m3{{4,10,3,12,13,1,5,7,8,8,11,7,11,3,1,6,12,3,2,0,15},{4,11,13,12,7,3,12,3,13,15,10,6,1,0,4,1,5,7,14,15,3}};

    vector<vector<int> > P=P_3_8_m2_unf;
    string filename = "H_3_8_m2_unf.txt";
    int r=r1*n2+r2*n1;
    int n=n1*n2+r1*r2;

    vector<vector<int> > a= initialize(m1,m2);
    vector<vector<int> > b=a;

    MD md(a,b);
    vector<vector<int> > Ha(r1);
    vector<vector<int> > Hb(r2);
    for(int i=0;i<r1;i++){
        vector<int> t(n1);
        for(int j=0;j<n1;j++)
            t[j]=1;
        Ha[i]=t;
    }
    for(int i=0;i<r2;i++){
        vector<int> t(n2);
        for(int j=0;j<n1;j++)
            t[j]=1;
        Hb[i]=t;
    }
    comcodeqsc Hqsc(Ha,Hb);

    vector<vector<int> > H=md.SC(P,Hqsc.H[0].elist,Hqsc.H[1].elist,Hqsc.n1,Hqsc.n2,Hqsc.r1,Hqsc.r2,L1,L2);
    /* cycle-(2g+4) list for A(B): H[0](H[1]).cyclelist[g]:
     * each row: a cycle
     * each entry: edge index of the cycle, label in H[0](H[1]).pcme */

    ofstream outfile(filename);
    for(int i=0;i<H.size();i++){
        for(int j=0;j<H[0].size();j++)
            outfile<<H[i][j]<<" ";
        outfile<<endl;
    }
    outfile.close();


    return 0;
}

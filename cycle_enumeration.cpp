//
// Created by Siyi Yang on 1/3/24.
//


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"
#include "MD.h"
#include "sc_qldpc_optimization.h"

using namespace std;

int main() {
    int r1 = 3;
    int r2 = 3;
    int n1 = 7;
    int n2 = 7;
    int m1 = 1;
    int m2 = 1;
    int L1 = 10;
    int L2 = 10;
    vector<vector<int> > P_uni_opt{{2, 1, 3, 8, 4, 8, 3, 3, 2, 0, 6, 1, 6, 6, 2, 5, 6, 8, 2, 0, 4, 1, 5, 7},
                                   {2, 2, 6, 5, 6, 3, 1, 0, 7, 6, 2, 0, 0, 4, 3, 8, 6, 0, 0, 7, 5, 8, 5, 3}};// 0,11,5113
    vector<vector<int> > P_uni_ini{{2, 3, 5, 3, 4, 0, 7, 0, 3, 6, 1, 6, 6, 0, 2, 8, 4, 8, 2, 7, 8, 5, 5, 1},
                                   {3, 3, 6, 6, 8, 3, 2, 5, 1, 4, 2, 0, 0, 4, 7, 8, 5, 1, 0, 7, 5, 8, 6, 2}};// 110,264,48142,
    vector<vector<int> > P_gd_opt{{8, 1, 2, 2, 6, 8, 6, 7, 0, 2, 4, 5, 3, 7, 1, 5, 8, 6, 8, 6, 2, 2, 1, 0},
                                  {8, 8, 2, 6, 6, 0, 3, 0, 7, 6, 6, 2, 1, 8, 2, 5, 3, 1, 8, 5, 2, 3, 8, 7}};// 0,11,5131, memort 2,2
    vector<vector<int> > P_gd_ini{{0, 1, 7, 2, 5, 8, 6, 8, 0, 2, 4, 8, 3, 6, 0, 5, 3, 6, 8, 6, 2, 2, 1, 0},
                                  {0, 5, 1, 6, 3, 3, 2, 0, 6, 6, 8, 0, 1, 8, 2, 5, 0, 4, 8, 8, 2, 6, 2, 7}};// 66,143,23120,
    vector<vector<int> > P3_gd_opt{{14, 12, 3, 3, 9,  11, 0, 2, 15, 6, 8, 12, 1, 14, 0, 7, 4, 0, 15, 10, 15, 3, 11, 12},
                                   {14, 9,  5, 5, 12, 0,  7, 3, 12, 3, 6, 2,  1, 15, 0, 8, 0, 1, 15, 14, 7,  2, 4,  15}};//memory 3,3
    vector<vector<int> > P2_uni_opt{{2, 5, 6, 8, 0, 6, 5, 6, 7, 2, 6, 5, 0, 6, 7, 0, 2, 1, 7, 8, 5},
                                    {2, 1, 3, 2, 0, 7, 6, 1, 8, 8, 6, 8, 0, 5, 6, 6, 2, 3, 0, 5, 1}};// 0,0,2316,
    vector<vector<int> > P2_uni_ini{{2, 0, 7, 1, 3, 4, 1, 6, 4, 2, 5, 1, 0, 6, 7, 3, 5, 4, 7, 8, 8},
                                    {6, 2, 8, 1, 0, 7, 4, 0, 8, 6, 7, 4, 4, 5, 3, 7, 2, 3, 1, 5, 1}};// 70,70,17708,
    vector<vector<int> > P2_gd_opt{{6, 2, 5, 1, 6, 0, 8, 4, 3, 8, 0, 1, 2, 2, 2, 6, 0, 8, 1, 4, 1},
                                   {6, 1, 0, 6, 7, 5, 2, 8, 3, 0, 2, 2, 4, 6, 1, 8, 3, 6, 8, 6, 0}};// 0,0,2466,
    vector<vector<int> > P2_gd_ini{{6, 2, 3, 0, 6, 7, 0, 2, 6, 8, 5, 6, 8, 7, 2, 0, 0, 8, 1, 4, 1},
                                   {2, 5, 7, 6, 7, 2, 2, 8, 0, 0, 1, 6, 4, 6, 1, 8, 3, 6, 8, 0, 0}};// 40,140,18710,
    vector<vector<int> > Pm1_3_7{{0, 3, 1, 3, 2, 1, 2, 1, 0, 2, 3, 1, 3, 0, 2, 3, 0, 0, 3, 2, 1},
                                 {1, 3, 3, 0, 0, 2, 1, 2, 1, 0, 2, 3, 3, 1, 3, 0, 2, 1, 3, 1, 2}};
    //vector<vector<int>> P_gd_opt({{5,8,0,3,8,5,0,0,1,3,5,4,6,6,2,4,2,7,2,7,1},{4,2,4,7,3,2,1,6,3,3,8,0,5,1,0,5,2,0,5,3,1}});// (3,3) 0,620,109722

    vector<vector<int> > P = Pm1_3_7;

    double w = 25;
    int r = r1 * n2 + r2 * n1;
    int n = n1 * n2 + r1 * r2;

    vector<vector<double> > we = get_weight(r1, r2, n1, n2);
    vector<vector<int> > a = initialize(m1, m2);
    vector<vector<int> > b = a;

    MD md(a, b);
    vector<vector<int> > Ha(r1);
    vector<vector<int> > Hb(r2);
    for (int i = 0; i < r1; i++) {
        vector<int> t(n1);
        for (int j = 0; j < n1; j++)
            t[j] = 1;
        Ha[i] = t;
    }
    for (int i = 0; i < r2; i++) {
        vector<int> t(n2);
        for (int j = 0; j < n1; j++)
            t[j] = 1;
        Hb[i] = t;
    }
    comcodeqsc Hqsc(Ha, Hb);

    vector<vector<int> > H = md.SC(P, Hqsc.H[0].elist, Hqsc.H[1].elist, Hqsc.n1, Hqsc.n2, Hqsc.r1, Hqsc.r2, L1, L2);
    //qCode qsc(H);

    return 0;
}
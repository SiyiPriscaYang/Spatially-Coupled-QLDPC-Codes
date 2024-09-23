//
// Created by Siyi Yang on 1/2/24.
//

#ifndef QLDPC_SC_OPTIMIZATION_H
#define QLDPC_SC_OPTIMIZATION_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"
#include "common_code.h"
#include "SC.h"

using namespace std;


class comcodesc: public comcode{
public:
    vector<vector<vector<int> > > cyclelist;// cyc_4,cyc_6,cyc_8. e.g., cyc_4: N_cyc_4*4
    vector<vector<vector<int> > > edge_cyclelist;// E_num*6*~ [e][0]: cycle-4 index [e][1]: sign while traversing the cycle-4
    comcodesc(){vector<vector<int> >  h(0);vector<vector<int> >  he(0);vector<vector<vector<int> > > vneigh(0);vector<vector<vector<int> > > cneigh(0);vector<vector<int> >  el(0);vector<vector<vector<int> > > cyclel(0);vector<vector<vector<int> > > ecyclel(0);pcm=h;pcme=he;vneighbor=vneigh;cneighbor=cneigh;elist=el;cyclelist=cyclel;edge_cyclelist=ecyclel;}
    void cycle();
    void cycle_RTC();
    void printcycle(){
        for(int i=0;i<cyclelist.size();i++){
            //print_matrix(cyclelist[i]);
            cout<<cyclelist[i].size()<<endl;
        }
    }
    vector<int> AO(vector<int> a,vector<int> P,int th1,int th2,vector<int> w);
    vector<int> AO_RTC(SC pattern,vector<int> P,int th1,int th2,int L);
};


void comcodesc::cycle(){
    int i1,i2,i3,i4,j1,j2,j3,j4,x1,x2,x3,x4,x5,x6,x7,e1,e2,e3,e4,e5,e6,e7,e8;
    int gamma=cneighbor.size();
    int kappa=vneighbor.size();
    vector<int> t1,t2,t3,t4,t5,t6,t7;
    vector<vector<vector<int> > > list(3);//3*Cyc(2g+4)_num*2g
    vector<vector<vector<int> > > ecyclist(elist.size());// E_num*6*~ [e][0]: cycle-4 index [e][1]: sign while traversing the cycle-4
    int n4=0;
    int n6=0;
    int n8=0;
    for(i1=0;i1<elist.size();i1++){
        vector<vector<int> >  t(6);
        ecyclist[i1]=t;
    }
    for(i1=0;i1<gamma;i1++){
        t1=cneighbor[i1][1];
        for(x1=0;x1<t1.size();x1++){
            j1=t1[x1];
            e1=cneighbor[i1][0][x1];//H[i1][j1]
            t2=vneighbor[j1][1];
            for(x2=0;x2<t2.size();x2++){
                i2=t2[x2];
                if(i2!=i1) {
                    e2 = vneighbor[j1][0][x2];//H[i2][j1]
                    t3 = cneighbor[i2][1];
                    for (x3 = 0; x3 < t3.size(); x3++) {
                        j2 = t3[x3];
                        if (j2 != j1) {
                            e3 = cneighbor[i2][0][x3];//H[i2][j2]
                            if (pcm[i1][j2])  {
                                e4 = pcme[i1][j2];//vneighbor[i2][0][x4];//H[i1][j2]
                                if (i2 > i1) {
                                    if (j2 > j1) {
                                        vector<int> t {e1, e2, e3, e4};//{i1,j1,i2,j2});
                                        list[0].push_back(t);
                                        ecyclist[e1][0].push_back(n4);
                                        ecyclist[e2][0].push_back(n4);
                                        ecyclist[e3][0].push_back(n4);
                                        ecyclist[e4][0].push_back(n4);
                                        ecyclist[e1][1].push_back(1);
                                        ecyclist[e2][1].push_back(-1);
                                        ecyclist[e3][1].push_back(1);
                                        ecyclist[e4][1].push_back(-1);
                                        n4++;

                                        for (x4 = 0; x4 < t1.size(); x4++) {
                                            j3 = t1[x4];
                                            if ((j3 != j1) && (j3 != j2) && pcm[i1][j3] && pcm[i2][j3]) {
                                                e5 = pcme[i1][j3];
                                                e6 = pcme[i2][j3];
                                                vector<int> tt {e5, e1, e2, e6, e5, e4, e3, e6};
                                                list[2].push_back(tt);
                                                ecyclist[e1][4].push_back(n8);
                                                ecyclist[e2][4].push_back(n8);
                                                ecyclist[e3][4].push_back(n8);
                                                ecyclist[e4][4].push_back(n8);
                                                ecyclist[e5][4].push_back(n8);
                                                ecyclist[e6][4].push_back(n8);
                                                ecyclist[e1][5].push_back(-1);
                                                ecyclist[e2][5].push_back(1);
                                                ecyclist[e3][5].push_back(1);
                                                ecyclist[e4][5].push_back(-1);
                                                ecyclist[e5][5].push_back(2);
                                                ecyclist[e6][5].push_back(-2);
                                                n8++;

                                                if (j3 > j1) {
                                                    for (x5 = 0; x5 < t1.size(); x5++) {
                                                        j4 = t1[x5];
                                                        if ((j4 > j1) && (j4 != j2) && (j4 != j3) && pcm[i1][j4] &&
                                                            pcm[i2][j4]) {
                                                            e7 = pcme[i1][j4];
                                                            e8 = pcme[i2][j4];
                                                            vector<int> ttt {e1, e2, e8, e7, e5, e6, e3, e4};
                                                            list[2].push_back(tt);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(1);
                                                            ecyclist[e2][5].push_back(-1);
                                                            ecyclist[e3][5].push_back(1);
                                                            ecyclist[e4][5].push_back(-1);
                                                            ecyclist[e5][5].push_back(1);
                                                            ecyclist[e6][5].push_back(-1);
                                                            ecyclist[e7][5].push_back(-1);
                                                            ecyclist[e8][5].push_back(1);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        for (x4 = 0; x4 < t2.size(); x4++) {
                                            i3 = t2[x4];
                                            if ((i3 != i1) && (i3 != i2) && pcm[i3][j1] && pcm[i3][j2]) {
                                                e5 = pcme[i3][j1];
                                                e6 = pcme[i3][j2];
                                                vector<int> tt {e5, e1, e4, e6, e5, e2, e3, e6};
                                                list[2].push_back(tt);
                                                ecyclist[e1][4].push_back(n8);
                                                ecyclist[e2][4].push_back(n8);
                                                ecyclist[e3][4].push_back(n8);
                                                ecyclist[e4][4].push_back(n8);
                                                ecyclist[e5][4].push_back(n8);
                                                ecyclist[e6][4].push_back(n8);
                                                ecyclist[e1][5].push_back(-1);
                                                ecyclist[e2][5].push_back(-1);
                                                ecyclist[e3][5].push_back(1);
                                                ecyclist[e4][5].push_back(1);
                                                ecyclist[e5][5].push_back(2);
                                                ecyclist[e6][5].push_back(-2);
                                                n8++;

                                                if (i3 > i1) {
                                                    for (x5 = 0; x5 < t2.size(); x5++) {
                                                        i4 = t2[x5];
                                                        if ((i4 > i1) && (i4 != i2) && (i4 != i3) && pcm[i4][j1] &&
                                                            pcm[i4][j2]) {
                                                            e7 = pcme[i4][j1];
                                                            e8 = pcme[i4][j2];
                                                            vector<int> ttt {e1, e2, e8, e7, e5, e6, e3, e4};
                                                            list[2].push_back(tt);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(1);
                                                            ecyclist[e2][5].push_back(-1);
                                                            ecyclist[e3][5].push_back(1);
                                                            ecyclist[e4][5].push_back(-1);
                                                            ecyclist[e5][5].push_back(1);
                                                            ecyclist[e6][5].push_back(-1);
                                                            ecyclist[e7][5].push_back(-1);
                                                            ecyclist[e8][5].push_back(1);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    for (x4 = 0; x4 < t2.size(); x4++) {
                                        i3 = t2[x4];
                                        if ((i3 != i1) && (i3 != i2)) {
                                            e5 = pcme[i3][j1];
                                            for (x5 = 0; x5 < t2.size(); x5++) {
                                                i4 = t2[x5];
                                                if ((i4 != i3) && (i4 != i2) && (i4 != i1)) {
                                                    e8 = pcme[i4][j1];
                                                    t4 = cneighbor[i3][1];
                                                    for (x6 = 0; x6 < t4.size(); x6++) {
                                                        j3 = t4[x6];
                                                        if ((j3 > j2) && (j3 != j1) && pcm[i4][j3]) {
                                                            e6 = pcme[i3][j3];
                                                            e7 = pcme[i4][j3];
                                                            vector<int> t {e2, e3, e4, e1, e5, e6, e7, e8};
                                                            list[2].push_back(t);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(-1);
                                                            ecyclist[e2][5].push_back(1);
                                                            ecyclist[e3][5].push_back(-1);
                                                            ecyclist[e4][5].push_back(1);
                                                            ecyclist[e5][5].push_back(1);
                                                            ecyclist[e6][5].push_back(-1);
                                                            ecyclist[e7][5].push_back(1);
                                                            ecyclist[e8][5].push_back(-1);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                if(j2>j1){
                                    for(x4=0;x4<t1.size();x4++){
                                        j3=t1[x4];
                                        if((j3!=j1)&&(j3!=j2)){
                                            e5=pcme[i1][j3];
                                            for(x5=0;x5<t1.size();x5++){
                                                j4=t1[x5];
                                                if((j4!=j1)&&(j4!=j2)&&(j4!=j3)){
                                                    e8=pcme[i1][j4];
                                                    t4=vneighbor[j3][1];
                                                    for(x6=0;x6<t4.size();x6++){
                                                        i3=t4[x6];
                                                        if((i3>i2)&&(i3!=i1)&&pcm[i3][j4]){
                                                            e6=pcme[i3][j3];
                                                            e7=pcme[i3][j4];
                                                            vector<int> t {e1,e5,e6,e7,e8,e4,e3,e2};
                                                            list[2].push_back(t);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(1);
                                                            ecyclist[e2][5].push_back(-1);
                                                            ecyclist[e3][5].push_back(1);
                                                            ecyclist[e4][5].push_back(-1);
                                                            ecyclist[e5][5].push_back(-1);
                                                            ecyclist[e6][5].push_back(1);
                                                            ecyclist[e7][5].push_back(-1);
                                                            ecyclist[e8][5].push_back(1);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            if(i2>i1){
                                t4 = vneighbor[j2][1];
                                for (x4 = 0; x4 < t4.size(); x4++) {
                                    i3 = t4[x4];
                                    if ((i3 != i2) && (i3 != i1)) {
                                        e4 = vneighbor[j2][0][x4];//H[i3][j2]
                                        t5 = cneighbor[i3][1];
                                        for (x5 = 0; x5 < t5.size(); x5++) {
                                            j3 = t5[x5];
                                            if ((j3 != j1) && (j3 != j2)) {
                                                e5 = cneighbor[i3][0][x5];//H[i3][j3]
                                                if (pcm[i1][j3]) {
                                                    if (i3 > i1) {
                                                        if (j3 < j1) {
                                                            e6 = pcme[i1][j3];
                                                            vector<int> tt(
                                                                    {e1, e2, e3, e4, e5, e6});//{i1,j1,i2,j2,i3,j3});
                                                            list[1].push_back(tt);
                                                            ecyclist[e1][2].push_back(n6);
                                                            ecyclist[e2][2].push_back(n6);
                                                            ecyclist[e3][2].push_back(n6);
                                                            ecyclist[e4][2].push_back(n6);
                                                            ecyclist[e5][2].push_back(n6);
                                                            ecyclist[e6][2].push_back(n6);
                                                            ecyclist[e1][3].push_back(1);
                                                            ecyclist[e2][3].push_back(-1);
                                                            ecyclist[e3][3].push_back(1);
                                                            ecyclist[e4][3].push_back(-1);
                                                            ecyclist[e5][3].push_back(1);
                                                            ecyclist[e6][3].push_back(-1);
                                                            n6++;
                                                        }

                                                        t6 = vneighbor[j3][1];
                                                        for (x6 = 0; x6 < t6.size(); x6++) {
                                                            i4 = t6[x6];
                                                            e6 = pcme[i4][j3];
                                                            if ((i4 > i1) && (i4 != i2) && (i4 != i3)) {
                                                                t7 = cneighbor[i4][1];
                                                                for (x7 = 0; x7 < t7.size(); x7++) {
                                                                    j4 = t7[x7];
                                                                    if ((j4 > j1) && (j4 != j2) && (j4 != j3)) {
                                                                        e7 = pcme[i4][j4];
                                                                        e8 = pcme[i1][j4];
                                                                        vector<int> tt {e1, e2, e3, e4, e5, e6, e7,
                                                                                        e8};//{i1,j1,i2,j2,i3,j3});
                                                                        list[2].push_back(tt);
                                                                        ecyclist[e1][4].push_back(n8);
                                                                        ecyclist[e2][4].push_back(n8);
                                                                        ecyclist[e3][4].push_back(n8);
                                                                        ecyclist[e4][4].push_back(n8);
                                                                        ecyclist[e5][4].push_back(n8);
                                                                        ecyclist[e6][4].push_back(n8);
                                                                        ecyclist[e7][4].push_back(n8);
                                                                        ecyclist[e1][5].push_back(1);
                                                                        ecyclist[e2][5].push_back(-1);
                                                                        ecyclist[e3][5].push_back(1);
                                                                        ecyclist[e4][5].push_back(-1);
                                                                        ecyclist[e5][5].push_back(1);
                                                                        ecyclist[e6][5].push_back(-1);
                                                                        ecyclist[e7][5].push_back(1);
                                                                        ecyclist[e8][5].push_back(-1);
                                                                        n8++;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                    if (pcm[i3][j1]) {
                                                        e6 = pcme[i1][j3];
                                                        e7 = pcme[i3][j1];
                                                        vector<int> tt {e1, e7, e4, e3, e2, e7, e5,
                                                                        e6};//{i1,j1,i2,j2,i3,j3});
                                                        list[2].push_back(tt);
                                                        ecyclist[e1][4].push_back(n8);
                                                        ecyclist[e2][4].push_back(n8);
                                                        ecyclist[e3][4].push_back(n8);
                                                        ecyclist[e4][4].push_back(n8);
                                                        ecyclist[e5][4].push_back(n8);
                                                        ecyclist[e6][4].push_back(n8);
                                                        ecyclist[e7][4].push_back(n8);
                                                        ecyclist[e1][5].push_back(1);
                                                        ecyclist[e2][5].push_back(1);
                                                        ecyclist[e3][5].push_back(-1);
                                                        ecyclist[e4][5].push_back(1);
                                                        ecyclist[e5][5].push_back(1);
                                                        ecyclist[e6][5].push_back(-1);
                                                        ecyclist[e7][5].push_back(-2);
                                                        n8++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cyclelist=list;
    edge_cyclelist=ecyclist;
}


void comcodesc::cycle_RTC(){
    int i1,i2,i3,i4,j1,j2,j3,j4,x1,x2,x3,x4,x5,x6,x7,e1,e2,e3,e4,e5,e6,e7,e8;
    int gamma=cneighbor.size();
    int kappa=vneighbor.size();
    vector<int> t1,t2,t3,t4,t5,t6,t7;
    vector<vector<vector<int> > > list(3);//3*Cyc(2g+4)_num*2g
    vector<vector<vector<int> > > ecyclist(elist.size());// E_num*6*~ [e][0]: cycle-4 index [e][1]: sign while traversing the cycle-4
    int n4=0;
    int n6=0;
    int n8=0;
    for(i1=0;i1<elist.size();i1++){
        vector<vector<int> >  t(6);
        ecyclist[i1]=t;
    }
    for(i1=0;i1<gamma;i1++){
        t1=cneighbor[i1][1];
        for(x1=0;x1<t1.size();x1++){
            j1=t1[x1];
            e1=cneighbor[i1][0][x1];//H[i1][j1]
            t2=vneighbor[j1][1];
            for(x2=0;x2<t2.size();x2++){
                i2=t2[x2];
                if(i2!=i1) {
                    e2 = vneighbor[j1][0][x2];//H[i2][j1]
                    t3 = cneighbor[i2][1];
                    for (x3 = 0; x3 < t3.size(); x3++) {
                        j2 = t3[x3];
                        if (j2 != j1) {
                            e3 = cneighbor[i2][0][x3];//H[i2][j2]
                            if (pcm[i1][j2])  {
                                e4 = pcme[i1][j2];//vneighbor[i2][0][x4];//H[i1][j2]
                                if (i2 > i1) {
                                    if (j2 > j1) {
                                        vector<int> t {e1, e2, e3, e4};//{i1,j1,i2,j2});
                                        list[0].push_back(t);
                                        ecyclist[e1][0].push_back(n4);
                                        ecyclist[e2][0].push_back(n4);
                                        ecyclist[e3][0].push_back(n4);
                                        ecyclist[e4][0].push_back(n4);
                                        ecyclist[e1][1].push_back(0);
                                        ecyclist[e2][1].push_back(1);
                                        ecyclist[e3][1].push_back(2);
                                        ecyclist[e4][1].push_back(3);
                                        n4++;

                                        for (x4 = 0; x4 < t1.size(); x4++) {
                                            j3 = t1[x4];
                                            if ((j3 != j1) && (j3 != j2) && pcm[i1][j3] && pcm[i2][j3]) {
                                                e5 = pcme[i1][j3];
                                                e6 = pcme[i2][j3];
                                                vector<int> tt {e5, e1, e2, e6, e5, e4, e3, e6};
                                                list[2].push_back(tt);
                                                ecyclist[e1][4].push_back(n8);
                                                ecyclist[e2][4].push_back(n8);
                                                ecyclist[e3][4].push_back(n8);
                                                ecyclist[e4][4].push_back(n8);
                                                ecyclist[e5][4].push_back(n8);
                                                ecyclist[e5][4].push_back(n8);
                                                ecyclist[e6][4].push_back(n8);
                                                ecyclist[e6][4].push_back(n8);
                                                ecyclist[e1][5].push_back(1);
                                                ecyclist[e2][5].push_back(2);
                                                ecyclist[e3][5].push_back(6);
                                                ecyclist[e4][5].push_back(5);
                                                ecyclist[e5][5].push_back(0);
                                                ecyclist[e5][5].push_back(4);
                                                ecyclist[e6][5].push_back(3);
                                                ecyclist[e6][5].push_back(7);
                                                n8++;

                                                if (j3 > j1) {
                                                    for (x5 = 0; x5 < t1.size(); x5++) {
                                                        j4 = t1[x5];
                                                        if ((j4 > j1) && (j4 != j2) && (j4 != j3) && pcm[i1][j4] &&
                                                            pcm[i2][j4]) {
                                                            e7 = pcme[i1][j4];
                                                            e8 = pcme[i2][j4];
                                                            vector<int> ttt {e1, e2, e8, e7, e5, e6, e3, e4};
                                                            list[2].push_back(tt);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(0);
                                                            ecyclist[e2][5].push_back(1);
                                                            ecyclist[e3][5].push_back(6);
                                                            ecyclist[e4][5].push_back(7);
                                                            ecyclist[e5][5].push_back(4);
                                                            ecyclist[e6][5].push_back(5);
                                                            ecyclist[e7][5].push_back(3);
                                                            ecyclist[e8][5].push_back(2);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        for (x4 = 0; x4 < t2.size(); x4++) {
                                            i3 = t2[x4];
                                            if ((i3 != i1) && (i3 != i2) && pcm[i3][j1] && pcm[i3][j2]) {
                                                e5 = pcme[i3][j1];
                                                e6 = pcme[i3][j2];
                                                vector<int> tt {e5, e1, e4, e6, e5, e2, e3, e6};
                                                list[2].push_back(tt);
                                                ecyclist[e1][4].push_back(n8);
                                                ecyclist[e2][4].push_back(n8);
                                                ecyclist[e3][4].push_back(n8);
                                                ecyclist[e4][4].push_back(n8);
                                                ecyclist[e5][4].push_back(n8);
                                                ecyclist[e6][4].push_back(n8);
                                                ecyclist[e5][4].push_back(n8);
                                                ecyclist[e6][4].push_back(n8);
                                                ecyclist[e1][5].push_back(1);
                                                ecyclist[e2][5].push_back(5);
                                                ecyclist[e3][5].push_back(6);
                                                ecyclist[e4][5].push_back(2);
                                                ecyclist[e5][5].push_back(0);
                                                ecyclist[e6][5].push_back(3);
                                                ecyclist[e5][5].push_back(4);
                                                ecyclist[e6][5].push_back(7);
                                                n8++;

                                                if (i3 > i1) {
                                                    for (x5 = 0; x5 < t2.size(); x5++) {
                                                        i4 = t2[x5];
                                                        if ((i4 > i1) && (i4 != i2) && (i4 != i3) && pcm[i4][j1] &&
                                                            pcm[i4][j2]) {
                                                            e7 = pcme[i4][j1];
                                                            e8 = pcme[i4][j2];
                                                            vector<int> ttt {e1, e2, e8, e7, e5, e6, e3, e4};
                                                            list[2].push_back(tt);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(0);
                                                            ecyclist[e2][5].push_back(1);
                                                            ecyclist[e3][5].push_back(6);
                                                            ecyclist[e4][5].push_back(7);
                                                            ecyclist[e5][5].push_back(4);
                                                            ecyclist[e6][5].push_back(5);
                                                            ecyclist[e7][5].push_back(3);
                                                            ecyclist[e8][5].push_back(2);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    for (x4 = 0; x4 < t2.size(); x4++) {
                                        i3 = t2[x4];
                                        if ((i3 != i1) && (i3 != i2)) {
                                            e5 = pcme[i3][j1];
                                            for (x5 = 0; x5 < t2.size(); x5++) {
                                                i4 = t2[x5];
                                                if ((i4 != i3) && (i4 != i2) && (i4 != i1)) {
                                                    e8 = pcme[i4][j1];
                                                    t4 = cneighbor[i3][1];
                                                    for (x6 = 0; x6 < t4.size(); x6++) {
                                                        j3 = t4[x6];
                                                        if ((j3 > j2) && (j3 != j1) && pcm[i4][j3]) {
                                                            e6 = pcme[i3][j3];
                                                            e7 = pcme[i4][j3];
                                                            vector<int> t {e2, e3, e4, e1, e5, e6, e7, e8};
                                                            list[2].push_back(t);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(3);
                                                            ecyclist[e2][5].push_back(0);
                                                            ecyclist[e3][5].push_back(1);
                                                            ecyclist[e4][5].push_back(2);
                                                            ecyclist[e5][5].push_back(4);
                                                            ecyclist[e6][5].push_back(5);
                                                            ecyclist[e7][5].push_back(6);
                                                            ecyclist[e8][5].push_back(7);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                if(j2>j1){
                                    for(x4=0;x4<t1.size();x4++){
                                        j3=t1[x4];
                                        if((j3!=j1)&&(j3!=j2)){
                                            e5=pcme[i1][j3];
                                            for(x5=0;x5<t1.size();x5++){
                                                j4=t1[x5];
                                                if((j4!=j1)&&(j4!=j2)&&(j4!=j3)){
                                                    e8=pcme[i1][j4];
                                                    t4=vneighbor[j3][1];
                                                    for(x6=0;x6<t4.size();x6++){
                                                        i3=t4[x6];
                                                        if((i3>i2)&&(i3!=i1)&&pcm[i3][j4]){
                                                            e6=pcme[i3][j3];
                                                            e7=pcme[i3][j4];
                                                            vector<int> t {e1,e5,e6,e7,e8,e4,e3,e2};
                                                            list[2].push_back(t);
                                                            ecyclist[e1][4].push_back(n8);
                                                            ecyclist[e2][4].push_back(n8);
                                                            ecyclist[e3][4].push_back(n8);
                                                            ecyclist[e4][4].push_back(n8);
                                                            ecyclist[e5][4].push_back(n8);
                                                            ecyclist[e6][4].push_back(n8);
                                                            ecyclist[e7][4].push_back(n8);
                                                            ecyclist[e8][4].push_back(n8);
                                                            ecyclist[e1][5].push_back(0);
                                                            ecyclist[e2][5].push_back(7);
                                                            ecyclist[e3][5].push_back(6);
                                                            ecyclist[e4][5].push_back(5);
                                                            ecyclist[e5][5].push_back(1);
                                                            ecyclist[e6][5].push_back(2);
                                                            ecyclist[e7][5].push_back(3);
                                                            ecyclist[e8][5].push_back(4);
                                                            n8++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            if(i2>i1){
                                t4 = vneighbor[j2][1];
                                for (x4 = 0; x4 < t4.size(); x4++) {
                                    i3 = t4[x4];
                                    if ((i3 != i2) && (i3 != i1)) {
                                        e4 = vneighbor[j2][0][x4];//H[i3][j2]
                                        t5 = cneighbor[i3][1];
                                        for (x5 = 0; x5 < t5.size(); x5++) {
                                            j3 = t5[x5];
                                            if ((j3 != j1) && (j3 != j2)) {
                                                e5 = cneighbor[i3][0][x5];//H[i3][j3]
                                                if (pcm[i1][j3]) {
                                                    if (i3 > i1) {
                                                        if (j3 < j1) {
                                                            e6 = pcme[i1][j3];
                                                            vector<int> tt(
                                                                    {e1, e2, e3, e4, e5, e6});//{i1,j1,i2,j2,i3,j3});
                                                            list[1].push_back(tt);
                                                            ecyclist[e1][2].push_back(n6);
                                                            ecyclist[e2][2].push_back(n6);
                                                            ecyclist[e3][2].push_back(n6);
                                                            ecyclist[e4][2].push_back(n6);
                                                            ecyclist[e5][2].push_back(n6);
                                                            ecyclist[e6][2].push_back(n6);
                                                            ecyclist[e1][3].push_back(0);
                                                            ecyclist[e2][3].push_back(1);
                                                            ecyclist[e3][3].push_back(2);
                                                            ecyclist[e4][3].push_back(3);
                                                            ecyclist[e5][3].push_back(4);
                                                            ecyclist[e6][3].push_back(5);
                                                            n6++;
                                                        }

                                                        t6 = vneighbor[j3][1];
                                                        for (x6 = 0; x6 < t6.size(); x6++) {
                                                            i4 = t6[x6];
                                                            e6 = pcme[i4][j3];
                                                            if ((i4 > i1) && (i4 != i2) && (i4 != i3)) {
                                                                t7 = cneighbor[i4][1];
                                                                for (x7 = 0; x7 < t7.size(); x7++) {
                                                                    j4 = t7[x7];
                                                                    if ((j4 > j1) && (j4 != j2) && (j4 != j3)) {
                                                                        e7 = pcme[i4][j4];
                                                                        e8 = pcme[i1][j4];
                                                                        vector<int> tt {e1, e2, e3, e4, e5, e6, e7,
                                                                                        e8};//{i1,j1,i2,j2,i3,j3});
                                                                        list[2].push_back(tt);
                                                                        ecyclist[e1][4].push_back(n8);
                                                                        ecyclist[e2][4].push_back(n8);
                                                                        ecyclist[e3][4].push_back(n8);
                                                                        ecyclist[e4][4].push_back(n8);
                                                                        ecyclist[e5][4].push_back(n8);
                                                                        ecyclist[e6][4].push_back(n8);
                                                                        ecyclist[e7][4].push_back(n8);
                                                                        ecyclist[e1][5].push_back(0);
                                                                        ecyclist[e2][5].push_back(1);
                                                                        ecyclist[e3][5].push_back(2);
                                                                        ecyclist[e4][5].push_back(3);
                                                                        ecyclist[e5][5].push_back(4);
                                                                        ecyclist[e6][5].push_back(5);
                                                                        ecyclist[e7][5].push_back(6);
                                                                        ecyclist[e8][5].push_back(7);
                                                                        n8++;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                    if (pcm[i3][j1]) {
                                                        e6 = pcme[i1][j3];
                                                        e7 = pcme[i3][j1];
                                                        vector<int> tt {e1, e7, e4, e3, e2, e7, e5,
                                                                        e6};//{i1,j1,i2,j2,i3,j3});
                                                        list[2].push_back(tt);
                                                        ecyclist[e1][4].push_back(n8);
                                                        ecyclist[e2][4].push_back(n8);
                                                        ecyclist[e3][4].push_back(n8);
                                                        ecyclist[e4][4].push_back(n8);
                                                        ecyclist[e5][4].push_back(n8);
                                                        ecyclist[e6][4].push_back(n8);
                                                        ecyclist[e7][4].push_back(n8);
                                                        ecyclist[e7][4].push_back(n8);
                                                        ecyclist[e1][5].push_back(0);
                                                        ecyclist[e2][5].push_back(4);
                                                        ecyclist[e3][5].push_back(3);
                                                        ecyclist[e4][5].push_back(2);
                                                        ecyclist[e5][5].push_back(6);
                                                        ecyclist[e6][5].push_back(7);
                                                        ecyclist[e7][5].push_back(1);
                                                        ecyclist[e7][5].push_back(5);
                                                        n8++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cyclelist=list;
    edge_cyclelist=ecyclist;
}


vector<int> comcodesc::AO(vector<int> a,vector<int> P,int th1,int th2,vector<int> w){
    vector<int> Pt(P);
    int i,j,k,l,e,s_old,s_new,v_diff,v_temp,v_opt,v_cur,diff_temp,diff_opt;
    bool start=true;
    bool nonstop=true;
    int g=cyclelist.size();
    vector<int> ncyc(g);
    vector<int> d_temp(g);
    vector<int> d_opt(g);
    vector<vector<int> >  sum(g);
    for(i=0;i<g;i++){
        vector<int> t(cyclelist[i].size());
        sum[i]=t;
    }
    int m= vec_max(a);
    vector<int> inda(m+1);
    for(k=0;k<a.size();k++)
        inda[a[k]]=k;
    vector<int> count(a.size());

    cout<<"entering while loop\n";

    int tnum;
    while(start||nonstop){
        if(start){
            start=false;
            for(i=0;i<g;i++){
                for(j=0;j<cyclelist[i].size();j++){
                    vector<int> t=cyclelist[i][j];
                    tnum=0;
                    for(k=0;k<=(i+1);k++)
                        tnum+=Pt[t[2*k]]-Pt[t[2*k+1]];
                    sum[i][j]=tnum;
                    if(!tnum)// other condition for other couplings
                        ncyc[i]++;
                }
            }
            print_int(ncyc);
        }

        nonstop=false;
        for(e=0;e<elist.size();e++){
            v_cur=v_opt=Pt[e];
            diff_opt=0;
            for(i=0;i<g;i++)
                d_opt[i]=0;
            for(k=0;k<a.size();k++){
                v_temp=a[k];
                v_diff=v_temp-v_cur;
                if(v_diff){
                    diff_temp=0;
                    for(i=0;i<g;i++){
                        d_temp[i]=0;
                        vector<int> t1=edge_cyclelist[e][2*i];
                        vector<int> t2=edge_cyclelist[e][2*i+1];
                        for(j=0;j<t1.size();j++){
                            s_new=s_old=sum[i][t1[j]];
                            s_new+=((t2[j]>0)? v_diff:(-v_diff));
                            d_temp[i]+=((s_new)? 0:1)-((s_old)? 0:1);// for non tail-biting SC codes. change to modulo L for tail-biting codes, add necessary modulo L parts of s_new/old for rotate TC code
                        }
                        diff_temp+=(w[i]*d_temp[i]);//  for non-tail-biting SC codes or Rotate TC, involve L and multiply further by the multiplicities
                    }
                    if((d_temp[0]<d_opt[0])||((d_temp[0]<=d_opt[0])&&(diff_temp<diff_opt))) {// the condition could be changed to ensure d[i] to be 0 if it is in the previous round
                        bool edit= edit_norm(count,inda[v_cur],k,th1,th2);
                        //cout<<"edit!\n";
                        if(edit){
                            v_opt = a[k];
                            diff_opt=diff_temp;
                            d_opt=d_temp;
                        }
                    }
                }
            }
            v_diff=v_opt-v_cur;
            if(v_diff) {
                nonstop = true;
                Pt[e] = v_opt;
                count[inda[v_opt]]++;
                count[inda[v_cur]]--;
                for(i=0;i<g;i++) {
                    vector<int> t1 = edge_cyclelist[e][2 * i];
                    vector<int> t2 = edge_cyclelist[e][2 * i + 1];
                    for (j = 0; j < t1.size(); j++) {
                        s_new=s_old=sum[i][t1[j]];
                        s_new+=((t2[j]>0)? v_diff:(-v_diff));
                        sum[i][t1[j]]=s_new;
                        ncyc[i]+=((s_new)? 0:1)-((s_old)? 0:1);
                    }
                }
            }
        }
    }

    print_int(ncyc);
    return Pt;
}


// unlike in non-RTC codes, P here should record the index instead of a (i instead of a_i)
vector<int> comcodesc::AO_RTC(SC pattern,vector<int> P,int th1,int th2,int L){
    vector<int> Pt(P);

    int i,j,k,l,e,s_old,s_new,v_diff,v_temp,v_opt,v_cur,diff_temp,diff_opt;
    bool start=true;
    bool nonstop=true;
    int g=cyclelist.size();
    vector<int> ncyc(elist.size());
    int Ncyc=0;
    for(i=0;i<g;i++){
        vector<int> t(cyclelist[i].size());
    }
    int m= pattern.m;
    vector<int> a=pattern.a;
    vector<int> inda=pattern.inda;
    vector<int> count(pattern.mt+1);

    vector<int> numcyc6(cyclelist[1].size());
    vector<int> numcyc6_temp=numcyc6;
    vector<int> numcyc6_opt=numcyc6;

    int eid;

    cout<<"entering while loop\n";

    int tnum;
    while(start||nonstop){
        if(start){
            start=false;
            for(j=0;j<cyclelist[1].size();j++){
                vector<int> t=cyclelist[i][j];
                tnum=pattern.count_RTC(Pt[t[0]],Pt[t[1]],Pt[t[2]],Pt[t[3]],Pt[t[4]],Pt[t[5]],L);
                numcyc6[j]=tnum;
                if(!tnum)
                    Ncyc+=tnum;
            }
            cout<<"Initial num of cycle 6="<<Ncyc<<endl;
        }

        nonstop=false;
        for(e=0;e<elist.size();e++){
            v_cur=v_opt=Pt[e];
            numcyc6_opt=numcyc6;
            diff_opt=0;
            for(k=0;k<a.size();k++){
                if(k!=v_cur){
                    diff_temp=0;
                    numcyc6_temp=numcyc6_opt;
                    vector<int> t1=edge_cyclelist[e][2];// cycle indices contain the edge
                    vector<int> t2=edge_cyclelist[e][3];// corresponding position within the cycle
                    for(j=0;j<t1.size();j++){
                        eid=t1[j];
                        vector<int> cyc_e=cyclelist[1][eid];
                        cyc_e[t2[j]]=k;
                        numcyc6_temp[eid]=pattern.count_RTC(Pt[cyc_e[0]],Pt[cyc_e[1]],Pt[cyc_e[2]],Pt[cyc_e[3]],Pt[cyc_e[4]],Pt[cyc_e[5]],L);
                        diff_temp+=numcyc6_temp[eid]-numcyc6[eid];// for non tail-biting SC codes. change to modulo L for tail-biting codes, add necessary modulo L parts of s_new/old for rotate TC code
                    }
                    if(diff_temp<diff_opt) {// the condition could be changed to ensure d[i] to be 0 if it is in the previous round
                        bool edit= edit_norm(count,v_cur,k,th1,th2);
                        if(edit){
                            v_opt = k;
                            diff_opt=diff_temp;
                            numcyc6_opt=numcyc6_temp;
                        }
                    }
                }
            }
            if(v_opt!=v_cur) {
                nonstop = true;
                numcyc6=numcyc6_opt;
                Pt[e] = v_opt;
                count[v_opt]++;
                count[v_cur]--;
                Ncyc+=diff_opt;
            }
        }
    }

    cout<<"Final number of cycle 6="<<Ncyc<<endl;

    return Pt;
}




#endif //QLDPC_SC_OPTIMIZATION_H

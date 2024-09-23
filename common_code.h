//
// Created by Siyi Yang on 1/2/24.
//

#ifndef QLDPC_COMMON_CODE_H
#define QLDPC_COMMON_CODE_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"

using namespace std;

class comcode{
public:
    vector<vector<int> >  pcm;// H
    vector<vector<int> >  pcme;// H with 1'th replaced by edge indices
    vector<vector<vector<int> > > vneighbor;// nv*2*~, [v][0]: ne; [v][1]: j
    vector<vector<vector<int> > > cneighbor;// nc*2*~, [c][0]: ne; [c][1]: i
    vector<vector<int> >  elist;// E_num*2  e=(i,j)
    comcode(){vector<vector<int> >  h(0);vector<vector<int> >  he(0);vector<vector<vector<int> > > vneigh(0);vector<vector<vector<int> > > cneigh(0);vector<vector<int> >  el(0);pcm=h;pcme=he;vneighbor=vneigh;cneighbor=cneigh;elist=el;}
    void init_code(vector<vector<int> >  H);
    vector<vector<int> >  v2h(vector<int> v);
    vector<int> h2v(vector<vector<int> >  h);
    vector<int> get_vsyndrome(vector<int> x);
    vector<int> get_csyndrome(vector<int> x);
};

vector<int> comcode::get_vsyndrome(vector<int> x){
    int i,j,k;
    int N=vneighbor.size();
    int M=x.size();
    vector<int> t;
    vector<int> y(N);
    for(i=0;i<N;i++){
        k=0;
        t=vneighbor[i][1];
        for(j=0;j<t.size();j++)
            k=k^x[t[j]];
        y[i]=k;
    }
    return y;
}

vector<int> comcode::get_csyndrome(vector<int> x){
    int i,j,k;
    int N=cneighbor.size();
    int M=x.size();
    vector<int> t;
    vector<int> y(N);
    for(i=0;i<N;i++){
        k=0;
        t=cneighbor[i][1];
        for(j=0;j<t.size();j++)
            k=k^x[t[j]];
        y[i]=k;
    }
    return y;
}


void comcode::init_code(vector<vector<int> >  H) {
    int i,j;
    int ne=0;
    int nc=H.size();
    int nv=H[0].size();
    vector<vector<int> >  Ht=H;
    vector<vector<vector<int> > > vneigh(nv);
    vector<vector<vector<int> > > cneigh(nc);
    vector<vector<int> >  et(0);
    for(i=0;i<nc;i++) {
        vector<vector<int> >  ct(2);
        for (j = 0; j < nv; j++)
            if (H[i][j]) {
                ct[0].push_back(ne);
                ct[1].push_back(j);
                vector<int> t {i,j};
                et.push_back(t);
                Ht[i][j]=ne;
                ne++;
            }
        cneigh[i]=ct;
    }

    for(j=0;j<nv;j++){
        vector<vector<int> >  vt(2);
        for(i=0;i<nc;i++)
            if(H[i][j]){
                vt[0].push_back(Ht[i][j]);
                vt[1].push_back(i);
            }
        vneigh[j]=vt;
    }

    pcm=H;
    pcme=Ht;
    vneighbor=vneigh;
    cneighbor=cneigh;
    elist=et;
}

vector<vector<int> >  comcode::v2h(vector<int> v){
    int i,e;
    vector<vector<int> >  h(cneighbor.size());
    for(i=0;i<h.size();i++){
        vector<int> t(vneighbor.size());
        h[i]=t;
    }
    for(e=0;e<v.size();e++)
        h[elist[e][0]][elist[e][1]]=v[e];
    return h;
}

vector<int> comcode::h2v(vector<vector<int> >  h){
    int e;
    vector<int> v(elist.size());
    for(e=0;e<v.size();e++)
        v[e]=h[elist[e][0]][elist[e][1]];
    return v;
}




#endif //QLDPC_COMMON_CODE_H

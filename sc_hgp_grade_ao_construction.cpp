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

int main(){
    int r1=3;
    int r2=3;
    int n1=7;
    int n2=7;
    int m1=2;
    int m2=2;
    double w=25;
    int L1=10;
    int L2=10;
    double step=0.02;

    int r=r1*n2+r2*n1;
    int n=n1*n2+r1*r2;

    vector<vector<double> >  we=get_weight(r1,r2,n1,n2);
    vector<vector<int> >  a= initialize(m1,m2);
    vector<vector<int> >  b=a;

    MD md(a,b);
    vector<vector<double> >  q=md.Q_GRADE(w,we[0],we[1],L1,L2,step);
    print_double(q[0]);
    print_double(q[1]);

    int lent=q[0].size();
    vector<double> p_uni_a(lent);
    double unit=1.0/((double) lent);
    for(int i=0;i<lent;i++)
        p_uni_a[i]=unit;
    vector<double> p_uni_b=p_uni_a;

    vector<int> L{L1,L2};
    int len1=r1*n1;
    int len2=r2*n2;

    vector<vector<int> >  Ha(r1);
    vector<vector<int> >  Hb(r2);
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

    vector<vector<int> >  coun_uni=p2count(p_uni_a,len1,(m1+1)*(m2+1)-1);
    vector<int> Pa_uni=coun_uni[1];
    print_int(Pa_uni);
    coun_uni=p2count(p_uni_b,len1,(m1+1)*(m2+1)-1);
    vector<int> Pb_uni=coun_uni[1];
    print_int(Pb_uni);

    vector<vector<int> >  coun=p2count(q[0],len1,(m1+1)*(m2+1)-1);
    vector<int> Pa=coun[1];
    print_int(Pa);
    coun=p2count(q[1],len2,(m1+1)*(m2+1)-1);
    vector<int> Pb=coun[1];
    print_int(Pb);

    vector<vector<int> >  P=Hqsc.AO(a,b,Pa,Pb,L,8,2,w);
    cout<<"finished P_gd ao\n";
    print_matrix(P);
    vector<vector<int> >  P1=Hqsc.AO(a,b,P[0],P[1],L,8,2,w);
    cout<<"check P_gd ao\n";
    print_matrix(P1);
    vector<vector<int> >  P2=Hqsc.AO(a,b,P1[0],P1[1],L,8,2,w);
    cout<<"check P1_gd ao\n";
    print_matrix(P2);
    vector<vector<int> >  P_uni=Hqsc.AO(a,b,Pa_uni,Pb_uni,L,8,2,w);
    cout<<"finished P_uni ao\n";
    print_matrix(P_uni);
    vector<vector<int> >  P1_uni=Hqsc.AO(a,b,P_uni[0],P_uni[1],L,8,2,w);
    cout<<"Check P_uni ao\n";
    print_matrix(P1_uni);
    vector<vector<int> >  P2_uni=Hqsc.AO(a,b,P1_uni[0],P1_uni[1],L,8,2,w);
    cout<<"Check P1_uni ao\n";
    print_matrix(P2_uni);

    return 0;
}
//
// Created by Siyi Yang on 1/2/24.
//

#ifndef QLDPC_MD_H
#define QLDPC_MD_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"

using namespace std;

class MD{
public:
    vector<vector<int> >  a;
    vector<vector<int> >  b;
    vector<vector<int> >  inda;//e.g. [0,1,0,0,2,0,3] for [0,1,4,6], inda[a[k]]=k
    vector<vector<int> >  indb;//e.g. [0,1,0,0,2,0,3] for [0,1,4,6], inda[a[k]]=k
    int m1;
    int m2;
    int mta;
    int mtb;
    MD(vector<vector<int> >  at,vector<vector<int> >  bt);

    //vector<int> count_freq(vector<int> P);
    vector<vector<vector<int> > > ao_init(vector<int> inia,vector<int> inib);
    vector<vector<double> >  Q_GRADE(double w,vector<double> wa, vector<double> wb,int L1,int L2,double step);
    vector<vector<int> >  SC(vector<vector<int> >  P,vector<vector<int> >  Ea,vector<vector<int> >  Eb,int n1,int n2,int r1,int r2,int L1,int L2);
    void print(){
        cout<<"a:"<<endl;
        print_matrix(a);
        cout<<"inda:"<<endl;
        print_matrix(inda);
        cout<<"b:"<<endl;
        print_matrix(b);
        cout<<"indb:"<<endl;
        print_matrix(indb);
        cout<<"m1="<<m1<<endl;
        cout<<"m2="<<m2<<endl;
        cout<<"mta="<<mta<<endl;
        cout<<"mtb="<<mtb<<endl;
    }
};

MD::MD(vector<vector<int> >  at,vector<vector<int> >  bt){
    a=at;
    b=bt;
    mta=a[0].size()-1;
    mtb=b[0].size()-1;
    int m1a=vec_max(a[0]);//a[0]: m1
    int m2a=vec_max(a[1]);//a[1]: m2
    int m1b=vec_max(b[0]);//a[0]: m1
    int m2b=vec_max(b[1]);//a[1]: m2
    m1=(m1a>m1b)? m1a:m1b;
    m2=(m2a>m2b)? m2a:m2b;
    vector<vector<int> >  ind1(m1+1);
    vector<vector<int> >  ind2(m1+1);
    for(int i=0;i<=m1;i++){
        vector<int> t(m2+1);
        ind1[i]=t;
        ind2[i]=t;
    }
    for(int k=0;k<a[0].size();k++){
        ind1[a[0][k]][a[1][k]]=k;
    }
    inda=ind1;
    for(int k=0;k<b[0].size();k++){
        ind2[b[0][k]][b[1][k]]=k;
    }
    indb=ind2;
}

vector<vector<int> >  MD::SC(vector<vector<int> >  P,vector<vector<int> >  Ea,vector<vector<int> >  Eb,int n1,int n2,int r1,int r2,int L1,int L2) {
    int i,ea, eb, ia, ib, ja, jb, j1, j2, ix, iz, l1, l2, v1, v2, cx, cz, m1at, m2at, m1bt, m2bt, p1, p2;

    int r1n2 = r1 * n2;
    int n1n2 = n1 * n2;
    int r = r1n2 + r2 * n1;
    int n = n1n2 + r1 * r2;
    int nL1 = n * L1;
    int rL1 = r * L1;
    int L1L2 = L1 * L2;
    int num_cn = rL1 * L2;
    int num_vn = nL1 * L2;

    vector<vector<int> >  H(num_cn);
    vector<int> t(num_vn);
    for (i = 0; i<num_cn; i++)
        H[i] = t;

    for(ea=0;ea<Ea.size();ea++) {
        ia = Ea[ea][0];
        ja = Ea[ea][1];
        p1 = P[0][ea];
        m1at = a[0][p1];
        m2at = a[1][p1];
        for(i=0;i<n2;i++){
            ix=i*r1+ia;
            j1=i*n1+ja;
            for(l1=0;l1<L1;l1++)
                for(l2=0;l2<L2;l2++) {
                    v1=l2*nL1+l1*n+j1;
                    cx=((l2+m2at)%L2)*rL1+((l1+m1at)%L1)*r+ix;
                    H[cx][v1]=1;
                }
        }
        for(i=0;i<r2;i++){
            iz=r1n2+i*n1+ja;
            j2=n1n2+i*r1+ia;
            for(l1=0;l1<L1;l1++)
                for(l2=0;l2<L2;l2++) {
                    v2=l2*nL1+l1*n+j2;
                    cz=((l2+m2-m2at)%L2)*rL1+((l1+m1-m1at)%L1)*r+iz;
                    H[cz][v2]=3;
                }
        }
    }

    for(eb=0;eb<Eb.size();eb++) {
        ib = Eb[eb][0];
        jb = Eb[eb][1];
        p2 = P[1][eb];
        m1bt = b[0][p2];
        m2bt = b[1][p2];
        for(i=0;i<n1;i++){
            j1=jb*n1+i;
            iz=r1n2+ib*n1+i;
            for(l1=0;l1<L1;l1++)
                for(l2=0;l2<L2;l2++) {
                    v1=l2*nL1+l1*n+j1;
                    cz=((l2+m2bt)%L2)*rL1+((l1+m1bt)%L1)*r+iz;
                    H[cz][v1]=3;
                }
        }
        for(i=0;i<r1;i++){
            j2=n1n2+ib*r1+i;
            ix=jb*r1+i;
            for(l1=0;l1<L1;l1++)
                for(l2=0;l2<L2;l2++) {
                    v2=l2*nL1+l1*n+j2;
                    cx=((l2+m2-m2bt)%L2)*rL1+((l1+m1-m1bt)%L1)*r+ix;
                    H[cx][v2]=1;
                }
        }
    }
    return H;
}

vector<vector<vector<int> > > MD::ao_init(vector<int> inia,vector<int> inib){
    vector<vector<int> >  pa(inia.size());
    vector<vector<int> >  pb(inib.size());
    for(int i=0;i<inia.size();i++) {
        vector<int> t {a[0][inia[i]],a[1][inia[i]]};
        pa[i]=t;
    }
    for(int i=0;i<inib.size();i++) {
        vector<int> t {b[0][inib[i]],b[1][inib[i]]};
        pb[i]=t;
    }
    vector<vector<vector<int> > > p(2);
    p[0]=pa;
    p[1]=pb;
    return p;
}

vector<vector<double> >  MD::Q_GRADE(double w,vector<double> wa, vector<double> wb,int L1,int L2,double step){
    int i,j,k;
    vector<vector<double> >  p(2);
    double v_prev=0.0;
    double v_cur=10.0;
    double eps=0.5;

    vector<double> pa(mta+1);
    vector<double> pb(mtb+1);
    for(i=0;i<=mta;i++)
        pa[i]=((double) 1)/((double) mta+1);
    for(i=0;i<=mtb;i++)
        pb[i]=((double) 1)/((double) mtb+1);

    vector<double> ta(mta+1);
    vector<double> tb(mtb+1);
    int count=0;
    bool is_first_round= true;
    vector<vector<double> >  fa,fb,fa_bar,fb_bar;
    vector<vector<double> >  ga,gb,ga2,gb2,ga3,gb3,ga4,gb4,gab;
    vector<vector<double> >  ha3,hb3,ha4,hb4,ha,hb;
    vector<double> grada,gradb;
    double f;

    double cga3=w*wa[0]*wb[2]+30*wa[0]*wb[3];
    double cgb3=w*wa[2]*wb[0]+30*wa[3]*wb[0];
    double cga4=wa[1]*wb[2];
    double cgb4=wa[2]*wb[1];
    double cgab=248.0;

    double cha3=6*cga3;
    double chb3=6*cgb3;
    double cha4=8*cga4;
    double chb4=8*cgb4;
    double chab=4*cgab;

    int m14=4*m1;
    int m24=4*m2;
    int m13=3*m1;
    int m23=3*m2;
    //cout<<"finished weight\n";

    while((abs(v_cur-v_prev)>eps)&&(count<100)){
        fa=new_matrix_index(pa,a,m1,m2);
        fa_bar= flip2(fa);
        ga=conv2(fa,fa_bar);
        ga2=conv2(ga,ga);
        ga3=conv2(ga2,ga2);
        ga4=conv2(ga3,ga3);
        ha3=conv2(ga2,fa);
        ha4=conv2(ga3,fa);
        //cout<<"finished a\n";

        fb=new_matrix_index(pb,b,m1,m2);
        fb_bar= flip2(fb);
        gb=conv2(fb,fb_bar);
        gb2=conv2(gb,gb);
        gb3=conv2(gb2,gb2);
        gb4=conv2(gb3,gb3);
        hb3=conv2(gb2,fb);
        hb4=conv2(gb3,fb);
        //cout<<"finished b\n";

        ha=conv2(gb2,conv2(ga,fa));
        hb=conv2(ga2,conv2(gb,fb));
        gab=conv2(ga2,gb2);
        //cout<<"finished ab\n";

        f= cgab*index_sum(gab,L1,L2,m14,m24);
        f+=cga3*index_sum(ga3,L1,L2,m13,m23);
        f+=cga4*index_sum(ga4,L1,L2,m14,m24);
        f+=cgb3*index_sum(gb3,L1,L2,m13,m23);
        f+=cgb4*index_sum(gb4,L1,L2,m14,m24);

        //cout<<"f="<<f<<endl;
        grada= mul(cha3,index_sum_grad(ha3,m1,m2,L1,L2,m13,m23,a));
        gradb=mul(chb3,index_sum_grad(hb3,m1,m2,L1,L2,m13,m23,b));
        grada=add(grada, index_sum_grad(ha4,m1,m2,L1,L2,m14,m24,a),cha4);
        gradb=add(gradb, index_sum_grad(hb4,m1,m2,L1,L2,m14,m24,b),chb4);
        grada=add(grada, index_sum_grad(ha,m1,m2,L1,L2,m14,m24,a),chab);
        gradb=add(gradb, index_sum_grad(hb,m1,m2,L1,L2,m14,m24,b),chab);
        grada= normalize(grada);
        gradb= normalize(gradb);
        //print_double(grada);
        //print_double(gradb);

        //cout<<"grada"<<endl;

        v_prev=v_cur;
        v_cur=f;

        if(is_first_round) {
            eps=v_cur/1e6;
            is_first_round= false;
        }
        else{
            if(v_cur>v_prev)
                step=step/2;
        }

        if(abs(v_cur-v_prev)>eps){
            //cout<<"update\n";
            pa=add(pa,grada,-step);
            //print_double(pa);
            pb=add(pb,gradb,-step);
            //print_double(pb);
        }
        count++;
    }

    p[0]=pa;
    p[1]=pb;

    return p;
}


#endif //QLDPC_MD_H

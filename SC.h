//
// Created by Siyi Yang on 1/2/24.
//

#ifndef QLDPC_SC_H
#define QLDPC_SC_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"

using namespace std;

class SC{
public:
    vector<int> a;//coupling pattern, e.g. [0,1,4,6]
    vector<int> b;//e.g. [0,0,2,3] for [0,1,4,6]
    vector<int> inda;//e.g. [0,1,0,0,2,0,3] for [0,1,4,6], inda[a[k]]=k
    int m;
    int mt;
    SC(vector<int> at);
    vector<int> count_freq(vector<int> P);
    vector<int> ao_init(vector<int> ini);
    vector<double> GRADE(int g,double step);
    vector<double> GRADE_RTC(int L,double step);// cycle 6 only
    int count_RTC(int i1,int i2,int i3,int i4,int i5,int i6,int L);
    void print(){
        print_int(a);
        print_int(b);
        print_int(inda);
        cout<<"m="<<m;
        cout<<"mt="<<mt;
    }
};

SC::SC(vector<int> at) {
    a=at;
    mt=a.size()-1;
    vector<int> t(a);
    for(int i=0;i<a.size();i++)
        t[i]=t[i]-i;
    b=t;
    m= vec_max(a);
    vector<int> ind(m+1);
    for(int k=0;k<a.size();k++)
        ind[a[k]]=k;
    inda=ind;
}

vector<int> SC::count_freq(vector<int> P){
    vector<int> count(mt+1);
    for(int i=0;i<P.size();i++)
        count[inda[P[i]]]++;
    return count;
}

vector<int> SC::ao_init(vector<int> ini){
    vector<int> P(ini);
    for(int i=0;i<P.size();i++)
        P[i]=a[ini[i]];
    return P;
}

// d girth/2
// L coupling length
vector<double> SC::GRADE(int d,double step){
    int i,j,k;
    double v_prev=0.0;
    double v_cur=1.0;
    double eps=0.5;
    vector<double> p(mt+1);
    double sum=0;
    for(i=0;i<=mt;i++) {
        //p[i] = ((double) 1) / ((double) mt + 1);
        p[i]=((double) rand()/(RAND_MAX));
        sum+=p[i];
    }
    for(i=0;i<=mt;i++)
        p[i]/=sum;
    vector<double> t(m+1);
    int count=0;
    bool is_first_round= true;
    vector<double> f(m+1);
    vector<double> f_bar(m+1);
    vector<double> h0,h,g,grad;

    while((abs(v_cur-v_prev)>eps)&&(count<100)){
        for(i=0;i<f.size();i++)
            f[i]=0;
        for(i=0;i<a.size();i++)
            f[a[i]]=p[i];
        f_bar= flip1(f);
        h0=conv1(f,f_bar);
        h=h0;
        g=mul(2*d,f);
        for(i=1;i<d;i++){
            h=conv1(h,h0);
            g=conv1(g,h0);
        }
        v_prev=v_cur;
        v_cur=h[d*m];
        grad= indexing(g,(d-1)*m,a);
        //grad= normalize(grad);
        if(is_first_round) {
            eps=v_cur/1e6;
            is_first_round= false;
            /*
            beta=2.0*v_cur/(mt+1);
            for(int tt=0;tt<=mt;tt++) {
                double xx=tanh(15.0 * p[tt]);
                v_cur += beta *(1.0+xx);
                grad[tt]+=15.0*beta*(1.0-xx*xx);
            }
             */
        }
        else{
            /*
            for(int tt=0;tt<=mt;tt++) {
                double xx=tanh(15.0 * p[tt]);
                v_cur += beta *(1.0+xx);
                grad[tt]+=15.0*beta*(1.0-xx*xx);
            }
             */
            if(v_cur>v_prev)
                step=step/2;
            else
                step=step;
        }
        //cout<<"v_cur="<<v_cur<<endl;
        grad= normalize(grad);
        p=add(p,grad,-step);
        count++;
    }
    cout<<v_cur<<endl;
    //cout<<"total "<<count<<" rounds\n";
    return p;
}


vector<double> SC::GRADE_RTC(int L,double step) {
    int i1,i2,i3,i4,i5,i6,st;
    vector<vector<int> >  solv(0);
    vector<double> w(0);
    int num_cyc;
    for(i1=0;i1<=mt;i1++){
        for(i2=0;i2<=mt;i2++){
            for(i3=0;i3<=mt;i3++) {
                for (i4 = 0; i4 <= mt; i4++) {
                    for (i5 = 0; i5 <= mt; i5++) {
                        for (i6 = 0; i6 <= mt; i6++) {
                            //st= count_RTC(i1,i2,i3,i4,i5,i6,L);
                            st= count_RTC(i1,i2,i3,i4,i5,i6,L);//+count_RTC(i2,i1,i6,i5,i4,i3,L)+count_RTC(i3,i4,i5,i6,i1,i2,L)+count_RTC(i4,i3,i2,i1,i6,i5,L)+count_RTC(i5,i6,i1,i2,i3,i4,L)+count_RTC(i6,i5,i4,i3,i2,i1,L);
                            if(st){
                                vector<int> solvt {i1,i2,i3,i4,i5,i6};
                                solv.push_back(solvt);
                                w.push_back((double) st);
                                //num_cyc+=st;
                                /*
                                int xxx=0;
                                for(int i=0;i<6;i++) {
                                    cout << a[solvt[i]] << ',';
                                    xxx+=((i%2)? a[solvt[i]]:(-a[solvt[i]]));
                                }
                                cout<<xxx<<','<<st<<endl;
                                 */
                            }
                        }
                    }
                }
            }
        }
    }

    //cout<<"num of cycle="<<num_cyc<<endl;

    int i,j,k;
    double v_prev=0.0;
    double v_cur=1.0;
    double eps=0.5;
    vector<double> p(mt+1);
    double sum=0;
    for(i=0;i<=mt;i++) {
        //p[i] = ((double) 1) / ((double) mt + 1);
        p[i]=((double) rand()/(RAND_MAX));
        sum+=p[i];
    }
    for(i=0;i<=mt;i++)
        p[i]/=sum;

    int count=0;
    bool is_first_round= true;
    double f,ft,gt;
    vector<double> g(mt+1);
    while((abs(v_cur-v_prev)>eps)&&(count<100)) {
        f=0.0;
        for(i=0;i<=mt;i++)
            g[i]=0.0;
        for(i=0;i<solv.size();i++){
            vector<int> ind=solv[i];
            ft=w[i];
            for(j=0;j<6;j++){
                ft*=p[ind[j]];
                gt=w[i];
                for(k=0;k<6;k++)
                    if(k!=j)
                        gt*=p[ind[k]];
                g[ind[j]]+=gt;
            }
            f+=ft;
        }

        v_prev=v_cur;
        v_cur=f;
        if(is_first_round) {
            eps=v_cur/1e6;
            is_first_round= false;
        }
        else{
            if(v_cur>v_prev)
                step=step/2;
            else
                step=step;
        }
        g= normalize(g);
        p=add(p,g,-step);
        count++;
    }

    cout<<v_cur<<endl;
    return p;
}


int SC::count_RTC(int i1,int i2,int i3,int i4,int i5,int i6,int L){
    int N=0;
    int sum=a[i1]-a[i2]+a[i3]-a[i4]+a[i5]-a[i6];
    int min45,max45,min23,max23,min16,max16,d65,d56,d12,d21,d34,d43,lb,ub,t;
    if(!(sum%L)){
        if(sum==(-2*L)) {
            t=L-a[i6]+a[i1]+i5-i2;
            if(t>0)
                N+=t;
            t=L-a[i2]+a[i3]+i1-i4;
            if(t>0)
                N+=t;
            t=L-a[i4]+a[i5]+i3-i6;
            if(t>0)
                N+=t;
            t=L-a[i1]+a[i6]+i2-i5;
            if(t>0)
                N+=t;
            t=L-a[i3]+a[i2]+i4-i1;
            if(t>0)
                N+=t;
            t=L-a[i5]+a[i4]+i6-i3;
            if(t>0)
                N+=t;
        }
        else{
            if(b[i1]>b[i6]){
                min16=b[i6];
                max16=b[i1];
            }
            else{
                min16=b[i1];
                max16=b[i6];
            }
            if(b[i2]>b[i3]){
                min23=b[i3];
                max23=b[i2];
            }
            else{
                min23=b[i2];
                max23=b[i3];
            }
            if(b[i4]>b[i5]){
                min45=b[i5];
                max45=b[i4];
            }
            else{
                min45=b[i4];
                max45=b[i5];
            }
            d65=a[i6]-a[i5];
            d12=a[i1]-a[i2];
            d34=a[i3]-a[i4];
            d56=-d65;
            d21=-d12;
            d43=-d34;

            if(sum==L){
                if(i1>i6){
                    t=min3(L-max45-d65,2*L-max23-d12,L-b[i6])-max3(-min45-d65,L-min23-d12,L-b[i1]);
                    if(t>0)
                        N+=t;
                }
                if(i1<i6){
                    t=min3(-b[i4]-d65,i2-a[i1],-b[i1])-max3(-b[i3]-d12,i5-a[i6],-b[i6]);
                    if(t>0)
                        N+=t;
                }

                if(i3>i2){
                    t=min3(L-max16-d21,2*L-max45-d34,L-b[i2])-max3(-min16-d21,L-min45-d34,L-b[i3]);
                    if(t>0)
                        N+=t;
                }
                if(i3<i2){
                    t=min3(-b[i6]-d21,i4-a[i3],-b[i3])-max3(-b[i5]-d34,i1-a[i2],-b[i2]);
                    if(t>0)
                        N+=t;
                }

                if(i5>i4){
                    t=min3(L-max23-d43,2*L-max16-d56,L-b[i4])-max3(-min23-d43,L-min16-d56,L-b[i5]);
                    if(t>0)
                        N+=t;
                }
                if(i5<i4){
                    t=min3(-b[i2]-d43,i6-a[i5],-b[i5])-max3(-b[i1]-d56,i3-a[i4],-b[i4]);
                    if(t>0)
                        N+=t;
                }
            }

            if(sum==(-L)){
                if(i1<i6){
                    t=min3(L-max45-d12,2*L-max23-d65,L-b[i1])-max3(-min45-d12,L-min23-d65,L-b[i6]);
                    if(t>0)
                        N+=t;
                }
                if(i1>i6){
                    t=min3(-b[i3]-d12,i5-a[i6],-b[i6])-max3(-b[i4]-d65,i2-a[i1],-b[i1]);
                    if(t>0)
                        N+=t;
                }

                if(i3<i2){
                    t=min3(L-max16-d34,2*L-max45-d21,L-b[i3])-max3(-min16-d34,L-min45-d21,L-b[i2]);
                    if(t>0)
                        N+=t;
                }
                if(i3>i2){
                    t=min3(-b[i5]-d34,i1-a[i2],-b[i2])-max3(-b[i6]-d21,i4-a[i3],-b[i3]);
                    if(t>0)
                        N+=t;
                }

                if(i5<i4){
                    t=min3(L-max23-d56,2*L-max16-d43,L-b[i5])-max3(-min23-d56,L-min16-d43,L-b[i4]);
                    if(t>0)
                        N+=t;
                }
                if(i5>i4){
                    t=min3(-b[i1]-d56,i3-a[i4],-b[i4])-max3(-b[i2]-d43,i6-a[i5],-b[i5]);
                    if(t>0)
                        N+=t;
                }
            }

            if(!sum){
                t=min3(L-max16,L-max23-d12,L-max45-d65)-max3(-min16,-min23-d12,-min45-d65);
                N+=((t>0)? t:0);
                if((i2<i3)&&(i5<i4)) {
                    t = min3(L - max16, L - b[i2] - d12, L - b[i5] - d65) -
                        max3(-min16, L - b[i3] - d12, L - b[i4] - d65);
                    N += ((t > 0) ? t : 0);
                }
                if((i4<i5)&&(i1<i6)) {
                    t = min3(L - max23, L - b[i4] - d34, L - b[i1] - d21) -
                        max3(-min23, L - b[i5] - d34, L - b[i6] - d21);
                    N += ((t > 0) ? t : 0);
                }
                if((i6<i1)&&(i3<i2)) {
                    t = min3(L - max45, L -b[i6]- d56, L - b[i3] - d43) - max3(-min45, L -b[i1]- d56, L - b[i2] - d43);
                    N += ((t > 0) ? t : 0);
                }
                if((i2>i3)&&(i5>i4)) {
                    t = min3(L - max16, -b[i3] - d12, -b[i4] - d65) - max3(-min16, -b[i2] - d12, -b[i5] - d65);
                    N += ((t > 0) ? t : 0);
                }
                if((i4>i5)&&(i1>i6)) {
                    t = min3(L - max23, -b[i5] - d34, -b[i6] - d21) - max3(-min23, -b[i4] - d34, -b[i1] - d21);
                    N += ((t > 0) ? t : 0);
                }
                if((i6>i1)&&(i3>i2)) {
                    t = min3(L - max45, -b[i1]- d56, - b[i2] - d43) - max3(-min45, -b[i6]- d56, - b[i3] - d43);
                    N += ((t > 0) ? t : 0);
                }
            }
        }
    }
    else
        N=0;
    return N;
}


#endif //QLDPC_SC_H

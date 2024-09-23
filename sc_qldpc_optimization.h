//
// Created by Siyi Yang on 1/2/24.
//

#ifndef QLDPC_SC_QLDPC_OPTIMIZATION_H
#define QLDPC_SC_QLDPC_OPTIMIZATION_H


#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>
#include "tools.h"
#include "common_code.h"
#include "sc_optimization.h"

using namespace std;


class comcodeqsc{
public:
    vector<comcodesc> H;
    int r1;
    int n1;
    int e1;
    int r2;
    int n2;
    int e2;
    comcodeqsc(vector<vector<int> >  Ha,vector<vector<int> >  Hb){
        comcodesc H0;
        H0.init_code(Ha);
        H0.cycle();
        comcodesc H1;
        H1.init_code(Hb);
        H1.cycle();
        vector<comcodesc> Ht(2);
        Ht[0]=H0;
        Ht[1]=H1;
        H=Ht;
        r1=H[0].pcm.size();
        n1=H[0].pcm[0].size();
        e1=H[0].elist.size();
        r2=H[1].pcm.size();
        n2=H[1].pcm[0].size();
        e2=H[1].elist.size();
        //cout<<"(r1,n1,e1)=("<<r1<<','<<n1<<','<<e1<<')'<<endl;
        //cout<<"(r2,n2,e2)=("<<r2<<','<<n2<<','<<e2<<')'<<endl;
        /*
        for(int r=0;r<2;r++)
            for(int i=0;i<3;i++)
                cout<<"H["<<r<<"].cyclelist["<<i<<"].size()="<<H[r].cyclelist[i].size()<<endl;
                */
    }
    void print(){
        //cout<<"start printing\n";
        for(int r=0;r<2;r++) {
            //cout<<"printing r="<<r<<endl;
            for (int i = 0; i < 3; i++) {
                //cout<<"printing i="<<i<<endl;
                cout << "H[" << r << "].cyclelist[" << i << "].size()=" << H[r].cyclelist[i].size() << endl;
            }
        }
    }
    vector<int> num_cycle(vector<vector<int> >  ncyc,vector<vector<vector<int> > > ncyc4,vector<int> L,vector<vector<int> >  m);
    vector<vector<int> >  AO(vector<vector<int> >  a,vector<vector<int> >  b,vector<int> Pa,vector<int> Pb,vector<int> L,int th1,int th2,double w);
    void  sample_list(vector<vector<int> >  a,vector<vector<int> >  b,vector<int> Pa,vector<int> Pb,vector<int> L);
};

vector<int> comcodeqsc::num_cycle(vector<vector<int> >  ncyc,vector<vector<vector<int> > > ncyc4,vector<int> L,vector<vector<int> >  m){
    vector<int> n_cur(3);
    int i,j,k;
    n_cur[0]=ncyc4[0][0][0]*(r2+n2)+ncyc4[1][0][0]*(r1+n1);
    n_cur[1]=ncyc[0][0]*(r2+n2)+ncyc[1][0]*(r1+n1);
    n_cur[2]=ncyc[0][1]*(r2+n2)+ncyc[1][1]*(r1+n1)+30*(ncyc[0][0]*e2+ncyc[1][0]*e1);
    k=2*ncyc4[0][0][0]*ncyc4[1][0][0];
    for(i=0;i<=m[0][0];i++)
        for(j=0;j<=m[0][1];j++)
            if((i||j)&&(i<=(L[0]-i))&&(j<=(L[1]-j))){
                if((i<(L[0]-i))&&(j<(L[1]-j)))
                    k+=(ncyc4[0][i][j]+ncyc4[0][L[0]-1-i][L[1]-1-j])*(ncyc4[1][i][j]+ncyc4[1][L[0]-1-i][L[1]-1-j]);
                else
                    k+=2*ncyc4[0][i][j]*ncyc4[1][i][j];
            }
    n_cur[2]+=(128*k);
    return n_cur;
}

void  sample_list(vector<vector<int> >  a,vector<vector<int> >  b,vector<int> Pa,vector<int> Pb,vector<int> L){
    vector<vector<int> > Pt(2);
    Pt[0] = Pa;
    Pt[1] = Pb;
    vector<vector<vector<int> > > ab(2);
    ab[0] = a;
    ab[1] = b;

    int r,i,j,l,k,sum;

    vector<vector<int> > ncyc(2);// (a,b)*(6,8)
    vector<vector<vector<int> > > ncyc4(2);// (a,b)*(L1,L2) : s4
    vector<vector<int> > m(2);// (a,b)*(m1,m2)

    for (r = 0; r < 2; r++) {
        vector<int> t(2);
        t[0] = t[1] = 0;
        ncyc[r] = t;
        m[r] = t;
        for (i = 0; i < 2; i++)
            m[r][i] = vec_max(ab[r][i]);
        vector<vector<int> > tt(L[0]);
        for (i = 0; i < L[0]; i++) {
            vector<int> ttt(L[1]);
            for (j = 0; j < L[1]; j++)
                ttt[i] = 0;
            tt[i] = ttt;
        }
        ncyc4[r] = tt;
    }


    vector<int> v_temp(2);//(a,b)*(s1,s2)

    for (r = 0; r < 2; r++) {
        for (i = 0; i < 3; i++) {
            for (j = 0; j < H[r].cyclelist[i].size(); j++) {
                vector<int> t = H[r].cyclelist[i][j];
                v_temp[0] = v_temp[1] = 0;
                for (l = 0; l < 2; l++) {
                    for (k = 0; k <= (i + 1); k++)
                        v_temp[l] += (ab[r][l][Pt[r][t[2 * k]]] - ab[r][l][Pt[r][t[2 * k + 1]]]);
                    v_temp[l] = ((v_temp[l] >= 0) ? (v_temp[l] % L[l]) : (L[l] - 1 -
                                                                          ((-1 - v_temp[l]) % L[l])));
                    //cout<<r<<','<<l<<','<<i<<','<<j<<':'<<sum[r][l][i][j]<<endl;
                }
                if (i) { // g=6,8
                    if (!(v_temp[0] || v_temp[1]))
                        ncyc[r][i - 1]++;
                } else// g=4
                    ncyc4[r][v_temp[0]][v_temp[1]]++;
            }
        }
    }

}

vector<vector<int> >  comcodeqsc::AO(vector<vector<int> >  a,vector<vector<int> >  b,vector<int> Pa,vector<int> Pb,vector<int> L,int th1,int th2,double w) {
    bool start = true;
    bool nonstop = true;
    int g, r, i, j, k, l, e, v_opt, v_cur, N_opt, N_temp, N_cur;
    vector<int> n_opt(3);
    vector<int> n_temp(3);
    vector<int> n_cur(3);

    vector<vector<int> > Pt(2);
    Pt[0] = Pa;
    Pt[1] = Pb;
    vector<vector<vector<int> > > ab(2);
    ab[0] = a;
    ab[1] = b;
    cout << "P_init:\n";
    print_matrix(Pt);
    /*
    cout<<"a:\n";
    print_matrix(ab[0]);
    cout<<"b:\n";
    print_matrix(ab[1]);
     */

    vector<int> v_diff(2);//(a,b)*(s1,s2)
    vector<int> v_temp(2);//(a,b)*(s1,s2)
    vector<int> s_new(2);
    vector<int> s_old(2);

    vector<vector<int> > ncyc(2);// (a,b)*(6,8)
    vector<vector<vector<int> > > ncyc4(2);// (a,b)*(L1,L2) : s4
    vector<vector<int> > m(2);// (a,b)*(m1,m2)

    for (r = 0; r < 2; r++) {
        vector<int> t(2);
        t[0] = t[1] = 0;
        ncyc[r] = t;
        m[r] = t;
        for (i = 0; i < 2; i++)
            m[r][i] = vec_max(ab[r][i]);
        vector<vector<int> > tt(L[0]);
        for (i = 0; i < L[0]; i++) {
            vector<int> ttt(L[1]);
            for (j = 0; j < L[1]; j++)
                ttt[i] = 0;
            tt[i] = ttt;
        }
        ncyc4[r] = tt;
    }
    //cout<<"m:\n";
    //print_matrix(m);

    vector<vector<int> > ncyc_temp = ncyc;// (a,b)*(6,8)
    vector<vector<vector<int> > > ncyc4_temp = ncyc4;// (a,b)*(L1,L2) : s4
    vector<vector<int> > ncyc_opt = ncyc;// (a,b)*(6,8)
    vector<vector<vector<int> > > ncyc4_opt = ncyc4;// (a,b)*(L1,L2) : s4

    vector<vector<vector<int> > > ind(2);
    vector<vector<int> > count(2);
    for (r = 0; r < 2; r++) {
        vector<vector<int> > indt(m[r][0] + 1);//m1a,m1b
        for (i = 0; i <= m[r][0]; i++) {
            vector<int> indtt(m[r][1] + 1);
            indt[i] = indtt;
        }
        ind[r] = indt;

        vector<int> countt(ab[r][0].size());
        count[r] = countt;
    }
    //cout<<"Finished ind and count\n";

    vector<vector<vector<vector<int> > > > sum(2);// 2*(0,1)*g*num_cycle_2g
    for (r = 0; r < 2; r++) {
        vector<vector<vector<int> > > tsum(2);
        for (l = 0; l < 2; l++) {
            vector<vector<int> > sumt(3);// g*~
            for (i = 0; i < 3; i++) {
                vector<int> t(H[r].cyclelist[i].size());
                //cout<<"H["<<r<<"].cyclelist["<<i<<"].size()="<<H[r].cyclelist[i].size()<<endl;
                sumt[i] = t;
            }
            tsum[l] = sumt;
        }
        sum[r] = tsum;
    }
    vector<vector<vector<vector<int> > > > sum_temp = sum;
    vector<vector<vector<vector<int> > > > sum_opt = sum;
    //cout<<"Finished sum\n";

    vector<vector<int> > wab(2);
    for (r = 0; r < 2; r++) {
        vector<int> wabt(3);
        for (i = 0; i < 3; i++)
            wabt[i] = H[r].cyclelist[i].size();
        wab[r] = wabt;
    }
    //cout<<"Finished wab\n";

    while (start || nonstop) {
        if (start) {
            start = false;
            for (r = 0; r < 2; r++) {
                for (i = 0; i < 3; i++) {
                    for (j = 0; j < H[r].cyclelist[i].size(); j++) {
                        vector<int> t = H[r].cyclelist[i][j];
                        v_temp[0] = v_temp[1] = 0;
                        for (l = 0; l < 2; l++) {
                            for (k = 0; k <= (i + 1); k++)
                                v_temp[l] += (ab[r][l][Pt[r][t[2 * k]]] - ab[r][l][Pt[r][t[2 * k + 1]]]);
                            v_temp[l] = ((v_temp[l] >= 0) ? (v_temp[l] % L[l]) : (L[l] - 1 -
                                                                                  ((-1 - v_temp[l]) % L[l])));
                            sum[r][l][i][j] = v_temp[l];
                            //cout<<r<<','<<l<<','<<i<<','<<j<<':'<<sum[r][l][i][j]<<endl;
                        }
                        if (i) { // g=6,8
                            if (!(v_temp[0] || v_temp[1]))
                                ncyc[r][i - 1]++;
                        } else// g=4
                            ncyc4[r][v_temp[0]][v_temp[1]]++;
                    }
                }
            }
            print_matrix(ncyc);
            cout << "compute the total start number\n";
            n_cur = num_cycle(ncyc, ncyc4, L, m);
            N_cur = w * n_cur[1] + n_cur[2];
            cout << "n_cur start:\n";
            print_int(n_cur);
            cout << "N_cur start=" << N_cur << endl;
        }
        nonstop = false;

        for (r = 0; r < 2; r++) {
            for (e = 0; e < H[r].elist.size(); e++) {
                v_cur = v_opt = Pt[r][e];// the index of the corresponding edge in ab[r]
                ncyc4_opt = ncyc4;
                ncyc_opt = ncyc;
                n_opt = n_cur;
                N_opt = N_cur;
                sum_opt = sum;
                for (k = 0; k < ab[r][0].size(); k++) {
                    ncyc4_temp = ncyc4;
                    ncyc_temp = ncyc;
                    sum_temp = sum;
                    for (i = 0; i < 3; i++) {
                        vector<int> t1 = H[r].edge_cyclelist[e][2 * i];
                        vector<int> t2 = H[r].edge_cyclelist[e][2 * i + 1];
                        for (j = 0; j < t1.size(); j++) {
                            for (l = 0; l < 2; l++) {
                                v_temp[l] = ab[r][l][k];
                                v_diff[l] = v_temp[l] - ab[r][l][v_cur];
                                s_new[l] = s_old[l] = sum[r][l][i][t1[j]];
                                s_new[l] += ((t2[j] > 0) ? v_diff[l] : (-v_diff[l]));
                                s_new[l] = (s_new[l] >= 0 ? (s_new[l] % L[l]) : (L[l] - 1 - ((-1 - s_new[l]) % L[l])));
                                sum_temp[r][l][i][t1[j]] = s_new[l];
                            }
                            if (i)
                                ncyc_temp[r][i - 1] += (((s_new[0] || s_new[1]) ? 0 : 1) -
                                                        ((s_old[0] || s_old[1]) ? 0 : 1));
                            else {
                                ncyc4_temp[r][s_old[0]][s_old[1]]--;
                                ncyc4_temp[r][s_new[0]][s_new[1]]++;
                            }
                        }
                    }
                    // compute the total number
                    n_temp = num_cycle(ncyc_temp, ncyc4_temp, L, m);
                    N_temp = w * n_temp[1] + n_temp[2];

                    if ((n_temp[0] < n_opt[0]) || ((n_temp[0] <= n_opt[0]) && (N_temp <
                                                                               N_opt))) {// the condition could be changed to ensure d[i] to be 0 if it is in the previous round
                        bool edit = edit_norm(count[r], v_cur, k, th1, th2);
                        //cout<<"edit!\n";
                        if (edit) {
                            v_opt = k;
                            //cout<<"(N_temp,N_opt)=("<<N_temp<<','<<N_opt<<')'<<endl;
                            //cout<<"(n4_temp,n4_opt)=("<<n_temp[0]<<','<<n_opt[0]<<')'<<endl;
                            ncyc4_opt = ncyc4_temp;
                            ncyc_opt = ncyc_temp;
                            n_opt = n_temp;
                            N_opt = N_temp;
                            sum_opt = sum_temp;
                        }
                    }
                }

                if (v_opt != v_cur) {
                    nonstop = true;
                    count[r][v_opt]++;
                    count[r][v_cur]--;
                    Pt[r][e] = v_opt;
                    ncyc = ncyc_opt;
                    ncyc4 = ncyc4_opt;
                    n_cur = n_opt;
                    N_cur = N_opt;
                    sum = sum_opt;
                    /*
                    cout<<"n_cur update:\n";
                    print_int(n_cur);
                    cout<<"N_cur update="<<N_cur<<endl;
                    cout<<"ncyc update:\n";
                    print_matrix(ncyc);
                     */
                    /*
                    cout<<"ncyc4[0]:\n";
                    print_matrix(ncyc4[0]);
                    cout<<"ncyc4[1]:\n";
                    print_matrix(ncyc4[1]);
                     */
                }
            }
        }
    }

    cout << "P_end:" << endl;
    print_matrix(Pt);
    return Pt;
}

#endif //QLDPC_SC_QLDPC_OPTIMIZATION_H

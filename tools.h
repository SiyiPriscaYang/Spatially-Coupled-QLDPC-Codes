//
// Created by Siyi Yang on 1/2/24.
//

#ifndef QLDPC_TOOLS_H
#define QLDPC_TOOLS_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;

double mexp(double x);
double mlog(double x);
double mllr(double x);

void print_int(vector<int>);
void print_double(vector<double>);
void print_matrix(vector<vector<int> >  H);
void print_matrix_double(vector<vector<double> >  H);

vector<int> is_in_row_space(vector<vector<int> > G,vector<int> x);
vector<vector<int> > generator(vector<vector<int> > H);

int vec_norm(vector<int> v);
bool edit_norm(vector<int> v, int ind1,int ind2,int th1,int th2);

int vec_max(vector<int> v);
int matrix_max(vector<vector<int> >  H);

double phi(double x);

vector<int> gen_rec(double p, int len);
vector<int> gen_qrec(vector<double> p, int len);
vector<double> get_lch(double p, int len);
vector<double> get_qlch_d1(vector<double> p);
vector<vector<double> >  get_qlch(vector<double> p, int len);

vector<double> F(vector<double> p);
vector<double> IF(vector<double> P);
int estimate(vector<double> llr);
double q2b(double x,double y);

vector<double> conv1(vector<double> x,vector<double> y);
vector<vector<double> >  conv2(vector<vector<double> >  x,vector<vector<double> >  y);
vector<double> flip1(vector<double> x);
vector<vector<double> >  flip2(vector<vector<double> >  x);

int min3(int x,int y,int z);
int max3(int x,int y,int z);

vector<double> indexing(vector<double> x,int shift,vector<int> ind);
vector<double> add(vector<double> x,vector<double> y,double scalar);
vector<double> mul(double alpha,vector<double> x);
vector<double> normalize(vector<double> x);

double index_sum(vector<vector<double> >  x, int L1,int L2,int s1,int s2);
vector<double> index_sum_grad(vector<vector<double> >  x,int m1,int m2,int L1,int L2,int s1,int s2,vector<vector<int> >  ind);
vector<vector<double> >  new_matrix_dim(int r,int n);
vector<vector<double> >  new_matrix_index(vector<double> x,vector<vector<int> >  a,int m1,int m2);

vector<vector<double> >  get_weight(int ir1,int ir2,int in1,int in2);
vector<vector<int> >  initialize(int m1,int m2);
vector<vector<int> >  p2count(vector<double> p,int len,int mt);

int generateSeed (int i);
int generateSeed (int i) {
    return rand() % i;
}

double mexp(double x) {
    if (x > 27.5) {
        return 5.32e11;
    } else {
        return exp(x);
    }
}

double mlog(double x) {
    if (x < 1e-12) {
        return -27.6;
    } else {
        return log(x);
    }
}

double mllr(double x) {
    double llr;
    if(x< 1e-12)
        llr=-27.6;
    else{
        if(x<(1-1e-12))
            llr=log(x/(1-x));
        else
            llr=27.6;
    }
    return llr;
}


// checked
void print_int(vector<int> x){
    for(int i=0;i<x.size();i++)
        cout<<x[i]<<",";
    cout<<endl;
}

// checked
void print_double(vector<double> x){
    for(int i=0;i<x.size();i++)
        cout<<x[i]<<",";
    cout<<endl;
}

// checked
void print_matrix(vector<vector<int> >  H){
    for(int i=0;i<H.size();i++)
        print_int(H[i]);
}

void print_matrix_double(vector<vector<double> >  H){
    for(int i=0;i<H.size();i++)
        print_double(H[i]);
}


vector<int> is_in_row_space(vector<vector<int> > G,vector<int> x) {
    int m = G.size()-1;
    int n = G[0].size();
    vector<int> ind_sys=G[m];
    vector<int> res = x;

    vector<int> lcomb(0);
    vector<int> t{-1};

    // subtract the rows of G corresponding to indices in ind_sys from res
    for (int i=0;i<ind_sys.size();i++) {
        if (x[ind_sys[i]]) {
            for (int j = 0; j < n; j++) {
                res[j] ^= G[i][j];
            }
            lcomb.push_back(i);
        }
    }

    // check if res is all zero
    for (int i = 0; i < n; i++) {
        if (res[i]) {
            return t;
        }
    }

    return lcomb;
}

vector<vector<int> > generator(vector<vector<int> > H){
    int n = H.size();
    int m = H[0].size();

    vector<vector<int> > G = H;
    vector<int> ind_sys(0);

    for (int j = 0, i = 0; j < m && i < n; j++) {
        int k = i;
        for (int p = i; p < n; p++) {
            if (G[p][j] ) {
                k = p;
                break;
            }
        }

        if (G[k][j] ) {
            // swap rows k and i
            swap(G[k], G[i]);

            // clear entries below the pivot
            for (int p = i + 1; p < n; p++) {
                if (G[p][j] == 1) {
                    for (int q = j; q < m; q++) {
                        G[p][q] = (G[p][q] + G[i][q]) % 2;
                    }
                }
            }

            ind_sys.push_back(j);
            i++;
        }
    }

    // create identity matrix
    for (int i = 0; i < ind_sys.size(); i++) {
        G[i][ind_sys[i]] = 1;
        for (int j = 0; j < i; j++) {
            if (G[j][ind_sys[i]] == 1) {
                for (int k = 0; k < m; k++) {
                    G[j][k] = (G[j][k] + G[i][k]) % 2;
                }
            }
        }
    }

    G.push_back(ind_sys);

    return G;
}


int vec_norm(vector<int> v){
    int n=0;
    for(int i=0;i<v.size();i++)
        n+=((v[i]>0)? v[i]:(-v[i]));
    return n;
}

bool edit_norm(vector<int> v, int ind1,int ind2,int th1,int th2){
    int n1=0;
    int n2=0;
    int nt=0;
    for(int i=0;i<v.size();i++) {
        if((i!=ind1)&&(i!=ind2))
            nt=((v[i]>0)? v[i]:(-v[i]));
        else {
            if(i==ind1)
                nt=((v[i]>1)? (v[i]-1):(1-v[i]));
            else
                nt=((v[i]>-1)? (v[i]+1):(-1-v[i]));
        }
        if(nt>n2)
            n2=nt;
        n1+=nt;
    }
    return (((n1<=th1)&&(n2<=th2))? true:false);
}

// checked
int matrix_max(vector<vector<int> >  H){
    int M=H[0][0];
    //cout<<H.size()<<",";
    //cout<<H[0].size()<<",";
    for(int i=0;i<H.size();i++)
        for(int j=0;j<H[0].size();j++)
            if(H[i][j]>M)
                M=H[i][j];
    return M;
}

int vec_max(vector<int> v){
    int M=v[0];
    for(int i=0;i<v.size();i++)
        if(v[i]>M)
            M=v[i];
    return M;
}

// checked
double phi(double x)
{
    double t=0;
    if(x>27.5)
        t=0;//2*exp(-x);
    else if(x<1e-12)
        t=28;//-log(x/2);
    else
        t=-log(tanh(x/2));
    return t;
}


// checked
vector<int> gen_rec(double p, int len){
    vector<int> x(len);
    double rand_num;
    for(int i=0;i<len;i++){
        rand_num=((double) rand()/(RAND_MAX));
        x[i]=(rand_num>p)? 0:1;
    }
    return x;
}

// checked
vector<double> get_lch(double p, int len){
    vector<double> lch(len);
    double llr;
    if(p==0)
        llr=22;
    else
        llr=log((1-p)/p);
    for(int i=0;i<len;i++)
        lch[i]=llr;//(y[i])? llr:-llr;
    return lch;
}

// checked
vector<int> gen_qrec(vector<double> p, int len){
    int i,j,t;
    double x;
    int plen=p.size();
    vector<int> y(len);
    vector<double> cdf(plen);
    x=1;
    for(i=0;i<plen;i++)
        x-=p[i];
    cdf[0]=x;
    for(i=0;i<(plen-1);i++)
        cdf[i+1]=p[i]+cdf[i];
    double rand_num;
    for(i=0;i<len;i++){
        rand_num=((double) rand()/(RAND_MAX));
        t=0;
        for(j=0;j<plen;j++)
            if(rand_num>cdf[j])
                t++;
        y[i]=t;
    }
    return y;
}

// checked
vector<vector<double> >  get_qlch(vector<double> p, int len){
    int t=p.size();
    int i;
    double x=1-(p[0]+p[1]+p[2]);
    vector<vector<double> >  lch(len);
    vector<double> llr(p.size());
    if(x>(1-3*exp(-22)))
        llr[0]=llr[1]=llr[2]=-22;
    else{
        for(i=0;i<t;i++)
            llr[i]=log(p[i]/x);
    }
    for(i=0;i<len;i++)
        lch[i]=llr;//(y[i])? llr:-llr;
    return lch;
}

vector<double> get_qlch_d1(vector<double> p) {
    int t = p.size();
    int i;
    double x = 1 - (p[0] + p[1] + p[2]);
    vector<double> llr(p.size());
    if (x > (1 - 3 * exp(-22)))
        llr[0] = llr[1] = llr[2] = -22;
    else {
        for (i = 0; i < t; i++)
            llr[i] = log(p[i] / x);
    }
    double lt = llr[1];
    llr[1] = llr[2];
    llr[2] = lt;
    return llr;
}


vector<double> F(vector<double> p)
{
    double a=1.0;
    double b=exp(p[0]);
    double c=exp(p[1]);
    double d=exp(p[2]);
    double s=a+b+c+d;
    a=a/s;
    b=b/s;
    c=c/s;
    d=d/s;
    vector<double> P(3);
    P[0]=a-b+c-d;
    P[1]=a+b-c-d;
    P[2]=a-b-c+d;
    return P;
}

vector<double> IF(vector<double> P)
{
    vector<double> p(3);
    double s=mlog(1.0+P[0]+P[1]+P[2]);
    p[0]=mlog(1.0-P[0]+P[1]-P[2])-s;
    p[1]=mlog(1.0+P[0]-P[1]-P[2])-s;
    p[2]=mlog(1.0-P[0]-P[1]+P[2])-s;
    return p;
}

//checked
int estimate(vector<double> llr){
    int i,t;
    int v=0;
    double m=0.0;
    for(i=0;i<3;i++){
        t=llr[i];
        if(t>m){
            v=i+1;
            m=t;
        }
    }
    return v;
}

double q2b(double x,double y){
    double l=0;
    if(y>x)
        l=((y-x)>27.5? y:(x+log(1+exp(y-x))));
    else
        l=((x-y)>27.5? x:(y+log(1+exp(x-y))));
    return l;
}


vector<double> conv1(vector<double> x,vector<double> y){
    int i,k;
    double tt;
    vector<double> t(x.size()+y.size()-1);
    for(k=0;k<(x.size()+y.size()-1);k++) {
        tt=0;
        for (i = 0; i < x.size(); i++) {
            if (((i +y.size())>k) && (i <= k)) {
                tt += (x[i] * y[k - i]);// 0<=k-i<y.size() k>=i>k-y.size()
            }
        }
        t[k]=tt;
    }
    return t;
}

vector<vector<double> >  conv2(vector<vector<double> >  x,vector<vector<double> >  y){
    int i1,k1,i2,k2;
    int n=x.size()+y.size()-1;
    int m=x[0].size()+y[0].size()-1;
    double ttt;
    vector<vector<double> >  t(n);
    for(i1=0;i1<n;i1++){
        vector<double> tt(m);
        t[i1]=tt;
    }
    for(k1=0;k1<n;k1++)
        for(k2=0;k2<m;k2++) {
            ttt=0;
            for (i1 = 0; i1 < x.size(); i1++)
                if (((i1 + y.size()) > k1) && (i1 <= k1))
                    for (i2 = 0; i2 < x[0].size(); i2++)
                        if (((i2 + y[0].size()) > k2) && (i2 <= k2))
                            ttt += (x[i1][i2] * y[k1 - i1][k2 - i2]);
            t[k1][k2]=ttt;
        }
    return t;
}

vector<double> flip1(vector<double> x){
    int l=x.size();
    vector<double> t(l);
    for(int i=0;i<l;i++)
        t[l-1-i]=x[i];
    return t;
}

vector<vector<double> >  flip2(vector<vector<double> >  x){
    int l=x.size();
    vector<vector<double> >  t(l);
    for(int i=0;i<l;i++){
        vector<double> tt=flip1(x[l-1-i]);
        t[i]=tt;
    }
    return t;
}



vector<double> indexing(vector<double> x,int shift,vector<int> ind){
    vector<double> y(ind.size());
    for(int i=0;i<ind.size();i++)
        y[i]=x[ind[i]+shift];
    return y;
}


vector<double> add(vector<double> x,vector<double> y,double scalar){
    vector<double> z(x.size());
    for(int i=0;i<x.size();i++)
        z[i]=x[i]+scalar*y[i];
    return z;
}

vector<double> mul(double alpha,vector<double> x){
    vector<double> y(x);
    for(int i=0;i<y.size();i++)
        y[i]*=alpha;
    return y;
}

vector<double> normalize(vector<double> x){
    int len=x.size();
    vector<double> y(len);
    double sum=0;
    int i;
    for(i=0;i<len;i++)
        sum+=x[i];
    sum/=((double) len);
    for(i=0;i<len;i++)
        y[i]=x[i]-sum;
    sum=0;
    for(i=0;i<len;i++)
        sum+=pow(y[i],2);
    sum= sqrt(sum);
    for(i=0;i<len;i++)
        y[i]/=sum;
    return y;
}



double index_sum(vector<vector<double> >  x,int L1,int L2,int s1,int s2){
    double sum=x[s1][s2];
    int i,j;
    for(i=L1;i<=s1;i+=L1)
        sum+=x[s1+i][s2]+x[s1-i][s2];
    for(j=L2;j<=s2;j+=L2)
        sum+=x[s1][s2+j]+x[s1][s2-j];
    for(i=L1;i<=s1;i+=L1)
        for(j=L2;j<=s2;j+=L2)
            sum+=x[s1+i][s2+j]+x[s1+i][s2-j]+x[s1-i][s2+j]+x[s1-i][s2-j];
    return sum;
}

vector<double> index_sum_grad(vector<vector<double> >  x,int m1,int m2,int L1,int L2,int s1,int s2,vector<vector<int> >  ind){
    vector<double> y(ind[0].size());
    double sum;
    int i,j,k,i0,j0,r1,r2;
    for(k=0;k<ind[0].size();k++){
        i0=ind[0][k];
        j0=ind[1][k];
        sum=0;
        r1=(s1-i0)%L1;
        r2=(s2-j0)%L2;
        r1=2*s1-r1-m1;
        r2=2*s2-r2-m2;
        for(i=r1;i>=0;i-=L1){
            for(j=r2;j>=0;j-=L2){
                sum+=x[i][j];
            }
        }
        y[k]=sum;
    }

    return y;
}

vector<vector<double> >  new_matrix_dim(int r,int n){
    vector<vector<double> >  t(r);
    for(int i=0;i<r;i++){
        vector<double> ti(n);
        t[i]=ti;
    }
    return t;
}

vector<vector<double> >  new_matrix_index(vector<double> x,vector<vector<int> >  a,int m1,int m2){
    vector<vector<double> >  t=new_matrix_dim(m1+1,m2+1);
    for(int i=0;i<a[0].size();i++)
        t[a[0][i]][a[1][i]]=x[i];
    return t;
}


vector<vector<int> >  p2count(vector<double> p,int len,int mt){
    vector<int> count(mt+1);
    vector<int> countt(mt+1);
    vector<int> ind(mt+1);
    vector<double> order_p(p);
    int i,j,k;
    double min,t,sum,dlen;
    for(i=0;i<=mt;i++)
        ind[i]=i;
    for(i=0;i<mt;i++){
        min=order_p[i];
        k=i;
        for(j=i+1;j<=mt;j++){
            t=order_p[j];
            if(t<min){
                k=j;
                min=t;
            }
        }
        j=ind[i];
        ind[i]=ind[k];
        ind[k]=j;
        order_p[k]=order_p[i];
        order_p[i]=min;
    }

    //print_int(ind);
    //print_double(order_p);
    vector<double> cdf(order_p);
    sum=0;
    dlen=(double) len;
    for(i=0;i<mt;i++){
        sum+=cdf[i];
        cdf[i]=sum*dlen;
    }
    cdf[mt]=dlen;
    //print_double(cdf);

    int isum=len;
    for(i=mt;i>0;i--){
        k=(int) round(cdf[i-1]);
        count[i]=isum-k;
        isum=k;
    }
    count[0]=isum;
    //print_int(count);

    for(i=0;i<=mt;i++)
        countt[ind[i]]=count[i];

    srand ( unsigned ( time(0) ) );

    vector<int> myvector;
    for (i=0; i<=mt; i++) {
        k=i;//a[i];
        for (j = 0; j < countt[i]; j++)
            myvector.push_back(k);
    }

    // shuffle all elements inside the vector
    random_shuffle(myvector.begin(), myvector.end());
    random_shuffle(myvector.begin(), myvector.end(), generateSeed);

    vector<vector<int> >  ini(2);
    ini[0]=countt;
    ini[1]=myvector;
    return ini;
}


vector<vector<double> >  get_weight(int ir1,int ir2,int in1,int in2){
    vector<vector<double> >  w(2);
    vector<double> wa(4);
    vector<double> wb(4);
    double r1=((double) ir1);
    double r2=((double) ir2);
    double n1=((double) in1);
    double n2=((double) in2);

    wa[0]=(2.0/3.0)*(r1-2)*(n1-2);
    wa[1]=(1.0/2.0)*(r1*(r1-3)+3)*(n1*(n1-3)+3);
    wa[2]=4.0*(r1+n1)/(r1*(r1-1)*n1*(n1-1));
    wa[3]=4.0/((r1-1)*(n1-1));
    wb[0]=(2.0/3.0)*(r2-2)*(n2-2);
    wb[1]=(1.0/2.0)*(r2*(r2-3)+3)*(n2*(n2-3)+3);
    wb[2]=4.0*(r2+n2)/(r2*(r2-1)*n2*(n2-1));
    wb[3]=4.0/((r2-1)*(n2-1));

    w[0]=wa;
    w[1]=wb;

    return w;
}

vector<vector<int> >  initialize(int m1,int m2){
    vector<vector<int> >  ind(2);
    vector<int> ind1((m1+1)*(m2+1));
    vector<int> ind2((m1+1)*(m2+1));
    int k=0;
    int i,j;
    for(i=0;i<=m1;i++)
        for(j=0;j<=m2;j++){
            ind1[k]=i;
            ind2[k]=j;
            k++;
        }
    ind[0]=ind1;
    ind[1]=ind2;

    return ind;
}


int min3(int x,int y,int z){
    int m=0;
    m=(y>x)? x:y;
    m=(z>m)? m:z;
    return m;
}

int max3(int x,int y,int z){
    int m=0;
    m=(y<x)? x:y;
    m=(z<m)? m:z;
    return m;
}

#endif //QLDPC_TOOLS_H

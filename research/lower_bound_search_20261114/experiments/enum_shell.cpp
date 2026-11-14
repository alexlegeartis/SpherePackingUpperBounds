#include <bits/stdc++.h>
using namespace std;
int n; double R[64][64],bound,shift[64]; long long z[64]; unsigned long long nodes=0,cap; ofstream out;
void go(int k,double q){
 if(++nodes>cap){cerr<<"CAPPED nodes "<<nodes<<endl;exit(2);}
 if(k<0){if(q>1e-6){for(int i=0;i<n;i++)out<<z[i]<<' ';out<<'\n';}return;}
 double a=-shift[k];for(int j=k+1;j<n;j++)a-=R[k][j]*(z[j]+shift[j])/R[k][k];
 double r=sqrt(max(0.,bound-q))/R[k][k];
 for(long long x=ceil(a-r-1e-8);x<=floor(a+r+1e-8);x++){
  z[k]=x;double v=R[k][k]*(x-a);if(q+v*v<=bound+1e-7)go(k-1,q+v*v);
 }
}
int main(int argc,char**argv){
 if(argc!=7&&argc!=8)return 1;n=atoi(argv[1]);bound=atof(argv[2]);cap=strtoull(argv[3],0,10);
 ifstream f(argv[4]);for(int i=0;i<n;i++)for(int j=0;j<n;j++)f>>R[i][j];if(!f)return 3;
 if(argc==8){ifstream s(argv[7]);for(int i=0;i<n;i++)s>>shift[i];if(!s)return 3;}
 out.open(argv[5]);auto t=chrono::steady_clock::now();go(n-1,0);out.close();
 ofstream log(argv[6]);log<<"FLOATING_COMPLETE nodes "<<nodes<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl;
}

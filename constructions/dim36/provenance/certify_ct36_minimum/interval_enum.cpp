#include <bits/stdc++.h>
using namespace std;
const int N=36; const double INF=numeric_limits<double>::infinity();
inline double dn(double x){return nextafter(x,-INF);} inline double up(double x){return nextafter(x,INF);}
struct I{double l,h;}; I mu[N][N],d[N]; long long Q[N][N],z[N]; unsigned long long nodes=0,leaves=0; bool found=false;chrono::steady_clock::time_point start;
void go(int k,double q){
 if(found)return; ++nodes;
 if(nodes%100000000==0){cerr<<"nodes "<<nodes<<" elapsed "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<" depth "<<k<<endl;}
 if(k<0){++leaves;__int128 s=0;for(int i=0;i<N;i++)for(int j=0;j<N;j++)s+=(__int128)z[i]*Q[i][j]*z[j];if(s>0&&s<160){found=true;cout<<"COUNTEREXAMPLE norm "<<(long long)s<<" vector ";for(auto x:z)cout<<x<<' ';cout<<endl;}return;}
 double cl=0,ch=0;
 for(int j=k+1;j<N;j++)if(z[j]){double a=mu[k][j].l*z[j],b=mu[k][j].h*z[j];cl=dn(cl+dn(min(a,b)));ch=up(ch+up(max(a,b)));}
 // center lies in [-ch,-cl]. Radius uses lower partial norm and lower d.
 double rad=up(sqrt(up(up(159-q)/d[k].l)));
 if(!isfinite(rad)||fabs(cl)>1e6||fabs(ch)>1e6||rad>1e6){cerr<<"ABORT unsafe bound"<<endl;exit(3);}
 long long lo=(long long)ceil(dn(-ch-rad)),hi=(long long)floor(up(-cl+rad));
 if(abs(lo)>1000000||abs(hi)>1000000){cerr<<"ABORT unsafe coordinate"<<endl;exit(3);}
 for(long long x=lo;x<=hi;x++){
  double a=dn(x+cl),b=up(x+ch);double near=(a>0?a:(b<0?-b:0));double add=dn(d[k].l*dn(near*near));double nq=dn(q+max(0.,add));
  if(nq<=159){z[k]=x;go(k-1,nq);if(found)return;}
 }
}
int main(int argc,char**argv){if(argc!=3)return 2;ifstream f(argv[1]),g(argv[2]);long long a,b;for(int i=0;i<N;i++){f>>a>>b;d[i]={ldexp((double)a,-40),ldexp((double)b,-40)};for(int j=0;j<N;j++){f>>a>>b;mu[i][j]={ldexp((double)a,-40),ldexp((double)b,-40)};}}for(auto&r:Q)for(auto&x:r)g>>x;if(!f||!g)return 2;for(auto x:d)if(!(x.l>0)||!isfinite(x.h))return 3;start=chrono::steady_clock::now();go(N-1,0);cout<<"COMPLETE "<<!found<<" nodes "<<nodes<<" leaves "<<leaves<<endl;}

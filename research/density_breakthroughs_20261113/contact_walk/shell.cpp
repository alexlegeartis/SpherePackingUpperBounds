#include <bits/stdc++.h>
using namespace std;
// Discovery oracle only: floating Cholesky pruning is not a minimum proof.
int n; vector<vector<double>> R; vector<vector<long long>> Q;
vector<long long> z; unsigned long long nodes=0, leaves=0, cap;
double bound; bool capped=false, first=false, found=false;
ofstream out;
void go(int k,double q) {
 if(capped||found)return;
 if(++nodes>cap){capped=true;return;}
 if(k<0){
  __int128 s=0;bool nonzero=false,positive=false;
  for(int i=0;i<n;i++){if(!nonzero&&z[i]){nonzero=true;positive=z[i]>0;}for(int j=0;j<n;j++)s+=(__int128)z[i]*Q[i][j]*z[j];}
  if(!nonzero||!positive||s>(long long)llround(bound))return;
  ++leaves;out<<(long long)s;for(auto x:z)out<<' '<<x;out<<'\n';
  if(first)found=true;
  return;
 }
 double a=0;for(int j=k+1;j<n;j++)a-=R[k][j]*z[j]/R[k][k];
 double rem=bound-q;if(rem< -1e-6)return;
 double rad=sqrt(max(0.,rem))/R[k][k];
 long long lo=ceil(a-rad-1e-7),hi=floor(a+rad+1e-7),mid=llround(a);
 if(abs(lo)>1000000||abs(hi)>1000000)throw runtime_error("unsafe range");
 for(long long s=0;s<=max(mid-lo,hi-mid);s++){
  for(int sg : (s?vector<int>{-1,1}:vector<int>{1})){
   long long x=mid+sg*s;if(x<lo||x>hi)continue;z[k]=x;
   double v=R[k][k]*(x-a);if(q+v*v<=bound+1e-6)go(k-1,q+v*v);
   if(capped||found)return;
  }
 }
}
int main(int argc,char**argv){
 if(argc!=8)return 2;
 n=stoi(argv[1]);bound=stod(argv[4]);cap=stoull(argv[5]);first=stoi(argv[7]);
 Q.assign(n,vector<long long>(n));R.assign(n,vector<double>(n));z.assign(n,0);
 ifstream qf(argv[2]),rf(argv[3]);for(auto&r:Q)for(auto&x:r)qf>>x;for(auto&r:R)for(auto&x:r)rf>>x;
 if(!qf||!rf)return 3;out.open(argv[6]);
 auto t=chrono::steady_clock::now();go(n-1,0);
 cout<<"{\"nodes\":"<<nodes<<",\"leaves\":"<<leaves<<",\"capped\":"<<(capped?"true":"false")<<",\"found\":"<<(found?"true":"false")<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<"}"<<endl;
}

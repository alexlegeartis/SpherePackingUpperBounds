#include <bits/stdc++.h>
using namespace std;
// Discovery only: floating point pruning and leaf norms are NOT certificates.
int n; vector<vector<double>> Q,R; vector<long long> z;
uint64_t nodes=0, cap, leaves=0, maxleaves; double bound; bool capped=false; ofstream out;
void rec(int k,double partial){
 if(capped)return;
 if(++nodes>cap){capped=true;return;}
 if(k<0){
  bool nz=false,pos=false;double q=0;
  for(int i=0;i<n;i++){if(!nz&&z[i]){nz=true;pos=z[i]>0;}for(int j=0;j<n;j++)q+=z[i]*Q[i][j]*z[j];}
  if(!nz||!pos||q>bound+1e-7)return;
  out<<setprecision(17)<<q;for(auto x:z)out<<' '<<x;out<<'\n';
  if(++leaves>=maxleaves)capped=true;
  return;
 }
 double c=0;for(int j=k+1;j<n;j++)c-=R[k][j]*z[j]/R[k][k];
 double rem=bound-partial;if(rem< -1e-7)return;
 double rad=sqrt(max(0.,rem))/R[k][k];
 long long lo=ceil(c-rad-1e-7),hi=floor(c+rad+1e-7),mid=llround(c);
 for(long long s=0;s<=max(mid-lo,hi-mid);s++)for(int sg=-1;sg<=1;sg+=2){
  if(s==0&&sg==1)continue;auto x=mid+sg*s;if(x<lo||x>hi)continue;
  z[k]=x;double a=R[k][k]*(x-c);if(partial+a*a<=bound+1e-7)rec(k-1,partial+a*a);
  if(capped)return;
 }
}
int main(int argc,char**argv){
 if(argc!=8)return 2;n=stoi(argv[1]);bound=stod(argv[4]);cap=stoull(argv[5]);maxleaves=stoull(argv[7]);
 Q.assign(n,vector<double>(n));R=Q;z.assign(n,0);ifstream qf(argv[2]),rf(argv[3]);
 for(auto&r:Q)for(auto&x:r)qf>>x;for(auto&r:R)for(auto&x:r)rf>>x;if(!qf||!rf)return 3;
 out.open(argv[6]);auto t=chrono::steady_clock::now();rec(n-1,0);
 cout<<"{\"nodes\":"<<nodes<<",\"leaves\":"<<leaves<<",\"capped\":"<<(capped?"true":"false")<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<"}"<<endl;
}

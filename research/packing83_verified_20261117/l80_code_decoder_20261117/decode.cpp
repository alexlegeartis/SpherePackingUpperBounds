#include <array>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <chrono>
using namespace std;
using V=array<int,8>;
vector<V>D;vector<int> dc;int J[8][8],lookup[256],X[10][8],Z[16][8],S[256][8],maskv[5][5],parity[32];
uint16_t cost[10][256][16],lo[10][256],lb[6][5][256];vector<int>ord[5];int sel[10],part[5],curz;uint64_t nodes=0,leaves=0,cap=100000000;constexpr int cut=640,inf=10000;
int floordiv(int a,int b){int q=a/b,r=a%b;return q-(r<0);}
pair<V,int> nearest(V a){V best;int bc=inf;for(int off:{0,8}){V z,res;int sum=0,c=0;for(int i=0;i<8;i++){z[i]=floordiv(a[i]-off+8,16);res[i]=a[i]-off-16*z[i];sum+=z[i];c+=res[i]*res[i];}if(sum%2){int k=0;for(int i=1;i<8;i++)if(abs(res[i])>abs(res[k]))k=i;c+=256-32*abs(res[k]);z[k]+=(res[k]>=0?1:-1);}if(c<bc){bc=c;for(int i=0;i<8;i++)best[i]=2*z[i]+off/8;}}return {best,bc};}
int carry(V y){int b=0;for(int j=0;j<8;j++){int v=0;for(int i=0;i<8;i++)v+=y[i]*J[i][j];if(v%4)throw runtime_error("nonintegral E8 coordinate");b|=((v/4%2+2)%2)<<j;}return lookup[b];}
V numerator(int i,int s){V a;for(int j=0;j<8;j++)a[j]=X[i][j]-4*S[s][j]-2*Z[curz][j];return a;}
void tables(){for(int i=0;i<10;i++)for(int s=0;s<256;s++){fill(cost[i][s],cost[i][s]+16,inf);auto a=numerator(i,s);auto [v,best]=nearest(a);if(best>256)throw runtime_error("local coverage premise failed");int cv=carry(v);V res;for(int j=0;j<8;j++)res[j]=a[j]-8*v[j];for(size_t k=0;k<D.size();k++){int q=0;for(int j=0;j<8;j++){int e=res[j]-8*D[k][j];q+=e*e;}if(q<=cut)cost[i][s][cv^dc[k]]=min<int>(cost[i][s][cv^dc[k]],q);}lo[i][s]=*min_element(cost[i][s],cost[i][s]+16);}
 for(int i=0;i<5;i++){ord[i].resize(256);for(int s=0;s<256;s++)ord[i][s]=s;sort(ord[i].begin(),ord[i].end(),[&](int a,int b){return lo[i][a]<lo[i][b];});}
 for(int k=0;k<=5;k++)for(int j=0;j<5;j++){int unknown=0;for(int i=k;i<5;i++)unknown|=maskv[i][j];for(int s=0;s<256;s++){int m=inf;for(int t=0;t<256;t++)if((t&~unknown)==(s&~unknown))m=min<int>(m,lo[5+j][t]);lb[k][j][s]=m;}}
}
V witness(int i,int s,int c,int want){V a=numerator(i,s);auto[v,nn]=nearest(a);int cv=carry(v);for(size_t k=0;k<D.size();k++)if((cv^dc[k])==c){V y;int q=0;for(int j=0;j<8;j++){y[j]=v[j]+D[k][j];int e=a[j]-8*y[j];q+=e*e;}if(q==want)return y;}throw runtime_error("missing local witness");}
bool leaf(){leaves++;for(int j=0;j<5;j++)sel[5+j]=part[j];int target=0;for(int h=0;h<4;h++){int w=0;for(int i=0;i<5;i++)w|=((sel[i]>>(4+h))&1)<<i;target|=parity[w]<<h;}
 int dp[11][16],prev[11][16],choice[11][16];for(auto &row:dp)fill(row,row+16,inf);dp[0][0]=0;for(int i=0;i<10;i++)for(int a=0;a<16;a++)for(int c=0;c<16;c++){int q=dp[i][a]+cost[i][sel[i]][c];if(q<dp[i+1][a^c]){dp[i+1][a^c]=q;prev[i+1][a^c]=a;choice[i+1][a^c]=c;}}
 if(dp[10][target]>cut)return false;int cs[10],v=target;for(int i=9;i>=0;i--){cs[i]=choice[i+1][v];v=prev[i+1][v];}
 cout<<"{\"status\":\"EXACT_WITNESS\",\"glue\":"<<curz<<",\"norm_numerator_over128\":"<<dp[10][target]<<",\"states\":[";for(int i=0;i<10;i++)cout<<(i?",":"")<<sel[i];cout<<"],\"carries\":[";for(int i=0;i<10;i++)cout<<(i?",":"")<<cs[i];cout<<"],\"local_points_times2\":[";for(int i=0;i<10;i++){auto y=witness(i,sel[i],cs[i],cost[i][sel[i]][cs[i]]);cout<<(i?",":"")<<"[";for(int j=0;j<8;j++)cout<<(j?",":"")<<y[j];cout<<"]";}cout<<"]}"<<endl;return true;
}
bool dfs(int k,int spent){if(++nodes>cap)throw runtime_error("NODE_CAP_NOT_CERTIFIED");int bound=spent;for(int i=k;i<5;i++)bound+=lo[i][ord[i][0]];for(int j=0;j<5;j++)bound+=lb[k][j][part[j]];if(bound>cut)return false;if(k==5)return leaf();for(int s:ord[k]){int q=spent+lo[k][s];if(q>cut)break;sel[k]=s;for(int j=0;j<5;j++)part[j]^=s&maskv[k][j];if(dfs(k+1,q))return true;for(int j=0;j<5;j++)part[j]^=s&maskv[k][j];}return false;}
int main(int argc,char**argv){try{if(argc<2)throw runtime_error("input");if(argc>2)cap=stoull(argv[2]);ifstream f(argv[1]);int n;f>>n;D.resize(n);dc.resize(n);for(auto&v:D)for(auto&x:v)f>>x;for(auto&x:dc)f>>x;for(auto&r:J)for(auto&x:r)f>>x;for(auto&x:lookup)f>>x;for(auto&r:X)for(auto&x:r)f>>x;for(auto&r:Z)for(auto&x:r)f>>x;for(auto&r:S)for(auto&x:r)f>>x;for(auto&r:maskv)for(auto&x:r)f>>x;for(auto&x:parity)f>>x;if(!f)throw runtime_error("input truncated");for(curz=0;curz<16;curz++){tables();fill(part,part+5,0);if(dfs(0,0))return 0;cerr<<"GLUE_COMPLETE "<<curz<<" nodes="<<nodes<<" leaves="<<leaves<<endl;}cout<<"{\"status\":\"EXHAUSTIVE_NO_VECTOR_PENDING_AUDIT\",\"nodes\":"<<nodes<<",\"leaves\":"<<leaves<<"}"<<endl;}catch(exception&e){cout<<"{\"status\":\"NOT_CERTIFIED\",\"reason\":\""<<e.what()<<"\",\"nodes\":"<<nodes<<"}"<<endl;return 20;}}

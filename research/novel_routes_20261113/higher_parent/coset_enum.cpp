// Floating discovery only: congruence-constrained enumeration with exact leaves.
// An empty traversal is NOT a formal certificate. A node-capped traversal is unresolved.
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <string>
using namespace std;
static constexpr int N=64;
long long Q[N][N],z[N],c[N],witness[N];double U[N][N],dd[N];
uint64_t nodes,cap;bool stopped,found;long long leafnorm;
bool walk(int i,double partial){
 if(++nodes>=cap){stopped=true;return true;}
 if(i<0){__int128 q=0;for(int a=0;a<N;a++)for(int b=0;b<N;b++)q+=(__int128)z[a]*Q[a][b]*z[b];if(q<=60){leafnorm=(long long)q;copy(z,z+N,witness);found=true;return true;}return false;}
 double center=0;for(int j=i+1;j<N;j++)center-=U[i][j]*z[j];
 double rem=60.000001-partial;if(rem<0)return false;
 double radius=sqrt(rem/dd[i]);
 long long lo=(long long)ceil((center-radius-c[i])/3.0-1e-9),hi=(long long)floor((center+radius-c[i])/3.0+1e-9);
 if(lo>hi)return false;
 long long mid=(long long)nearbyint((center-c[i])/3.0);mid=max(lo,min(hi,mid));
 for(long long step=0;step<=max(mid-lo,hi-mid);step++){
  long long a=mid+step;
  if(a<=hi){z[i]=c[i]+3*a;double f=z[i]-center;if(walk(i-1,partial+dd[i]*f*f))return true;}
  if(step){a=mid-step;if(a>=lo){z[i]=c[i]+3*a;double f=z[i]-center;if(walk(i-1,partial+dd[i]*f*f))return true;}}
 }
 return false;
}
int main(int argc,char**argv){
 if(argc<4)return 2;ifstream f(argv[1]);for(int i=0;i<N;i++)for(int j=0;j<N;j++)f>>Q[i][j];if(!f)return 3;
 long double L[N][N]={},D[N]={};for(int i=0;i<N;i++){L[i][i]=1;for(int j=0;j<i;j++){long double v=Q[i][j];for(int k=0;k<j;k++)v-=L[i][k]*D[k]*L[j][k];L[i][j]=v/D[j];}long double d=Q[i][i];for(int k=0;k<i;k++)d-=L[i][k]*L[i][k]*D[k];if(d<=0)return 4;D[i]=d;dd[i]=(double)d;for(int j=0;j<i;j++)U[j][i]=(double)L[i][j];}
 cap=stoull(argv[3]);ifstream in(argv[2]);int id;
 while(in>>id){for(int i=0;i<N;i++)in>>c[i];if(!in)return 5;nodes=0;stopped=false;found=false;leafnorm=-1;auto t=chrono::steady_clock::now();walk(N-1,0);double seconds=chrono::duration<double>(chrono::steady_clock::now()-t).count();cout<<"{\"id\":"<<id<<",\"status\":\""<<(found?"exact_reject":stopped?"node_cap":"empty_float_candidate")<<"\",\"nodes\":"<<nodes<<",\"seconds\":"<<seconds<<",\"norm\":"<<leafnorm<<",\"witness\":[";if(found)for(int i=0;i<N;i++)cout<<(i?",":"")<<witness[i];cout<<"]}"<<endl;
 }
}

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
long long Q[N][N],z[N],c[N],witness[N];double U[N][N],dd[N],zd[N];
uint64_t nodes,cap;bool stopped,found;long long leafnorm;
bool walk(int i,double partial){
 if(++nodes>=cap){stopped=true;return true;}
 if(i<0){__int128 q=0;for(int a=0;a<N;a++)for(int b=0;b<N;b++)q+=(__int128)z[a]*Q[a][b]*z[b];if(q<=60){leafnorm=(long long)q;copy(z,z+N,witness);found=true;return true;}return false;}
 double c0=0,c1=0,c2=0,c3=0;int j=i+1;
 for(;j+3<N;j+=4){c0+=U[i][j]*zd[j];c1+=U[i][j+1]*zd[j+1];c2+=U[i][j+2]*zd[j+2];c3+=U[i][j+3]*zd[j+3];}
 for(;j<N;j++)c0+=U[i][j]*zd[j];double center=-((c0+c1)+(c2+c3));
 double rem=60.000001-partial;if(rem<0)return false;
 long long mid=(long long)nearbyint((center-c[i])/3.0);long long base=c[i]+3*mid;double delta=base-center;
 if(dd[i]*delta*delta>rem)return false;
 z[i]=base;zd[i]=(double)base;if(walk(i-1,partial+dd[i]*delta*delta))return true;
 int sign=delta>=0?-1:1;
 for(long long step=1;;step++){
  double v=delta+3*step*sign,vv=delta-3*step*sign;bool one=false;
  if(dd[i]*v*v<=rem){one=true;z[i]=base+3*step*sign;zd[i]=(double)z[i];if(walk(i-1,partial+dd[i]*v*v))return true;}
  if(dd[i]*vv*vv<=rem){one=true;z[i]=base-3*step*sign;zd[i]=(double)z[i];if(walk(i-1,partial+dd[i]*vv*vv))return true;}
  if(!one)break;
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

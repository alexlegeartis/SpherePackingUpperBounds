#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
struct I{double l,h;};
const double INF=numeric_limits<double>::infinity();
double dn(double x){if(!isfinite(x))throw runtime_error("nonfinite arithmetic");return nextafter(x,-INF);}
double up(double x){if(!isfinite(x))throw runtime_error("nonfinite arithmetic");return nextafter(x,INF);}
I add(I a,I b){return {dn(a.l+b.l),up(a.h+b.h)};}
I mul(I a,I b){double v[4]={a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h};return {dn(*min_element(v,v+4)),up(*max_element(v,v+4))};}
int n;long long threshold;vector<long long> r,z;vector<vector<long long>> G;vector<I>d,x;vector<vector<I>> mu;uint64_t nodes=0,cap;bool fast=false;bool found=false;vector<long long>witness;__int128 witnessq;
string decimal(__int128 q){if(q==0)return "0";bool neg=q<0;if(neg)q=-q;string s;while(q){s.push_back('0'+q%10);q/=10;}if(neg)s.push_back('-');reverse(s.begin(),s.end());return s;}
void leaf(){vector<long long> v(n);for(int i=0;i<n;i++){v[i]=r[i]-2*z[i];if(llabs(v[i])>(1LL<<22))throw runtime_error("exact accumulator range exceeded");}
 __int128 q=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)q+=(__int128)v[i]*G[i][j]*v[j];
 if(q<=threshold){found=true;witness=v;witnessq=q;}
}
void dfs(int i,double partial,double R){if(found)return;if(++nodes>cap)throw runtime_error("node cap");if(i<0){leaf();return;}
 I c; if(fast){double v=double(r[i])/2;for(int j=i+1;j<n;j++)v+=mu[j][i].l*(double(r[j])/2-double(z[j]));c={dn(v-ldexp(1.0,-10)),up(v+ldexp(1.0,-10))};}else{c=x[i];for(int j=i+1;j<n;j++){I e={dn(x[j].l-double(z[j])),up(x[j].h-double(z[j]))};c=add(c,mul(mu[j][i],e));}}
 double rem=up(R-partial);if(rem<0)return;double rad=up(sqrt(max(0.0,up(rem/d[i].l))));
 double l=ceil(dn(c.l-rad)),h=floor(up(c.h+rad));if(l>h)return;if(l<-(1LL<<20)||h>(1LL<<20))throw runtime_error("coordinate range exceeded");
 long long lo=(long long)l,hi=(long long)h,mid=(long long)llround((c.l+c.h)/2);mid=max(lo,min(hi,mid));long long maxoff=max(mid-lo,hi-mid);
 auto visit=[&](long long q){z[i]=q;I e={dn(double(q)-c.h),up(double(q)-c.l)};double sq=0;if(e.l>0)sq=max(0.0,dn(e.l*e.l));else if(e.h<0)sq=max(0.0,dn(e.h*e.h));double contrib=max(0.0,dn(d[i].l*sq));double np=max(0.0,dn(partial+contrib));if(np<=R)dfs(i-1,np,R);};
 visit(mid);for(long long off=1;off<=maxoff&&!found;off++){if(mid-off>=lo)visit(mid-off);if(mid+off<=hi&&!found)visit(mid+off);}
}
int main(int argc,char**argv){auto start=chrono::steady_clock::now();try{if(argc<3)throw runtime_error("usage: verifier input node_cap");ifstream f(argv[1]);f>>n>>threshold;if(n<1||n>110||threshold<0||threshold>1000000)throw runtime_error("invalid dimensions/threshold");bool bits=(argc>4&&string(argv[4])=="bits");auto readI=[&](I&v){if(bits){uint64_t l,h;f>>l>>h;static_assert(sizeof(double)==8);memcpy(&v.l,&l,8);memcpy(&v.h,&h,8);}else f>>v.l>>v.h;};cap=stoull(argv[2]);fast=(argc>3&&string(argv[3])=="fast");if(!numeric_limits<double>::is_iec559||fegetround()!=FE_TONEAREST)throw runtime_error("IEEE nearest arithmetic required");r.resize(n);z.resize(n);x.resize(n);d.resize(n);G.assign(n,vector<long long>(n));mu.assign(n,vector<I>(n));for(auto &v:r)f>>v;for(auto&row:G)for(auto&v:row){f>>v;if(llabs(v)>(1LL<<30))throw runtime_error("Gram range");}for(auto&v:d){readI(v);if(!(v.l>0&&v.h>=v.l))throw runtime_error("bad diagonal");}for(auto&row:mu)for(auto&v:row){readI(v);if(!(isfinite(v.l)&&isfinite(v.h)&&v.l<=v.h))throw runtime_error("bad interval");}if(!f)throw runtime_error("incomplete input");if(fast){for(auto v:r)if(v<0||v>2)throw runtime_error("fast target restriction");for(auto row:mu)for(auto v:row)if(v.l < -16 || v.h > 16 || up(v.h-v.l)>ldexp(1.0,-47))throw runtime_error("fast coefficient restriction");}for(int i=0;i<n;i++){if(llabs(r[i])>1000000)throw runtime_error("target range");x[i]={dn(double(r[i])/2),up(double(r[i])/2)};for(int j=0;j<n;j++)if(G[i][j]!=G[j][i])throw runtime_error("asymmetric Gram");}
 double R=up(double(threshold)/4);dfs(n-1,0,R);double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();cout<<setprecision(17)<<"{\"status\":\""<<(found?"EXACT_WITNESS":"EXHAUSTIVE_NO_VECTOR")<<"\",\"nodes\":"<<nodes<<",\"seconds\":"<<seconds;if(found){cout<<",\"norm_numerator_over4\":"<<decimal(witnessq)<<",\"residual_numerator\":[";for(int i=0;i<n;i++){if(i)cout<<',';cout<<witness[i];}cout<<']';}cout<<"}\n";return found?10:0;
 }catch(const exception&e){cout<<"{\"status\":\"NOT_CERTIFIED\",\"reason\":\""<<e.what()<<"\",\"nodes\":"<<nodes<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";return 20;}}

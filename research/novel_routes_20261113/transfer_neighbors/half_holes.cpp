#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
// Discovery only: floating Cholesky pruning, exact integer leaf norms.
// Every negative result is heuristic unless a short witness was returned.
using V=std::array<int,32>;
int G[32][32]; double R[32][32]; V z,bestz; uint32_t mask;
long long nodes,totalnodes; long long cap=2000000;
int best; bool capped;
bool collect=false;std::vector<V> contacts;
long long exactnorm(const V&x){long long s=0;for(int i=0;i<32;i++)for(int j=0;j<32;j++)s+=(long long)x[i]*G[i][j]*x[j];return s;}
void rec(int k,double used){
 if(++nodes>cap){capped=true;return;}
 if(k<0){auto q=exactnorm(z);if(collect){if(q<best)contacts.push_back(z);}else if(q<best){best=q;bestz=z;}return;}
 double tail=0;for(int j=k+1;j<32;j++)tail+=R[k][j]*z[j];
 double center=-tail/R[k][k];int parity=(mask>>k)&1;
 double avail=best-1e-7-used;if(avail<0)return;
 double rad=std::sqrt(avail)/R[k][k];
 int lo=(int)std::ceil((center-rad-parity)/2),hi=(int)std::floor((center+rad-parity)/2);
 int mid=(int)std::llround((center-parity)/2);
 for(int step=0;step<=std::max(mid-lo,hi-mid);step++)for(int side=0;side<(step==0?1:2);side++){
  int t=mid+(side?-step:step);if(t<lo||t>hi)continue;
  z[k]=2*t+parity;double val=R[k][k]*z[k]+tail;
  if(used+val*val<best-1e-7)rec(k-1,used+val*val);
  if(capped)return;
 }
}
int cvp(uint32_t m){mask=m;nodes=0;capped=false;
 for(int k=31;k>=0;k--){double t=0;for(int j=k+1;j<32;j++)t+=R[k][j]*z[j];int p=(m>>k)&1;z[k]=p+2*(int)std::llround((-t/R[k][k]-p)/2);}
 best=exactnorm(z);bestz=z;rec(31,0);totalnodes+=nodes;return best;
}
int main(int argc,char**argv){
 if(argc<5)return 1;std::ifstream f(argv[1]);for(auto&r:G)for(auto &v:r)f>>v;
 std::ifstream r(argv[2]);for(auto&a:R)for(auto&v:a)r>>v;
 double seconds=std::stod(argv[3]);std::mt19937 rng(std::stoul(argv[4]));
 auto start=std::chrono::steady_clock::now();long long calls=0,caps=0;int hist[200]={0};int record=0;uint32_t bm=0,current=rng();int cq=0;V saved;
 while(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<seconds){
  uint32_t cand=(calls%64==0)?rng():current^(1u<<(rng()%32));
  int q=cvp(cand);calls++;caps+=capped;hist[std::min(q,199)]++;
  if(q>record){record=q;bm=cand;saved=bestz;std::cerr<<"record "<<record<<" mask "<<bm<<" calls "<<calls<<" cap "<<capped<<"\n";}
  if(calls%64==0||q>=cq||rng()%50==0){current=cand;cq=q;}
  if(q>16&&!capped){std::cerr<<"CANDIDATE min "<<q<<"\n";break;}
 }
 std::cout<<"{\"calls\":"<<calls<<",\"nodes\":"<<totalnodes<<",\"caps\":"<<caps<<",\"best_found_norm\":"<<record<<",\"parity\":"<<bm<<",\"witness\":[";
 for(int i=0;i<32;i++)std::cout<<(i?",":"")<<saved[i];std::cout<<"],\"hist\":{";
 bool first=true;for(int i=0;i<200;i++)if(hist[i]){std::cout<<(first?"":",")<<"\""<<i<<"\":"<<hist[i];first=false;}
 std::cout<<"}}\n";
 if(argc>5){collect=true;mask=bm;best=record+1;nodes=0;capped=false;rec(31,0);std::ofstream out(argv[5]);for(auto&w:contacts){for(int i=0;i<32;i++)out<<(i?" ":"")<<w[i];out<<"\n";}std::cerr<<"contacts "<<contacts.size()<<" cap "<<capped<<"\n";}
}

// Heuristic search for ternary parity columns eliminating all bad outer supports.
// All support tests and final finite checks are exact; search nonexistence is NOT.
// Usage: search_parity supports.txt rank seconds seed result.json
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>
using Clock=std::chrono::steady_clock;
struct Hash{size_t operator()(uint64_t x)const{x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}};
int r,q,group,cols[13][3],plusTable[13][27],minusTable[13][27];
std::unordered_set<uint64_t,Hash> good;
int add(int a,int c,int sign){int ret=0,m=1;for(int j=0;j<r;j++){ret+=((a%3+sign*cols[c][j]+6)%3)*m;a/=3;m*=3;}return ret;}
uint32_t step(uint32_t bits,int c){uint32_t out=0;while(bits){int i=__builtin_ctz(bits);bits&=bits-1;out|=1U<<plusTable[c][i];out|=1U<<minusTable[c][i];}return out;}
void patterns(int left,int first,uint64_t key,uint32_t sums){
 if(!left){if(!(sums&1))good.insert(key);return;}
 for(int c=first;c<q;c++)patterns(left-1,c,key+(1ULL<<(4*c)),step(sums,c));
}
int main(int argc,char**argv){
 if(argc!=6)return 2;
 std::ifstream in(argv[1]);int n,m;in>>n>>m;r=std::stoi(argv[2]);double seconds=std::stod(argv[3]);uint64_t seed=std::stoull(argv[4]);
 if(r<2||r>3)return 3;
 group=1;for(int j=0;j<r;j++)group*=3;q=0;
 for(int x=1;x<group;x++){int t=x,first=0,d[3]={};for(int j=0;j<r;j++){d[j]=t%3;t/=3;if(!first&&d[j])first=d[j];}if(first!=1)continue;for(int j=0;j<r;j++)cols[q][j]=d[j];q++;}
 for(int c=0;c<q;c++)for(int a=0;a<group;a++){plusTable[c][a]=add(a,c,1);minusTable[c][a]=add(a,c,-1);}
 good.reserve(100000);patterns(6,0,0,1);patterns(7,0,0,1);
 std::vector<std::vector<int>> supports(m),incident(n);
 for(int a=0;a<m;a++){uint32_t s;in>>s;for(int j=0;j<n;j++)if(s&(1U<<j)){supports[a].push_back(j);incident[j].push_back(a);}}
 std::vector<int> h(n),besth(n);std::vector<uint64_t> keys(m);std::mt19937_64 rng(seed);
 std::uniform_real_distribution<double> unit(0,1);int best=m+1;uint64_t calls=0;auto start=Clock::now();
 auto elapsed=[&](){return std::chrono::duration<double>(Clock::now()-start).count();};
 auto valid=[&](uint64_t k){return good.count(k)!=0;};
 auto save=[&](int restart){
  std::ofstream f(argv[5]);f<<"{\"length\":"<<n<<",\"rank_limit\":"<<r<<",\"seed\":"<<seed<<",\"bad_supports\":"<<m<<",\"uncovered_supports\":"<<best<<",\"moves\":"<<calls<<",\"seconds\":"<<elapsed()<<",\"restart\":"<<restart<<",\"columns\":[";
  for(int j=0;j<n;j++){if(j)f<<",";f<<"[";for(int k=0;k<r;k++)f<<(k?",":"")<<cols[besth[j]][k];f<<"]";}f<<"],\"scope\":\"heuristic search, exact finite support objective; failure is not a nonexistence proof\"}\n";
  std::cout<<"best "<<best<<"/"<<m<<" moves "<<calls<<" time "<<elapsed()<<std::endl;
 };
 for(int restart=0;elapsed()<seconds;restart++){
  if(restart%4==1&&best<=m){h=besth;for(int t=0;t<1+restart%7;t++)h[rng()%n]=rng()%q;}
  else for(int j=0;j<n;j++)h[j]=(restart%3==0?int(rng()%r):int(rng()%q));
  int cost=0;for(int a=0;a<m;a++){keys[a]=0;for(int j:supports[a])keys[a]+=1ULL<<(4*h[j]);cost+=!valid(keys[a]);}
  for(int t=0;t<6000&&elapsed()<seconds;t++){
   if(cost<best){best=cost;besth=h;save(restart);if(!best)return 0;}
   int j=rng()%n,c=rng()%q;if(c==h[j])continue;
   uint64_t delta=(1ULL<<(4*c))-(1ULL<<(4*h[j]));int diff=0;
   for(int a:incident[j])diff+=int(valid(keys[a]))-int(valid(keys[a]+delta));
   double temp=std::max(0.15,0.015*m*(1.0-double(t)/6000));
   if(diff<=0||unit(rng)<std::exp(-double(diff)/temp)){cost+=diff;h[j]=c;for(int a:incident[j])keys[a]+=delta;}
   calls++;
  }
 }
 save(-1);return 0;
}

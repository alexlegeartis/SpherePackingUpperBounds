#include <array>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <stdexcept>
using Row=std::array<int,8>;
constexpr int MAX=1000000,CUT=640;
std::vector<Row> shell; std::vector<int> shellcarry;
int inverse8[8][8],labels[256],target[10][8],glues[16][8],representative[256][8],masks[5][5],twist[32];
int costs[10][256][16],cheap[10][256],relaxation[6][5][256],chosen[10],partial[5],order[5]={4,3,2,1,0};
std::array<std::vector<int>,5> ordered;
uint64_t nodes=0,leaves=0;int glue;
int floor16(int a){int q=a/16;return q-((a%16)<0);}
int carry(const Row&r){int mask=0;for(int j=0;j<8;j++){int t=0;for(int i=0;i<8;i++)t+=r[i]*inverse8[i][j];if(t%4)throw std::runtime_error("bad E8 coordinate");mask|=((t/4%2+2)%2)<<j;}return labels[mask];}
Row closest(const Row&a){
 int best=MAX;Row answer;
 // Exhaust all lower/upper coordinate choices in each E8 branch, unlike a least-cost-flip routine.
 for(int offset:{0,8}){Row lower;for(int j=0;j<8;j++)lower[j]=floor16(a[j]-offset);
  for(int word=0;word<256;word++){int parity=0,cost=0;Row r;for(int j=0;j<8;j++){int v=lower[j]+((word>>j)&1);parity+=v;r[j]=2*v+offset/8;int e=a[j]-8*r[j];cost+=e*e;}
   if(parity%2==0 && cost<best){best=cost;answer=r;}
  }
 }
 if(best>256)throw std::runtime_error("nearest distance exceeds covering premise");return answer;
}
void build(){
 for(int i=0;i<10;i++)for(int word=0;word<256;word++){
  Row a;for(int j=0;j<8;j++)a[j]=target[i][j]-4*representative[word][j]-2*glues[glue][j];
  Row center=closest(a);int cl=carry(center);std::fill(costs[i][word],costs[i][word]+16,MAX);
  for(size_t k=0;k<shell.size();k++){int norm=0;for(int j=0;j<8;j++){int e=a[j]-8*(center[j]+shell[k][j]);norm+=e*e;if(norm>CUT)break;}
   if(norm<=CUT){int c=cl^shellcarry[k];costs[i][word][c]=std::min(costs[i][word][c],norm);}
  }
  cheap[i][word]=*std::min_element(costs[i][word],costs[i][word]+16);
 }
 for(int k=0;k<5;k++){ordered[k].resize(256);for(int s=0;s<256;s++)ordered[k][s]=s;int row=order[k];std::sort(ordered[k].begin(),ordered[k].end(),[row](int a,int b){return cheap[row][a]!=cheap[row][b]?cheap[row][a]<cheap[row][b]:a>b;});}
 for(int depth=0;depth<=5;depth++)for(int j=0;j<5;j++){
  int freebits=0;for(int k=depth;k<5;k++)freebits|=masks[order[k]][j];
  for(int s=0;s<256;s++){int lower=MAX;for(int t=0;t<256;t++)if(((t^s)&(~freebits&255))==0)lower=std::min(lower,cheap[j+5][t]);relaxation[depth][j][s]=lower;}
 }
}
bool complete(){
 leaves++;for(int j=0;j<5;j++)chosen[j+5]=partial[j];int wanted=0;
 for(int h=0;h<4;h++){int word=0;for(int i=0;i<5;i++)word|=((chosen[i]>>(h+4))&1)<<i;wanted|=twist[word]<<h;}
 std::array<int,16> old,next;old.fill(MAX);old[0]=0;
 for(int row=0;row<10;row++){next.fill(MAX);for(int out=0;out<16;out++)for(int c=0;c<16;c++)next[out]=std::min(next[out],old[out^c]+costs[row][chosen[row]][c]);old=next;}
 if(old[wanted]<=CUT){std::cout<<"COUNTEREXAMPLE glue="<<glue<<" normnumerator="<<old[wanted]<<std::endl;return true;}return false;
}
bool walk(int depth,int used){
 if(++nodes>20000000000ULL)throw std::runtime_error("node limit, incomplete");
 int bound=used;for(int k=depth;k<5;k++)bound+=cheap[order[k]][ordered[k][0]];for(int j=0;j<5;j++)bound+=relaxation[depth][j][partial[j]];
 if(bound>CUT)return false;if(depth==5)return complete();int row=order[depth];
 for(int s:ordered[depth]){int next=used+cheap[row][s];if(next>CUT)break;chosen[row]=s;for(int j=0;j<5;j++)partial[j]^=s&masks[row][j];bool found=walk(depth+1,next);for(int j=0;j<5;j++)partial[j]^=s&masks[row][j];if(found)return true;}return false;
}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("input required");std::ifstream in(argv[1]);int size;in>>size;if(size!=9121)throw std::runtime_error("shell size");shell.resize(size);shellcarry.resize(size);
 for(auto&r:shell)for(int&v:r){in>>v;if(v<-4||v>4)throw std::runtime_error("shell coordinate range");}
 for(int&v:shellcarry){in>>v;if(v<0||v>=16)throw std::runtime_error("carry range");}
 for(auto&r:inverse8)for(int&v:r)in>>v;for(int&v:labels)in>>v;for(auto&r:target)for(int&v:r)in>>v;for(auto&r:glues)for(int&v:r)in>>v;for(auto&r:representative)for(int&v:r)in>>v;for(auto&r:masks)for(int&v:r)in>>v;for(int&v:twist)in>>v;if(!in)throw std::runtime_error("truncated input");int extra;if(in>>extra)throw std::runtime_error("extra input");
 for(auto&r:target)for(int v:r)if(v<-10000||v>10000)throw std::runtime_error("target range");for(auto&r:glues)for(int v:r)if(v<-1000||v>1000)throw std::runtime_error("glue range");for(auto&r:representative)for(int v:r)if(v<-1000||v>1000)throw std::runtime_error("representative range");
 for(glue=0;glue<16;glue++){build();std::fill(partial,partial+5,0);if(walk(0,0))return 2;std::cout<<"COMPLETE_GLUE "<<glue<<" nodes "<<nodes<<" leaves "<<leaves<<std::endl;}
 std::cout<<"EXHAUSTIVE_NO_VECTOR independent_reverse_rows nodes="<<nodes<<" leaves="<<leaves<<std::endl;
 }catch(const std::exception&e){std::cerr<<"INCOMPLETE "<<e.what()<<std::endl;return 1;}return 0;}

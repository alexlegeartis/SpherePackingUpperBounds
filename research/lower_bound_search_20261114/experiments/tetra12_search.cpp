#include <bits/stdc++.h>
using namespace std;
int n,target; vector<uint32_t> blocks; vector<vector<int>> nexts; vector<unsigned char> bad6,bad9;
unsigned long long edges=0,nodes=0,triple_checks=0; vector<int> answer;
array<uint32_t,4> tris(uint32_t a){array<uint32_t,4>r;uint32_t b=a;int i=0;while(b){uint32_t x=b&-b;b-=x;r[i++]=a^x;}return r;}
bool pairgood(uint32_t a,uint32_t b){if(a&b)return false;auto x=tris(a),y=tris(b);for(auto u:x)for(auto v:y)if(bad6[u|v])return false;return true;}
bool triplegood(uint32_t a,uint32_t b,uint32_t c){++triple_checks;auto x=tris(a),y=tris(b),z=tris(c);for(auto u:x)for(auto v:y)for(auto w:z)if(bad9[u|v|w])return false;return true;}
bool go(vector<int>&chosen,const vector<int>&candidates,uint32_t used){
 ++nodes;if((int)chosen.size()==target){answer=chosen;return true;}if((int)chosen.size()+(int)candidates.size()<target)return false;
 for(int k=0;k<(int)candidates.size();++k){int a=candidates[k];if((int)chosen.size()+(int)candidates.size()-k<target)break;bool okay=true;
  for(int i=0;i<(int)chosen.size()&&okay;++i)for(int j=i+1;j<(int)chosen.size();++j)if(!triplegood(blocks[a],blocks[chosen[i]],blocks[chosen[j]])){okay=false;break;}
  if(!okay)continue;
  vector<int> rest;set_intersection(candidates.begin()+k+1,candidates.end(),nexts[a].begin(),nexts[a].end(),back_inserter(rest));chosen.push_back(a);if(go(chosen,rest,used|blocks[a]))return true;chosen.pop_back();
 }return false;
}
int main(int argc,char**argv){if(argc!=3)return 1;ifstream f(argv[1]);int s6,s7,s9;f>>n>>target>>s6>>s7>>s9;if(!f||n>20)return 2;bad6.assign(1u<<n,0);bad9=bad6;
 for(int i=0;i<s6;++i){uint32_t s;f>>s;bad6[s]=1;}
 for(int i=0;i<s7;++i){uint32_t s;f>>s;uint32_t t=s;while(t){uint32_t b=t&-t;t-=b;bad6[s^b]=1;}}
 for(int i=0;i<s9;++i){uint32_t s;f>>s;bad9[s]=1;}
 if(!f)return 3;
 for(uint32_t a=0;a<(1u<<n);++a)if(__builtin_popcount(a)==4)blocks.push_back(a);nexts.resize(blocks.size());
 for(int i=0;i<(int)blocks.size();++i)for(int j=i+1;j<(int)blocks.size();++j)if(pairgood(blocks[i],blocks[j])){nexts[i].push_back(j);++edges;}
 vector<int> candidates(blocks.size());iota(candidates.begin(),candidates.end(),0);vector<int> chosen;bool sat=go(chosen,candidates,0);ofstream out(argv[2]);out<<(sat?"SAT":"UNSAT")<<" n "<<n<<" tetracodes "<<target<<" blocks "<<blocks.size()<<" edges "<<edges<<" nodes "<<nodes<<" triple_checks "<<triple_checks<<"\n";if(sat)for(int i:answer)out<<blocks[i]<<"\n";cout<<(sat?"SAT":"UNSAT")<<" edges "<<edges<<" nodes "<<nodes<<"\n";return 0;
}

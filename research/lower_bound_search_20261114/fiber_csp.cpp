// Exhaustive search for a common rank-at-most-three ternary fiber check.
// Projective columns identify h and -h, which have identical signed sums.
#include <bits/stdc++.h>
using namespace std;
constexpr uint32_t ALL=(1u<<27)-1; constexpr unsigned COLORS=(1u<<13)-1;
int n,m; vector<uint32_t> supports,reach; vector<vector<int>> incident;
vector<array<int,3>> points; unsigned dom[32];int choice[32];
unsigned long long nodes=0,cap; bool capped=false;
uint32_t digitlow[3],digithigh[3];
uint32_t translate(uint32_t a,int color,int sign){
#ifdef INDEPENDENT_STEPS
 uint32_t result=0;
 for(int x=0;x<27;x++)if(a>>x&1){
  int y=0;for(int d=0,k=1;d<3;d++,k*=3)y+=((x/k%3+sign*points[color][d]+3)%3)*k;
  result|=1u<<y;
 }return result;
#else
 for(int d=0,step=1;d<3;d++,step*=3){
  int t=(points[color][d]*sign+3)%3;
  for(int j=0;j<t;j++) a=((a&digitlow[d])<<step)|((a&digithigh[d])>>(2*step));
 }return a;
#endif
}
uint32_t extend(uint32_t a,int color){return translate(a,color,1)|translate(a,color,-1);}
bool dfs(uint32_t left){
 if(++nodes>cap){capped=true;return false;}
 if(!left)return true;
 int v=-1,best=99;
 for(int j=0;j<n;j++)if(left>>j&1){int c=__builtin_popcount(dom[j]);if(c<best){best=c;v=j;}}
 if(best==0)return false;
 for(unsigned avail=dom[v];avail;avail&=avail-1){
  int col=__builtin_ctz(avail);uint32_t rest=left^(1u<<v);bool ok=true;
  vector<pair<int,uint32_t>> oldreach;vector<pair<int,unsigned>> olddom;
  for(int j:incident[v]){
   oldreach.emplace_back(j,reach[j]);reach[j]=extend(reach[j],col);
   uint32_t remain=supports[j]&rest;
   if(!remain){if(reach[j]&1){ok=false;break;}}
   else if(reach[j]==ALL){ok=false;break;}
   else if(!(remain&(remain-1))){
    int u=__builtin_ctz(remain);unsigned allowed=0;
    for(int c=0;c<13;c++)if(dom[u]>>c&1){if(!(extend(reach[j],c)&1))allowed|=1u<<c;}
    if(allowed!=dom[u]){olddom.emplace_back(u,dom[u]);dom[u]=allowed;}
    if(!allowed){ok=false;break;}
   }
  }
  choice[v]=col;
  if(ok && dfs(rest))return true;
  for(auto it=olddom.rbegin();it!=olddom.rend();it++)dom[it->first]=it->second;
  for(auto [j,a]:oldreach)reach[j]=a;
  choice[v]=-1;
  if(capped)return false;
 }
 return false;
}
int main(int argc,char**argv){
 if(argc!=3&&argc!=4)return 1;ifstream f(argv[1]);f>>n>>m;if(n>30||n<2||m<=0)return 1;
 cap=strtoull(argv[2],0,10);supports.resize(m);reach.assign(m,1);incident.resize(n);
 for(int j=0;j<m;j++){f>>supports[j];for(int v=0;v<n;v++)if(supports[j]>>v&1)incident[v].push_back(j);}
 if(!f)return 1;
 points={array<int,3>{1,0,0},array<int,3>{0,1,0},array<int,3>{0,0,1}};
 for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int c=0;c<3;c++){
  array<int,3> v{a,b,c};int first=a?a:(b?b:c);if(first!=1)continue;
  if(find(points.begin(),points.end(),v)==points.end())points.push_back(v);
 }
 if(points.size()!=13)return 2;
 for(int x=0;x<27;x++)for(int d=0,k=1;d<3;d++,k*=3){
  if(x/k%3==2)digithigh[d]|=1u<<x;else digitlow[d]|=1u<<x;
 }
 // Self-check each bit permutation against independent coordinate arithmetic.
 for(int x=0;x<27;x++)for(int c=0;c<13;c++)for(int sg:{-1,1}){
  int y=0;for(int d=0,k=1;d<3;d++,k*=3)y+=((x/k%3+sg*points[c][d]+3)%3)*k;
  if(translate(1u<<x,c,sg)!=(1u<<y))return 3;
 }
 fill(dom,dom+32,COLORS);fill(choice,choice+32,-1);
 // GL(3,3), and independent sign changes, normalize columns 0 and 1.
 if(argc!=4){dom[0]=1;dom[1]=3;}
 auto t=chrono::steady_clock::now();bool yes=dfs((1u<<n)-1);
 cout<<(yes?"SAT":(capped?"CAPPED":"UNSAT"))<<" nodes "<<nodes<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl;
 if(yes)for(int i=0;i<n;i++)cout<<points[choice[i]][0]<<' '<<points[choice[i]][1]<<' '<<points[choice[i]][2]<<'\n';
}

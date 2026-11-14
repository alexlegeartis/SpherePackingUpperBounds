#include <bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){
 if(argc!=5)return 1;ifstream f(argv[1]);int n,k,s6,s7,s9;f>>n>>k>>s6>>s7>>s9;if(!f||n>20||k>10)return 2;
 vector<unsigned char> cost(1u<<n,0);
 for(int i=0;i<s6;++i){unsigned s;f>>s;cost[s]=10;}
 for(int i=0;i<s7;++i){unsigned s;f>>s;for(unsigned t=s;t;){unsigned b=t&-t;t-=b;cost[s^b]=max<int>(cost[s^b],1);}}
 for(int i=0;i<s9;++i){unsigned s;f>>s;cost[s]=max<int>(cost[s],1);}
 vector<vector<int>> G(k,vector<int>(n));for(auto &r:G)for(auto &a:r)f>>a;if(!f)return 3;
 vector<vector<unsigned char>> X;
 int total=1;for(int i=0;i<k;++i)total*=3;
 for(int q=1;q<total;++q){int a=q;vector<int> coeff(k);int highest=0;for(int i=0;i<k;++i){coeff[i]=a%3;if(coeff[i])highest=coeff[i];a/=3;}if(highest!=1)continue;vector<unsigned char> row(n);for(int i=0;i<n;++i){int v=0;for(int j=0;j<k;++j)v+=coeff[j]*G[j][i];row[i]=v%3;}X.push_back(row);}
 vector<unsigned> masks(X.size());int energy=0;for(int j=0;j<(int)X.size();++j){for(int i=0;i<n;++i)if(X[j][i])masks[j]|=1u<<i;energy+=cost[masks[j]];}
 vector<unsigned> proposal(masks.size());vector<array<unsigned char,4>> values(masks.size());
 double seconds=atof(argv[3]);uint64_t cap=strtoull(argv[4],0,10);mt19937_64 rng(20261114+n);uniform_real_distribution<double> unif(0,1);auto start=chrono::steady_clock::now();int best=energy;vector<vector<int>> bestG=G;unsigned long long accepted=0;
 auto save=[&](uint64_t step){ofstream out(argv[2]);out<<n<<" "<<k<<" "<<best<<" "<<step<<"\n";for(auto&r:bestG){for(int i=0;i<n;++i){if(i)out<<" ";out<<r[i];}out<<"\n";}cout<<"best "<<best<<" step "<<step<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;};save(0);
 for(uint64_t step=1;step<=cap;++step){if(step%128==0&&chrono::duration<double>(chrono::steady_clock::now()-start).count()>seconds)break;
  array<int,4> pos,sign;unsigned used=0;for(int i=0;i<4;++i){do{pos[i]=rng()%n;}while(used>>pos[i]&1);used|=1u<<pos[i];sign[i]=1+rng()%2;}
  int change=0;for(int j=0;j<(int)X.size();++j){int dot=0;for(int i=0;i<4;++i)dot+=X[j][pos[i]]*sign[i];dot%=3;unsigned m=masks[j]&~used;
   for(int i=0;i<4;++i){int v=(X[j][pos[i]]+dot*sign[i])%3;values[j][i]=v;if(v)m|=1u<<pos[i];}proposal[j]=m;change+=(int)cost[m]-(int)cost[masks[j]];
  }
  double phase=(step%24000)/24000.;double temperature=step<250?1e9:8*pow(0.025,phase);
  if(change<=0||unif(rng)<exp(-change/temperature)){
   energy+=change;++accepted;for(int j=0;j<(int)X.size();++j)for(int i=0;i<4;++i)X[j][pos[i]]=values[j][i];masks.swap(proposal);
   for(auto&r:G){int dot=0;for(int i=0;i<4;++i)dot+=r[pos[i]]*sign[i];dot%=3;for(int i=0;i<4;++i)r[pos[i]]=(r[pos[i]]+dot*sign[i])%3;}
   if(energy<best){best=energy;bestG=G;save(step);if(!best)break;}
  }
 }
 cout<<"FINISHED best "<<best<<" accepted "<<accepted<<endl;return 0;
}

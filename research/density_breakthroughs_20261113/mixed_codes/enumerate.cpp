// Exact exhaustive GF(4) codeword enumeration using its binary Gray order.
// Input: length n, binary rank b, b packed generator rows (two bits/coordinate).
// Output: all distinct supports of weight six or seven, plus weight counts.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
int main(int argc,char**argv){
 if(argc!=3) return 2;
 std::ifstream in(argv[1]);int n,b;in>>n>>b;
 if(n>31||b>32) return 3;
 std::vector<uint64_t> g(b),hist(n+1);for(auto&v:g)in>>v;
 std::set<uint32_t> bad;uint64_t word=0;
 uint64_t total=uint64_t(1)<<b;
 for(uint64_t i=0;i<total;i++){
   if(i) word^=g[__builtin_ctzll(i)];
   uint64_t bits=(word|(word>>1))&0x5555555555555555ULL;
   int w=__builtin_popcountll(bits);hist[w]++;
   if(w==6||w==7){uint32_t s=0;for(int j=0;j<n;j++)if(bits&(uint64_t(1)<<(2*j)))s|=1U<<j;bad.insert(s);}
 }
 if(hist[0]!=1) return 4;
 for(int w=1;w<6;w++)if(hist[w]) return 5;
 std::ofstream out(argv[2]);out<<n<<" "<<bad.size()<<"\n";for(auto s:bad)out<<s<<"\n";
 std::cout<<"{\"length\":"<<n<<",\"binary_dimension\":"<<b<<",\"total\":"<<total<<",\"bad_supports\":"<<bad.size()<<",\"weight_enumerator\":[";
 for(int i=0;i<=n;i++)std::cout<<(i?",":"")<<hist[i];std::cout<<"]}\n";
}

#include <cmath>
#include <cfenv>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <stdexcept>
#include <chrono>
using I=long long; using F=long double;
int n;I bound;std::vector<I>r,z;std::vector<std::vector<I>>G;std::vector<F>d;std::vector<std::vector<F>>m;
uint64_t nodes=0,leaves=0;constexpr uint64_t limit=6000000000ULL;
F below(F a){if(!std::isfinite(a))throw std::runtime_error("nonfinite");return std::nextafter(a,-std::numeric_limits<F>::infinity());}
F above(F a){if(!std::isfinite(a))throw std::runtime_error("nonfinite");return std::nextafter(a,std::numeric_limits<F>::infinity());}
bool visit(int row,F spent){
 if(++nodes>limit)throw std::runtime_error("unfinished node limit");
 if(row<0){++leaves;__int128 val=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++){I u=r[i]-2*z[i],v=r[j]-2*z[j];if(std::abs(u)>(1LL<<22)||std::abs(v)>(1LL<<22))throw std::runtime_error("leaf overflow range");val+=(__int128)u*G[i][j]*v;}if(val<=bound){std::cout<<"FOUND ";for(auto a:z)std::cout<<a<<' ';std::cout<<'\n';return true;}return false;}
 F center=F(r[row])/2;for(int j=row+1;j<n;j++)center+=m[j][row]*(F(r[j])/2-F(z[j]));
 // Input coefficient error<2^-47, residual<2^21, n<=110: total<2^-19.
 // 64-bit significand long-double rounding contributes<2^-24. 2^-16 is conservative.
 F clo=below(center-std::ldexp(F(1),-16)),chi=above(center+std::ldexp(F(1),-16));
 F remain=above(F(bound)/4-spent);if(remain<0)return false;
 F rad=above(std::sqrt(above(remain/d[row])));F l=std::ceil(below(clo-rad)),h=std::floor(above(chi+rad));
 if(l>h)return false;if(l<-(1LL<<20)||h>(1LL<<20))throw std::runtime_error("coordinate range");
 for(I q=I(l);q<=I(h);q++){
  z[row]=q;F gap=0;if(F(q)<clo)gap=std::max(F(0),below(clo-F(q)));else if(F(q)>chi)gap=std::max(F(0),below(F(q)-chi));
  F term=std::max(F(0),below(d[row]*std::max(F(0),below(gap*gap))));F child=std::max(F(0),below(spent+term));
  if(child<=F(bound)/4 && visit(row-1,child))return true;
 }
 return false;
}
int main(int argc,char**argv){auto start=std::chrono::steady_clock::now();try{
 if(argc!=2)throw std::runtime_error("input required");if(std::numeric_limits<F>::digits<64||!std::numeric_limits<double>::is_iec559||std::fegetround()!=FE_TONEAREST)throw std::runtime_error("unsupported arithmetic");
 std::ifstream in(argv[1]);in>>n>>bound;if(n!=56||bound!=24)throw std::runtime_error("wrong case");r.resize(n);z.resize(n);G.assign(n,std::vector<I>(n));m.assign(n,std::vector<F>(n));d.resize(n);
 for(auto &a:r){in>>a;if(a<0||a>1)throw std::runtime_error("target");}for(auto &row:G)for(auto&a:row){in>>a;if(std::abs(a)>(1LL<<30))throw std::runtime_error("Gram");}
 auto interval=[&](){uint64_t a,b;in>>a>>b;double lo,hi;std::memcpy(&lo,&a,8);std::memcpy(&hi,&b,8);if(!std::isfinite(lo)||!std::isfinite(hi)||lo>hi)throw std::runtime_error("input interval");return std::pair<F,F>(F(lo),F(hi));};
 for(auto&a:d){auto x=interval();a=x.first;if(a<=0)throw std::runtime_error("diagonal");}
 for(auto&row:m)for(auto&a:row){auto x=interval();if(x.first < -16||x.second>16||x.second-x.first>std::ldexp(F(1),-47))throw std::runtime_error("coefficient");a=x.first;}
 if(!in)throw std::runtime_error("truncated input");bool found=visit(n-1,0);
 std::cout<<std::setprecision(18)<<"{\"status\":\""<<(found?"VECTOR":"EXHAUSTIVE_NO_VECTOR")<<"\",\"nodes\":"<<nodes<<",\"leaves\":"<<leaves<<",\"long_double_digits\":"<<std::numeric_limits<F>::digits<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";return found?10:0;
 }catch(const std::exception&e){std::cout<<"NOT_CERTIFIED "<<e.what()<<" nodes="<<nodes<<'\n';return 20;}}

#pragma once
struct Pair {
 long long n=0;long double x=0,y=0,xx=0,yy=0,xy=0;
 void add(double a,double b){++n;x+=a;y+=b;xx+=a*a;yy+=b*b;xy+=a*b;}
 void merge(const Pair&b){n+=b.n;x+=b.x;y+=b.y;xx+=b.xx;yy+=b.yy;xy+=b.xy;}
 Pair operator-(const Pair&b)const{return {n-b.n,x-b.x,y-b.y,xx-b.xx,yy-b.yy,xy-b.xy};}
 double corr()const{if(n<=2)return NaN;long double a=xx-x*x/n,b=yy-y*y/n;return a>0&&b>0?double((xy-x*y/n)/std::sqrt(a*b)):NaN;}
 double cov()const{return n>1?double((xy-x*y/n)/(n-1)):NaN;}
};
// Eight sufficient statistics: temperature, two crossing positions, sx, sy, before, after, angle.
// Per-channel/run centering and block score sums implement the same fixed-effect OLS as before.
struct Sums {
 long long n=0;long double s[8]={},q[8][8]={},centered[8][8]={};double lo=1e100,hi=-1e100;
 void add(const Hit&h){double a[8]={h.temperature,h.longitudinal,h.transverse,h.sx,h.sy,h.before,h.after,std::hypot(h.sx,h.sy)};long double delta[8];for(int i=0;i<8;++i)delta[i]=n?a[i]-s[i]/n:0;long double weight=static_cast<long double>(n)/(n+1);++n;for(int i=0;i<8;++i){s[i]+=a[i];for(int j=0;j<8;++j){q[i][j]+=a[i]*a[j];centered[i][j]+=delta[i]*delta[j]*weight;}}lo=std::min(lo,a[0]);hi=std::max(hi,a[0]);}
 void merge(const Sums&b){if(!b.n)return;if(!n){*this=b;return;}long double delta[8];for(int i=0;i<8;++i)delta[i]=b.s[i]/b.n-s[i]/n;long double weight=static_cast<long double>(n)*b.n/(n+b.n);for(int i=0;i<8;++i)for(int j=0;j<8;++j)centered[i][j]+=b.centered[i][j]+delta[i]*delta[j]*weight;n+=b.n;for(int i=0;i<8;++i){s[i]+=b.s[i];for(int j=0;j<8;++j)q[i][j]+=b.q[i][j];}lo=std::min(lo,b.lo);hi=std::max(hi,b.hi);}
};
using RegKey=std::array<int,5>; // sample group id stage model
using RegBlockKey=std::array<Long64_t,5>; // sample cell run file block
using MeanKey=std::array<int,3>; // sample cell run (-1 pooled)
using ScoreKey=std::array<Long64_t,8>; // RegKey + block
struct Score {long double u=0;void merge(const Score&b){u+=b.u;}};
struct Regression {
 std::string status="insufficient_statistics";long long n=0,fe=0,blocks=0;int rank=0;std::set<int>runs;
 long double xx[6][6]={},angleSum=0,angleQ=0,temperatureAngle=0,clusterQ=0;
 double range=0,beta=NaN,se=NaN,cluster=NaN,tx=NaN,ty=NaN,ta=NaN;
 std::vector<int>cols;double norm[5]={},coef[5]={},influence[5]={};bool fitted=false;
 void add(const Sums&a,int stage){int ids[6]={0,1,2,3,4,5+stage};n+=a.n;++fe;range=std::max(range,a.hi-a.lo);angleSum+=a.s[7];angleQ+=a.q[7][7];temperatureAngle+=a.centered[0][7];for(int j=0;j<6;++j)for(int k=0;k<6;++k)xx[j][k]+=a.centered[ids[j]][ids[k]];}
 void solve(int model,const Config&c){
  if(n<c.integer("min-slope-hits",40))return;
  if(range<c.number("min-temperature-span",.5)){status=model?"run_temperature_confounded":"insufficient_temperature_span";return;}
  for(int j=0;j<5;++j)norm[j]=double(std::max(0.L,xx[j][j]));
  auto corr=[&](int j){return norm[0]>0&&norm[j]>0?double(xx[0][j]/std::sqrt(norm[0]*norm[j])):NaN;};tx=corr(1);ty=corr(2);long double av=angleQ-angleSum*angleSum/n;ta=norm[0]>0&&av>0?double(temperatureAngle/std::sqrt(norm[0]*av)):NaN;
  if(norm[0]<1e-12){status="no_within_temperature_variance";return;}
  cols={0};for(int j=1;j<5;++j)if(norm[j]>1e-12)cols.push_back(j);int m=cols.size();TMatrixD mat(m,m);TVectorD rhs(m);
  for(int j=0;j<m;++j){rhs[j]=xx[cols[j]][5]/std::sqrt(norm[cols[j]]);for(int k=0;k<m;++k)mat(j,k)=xx[cols[j]][cols[k]]/std::sqrt(norm[cols[j]]*norm[cols[k]]);}
  TDecompSVD svd(mat);auto sig=svd.GetSig();for(int j=0;j<m;++j)rank+=sig[j]>sig[0]*1e-8;if(rank<m){status="temperature_position_angle_degenerate";return;}
  Bool_t ok;TMatrixD inv=svd.Invert(ok);if(!ok){status="singular_fit";return;}TVectorD b=inv*rhs;
  for(int j=0;j<m;++j){coef[cols[j]]=b[j]/std::sqrt(norm[cols[j]]);influence[cols[j]]=inv(0,j)/std::sqrt(norm[cols[j]]);}
  long long dof=n-fe-m;if(dof<=0){status="no_residual_dof";return;}
  long double rss=xx[5][5];for(int j=0;j<5;++j){rss-=2*coef[j]*xx[j][5];for(int k=0;k<5;++k)rss+=coef[j]*coef[k]*xx[j][k];}
  beta=coef[0];se=std::sqrt(std::max(0.L,inv(0,0)*rss/dof/norm[0]));fitted=true;
 }
 void finish(){if(!fitted)return;if(blocks>=4){cluster=std::sqrt(clusterQ*blocks/(blocks-1.L)/norm[0]);status="diagnostic_fit_not_calibration";}else status="fit_insufficient_independent_blocks";}
};
using CorrKey=std::array<int,6>;
using CorrBlockKey=std::array<Long64_t,9>;
using ResidualKey=std::array<Long64_t,13>;
struct Correlation {Pair total;long long blocks=0;Moment jack;};
struct ResidualSummary {
 Summary sum;long long blocks=0;long double ss=0,sn=0,nn=0;
 void add(const Summary&b){sum.merge(b);++blocks;ss+=b.value.s*b.value.s;sn+=b.value.s*b.value.n;nn+=static_cast<long double>(b.value.n)*b.value.n;}
};
struct RunSummary {Moment temp,value;void merge(const RunSummary&b){mergeMoment(temp,b.temp);mergeMoment(value,b.value);}};
struct Validation {Moment before,after,hits;void merge(const Validation&b){mergeMoment(before,b.before);mergeMoment(after,b.after);mergeMoment(hits,b.hits);}};
std::array<double,3>performance(const Moment&e,const Moment&h){return {e.mean(),e.mean()>0?e.sd()/e.mean():NaN,h.mean()};}

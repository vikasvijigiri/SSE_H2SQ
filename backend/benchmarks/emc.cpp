// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
// Independent equilibrium classical Monte Carlo. No SSE Geometry or weights.
#include "../include/read_input.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <numeric>

struct ClassicalLattice {
 int lx,ly,n;std::vector<std::array<int,4>> squares;std::vector<std::array<int,2>> pairs;
 explicit ClassicalLattice(int x,int y):lx(x),ly(y),n(x*y){
  std::map<std::pair<int,int>,int> index;
  for(int iy=0;iy<ly;++iy)for(int ix=0;ix<lx;++ix)if((ix+iy)%2==0){index[{ix,iy}]=squares.size();squares.push_back({ix+lx*iy,ix+lx*((iy+1)%ly),(ix+1)%lx+lx*((iy+1)%ly),(ix+1)%lx+lx*iy});}
  for(const auto& entry:index){auto [ix,iy]=entry.first;for(int dx:{-1,1})pairs.push_back({entry.second,index.at({(ix+dx+lx)%lx,(iy+1)%ly})});}
 }
 std::array<double,15> measure(const std::vector<int>& s,const Parameters& p)const{
  std::array<double,15> a{};std::vector<std::array<double,2>> dipoles;double m=0,x=0,y=0,px=0,py=0,energy=0.;
  for(int i=0;i<n;++i){int ix=i%lx,iy=i/lx;m+=.5*s[i]*((ix+iy)%2?-1:1)/n;x+=.5*s[i]*(ix%2?-1:1)/n;y+=.5*s[i]*(iy%2?-1:1)/n;}
  for(const auto& q:squares){int w=s[q[0]],b=s[q[1]],c=s[q[2]],d=s[q[3]];double vx=(w-b-c+d)/4.,vy=(w+b-c-d)/4.;
   energy+=-p.J0*w*b*c*d+p.J1*(w*c+b*d);dipoles.push_back({vx,vy});double sign=(q[0]%lx)%2?-1.:1.;px+=sign*vx/squares.size();py+=sign*vy/squares.size();
   a[7]+=double(w*b*c*d)/squares.size();a[8]+=double(w+b+c+d==0)/squares.size();a[14]+=.5*(vx*vx+vy*vy)/squares.size();
  }
  for(const auto& pair:pairs)energy+=p.J2*(dipoles[pair[0]][0]*dipoles[pair[1]][0]+dipoles[pair[0]][1]*dipoles[pair[1]][1]);
  a[0]=a[1]=energy/n;a[2]=a[0]*a[0];a[3]=a[2]*a[2];a[4]=m*m;a[5]=a[4]*a[4];a[6]=0.;a[9]=(4*a[8]-1)/3;a[10]=x*x+y*y;a[11]=a[10]*a[10];a[12]=px*px+py*py;a[13]=a[12]*a[12];return a;
 }
 double energy(const std::vector<int>& s,const Parameters& p)const{return n*measure(s,p)[0];}
};

int main(int argc,char** argv){
 if(argc!=5){std::cerr<<"input trace seed run|enumerate\n";return 2;}
 Parameters p;read_parameters(argv[1],p);if(p.hx!=0 || p.Lz!=1){std::cerr<<"Classical EMC requires h=0 and one layer\n";return 2;}
 ClassicalLattice lattice(p.Lx,p.Ly);RNG rng(std::stoull(argv[3]));
 std::ofstream out(argv[2]);out<<std::setprecision(17)<<"E,E_count,E2,E4,mz2,mz4,mx,plaquette,ice,rho,stripe2,stripe4,polar2,polar4,ct0\n";
 if(std::string(argv[4])=="enumerate"){
  if(p.Ns>20)return 2;
  std::vector<int> s(p.Ns);
  for(unsigned state=0;state<(1u<<p.Ns);++state){for(int i=0;i<p.Ns;++i)s[i]=(state>>i)&1?1:-1;auto a=lattice.measure(s,p);for(int k=0;k<15;++k)out<<(k?",":"")<<a[k];out<<'\n';}return 0;
 }
 double hot=std::min(p.beta,.1/std::max({1.,std::abs(p.J0),std::abs(p.J1),std::abs(p.J2)}));
 int replicas=std::max(2,int(std::ceil(std::log(p.beta/hot)/std::log(1.5)))+1);
 std::vector<double> betas(replicas),energies(replicas);std::vector<std::vector<int>> spins(replicas,std::vector<int>(p.Ns));std::vector<int> labels(replicas),phase(replicas);
 for(int r=0;r<replicas;++r){betas[r]=hot*std::pow(p.beta/hot,double(r)/(replicas-1));labels[r]=r;for(auto& s:spins[r])s=rng.uniform()<.5?-1:1;energies[r]=lattice.energy(spins[r],p);}
 std::vector<long long> attempted(replicas-1),accepted(replicas-1);long long trips=0;
 auto start=std::chrono::steady_clock::now();
 for(int step=-p.n_therm;step<p.n_measure;++step){
  if(step==0){std::fill(phase.begin(),phase.end(),0);std::fill(attempted.begin(),attempted.end(),0);std::fill(accepted.begin(),accepted.end(),0);trips=0;}
  for(int r=0;r<replicas;++r){
   for(int k=0;k<p.Ns;++k){int site=int(rng.uniform()*p.Ns);spins[r][site]*=-1;double next=lattice.energy(spins[r],p);if(std::log(std::max(rng.uniform(),1e-300))<-betas[r]*(next-energies[r]))energies[r]=next;else spins[r][site]*=-1;}
   unsigned mask=unsigned(rng.uniform()<.5)+2u*unsigned(rng.uniform()<.5);
   for(int i=0;i<p.Ns;++i)if(mask&(1u<<((i%p.Lx+i/p.Lx)%2)))spins[r][i]*=-1;
   double next=lattice.energy(spins[r],p);
   if(std::log(std::max(rng.uniform(),1e-300))<-betas[r]*(next-energies[r]))energies[r]=next;
   else for(int i=0;i<p.Ns;++i)if(mask&(1u<<((i%p.Lx+i/p.Lx)%2)))spins[r][i]*=-1;
  }
  for(int r=(step+p.n_therm)%2;r+1<replicas;r+=2){++attempted[r];if(std::log(std::max(rng.uniform(),1e-300))<(betas[r]-betas[r+1])*(energies[r]-energies[r+1])){std::swap(spins[r],spins[r+1]);std::swap(energies[r],energies[r+1]);std::swap(labels[r],labels[r+1]);++accepted[r];}}
  if(step>=0){if(phase[labels.back()]==2)++trips;phase[labels.back()]=1;if(phase[labels.front()]==1)phase[labels.front()]=2;auto a=lattice.measure(spins.back(),p);for(int k=0;k<15;++k)out<<(k?",":"")<<a[k];out<<'\n';}
 }
 std::cout<<std::setprecision(17)<<"{\"ensemble\":\"independent classical Boltzmann equilibrium MC\",\"replicas\":"<<replicas<<",\"round_trips\":"<<trips<<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"swap_attempts\":[";
 for(int r=0;r+1<replicas;++r)std::cout<<(r?",":"")<<attempted[r];
 std::cout<<"],\"swap_accepted\":[";
 for(int r=0;r+1<replicas;++r)std::cout<<(r?",":"")<<accepted[r];
 std::cout<<"]}";
}

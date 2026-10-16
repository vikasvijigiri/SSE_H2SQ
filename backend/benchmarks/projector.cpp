// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
// <+| exp(-T H) |+> with free time endpoints. Measurements insert at T/2.
#include "../include/updates.hpp"
#include "../include/read_input.hpp"
#include "../include/observables.hpp"
#include "../include/projector_estimators.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <numeric>

struct SSEKernelContract {
 static std::vector<int> initial_components(const SSEConfig& cfg,const SSEUpdates& u){std::vector<int> result(cfg.spins.spin.size());for(size_t s=0;s<result.size();++s)result[s]=u.first_leg[s]<0?u.n_legs+int(s):u.m_cluster_id[u.first_leg[s]];return result;}
 static std::array<int,2> field_components(const SSEUpdates& u,int rank){int base=u.vertex_base_vec[rank];return {u.m_cluster_id[base],u.m_cluster_id[base+1]};}
 static bool valid(const SSEConfig& cfg,const Parameters& p,const Geometry& g,SSEUpdates& u){
  SSEConfig state;state.spins=cfg.spins;int count=0;
  for(const auto& op:cfg.op_string){
   if(op.type==OpType::Identity)continue;
   ++count;
   if(op.type==OpType::Diagonal){if(!(u.diagonal_weight(int(op.subtype),op.index,state,p,g)>0))return false;}
   else if(op.subtype==OpSubtype::TransverseX_offdiag)state.spins.spin[op.index]*=-1;
   else return false;
  }
  return count==cfg.n_ops;
 }
};

struct MidpointBinomial {
 std::vector<double> mass,cdf,event;
 static std::vector<double> probabilities(int n){
  std::vector<double> a(n+1);int mode=n/2;a[mode]=1.;
  for(int k=mode;k<n;++k)a[k+1]=a[k]*(n-k)/(k+1.);
  for(int k=mode;k>0;--k)a[k-1]=a[k]*k/(n-k+1.);
  double z=std::accumulate(a.begin(),a.end(),0.);for(auto& x:a)x/=z;return a;
 }
 explicit MidpointBinomial(int n):mass(probabilities(n)),cdf(n+1),event(n?probabilities(n-1):std::vector<double>{}){std::partial_sum(mass.begin(),mass.end(),cdf.begin());}
 double interval(int first,int end)const{return cdf[end-1]-(first?cdf[first-1]:0.);}
};

// Order: diagonal energy/spin, mz2,mz4,plaquette,ice,rho,stripe2,stripe4,polar2,polar4,ct0.
std::array<double,11> diagonal_measure(const std::vector<int>& s,const Parameters& p,const Geometry& g){
 std::array<double,11> a{};double m=0,x=0,y=0,px=0,py=0;
 for(int i=0;i<p.Ns;++i){int ix=i%p.Lx,iy=(i/p.Lx)%p.Ly;m+=.5*s[i]*((ix+iy)%2?-1:1)/p.Ns;x+=.5*s[i]*(ix%2?-1:1)/p.Ns;y+=.5*s[i]*(iy%2?-1:1)/p.Ns;}
 a[1]=m*m;a[2]=m*m*m*m;a[6]=x*x+y*y;a[7]=a[6]*a[6];
 double np=g.blackPlaquettes.size();
 for(const auto& q:g.blackPlaquettes){int w=s[q.sites[0]],x0=s[q.sites[1]],y0=s[q.sites[2]],z=s[q.sites[3]];
  double vx=(w-x0-y0+z)/4.,vy=(w+x0-y0-z)/4.;int sign=(q.sites[0]%p.Lx)%2?-1:1;
  a[0]+=-p.J0*w*x0*y0*z+p.J1*(w*y0+x0*z);a[3]+=double(w*x0*y0*z)/np;a[4]+=double(w+x0+y0+z==0)/np;
  px+=sign*vx/np;py+=sign*vy/np;a[10]+=.5*(vx*vx+vy*vy)/np;
 }
 for(const auto& bond:g.j2bonds)a[0]+=p.J2*bond.J_val*s[bond.s1]*s[bond.s2];
 a[0]+=p.J2*g.j2_constant_term;
 for(const auto& bond:g.j3bonds)a[0]+=p.J3*s[bond[0]]*s[bond[1]];
 a[0]/=p.Ns;a[5]=(4*a[4]-1)/3;a[8]=px*px+py*py;a[9]=a[8]*a[8];return a;
}

int main(int argc,char** argv){
 if(argc!=4 && argc!=5){std::cerr<<"input trace seed [orbit-map]\n";return 2;}
 Parameters p;read_parameters(argv[1],p);if(p.hx==0){std::cerr<<"Projector requires h != 0 for unique ground validation\n";return 2;}
 Geometry g;g.build(p.Lx,p.Ly,p.Lz);RNG rng(std::stoull(argv[3])),measurement_rng(std::stoull(argv[3])^0x9e3779b97f4a7c15ULL);
 ProjectorDiagonal conditional(p,g);
 std::vector<int> orbit_map;std::vector<double> histogram;
 if(argc==5){
  if(p.Ns>16){std::cerr<<"Orbit histograms are limited to 16 spins\n";return 2;}
  orbit_map.resize(1u<<p.Ns);std::ifstream map_file(argv[4],std::ios::binary);
  map_file.read(reinterpret_cast<char*>(orbit_map.data()),orbit_map.size()*sizeof(int));
  if(!map_file || *std::min_element(orbit_map.begin(),orbit_map.end())<0 || *std::max_element(orbit_map.begin(),orbit_map.end())>=int(orbit_map.size()))return 2;
  histogram.resize(*std::max_element(orbit_map.begin(),orbit_map.end())+1);
 }
 ProjectorHistogram histogram_estimator(orbit_map,histogram);
 SSEConfig cfg;cfg.n_ops=0;cfg.max_ops=20;cfg.op_string.resize(20);cfg.spins.spin.resize(p.Ns);
 for(auto& s:cfg.spins.spin)s=rng.uniform()<.5?-1:1;
 SSEUpdates u;u.open_boundaries=true;long long checks=0,open_ends=0,saturations=0;
 for(int i=0;i<p.n_therm;++i){u.sweep(cfg,p,g,rng);u.adjust_operator_string(cfg);if(i%1000==0){++checks;if(!SSEKernelContract::valid(cfg,p,g,u))return 3;}}
 u.reset_symmetry_measurements();u.diagonal_seconds=u.cluster_seconds=u.gauge_seconds=u.canonical_seconds=0.;
 std::ofstream out(argv[2]);out<<std::setprecision(17)<<"E,E_count,E2,E4,mz2,mz4,mx,plaquette,ice,rho,stripe2,stripe4,polar2,polar4,ct0,raw_mx,raw_mz2,raw_mz4,raw_ice\n";
 std::map<int,MidpointBinomial> bins;auto started=std::chrono::steady_clock::now();
 for(int i=0;i<p.n_measure;++i){
  u.sweep(cfg,p,g,rng);if(i%1000==0){++checks;if(!SSEKernelContract::valid(cfg,p,g,u))return 3;}
  if(!bins.count(cfg.n_ops)){if(bins.size()>512)bins.clear();bins.emplace(cfg.n_ops,MidpointBinomial(cfg.n_ops));}
  const auto& bin=bins.at(cfg.n_ops);auto spin=cfg.spins.spin;auto components=SSEKernelContract::initial_components(cfg,u);int rank=0,begin=0;double density=0.,raw_density=0.;std::array<double,11> mid{},raw_mid{};
  int sampled_rank=std::min(cfg.n_ops,int(std::lower_bound(bin.cdf.begin(),bin.cdf.end(),measurement_rng.uniform())-bin.cdf.begin()));
  auto add=[&](int end){double w=bin.interval(begin,end);if(w>0){auto v=conditional.measure(spin,components);auto raw=diagonal_measure(spin,p,g);for(size_t k=0;k<mid.size();++k){mid[k]+=w*v[k];raw_mid[k]+=w*raw[k];}}
   if(!histogram.empty() && sampled_rank>=begin && sampled_rank<end)histogram_estimator.add(spin,components);
  };
  for(const auto& op:cfg.op_string)if(op.type!=OpType::Identity){
   if(op.subtype==OpSubtype::TransverseX_diag || op.subtype==OpSubtype::TransverseX_offdiag){
    add(rank+1);auto cid=SSEKernelContract::field_components(u,rank);bool off=op.type==OpType::OffDiagonal;
    double intensity=cfg.n_ops*bin.event[rank]/p.beta;density+=intensity*(cid[0]==cid[1]?double(off):.5);raw_density+=intensity*double(off);
    if(off)spin[op.index]*=-1;
    components[op.index]=cid[1];begin=rank+1;
   }
   ++rank;
  }
  add(rank+1);open_ends+=(spin!=cfg.spins.spin);saturations+=(cfg.n_ops==cfg.max_ops);
  double mx=density/(p.hx*p.Ns),energy=mid[0]-p.hx*mx;
  Observables obs;obs.accumulate(cfg,p,g);
  out<<energy<<','<<obs.E<<','<<obs.E2<<','<<obs.E4<<','<<mid[1]<<','<<mid[2]<<','<<mx;
  for(int k=3;k<11;++k)out<<','<<mid[k];
  out<<','<<raw_density/(p.hx*p.Ns)<<','<<raw_mid[1]<<','<<raw_mid[2]<<','<<raw_mid[4];
  out<<'\n';
 }
 Observables obs;obs.accumulate(cfg,p,g);double shift=obs.E*p.Ns+cfg.n_ops/p.beta;
 if(!histogram.empty()){std::ofstream hist_file(std::string(argv[2])+".hist.bin",std::ios::binary);hist_file.write(reinterpret_cast<const char*>(histogram.data()),histogram.size()*sizeof(double));if(!hist_file)return 4;}
 std::cout<<std::setprecision(17)<<"{\"ensemble\":\"open-boundary two-sided projector\",\"projection_total\":"<<p.beta<<",\"projection_per_side\":"<<p.beta/2<<",\"operator_cutoff\":"<<cfg.max_ops<<",\"norm_mean_bound\":"<<u.site_symmetry->norm_mean<<",\"energy_shift_total\":"<<shift<<",\"inverse_temperature\":"<<p.beta<<",\"spins\":"<<p.Ns<<",\"local_skip_probability\":"<<u.local_skip_probability<<",\"configuration_checks\":"<<checks<<",\"open_endpoint_rows\":"<<open_ends<<",\"cutoff_saturations\":"<<saturations<<",\"wall_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<",\"site_symmetry_orbit_counts\":[";
 for(int k=0;k<4;++k)std::cout<<(k?",":"")<<u.site_symmetry->orbit_counts[k];
 std::cout<<"],\"histogram_orbits\":"<<histogram.size()<<",\"histogram_rows\":"<<std::accumulate(histogram.begin(),histogram.end(),0.)
 <<",\"diagonal_skips\":"<<u.diagonal_skips<<",\"canonical_skips\":"<<u.canonical_skips
 <<",\"update_seconds\":{\"diagonal\":"<<u.diagonal_seconds<<",\"cluster\":"<<u.cluster_seconds
 <<",\"gauge\":"<<u.gauge_seconds<<",\"canonical\":"<<u.canonical_seconds<<",\"symmetry\":"<<u.symmetry_seconds<<"}}";
}

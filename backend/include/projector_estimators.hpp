// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#ifndef SSE_PROJECTOR_ESTIMATORS_HPP
#define SSE_PROJECTOR_ESTIMATORS_HPP
#include "variables.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

// At a fixed graph and its relative spin labels every connected component
// has an independent, uniform orientation. This diagonal conditional average
// remains valid with frustrated relative signs. It does NOT use a bridge-X
// insertion formula. Global even/odd spin orbits are also averaged exactly.
class ProjectorDiagonal {
 const Parameters& p;const Geometry& g;
 std::vector<std::array<double,6>> weights;
public:
 ProjectorDiagonal(const Parameters& parameters,const Geometry& geometry):p(parameters),g(geometry),weights(p.Ns){
  for(int i=0;i<p.Ns;++i){int x=i%p.Lx,y=(i/p.Lx)%p.Ly;weights[i]={((x+y)%2?-1.:1.)/(2*p.Ns),1./(2*p.Ns),(x%2?-1.:1.)/(2*p.Ns),(y%2?-1.:1.)/(2*p.Ns),0.,0.};}
  const double cx[4]={.25,-.25,-.25,.25},cy[4]={.25,.25,-.25,-.25};
  for(const auto& q:g.blackPlaquettes){double sign=(q.sites[0]%p.Lx)%2?-1.:1.;for(int k=0;k<4;++k){weights[q.sites[k]][4]+=sign*cx[k]/g.blackPlaquettes.size();weights[q.sites[k]][5]+=sign*cy[k]/g.blackPlaquettes.size();}}
 }
 std::array<double,11> measure(const std::vector<int>& spin,const std::vector<int>& component)const{
  std::vector<int> ids;std::vector<std::array<double,6>> coefficients;
  for(int site=0;site<p.Ns;++site){auto found=std::find(ids.begin(),ids.end(),component[site]);size_t c=found-ids.begin();if(c==ids.size()){ids.push_back(component[site]);coefficients.push_back({});}for(int mode=0;mode<6;++mode)coefficients[c][mode]+=spin[site]*weights[site][mode];}
  std::array<double,6> second{},fourth{};
  for(const auto& c:coefficients)for(int k=0;k<6;++k){second[k]+=c[k]*c[k];fourth[k]+=c[k]*c[k]*c[k]*c[k];}
  for(int k=0;k<6;++k)fourth[k]=3*second[k]*second[k]-2*fourth[k];
  auto vector_fourth=[&](int a,int b){double cross=0,correction=0;for(const auto& c:coefficients){cross+=c[a]*c[b];double r=c[a]*c[a]+c[b]*c[b];correction+=r*r;}return std::pow(second[a]+second[b],2)+2*(second[a]*second[a]+second[b]*second[b]+2*cross*cross)-2*correction;};
  auto pair=[&](int a,int b){return component[a]==component[b]?double(spin[a]*spin[b]):0.;};
  std::array<double,11> out{};out[1]=.5*(second[0]+second[1]);out[2]=.5*(fourth[0]+fourth[1]);out[6]=second[2]+second[3];out[7]=vector_fourth(2,3);out[8]=second[4]+second[5];out[9]=vector_fourth(4,5);
  double np=g.blackPlaquettes.size();
  for(const auto& q:g.blackPlaquettes){int a=q.sites[0],b=q.sites[1],c=q.sites[2],d=q.sites[3];
   bool even=(component[a]==component[b] && component[c]==component[d]) || (component[a]==component[c] && component[b]==component[d]) || (component[a]==component[d] && component[b]==component[c]);
   double product=even?double(spin[a]*spin[b]*spin[c]*spin[d]):0.,diagonals=pair(a,c)+pair(b,d);
   out[0]+=-p.J0*product+p.J1*diagonals;out[3]+=product/np;out[4]+=(3-diagonals+3*product)/(8*np);out[10]+=(.25-.125*diagonals)/np;
  }
  for(const auto& bond:g.j2bonds)out[0]+=p.J2*bond.J_val*pair(bond.s1,bond.s2);
  out[0]+=p.J2*g.j2_constant_term;for(const auto& bond:g.j3bonds)out[0]+=p.J3*pair(bond[0],bond[1]);
  out[0]/=p.Ns;out[5]=(4*out[4]-1)/3;return out;
 }
};

// Rao-Blackwell histogram at one exactly sampled central-time index. Only
// components meeting that time slice need enumeration. No ED probabilities
// enter this calculation: the supplied map contains geometry orbits only.
class ProjectorHistogram {
 const std::vector<int>& map;std::vector<double>& counts;std::vector<double> uniform;
 std::map<std::vector<unsigned>,std::vector<std::pair<int,double>>> cache;
public:
 ProjectorHistogram(const std::vector<int>& mapping,std::vector<double>& output):map(mapping),counts(output),uniform(output.size()){
  for(int orbit:map)uniform[orbit]+=1./map.size();
 }
 void add(const std::vector<int>& spin,const std::vector<int>& component){
 std::vector<int> ids;std::vector<unsigned> masks;unsigned state=0;
 for(size_t i=0;i<spin.size();++i){if(spin[i]>0)state|=1u<<i;auto found=std::find(ids.begin(),ids.end(),component[i]);size_t c=found-ids.begin();if(c==ids.size()){ids.push_back(component[i]);masks.push_back(0);}masks[c]|=1u<<i;}
 if(ids.size()==spin.size()){
  for(size_t i=0;i<counts.size();++i)counts[i]+=uniform[i];
  return;
 }
 unsigned canonical=state;for(unsigned mask:masks)if(canonical&(mask&(~mask+1)))canonical^=mask;
 std::vector<unsigned> key=masks;key.push_back(canonical);
 auto found=cache.end();if(masks.size()>=10)found=cache.find(key);
 if(found!=cache.end()){for(const auto& item:found->second)counts[item.first]+=item.second;return;}
 unsigned alternatives=1u<<masks.size();double weight=1./alternatives;state=canonical;
 std::vector<int> tally;if(masks.size()>=10)tally.resize(counts.size());
 for(unsigned i=0;i<alternatives;++i){if(i){unsigned changed=i&(~i+1);int bit=0;while((changed>>bit)!=1u)++bit;state^=masks[bit];}if(tally.empty())counts[map[state]]+=weight;else ++tally[map[state]];}
 if(!tally.empty()){if(cache.size()>=1024)cache.clear();std::vector<std::pair<int,double>> distribution;for(size_t i=0;i<tally.size();++i)if(tally[i]){double probability=tally[i]*weight;counts[i]+=probability;distribution.emplace_back(int(i),probability);}cache.emplace(std::move(key),std::move(distribution));}
 }
};
#endif

// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#ifndef SSE_SITE_SYMMETRY_HPP
#define SSE_SITE_SYMMETRY_HPP
#include "variables.hpp"
#include <stdexcept>
#include <vector>
#include <array>
#include <cmath>
#include <limits>

// Exact heat bath over the Z2 x Z2 orbit of the two SPIN sublattices.
// These are site parity (x+y)%2, distinct from the A/B molecular convention.
// Every retained J1/J2/J3 bond has both ends on one spin sublattice and every
// J0 plaquette has two sites on each. Transverse flips commute with this move.
// The SSE operator string and every local matrix element remain unchanged.
class SSESiteSymmetry {
    std::vector<int> parity;
public:
    std::array<long long,4> orbit_counts{};
    int lx,ly,lz;
    std::array<double,6> parameter_key;
    double norm_mean;
    int cutoff_floor;
    static double expansion_bound(const Parameters& p,const Geometry& g){
        double bound=2*g.blackPlaquettes.size()*(std::abs(p.J0)+2*std::abs(p.J1))+2*p.Ns*std::abs(p.hx);
        for(const auto& b:g.j2bonds)bound+=2*std::abs(b.J_val*p.J2);
        bound+=2*g.j3bonds.size()*std::abs(p.J3);
        return p.beta*bound*(1+1024*std::numeric_limits<double>::epsilon());
    }
    SSESiteSymmetry(const Parameters& p,const Geometry& g):parity(p.Ns),lx(p.Lx),ly(p.Ly),lz(p.Lz){
        parameter_key={p.beta,p.J0,p.J1,p.J2,p.J3,p.hx};
        norm_mean=expansion_bound(p,g);
        double floor=std::ceil(norm_mean+12*std::sqrt(norm_mean+1)+32);
        if(!std::isfinite(floor) || floor>std::numeric_limits<int>::max())throw std::runtime_error("Requested SSE cutoff exceeds integer capacity.");
        cutoff_floor=int(floor);
        if(p.Lx%2 || p.Ly%2)throw std::runtime_error("Site-sublattice symmetry requires even periodic dimensions.");
        for(int s=0;s<p.Ns;++s)parity[s]=(s%p.Lx+(s/p.Lx)%p.Ly)%2;
        for(const auto& q:g.blackPlaquettes){
            int count=0;for(int s:q.sites)count+=parity[s];
            if(count!=2 || parity[q.sites[0]]!=parity[q.sites[2]] || parity[q.sites[1]]!=parity[q.sites[3]])
                throw std::runtime_error("Plaquette does not preserve independent site-sublattice flips.");
        }
        for(const auto& b:g.j2bonds)if(parity[b.s1]!=parity[b.s2])
            throw std::runtime_error("J2 bond breaks independent site-sublattice flips.");
        for(const auto& b:g.j3bonds)if(parity[b[0]]!=parity[b[1]])
            throw std::runtime_error("J3 bond breaks independent site-sublattice flips.");
    }
    bool matches(const Parameters& p) const{return lx==p.Lx && ly==p.Ly && lz==p.Lz && parameter_key==std::array<double,6>{p.beta,p.J0,p.J1,p.J2,p.J3,p.hx};}
    void ensure_cutoff(SSEConfig& cfg) const {
        if(cfg.max_ops<cutoff_floor){cfg.op_string.resize(cutoff_floor);cfg.max_ops=cutoff_floor;}
    }
    void apply_mask(SSEConfig& cfg,unsigned mask){
        if(mask>3 || cfg.spins.spin.size()!=parity.size())throw std::runtime_error("Invalid site-symmetry orbit/configuration.");
        ++orbit_counts[mask];
        for(size_t s=0;s<parity.size();++s)if(mask&(1u<<parity[s]))cfg.spins.spin[s]*=-1;
    }
    void sweep(SSEConfig& cfg,RNG& rng){
        unsigned mask=unsigned(rng.uniform()<.5)+2u*unsigned(rng.uniform()<.5);
        apply_mask(cfg,mask);
    }
};
#endif

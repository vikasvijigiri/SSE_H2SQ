// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#ifndef SSE_TEMPERING_HPP
#define SSE_TEMPERING_HPP
#include "updates.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

// Identical couplings and equal padding lengths give the exact SSE exchange
// ratio (beta_i/beta_j)^(n_j-n_i). Local weights and constant shifts cancel.
// Padding adapts only during thermalization, then stays fixed in measurements.
class SSETempering {
public:
    static double swap_log_ratio(double beta_a,double beta_b,int n_a,int n_b) {
        return (n_b-n_a)*std::log(beta_a/beta_b);
    }
    std::vector<Parameters> parameters;
    std::vector<SSEConfig> warm_configs;
    std::vector<SSEUpdates> warm_updates;
    std::vector<SSEConfig*> configs;
    std::vector<SSEUpdates*> updates;
    std::vector<long long> attempts, accepted;
    long long round_trips=0;
private:
    std::vector<int> labels, journey;
    long long steps=0;
public:
    SSETempering(const Parameters& target, SSEConfig& cold, SSEUpdates& cold_update, RNG& rng) {
        double scale=std::max({1.,std::abs(target.J0),std::abs(target.J1),std::abs(target.J2),std::abs(target.J3),std::abs(target.hx)});
        double beta=std::min(.25/scale,target.beta);
        while(beta<target.beta*(1-1e-12)) {
            parameters.push_back(target);parameters.back().beta=beta;parameters.back().T=1/beta;
            beta=std::min(target.beta,beta*1.6);
        }
        parameters.push_back(target);
        warm_configs.resize(parameters.size()-1);warm_updates.resize(warm_configs.size());
        for(size_t i=0;i<warm_configs.size();++i){
            auto& c=warm_configs[i];c.max_ops=20;c.n_ops=0;c.op_string.resize(20);c.spins.spin.resize(target.Ns);
            for(auto& spin:c.spins.spin)spin=rng.uniform()<.5?-1:1;
            configs.push_back(&c);updates.push_back(&warm_updates[i]);
        }
        configs.push_back(&cold);updates.push_back(&cold_update);
        for(auto* u:updates){u->spatial_premarking=true;u->canonical_enabled=true;}
        labels.resize(configs.size());std::iota(labels.begin(),labels.end(),0);
        journey.assign(configs.size(),0);attempts.assign(configs.size()-1,0);accepted=attempts;
    }
    void sweep(const Geometry& g,RNG& rng,bool thermal) {
        for(size_t i=0;i<configs.size();++i){
            updates[i]->sweep(*configs[i],parameters[i],g,rng);
            if(thermal)updates[i]->adjust_operator_string(*configs[i]);
        }
        if(thermal){
            int common=0;for(auto* c:configs)common=std::max(common,c->max_ops);
            for(auto* c:configs){c->op_string.resize(common);c->max_ops=common;}
        }
        // Alternate disjoint even/odd neighboring exchanges every ten sweeps.
        if(++steps%10==0){
            size_t parity=(steps/10)%2;
            for(size_t i=parity;i+1<configs.size();i+=2){
                ++attempts[i];
                double ratio=swap_log_ratio(parameters[i].beta,parameters[i+1].beta,configs[i]->n_ops,configs[i+1]->n_ops);
                if(std::log(std::max(rng.uniform(),1e-300))<std::min(0.,ratio)){
                    std::swap(*configs[i],*configs[i+1]);std::swap(labels[i],labels[i+1]);++accepted[i];
                }
            }
        }
        if(configs.size()>1){
            int hot=labels.front(),cold=labels.back();
            if(journey[hot]==2)++round_trips;
            journey[hot]=1;
            if(journey[cold]==1)journey[cold]=2;
        }
    }
    void reset_measurements(){
        std::fill(attempts.begin(),attempts.end(),0);std::fill(accepted.begin(),accepted.end(),0);
        std::fill(journey.begin(),journey.end(),0);round_trips=0;
        for(auto* u:updates){u->diagonal_seconds=u->cluster_seconds=u->gauge_seconds=u->canonical_seconds=0;u->canonical_attempts=u->canonical_accepted=0;}
    }
};
#endif

// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#include "../include/variables.hpp"
#include "../include/updates.hpp"
#include "../include/observables.hpp"
#include "../include/read_input.hpp"
#include "../include/tempering.hpp"
#include <memory>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <cmath>

struct SSEKernelContract {
    static bool valid(const SSEConfig& cfg, const Parameters& p, const Geometry& g, SSEUpdates& u) {
        SSEConfig propagated;
        propagated.spins = cfg.spins;
        int count = 0;
        for (const auto& op : cfg.op_string) {
            if (op.type == OpType::Identity) continue;
            ++count;
            if (op.type == OpType::Diagonal) {
                if (!(u.diagonal_weight(int(op.subtype), op.index, propagated, p, g) > 0)) return false;
            } else {
                if (op.subtype != OpSubtype::TransverseX_offdiag) return false;
                propagated.spins.spin[op.index] *= -1;
            }
        }
        return count == cfg.n_ops && propagated.spins.spin == cfg.spins.spin;
    }
};

struct DiagonalMeasurement {
    double m=0, stripe2=0, polar2=0, parity=0, ice=0, px=0, py=0;
    std::vector<double> molecular_x;
};

DiagonalMeasurement measure(const std::vector<int>& spin, const Parameters& p, const Geometry& g) {
    DiagonalMeasurement r;
    double sx=0, sy=0;
    for(int i=0;i<p.Ns;++i) {
        int x=i%p.Lx, y=(i/p.Lx)%p.Ly;
        r.m += .5*spin[i]*((x+y)%2?-1:1)/p.Ns;
        sx += .5*spin[i]*(x%2?-1:1)/p.Ns;
        sy += .5*spin[i]*(y%2?-1:1)/p.Ns;
    }
    r.stripe2=sx*sx+sy*sy;
    for(const auto& q:g.blackPlaquettes) {
        int a=spin[q.sites[0]], b=spin[q.sites[1]], c=spin[q.sites[2]], d=spin[q.sites[3]];
        double px=(a-b-c+d)/4., py=(a+b-c-d)/4.;
        int sign=(q.sites[0]%p.Lx)%2?-1:1;
        r.molecular_x.push_back(px);
        r.px+=sign*px/g.blackPlaquettes.size(); r.py+=sign*py/g.blackPlaquettes.size();
        r.parity+=double(a*b*c*d)/g.blackPlaquettes.size();
        r.ice+=double(a+b+c+d==0)/g.blackPlaquettes.size();
    }
    r.polar2=r.px*r.px+r.py*r.py;
    return r;
}

int main(int argc,char** argv) {
    if(argc!=5){std::cerr<<"input trace algorithm seed\n";return 2;}
    Parameters p;read_parameters(argv[1],p);
    Geometry g;g.build(p.Lx,p.Ly,p.Lz);
    RNG rng(std::stoull(argv[4])), measurement_rng(std::stoull(argv[4]) ^ 0x9e3779b97f4a7c15ULL);
    SSEConfig cfg;cfg.max_ops=20;cfg.n_ops=0;cfg.op_string.resize(20);cfg.spins.spin.resize(p.Ns);
    for(auto& spin:cfg.spins.spin)spin=rng.uniform()<.5?-1:1;
    SSEUpdates update;
    update.spatial_premarking=std::string(argv[3])!="vertex";
    update.canonical_enabled=std::string(argv[3])=="canonical";
    std::unique_ptr<SSETempering> tempering;
    if(std::string(argv[3])=="tempered")tempering=std::make_unique<SSETempering>(p,cfg,update,rng);
    long long configuration_checks=0;
    for(int step=0;step<p.n_therm;++step){
        if(tempering)tempering->sweep(g,rng,true);
        else {update.sweep(cfg,p,g,rng);update.adjust_operator_string(cfg);}
        if(step%1000==0){++configuration_checks;if(!SSEKernelContract::valid(cfg,p,g,update)){std::cerr<<"Invalid thermalized SSE configuration";return 3;}}
    }
    update.diagonal_seconds=update.cluster_seconds=update.gauge_seconds=update.canonical_seconds=0;
    update.canonical_attempts=update.canonical_accepted=0;
    if(tempering)tempering->reset_measurements();
    std::ofstream out(argv[2]);out<<std::setprecision(17);
    out<<"E,E2,E4,mz,mz2,mz4,mx,plaquette,ice,rho,stripe2,stripe4,polar2,polar4,bar_m,bar_m2,bar_px,bar_py,bar_p2,nx,nxf2,ct0,ct1,ct2,ct3,ct4\n";
    long long non_global=0,saturated=0;
    auto start=std::chrono::steady_clock::now();
    for(int step=0;step<p.n_measure;++step) {
        auto before=cfg.spins.spin;
        if(tempering)tempering->sweep(g,rng,false);
        else update.sweep(cfg,p,g,rng);
        if(step%1000==0){++configuration_checks;if(!SSEKernelContract::valid(cfg,p,g,update)){std::cerr<<"Invalid measured SSE configuration";return 3;}}
        int flips=0;for(int i=0;i<p.Ns;++i)flips+=(before[i]!=cfg.spins.spin[i]);
        non_global+=(flips>0 && flips<p.Ns);saturated+=(cfg.n_ops==cfg.max_ops);
        Observables obs;obs.accumulate(cfg,p,g);
        auto base=measure(cfg.spins.spin,p,g);
        std::vector<double> times;
        times.reserve(cfg.n_ops);
        if(p.hx!=0)for(int i=0;i<cfg.n_ops;++i)times.push_back(p.beta*measurement_rng.uniform());
        std::sort(times.begin(),times.end());
        std::vector<std::pair<double,int>> events;
        int active=0, nx=0;
        for(const auto& op:cfg.op_string)if(op.type!=OpType::Identity){
            if(op.type==OpType::OffDiagonal){++nx;events.push_back({times[active],op.index});}
            ++active;
        }
        auto spin=cfg.spins.spin;auto current=base;
        double last=0,bar_m=0,bar_px=0,bar_py=0;
        for(auto event:events){
            double dt=(event.first-last)/p.beta;
            bar_m+=dt*current.m;bar_px+=dt*current.px;bar_py+=dt*current.py;
            spin[event.second]*=-1;current=measure(spin,p,g);last=event.first;
        }
        double dt=(p.beta-last)/p.beta;
        bar_m+=dt*current.m;bar_px+=dt*current.px;bar_py+=dt*current.py;
        double ct[5]={};double fractions[5]={0,.125,.25,.375,.5};
        spin=cfg.spins.spin;size_t next=0;
        for(int t=0;t<5;++t){
            while(next<events.size() && events[next].first<=fractions[t]*p.beta){spin[events[next].second]*=-1;++next;}
            auto at_time=measure(spin,p,g);
            for(size_t q=0;q<base.molecular_x.size();++q)ct[t]+=base.molecular_x[q]*at_time.molecular_x[q]/base.molecular_x.size();
        }
        out<<obs.E<<','<<obs.E2<<','<<obs.E4<<','<<base.m<<','<<base.m*base.m<<','<<std::pow(base.m,4)<<','<<obs.mx<<','<<base.parity<<','<<base.ice<<','<<(4*base.ice-1)/3<<','<<base.stripe2<<','<<base.stripe2*base.stripe2<<','<<base.polar2<<','<<base.polar2*base.polar2<<','<<bar_m<<','<<bar_m*bar_m<<','<<bar_px<<','<<bar_py<<','<<bar_px*bar_px+bar_py*bar_py<<','<<nx<<','<<nx*(nx-1.);
        for(auto value:ct)out<<','<<value;
        out<<'\n';
    }
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if(tempering)for(const auto& u:tempering->warm_updates){
        update.diagonal_seconds+=u.diagonal_seconds;update.cluster_seconds+=u.cluster_seconds;
        update.gauge_seconds+=u.gauge_seconds;update.canonical_seconds+=u.canonical_seconds;
        update.canonical_attempts+=u.canonical_attempts;update.canonical_accepted+=u.canonical_accepted;
    }
    std::cout<<std::setprecision(17)<<"{\"replicas\":"<<(tempering?tempering->configs.size():1)<<",\"round_trips\":"<<(tempering?tempering->round_trips:0)<<",\"beta_ladder\":[";
    if(tempering)for(size_t i=0;i<tempering->parameters.size();++i)std::cout<<(i?",":"")<<tempering->parameters[i].beta;
    std::cout<<"],\"swap_attempts\":[";
    if(tempering)for(size_t i=0;i<tempering->attempts.size();++i)std::cout<<(i?",":"")<<tempering->attempts[i];
    std::cout<<"],\"swap_accepted\":[";
    if(tempering)for(size_t i=0;i<tempering->accepted.size();++i)std::cout<<(i?",":"")<<tempering->accepted[i];
    std::cout<<"],";
    std::cout<<std::setprecision(17)<<"\"measurement_seconds\":"<<elapsed<<",\"configuration_checks\":"<<configuration_checks<<",\"non_global_moves\":"<<non_global<<",\"cutoff_saturations\":"<<saturated<<",\"update_seconds\":{\"diagonal\":"<<update.diagonal_seconds<<",\"cluster\":"<<update.cluster_seconds<<",\"gauge\":"<<update.gauge_seconds<<",\"canonical\":"<<update.canonical_seconds<<"},\"canonical_attempts\":"<<update.canonical_attempts<<",\"canonical_accepted\":"<<update.canonical_accepted<<"}\n";
}

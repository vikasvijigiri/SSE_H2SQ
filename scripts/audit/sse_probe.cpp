// Audit-only access to cluster instrumentation; production sources are unchanged.
#include <bits/stdc++.h>
#define private public
#include "../../backend/include/variables.hpp"
#include "../../backend/include/updates.hpp"
#undef private
#include "../../backend/include/observables.hpp"

double physical(const SSEConfig& cfg, const Parameters& prm, const Geometry& g) {
    double e=0;
    for (const auto& p:g.blackPlaquettes) {
        int a=cfg.spins[p.sites[0]]*cfg.spins[p.sites[2]];
        int b=cfg.spins[p.sites[1]]*cfg.spins[p.sites[3]];
        e-=prm.J0*a*b;
        e+=std::abs(prm.J1)*(1-a-b-a*b);
    }
    for (const auto& b:g.j2bonds) e+=prm.J2*b.J_val*cfg.spins[b.s1]*cfg.spins[b.s2];
    e+=prm.J2*g.j2_constant_term;
    for (const auto& b:g.j3bonds) e+=prm.J3*cfg.spins[b[0]]*cfg.spins[b[1]];
    return e/prm.Ns;
}

int main(int argc,char** argv) {
    std::cout<<std::setprecision(17);
    if (argc>=5 && std::string(argv[1])=="geometry") {
        Geometry g; g.build(std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]));
        std::cout<<"{\"plaquettes\":[";
        for(size_t i=0;i<g.blackPlaquettes.size();++i) {
            if(i)std::cout<<",";
            auto s=g.blackPlaquettes[i].sites;
            std::cout<<"["<<s[0]<<","<<s[1]<<","<<s[2]<<","<<s[3]<<"]";
        }
        std::cout<<"],\"j2bonds\":[";
        for(size_t i=0;i<g.j2bonds.size();++i) {
            if(i)std::cout<<",";
            auto b=g.j2bonds[i];std::cout<<"["<<b.s1<<","<<b.s2<<","<<b.J_val<<"]";
        }
        std::cout<<"],\"j2constant\":"<<g.j2_constant_term<<",\"j3bonds\":[";
        for(size_t i=0;i<g.j3bonds.size();++i){if(i)std::cout<<",";auto b=g.j3bonds[i];std::cout<<"["<<b[0]<<","<<b[1]<<"]";}
        std::cout<<"]}\n";
        return 0;
    }
    if(argc!=11 && argc!=13) return 2;
    Parameters prm{};
    prm.Lx=argc==13?std::stoi(argv[11]):4;prm.Ly=argc==13?std::stoi(argv[12]):4;prm.Lz=1;prm.Ns=prm.Lx*prm.Ly;prm.Np=prm.Ns/2;
    prm.beta=std::stod(argv[1]);prm.J0=std::stod(argv[2]);prm.J1=std::stod(argv[3]);prm.J2=std::stod(argv[4]);prm.hx=std::stod(argv[5]);
    prm.n_therm=std::stoi(argv[6]);prm.n_measure=std::stoi(argv[7]);
    RNG rng;rng.gen.seed(std::stoull(argv[8]));
    Geometry geom;geom.build(prm.Lx,prm.Ly,1);
    SSEConfig cfg;cfg.max_ops=20;cfg.n_ops=0;cfg.op_string.resize(20);cfg.spins.spin.resize(prm.Ns);
    std::string initial=argv[9];
    for(int i=0;i<prm.Ns;++i) cfg.spins.spin[i]=initial=="cold"?1:initial=="stripe"?(i%prm.Lx<prm.Lx/2?1:-1):(rng.uniform()<0.5?-1:1);
    SSEUpdates update;Observables obs;obs.reset();
    long long changed=0,non_global=0,invalid=0,saturated=0;
    double largest_sum=0,physical_sum=0;
    std::unordered_set<unsigned long long> audit_states;
    std::ofstream trace(argv[10]);trace<<"step,energy,m2,n_ops,largest_cluster,flip_count\n";
    for(int step=-prm.n_therm;step<prm.n_measure;++step) {
        auto before=cfg.spins.spin;
        update.sweep(cfg,prm,geom,rng);
        if(step<0){update.adjust_operator_string(cfg);continue;}
        int flips=0;for(int i=0;i<prm.Ns;++i)flips+=(cfg.spins.spin[i]!=before[i]);
        changed+=(flips>0);non_global+=(flips>0 && flips<prm.Ns);
        std::vector<int> sizes(update.m_cluster_flipped.size(),0);
        for(auto cid:update.m_cluster_id)if(cid>=0)++sizes[cid];
        double largest=sizes.empty()?0:double(*std::max_element(sizes.begin(),sizes.end()))/update.n_legs;
        largest_sum+=largest;saturated+=(cfg.n_ops==cfg.max_ops);
        auto spin=cfg.spins.spin;
        for(const auto& op:cfg.op_string) {
            if(op.type==OpType::OffDiagonal)spin[op.index]*=-1;
            else if(op.type==OpType::Diagonal) {
                SSEConfig view;view.spins.spin=spin;
                if(update.diagonal_weight(int(op.subtype),op.index,view,prm,geom)<=0)++invalid;
            }
        }
        unsigned long long audit_mask=0;
        for(int i=0;i<prm.Ns;++i) if(cfg.spins.spin[i]>0) audit_mask |= (1ULL<<i);
        audit_states.insert(audit_mask);
        obs.accumulate(cfg,prm,geom);
        double mz=0;for(int i=0;i<prm.Ns;++i)mz+=0.5*cfg.spins.spin[i]*(((i%prm.Lx+i/prm.Lx)%2)?-1:1);mz/=prm.Ns;
        double e=physical(cfg,prm,geom);physical_sum+=e;
        trace<<step<<","<<e<<","<<mz*mz<<","<<cfg.n_ops<<","<<largest<<","<<flips<<"\n";
    }
    obs.normalize();
    std::cout<<"{\"E\":"<<obs.E<<",\"E2\":"<<obs.E2<<",\"E4\":"<<obs.E4<<",\"m2\":"<<obs.mz2<<",\"m4\":"<<obs.mz4<<",\"Mx\":"<<obs.mx<<",\"diagonal_energy\":"<<physical_sum/prm.n_measure<<",\"largest_cluster\":"<<largest_sum/prm.n_measure<<",\"non_global_moves\":"<<non_global<<",\"changed_moves\":"<<changed<<",\"invalid_vertices\":"<<invalid<<",\"cutoff_saturations\":"<<saturated<<",\"unique_states\":"<<audit_states.size()<<",\"max_ops\":"<<cfg.max_ops<<"}\n";
}

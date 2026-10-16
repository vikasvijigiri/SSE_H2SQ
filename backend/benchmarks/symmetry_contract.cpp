// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
// Independently verify the exact Z2 x Z2 site-sublattice symmetry of every
// shifted diagonal matrix element, not only the summed physical energy.
#include "../include/updates.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

struct SSEKernelContract {
    static std::vector<double> weights(SSEUpdates& update,const SSEConfig& cfg,const Parameters& p,const Geometry& g){
        std::vector<double> result;
        for(int type=0;type<6;++type){
            int count=type==0?g.blackPlaquettes.size():type==1?2*g.blackPlaquettes.size():type==2?g.j2bonds.size():type==3?g.j3bonds.size():p.Ns;
            for(int index=0;index<count;++index)result.push_back(update.diagonal_weight(type,index,cfg,p,g));
        }
        return result;
    }
};

int main(int argc,char** argv){
    if(argc!=3)return 2;
    Parameters p{};p.Lx=std::stoi(argv[1]);p.Ly=std::stoi(argv[2]);p.Lz=1;p.Ns=p.Lx*p.Ly;
    Geometry g;g.build(p.Lx,p.Ly,1);SSEUpdates update;SSEConfig cfg;cfg.spins.spin.resize(p.Ns);
    SSESiteSymmetry symmetry(p,g);
    Geometry broken=g;broken.j2bonds.push_back({0,1,1.});
    bool rejected=false;try{SSESiteSymmetry invalid(p,broken);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)return 3;
    double worst=0.;
    for(int sign1:{-1,1})for(int sign2:{-1,1}){
        p.J0=1.;p.J1=sign1*.3;p.J2=sign2*.047;p.hx=-.23;
        for(int state=0;state<(1<<p.Ns);++state){
            for(int site=0;site<p.Ns;++site)cfg.spins.spin[site]=(state>>site&1)?1:-1;
            auto before=SSEKernelContract::weights(update,cfg,p,g);
            for(int parity=0;parity<2;++parity){
                symmetry.apply_mask(cfg,1u<<parity);
                auto after=SSEKernelContract::weights(update,cfg,p,g);
                for(size_t i=0;i<before.size();++i)worst=std::max(worst,std::abs(before[i]-after[i]));
                symmetry.apply_mask(cfg,1u<<parity);
            }
        }
    }
    std::cout<<worst<<'\n';return worst<1e-13?0:1;
}

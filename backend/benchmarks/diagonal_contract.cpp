// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
// Cross-language contract: compare the production shifted operator sum with
// independently constructed molecular-dipole Hamiltonians for EVERY basis state.
#include "../include/variables.hpp"
#include "../include/updates.hpp"
#include "../include/observables.hpp"
#include "../include/tempering.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
struct SSEKernelContract {
    static double weight(SSEUpdates& u, OpSubtype type, int index, const SSEConfig& c, const Parameters& p, const Geometry& g) {
        return u.diagonal_weight(int(type), index, c, p, g);
    }
};
int main(int argc,char** argv) {
    if(argc!=6)return 2;
    if(std::string(argv[1])=="swap"){
        std::cout<<std::setprecision(17)<<SSETempering::swap_log_ratio(std::stod(argv[2]),std::stod(argv[3]),std::stoi(argv[4]),std::stoi(argv[5]));
        return 0;
    }
    Parameters p;p.Lx=std::stoi(argv[1]);p.Ly=std::stoi(argv[2]);p.Lz=1;
    p.Ns=p.Lx*p.Ly;p.beta=1;p.J0=std::stod(argv[3]);p.J1=std::stod(argv[4]);p.J2=std::stod(argv[5]);p.J3=p.hx=0;
    Geometry g;g.build(p.Lx,p.Ly,1);
    double shift=g.blackPlaquettes.size()*(std::abs(p.J0)+2*std::abs(p.J1))+g.j2_constant_term*p.J2;
    for(const auto& b:g.j2bonds)shift+=std::abs(b.J_val*p.J2);
    SSEUpdates u;SSEConfig c;c.spins.spin.resize(p.Ns);
    std::cout<<std::setprecision(17);
    for(int state=0;state<(1<<p.Ns);++state){
        for(int site=0;site<p.Ns;++site)c.spins.spin[site]=((state>>site)&1)?1:-1;
        double energy=shift;
        for(size_t i=0;i<g.blackPlaquettes.size();++i){
            energy-=SSEKernelContract::weight(u,OpSubtype::J0_Plaquette,i,c,p,g);
            energy-=SSEKernelContract::weight(u,OpSubtype::J1_Plaquette,2*i,c,p,g);
            energy-=SSEKernelContract::weight(u,OpSubtype::J1_Plaquette,2*i+1,c,p,g);
        }
        for(size_t i=0;i<g.j2bonds.size();++i)energy-=SSEKernelContract::weight(u,OpSubtype::J2_Dipole,i,c,p,g);
        std::cout<<energy<<'\n';
    }
}

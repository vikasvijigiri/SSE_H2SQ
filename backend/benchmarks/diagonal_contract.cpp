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
#include <algorithm>
#include <stdexcept>
struct SSEKernelContract {
    static double weight(SSEUpdates& u, OpSubtype type, int index, const SSEConfig& c, const Parameters& p, const Geometry& g) {
        return u.diagonal_weight(int(type), index, c, p, g);
    }
};
int main(int argc,char** argv) {
    if(argc==5 && std::string(argv[1])=="geometry") {
        const int lx=std::stoi(argv[2]),ly=std::stoi(argv[3]),lz=std::stoi(argv[4]);
        Geometry g;g.build(lx,ly,lz);
        const int ns=lx*ly*lz;
        auto xyz=[&](int i){return std::array<int,3>{i%lx,(i/lx)%ly,i/(lx*ly)};};
        if(g.blackPlaquettes.size()!=static_cast<size_t>(ns/2) ||
           g.whitePlaquettes.size()!=static_cast<size_t>(ns/2) ||
           g.j2pairs.size()!=static_cast<size_t>(ns) ||
           g.j3bonds.size()!=static_cast<size_t>(lz>1?ns:0) ||
           g.siteToPlaquette.size()!=static_cast<size_t>(ns)) return 3;
        for(const auto& p:g.blackPlaquettes) {
            const auto c=xyz(p.sites[0]);
            const int x=c[0],y=c[1],z=c[2],off=lx*ly*z;
            const std::array<int,4> expected={y*lx+x+off,((y+1)%ly)*lx+x+off,
                ((y+1)%ly)*lx+(x+1)%lx+off,y*lx+(x+1)%lx+off};
            if(p.sites!=expected || ((x+y)%2)!=0) {
                std::cerr<<"bad plaquette: xyz="<<x<<','<<y<<','<<z<<" sites=";
                for(int s:p.sites)std::cerr<<s<<',';
                std::cerr<<" expected=";
                for(int s:expected)std::cerr<<s<<',';
                std::cerr<<'\n';
                return 4;
            }
        }
        std::vector<int> outgoing(g.blackPlaquettes.size(),0);
        for(const auto& pair:g.j2pairs) {
            if(pair[0]<0 || pair[1]<0 || pair[0]>=static_cast<int>(g.blackPlaquettes.size()) ||
               pair[1]>=static_cast<int>(g.blackPlaquettes.size())) return 5;
            ++outgoing[pair[0]];
            const auto a=xyz(g.blackPlaquettes[pair[0]].sites[0]);
            const auto b=xyz(g.blackPlaquettes[pair[1]].sites[0]);
            const int dx=(b[0]-a[0]+lx)%lx,dy=(b[1]-a[1]+ly)%ly;
            if(a[2]!=b[2] || dy!=1 || (dx!=1 && dx!=lx-1)) return 6;
        }
        if(std::any_of(outgoing.begin(),outgoing.end(),[](int n){return n!=2;})) return 9;
        for(const auto& bond:g.j3bonds) {
            const auto a=xyz(bond[0]),b=xyz(bond[1]);
            if(a[0]!=b[0] || a[1]!=b[1] || b[2]!=(a[2]+1)%lz) return 7;
        }
        for(int site=0;site<ns;++site) {
            const auto p=g.siteToPlaquette[site];
            if(p[0]<0 || p[1]<0 || p[0]==p[1] ||
               p[0]>=static_cast<int>(g.blackPlaquettes.size()) ||
               p[1]>=static_cast<int>(g.blackPlaquettes.size())) return 10;
            for(int q:p) {
                const auto& corners=g.blackPlaquettes[q].sites;
                if(std::find(corners.begin(),corners.end(),site)==corners.end()) return 11;
            }
        }
        for(const auto& bond:g.j2bonds)
            if(bond.s1<0 || bond.s2<0 || bond.s1>=ns || bond.s2>=ns || bond.s1==bond.s2) return 8;
        std::cout<<ns<<' '<<g.blackPlaquettes.size()<<' '<<g.whitePlaquettes.size()<<' '
                 <<g.j2pairs.size()<<' '<<g.j2bonds.size()<<' '<<g.j3bonds.size()<<'\n';
        return 0;
    }
    if(argc==9 && std::string(argv[1])=="stack") {
        Parameters p;p.Lx=std::stoi(argv[2]);p.Ly=std::stoi(argv[3]);p.Lz=std::stoi(argv[4]);
        p.Ns=p.Lx*p.Ly*p.Lz;p.beta=1.;p.J0=std::stod(argv[5]);p.J1=std::stod(argv[6]);
        p.J2=std::stod(argv[7]);p.J3=std::stod(argv[8]);p.hx=0.;
        if(p.Ns>16)return 9;
        Geometry g;g.build(p.Lx,p.Ly,p.Lz);
        double shift=g.blackPlaquettes.size()*(std::abs(p.J0)+2*std::abs(p.J1))+
                     g.j2_constant_term*p.J2+g.j3bonds.size()*std::abs(p.J3);
        for(const auto& b:g.j2bonds)shift+=std::abs(b.J_val*p.J2);
        SSEUpdates u;SSEConfig c;c.spins.spin.resize(p.Ns);
        std::cout<<std::setprecision(17);
        for(int state=0;state<(1<<p.Ns);++state) {
            for(int site=0;site<p.Ns;++site)c.spins.spin[site]=((state>>site)&1)?1:-1;
            double energy=shift;
            for(size_t i=0;i<g.blackPlaquettes.size();++i) {
                energy-=SSEKernelContract::weight(u,OpSubtype::J0_Plaquette,i,c,p,g);
                energy-=SSEKernelContract::weight(u,OpSubtype::J1_Plaquette,2*i,c,p,g);
                energy-=SSEKernelContract::weight(u,OpSubtype::J1_Plaquette,2*i+1,c,p,g);
            }
            for(size_t i=0;i<g.j2bonds.size();++i)energy-=SSEKernelContract::weight(u,OpSubtype::J2_Dipole,i,c,p,g);
            for(size_t i=0;i<g.j3bonds.size();++i)energy-=SSEKernelContract::weight(u,OpSubtype::J3_Inter,i,c,p,g);
            std::cout<<energy<<'\n';
        }
        return 0;
    }
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

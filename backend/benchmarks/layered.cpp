// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
// Small-stack SSE trace for the exact uniform-J3 sign-mapping benchmark.
#include "../include/variables.hpp"
#include "../include/updates.hpp"
#include "../include/observables.hpp"
#include "../include/read_input.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct LayerMeasure {
    double ice=0., rho4=0., plaquette=0.;
    double spin_q0=0., spin_qpi=0.;
    double dipole_q0=0., dipole_qpi=0.;
    std::vector<double> molecular_x;
};

static LayerMeasure measure(const std::vector<int>& spin,
                            const Parameters& p,
                            const Geometry& g) {
    LayerMeasure m;
    const int area=p.Lx*p.Ly;
    const int per_layer=std::max(1,area/2);
    std::vector<double> layer_px(p.Lz,0.), layer_py(p.Lz,0.);
    for(int z=0;z<p.Lz;++z) {
        double staggered=0.;
        for(int y=0;y<p.Ly;++y) for(int x=0;x<p.Lx;++x) {
            const int s=spin[z*area+x+p.Lx*y];
            staggered+=((x+y)%2 ? -1. : 1.)*s;
        }
        const double mz=staggered/double(area);
        m.spin_q0+=mz/double(p.Lz);
        m.spin_qpi+=(z%2 ? -mz : mz)/double(p.Lz);
    }
    for(const auto& q:g.blackPlaquettes) {
        const int z=q.sites[0]/area;
        const int a=spin[q.sites[0]], b=spin[q.sites[1]],
                  c=spin[q.sites[2]], d=spin[q.sites[3]];
        const double px=(a-b-c+d)/4., py=(a+b-c-d)/4.;
        const double eta=((q.sites[0]%p.Lx)%2 ? -1. : 1.);
        const int sum=a+b+c+d;
        const bool two_up_two_down=(sum==0);
        const bool polarized_ice=two_up_two_down && (px*px+py*py>.5);
        m.ice+=two_up_two_down ? 1. : 0.;
        m.rho4+=polarized_ice ? 1. : -1./3.;
        m.plaquette+=double(a*b*c*d);
        layer_px[z]+=eta*px/double(per_layer);
        layer_py[z]+=eta*py/double(per_layer);
        m.molecular_x.push_back(px);
    }
    const double np=double(g.blackPlaquettes.size());
    m.ice/=np; m.rho4/=np; m.plaquette/=np;
    for(int z=0;z<p.Lz;++z) {
        m.dipole_q0+=(layer_px[z]*layer_px[z]+layer_py[z]*layer_py[z]);
        const double sign=(z%2 ? -1. : 1.);
        m.dipole_qpi+=sign*(layer_px[z]*layer_px[z]+layer_py[z]*layer_py[z]);
    }
    // The sign-sensitive structure factors are squares of the layer-summed
    // physical dipole.  Keep both Cartesian components before squaring.
    double p0x=0.,p0y=0.,ppix=0.,ppiy=0.;
    for(int z=0;z<p.Lz;++z) {
        p0x+=layer_px[z]/p.Lz; p0y+=layer_py[z]/p.Lz;
        const double sign=(z%2 ? -1. : 1.);
        ppix+=sign*layer_px[z]/p.Lz; ppiy+=sign*layer_py[z]/p.Lz;
    }
    m.dipole_q0=p0x*p0x+p0y*p0y;
    m.dipole_qpi=ppix*ppix+ppiy*ppiy;
    return m;
}

#ifndef LAYERED_MEASURE_ONLY
int main(int argc,char** argv) {
    if(argc!=5) { std::cerr<<"input trace algorithm seed\n"; return 2; }
    Parameters p; read_parameters(argv[1],p);
    if(p.Lz<2 || p.Lz%2) { std::cerr<<"Use an even periodic stack with Lz >= 2.\n"; return 2; }
    Geometry g; g.build(p.Lx,p.Ly,p.Lz);
    RNG rng(std::stoull(argv[4]));
    SSEConfig cfg; cfg.max_ops=20; cfg.n_ops=0; cfg.op_string.resize(20);
    cfg.spins.spin.resize(p.Ns);
    for(auto& s:cfg.spins.spin) s=rng.uniform()<.5 ? -1 : 1;
    SSEUpdates update;
    const std::string algorithm=argv[3];
    update.spatial_premarking=(algorithm!="vertex");
    update.canonical_enabled=(algorithm=="canonical" || algorithm=="symmetry");
    update.symmetry_enabled=(algorithm=="symmetry");
    update.local_skip_probability=update.symmetry_enabled ? .01 : 0.;
    for(int step=0;step<p.n_therm;++step) {
        update.sweep(cfg,p,g,rng); update.adjust_operator_string(cfg);
    }
    std::ofstream out(argv[2]);
    if(!out) { std::cerr<<"Cannot create trace file.\n"; return 2; }
    out<<"E,E2,ice,rho4,plaquette,spin_q0_sq,spin_qpi_sq,dipole_q0_sq,dipole_qpi_sq\n";
    out<<std::setprecision(17);
    for(int step=0;step<p.n_measure;++step) {
        update.sweep(cfg,p,g,rng);
        Observables energy; energy.accumulate(cfg,p,g);
        const int count=std::max(1,cfg.n_ops);
        std::vector<int> propagated=cfg.spins.spin;
        LayerMeasure state=measure(propagated,p,g);
        double ice_sum=0.,rho_sum=0.,plaq_sum=0.,sq0_sum=0.,sqpi_sum=0.,dq0_sum=0.,dqpi_sum=0.;
        int active=0,segment_begin=0;
        auto add=[&](const LayerMeasure& v,int copies) {
            ice_sum+=copies*v.ice; rho_sum+=copies*v.rho4; plaq_sum+=copies*v.plaquette;
            sq0_sum+=copies*v.spin_q0*v.spin_q0; sqpi_sum+=copies*v.spin_qpi*v.spin_qpi;
            dq0_sum+=copies*v.dipole_q0; dqpi_sum+=copies*v.dipole_qpi;
        };
        for(const auto& op:cfg.op_string) if(op.type!=OpType::Identity) {
            if(op.type==OpType::OffDiagonal) {
                add(state,active+1-segment_begin);
                segment_begin=active+1;
                propagated[op.index]*=-1;
                state=measure(propagated,p,g);
            }
            ++active;
        }
        if(active==0) add(state,1);
        else if(segment_begin<active) add(state,active-segment_begin);
        const double denom=double(count);
        out<<energy.E<<','<<energy.E2<<','<<ice_sum/denom<<','<<rho_sum/denom<<','<<plaq_sum/denom<<','
           <<sq0_sum/denom<<','<<sqpi_sum/denom<<','<<dq0_sum/denom<<','<<dqpi_sum/denom<<'\n';
    }
    return 0;
}
#endif

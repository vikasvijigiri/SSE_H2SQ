// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#define LAYERED_MEASURE_ONLY
#include "layered.cpp"

int main() {
    Parameters p; p.Lx=2; p.Ly=2; p.Lz=2; p.Ns=8; p.Np=4;
    Geometry g; g.build(p.Lx,p.Ly,p.Lz);
    std::cout<<std::setprecision(17);
    for(int basis=0;basis<(1<<p.Ns);++basis) {
        std::vector<int> spins(p.Ns);
        for(int i=0;i<p.Ns;++i) spins[i]=((basis>>i)&1) ? 1 : -1;
        const auto value=measure(spins,p,g);
        std::cout<<basis<<','<<value.ice<<','<<value.rho4<<','<<value.plaquette<<','
                 <<value.spin_q0*value.spin_q0<<','<<value.spin_qpi*value.spin_qpi<<','
                 <<value.dipole_q0<<','<<value.dipole_qpi<<'\n';
    }
    return 0;
}

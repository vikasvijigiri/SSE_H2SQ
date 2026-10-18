#include "../include/variables.hpp"

#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <iostream>
#include <algorithm>
#include <stdexcept>

// --------------------------------------------------
// Utility
// --------------------------------------------------
static inline int mod(int a, int b)
{
    return (a % b + b) % b;
}

// --------------------------------------------------
// Build lattice geometry
// --------------------------------------------------

void Geometry::build(int Lx, int Ly, int Lz)
{
    // The periodic checkerboard coloring used by J0/J1/J2 closes only when
    // both in-plane periods are even. Reject invalid cells before allocating
    // or silently constructing a lattice with missing/miscolored bonds.
    if (Lx < 2 || Ly < 2 || Lz < 1 || (Lx % 2) != 0 || (Ly % 2) != 0) {
        throw std::invalid_argument(
            "Periodic checkerboard geometry requires even Lx,Ly >= 2 and Lz >= 1.");
    }

    blackPlaquettes.clear();
    whitePlaquettes.clear();
    neighbors.clear();
    siteToPlaquette.clear();
    j2pairs.clear();
    j3bonds.clear();
    j2_constant_term = 0.0;

    const int Ns = Lx * Ly * Lz;

    neighbors.resize(Ns);
    siteToPlaquette.assign(Ns, std::array<int, 2>{-1, -1});

    // Grid to map (x,y,z) -> {color, index}
    // Only XY plaquettes for now (8*4=32 logic)
    struct PlaqInfo { int color; int idx; };
    std::vector<std::vector<std::vector<PlaqInfo>>> grid(
        Lz, std::vector<std::vector<PlaqInfo>>(
            Ly, std::vector<PlaqInfo>(Lx, {-1, -1})));

    for (int z = 0; z < Lz; ++z) {
        int z_offset = Lx * Ly * z;

        for (int y = 0; y < Ly; ++y) {
            for (int x = 0; x < Lx; ++x) {
                int i = x + Lx * y + z_offset;

                // Neighbors (cubic, coordination 6)
                neighbors[i].sites = {
                    mod(x + 1, Lx) + Lx * y + z_offset,
                    mod(x - 1, Lx) + Lx * y + z_offset,
                    x + Lx * mod(y + 1, Ly) + z_offset,
                    x + Lx * mod(y - 1, Ly) + z_offset,
                    x + Lx * y + Lx * Ly * mod(z + 1, Lz),
                    x + Lx * y + Lx * Ly * mod(z - 1, Lz)
                };

                // J3 bonds (Vertical)
                if (Lz > 1) j3bonds.push_back({i, (int)neighbors[i].sites[4]});

                // XY Plaquette at (x,y,z)
                Plaquette p;
                p.sites = {
                    i,
                    x + Lx * mod(y + 1, Ly) + z_offset,
                    mod(x + 1, Lx) + Lx * mod(y + 1, Ly) + z_offset,
                    mod(x + 1, Lx) + Lx * y + z_offset
                };

                if ((x + y) % 2 == 0) {
                    grid[z][y][x] = {0, (int)blackPlaquettes.size()};
                    blackPlaquettes.push_back(p);
                } else {
                    grid[z][y][x] = {1, (int)whitePlaquettes.size()};
                    whitePlaquettes.push_back(p);
                }
            }
        }
    }

    // Site -> black plaquettes mapping
    std::vector<int> counts(Ns, 0);
    for (int p_idx = 0; p_idx < (int)blackPlaquettes.size(); ++p_idx) {
        for (int s : blackPlaquettes[p_idx].sites) {
            if (counts[s] < 2) {
                siteToPlaquette[s][counts[s]] = p_idx;
                counts[s]++;
            }
        }
    }

    // J2 pairs (Nearest Neighbor Black-Black Plaquettes)
    // On the Checkerboard, "nearest" black plaquettes are diagonal neighbors.
    // e.g. (x,y) and (x+1,y+1) or (x+1,y-1).
    // We iterate all black plaquettes and find their unique neighbors to avoid double counting.
    // Neighbors of (x,y) are (x+1, y+1), (x+1, y-1), (x-1, y+1), (x-1, y-1).
    // To count each pair once, we can just look "forward" in x?
    // Let's look at (x+1, y+1) and (x+1, y-1).
    
    // J2 pairs (Nearest Neighbor Black-Black Plaquettes)
    // User definition: Interaction with Top Right (TR) and Top Left (TL) neighbors.
    // Sum_P [ P.TR + P.TL ] covers all diagonal bonds exactly once.
    // TR: (x+1, y+1)
    // TL: (x-1, y+1)
    
    // We also store neighbor indices for printing the requested table.
    // Use a temporary structure or just recalculate during print. 
    // Let's recalculate during print to keep Geometry clean, but use strictly this logic for building J2Pairs.
    
    for (int z = 0; z < Lz; ++z) {
        for (int y = 0; y < Ly; ++y) {
            for (int x = 0; x < Lx; ++x) {
                if (grid[z][y][x].color == 0) { // Black
                    int p_idx = grid[z][y][x].idx;
                    
                    // 1. Top Right (x+1, y+1)
                    int nx_tr = mod(x + 1, Lx);
                    int ny_tr = mod(y + 1, Ly);
                    if (grid[z][ny_tr][nx_tr].color == 0) {
                        int tr_idx = grid[z][ny_tr][nx_tr].idx;
                        j2pairs.push_back({p_idx, tr_idx});
                    }

                    // 2. Top Left (x-1, y+1)
                    int nx_tl = mod(x - 1, Lx);
                    int ny_tl = mod(y + 1, Ly);
                    if (grid[z][ny_tl][nx_tl].color == 0) {
                        int tl_idx = grid[z][ny_tl][nx_tl].idx;
                        j2pairs.push_back({p_idx, tl_idx});
                    }
                }
            }
        }
    }
    
    // Sort the periodic bond multiset; retain directional multiplicity.
    std::sort(j2pairs.begin(), j2pairs.end());
    // Keep periodic parallel bonds when a dimension is two.

    // --------------------------------------------------
    // Expand J2 pairs into Bonds
    // P_A . P_B = (PxA PxB + PyA PyB)
    // PxA = sum_i c_xi Si. c_x = {0.25, 0.25, -0.25, -0.25} for {0,1,2,3}
    // PyA = sum_i c_yi Si. c_y = {0.25, -0.25, -0.25, 0.25}
    // --------------------------------------------------
    
    // Temporary map to accumulate J_ij
    // Key: pair<int,int> (s1 < s2), Value: J
    std::map<std::pair<int, int>, double> bond_map;

    for (const auto& pair : j2pairs) {
        const auto& pA = blackPlaquettes[pair[0]];
        const auto& pB = blackPlaquettes[pair[1]];

        // Sites
        int sA[4] = { pA.sites[0], pA.sites[1], pA.sites[2], pA.sites[3] };
        int sB[4] = { pB.sites[0], pB.sites[1], pB.sites[2], pB.sites[3] };

        // Coefficients corrected to match ED definition
        // Px (Columnar/Vertical): (BL + TL) - (TR + BR) -> { +, +, -, - }
        double cx[4] = { 0.25, 0.25, -0.25, -0.25 };
        
        // Py: (BL + BR) - (TL + TR), cyclic ordering.
        double cy[4] = { 0.25, -0.25, -0.25, 0.25 };

        // Iterate all 16 pairs (i in A, j in B)
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                int site1 = sA[i];
                int site2 = sB[j];

                // Effective coupling contribution from this pair of plaquettes
                // Term = J2 * (cx_i*cx_j + cy_i*cy_j) * S_i * S_j
                double term = (cx[i] * cx[j] + cy[i] * cy[j]);
                
                // If sites are the same, this is an identity operator contribution
                // C * S_i * S_i = C * 1 = C
                if (site1 == site2) {
                    j2_constant_term += term;
                    continue;
                }
                
                // Store sorted
                int u = std::min(site1, site2);
                int v = std::max(site1, site2);
                
                bond_map[{u, v}] += term;
            }
        }
    }

    // Convert map to vector
    j2bonds.clear();
    for (const auto& kv : bond_map) {
        if (std::abs(kv.second) > 1e-10) { // Filter zero couplings
            J2Bond b;
            b.s1 = kv.first.first;
            b.s2 = kv.first.second;
            b.J_val = kv.second; // This is the coeff C_kl. Total Hamiltonian H += J2 * C_kl * S_k * S_l
            j2bonds.push_back(b);
        }
    }
}

// ======================================================
// Print black plaquettes
// ======================================================
void Geometry::print_black_plaquettes() const
{
    // Reconstruct grid for printing neighbors (since we don't store the grid permanently)
    // Or just use the same logic as build:
    // Plaquette p is at (x,y) if (x+y)%2==0.
    // We need to map linear index p back to (x,y).
    // This is hard without storing the grid.
    
    // Easier: Store the neighbors explicitly during build or just iterate grid again here.
    // Since we didn't change the struct to store neighbors, we will rebuild the grid locally.
    
    // We need Lx, Ly, Lz from somewhere? 
    // They are not stored in Geometry directly, passed to build.
    // We should assume standard sizes or pass them?
    // Wait, Geometry doesn't know Lx, Ly.
    // But print_black_plaquettes() is called from main where we have prm.
    // The signature is const void. 
    
    // Actually, siteToPlaquette map can help?
    // A plaquette P contains sites (s0, s1, s2, s3).
    // s0 = (x,y). s2 = (x+1, y+1).
    // TR neighbor shares site s2? 
    // TR neighbor of P(x,y) is P(x+1, y+1).
    // P(x+1, y+1) has BL corner at (x+1, y+1) aka s2 of P.
    // So TR neighbor contains site s2 at index 0?
    // Let's use site connectivity.
    
    std::cout << "\n=======================================================\n";
    std::cout << "Black Plaquette Neighbor Table (Reference for J2)\n";
    std::cout << "Plaquette ID | Top Right Neighbor | Top Left Neighbor\n";
    std::cout << "-------------------------------------------------------\n";

    for (int p_idx = 0; p_idx < (int)blackPlaquettes.size(); ++p_idx) {
        const auto& p = blackPlaquettes[p_idx];
        // Sites: s0(BL), s1(TL), s2(TR), s3(BR)
        
        // TR Neighbor: P(x+1, y+1). Its BL corner is (x+1, y+1).
        // This is site p.sites[2] !
        // So the TR neighbor is the black plaquette whose site 0 is p.sites[2].
        // Or strictly, site 0? Not necessarily.
        // It's the black plaquette containing site p.sites[2].
        // Wait, p.sites[2] is shared by 4 plaquettes (2 black, 2 white).
        // P(x,y) [Black], P(x+1,y) [White], P(x,y+1) [White], P(x+1,y+1) [Black].
        // So p.sites[2] is site (x+1,y+1).
        // It is BL corner of P(x+1,y+1).
        // It is TL corner of P(x+1,y).
        // It is BR corner of P(x,y+1).
        // It is TR corner of P(x,y).
        
        // So we look for a black plaquette Q != P that contains p.sites[2].
        int tr_neighbor = -1;
        int s_tr = p.sites[2];
        
        // Use siteToPlaquette lookup
        if (s_tr < (int)siteToPlaquette.size()) {
             // siteToPlaquette stores up to 2 black plaquettes per site.
             // One is P. The other must be TR (or BL of neighbor).
             int cand1 = siteToPlaquette[s_tr][0];
             int cand2 = siteToPlaquette[s_tr][1];
             if (cand1 != -1 && cand1 != p_idx) tr_neighbor = cand1;
             else if (cand2 != -1 && cand2 != p_idx) tr_neighbor = cand2;
        }

        // TL Neighbor: P(x-1, y+1).
        // Sites of P(x,y): (x,y){0} ...
        // P(x-1, y+1) has BR corner at (x, y+1).
        // Site (x, y+1) is P.sites[1] (TL corner of P).
        // So TL neighbor shares site P.sites[1].
        int tl_neighbor = -1;
        int s_tl = p.sites[1];
        
        if (s_tl < (int)siteToPlaquette.size()) {
             int cand1 = siteToPlaquette[s_tl][0];
             int cand2 = siteToPlaquette[s_tl][1];
             if (cand1 != -1 && cand1 != p_idx) tl_neighbor = cand1;
             else if (cand2 != -1 && cand2 != p_idx) tl_neighbor = cand2;
        }

        std::cout << p_idx << " | " << tr_neighbor << " | " << tl_neighbor << "\n";
    }
    std::cout << "=======================================================\n";
    std::cout << std::flush;
}

void Geometry::print_j2_bonds() const
{
    std::cout << "\n------------------------------------------------------\n";
    std::cout << "J2 Bond Decomposition\n";
    std::cout << "------------------------------------------------------\n";
    std::cout << "Total Bonds: " << j2bonds.size() << "\n";
    std::cout << "Constant Shift (from identity terms): " << j2_constant_term << "\n\n";
    
    std::cout << "Format: Bond Index : S1 - S2 : Coefficient\n";
    
    int idx = 0;
    for (const auto& bond : j2bonds) {
        std::cout << "Bond " << idx++ << " : " 
                  << bond.s1 << " - " << bond.s2 << " : " 
                  << bond.J_val << "\n";
    }
    


    std::cout << "------------------------------------------------------\n\n";
    std::cout << std::flush;
}

void Parameters::print_parameters()
{
    std::cout << "\n================ Simulation Parameters ================\n";

    std::cout << "Lattice:\n";
    std::cout << "  Lx            = " << Lx << "\n";
    std::cout << "  Ly            = " << Ly << "\n";
    std::cout << "  Lz            = " << Lz << "\n";

    std::cout << "\nTemperature:\n";
    std::cout << "  beta          = " << beta << "\n";

    std::cout << "\nCouplings:\n";
    std::cout << "  J0 (plaquette)= " << J0 << "\n";
    std::cout << "  J1            = " << J1 << "\n";
    std::cout << "  J2            = " << J2 << "\n";
    std::cout << "  J3            = " << J3 << "\n";
    std::cout << "  hx            = " << hx << "\n";

    std::cout << "\nMonte Carlo:\n";
    std::cout << "  n_therm       = " << n_therm << "\n";
    std::cout << "  n_measure     = " << n_measure << "\n";

    std::cout << "=======================================================\n\n";
}

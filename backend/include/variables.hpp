#ifndef SSE_VARIABLES_HPP
#define SSE_VARIABLES_HPP

#include <vector>
#include <array>
#include <random>
#include <chrono>
#include <string>

// ======================================================
// Geometry data (NO Monte Carlo logic here)
// ======================================================

struct Plaquette {
    std::array<int, 4> sites;   // 4 sites forming plaquette
};

struct SiteNeighbors {
    std::array<int, 6> sites;   // coordination = 6 (your lattice)
};

// ======================================================
// Spin configuration (σ^z basis)
// ======================================================

struct SpinConfig {
    std::vector<int> spin;   // +1 / -1

    inline int operator[](int i) const { return spin[i]; }
    inline int& operator[](int i) { return spin[i]; }

    inline void flip(int i) { spin[i] *= -1; }
};

// ======================================================
// SSE operator definitions
// ======================================================

enum class OpType : int {
    Identity    = 0,
    Diagonal    = 1,
    OffDiagonal = 2
};

// Operator subtype identifiers
enum class OpSubtype : int {
    None = -1,
    J0_Plaquette = 0,
    J1_Plaquette = 1,
    J2_Dipole    = 2,
    J3_Inter     = 3,
    TransverseX_diag = 4,
    TransverseX_offdiag = 5
};


struct SSEOperator {
    OpType type;
    OpSubtype subtype;
    int index;     // site / plaquette / bond index

    SSEOperator()
        : type(OpType::Identity), subtype(OpSubtype::J0_Plaquette), index(-1) {}
};

// ======================================================
// Complete SSE configuration
// ======================================================

struct SSEConfig {
    SpinConfig spins;                    // |α⟩ state
    std::vector<SSEOperator> op_string;  // operator sequence
    int n_ops;                           // expansion order
    int max_ops;                         // allocated length
};

// ======================================================
// Simulation parameters
// ======================================================

struct Parameters {
    int Lx, Ly, Lz;
    int Ns;      // number of sites
    int Np;      // number of plaquettes

    double beta;
    double T;

    // Couplings
    double J0, J1, J2, J3;
    double hx;

    // MC control
    int n_therm;
    int n_measure;

    
    void print_parameters();
};

// ======================================================
// Lattice & connectivity container
// ======================================================

struct Geometry {
    std::vector<Plaquette> blackPlaquettes;
    std::vector<Plaquette> whitePlaquettes;
    std::vector<SiteNeighbors> neighbors;

    std::vector<std::array<int,2>> siteToPlaquette;

    // New for 3D and J2/J3
    std::vector<std::array<int, 2>> j2pairs; // Pairs of interacting plaquettes (index in black, index in white)
    
    struct J2Bond {
        int s1, s2;
        double J_val; // effective coupling magnitude (signed)
    };
    std::vector<J2Bond> j2bonds;
    double j2_constant_term = 0.0; // Sum of identity terms from expansion

    std::vector<std::array<int, 2>> j3bonds; // Interlayer bonds (site_L, site_{L+1})

    void print_black_plaquettes() const;
    void print_j2_bonds() const;

    void build(int Lx, int Ly, int Lz);

};

// ======================================================
// Random number generator
// ======================================================

class RNG {
private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dist;

public:
    explicit RNG(unsigned long long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count())
    : gen(seed),
      dist(0.0, 1.0) {}

    inline double uniform() { return dist(gen); }
};

#endif

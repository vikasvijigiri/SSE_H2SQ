#ifndef SSE_UPDATES_HPP
#define SSE_UPDATES_HPP

#include <vector>
#include <optional>
#include "site_symmetry.hpp"
#include "variables.hpp"

// ======================================================
// SSE Update Engine
// ======================================================

class SSEUpdates {
    friend struct SSEKernelContract;
public:
    SSEUpdates() = default;
    bool spatial_premarking = true;
    bool canonical_enabled = true;
    bool symmetry_enabled = true;
    double local_skip_probability = .01;
    long long diagonal_skips=0,canonical_skips=0;
    double symmetry_seconds = 0.;
    std::optional<SSESiteSymmetry> site_symmetry;
    void reset_symmetry_measurements(){symmetry_seconds=0.;diagonal_skips=canonical_skips=0;if(site_symmetry)site_symmetry->orbit_counts.fill(0);}
    double diagonal_seconds = 0, cluster_seconds = 0, gauge_seconds = 0, canonical_seconds = 0;
    long long canonical_attempts = 0, canonical_accepted = 0;
    void canonical_update(SSEConfig&, const Parameters&, const Geometry&, RNG&);

    // --------------------------------------------------
    // Top-level Monte Carlo sweep
    // --------------------------------------------------
    void sweep(SSEConfig& cfg,
               const Parameters& prm,
               const Geometry& geom,
               RNG& rng);

    // --------------------------------------------------
    // Core SSE updates
    // --------------------------------------------------
    void adjust_operator_string(SSEConfig& cfg);

    void check_spin_consistency(const SSEConfig& cfg,
                                const Geometry& geom);

    // Diagonal insertion / removal
    void diagonal_update(SSEConfig& cfg,
                         const Parameters& prm,
                         const Geometry& geom,
                         RNG& rng);

    // Quantum loop-cluster update (Melko-style SSE BFS cluster algorithm).
    //
    // For each non-identity operator vertex in the SSE string, internal
    // "partner" links encode which pairs of legs must flip together to
    // preserve that vertex's weight:
    //
    //   J0 plaquette  : randomly decomposes σ1σ2σ3σ4 into two 2-body
    //                   bonds chosen from the three pairings
    //                   {(0,1)&(2,3)}, {(0,2)&(1,3)}, {(0,3)&(1,2)}.
    //                   Each pair forms an independent bond in the cluster
    //                   graph — the two bonds CAN end up in different
    //                   clusters.
    //   J1 bond       : independent two-site operators for (0,2) and (1,3).
    //                   Together these represent J1 = σ0σ2 + σ1σ3 structure.
    //   J2 / J3 bond  : both sites in the pair are bonded (flip together).
    //   hx (transv.)  : no internal bond; instead, when a cluster boundary
    //                   straddles an hx vertex the operator toggles between
    //                   hx_diag ↔ hx_offdiag to maintain spin consistency.
    //
    // BFS identifies connected components; each cluster is flipped with
    // probability 1/2.
    void loop_cluster_update(SSEConfig& cfg,
                             const Parameters& prm,
                             const Geometry& geom,
                             RNG& rng);

    // --------------------------------------------------
    // Optional: Self-learning Monte Carlo (global move)
    // --------------------------------------------------
    void slmc_update(SSEConfig& cfg,
                     const Parameters& prm,
                     const Geometry& geom,
                     RNG& rng);

    // Supplementary pure-J0 row/column moves; these alone do not span all sectors.
    // Only active for J0-only model (J1=J2=J3=hx=0).
    void gauge_update(SSEConfig& cfg,
                      const Parameters& prm,
                      const Geometry& geom,
                      RNG& rng);

private:
    // ==================================================
    // Internal helpers
    // ==================================================

    // Build SSE imaginary-time vertex list and partner/time links
    void build_vertex_list(SSEConfig& cfg,
                           const Geometry& geom,
                           RNG& rng);

    void check_integrity();

    // Propagate t=0 boundary spins from vertex_spin after cluster update
    void propagate_spins(SSEConfig& cfg);

    // Diagonal operator weights
    double diagonal_weight(int subtype,
                           int index,
                           const SSEConfig& cfg,
                           const Parameters& prm,
                           const Geometry& geom);

    // ==================================================
    // Internal storage (reused every sweep — no per-sweep allocation)
    // ==================================================

    std::vector<int> vertex_op_index;       // indices of non-identity ops
    std::vector<int> vertex_decomp;         // J0 pair decomposition type per op
    std::vector<int> linked_leg;            // imaginary-time links between legs
    std::vector<int> linked_leg_to_op_index;
    std::vector<int> vertex_partner;        // within-vertex pair-bond links
    std::vector<int> vertex_spin;           // spin value on each leg
    std::vector<int> leg_site;
    std::vector<int> base_vertex_leg;
    std::vector<int> first_leg;             // per site: first leg in imaginary time
    std::vector<int> last_leg;              // per site: last leg in imaginary time

    int n_legs = 0;

    // BFS scratch (member to avoid per-call allocation)
    std::vector<int>  m_cluster_id;
    std::vector<bool> m_cluster_flipped;
    std::vector<int>  m_bfs_queue;

    // Per-vertex base leg index (parallel to vertex_op_index)
    std::vector<int> vertex_base_vec;

    // Cached active subtypes — rebuilt each sweep
    std::vector<OpSubtype> m_active_subtypes;
    std::vector<int>       m_n_diags;
    int                    m_N_active = 0;
    int                    m_n_diag_by_sub[6] = {};

    bool use_slmc = false;
};

#endif // SSE_UPDATES_HPP

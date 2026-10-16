#include "../include/updates.hpp"
#include "../include/variables.hpp"
#include "../include/sse_debug.hpp"  // no-op unless SSE_DEBUG is defined

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>

// ======================================================
// Top-level sweep
// ======================================================
void SSEUpdates::sweep(SSEConfig &cfg, const Parameters &prm,
                       const Geometry &geom, RNG &rng) {
  if(symmetry_enabled){
    if(!site_symmetry || !site_symmetry->matches(prm))site_symmetry.emplace(prm,geom);
    // This fixed model-based floor is set before sampling; padding is then
    // allowed to grow during thermalization only. Production keeps M fixed.
    site_symmetry->ensure_cutoff(cfg);
  }
  auto start = std::chrono::steady_clock::now();
  diagonal_update(cfg, prm, geom, rng);
  auto diagonal_end = std::chrono::steady_clock::now();
#ifdef SSE_DEBUG
  dbg_print_op_string(cfg, prm, geom);
#endif
  for (int i = 0; i < 2; ++i) loop_cluster_update(cfg, prm, geom, rng);
  auto cluster_end = std::chrono::steady_clock::now();
  gauge_update(cfg, prm, geom, rng);
  auto gauge_end = std::chrono::steady_clock::now();
  if (canonical_enabled) canonical_update(cfg, prm, geom, rng);
  auto canonical_end = std::chrono::steady_clock::now();
  if(symmetry_enabled){
    if(!site_symmetry || !site_symmetry->matches(prm))site_symmetry.emplace(prm,geom);
    site_symmetry->sweep(cfg,rng);
  }
  auto end = std::chrono::steady_clock::now();
  diagonal_seconds += std::chrono::duration<double>(diagonal_end-start).count();
  cluster_seconds += std::chrono::duration<double>(cluster_end-diagonal_end).count();
  gauge_seconds += std::chrono::duration<double>(gauge_end-cluster_end).count();
  canonical_seconds += std::chrono::duration<double>(canonical_end-gauge_end).count();
  if(symmetry_enabled)symmetry_seconds += std::chrono::duration<double>(end-canonical_end).count();
  if (use_slmc) slmc_update(cfg, prm, geom, rng);
}

// ======================================================
// 1. Diagonal update
// ======================================================
void SSEUpdates::diagonal_update(SSEConfig &cfg, const Parameters &prm,
                                 const Geometry &geom, RNG &rng) {
  // Rebuild active-subtype cache — member vectors avoid per-sweep allocation
  m_active_subtypes.clear();
  m_n_diags.clear();
  std::fill(std::begin(m_n_diag_by_sub), std::end(m_n_diag_by_sub), 0);

  auto add_type = [&](OpSubtype sub, int nd) {
    if (nd <= 0) return;
    m_active_subtypes.push_back(sub);
    m_n_diags.push_back(nd);
    m_n_diag_by_sub[static_cast<int>(sub)] = nd;
  };

  if (prm.J0 != 0.0)
    add_type(OpSubtype::J0_Plaquette,
             static_cast<int>(geom.blackPlaquettes.size()));
  if (prm.J1 != 0.0)
    add_type(OpSubtype::J1_Plaquette,
             2 * static_cast<int>(geom.blackPlaquettes.size()));
  if (prm.J2 != 0.0)
    add_type(OpSubtype::J2_Dipole, static_cast<int>(geom.j2bonds.size()));
  if (prm.J3 != 0.0)
    add_type(OpSubtype::J3_Inter, static_cast<int>(geom.j3bonds.size()));
  if (prm.hx != 0.0)
    add_type(OpSubtype::TransverseX_diag, prm.Ns);

  m_N_active = static_cast<int>(m_active_subtypes.size());

  for (int p = 0; p < cfg.max_ops; ++p) {
    SSEOperator &op = cfg.op_string[p];

    if (op.type == OpType::Identity) {
      if(local_skip_probability>0 && rng.uniform()<local_skip_probability){++diagonal_skips;continue;}
      if (m_N_active > 0) {
        int type_idx = int(rng.uniform() * m_N_active);
        OpSubtype sub = m_active_subtypes[type_idx];
        int n_diag = m_n_diags[type_idx];
        int index = int(rng.uniform() * n_diag);

        double w =
            diagonal_weight(static_cast<int>(sub), index, cfg, prm, geom);
        if (w > 0.0) {
          double prob = m_N_active * prm.beta * n_diag * w /
                        double(cfg.max_ops - cfg.n_ops);
          if (rng.uniform() < prob) {
            op.type = OpType::Diagonal;
            op.subtype = sub;
            op.index = index;
            ++cfg.n_ops;
          }
        }
      }
    } else if (op.type == OpType::Diagonal) {
      if(local_skip_probability>0 && rng.uniform()<local_skip_probability){++diagonal_skips;continue;}
      OpSubtype sub = op.subtype;
      int n_diag = m_n_diag_by_sub[static_cast<int>(sub)]; // O(1) array lookup

      double w =
          diagonal_weight(static_cast<int>(sub), op.index, cfg, prm, geom);
      if (w > 0.0) {
        double prob = double(cfg.max_ops - cfg.n_ops + 1) /
                      (m_N_active * prm.beta * n_diag * w);
        if (rng.uniform() < prob) {
          op.type = OpType::Identity;
          op.subtype = OpSubtype::None;
          op.index = -1;
          --cfg.n_ops;
        }
      }
    } else if (op.subtype == OpSubtype::TransverseX_offdiag) {
      cfg.spins.spin[op.index] *= -1;
    }
  }
}

// ======================================================
// 3. Diagonal operator weight
// ======================================================
double SSEUpdates::diagonal_weight(int subtype, int index, const SSEConfig &cfg,
                                   const Parameters &prm,
                                   const Geometry &geom) {
  switch (static_cast<OpSubtype>(subtype)) {

  case OpSubtype::J0_Plaquette: {
    const auto &p = geom.blackPlaquettes[index];
    int prod = 1;
    for (int s : p.sites)
      prod *= cfg.spins.spin[s];
    return std::abs(prm.J0) + prm.J0 * prod;
  }

  case OpSubtype::J1_Plaquette: {
    const auto &p = geom.blackPlaquettes[index / 2];
    int a = (index % 2 == 0) ? 0 : 1;
    int b = a + 2;
    return std::abs(prm.J1) - prm.J1 * cfg.spins.spin[p.sites[a]] * cfg.spins.spin[p.sites[b]];
  }

  case OpSubtype::J2_Dipole: {
    const auto &bond = geom.j2bonds[index];
    double J_local = bond.J_val * prm.J2;
    return std::abs(J_local) -
           J_local * cfg.spins.spin[bond.s1] * cfg.spins.spin[bond.s2];
  }

  case OpSubtype::J3_Inter: {
    const auto &bond = geom.j3bonds[index];
    return std::abs(prm.J3) - prm.J3 * cfg.spins.spin[bond[0]] * cfg.spins.spin[bond[1]];
  }

  case OpSubtype::TransverseX_diag:
    return std::abs(prm.hx);

  default:
    return 0.0;
  }
}

// ======================================================
// 2. Quantum loop-cluster update (Melko-style BFS cluster algorithm)
// ======================================================
void SSEUpdates::loop_cluster_update(SSEConfig &cfg, const Parameters &/*prm*/,
                                      const Geometry &geom, RNG &rng) {
  build_vertex_list(cfg, geom, rng);
#ifndef NDEBUG
  check_integrity();
#endif

#ifdef SSE_DEBUG
  static int dbg_pass = 0;
  ++dbg_pass;
  constexpr int DBG_MAX_CLUSTERS = 2;
  dbg_loop_pass_hdr(dbg_pass, n_legs);
  dbg_print_vertex_list(
      vertex_op_index, linked_leg, vertex_partner, vertex_spin,
      leg_site, vertex_base_vec, vertex_decomp,
      cfg.op_string, first_leg, last_leg,
      n_legs, (int)cfg.spins.spin.size(), dbg_pass);
#endif

  // --- 1. Identify clusters via BFS (member scratch — no per-call allocation)
  // ---
  m_cluster_id.assign(n_legs, -1);
  m_cluster_flipped.clear();
  m_cluster_flipped.reserve(n_legs);
  int n_clusters = 0;

#ifdef SSE_DEBUG
  std::vector<int> dbg_spin_before; // saved just before spin update loop
#endif

  for (int i = 0; i < n_legs; ++i) {
    if (m_cluster_id[i] != -1)
      continue;

    const int cid = n_clusters++;

#ifdef SSE_DEBUG
    const double dbg_r = rng.uniform();
    const bool do_flip = (dbg_r < 0.5);
    if (cid < DBG_MAX_CLUSTERS)
      dbg_loop_cluster_start(cid, i, dbg_r, do_flip, leg_site);
#else
    const bool do_flip = (rng.uniform() < 0.5);
#endif

    m_cluster_flipped.push_back(do_flip);

    m_cluster_id[i] = cid;
    m_bfs_queue.clear();
    m_bfs_queue.push_back(i);

    for (size_t head = 0; head < m_bfs_queue.size(); ++head) {
      const int u = m_bfs_queue[head];

      const int v1 = linked_leg[u];
      if (v1 != -1 && m_cluster_id[v1] == -1) {
        m_cluster_id[v1] = cid;
        m_bfs_queue.push_back(v1);
#ifdef SSE_DEBUG
        if (cid < DBG_MAX_CLUSTERS)
          dbg_loop_bfs_expand(u, v1, "linked ", leg_site);
#endif
      }
      const int v2 = vertex_partner[u];
      if (v2 != -1 && m_cluster_id[v2] == -1) {
        m_cluster_id[v2] = cid;
        m_bfs_queue.push_back(v2);
#ifdef SSE_DEBUG
        if (cid < DBG_MAX_CLUSTERS)
          dbg_loop_bfs_expand(u, v2, "partner", leg_site);
#endif
      }
    }

#ifdef SSE_DEBUG
    if (cid < DBG_MAX_CLUSTERS)
      dbg_loop_cluster_end(cid, (int)m_bfs_queue.size());
    else if (cid == DBG_MAX_CLUSTERS)
      std::cout << PB << "  " << AGY << ADIM
                << "  ... remaining clusters not shown (only first "
                << DBG_MAX_CLUSTERS << " printed) ..." << A0 << "\n";
#endif
  }

  // --- 2. Update t=0 boundary spins ---
  const int N = static_cast<int>(cfg.spins.spin.size());

#ifdef SSE_DEBUG
  dbg_spin_before = cfg.spins.spin;
#endif

  for (int s = 0; s < N; ++s) {
    if (first_leg[s] != -1) {
      if (m_cluster_flipped[m_cluster_id[first_leg[s]]])
        cfg.spins.spin[s] *= -1;
    } else {
      if (rng.uniform() < 0.5)
        cfg.spins.spin[s] *= -1;
    }
  }

#ifdef SSE_DEBUG
  dbg_final_spins(dbg_spin_before, cfg.spins.spin);
#endif

  // --- 3. Toggle TransverseX operators that straddle a cluster boundary ---
  // vertex_base_vec[vp] pre-computed in build_vertex_list: no base-tracking
  // needed here.
#ifdef SSE_DEBUG
  bool dbg_had_toggle = false;
#endif

  for (int vp = 0; vp < static_cast<int>(vertex_op_index.size()); ++vp) {
    const int p = vertex_op_index[vp];
    auto &op = cfg.op_string[p];

    if (op.subtype != OpSubtype::TransverseX_diag &&
        op.subtype != OpSubtype::TransverseX_offdiag)
      continue;

    const int base = vertex_base_vec[vp];
    const bool flip_in = m_cluster_flipped[m_cluster_id[base]];
    const bool flip_out = m_cluster_flipped[m_cluster_id[base + 1]];

    if (flip_in != flip_out) {
#ifdef SSE_DEBUG
      if (!dbg_had_toggle) { dbg_loop_toggle_section(); dbg_had_toggle = true; }
      const char* from_s = (op.subtype == OpSubtype::TransverseX_diag) ? "hx_diag " : "hx_offD ";
      const char* to_s   = (op.subtype == OpSubtype::TransverseX_diag) ? "hx_offD " : "hx_diag ";
      dbg_loop_op_toggle(p, op.index, from_s, to_s);
#endif
      if (op.subtype == OpSubtype::TransverseX_diag) {
        op.subtype = OpSubtype::TransverseX_offdiag;
        op.type = OpType::OffDiagonal;
      } else {
        op.subtype = OpSubtype::TransverseX_diag;
        op.type = OpType::Diagonal;
      }
    }
  }

#ifdef SSE_DEBUG
  {
    int dbg_n_flipped = 0;
    for (bool f : m_cluster_flipped) if (f) ++dbg_n_flipped;
    dbg_loop_footer(n_clusters, dbg_n_flipped);
  }
#endif

#ifndef NDEBUG
  check_spin_consistency(cfg, geom);
#endif
}

// ======================================================
// 4. Build vertex list (SSE imaginary-time linking)
// ======================================================
void SSEUpdates::build_vertex_list(SSEConfig &cfg, const Geometry &geom,
                                   RNG &rng) {
  vertex_op_index.clear();
  linked_leg.clear();
  vertex_spin.clear();
  leg_site.clear();
  linked_leg_to_op_index.clear();
  base_vertex_leg.clear();
  vertex_partner.clear();
  vertex_base_vec.clear();

  const int N = static_cast<int>(cfg.spins.spin.size());
  const int M = cfg.max_ops;

  for (int p = 0; p < M; ++p)
    if (cfg.op_string[p].type != OpType::Identity)
      vertex_op_index.push_back(p);

  const int max_legs = cfg.n_ops * 8;
  linked_leg.assign(max_legs, -1);
  vertex_partner.assign(max_legs, -1);
  vertex_spin.resize(max_legs);
  leg_site.resize(max_legs);
  linked_leg_to_op_index.resize(max_legs);
  base_vertex_leg.resize(max_legs);

  last_leg.assign(N, -1);
  first_leg.assign(N, -1);

  vertex_base_vec.reserve(vertex_op_index.size());
  std::vector<int> spatial_motifs(geom.blackPlaquettes.size());
  for (auto& motif : spatial_motifs) motif = int(rng.uniform() * 3.0);

  std::vector<int> spin = cfg.spins.spin;

  int leg_counter = 0;

  if (static_cast<int>(vertex_decomp.size()) < M)
    vertex_decomp.resize(M);

  for (int vp = 0; vp < static_cast<int>(vertex_op_index.size()); ++vp) {
    const int p = vertex_op_index[vp];
    const auto &op = cfg.op_string[p];
    int decomp_type = 0;

    // Stack-allocated buffer — avoids one heap allocation per vertex in the hot
    // loop
    int sites_buf[4];
    int n_sites = 0;
    int n_op_legs = 0;

    if (op.subtype == OpSubtype::TransverseX_diag ||
        op.subtype == OpSubtype::TransverseX_offdiag) {
      sites_buf[0] = op.index;
      n_sites = 1;
      n_op_legs = 2;
    } else if (op.subtype == OpSubtype::J0_Plaquette) {
      const auto &plaq = geom.blackPlaquettes[op.index];
      sites_buf[0] = plaq.sites[0];
      sites_buf[1] = plaq.sites[1];
      sites_buf[2] = plaq.sites[2];
      sites_buf[3] = plaq.sites[3];
      n_sites = 4;
      n_op_legs = 8;
      decomp_type = spatial_premarking ? spatial_motifs[op.index] : int(rng.uniform() * 3.0);
    } else if (op.subtype == OpSubtype::J1_Plaquette) {
      const auto &plaq = geom.blackPlaquettes[op.index / 2];
      int a = (op.index % 2 == 0) ? 0 : 1;
      sites_buf[0] = plaq.sites[a];
      sites_buf[1] = plaq.sites[a + 2];
      n_sites = 2;
      n_op_legs = 4;
    } else if (op.subtype == OpSubtype::J2_Dipole) {
      const auto &bond = geom.j2bonds[op.index];
      sites_buf[0] = bond.s1;
      sites_buf[1] = bond.s2;
      n_sites = 2;
      n_op_legs = 4;
    } else if (op.subtype == OpSubtype::J3_Inter) {
      const auto &bond = geom.j3bonds[op.index];
      sites_buf[0] = bond[0];
      sites_buf[1] = bond[1];
      n_sites = 2;
      n_op_legs = 4;
    }

    vertex_decomp[p] = decomp_type;
    vertex_base_vec.push_back(
        leg_counter); // store base for directed_loop_update

    if (n_op_legs > 0) {
      const int base = leg_counter;
      leg_counter += n_op_legs;

      if (op.subtype == OpSubtype::TransverseX_diag ||
          op.subtype == OpSubtype::TransverseX_offdiag) {

        const int site = sites_buf[0];
        const int in_leg = base;
        const int out_leg = base + 1;

        vertex_spin[in_leg]  = spin[site];
        // Out-leg: spin flips for off-diagonal; propagate through imaginary time
        vertex_spin[out_leg] = (op.subtype == OpSubtype::TransverseX_offdiag)
                               ? -spin[site] : spin[site];
        if (op.subtype == OpSubtype::TransverseX_offdiag) spin[site] *= -1;

        leg_site[in_leg] = site;
        leg_site[out_leg] = site;

        linked_leg_to_op_index[in_leg] = p;
        linked_leg_to_op_index[out_leg] = p;

        base_vertex_leg[in_leg] = base;
        base_vertex_leg[out_leg] = base;

        if (last_leg[site] != -1) {
          linked_leg[last_leg[site]] = in_leg;
          linked_leg[in_leg] = last_leg[site];
        } else {
          first_leg[site] = in_leg;
        }
        last_leg[site] = out_leg;
      } else {
        // J0, J1, J2, J3: time-link all site legs
        for (int i = 0; i < n_sites; ++i) {
          const int site = sites_buf[i];
          const int l_in = base + 2 * i;
          const int l_out = l_in + 1;

          // Populate vertex_spin for all diagonal-op legs (in == out for diagonal)
          vertex_spin[l_in]  = spin[site];
          vertex_spin[l_out] = spin[site];

          leg_site[l_in] = site;
          leg_site[l_out] = site;

          linked_leg_to_op_index[l_in] = p;
          linked_leg_to_op_index[l_out] = p;

          base_vertex_leg[l_in] = base;
          base_vertex_leg[l_out] = base;

          if (last_leg[site] != -1) {
            linked_leg[last_leg[site]] = l_in;
            linked_leg[l_in] = last_leg[site];
          } else {
            first_leg[site] = l_in;
          }
          last_leg[site] = l_out;
        }

        // Internal vertex partner links (4-cycle per pair).
        //
        // For each paired pair of sites (A, B) in a diagonal vertex, the
        // 4 legs {in_A, in_B, out_B, out_A} form a directed cycle:
        //   in_A → in_B → out_B → out_A → in_A
        //
        // This forces all 4 legs into the same BFS cluster, which guarantees:
        //  (a) Both sites A and B are always flipped together (preserves the
        //      diagonal coupling σ_A σ_B).
        //  (b) Same-site in_A and out_A are always co-clustered, so the
        //      imaginary-time chain for hx operators passing through this
        //      vertex is never "split" — critical for hx parity correctness.
        //
        // vertex_partner stores a single "next-in-cycle" pointer per leg.
        // BFS follows all partner edges transitively, so the full 4-leg set
        // is visited in one connected component regardless of entry point.
        auto set_pair_cycle = [&](int liA, int liB) {
          // 4-cycle: liA(in_A) → liB(in_B) → liB+1(out_B) → liA+1(out_A) → liA
          vertex_partner[liA]     = liB;
          vertex_partner[liB]     = liB + 1;
          vertex_partner[liB + 1] = liA + 1;
          vertex_partner[liA + 1] = liA;
        };

        if (op.subtype == OpSubtype::J0_Plaquette) {
          const int type = vertex_decomp[p];
          int p1A, p1B, p2A, p2B;
          if (type == 0) {
            p1A = 0; p1B = 1; p2A = 2; p2B = 3;
          } else if (type == 1) {
            p1A = 0; p1B = 2; p2A = 1; p2B = 3;
          } else {
            p1A = 0; p1B = 3; p2A = 1; p2B = 2;
          }
          set_pair_cycle(base + 2 * p1A, base + 2 * p1B);
          set_pair_cycle(base + 2 * p2A, base + 2 * p2B);
        } else if (op.subtype == OpSubtype::J1_Plaquette) {
          // Fixed pairing: sites (0,2) and (1,3)
          set_pair_cycle(base, base + 2);
        } else {
          // J2/J3: two-site bond — pair sites 0 and 1
          set_pair_cycle(base + 0, base + 2);
        }
      }
    }
  }

  // Trim to actual usage
  linked_leg.resize(leg_counter);
  vertex_partner.resize(leg_counter);
  vertex_spin.resize(leg_counter);
  leg_site.resize(leg_counter);
  linked_leg_to_op_index.resize(leg_counter);
  base_vertex_leg.resize(leg_counter);

  // Close periodic boundary conditions
  for (int s = 0; s < N; ++s) {
    if (first_leg[s] != -1) {
      linked_leg[last_leg[s]] = first_leg[s];
      linked_leg[first_leg[s]] = last_leg[s];
    }
  }
  n_legs = leg_counter;
}

// Check linked list basic integrity
void SSEUpdates::check_integrity() {
  const int n = static_cast<int>(linked_leg.size());
  for (int l = 0; l < n; ++l) {
    int next = linked_leg[l];
    if (next < 0 || next >= n) {
      std::cerr << "Integrity Error: Leg " << l << " links to invalid " << next
                << "\n";
      std::abort();
    }
    if (linked_leg[next] != l) {
      std::cerr << "Integrity Error: Asymmetry at leg " << l << " -> " << next
                << " -> " << linked_leg[next] << "\n";
      std::abort();
    }
  }
}

// =====================================================
// Propagate spins after cluster update
// ======================================================
void SSEUpdates::propagate_spins(SSEConfig &cfg) {
  const int N = static_cast<int>(cfg.spins.spin.size());
  for (int site = 0; site < N; ++site) {
    const int leg0 = first_leg[site];
    if (leg0 == -1)
      continue;
    cfg.spins.spin[site] = vertex_spin[leg0];
  }
}

void SSEUpdates::check_spin_consistency(const SSEConfig &cfg,
                                        const Geometry & /*geom*/) {
  std::vector<int> current_spins = cfg.spins.spin;
  const std::vector<int> &initial_spins = cfg.spins.spin;

  for (int p = 0; p < cfg.max_ops; ++p) {
    const auto &op = cfg.op_string[p];
    if (op.type == OpType::OffDiagonal &&
        op.subtype == OpSubtype::TransverseX_offdiag && op.index >= 0 &&
        op.index < static_cast<int>(current_spins.size())) {
      current_spins[op.index] *= -1;
    }
  }

  for (size_t i = 0; i < current_spins.size(); ++i) {
    if (current_spins[i] != initial_spins[i]) {
      std::cerr << "Spin consistency check failed at site " << i << "\n";
      std::exit(EXIT_FAILURE);
    }
  }
}

void SSEUpdates::adjust_operator_string(SSEConfig &cfg) {
  // Sandvik scaling: new_L = n + n/3
  int new_L = cfg.n_ops + cfg.n_ops / 3;
  if (new_L < 20)
    new_L = 20;

  if (new_L > cfg.max_ops) {
    cfg.op_string.resize(new_L);
    for (int p = cfg.max_ops; p < new_L; ++p) {
      cfg.op_string[p].type = OpType::Identity;
      cfg.op_string[p].subtype = OpSubtype::None;
      cfg.op_string[p].index = -1;
    }
    cfg.max_ops = new_L;
  }
}

// ======================================================
// 9. SLMC global proposal (skeleton)
// ======================================================
void SSEUpdates::slmc_update(SSEConfig &cfg, const Parameters & /*prm*/,
                             const Geometry & /*geom*/, RNG &rng) {
  (void)cfg;
  (void)rng;
  throw std::logic_error("SLMC is not implemented; use the validated SSE baseline.");
}

// ======================================================
// 10. Gauge-sector update (row / column flips)
// ======================================================
// Flipping every spin in row r touches exactly the 2 sites of each black
// plaquette that straddles row r (its bottom-row sites and its top-row sites
// each come in pairs within that plaquette). Flipping both → prod unchanged.
// Same argument holds for column flips.  Both moves always satisfy detailed
// balance (w_new = w_old) so they are accepted unconditionally.
// These are supplementary pure-J0 sector moves. Row/column flips alone do
// not connect the complete ground-state manifold.
void SSEUpdates::gauge_update(SSEConfig &cfg, const Parameters &prm,
                              const Geometry &/*geom*/, RNG &rng) {
  // Only valid when J0 is the only interaction. J1 plaquette weights use
  // diagonal spin products s[0]*s[2] + s[1]*s[3]: flipping a full row negates
  // both products (sign flip), so row flips invalidate J1 weights. J2/J3/hx
  // add bonds that cross rows/columns and are similarly broken by row/col flips.
  if (prm.J1 != 0.0 || prm.J2 != 0.0 || prm.J3 != 0.0 || prm.hx != 0.0)
    return;

  const int Lx = prm.Lx;
  const int Ly = prm.Ly;

  // Attempt an independent row flip for each row
  for (int r = 0; r < Ly; ++r) {
    if (rng.uniform() < 0.5) {
      for (int x = 0; x < Lx; ++x)
        cfg.spins.spin[x + Lx * r] *= -1;
    }
  }

  // Attempt an independent column flip for each column
  for (int c = 0; c < Lx; ++c) {
    if (rng.uniform() < 0.5) {
      for (int y = 0; y < Ly; ++y)
        cfg.spins.spin[c + Lx * y] *= -1;
    }
  }
}


// Symmetric whole-worldline proposal. Transverse-X events commute with the
// boundary spin inversion; only diagonal matrix elements change. Zero weights
// reject the move. The ratio uses every affected propagated vertex.
void SSEUpdates::canonical_update(SSEConfig& cfg, const Parameters& prm, const Geometry& geom, RNG& rng) {
  if(local_skip_probability>0 && rng.uniform()<local_skip_probability){++canonical_skips;return;}
  ++canonical_attempts;
  int site = int(rng.uniform() * prm.Ns);
  SSEConfig old_state, new_state;
  old_state.spins = cfg.spins;
  new_state.spins = cfg.spins;
  new_state.spins.spin[site] *= -1;
  double log_ratio = 0.0;
  for (const auto& op : cfg.op_string) {
    if (op.type == OpType::OffDiagonal) {
      old_state.spins.spin[op.index] *= -1;
      new_state.spins.spin[op.index] *= -1;
    } else if (op.type == OpType::Diagonal) {
      double before = diagonal_weight(int(op.subtype),op.index,old_state,prm,geom);
      double after = diagonal_weight(int(op.subtype),op.index,new_state,prm,geom);
      if (after <= 0 || before <= 0) return;
      log_ratio += std::log(after/before);
    }
  }
  if (std::log(std::max(rng.uniform(),1e-300)) < std::min(0.0,log_ratio)) {
    cfg.spins.spin[site] *= -1;
    ++canonical_accepted;
  }
}

#pragma once
// ══════════════════════════════════════════════════════════════════════════════
// sse_debug.hpp — SSE diagnostic / debug output
//
// Enable: compile with -DSSE_DEBUG  (make debug  from backend/build/)
//
// Every output line carries a [DBG-XXX] prefix so the web-UI terminal can
// colour-code it by tag.  ANSI escape codes are ALSO embedded so the output
// looks rich in a native terminal (VS Code, Windows Terminal, etc.).
// The web UI strips ANSI via:  str.replace(/\x1b\[[0-9;]*m/g, '')
// ══════════════════════════════════════════════════════════════════════════════
#ifdef SSE_DEBUG

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "variables.hpp"

// ── ANSI colours ─────────────────────────────────────────────────────────────
static const char* const A0    = "\033[0m";       // reset
static const char* const ABLD  = "\033[1m";       // bold
static const char* const ADIM  = "\033[2m";       // dim

static const char* const AW    = "\033[97m";      // bright white
static const char* const AGY   = "\033[90m";      // gray
static const char* const ACY   = "\033[96m";      // cyan
static const char* const AYL   = "\033[93m";      // yellow
static const char* const AGN   = "\033[92m";      // green
static const char* const ARD   = "\033[91m";      // red
static const char* const AMG   = "\033[95m";      // magenta
static const char* const ABL   = "\033[94m";      // blue
static const char* const AOR   = "\033[38;5;208m";// orange

// Dark background panels (256-colour; visible in native terminal only)
static const char* const BGI   = "\033[48;5;17m"; // navy   → init
static const char* const BGO   = "\033[48;5;22m"; // green  → op-string
static const char* const BGL   = "\033[48;5;53m"; // violet → linked list
static const char* const BGB   = "\033[48;5;88m"; // red    → loop

// ── UI tag prefixes  (web terminal colours on these) ─────────────────────────
static const char* const PI    = "[DBG-INIT] ";
static const char* const PO    = "[DBG-OPS]  ";
static const char* const PL    = "[DBG-LIST] ";
static const char* const PB    = "[DBG-LOOP] ";
static const char* const PW    = "[DBG-WARN] ";

// ── Tiny helpers ──────────────────────────────────────────────────────────────
static inline std::string spin_icon(int s) {
    if (s > 0) return std::string(AGN) + "UP" + A0;
    if (s < 0) return std::string(ARD) + "DN" + A0;
    return std::string(AYL) + "?" + A0;   // 0 = not set
}

static inline std::string ok_fail(bool b) {
    return b ? (std::string(AGN) + ABLD + "OK" + A0)
             : (std::string(ARD) + ABLD + "NO" + A0);
}

static inline void section_hdr(const char* pfx, const char* bg, const char* title) {
    std::string t = title;
    int pad = std::max(0, 62 - (int)t.size());
    std::cout << pfx << ABLD << bg << AW << "  -- " << t << " ";
    for (int k = 0; k < pad; ++k) std::cout << "-";
    std::cout << A0 << "\n";
}
static inline void section_ftr(const char* pfx, const char* bg) {
    std::cout << pfx << ABLD << bg << AW << "  ";
    for (int k = 0; k < 68; ++k) std::cout << "-";
    std::cout << A0 << "\n" << std::flush;
}

static inline const char* sub_clr(OpSubtype s) {
    switch (s) {
        case OpSubtype::J0_Plaquette:        return AGN;
        case OpSubtype::J1_Plaquette:        return AYL;
        case OpSubtype::J2_Dipole:           return ACY;
        case OpSubtype::J3_Inter:            return ABL;
        case OpSubtype::TransverseX_diag:    return AOR;
        case OpSubtype::TransverseX_offdiag: return AMG;
        default:                             return AGY;
    }
}
static inline const char* sub_nm(OpSubtype s) {
    switch (s) {
        case OpSubtype::J0_Plaquette:        return "J0_plaq ";
        case OpSubtype::J1_Plaquette:        return "J1_bond ";
        case OpSubtype::J2_Dipole:           return "J2_dip  ";
        case OpSubtype::J3_Inter:            return "J3_inter";
        case OpSubtype::TransverseX_diag:    return "hx_diag ";
        case OpSubtype::TransverseX_offdiag: return "hx_offD ";
        default:                             return "None    ";
    }
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// 0. INIT CHECK
//    Call once after initializing cfg / prm / geom and before the sweep loop.
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
inline void dbg_check_init(const SSEConfig&   cfg,
                           const Parameters&  prm,
                           const Geometry&    geom)
{
    section_hdr(PI, BGI, "0. INIT CHECK");

    // ── spin vector ──────────────────────────────────────────────────────────
    int n_up = 0, n_dn = 0, n_bad = 0;
    for (int s : cfg.spins.spin) {
        if      (s ==  1) ++n_up;
        else if (s == -1) ++n_dn;
        else               ++n_bad;
    }
    bool spins_ok = (n_bad == 0 && (int)cfg.spins.spin.size() == prm.Ns);
    std::cout << PI << "  Spins     "
              << AW << cfg.spins.spin.size() << A0 << " total  "
              << AGN << "(+)" << A0 << n_up << "  "
              << ARD << "(-)" << A0 << n_dn;
    if (n_bad) std::cout << "  " << ARD << ABLD << "BAD=" << n_bad << A0;
    std::cout << "  (expected Ns=" << prm.Ns << ")  " << ok_fail(spins_ok) << "\n";
    std::cout << PI << "  Spin vals: ";
    for (int i = 0; i < (int)cfg.spins.spin.size(); ++i)
        std::cout << "s" << i << "=" << spin_icon(cfg.spins.spin[i]) << " ";
    std::cout << "\n";

    // ── op_string ─────────────────────────────────────────────────────────────
    bool all_id = true;
    for (int p = 0; p < cfg.max_ops; ++p)
        if (cfg.op_string[p].type != OpType::Identity) { all_id = false; break; }
    // Note: default constructor sets subtype=J0_Plaquette, not None
    bool subtype_warn = !cfg.op_string.empty() && (cfg.op_string[0].subtype != OpSubtype::None);
    std::cout << PI << "  op_string "
              << AW << cfg.op_string.size() << A0 << " slots  "
              << "n_ops=" << cfg.n_ops << "  max_ops=" << cfg.max_ops
              << "  all_ID=" << ok_fail(all_id);
    if (subtype_warn)
        std::cout << "  " << AYL << "(!) default subtype=J0_Plaquette (should be None)" << A0;
    std::cout << "\n";

    // ── parameters ────────────────────────────────────────────────────────────
    bool prm_ok = true;
    auto ck = [&](double v, const char* nm) {
        if (std::isnan(v) || std::isinf(v)) {
            std::cout << PW << "  " << ARD << nm << " is NaN/Inf!" << A0 << "\n";
            prm_ok = false;
        }
    };
    ck(prm.J0,"J0"); ck(prm.J1,"J1"); ck(prm.J2,"J2"); ck(prm.J3,"J3");
    ck(prm.hx,"hx"); ck(prm.beta,"beta");
    bool ns_ok = (prm.Ns == prm.Lx * prm.Ly * prm.Lz);
    if (!ns_ok)
        std::cout << PW << "  " << ARD << "Ns=" << prm.Ns
                  << " != Lx*Ly*Lz=" << (prm.Lx*prm.Ly*prm.Lz) << A0 << "\n";
    std::cout << PI << "  Params    "
              << "J0=" << AYL << prm.J0 << A0
              << " J1=" << AYL << prm.J1 << A0
              << " J2=" << AYL << prm.J2 << A0
              << " J3=" << AYL << prm.J3 << A0
              << " hx=" << AYL << prm.hx << A0
              << " beta="  << ACY << prm.beta << A0
              << " Ns=" << prm.Ns << " Np=" << prm.Np
              << "  " << ok_fail(prm_ok && ns_ok) << "\n";

    // ── geometry ──────────────────────────────────────────────────────────────
    bool geom_ok = !geom.blackPlaquettes.empty();
    std::cout << PI << "  Geometry  "
              << "blackPlaqs=" << AW << geom.blackPlaquettes.size() << A0
              << "  j2bonds="  << AW << geom.j2bonds.size()         << A0
              << "  j3bonds="  << AW << geom.j3bonds.size()         << A0
              << "  " << ok_fail(geom_ok) << "\n";

    // ── known issues ──────────────────────────────────────────────────────────
    std::cout << PW << "  " << AYL
              << "vertex_spin was NOT set for J0/J1/J2/J3 legs in original build_vertex_list"
              << A0 << " - " << AGN << "FIXED in debug build" << A0 << "\n";
    std::cout << PW << "  " << AGY
              << "choose_exit_leg() is never called from directed_loop_update() "
              << "— current update is Wolff BFS cluster (not true directed loop)"
              << A0 << "\n";

    section_ftr(PI, BGI);
    std::cout << "\n";
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// 1. OPERATOR STRING  (call after diagonal_update)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
inline void dbg_print_op_string(const SSEConfig&  cfg,
                                const Parameters& prm,
                                const Geometry&   geom)
{
    section_hdr(PO, BGO, "1. OPERATOR STRING  (after diagonal_update)");
    std::cout << PO << "  n_ops=" << AW << cfg.n_ops << A0
              << "  max_ops=" << cfg.max_ops << "\n";
    std::cout << PO << "  t=0 spins: ";
    for (int i = 0; i < (int)cfg.spins.spin.size(); ++i)
        std::cout << "s" << i << "=" << spin_icon(cfg.spins.spin[i]) << " ";
    std::cout << "\n\n";

    int id_run = 0;
    for (int p = 0; p < cfg.max_ops; ++p) {
        const auto& op = cfg.op_string[p];

        if (op.type == OpType::Identity) { ++id_run; continue; }

        if (id_run > 0) {
            std::cout << PO << "  " << AGY << ADIM
                      << "  ... " << id_run << " Identity slot"
                      << (id_run > 1 ? "s" : "") << " ···" << A0 << "\n";
            id_run = 0;
        }

        const char* sc = sub_clr(op.subtype);

        if (op.type == OpType::Diagonal) {
            std::cout << PO << "  p=" << std::setw(3) << p
                      << "  " << sc << ABLD
                      << "[DIA " << sub_nm(op.subtype)
                      << " idx=" << std::setw(3) << op.index << "]" << A0;

            if (op.subtype == OpSubtype::J0_Plaquette
                && op.index >= 0 && op.index < (int)geom.blackPlaquettes.size()) {
                const auto& pl = geom.blackPlaquettes[op.index];
                int prod = 1;
                for (int s : pl.sites) prod *= cfg.spins.spin[s];
                double w = std::abs(prod + 1) * std::abs(prm.J0);
                std::cout << "  sites=[" << pl.sites[0] << "," << pl.sites[1]
                          << "," << pl.sites[2] << "," << pl.sites[3] << "]"
                          << "  sprod=" << (prod > 0 ? std::string(AGN)+"+1"+A0 : std::string(ARD)+"-1"+A0)
                          << "  " << ACY << "w=" << w << A0;

            } else if (op.subtype == OpSubtype::J1_Plaquette && op.index >= 0 && op.index/2 < (int)geom.blackPlaquettes.size()) {
                const auto& pl = geom.blackPlaquettes[op.index/2];
                int a=op.index%2, b=a+2;
                std::cout << " bond=" << pl.sites[a] << "-" << pl.sites[b];
            } else if (op.subtype == OpSubtype::J2_Dipole
                       && op.index >= 0 && op.index < (int)geom.j2bonds.size()) {
                const auto& b = geom.j2bonds[op.index];
                double Jl = b.J_val * prm.J2;
                double w  = std::abs(Jl) - Jl * cfg.spins.spin[b.s1] * cfg.spins.spin[b.s2];
                std::cout << "  sites=[" << b.s1 << "," << b.s2 << "]"
                          << "  J_val=" << b.J_val
                          << "  " << spin_icon(cfg.spins.spin[b.s1])
                          << spin_icon(cfg.spins.spin[b.s2])
                          << "  " << ACY << "w=" << w << A0;

            } else if (op.subtype == OpSubtype::J3_Inter
                       && op.index >= 0 && op.index < (int)geom.j3bonds.size()) {
                const auto& b = geom.j3bonds[op.index];
                int s0 = cfg.spins.spin[b[0]], s1 = cfg.spins.spin[b[1]];
                double w = (s0 * s1 == -1) ? 2.0 * std::abs(prm.J3) : 0.0;
                std::cout << "  sites=[" << b[0] << "," << b[1] << "]"
                          << "  " << spin_icon(s0) << spin_icon(s1)
                          << "  " << ACY << "w=" << w << A0;

            } else if (op.subtype == OpSubtype::TransverseX_diag && op.index >= 0) {
                std::cout << "  site=" << op.index
                          << "  s=" << spin_icon(cfg.spins.spin[op.index])
                          << "  " << ACY << "w=" << std::abs(prm.hx) << A0;
            }
        } else { // OffDiagonal
            std::cout << PO << "  p=" << std::setw(3) << p
                      << "  " << sc << ABLD
                      << "[OFF " << sub_nm(op.subtype)
                      << " idx=" << std::setw(3) << op.index << "]" << A0;
            if (op.index >= 0 && op.index < (int)cfg.spins.spin.size())
                std::cout << "  site=" << op.index
                          << "  s_in=" << spin_icon(cfg.spins.spin[op.index])
                          << "  " << AMG << "(spin flips here)" << A0;
        }
        std::cout << "\n";
    }
    if (id_run > 0) {
        std::cout << PO << "  " << AGY << ADIM
                  << "  ··· " << id_run << " Identity slot"
                  << (id_run > 1 ? "s" : "") << " ···" << A0 << "\n";
    }
    section_ftr(PO, BGO);
    std::cout << "\n";
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// 3. VERTEX LINKED LIST  (call after build_vertex_list)
//    All vectors are passed by const-ref so we don't need friend access.
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
inline void dbg_print_vertex_list(
    const std::vector<int>&         vertex_op_index,
    const std::vector<int>&         linked_leg,
    const std::vector<int>&         vertex_partner,
    const std::vector<int>&         vertex_spin,
    const std::vector<int>&         leg_site,
    const std::vector<int>&         vertex_base_vec,  // per-vertex (size=n_vtx)
    const std::vector<int>&         vertex_decomp,    // per op-slot  (size=max_ops)
    const std::vector<SSEOperator>& op_string,
    const std::vector<int>&         first_leg,
    const std::vector<int>&         last_leg,
    int                             n_legs,
    int                             n_sites,
    int                             pass_idx)
{
    std::string t = "3. VERTEX LINKED LIST  (directed_loop pass "
                  + std::to_string(pass_idx) + ")";
    section_hdr(PL, BGL, t.c_str());

    int n_vtx = (int)vertex_op_index.size();
    std::cout << PL << "  n_vertices=" << AW << n_vtx << A0
              << "  n_legs=" << AW << n_legs << A0 << "\n";

    // ── imaginary-time boundaries ─────────────────────────────────────────────
    std::cout << PL << "\n"
              << PL << "  " << ABLD << "Periodic boundaries (first/last per site):" << A0 << "\n";
    for (int s = 0; s < n_sites; ++s) {
        if (first_leg[s] == -1) {
            std::cout << PL << "    site " << std::setw(3) << s
                      << ": " << AGY << "no operators on this site" << A0 << "\n";
        } else {
            int fl = first_leg[s], ll = last_leg[s];
            bool pbc_ok = (linked_leg[ll] == fl) && (linked_leg[fl] == ll);
            std::cout << PL << "    site " << std::setw(3) << s
                      << ": first=leg" << std::setw(4) << fl
                      << "  last=leg"  << std::setw(4) << ll
                      << "  " << ABL << "linked[last]=" << linked_leg[ll]
                      << " linked[first]=" << linked_leg[fl] << A0
                      << "  PBC:" << ok_fail(pbc_ok) << "\n";
        }
    }
    std::cout << PL << "\n";

    // ── per-vertex detail ─────────────────────────────────────────────────────
    for (int vi = 0; vi < n_vtx; ++vi) {
        int p        = vertex_op_index[vi];
        const auto& op = op_string[p];
        int base     = vertex_base_vec[vi];

        int n_vlegs = 2;
        if (op.subtype == OpSubtype::J0_Plaquette)
            n_vlegs = 8;
        else if (op.subtype == OpSubtype::J1_Plaquette || op.subtype == OpSubtype::J2_Dipole || op.subtype == OpSubtype::J3_Inter)
            n_vlegs = 4;

        const char* sc = sub_clr(op.subtype);
        std::cout << PL << "  " << ABLD << sc
                  << "Vtx " << std::setw(2) << vi
                  << "  op[" << std::setw(3) << p << "]"
                  << "  " << (op.type == OpType::Diagonal ? "DIA" : "OFF")
                  << " " << sub_nm(op.subtype)
                  << " idx=" << op.index;
        if ((op.subtype == OpSubtype::J0_Plaquette || op.subtype == OpSubtype::J1_Plaquette)
            && p < (int)vertex_decomp.size())
            std::cout << "  decomp=" << ACY << vertex_decomp[p] << A0;
        std::cout << A0 << "\n";

        for (int l = 0; l < n_vlegs; ++l) {
            int leg     = base + l;
            int site    = (leg < (int)leg_site.size())      ? leg_site[leg]      : -1;
            int linked  = (leg < (int)linked_leg.size())    ? linked_leg[leg]    : -9;
            int partner = (leg < (int)vertex_partner.size())? vertex_partner[leg] : -9;
            int sv      = (leg < (int)vertex_spin.size())   ? vertex_spin[leg]   :  0;
            bool is_in  = (l % 2 == 0);

            std::cout << PL << "    leg" << std::setw(4) << leg
                      << (is_in ? " [in ] " : " [out] ")
                      << " site=" << std::setw(3) << site
                      << "  s=" << spin_icon(sv)
                      << "  <t> " << ABL << "leg" << std::setw(4) << linked << A0;
            if (partner >= 0)
                std::cout << "  partner→ " << AOR << "leg" << std::setw(4) << partner << A0;
            std::cout << "\n";
        }
    }

    section_ftr(PL, BGL);
    std::cout << "\n";
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// 4. DIRECTED LOOP trace helpers
//    Call these at the indicated points inside directed_loop_update().
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

static inline void dbg_loop_pass_hdr(int pass, int n_legs) {
    std::string t = "4. DIRECTED LOOP BFS  pass " + std::to_string(pass)
                  + "  n_legs=" + std::to_string(n_legs);
    section_hdr(PB, BGB, t.c_str());
    std::cout << PB << "  " << AYL
              << "(BFS Wolff-cluster — every cluster flipped independently with p=0.5)"
              << A0 << "\n\n";
}

static inline void dbg_loop_cluster_start(int cid, int start_leg, double r, bool flip,
                                          const std::vector<int>& leg_site) {
    int site = (start_leg < (int)leg_site.size()) ? leg_site[start_leg] : -1;
    std::cout << PB << "  " << ABLD
              << (flip ? AGN : AGY) << "Cluster " << cid << A0
              << "  start=leg" << std::setw(4) << start_leg
              << " (site " << site << ")"
              << "  r=" << ACY << std::fixed << std::setprecision(4) << r << A0
              << "  --> " << (flip
                  ? (std::string(ABLD) + AGN + "FLIP" + A0)
                  : (std::string(AGY)  + "KEEP" + A0))
              << "\n";
}

static inline void dbg_loop_bfs_expand(int from_leg, int to_leg,
                                       const char* via,
                                       const std::vector<int>& leg_site) {
    int fs = (from_leg < (int)leg_site.size()) ? leg_site[from_leg] : -1;
    int ts = (to_leg   < (int)leg_site.size()) ? leg_site[to_leg]   : -1;
    std::cout << PB << "    " << AGY
              << "leg" << std::setw(4) << from_leg << "(s" << fs << ")"
              << " --" << via << "> "
              << "leg" << std::setw(4) << to_leg   << "(s" << ts << ")"
              << A0 << "\n";
}

static inline void dbg_loop_cluster_end(int cid, int n_members) {
    std::cout << PB << "    " << AGY << ADIM
              << "[cluster " << cid << ": " << n_members << " leg"
              << (n_members != 1 ? "s" : "") << " total]"
              << A0 << "\n";
}

static inline void dbg_final_spins(const std::vector<int>& before,
                                   const std::vector<int>& after) {
    std::cout << PB << "\n"
              << PB << "  " << ABLD << AW
              << "-- Final spins after cluster flips --" << A0 << "\n";
    int n_flip = 0;
    for (int s = 0; s < (int)after.size(); ++s) {
        bool changed = (before[s] != after[s]);
        if (changed) ++n_flip;
        std::cout << PB << "    site " << std::setw(3) << s << ":  "
                  << spin_icon(before[s]) << " --> " << spin_icon(after[s])
                  << "  " << (changed
                      ? (std::string(ABLD) + AGN + "FLIPPED" + A0)
                      : (std::string(AGY)  + "kept   " + A0))
                  << "\n";
    }
    std::cout << PB << "  total: "
              << AGN << n_flip << " site" << (n_flip != 1 ? "s" : "") << " flipped" << A0
              << "  " << AGY << (int)after.size() - n_flip << " kept" << A0 << "\n";
}

static inline void dbg_loop_toggle_section() {
    std::cout << PB << "\n"
              << PB << "  " << ABLD << AW
              << "-- TransverseX operator toggles --" << A0 << "\n";
}

static inline void dbg_loop_op_toggle(int p, int site,
                                      const char* from_s, const char* to_s) {
    std::cout << PB << "    op[" << std::setw(3) << p << "]"
              << " site=" << site
              << "  " << AOR << from_s << A0 << " -> " << AMG << to_s << A0
              << "  " << AMG << "(cluster-boundary crossing)" << A0 << "\n";
}

static inline void dbg_loop_footer(int n_clusters, int n_flipped) {
    std::cout << PB << "\n"
              << PB << "  Summary: clusters=" << AW << n_clusters << A0
              << "  flipped=" << AGN << n_flipped << A0
              << "  kept=" << AGY << (n_clusters - n_flipped) << A0 << "\n";
    section_ftr(PB, BGB);
    std::cout << "\n";
}

#endif // SSE_DEBUG

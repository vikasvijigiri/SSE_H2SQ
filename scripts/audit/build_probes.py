"""Reproduce historical defect-isolation experiments in isolated executables."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import uuid

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / ".cache/physics-audit"
BASELINE_REVISION = "f313315ad1b3bd0c1d7f632ba078cbae5061a8b7"


def baseline_source(path):
    return subprocess.check_output(["git", "show", f"{BASELINE_REVISION}:{path}"], cwd=ROOT).decode("utf-8")


def prepare_baseline():
    paths = subprocess.check_output(["git", "ls-tree", "-r", "--name-only", BASELINE_REVISION, "backend/include"], cwd=ROOT, text=True).splitlines()
    paths += ["backend/src/variables.cpp", "backend/src/observables.cpp", "backend/Analysis/ED_analysis/exact_diagonalzation_H2SQ.py"]
    for path in paths:
        destination = OUT / "historical-baseline" / Path(path).relative_to("backend")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(baseline_source(path), encoding="utf-8")

def replace_once(source, old, new):
    if source.count(old) != 1:
        raise RuntimeError("Production source changed; review audit transformation: " + old[:80])
    return source.replace(old, new, 1)

def transformed(variant):
    source = baseline_source("backend/src/updates.cpp")
    source = source.replace('#include "../include/', '#include "historical-baseline/include/')
    if variant != "baseline":
        source = replace_once(source, '  vertex_base_vec.reserve(vertex_op_index.size());',
            '  vertex_base_vec.reserve(vertex_op_index.size());\n  std::vector<int> audit_motifs(geom.blackPlaquettes.size());\n  for (auto& motif : audit_motifs) motif = int(rng.uniform() * 3.0);')
        source = replace_once(source, 'decomp_type = static_cast<int>(rng.uniform() * 3.0);', 'decomp_type = audit_motifs[op.index];')
    if variant == "intended_combined":
        source = replace_once(source, 'return (prod_sum == 2 ? 4.0 * std::abs(prm.J1) : 0.0);', 'return 2.0 * std::abs(prm.J1) - prm.J1 * prod_sum;')
    if variant == "intended_bonds":
        source = replace_once(source, 'add_type(OpSubtype::J1_Plaquette,\n             static_cast<int>(geom.blackPlaquettes.size()));', 'add_type(OpSubtype::J1_Plaquette,\n             2 * static_cast<int>(geom.blackPlaquettes.size()));')
        start = source.index('    const auto &p = geom.blackPlaquettes[index];',source.index('case OpSubtype::J1_Plaquette:'))
        end = source.index('\n  }',start)
        source = source[:start] + '''    const auto &p = geom.blackPlaquettes[index / 2];
    int a = (index % 2 == 0) ? 0 : 1;
    int b = a + 2;
    return std::abs(prm.J1) - prm.J1 * cfg.spins.spin[p.sites[a]] * cfg.spins.spin[p.sites[b]];''' + source[end:]
        start = source.index('      const auto &plaq =',source.index('} else if (op.subtype == OpSubtype::J1_Plaquette)'))
        end = source.index('\n    } else if',start)
        source = source[:start] + '''      const auto &plaq = geom.blackPlaquettes[op.index / 2];
      int a = (op.index % 2 == 0) ? 0 : 1;
      sites_buf[0] = plaq.sites[a];
      sites_buf[1] = plaq.sites[a + 2];
      n_sites = 2;
      n_op_legs = 4;''' + source[end:]
        source = replace_once(source, '          set_pair_cycle(base + 0, base + 4); // 2*0, 2*2\n          set_pair_cycle(base + 2, base + 6); // 2*1, 2*3', '          set_pair_cycle(base, base + 2);')
    return source

def build(variant):
    OUT.mkdir(parents=True, exist_ok=True)
    prepare_baseline()
    update = OUT / (variant + "_updates.cpp")
    update.write_text(transformed(variant), encoding="utf-8")
    probe = (ROOT / "scripts/audit/sse_probe.cpp").read_text(encoding="utf-8").replace('../../backend/include/', 'historical-baseline/include/')
    if variant.startswith("intended"):
        probe = replace_once(probe, 'e+=std::abs(prm.J1)*(1-a-b-a*b);', 'e+=prm.J1*(a+b);')
    probe_path = OUT / (variant + "_probe.cpp")
    probe_path.write_text(probe, encoding="utf-8")
    names = {"baseline": "sse_probe", "premarked": "sse_probe_premarked", "intended_combined": "sse_probe_intended", "intended_bonds": "sse_probe_pair_bonds"}
    executable = OUT / (names[variant] + (".exe" if os.name == "nt" else ""))
    command = ["g++", "-O2", "-std=c++17", str(probe_path), str(update), str(OUT / "historical-baseline/src/variables.cpp"), str(OUT / "historical-baseline/src/observables.cpp")]
    if os.name == "nt":
        command += ["-static-libgcc", "-static-libstdc++", "-Wl,-Bstatic,--whole-archive", "-lwinpthread", "-Wl,--no-whole-archive,-Bdynamic"]
    staged_executable = executable.with_name(executable.stem + "-" + uuid.uuid4().hex + executable.suffix)
    command += ["-o", str(staged_executable)]
    subprocess.run(command, check=True)
    staged_executable.replace(executable)
    provenance = {"variant": variant, "command": command, "historical_baseline_revision": BASELINE_REVISION, "git_revision": subprocess.check_output(["git", "rev-parse", "HEAD"],cwd=ROOT,text=True).strip(), "updates_sha256": hashlib.sha256(update.read_bytes()).hexdigest(), "probe_sha256": hashlib.sha256(probe_path.read_bytes()).hexdigest()}
    (OUT / (variant + "_build.json")).write_text(json.dumps(provenance, indent=2))
    print(executable)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant", choices=["baseline", "premarked", "intended_combined", "intended_bonds"])
    build(parser.parse_args().variant)

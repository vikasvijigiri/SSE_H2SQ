Lx, Ly = 4, 4
cell = 80
margin = 55
r = 9
W = Lx * cell + 2 * margin
H = Ly * cell + 2 * margin

def sx(i): return margin + i * cell
def sy(j): return margin + j * cell

# Build plaquette map: SVG cell (ci,cj) -> blackPlaquette index
plaq = {}
pid = 0
for y in range(Ly):
    for x in range(Lx):
        if (x + y) % 2 == 0:
            plaq[(x, (Ly - 2 - y) % Ly)] = pid
            pid += 1

BG          = "#ffffff"
CELL_FILL   = "#1a1a2e"   # dark navy — strong contrast against white bg
HATCH_COL   = "#ffffff"   # pure white — maximum contrast against dark cell fill
BORDER_COL  = "#3a3a90"   # blue border
LABEL_COL   = "#ffffff"   # white labels on dark cell
BOND_COL    = "#222222"   # near-black bonds, thick, clearly visible on white bg
NODE_FILL   = "#1a50a0"   # blue nodes
NODE_GHOST  = "#aaaaaa"   # gray ghost nodes

def crosshatch(out, x0, y0, sz, spacing, color, sw):
    """Draw explicit NW-SE and SW-NE diagonal lines clipped to square cell."""
    for d in range(-(sz - spacing), sz, spacing):
        # NW-SE (slope +1)
        if d >= 0:
            ax, ay, bx, by = x0+d, y0, x0+sz, y0+sz-d
        else:
            ax, ay, bx, by = x0, y0-d, x0+sz+d, y0+sz
        out.append(f'  <line x1="{ax}" y1="{ay}" x2="{bx}" y2="{by}" stroke="{color}" stroke-width="{sw}" stroke-linecap="square"/>')
        # SW-NE (slope -1)
        if d >= 0:
            ax, ay, bx, by = x0+d, y0+sz, x0+sz, y0+d
        else:
            ax, ay, bx, by = x0, y0+sz+d, x0+sz+d, y0
        out.append(f'  <line x1="{ax}" y1="{ay}" x2="{bx}" y2="{by}" stroke="{color}" stroke-width="{sw}" stroke-linecap="square"/>')

out = []
out.append(f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg">')
out.append(f'  <rect width="{W}" height="{H}" fill="{BG}"/>')

# Black plaquettes: solid fill + crosshatch + border + label
for (ci, cj), p in plaq.items():
    x0, y0 = sx(ci), sy(cj)
    cx, cy = x0 + cell // 2, y0 + cell // 2
    out.append(f'  <rect x="{x0}" y="{y0}" width="{cell}" height="{cell}" fill="{CELL_FILL}"/>')
    crosshatch(out, x0, y0, cell, 12, HATCH_COL, 2.5)
    out.append(f'  <rect x="{x0}" y="{y0}" width="{cell}" height="{cell}" fill="none" stroke="{BORDER_COL}" stroke-width="2"/>')
    out.append(f'  <text x="{cx}" y="{cy+1}" text-anchor="middle" dominant-baseline="middle" '
               f'fill="{LABEL_COL}" font-size="13" font-family="monospace" font-weight="bold">P{p}</text>')

# Bonds (drawn after hatch so they're on top)
for nj in range(Ly + 1):
    for ni in range(Lx):
        out.append(f'  <line x1="{sx(ni)}" y1="{sy(nj)}" x2="{sx(ni+1)}" y2="{sy(nj)}" stroke="{BOND_COL}" stroke-width="2.5"/>')
for ni in range(Lx + 1):
    for nj in range(Ly):
        out.append(f'  <line x1="{sx(ni)}" y1="{sy(nj)}" x2="{sx(ni)}" y2="{sy(nj+1)}" stroke="{BOND_COL}" stroke-width="2.5"/>')

# Site nodes
for ni in range(Lx + 1):
    for nj in range(Ly + 1):
        px, py = sx(ni), sy(nj)
        lat_y = (Ly - 1 - nj) % Ly
        site  = (ni % Lx) + Lx * lat_y
        periodic = (ni == Lx or nj == Ly)
        fc = NODE_GHOST if periodic else NODE_FILL
        sc = "#888888" if periodic else "#0a2870"
        dash = 'stroke-dasharray="3,2"' if periodic else ''
        out.append(f'  <circle cx="{px}" cy="{py}" r="{r}" fill="{fc}" stroke="{sc}" stroke-width="1.5" {dash}/>')
        label_y = py - r - 4 if nj > 0 else py + r + 13
        out.append(f'  <text x="{px}" y="{label_y}" text-anchor="middle" fill="{fc}" font-size="9" font-family="monospace">{site}</text>')

# Title
out.append(f'  <text x="{W//2}" y="22" text-anchor="middle" fill="#222266" font-size="11" font-family="monospace">'
           f'H2SQ {Lx}x{Ly}  |  Black Plaquettes P0..P{pid-1}</text>')
out.append(f'  <text x="14" y="{H//2}" text-anchor="middle" fill="#444444" font-size="9" font-family="monospace" '
           f'transform="rotate(-90,14,{H//2})">y (lattice)</text>')
out.append(f'  <text x="{W//2}" y="{H-8}" text-anchor="middle" fill="#444444" font-size="9" font-family="monospace">x (lattice)</text>')
out.append('</svg>')

path = __import__('pathlib').Path(__file__).resolve().with_name('lattice_plaquettes.svg')
with open(path, 'w') as f:
    f.write('\n'.join(out))
print(f'Written to {path}')

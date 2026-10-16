import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import { Play, Square, Sparkles, TableProperties, Grid3x3, Terminal, Copy, Check, Bug } from 'lucide-react';
import { motion, AnimatePresence } from 'framer-motion';

const hostname = typeof window !== 'undefined' && window.location.hostname ? window.location.hostname : '127.0.0.1';
const actualHost = hostname === 'localhost' ? '127.0.0.1' : hostname;
const API_URL = `http://${actualHost}:8000/api`;
const WS_URL = `ws://${actualHost}:8000/ws`;

const InputField = ({ 
  label, 
  value, 
  onChange, 
  type = 'number', 
  step = '1' 
}: { 
  label: string; 
  value: any; 
  onChange: (val: any) => void; 
  type?: string; 
  step?: string; 
}) => (
  <div className="flex flex-col gap-1.5">
    <label className="text-[11px] font-semibold tracking-wide text-slate-400/80 uppercase">{label}</label>
    <input 
      type={type} 
      step={step}
      value={value}
      onChange={(e) => {
        const val = e.target.value;
        onChange(type === 'number' ? (val === '' ? '' : parseFloat(val)) : val);
      }}
      className="bg-slate-950/60 border border-slate-800/80 rounded-lg px-3 py-2 text-xs text-white focus:border-indigo-500/80 focus:ring-1 focus:ring-indigo-500/30 outline-none transition-all duration-200"
    />
  </div>
);

const SpinSVG = ({ spins, step, Lx, Ly, Lz = 1, cellSize, pad }: { spins: number[]; step: number; Lx: number; Ly: number; Lz?: number; cellSize: number; pad: number }) => {
  const plaquettes = [] as React.ReactElement[];
  const bonds      = [] as React.ReactElement[];
  const circles    = [] as React.ReactElement[];

  const actualLz = Lz > 1 ? Lz : Math.max(1, Math.floor(spins.length / (Lx * Ly)));

  const layerSpacing = Lx * cellSize + pad * 2;
  const svgW = actualLz * layerSpacing;
  const svgH = Ly * cellSize + 2 * pad;
  const r = cellSize * 0.38;

  // Plaquette index map: same enumeration as C++ blackPlaquettes
  // for y in 0..Ly-1, for x in 0..Lx-1, if (x+y)%2==0 → push
  const pidMap = new Map<string, number>();
  { let pid = 0;
    for (let y = 0; y < Ly; y++)
      for (let x = 0; x < Lx; x++)
        if ((x + y) % 2 === 0)
          pidMap.set(`${x},${y}`, pid++);
  }

  for (let z = 0; z < actualLz; z++) {
    const xOffset = z * layerSpacing;

    // Black plaquettes: fill + X mark + index label.
    // The SVG is sized so that plaquettes at x=Lx-1 and y=Ly-1 (periodic
    // boundary) also fit exactly at the edge — no clipping needed.
    for (let y = 0; y < Ly; y++) {
      for (let x = 0; x < Lx; x++) {
        if ((x + y) % 2 !== 0) continue;
        const pid = pidMap.get(`${x},${y}`) ?? -1;
        const px  = xOffset + pad + x * cellSize + cellSize / 2;
        const py  = pad + y * cellSize + cellSize / 2;
        const xcx = px + cellSize / 2;
        const ycy = py + cellSize / 2;
        const xp  = cellSize * 0.18;          // padding from corner for X arms

        plaquettes.push(
          <rect key={`plaq-${z}-${x}-${y}`}
            x={px} y={py} width={cellSize} height={cellSize}
            fill="rgba(80,80,200,0.20)" rx={2} />
        );
        {/* X mark — two diagonal lines, white semi-transparent */}
        plaquettes.push(
          <line key={`xa-${z}-${x}-${y}`}
            x1={px+xp} y1={py+xp} x2={px+cellSize-xp} y2={py+cellSize-xp}
            stroke="rgba(255,255,255,0.60)" strokeWidth={1.5} strokeLinecap="round" />
        );
        plaquettes.push(
          <line key={`xb-${z}-${x}-${y}`}
            x1={px+cellSize-xp} y1={py+xp} x2={px+xp} y2={py+cellSize-xp}
            stroke="rgba(255,255,255,0.60)" strokeWidth={1.5} strokeLinecap="round" />
        );
        {/* Plaquette index label */}
        {pid >= 0 && cellSize >= 18 && plaquettes.push(
          <text key={`plabel-${z}-${x}-${y}`}
            x={xcx} y={ycy}
            textAnchor="middle" dominantBaseline="middle"
            fontSize={Math.max(6, Math.floor(cellSize * 0.28))}
            fill="rgba(200,210,255,0.90)" fontFamily="monospace" fontWeight="bold"
            style={{ pointerEvents: 'none', userSelect: 'none' as const }}>
            P{pid}
          </text>
        )}
      }
    }

    // Bonds — white/light so they're visible on the dark background
    for (let y = 0; y < Ly; y++) {
      for (let x = 0; x < Lx; x++) {
        const cx = xOffset + pad + x * cellSize + cellSize / 2;
        const cy = pad + y * cellSize + cellSize / 2;
        if (x < Lx - 1) {
          const nx = xOffset + pad + (x + 1) * cellSize + cellSize / 2;
          bonds.push(<line key={`h-${z}-${x}-${y}`} x1={cx} y1={cy} x2={nx} y2={cy} stroke="#cbd5e1" strokeWidth={2} />);
        }
        if (y < Ly - 1) {
          const ny = pad + (y + 1) * cellSize + cellSize / 2;
          bonds.push(<line key={`v-${z}-${x}-${y}`} x1={cx} y1={cy} x2={cx} y2={ny} stroke="#cbd5e1" strokeWidth={2} />);
        }
      }
    }

    // Interlayer bonds
    if (z > 0) {
      for (let y = 0; y < Ly; y++) {
        for (let x = 0; x < Lx; x++) {
          const cx = xOffset + pad + x * cellSize + cellSize / 2;
          const cy = pad + y * cellSize + cellSize / 2;
          const prevCx = (z - 1) * layerSpacing + pad + x * cellSize + cellSize / 2;
          bonds.push(<line key={`il-${z}-${x}-${y}`} x1={prevCx} y1={cy} x2={cx} y2={cy} stroke="#94a3b8" strokeWidth={1.5} strokeDasharray="4,4" opacity={0.7} />);
        }
      }
    }

    // Spins
    for (let y = 0; y < Ly; y++) {
      for (let x = 0; x < Lx; x++) {
        const site = z * (Lx * Ly) + y * Lx + x;
        if (site >= spins.length) continue;
        const cx = xOffset + pad + x * cellSize + cellSize / 2;
        const cy = pad + y * cellSize + cellSize / 2;
        const isUp = spins[site] > 0;
        circles.push(
          <g key={`spin-${site}`}>
            <circle cx={cx} cy={cy} r={r + 1} fill={isUp ? '#1d4ed8' : '#b91c1c'} opacity={0.35} />
            <circle cx={cx} cy={cy} r={r} fill={isUp ? '#3b82f6' : '#ef4444'} stroke={isUp ? '#93c5fd' : '#fca5a5'} strokeWidth={1} />
            <text x={cx} y={cy + (cellSize > 18 ? 3.5 : 2.5)} textAnchor="middle" fontSize={cellSize > 18 ? 8 : 6} fontWeight="bold" fill="white" style={{ pointerEvents: 'none', userSelect: 'none' }}>
              {isUp ? '↑' : '↓'}
            </text>
          </g>
        );
      }
    }
  }

  return (
    <div className="flex flex-col items-center gap-1.5">
      <div className="text-[9px] font-mono font-bold text-slate-500 uppercase tracking-wider flex items-center justify-between w-full px-2">
        <span>Step {step.toLocaleString()}</span>
        {actualLz > 1 && <span className="text-slate-600">({actualLz} Layers)</span>}
      </div>
      <svg width={svgW} height={svgH} className="rounded-lg bg-slate-950/60 border border-slate-800/40">
        {plaquettes}
        {bonds}
        {circles}
      </svg>
    </div>
  );
};

const LatticeViewer = ({
  spinConfigs, Lx, Ly, Lz = 1
}: {
  spinConfigs: { step: number; spins: number[] }[];
  Lx: number;
  Ly: number;
  Lz?: number;
}) => {
  if (!spinConfigs || spinConfigs.length === 0) return null;

  // Fixed cell size so every site is readable, the SVG expands to fit.
  // cellSize target: ~28px for small lattices, scales down gracefully for larger ones.
  const maxDim = Math.max(Lx, Ly);
  const cellSize = Math.min(36, Math.max(12, Math.floor(360 / (maxDim * Math.max(1, Lz)))));
  const pad = cellSize * 0.5;
  const configs = spinConfigs.slice(0, 4);

  return (
    <div className="bg-[#0e0f1e]/60 border border-slate-800/60 rounded-2xl p-4 backdrop-blur-xl shadow-2xl relative overflow-x-auto custom-scrollbar">
      <div className="absolute top-0 left-0 right-0 h-[1px] bg-gradient-to-r from-transparent via-cyan-500/40 to-transparent" />
      <div className="flex items-center justify-between mb-4 flex-wrap gap-2">
        <h3 className="text-xs font-bold tracking-wide text-white uppercase flex items-center gap-1.5 min-w-max">
          <Grid3x3 className="w-3.5 h-3.5 text-cyan-400" />
          Spin Configuration Snapshots
          <span className="text-[9px] text-slate-500 font-normal normal-case tracking-normal ml-2">
            ({Lx}×{Ly}{Lz > 1 ? `×${Lz}` : ''} LATTICE · {spinConfigs.length} SNAPSHOTS DURING MEASUREMENT · BLUE=↑ RED=↓)
          </span>
        </h3>
      </div>
      
      <div className="flex items-center gap-4 w-max">
        {configs.map((config, idx) => (
          <SpinSVG 
            key={idx} 
            spins={config.spins} 
            step={config.step} 
            Lx={Lx} 
            Ly={Ly}
            Lz={Lz}
            cellSize={cellSize} 
            pad={pad} 
          />
        ))}
      </div>
    </div>
  );
};

const AlgorithmChecksPanel = ({ summary }: { summary: Record<string, any> }) => {
  const checkNames = Object.keys(summary).filter(k => !k.startsWith('_'));
  if (checkNames.length === 0) return null;
  const total = summary._total || 0;
  const failed = summary._failed || 0;
  const allPassed = summary._all_passed === true;

  return (
    <div className="bg-[#0e0f1e]/60 border border-slate-800/60 rounded-2xl p-4 backdrop-blur-xl shadow-2xl relative overflow-hidden">
      <div className="absolute top-0 left-0 right-0 h-[1px] bg-gradient-to-r from-transparent via-emerald-500/30 to-transparent" />
      <div className="flex items-center justify-between mb-3">
        <h3 className="text-xs font-bold tracking-wide text-white uppercase flex items-center gap-1.5">
          <span className="text-emerald-400 text-sm">⊕</span>
          Algorithm Checks
          <span className="text-[9px] font-normal text-slate-500 ml-1 lowercase">(n_ops · weights · τ-periodicity)</span>
        </h3>
        <span className={`text-[9px] font-bold px-2 py-0.5 rounded-full border ${
          allPassed
            ? 'text-emerald-400 border-emerald-500/30 bg-emerald-950/40'
            : 'text-red-400 border-red-500/30 bg-red-950/40'
        }`}>
          {total - failed}/{total} PASS
        </span>
      </div>
      <div className="flex items-center gap-2 flex-wrap">
        {checkNames.map(name => {
          const passed = summary[name] === 'PASS';
          return (
            <div key={name} className={`flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg border text-[10px] font-bold ${
              passed
                ? 'text-emerald-400 border-emerald-500/20 bg-emerald-950/30'
                : 'text-red-400 border-red-500/20 bg-red-950/30'
            }`}>
              <span>{passed ? '✓' : '✗'}</span>
              <span className="font-mono">{name}</span>
            </div>
          );
        })}
        {summary._energy_inst !== undefined && (
          <div className="flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg border text-[10px] border-slate-700/50 bg-slate-900/30 text-slate-400">
            <span className="text-slate-500">E_inst/site</span>
            <span className="font-mono text-slate-300">{Number(summary._energy_inst).toFixed(4)}</span>
          </div>
        )}
      </div>
      {!allPassed && (
        <p className="mt-2 text-[9px] text-red-400/80">
          One or more SSE algorithm checks failed — results may be unreliable.
        </p>
      )}
    </div>
  );
};

// ─── ANSI strip + terminal colour coding ─────────────────────────────────────
// Strip ANSI escape sequences (added by C++ debug output) before display.
// The web terminal doesn't interpret them natively, so we strip and let the
// prefix-based CSS colouring below provide the visual hierarchy.
const stripAnsi = (s: string) => s.replace(/\x1b\[[0-9;]*m/g, '');

const termLineClass = (line: string): string => {
  if (!line.trim()) return 'text-transparent select-none';
  // Debug sections — matched BEFORE the generic '[' check below
  if (line.startsWith('[DBG-INIT]')) return 'text-blue-400';
  if (line.startsWith('[DBG-OPS]'))  return 'text-emerald-400';
  if (line.startsWith('[DBG-LIST]')) return 'text-violet-400';
  if (line.startsWith('[DBG-LOOP]')) return 'text-rose-400';
  if (line.startsWith('[DBG-WARN]')) return 'text-amber-400';
  // Standard SSE output
  if (line.startsWith('[PROGRESS]'))           return 'text-slate-600';
  if (line.includes('[SSE-CHECK] PASS'))       return 'text-emerald-400';
  if (line.includes('[SSE-CHECK] FAIL'))       return 'text-red-400';
  if (line.startsWith('[SSE-CHECKS-SUMMARY]')) return 'text-cyan-300';
  if (line.startsWith('['))                    return 'text-amber-300';
  if (line.startsWith('#'))                    return 'text-slate-500';
  return 'text-slate-300';
};

const TerminalPanel = ({ lines, isLive }: { lines: string[]; isLive: boolean }) => {
  const scrollRef = React.useRef<HTMLDivElement>(null);
  const [autoScroll, setAutoScroll] = React.useState(true);
  const [copied, setCopied] = React.useState(false);

  const handleCopy = () => {
    navigator.clipboard?.writeText(lines.join('\n')).then(() => {
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    });
  };

  React.useEffect(() => {
    if (autoScroll && scrollRef.current)
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
  }, [lines, autoScroll]);

  return (
    <div className="bg-[#06060a] border border-slate-800/60 rounded-2xl overflow-hidden shadow-2xl relative">
      {/* top accent */}
      <div className="absolute top-0 left-0 right-0 h-[1px] bg-gradient-to-r from-transparent via-emerald-500/25 to-transparent" />

      {/* title bar */}
      <div className="flex items-center justify-between px-4 py-2.5 border-b border-slate-800/50 bg-slate-950/60">
        <div className="flex items-center gap-3">
          {/* macOS traffic lights */}
          <div className="flex gap-1.5">
            <div className="w-2.5 h-2.5 rounded-full bg-red-500/60" />
            <div className="w-2.5 h-2.5 rounded-full bg-amber-500/60" />
            <div className="w-2.5 h-2.5 rounded-full bg-emerald-500/60" />
          </div>
          <Terminal className="w-3 h-3 text-slate-500" />
          <span className="text-[10px] font-bold text-slate-400 tracking-widest uppercase font-mono">
            sse_sim_test — stdout
          </span>
          {isLive && (
            <span className="flex items-center gap-1.5 text-[9px] font-mono text-emerald-400">
              <span className="relative flex h-1.5 w-1.5">
                <span className="animate-ping absolute inline-flex h-full w-full rounded-full bg-emerald-400 opacity-75" />
                <span className="relative inline-flex rounded-full h-1.5 w-1.5 bg-emerald-400" />
              </span>
              LIVE
            </span>
          )}
        </div>

        <div className="flex items-center gap-2">
          <span className="text-[9px] text-slate-600 font-mono">{lines.length} lines</span>
          <button
            onClick={() => setAutoScroll(a => !a)}
            title="Toggle auto-scroll"
            className={`text-[8px] font-bold px-1.5 py-0.5 rounded border transition-colors ${
              autoScroll
                ? 'text-emerald-400 border-emerald-500/30 bg-emerald-950/30'
                : 'text-slate-600 border-slate-700/30 bg-slate-900/30 hover:text-slate-300'
            }`}
          >
            AUTO↓
          </button>
          {lines.length > 0 && (
            <button
              onClick={handleCopy}
              title="Copy all output"
              className={`flex items-center gap-1 text-[8px] font-bold px-1.5 py-0.5 rounded border transition-all duration-200 ${
                copied
                  ? 'text-emerald-400 border-emerald-500/40 bg-emerald-950/30'
                  : 'text-slate-500 border-slate-700/30 bg-slate-900/30 hover:text-slate-300'
              }`}
            >
              {copied ? (
                <>
                  <Check className="w-2.5 h-2.5" />
                  COPIED
                </>
              ) : (
                <>
                  <Copy className="w-2.5 h-2.5" />
                  COPY
                </>
              )}
            </button>
          )}
        </div>
      </div>

      {/* body */}
      <div
        ref={scrollRef}
        onScroll={() => {
          if (!scrollRef.current) return;
          const { scrollTop, scrollHeight, clientHeight } = scrollRef.current;
          setAutoScroll(scrollHeight - scrollTop - clientHeight < 24);
        }}
        className="h-[calc(100vh-155px)] overflow-y-auto p-3 custom-scrollbar"
        style={{ fontFamily: "'Cascadia Code','Fira Code','Consolas','Courier New',monospace" }}
      >
        {lines.length === 0 ? (
          <div className="flex items-center gap-2 justify-center h-full text-slate-700 text-[11px] font-mono">
            <span>Waiting for output</span>
            <span className="animate-pulse text-emerald-700">▊</span>
          </div>
        ) : (
          lines.map((line, i) => (
            <div key={i} className={`flex items-baseline text-[10px] leading-[1.65] ${termLineClass(stripAnsi(line))}`}>
              <span className="text-slate-800 mr-3 select-none tabular-nums w-7 text-right flex-shrink-0 font-mono">
                {i + 1}
              </span>
              <span className="break-all whitespace-pre-wrap">{stripAnsi(line) || ' '}</span>
            </div>
          ))
        )}
        {isLive && (
          <div className="flex items-baseline text-[10px] leading-[1.65] text-emerald-600 mt-0.5">
            <span className="text-slate-800 mr-3 select-none w-7 text-right flex-shrink-0 font-mono">$</span>
            <span className="animate-pulse">▊</span>
          </div>
        )}
      </div>
    </div>
  );
};

export default function SimulationDashboard({ isServerReady = true }: { isServerReady?: boolean }) {
  const getSavedParams = (key: string, defaultParams: any) => {
    if (typeof window !== 'undefined') {
      const saved = localStorage.getItem(key);
      if (saved) {
        try {
          return { ...defaultParams, ...JSON.parse(saved) };
        } catch (e) {}
      }
    }
    return defaultParams;
  };

  const defaultEdParams = { Lx: 4, Ly: 4, Nl: 1, J0: 1.0, J1: 0.3, J2: 0.025, h: 0.0, beta: 1.0, n: 20, sparse: true };
  const defaultQmcParams = { Lx: 4, Ly: 4, Lz: 1, J0: 1.0, J1: 0.3, J2: 0.025, J3: 0.0, hx: 0.0, beta: 1.0, n_therm: 50000, n_measure: 50000 };

  const [edParams, setEdParams] = useState(() => getSavedParams('edParams_v1', defaultEdParams));
  const [qmcParams, setQmcParams] = useState(() => getSavedParams('qmcParams_v1', defaultQmcParams));

  useEffect(() => {
    localStorage.setItem('edParams_v1', JSON.stringify(edParams));
  }, [edParams]);

  useEffect(() => {
    localStorage.setItem('qmcParams_v1', JSON.stringify(qmcParams));
  }, [qmcParams]);

  const [simType, setSimType] = useState<string>(() => localStorage.getItem('simType') ?? 'ED');

  type PhaseKey = 'therm' | 'measure' | 'hamiltonian' | 'diagonalization' | 'observables' | 'compile';
  type PhaseState = {
    eta: number;        // original calibrated estimate (s)
    deadline: number;   // epoch ms when this phase will finish
    start: number;      // epoch ms when this phase started
    displayEta: number; // live countdown (s), updated by interval
    progress: number;   // 0–98 synthetic from elapsed / eta
    done: boolean;
  };
  type Job = {
    id: string;
    type: string;
    params: any;
    startTime: number;
    currentPhase: PhaseKey | null;
    phases: Partial<Record<PhaseKey, PhaseState>>;
    checks: Record<string, { total: number; passed: number }>;
  };
  const [runningJobs, setRunningJobs] = useState<Job[]>([]);
  const wsRefs = useRef<Map<string, WebSocket>>(new Map());
  const cancelledJobs = useRef<Set<string>>(new Set());

  // ─── History State ──────────────────────────────────────────────────────────
  // @ts-ignore
  const [historyData, setHistoryData] = useState<any[]>([]);
  // @ts-ignore
  const [filteredHistory, setFilteredHistory] = useState<any[]>([]);
  // @ts-ignore
  const [startDate, setStartDate] = useState<string>('');
  // @ts-ignore
  const [endDate, setEndDate] = useState<string>('');
  // @ts-ignore
  const [historyError, setHistoryError] = useState<string | null>(null);

  const fetchHistory = async () => {
    try {
      setHistoryError(null);
      const res = await axios.get(`${API_URL}/history`);
      if (res.data && res.data.status === 'success') {
        const data = res.data.data;
        setHistoryData(data);
        setFilteredHistory(data);
        // Always use the latest from the server — keeps comparison table
        // and spin config snapshots in sync with the most recent QMC run.
        setResults(data.slice(0, 5));
        setSelectedIndex(0);
      } else {
        setHistoryError("Invalid status from server: " + JSON.stringify(res.data));
      }
    } catch (err: any) {
      console.error("Failed to load history:", err);
      setHistoryError(err.message || "Unknown network error");
    }
  };

  useEffect(() => {
    fetchHistory();
  }, []);

  useEffect(() => {
    let filtered = historyData;
    if (startDate) {
      filtered = filtered.filter(item => new Date(item.timestamp) >= new Date(startDate));
    }
    if (endDate) {
      const end = new Date(endDate);
      end.setHours(23, 59, 59, 999);
      filtered = filtered.filter(item => new Date(item.timestamp) <= end);
    }
    setFilteredHistory(filtered);
  }, [historyData, startDate, endDate]);

  // @ts-ignore
  const loadHistoryRun = (historyItem: any) => {
    // If it's already at index 0, do nothing. Else prepend it.
    setResults(prev => [historyItem, ...prev].slice(0, 5));
    setSelectedIndex(0);
  };

  // Results always come from backend history — never from localStorage.
  // This ensures the comparison table and spin config snapshots are always
  // consistent with the same (latest) QMC run.
  const [results, setResults] = useState<any[]>([]);
  const [selectedIndex, setSelectedIndex] = useState(0);

  // ─── Terminal / view state ────────────────────────────────────────────────
  const [terminalLines, setTerminalLines] = useState<string[]>([]);
  const [showTerminal, setShowTerminal] = useState(false);

  useEffect(() => { localStorage.setItem('simType', simType); }, [simType]);

  // ── Ticker: smooth countdown every 200ms ─────────────────────────────────
  // The displayed ETA counts down from the server's deadline estimate.
  // Between server messages it decrements naturally so the UI never freezes.
  // Progress bar uses elapsed / eta (capped at 98% until server says done).
  useEffect(() => {
    const id = setInterval(() => {
      const now = Date.now();
      setRunningJobs(prev => prev.map(job => {
        const updated: typeof job.phases = {};
        let dirty = false;
        for (const [k, ps] of Object.entries(job.phases) as [PhaseKey, PhaseState][]) {
          if (ps.done) { updated[k] = ps; continue; }
          // Smooth countdown: deadline - now (never goes below 0)
          const displayEta = Math.max(0, (ps.deadline - now) / 1000);
          // Progress: use server-reported progress for QMC phases; estimate synthetically for ED phases
          const isQmcPhase = k === 'therm' || k === 'measure';
          const progress   = isQmcPhase
            ? ps.progress
            : (ps.eta > 0 ? Math.min(98, ((now - ps.start) / 1000) / ps.eta * 100) : 0);
          updated[k] = { ...ps, displayEta, progress };
          dirty = true;
        }
        return dirty ? { ...job, phases: updated } : job;
      }));
    }, 200);  // 200ms for silky-smooth countdown
    return () => clearInterval(id);
  }, []);

  const handleStop = async (jobId: string) => {
    cancelledJobs.current.add(jobId);
    try { await axios.post(`${API_URL}/cancel/${jobId}`); } catch { /* ignore */ }
    const ws = wsRefs.current.get(jobId);
    if (ws) { ws.close(); wsRefs.current.delete(jobId); }
    setRunningJobs(prev => prev.filter(j => j.id !== jobId));
  };

  const handleRun = async (overrideType?: 'ED' | 'QMC', overrideParams?: any) => {
    const activeSimType = overrideType || simType;
    const currentParams = overrideParams || (activeSimType === 'ED' ? edParams : qmcParams);
    
    // 1. Prepare job payload and UI state
    const jobId = Math.random().toString(36).substring(7);
    const newJob: Job = {
      id: jobId,
      type: activeSimType,
      params: { ...currentParams },
      startTime: Date.now(),
      currentPhase: null,
      phases: {},
      checks: {},
    };
    setRunningJobs(prev => [...prev, newJob]);
    setTerminalLines([]);  // clear terminal output for the new job
    // 2. Open WebSocket for progress updates (best-effort — don't block the HTTP call)
    const ws = new WebSocket(`${WS_URL}/${jobId}`);
    wsRefs.current.set(jobId, ws);

    ws.onerror = (e) => {
      console.warn('[WS] Connection error for job', jobId, e);
      // Don't cancel the job — the HTTP POST still runs and returns results.
    };

    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);

      // Raw terminal output — handled separately from job state
      if (data.type === 'stdout') {
        setTerminalLines(prev => [...prev, data.line as string].slice(-1000));
        return;
      }

      setRunningJobs(prev => prev.map(job => {
        if (job.id !== jobId) return job;
        
        if (data.type === 'time_estimate' || data.type === 'progress') {
          const phaseName = data.phase as PhaseKey;
          const now = Date.now();
          const updatedPhases: typeof job.phases = { ...job.phases };
          
          // Mark previous phase done ONLY if it's a different phase
          if (job.currentPhase && job.currentPhase !== phaseName && updatedPhases[job.currentPhase]) {
            updatedPhases[job.currentPhase] = {
              ...updatedPhases[job.currentPhase]!,
              done: true, progress: 100, displayEta: 0,
            };
          }
          
          if (data.type === 'progress') {
            const pct = (data.current / data.total) * 100;
            const prevPs = updatedPhases[phaseName];

            // ── tqdm-style EMA blend on incoming server estimate ──
            // alpha=0.3 (tqdm default): responsive but stable.
            // We blend into the *remaining* ETA to avoid jumps.
            const EMA_ALPHA = 0.3;
            const serverEta = data.eta as number;
            let blendedEta: number;
            if (!prevPs || prevPs.eta <= 0) {
              blendedEta = serverEta;  // warm-start: accept first estimate directly
            } else {
              // Count down prevPs by elapsed time, then blend toward server value
              const prevRemaining = Math.max(0, (prevPs.deadline - now) / 1000);
              blendedEta = EMA_ALPHA * serverEta + (1 - EMA_ALPHA) * prevRemaining;
            }

            updatedPhases[phaseName] = {
              eta:        blendedEta,
              deadline:   now + blendedEta * 1000,   // deadline drives the countdown
              start:      prevPs?.start || now,
              displayEta: blendedEta,
              progress:   pct,
              done:       false,
            };
          } else {
            // Original time_estimate logic
            updatedPhases[phaseName] = {
              eta: data.eta,
              deadline: now + data.eta * 1000,
              start: now,
              displayEta: data.eta,
              progress: 0,
              done: false,
            };
          }
          return { ...job, currentPhase: phaseName, phases: updatedPhases };
        } else if (data.type === 'check_result') {
          const name = data.name as string;
          const prev_c = job.checks[name] || { total: 0, passed: 0 };
          return {
            ...job,
            checks: {
              ...job.checks,
              [name]: { total: prev_c.total + 1, passed: prev_c.passed + (data.passed ? 1 : 0) },
            },
          };
        } else if (data.type === 'complete') {
          // Mark the final active phase done.
          const updatedPhases: typeof job.phases = { ...job.phases };
          if (job.currentPhase && updatedPhases[job.currentPhase]) {
            updatedPhases[job.currentPhase] = {
              ...updatedPhases[job.currentPhase]!,
              done: true, progress: 100, displayEta: 0,
            };
          }
          return { ...job, currentPhase: null, phases: updatedPhases };
        }
        return job;
      }));
    };

    // 3. Launch API request asynchronously (don't block the UI)
    try {
      const endpoint = activeSimType === 'ED' ? '/run-ed' : '/run-qmc';
      const payload = { ...currentParams, client_id: jobId };
      const res = await axios.post(`${API_URL}${endpoint}`, payload);
      const runData = res.data.data;
      
      ws.close();
      wsRefs.current.delete(jobId);
      setRunningJobs(prev => prev.filter(j => j.id !== jobId));

      // Ignore result if the user clicked Stop while the job was running
      if (cancelledJobs.current.has(jobId)) {
        cancelledJobs.current.delete(jobId);
        return;
      }

      if (runData && runData.length > 0) {
        let { ed, qmc } = runData[0];

        if (ed && qmc && !edQmcParamsMatch(ed, qmc)) {
          // Params don't match — keep only the primary simulation result
          if (activeSimType === 'ED') qmc = null;
          else ed = null;
        }

        // Immediately show new result at top, then refresh from backend
        setResults(prev => [{ ed, qmc }, ...prev].slice(0, 5));
        setSelectedIndex(0);
        // Sync with backend so spin configs + comparison table stay consistent
        await fetchHistory();
      }
    } catch (error) {
      console.error("Simulation error:", error);
      ws.close();
      wsRefs.current.delete(jobId);
      cancelledJobs.current.delete(jobId);
      setRunningJobs(prev => prev.filter(j => j.id !== jobId));
    }
  };

  const formatScientific = (v: number, precision: number): string => {
    if (v === 0) return '0';
    const str = v.toExponential(precision);
    const [base, expStr] = str.split('e');
    const exp = parseInt(expStr, 10);
    if (exp === 0) return base;
    const supMap: Record<string, string> = {
      '-': '⁻', '+': '⁺', '0': '⁰', '1': '¹', '2': '²', 
      '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'
    };
    const unicodeExp = exp.toString().split('').map(c => supMap[c]).join('');
    return `${base} × 10${unicodeExp}`;
  };

  const formatObsVal = (val: number | null | undefined, formatType: 'float' | 'exp') => {
    if (val === null || val === undefined) return '--';
    return formatType === 'float' ? val.toFixed(6) : formatScientific(val, 4);
  };

  const formatEta = (val: any) => {
    if (val === undefined || val === null || isNaN(val) || !isFinite(val)) return '...';
    if (val <= 0) return 'Finishing...';
    if (val >= 60) {
      const m = Math.floor(val / 60);
      const s = Math.floor(val % 60);
      return `${m}m ${s.toString().padStart(2, '0')}s`;
    }
    return `${Number(val).toFixed(1)}s`;
  };

  // Returns true only when the two records were computed for IDENTICAL physics.
  // Used in handleRun (live results) and getObservableRows (display / history).
  const edQmcParamsMatch = (ed: any, qmc: any): boolean => {
    if (!ed || !qmc) return false;
    const edNl  = ed.Nl  ?? ed.Lz ?? 1;
    const qmcLz = qmc.Lz ?? qmc.Nl ?? 1;
    const edH   = ed.h   ?? ed.hx ?? 0;
    const qmcHx = qmc.hx ?? qmc.h ?? 0;
    const qmcJ3 = qmc.J3 ?? 0;
    return (
      ed.Lx === qmc.Lx &&
      ed.Ly === qmc.Ly &&
      edNl  === qmcLz &&
      Math.abs(ed.beta  - qmc.beta) < 1e-6 &&
      Math.abs(ed.J0   - qmc.J0)   < 1e-6 &&
      Math.abs(ed.J1   - qmc.J1)   < 1e-6 &&
      Math.abs(ed.J2   - qmc.J2)   < 1e-6 &&
      Math.abs(edH     - qmcHx)    < 1e-6 &&
      Math.abs(qmcJ3)              < 1e-6
    );
  };

  const getObservableRows = () => {
    if (!results || results.length === 0 || !results[selectedIndex]) return [];
    const item = results[selectedIndex];
    const ed  = item.ed;
    const qmc = item.qmc;

    // Guard: never compare ED and QMC from different parameter sets.
    // If both are present but don't match, show values for each column but
    // suppress the Δ and Agreement columns (they'd be meaningless).
    const bothPresent = !!ed && !!qmc;
    const matched     = bothPresent ? edQmcParamsMatch(ed, qmc) : true;

    const observablesList = [
      { name: "Energy ⟨E⟩",          edKey: "enrg",        qmcKey: "enrg",        format: "float" },
      { name: "⟨E²⟩",                edKey: "enrg2",       qmcKey: "enrg2",       format: "float" },
      { name: "⟨E⁴⟩",                edKey: "enrg4",       qmcKey: "enrg4",       format: "float" },
      { name: "Stag. Mag² ⟨m²⟩",    edKey: "SMag_square", qmcKey: "SMag_square", format: "exp"   },
      { name: "Stag. Mag⁴ ⟨m⁴⟩",    edKey: "SMag_four",   qmcKey: "SMag_four",   format: "exp"   },
      { name: "Transverse Mag ⟨Mₓ⟩", edKey: "Mag_x",       qmcKey: "Mag_x",       format: "exp"   },
    ];

    return observablesList.map(obs => {
      const edVal  = ed  ? ed[obs.edKey]   : null;
      const qmcVal = qmc ? qmc[obs.qmcKey] : null;

      // Only compute Δ / agreement when both values exist AND params match
      const canCompare = matched && edVal != null && qmcVal != null;
      const delta = canCompare ? Math.abs(qmcVal - edVal) : null;
      const agreement = (canCompare && edVal !== 0)
        ? Math.max(0, 1 - Math.abs(qmcVal - edVal) / Math.abs(edVal)) * 100
        : null;

      return { name: obs.name, edVal, qmcVal, delta, agreement, format: obs.format, matched };
    });
  };



  return (
    <>
    <div className="grid grid-cols-1 lg:grid-cols-12 gap-4 items-start">

      {/* ── Left: Parameters ── */}
      <div className="lg:col-span-4">
        <div className="bg-[#0e0f1e]/60 border border-slate-800/60 rounded-2xl p-4 backdrop-blur-xl shadow-2xl relative overflow-hidden">
          <div className="absolute top-0 left-0 right-0 h-[1px] bg-gradient-to-r from-transparent via-indigo-500/50 to-transparent" />

          <div className="flex items-center justify-between mb-4">
            <h2 className="text-xs font-bold tracking-wide text-white uppercase flex items-center gap-1.5">
              <Sparkles className="w-3.5 h-3.5 text-indigo-400" />
              Parameters
            </h2>
            <div className="bg-slate-950/80 p-0.5 border border-slate-800/80 rounded-lg flex">
              {['ED', 'QMC'].map(type => (
                <button
                  key={type}
                  onClick={() => setSimType(type)}
                  className={`px-3 py-1 text-[10px] font-bold rounded-md transition-all duration-200 ${simType === type ? 'bg-indigo-500/20 text-indigo-300 border border-indigo-500/30' : 'text-slate-500 hover:text-slate-300 border border-transparent'}`}
                >
                  {type}
                </button>
              ))}
            </div>
          </div>

          {simType === 'ED' ? (
            <div className="space-y-2.5">
              <div className="grid grid-cols-3 gap-2">
                <InputField label="Lx" value={edParams.Lx} onChange={(val) => setEdParams({ ...edParams, Lx: val })} />
                <InputField label="Ly" value={edParams.Ly} onChange={(val) => setEdParams({ ...edParams, Ly: val })} />
                <InputField label="Layers" value={edParams.Nl} onChange={(val) => setEdParams({ ...edParams, Nl: val })} />
              </div>
              <div className="grid grid-cols-3 gap-2 pt-2.5 border-t border-slate-800/30">
                <InputField label="J₀" value={edParams.J0} step="0.1" onChange={(val) => setEdParams({ ...edParams, J0: val })} />
                <InputField label="J₁" value={edParams.J1} step="0.1" onChange={(val) => setEdParams({ ...edParams, J1: val })} />
                <InputField label="J₂" value={edParams.J2} step="0.01" onChange={(val) => setEdParams({ ...edParams, J2: val })} />
              </div>
              <div className="grid grid-cols-4 gap-2 pt-2.5 border-t border-slate-800/30">
                <InputField label="Field h" value={edParams.h} step="0.1" onChange={(val) => setEdParams({ ...edParams, h: val })} />
                <InputField label="β = 1/T" value={edParams.beta} step="0.5" onChange={(val) => setEdParams({ ...edParams, beta: val })} />
                <InputField label="Basis n" value={edParams.n} onChange={(val) => setEdParams({ ...edParams, n: val })} />
                <div className="flex flex-col gap-1.5">
                  <label className="text-[11px] font-semibold tracking-wide text-slate-400/80 uppercase">Sparse</label>
                  <div className="flex items-center h-[30px]">
                    <input
                      type="checkbox"
                      checked={edParams.sparse}
                      onChange={(e) => setEdParams({ ...edParams, sparse: e.target.checked })}
                      className="w-4 h-4 rounded border-slate-800 bg-slate-950 outline-none cursor-pointer accent-indigo-500"
                    />
                  </div>
                </div>
              </div>
            </div>
          ) : (
            <div className="space-y-2.5">
              <div className="grid grid-cols-3 gap-2">
                <InputField label="Lx" value={qmcParams.Lx} onChange={(val) => setQmcParams({ ...qmcParams, Lx: val })} />
                <InputField label="Ly" value={qmcParams.Ly} onChange={(val) => setQmcParams({ ...qmcParams, Ly: val })} />
                <InputField label="Layers" value={qmcParams.Lz} onChange={(val) => setQmcParams({ ...qmcParams, Lz: val })} />
              </div>
              <div className="grid grid-cols-4 gap-2 pt-2.5 border-t border-slate-800/30">
                <InputField label="J₀" value={qmcParams.J0} step="0.1" onChange={(val) => setQmcParams({ ...qmcParams, J0: val })} />
                <InputField label="J₁" value={qmcParams.J1} step="0.1" onChange={(val) => setQmcParams({ ...qmcParams, J1: val })} />
                <InputField label="J₂" value={qmcParams.J2} step="0.01" onChange={(val) => setQmcParams({ ...qmcParams, J2: val })} />
                <InputField label="J₃" value={qmcParams.J3} step="0.01" onChange={(val) => setQmcParams({ ...qmcParams, J3: val })} />
              </div>
              <div className="grid grid-cols-2 gap-2 pt-2.5 border-t border-slate-800/30">
                <InputField label="Field hₓ" value={qmcParams.hx} step="0.1" onChange={(val) => setQmcParams({ ...qmcParams, hx: val })} />
                <InputField label="β = 1/T" value={qmcParams.beta} step="0.5" onChange={(val) => setQmcParams({ ...qmcParams, beta: val })} />
              </div>
              <div className="grid grid-cols-2 gap-2">
                <InputField label="Therm sweeps" value={qmcParams.n_therm} step="5000" onChange={(val) => setQmcParams({ ...qmcParams, n_therm: val })} />
                <InputField label="Meas sweeps" value={qmcParams.n_measure} step="5000" onChange={(val) => setQmcParams({ ...qmcParams, n_measure: val })} />
              </div>
            </div>
          )}

          <button
            onClick={() => handleRun()}
            disabled={!isServerReady}
            className={`mt-4 w-full text-white text-xs font-bold py-2.5 rounded-xl flex items-center justify-center gap-2 transition-all ${
              !isServerReady 
                ? "bg-slate-800 text-slate-500 cursor-not-allowed" 
                : "bg-gradient-to-r from-indigo-500 via-purple-500 to-indigo-600 shadow-lg shadow-indigo-500/20 active:scale-[0.98] cursor-pointer hover:brightness-110"
            }`}
          >
            {runningJobs.length > 0 ? (
              <>
                <span className="relative flex h-2 w-2">
                  <span className="animate-ping absolute inline-flex h-full w-full rounded-full bg-white opacity-60" />
                  <span className="relative inline-flex rounded-full h-2 w-2 bg-white" />
                </span>
                RUNNING — QUEUE {simType}
              </>
            ) : (
              <>
                <Play className="w-3.5 h-3.5 fill-current" />
                LAUNCH {simType}
              </>
            )}
          </button>

          {simType === 'QMC' && (
            <button
              onClick={() => {
                setShowTerminal(true);
                handleRun('QMC', { ...qmcParams, debug_mode: true, n_therm: 1, n_measure: 0 });
              }}
              disabled={!isServerReady}
              className={`mt-1.5 w-full text-xs font-bold py-2 rounded-xl flex items-center justify-center gap-2 transition-all border ${
                !isServerReady
                  ? 'border-slate-800 text-slate-600 cursor-not-allowed'
                  : 'border-amber-500/40 text-amber-400 hover:border-amber-400 hover:bg-amber-500/10 active:scale-[0.98] cursor-pointer'
              }`}
            >
              <Bug className="w-3 h-3" />
              DEBUG RUN
            </button>
          )}

          {/* ── View toggle ── */}
          <div className="mt-3 pt-3 border-t border-slate-800/40 flex items-center justify-center gap-3">
            <button
              onClick={() => setShowTerminal(false)}
              className={`flex items-center gap-1.5 text-[10px] font-bold uppercase tracking-widest transition-colors duration-150 ${!showTerminal ? 'text-white' : 'text-slate-500 hover:text-slate-400'}`}
            >
              <TableProperties className="w-3 h-3" />
              Results
            </button>
            <button
              onClick={() => setShowTerminal(t => !t)}
              aria-label="Toggle terminal view"
              className="relative w-11 h-[22px] rounded-full border border-slate-700/60 bg-slate-900 focus:outline-none"
            >
              <motion.div
                animate={{ x: showTerminal ? 22 : 2 }}
                transition={{ type: 'spring', stiffness: 500, damping: 32 }}
                className={`absolute top-[3px] w-4 h-4 rounded-full shadow-md transition-colors duration-200 ${
                  showTerminal
                    ? 'bg-gradient-to-tr from-emerald-500 to-cyan-400'
                    : 'bg-gradient-to-tr from-indigo-500 to-purple-500'
                }`}
              />
            </button>
            <button
              onClick={() => setShowTerminal(true)}
              className={`flex items-center gap-1.5 text-[10px] font-bold uppercase tracking-widest transition-colors duration-150 ${showTerminal ? 'text-white' : 'text-slate-500 hover:text-slate-400'}`}
            >
              <Terminal className="w-3 h-3" />
              Terminal
              {runningJobs.length > 0 && (
                <span className="relative flex h-1.5 w-1.5 ml-0.5">
                  <span className="animate-ping absolute inline-flex h-full w-full rounded-full bg-emerald-400 opacity-75" />
                  <span className="relative inline-flex rounded-full h-1.5 w-1.5 bg-emerald-400" />
                </span>
              )}
            </button>
          </div>

          <AnimatePresence>
            {runningJobs.length > 0 && (
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: 'auto' }}
                exit={{ opacity: 0, height: 0 }}
                className="mt-3 space-y-2 border-t border-slate-800/40 pt-3"
              >
                <span className="text-[9px] font-bold text-slate-600 uppercase tracking-widest">Active Jobs</span>
                {runningJobs.map(job => {
                  // Two phase rows differ by sim type
                  const phaseRows: { key: PhaseKey; label: string; barCls: string }[] =
                    job.type === 'ED'
                      ? [
                          { key: 'hamiltonian',     label: 'BUILD H', barCls: 'from-amber-500 to-orange-400'   },
                          { key: 'diagonalization', label: 'DIAG',    barCls: 'from-cyan-500 to-blue-400'      },
                          { key: 'observables',     label: 'OBS',     barCls: 'from-violet-500 to-purple-400'  },
                        ]
                      : [
                          // Compile row only appears when the backend signals a rebuild is needed
                          ...(job.phases['compile'] || job.currentPhase === 'compile'
                            ? [{ key: 'compile' as PhaseKey, label: 'C++ BUILD', barCls: 'from-rose-500 to-pink-400' }]
                            : []),
                          { key: 'therm',   label: 'THERM',   barCls: 'from-amber-500 to-orange-400' },
                          { key: 'measure', label: 'MEASURE', barCls: 'from-cyan-500 to-blue-400'    },
                        ];

                  return (
                    <div key={job.id} className="bg-slate-900/50 rounded-xl p-2.5 border border-slate-800/50 space-y-1.5">
                      {/* Job header */}
                      <div className="flex items-center gap-1.5">
                        <span className="relative flex h-1.5 w-1.5">
                          <span className="animate-ping absolute inline-flex h-full w-full rounded-full bg-indigo-400 opacity-75" />
                          <span className="relative inline-flex rounded-full h-1.5 w-1.5 bg-indigo-500" />
                        </span>
                        <span className="text-[10px] font-extrabold text-indigo-400">{job.type}</span>
                        <span className="text-[8px] text-slate-600 font-mono">β={job.params.beta} · {job.params.Lx}×{job.params.Ly}</span>
                        <button
                          onClick={() => handleStop(job.id)}
                          title="Stop simulation"
                          className="ml-auto flex items-center gap-0.5 px-1.5 py-0.5 rounded text-[8px] font-bold text-red-400 border border-red-500/30 bg-red-950/30 hover:bg-red-500/20 transition-colors"
                        >
                          <Square className="w-2 h-2 fill-current" />
                          STOP
                        </button>
                      </div>

                      {/* One row per phase */}
                      {phaseRows.map(({ key, label, barCls }) => {
                        const ps = job.phases[key];
                        const isActive = job.currentPhase === key;
                        const isDone   = ps?.done ?? false;
                        const isPending = !ps && !isDone;

                        let etaText = 'waiting';
                        if (isDone)    etaText = '✓ done';
                        else if (ps && !isDone) etaText = formatEta(ps.displayEta);

                        const progress = isDone ? 100 : (ps?.progress ?? 0);

                        const labelCls = isDone
                          ? 'text-emerald-400 border-emerald-500/30 bg-emerald-950/40'
                          : isActive
                          ? (key === 'compile'
                              ? 'text-rose-400 border-rose-500/25 bg-rose-500/10'
                              : key === 'therm' || key === 'hamiltonian'
                              ? 'text-amber-400 border-amber-500/25 bg-amber-500/10'
                              : 'text-cyan-400 border-cyan-500/25 bg-cyan-500/10')
                          : 'text-slate-600 border-slate-700/30 bg-slate-800/30';

                        const barGrad = isDone
                          ? 'from-emerald-500 to-emerald-400'
                          : isPending
                          ? 'from-slate-700 to-slate-700'
                          : barCls;

                        return (
                          <div key={key}>
                            <div className="flex items-center justify-between mb-0.5">
                              <span className={`text-[8px] font-bold px-1.5 py-0.5 rounded border ${labelCls}`}>{label}</span>
                              <span className={`text-[9px] font-mono ${isDone ? 'text-emerald-400' : isActive ? 'text-slate-300' : 'text-slate-600'}`}>
                                {etaText}
                              </span>
                            </div>
                            <div className="w-full bg-slate-950/80 rounded-full h-[3px] overflow-hidden">
                              <motion.div
                                className={`bg-gradient-to-r ${barGrad} h-[3px] rounded-full`}
                                initial={{ width: 0 }}
                                animate={{ width: `${progress}%` }}
                                transition={{ ease: 'linear' }}
                              />
                            </div>
                          </div>
                        );
                      })}
                      {/* Live algorithm checks (appear as first [SSE-CHECK] lines arrive) */}
                      {Object.keys(job.checks).length > 0 && (
                        <div className="flex items-center gap-1.5 pt-1.5 border-t border-slate-800/40 mt-0.5 flex-wrap">
                          <span className="text-[7px] font-bold text-slate-600 uppercase tracking-widest">Checks</span>
                          {Object.entries(job.checks).map(([name, { total, passed }]) => {
                            const ok = passed === total;
                            return (
                              <span key={name} className={`text-[7px] font-bold px-1 py-0.5 rounded border ${
                                ok ? 'text-emerald-400 border-emerald-500/30 bg-emerald-950/30'
                                   : 'text-red-400 border-red-500/30 bg-red-950/30'
                              }`}>
                                {ok ? '✓' : '✗'} {name}
                              </span>
                            );
                          })}
                        </div>
                      )}
                    </div>
                  );
                })}
              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </div>

      {/* ── Right: Analytics + Results ── */}
      <div className="lg:col-span-8">

        {!showTerminal ? (<>
        {/* Results Table */}
        <div className="bg-[#0e0f1e]/60 border border-slate-800/60 rounded-2xl p-4 backdrop-blur-xl shadow-2xl relative overflow-hidden">
          <div className="flex items-center justify-between mb-3 gap-2 flex-wrap">
            <h3 className="text-xs font-bold tracking-wide text-white uppercase flex items-center gap-1.5">
              <TableProperties className="w-3.5 h-3.5 text-emerald-400" />
              {selectedIndex === 0 ? 'Latest Run' : 'Historical Run'}
            </h3>
            <div className="flex items-center gap-1.5 flex-wrap">
              {/* Params for the selected run */}
              {results.length > 0 && results[selectedIndex] && (() => {
                const src = results[selectedIndex].ed ?? results[selectedIndex].qmc;
                if (!src) return null;
                const pills = [
                  { k: 'L',  v: `${src.Lx}×${src.Ly}${src.Lz ? '×' + src.Lz : (src.Nl ? '×' + src.Nl : '')}` },
                  { k: 'β',  v: String(src.beta) },
                  { k: 'J₀', v: String(src.J0) },
                  { k: 'J₁', v: String(src.J1) },
                  { k: 'J₂', v: String(src.J2) },
                  ...(src.h  != null && src.h  !== 0 ? [{ k: 'h',  v: String(src.h)  }] : []),
                  ...(src.hx != null && src.hx !== 0 ? [{ k: 'hₓ', v: String(src.hx) }] : []),
                  ...(src.J3 != null && src.J3 !== 0 ? [{ k: 'J₃', v: String(src.J3) }] : []),
                ];
                return (
                  <div className="flex items-center gap-1 mr-2 flex-wrap">
                    {pills.map(({ k, v }) => (
                      <span key={k} className="inline-flex items-center gap-0.5 text-[9px] font-mono bg-slate-900/60 border border-slate-800/60 rounded px-1.5 py-0.5">
                        <span className="text-slate-500">{k}=</span>
                        <span className="text-slate-300">{v}</span>
                      </span>
                    ))}
                  </div>
                );
              })()}

              {/* Source badges — show which simulations contributed */}
              {results.length > 0 && (() => {
                const r = results[selectedIndex];
                const hasED  = !!r?.ed;
                const hasQMC = !!r?.qmc;
                const fmtTime = (ts?: string) => ts ? new Date(ts).toLocaleString([], { month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' }) : '';
                return (
                  <>
                    <span className={`text-[8px] font-bold px-1.5 py-0.5 rounded border flex items-center gap-1.5 ${
                      hasED
                        ? 'text-indigo-400 border-indigo-500/30 bg-indigo-950/40'
                        : 'text-slate-600 border-slate-700/30 bg-slate-900/30 line-through'
                    }`}>
                      <span>ED</span>
                      {hasED && r.ed?.timestamp && <span className="font-mono font-normal text-[7px] opacity-70 tracking-wider">{fmtTime(r.ed.timestamp)}</span>}
                    </span>
                    <span className={`text-[8px] font-bold px-1.5 py-0.5 rounded border flex items-center gap-1.5 ${
                      hasQMC
                        ? 'text-purple-400 border-purple-500/30 bg-purple-950/40'
                        : 'text-slate-600 border-slate-700/30 bg-slate-900/30 line-through'
                    }`}>
                      <span>QMC</span>
                      {hasQMC && r.qmc?.timestamp && <span className="font-mono font-normal text-[7px] opacity-70 tracking-wider">{fmtTime(r.qmc.timestamp)}</span>}
                    </span>
                  </>
                );
              })()}
              {results.length > 1 && (
                <select
                  value={selectedIndex}
                  onChange={(e) => setSelectedIndex(Number(e.target.value))}
                  className="bg-slate-900 border border-slate-700 text-[10px] text-slate-300 rounded px-2 py-1 outline-none"
                >
                  {results.map((r, idx) => {
                    const ts = r.qmc?.timestamp || r.ed?.timestamp;
                    const dateStr = ts ? new Date(ts).toLocaleString([], { month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' }) : 'Unknown Date';
                    return (
                      <option key={idx} value={idx}>
                        {idx === 0 ? 'Latest' : `Run −${idx}`} · {dateStr}
                      </option>
                    );
                  })}
                </select>
              )}
            </div>
          </div>

          {results.length > 0 ? (
            <div className="overflow-x-auto">
              {/* Missing ED or Mismatch warning */}
              {(() => {
                const r = results[selectedIndex];
                if (!r) return null;

                if (r.ed && r.qmc && !edQmcParamsMatch(r.ed, r.qmc)) {
                  return (
                    <div className="mb-2 px-3 py-1.5 rounded-lg bg-amber-500/10 border border-amber-500/30 text-amber-400 text-[10px] font-semibold flex items-center gap-2">
                      ⚠ ED and QMC results have different parameters — Δ and Agreement are hidden
                    </div>
                  );
                }
                if (r.qmc && !r.ed) {
                  const hasJ3 = r.qmc.J3 && Math.abs(r.qmc.J3) > 1e-6;
                  const ns = (r.qmc.Lx ?? 2) * (r.qmc.Ly ?? 2) * (r.qmc.Lz ?? r.qmc.Nl ?? 1);
                  const edFeasible = ns <= 16;
                  if (hasJ3) {
                    return (
                      <div className="mb-2 px-3 py-1.5 rounded-lg bg-slate-800/50 border border-slate-700/50 text-slate-400 text-[10px] flex items-center gap-2">
                        <span>ℹ️ Exact Diagonalization (ED) is not available for models with J₃ ≠ 0.</span>
                      </div>
                    );
                  } else if (!edFeasible) {
                    return (
                      <div className="mb-2 px-3 py-1.5 rounded-lg bg-slate-800/50 border border-slate-700/50 text-slate-400 text-[10px] flex items-center gap-2">
                        <span>ℹ️ ED is limited to ≤ 16 spins. This {r.qmc.Lx}×{r.qmc.Ly} lattice ({ns} spins) is too large for exact diagonalization.</span>
                      </div>
                    );
                  } else {
                    return (
                      <div className="mb-3 px-4 py-3 rounded-xl bg-indigo-500/10 border border-indigo-500/40 flex items-center justify-between gap-4 shadow-lg shadow-indigo-900/10">
                        <div className="flex items-center gap-3">
                          <div className="w-8 h-8 rounded-lg bg-indigo-500/20 border border-indigo-500/30 flex items-center justify-center flex-shrink-0">
                            <Sparkles className="w-4 h-4 text-indigo-400" />
                          </div>
                          <div>
                            <p className="text-[11px] font-semibold text-indigo-300">No ED comparison cached for these parameters</p>
                            <p className="text-[9px] text-indigo-400/70 mt-0.5">
                              {r.qmc.Lx}×{r.qmc.Ly}, β={r.qmc.beta}, J₀={r.qmc.J0}, J₁={r.qmc.J1}, J₂={r.qmc.J2}
                              {(r.qmc.hx ?? 0) !== 0 ? `, hₓ=${r.qmc.hx}` : ''} — {ns} spins (ED feasible)
                            </p>
                          </div>
                        </div>
                        <button
                          onClick={() => {
                            const edP = {
                              Lx: r.qmc.Lx, Ly: r.qmc.Ly, Nl: r.qmc.Lz ?? r.qmc.Nl ?? 1,
                              J0: r.qmc.J0, J1: r.qmc.J1, J2: r.qmc.J2,
                              h: r.qmc.hx ?? r.qmc.h ?? 0,
                              beta: r.qmc.beta,
                              n: 16, sparse: false,
                            };
                            setEdParams(edP as any);
                            setSimType('ED');
                            handleRun('ED', edP);
                          }}
                          className="flex-shrink-0 bg-indigo-600 hover:bg-indigo-500 active:scale-95 text-white text-[10px] font-bold py-1.5 px-4 rounded-lg shadow-md transition-all whitespace-nowrap"
                        >
                          ⚡ Run ED Now
                        </button>
                      </div>
                    );
                  }
                }
                return null;
              })()}
              <table className="w-full text-left border-collapse">
                <thead>
                  <tr className="border-b border-slate-800/60 text-[9px] font-bold text-slate-500 uppercase tracking-wider">
                    <th className="py-2 px-2.5">Observable</th>
                    <th className={`py-2 px-2.5 ${results[selectedIndex]?.ed ? 'text-indigo-300/80' : 'text-slate-700'}`}>Exact ED</th>
                    <th className={`py-2 px-2.5 ${results[selectedIndex]?.qmc ? 'text-purple-300/80' : 'text-slate-700'}`}>QMC SSE</th>
                    <th className="py-2 px-2.5">Δ</th>
                    <th className="py-2 px-2.5 text-right">Agreement</th>
                  </tr>
                </thead>
                <tbody className="text-xs font-mono">
                  {getObservableRows().map((row: any, i: number) => {
                    const agreementColor = row.agreement >= 99.9 ? 'text-emerald-400' : row.agreement >= 99.0 ? 'text-indigo-400' : 'text-amber-400';
                    return (
                      <tr key={i} className="border-b border-slate-800/20 hover:bg-slate-900/30 transition-colors">
                        <td className="py-1.5 px-2.5 font-semibold text-slate-300 font-sans text-[11px]">{row.name}</td>
                        <td className="py-1.5 px-2.5 text-indigo-400">{formatObsVal(row.edVal, row.format)}</td>
                        <td className="py-1.5 px-2.5 text-purple-400">{formatObsVal(row.qmcVal, row.format)}</td>
                        <td className="py-1.5 px-2.5 text-slate-500 text-[10px]">{row.delta !== null ? formatScientific(row.delta, 4) : '--'}</td>
                        <td className={`py-1.5 px-2.5 text-right font-bold text-[11px] ${row.agreement !== null ? agreementColor : 'text-slate-600'}`}>
                          {row.agreement !== null ? `${row.agreement.toFixed(3)}%` : '--'}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          ) : (
            <div className="py-8 flex flex-col items-center gap-2 text-center">
              <div className="w-9 h-9 rounded-full bg-slate-900/80 border border-slate-800/60 flex items-center justify-center">
                <TableProperties className="w-3.5 h-3.5 text-slate-600" />
              </div>
              <div>
                <p className="text-xs font-semibold text-slate-500">No results yet</p>
                <p className="text-[10px] text-slate-600 mt-0.5">Launch a simulation to see observables</p>
              </div>
            </div>
          )}
        </div>

        {/* Algorithm Checks Panel (QMC only — appears when SSE_DEBUG checks ran) */}
        {results.length > 0 && (() => {
          const summary = results[selectedIndex]?.qmc?.checks_summary;
          if (!summary) return null;
          const checkNames = Object.keys(summary).filter(k => !k.startsWith('_'));
          if (checkNames.length === 0) return null;
          return <AlgorithmChecksPanel summary={summary} />;
        })()}

        {/* Lattice snapshots (QMC only) */}
        {results.length > 0 && results[selectedIndex]?.qmc?.spin_configs?.length > 0 && (
          <LatticeViewer
            spinConfigs={results[selectedIndex].qmc.spin_configs}
            Lx={results[selectedIndex].qmc.Lx}
            Ly={results[selectedIndex].qmc.Ly}
            Lz={results[selectedIndex].qmc.Lz ?? results[selectedIndex].qmc.Nl ?? 1}
          />
        )}

        {/* Run History UI has been removed as requested */}
        </>) : (
          <TerminalPanel lines={terminalLines} isLive={runningJobs.length > 0} />
        )}

      </div>

    </div>
    </>
  );
}

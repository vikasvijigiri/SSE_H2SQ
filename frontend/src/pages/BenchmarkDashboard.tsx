// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
import { useEffect, useState } from 'react';
import axios from 'axios';
import { ComposedChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ErrorBar, ReferenceLine, BarChart, Bar } from 'recharts';
import { Play, Square, Download, Activity, ChevronDown, FileText } from 'lucide-react';

const API = 'http://127.0.0.1:8000/api/benchmarks';
const colors: Record<string, string> = { vertex: '#0072B2', premarked: '#D55E00', canonical: '#009E73', tempered: '#CC79A7' };
const names: Record<string, string> = { vertex: 'Per-vertex pairing', premarked: 'Spatial premarking', canonical: 'Premarking + worldline MH', tempered: 'Tempered premarking + worldline MH' };
const axisNames: Record<string, string> = { beta: 'Inverse temperature β', J1: 'Diagonal coupling J₁', J2: 'Dipole coupling J₂', hx: 'Transverse field hₓ' };
const defaultAxes = { beta: '0.5, 2, 8', J1: '0, 0.15, 0.3', J2: '0, 0.025, 0.05', hx: '0, 0.2, 0.5' };
interface Observation { mean: number | null; error: number | null; tau_int?: number | null; tau_error?: number | null; ess?: number | null; ess_per_second?: number | null; rhat?: number | null; difference?: number | null; pull?: number | null; status: string; flags?: string[]; acf?: { lag: number; value: number }[] }
interface Method { validation_status?: string; observables: Record<string, Observation>; update_seconds: Record<string, number>; runs: { seed: number; non_global_moves: number; cutoff_saturations: number; canonical_attempts: number; canonical_accepted: number; replicas?: number; round_trips?: number; beta_ladder?: number[]; swap_attempts?: number[]; swap_accepted?: number[] }[]; block_size: number; retained_samples: number; total_seconds: number }
interface Point { axis: string; x: number; point_id: string; params: Record<string, number>; ed: { values: Record<string, number | null>; numerical_residual: number; reference: string; dimension: number }; methods: Record<string, Method>; reused?: boolean }
interface Job { id: string; created: string; status: string; request: Record<string, unknown>; completed_runs: number; total_runs: number; message: string; points: Point[]; observables: Record<string, [string, string]>; definitions: Record<string, string>; provenance: { model_version: string; revision: string; compiler: string; source_sha256: Record<string, string> }; citations: { title: string; url: string }[] }
type JobRow = Pick<Job, 'id' | 'created' | 'status' | 'completed_runs' | 'total_runs' | 'request'>;
const fmt = (value: number | null | undefined, digits = 5) => value == null ? 'Unresolved' : Math.abs(value) > 0 && (Math.abs(value) < .001 || Math.abs(value) > 10000) ? value.toExponential(3) : value.toPrecision(digits);
const box = 'rounded-2xl border border-slate-800 bg-[#101222] p-5';
const input = 'w-full rounded-lg border border-slate-700 bg-slate-950 px-3 py-2 text-sm text-slate-100';
const chartStyle = { fontFamily: "'Times New Roman', serif", color: '#16202a' };
const errorMessage = (error: unknown) => axios.isAxiosError(error) ? typeof error.response?.data?.detail === 'string' ? error.response.data.detail : JSON.stringify(error.response?.data?.detail || error.message) : String(error);

function ScientificPanel({ pointRows, job, axis, observable, quantity = 'mean' }: { pointRows: Point[]; job: Job; axis: string; observable: string; quantity?: 'mean' | 'tau_int' | 'ess_per_second' | 'difference' }) {
  const data = pointRows.map(point => {
    const row: Record<string, number | null> = { x: point.x, ED_error: 0, ED: quantity === 'mean' ? point.ed.values[observable] : quantity === 'difference' ? 0 : null };
    for (const method of Object.keys(names)) {
      const value = point.methods[method]?.observables[observable];
      row[method] = value?.[quantity] ?? null;
      row[`${method}_error`] = ['mean', 'difference'].includes(quantity) ? value?.error ?? null : quantity === 'tau_int' ? value?.tau_error ?? null : null;
    }
    return row;
  });
  const ylabel = quantity === 'mean' ? job.observables[observable]?.[1] : quantity === 'difference' ? 'QMC ? ED' : quantity === 'tau_int' ? 'τint (complete sweeps)' : 'Effective samples / second';
  return <div className="rounded-xl bg-white p-4" style={chartStyle}>
    <div className="flex items-start justify-between gap-3 mb-1">
      <h3 className="font-semibold text-base">{quantity === 'mean' ? job.observables[observable]?.[0] : ylabel}</h3>
      {<div className="flex gap-2 text-xs font-sans text-blue-700">{['svg', 'pdf', 'png'].map(format => <a key={format} href={`${API}/${job.id}/export?format=${format}&axis=${axis}&observable=${observable}&quantity=${quantity}`}>{format.toUpperCase()}</a>)}</div>}
    </div>
    <ResponsiveContainer width="100%" height={300}>
      <ComposedChart data={data} margin={{ top: 12, right: 18, left: 12, bottom: 25 }}>
        <CartesianGrid stroke="#e6e9ed" vertical={false} />
        <XAxis dataKey="x" type="number" domain={['dataMin', 'dataMax']} label={{ value: axisNames[axis], position: 'bottom', offset: 6 }} tick={{ fontSize: 13 }} />
        <YAxis width={72} domain={['auto', 'auto']} tickFormatter={value => fmt(value, 3)} tick={{ fontSize: 12 }} />
        <Tooltip formatter={(value, name) => [typeof value === 'number' ? fmt(value, 7) : 'Unavailable', name]} labelFormatter={value => `${axis} = ${value}`} />
        <Legend verticalAlign="top" wrapperStyle={{ fontSize: 11 }} />
        {['mean', 'difference'].includes(quantity) && <Line dataKey="ED" name="Complete ED · SE = 0" stroke="#111827" strokeDasharray="5 4" strokeWidth={1.5} dot={{ r: 4, fill: '#111827' }} isAnimationActive={false} connectNulls={false}><ErrorBar dataKey="ED_error" width={4} /></Line>}
        {Object.keys(names).map(method => <Line key={method} dataKey={method} name={names[method]} stroke={colors[method]} strokeWidth={1.5} dot={{ r: 3.5 }} isAnimationActive={false} connectNulls={false}>{quantity !== 'ess_per_second' && <ErrorBar dataKey={`${method}_error`} width={5} strokeWidth={1.4} />}</Line>)}
      </ComposedChart>
    </ResponsiveContainer>
    <p className="text-xs text-slate-500 font-sans">{ylabel} · {['mean', 'difference'].includes(quantity) ? 'QMC bars: 1 SE. ED sampling error: zero. Missing QMC bars require review.' : quantity === 'tau_int' ? 'τint = ½ + Σρ(t). Uncertainty is an asymptotic estimate.' : 'Uses summed chain wall time including thermalization/I/O. Compare validated points; unseen sectors can give misleadingly short ?.'}</p>
  </div>;
}

export default function BenchmarkDashboard() {
  const [jobs, setJobs] = useState<JobRow[]>([]);
  const [job, setJob] = useState<Job | null>(null);
  const [selected, setSelected] = useState('');
  const [axis, setAxis] = useState('beta');
  const [pointIndex, setPointIndex] = useState(0);
  const [observable, setObservable] = useState('E');
  const [panels, setPanels] = useState(['E', 'mx', 'stripe2', 'rho', 'cv', 'polar2']);
  const [axes, setAxes] = useState<Record<string, string>>(defaultAxes);
  const [methods, setMethods] = useState(Object.keys(names));
  const [parameters, setParameters] = useState<Record<string, number>>({ Lx: 4, Ly: 2, J0: 1, beta: 2, J1: .3, J2: .025, hx: .2, chains: 4, n_therm: 50000, n_measure: 105000, seed: 20261016 });
  const [error, setError] = useState('');
  const [busy, setBusy] = useState(false);
  const [showSetup, setShowSetup] = useState(true);

  useEffect(() => {
    let mounted = true; let loaded = false;
    const refresh = async () => {
      try { const { data } = await axios.get<JobRow[]>(API); if (mounted && !loaded) { setShowSetup(data.length === 0); loaded = true; } if (mounted) { setJobs(data); setSelected(previous => previous || data.find(j => Object.keys(j.request.axes as object).length > 1 && ['complete', 'running', 'queued'].includes(j.status))?.id || data[0]?.id || ''); } }
      catch (error) { if (mounted) setError(errorMessage(error)); }
    };
    void refresh(); const timer = setInterval(refresh, 6000);
    return () => { mounted = false; clearInterval(timer); };
  }, []);
  useEffect(() => {
    if (!selected) return;
    let mounted = true;
    const refresh = async () => {
      try { const { data } = await axios.get<Job>(`${API}/${selected}`); if (mounted) { setJob(data); setError(''); setAxis(previous => Object.keys(data.request.axes as object).includes(previous) ? previous : Object.keys(data.request.axes as object)[0]); } }
      catch (error) { if (mounted) setError(errorMessage(error)); }
    };
    void refresh(); const timer = setInterval(refresh, 3000);
    return () => { mounted = false; clearInterval(timer); };
  }, [selected]);

  async function start() {
    setBusy(true); setError('');
    try {
      const parsedAxes: Record<string, number[]> = {};
      for (const [key, value] of Object.entries(axes)) if (value.trim()) {
        const values = value.split(',').map(v => Number(v.trim()));
        if (values.some(v => !Number.isFinite(v)) || value.split(',').some(v => !v.trim())) throw new Error('Enter comma-separated numeric sweep values; leave an axis empty to omit it.');
        parsedAxes[key] = values;
      }
      const { data } = await axios.post<Job>(API, { ...parameters, axes: parsedAxes, algorithms: methods });
      setJob(data); setSelected(data.id); setShowSetup(false); setPointIndex(0);
      setJobs(previous => [data, ...previous]);
    } catch (error) { setError(errorMessage(error)); } finally { setBusy(false); }
  }
  function preset(classical = false) {
    setParameters({ Lx: 4, Ly: classical ? 4 : 2, J0: 1, beta: classical ? 10 : 2, J1: .3, J2: classical ? 0 : .025, hx: classical ? 0 : .2, chains: 4, n_therm: 50000, n_measure: 105000, seed: 20261016 });
    setAxes(classical ? { beta: '10', J1: '', J2: '', hx: '' } : defaultAxes);
    setMethods(classical ? ['premarked', 'canonical'] : Object.keys(names));
  }
  const comparison = job && job.status === 'complete' ? jobs.find(candidate => candidate.id !== job.id && candidate.status === 'complete' &&
    ['Lx', 'Ly', 'J0', 'beta', 'J1', 'J2', 'hx', 'seed', 'chains', 'n_therm', 'n_measure'].every(key => candidate.request[key] === job.request[key]) &&
    (candidate.request.algorithms as string[]).every(method => !(job.request.algorithms as string[]).includes(method))) : null;
  async function compare() {
    if (!job || !comparison) return;
    try { const { data } = await axios.post<Job>(`${API}/compare`, [job.id, comparison.id]); setJobs(previous => [data, ...previous]); setSelected(data.id); setJob(data); }
    catch (error) { setError(errorMessage(error)); }
  }
  const rows = job?.points.filter(p => p.axis === axis).sort((a, b) => a.x - b.x) || [];
  const point = rows[Math.min(pointIndex, Math.max(0, rows.length - 1))];
  const running = jobs.some(j => ['running', 'queued'].includes(j.status)) || (job && ['running', 'queued'].includes(job.status));
  const acf = point ? Array.from(new Set(Object.values(point.methods).flatMap(m => m.observables[observable]?.acf?.map(v => v.lag) || []))).sort((a, b) => a - b).map(lag => ({ lag, ...Object.fromEntries(Object.entries(point.methods).map(([method, value]) => [method, value.observables[observable]?.acf?.find(v => v.lag === lag)?.value ?? null])) })) : [];
  const stageCosts = point ? ['diagonal', 'cluster', 'gauge', 'canonical'].map(stage => ({ stage, ...Object.fromEntries(Object.entries(point.methods).map(([method, value]) => [method, 1e6 * value.update_seconds[stage] / (Number(job?.request.chains) * Number(job?.request.n_measure))])) })) : [];
  const ctRows = point ? [0, .125, .25, .375, .5].map((fraction, i) => ({ fraction, ED: point.ed.values[`ct${i}`], ...Object.fromEntries(Object.entries(point.methods).flatMap(([method, value]) => [[method, value.observables[`ct${i}`]?.mean], [`${method}_error`, value.observables[`ct${i}`]?.error]])) })) : [];

  return <main className="max-w-[1500px] mx-auto px-5 py-8 relative z-10 space-y-6">
    <div className="flex flex-wrap justify-between gap-4 items-end">
      <div><p className="text-xs tracking-widest text-cyan-400 uppercase mb-2">Finite-system validation</p><h2 className="text-3xl font-semibold">ED × QMC benchmark laboratory</h2><p className="text-slate-400 mt-2 text-sm">Parameter sweeps, measured uncertainty, and update efficiency in one reproducible experiment.</p></div>
      <button className="border border-slate-700 rounded-lg px-4 py-2 text-sm flex gap-2 items-center" onClick={() => setShowSetup(!showSetup)}>Experiment settings <ChevronDown size={16} /></button>
    </div>
    {error && <div role="alert" className="p-4 border border-red-700 rounded-xl bg-red-950/40 text-red-200 break-words">{error}</div>}
    {showSetup && <section className={box}>
      <div className="flex flex-wrap justify-between gap-3 mb-5"><h3 className="font-semibold">Design an experiment</h3><div className="flex gap-3 text-xs"><button className="text-cyan-300" onClick={() => preset()}>Recommended quantum suite · 4×2</button><button className="text-purple-300" onClick={() => preset(true)}>Screenshot reference · 4×4, h = 0</button></div></div>
      <div className="grid grid-cols-3 sm:grid-cols-6 xl:grid-cols-12 gap-3">{Object.entries(parameters).map(([key, value]) => <label key={key} className="text-xs text-slate-400 space-y-1 block">{key}<input aria-label={key} className={input} type="number" step={['J0', 'J1', 'J2', 'hx', 'beta'].includes(key) ? 'any' : '1'} value={value} onChange={event => setParameters(previous => ({ ...previous, [key]: Number(event.target.value) }))} /></label>)}</div>
      <div className="grid sm:grid-cols-2 xl:grid-cols-4 gap-4 mt-5">{Object.entries(axes).map(([key, value]) => <label className="block text-xs text-slate-400 space-y-2" key={key}>{axisNames[key]}<input className={input} aria-label={`${key} sweep values`} value={value} onChange={event => setAxes(previous => ({ ...previous, [key]: event.target.value }))} /></label>)}</div>
      <div className="flex flex-wrap gap-5 mt-5">{Object.entries(names).map(([key, value]) => <label key={key} className="text-sm flex items-center gap-2"><input type="checkbox" checked={methods.includes(key)} onChange={() => setMethods(previous => previous.includes(key) ? previous.filter(m => m !== key) : [...previous, key])} />{value}</label>)}</div>
      <div className="flex flex-wrap justify-between gap-4 mt-5 items-center"><p className="text-xs text-slate-400 max-w-3xl">Each axis varies independently around the base parameters. Identical points share the same measurements. Four independent seeds are recommended. Complete quantum ED: ≤10 spins; zero-field enumeration: ≤16 spins. One layer, J₃ = 0.</p><button disabled={busy || !!running} className="bg-cyan-400 text-slate-950 rounded-lg px-5 py-2.5 font-semibold flex gap-2 items-center disabled:opacity-40" onClick={() => void start()}><Play size={16} />{busy ? 'Starting…' : 'Run benchmark'}</button></div>
    </section>}
    <section className={`${box} flex flex-wrap gap-5 items-center justify-between`}>
      <div className="min-w-64 flex-1"><label className="text-xs text-slate-400 block mb-2">Saved local experiments</label><select aria-label="Saved experiment" className={input} value={selected} onChange={event => { setSelected(event.target.value); setPointIndex(0); }}>{!jobs.length && <option value="">No experiments yet</option>}{jobs.map(j => <option key={j.id} value={j.id}>{new Date(j.created).toLocaleString()} · {String(j.request.Lx)}×{String(j.request.Ly)} · {j.status} · {j.completed_runs}/{j.total_runs} chains</option>)}</select></div>
      {job && <div className="flex-1 min-w-64"><div className="flex justify-between text-sm mb-2"><span className={job.status === 'complete' ? 'text-emerald-400' : 'text-amber-300'}>{job.status.toUpperCase()}</span><span>{job.completed_runs} / {job.total_runs} chains</span></div><div className="h-1.5 rounded-full bg-slate-800 overflow-hidden"><div className="h-full bg-cyan-400 transition-all" style={{ width: `${100 * job.completed_runs / Math.max(1, job.total_runs)}%` }} /></div><p className="text-xs text-slate-400 mt-2">{job.message}</p></div>}
      {job && <div className="flex gap-3 items-center">{comparison && <button className="text-purple-300 text-xs" onClick={() => void compare()}>Compare saved schedules</button>}{['csv', 'json'].map(format => <a className="flex items-center gap-1 text-xs text-cyan-300" key={format} href={`${API}/${job.id}/export?format=${format}`}><Download size={14} />{format.toUpperCase()}</a>)}{['failed', 'cancelled', 'interrupted'].includes(job.status) && <button className="text-cyan-300 text-xs" disabled={!!running} onClick={() => void axios.post<Job>(`${API}/${job.id}/resume`).then(response => setJob(response.data)).catch(error => setError(errorMessage(error)))}>Resume stored chains</button>}{['running', 'queued'].includes(job.status) && <button className="text-red-300 text-xs flex gap-1" onClick={() => void axios.post(`${API}/${job.id}/cancel`).catch(error => setError(errorMessage(error)))}><Square size={13} />Cancel</button>}</div>}
    </section>
    {!job && <section className={`${box} text-center py-12`}><Activity className="mx-auto text-cyan-400 mb-4" /><p>Start a benchmark to populate ED comparisons and measured QMC diagnostics.</p></section>}
    {job && <>
      <section className="rounded-xl bg-cyan-950/30 border border-cyan-900/70 p-4 text-sm text-cyan-100 leading-relaxed">H = −J₀ ΣAₚ + J₁ Σ(σ₀σ₂ + σ₁σ₃) − J₂ ΣPₐ·Pᵦ − hₓ Σσˣ. The B-sublattice dipole sign is folded into the implemented +J₂ term. SSE constant shifts are removed before comparing energy moments. ED has zero sampling uncertainty; QMC bars show one standard error. Flagged points need review before publication.</section>
      <div className="flex flex-wrap gap-2">{Object.keys(job.request.axes as object).map(key => <button key={key} onClick={() => { setAxis(key); setPointIndex(0); }} className={`px-5 py-2 rounded-lg border text-sm ${axis === key ? 'border-cyan-500 bg-cyan-950 text-cyan-200' : 'border-slate-800 text-slate-400'}`}>{axisNames[key]}</button>)}</div>
      <section className={box}><h3 className="font-semibold mb-3">Observable panels</h3><div className="flex flex-wrap gap-x-5 gap-y-2">{Object.entries(job.observables).filter(([key]) => !key.startsWith('ct')).map(([key, [name]]) => <label key={key} className="text-xs text-slate-300 flex gap-2 items-center"><input type="checkbox" checked={panels.includes(key)} onChange={() => setPanels(previous => previous.includes(key) ? previous.filter(p => p !== key) : [...previous, key])} />{name}</label>)}</div></section>
      <p className="text-xs text-slate-400">{String(job.request.Lx)}×{String(job.request.Ly)} periodic · Fixed: {['beta', 'J0', 'J1', 'J2', 'hx'].filter(k => k !== axis).map(k => `${k}=${job.request[k]}`).join(' · ')} · ED reference includes every basis state.</p>
      <div className="grid lg:grid-cols-2 gap-5">{panels.map(key => <ScientificPanel key={key} job={job} axis={axis} observable={key} pointRows={rows} />)}</div>
      <section className={box}>
        <div className="flex flex-wrap justify-between gap-4 mb-5"><div><h3 className="font-semibold">Mixing and update efficiency</h3><p className="text-xs text-slate-400 mt-1">Autocorrelation belongs to a complete sweep schedule. The cost of every update stage is measured separately.</p></div><div className="flex gap-3"><select aria-label="Diagnostic observable" className={input} value={observable} onChange={event => setObservable(event.target.value)}>{Object.entries(job.observables).map(([key, [name]]) => <option key={key} value={key}>{name}</option>)}</select><select aria-label="Diagnostic parameter point" className={input} value={Math.min(pointIndex, Math.max(0, rows.length - 1))} onChange={event => setPointIndex(Number(event.target.value))}>{rows.map((p, index) => <option key={index} value={index}>{axis}={p.x}</option>)}</select></div></div>
        <div className="grid lg:grid-cols-2 gap-5"><ScientificPanel job={job} axis={axis} observable={observable} pointRows={rows} quantity="tau_int" /><ScientificPanel job={job} axis={axis} observable={observable} pointRows={rows} quantity="ess_per_second" /></div>
        {point && <><div className="mt-5"><ScientificPanel job={job} axis={axis} observable={observable} pointRows={rows} quantity="difference" /></div><div className="grid lg:grid-cols-2 gap-5 mt-5">
          <div className="bg-white rounded-xl p-4" style={chartStyle}><h4 className="font-semibold mb-3">Autocorrelation · {axis} = {point.x}</h4><ResponsiveContainer width="100%" height={260}><ComposedChart data={acf} margin={{ bottom: 20, right: 20 }}><CartesianGrid vertical={false} stroke="#eee" /><XAxis dataKey="lag" type="number" label={{ value: 'Lag (complete sweeps)', position: 'bottom' }} /><YAxis domain={[-.1, 1]} /><Tooltip /><ReferenceLine y={0} stroke="#999" /><Legend verticalAlign="top" />{Object.keys(names).map(method => <Line key={method} dataKey={method} name={names[method]} stroke={colors[method]} dot={false} isAnimationActive={false} />)}</ComposedChart></ResponsiveContainer></div>
          <div className="bg-white rounded-xl p-4" style={chartStyle}><h4 className="font-semibold mb-3">Update stage cost · μs / sweep</h4><ResponsiveContainer width="100%" height={260}><BarChart data={stageCosts}><CartesianGrid vertical={false} stroke="#eee" /><XAxis dataKey="stage" /><YAxis /><Tooltip /><Legend verticalAlign="top" />{Object.keys(names).map(method => <Bar key={method} dataKey={method} name={names[method]} fill={colors[method]} />)}</BarChart></ResponsiveContainer></div>
        </div>
        <div className="overflow-x-auto mt-5"><table className="w-full text-xs text-left"><thead className="text-slate-400 border-b border-slate-700"><tr>{['Schedule', 'QMC ± SE', 'ED', 'Δ / SE', 'τint ± estimated SE', 'ESS', 'R̂', 'Diagnostics / trace'].map(label => <th key={label} className="p-3 font-medium">{label}</th>)}</tr></thead><tbody>{Object.entries(point.methods).map(([method, summary]) => { const v = summary.observables[observable]; return <tr key={method} className="border-b border-slate-800"><td className="p-3" style={{ color: colors[method] }}>{names[method]}{summary.validation_status === 'review' && <p className="text-amber-300 mt-1">Schedule has flagged observables</p>}</td><td className="p-3 font-mono">{fmt(v?.mean)} ± {fmt(v?.error, 3)}</td><td className="p-3 font-mono">{fmt(point.ed.values[observable])}</td><td className="p-3">{fmt(v?.pull, 3)}</td><td className="p-3">{fmt(v?.tau_int, 4)} ± {fmt(v?.tau_error, 2)}</td><td className="p-3">{fmt(v?.ess, 4)}</td><td className="p-3">{fmt(v?.rhat, 4)}</td><td className="p-3 max-w-72"><span className={v?.status === 'checks passed' || v?.status === 'analytic zero' ? 'text-emerald-400' : 'text-amber-300'}>{v?.status}</span><p className="text-slate-400 mt-1">{v?.flags?.join('; ')}</p><div className="flex gap-2 mt-2">{summary.runs.map((run, chain) => <a className="text-cyan-300 flex gap-1 items-center" key={chain} href={`${API}/${job.id}/traces/${point.point_id}/${method}/${chain}`} title={`Seed ${run.seed}. NumPy structured array in NPZ; allow_pickle=False.`}><FileText size={12} />{chain + 1}</a>)}</div></td></tr>; })}</tbody></table></div>
        <div className="text-xs text-slate-400 flex flex-wrap gap-x-6 gap-y-2 mt-4"><span>ED: {point.ed.reference}, dim = {point.ed.dimension}, max eigenpair residual = {fmt(point.ed.numerical_residual, 3)}.</span>{Object.entries(point.methods).map(([method, m]) => <span key={method}>{names[method]}: block {m.block_size} sweeps; {m.runs.reduce((n, run) => n + run.non_global_moves, 0).toLocaleString()} non-global changes; {m.runs.reduce((n, run) => n + run.cutoff_saturations, 0)} cutoff saturations.{method === 'tempered' && ` ${m.runs[0]?.replicas} temperatures; ${m.runs.reduce((n, run) => n + (run.round_trips || 0), 0)} measured round trips.`}</span>)}</div>
        </>}
      </section>
      {point && <section className={box}><h3 className="font-semibold mb-4">Imaginary-time dipole correlation · {axis} = {point.x}</h3><div className="bg-white rounded-xl p-4" style={chartStyle}><ResponsiveContainer width="100%" height={320}><ComposedChart data={ctRows} margin={{ bottom: 25, right: 20 }}><CartesianGrid vertical={false} stroke="#eee" /><XAxis dataKey="fraction" type="number" label={{ value: 'Imaginary time τ / β', position: 'bottom', offset: 5 }} /><YAxis /><Tooltip /><Legend verticalAlign="top" /><Line dataKey="ED" name="Complete ED · SE = 0" stroke="#111827" strokeDasharray="5 4" isAnimationActive={false} />{Object.keys(names).map(method => <Line key={method} dataKey={method} name={names[method]} stroke={colors[method]} isAnimationActive={false}><ErrorBar dataKey={`${method}_error`} width={5} /></Line>)}</ComposedChart></ResponsiveContainer><p className="text-xs font-sans text-slate-500">c(τ) = Nₚ⁻¹ Σₚ ⟨Pₓ(p,0)Pₓ(p,τ)⟩. Markov-chain lag and imaginary time are different axes.</p></div></section>}
      <details className={box}><summary className="cursor-pointer font-semibold">Definitions, provenance and research sources</summary><div className="mt-4 space-y-3 text-sm text-slate-400">{Object.entries(job.definitions).map(([key, value]) => <p key={key}><strong className="text-slate-200">{key}: </strong>{value}</p>)}<p>Model: {job.provenance.model_version} · source revision: {job.provenance.revision} · {job.provenance.compiler}. Exact source hashes and seeds are included in JSON and raw traces.</p><div className="flex flex-wrap gap-4">{job.citations.map(source => <a className="text-cyan-300 underline" href={source.url} target="_blank" rel="noreferrer" key={source.url}>{source.title}</a>)}</div><p>Passing these diagnostics validates this finite-system experiment. It does not prove ergodicity or ground-state convergence in every parameter regime. Extend β, lattice size and run length before making a thermodynamic claim.</p></div></details>
    </>}
  </main>;
}

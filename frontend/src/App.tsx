import { useState, useEffect } from 'react';
import { Cpu, Sparkles, FlaskConical, BarChart2, ServerOff } from 'lucide-react';
import axios from 'axios';
import SimulationDashboard from './pages/SimulationDashboard.tsx';
import PhysicsPlots from './pages/PhysicsPlots.tsx';
import BenchmarkDashboard from './pages/BenchmarkDashboard.tsx';

type Page = 'simulation' | 'plots' | 'benchmarks';

function Header({ page, setPage, isServerReady }: { page: Page; setPage: (p: Page) => void, isServerReady: boolean }) {
  return (
    <header className="border-b border-slate-800/60 bg-[#080914]/40 backdrop-blur-xl sticky top-0 z-50">
      <div className="max-w-7xl mx-auto px-6 py-4 flex justify-between items-center">
        <div className="flex items-center gap-3">
          <div className="w-10 h-10 rounded-xl bg-gradient-to-tr from-indigo-500 via-purple-500 to-cyan-400 flex items-center justify-center shadow-lg shadow-indigo-500/25">
            <Cpu className="w-5 h-5 text-white" />
          </div>
          <div>
            <h1 className="text-base font-bold tracking-tight bg-gradient-to-r from-white via-slate-100 to-indigo-200 bg-clip-text text-transparent">
              SSE H2SQ Simulation Hub
            </h1>
            <p className="text-[9px] text-indigo-400/80 uppercase tracking-widest font-bold flex items-center gap-1">
              <Sparkles className="w-2.5 h-2.5" /> High Performance Computing
            </p>
          </div>
        </div>

        <div className="flex items-center gap-4">
          {/* Page navigation */}
          <nav className="bg-slate-950/80 p-0.5 border border-slate-800/80 rounded-lg flex">
            <button
              onClick={() => setPage('simulation')}
              className={`flex items-center gap-1.5 px-3 py-1.5 text-[10px] font-bold rounded-md transition-all duration-200 ${page === 'simulation'
                  ? 'bg-indigo-500/20 text-indigo-300 border border-indigo-500/30'
                  : 'text-slate-500 hover:text-slate-300 border border-transparent'
                }`}
            >
              <FlaskConical className="w-3 h-3" />
              Simulation
            </button>
            <button
              onClick={() => setPage('plots')}
              className={`flex items-center gap-1.5 px-3 py-1.5 text-[10px] font-bold rounded-md transition-all duration-200 ${page === 'plots'
                  ? 'bg-purple-500/20 text-purple-300 border border-purple-500/30'
                  : 'text-slate-500 hover:text-slate-300 border border-transparent'
                }`}
            >
              <BarChart2 className="w-3 h-3" />
              Publication Plots
            </button>
            <button onClick={() => setPage('benchmarks')} className={`flex items-center gap-1.5 px-3 py-1.5 text-[10px] font-bold rounded-md border ${page === 'benchmarks' ? 'bg-cyan-500/20 text-cyan-300 border-cyan-500/30' : 'text-slate-500 border-transparent hover:text-slate-300'}`}><BarChart2 className="w-3 h-3" />ED / QMC Benchmarks</button>
          </nav>

          {isServerReady ? (
            <div className="px-3 py-1 rounded-full bg-emerald-500/10 border border-emerald-500/20 text-emerald-400 text-xs font-semibold flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-emerald-400 animate-ping"></span>
              Server Ready
            </div>
          ) : (
            <div className="px-3 py-1 rounded-full bg-red-500/10 border border-red-500/20 text-red-400 text-xs font-semibold flex items-center gap-2">
              <ServerOff className="w-3 h-3" />
              Server Offline
            </div>
          )}
        </div>
      </div>
    </header>
  );
}

function App() {
  const [page, setPage] = useState<Page>(window.location.hash === '#benchmarks' ? 'benchmarks' : window.location.hash === '#plots' ? 'plots' : 'simulation');
  const [isServerReady, setIsServerReady] = useState(false);

  useEffect(() => { window.history.replaceState(null, '', `#${page}`); }, [page]);

  useEffect(() => {
    const checkHealth = async () => {
      try {
        await axios.get('http://127.0.0.1:8000/', { timeout: 2000 });
        setIsServerReady(true);
      } catch (error) {
        setIsServerReady(false);
      }
    };

    // Initial check
    checkHealth();

    // Poll every 5 seconds
    const interval = setInterval(checkHealth, 5000);
    return () => clearInterval(interval);
  }, []);

  return (
    <div className="min-h-screen text-slate-100 selection:bg-indigo-500/30 relative">
      <div className="glow-bg"></div>
      <div className="glow-bg-2"></div>

      <Header page={page} setPage={setPage} isServerReady={isServerReady} />

      {page === 'simulation' ? (
        <main className="max-w-7xl mx-auto px-6 py-8 relative z-10">
          <SimulationDashboard isServerReady={isServerReady} />
        </main>
      ) : page === 'benchmarks' ? (
        <BenchmarkDashboard />
      ) : (
        <PhysicsPlots />
      )}
    </div>
  );
}

export default App;

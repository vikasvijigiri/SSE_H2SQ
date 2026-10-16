"""Copyright (c) 2026 Vikas Vijigiri. All rights reserved.

Local services: python scripts/deploy_all.py [up|status|stop|restart]
The default 'up' prepares dependencies, builds QMC and starts both services.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
import uuid

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / ".cache" / "services"
STATE = RUNTIME / "state.json"
PYTHON = ROOT / ".venv" / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
URLS = {"backend": "http://127.0.0.1:8000/", "frontend": "http://127.0.0.1:5173/"}


def source_version():
    digest = hashlib.sha256()
    paths = [Path(__file__).resolve(), ROOT / "backend/build/Makefile", ROOT / "backend/requirements.txt", ROOT / "frontend/package-lock.json"]
    for folder in ("backend/src", "backend/include", "backend/benchmarks", "backend/app", "frontend/src"):
        paths.extend(p for p in (ROOT / folder).rglob("*") if p.is_file() and p.suffix in {".py", ".cpp", ".hpp", ".ts", ".tsx", ".css"})
    paths.extend(p for p in (ROOT / "frontend").iterdir() if p.suffix in {".json", ".js", ".ts", ".html"})
    for path in sorted(set(paths)):
        digest.update(path.relative_to(ROOT).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def read_state():
    try:
        return json.loads(STATE.read_text())
    except (OSError, ValueError):
        return {}


def healthy(name):
    try:
        with urllib.request.urlopen(URLS[name], timeout=2) as response:
            if name == "backend":
                return json.load(response).get("message") == "Welcome to the SSE_H2SQ Simulation API"
            return response.status == 200 and b"/src/main.tsx" in response.read()
    except (OSError, ValueError):
        return False


def status():
    state = read_state()
    results = {name: healthy(name) for name in URLS}
    for name, ready in results.items():
        print(f"{name}: {'ready' if ready else 'offline'} {URLS[name]}")
    print(f"Logs: {RUNTIME}")
    return bool(state) and all(results.values())


def run(command, cwd=ROOT):
    print("Running:", subprocess.list2cmdline([str(x) for x in command]), flush=True)
    subprocess.run([str(x) for x in command], cwd=cwd, check=True)


def prepare():
    node = shutil.which("node")
    npm = shutil.which("npm.cmd" if os.name == "nt" else "npm")
    if not node or not npm or not shutil.which("make") or not shutil.which("g++"):
        raise RuntimeError("Install Node.js/npm, GNU Make and g++ before deploying.")
    if not PYTHON.exists():
        run([sys.executable, "-m", "venv", ROOT / ".venv"])
    requirements = ROOT / "backend/requirements.txt"
    stamp = RUNTIME / "requirements.installed"
    if not stamp.exists() or stamp.read_bytes() != requirements.read_bytes():
        run([PYTHON, "-m", "pip", "install", "-r", requirements])
        stamp.write_bytes(requirements.read_bytes())
    lock = ROOT / "frontend/package-lock.json"
    stamp = RUNTIME / "frontend.installed"
    vite = ROOT / "frontend/node_modules/vite/bin/vite.js"
    if not vite.exists() or not stamp.exists() or stamp.read_bytes() != lock.read_bytes():
        if os.name == "nt":
            npm_cli = Path(npm).parent / "node_modules/npm/bin/npm-cli.js"
            if not npm_cli.exists():
                raise RuntimeError(f"Cannot locate npm CLI: {npm_cli}")
            run([node, npm_cli, "ci"], ROOT / "frontend")
        else:
            run([npm, "ci"], ROOT / "frontend")
        stamp.write_bytes(lock.read_bytes())
    run(["make", "-C", ROOT / "backend/build", "all", "benchmark"])
    run(["make", "-C", ROOT / "backend/build", "debugbin"])


def worker(token):
    """Own child handles; stop files never kill unrelated processes by PID."""
    processes = []
    logs = []
    stop_file = RUNTIME / f"{token}.stop"
    heartbeat = RUNTIME / f"{token}.heartbeat"
    commands = {
        "backend": [str(PYTHON), "-m", "uvicorn", "app.main:app", "--app-dir", str(ROOT / "backend"), "--host", "127.0.0.1", "--port", "8000"],
        "frontend": [shutil.which("node"), str(ROOT / "frontend/node_modules/vite/bin/vite.js"), "--host", "127.0.0.1", "--port", "5173", "--strictPort"],
    }
    scientific_environment=os.environ.copy()
    # QMC chains already run in separate processes; avoid nested BLAS pools.
    for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
        scientific_environment.setdefault(variable,"1")
    try:
        for name, command in commands.items():
            log = (RUNTIME / f"{name}.log").open("ab")
            logs.append(log)
            processes.append(subprocess.Popen(command, cwd=ROOT / ("frontend" if name == "frontend" else "backend"), env=scientific_environment, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL))
        state = {"token": token, "supervisor": os.getpid(), "children": [p.pid for p in processes], "source_version": source_version()}
        heartbeat.touch()
        temporary = STATE.with_suffix(".tmp")
        temporary.write_text(json.dumps(state))
        temporary.replace(STATE)
        restarts = {name: 0 for name in commands}
        while not stop_file.exists():
            for index, (name, command) in enumerate(commands.items()):
                code = processes[index].poll()
                if code is not None:
                    restarts[name] += 1
                    print(f"{name} exited with code {code}; restarting (attempt {restarts[name]}).", flush=True)
                    time.sleep(min(10, .5 * 2**min(restarts[name]-1, 5)))
                    if stop_file.exists():
                        break
                    processes[index] = subprocess.Popen(command, cwd=ROOT / ("frontend" if name == "frontend" else "backend"), env=scientific_environment, stdout=logs[index], stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
                    state["children"] = [process.pid for process in processes]
                    for attempt in range(20):
                        try:
                            temporary.write_text(json.dumps(state))
                            temporary.replace(STATE)
                            break
                        except PermissionError:
                            if attempt == 19:
                                raise
                            time.sleep(.05)
            try:
                heartbeat.touch()
            except PermissionError:
                pass  # A transient Windows scanner lock must not stop services.
            time.sleep(0.5)
    finally:
        for process in processes:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
        for log in logs:
            log.close()
        if read_state().get("token") == token:
            STATE.unlink(missing_ok=True)
        stop_file.unlink(missing_ok=True)
        heartbeat.unlink(missing_ok=True)


def stop():
    state = read_state()
    if not state:
        print("No managed services are running.")
        return
    heartbeat = RUNTIME / f"{state['token']}.heartbeat"
    if not heartbeat.exists() or time.time() - heartbeat.stat().st_mtime > 30:
        STATE.unlink(missing_ok=True)
        heartbeat.unlink(missing_ok=True)
        print("Cleared expired supervisor state; no processes were killed.")
        return
    (RUNTIME / f"{state['token']}.stop").touch()
    deadline = time.monotonic() + 30
    while read_state().get("token") == state["token"] and time.monotonic() < deadline:
        time.sleep(0.5)
    if read_state().get("token") == state["token"]:
        raise RuntimeError("Supervisor did not respond. Inspect service logs; no unrelated processes were killed.")
    print("Stopped frontend and backend.")


def up():
    if read_state():
        if status() and read_state().get("source_version") == source_version():
            print("Both services are already running.")
            return
        print("Restarting services to deploy current source.")
        stop()
    for port in (8000, 5173):
        with socket.socket() as connection:
            if connection.connect_ex(("127.0.0.1", port)) == 0:
                raise RuntimeError(f"Port {port} is occupied by an unmanaged service. Stop it before deploying.")
    prepare()
    token = uuid.uuid4().hex
    options = {"creationflags": subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {"start_new_session": True}
    with (RUNTIME / "supervisor.log").open("ab") as log:
        process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "_worker", "--token", token], cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, **options)
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(f"Services exited during startup. Inspect logs in {RUNTIME}")
        if read_state().get("token") == token and all(healthy(name) for name in URLS):
            status()
            return
        time.sleep(1)
    stop()
    raise RuntimeError(f"Startup timed out. Inspect logs in {RUNTIME}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", nargs="?", default="up", choices=["up", "status", "stop", "restart", "_worker"])
    parser.add_argument("--token", help=argparse.SUPPRESS)
    args = parser.parse_args()
    RUNTIME.mkdir(parents=True, exist_ok=True)
    if args.action == "_worker":
        worker(args.token)
        return
    if args.action == "status":
        raise SystemExit(0 if status() else 1)
    lock = RUNTIME / "operation.lock"
    try:
        lock.mkdir()
    except FileExistsError:
        raise RuntimeError(f"Another deployment operation is in progress: {lock}")
    try:
        if args.action in {"stop", "restart"}:
            stop()
        if args.action in {"up", "restart"}:
            up()
    finally:
        lock.rmdir()


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Deployment failed: {error}", file=sys.stderr)
        raise SystemExit(1)

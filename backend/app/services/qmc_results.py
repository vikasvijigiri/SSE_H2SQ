"""Persistent QMC runs and delete-one-run jackknife estimates."""
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import uuid
from contextlib import contextmanager

DATABASE = Path(__file__).resolve().parents[2] / "files" / "qmc_runs.sqlite3"
OBSERVABLES = ("enrg", "enrg2", "enrg4", "SMag_square", "SMag_four", "Mag_x")


@contextmanager
def connect():
    DATABASE.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(DATABASE, timeout=30)
    try:
        connection.execute("CREATE TABLE IF NOT EXISTS runs (id TEXT PRIMARY KEY, parameter_key TEXT NOT NULL, record TEXT NOT NULL, raw_output TEXT)")
        connection.execute("CREATE INDEX IF NOT EXISTS runs_parameters ON runs(parameter_key)")
        with connection:
            yield connection
    finally:
        connection.close()


def save_run(key, record, raw_output=None, legacy=False):
    result = dict(record)
    identity = result.get("run_id")
    if not identity:
        identity = ("legacy-" + hashlib.sha256((key + json.dumps(record, sort_keys=True)).encode()).hexdigest()) if legacy else uuid.uuid4().hex
    result["run_id"] = identity
    if raw_output is not None:
        result["result_file_url"] = f"/api/qmc-runs/{identity}/result"
    with connect() as connection:
        connection.execute("INSERT OR IGNORE INTO runs VALUES (?, ?, ?, ?)", (identity, key, json.dumps(result), raw_output))
    return result


def import_legacy(key, record):
    if record:
        save_run(key, record, legacy=True)


def jackknife_summary(key):
    with connect() as connection:
        records = [json.loads(row[0]) for row in connection.execute("SELECT record FROM runs WHERE parameter_key = ?", (key,))]
    values = {}
    for observable in OBSERVABLES:
        samples = [float(record[observable]) for record in records if isinstance(record.get(observable), (int, float)) and math.isfinite(record[observable])]
        n = len(samples)
        if not n:
            continue
        mean = math.fsum(samples) / n
        # For a linear mean, the bias-corrected jackknife estimate equals the mean.
        # Delete-one-run SE simplifies to sqrt(sum((x-mean)^2)/(n*(n-1))).
        error = math.sqrt(math.fsum((x - mean) ** 2 for x in samples) / (n * (n - 1))) if n > 1 else None
        values[observable] = {"mean": mean, "standard_error": error, "n_runs": n}
    return {"n_runs": len(records), "observables": values, "method": "equal-weight delete-one-run jackknife; between-run standard error"}


def result_text(run_id):
    with connect() as connection:
        row = connection.execute("SELECT raw_output FROM runs WHERE id = ?", (run_id,)).fetchone()
    return row[0] if row else None

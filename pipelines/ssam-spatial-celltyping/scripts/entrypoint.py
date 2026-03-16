#!/opt/conda/bin/python3
"""
Entry point: runs SSAM pipeline without invoking any system shell.
Reads config.yaml and calls run_ssam.py directly via subprocess with
the conda Python interpreter — bypasses /bin/sh and /usr/bin/bash entirely.
"""
import sys
import os
import subprocess
import yaml

CONFIG_PATH = "/pipeline/config.yaml"
SCRIPT_PATH = "/pipeline/scripts/run_ssam.py"
PYTHON      = "/opt/conda/bin/python3"

def main():
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)

    input_mode  = cfg.get("input_mode",  "zarr")
    zarr        = cfg.get("zarr",        "")
    coordinates = cfg.get("coordinates", "")
    signatures  = cfg.get("signatures",  "")
    bandwidth   = str(cfg.get("bandwidth",  2.5))
    threshold   = str(cfg.get("threshold",  0.2))
    map_width   = str(cfg.get("map_width",  2000))
    threads     = str(cfg.get("threads",    4))
    output_dir  = "/output"

    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        PYTHON, SCRIPT_PATH,
        "--input-mode", input_mode,
        "--output-dir", output_dir,
        "--bandwidth",  bandwidth,
        "--threshold",  threshold,
        "--map-width",  map_width,
        "--threads",    threads,
    ]
    if input_mode == "zarr":
        cmd += ["--zarr", zarr]
    else:
        cmd += ["--coordinates", coordinates, "--signatures", signatures]

    log_path = os.path.join(output_dir, "ssam.log")
    print(f"[entrypoint] Running: {' '.join(cmd)}", flush=True)

    with open(log_path, "w") as logf:
        result = subprocess.run(cmd, stdout=logf, stderr=logf)

    with open(log_path) as logf:
        print(logf.read(), flush=True)

    if result.returncode != 0:
        sys.exit(f"[entrypoint] FAILED (exit {result.returncode}). See {log_path}")

    print("[entrypoint] Done.", flush=True)

if __name__ == "__main__":
    main()


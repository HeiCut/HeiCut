#!/usr/bin/env python3
import csv
import os
from pathlib import Path

# Relative paths of *raw* CSVs you want to convert (from this file’s location)
INPUT_CSVS = [
    "../generated/all_results/ilp/all_results.csv",
    "../generated/all_results/kernelizer_IT0/all_results.csv",
    "../generated/all_results/kernelizer_IT1/all_results.csv",
    "../generated/all_results/submodular_tight_single/all_results.csv"
]
# Where the per-algorithm CSVs should land (R script expects these)
# Make this point to the folder that `instances.R` reads from
OUTPUT_DIR = Path("../generated/all_results")
# ------------------------------------------------

def convert_one(input_csv_path: Path, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)

    with input_csv_path.open(mode='r', newline='') as infile:
        reader = csv.reader(infile, delimiter=';')
        rows = list(reader)

    algorithms = {}
    total_failed_count = 0

    for row in rows:
        if len(row) != 7:
            print(f"[{input_csv_path}] Skipping malformed row: {row}")
            continue

        mincut, time, peak_memory, is_exact, hypergraph, seed, algorithm = [col.strip() for col in row]

        try:
            # treat empty mincut as "failed" 
            is_exact_flag = 0 if mincut.strip() else 1
            # is_exact_flag = 0 if int(is_exact) == 1 else 1  # original alt
        except ValueError:
            print(f"[{input_csv_path}] Skipping row due to invalid is_exact value: {row}")
            continue

        try:
            if mincut:
                mincut = str(int(mincut) + 1)
        except ValueError:
            print(f"[{input_csv_path}] Skipping row due to invalid mincut value: {row}")
            continue

        if is_exact_flag == 1:
            total_failed_count += 1
            time, peak_memory, mincut = '0', '0', '0'

        algorithms.setdefault(algorithm, []).append({
            'algorithm_name': algorithm,
            'graph_name': hypergraph,
            'runtime': time,
            'memory': peak_memory,
            'mincut': mincut,
            'k_value': seed,
            'Failed': is_exact_flag
        })

    print(f"[{input_csv_path}] Total failed = {total_failed_count}")

    # Write one CSV per algorithm into OUTPUT_DIR
    for algorithm, entries in algorithms.items():
        # NOTE: instances.R uses base names like "SMHM_ilp", "kernelizer_IT0", etc.
        # If your input 'algorithm' values already match those, this will line up.
        out_path = out_dir / f"{algorithm}.csv"
        with out_path.open(mode='w', newline='') as outfile:
            writer = csv.DictWriter(
                outfile,
                fieldnames=['algorithm_name', 'graph_name', 'runtime', 'memory', 'mincut', 'k_value', 'Failed']
            )
            writer.writeheader()
            writer.writerows(entries)
        print(f"  -> wrote {out_path}")

def main():
    here = Path(__file__).resolve().parent
    for rel in INPUT_CSVS:
        input_csv = (here / rel).resolve()
        if not input_csv.exists():
            print(f"[WARN] Missing input: {rel}")
            continue
        convert_one(input_csv, OUTPUT_DIR)

if __name__ == "__main__":
    main()

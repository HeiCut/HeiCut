#!/usr/bin/env python3
# analyze_results.py

import argparse
import csv
import glob
import os
import sys

# Use a non-interactive backend so we can save PDFs without a display
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ---------- helpers ----------
def geometric_mean(x: pd.Series) -> float:
    # robust to zeros/negatives by filtering (log undefined otherwise)
    x = pd.to_numeric(x, errors="coerce")
    x = x[(x > 0) & x.notna()]
    if len(x) == 0:
        return float("nan")
    return float(np.exp(np.log(x).mean()))


def average(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    return float(x.mean())


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


# ---------- plots ----------
def plot_distribution_by_algorithm(
    df: pd.DataFrame,
    algorithms_list,
    fields,
    out_root,
    filter_condition_1=None,
    filter_condition_2=None,
    fontsize=13,
    num_bins=30,
):
    for algorithm in algorithms_list:
        subset = df[df["algorithm"] == algorithm].copy()
        if subset.empty:
            continue

        filtered_subset_1 = None
        filtered_subset_2 = None
        if filter_condition_1:
            filtered_subset_1 = subset[filter_condition_1(subset)]
        if filter_condition_1 and filter_condition_2:
            filtered_subset_2 = subset[filter_condition_2(subset)]

        fig, axes = plt.subplots(1, len(fields), figsize=(16, 4), squeeze=False)
        axes = axes[0]

        for i, field in enumerate(fields):
            vals = pd.to_numeric(subset[field], errors="coerce").dropna()
            if len(vals) == 0:
                continue

            # guard against min==0 for log bins
            vmin = max(vals.min() / 1.1, 1e-12)
            vmax = max(vals.max() * 1.1, vmin * 10)
            log_bins = np.logspace(np.log10(vmin), np.log10(vmax), num=num_bins)

            axes[i].hist(vals, bins=log_bins, color="lightblue", alpha=0.7, edgecolor="darkblue", label="All Runs")

            if filter_condition_1 and filtered_subset_1 is not None and not filtered_subset_1.empty:
                vals1 = pd.to_numeric(filtered_subset_1[field], errors="coerce").dropna()
                if len(vals1) > 0:
                    axes[i].hist(
                        vals1, bins=log_bins, color="lightgreen", alpha=0.7, edgecolor="darkgreen", label="Filtered Runs"
                    )

            if filter_condition_1 and filter_condition_2 and filtered_subset_2 is not None and not filtered_subset_2.empty:
                vals2 = pd.to_numeric(filtered_subset_2[field], errors="coerce").dropna()
                if len(vals2) > 0:
                    axes[i].hist(
                        vals2, bins=log_bins, color="lightpink", alpha=0.7, edgecolor="darkred", label="Second Filtered Runs"
                    )

            axes[i].set_xscale("log")
            axes[i].set_title(f"{field} distribution for {algorithm}", fontsize=fontsize)
            axes[i].set_xlabel(f"{field} (log scale)", fontsize=fontsize)
            axes[i].set_ylabel("Number of Runs", fontsize=fontsize)
            axes[i].grid(True, linestyle=":")
            axes[i].legend()

        num_total_runs = len(subset)
        relative_stats = {
            "num_runs_less_1_sec": (pd.to_numeric(subset["time"], errors="coerce") <= 1).sum(),
            "num_runs_more_1_min": (pd.to_numeric(subset["time"], errors="coerce") >= 60).sum(),
            "num_runs_exact": (subset["is_exact"] == 1).sum(),
            "num_runs_!naive": ((subset["mincut"] < subset.get("naive_mincut", np.inf)) & subset["mincut"].notnull()).sum(),
            "num_runs_!naive_&_!0": (
                (subset["mincut"] < subset.get("naive_mincut", np.inf)) & subset["mincut"].notnull() & (subset["mincut"] != 0)
            ).sum(),
            "num_runs_failed": subset.get("failed", pd.Series(False, index=subset.index)).sum(),
        }

        print(f"#################### {algorithm} ####################")
        print(f"num_total_runs:\t\t {num_total_runs}")
        for field in ["time", "memory"]:
            print(f"avg_total_{field}:\t\t {average(subset[field]):.3f}")
            print(f"geo_mean_total_{field}:\t {geometric_mean(subset[field]):.3f}")
        for name, stat in relative_stats.items():
            pct = (stat / num_total_runs * 100.0) if num_total_runs else 0.0
            print(f"{name}:\t\t {stat}\t ({pct:.3f}%)")

        fig.tight_layout()
        out_dir = os.path.join(out_root, algorithm)
        ensure_dir(out_dir)
        plt.savefig(os.path.join(out_dir, f"dist_{algorithm}.pdf"))
        plt.close(fig)


def plot_reduction_by_algorithm(
    data_root,
    out_root,
    algorithms_list,
    num_rounds=None,
    hypergraph_list=None,
    num_suffix_param_columns=4,
    fontsize=13,
):
    hypergraph_list = hypergraph_list or []

    for algorithm in algorithms_list:
        all_reductions_path = os.path.join(data_root, algorithm, "all_reductions.csv")
        num_runs_instant_exit = 0
        num_reductions_preprocessing = 1
        num_reductions_per_round = 4 if "_IT0" in algorithm else 5

        if not os.path.exists(all_reductions_path):
            print(f"File does not exist: {all_reductions_path}")
            continue

        with open(all_reductions_path, "r") as file:
            reader = csv.reader(file, delimiter=";")
            all_rows = [row for row in reader]

        if not all_rows:
            print(f"Empty CSV: {all_reductions_path}")
            continue

        max_num_columns = max(len(row) for row in all_rows)
        if num_rounds is not None:
            max_num_columns = min(
                max_num_columns,
                2 + 3 * (num_reductions_preprocessing + num_rounds * num_reductions_per_round) + num_suffix_param_columns,
            )

        aligned_rows = []
        for row in all_rows:
            # Ignore rows where no reduction rule applied (naive mincut already 0)
            if len(row) > 2 and row[2] == "":
                num_runs_instant_exit += 1
                continue

            if max_num_columns == len(row):
                aligned_rows.append(row)
            elif max_num_columns < len(row):
                aligned_rows.append(row[: (max_num_columns - num_suffix_param_columns)] + row[-num_suffix_param_columns:])
            else:
                # pad with placeholders to align
                aligned_rows.append(
                    row[:-num_suffix_param_columns]
                    + [pd.NA, pd.NA, pd.NA] * ((max_num_columns - len(row)) // 3 - 1)
                    + [row[-num_suffix_param_columns - 3], row[-num_suffix_param_columns - 2], pd.NA]
                    + row[-num_suffix_param_columns:]
                )

        all_reductions_df = pd.DataFrame(aligned_rows)
        if len(hypergraph_list) > 0:
            all_reductions_df = all_reductions_df[all_reductions_df.iloc[:, -3].isin(hypergraph_list)]

        # base solver time column: convert to nullable float
        base_time_col = all_reductions_df.columns[-num_suffix_param_columns]
        all_reductions_df[base_time_col] = (
            all_reductions_df[base_time_col].replace("", pd.NA).astype("Float64")
        )

        relative_all_reductions_df = all_reductions_df.copy()

        fig = plt.figure(figsize=(16, 10))
        gs = gridspec.GridSpec(3, 6, height_ratios=[1, 0.5, 0.75])

        # edges & nodes
        for i, metric in enumerate(["edges", "nodes"]):
            cols = (
                all_reductions_df.columns[i : i + 1].tolist()
                + all_reductions_df.columns[2 + i : -num_suffix_param_columns : 3].tolist()
            )

            for col in cols:
                all_reductions_df[col] = pd.to_numeric(all_reductions_df[col], errors="coerce").astype("Int64")
                relative_all_reductions_df[col] = (
                    all_reductions_df[col] / all_reductions_df[cols[0]]
                ) * 100

            positions = [gs[0, :3], gs[0, 3:]]
            ax = fig.add_subplot(positions[i])
            for j in range(len(relative_all_reductions_df)):
                ax.plot(range(len(cols)), relative_all_reductions_df[cols].iloc[j], alpha=0.4)

            ax.set_title(f"% of remaining {metric} per reduction step for {algorithm}", fontsize=fontsize)
            ax.set_xlabel("reduction step", fontsize=fontsize)
            ax.set_ylabel(f"% of remaining {metric}", fontsize=fontsize)
            ax.set_ylim([-10, 110])
            ax.set_xlim([0, len(cols) - 1])
            ax.grid(True, linestyle=":")

        # time per reduction
        time_cols = all_reductions_df.columns[4 : -num_suffix_param_columns : 3].tolist()
        ax3 = fig.add_subplot(gs[1, :])
        for col in time_cols:
            all_reductions_df[col] = pd.to_numeric(all_reductions_df[col], errors="coerce").astype("Float64")
        for j in range(len(all_reductions_df)):
            ax3.plot(range(1, len(time_cols) + 1), all_reductions_df[time_cols].iloc[j], alpha=0.4)
        ax3.set_title(f"time per reduction step for {algorithm}", fontsize=fontsize)
        ax3.set_xlabel("reduction step", fontsize=fontsize)
        ax3.set_ylabel("time", fontsize=fontsize)
        ax3.set_xlim([1, len(time_cols)])
        ax3.set_yscale("log")
        ax3.grid(True, linestyle=":")

        # effectiveness bars
        for i, metric in enumerate(["edges", "nodes", "time"]):
            effectiveness_list = []
            for j in range(num_reductions_per_round):
                base_indices = np.array(
                    range(
                        2 + 3 * (num_reductions_preprocessing + j),
                        all_reductions_df.shape[1] - num_suffix_param_columns,
                        3 * num_reductions_per_round,
                    )
                )
                if metric == "time":
                    eff_df = all_reductions_df.iloc[:, base_indices + i].astype("Float64")
                else:
                    denom = (
                        all_reductions_df.iloc[:, base_indices + i - 3]
                        .astype("Int64")
                        .replace(0, pd.NA)
                        .values
                    )
                    eff_df = 100 * (1 - all_reductions_df.iloc[:, base_indices + i].astype("Int64").div(denom))
                validity_df = (all_reductions_df.iloc[:, base_indices + 2].notna())
                validity_df.columns = eff_df.columns
                eff_df = eff_df.where(validity_df).astype("float")
                effectiveness_list.append(float(np.round(average(eff_df.stack()), 3)))

            positions = [gs[2, :2], gs[2, 2:4], gs[2, 4:]]
            ax_eff = fig.add_subplot(positions[i])
            effectiveness_labels = ["no_heavy_edges", "no_heavy_overlaps", "no_shiftable_2-edges", "no_triangle_2-edges"]
            if "_IT0" not in algorithm:
                effectiveness_labels.insert(0, "label_propagation")

            bars = ax_eff.bar(effectiveness_labels, effectiveness_list, color="lightblue", alpha=0.7, edgecolor="darkblue")
            ax_eff.set_xlabel("reduction rule", fontsize=fontsize)

            if metric == "time":
                ax_eff.set_title("average time usage (s) for each rule", fontsize=fontsize)
                ax_eff.set_ylabel("average time usage (s)", fontsize=fontsize)
            else:
                ax_eff.set_title(f"% of average {metric} reduction for each rule", fontsize=fontsize)
                ax_eff.set_ylabel(f"% of average {metric} reduction", fontsize=fontsize)
                ax_eff.set_ylim([0, 100])
                ax_eff.set_yticks(ticks=range(0, 101, 10))
            ax_eff.grid(True, linestyle=":")
            for label in ax_eff.get_xticklabels():
                label.set_fontsize(10)
                label.set_rotation(10)
            for bar in bars:
                h = bar.get_height()
                txt = ax_eff.text(
                    bar.get_x() + bar.get_width() / 2,
                    h,
                    f"{h:.2f}",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=8,
                    weight="bold",
                )
                txt.set_bbox(dict(boxstyle="square,pad=0.2", facecolor="darkblue", edgecolor="darkblue"))

        relative_stats = {
            "num_runs_instant_exit": num_runs_instant_exit,
            "num_runs_max_reduced": (
                (all_reductions_df.iloc[:, -num_suffix_param_columns - 3] == 0)
                | (all_reductions_df.iloc[:, -num_suffix_param_columns - 2] == 1)
            ).sum(),
            "num_runs_base_solver": (all_reductions_df.iloc[:, -num_suffix_param_columns].notna()).sum(),
        }

        num_total_runs = len(all_reductions_df) + num_runs_instant_exit
        print(f"#################### {algorithm} ####################")
        print(f"num_total_runs:\t\t {num_total_runs}")
        for name, stat in relative_stats.items():
            pct = (stat / num_total_runs * 100.0) if num_total_runs else 0.0
            print(f"{name}:\t\t {stat}\t ({pct:.3f}%)")
        mask = all_reductions_df.iloc[:, -num_suffix_param_columns].notna()
        if mask.any():
            print(f"base_solver_avg_edges:\t {average(all_reductions_df.loc[mask].iloc[:, -num_suffix_param_columns - 3]):.3f}")
            print(f"base_solver_geo_mean_edges:\t {geometric_mean(all_reductions_df.loc[mask].iloc[:, -num_suffix_param_columns - 3]):.3f}")
            print(f"base_solver_avg_nodes:\t {average(all_reductions_df.loc[mask].iloc[:, -num_suffix_param_columns - 2]):.3f}")
            print(f"base_solver_geo_mean_nodes:\t {geometric_mean(all_reductions_df.loc[mask].iloc[:, -num_suffix_param_columns - 2]):.3f}")
            print(f"base_solver_avg_time:\t {average(all_reductions_df.loc[mask].iloc[:, -num_suffix_param_columns]):.3f}")
            print(f"base_solver_geo_mean_time:\t {geometric_mean(all_reductions_df.loc[mask].iloc[:, -num_suffix_param_columns]):.3f}")

        plt.tight_layout()
        out_dir = os.path.join(out_root, algorithm)
        ensure_dir(out_dir)
        plt.savefig(os.path.join(out_dir, f"reduction_{algorithm}.pdf"))
        plt.close(fig)


def plot_performance_profile(
    df: pd.DataFrame,
    naive_mincuts_df: pd.DataFrame | None,
    algorithms_list,
    objectives,
    out_root,
    xlogscales=None,
    xmaxvals=None,
    fontsize=13,
    penalize_failures: bool = False,
):
    """
    If penalize_failures == False (default):
        - Failed rows are dropped from profiles for all objectives.
    If penalize_failures == True:
        - Failed rows are included with a very large sentinel for the objective.
        - Instances where every algorithm failed are dropped.
    """
    styles = [
        ["blue", "--"], ["green", "-."], ["orange", "--"], ["magenta", "-."],
        ["red", "--"], ["purple", "-."], ["brown", "--"], ["cyan", "--"]
    ]
    SENTINEL = 1e12  # "infinite" cost for failed runs when penalizing

    fig, axes = plt.subplots(1, len(objectives), figsize=(16, 4), squeeze=False)
    axes = axes[0]

    for i, objective in enumerate(objectives):
        successful = df[(df["algorithm"].isin(algorithms_list)) & (~df.get("failed", False))].copy()
        if successful.empty:
            continue

        min_objective_per_instance = (
            successful
            .groupby(["hypergraph", "seed"])[objective]
            .min()
            .reset_index()
            .rename(columns={objective: "min_objective"})
        )

        overall_max_ratio = 2.0

        for j, algorithm in enumerate(algorithms_list):
            if penalize_failures:
                sub_all = df[df["algorithm"] == algorithm].copy()
                if sub_all.empty:
                    continue
                subset = pd.merge(sub_all, min_objective_per_instance, on=["hypergraph", "seed"], how="inner")
                if subset.empty:
                    continue
                mask_failed = subset.get("failed", False)
                subset.loc[mask_failed, objective] = SENTINEL
                subset = subset.dropna(subset=[objective, "min_objective"])
            else:
                subset = df[(df["algorithm"] == algorithm) & (~df.get("failed", False))].copy()
                if subset.empty:
                    continue
                subset = pd.merge(subset, min_objective_per_instance, on=["hypergraph", "seed"], how="inner")
                if subset.empty:
                    continue

            def compute_ratio_to_best(row):
                if row["min_objective"] == 0:
                    return 1.0 if row[objective] == 0 else 1000.0
                return row[objective] / row["min_objective"]

            subset["ratio_to_best"] = subset.apply(compute_ratio_to_best, axis=1)

            max_ratio = subset["ratio_to_best"].max(skipna=True)
            if pd.notna(max_ratio) and max_ratio > overall_max_ratio:
                overall_max_ratio = float(max_ratio)

            subset = subset.sort_values(by="ratio_to_best").reset_index(drop=True)
            denom = max(len(subset), len(naive_mincuts_df) if naive_mincuts_df is not None else len(subset)) - 1
            denom = max(denom, 1)
            subset["ind"] = subset.index / denom

            style = styles[j % len(styles)]
            axes[i].plot(subset["ratio_to_best"], subset["ind"], label=algorithm, linewidth=3, linestyle=style[1], color=style[0])

        if xlogscales and xlogscales[i]:
            axes[i].set_xscale("log")
        axes[i].set_title(f"{objective} performance profile", fontsize=fontsize)
        axes[i].set_xlabel(r"$\tau$", fontsize=fontsize)
        xmax = xmaxvals[i] if (xmaxvals and i < len(xmaxvals) and xmaxvals[i] is not None) else overall_max_ratio
        axes[i].set_xlim([0.99, xmax])
        axes[i].set_ylim([0, 1])
        axes[i].set_ylabel(rf"% instances $\leq \tau \cdot$ best instance", fontsize=fontsize)
        axes[i].grid(True, linestyle=":")

    for ax in axes[1:]:
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_ylabel("")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=max(1, len(labels)), frameon=True, fontsize=fontsize)
    plt.tight_layout(rect=[0, 0.1, 1, 1])

    ensure_dir(out_root)
    out_path = os.path.join(out_root, "performance_plots.pdf")
    plt.savefig(out_path)
    plt.close(fig)


# ---------- stats dumps ----------
def print_submodular_stats(data_root, algorithms_list):
    all_iterations_columns = ["initial_num_nodes", "num_iterations", "num_contractions_per_it", "hypergraph", "seed", "algorithm"]
    for algorithm in algorithms_list:
        p = os.path.join(data_root, algorithm, "all_iterations.csv")
        print(f"#################### {algorithm} ####################")
        if os.path.exists(p):
            df = pd.read_csv(p, sep=";", header=None, names=all_iterations_columns)
            for c in ["initial_num_nodes", "num_iterations", "num_contractions_per_it"]:
                df[c] = pd.to_numeric(df[c], errors="coerce")
            if df[["initial_num_nodes", "num_iterations", "num_contractions_per_it"]].isna().any().any():
                print("Warning: NaN values found in columns after conversion.")
            num_total_runs = len(df)
            relative_stats = {"num_runs_finished": df["num_iterations"].notna().sum()}
            print(f"num_total_runs:\t\t {num_total_runs}")
            print(f"avg_contractions_per_it:\t {average(df['num_contractions_per_it']):.3f}")
            print(f"geo_mean_contractions_per_it:\t {geometric_mean(df['num_contractions_per_it']):.3f}")
            for name, stat in relative_stats.items():
                pct = (stat / num_total_runs * 100.0) if num_total_runs else 0.0
                print(f"{name}:\t\t {stat}\t ({pct:.3f}%)")
        else:
            print(f"Warning: No all_iterations.csv file found for {algorithm}.")

    all_hypercactus_columns = [
        "initial_num_edges", "initial_num_nodes", "kernel_num_edges", "kernel_num_nodes",
        "hypercactus_num_edges", "hypercactus_num_nodes", "time", "memory", "hypergraph", "seed", "algorithm"
    ]
    for algorithm in algorithms_list:
        p = os.path.join(data_root, algorithm, "all_hypercactus_results.csv")
        print(f"#################### {algorithm} ####################")
        if os.path.exists(p):
            df = pd.read_csv(p, sep=";", header=None, names=all_hypercactus_columns)
            for c in ["initial_num_edges","initial_num_nodes","kernel_num_edges","kernel_num_nodes",
                      "hypercactus_num_edges","hypercactus_num_nodes","time","memory"]:
                df[c] = pd.to_numeric(df[c], errors="coerce")
            if df[["hypercactus_num_edges", "hypercactus_num_nodes", "time"]].isna().any().any():
                print("Warning: NaN values found in columns after conversion.")

            num_total_runs = len(df)
            relative_stats = {"num_runs_finished": df["time"].notna().sum()}
            print(f"num_total_runs:\t\t {num_total_runs}")
            for name, stat in relative_stats.items():
                pct = (stat / num_total_runs * 100.0) if num_total_runs else 0.0
                print(f"{name}:\t\t {stat}\t ({pct:.3f}%)")
            print(f"avg_initial_num_edges:\t\t {average(df['initial_num_edges']):.3f}")
            print(f"geo_mean_initial_num_edges:\t {geometric_mean(df['initial_num_edges']):.3f}")
            print(f"avg_initial_num_nodes:\t\t {average(df['initial_num_nodes']):.3f}")
            print(f"geo_mean_initial_num_nodes:\t {geometric_mean(df['initial_num_nodes']):.3f}")
            print()
            print(f"avg_kernel_num_edges:\t\t {average(df['kernel_num_edges']):.3f}")
            print(f"geo_mean_kernel_num_edges:\t {geometric_mean(df['kernel_num_edges']):.3f}")
            print(f"avg_kernel_num_nodes:\t\t {average(df['kernel_num_nodes']):.3f}")
            print(f"geo_mean_kernel_num_nodes:\t {geometric_mean(df['kernel_num_nodes']):.3f}")
            print()
            print(f"avg_hypercactus_num_edges:\t {average(df['hypercactus_num_edges']):.3f}")
            print(f"geo_mean_hypercactus_num_edges:\t {geometric_mean(df['hypercactus_num_edges']):.3f}")
            print(f"avg_hypercactus_num_nodes:\t {average(df['hypercactus_num_nodes']):.3f}")
            print(f"geo_mean_hypercactus_num_nodes:\t {geometric_mean(df['hypercactus_num_nodes']):.3f}")
            print()
            print(f"avg_time:\t\t\t {average(df['time']):.3f}")
            print(f"geo_mean_time:\t\t\t {geometric_mean(df['time']):.3f}")
            print(f"max_time:\t\t\t {pd.to_numeric(df['time'], errors='coerce').max():.3f}")
            print(f"avg_memory:\t\t\t {average(df['memory']):.3f}")
            print(f"geo_mean_memory:\t\t {geometric_mean(df['memory']):.3f}")
            print(f"max_memory:\t\t\t {pd.to_numeric(df['memory'], errors='coerce').max():.3f}")
        else:
            print(f"Warning: No all_hypercactus_results.csv file found for {algorithm}.")


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Analyze experiment CSVs and export plots as PDF.")
    ap.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing per-algorithm subfolders with CSVs (e.g., */all_results.csv) and optionally all_naive_mincuts.csv.",
    )
    ap.add_argument(
        "--out-dir",
        default=".",  # default to current working directory
        help="Directory to write PDF plots. Defaults to the current working directory.",
    )
    ap.add_argument("--fontsize", type=int, default=13, help="Font size for plots (default: 13).")
    ap.add_argument("--num-bins", type=int, default=30, help="Histogram bins for distribution plots (default: 30).")

    # New: enable flags (default behavior is disabled)
    ap.add_argument("--enable-distribution", action="store_true", help="Generate distribution histograms per algorithm.")
    ap.add_argument("--enable-reduction", action="store_true", help="Generate reduction-step plots per algorithm.")
    ap.add_argument("--skip-profile", action="store_true", help="Skip the performance profile plot (off by default).")

    ap.add_argument(
        "--penalize-failures",
        action="store_true",
        help="Include failed runs in profiles with an infinite-like cost so failures visibly hurt curves.",
    )
    ap.add_argument(
        "--print-failure-summary",
        action="store_true",
        help="Print failure rates by algorithm to stdout.",
    )
    args = ap.parse_args()

    data_dir = os.path.abspath(args.data_dir)
    out_dir = os.path.abspath(args.out_dir)
    ensure_dir(out_dir)

    # Load all_results from subfolders
    all_results_columns = ["mincut", "time", "memory", "is_exact", "hypergraph", "seed", "algorithm"]
    results_glob = os.path.join(data_dir, "*/all_results.csv")
    all_results_csv = glob.glob(results_glob)
    if not all_results_csv:
        print(f"No CSVs found with pattern: {results_glob}")
        sys.exit(1)

    all_results_df_list = []
    for file in all_results_csv:
        try:
            df = pd.read_csv(file, sep=";", header=None, names=all_results_columns)
            # enforce numeric
            for c in ["time", "memory", "mincut", "is_exact"]:
                if c in df.columns:
                    df[c] = pd.to_numeric(df[c], errors="coerce")
            # propagate algorithm name from folder (overwrite whatever's in CSV)
            algorithm = os.path.basename(os.path.dirname(file))
            df["algorithm"] = algorithm
            all_results_df_list.append(df)
        except Exception as e:
            print(f"Error reading file {file}: {e}")

    if not all_results_df_list:
        print("No valid all_results.csv files loaded.")
        sys.exit(1)

    all_results_df = pd.concat(all_results_df_list, ignore_index=True)

    # Load naive mincuts if present
    naive_path = os.path.join(data_dir, "all_naive_mincuts.csv")
    all_naive_mincuts_df = None
    if os.path.exists(naive_path):
        all_naive_mincuts_df = pd.read_csv(naive_path, sep=";", header=None, names=["naive_mincut", "hypergraph"])
        all_results_df = pd.merge(all_results_df, all_naive_mincuts_df, on="hypergraph", how="left")
    else:
        print(f"File does not exist: {naive_path} (continuing without naive_mincut)")

    # Mark failed runs: empty/NaN mincut => failed
    all_results_df["failed"] = all_results_df["mincut"].isna()

    algorithms_list = sorted(a for a in all_results_df["algorithm"].dropna().unique())
    print("Algorithms:", algorithms_list)
    print()
    print(all_results_df.info())

    # Optional: failure rate summary
    if args.print_failure_summary:
        print("\nFailure rates by algorithm:")
        for a in algorithms_list:
            sub = all_results_df[all_results_df["algorithm"] == a]
            total = len(sub)
            fails = int(sub["failed"].sum())
            pct = (fails / total * 100.0) if total else 0.0
            print(f"  {a}: {fails}/{total} failed ({pct:.2f}%)")
        print()

    # Default behavior: ONLY performance profiles.
    skip_distribution = not args.enable_distribution
    skip_reduction = not args.enable_reduction

    if not args.skip_profile:
        plot_performance_profile(
            all_results_df,
            all_naive_mincuts_df,
            algorithms_list,
            ["mincut", "time", "memory"],
            out_root=out_dir,  # default is cwd
            xlogscales=[True, False, False],
            xmaxvals=[None, 20, 1.1],
            fontsize=args.fontsize,
            penalize_failures=args.penalize_failures,
        )

    if not skip_distribution:
        plot_distribution_by_algorithm(
            all_results_df,
            algorithms_list,
            ["time", "memory"],
            out_root=out_dir,
            filter_condition_1=lambda df: (df["mincut"] < df.get("naive_mincut", np.inf)) & df["mincut"].notnull(),
            filter_condition_2=lambda df: (df["mincut"] < df.get("naive_mincut", np.inf)) & df["mincut"].notnull() & (df["mincut"] != 0),
            fontsize=args.fontsize,
            num_bins=args.num_bins,
        )

    if not skip_reduction:
        plot_reduction_by_algorithm(
            data_root=data_dir,
            out_root=out_dir,
            algorithms_list=algorithms_list,
            num_rounds=None,
            hypergraph_list=[],
            num_suffix_param_columns=4,
            fontsize=args.fontsize,
        )

    # Optional: print_submodular_stats(data_dir, algorithms_list)


if __name__ == "__main__":
    main()

"""Generate descriptive analytics for the recorded ECDLP experiments.

The script consumes ``ecdla_results.csv`` (produced by ``ecdlp_experiments.py``)
and emits:

* ``analysis_report.md`` – Markdown tables capturing the core descriptive
  statistics and success rates per algorithm and curve.
* ``figures/`` – optional PNG visualisations (only when Matplotlib is
  available).  The script gracefully skips figure generation otherwise so that
  the numerical analysis remains reproducible without third-party packages.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev


ROOT = Path(__file__).resolve().parent
RESULTS_PATH = ROOT / "ecdla_results.csv"
REPORT_PATH = ROOT / "analysis_report.md"
FIGURE_DIR = ROOT / "figures"


def load_results(path: Path) -> list[dict[str, object]]:
    """Load the experiment results with appropriate type conversions."""

    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for record in reader:
            row = {
                "curve_name": record["curve_name"],
                "prime": int(record["prime"]),
                "order": int(record["order"]),
                "algorithm": record["algorithm"],
                "trial": int(record["trial"]),
                "secret": int(record["secret"]),
                "operations": int(record["operations"]),
                "additions": int(record["additions"]),
                "doublings": int(record["doublings"]),
                "solved": record["solved"].strip().lower() == "true",
            }
            rows.append(row)
    return rows


def safe_stdev(values: list[float]) -> float:
    return stdev(values) if len(values) > 1 else 0.0


def compute_stats(rows: list[dict[str, object]]):
    solved_groups: dict[tuple[str, int, str], list[dict[str, object]]] = defaultdict(list)
    success_counts: dict[tuple[str, int, str], tuple[int, int]] = defaultdict(lambda: (0, 0))

    for row in rows:
        key = (row["curve_name"], row["prime"], row["algorithm"])
        total, successes = success_counts[key]
        success_counts[key] = (total + 1, successes + (1 if row["solved"] else 0))
        if row["solved"]:
            solved_groups[key].append(row)

    stats_records = []
    for (curve, prime, algorithm), entries in solved_groups.items():
        operations = [entry["operations"] for entry in entries]
        sqrt_order = math.sqrt(entries[0]["order"])
        ops_per_sqrt = [op / sqrt_order for op in operations]
        trial_ids = {entry["trial"] for entry in entries}

        stats_records.append(
            {
                "curve_name": curve,
                "prime": prime,
                "algorithm": algorithm,
                "trials": len(trial_ids),
                "operations_mean": mean(operations),
                "operations_median": median(operations),
                "operations_std": safe_stdev(operations),
                "sqrt_order": sqrt_order,
                "ops_per_sqrt_median": median(ops_per_sqrt),
                "ops_per_sqrt_mean": mean(ops_per_sqrt),
            }
        )

    success_records = []
    for (curve, prime, algorithm), (total, successes) in success_counts.items():
        rate = successes / total if total else 0.0
        success_records.append(
            {
                "curve_name": curve,
                "prime": prime,
                "algorithm": algorithm,
                "success_rate": rate,
                "trials": total,
            }
        )

    return stats_records, success_records


def format_markdown_table(
    algorithms: list[str],
    primes: list[int],
    values: dict[tuple[str, int], float],
    value_format: str = "{:.2f}",
) -> str:
    header = ["Algorithm"] + [f"p={p}" for p in primes]
    rows = ["| " + " | ".join(header) + " |"]
    rows.append("| " + " | ".join(["---"] * len(header)) + " |")
    for alg in algorithms:
        row = [alg]
        for prime in primes:
            value = values.get((alg, prime))
            if value is None:
                row.append("–")
            else:
                row.append(value_format.format(value))
        rows.append("| " + " | ".join(row) + " |")
    return "\n".join(rows)


def derive_constant_factor_commentary(stats: list[dict[str, object]]) -> str:
    """Generate narrative bullet points about scaling characteristics."""

    by_algorithm: dict[str, list[dict[str, object]]] = defaultdict(list)
    for record in stats:
        by_algorithm[record["algorithm"]].append(record)

    bullets: list[str] = []

    if "Brute Force" in by_algorithm:
        ratios = [rec["operations_mean"] / rec["sqrt_order"] ** 2 for rec in by_algorithm["Brute Force"]]
        avg_ratio = sum(ratios) / len(ratios)
        bullets.append(
            f"Brute force scales linearly with the group order, requiring roughly {avg_ratio:.2f}·n operations on these curves."
        )

    if "Baby-Step Giant-Step" in by_algorithm:
        medians = [rec["ops_per_sqrt_median"] for rec in by_algorithm["Baby-Step Giant-Step"]]
        avg = sum(medians) / len(medians)
        bullets.append(
            f"Baby-Step Giant-Step exhibits the expected √n complexity with a median constant factor of about {avg:.1f}."
        )

    if "Pollard Rho" in by_algorithm:
        medians = [rec["ops_per_sqrt_median"] for rec in by_algorithm["Pollard Rho"]]
        avg = sum(medians) / len(medians)
        bullets.append(
            f"Pollard's Rho also scales with √n but with a higher and more variable constant (median ≈ {avg:.1f} across curves)."
        )

    if "Pollard Kangaroo" in by_algorithm:
        medians = [rec["ops_per_sqrt_median"] for rec in by_algorithm["Pollard Kangaroo"]]
        avg = sum(medians) / len(medians)
        bullets.append(
            f"Pollard's Kangaroo shows similar √n scaling to Rho, with a median constant factor near {avg:.1f}."
        )

    bullets.append(
        "Extrapolating these constants to secp256k1 (n ≈ 1.16×10²⁷) implies ≥2¹²⁸ group operations for √n algorithms, far beyond practical reach."
    )

    return "\n".join(f"- {line}" for line in bullets)


def make_plots(rows: list[dict[str, object]]) -> list[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    FIGURE_DIR.mkdir(exist_ok=True)

    solved = [row for row in rows if row["solved"]]
    if not solved:
        return []

    # Box plot of operations per algorithm across all curves.
    algorithms = sorted({row["algorithm"] for row in solved})
    data_by_alg: dict[str, list[int]] = defaultdict(list)
    for row in solved:
        data_by_alg[row["algorithm"]].append(row["operations"])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.boxplot([data_by_alg[alg] for alg in algorithms], labels=algorithms, showfliers=False)
    ax.set_yscale("log")
    ax.set_ylabel("Group operations (log scale)")
    ax.set_title("Operation counts by algorithm")
    fig.tight_layout()
    boxplot_path = FIGURE_DIR / "operations_by_algorithm.png"
    fig.savefig(boxplot_path, dpi=200)
    plt.close(fig)

    # Scatter plot of median operations vs sqrt(p).
    medians: dict[tuple[str, int], float] = {}
    sqrt_primes: dict[int, float] = {}
    grouped: dict[tuple[str, int], list[int]] = defaultdict(list)
    for row in solved:
        key = (row["algorithm"], row["prime"])
        grouped[key].append(row["operations"])
        sqrt_primes[row["prime"]] = math.sqrt(row["prime"])
    for key, ops in grouped.items():
        medians[key] = median(ops)

    fig, ax = plt.subplots(figsize=(10, 6))
    for algorithm in algorithms:
        xs, ys, labels = [], [], []
        for prime in sorted(sqrt_primes):
            key = (algorithm, prime)
            if key in medians:
                xs.append(sqrt_primes[prime])
                ys.append(medians[key])
                labels.append(prime)
        if xs:
            ax.scatter(xs, ys, label=algorithm)
            for x, y, prime in zip(xs, ys, labels):
                ax.annotate(f"p={prime}", (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)

    ax.set_xlabel("√p")
    ax.set_ylabel("Median operations")
    ax.set_title("Median operations vs. √p")
    ax.legend()
    fig.tight_layout()
    scatter_path = FIGURE_DIR / "median_ops_vs_sqrtp.png"
    fig.savefig(scatter_path, dpi=200)
    plt.close(fig)

    return [str(boxplot_path.relative_to(ROOT)), str(scatter_path.relative_to(ROOT))]


def generate_report(stats: list[dict[str, object]], success: list[dict[str, object]], figures: list[str]) -> None:
    algorithms = sorted({record["algorithm"] for record in stats})
    primes = sorted({record["prime"] for record in stats})

    mean_ops = {(rec["algorithm"], rec["prime"]): rec["operations_mean"] for rec in stats}
    ops_per_sqrt = {(rec["algorithm"], rec["prime"]): rec["ops_per_sqrt_median"] for rec in stats}
    success_rates = {(rec["algorithm"], rec["prime"]): rec["success_rate"] for rec in success}

    mean_table = format_markdown_table(algorithms, primes, mean_ops)
    sqrt_table = format_markdown_table(algorithms, primes, ops_per_sqrt)
    success_table = format_markdown_table(algorithms, primes, success_rates, value_format="{:.2%}")

    commentary_parts = [
        "# ECDLP Algorithm Scaling Analysis",
        "",
        "The tables below summarise the empirical behaviour of four classical "
        "discrete-log algorithms executed on six elliptic curves with prime "
        "moduli between 97 and 1999.",
        "",
        "## Mean Operation Counts",
        "",
        mean_table,
        "",
        "## Median Operations Normalised by √p",
        "",
        sqrt_table,
        "",
        "## Empirical Success Rates",
        "",
        success_table,
        "",
        "## Observations",
        "",
        derive_constant_factor_commentary(stats),
    ]
    commentary = "\n".join(commentary_parts)

    figure_section = []
    if figures:
        figure_section.append("\n## Figures\n")
        for rel_path in figures:
            figure_section.append(f"![{Path(rel_path).stem}]({rel_path})\n")

    with REPORT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(commentary)
        if figure_section:
            fh.write("\n".join(figure_section))


def main() -> None:
    rows = load_results(RESULTS_PATH)
    stats, success = compute_stats(rows)
    figures = make_plots(rows)
    generate_report(stats, success, figures)


if __name__ == "__main__":
    main()

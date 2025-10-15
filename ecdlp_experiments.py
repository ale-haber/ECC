import math
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

try:  # Optional plotting dependencies
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:  # pragma: no cover - graceful degradation when plotting libs missing
    HAS_PLOTTING = False
    plt = None  # type: ignore
    sns = None  # type: ignore

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:  # pragma: no cover - optional analytics dependency
    HAS_PANDAS = False
    pd = None  # type: ignore

import csv
import statistics


Point = Optional[Tuple[int, int]]


@dataclass
class OperationCounter:
    additions: int = 0
    doublings: int = 0

    def record_add(self, is_double: bool = False) -> None:
        if is_double:
            self.doublings += 1
        else:
            self.additions += 1

    @property
    def total(self) -> int:
        return self.additions + self.doublings

    def reset(self) -> None:
        self.additions = 0
        self.doublings = 0


class EllipticCurve:
    def __init__(self, p: int, a: int, b: int):
        self.p = p
        self.a = a
        self.b = b
        self.infinity: Point = None

    def is_on_curve(self, P: Point) -> bool:
        if P is None:
            return True
        x, y = P
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

    def negate(self, P: Point) -> Point:
        if P is None:
            return None
        x, y = P
        return (x, (-y) % self.p)

    def add(self, P: Point, Q: Point, counter: Optional[OperationCounter] = None) -> Point:
        if P is None:
            return Q
        if Q is None:
            return P
        if P == self.negate(Q):
            return None

        if P == Q:
            return self._double(P, counter)

        px, py = P
        qx, qy = Q
        if px == qx:
            return None

        m = ((qy - py) * pow(qx - px, -1, self.p)) % self.p
        rx = (m * m - px - qx) % self.p
        ry = (m * (px - rx) - py) % self.p
        if counter is not None:
            counter.record_add(is_double=False)
        return (rx, ry)

    def _double(self, P: Point, counter: Optional[OperationCounter] = None) -> Point:
        if P is None:
            return None
        px, py = P
        if py == 0:
            return None
        m = ((3 * px * px + self.a) * pow(2 * py, -1, self.p)) % self.p
        rx = (m * m - 2 * px) % self.p
        ry = (m * (px - rx) - py) % self.p
        if counter is not None:
            counter.record_add(is_double=True)
        return (rx, ry)

    def scalar_mult(self, k: int, P: Point, counter: Optional[OperationCounter] = None) -> Point:
        if k == 0 or P is None:
            return None
        result = None
        addend = P
        kk = k
        while kk > 0:
            if kk & 1:
                result = self.add(result, addend, counter)
            addend = self.add(addend, addend, counter)
            kk >>= 1
        return result

    def enumerate_points(self) -> List[Point]:
        points: List[Point] = [None]
        for x in range(self.p):
            rhs = (x ** 3 + self.a * x + self.b) % self.p
            y_sq = rhs % self.p
            ys = self._sqrt_mod(y_sq)
            for y in ys:
                points.append((x, y))
        return points

    def _sqrt_mod(self, a: int) -> List[int]:
        # Tonelli-Shanks for prime p (p odd)
        if a == 0:
            return [0]
        if self.p % 4 == 3:
            r = pow(a, (self.p + 1) // 4, self.p)
            if (r * r) % self.p == a % self.p:
                return [r, (-r) % self.p]
            return []
        # General Tonelli-Shanks
        q = self.p - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1
        if s == 1:
            r = pow(a, (self.p + 1) // 4, self.p)
            return [r, (-r) % self.p]
        for z in range(2, self.p):
            if pow(z, (self.p - 1) // 2, self.p) == self.p - 1:
                c = pow(z, q, self.p)
                break
        else:
            return []
        r = pow(a, (q + 1) // 2, self.p)
        t = pow(a, q, self.p)
        m = s
        while t != 1:
            temp = t
            i = 0
            while temp != 1:
                temp = pow(temp, 2, self.p)
                i += 1
                if i == m:
                    return []
            b = pow(c, 1 << (m - i - 1), self.p)
            r = (r * b) % self.p
            t = (t * b * b) % self.p
            c = (b * b) % self.p
            m = i
        return [r, (-r) % self.p]


@dataclass
class CurveInfo:
    name: str
    curve: EllipticCurve
    generator: Point
    order: int


@dataclass
class ExperimentResult:
    curve_name: str
    prime: int
    order: int
    algorithm: str
    trial: int
    secret: int
    operations: int
    additions: int
    doublings: int
    solved: bool


# ----- Utility helpers -----

def point_order(curve: EllipticCurve, P: Point) -> int:
    counter = OperationCounter()
    current = P
    for k in range(1, curve.p * 2):
        if current is None:
            return k
        current = curve.add(current, P, counter)
    raise ValueError("Order computation failed")


def find_generator(curve: EllipticCurve) -> Tuple[Point, int]:
    points = curve.enumerate_points()
    non_inf = [P for P in points if P is not None]
    for P in non_inf:
        ord_p = point_order(curve, P)
        if ord_p > 20:  # heuristic threshold
            return P, ord_p
    raise ValueError("No suitable generator found")


# ----- Algorithm implementations -----

def brute_force(curve: EllipticCurve, G: Point, Q: Point, order: int) -> Tuple[int, OperationCounter]:
    counter = OperationCounter()
    current = None
    for k in range(order):
        if current == Q:
            return k, counter
        current = curve.add(current, G, counter)
    raise ValueError("Log not found by brute force")


def baby_step_giant_step(curve: EllipticCurve, G: Point, Q: Point, order: int) -> Tuple[int, OperationCounter]:
    m = math.isqrt(order) + 1
    counter = OperationCounter()
    baby_steps: Dict[Point, int] = {}
    current = None
    for j in range(m):
        baby_steps[current] = j
        current = curve.add(current, G, counter)
    mG = curve.scalar_mult(m, G, counter)
    gamma = Q
    neg_mG = curve.negate(mG)
    for i in range(m + 1):
        if gamma in baby_steps:
            j = baby_steps[gamma]
            return (i * m + j) % order, counter
        gamma = curve.add(gamma, neg_mG, counter)
    raise ValueError("Log not found by BSGS")


def pollards_rho(curve: EllipticCurve, G: Point, Q: Point, order: int, seed: Optional[int] = None) -> Tuple[int, OperationCounter]:
    rng = random.Random(seed)
    counter = OperationCounter()

    def partition(P: Point) -> int:
        if P is None:
            return 0
        return P[0] % 3

    def step(P: Point, a: int, b: int) -> Tuple[Point, int, int]:
        subset = partition(P)
        if subset == 0:
            P = curve.add(P, G, counter)
            a = (a + 1) % order
        elif subset == 1:
            P = curve.add(P, Q, counter)
            b = (b + 1) % order
        else:
            P = curve.add(P, P, counter)
            a = (2 * a) % order
            b = (2 * b) % order
        return P, a, b

    for attempt in range(5):
        a = rng.randrange(order)
        b = rng.randrange(order)
        X = curve.add(curve.scalar_mult(a, G, counter), curve.scalar_mult(b, Q, counter), counter)
        Y, c, d = X, a, b
        while True:
            X, a, b = step(X, a, b)
            Y, c, d = step(*step(Y, c, d))
            if X == Y:
                if b == d:
                    break
                numerator = (a - c) % order
                denominator = (d - b) % order
                try:
                    inv = pow(denominator, -1, order)
                except ValueError:
                    break
                k = (numerator * inv) % order
                if curve.scalar_mult(k, G) == Q:
                    return k, counter
                else:
                    break
    raise ValueError("Pollard's Rho failed to find the log")


def pollards_kangaroo(curve: EllipticCurve, G: Point, Q: Point, order: int) -> Tuple[int, OperationCounter]:
    counter = OperationCounter()
    N = order - 1
    step_sizes = [1, 2, 3, 5, 8, 13, 21, 34]
    traps: Dict[Point, int] = {}

    precomputed: List[Point] = []
    for s in step_sizes:
        precomputed.append(curve.scalar_mult(s, G, counter))

    def f(P: Point) -> int:
        if P is None:
            return 0
        return P[0] % len(step_sizes)

    tame = curve.scalar_mult(N // 2, G, counter)
    d_t = N // 2
    max_tame_steps = 3 * len(step_sizes) * math.isqrt(order)
    for _ in range(max_tame_steps):
        idx = f(tame)
        tame = curve.add(tame, precomputed[idx], counter)
        d_t = (d_t + step_sizes[idx]) % order
        traps[tame] = d_t

    wild = Q
    d_w = 0
    max_wild_steps = 2 * order
    for _ in range(max_wild_steps):
        idx = f(wild)
        wild = curve.add(wild, precomputed[idx], counter)
        d_w = (d_w + step_sizes[idx]) % order
        if wild in traps:
            k = (traps[wild] - d_w) % order
            if curve.scalar_mult(k, G) == Q:
                return k, counter
    raise ValueError("Pollard's Kangaroo failed to find the log")


# ----- Experiment driver -----

def build_curves() -> List[CurveInfo]:
    primes = [97, 499, 911, 1237, 1613, 1999]
    curves: List[CurveInfo] = []
    for p in primes:
        curve = EllipticCurve(p, a=2, b=3)
        if (4 * curve.a ** 3 + 27 * curve.b ** 2) % p == 0:
            curve = EllipticCurve(p, a=5, b=7)
        G, order = find_generator(curve)
        curves.append(CurveInfo(name=f"p={p}", curve=curve, generator=G, order=order))
    return curves


def run_experiments(trials_per_curve: int = 25, seed: int = 1234) -> List[ExperimentResult]:
    rng = random.Random(seed)
    records: List[ExperimentResult] = []
    curves = build_curves()
    algorithms = {
        "Brute Force": brute_force,
        "Baby-Step Giant-Step": baby_step_giant_step,
        "Pollard Rho": pollards_rho,
        "Pollard Kangaroo": pollards_kangaroo,
    }

    for info in curves:
        curve, G, order = info.curve, info.generator, info.order
        for trial in range(1, trials_per_curve + 1):
            secret = rng.randrange(1, order)
            Q = curve.scalar_mult(secret, G)
            for name, algo in algorithms.items():
                try:
                    k, counter = algo(curve, G, Q, order)
                except ValueError:
                    k = None
                    counter = OperationCounter()
                records.append(
                    ExperimentResult(
                        curve_name=info.name,
                        prime=curve.p,
                        order=order,
                        algorithm=name,
                        trial=trial,
                        secret=secret,
                        operations=counter.total,
                        additions=counter.additions,
                        doublings=counter.doublings,
                        solved=k is not None,
                    )
                )
    return records


def summarize_results(records: List[ExperimentResult]):
    if HAS_PANDAS:
        df = pd.DataFrame([r.__dict__ for r in records])
        summary = (
            df.groupby(["curve_name", "algorithm"])
            .agg(
                mean_ops=("operations", "mean"),
                median_ops=("operations", "median"),
                std_ops=("operations", "std"),
                success_rate=("solved", "mean"),
            )
            .reset_index()
        )
        return summary

    grouped: Dict[Tuple[str, str], List[ExperimentResult]] = {}
    for rec in records:
        grouped.setdefault((rec.curve_name, rec.algorithm), []).append(rec)

    summary_rows = []
    for (curve_name, algorithm), recs in grouped.items():
        ops = [r.operations for r in recs]
        summary_rows.append(
            {
                "curve_name": curve_name,
                "algorithm": algorithm,
                "mean_ops": statistics.mean(ops),
                "median_ops": statistics.median(ops),
                "std_ops": statistics.pstdev(ops) if len(ops) > 1 else 0.0,
                "success_rate": sum(1 for r in recs if r.solved) / len(recs),
            }
        )
    return summary_rows


def plot_results(records: List[ExperimentResult]) -> None:
    if not HAS_PLOTTING or not HAS_PANDAS:
        print("Plotting libraries not available; skipping visualization generation.")
        return

    df = pd.DataFrame([r.__dict__ for r in records])
    sns.set_theme(style="whitegrid")

    plt.figure(figsize=(14, 8))
    sns.boxplot(
        data=df,
        x="curve_name",
        y="operations",
        hue="algorithm",
        showfliers=False,
    )
    plt.title("Distribution of Group Operations by Algorithm and Curve")
    plt.ylabel("Group operations (additions + doublings)")
    plt.xlabel("Curve (prime modulus)")
    plt.tight_layout()
    plt.savefig("fig_boxplot_operations.png", dpi=200)

    plt.figure(figsize=(12, 7))
    sns.lineplot(
        data=df.groupby(["curve_name", "algorithm"]).operations.mean().reset_index(),
        x="curve_name",
        y="operations",
        hue="algorithm",
        marker="o",
    )
    plt.title("Average Operations vs Curve Size")
    plt.ylabel("Average operations")
    plt.xlabel("Curve (prime modulus)")
    plt.tight_layout()
    plt.savefig("fig_mean_operations.png", dpi=200)

    df_with_sqrt = df.copy()
    df_with_sqrt["sqrt_order"] = df_with_sqrt["order"].apply(math.isqrt)

    plt.figure(figsize=(12, 7))
    sns.scatterplot(
        data=df_with_sqrt,
        x="sqrt_order",
        y="operations",
        hue="algorithm",
        alpha=0.6,
    )
    plt.title("Operations vs sqrt(order)")
    plt.xlabel("sqrt(order)")
    plt.ylabel("Operations")
    plt.tight_layout()
    plt.savefig("fig_scatter_sqrt.png", dpi=200)

    plt.figure(figsize=(14, 8))
    sns.violinplot(
        data=df,
        x="algorithm",
        y="operations",
        inner="quartile",
    )
    plt.title("Overall Operation Count Distributions")
    plt.ylabel("Operations")
    plt.xlabel("Algorithm")
    plt.tight_layout()
    plt.savefig("fig_violin_algorithms.png", dpi=200)


def write_results_csv(records: List[ExperimentResult], path: str) -> None:
    fieldnames = [
        "curve_name",
        "prime",
        "order",
        "algorithm",
        "trial",
        "secret",
        "operations",
        "additions",
        "doublings",
        "solved",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow({
                "curve_name": rec.curve_name,
                "prime": rec.prime,
                "order": rec.order,
                "algorithm": rec.algorithm,
                "trial": rec.trial,
                "secret": rec.secret,
                "operations": rec.operations,
                "additions": rec.additions,
                "doublings": rec.doublings,
                "solved": rec.solved,
            })


def write_summary(summary, path: str) -> None:
    fieldnames = ["curve_name", "algorithm", "mean_ops", "median_ops", "std_ops", "success_rate"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary:
            writer.writerow(row if isinstance(row, dict) else row._asdict())


def main() -> None:
    records = run_experiments(trials_per_curve=25, seed=2024)
    write_results_csv(records, "ecdla_results.csv")
    summary = summarize_results(records)
    if HAS_PANDAS:
        summary_df = summary
        summary_df.to_csv("ecdla_summary.csv", index=False)
        print("Summary statistics:\n", summary_df)
    else:
        write_summary(summary, "ecdla_summary.csv")
        print("Summary statistics:")
        for row in summary:
            print(row)
    plot_results(records)
    if HAS_PLOTTING and HAS_PANDAS:
        print("Saved detailed results to ecdla_results.csv and generated visualizations.")
    else:
        print("Saved detailed results to ecdla_results.csv; visualization step skipped due to missing libraries.")


if __name__ == "__main__":
    main()
